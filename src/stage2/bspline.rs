//! Stage 2.7: NURBS fallback — fit B-spline surfaces to freeform regions.
//!
//! After stages 2.1–2.6, faces that fit no analytic surface remain as
//! single-face planar selections. This stage gathers connected components of
//! those faces and fits each one with a `Geom_BSplineSurface`:
//!
//! 1. Parameterize the region onto the unit square with Floater mean-value
//!    coordinates (boundary chord-length-pinned to the square perimeter,
//!    interior solved by Gauss–Seidel).
//! 2. Resample a regular UV grid of 3D points from the parameterized mesh
//!    (piecewise-linear, i.e. chordal, between vertices; hole regions filled
//!    with a Laplacian membrane so the approximator sees smooth data).
//! 3. Fit with `GeomAPI_PointsToBSplineSurface` and validate every region
//!    vertex and face centroid against `--surface-tolerance`.
//!
//! The grid is oriented so the fitted surface's natural normal (d1u × d1v)
//! matches the mesh's outward normals: B-spline faces are always
//! `concave = false`/FORWARD, and downstream orientation logic needs no
//! per-face sign heuristics (DEVELOPMENT_PLAN §2.7 implementation notes).
//!
//! The mesh styles this must handle range from scripted regular grids to
//! CAD-export STRIP meshes (full-length quads with no interior vertices along
//! a ruling, e.g. tessellated extrusions): the parameterization is driven by
//! the mesh boundary, so a 2-ring strip region works the same way as a dense
//! grid — the ruling direction simply has no interior unknowns to solve.

use crate::config::Config;
use crate::stage1::ConnectedMesh;
use crate::viz::{self, VizAction, VizSender};
use opencascade_sys::{
    approx, geom, geom_abs, geom_api, gp, shape_analysis, t_col_std, t_colgp, OwnedPtr,
};
use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Instant;

use super::{
    face_centroid, planar, ConicalHypothesis, PlanarHypothesis, SelectedSurface,
    SphericalHypothesis, ToroidalHypothesis,
};

/// Minimum number of mesh faces for a B-spline fit attempt. A bicubic patch
/// has many DOF; below this the region is noise (and DEVELOPMENT_PLAN §2.7
/// specifies ≥ 4).
pub(crate) const MIN_BSPLINE_FACES: usize = 4;

/// Selections up to this many faces are candidates for absorption into a
/// freeform region they connect to smoothly. Freeform surfaces bait the
/// analytic detectors locally: an inflection produces short runs of
/// accidentally-coplanar facets fitted as tiny planes (bspline_wave),
/// bands of near-constant curvature fit as narrow cylinder ribbons
/// (bspline_saddle: 16 six-face cylinders carving the patch apart), and a
/// Gaussian bump peak fits 9–12-face cones (hill_with_hole). Left
/// selected, they punch holes in the freeform region and leave fake
/// feature edges. Real small features are protected twice: a component
/// with no single-face planar member is never fitted, and absorbed
/// non-single faces must not outnumber the single-face leftovers (see
/// `deduce_bspline_surfaces`).
const SMALL_SELECTION_ABSORB_FACES: usize = 16;

/// Precision passed to ShapeAnalysis_Surface projections (mm). This is the
/// convergence target, not an acceptance gate.
pub(crate) const PROJ_PRECISION: f64 = 1e-7;

/// Above this projection distance (mm), a seeded NextValueOfUV solution is
/// considered lost and ShapeAnalysis_Surface internally re-runs the global
/// ValueOfUV search. Chordal queries legitimately sit up to ~surface
/// tolerance off the fit, so this must sit well above that.
pub(crate) const PROJ_MAX_PRECISION: f64 = 2.0;

/// A seeded projection whose gap exceeds this (mm) is re-verified with a
/// fresh global ValueOfUV and the better solution wins. NextValueOfUV's
/// Newton step can converge to a nearby local minimum on a many-knot
/// B-spline (observed 1.4 mm bogus gaps on points lying exactly on the
/// surface), and its own internal fallback only fires above
/// PROJ_MAX_PRECISION. Genuine query points sit within chordal sag of the
/// fit, so rechecks beyond this threshold are rare in the hot paths.
pub(crate) const PROJ_GAP_RECHECK: f64 = 0.02;

// ---------------------------------------------------------------------------
// Hypothesis
// ---------------------------------------------------------------------------

/// Cached OCCT projector for a fitted B-spline surface.
/// `ShapeAnalysis_Surface` caches internally (hence `&mut` methods), and
/// seeding each projection with the previous solution (`NextValueOfUV`) is
/// the B-spline analog of the analytic charts' branch continuity.
struct Projector {
    sas: OwnedPtr<shape_analysis::Surface>,
    seed: Option<(f64, f64)>,
}

/// A B-spline surface hypothesis fitted to a freeform region of mesh faces.
///
/// Unlike the analytic hypotheses this owns the actual OCCT surface: fitted
/// geometry is not reconstructible from a few parameters, and stage 3 must
/// use *this exact* `Geom_BSplineSurface` object (charts wrap it via
/// ShapeAnalysis_Surface, so a rebuilt copy would silently disagree).
pub struct BSplineHypothesis {
    /// The fitted surface. Stage 3 mints type-erased handles to this same
    /// object via `to_handle_surface()`.
    pub surface: OwnedPtr<geom::HandleGeomBSplineSurface>,
    /// Cached projector for distance/normal queries (interior-mutable
    /// because ShapeAnalysis_Surface caches internally).
    projector: RefCell<Projector>,
    /// Set of mesh face indices covered by this surface.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this surface.
    pub vertices: Vec<usize>,
    /// Maximum distance from any region vertex to the fitted surface.
    pub error_max: f64,
    /// Maximum distance from any face centroid to the fitted surface.
    pub centroid_error_max: f64,
    /// Sum of vertex-to-surface distances.
    pub error_abs_sum: f64,
    /// Fit shape summary for logging.
    pub u_degree: i32,
    pub v_degree: i32,
    pub nb_u_poles: i32,
    pub nb_v_poles: i32,
    /// True for a closed freeform band fitted periodically around the band:
    /// the surface is U-periodic with its U knots rescaled to exactly
    /// [0, 2π] (the period stage 3.4's wrapped-face machinery hard-codes),
    /// and the U=0 iso is the gauge seam at the region's cut path.
    pub u_periodic: bool,
}

impl std::fmt::Debug for BSplineHypothesis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BSplineHypothesis")
            .field("surface", &"<HandleGeomBSplineSurface>")
            .field("faces", &self.faces.len())
            .field("vertices", &self.vertices.len())
            .field("error_max", &self.error_max)
            .field("centroid_error_max", &self.centroid_error_max)
            .field("degree", &(self.u_degree, self.v_degree))
            .field("poles", &(self.nb_u_poles, self.nb_v_poles))
            .field("u_periodic", &self.u_periodic)
            .finish()
    }
}

impl BSplineHypothesis {
    fn new(
        surface: OwnedPtr<geom::HandleGeomBSplineSurface>,
        faces: Vec<usize>,
        vertices: Vec<usize>,
    ) -> Self {
        let bs = surface.get();
        let (u_degree, v_degree) = (bs.u_degree(), bs.v_degree());
        let (nb_u_poles, nb_v_poles) = (bs.nb_u_poles(), bs.nb_v_poles());
        let erased = surface.to_handle_surface();
        let sas = shape_analysis::Surface::new_handlegeomsurface(&erased);
        BSplineHypothesis {
            surface,
            projector: RefCell::new(Projector { sas, seed: None }),
            faces,
            vertices,
            error_max: 0.0,
            centroid_error_max: 0.0,
            error_abs_sum: 0.0,
            u_degree,
            v_degree,
            nb_u_poles,
            nb_v_poles,
            u_periodic: false,
        }
    }

    /// UV of the surface point nearest to `p` (robust OCCT projection,
    /// seeded from the previous query), plus the projection gap.
    /// A seeded solution with a suspicious gap is re-verified with a fresh
    /// global search — NextValueOfUV's Newton step can settle in a nearby
    /// local minimum without triggering its internal fallback.
    fn point_uv_gap(&self, p: &[f64; 3]) -> (f64, f64, f64) {
        let pnt = gp::Pnt::new_real3(p[0], p[1], p[2]);
        let mut proj = self.projector.borrow_mut();
        let mut best = match proj.seed {
            Some((su, sv)) => {
                let prev = gp::Pnt2d::new_real2(su, sv);
                let uv = proj
                    .sas
                    .next_value_of_uv(&prev, &pnt, PROJ_PRECISION, PROJ_MAX_PRECISION);
                (uv.x(), uv.y(), proj.sas.gap())
            }
            None => {
                let uv = proj.sas.value_of_uv(&pnt, PROJ_PRECISION);
                (uv.x(), uv.y(), proj.sas.gap())
            }
        };
        if proj.seed.is_some() && best.2 > PROJ_GAP_RECHECK {
            let uv = proj.sas.value_of_uv(&pnt, PROJ_PRECISION);
            let gap = proj.sas.gap();
            if gap < best.2 {
                best = (uv.x(), uv.y(), gap);
            }
        }
        proj.seed = Some((best.0, best.1));
        best
    }

    /// UV of the surface point nearest to `p`.
    pub fn point_uv(&self, p: &[f64; 3]) -> (f64, f64) {
        let (u, v, _) = self.point_uv_gap(p);
        (u, v)
    }

    /// Distance from `p` to the fitted surface.
    pub fn distance(&self, p: &[f64; 3]) -> f64 {
        self.point_uv_gap(p).2
    }

    /// Distance from `p` to the fitted surface, seeding the projection at a
    /// known parameter-space location (normalized [0,1]² fractions of the
    /// surface's knot range). Used during fit validation, where the
    /// parameterization that produced the fit gives every mesh vertex a
    /// trustworthy starting point — sequential seeding breaks at mesh row
    /// wraps and neither NextValueOfUV nor ValueOfUV recovers reliably on
    /// many-knot surfaces.
    pub fn distance_seeded(&self, p: &[f64; 3], seed_frac: (f64, f64)) -> f64 {
        let (mut u0, mut u1, mut v0, mut v1) = (0.0, 0.0, 0.0, 0.0);
        self.surface.get().bounds(&mut u0, &mut u1, &mut v0, &mut v1);
        let seed = (
            u0 + seed_frac.0.clamp(0.0, 1.0) * (u1 - u0),
            v0 + seed_frac.1.clamp(0.0, 1.0) * (v1 - v0),
        );
        self.projector.borrow_mut().seed = Some(seed);
        self.point_uv_gap(p).2
    }

    /// Nearest point on the fitted surface.
    pub fn project_point(&self, p: &[f64; 3]) -> [f64; 3] {
        let (u, v) = self.point_uv(p);
        let mut proj = self.projector.borrow_mut();
        let q = proj.sas.value_real2(u, v);
        [q.x(), q.y(), q.z()]
    }

    /// Outward surface normal at the point of the surface nearest to `p`.
    /// The fit orients the natural normal (d1u × d1v) outward, so no sign
    /// logic is needed (`concave` is always false for B-spline faces).
    pub fn normal_at(&self, p: &[f64; 3]) -> Option<[f64; 3]> {
        let (u, v) = self.point_uv(p);
        let mut pt = gp::Pnt::new_real3(0.0, 0.0, 0.0);
        let mut d1u = gp::Vec::new();
        let mut d1v = gp::Vec::new();
        self.surface.get().d1(u, v, &mut pt, &mut d1u, &mut d1v);
        let n = [
            d1u.y() * d1v.z() - d1u.z() * d1v.y(),
            d1u.z() * d1v.x() - d1u.x() * d1v.z(),
            d1u.x() * d1v.y() - d1u.y() * d1v.x(),
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len < 1e-15 {
            return None;
        }
        Some([n[0] / len, n[1] / len, n[2] / len])
    }
}

// ---------------------------------------------------------------------------
// Region extraction
// ---------------------------------------------------------------------------

/// Sentinel "outside selection" id for boundary edges created by a seam
/// cut: they were interior mesh edges, so no outside selection exists. A
/// distinct value from usize::MAX (open mesh boundary) so corner detection
/// sees transitions exactly at the cut endpoints.
const CUT_OUTSIDE: usize = usize::MAX - 1;

/// A connected freeform region: faces, triangulated geometry, and ordered
/// boundary loops (region on the left of each directed boundary edge, i.e.
/// counter-clockwise seen from outside the solid). Vertices are indexed
/// LOCALLY: a seam cut duplicates the cut-path vertices, so one global mesh
/// vertex can appear as two local vertices.
struct Region {
    faces: Vec<usize>,
    /// Region triangulation (quads split along the v0–v2 diagonal), local
    /// vertex indices, mesh winding preserved.
    tris: Vec<[usize; 3]>,
    /// Region face index (into `faces`) of each triangle.
    tri_face: Vec<usize>,
    /// Local → global mesh vertex index (seam-cut clones repeat one).
    verts: Vec<usize>,
    /// Local vertex positions (clones repeat their original's).
    positions: Vec<[f64; 3]>,
    /// Global → one local index (originals win over clones), for seeding.
    global_to_local: HashMap<usize, usize>,
    /// Outside selection id per undirected GLOBAL vertex pair of the
    /// original region boundary (global keys are clone-transparent).
    outside: HashMap<(usize, usize), usize>,
    /// Boundary loops as ordered (local vertex, outside_surface_id) pairs:
    /// the directed boundary edge starts at `vertex` and is shared with a
    /// face of selected surface `outside_surface_id` (or CUT_OUTSIDE).
    loops: Vec<Vec<(usize, usize)>>,
}

fn build_region(
    mesh: &ConnectedMesh,
    faces: Vec<usize>,
    face_to_selected: &[usize],
    in_region: &[bool],
) -> Result<Region, String> {
    let mut vert_set = HashSet::new();
    for &fi in &faces {
        let face = &mesh.faces[fi];
        for &v in &face.vertex_indices[..face.vertex_count as usize] {
            vert_set.insert(v);
        }
    }
    let mut verts: Vec<usize> = vert_set.into_iter().collect();
    verts.sort_unstable();
    let global_to_local: HashMap<usize, usize> =
        verts.iter().enumerate().map(|(li, &vi)| (vi, li)).collect();
    let positions: Vec<[f64; 3]> = verts
        .iter()
        .map(|&vi| {
            let v = &mesh.vertices[vi];
            [v.x, v.y, v.z]
        })
        .collect();

    let mut tris = Vec::with_capacity(faces.len() * 2);
    let mut tri_face = Vec::with_capacity(faces.len() * 2);
    let mut outside: HashMap<(usize, usize), usize> = HashMap::new();
    for (rf, &fi) in faces.iter().enumerate() {
        let face = &mesh.faces[fi];
        let vc = face.vertex_count as usize;
        let vi = &face.vertex_indices;
        let l = |g: usize| global_to_local[&g];
        tris.push([l(vi[0]), l(vi[1]), l(vi[2])]);
        tri_face.push(rf);
        if vc == 4 {
            tris.push([l(vi[0]), l(vi[2]), l(vi[3])]);
            tri_face.push(rf);
        }
        for k in 0..vc {
            let ni = face.neighbors[k];
            if ni < 0 || !in_region[ni as usize] {
                let a = vi[k];
                let b = vi[(k + 1) % vc];
                let outside_id = if ni >= 0 {
                    face_to_selected[ni as usize]
                } else {
                    usize::MAX
                };
                outside.insert((a.min(b), a.max(b)), outside_id);
            }
        }
    }

    let mut region = Region {
        faces,
        tris,
        tri_face,
        verts,
        positions,
        global_to_local,
        outside,
        loops: Vec::new(),
    };
    region.loops = chain_boundary_loops(&region)?;
    if region.loops.is_empty() {
        return Err("region has no boundary loop".to_string());
    }
    Ok(region)
}

/// Chain the triangulation's boundary edges (edges with exactly one
/// incident triangle) into ordered loops, region on the left. Outside ids
/// come from the region's global-pair map; edges absent from it are
/// seam-cut edges.
fn chain_boundary_loops(region: &Region) -> Result<Vec<Vec<(usize, usize)>>, String> {
    // Count incidences per undirected local edge; remember the directed
    // form each edge takes in its (unique, for boundary edges) triangle.
    let mut count: HashMap<(usize, usize), (u32, (usize, usize))> = HashMap::new();
    for tri in &region.tris {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            let e = count.entry((a.min(b), a.max(b))).or_insert((0, (a, b)));
            e.0 += 1;
            e.1 = (a, b);
        }
    }
    let mut boundary_next: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    let mut n_boundary_edges = 0usize;
    for (&key, &(c, (a, b))) in &count {
        if c == 1 {
            let ga = region.verts[key.0];
            let gb = region.verts[key.1];
            let outside_id = region
                .outside
                .get(&(ga.min(gb), ga.max(gb)))
                .copied()
                .unwrap_or(CUT_OUTSIDE);
            boundary_next.entry(a).or_default().push((b, outside_id));
            n_boundary_edges += 1;
        }
    }
    // Deterministic chaining: smallest start vertex, smallest next first.
    for outs in boundary_next.values_mut() {
        outs.sort_unstable();
    }
    let mut loops = Vec::new();
    let mut consumed = 0usize;
    while consumed < n_boundary_edges {
        let start = *boundary_next
            .iter()
            .filter(|(_, outs)| !outs.is_empty())
            .map(|(a, _)| a)
            .min()
            .ok_or("boundary chaining lost edges")?;
        let mut lp: Vec<(usize, usize)> = Vec::new();
        let mut cur = start;
        loop {
            let outs = boundary_next
                .get_mut(&cur)
                .filter(|o| !o.is_empty())
                .ok_or_else(|| format!("open boundary chain at local vertex {cur}"))?;
            // A vertex can have two outgoing boundary edges when the region
            // touches itself at a vertex; taking any unconsumed edge still
            // partitions the boundary into valid loops.
            let (next, outside_id) = outs.remove(0);
            lp.push((cur, outside_id));
            consumed += 1;
            cur = next;
            if cur == start {
                break;
            }
            if lp.len() > n_boundary_edges {
                return Err("boundary chaining did not close".to_string());
            }
        }
        if lp.len() >= 3 {
            loops.push(lp);
        }
    }
    Ok(loops)
}

/// Squared bounding-box diagonal of each boundary loop.
fn loop_diags_sq(region: &Region) -> Vec<f64> {
    region
        .loops
        .iter()
        .map(|lp| {
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for &(li, _) in lp {
                let p = region.positions[li];
                for d in 0..3 {
                    lo[d] = lo[d].min(p[d]);
                    hi[d] = hi[d].max(p[d]);
                }
            }
            (0..3).map(|d| (hi[d] - lo[d]).powi(2)).sum::<f64>()
        })
        .collect()
}

/// Index (into `loops`) of the outer boundary loop: the one with the largest
/// 3D bounding-box diagonal. Other loops are holes and float freely in the
/// parameterization.
fn outer_loop_index(region: &Region) -> usize {
    loop_diags_sq(region)
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// True when the region is a closed band (annulus) rather than a disk with
/// holes: it has a second boundary loop of extent comparable to the largest
/// (≥ half its bounding-box diagonal). A small hole floats freely in the
/// disk parameterization; a comparable second loop means the region wraps
/// around (e.g. the wall of a closed prism) and must be seam-cut first.
fn is_closed_band(region: &Region) -> bool {
    let mut d = loop_diags_sq(region);
    if d.len() < 2 {
        return false;
    }
    d.sort_by(|a, b| b.partial_cmp(a).unwrap());
    d[1] >= 0.25 * d[0]
}

/// Cut a closed-band region into a disk along a shortest vertex path
/// between its two largest boundary loops. The path vertices are
/// duplicated: triangles on the path's right side (region winding) switch
/// to the clones, and the boundary is re-chained — the two loops merge
/// into one, with the two cut copies flagged CUT_OUTSIDE. Clones keep
/// their original's global index and position, so the sampled seam columns
/// at the parameter extremes carry identical 3D data and the fitted
/// surface closes exactly.
fn cut_region_seam(
    mut region: Region,
    mesh: &ConnectedMesh,
    verbosity: u8,
) -> Result<Region, String> {
    // The two largest loops.
    let diags = loop_diags_sq(&region);
    let mut order: Vec<usize> = (0..diags.len()).collect();
    order.sort_by(|&a, &b| diags[b].partial_cmp(&diags[a]).unwrap());
    let (la, lb) = (order[0], order[1]);
    let on_a: HashSet<usize> = region.loops[la].iter().map(|&(li, _)| li).collect();
    let on_b: HashSet<usize> = region.loops[lb].iter().map(|&(li, _)| li).collect();

    // Vertex adjacency over the region faces' PERIMETER edges — never quad
    // diagonals: a path along a diagonal would split that quad's two
    // triangles onto opposite sides of the cut, and the face partition
    // downstream must follow the cut exactly.
    let nv = region.verts.len();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for &fi in &region.faces {
        let face = &mesh.faces[fi];
        let vc = face.vertex_count as usize;
        for k in 0..vc {
            let a = region.global_to_local[&face.vertex_indices[k]];
            let b = region.global_to_local[&face.vertex_indices[(k + 1) % vc]];
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    for l in adj.iter_mut() {
        l.sort_unstable();
        l.dedup();
    }

    // Multi-source BFS (hop count) from loop A to loop B; deterministic
    // frontier order.
    let mut parent = vec![usize::MAX; nv];
    let mut seen = vec![false; nv];
    let mut frontier: Vec<usize> = {
        let mut f: Vec<usize> = on_a.iter().copied().collect();
        f.sort_unstable();
        f
    };
    for &s in &frontier {
        seen[s] = true;
    }
    let mut hit = None;
    'bfs: while !frontier.is_empty() {
        let mut next = Vec::new();
        for &v in &frontier {
            for &w in &adj[v] {
                if !seen[w] {
                    seen[w] = true;
                    parent[w] = v;
                    if on_b.contains(&w) {
                        hit = Some(w);
                        break 'bfs;
                    }
                    next.push(w);
                }
            }
        }
        frontier = next;
    }
    let Some(end) = hit else {
        return Err("closed band: no path between boundary loops".to_string());
    };
    let mut path = vec![end];
    while parent[*path.last().unwrap()] != usize::MAX {
        path.push(parent[*path.last().unwrap()]);
    }
    path.reverse(); // loop A → loop B
    if verbosity >= 2 {
        eprintln!(
            "  Stage 2.7: closed band — seam cut along {} vertices from \
             [{:.3},{:.3},{:.3}] to [{:.3},{:.3},{:.3}]",
            path.len(),
            region.positions[path[0]][0],
            region.positions[path[0]][1],
            region.positions[path[0]][2],
            region.positions[path[path.len() - 1]][0],
            region.positions[path[path.len() - 1]][1],
            region.positions[path[path.len() - 1]][2],
        );
    }

    // Classify the triangles flanking each path edge: with consistent mesh
    // winding, the triangle whose winding contains pi→pi+1 is on the LEFT
    // of the directed path; the one containing pi+1→pi is on the RIGHT.
    // Right-side fans switch to the clones.
    let path_edge: HashSet<(usize, usize)> = path
        .windows(2)
        .map(|w| (w[0].min(w[1]), w[0].max(w[1])))
        .collect();
    let mut right_seed: HashMap<usize, Vec<usize>> = HashMap::new(); // path vtx → tri ids
    for (ti, tri) in region.tris.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            if path_edge.contains(&(a.min(b), a.max(b))) {
                // Directed path order: does (a, b) match path direction?
                let pa = path.iter().position(|&p| p == a).unwrap();
                let along = pa + 1 < path.len() && path[pa + 1] == b;
                if !along {
                    // Contains the reversed edge → right side.
                    right_seed.entry(a).or_default().push(ti);
                    right_seed.entry(b).or_default().push(ti);
                }
            }
        }
    }

    // Per path vertex: flood its star from the right-side seeds without
    // crossing path edges; those triangles get the clone.
    let mut vert_tris: Vec<Vec<usize>> = vec![Vec::new(); nv];
    for (ti, tri) in region.tris.iter().enumerate() {
        for &v in tri {
            vert_tris[v].push(ti);
        }
    }
    for (pi, &pv) in path.iter().enumerate() {
        let clone = region.verts.len();
        region.verts.push(region.verts[pv]);
        region.positions.push(region.positions[pv]);
        vert_tris.push(Vec::new());

        let star: Vec<usize> = vert_tris[pv].clone();
        let seeds: Vec<usize> = right_seed.get(&pv).cloned().unwrap_or_default();
        if seeds.is_empty() {
            return Err(format!("seam cut: path vertex {pi} has no right-side triangle"));
        }
        // Star adjacency: two star triangles are connected iff they share
        // an edge incident to pv that is not a path edge.
        let mut in_right: HashSet<usize> = seeds.iter().copied().collect();
        let mut stack = seeds;
        while let Some(ti) = stack.pop() {
            let tri = region.tris[ti];
            for k in 0..3 {
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                if a != pv && b != pv {
                    continue; // not incident to pv
                }
                let key = (a.min(b), a.max(b));
                if path_edge.contains(&key) {
                    continue; // never cross the cut
                }
                for &tj in &vert_tris[pv] {
                    if tj != ti && !in_right.contains(&tj) {
                        let tjt = region.tris[tj];
                        let has_edge = (0..3).any(|m| {
                            let (c, d) = (tjt[m], tjt[(m + 1) % 3]);
                            (c.min(d), c.max(d)) == key
                        });
                        if has_edge {
                            in_right.insert(tj);
                            stack.push(tj);
                        }
                    }
                }
            }
        }
        if in_right.len() >= star.len() {
            return Err(format!(
                "seam cut: path vertex {pi} star did not split ({} of {} triangles on the right)",
                in_right.len(),
                star.len(),
            ));
        }
        for &ti in &in_right {
            for v in region.tris[ti].iter_mut() {
                if *v == pv {
                    *v = clone;
                }
            }
            vert_tris[clone].push(ti);
        }
        vert_tris[pv].retain(|ti| !in_right.contains(ti));
    }

    region.loops = chain_boundary_loops(&region)?;
    if region.loops.len() != diags.len() - 1 {
        return Err(format!(
            "seam cut did not merge the loops ({} before, {} after)",
            diags.len(),
            region.loops.len(),
        ));
    }
    Ok(region)
}

// ---------------------------------------------------------------------------
// Boundary corners and square mapping
// ---------------------------------------------------------------------------

/// Choose 4 corner positions (indices into the loop) to pin to the UV square
/// corners. Primary signal: vertices where the adjacent *outside* surface
/// changes — region corners in the B-rep sense. When the boundary offers a
/// different number of such transitions (e.g. a tangent bump inside one
/// plane has zero), fall back to arc-length quartiles so the mapping stays
/// well-formed; the corners are then a gauge choice, not a feature claim.
fn choose_corners(
    lp: &[(usize, usize)],
    positions: &[[f64; 3]],
) -> Vec<usize> {
    let n = lp.len();
    // Cumulative chord length at each loop position.
    let mut cum = vec![0.0_f64; n + 1];
    for i in 0..n {
        let a = positions[lp[i].0];
        let b = positions[lp[(i + 1) % n].0];
        let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
        cum[i + 1] = cum[i] + d;
    }
    let total = cum[n];

    // Outside-surface transitions: vertex i sits between boundary edge i-1
    // and edge i; it is a corner candidate when their outside ids differ.
    let mut candidates: Vec<usize> = (0..n)
        .filter(|&i| lp[(i + n - 1) % n].1 != lp[i].1)
        .collect();

    if candidates.len() == 4 {
        return candidates;
    }
    if candidates.is_empty() {
        candidates.push(0);
    }
    // Anchor at the first candidate, then pick the loop positions closest to
    // the quarter arc-length marks (from remaining candidates when there are
    // too many, from any loop vertex when there are too few).
    let anchor = candidates[0];
    let mut chosen = vec![anchor];
    let pool: Vec<usize> = if candidates.len() > 4 {
        candidates.clone()
    } else {
        (0..n).collect()
    };
    for k in 1..4 {
        let target = (cum[anchor] + total * k as f64 / 4.0) % total;
        let best = pool
            .iter()
            .copied()
            .filter(|i| !chosen.contains(i))
            .min_by(|&a, &b| {
                let da = arc_dist(cum[a], target, total);
                let db = arc_dist(cum[b], target, total);
                da.partial_cmp(&db).unwrap()
            });
        if let Some(i) = best {
            chosen.push(i);
        }
    }
    chosen.sort_unstable();
    chosen.dedup();
    chosen
}

fn arc_dist(a: f64, b: f64, total: f64) -> f64 {
    let d = (a - b).rem_euclid(total);
    d.min(total - d)
}

/// Pin the outer loop to the unit-square perimeter: corners at the square
/// corners, each side spaced by cumulative chord length. The loop is CCW
/// around the region, and the square is traced CCW in UV, so the flattening
/// preserves orientation (the fitted normal d1u × d1v points outward).
fn map_outer_loop_to_square(
    lp: &[(usize, usize)],
    corners: &[usize],
    positions: &[[f64; 3]],
    uv: &mut [(f64, f64)],
) {
    let n = lp.len();
    for s in 0..4 {
        let c0 = corners[s];
        let c1 = corners[(s + 1) % 4];
        // Loop positions c0..=c1 (cyclic).
        let len = (c1 + n - c0) % n; // edges on this side
        let mut side_cum = vec![0.0_f64; len + 1];
        for e in 0..len {
            let i = (c0 + e) % n;
            let a = positions[lp[i].0];
            let b = positions[lp[(i + 1) % n].0];
            let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2))
                .sqrt();
            side_cum[e + 1] = side_cum[e] + d;
        }
        let side_total = side_cum[len].max(1e-30);
        for (e, &cum_e) in side_cum[..len].iter().enumerate() {
            let i = (c0 + e) % n;
            let t = cum_e / side_total;
            uv[lp[i].0] = match s {
                0 => (t, 0.0),
                1 => (1.0, t),
                2 => (1.0 - t, 1.0),
                _ => (0.0, 1.0 - t),
            };
        }
    }
}

// ---------------------------------------------------------------------------
// Floater mean-value parameterization
// ---------------------------------------------------------------------------

/// Solve for interior (and hole-boundary) vertex UVs with mean-value
/// weights: each unknown vertex is the convex combination of its neighbors,
/// w_ij = Σ tan(γ/2)/|p_i − p_j| over the triangles flanking edge (i,j).
/// Plain Gauss–Seidel; the systems here are at most a few thousand unknowns.
fn solve_interior_uv(
    region: &Region,
    pinned: &[bool],
    uv: &mut [(f64, f64)],
    verbosity: u8,
) -> Result<(), String> {
    let nv = region.verts.len();
    // Accumulate weights in BTreeMaps for deterministic iteration order.
    let mut weights: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); nv];
    for tri in &region.tris {
        for k in 0..3 {
            let i = tri[k];
            if pinned[i] {
                continue;
            }
            let pi = region.positions[tri[k]];
            let pa = region.positions[tri[(k + 1) % 3]];
            let pb = region.positions[tri[(k + 2) % 3]];
            let ea = [pa[0] - pi[0], pa[1] - pi[1], pa[2] - pi[2]];
            let eb = [pb[0] - pi[0], pb[1] - pi[1], pb[2] - pi[2]];
            let la = (ea[0] * ea[0] + ea[1] * ea[1] + ea[2] * ea[2]).sqrt();
            let lb = (eb[0] * eb[0] + eb[1] * eb[1] + eb[2] * eb[2]).sqrt();
            if la < 1e-30 || lb < 1e-30 {
                return Err("degenerate triangle in region".to_string());
            }
            let cosg = ((ea[0] * eb[0] + ea[1] * eb[1] + ea[2] * eb[2]) / (la * lb))
                .clamp(-1.0, 1.0);
            let sing = (1.0 - cosg * cosg).sqrt();
            // tan(γ/2), stable form; clamp near-degenerate angles.
            let tan_half = if 1.0 + cosg < 1e-12 { 1e6 } else { sing / (1.0 + cosg) };
            let a = tri[(k + 1) % 3];
            let b = tri[(k + 2) % 3];
            *weights[i].entry(a).or_insert(0.0) += tan_half / la;
            *weights[i].entry(b).or_insert(0.0) += tan_half / lb;
        }
    }

    let unknowns: Vec<usize> = (0..nv).filter(|&i| !pinned[i]).collect();
    if unknowns.is_empty() {
        return Ok(());
    }
    for &i in &unknowns {
        if weights[i].is_empty() {
            return Err("isolated interior vertex in region".to_string());
        }
        // Start unknowns at the square center.
        uv[i] = (0.5, 0.5);
    }
    let mut max_delta = f64::INFINITY;
    let mut iters = 0usize;
    for _iter in 0..30_000 {
        iters += 1;
        max_delta = 0.0;
        for &i in &unknowns {
            let mut su = 0.0;
            let mut sv = 0.0;
            let mut sw = 0.0;
            for (&j, &w) in &weights[i] {
                su += w * uv[j].0;
                sv += w * uv[j].1;
                sw += w;
            }
            let nu = su / sw;
            let nv2 = sv / sw;
            max_delta = max_delta.max((nu - uv[i].0).abs()).max((nv2 - uv[i].1).abs());
            uv[i] = (nu, nv2);
        }
        if max_delta < 1e-10 {
            break;
        }
    }
    if verbosity >= 3 {
        eprintln!(
            "  Stage 2.7: mean-value solve: {} unknowns, {iters} iterations, final delta {max_delta:.2e}",
            unknowns.len(),
        );
    }
    if !max_delta.is_finite() {
        return Err("parameterization solve diverged".to_string());
    }
    // An unconverged map silently corrupts every sample; Floater
    // bijectivity holds only for the exact solution.
    if max_delta > 1e-8 {
        return Err(format!(
            "parameterization solve did not converge (final delta {max_delta:.2e})"
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Grid resampling
// ---------------------------------------------------------------------------

/// Piecewise-linear sampler over the flattened region triangulation, with a
/// uniform UV bin grid for point location.
struct UvSampler<'a> {
    region: &'a Region,
    uv: &'a [(f64, f64)],
    bins: Vec<Vec<usize>>,
    nbins: usize,
}

impl<'a> UvSampler<'a> {
    fn new(region: &'a Region, uv: &'a [(f64, f64)]) -> Self {
        let nbins = 64usize;
        let mut bins = vec![Vec::new(); nbins * nbins];
        for (ti, tri) in region.tris.iter().enumerate() {
            let uvs = tri.map(|li| uv[li]);
            let (mut lo_u, mut hi_u) = (f64::INFINITY, f64::NEG_INFINITY);
            let (mut lo_v, mut hi_v) = (f64::INFINITY, f64::NEG_INFINITY);
            for (u, v) in uvs {
                lo_u = lo_u.min(u);
                hi_u = hi_u.max(u);
                lo_v = lo_v.min(v);
                hi_v = hi_v.max(v);
            }
            let bin = |x: f64| ((x * nbins as f64) as isize).clamp(0, nbins as isize - 1) as usize;
            for bu in bin(lo_u - 1e-9)..=bin(hi_u + 1e-9) {
                for bv in bin(lo_v - 1e-9)..=bin(hi_v + 1e-9) {
                    bins[bu * nbins + bv].push(ti);
                }
            }
        }
        UvSampler { region, uv, bins, nbins }
    }

    /// 3D point at parameter (u, v), or None when (u, v) is not covered by
    /// the flattened region (a genuine hole). Numeric misses — the Floater
    /// map squeezes triangles near the square corners into slivers that
    /// point-in-triangle tests fall between — are recovered by clamping to
    /// the nearest triangle within a one-bin UV distance; feeding such
    /// nodes to the hole-fill membrane instead extrapolates the surface and
    /// poisons the fit (observed 1.8 mm grid garbage on the saddle).
    fn sample(&self, u: f64, v: f64) -> Option<[f64; 3]> {
        let bin =
            |x: f64| ((x * self.nbins as f64) as isize).clamp(0, self.nbins as isize - 1) as usize;
        let (bu, bv) = (bin(u), bin(v));
        let cands = &self.bins[bu * self.nbins + bv];
        let mut best: Option<(f64, [f64; 3])> = None;
        for &ti in cands {
            let (worst, pt) = self.barycentric(ti, u, v);
            if best.as_ref().is_none_or(|(bw, _)| worst > *bw) {
                best = Some((worst, pt));
            }
        }
        if let Some((worst, pt)) = best {
            if worst >= -1e-9 {
                return Some(pt);
            }
        }
        // Clamped fallback over the 3×3 bin neighborhood.
        let mut nearest: Option<(f64, [f64; 3])> = None;
        for du in -1_isize..=1 {
            for dv in -1_isize..=1 {
                let (nu, nv) = (bu as isize + du, bv as isize + dv);
                if nu < 0 || nv < 0 || nu >= self.nbins as isize || nv >= self.nbins as isize {
                    continue;
                }
                for &ti in &self.bins[nu as usize * self.nbins + nv as usize] {
                    if let Some((d2, pt)) = self.clamp_to_triangle(ti, u, v) {
                        if nearest.as_ref().is_none_or(|(bd, _)| d2 < *bd) {
                            nearest = Some((d2, pt));
                        }
                    }
                }
            }
        }
        let max_dist = 1.0 / self.nbins as f64;
        match nearest {
            Some((d2, pt)) if d2.sqrt() <= max_dist => Some(pt),
            _ => None,
        }
    }

    /// Barycentric coordinates of (u, v) in flattened triangle `ti` and the
    /// corresponding (unclamped) 3D interpolation. Returns the worst
    /// coordinate (≥ 0 means inside).
    fn barycentric(&self, ti: usize, u: f64, v: f64) -> (f64, [f64; 3]) {
        let tri = self.region.tris[ti];
        let [a, b, c] = tri.map(|li| self.uv[li]);
        let det = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
        if det.abs() < 1e-30 {
            return (f64::NEG_INFINITY, [0.0; 3]);
        }
        let l1 = ((u - a.0) * (c.1 - a.1) - (c.0 - a.0) * (v - a.1)) / det;
        let l2 = ((b.0 - a.0) * (v - a.1) - (u - a.0) * (b.1 - a.1)) / det;
        let l0 = 1.0 - l1 - l2;
        let p3 = tri.map(|li| self.region.positions[li]);
        let pt = [
            l0 * p3[0][0] + l1 * p3[1][0] + l2 * p3[2][0],
            l0 * p3[0][1] + l1 * p3[1][1] + l2 * p3[2][1],
            l0 * p3[0][2] + l1 * p3[1][2] + l2 * p3[2][2],
        ];
        (l0.min(l1).min(l2), pt)
    }

    /// Closest point of flattened triangle `ti` to (u, v): squared UV
    /// distance and the 3D point at the clamped location.
    fn clamp_to_triangle(&self, ti: usize, u: f64, v: f64) -> Option<(f64, [f64; 3])> {
        let (worst, pt) = self.barycentric(ti, u, v);
        if worst == f64::NEG_INFINITY {
            return None;
        }
        if worst >= 0.0 {
            return Some((0.0, pt));
        }
        // Outside: closest point lies on one of the edges.
        let tri = self.region.tris[ti];
        let uvs = tri.map(|li| self.uv[li]);
        let p3 = tri.map(|li| self.region.positions[li]);
        let mut best: Option<(f64, [f64; 3])> = None;
        for e in 0..3 {
            let (a, b) = (uvs[e], uvs[(e + 1) % 3]);
            let (pa, pb) = (p3[e], p3[(e + 1) % 3]);
            let (ex, ey) = (b.0 - a.0, b.1 - a.1);
            let len2 = ex * ex + ey * ey;
            let t = if len2 < 1e-30 {
                0.0
            } else {
                (((u - a.0) * ex + (v - a.1) * ey) / len2).clamp(0.0, 1.0)
            };
            let (qx, qy) = (a.0 + t * ex, a.1 + t * ey);
            let d2 = (u - qx) * (u - qx) + (v - qy) * (v - qy);
            if best.as_ref().is_none_or(|(bd, _)| d2 < *bd) {
                let pt = [
                    pa[0] + t * (pb[0] - pa[0]),
                    pa[1] + t * (pb[1] - pa[1]),
                    pa[2] + t * (pb[2] - pa[2]),
                ];
                best = Some((d2, pt));
            }
        }
        best
    }
}

/// Monotone reparameterizations of the two UV axes that equalize 3D arc
/// spacing. The Floater map compresses high-curvature areas, so sampling a
/// UV-uniform grid gives wildly uneven 3D point spacing; OCCT parameterizes
/// the fit by chord length *averaged across rows*, and locally
/// near-coincident stations make the interpolation system oscillate between
/// nodes (observed ~1 mm bulges on the saddle with all node residuals at
/// 1e-2). Sampling at arc-length-equalized stations keeps the fit's
/// parameterization consistent with the data everywhere.
struct AxisMaps {
    /// (uv, arc_fraction) samples of each axis's averaged cumulative arc
    /// profile G — strictly monotone.
    u_profile: Vec<(f64, f64)>,
    v_profile: Vec<(f64, f64)>,
}

impl AxisMaps {
    fn build(sampler: &UvSampler) -> AxisMaps {
        AxisMaps {
            u_profile: axis_profile(sampler, true),
            v_profile: axis_profile(sampler, false),
        }
    }

    /// Arc-length-uniform sample stations for `n` grid nodes along an axis:
    /// uv values whose arc fractions are uniform.
    fn stations(profile: &[(f64, f64)], n: usize) -> Vec<f64> {
        let mut out = Vec::with_capacity(n);
        let mut k = 0usize;
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            while k + 1 < profile.len() && profile[k + 1].1 < t {
                k += 1;
            }
            let (u0, g0) = profile[k];
            let (u1, g1) = profile[(k + 1).min(profile.len() - 1)];
            let s = if g1 > g0 { (t - g0) / (g1 - g0) } else { 0.0 };
            out.push((u0 + s.clamp(0.0, 1.0) * (u1 - u0)).clamp(0.0, 1.0));
        }
        out
    }

    /// Arc fraction G(uv) along an axis (for projection seeding: the fit's
    /// parameter is ~linear in arc fraction, not in raw UV).
    fn arc_fraction(profile: &[(f64, f64)], uv: f64) -> f64 {
        let mut k = 0usize;
        while k + 1 < profile.len() && profile[k + 1].0 < uv {
            k += 1;
        }
        let (u0, g0) = profile[k];
        let (u1, g1) = profile[(k + 1).min(profile.len() - 1)];
        let s = if u1 > u0 { (uv - u0) / (u1 - u0) } else { 0.0 };
        (g0 + s.clamp(0.0, 1.0) * (g1 - g0)).clamp(0.0, 1.0)
    }
}

/// Cumulative 3D arc profile of one UV axis, averaged over several isolines
/// of the other axis, normalized to [0, 1] and made strictly monotone. A
/// 15% blend with the identity bounds the distortion correction, keeping
/// stations well-defined where probe isolines cross data gaps.
fn axis_profile(sampler: &UvSampler, u_axis: bool) -> Vec<(f64, f64)> {
    const PROBES: usize = 9;
    const FINE: usize = 257;
    // Per-station segment lengths, summed over the probe isolines.
    let mut segs = vec![0.0_f64; FINE];
    for pi in 0..PROBES {
        let w = pi as f64 / (PROBES - 1) as f64;
        let mut prev: Option<[f64; 3]> = None;
        for (k, s) in segs.iter_mut().enumerate() {
            let t = k as f64 / (FINE - 1) as f64;
            let p = if u_axis {
                sampler.sample(t, w)
            } else {
                sampler.sample(w, t)
            };
            if let Some(p) = p {
                if let Some(q) = prev {
                    *s += ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)
                        + (p[2] - q[2]).powi(2))
                    .sqrt();
                }
                prev = Some(p);
            }
        }
    }
    let mut profile = Vec::with_capacity(FINE);
    let mut acc = 0.0_f64;
    for (k, s) in segs.iter().enumerate() {
        if k > 0 {
            acc += s;
        }
        profile.push((k as f64 / (FINE - 1) as f64, acc));
    }
    let total = acc.max(1e-30);
    for (u, g) in profile.iter_mut() {
        *g = 0.85 * (*g / total) + 0.15 * *u;
    }
    profile
}

/// Fill grid nodes not covered by the region (hole interiors) with a
/// Laplacian membrane over the grid graph, so the B-spline approximator
/// sees smooth data there. The membrane deviates from the true surface
/// inside the hole, but no mesh exists there to validate against and the
/// face is trimmed by the hole's boundary loop downstream; what matters is
/// that the data stays smooth *near* the hole rim.
fn fill_missing_grid_nodes(
    grid: &mut [Option<[f64; 3]>],
    n_u: usize,
    n_v: usize,
) -> Result<(), String> {
    let missing: Vec<usize> = (0..n_u * n_v).filter(|&i| grid[i].is_none()).collect();
    if missing.is_empty() {
        return Ok(());
    }
    if missing.len() == n_u * n_v {
        return Err("no grid node hit the region".to_string());
    }
    // Initialize each missing node from the nearest valid node in its row
    // (then column), then relax.
    let idx = |i: usize, j: usize| i * n_v + j;
    let mut values: Vec<[f64; 3]> = Vec::with_capacity(grid.len());
    for i in 0..n_u {
        for j in 0..n_v {
            let g = &grid[idx(i, j)];
            if let Some(p) = g {
                values.push(*p);
                continue;
            }
            let mut init = None;
            for r in 1..n_u.max(n_v) {
                for (ii, jj) in [
                    (i as isize - r as isize, j as isize),
                    (i as isize + r as isize, j as isize),
                    (i as isize, j as isize - r as isize),
                    (i as isize, j as isize + r as isize),
                ] {
                    if ii >= 0 && (ii as usize) < n_u && jj >= 0 && (jj as usize) < n_v {
                        if let Some(p) = grid[idx(ii as usize, jj as usize)] {
                            init = Some(p);
                            break;
                        }
                    }
                }
                if init.is_some() {
                    break;
                }
            }
            values.push(init.ok_or("hole fill found no valid node")?);
        }
    }
    // Gauss–Seidel relaxation of the missing nodes only.
    for _ in 0..2000 {
        let mut max_delta = 0.0_f64;
        for &m in &missing {
            let (i, j) = (m / n_v, m % n_v);
            let mut acc = [0.0_f64; 3];
            let mut cnt = 0.0;
            for (ii, jj) in [
                (i.wrapping_sub(1), j),
                (i + 1, j),
                (i, j.wrapping_sub(1)),
                (i, j + 1),
            ] {
                if ii < n_u && jj < n_v {
                    let p = values[idx(ii, jj)];
                    acc = [acc[0] + p[0], acc[1] + p[1], acc[2] + p[2]];
                    cnt += 1.0;
                }
            }
            let nu = [acc[0] / cnt, acc[1] / cnt, acc[2] / cnt];
            let old = values[m];
            max_delta = max_delta
                .max((nu[0] - old[0]).abs())
                .max((nu[1] - old[1]).abs())
                .max((nu[2] - old[2]).abs());
            values[m] = nu;
        }
        if max_delta < 1e-9 {
            break;
        }
    }
    for (g, v) in grid.iter_mut().zip(values) {
        *g = Some(v);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Fitting
// ---------------------------------------------------------------------------

/// Fit one freeform region with a B-spline surface. Returns the validated
/// hypothesis or a reason string when the region can't be fitted (the caller
/// leaves it as single-face planar surfaces).
fn fit_component(
    mesh: &ConnectedMesh,
    component: Vec<usize>,
    face_to_selected: &[usize],
    in_region: &[bool],
    config: &Config,
) -> Result<Vec<BSplineHypothesis>, String> {
    let region = build_region(mesh, component, face_to_selected, in_region)?;
    // A closed band (e.g. the wall of a closed prism) cannot be flattened
    // to a disk; cut it along a seam path first. The band is then fitted
    // PERIODICALLY around the cut (one surface, one face, a real seam edge
    // in stage 3.4 — the form CAD kernels write). The clamped fit + seam
    // split into two faces remains as a fallback for closed bands the
    // periodic fit cannot validate.
    let cut = is_closed_band(&region);
    let region = if cut {
        cut_region_seam(region, mesh, config.verbosity)?
    } else {
        region
    };
    let outer = outer_loop_index(&region);
    if config.verbosity >= 2 {
        let loop_lens: Vec<usize> = region.loops.iter().map(|l| l.len()).collect();
        eprintln!(
            "  Stage 2.7: region {} faces, {} verts, {} tris, boundary loops {:?} (outer = #{outer}{})",
            region.faces.len(),
            region.verts.len(),
            region.tris.len(),
            loop_lens,
            if cut { ", seam-cut" } else { "" },
        );
    }
    if config.verbosity >= 3 {
        let pts: Vec<[f64; 3]> = region.loops[outer]
            .iter()
            .step_by((region.loops[outer].len() / 24).max(1))
            .map(|&(li, _)| region.positions[li])
            .collect();
        eprintln!("  Stage 2.7: outer loop samples: {pts:.3?}");
    }

    // Parameterize: outer loop pinned to the unit square, everything else
    // (interior and hole-boundary vertices) solved with mean-value weights.
    let nv = region.verts.len();
    let mut uv = vec![(0.0_f64, 0.0_f64); nv];
    let mut pinned = vec![false; nv];
    let lp = &region.loops[outer];
    let corners = choose_corners(lp, &region.positions);
    if corners.len() != 4 {
        return Err(format!(
            "could not pin 4 boundary corners (found {})",
            corners.len()
        ));
    }
    map_outer_loop_to_square(lp, &corners, &region.positions, &mut uv);
    for &(li, _) in lp {
        pinned[li] = true;
    }
    solve_interior_uv(&region, &pinned, &mut uv, config.verbosity)?;

    // Grid resolution baseline: the mesh's own boundary density (sides of
    // the square). Oversampled to avoid aliasing between vertices; the
    // retry ladder below raises density when validation fails (the Floater
    // map compresses high-curvature areas, so a uniform UV grid can
    // undersample them in 3D).
    let n = lp.len();
    let side_len = |s: usize| (corners[(s + 1) % 4] + n - corners[s]) % n + 1;
    let base_u = side_len(0).max(side_len(2));
    let base_v = side_len(1).max(side_len(3));

    if config.verbosity >= 2 {
        let corner_pts: Vec<[f64; 3]> = corners
            .iter()
            .map(|&c| region.positions[lp[c].0])
            .collect();
        eprintln!(
            "  Stage 2.7: sides {}x{} verts, corners at {:?}",
            base_u, base_v, corner_pts,
        );
    }

    // Pinned positions of the boundary vertices along each pair of opposite
    // sides: the base grid stations (see `axis_stations`).
    let side_coords = |s: usize, take_u: bool| -> Vec<f64> {
        let n = lp.len();
        let len = (corners[(s + 1) % 4] + n - corners[s]) % n;
        let mut out: Vec<f64> = (0..=len)
            .map(|e| {
                let li = lp[(corners[s] + e) % n].0;
                if take_u { uv[li].0 } else { uv[li].1 }
            })
            .collect();
        out.sort_by(|a, b| a.partial_cmp(b).unwrap());
        out
    };
    let (u_side_a, u_side_b) = (side_coords(0, true), side_coords(2, true));
    let (v_side_a, v_side_b) = (side_coords(1, false), side_coords(3, false));

    let sampler = UvSampler::new(&region, &uv);
    let maps = AxisMaps::build(&sampler);
    // Fit modes, in order of preference: a cut band tries the periodic fit
    // around the band first, then falls back to the clamped fit that gets
    // seam-split into two faces sharing the surface.
    let modes: Vec<Option<bool>> = if cut {
        vec![Some(cut_around_axis(&region, &uv)?), None]
    } else {
        vec![None]
    };
    let mut last_err = String::new();
    for around in modes {
        for (density, cap) in [(1.5, 121), (2.5, 181), (3.5, 241)] {
            // Interior data density sets the enrichment target: an axis's
            // boundary can be much coarser than the region it bounds (flat
            // margins tessellate coarsely — an 8-vertex side on a 946-vertex
            // region).
            let target = ((region.verts.len() as f64).sqrt() * density).round() as usize;
            let target = target.clamp(8, cap);
            let u_stations = axis_stations(&u_side_a, &u_side_b, target, &maps.u_profile);
            let v_stations = axis_stations(&v_side_a, &v_side_b, target, &maps.v_profile);
            match fit_grid(
                mesh, &region, &sampler, &maps, &u_stations, &v_stations, around, config,
            ) {
                Ok(hyp) => {
                    return if cut && around.is_none() {
                        if config.verbosity >= 1 {
                            eprintln!(
                                "  Stage 2.7: periodic fit failed for the closed band \
                                 ({last_err}); falling back to the clamped seam-split"
                            );
                        }
                        split_cut_hypothesis(hyp, mesh, &region, &uv)
                    } else {
                        Ok(vec![hyp])
                    };
                }
                Err(e) => {
                    if config.verbosity >= 2 {
                        eprintln!(
                            "  Stage 2.7: {}x{}{} grid attempt failed: {e}",
                            u_stations.len(),
                            v_stations.len(),
                            if around.is_some() { " periodic" } else { "" },
                        );
                    }
                    // A collapsed parameterization won't improve with density.
                    if e.contains("parameterization covers") {
                        return Err(e);
                    }
                    last_err = e;
                }
            }
        }
    }
    Err(last_err)
}

/// Which unit-square axis runs around a seam-cut band (from cut copy to cut
/// copy): the cut vertices sit at that axis's parameter extremes.
fn cut_around_axis(region: &Region, uv: &[(f64, f64)]) -> Result<bool, String> {
    let cut_uv = region
        .loops
        .iter()
        .flatten()
        .find(|&&(_, out)| out == CUT_OUTSIDE)
        .map(|&(li, _)| uv[li])
        .ok_or("seam-cut region has no cut boundary")?;
    Ok(cut_uv.0.min(1.0 - cut_uv.0) < cut_uv.1.min(1.0 - cut_uv.1))
}

/// FALLBACK for closed bands the periodic fit cannot validate: split the
/// clamped fitted hypothesis of a seam-cut band into two halves along the
/// around-the-band axis. Both halves share the SAME fitted surface object.
/// Each half spans only part of the surface's parameter range, so its two
/// cut edges are geometrically distinct and no periodic-seam machinery is
/// needed — but those two edges are reconstruction artifacts a reference
/// model (one periodic face with a real seam) does not have, so the
/// stage 3.1 edge compare will flag them.
fn split_cut_hypothesis(
    hyp: BSplineHypothesis,
    mesh: &ConnectedMesh,
    region: &Region,
    uv: &[(f64, f64)],
) -> Result<Vec<BSplineHypothesis>, String> {
    // The around axis runs from cut copy to cut copy.
    let around_u = cut_around_axis(region, uv)?;

    // Mean around-coordinate per region face, from its triangles.
    let mut face_sum = vec![(0.0_f64, 0usize); region.faces.len()];
    for (ti, tri) in region.tris.iter().enumerate() {
        let rf = region.tri_face[ti];
        for &li in tri {
            face_sum[rf].0 += if around_u { uv[li].0 } else { uv[li].1 };
            face_sum[rf].1 += 1;
        }
    }
    let mut faces_a = Vec::new();
    let mut faces_b = Vec::new();
    for (rf, &fi) in region.faces.iter().enumerate() {
        let mean = face_sum[rf].0 / face_sum[rf].1.max(1) as f64;
        if std::env::var_os("BREPPER_DEBUG_SPLIT").is_some() && !(0.05..=0.95).contains(&mean) {
            let corners: Vec<(usize, f64, f64)> = region
                .tris
                .iter()
                .enumerate()
                .filter(|(ti, _)| region.tri_face[*ti] == rf)
                .flat_map(|(_, tri)| tri.iter().map(|&li| (li, uv[li].0, uv[li].1)))
                .collect();
            eprintln!("      split face {fi}: mean={mean:.4} corners={corners:.4?}");
        }
        if mean < 0.5 {
            faces_a.push(fi);
        } else {
            faces_b.push(fi);
        }
    }
    if faces_a.is_empty() || faces_b.is_empty() {
        return Err("seam split left one half empty".to_string());
    }
    let verts_of = |faces: &[usize]| -> Vec<usize> {
        let mut vs: Vec<usize> = faces
            .iter()
            .flat_map(|&fi| {
                let f = &mesh.faces[fi];
                f.vertex_indices[..f.vertex_count as usize].to_vec()
            })
            .collect();
        vs.sort_unstable();
        vs.dedup();
        vs
    };

    let surface_b = hyp
        .surface
        .to_handle_surface()
        .downcast_to_b_spline_surface()
        .ok_or("could not mint a second handle to the fitted surface")?;
    let mut hyp_b = BSplineHypothesis::new(surface_b, faces_b.clone(), verts_of(&faces_b));
    hyp_b.error_max = hyp.error_max;
    hyp_b.centroid_error_max = hyp.centroid_error_max;
    hyp_b.error_abs_sum = hyp.error_abs_sum;

    let mut hyp_a = hyp;
    hyp_a.vertices = verts_of(&faces_a);
    hyp_a.faces = faces_a;

    Ok(vec![hyp_a, hyp_b])
}

/// Grid stations for one axis. The base stations are the boundary
/// vertices' own pinned positions (averaged between the two opposite
/// sides — structured or not, the boundary is where the region has exact
/// on-surface data at exact parameters): stations *between* the vertices
/// of a coarse profile sample chord interiors and pin the interpolant to
/// the chordal prism (a 1.7% volume deficit on the arch, the mesh's own
/// chordal error). The base is then evenly subdivided up to the interior-
/// density target when the boundary is coarser than the region it bounds
/// (flat margins tessellate coarsely: hill_with_hole has 8-vertex sides on
/// a 946-vertex region, and an 8×8 grid cannot carry the bump). Sides with
/// mismatched vertex counts fall back to arc-length-uniform stations.
fn axis_stations(
    side_a: &[f64],
    side_b: &[f64],
    target: usize,
    profile: &[(f64, f64)],
) -> Vec<f64> {
    if side_a.len() != side_b.len() || side_a.len() < 2 {
        return AxisMaps::stations(profile, target);
    }
    let base: Vec<f64> = side_a
        .iter()
        .zip(side_b)
        .map(|(a, b)| ((a + b) / 2.0).clamp(0.0, 1.0))
        .collect();
    let k = (target.saturating_sub(1)).div_ceil(base.len() - 1).max(1);
    let mut out = Vec::with_capacity((base.len() - 1) * k + 1);
    for w in base.windows(2) {
        for step in 0..k {
            out.push(w[0] + (w[1] - w[0]) * step as f64 / k as f64);
        }
    }
    out.push(*base.last().unwrap());
    out
}

/// One fit attempt at fixed grid stations: resample, orient, fit,
/// validate against `--surface-tolerance`.
///
/// `around` selects the fit mode: None = clamped fit of an open (disk)
/// region; Some(around_u) = periodic fit of a seam-cut closed band whose
/// around-the-band unit-square axis is u (true) or v (false).
#[allow(clippy::too_many_arguments)]
fn fit_grid(
    mesh: &ConnectedMesh,
    region: &Region,
    sampler: &UvSampler,
    maps: &AxisMaps,
    u_stations: &[f64],
    v_stations: &[f64],
    around: Option<bool>,
    config: &Config,
) -> Result<BSplineHypothesis, String> {
    let (n_u, n_v) = (u_stations.len(), v_stations.len());
    let mut grid: Vec<Option<[f64; 3]>> = Vec::with_capacity(n_u * n_v);
    for &su in u_stations {
        for &sv in v_stations {
            grid.push(sampler.sample(su, sv));
        }
    }
    let covered = grid.iter().filter(|g| g.is_some()).count();
    if config.verbosity >= 3 && covered < grid.len() {
        let mut missing_uvs: Vec<(f64, f64)> = Vec::new();
        for i in 0..n_u {
            for j in 0..n_v {
                if grid[i * n_v + j].is_none() {
                    missing_uvs.push((u_stations[i], v_stations[j]));
                }
            }
        }
        eprintln!(
            "  Stage 2.7: {}/{} grid nodes missing, first few at {:?}",
            grid.len() - covered,
            grid.len(),
            &missing_uvs[..missing_uvs.len().min(8)],
        );
    }
    // A mostly-empty grid means the flattening collapsed (e.g. a closed
    // periodic wall parameterized as an annulus) — not fittable this way.
    if covered * 2 < grid.len() {
        return Err(format!(
            "parameterization covers only {covered}/{} grid nodes",
            grid.len()
        ));
    }
    fill_missing_grid_nodes(&mut grid, n_u, n_v)?;
    let grid: Vec<[f64; 3]> = grid.into_iter().map(Option::unwrap).collect();

    // Orient the grid so the fitted surface's natural normal (d1u × d1v)
    // points outward, matching the mesh winding normals. Fixing the
    // convention at fit time (one row flip) beats recovering it later with
    // per-face heuristics.
    let idx = |i: usize, j: usize| i * n_v + j;
    let (ci, cj) = (n_u / 2, n_v / 2);
    let du = sub3(grid[idx(ci + 1, cj)], grid[idx(ci - 1, cj)]);
    let dv = sub3(grid[idx(ci, cj + 1)], grid[idx(ci, cj - 1)]);
    let grid_n = cross3(du, dv);
    let center = grid[idx(ci, cj)];
    let mesh_n = region_normal_near(mesh, region, center);
    let flip = dot3(grid_n, mesh_n) < 0.0;

    // Array layout: the fitted surface's U is the array's FIRST index. A
    // periodic fit needs the around-the-band axis there (transpose when the
    // around axis is v), with its closing station — the second cut copy,
    // whose column duplicates the first — omitted: OCCT's periodic
    // Interpolate appends the wrap-around copy of the first row itself
    // (probed, see DEVELOPMENT_PLAN §2.7). Orientation: transposing negates
    // d1u×d1v, so the outward flip reverses the SECOND index and accounts
    // for the transpose.
    let (n_a, n_b, transpose) = match around {
        Some(false) => (n_v - 1, n_u, true),
        Some(true) => (n_u - 1, n_v, false),
        None => (n_u, n_v, false),
    };
    let flip_b = flip != transpose;
    let mut pts = t_colgp::Array2OfPnt::new_int4(1, n_a as i32, 1, n_b as i32);
    for a in 0..n_a {
        for b in 0..n_b {
            let bb = if flip_b { n_b - 1 - b } else { b };
            let (i, j) = if transpose { (bb, a) } else { (a, bb) };
            let p = grid[idx(i, j)];
            let pnt = gp::Pnt::new_real3(p[0], p[1], p[2]);
            pts.set_value(a as i32 + 1, b as i32 + 1, &pnt);
        }
    }
    // Interpolation first: the grid is sampled from an exact-on-surface
    // mesh, so a C2 cubic surface through every node is the most faithful
    // reconstruction (OCCT's least-squares approximator reports IsDone
    // while overshooting its Tol3D by millimeters near strongly curved
    // patch boundaries — measured on the saddle model). The approximator
    // remains as a fallback for data interpolation can't handle — except on
    // the periodic path, which is interpolation-only (the approximator
    // binding carries no periodic flag).
    let approx_tols = [
        (config.surface_tolerance_mm / 8.0).max(1e-4),
        (config.surface_tolerance_mm / 2.0).max(1e-4),
    ];
    let methods: &[Option<f64>] = if around.is_some() {
        &[None]
    } else {
        &[None, Some(approx_tols[0]), Some(approx_tols[1])]
    };
    let mut last_err = String::new();
    for &method in methods {
        let surface = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let fit = match method {
                None => {
                    let mut fit = geom_api::PointsToBSplineSurface::new();
                    fit.interpolate_array2ofpnt_parametrizationtype_bool(
                        &pts,
                        approx::ParametrizationType::Chordlength,
                        around.is_some(),
                    );
                    fit
                }
                Some(tol3d) => {
                    geom_api::PointsToBSplineSurface::new_array2ofpnt_parametrizationtype_int2_shape_real(
                        &pts,
                        approx::ParametrizationType::Chordlength,
                        3,
                        8,
                        geom_abs::Shape::C2,
                        tol3d,
                    )
                }
            };
            if fit.is_done() {
                fit.surface().to_handle_surface().downcast_to_b_spline_surface()
            } else {
                None
            }
        }));
        let mut surface = match surface {
            Ok(Some(s)) => s,
            Ok(None) => {
                last_err = "GeomAPI_PointsToBSplineSurface did not converge".to_string();
                continue;
            }
            Err(_) => {
                last_err = "GeomAPI_PointsToBSplineSurface raised an exception".to_string();
                continue;
            }
        };
        if around.is_some() {
            // Rescale the periodic (around) axis's knots to exactly [0, 2π]
            // — an exact affine reparameterization (probed: geometry drift
            // ~1e-15). Stage 3.4's wrapped-face machinery (branch unwrap,
            // winding, pcurve shifts, seam placement) hard-codes that
            // period.
            let bs = surface.get_mut();
            if !bs.is_u_periodic() {
                last_err = "periodic interpolation returned a non-periodic surface".to_string();
                continue;
            }
            let nk = bs.nb_u_knots();
            let (k0, k1) = (bs.u_knot(1), bs.u_knot(nk));
            let mut knots = t_col_std::Array1OfReal::new_int2(1, nk);
            for k in 1..=nk {
                let scaled = (bs.u_knot(k) - k0) / (k1 - k0) * std::f64::consts::TAU;
                knots.set_value(k, &scaled);
            }
            bs.set_u_knots(&knots);
        }
        match validate_fit(
            mesh, region, sampler.uv, maps, around, flip_b, surface, n_u, n_v, &grid, config,
        ) {
            Ok(hyp) => return Ok(hyp),
            Err(e) => last_err = e,
        }
    }
    Err(last_err)
}

/// Validate a fitted surface: every region vertex and face centroid must
/// lie within surface tolerance. Each projection is seeded from the
/// vertex's own Floater UV (`uv`, indexed like `region.verts`), the one
/// reliable way to find the true nearest point on a many-knot surface.
#[allow(clippy::too_many_arguments)]
fn validate_fit(
    mesh: &ConnectedMesh,
    region: &Region,
    uv: &[(f64, f64)],
    maps: &AxisMaps,
    around: Option<bool>,
    flip_b: bool,
    surface: OwnedPtr<geom::HandleGeomBSplineSurface>,
    n_u: usize,
    n_v: usize,
    grid: &[[f64; 3]],
    config: &Config,
) -> Result<BSplineHypothesis, String> {
    // Hypothesis vertices are unique GLOBAL mesh vertices (region.verts may
    // repeat globals for seam-cut clones).
    let mut hyp_verts = region.verts.clone();
    hyp_verts.sort_unstable();
    hyp_verts.dedup();
    let mut hyp = BSplineHypothesis::new(surface, region.faces.clone(), hyp_verts);
    hyp.u_periodic = around.is_some();
    let tol = config.surface_tolerance_mm;
    // The fit's parameter is ~linear in arc fraction (the grid stations
    // equalize arc). Map the Floater axes onto the fitted surface's axes:
    // transposed when the band's around axis is v (the surface's U is
    // always the array's first index), with the second axis reversed when
    // the grid was flipped for outward orientation.
    let seed_of = |li: usize| -> (f64, f64) {
        let (u, v) = uv[li];
        let gu = AxisMaps::arc_fraction(&maps.u_profile, u);
        let gv = AxisMaps::arc_fraction(&maps.v_profile, v);
        let (ga, gb) = if around == Some(false) { (gv, gu) } else { (gu, gv) };
        (ga, if flip_b { 1.0 - gb } else { gb })
    };
    if config.verbosity >= 3 {
        // Probe: the fit's residual against its own input grid. Large gaps
        // here mean the OCCT fit itself missed the data; small gaps with
        // failing vertex validation implicate the sampler or the
        // validation projection instead.
        let mut max_gap = 0.0_f64;
        let mut max_at = [0.0_f64; 3];
        for p in grid.iter().step_by(3) {
            let d = hyp.distance(p);
            if d > max_gap {
                max_gap = d;
                max_at = *p;
            }
        }
        eprintln!(
            "  Stage 2.7: {n_u}x{n_v} fit (deg {}x{}, {}x{} poles) residual vs own grid: \
             max {max_gap:.3e} mm at [{:.3},{:.3},{:.3}]",
            hyp.u_degree, hyp.v_degree, hyp.nb_u_poles, hyp.nb_v_poles,
            max_at[0], max_at[1], max_at[2],
        );
    }
    for (li, &vi) in region.verts.iter().enumerate() {
        let p = region.positions[li];
        let d = hyp.distance_seeded(&p, seed_of(li));
        hyp.error_max = hyp.error_max.max(d);
        hyp.error_abs_sum += d;
        if hyp.error_max > tol {
            if config.verbosity >= 3 {
                // Fold detector: sampling the flattening at this vertex's
                // own UV must reproduce the vertex; anything else means an
                // overlapping (folded) triangle wins the point location.
                let (u, vv) = uv[li];
                if let Some(smp) = UvSampler::new(region, uv).sample(u, vv) {
                    eprintln!(
                        "  Stage 2.7: vertex {vi} uv=({u:.4},{vv:.4}); sampler there gives \
                         [{:.3},{:.3},{:.3}] (vertex is [{:.3},{:.3},{:.3}])",
                        smp[0], smp[1], smp[2], p[0], p[1], p[2],
                    );
                }
            }
            return Err(format!(
                "fit misses vertex {vi} [{:.3},{:.3},{:.3}] by {:.3e} mm (tolerance {tol:.3e})",
                p[0], p[1], p[2], hyp.error_max
            ));
        }
    }
    for &fi in &region.faces {
        let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
        let seed_li = region.global_to_local[&mesh.faces[fi].vertex_indices[0]];
        let d = hyp.distance_seeded(&c, seed_of(seed_li));
        hyp.centroid_error_max = hyp.centroid_error_max.max(d);
        if hyp.centroid_error_max > tol {
            return Err(format!(
                "fit misses face {fi} centroid by {:.3e} mm (tolerance {tol:.3e})",
                hyp.centroid_error_max
            ));
        }
    }
    Ok(hyp)
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Outward mesh normal of the region face whose centroid is nearest to `p`.
fn region_normal_near(mesh: &ConnectedMesh, region: &Region, p: [f64; 3]) -> [f64; 3] {
    let mut best = [0.0, 0.0, 1.0];
    let mut best_d = f64::INFINITY;
    for &fi in &region.faces {
        let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
        let d = (c[0] - p[0]).powi(2) + (c[1] - p[1]).powi(2) + (c[2] - p[2]).powi(2);
        if d < best_d {
            if let Some(n) = mesh.faces[fi].normal {
                best_d = d;
                best = n;
            }
        }
    }
    best
}

// ---------------------------------------------------------------------------
// Stage 2.7 driver
// ---------------------------------------------------------------------------

/// Stage 2.7 entry: gather connected components of single-face planar
/// selections and fit each with a B-spline surface (or, for regions that are
/// planar to fitting exactness, a multi-face plane). Successful fits replace
/// the component's single-face planar selections in `selected_surfaces`.
#[allow(clippy::too_many_arguments)]
pub(super) fn deduce_bspline_surfaces(
    mesh: &ConnectedMesh,
    selected_surfaces: &mut Vec<SelectedSurface>,
    planar_hypotheses: &mut Vec<PlanarHypothesis>,
    conical_hypotheses: &[ConicalHypothesis],
    spherical_hypotheses: &[SphericalHypothesis],
    toroidal_hypotheses: &[ToroidalHypothesis],
    bspline_hypotheses: &mut Vec<BSplineHypothesis>,
    viz: Option<&VizSender>,
    config: &Config,
) {
    let t = Instant::now();
    let num_faces = mesh.faces.len();

    // Face → selected surface index, and the single-face planar face set.
    let mut face_to_selected = vec![usize::MAX; num_faces];
    for (si, s) in selected_surfaces.iter().enumerate() {
        let faces: &[usize] = match s {
            SelectedSurface::Planar(i) => &planar_hypotheses[*i].faces,
            SelectedSurface::Conical(i) => &conical_hypotheses[*i].faces,
            SelectedSurface::Spherical(i) => &spherical_hypotheses[*i].faces,
            SelectedSurface::Toroidal(i) => &toroidal_hypotheses[*i].faces,
            SelectedSurface::BSpline(i) => &bspline_hypotheses[*i].faces,
        };
        for &fi in faces {
            face_to_selected[fi] = si;
        }
    }
    let selection_faces = |si: usize| -> &[usize] {
        match selected_surfaces[si] {
            SelectedSurface::Planar(i) => &planar_hypotheses[i].faces,
            SelectedSurface::Conical(i) => &conical_hypotheses[i].faces,
            SelectedSurface::Spherical(i) => &spherical_hypotheses[i].faces,
            SelectedSurface::Toroidal(i) => &toroidal_hypotheses[i].faces,
            SelectedSurface::BSpline(i) => &bspline_hypotheses[i].faces,
        }
    };
    let is_single_planar = |si: usize| -> bool {
        matches!(
            selected_surfaces[si],
            SelectedSurface::Planar(pi) if planar_hypotheses[pi].faces.len() == 1
        )
    };
    let single: Vec<bool> = (0..num_faces)
        .map(|fi| face_to_selected[fi] != usize::MAX && is_single_planar(face_to_selected[fi]))
        .collect();
    let candidate: Vec<bool> = (0..num_faces)
        .map(|fi| {
            face_to_selected[fi] != usize::MAX
                && selection_faces(face_to_selected[fi]).len() <= SMALL_SELECTION_ABSORB_FACES
        })
        .collect();

    // Connected components over mesh adjacency, lowest face index first.
    // Components grow only across SMOOTH edges: a freeform surface varies
    // by at most the tessellation's angular step between neighbors (the
    // same angular tolerance the analytic BFS gathering uses), while a
    // sharp crease marks a real B-rep edge between distinct surfaces. This
    // keeps legitimately planar orphans (a single-quad bottom face, a 45°
    // chamfer strip) from gluing freeform regions — or entire chamfered
    // cubes — into impossible components.
    let smooth_cos = config.angular_tolerance_rad.min(std::f64::consts::PI).cos();
    let smooth_edge = |a: usize, b: usize| -> bool {
        match (mesh.faces[a].normal, mesh.faces[b].normal) {
            (Some(na), Some(nb)) => {
                na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2] >= smooth_cos
            }
            _ => false,
        }
    };
    let mut visited = vec![false; num_faces];
    let mut components: Vec<Vec<usize>> = Vec::new();
    // Planes already tangent-merged into an earlier component (a plane
    // must not be consumed by two components).
    let mut merged_planes: HashSet<usize> = HashSet::new();
    for start in 0..num_faces {
        if !candidate[start] || visited[start] {
            continue;
        }
        let mut comp = Vec::new();
        let mut stack = vec![start];
        visited[start] = true;
        while let Some(fi) = stack.pop() {
            comp.push(fi);
            let face = &mesh.faces[fi];
            for k in 0..face.vertex_count as usize {
                let ni = face.neighbors[k];
                if ni >= 0 {
                    let ni = ni as usize;
                    if candidate[ni] && !visited[ni] && smooth_edge(fi, ni) {
                        visited[ni] = true;
                        stack.push(ni);
                    }
                }
            }
        }
        comp.sort_unstable();
        // A component with no single-face planar member is a real small
        // feature caught by the widened candidate set — leave it.
        let n_single = comp.iter().filter(|&&fi| single[fi]).count();
        if config.verbosity >= 2 {
            eprintln!(
                "  Stage 2.7: component of {} faces ({} single-face planar, {} absorbed)",
                comp.len(),
                n_single,
                comp.len() - n_single,
            );
        }
        if n_single == 0 {
            continue;
        }
        // Small analytic/planar selections are only *evidence of accident*
        // when the unexplained (single-face) area dominates. When absorbed
        // faces would outnumber the single-face leftovers, trust the
        // analytic fits and keep just the leftovers.
        if comp.len() - n_single > n_single {
            comp.retain(|&fi| single[fi]);
        }
        // Merge in multi-face planar selections whose ENTIRE shared
        // boundary with the component is smooth: a plane tangent to a
        // freeform region along its whole junction is the flat part of the
        // same surface (B-spline local support makes exactly-planar
        // margins, e.g. bspline_blister), and keeping it separate invents a
        // feature edge the reference doesn't have. Full-boundary tangency
        // is the discriminator — real planar neighbors (walls, bottoms)
        // meet freeform regions across at least some sharp edges. Analytic
        // models are untouched: they have no freeform components at all.
        loop {
            let comp_set: HashSet<usize> = comp.iter().copied().collect();
            let mut boundary_planes: Vec<usize> = Vec::new();
            for &fi in &comp {
                let face = &mesh.faces[fi];
                for k in 0..face.vertex_count as usize {
                    let ni = face.neighbors[k];
                    if ni < 0 {
                        continue;
                    }
                    let ni = ni as usize;
                    if comp_set.contains(&ni) {
                        continue;
                    }
                    let si = face_to_selected[ni];
                    if si != usize::MAX
                        && !merged_planes.contains(&si)
                        && matches!(selected_surfaces[si], SelectedSurface::Planar(_))
                        && selection_faces(si).len() > 1
                    {
                        boundary_planes.push(si);
                    }
                }
            }
            boundary_planes.sort_unstable();
            boundary_planes.dedup();
            let mut grew = false;
            for si in boundary_planes {
                let pfaces = selection_faces(si);
                let all_smooth = pfaces.iter().all(|&pf| {
                    let face = &mesh.faces[pf];
                    (0..face.vertex_count as usize).all(|k| {
                        let ni = face.neighbors[k];
                        if ni < 0 {
                            return true;
                        }
                        let ni = ni as usize;
                        !comp_set.contains(&ni) || smooth_edge(pf, ni)
                    })
                });
                if all_smooth {
                    comp.extend_from_slice(pfaces);
                    merged_planes.insert(si);
                    grew = true;
                }
            }
            if !grew {
                break;
            }
            comp.sort_unstable();
            comp.dedup();
        }

        // Defensive: never consume a selection partially (its remaining
        // faces would end up covered by no surface). Same-selection faces
        // are mutually smooth-adjacent in practice, so a touched selection
        // is normally entirely inside the component.
        let comp_set: HashSet<usize> = comp.iter().copied().collect();
        let mut selections: Vec<usize> = comp.iter().map(|&fi| face_to_selected[fi]).collect();
        selections.sort_unstable();
        selections.dedup();
        for si in selections {
            let faces = selection_faces(si);
            if !faces.iter().all(|fi| comp_set.contains(fi)) {
                comp.retain(|fi| !faces.contains(fi));
            }
        }
        components.push(comp);
    }

    let mut fitted = 0usize;
    let mut planarized = 0usize;
    let mut skipped_small = 0usize;
    let mut failed = 0usize;
    let mut covered_faces = 0usize;
    // Selected-surface indices consumed by fits, removed in one pass at the
    // end: `face_to_selected` holds positions in the *original*
    // `selected_surfaces` order, and new surfaces are only appended, so the
    // recorded indices stay valid until the final retain.
    let mut consumed_selected: HashSet<usize> = HashSet::new();

    for comp in components {
        if comp.len() < MIN_BSPLINE_FACES {
            skipped_small += 1;
            continue;
        }

        let comp_len = comp.len();
        let vert_set: HashSet<usize> = comp
            .iter()
            .flat_map(|&fi| {
                let f = &mesh.faces[fi];
                f.vertex_indices[..f.vertex_count as usize].to_vec()
            })
            .collect();

        // Nearly planar regions get a plane, not a B-spline that overfits
        // noise. "Nearly" means planar to fitting exactness (the absorb
        // pass's 10× vertex tolerance), not to chordal surface tolerance —
        // a shallow dome must not be flattened. The area-weighted normal of
        // a closed or symmetric region cancels to ~zero and normalizes to
        // NaN, which f64::max would silently ignore — treat any non-finite
        // plane as non-planar.
        let (pn, pd) = planar::fit_plane(&comp, &vert_set, &mesh.faces, &mesh.vertices);
        let plane_tol = 10.0 * config.vertex_tolerance_mm;
        let plane_err = if pn.iter().all(|c| c.is_finite()) && pd.is_finite() {
            vert_set
                .iter()
                .map(|&vi| {
                    let v = &mesh.vertices[vi];
                    (pn[0] * v.x + pn[1] * v.y + pn[2] * v.z - pd).abs()
                })
                .fold(0.0_f64, f64::max)
        } else {
            f64::INFINITY
        };
        if plane_err <= plane_tol {
            let mut verts: Vec<usize> = vert_set.into_iter().collect();
            verts.sort_unstable();
            let pi = planar_hypotheses.len();
            planar_hypotheses.push(PlanarHypothesis {
                normal: pn,
                distance: pd,
                faces: comp.clone(),
                vertices: verts,
                error_max: plane_err,
                error_min: -plane_err,
                error_abs_sum: 0.0,
            });
            consumed_selected.extend(comp.iter().map(|&fi| face_to_selected[fi]));
            selected_surfaces.push(SelectedSurface::Planar(pi));
            planarized += 1;
            covered_faces += comp_len;
            if config.verbosity >= 1 {
                eprintln!(
                    "  Stage 2.7: {comp_len}-face region is planar to {plane_err:.2e} mm — \
                     selected a plane instead of a B-spline"
                );
            }
            continue;
        }

        let in_region: Vec<bool> = {
            let comp_set: HashSet<usize> = comp.iter().copied().collect();
            (0..num_faces).map(|fi| comp_set.contains(&fi)).collect()
        };
        match fit_component(mesh, comp.clone(), &face_to_selected, &in_region, config) {
            Ok(hyps) => {
                let split = hyps.len() == 2;
                if config.verbosity >= 1 {
                    let hyp = &hyps[0];
                    eprintln!(
                        "  Stage 2.7: fitted B-spline (deg {}×{}, {}×{} poles) to {} faces, \
                         vtx_err_max={:.2e}, cen_err_max={:.2e}{}",
                        hyp.u_degree, hyp.v_degree, hyp.nb_u_poles, hyp.nb_v_poles,
                        comp_len, hyp.error_max, hyp.centroid_error_max,
                        if split {
                            " (closed band: seam-split fallback, two faces)"
                        } else if hyp.u_periodic {
                            " (closed band: one periodic face)"
                        } else {
                            ""
                        },
                    );
                }
                if let Some(viz_sender) = viz {
                    let hyp = &hyps[0];
                    let mut overlay = viz::VizOverlay::new();
                    overlay.face_highlights.push(viz::FaceHighlight {
                        face_indices: comp.clone(),
                        color: [0.0, 0.8, 0.0, 1.0],
                    });
                    overlay.status_text = format!(
                        "Stage 2.7: fitted B-spline (deg {}×{}, {}×{} poles) to {} faces \
                         [space=next]",
                        hyp.u_degree, hyp.v_degree, hyp.nb_u_poles, hyp.nb_v_poles, comp_len,
                    );
                    // On Quit just stop showing overlays; the fit results
                    // must still be committed consistently below.
                    let _ = viz_sender.show_and_wait(overlay) == VizAction::Quit;
                }
                for hyp in hyps {
                    let bi = bspline_hypotheses.len();
                    bspline_hypotheses.push(hyp);
                    selected_surfaces.push(SelectedSurface::BSpline(bi));
                }
                consumed_selected.extend(comp.iter().map(|&fi| face_to_selected[fi]));
                fitted += 1;
                covered_faces += comp_len;
            }
            Err(reason) => {
                failed += 1;
                if !config.quiet {
                    eprintln!(
                        "  Stage 2.7 WARNING: {comp_len}-face freeform region left as \
                         single-face planes: {reason}"
                    );
                }
            }
        }
    }

    // Drop the consumed single-face planar selections in one pass (appended
    // surfaces sit past the original length and are never in the set).
    if !consumed_selected.is_empty() {
        let mut si = 0usize;
        selected_surfaces.retain(|_| {
            let keep = !consumed_selected.contains(&si);
            si += 1;
            keep
        });
    }

    if config.verbosity >= 1 || fitted + planarized + failed > 0 {
        eprintln!(
            "Stage 2.7 ({:.3}s): Fitted {} B-spline surface(s) covering {} faces \
             ({} region(s) planarized, {} failed, {} below {} faces)",
            t.elapsed().as_secs_f64(),
            fitted,
            covered_faces,
            planarized,
            failed,
            skipped_small,
            MIN_BSPLINE_FACES,
        );
    }
}

