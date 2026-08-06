//! `--invariant-report`: measure the DEVELOPMENT_PLAN "Invariants" section
//! against ground truth. **Diagnostics only — nothing downstream reads this.**
//!
//! The proposed invariants are statements of the form "every facet of a true
//! surface satisfies X" and "the discovered surface does not depend on where
//! you started". Both are only measurable against a reference that says which
//! facets belong to which surface, which is exactly what
//! [`crate::reference::label_faces`] produces for free from the `--compare`
//! distance query. So this module is a second consumer of that ground truth,
//! run after stage 2.2 (where the jet and the smooth regions exist) and before
//! 2.3 (whose decomposition it deliberately does not use — a fragile invariant
//! must be visible on models whose decomposition fails outright).
//!
//! Five families, matching the section's five headings:
//!
//! - **F (Fitting Requirement)** — with the reference's *own exact* surface
//!   parameters, how far off it do the facets' vertices, edge midpoints, facet
//!   centroids and facet interiors actually land, and how large is the
//!   within-surface dihedral? This asks whether the requirement is satisfiable
//!   at all before asking whether any algorithm satisfies it. It also measures
//!   the two probe choices the section makes (edge centres for developables,
//!   facet centres for spheres) against the true whole-facet maximum.
//! - **A (angular separability)** — the largest dihedral *inside* a reference
//!   face against the smallest dihedral *across* its boundary. Where the first
//!   exceeds the second, no angular tolerance can separate that pair.
//! - **I (Surface Identification)** — fit all four analytic types to the true
//!   facet set of each reference face and report each one's worst vertex
//!   residual, in millimetres and in float32 quanta. This is the distribution
//!   behind the section's open question about a cone that fits a cylinder.
//! - **S (Seed Independence)** — run the production grow (`compete::grow`) from
//!   many different seeds inside one reference face and count how many distinct
//!   face sets come back, how many equal the truth, and the worst Jaccard.
//!   Also runs the full multi-type competition (`compete::compete_from_seed`)
//!   per seed, which is what the pipeline actually does.
//! - **G (Growth Independence)** — the flood re-run with the neighbour
//!   iteration reversed (a direct test of traversal-order independence), the
//!   refit fixpoint started from the whole true face set instead of the seed
//!   patch (a direct test of start independence), and how often the refit loop
//!   hits its round cap without converging.
//!
//! Output is one TSV row per reference face on **stdout** (prefix `INV`) plus a
//! per-model summary on stderr, so a shell harness can concatenate models.

use std::collections::{HashMap, HashSet};

use crate::config::Config;
use crate::reference::{FaceLabel, RefSurface, ReferenceModel};
use crate::stage1::{float32_position_bound_mm, ConnectedMesh};

use super::compete;
use super::conical::ConicalSurface;
use super::edges::{EdgeClass, EdgeClassification, JetField};
use super::planar::PlanarSurface;
use super::spherical::SphericalSurface;
use super::surface::SurfaceType;
use super::toroidal::ToroidalSurface;
use super::{
    face_area, face_centroid, MIN_CYLINDER_FACES, MIN_SPHERE_FACES, MIN_TORUS_FACES,
};

/// How many seeds to try per reference face in the seed-independence test.
/// Seeds are taken at a fixed stride over the face's sorted facet list, so the
/// sample is deterministic and spread over the surface rather than clustered.
/// Override with `BREPPER_INV_SEEDS`; the per-seed competition costs about as
/// much as one stage-2.3 seed, which on a 10k-face smooth region is ~0.4 s.
const DEFAULT_SEEDS_PER_FACE: usize = 12;

fn seeds_per_face() -> usize {
    std::env::var("BREPPER_INV_SEEDS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(DEFAULT_SEEDS_PER_FACE)
}

/// `BREPPER_INV_NO_COMPETE=1` drops the four-type competition per seed, which
/// is by far the most expensive measurement, keeping the single-type grow, the
/// traversal-order test and the oracle start.
fn compete_enabled() -> bool {
    std::env::var("BREPPER_INV_NO_COMPETE").is_err()
}

/// Reference faces with fewer facets than this are reported but not
/// seed-tested: with 3 facets there is nothing for a seed to be independent of.
const MIN_FACETS_FOR_SEED_TEST: usize = 6;

/// Barycentric sampling order for the whole-facet maximum. Order n gives
/// (n-1)(n-2)/2 strictly interior points per triangle — 21 at n = 8.
const INTERIOR_ORDER: usize = 8;

// ---------------------------------------------------------------------------
// Exact distances to a reference surface
// ---------------------------------------------------------------------------

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n > 0.0 {
        [a[0] / n, a[1] / n, a[2] / n]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Perpendicular distance from `p` to a reference surface, or `None` for the
/// freeform types (no closed form; measuring those needs OCCT and is out of
/// scope for a per-sample-point inner loop).
fn ref_distance(s: &RefSurface, p: [f64; 3]) -> Option<f64> {
    Some(match *s {
        RefSurface::Plane { point, normal } => dot(sub(p, point), unit(normal)).abs(),
        RefSurface::Cylinder {
            axis_pos,
            axis_dir,
            radius,
        } => {
            let d = unit(axis_dir);
            let v = sub(p, axis_pos);
            let h = dot(v, d);
            let rho = norm(sub(v, [d[0] * h, d[1] * h, d[2] * h]));
            (rho - radius).abs()
        }
        RefSurface::Cone {
            apex,
            axis_dir,
            half_angle,
        } => {
            // Distance to the double cone's nearer nappe, measured along the
            // surface normal: rho*cos(a) - h*sin(a) with h along the opening
            // direction. OCCT's semi-angle is signed with the axis direction.
            let d = unit(axis_dir);
            let v = sub(p, apex);
            let h = dot(v, d);
            let rho = norm(sub(v, [d[0] * h, d[1] * h, d[2] * h]));
            (rho * half_angle.cos() - h * half_angle.sin()).abs()
        }
        RefSurface::Sphere { centre, radius } => (norm(sub(p, centre)) - radius).abs(),
        RefSurface::Torus {
            centre,
            axis_dir,
            major,
            minor,
        } => {
            let d = unit(axis_dir);
            let v = sub(p, centre);
            let z = dot(v, d);
            let rho = norm(sub(v, [d[0] * z, d[1] * z, d[2] * z]));
            (((rho - major).powi(2) + z * z).sqrt() - minor).abs()
        }
        RefSurface::Other { .. } => return None,
    })
}

// ---------------------------------------------------------------------------
// Facet geometry helpers
// ---------------------------------------------------------------------------

fn vertex_pos(mesh: &ConnectedMesh, vi: usize) -> [f64; 3] {
    let v = &mesh.vertices[vi];
    [v.x, v.y, v.z]
}

fn facet_points(mesh: &ConnectedMesh, fi: usize) -> Vec<[f64; 3]> {
    let f = &mesh.faces[fi];
    (0..f.vertex_count as usize)
        .map(|k| vertex_pos(mesh, f.vertex_indices[k]))
        .collect()
}

/// Unit facet normal — the mesh's own stored normal where there is one, so
/// this measures the *same* quantity `edges::classify_edges` gates the sharp
/// partition on. Recomputed from the first three vertices only as a fallback.
fn facet_normal(mesh: &ConnectedMesh, fi: usize) -> [f64; 3] {
    if let Some(n) = mesh.faces[fi].normal {
        return unit(n);
    }
    let p = facet_points(mesh, fi);
    if p.len() < 3 {
        return [0.0, 0.0, 1.0];
    }
    let a = sub(p[1], p[0]);
    let b = sub(p[2], p[0]);
    unit([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ])
}

/// Dihedral angle between two adjacent facets, in degrees.
fn dihedral_deg(mesh: &ConnectedMesh, a: usize, b: usize) -> f64 {
    let (na, nb) = (facet_normal(mesh, a), facet_normal(mesh, b));
    // Angle via the cross product, which stays linear in the error near 0°
    // (see DEVELOPMENT_PLAN: acos amplifies rounding near an argument of 1).
    let c = [
        na[1] * nb[2] - na[2] * nb[1],
        na[2] * nb[0] - na[0] * nb[2],
        na[0] * nb[1] - na[1] * nb[0],
    ];
    norm(c).atan2(dot(na, nb)).to_degrees().abs()
}

/// The three probes the Invariants section names, plus the true whole-facet
/// maximum from a dense interior sample: (vertex, edge midpoint, centroid,
/// interior).
fn facet_probes(
    mesh: &ConnectedMesh,
    fi: usize,
    s: &RefSurface,
) -> Option<(f64, f64, f64, f64, f64)> {
    let p = facet_points(mesh, fi);
    let n = p.len();
    if n < 3 {
        return None;
    }
    let mut v_max = 0.0_f64;
    for &q in &p {
        v_max = v_max.max(ref_distance(s, q)?);
    }
    // Boundary edge midpoints, PLUS the diagonals of a fused quad. The
    // diagonals matter and the section's text does not mention them: the
    // whole-facet maximum is bounded by 4/3 of the largest midpoint over the
    // edges of a TRIANGULATION, and a quad's triangulation needs a diagonal.
    // Left out, a square quad on a sphere reads exactly half the truth (its
    // diagonal midpoint IS the quad's centre, where the sag is 2x the side
    // midpoints') — measured at 1.999 on `sphere_25_fine` before this was
    // added.
    let mut e_max = 0.0_f64;
    let mid = |a: [f64; 3], b: [f64; 3]| -> Option<f64> {
        ref_distance(s, [
            0.5 * (a[0] + b[0]),
            0.5 * (a[1] + b[1]),
            0.5 * (a[2] + b[2]),
        ])
    };
    for k in 0..n {
        e_max = e_max.max(mid(p[k], p[(k + 1) % n])?);
    }
    let mut d_max = 0.0_f64;
    if n == 4 {
        d_max = d_max.max(mid(p[0], p[2])?).max(mid(p[1], p[3])?);
    }
    let e_max_with_diag = e_max.max(d_max);
    let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
    let c_max = ref_distance(s, c)?;

    // Interior: barycentric grid over the fan triangulation.
    let mut i_max = 0.0_f64;
    let tris: Vec<[usize; 3]> = (1..n - 1).map(|k| [0, k, k + 1]).collect();
    let m = INTERIOR_ORDER;
    for t in tris {
        for i in 1..m {
            for j in 1..(m - i) {
                let k = m - i - j;
                if k == 0 {
                    continue;
                }
                let (a, b, cc) = (i as f64 / m as f64, j as f64 / m as f64, k as f64 / m as f64);
                let q = [
                    a * p[t[0]][0] + b * p[t[1]][0] + cc * p[t[2]][0],
                    a * p[t[0]][1] + b * p[t[1]][1] + cc * p[t[2]][1],
                    a * p[t[0]][2] + b * p[t[1]][2] + cc * p[t[2]][2],
                ];
                i_max = i_max.max(ref_distance(s, q)?);
            }
        }
    }
    Some((v_max, e_max, c_max, i_max, e_max_with_diag))
}

// ---------------------------------------------------------------------------
// Seed-rule probe: are the section's stated seed requirements the right shape?
// ---------------------------------------------------------------------------

/// How many 6-vertex chains to sample per developable reference face.
const MAX_SEED_PAIRS: usize = 200;

/// The two readings of "no more than three of which are coplanar", plus the
/// closed-form axis error the seed actually yields.
struct SeedRuleStats {
    pairs: usize,
    /// Union of the pair has exactly 6 distinct vertices (the section's rule).
    six_verts: usize,
    /// All 6 lie in a plane — the genuinely degenerate case (two rulings).
    all_coplanar: usize,
    /// SOME four of them lie in a plane — the literal reading of the rule,
    /// which a fused quad satisfies by construction.
    four_coplanar: usize,
    /// Axis error (degrees) of `initialize::cylinder` on the pair alone.
    ax_err: Vec<f64>,
    /// Same, restricted to chains the section's rule would admit under the
    /// "not all six coplanar" reading.
    ax_err_admitted: Vec<f64>,
    /// Axis error of the Phase-1.1 initializer over the WHOLE face — the
    /// selection rule, for contrast with a bare minimal seed.
    selected_err: f64,
}

/// Largest and smallest spread of a point set about its centroid, and the
/// out-of-plane spread — enough to answer "are these coplanar?" scale-free.
fn planarity(points: &[[f64; 3]]) -> f64 {
    let n = points.len();
    if n < 4 {
        return 0.0;
    }
    let inv = 1.0 / n as f64;
    let mut c = [0.0_f64; 3];
    for p in points {
        for k in 0..3 {
            c[k] += p[k] * inv;
        }
    }
    // 3x3 covariance; the smallest eigenvalue's sqrt is the out-of-plane RMS.
    let mut m = [[0.0_f64; 3]; 3];
    let mut extent = 0.0_f64;
    for p in points {
        let d = sub(*p, c);
        extent = extent.max(norm(d));
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += d[i] * d[j] * inv;
            }
        }
    }
    if extent <= 0.0 {
        return 0.0;
    }
    // Smallest eigenvalue of a symmetric 3x3 by the closed form used elsewhere
    // in stage 2 would be overkill here: power-iterate the inverse via the
    // characteristic polynomial's smallest root, computed by bisection on the
    // determinant, is also overkill. Use the analytic symmetric eigenvalues.
    let lam = sym3_eigenvalues(&m);
    (lam[0].max(0.0)).sqrt() / extent
}

/// Ascending eigenvalues of a symmetric 3x3 matrix (closed form).
#[allow(clippy::needless_range_loop)] // matrix code: i/j index several arrays at once
fn sym3_eigenvalues(m: &[[f64; 3]; 3]) -> [f64; 3] {
    let p1 = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];
    let tr = m[0][0] + m[1][1] + m[2][2];
    if p1 <= 0.0 {
        let mut e = [m[0][0], m[1][1], m[2][2]];
        e.sort_by(|a, b| a.partial_cmp(b).unwrap());
        return e;
    }
    let q = tr / 3.0;
    let p2 = (m[0][0] - q).powi(2) + (m[1][1] - q).powi(2) + (m[2][2] - q).powi(2) + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();
    let mut b = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            b[i][j] = (m[i][j] - if i == j { q } else { 0.0 }) / p;
        }
    }
    let det = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
        - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
        + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
    let r = (det / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let e1 = q + 2.0 * p * phi.cos();
    let e3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).cos();
    let e2 = tr - e1 - e3;
    let mut e = [e1, e2, e3];
    e.sort_by(|a, b| a.partial_cmp(b).unwrap());
    e
}

/// Do some four of `points` lie in a plane to within `rel` of their extent?
fn any_four_coplanar(points: &[[f64; 3]], rel: f64) -> bool {
    let n = points.len();
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                for l in (k + 1)..n {
                    let quad = [points[i], points[j], points[k], points[l]];
                    if planarity(&quad) <= rel {
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// Enumerate facet **chains** (paths in the dual graph, the section's
/// "non-contiguous facet chain") of up to `MAX_CHAIN` facets inside one
/// reference face, and keep those with exactly 6 distinct vertices — the
/// section's stated cylinder/cone seed.
const MAX_CHAIN: usize = 4;

fn chains_with_six_vertices(
    mesh: &ConnectedMesh,
    faces: &[usize],
    member: &HashSet<usize>,
    cap: usize,
) -> (Vec<Vec<usize>>, usize) {
    let mut out: Vec<Vec<usize>> = Vec::new();
    // Distinct facet SETS considered. The DFS walks each undirected chain once
    // per direction, and the returned `out` is deduplicated by face set, so
    // counting raw paths here would put an ordered count over a deduplicated
    // one — a ratio formed that way reads about half the truth.
    let mut considered: HashSet<Vec<usize>> = HashSet::new();
    for &start in faces {
        if out.len() >= cap {
            break;
        }
        // Depth-first enumeration of simple paths from `start`.
        let mut stack: Vec<Vec<usize>> = vec![vec![start]];
        while let Some(path) = stack.pop() {
            if out.len() >= cap {
                break;
            }
            if path.len() >= 2 {
                let mut key = path.clone();
                key.sort_unstable();
                let fresh = considered.insert(key);
                if fresh && vertices_of(mesh, &path).len() == 6 {
                    let mut p = path.clone();
                    p.sort_unstable();
                    out.push(p);
                }
            }
            if path.len() >= MAX_CHAIN {
                continue;
            }
            let last = *path.last().unwrap();
            let vc = mesh.faces[last].vertex_count as usize;
            for &ni in &mesh.faces[last].neighbors[..vc] {
                if ni < 0 {
                    continue;
                }
                let nj = ni as usize;
                if !member.contains(&nj) || path.contains(&nj) {
                    continue;
                }
                let mut next = path.clone();
                next.push(nj);
                stack.push(next);
            }
        }
    }
    out.sort();
    out.dedup();
    (out, considered.len())
}

// ---------------------------------------------------------------------------
// The clean grow operator, and its basin of attraction
// ---------------------------------------------------------------------------

/// Refit-and-reflood rounds before the clean operator gives up.
const MAX_CLEAN_ROUNDS: usize = 12;

/// What happened when a seed SET was grown under the clean operator.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Basin {
    /// Converged to exactly the reference face. The seed is in the basin.
    Exact,
    /// Converged to a strict subset — the fit stalled inside its own surface.
    Short,
    /// Converged to a strict superset — the fit leaked past the true boundary.
    Leak,
    /// Converged to something that is neither.
    Other,
    /// No fit of this type exists on the seed set (or on some later support).
    NoFit,
    /// The fit of the seed set does not validate any facet of the seed set.
    SeedRejected,
    /// Still moving at the round cap.
    NoConverge,
}

impl Basin {
    fn code(self) -> u8 {
        match self {
            Basin::Exact => 0,
            Basin::Short => 1,
            Basin::Leak => 2,
            Basin::Other => 3,
            Basin::NoFit => 4,
            Basin::SeedRejected => 5,
            Basin::NoConverge => 6,
        }
    }
}

/// **The clean operator.** `F₀ = S`; `F_{k+1}` = the union of the validated
/// components containing the seed set, under a fresh least-squares fit of
/// `F_k`; iterate to a fixpoint.
///
/// No `seed_patch`, no `condition_patch`, no trimming, no jet — the seed set is
/// the entire initial fit support and validation is the only thing that decides
/// membership. This is the baseline the Invariants section is trying to
/// establish, so the question "which seed sets are valid" becomes "which seed
/// sets lie in the basin of attraction of their own surface", which is what
/// this measures.
///
/// Growth is confined to the smooth region, as the pipeline does: the sharp
/// partition is the one thing measured exact, so leaving it in place isolates
/// the seeding question from the partitioning one.
fn grow_clean<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    seed: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
) -> (Basin, Vec<usize>, usize, bool, Option<S::Fit>) {
    let vtol = config.vertex_tolerance_mm;
    let stol = config.surface_tolerance_mm;
    let mut cur: Vec<usize> = seed.to_vec();
    cur.sort_unstable();
    let mut self_ok = false;
    for round in 0..MAX_CLEAN_ROUNDS {
        let verts = vertices_of(mesh, &cur);
        let Some(mut fit) = S::fit_surface(mesh, &cur, &verts, 0) else {
            return (Basin::NoFit, cur, round, self_ok, None);
        };
        S::determine_convexity(mesh, &cur, &mut fit);
        if round == 0 {
            // Does the seed set's own fit explain the whole seed set? A
            // ground-truth-free predicate, so it could gate admission in
            // production.
            self_ok = seed
                .iter()
                .all(|&fi| compete::face_validates::<S>(mesh, fi, &fit, centroids, vtol, stol));
        }
        let mut next: HashSet<usize> = HashSet::new();
        for &st in seed {
            if next.contains(&st) {
                continue;
            }
            let comp = compete::flood::<S>(mesh, in_region, st, &fit, centroids, vtol, stol);
            next.extend(comp);
        }
        if next.is_empty() {
            return (Basin::SeedRejected, cur, round, self_ok, Some(fit));
        }
        let mut nv: Vec<usize> = next.into_iter().collect();
        nv.sort_unstable();
        if nv == cur {
            return (Basin::Exact, cur, round, self_ok, Some(fit)); // caller re-classifies
        }
        cur = nv;
    }
    (Basin::NoConverge, cur, MAX_CLEAN_ROUNDS, self_ok, None)
}

/// Classify a converged result against the reference face.
fn classify(result: &[usize], truth: &HashSet<usize>) -> Basin {
    let inside = result.iter().filter(|f| truth.contains(f)).count();
    if inside == result.len() && result.len() == truth.len() {
        Basin::Exact
    } else if inside == result.len() {
        Basin::Short
    } else if inside == truth.len() {
        Basin::Leak
    } else {
        Basin::Other
    }
}

// ---------------------------------------------------------------------------
// The DEVELOPMENT_PLAN "Surface Seed Requirements", as literal generators
// ---------------------------------------------------------------------------

/// Vertex -> the facets of one surface incident on it.
fn vertex_star_map(mesh: &ConnectedMesh, faces: &[usize]) -> HashMap<usize, Vec<usize>> {
    let mut m: HashMap<usize, Vec<usize>> = HashMap::new();
    for &fi in faces {
        let vc = mesh.faces[fi].vertex_count as usize;
        for &vi in &mesh.faces[fi].vertex_indices[..vc] {
            m.entry(vi).or_default().push(fi);
        }
    }
    for v in m.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    m
}

/// Is every facet of the mesh incident on `vi` inside `member`? A complete star
/// means the vertex is interior to the surface, so the ring closes.
fn star_is_complete(mesh: &ConnectedMesh, vi: usize, star: &[usize], member: &HashSet<usize>) -> bool {
    // Walk out from the star: any mesh facet sharing `vi` with a star member
    // that is not itself in `member` means the ring is cut by the boundary.
    for &fi in star {
        let vc = mesh.faces[fi].vertex_count as usize;
        for &ni in &mesh.faces[fi].neighbors[..vc] {
            if ni < 0 {
                continue;
            }
            let nj = ni as usize;
            let nvc = mesh.faces[nj].vertex_count as usize;
            if mesh.faces[nj].vertex_indices[..nvc].contains(&vi) && !member.contains(&nj) {
                return false;
            }
        }
    }
    true
}

/// Connected facet subsets of size 2..=`max_n` containing `anchor`, confined to
/// `member`. The section calls these "chains"; a subset that happens not to be
/// a simple path is reported as such rather than silently admitted.
fn connected_subsets_containing(
    mesh: &ConnectedMesh,
    member: &HashSet<usize>,
    anchor: usize,
    max_n: usize,
    cap: usize,
) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut frontier: Vec<Vec<usize>> = vec![vec![anchor]];
    for _ in 1..max_n {
        let mut next: Vec<Vec<usize>> = Vec::new();
        for set in &frontier {
            for &fi in set {
                let vc = mesh.faces[fi].vertex_count as usize;
                for &ni in &mesh.faces[fi].neighbors[..vc] {
                    if ni < 0 {
                        continue;
                    }
                    let nj = ni as usize;
                    if !member.contains(&nj) || set.contains(&nj) {
                        continue;
                    }
                    let mut grown = set.clone();
                    grown.push(nj);
                    grown.sort_unstable();
                    if seen.insert(grown.clone()) {
                        next.push(grown.clone());
                        out.push(grown);
                        if out.len() >= cap {
                            return out;
                        }
                    }
                }
            }
        }
        if next.is_empty() {
            break;
        }
        frontier = next;
    }
    out
}

/// Is `set` a simple path in the dual graph (every member adjacent to at most
/// two others within the set, and connected)?
fn is_chain(mesh: &ConnectedMesh, set: &[usize]) -> bool {
    if !is_connected(mesh, set) {
        return false;
    }
    let m: HashSet<usize> = set.iter().copied().collect();
    let mut ends = 0;
    for &fi in set {
        let vc = mesh.faces[fi].vertex_count as usize;
        let deg = mesh.faces[fi].neighbors[..vc]
            .iter()
            .filter(|&&ni| ni >= 0 && m.contains(&(ni as usize)))
            .count();
        if deg > 2 {
            return false;
        }
        if deg <= 1 {
            ends += 1;
        }
    }
    set.len() == 1 || ends == 2
}

/// The seed sets the section's **Surface Seed Requirements** admit that contain
/// `anchor`, all confined to one true surface.
///
/// This is the section read literally, with one addition the section is missing
/// (a torus rule, see below) and one variant it needs measured both ways (the
/// coplanarity clause).
///
/// - **Plane** — the facet itself.
/// - **Cylinder** — a connected facet chain of at most `max_chain` facets with
///   **exactly 6 distinct vertices**. Measured, exactly three shapes reach 6:
///   two quads, four triangles, or two triangles and a quad.
/// - **Cone** — the same, plus three triangular facets sharing one vertex (the
///   candidate apex).
/// - **Sphere** — the anchor plus every facet sharing an edge with it, and the
///   same ring taken about each of its neighbours (both contain the anchor).
/// - **Torus** — the **vertex star**: every facet incident on one vertex. This
///   is the rule the section omits, and the reason it has to be a star rather
///   than a chain is identifiability, not taste: a torus's normal lines cross
///   the axis at a point that slides with the *tube* angle, so a seed confined
///   to one tube angle cannot pin the axis — which is exactly what a chain of
///   two quads along the major direction is. A star wraps a grid node, so it
///   spans both parametric directions by construction.
#[allow(clippy::too_many_arguments)]
fn seed_family(
    mesh: &ConnectedMesh,
    member: &HashSet<usize>,
    vstar: &HashMap<usize, Vec<usize>>,
    kind: Kind,
    anchor: usize,
    max_chain: usize,
    require_coplanar_clause: bool,
    spec_only: bool,
    cap: usize,
) -> Vec<(&'static str, Vec<usize>)> {
    let mut out: Vec<(&'static str, Vec<usize>)> = Vec::new();
    let copl_rel = 1.0e-6;
    match kind {
        Kind::Plane => out.push(("facet", vec![anchor])),
        Kind::Cylinder | Kind::Cone => {
            for set in connected_subsets_containing(mesh, member, anchor, max_chain, cap * 8) {
                let vs = vertices_of(mesh, &set);
                if vs.len() != 6 || !is_chain(mesh, &set) {
                    continue;
                }
                if require_coplanar_clause {
                    let pts: Vec<[f64; 3]> = vs.iter().map(|&v| vertex_pos(mesh, v)).collect();
                    if any_four_coplanar(&pts, copl_rel) {
                        continue;
                    }
                }
                out.push(("chain6", set));
                if out.len() >= cap {
                    break;
                }
            }
            if kind == Kind::Cone {
                // Three triangles sharing a vertex — the apex alternative.
                let vc = mesh.faces[anchor].vertex_count as usize;
                for &vi in &mesh.faces[anchor].vertex_indices[..vc] {
                    let star = &vstar[&vi];
                    let tris: Vec<usize> = star
                        .iter()
                        .copied()
                        .filter(|&f| mesh.faces[f].vertex_count == 3)
                        .collect();
                    if tris.len() >= 3 && tris.contains(&anchor) {
                        let mut pick: Vec<usize> =
                            tris.iter().copied().filter(|&f| f != anchor).take(2).collect();
                        pick.push(anchor);
                        pick.sort_unstable();
                        out.push(("apex3", pick));
                    }
                }
            }
        }
        Kind::Sphere => {
            let ring = |f: usize| -> Vec<usize> {
                let mut r = vec![f];
                let vc = mesh.faces[f].vertex_count as usize;
                for &ni in &mesh.faces[f].neighbors[..vc] {
                    if ni >= 0 && member.contains(&(ni as usize)) {
                        r.push(ni as usize);
                    }
                }
                r.sort_unstable();
                r
            };
            out.push(("edge-ring", ring(anchor)));
            let vc = mesh.faces[anchor].vertex_count as usize;
            for &ni in &mesh.faces[anchor].neighbors[..vc] {
                if ni >= 0 && member.contains(&(ni as usize)) {
                    let r = ring(ni as usize);
                    if r.contains(&anchor) {
                        out.push(("edge-ring", r));
                    }
                }
                if out.len() >= cap {
                    break;
                }
            }
        }
        Kind::Torus => {
            let vc = mesh.faces[anchor].vertex_count as usize;
            // Single vertex star — one closed ring of facets about one vertex.
            for &vi in &mesh.faces[anchor].vertex_indices[..vc] {
                let star = vstar[&vi].clone();
                if star.len() < 4 || !star.contains(&anchor) {
                    continue;
                }
                let name = if star_is_complete(mesh, vi, &star, member) {
                    "star"
                } else {
                    "star-cut"
                };
                out.push((name, star));
            }
            // The union of the stars of the anchor's OWN vertices — its full
            // vertex neighbourhood. A 7-DOF torus fitted to a single star of
            // four quads has nine vertices against seven unknowns, which is
            // determining only just; this is the next widening that is still
            // "rings of facets about a vertex".
            let mut wide: Vec<usize> = Vec::new();
            let mut all_complete = true;
            for &vi in &mesh.faces[anchor].vertex_indices[..vc] {
                if let Some(st) = vstar.get(&vi) {
                    if !star_is_complete(mesh, vi, st, member) {
                        all_complete = false;
                    }
                    wide.extend(st.iter().copied());
                }
            }
            wide.sort_unstable();
            wide.dedup();
            if wide.len() >= 4 && wide.contains(&anchor) {
                out.push((if all_complete { "star2" } else { "star2-cut" }, wide));
            }
            out.truncate(cap.max(4));
        }
        Kind::Other => {}
    }
    // Fallback ladder, tried only after the section's own patterns: plain
    // connected balls of increasing size about the anchor. Because the family is
    // tried smallest-first and short-circuits on the first success, the winning
    // name reports the LEAST support that actually works — the section's own
    // rule where that suffices, and otherwise the size it takes.
    if kind != Kind::Other && !spec_only {
        for &k in &[8usize, 12, 16, 24, 32, 48] {
            let mut ball = vec![anchor];
            let mut seen: HashSet<usize> = HashSet::from([anchor]);
            let mut frontier = vec![anchor];
            while ball.len() < k && !frontier.is_empty() {
                let mut ring = Vec::new();
                for &fi in &frontier {
                    let fvc = mesh.faces[fi].vertex_count as usize;
                    for &ni in &mesh.faces[fi].neighbors[..fvc] {
                        if ni < 0 || ball.len() >= k {
                            continue;
                        }
                        let nj = ni as usize;
                        if member.contains(&nj) && seen.insert(nj) {
                            ball.push(nj);
                            ring.push(nj);
                        }
                    }
                }
                if ring.is_empty() {
                    break;
                }
                frontier = ring;
            }
            if ball.len() < k {
                break; // the surface itself is smaller than this rung
            }
            ball.sort_unstable();
            out.push((
                match k {
                    8 => "ball8",
                    12 => "ball12",
                    16 => "ball16",
                    24 => "ball24",
                    32 => "ball32",
                    _ => "ball48",
                },
                ball,
            ));
        }
    }

    // The section's own patterns first, then the ladder by size, so the winner
    // is always the least admissible support.
    let rank = |n: &str| match n {
        "facet" | "chain6" | "apex3" | "edge-ring" | "star" | "star-cut" | "star2"
        | "star2-cut" => 0,
        _ => 1,
    };
    out.sort_by(|a, b| {
        rank(a.0)
            .cmp(&rank(b.0))
            .then_with(|| a.1.len().cmp(&b.1.len()))
            .then_with(|| a.1.cmp(&b.1))
    });
    out.dedup_by(|a, b| a.1 == b.1);
    out
}

// ---------------------------------------------------------------------------
// `--viz-seeds`: watch seed sets grow, and see why the failures fail
// ---------------------------------------------------------------------------

/// Seed green; covered blue; overshoot red; missed yellow; blocked orange.
const C_SEED: [f32; 4] = [0.15, 0.85, 0.25, 1.0];
const C_COVERED: [f32; 4] = [0.20, 0.45, 0.95, 1.0];
const C_OVERSHOOT: [f32; 4] = [0.95, 0.15, 0.15, 1.0];
const C_MISSED: [f32; 4] = [0.95, 0.90, 0.15, 1.0];
const C_BLOCKED: [f32; 4] = [1.00, 0.55, 0.05, 1.0];
const C_ANCHOR: [f32; 4] = [1.00, 1.00, 1.00, 1.0];

/// Why one facet passed or failed `face_validates`, clause by clause.
///
/// `face_validates` is a conjunction of three tests — vertex distance, centroid
/// distance, and orientation coherence — and a diagnosis that reports only the
/// first is actively misleading: a facet can sit exactly on the surface (0.00x
/// vertex tolerance) and still be refused for the other two. Both of those
/// happen in the corpus.
#[derive(Clone, Copy)]
struct FacetVerdict {
    /// Worst vertex distance, in units of `--vertex-tolerance`.
    vratio: f64,
    /// Centroid distance, in units of `--surface-tolerance`.
    cratio: f64,
    /// Does the facet normal agree in sign with the surface normal?
    orient: bool,
}

impl FacetVerdict {
    fn ok(&self) -> bool {
        self.orient && self.vratio <= 1.0 && self.cratio <= 1.0
    }
    /// The clause that refused this facet, and by how much.
    fn cause(&self) -> String {
        if !self.orient {
            "ORIENTATION (facet normal opposes the surface)".to_string()
        } else if self.vratio > 1.0 && self.cratio > 1.0 {
            format!(
                "VERTEX {:.1}x vtol and CENTROID {:.1}x stol",
                self.vratio, self.cratio
            )
        } else if self.vratio > 1.0 {
            format!("VERTEX {:.1}x vtol", self.vratio)
        } else if self.cratio > 1.0 {
            format!("CENTROID {:.1}x stol", self.cratio)
        } else {
            "nothing (it passes)".to_string()
        }
    }
}

fn facet_verdict<S: SurfaceType>(
    mesh: &ConnectedMesh,
    fi: usize,
    fit: &S::Fit,
    centroids: &[[f64; 3]],
    vtol: f64,
    stol: f64,
) -> FacetVerdict {
    let vratio = facet_points(mesh, fi)
        .iter()
        .map(|p| S::vertex_distance(fit, p))
        .fold(0.0_f64, f64::max)
        / vtol.max(1e-30);
    FacetVerdict {
        vratio,
        cratio: S::centroid_distance(fit, &centroids[fi]) / stol.max(1e-30),
        orient: S::convexity_check(fit, mesh, fi),
    }
}

/// Summarise a set of facet verdicts into one dominant cause.
fn dominant_cause(vs: &[(usize, FacetVerdict)]) -> String {
    let bad: Vec<&(usize, FacetVerdict)> = vs.iter().filter(|(_, v)| !v.ok()).collect();
    if bad.is_empty() {
        return "nothing — every facet passes".to_string();
    }
    let (orient, vtx, cent) = (
        bad.iter().filter(|(_, v)| !v.orient).count(),
        bad.iter().filter(|(_, v)| v.orient && v.vratio > 1.0).count(),
        bad.iter()
            .filter(|(_, v)| v.orient && v.vratio <= 1.0 && v.cratio > 1.0)
            .count(),
    );
    let mut parts = Vec::new();
    if orient > 0 {
        parts.push(format!("{orient} on ORIENTATION"));
    }
    if vtx > 0 {
        let worst = bad
            .iter()
            .filter(|(_, v)| v.orient && v.vratio > 1.0)
            .map(|(_, v)| v.vratio)
            .fold(0.0_f64, f64::max);
        let closest = bad
            .iter()
            .filter(|(_, v)| v.orient && v.vratio > 1.0)
            .map(|(_, v)| v.vratio)
            .fold(f64::INFINITY, f64::min);
        parts.push(format!(
            "{vtx} on VERTEX ({closest:.1}x-{worst:.1}x vtol)"
        ));
    }
    if cent > 0 {
        let worst = bad
            .iter()
            .filter(|(_, v)| v.orient && v.vratio <= 1.0 && v.cratio > 1.0)
            .map(|(_, v)| v.cratio)
            .fold(0.0_f64, f64::max);
        parts.push(format!("{cent} on CENTROID (up to {worst:.1}x stol)"));
    }
    parts.join(", ")
}

/// One seed set to show, with everything needed to explain its outcome.
struct SeedCase {
    shape: &'static str,
    seed: Vec<usize>,
    verdict: Basin,
    result: Vec<usize>,
    rounds: usize,
    self_ok: bool,
    /// Truth facets the growth never reached.
    missed: Vec<usize>,
    /// Result facets outside the truth.
    overshoot: Vec<usize>,
    /// Missed truth facets that touch the result — the ones that actually
    /// stopped it — with the clause-by-clause verdict for each.
    blocked: Vec<(usize, FacetVerdict)>,
    /// Missed truth facets with NO neighbour in the result: not refused by
    /// validation at all, simply not reachable from the seed without leaving
    /// the surface. A connectivity fact about the reference face, not a
    /// seeding failure.
    unreachable: Vec<usize>,
    /// Missed truth facets that lie OUTSIDE the smooth region the seed is in.
    /// The flood is confined to one region, so these are unreachable no matter
    /// what the seed is or how well it fits — stage 2.1's sharp-edge partition
    /// cut this true surface in two before any fitting happened.
    outside_region: Vec<usize>,
    /// Is the reference face even connected through its own facets?
    truth_connected: bool,
    /// Clause-by-clause verdict of every seed facet against the seed set's OWN
    /// first fit. This is the admission quantity.
    seed_verdicts: Vec<(usize, FacetVerdict)>,
}

impl SeedCase {
    fn covers(&self) -> bool {
        self.missed.is_empty() && !matches!(self.verdict, Basin::NoFit | Basin::SeedRejected)
    }
    /// A one-line diagnosis, in the vocabulary of what went wrong.
    fn why(&self, n_truth: usize) -> String {
        match self.verdict {
            Basin::NoFit => format!(
                "NO FIT — the fitter declines a {}-facet support outright (under-determined); \
                 nothing to validate or grow",
                self.seed.len()
            ),
            Basin::SeedRejected if self.rounds == 0 => format!(
                "SEED REJECTED AT ROUND 0 — the seed set's own fit refuses its own facets, so \
                 nothing floods: {}",
                dominant_cause(&self.seed_verdicts)
            ),
            // Rejected on a LATER round is a different failure and used to be
            // reported as the same one: the seed's own fit was fine, growth
            // happened, and then the refit on the grown support drifted far
            // enough that it no longer validates the seed it came from. The
            // iteration diverged rather than the seed being bad.
            Basin::SeedRejected => format!(
                "REFIT DRIFTED OFF ITS OWN SEED after {} round(s) — the seed's own fit was fine; \
                 the fit of the grown support no longer validates the seed",
                self.rounds
            ),
            Basin::Short => {
                let refused = self.blocked.iter().filter(|(_, v)| !v.ok()).count();
                if refused == 0 && !self.outside_region.is_empty() {
                    format!(
                        "CUT BY THE PARTITION — reached {}/{} facets; the {} missing one(s) are \
                         in a DIFFERENT smooth region, so no seed can reach them. Stage 2.1 called \
                         an edge inside this true surface sharp. Not a seeding failure.",
                        self.result.len(),
                        n_truth,
                        self.outside_region.len(),
                    )
                } else if refused == 0 && !self.unreachable.is_empty() {
                    format!(
                        "STALLED — reached {}/{} facets, and the {} missing one(s) are NOT REFUSED: \
                         they have no mesh neighbour in the grown set{}. The reference face is \
                         {}connected through its own facets.",
                        self.result.len(),
                        n_truth,
                        self.unreachable.len(),
                        if self.truth_connected {
                            ", so growth would have to leave the surface and come back"
                        } else {
                            ""
                        },
                        if self.truth_connected { "" } else { "NOT " },
                    )
                } else {
                    format!(
                        "STALLED — reached {}/{} facets; {} true facet(s) on the frontier were \
                         refused: {}{}",
                        self.result.len(),
                        n_truth,
                        refused,
                        dominant_cause(&self.blocked),
                        match (self.unreachable.len(), self.outside_region.len()) {
                            (0, 0) => String::new(),
                            (u, 0) => format!(" (+{u} unreachable)"),
                            (0, o) => format!(" (+{o} cut into another region)"),
                            (u, o) => format!(" (+{u} unreachable, +{o} in another region)"),
                        },
                    )
                }
            }
            Basin::Leak => format!(
                "COVERS + OVERSHOOT — {}/{} true facets plus {} beyond the true boundary",
                self.result.len() - self.overshoot.len(),
                n_truth,
                self.overshoot.len()
            ),
            Basin::Other => format!(
                "MIXED — {} true facets missed AND {} claimed beyond the boundary",
                self.missed.len(),
                self.overshoot.len()
            ),
            Basin::NoConverge => format!(
                "NO FIXPOINT — still moving after {} refit rounds ({} facets)",
                self.rounds,
                self.result.len()
            ),
            Basin::Exact => format!("EXACT — {} facets, {} refit round(s)", n_truth, self.rounds),
        }
    }
}

/// Why did a seed set's own fit refuse its own facets?
///
/// The question this answers is whether the failure is a shortage of
/// *information* or a failure of the *fitter*: it prints the residual the
/// reference's own exact parameters achieve on the very same facets — the floor
/// any fit of that type could reach — beside the residual the fitter actually
/// achieved, and the two parameter sets side by side. Gated by
/// `BREPPER_INV_FITDUMP=1`.
fn fit_dump(
    mesh: &ConnectedMesh,
    seed: &[usize],
    ref_surf: &RefSurface,
    kind: Kind,
    vtol: f64,
    quantum: f64,
    label: &str,
) {
    let verts = vertices_of(mesh, seed);
    let truth_max = verts
        .iter()
        .filter_map(|&v| ref_distance(ref_surf, vertex_pos(mesh, v)))
        .fold(0.0_f64, f64::max);
    eprintln!(
        "  [fit] {label}: {} facets, {} vertices | the REFERENCE's own parameters fit these \
         vertices to {:.3e} mm = {:.2}x vtol = {:.2} quanta",
        seed.len(),
        verts.len(),
        truth_max,
        truth_max / vtol,
        truth_max / quantum,
    );
    match kind {
        Kind::Torus => {
            let RefSurface::Torus {
                centre,
                axis_dir,
                major,
                minor,
            } = *ref_surf
            else {
                return;
            };
            match ToroidalSurface::fit_surface(mesh, seed, &verts, 0) {
                Some(f) => {
                    let fitted_max = verts
                        .iter()
                        .map(|&v| ToroidalSurface::vertex_distance(&f, &vertex_pos(mesh, v)))
                        .fold(0.0_f64, f64::max);
                    let d = unit(axis_dir);
                    let fd = unit(f.axis_direction);
                    let axis_err = norm([
                        fd[1] * d[2] - fd[2] * d[1],
                        fd[2] * d[0] - fd[0] * d[2],
                        fd[0] * d[1] - fd[1] * d[0],
                    ])
                    .clamp(0.0, 1.0)
                    .asin()
                    .to_degrees();
                    eprintln!(
                        "  [fit]   truth: R {major:.6} r {minor:.6} centre ({:.4},{:.4},{:.4})",
                        centre[0], centre[1], centre[2]
                    );
                    eprintln!(
                        "  [fit]   fitted: R {:.6} r {:.6} centre ({:.4},{:.4},{:.4})  \
                         axis off by {axis_err:.4}deg  -> residual {:.3e} mm = {:.1}x vtol",
                        f.major_radius,
                        f.minor_radius,
                        f.center[0],
                        f.center[1],
                        f.center[2],
                        fitted_max,
                        fitted_max / vtol,
                    );
                    eprintln!(
                        "  [fit]   R error {:+.4} mm ({:+.1}%), r error {:+.6} mm ({:+.2}%)",
                        f.major_radius - major,
                        100.0 * (f.major_radius - major) / major.max(1e-30),
                        f.minor_radius - minor,
                        100.0 * (f.minor_radius - minor) / minor.max(1e-30),
                    );
                    // Two quantities that say WHICH kind of failure this is.
                    //
                    // `R - r` is the fitted torus's inner-equator radius. When
                    // it equals the TRUE minor radius, the fitter has presented
                    // the saddle at its own inner equator to the patch and
                    // matched the fillet's tube curvature with the wrong
                    // principal direction — which forces its axis 90 degrees off.
                    //
                    // The sweep of the patch about the FITTED axis says whether
                    // R was identifiable there at all: below a few degrees the
                    // objective has a valley along R and the optimizer is free
                    // to wander down it, which is why R varies wildly between
                    // seeds while `R - r` stays pinned.
                    let sweep = {
                        let d = unit(f.axis_direction);
                        let (bu, bv) = super::perpendicular_basis(&d);
                        let mut ang: Vec<f64> = verts
                            .iter()
                            .filter_map(|&vi| {
                                let w = sub(vertex_pos(mesh, vi), f.center);
                                let z = dot(w, d);
                                let rad = sub(w, [d[0] * z, d[1] * z, d[2] * z]);
                                let (x, y) = (dot(rad, bu), dot(rad, bv));
                                (x.hypot(y) > 1e-12).then(|| y.atan2(x))
                            })
                            .collect();
                        if ang.len() < 2 {
                            0.0
                        } else {
                            ang.sort_by(|a, b| a.partial_cmp(b).unwrap());
                            let mut gap = ang[0] + std::f64::consts::TAU - ang[ang.len() - 1];
                            for w in ang.windows(2) {
                                gap = gap.max(w[1] - w[0]);
                            }
                            (std::f64::consts::TAU - gap).max(0.0).to_degrees()
                        }
                    };
                    let sse = |g: &super::toroidal::TorusFit| -> f64 {
                        verts
                            .iter()
                            .map(|&v| {
                                ToroidalSurface::vertex_distance(g, &vertex_pos(mesh, v)).powi(2)
                            })
                            .sum::<f64>()
                    };
                    let truth_fit = super::toroidal::TorusFit {
                        center: centre,
                        axis_direction: axis_dir,
                        major_radius: major,
                        minor_radius: minor,
                        convex: None,
                    };
                    eprintln!(
                        "  [fit]   fitted R-r = {:.4} (true r = {minor:.4}) | patch sweeps \
                         {sweep:.2}deg about the FITTED axis | sum-sq residual: truth {:.3e} vs \
                         fitted {:.3e} ({:.0}x worse)",
                        f.major_radius - f.minor_radius,
                        sse(&truth_fit),
                        sse(&f),
                        sse(&f) / sse(&truth_fit).max(1e-300),
                    );
                    // Would Phase 1.1's closed-form initializer (built, unit
                    // tested, never wired into the pipeline) have started the
                    // solve in the right basin on this very seed set?
                    match super::initialize::torus(seed, mesh) {
                        Some(super::initialize::Initial::Torus {
                            centre: ic,
                            axis_dir: id,
                            major: im,
                            minor: ir,
                        }) => {
                            let idu = unit(id);
                            let ierr = norm([
                                idu[1] * d[2] - idu[2] * d[1],
                                idu[2] * d[0] - idu[0] * d[2],
                                idu[0] * d[1] - idu[1] * d[0],
                            ])
                            .clamp(0.0, 1.0)
                            .asin()
                            .to_degrees();
                            let ifit = super::toroidal::TorusFit {
                                center: ic,
                                axis_direction: id,
                                major_radius: im,
                                minor_radius: ir,
                                convex: None,
                            };
                            let imax = verts
                                .iter()
                                .map(|&v| {
                                    ToroidalSurface::vertex_distance(&ifit, &vertex_pos(mesh, v))
                                })
                                .fold(0.0_f64, f64::max);
                            eprintln!(
                                "  [fit]   Phase-1.1 closed form: R {im:.6} r {ir:.6} axis off \
                                 {ierr:.4}deg -> residual {imax:.3e} mm = {:.1}x vtol",
                                imax / vtol
                            );
                        }
                        _ => eprintln!("  [fit]   Phase-1.1 closed form: declined"),
                    }
                }
                None => eprintln!("  [fit]   fitted: none (fit_torus declined)"),
            }
        }
        Kind::Cone => {
            let RefSurface::Cone {
                apex,
                axis_dir,
                half_angle,
            } = *ref_surf
            else {
                return;
            };
            let extent = {
                let mut e = 0.0_f64;
                let c = vertex_pos(mesh, verts[0]);
                for &v in &verts {
                    e = e.max(norm(sub(vertex_pos(mesh, v), c)));
                }
                e.max(1e-30)
            };
            let d = unit(axis_dir);
            eprintln!(
                "  [fit]   truth: half-angle {:.4}deg apex ({:.4},{:.4},{:.4})",
                half_angle.to_degrees(),
                apex[0],
                apex[1],
                apex[2]
            );
            let report = |tag: &str, fapex: [f64; 3], fdir: [f64; 3], fhalf: f64| {
                let fd = unit(fdir);
                let ax = norm([
                    fd[1] * d[2] - fd[2] * d[1],
                    fd[2] * d[0] - fd[0] * d[2],
                    fd[0] * d[1] - fd[1] * d[0],
                ])
                .clamp(0.0, 1.0)
                .asin()
                .to_degrees();
                let f = super::conical::ConeFit {
                    apex: fapex,
                    axis_direction: fdir,
                    half_angle: fhalf,
                    radius: 0.0,
                    convex: None,
                };
                let (mx, sse) = verts.iter().fold((0.0_f64, 0.0_f64), |(m, s2), &v| {
                    let e = ConicalSurface::vertex_distance(&f, &vertex_pos(mesh, v));
                    (m.max(e), s2 + e * e)
                });
                eprintln!(
                    "  [fit]   {tag}: half-angle {:.4}deg ({:+.4} off) axis {ax:.4}deg off, \
                     apex ({:.3},{:.3},{:.3}) = {:.2}x the patch extent away -> residual \
                     {mx:.3e} mm = {:.1}x vtol, sse {sse:.3e}",
                    fhalf.to_degrees(),
                    (fhalf - half_angle).to_degrees(),
                    fapex[0],
                    fapex[1],
                    fapex[2],
                    norm(sub(fapex, apex)) / extent,
                    mx / vtol,
                );
            };
            let truth_sse: f64 = {
                let tf = super::conical::ConeFit {
                    apex,
                    axis_direction: axis_dir,
                    half_angle,
                    radius: 0.0,
                    convex: None,
                };
                verts
                    .iter()
                    .map(|&v| ConicalSurface::vertex_distance(&tf, &vertex_pos(mesh, v)).powi(2))
                    .sum()
            };
            // --- the two-circle construction, and why it is or is not enough ---
            //
            // Two ruling-aligned quads on a cone band give 6 vertices in two
            // rings of 3, one on each rim. Three points in space determine a
            // unique circle, so each rim gives a centre, a radius and a normal;
            // the line through the two centres is the axis and the two radii at
            // two heights give the half-angle by similar triangles. Exact
            // arithmetic, no fitting. What follows is that construction carried
            // out on this actual seed, with its conditioning.
            {
                let mut rings: Vec<(f64, Vec<[f64; 3]>)> = Vec::new();
                for &vi in &verts {
                    let p = vertex_pos(mesh, vi);
                    let z = dot(sub(p, apex), d);
                    match rings.iter_mut().find(|(rz, _)| (rz - z).abs() < 1e-6) {
                        Some((_, v)) => v.push(p),
                        None => rings.push((z, vec![p])),
                    }
                }
                rings.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                eprintln!(
                    "  [fit]   two-circle check: {} vertices fall into {} ring(s) along the true \
                     axis",
                    verts.len(),
                    rings.len()
                );
                // Circumcircle of three points, and the arc they span.
                let circum = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| -> Option<([f64; 3], f64)> {
                    let ab = sub(b, a);
                    let ac = sub(c, a);
                    let n = [
                        ab[1] * ac[2] - ab[2] * ac[1],
                        ab[2] * ac[0] - ab[0] * ac[2],
                        ab[0] * ac[1] - ab[1] * ac[0],
                    ];
                    let den = 2.0 * dot(n, n);
                    if den <= 0.0 {
                        return None;
                    }
                    let t1 = [
                        n[1] * ab[2] - n[2] * ab[1],
                        n[2] * ab[0] - n[0] * ab[2],
                        n[0] * ab[1] - n[1] * ab[0],
                    ];
                    let t2 = [
                        ac[1] * n[2] - ac[2] * n[1],
                        ac[2] * n[0] - ac[0] * n[2],
                        ac[0] * n[1] - ac[1] * n[0],
                    ];
                    let (ac2, ab2) = (dot(ac, ac), dot(ab, ab));
                    let ctr = [
                        a[0] + (t1[0] * ac2 + t2[0] * ab2) / den,
                        a[1] + (t1[1] * ac2 + t2[1] * ab2) / den,
                        a[2] + (t1[2] * ac2 + t2[2] * ab2) / den,
                    ];
                    Some((ctr, norm(sub(a, ctr))))
                };
                let mut built: Vec<([f64; 3], f64, f64, f64)> = Vec::new(); // ctr, r, r_true, arc
                for (z, pts) in &rings {
                    let r_true = z * half_angle.tan();
                    let (bu, bv) = super::perpendicular_basis(&d);
                    let mut ang: Vec<f64> = pts
                        .iter()
                        .map(|&p| {
                            let w = sub(p, apex);
                            let zz = dot(w, d);
                            let rad = sub(w, [d[0] * zz, d[1] * zz, d[2] * zz]);
                            dot(rad, bv).atan2(dot(rad, bu))
                        })
                        .collect();
                    ang.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    let arc = (ang[ang.len() - 1] - ang[0]).to_degrees();
                    if pts.len() >= 3 {
                        if let Some((ctr, r)) = circum(pts[0], pts[1], pts[2]) {
                            eprintln!(
                                "  [fit]     ring at z {z:8.4}: {} pts spanning {arc:6.3}deg | \
                                 true radius {r_true:8.5} | circumcircle radius {r:8.5} \
                                 (err {:+.3e} mm = {:.1}x vtol)",
                                pts.len(),
                                r - r_true,
                                (r - r_true).abs() / vtol,
                            );
                            built.push((ctr, r, r_true, arc));
                        }
                    } else {
                        eprintln!(
                            "  [fit]     ring at z {z:8.4}: {} pts spanning {arc:6.3}deg | \
                             true radius {r_true:8.5} | fewer than 3 points, no circle",
                            pts.len()
                        );
                    }
                }
                if built.len() >= 2 {
                    let (c0, r0, _, arc0) = built[0];
                    let (c1, r1, _, _) = built[built.len() - 1];
                    let axis = unit(sub(c1, c0));
                    let ax_err = norm([
                        axis[1] * d[2] - axis[2] * d[1],
                        axis[2] * d[0] - axis[0] * d[2],
                        axis[0] * d[1] - axis[1] * d[0],
                    ])
                    .clamp(0.0, 1.0)
                    .asin()
                    .to_degrees();
                    let h = norm(sub(c1, c0));
                    let ha = ((r1 - r0) / h).atan().to_degrees();
                    eprintln!(
                        "  [fit]   two-circle result: axis {ax_err:.4}deg off, half-angle \
                         {ha:.4}deg (true {:.4}, err {:+.4}deg)",
                        half_angle.to_degrees(),
                        ha - half_angle.to_degrees(),
                    );
                    // The sensitivity that governs it: a circle through points
                    // spanning half-arc `a` recovers its radius with error
                    // ~2*eps/a^2, and the half-angle takes the DIFFERENCE of two
                    // such radii over the height, so the error is divided by h.
                    let a_rad = (arc0 / 2.0).to_radians();
                    let eps = truth_max.max(1e-12);
                    let dr = 2.0 * eps / (a_rad * a_rad).max(1e-30);
                    eprintln!(
                        "  [fit]   conditioning: half-arc {:.3}deg, mesh error {eps:.2e} mm -> \
                         radius uncertainty ~2e/a^2 = {dr:.3e} mm ({:.0}x vtol); over height \
                         {h:.3} that is a half-angle uncertainty of ~{:.3}deg",
                        arc0 / 2.0,
                        dr / vtol,
                        (2.0 * dr / h).atan().to_degrees(),
                    );
                }
            }
            report("truth-as-fit", apex, axis_dir, half_angle);
            // Verbosity 4 so the fitter's own rejection reasons (half-angle
            // bounds, apex distance, `axis_normals_consistent`) reach stderr —
            // the dump is where those belong.
            match ConicalSurface::fit_surface(mesh, seed, &verts, 4) {
                Some(f) => {
                    if f.half_angle.abs() < 1e-9 {
                        eprintln!(
                            "  [fit]   fitted: a CYLINDER of radius {:.6} (the cone rung was not \
                             taken) -> residual {:.3e} mm = {:.1}x vtol",
                            f.radius,
                            verts
                                .iter()
                                .map(|&v| ConicalSurface::vertex_distance(&f, &vertex_pos(mesh, v)))
                                .fold(0.0_f64, f64::max),
                            verts
                                .iter()
                                .map(|&v| ConicalSurface::vertex_distance(&f, &vertex_pos(mesh, v)))
                                .fold(0.0_f64, f64::max)
                                / vtol,
                        );
                    } else {
                        report("fitted", f.apex, f.axis_direction, f.half_angle);
                        let fsse: f64 = verts
                            .iter()
                            .map(|&v| {
                                ConicalSurface::vertex_distance(&f, &vertex_pos(mesh, v)).powi(2)
                            })
                            .sum();
                        // The distinction that decides whether this is fixable.
                        if fsse < truth_sse {
                            eprintln!(
                                "  [fit]   VERDICT: NOT IDENTIFIABLE — the fitted cone explains \
                                 these {} vertices {:.0}x BETTER than the true one \
                                 ({fsse:.3e} vs {truth_sse:.3e}). A second exact cone passes \
                                 through them; no optimizer can prefer the truth.",
                                verts.len(),
                                truth_sse / fsse.max(1e-300),
                            );
                        } else {
                            eprintln!(
                                "  [fit]   VERDICT: CONVERGENCE — the truth explains them \
                                 {:.0}x better ({truth_sse:.3e} vs {fsse:.3e}), so the minimum \
                                 is at the truth and the fit failed to reach it.",
                                fsse / truth_sse.max(1e-300),
                            );
                        }
                    }
                }
                None => eprintln!("  [fit]   fitted: none (ConicalSurface declined)"),
            }
            match super::initialize::cone(seed, mesh, vtol) {
                Some(super::initialize::Initial::Cone {
                    apex: ia,
                    axis_dir: id,
                    half_angle: ih,
                }) => report("Phase-1.1 closed form", ia, id, ih),
                _ => eprintln!("  [fit]   Phase-1.1 closed form: declined"),
            }
        }
        Kind::Cylinder => {
            if let Some(f) = ConicalSurface::fit_surface(mesh, seed, &verts, 0) {
                let fitted_max = verts
                    .iter()
                    .map(|&v| ConicalSurface::vertex_distance(&f, &vertex_pos(mesh, v)))
                    .fold(0.0_f64, f64::max);
                eprintln!(
                    "  [fit]   fitted: half-angle {:.4}deg radius {:.6} -> residual {:.3e} mm \
                     = {:.1}x vtol",
                    f.half_angle.to_degrees(),
                    f.radius,
                    fitted_max,
                    fitted_max / vtol,
                );
            }
        }
        Kind::Sphere => {
            if let Some(f) = SphericalSurface::fit_surface(mesh, seed, &verts, 0) {
                let fitted_max = verts
                    .iter()
                    .map(|&v| SphericalSurface::vertex_distance(&f, &vertex_pos(mesh, v)))
                    .fold(0.0_f64, f64::max);
                eprintln!(
                    "  [fit]   fitted: r {:.6} -> residual {:.3e} mm = {:.1}x vtol",
                    f.radius,
                    fitted_max,
                    fitted_max / vtol
                );
            }
        }
        _ => {}
    }
}

/// Grow one seed set and work out, in detail, what happened to it.
#[allow(clippy::too_many_arguments)]
fn diagnose<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    seed: &[usize],
    truth: &HashSet<usize>,
    centroids: &[[f64; 3]],
    config: &Config,
    shape: &'static str,
    ref_surf: Option<(&RefSurface, Kind, f64, &str)>,
) -> SeedCase {
    let (verdict, result, rounds, self_ok, fit) =
        grow_clean::<S>(mesh, in_region, seed, centroids, config);
    analyse::<S>(
        mesh, in_region, seed, truth, centroids, config, shape, ref_surf, verdict, result, rounds,
        self_ok, fit,
    )
}

/// Work out, in detail, what happened to a grown seed set.
#[allow(clippy::too_many_arguments)]
fn analyse<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    seed: &[usize],
    truth: &HashSet<usize>,
    centroids: &[[f64; 3]],
    config: &Config,
    shape: &'static str,
    ref_surf: Option<(&RefSurface, Kind, f64, &str)>,
    mut verdict: Basin,
    result: Vec<usize>,
    rounds: usize,
    self_ok: bool,
    fit: Option<S::Fit>,
) -> SeedCase {
    let vtol = config.vertex_tolerance_mm;
    let stol = config.surface_tolerance_mm;
    if verdict == Basin::Exact {
        verdict = classify(&result, truth);
    }
    let rset: HashSet<usize> = result.iter().copied().collect();
    let mut missed: Vec<usize> = truth.difference(&rset).copied().collect();
    let mut overshoot: Vec<usize> = rset.difference(truth).copied().collect();
    missed.sort_unstable();
    overshoot.sort_unstable();

    // Split the missed facets by CAUSE: refused by validation at the frontier,
    // not mesh-reachable, or cut into a different smooth region by stage 2.1
    // (in which case no seed can reach them at all).
    let mut blocked: Vec<(usize, FacetVerdict)> = Vec::new();
    let mut unreachable: Vec<usize> = Vec::new();
    let mut outside_region: Vec<usize> = Vec::new();
    for &m in &missed {
        if !in_region[m] {
            outside_region.push(m);
            continue;
        }
        let vc = mesh.faces[m].vertex_count as usize;
        let touches = mesh.faces[m].neighbors[..vc]
            .iter()
            .any(|&ni| ni >= 0 && rset.contains(&(ni as usize)));
        match (touches, fit.as_ref()) {
            (true, Some(f)) => {
                blocked.push((m, facet_verdict::<S>(mesh, m, f, centroids, vtol, stol)))
            }
            (true, None) => {}
            (false, _) => unreachable.push(m),
        }
    }

    // Every seed facet against the seed set's OWN first fit — the admission
    // quantity, clause by clause.
    let verts = vertices_of(mesh, seed);
    let mut seed_verdicts: Vec<(usize, FacetVerdict)> = Vec::new();
    if let Some(mut f0) = S::fit_surface(mesh, seed, &verts, 0) {
        S::determine_convexity(mesh, seed, &mut f0);
        for &fi in seed {
            seed_verdicts.push((fi, facet_verdict::<S>(mesh, fi, &f0, centroids, vtol, stol)));
        }
    }

    if let Some((rs, kind, quantum, label)) = ref_surf {
        if !self_ok && std::env::var("BREPPER_INV_FITDUMP").is_ok() {
            fit_dump(mesh, seed, rs, kind, vtol, quantum, label);
        }
    }

    SeedCase {
        shape,
        seed: seed.to_vec(),
        verdict,
        result,
        rounds,
        self_ok,
        missed,
        overshoot,
        blocked,
        unreachable,
        outside_region,
        truth_connected: {
            let mut t: Vec<usize> = truth.iter().copied().collect();
            t.sort_unstable();
            is_connected(mesh, &t)
        },
        seed_verdicts,
    }
}

/// **The invariant, as stated:** for every facet of a surface, is there at
/// least one conforming seed set containing it that grows to a superset of the
/// surface?
///
/// Existential in the seed, universal in the facet — which is the right shape
/// for an algorithm that gets to *choose* its seeds. It is much weaker than
/// "every conforming seed works", and it is the property that actually has to
/// hold for a surface to be reachable at all.
///
/// Returns one outcome per sampled anchor facet.
struct AnchorOutcome {
    anchor: usize,
    /// Conforming seed sets containing this facet.
    attempts: Vec<(&'static str, Vec<usize>)>,
    /// The first one that grew to a superset, if any.
    winner: Option<usize>,
    /// Diagnosis of each attempt, in order (only filled up to the winner).
    cases: Vec<SeedCase>,
}

#[allow(clippy::too_many_arguments)]
fn probe_anchor<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    member: &HashSet<usize>,
    vstar: &HashMap<usize, Vec<usize>>,
    truth: &HashSet<usize>,
    centroids: &[[f64; 3]],
    config: &Config,
    kind: Kind,
    anchor: usize,
    ref_surf: &RefSurface,
    quantum: f64,
    label: &str,
) -> AnchorOutcome {
    let max_chain = config.viz_seed_sizes.iter().copied().max().unwrap_or(4);
    let attempts = seed_family(
        mesh,
        member,
        vstar,
        kind,
        anchor,
        max_chain,
        config.viz_seed_coplanar_clause,
        config.viz_seed_spec_only,
        config.viz_seed_max_attempts.max(4),
    );
    let mut cases = Vec::new();
    let mut winner = None;
    for (i, (shape, set)) in attempts.iter().enumerate() {
        let lbl = format!("{label} facet {anchor} seed {shape} x{}", set.len());
        let rs = Some((ref_surf, kind, quantum, lbl.as_str()));
        let case = diagnose::<S>(mesh, in_region, set, truth, centroids, config, shape, rs);
        let covers = case.covers();
        cases.push(case);
        if covers {
            winner = Some(i);
            break;
        }
    }
    AnchorOutcome {
        anchor,
        attempts,
        winner,
        cases,
    }
}

/// Interactive walk, one **anchor facet** at a time. A facet is shown only when
/// *every* conforming seed set containing it fails, and then each of its
/// attempts is shown in turn with its own diagnosis.
#[allow(clippy::too_many_arguments)]
fn viz_face<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    faces: &[usize],
    truth: &HashSet<usize>,
    centroids: &[[f64; 3]],
    config: &Config,
    viz: Option<&crate::viz::VizSender>,
    model: &str,
    ref_idx: usize,
    kind_name: &str,
    kind: Kind,
    ref_surf: &RefSurface,
    quantum: f64,
) -> bool {
    let member: HashSet<usize> = faces.iter().copied().collect();
    let vstar = vertex_star_map(mesh, faces);
    let n_sample = config.viz_seed_samples.max(1);
    let stride = (faces.len() / n_sample).max(1);
    let anchors: Vec<usize> = faces.iter().copied().step_by(stride).take(n_sample).collect();
    let show = config.viz_seed_show.as_str();

    for a in anchors {
        let out = probe_anchor::<S>(
            mesh, in_region, &member, &vstar, truth, centroids, config, kind, a, ref_surf,
            quantum, &format!("{model} ref#{ref_idx}"),
        );
        let facet_ok = out.winner.is_some();
        let wanted = match show {
            "all" => true,
            "covering" => facet_ok,
            _ => !facet_ok,
        };
        if !wanted {
            continue;
        }
        // One machine row per anchor facet: the invariant, as stated.
        println!(
            "ANCHOR\t{model}\t{ref_idx}\t{kind_name}\t{}\t{}\t{}\t{}\t{}\t{}",
            faces.len(),
            out.anchor,
            out.attempts.len(),
            u8::from(out.winner.is_some()),
            out.winner.map(|w| out.attempts[w].0).unwrap_or("-"),
            out.cases
                .last()
                .map(|c| match c.verdict {
                    Basin::Exact => "exact",
                    Basin::Leak => "leak",
                    Basin::Short if !c.outside_region.is_empty() => "partition-cut",
                    Basin::Short => "stalled",
                    Basin::NoFit => "no-fit",
                    Basin::SeedRejected => "seed-rejected",
                    Basin::NoConverge => "no-fixpoint",
                    Basin::Other => "mixed",
                })
                .unwrap_or("no-conforming-seed"),
        );
        if out.attempts.is_empty() {
            let status = format!(
                "{model} ref#{ref_idx} {kind_name} ({} facets) | facet {a} | NO CONFORMING SEED \
                 SET EXISTS containing this facet under the section's rule for a {kind_name}",
                faces.len()
            );
            if !emit(viz, mesh, None, truth, &status) {
                return false;
            }
            continue;
        }
        for (i, case) in out.cases.iter().enumerate() {
            let verdict = if facet_ok && i == out.cases.len() - 1 {
                "FACET OK — this seed set covers".to_string()
            } else if facet_ok {
                case.why(truth.len())
            } else {
                format!("FACET FAILS ({}/{} attempts) — {}", i + 1, out.attempts.len(), case.why(truth.len()))
            };
            let seed_worst = case
                .seed_verdicts
                .iter()
                .map(|(_, v)| v.vratio)
                .fold(0.0_f64, f64::max);
            let status = format!(
                "{model} ref#{ref_idx} {kind_name} ({} facets) | facet {a} | seed {} x{} | \
                 own fit: vertex {:.2}x vtol{} | {verdict}",
                faces.len(),
                case.shape,
                case.seed.len(),
                seed_worst,
                if case.seed_verdicts.is_empty() {
                    ", no fit exists".to_string()
                } else if case.self_ok {
                    String::new()
                } else {
                    format!(", ADMISSION FAILS on {}", dominant_cause(&case.seed_verdicts))
                },
            );
            if !emit(viz, mesh, Some((a, case)), truth, &status) {
                return false;
            }
        }
    }
    true
}

/// Draw one step (or print it, headless). Returns false when the user quits.
fn emit(
    viz: Option<&crate::viz::VizSender>,
    mesh: &ConnectedMesh,
    case: Option<(usize, &SeedCase)>,
    truth: &HashSet<usize>,
    status: &str,
) -> bool {
    let Some(viz) = viz else {
        let detail = case
            .and_then(|(_, c)| c.blocked.iter().find(|(_, v)| !v.ok()).map(|(_, v)| v.cause()))
            .map(|d| format!("  [blocked by {d}]"))
            .unwrap_or_default();
        eprintln!("  [seed] {status}{detail}");
        return true;
    };
    let mut ov = crate::viz::VizOverlay::new();
    let anchor = match case {
        Some((a, c)) => {
            let blocked_only: Vec<usize> = c
                .blocked
                .iter()
                .filter(|(_, v)| !v.ok())
                .map(|&(f, _)| f)
                .collect();
            let covered: Vec<usize> = c
                .result
                .iter()
                .copied()
                .filter(|f| truth.contains(f) && !c.seed.contains(f))
                .collect();
            let missed_only: Vec<usize> = c
                .missed
                .iter()
                .copied()
                .filter(|f| !blocked_only.contains(f))
                .collect();
            for (idxs, color) in [
                (covered, C_COVERED),
                (c.overshoot.clone(), C_OVERSHOOT),
                (missed_only, C_MISSED),
                (blocked_only, C_BLOCKED),
                (c.seed.clone(), C_SEED),
                (vec![a], C_ANCHOR),
            ] {
                if !idxs.is_empty() {
                    ov.face_highlights.push(crate::viz::FaceHighlight {
                        face_indices: idxs,
                        color,
                    });
                }
            }
            a
        }
        None => *truth.iter().next().unwrap_or(&0),
    };
    let c = face_centroid(&mesh.faces[anchor], &mesh.vertices);
    ov.focus_point = Some([c[0] as f32, c[1] as f32, c[2] as f32]);
    let n = facet_normal(mesh, anchor);
    ov.focus_normal = Some([n[0] as f32, n[1] as f32, n[2] as f32]);
    ov.status_text = status.to_string();
    !matches!(viz.show_and_wait(ov), crate::viz::VizAction::Quit)
}

/// Sample seed SETS inside one reference face and report, for each, whether the/// Sample seed SETS inside one reference face and report, for each, whether the
/// clean operator grows it to exactly that face.
///
/// Three shapes, because the answer may depend on more than size:
/// - `ball`   — a BFS neighbourhood of `k` facets, the compact seed.
/// - `split` — two BFS balls of `k/2` at the two most distant starts, which is
///   the well-conditioned shape for an axis (Phase 1.0's "angularly most
///   distant pair") and is *disconnected*.
/// - `strip` — `k` facets taken in index order from one start, which on a strip
///   tessellation walks the ruling (Phase 1.2's measured rank-0 failure).
#[allow(clippy::too_many_arguments)]
fn basin_probe<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    faces: &[usize],
    truth: &HashSet<usize>,
    centroids: &[[f64; 3]],
    config: &Config,
    truth_axis: Option<([f64; 3], [f64; 3])>,
    ref_of: &[Option<usize>],
    ref_surfaces: &[RefSurface],
    model: &str,
    ref_idx: usize,
    kind: &str,
) {
    // Reference faces adjacent to this one across a mesh edge — the only
    // places the section's second clause allows an overshoot to land.
    let mut neighbours: HashSet<usize> = HashSet::new();
    for &fi in faces {
        let vc = mesh.faces[fi].vertex_count as usize;
        for &ni in &mesh.faces[fi].neighbors[..vc] {
            if ni < 0 {
                continue;
            }
            if let Some(o) = ref_of[ni as usize] {
                if o != ref_idx {
                    neighbours.insert(o);
                }
            }
        }
    }
    const SIZES: [usize; 8] = [2, 3, 4, 6, 8, 12, 16, 24];
    const STARTS: usize = 5;
    let member: HashSet<usize> = faces.iter().copied().collect();
    let ball = |start: usize, k: usize| -> Vec<usize> {
        let mut out = vec![start];
        let mut seen: HashSet<usize> = HashSet::from([start]);
        let mut frontier = vec![start];
        while out.len() < k && !frontier.is_empty() {
            let mut ring = Vec::new();
            for &fi in &frontier {
                let vc = mesh.faces[fi].vertex_count as usize;
                for &ni in &mesh.faces[fi].neighbors[..vc] {
                    if ni < 0 || out.len() >= k {
                        continue;
                    }
                    let nj = ni as usize;
                    if member.contains(&nj) && seen.insert(nj) {
                        out.push(nj);
                        ring.push(nj);
                    }
                }
            }
            if ring.is_empty() {
                break;
            }
            frontier = ring;
        }
        out.sort_unstable();
        out
    };
    // Geometric spread of a seed set about the true axis, for the analysis.
    let spread = |set: &[usize]| -> (f64, f64) {
        let Some((tpos, tdir0)) = truth_axis else {
            return (f64::NAN, f64::NAN);
        };
        let tdir = unit(tdir0);
        let (bu, bv) = super::perpendicular_basis(&tdir);
        let mut angles: Vec<f64> = Vec::new();
        let (mut zmin, mut zmax, mut rmax) = (f64::INFINITY, f64::NEG_INFINITY, 0.0_f64);
        for &fi in set {
            for p in facet_points(mesh, fi) {
                let w = sub(p, tpos);
                let z = dot(w, tdir);
                let rad = sub(w, [tdir[0] * z, tdir[1] * z, tdir[2] * z]);
                let (x, y) = (dot(rad, bu), dot(rad, bv));
                let rr = x.hypot(y);
                if rr > 1e-12 {
                    angles.push(y.atan2(x));
                }
                zmin = zmin.min(z);
                zmax = zmax.max(z);
                rmax = rmax.max(rr);
            }
        }
        if angles.len() < 2 || rmax <= 0.0 {
            return (0.0, 0.0);
        }
        angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut gap = angles[0] + std::f64::consts::TAU - angles[angles.len() - 1];
        for w in angles.windows(2) {
            gap = gap.max(w[1] - w[0]);
        }
        (
            (std::f64::consts::TAU - gap).max(0.0).to_degrees(),
            (zmax - zmin) / rmax,
        )
    };

    let stride = (faces.len() / STARTS).max(1);
    let starts: Vec<usize> = faces.iter().copied().step_by(stride).take(STARTS).collect();
    let far = *faces.last().unwrap();
    for &start in &starts {
        for &k in &SIZES {
            if k > faces.len() {
                continue;
            }
            let mut cases: Vec<(&str, Vec<usize>)> = vec![("ball", ball(start, k))];
            if k >= 4 {
                let mut sp = ball(start, k / 2);
                sp.extend(ball(far, k - k / 2));
                sp.sort_unstable();
                sp.dedup();
                cases.push(("split", sp));
            }
            {
                let i0 = faces.iter().position(|&f| f == start).unwrap_or(0);
                let strip: Vec<usize> =
                    faces.iter().cycle().skip(i0).take(k).copied().collect();
                let mut st = strip;
                st.sort_unstable();
                st.dedup();
                cases.push(("strip", st));
            }
            for (shape, set) in cases {
                if set.len() < 2 {
                    continue;
                }
                let connected = is_connected(mesh, &set);
                let (sweep, axial) = spread(&set);
                let (mut verdict, result, rounds, self_ok, fit) =
                    grow_clean::<S>(mesh, in_region, &set, centroids, config);
                if verdict == Basin::Exact {
                    verdict = classify(&result, truth);
                }
                let inside = result.iter().filter(|f| truth.contains(f)).count();
                // --- the section's SECOND clause, scored directly ---
                // Every overshoot facet should (a) belong to a reference face
                // adjacent to this one, and (b) be explained strictly better by
                // that face's own surface than by the fit that claimed it.
                let over: Vec<usize> =
                    result.iter().copied().filter(|f| !truth.contains(f)).collect();
                let mut over_refs: HashSet<usize> = HashSet::new();
                let mut over_adj = true;
                let mut over_better = 0usize;
                let mut over_scored = 0usize;
                for &g in &over {
                    // An unlabelled (seam-ambiguous) facet is not evidence
                    // either way, so it is not counted against the clause.
                    if let Some(o) = ref_of[g] {
                        {
                            over_refs.insert(o);
                            if !neighbours.contains(&o) {
                                over_adj = false;
                            }
                            if let Some(f) = fit.as_ref() {
                                let mut d_fit = 0.0_f64;
                                let mut d_ref = 0.0_f64;
                                let mut ok = true;
                                for p in facet_points(mesh, g) {
                                    d_fit = d_fit.max(S::vertex_distance(f, &p));
                                    match ref_distance(&ref_surfaces[o], p) {
                                        Some(d) => d_ref = d_ref.max(d),
                                        None => ok = false,
                                    }
                                }
                                if ok {
                                    over_scored += 1;
                                    if d_ref < d_fit {
                                        over_better += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                let short_n = truth.len() - inside;
                let jacc = if result.is_empty() && truth.is_empty() {
                    1.0
                } else {
                    inside as f64 / (result.len() + truth.len() - inside).max(1) as f64
                };
                println!(
                    "BASIN\t{model}\t{ref_idx}\t{kind}\t{}\t{shape}\t{}\t{}\t{}\t{}\t{}\t\
                     {sweep:.3}\t{axial:.5}\t{jacc:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    faces.len(),
                    set.len(),
                    u8::from(connected),
                    u8::from(self_ok),
                    verdict.code(),
                    rounds,
                    result.len(),
                    over.len(),
                    short_n,
                    over_refs.len(),
                    u8::from(over_adj),
                    over_better,
                    over_scored,
                );
            }
        }
    }
}

/// Is a facet set connected through the mesh's face adjacency?
fn is_connected(mesh: &ConnectedMesh, set: &[usize]) -> bool {
    if set.is_empty() {
        return true;
    }
    let member: HashSet<usize> = set.iter().copied().collect();
    let mut seen: HashSet<usize> = HashSet::from([set[0]]);
    let mut stack = vec![set[0]];
    while let Some(fi) = stack.pop() {
        let vc = mesh.faces[fi].vertex_count as usize;
        for &ni in &mesh.faces[fi].neighbors[..vc] {
            if ni < 0 {
                continue;
            }
            let nj = ni as usize;
            if member.contains(&nj) && seen.insert(nj) {
                stack.push(nj);
            }
        }
    }
    seen.len() == set.len()
}

/// How large does a fit support have to be before it determines a developable?
///
/// The Invariants section legislates this with a vertex count ("exactly 6"),
/// and `compete.rs` answers it with `SEED_PATCH_CAP` = 32 facets at
/// `MAX_PATCH_DEPTH` = 3 hops. Neither is measured. This grows a BFS
/// neighbourhood inside one reference face to a ladder of sizes and reports the
/// closed-form axis error at each, so the question is answered by a curve
/// rather than by a constant.
///
/// Starts are taken at a stride over the face so the sample is spread rather
/// than clustered, and growth is confined to the reference face — this measures
/// how much SUPPORT a fit needs, not whether the patch strays.
fn patch_size_probe(
    mesh: &ConnectedMesh,
    faces: &[usize],
    truth_axis: ([f64; 3], [f64; 3]),
    is_cone: bool,
    model: &str,
    ref_idx: usize,
    kind: &str,
) {
    const SIZES: [usize; 10] = [2, 3, 4, 6, 8, 12, 16, 24, 32, 64];
    const STARTS: usize = 8;
    let member: HashSet<usize> = faces.iter().copied().collect();
    let (_, want_dir) = truth_axis;
    let want = unit(want_dir);
    let ang = |d: [f64; 3]| -> f64 {
        let d = unit(d);
        let c = [
            d[1] * want[2] - d[2] * want[1],
            d[2] * want[0] - d[0] * want[2],
            d[0] * want[1] - d[1] * want[0],
        ];
        norm(c).clamp(0.0, 1.0).asin().to_degrees()
    };
    let stride = (faces.len() / STARTS).max(1);
    for &start in faces.iter().step_by(stride).take(STARTS) {
        // One BFS from the start, recording the axis error each time the
        // accumulated patch reaches the next size on the ladder.
        let mut patch = vec![start];
        let mut seen: HashSet<usize> = HashSet::from([start]);
        let mut frontier = vec![start];
        let mut next_size = 0usize;
        while next_size < SIZES.len() {
            while next_size < SIZES.len() && patch.len() >= SIZES[next_size] {
                let target = SIZES[next_size];
                let sub: Vec<usize> = patch[..target.min(patch.len())].to_vec();
                let got = if is_cone {
                    super::initialize::cone(&sub, mesh, 5e-5)
                } else {
                    super::initialize::cylinder(&sub, mesh)
                };
                if let Some(
                    super::initialize::Initial::Cylinder { axis_dir, .. }
                    | super::initialize::Initial::Cone { axis_dir, .. },
                ) = got
                {
                    println!(
                        "PATCHSIZE\t{model}\t{ref_idx}\t{kind}\t{target}\t{:.6}",
                        ang(axis_dir)
                    );
                }
                next_size += 1;
            }
            // Expand one ring.
            let mut ring = Vec::new();
            for &fi in &frontier {
                let vc = mesh.faces[fi].vertex_count as usize;
                for &ni in &mesh.faces[fi].neighbors[..vc] {
                    if ni < 0 {
                        continue;
                    }
                    let nj = ni as usize;
                    if member.contains(&nj) && seen.insert(nj) {
                        patch.push(nj);
                        ring.push(nj);
                    }
                }
            }
            if ring.is_empty() {
                break;
            }
            frontier = ring;
        }
    }
}

/// Measure the section's cylinder/cone seed rule on one developable reference
/// face: enumerate the chains it admits, test both readings of its coplanarity
/// clause, and fit each chain's closed-form axis against ground truth.
fn seed_rule_probe(
    mesh: &ConnectedMesh,
    faces: &[usize],
    truth_axis: ([f64; 3], [f64; 3]),
    is_cone: bool,
) -> SeedRuleStats {
    let member: HashSet<usize> = faces.iter().copied().collect();
    let mut st = SeedRuleStats {
        pairs: 0,
        six_verts: 0,
        all_coplanar: 0,
        four_coplanar: 0,
        ax_err: Vec::new(),
        ax_err_admitted: Vec::new(),
        selected_err: f64::NAN,
    };
    let (_, want_dir) = truth_axis;
    let want = unit(want_dir);
    let ang = |d: [f64; 3]| -> f64 {
        let d = unit(d);
        let c = [
            d[1] * want[2] - d[2] * want[1],
            d[2] * want[0] - d[0] * want[2],
            d[0] * want[1] - d[1] * want[0],
        ];
        // |sin| form: linear in the error near 0, unlike acos.
        norm(c).clamp(0.0, 1.0).asin().to_degrees()
    };

    // The selection rule Phase 1.1 actually built: widest normal pair over the
    // whole face, re-picked among ruling-aligned facets. The contrast between
    // this and a bare chain is the point of the measurement.
    st.selected_err = match if is_cone {
        super::initialize::cone(faces, mesh, 5e-5)
    } else {
        super::initialize::cylinder(faces, mesh)
    } {
        Some(super::initialize::Initial::Cylinder { axis_dir, .. })
        | Some(super::initialize::Initial::Cone { axis_dir, .. }) => ang(axis_dir),
        _ => f64::NAN,
    };

    // "Coplanar" needs a scale-free bar. 1e-6 of the point set's own extent is
    // five orders below a tessellation step and is a validity floor, not a
    // tuning knob: a fused quad is coplanar to rounding, distinct rulings are
    // separated by the exporter's angular setting.
    let copl_rel = 1.0e-6;
    let (chains, considered) = chains_with_six_vertices(mesh, faces, &member, MAX_SEED_PAIRS);
    st.pairs = considered;
    st.six_verts = chains.len();
    for chain in &chains {
        let vs = vertices_of(mesh, chain);
        let pts: Vec<[f64; 3]> = vs.iter().map(|&v| vertex_pos(mesh, v)).collect();
        let all_cop = planarity(&pts) <= copl_rel;
        if all_cop {
            st.all_coplanar += 1;
        }
        if any_four_coplanar(&pts, copl_rel) {
            st.four_coplanar += 1;
        }
        let got = if is_cone {
            super::initialize::cone(chain, mesh, 5e-5)
        } else {
            super::initialize::cylinder(chain, mesh)
        };
        if let Some(super::initialize::Initial::Cylinder { axis_dir, .. })
        | Some(super::initialize::Initial::Cone { axis_dir, .. }) = got
        {
            let err = ang(axis_dir);
            st.ax_err.push(err);
            if !all_cop {
                st.ax_err_admitted.push(err);
            }
        }
    }
    st
}

fn quantile(v: &mut [f64], q: f64) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let i = ((v.len() - 1) as f64 * q).round() as usize;
    v[i]
}

// ---------------------------------------------------------------------------
// Per-reference-face measurement
// ---------------------------------------------------------------------------

/// Which production surface type stands in for a reference type.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Other,
}

fn kind_of(s: &RefSurface) -> Kind {
    match s {
        RefSurface::Plane { .. } => Kind::Plane,
        RefSurface::Cylinder { .. } => Kind::Cylinder,
        RefSurface::Cone { .. } => Kind::Cone,
        RefSurface::Sphere { .. } => Kind::Sphere,
        RefSurface::Torus { .. } => Kind::Torus,
        RefSurface::Other { .. } => Kind::Other,
    }
}

/// Worst vertex residual of a freshly fitted surface of type `S` over `faces`.
/// `None` when the fitter declines the face set.
fn fit_residual<S: SurfaceType>(mesh: &ConnectedMesh, faces: &[usize]) -> Option<f64> {
    let verts = vertices_of(mesh, faces);
    let fit = S::fit_surface(mesh, faces, &verts, 0)?;
    let mut worst = 0.0_f64;
    for &vi in &verts {
        worst = worst.max(S::vertex_distance(&fit, &vertex_pos(mesh, vi)));
    }
    Some(worst)
}

fn vertices_of(mesh: &ConnectedMesh, faces: &[usize]) -> Vec<usize> {
    let mut v: Vec<usize> = faces
        .iter()
        .flat_map(|&fi| {
            let vc = mesh.faces[fi].vertex_count as usize;
            mesh.faces[fi].vertex_indices[..vc].iter().copied()
        })
        .collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// A flood identical to `compete::flood` except that the neighbour list is
/// walked in reverse. Used only to test traversal-order independence.
fn flood_reversed<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    seed: usize,
    fit: &S::Fit,
    centroids: &[[f64; 3]],
    vtol: f64,
    stol: f64,
) -> Vec<usize> {
    if !compete::face_validates::<S>(mesh, seed, fit, centroids, vtol, stol) {
        return Vec::new();
    }
    let mut seen: HashSet<usize> = HashSet::new();
    seen.insert(seed);
    let mut out = Vec::new();
    let mut queue = std::collections::VecDeque::from([seed]);
    while let Some(fi) = queue.pop_front() {
        out.push(fi);
        let vc = mesh.faces[fi].vertex_count as usize;
        for &ni in mesh.faces[fi].neighbors[..vc].iter().rev() {
            if ni < 0 {
                continue;
            }
            let ni = ni as usize;
            if in_region[ni]
                && !seen.contains(&ni)
                && compete::face_validates::<S>(mesh, ni, fit, centroids, vtol, stol)
            {
                seen.insert(ni);
                queue.push_back(ni);
            }
        }
    }
    out.sort_unstable();
    out
}

/// One reference face's measurements.
struct Row {
    idx: usize,
    kind: &'static str,
    facets: usize,
    area: f64,
    // F — fitting requirement against the reference's own parameters
    v_max: f64,
    e_max: f64,
    c_max: f64,
    i_max: f64,
    /// Edge-midpoint max INCLUDING quad diagonals.
    ed_max: f64,
    // A — angular separability
    dih_in_max: f64,
    dih_out_min: f64,
    /// Edges the pipeline's OWN classifier calls sharp with both faces inside
    /// this one reference face — i.e. places where stage 2.1 cuts a true
    /// surface in half before any fitting happens.
    sharp_in: usize,
    // I — identification: worst vertex residual per candidate type, mm
    res: [f64; 4], // plane, sphere, cyl/cone, torus
    // S — seed independence
    seeds: usize,
    grew: usize,
    exact: usize,
    distinct: usize,
    jacc_min: f64,
    /// Worst fraction of a grown set that lies OUTSIDE the true face set
    /// (leak across a junction).
    over_max: f64,
    /// Worst fraction of the true face set a grown set MISSED (stopping short).
    under_max: f64,
    /// The same three counts restricted to seeds whose every mesh neighbour is
    /// in the same reference face — the section's invariant with an INTERIOR
    /// precondition, which is the version that might actually be true.
    seeds_int: usize,
    grew_int: usize,
    exact_int: usize,
    distinct_int: usize,
    /// How many reference faces share this one's surface parameters — >1 means
    /// the reference split a periodic surface at its seam.
    seam_group: usize,
    /// Seeds whose grown set equals the UNION of that seam group.
    exact_union: usize,
    // competition-level (all four types compete, as the pipeline does)
    c_distinct_types: usize,
    c_correct_type: usize,
    /// `type:count` histogram of what the competition picked, for the seeds
    /// that produced anything.
    c_hist: String,
    // G — growth independence
    order_diffs: usize,
    oracle_same: bool,
}

/// Jaccard index of two sorted face sets.
fn jaccard(a: &[usize], b: &HashSet<usize>) -> f64 {
    if a.is_empty() && b.is_empty() {
        return 1.0;
    }
    let inter = a.iter().filter(|f| b.contains(f)).count();
    let union = a.len() + b.len() - inter;
    if union == 0 {
        1.0
    } else {
        inter as f64 / union as f64
    }
}

/// Run the whole report. `jet`, `face_region` and `regions` come from stage 2.2.
#[allow(clippy::too_many_arguments)]
pub fn report(
    config: &Config,
    mesh: &ConnectedMesh,
    model: &ReferenceModel,
    labels: &[FaceLabel],
    jet: &[JetField],
    edges: &EdgeClassification,
    face_region: &[usize],
    regions: &[super::regions::Region],
    viz: Option<&crate::viz::VizSender>,
) {
    let model_name = config
        .input_stl
        .rsplit('/')
        .next()
        .unwrap_or(&config.input_stl)
        .trim_end_matches(".stl")
        .to_string();
    let vtol = config.vertex_tolerance_mm;
    let stol = config.surface_tolerance_mm;
    let quantum = float32_position_bound_mm(mesh);
    let centroids: Vec<[f64; 3]> = (0..mesh.faces.len())
        .map(|fi| face_centroid(&mesh.faces[fi], &mesh.vertices))
        .collect();

    // Group facets by reference face.
    let mut by_ref: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut unlabeled = 0usize;
    for (fi, lab) in labels.iter().enumerate() {
        match lab.face() {
            Some(r) => by_ref.entry(r).or_default().push(fi),
            None => unlabeled += 1,
        }
    }
    for v in by_ref.values_mut() {
        v.sort_unstable();
    }

    // Reverse lookup for the across-boundary dihedral.
    let mut ref_of: Vec<Option<usize>> = vec![None; mesh.faces.len()];
    for (&r, faces) in &by_ref {
        for &fi in faces {
            ref_of[fi] = Some(r);
        }
    }

    let mut keys: Vec<usize> = by_ref.keys().copied().collect();
    keys.sort_unstable();

    // A periodic surface is normally split across several reference faces at
    // its seam (OCCT gives both halves the same underlying geometry), and a
    // fit that spans the seam is RIGHT even though it does not equal any one
    // reference face. Group faces by identical surface parameters so the
    // seed-independence test can score against the union as well.
    let mut union_of: HashMap<usize, Vec<usize>> = HashMap::new();
    for &a in &keys {
        let mut grp = Vec::new();
        for &b in &keys {
            if model.surfaces[a] == model.surfaces[b] {
                grp.push(b);
            }
        }
        union_of.insert(a, grp);
    }

    println!(
        "#INV\tmodel\tref\tkind\tfacets\tarea\tvmax\temax\tcmax\timax\tedmax\tdih_in\tdih_out\tsharp_in\t\
         res_plane\tres_sphere\tres_cyl\tres_torus\tseeds\tgrew\texact\tdistinct\tjacc_min\tover_max\tunder_max\tseeds_int\tgrew_int\texact_int\tdistinct_int\tseam_group\texact_union\t\
         ctypes\tccorrect\tchist\torder_diffs\toracle\tvtol\tstol\tquantum"
    );

    println!(
        "#BASIN\tmodel\tref\tkind\tface_facets\tshape\tseed_facets\tconnected\tself_ok\t\
         verdict\trounds\tsweep_deg\taxial_rel\tjaccard\tresult_facets\tover_n\tshort_n\t\
         over_refs\tover_adjacent\tover_better\tover_scored"
    );
    println!(
        "#ANCHOR\tmodel\tref\tkind\tface_facets\tanchor\tattempts\tfacet_ok\twinner_shape\t\
         last_verdict"
    );
    println!("#PATCHSIZE\tmodel\tref\tkind\tfacets\taxis_err_deg");
    println!(
        "#SEEDRULE\tmodel\tref\tkind\tpairs\tsix_verts\tall_coplanar\tfour_coplanar\t\
         ax_p50\tax_p95\tax_max\tn_admitted\tadm_p50\tadm_max\tselected_err"
    );

    let mut seedrule_totals = (0usize, 0usize, 0usize, 0usize, 0usize);
    let mut seedrule_worst = 0.0_f64;
    let mut seedrule_worst_adm = 0.0_f64;
    let mut seedrule_selected_worst = 0.0_f64;
    let mut viz_quit = false;
    let mut rows: Vec<Row> = Vec::new();
    for &r in &keys {
        let faces = &by_ref[&r];
        let surf = &model.surfaces[r];
        let kind = kind_of(surf);
        let area: f64 = faces
            .iter()
            .map(|&fi| face_area(&mesh.faces[fi], &mesh.vertices))
            .sum();

        // --- F: fitting requirement against ground truth -------------------
        let (mut v_max, mut e_max, mut c_max, mut i_max) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
        let mut ed_max = 0.0_f64;
        let mut probes_ok = true;
        for &fi in faces {
            match facet_probes(mesh, fi, surf) {
                Some((v, e, c, i, ed)) => {
                    v_max = v_max.max(v);
                    e_max = e_max.max(e);
                    c_max = c_max.max(c);
                    i_max = i_max.max(i);
                    ed_max = ed_max.max(ed);
                }
                None => {
                    probes_ok = false;
                    break;
                }
            }
        }
        if !probes_ok {
            v_max = f64::NAN;
            e_max = f64::NAN;
            c_max = f64::NAN;
            i_max = f64::NAN;
            ed_max = f64::NAN;
        }

        // --- A: angular separability ---------------------------------------
        let mut dih_in_max = 0.0_f64;
        let mut dih_out_min = f64::INFINITY;
        let mut sharp_in = 0usize;
        for &fi in faces {
            let vc = mesh.faces[fi].vertex_count as usize;
            for (e_idx, &ni) in mesh.faces[fi].neighbors[..vc].iter().enumerate() {
                if ni < 0 {
                    continue;
                }
                let ni = ni as usize;
                let d = dihedral_deg(mesh, fi, ni);
                match ref_of[ni] {
                    Some(o) if o == r => {
                        dih_in_max = dih_in_max.max(d);
                        if fi < ni
                            && edges.class[fi][e_idx] == EdgeClass::Sharp
                        {
                            sharp_in += 1;
                        }
                    }
                    Some(_) => dih_out_min = dih_out_min.min(d),
                    None => {}
                }
            }
        }

        // --- I: all four types on the true support -------------------------
        let res = [
            fit_residual::<PlanarSurface>(mesh, faces).unwrap_or(f64::NAN),
            fit_residual::<SphericalSurface>(mesh, faces).unwrap_or(f64::NAN),
            fit_residual::<ConicalSurface>(mesh, faces).unwrap_or(f64::NAN),
            fit_residual::<ToroidalSurface>(mesh, faces).unwrap_or(f64::NAN),
        ];

        // --- E: the section's stated seed rule, on developables -------------
        if matches!(kind, Kind::Cylinder | Kind::Cone) {
            if let Some(ax) = surf.axis() {
                if std::env::var("BREPPER_INV_PATCHSIZE").is_ok() {
                patch_size_probe(
                    mesh,
                    faces,
                    ax,
                    kind == Kind::Cone,
                    &model_name,
                    r,
                    surf.type_name(),
                );
                }
                let st = seed_rule_probe(mesh, faces, ax, kind == Kind::Cone);
                let p50 = quantile(&mut st.ax_err.clone(), 0.50);
                let p95 = quantile(&mut st.ax_err.clone(), 0.95);
                let amax = st.ax_err.iter().cloned().fold(0.0_f64, f64::max);
                let ap50 = quantile(&mut st.ax_err_admitted.clone(), 0.50);
                let amax2 = st
                    .ax_err_admitted
                    .iter()
                    .cloned()
                    .fold(0.0_f64, f64::max);
                println!(
                    "SEEDRULE\t{model}\t{r}\t{k}\t{pairs}\t{six}\t{allc}\t{fourc}\t\
                     {p50:.6}\t{p95:.6}\t{amax:.6}\t{n_adm}\t{ap50:.6}\t{amax2:.6}\t{sel:.6}",
                    model = model_name,
                    r = r,
                    k = surf.type_name(),
                    pairs = st.pairs,
                    six = st.six_verts,
                    allc = st.all_coplanar,
                    fourc = st.four_coplanar,
                    n_adm = st.ax_err_admitted.len(),
                    sel = st.selected_err,
                );
                seedrule_totals.0 += st.pairs;
                seedrule_totals.1 += st.six_verts;
                seedrule_totals.2 += st.all_coplanar;
                seedrule_totals.3 += st.four_coplanar;
                seedrule_totals.4 += st.ax_err.len();
                seedrule_worst = seedrule_worst.max(amax);
                        seedrule_worst_adm = seedrule_worst_adm.max(amax2);
                if st.selected_err.is_finite() {
                    seedrule_selected_worst = seedrule_selected_worst.max(st.selected_err);
                }
            }
        }

        // --- S/G: seed and growth independence -----------------------------
        let mut row = Row {
            idx: r,
            kind: surf.type_name(),
            facets: faces.len(),
            area,
            v_max,
            e_max,
            c_max,
            i_max,
            ed_max,
            dih_in_max,
            dih_out_min,
            sharp_in,
            res,
            seeds: 0,
            grew: 0,
            exact: 0,
            distinct: 0,
            jacc_min: f64::NAN,
            over_max: 0.0,
            under_max: 0.0,
            seeds_int: 0,
            grew_int: 0,
            exact_int: 0,
            distinct_int: 0,
            seam_group: 1,
            exact_union: 0,
            c_distinct_types: 0,
            c_correct_type: 0,
            c_hist: String::new(),
            order_diffs: 0,
            oracle_same: false,
        };

        if faces.len() >= MIN_FACETS_FOR_SEED_TEST && kind != Kind::Other {
            let truth: HashSet<usize> = faces.iter().copied().collect();

            // `--viz-seeds`: walk seed sets interactively (failing ones by
            // default), showing what each grew to and why it stopped.
            if config.viz_seeds && !viz_quit {
                {
                    let vs = viz;
                    // By default the flood is limited ONLY by validation: the
                    // sharp-edge partition is disabled here so that a seed
                    // cannot fail merely because stage 2.1 cut its surface in
                    // two. `--viz-seed-confine-region` restores the confinement.
                    let mut inreg = vec![true; mesh.faces.len()];
                    if config.viz_seed_confine_region {
                        inreg = vec![false; mesh.faces.len()];
                        let reg0 = face_region[faces[0]];
                        if reg0 < regions.len() {
                            for &fi in &regions[reg0].faces {
                                inreg[fi] = true;
                            }
                        }
                    }
                    let keep = match kind {
                        Kind::Plane => viz_face::<PlanarSurface>(
                            mesh, &inreg, faces, &truth, &centroids, config, vs, &model_name, r,
                            surf.type_name(), kind, surf, quantum,
                        ),
                        Kind::Cylinder | Kind::Cone => viz_face::<ConicalSurface>(
                            mesh, &inreg, faces, &truth, &centroids, config, vs, &model_name, r,
                            surf.type_name(), kind, surf, quantum,
                        ),
                        Kind::Sphere => viz_face::<SphericalSurface>(
                            mesh, &inreg, faces, &truth, &centroids, config, vs, &model_name, r,
                            surf.type_name(), kind, surf, quantum,
                        ),
                        Kind::Torus => viz_face::<ToroidalSurface>(
                            mesh, &inreg, faces, &truth, &centroids, config, vs, &model_name, r,
                            surf.type_name(), kind, surf, quantum,
                        ),
                        Kind::Other => true,
                    };
                    if !keep {
                        viz_quit = true;
                    }
                }
            }

            // Basin-of-attraction sweep for the CLEAN operator (seed set = the
            // whole initial fit support, no seed_patch/condition_patch).
            if std::env::var("BREPPER_INV_BASIN").is_ok() {
                let reg0 = face_region[faces[0]];
                let mut inreg = vec![false; mesh.faces.len()];
                if reg0 < regions.len() {
                    for &fi in &regions[reg0].faces {
                        inreg[fi] = true;
                    }
                }
                let ax = surf.axis();
                match kind {
                    Kind::Plane => basin_probe::<PlanarSurface>(
                        mesh, &inreg, faces, &truth, &centroids, config, ax, &ref_of,
                        &model.surfaces, &model_name, r, surf.type_name(),
                    ),
                    Kind::Cylinder | Kind::Cone => basin_probe::<ConicalSurface>(
                        mesh, &inreg, faces, &truth, &centroids, config, ax, &ref_of,
                        &model.surfaces, &model_name, r, surf.type_name(),
                    ),
                    Kind::Sphere => basin_probe::<SphericalSurface>(
                        mesh, &inreg, faces, &truth, &centroids, config, ax, &ref_of,
                        &model.surfaces, &model_name, r, surf.type_name(),
                    ),
                    Kind::Torus => basin_probe::<ToroidalSurface>(
                        mesh, &inreg, faces, &truth, &centroids, config, ax, &ref_of,
                        &model.surfaces, &model_name, r, surf.type_name(),
                    ),
                    Kind::Other => {}
                }
            }
            let truth_union: HashSet<usize> = union_of[&r]
                .iter()
                .flat_map(|g| by_ref[g].iter().copied())
                .collect();
            row.seam_group = union_of[&r].len();
            // The smooth region the face lives in — the arena competition uses.
            let reg = face_region[faces[0]];
            let mut in_region = vec![false; mesh.faces.len()];
            if reg < regions.len() {
                for &fi in &regions[reg].faces {
                    in_region[fi] = true;
                }
            }

            let nseeds = seeds_per_face();
            let stride = (faces.len() / nseeds).max(1);
            let seeds: Vec<usize> = faces.iter().copied().step_by(stride).take(nseeds).collect();
            row.seeds = seeds.len();

            let mut outcomes: HashSet<Vec<usize>> = HashSet::new();
            let mut outcomes_int: HashSet<Vec<usize>> = HashSet::new();
            // A seed is INTERIOR when every mesh neighbour it has is labelled
            // to the same reference face — so nothing about its own
            // neighbourhood is ambiguous.
            let is_interior = |fi: usize| -> bool {
                let vc = mesh.faces[fi].vertex_count as usize;
                mesh.faces[fi].neighbors[..vc].iter().all(|&ni| {
                    ni >= 0 && ref_of[ni as usize] == Some(r)
                })
            };
            let mut jacc_min = f64::INFINITY;
            let do_compete = compete_enabled();
            let mut ctypes: HashSet<&'static str> = HashSet::new();
            let mut chist: HashMap<&'static str, usize> = HashMap::new();
            let want_type = match kind {
                Kind::Plane => "plane",
                Kind::Cylinder => "cyl",
                Kind::Cone => "cone",
                Kind::Sphere => "sphere",
                Kind::Torus => "torus",
                Kind::Other => "",
            };
            let flat_eps = {
                // Mirror compete_region's region-scale flatness floor.
                let mut min = [f64::INFINITY; 3];
                let mut max = [f64::NEG_INFINITY; 3];
                if reg < regions.len() {
                    for &vi in &regions[reg].vertices {
                        let v = &mesh.vertices[vi];
                        for (k, c) in [v.x, v.y, v.z].into_iter().enumerate() {
                            min[k] = min[k].min(c);
                            max[k] = max[k].max(c);
                        }
                    }
                }
                let extent = (max[0] - min[0]).max(max[1] - min[1]).max(max[2] - min[2]);
                if extent.is_finite() && extent > 0.0 {
                    1.0 / (2.0 * extent)
                } else {
                    1e-6
                }
            };

            for &seed in &seeds {
                // Single-type grow: the invariant as the section states it.
                let grown: Option<Vec<usize>> = match kind {
                    Kind::Plane => compete::grow::<PlanarSurface>(
                        mesh, &in_region, seed, jet, &centroids, config, 1,
                    )
                    .map(|(_, _, f)| f),
                    Kind::Cylinder | Kind::Cone => compete::grow::<ConicalSurface>(
                        mesh, &in_region, seed, jet, &centroids, config, MIN_CYLINDER_FACES,
                    )
                    .map(|(_, _, f)| f),
                    Kind::Sphere => compete::grow::<SphericalSurface>(
                        mesh, &in_region, seed, jet, &centroids, config, MIN_SPHERE_FACES,
                    )
                    .map(|(_, _, f)| f),
                    Kind::Torus => compete::grow::<ToroidalSurface>(
                        mesh, &in_region, seed, jet, &centroids, config, MIN_TORUS_FACES,
                    )
                    .map(|(_, _, f)| f),
                    Kind::Other => None,
                };
                let interior = is_interior(seed);
                if interior {
                    row.seeds_int += 1;
                }
                if let Some(g) = grown {
                    row.grew += 1;
                    if interior {
                        row.grew_int += 1;
                        if g.len() == truth.len() && g.iter().all(|f| truth.contains(f)) {
                            row.exact_int += 1;
                        }
                        outcomes_int.insert(g.clone());
                    }
                    if g.len() == truth.len() && g.iter().all(|f| truth.contains(f)) {
                        row.exact += 1;
                    }
                    if g.len() == truth_union.len() && g.iter().all(|f| truth_union.contains(f)) {
                        row.exact_union += 1;
                    }
                    jacc_min = jacc_min.min(jaccard(&g, &truth));
                    let inside = g.iter().filter(|f| truth.contains(f)).count();
                    row.over_max = row
                        .over_max
                        .max((g.len() - inside) as f64 / g.len().max(1) as f64);
                    row.under_max = row
                        .under_max
                        .max((truth.len() - inside) as f64 / truth.len().max(1) as f64);
                    outcomes.insert(g);
                }

                // Competition-level: what the pipeline would actually pick.
                if !do_compete {
                    continue;
                }
                if let Some((hyp, _, _)) = compete::compete_from_seed(
                    mesh, &in_region, seed, jet, &centroids, config, flat_eps, false,
                ) {
                    ctypes.insert(hyp.kind_name());
                    *chist.entry(hyp.kind_name()).or_insert(0) += 1;
                    if hyp.kind_name() == want_type {
                        row.c_correct_type += 1;
                    }
                }
            }
            row.distinct = outcomes.len();
            row.distinct_int = outcomes_int.len();
            row.jacc_min = if jacc_min.is_finite() { jacc_min } else { f64::NAN };
            row.c_distinct_types = ctypes.len();
            let mut hv: Vec<(&str, usize)> = chist.into_iter().collect();
            hv.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
            row.c_hist = hv
                .iter()
                .map(|(k, v)| format!("{k}:{v}"))
                .collect::<Vec<_>>()
                .join(",");

            // --- G: traversal-order independence of the flood ---------------
            // Fit the true support once, then flood from each seed in both
            // neighbour orders and compare the resulting sets.
            let verts = vertices_of(mesh, faces);
            macro_rules! order_test {
                ($S:ty) => {{
                    if let Some(mut fit) = <$S>::fit_surface(mesh, faces, &verts, 0) {
                        <$S>::determine_convexity(mesh, faces, &mut fit);
                        let mut diffs = 0;
                        for &seed in &seeds {
                            let a = compete::flood::<$S>(
                                mesh, &in_region, seed, &fit, &centroids, vtol, stol,
                            );
                            let b = flood_reversed::<$S>(
                                mesh, &in_region, seed, &fit, &centroids, vtol, stol,
                            );
                            if a != b {
                                diffs += 1;
                            }
                        }
                        row.order_diffs = diffs;
                    }
                }};
            }
            match kind {
                Kind::Plane => order_test!(PlanarSurface),
                Kind::Cylinder | Kind::Cone => order_test!(ConicalSurface),
                Kind::Sphere => order_test!(SphericalSurface),
                Kind::Torus => order_test!(ToroidalSurface),
                Kind::Other => {}
            }

            // --- G: start independence of the refit fixpoint ----------------
            // Feed `grow_from` the ENTIRE true face set as its initial patch —
            // the most different start possible from a 32-face seed patch — and
            // ask whether the fixpoint is the same one the seed reached.
            let seed0 = seeds[0];
            let oracle: Option<Vec<usize>> = match kind {
                Kind::Plane => compete::grow_from::<PlanarSurface>(
                    mesh, &in_region, seed0, faces.clone(), &centroids, config, 1,
                )
                .map(|(_, _, f)| f),
                Kind::Cylinder | Kind::Cone => compete::grow_from::<ConicalSurface>(
                    mesh, &in_region, seed0, faces.clone(), &centroids, config,
                    MIN_CYLINDER_FACES,
                )
                .map(|(_, _, f)| f),
                Kind::Sphere => compete::grow_from::<SphericalSurface>(
                    mesh, &in_region, seed0, faces.clone(), &centroids, config, MIN_SPHERE_FACES,
                )
                .map(|(_, _, f)| f),
                Kind::Torus => compete::grow_from::<ToroidalSurface>(
                    mesh, &in_region, seed0, faces.clone(), &centroids, config, MIN_TORUS_FACES,
                )
                .map(|(_, _, f)| f),
                Kind::Other => None,
            };
            row.oracle_same = match oracle {
                Some(o) => o.len() == truth.len() && o.iter().all(|f| truth.contains(f)),
                None => false,
            };
        }

        println!(
            "INV\t{model}\t{idx}\t{kind}\t{facets}\t{area:.6}\t{vmax:.6e}\t{emax:.6e}\t\
             {cmax:.6e}\t{imax:.6e}\t{edmax:.6e}\t{din:.4}\t{dout:.4}\t{sharpin}\t{r0:.6e}\t{r1:.6e}\t{r2:.6e}\t\
             {r3:.6e}\t{seeds}\t{grew}\t{exact}\t{distinct}\t{jacc:.4}\t{over:.4}\t{under:.4}\t{si}\t{gi}\t{ei}\t{di}\t{sg}\t{eu}\t{ct}\t{cc}\t{chist}\t{od}\t\
             {oracle}\t{vtol:.6e}\t{stol:.6e}\t{q:.6e}",
            model = model_name,
            idx = row.idx,
            kind = row.kind,
            facets = row.facets,
            area = row.area,
            vmax = row.v_max,
            emax = row.e_max,
            cmax = row.c_max,
            imax = row.i_max,
            edmax = row.ed_max,
            din = row.dih_in_max,
            sharpin = row.sharp_in,
            dout = if row.dih_out_min.is_finite() { row.dih_out_min } else { -1.0 },
            r0 = row.res[0],
            r1 = row.res[1],
            r2 = row.res[2],
            r3 = row.res[3],
            seeds = row.seeds,
            grew = row.grew,
            exact = row.exact,
            distinct = row.distinct,
            jacc = row.jacc_min,
            over = row.over_max,
            under = row.under_max,
            si = row.seeds_int,
            gi = row.grew_int,
            ei = row.exact_int,
            di = row.distinct_int,
            sg = row.seam_group,
            eu = row.exact_union,
            ct = row.c_distinct_types,
            cc = row.c_correct_type,
            chist = if row.c_hist.is_empty() { "-" } else { &row.c_hist },
            od = row.order_diffs,
            oracle = if row.oracle_same { 1 } else { 0 },
            vtol = vtol,
            stol = stol,
            q = quantum,
        );
        rows.push(row);
    }

    if seedrule_totals.0 > 0 {
        eprintln!(
            "    E  seed rule on developables: {} distinct chains (2\u{2013}4 facets) — {} have exactly 6 vertices \
             ({:.0}% of chains); of THOSE, {} are wholly coplanar ({:.1}%) and {} have SOME four coplanar ({:.0}%)",
            seedrule_totals.0,
            seedrule_totals.1,
            pct(seedrule_totals.1, seedrule_totals.0),
            seedrule_totals.2,
            pct(seedrule_totals.2, seedrule_totals.1),
            seedrule_totals.3,
            pct(seedrule_totals.3, seedrule_totals.1),
        );
        eprintln!(
            "    E  closed-form axis from a bare pair: worst {seedrule_worst:.4}\u{00b0} over {} \
             solvable chains; worst {seedrule_worst_adm:.4}\u{00b0} among the non-coplanar ones — \
             against {seedrule_selected_worst:.4}\u{00b0} for the whole-face SELECTED pair",
            seedrule_totals.4,
        );
    }
    summarize(&model_name, &rows, unlabeled, mesh.faces.len(), vtol, stol, quantum);
}

fn pct(a: usize, b: usize) -> f64 {
    if b == 0 {
        0.0
    } else {
        100.0 * a as f64 / b as f64
    }
}

#[allow(clippy::too_many_arguments)]
fn summarize(
    model: &str,
    rows: &[Row],
    unlabeled: usize,
    total_faces: usize,
    vtol: f64,
    stol: f64,
    quantum: f64,
) {
    let analytic: Vec<&Row> = rows.iter().filter(|r| r.v_max.is_finite()).collect();
    eprintln!("  Invariant report ({model}):");
    eprintln!(
        "    reference faces {} ({} analytic, {} freeform/unmeasurable); \
         {unlabeled} of {total_faces} facets unlabeled (seam-ambiguous)",
        rows.len(),
        analytic.len(),
        rows.len() - analytic.len(),
    );
    eprintln!(
        "    tolerances: vtol {vtol:.3e}  stol {stol:.3e}  f32 bound {quantum:.3e} mm"
    );

    if !analytic.is_empty() {
        let worst_v = analytic.iter().map(|r| r.v_max).fold(0.0_f64, f64::max);
        let worst_e = analytic.iter().map(|r| r.e_max).fold(0.0_f64, f64::max);
        let worst_c = analytic.iter().map(|r| r.c_max).fold(0.0_f64, f64::max);
        let worst_i = analytic.iter().map(|r| r.i_max).fold(0.0_f64, f64::max);
        let v_fail = analytic.iter().filter(|r| r.v_max > vtol).count();
        let e_fail = analytic.iter().filter(|r| r.e_max > stol).count();
        let c_fail = analytic.iter().filter(|r| r.c_max > stol).count();
        eprintln!(
            "    F  vertex max {worst_v:.3e} ({} q) — {v_fail} face(s) over vtol",
            (worst_v / quantum).round() as i64
        );
        eprintln!(
            "    F  edge-mid max {worst_e:.3e} ({e_fail} over stol); centroid max {worst_c:.3e} \
             ({c_fail} over stol); interior max {worst_i:.3e}"
        );
        // The 4/3 prediction: the whole-facet max lies in [edge-mid max, 4/3 × it].
        let ratios: Vec<f64> = analytic
            .iter()
            .filter(|r| r.ed_max > 0.0)
            .map(|r| r.i_max.max(r.c_max).max(r.ed_max) / r.ed_max)
            .collect();
        if !ratios.is_empty() {
            let hi = ratios.iter().cloned().fold(0.0_f64, f64::max);
            let over = ratios.iter().filter(|&&x| x > 4.0 / 3.0 + 1e-9).count();
            eprintln!(
                "    F  facet-max / edge-mid-max: worst {hi:.4} (predicted \u{2264} 1.3333); \
                 {over} face(s) over"
            );
        }
        let cratios: Vec<f64> = analytic
            .iter()
            .filter(|r| r.c_max > 0.0)
            .map(|r| r.i_max.max(r.e_max).max(r.c_max) / r.c_max)
            .collect();
        if !cratios.is_empty() {
            let hi = cratios.iter().cloned().fold(0.0_f64, f64::max);
            eprintln!("    F  facet-max / centroid: worst {hi:.4} (centroid under-reads by this)");
        }
    }

    let split = rows.iter().filter(|r| r.sharp_in > 0).count();
    let split_edges: usize = rows.iter().map(|r| r.sharp_in).sum();
    eprintln!(
        "    A  stage 2.1 calls {split_edges} edge(s) SHARP strictly inside a reference face, \
         cutting {split} true surface(s) before any fitting happens"
    );
    let sep: Vec<&Row> = rows
        .iter()
        .filter(|r| r.dih_out_min.is_finite() && r.dih_in_max > 0.0)
        .collect();
    let inseparable = sep.iter().filter(|r| r.dih_in_max >= r.dih_out_min).count();
    if !sep.is_empty() {
        let worst = sep
            .iter()
            .map(|r| r.dih_in_max / r.dih_out_min.max(1e-12))
            .fold(0.0_f64, f64::max);
        eprintln!(
            "    A  {inseparable} of {} faces have an interior dihedral \u{2265} their smallest \
             boundary dihedral (worst ratio {worst:.1}\u{00d7}) — no angular tolerance separates those",
            sep.len()
        );
    }

    // Identification: how often does a WRONG type also validate the true support?
    let names = ["plane", "sphere", "cyl/cone", "torus"];
    let want = |r: &Row| match r.kind {
        "plane" => 0usize,
        "sphere" => 1,
        "cylinder" | "cone" => 2,
        "torus" => 3,
        _ => usize::MAX,
    };
    let mut impostors = vec![0usize; 4];
    let mut own_fail = 0usize;
    let mut n_typed = 0usize;
    for r in rows {
        let w = want(r);
        if w == usize::MAX {
            continue;
        }
        n_typed += 1;
        if r.res[w].is_nan() || r.res[w] > vtol {
            own_fail += 1;
        }
        for (k, imp) in impostors.iter_mut().enumerate() {
            if k != w && r.res[k].is_finite() && r.res[k] <= vtol {
                *imp += 1;
            }
        }
    }
    if n_typed > 0 {
        eprintln!(
            "    I  {own_fail} of {n_typed} reference faces are NOT explained at vtol by a fit of \
             their OWN type"
        );
        eprintln!(
            "    I  wrong types that also validate the true support: {}",
            names
                .iter()
                .zip(&impostors)
                .filter(|(_, &c)| c > 0)
                .map(|(n, c)| format!("{n} {c}"))
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    let tested: Vec<&Row> = rows.iter().filter(|r| r.seeds > 0).collect();
    if !tested.is_empty() {
        let seeds: usize = tested.iter().map(|r| r.seeds).sum();
        let grew: usize = tested.iter().map(|r| r.grew).sum();
        let exact: usize = tested.iter().map(|r| r.exact).sum();
        let stable = tested.iter().filter(|r| r.distinct <= 1).count();
        let ccorrect: usize = tested.iter().map(|r| r.c_correct_type).sum();
        let ctype_stable = tested.iter().filter(|r| r.c_distinct_types <= 1).count();
        let worst_j = tested
            .iter()
            .filter(|r| r.jacc_min.is_finite())
            .map(|r| r.jacc_min)
            .fold(1.0_f64, f64::min);
        let exact_u: usize = tested.iter().map(|r| r.exact_union).sum();
        let seamed = tested.iter().filter(|r| r.seam_group > 1).count();
        eprintln!(
            "    S  {} faces seed-tested ({seamed} of them a seam-split half), {seeds} seeds: \
             {grew} grew ({:.1}%), {exact} recovered the true face set ({:.1}%), {exact_u} \
             recovered its whole seam group ({:.1}%)",
            tested.len(),
            pct(grew, seeds),
            pct(exact, seeds),
            pct(exact_u, seeds),
        );
        let worst_over = tested.iter().map(|r| r.over_max).fold(0.0_f64, f64::max);
        let worst_under = tested.iter().map(|r| r.under_max).fold(0.0_f64, f64::max);
        eprintln!(
            "    S  {stable} of {} faces gave ONE outcome across all seeds; worst Jaccard vs truth \
             {worst_j:.3} (worst leak-out {:.1}%, worst stop-short {:.1}%)",
            tested.len(),
            100.0 * worst_over,
            100.0 * worst_under,
        );
        eprintln!(
            "    S  competition: {ccorrect} of {seeds} seeds picked the reference's own type \
             ({:.1}%); {ctype_stable} of {} faces agreed on a single type",
            pct(ccorrect, seeds),
            tested.len()
        );
        let si: usize = tested.iter().map(|r| r.seeds_int).sum();
        let gi: usize = tested.iter().map(|r| r.grew_int).sum();
        let ei: usize = tested.iter().map(|r| r.exact_int).sum();
        let stable_int = tested
            .iter()
            .filter(|r| r.seeds_int > 0 && r.distinct_int <= 1 && r.grew_int == r.seeds_int)
            .count();
        let with_int = tested.iter().filter(|r| r.seeds_int > 0).count();
        eprintln!(
            "    S  INTERIOR seeds only ({si} of {seeds}): {gi} grew ({:.1}%), {ei} exact \
             ({:.1}%); {stable_int} of {with_int} faces fully agreed",
            pct(gi, si),
            pct(ei, si),
        );
        let od: usize = tested.iter().map(|r| r.order_diffs).sum();
        let oracle = tested.iter().filter(|r| r.oracle_same).count();
        eprintln!(
            "    G  flood traversal order changed the result on {od} of {seeds} seeds; \
             oracle-start fixpoint equals truth on {oracle} of {} faces",
            tested.len()
        );
    }
}
