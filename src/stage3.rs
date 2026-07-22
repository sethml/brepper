//! Stage 3: Surface Reconstruction
//!
//! 3.1: Create OCCT surface objects from hypotheses and build adjacency graph
//! 3.2: Detect and create tangency relationships
//! 3.3: Create OCCT edge wires (surface-surface intersections)
//! 3.4: Create OCCT faces (surface bounded by wires)
//! 3.5: Construct shells from connected faces
//! 3.6: Construct solids from shells

use crate::config::Config;
use crate::stage1::{ConnectedMesh, MeshVertex};
use crate::stage2::{SelectedSurface, Stage2Output};
use crate::stage1;
use opencascade_sys::{
    b_rep, b_rep_adaptor, b_rep_bnd_lib, b_rep_builder_api, b_rep_check, b_rep_extrema,
    b_rep_g_prop, b_rep_lib,
    bnd, convert, extrema,
    geom, geom2d, geom2d_api, geom_abs, geom_api, geom_convert, g_prop, gp, message,
    shape_analysis, shape_fix,
    n_collection, top_abs, top_exp, top_loc, topo_ds, OwnedPtr,
};
use std::cell::RefCell;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::time::Instant;
use std::error::Error;
use std::fmt::{Display, Formatter};

// ---------------------------------------------------------------------------
// Data structures — Stage 3 output
// ---------------------------------------------------------------------------

/// Descriptor for a reconstructed B-Rep face, corresponding to one selected surface hypothesis.
pub struct FaceDescriptor {
    /// Index of the selected surface (into Stage2Output.selected_surfaces).
    pub selected_surface_idx: usize,
    /// OCCT surface object (Geom_Plane, Geom_CylindricalSurface, or Geom_SphericalSurface),
    /// type-erased as HandleGeomSurface.
    pub surface: OwnedPtr<geom::Handle_Geom_Surface>,
    // The TopoDS_Face objects built from these descriptors live in
    // Stage3Output.faces (parallel to face_descriptors), populated in 3.4.
    /// Indices of adjacent FaceDescriptors, ordered topologically around this face's boundary.
    /// A face connecting to this one via multiple edges will appear multiple times.
    pub adjacent_faces: Vec<usize>,
    /// ReconEdge indices, one per adjacent face. `edge_indices[i]` is the edge between
    /// this face and `adjacent_faces[i]`.
    pub edge_indices: Vec<usize>,
    /// BRepVertex indices at corners. `vertex_indices[i]` is the vertex between
    /// `edge_indices[i]` and `edge_indices[(i+1) % N]`.
    pub vertex_indices: Vec<usize>,
}

impl std::fmt::Debug for FaceDescriptor {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FaceDescriptor")
            .field("selected_surface_idx", &self.selected_surface_idx)
            .field("surface", &"<HandleGeomSurface>")
            .field("adjacent_faces", &self.adjacent_faces)
            .field("edge_indices", &self.edge_indices)
            .field("vertex_indices", &self.vertex_indices)
            .finish()
    }
}

/// Descriptor for an edge (boundary segment) between two adjacent reconstructed faces.
pub struct ReconEdge {
    /// Indices of the two adjacent FaceDescriptors.
    pub face_indices: [usize; 2],
    /// Indices of the two BRepVertices at the endpoints.
    /// For closed-loop edges (no vertices), both are `usize::MAX`.
    pub vertex_indices: [usize; 2],
    /// 3D intersection curve, trimmed to vertex endpoints. Populated in stage 3.3.
    pub curve_3d: Option<OwnedPtr<geom::Handle_Geom_Curve>>,
    /// Whether the adjacent surfaces are tangent along this edge.
    pub tangent: bool,
    /// Mesh vertex indices along this boundary, ordered along the boundary.
    pub mesh_boundary_vertices: Vec<usize>,
}

impl std::fmt::Debug for ReconEdge {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ReconEdge")
            .field("face_indices", &self.face_indices)
            .field("vertex_indices", &self.vertex_indices)
            .field("curve_3d", &self.curve_3d.as_ref().map(|_| "<HandleGeomCurve>"))
            .field("tangent", &self.tangent)
            .field("mesh_boundary_vertices", &self.mesh_boundary_vertices)
            .finish()
    }
}

/// Vertex where three or more reconstructed faces meet.
#[derive(Debug)]
pub struct BRepVertex {
    /// Indices of adjacent FaceDescriptors, in topological order around the vertex.
    pub adjacent_faces: Vec<usize>,
    /// Indices of adjacent ReconEdges, in topological order.
    pub adjacent_edges: Vec<usize>,
    /// 3D point location of this vertex.
    pub point: [f64; 3],
}

/// The output of Stage 3: the fully reconstructed B-Rep topology.
pub struct Stage3Output {
    /// Stage 2 output, preserved for access to hypotheses and mesh.
    pub stage2: Stage2Output,
    /// Reconstructed face descriptors, one per selected surface.
    pub face_descriptors: Vec<FaceDescriptor>,
    /// Edges connecting adjacent faces.
    pub edges: Vec<ReconEdge>,
    /// Vertices where faces and edges meet.
    pub vertices: Vec<BRepVertex>,
    /// OCCT face objects created in stage 3.4 (one per face descriptor).
    pub faces: Vec<OwnedPtr<topo_ds::Face>>,
    /// Faces whose natural surface normal points into the material (their
    /// shell orientation is Reversed).
    pub concave_faces: Vec<bool>,
    /// OCCT shell(s) assembled in stage 3.5 from the shared-edge faces.
    pub shells: Vec<OwnedPtr<topo_ds::Shell>>,
    /// OCCT solid(s) created in stage 3.6.
    pub solids: Vec<OwnedPtr<topo_ds::Solid>>,
}

impl std::fmt::Debug for Stage3Output {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Stage3Output")
            .field("face_descriptors", &self.face_descriptors.len())
            .field("edges", &self.edges.len())
            .field("vertices", &self.vertices.len())
            .field("faces", &self.faces.len())
            .field("shells", &self.shells.len())
            .field("solids", &self.solids.len())
            .finish()
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Stage3Error {
    AdjacencyError(String),
    Compare(Stage3CompareError),
}

impl Display for Stage3Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage3Error::AdjacencyError(msg) => write!(f, "stage 3.1 adjacency error: {msg}"),
            Stage3Error::Compare(e) => write!(f, "stage 3.{} compare: {e}", e.substage),
        }
    }
}

impl Error for Stage3Error {}

impl From<Stage3CompareError> for Stage3Error {
    fn from(e: Stage3CompareError) -> Self {
        Stage3Error::Compare(e)
    }
}

#[derive(Debug)]
pub struct Stage3CompareError {
    pub substage: u8,
    pub check_type: &'static str,
    pub element_index: usize,
    pub max_distance: f64,
    pub tolerance: f64,
}

impl Display for Stage3CompareError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if self.check_type == "solid volume" {
            return write!(
                f,
                "solid {} volume differs from STEP by {:.2e} relative (tolerance {:.0e})",
                self.element_index, self.max_distance, self.tolerance
            );
        }
        write!(
            f,
            "{} {} is {:.6e} mm from nearest STEP {} (tolerance: {:.6e} mm)",
            self.check_type, self.element_index, self.max_distance,
            if self.check_type == "vertex" { "vertex" } else { "edge" },
            self.tolerance
        )
    }
}

impl Error for Stage3CompareError {}

/// Build a compound containing all edges from the reference STEP shape.
fn build_edge_compound(compare_shape: &topo_ds::Shape) -> OwnedPtr<topo_ds::Compound> {
    let builder = topo_ds::Builder::new();
    let mut compound = topo_ds::Compound::new();
    builder.make_compound(&mut compound);
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        compare_shape,
        top_abs::ShapeEnum::Edge,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        builder.add(compound.as_shape_mut(), explorer.value());
        explorer.next();
    }
    compound
}

/// Build a compound containing all vertices from the reference STEP shape.
fn build_vertex_compound(compare_shape: &topo_ds::Shape) -> OwnedPtr<topo_ds::Compound> {
    let builder = topo_ds::Builder::new();
    let mut compound = topo_ds::Compound::new();
    builder.make_compound(&mut compound);
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        compare_shape,
        top_abs::ShapeEnum::Vertex,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        builder.add(compound.as_shape_mut(), explorer.value());
        explorer.next();
    }
    compound
}

/// Compute minimum distance from a point to a shape.
fn min_distance_to_shape(pt: &gp::Pnt, shape: &topo_ds::Shape) -> f64 {
    stage1::min_distance_to_shape(pt, shape)
}

/// Compare stage 3.1 adjacency graph against reference STEP shape.
///
/// For each ReconEdge, samples mesh boundary vertices along the edge and checks
/// that each is within vertex_tolerance of a STEP edge.
/// For each BRepVertex, checks that it is within vertex_tolerance of a STEP vertex.
fn compare_adjacency_to_step(
    edges: &[ReconEdge],
    vertices: &[BRepVertex],
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), Stage3CompareError> {
    let mesh = &stage2.mesh;
    let compare_shape = config.compare_shape.as_ref().unwrap();

    // Build compounds of edges and vertices from the reference STEP shape
    let edge_compound = build_edge_compound(compare_shape);
    let vertex_compound = build_vertex_compound(compare_shape);

    // Check ReconEdge boundary vertices against STEP edges
    let mut edge_overall_max = 0.0_f64;
    for (ei, edge) in edges.iter().enumerate() {
        let mut max_dist = 0.0_f64;
        let boundary = &edge.mesh_boundary_vertices;

        // Sample mesh boundary vertices: for short edges check all, for longer ones sample
        let step = if boundary.len() <= 10 { 1 } else { boundary.len() / 10 };
        let mut sample_indices: Vec<usize> = (0..boundary.len()).step_by(step).collect();
        // Always include first and last
        if let Some(&last_idx) = sample_indices.last() {
            if last_idx != boundary.len() - 1 {
                sample_indices.push(boundary.len() - 1);
            }
        }

        let mut worst_vi = usize::MAX;
        for &si in &sample_indices {
            let vi = boundary[si];
            let v = &mesh.vertices[vi];
            let pt = gp::Pnt::new_real3(v.x, v.y, v.z);
            let d = min_distance_to_shape(&pt, edge_compound.as_shape());
            if d > max_dist {
                max_dist = d;
                worst_vi = vi;
            }
        }

        edge_overall_max = edge_overall_max.max(max_dist);
        // Boundary vertices of tangent edges are allowed the chordal surface
        // tolerance: an imprecise tessellation can place the region boundary
        // one row off the tangent junction (see stage3() comment).
        let edge_tol = if edge.tangent {
            config.surface_tolerance_mm
        } else {
            config.vertex_tolerance_mm
        };
        if max_dist > edge_tol {
            if !config.quiet && worst_vi != usize::MAX {
                let v = &mesh.vertices[worst_vi];
                eprintln!(
                    "  Compare 3.1: edge {ei} worst boundary vertex {worst_vi} at ({:.4}, {:.4}, {:.4}), dist {max_dist:.6e}",
                    v.x, v.y, v.z,
                );
            }
            return Err(Stage3CompareError {
                substage: 1,
                check_type: "edge",
                element_index: ei,
                max_distance: max_dist,
                tolerance: edge_tol,
            });
        }
    }

    // Check BRepVertices against STEP vertices. Vertices adjacent to any
    // TANGENT edge get the chordal surface tolerance for the same reason as
    // tangent-edge boundary vertices: when the tessellation lacks a precise
    // point at a corner where surfaces meet tangentially (e.g. triple points
    // of a fillet network), the mesh corner vertex can sit a fraction of a
    // facet off the true corner — harmless, since corner positions get
    // refined against the reconstructed surfaces downstream.
    let mut vert_overall_max = 0.0_f64;
    for (vi, vert) in vertices.iter().enumerate() {
        let pt = gp::Pnt::new_real3(vert.point[0], vert.point[1], vert.point[2]);
        let d = min_distance_to_shape(&pt, vertex_compound.as_shape());

        vert_overall_max = vert_overall_max.max(d);
        let near_tangent = vert
            .adjacent_edges
            .iter()
            .any(|&ei| edges.get(ei).is_some_and(|e| e.tangent));
        let vert_tol = if near_tangent {
            config.surface_tolerance_mm
        } else {
            config.vertex_tolerance_mm
        };
        if d > vert_tol {
            return Err(Stage3CompareError {
                substage: 1,
                check_type: "vertex",
                element_index: vi,
                max_distance: d,
                tolerance: vert_tol,
            });
        }
    }

    if !config.quiet {
        eprintln!(
            "  Compare 3.1: edge boundary max dist {:.6e} mm, vertex max dist {:.6e} mm",
            edge_overall_max, vert_overall_max
        );
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Stage 3.1: Create OCCT surfaces and build adjacency graph
// ---------------------------------------------------------------------------

/// Create an OCCT `HandleGeomSurface` for a selected surface hypothesis.
fn create_occt_surface(
    surface: &SelectedSurface,
    output: &Stage2Output,
) -> OwnedPtr<geom::Handle_Geom_Surface> {
    match surface {
        SelectedSurface::Planar(idx) => {
            let hyp = &output.planar_hypotheses[*idx];
            let origin = gp::Pnt::new_real3(
                hyp.normal[0] * hyp.distance,
                hyp.normal[1] * hyp.distance,
                hyp.normal[2] * hyp.distance,
            );
            let normal = gp::Dir::new_real3(hyp.normal[0], hyp.normal[1], hyp.normal[2]);
            let plane = geom::Plane::new_pnt_dir(&origin, &normal);
            geom::Plane::to_handle(plane).to_handle_surface()
        }
        SelectedSurface::Conical(idx) => {
            let hyp = &output.conical_hypotheses[*idx];
            let dir = gp::Dir::new_real3(
                hyp.axis_direction[0],
                hyp.axis_direction[1],
                hyp.axis_direction[2],
            );
            if hyp.is_cylinder() {
                let origin = gp::Pnt::new_real3(hyp.apex[0], hyp.apex[1], hyp.apex[2]);
                let ax3 = gp::Ax3::new_pnt_dir(&origin, &dir);
                let cyl = geom::CylindricalSurface::new_ax3_real(&ax3, hyp.radius_at_apex);
                geom::CylindricalSurface::to_handle(cyl).to_handle_surface()
            } else {
                let origin = gp::Pnt::new_real3(hyp.apex[0], hyp.apex[1], hyp.apex[2]);
                let ax3 = gp::Ax3::new_pnt_dir(&origin, &dir);
                let cone = geom::ConicalSurface::new_ax3_real2(&ax3, hyp.half_angle, 0.0);
                geom::ConicalSurface::to_handle(cone).to_handle_surface()
            }
        }
        SelectedSurface::Spherical(idx) => {
            let hyp = &output.spherical_hypotheses[*idx];
            let origin = gp::Pnt::new_real3(hyp.center[0], hyp.center[1], hyp.center[2]);
            // Sphere axis direction is arbitrary; use Z-up
            let dir = gp::Dir::new_real3(0.0, 0.0, 1.0);
            let ax3 = gp::Ax3::new_pnt_dir(&origin, &dir);
            let sphere = geom::SphericalSurface::new_ax3_real(&ax3, hyp.radius);
            geom::SphericalSurface::to_handle(sphere).to_handle_surface()
        }
        SelectedSurface::Toroidal(idx) => {
            let hyp = &output.toroidal_hypotheses[*idx];
            let origin = gp::Pnt::new_real3(hyp.center[0], hyp.center[1], hyp.center[2]);
            let dir = gp::Dir::new_real3(
                hyp.axis_direction[0],
                hyp.axis_direction[1],
                hyp.axis_direction[2],
            );
            let ax3 = gp::Ax3::new_pnt_dir(&origin, &dir);
            let torus = geom::ToroidalSurface::new_ax3_real2(
                &ax3,
                hyp.major_radius,
                hyp.minor_radius,
            );
            geom::ToroidalSurface::to_handle(torus).to_handle_surface()
        }
        SelectedSurface::BSpline(idx) => {
            // Not a reconstruction: a fresh handle to the exact fitted
            // Geom_BSplineSurface object owned by the stage-2.7 hypothesis.
            output.bspline_hypotheses[*idx].surface.to_handle_surface()
        }
    }
}

/// Get the faces belonging to a selected surface.
fn surface_faces<'a>(surface: &SelectedSurface, output: &'a Stage2Output) -> &'a [usize] {
    match surface {
        SelectedSurface::Planar(idx) => &output.planar_hypotheses[*idx].faces,
        SelectedSurface::Conical(idx) => &output.conical_hypotheses[*idx].faces,
        SelectedSurface::Spherical(idx) => &output.spherical_hypotheses[*idx].faces,
        SelectedSurface::Toroidal(idx) => &output.toroidal_hypotheses[*idx].faces,
        SelectedSurface::BSpline(idx) => &output.bspline_hypotheses[*idx].faces,
    }
}
/// Canonical (undirected) representation of a mesh edge for hashing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct UndirectedEdge(usize, usize);

impl UndirectedEdge {
    fn new(a: usize, b: usize) -> Self {
        if a <= b {
            UndirectedEdge(a, b)
        } else {
            UndirectedEdge(b, a)
        }
    }
}

/// An ordered pair of selected surface indices (smaller first).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct SurfacePair(usize, usize);

impl SurfacePair {
    fn new(a: usize, b: usize) -> Self {
        if a <= b {
            SurfacePair(a, b)
        } else {
            SurfacePair(b, a)
        }
    }
}

/// Build the mapping from mesh face index to selected surface index.
fn build_face_to_surface_map(output: &Stage2Output) -> Vec<usize> {
    let num_faces = output.mesh.faces.len();
    let mut face_to_surface = vec![usize::MAX; num_faces];
    for (si, surface) in output.selected_surfaces.iter().enumerate() {
        for &fi in surface_faces(surface, output) {
            debug_assert_eq!(
                face_to_surface[fi],
                usize::MAX,
                "face {fi} assigned to multiple selected surfaces"
            );
            face_to_surface[fi] = si;
        }
    }
    face_to_surface
}

/// Collect all mesh edges that lie on boundaries between different selected surfaces.
/// Returns a map from surface pair to ordered lists of boundary mesh edges.
fn collect_boundary_edges(
    mesh: &ConnectedMesh,
    face_to_surface: &[usize],
) -> HashMap<SurfacePair, Vec<UndirectedEdge>> {
    let mut boundary_edges: HashMap<SurfacePair, Vec<UndirectedEdge>> = HashMap::new();

    for (fi, face) in mesh.faces.iter().enumerate() {
        let si = face_to_surface[fi];
        let vc = face.vertex_count as usize;
        for ei in 0..vc {
            let ni = face.neighbors[ei];
            if ni < 0 {
                continue;
            }
            let ni = ni as usize;
            let sn = face_to_surface[ni];
            if si == sn {
                continue;
            }
            // This is a boundary edge between surfaces si and sn.
            // Only record each undirected edge once (when fi < ni).
            if fi < ni {
                let v0 = face.vertex_indices[ei];
                let v1 = face.vertex_indices[(ei + 1) % vc];
                let pair = SurfacePair::new(si, sn);
                boundary_edges
                    .entry(pair)
                    .or_default()
                    .push(UndirectedEdge::new(v0, v1));
            }
        }
    }

    boundary_edges
}

/// Find mesh vertices where 3+ selected surfaces meet.
/// Returns a map from mesh vertex index to the set of selected surface indices meeting there.
fn find_corner_vertices(
    mesh: &ConnectedMesh,
    face_to_surface: &[usize],
) -> HashMap<usize, BTreeSet<usize>> {
    let num_verts = mesh.vertices.len();
    let mut vertex_surfaces: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); num_verts];

    for (fi, face) in mesh.faces.iter().enumerate() {
        let si = face_to_surface[fi];
        let vc = face.vertex_count as usize;
        for vi in 0..vc {
            vertex_surfaces[face.vertex_indices[vi]].insert(si);
        }
    }

    let mut corners = HashMap::new();
    for (vi, surfaces) in vertex_surfaces.into_iter().enumerate() {
        if surfaces.len() >= 3 {
            corners.insert(vi, surfaces);
        }
    }
    corners
}

/// Chain boundary edges into contiguous sequences (boundary segments).
/// Each boundary segment is a connected path of edges between the same pair of surfaces.
/// Returns vectors of ordered mesh vertex indices along each chain.
fn chain_boundary_edges(edges: &[UndirectedEdge]) -> Vec<Vec<usize>> {
    if edges.is_empty() {
        return Vec::new();
    }

    // Build adjacency: vertex -> list of connected vertices via boundary edges
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for &UndirectedEdge(a, b) in edges {
        adj.entry(a).or_default().push(b);
        adj.entry(b).or_default().push(a);
    }

    let mut visited_edges: HashSet<UndirectedEdge> = HashSet::new();
    let mut chains = Vec::new();

    // Start chains from endpoints (degree 1 vertices) first, then from any unvisited edge.
    // Endpoints are vertices with exactly 1 boundary edge connection.
    let mut start_vertices: Vec<usize> = adj
        .iter()
        .filter(|(_, neighbors)| neighbors.len() == 1)
        .map(|(&v, _)| v)
        .collect();
    start_vertices.sort_unstable(); // deterministic ordering

    // Also consider all vertices as potential starts for closed loops
    let mut all_vertices: Vec<usize> = adj.keys().copied().collect();
    all_vertices.sort_unstable();

    // Try endpoints first, then any vertex
    for &start in start_vertices.iter().chain(all_vertices.iter()) {
        // Try to extend a chain from this vertex
        let has_unvisited = adj
            .get(&start)
            .is_some_and(|nbrs| {
                nbrs.iter()
                    .any(|&n| !visited_edges.contains(&UndirectedEdge::new(start, n)))
            });
        if !has_unvisited {
            continue;
        }

        let mut chain = vec![start];
        let mut current = start;

        loop {
            let next = adj.get(&current).and_then(|nbrs| {
                nbrs.iter()
                    .copied()
                    .find(|&n| !visited_edges.contains(&UndirectedEdge::new(current, n)))
            });

            match next {
                Some(n) => {
                    visited_edges.insert(UndirectedEdge::new(current, n));
                    if n == start && chain.len() > 1 {
                        // Closed loop — don't repeat the start vertex
                        break;
                    }
                    chain.push(n);
                    current = n;
                }
                None => break,
            }
        }

        if chain.len() >= 2 {
            chains.push(chain);
        }
    }

    chains
}

/// Compute the average position of a mesh vertex projected onto all adjacent surfaces.
/// For now, just use the raw mesh vertex position.
fn vertex_position(vi: usize, vertices: &[MeshVertex]) -> [f64; 3] {
    let v = &vertices[vi];
    [v.x, v.y, v.z]
}

/// Build stage 3.1 output: OCCT surfaces and adjacency graph.
/// Project a point onto a selected surface's fitted geometry (analytic
/// nearest-point for plane/cylinder/cone/sphere/torus).
fn project_point_to_selected(
    p: [f64; 3],
    surface: &SelectedSurface,
    stage2: &Stage2Output,
) -> [f64; 3] {
    match surface {
        SelectedSurface::Planar(i) => {
            let h = &stage2.planar_hypotheses[*i];
            let d = h.normal[0] * p[0] + h.normal[1] * p[1] + h.normal[2] * p[2] - h.distance;
            [
                p[0] - d * h.normal[0],
                p[1] - d * h.normal[1],
                p[2] - d * h.normal[2],
            ]
        }
        SelectedSurface::Conical(i) => {
            let h = &stage2.conical_hypotheses[*i];
            let d = [p[0] - h.apex[0], p[1] - h.apex[1], p[2] - h.apex[2]];
            let a = h.axis_direction;
            let ht = d[0] * a[0] + d[1] * a[1] + d[2] * a[2];
            let radial = [d[0] - ht * a[0], d[1] - ht * a[1], d[2] - ht * a[2]];
            let rd = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            let surface_r = h.radius_at_apex + ht * h.half_angle.tan();
            if rd > 1e-15 {
                let sc = surface_r / rd;
                [
                    h.apex[0] + ht * a[0] + radial[0] * sc,
                    h.apex[1] + ht * a[1] + radial[1] * sc,
                    h.apex[2] + ht * a[2] + radial[2] * sc,
                ]
            } else {
                p
            }
        }
        SelectedSurface::Spherical(i) => {
            let h = &stage2.spherical_hypotheses[*i];
            let d = [p[0] - h.center[0], p[1] - h.center[1], p[2] - h.center[2]];
            let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if dist > 1e-15 {
                let sc = h.radius / dist;
                [
                    h.center[0] + d[0] * sc,
                    h.center[1] + d[1] * sc,
                    h.center[2] + d[2] * sc,
                ]
            } else {
                p
            }
        }
        SelectedSurface::Toroidal(i) => {
            let h = &stage2.toroidal_hypotheses[*i];
            let d = [p[0] - h.center[0], p[1] - h.center[1], p[2] - h.center[2]];
            let a = h.axis_direction;
            let ax = d[0] * a[0] + d[1] * a[1] + d[2] * a[2];
            let radial = [d[0] - ax * a[0], d[1] - ax * a[1], d[2] - ax * a[2]];
            let rd = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            if rd > 1e-15 {
                let tc = [
                    h.center[0] + h.major_radius * radial[0] / rd,
                    h.center[1] + h.major_radius * radial[1] / rd,
                    h.center[2] + h.major_radius * radial[2] / rd,
                ];
                let tv = [p[0] - tc[0], p[1] - tc[1], p[2] - tc[2]];
                let td = (tv[0] * tv[0] + tv[1] * tv[1] + tv[2] * tv[2]).sqrt();
                if td > 1e-15 {
                    let sc = h.minor_radius / td;
                    [tc[0] + tv[0] * sc, tc[1] + tv[1] * sc, tc[2] + tv[2] * sc]
                } else {
                    p
                }
            } else {
                p
            }
        }
        SelectedSurface::BSpline(i) => stage2.bspline_hypotheses[*i].project_point(&p),
    }
}

/// Refine corner vertex positions by cyclic projection onto the adjacent
/// fitted surfaces (DEVELOPMENT_PLAN BRepVertex TODO). The raw mesh vertex
/// at a corner is only as accurate as the tessellation puts it — where a
/// fillet terminates, the mesh may have no vertex exactly at the CAD corner.
/// Alternating projections converge to a point on the intersection of the
/// adjacent surfaces. The refinement is kept only if it stays within
/// 2× surface tolerance of the mesh vertex (safety against divergence on
/// near-parallel surface configurations).
fn refine_corner_vertices(
    vertices: &mut [BRepVertex],
    edges: &[ReconEdge],
    face_descriptors: &[FaceDescriptor],
    stage2: &Stage2Output,
    config: &Config,
) {
    // Two adjacent corners must never refine onto the same point: the edge
    // between them collapses to zero length and MakeEdge rejects it
    // (DifferentPointsOnClosedCurve; even its straight-line fallback needs
    // distinct points). Adjacent surface sets that pin the same intersection
    // point arise when fits are coarser than the mesh — observed on 3DBenchy
    // with --skip-stage forcing planar-only fits, where near-coplanar facet
    // planes gave both endpoints of edge 5715 the same projection fixpoint
    // (4e-14 mm apart). Below this separation from any edge-adjacent corner,
    // the candidate is discarded and the corner keeps the raw mesh position
    // (mesh vertices themselves are always ≥ the weld tolerance apart).
    let min_sep = config.vertex_tolerance_mm.max(1e-6);

    for vi in 0..vertices.len() {
        let vert = &vertices[vi];
        if vert.adjacent_faces.len() < 2 {
            continue;
        }
        // Only refine SHARP corners: where any adjacent surface pair is
        // tangent, the surfaces intersect along a curve rather than an
        // isolated point, and cyclic projection slides along the tangent
        // direction to an arbitrary point on it (observed drifting 0.57 mm).
        // Tangent-adjacent corners keep the raw mesh vertex (and the relaxed
        // chordal compare tolerance).
        let near_tangent = vert
            .adjacent_edges
            .iter()
            .any(|&ei| edges.get(ei).is_some_and(|e| e.tangent));
        if near_tangent {
            continue;
        }
        // B-spline-adjacent corners also keep the raw mesh vertex: the
        // fitted surface only approximates the mesh (to fit tolerance,
        // ~1e-2 mm), so cyclic projection onto it can only degrade a mesh
        // corner that sits exactly on the true feature curve.
        let near_bspline = vert.adjacent_faces.iter().any(|&fi| {
            matches!(
                stage2.selected_surfaces[face_descriptors[fi].selected_surface_idx],
                SelectedSurface::BSpline(_)
            )
        });
        if near_bspline {
            continue;
        }
        let original = vert.point;
        let mut p = original;
        for _iter in 0..30 {
            let prev = p;
            for &fi in &vert.adjacent_faces {
                let surface =
                    &stage2.selected_surfaces[face_descriptors[fi].selected_surface_idx];
                p = project_point_to_selected(p, surface, stage2);
            }
            let dm = (p[0] - prev[0]).powi(2) + (p[1] - prev[1]).powi(2) + (p[2] - prev[2]).powi(2);
            if dm < 1e-20 {
                break;
            }
        }
        let moved = ((p[0] - original[0]).powi(2)
            + (p[1] - original[1]).powi(2)
            + (p[2] - original[2]).powi(2))
        .sqrt();
        if moved > config.surface_tolerance_mm * 2.0 {
            continue;
        }
        let collapses_adjacent_edge = vertices[vi].adjacent_edges.iter().any(|&ei| {
            edges.get(ei).is_some_and(|e| e.vertex_indices.iter().any(|&ui| {
                ui != usize::MAX && ui != vi && {
                    let q = vertices[ui].point;
                    (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)
                        < min_sep * min_sep
                }
            }))
        });
        if collapses_adjacent_edge {
            if config.verbosity >= 1 {
                eprintln!(
                    "  Corner refinement: vertex {vi} kept at mesh position — refined \
                     point within {min_sep:.1e} mm of an adjacent corner"
                );
            }
            continue;
        }
        vertices[vi].point = p;
    }
}

fn build_surfaces_and_adjacency(
    input: Stage2Output,
    config: &Config,
) -> Result<Stage3Output, Stage3Error> {
    let t = Instant::now();
    let num_surfaces = input.selected_surfaces.len();

    // 1. Build face-to-surface mapping
    let face_to_surface = build_face_to_surface_map(&input);

    // Verify all faces are assigned
    for (fi, &si) in face_to_surface.iter().enumerate() {
        if si == usize::MAX {
            return Err(Stage3Error::AdjacencyError(format!(
                "mesh face {fi} not assigned to any selected surface"
            )));
        }
    }

    // 2. Create OCCT surface objects
    let mut face_descriptors: Vec<FaceDescriptor> = Vec::with_capacity(num_surfaces);
    for (si, surface) in input.selected_surfaces.iter().enumerate() {
        let occt_surface = create_occt_surface(surface, &input);
        face_descriptors.push(FaceDescriptor {
            selected_surface_idx: si,
            surface: occt_surface,
            adjacent_faces: Vec::new(),
            edge_indices: Vec::new(),
            vertex_indices: Vec::new(),
        });
    }

    // 3. Collect boundary edges between surfaces and find corner vertices
    let boundary_edges = collect_boundary_edges(&input.mesh, &face_to_surface);
    let corner_vertices = find_corner_vertices(&input.mesh, &face_to_surface);

    // 4. Create BRepVertices for corners
    // Map from mesh vertex index to BRepVertex index
    let mut mesh_vert_to_brep_vert: HashMap<usize, usize> = HashMap::new();
    let mut brep_vertices: Vec<BRepVertex> = Vec::new();

    // Sort corner vertices for deterministic output
    let mut corner_list: Vec<(usize, BTreeSet<usize>)> = corner_vertices.into_iter().collect();
    corner_list.sort_by_key(|&(vi, _)| vi);

    for (vi, surfaces) in &corner_list {
        let brep_idx = brep_vertices.len();
        mesh_vert_to_brep_vert.insert(*vi, brep_idx);
        brep_vertices.push(BRepVertex {
            adjacent_faces: surfaces.iter().copied().collect(),
            adjacent_edges: Vec::new(), // filled below
            point: vertex_position(*vi, &input.mesh.vertices),
        });
    }

    // 5. Chain boundary edges into ReconEdges
    let mut recon_edges: Vec<ReconEdge> = Vec::new();

    // Process surface pairs in deterministic order
    let mut sorted_pairs: Vec<(SurfacePair, Vec<UndirectedEdge>)> =
        boundary_edges.into_iter().collect();
    sorted_pairs.sort_by_key(|(pair, _)| *pair);

    // Map from surface pair to list of ReconEdge indices for that pair
    let mut pair_to_edges: HashMap<SurfacePair, Vec<usize>> = HashMap::new();

    for (pair, edges) in sorted_pairs {
        let chains = chain_boundary_edges(&edges);

        for chain in chains {
            let edge_idx = recon_edges.len();

            // Determine endpoint BRepVertices
            let start_vi = chain.first().copied().unwrap();
            let end_vi = chain.last().copied().unwrap();
            let start_brep = mesh_vert_to_brep_vert
                .get(&start_vi)
                .copied()
                .unwrap_or(usize::MAX);
            let end_brep = mesh_vert_to_brep_vert
                .get(&end_vi)
                .copied()
                .unwrap_or(usize::MAX);

            // Register this edge with its BRepVertices
            if start_brep != usize::MAX {
                brep_vertices[start_brep].adjacent_edges.push(edge_idx);
            }
            if end_brep != usize::MAX && end_brep != start_brep {
                brep_vertices[end_brep].adjacent_edges.push(edge_idx);
            }

            recon_edges.push(ReconEdge {
                face_indices: [pair.0, pair.1],
                vertex_indices: [start_brep, end_brep],
                curve_3d: None, // populated in stage 3.3
                tangent: false, // determined in stage 3.2
                mesh_boundary_vertices: chain,
            });

            pair_to_edges.entry(pair).or_default().push(edge_idx);
        }
    }

    // 6. Build adjacency lists for each FaceDescriptor
    // For each face descriptor, collect all edges that touch it and the adjacent faces.
    let mut face_edges: Vec<Vec<(usize, usize)>> = vec![Vec::new(); num_surfaces];
    // face_edges[si] = list of (edge_idx, other_surface_idx)
    for (ei, edge) in recon_edges.iter().enumerate() {
        let [f0, f1] = edge.face_indices;
        face_edges[f0].push((ei, f1));
        face_edges[f1].push((ei, f0));
    }

    // For each face, order the edges topologically around the boundary using BRepVertex connections.
    // For simple cases (non-self-adjacent faces), the ordering follows the mesh boundary.
    for si in 0..num_surfaces {
        let edges_and_adj = &face_edges[si];
        if edges_and_adj.is_empty() {
            continue;
        }

        // Build a traversal order using edge-vertex-edge chains.
        // Each edge has two endpoint BRepVertices. We walk: edge -> end vertex -> next edge
        // sharing that vertex (and this face).
        let num_edges = edges_and_adj.len();

        if num_edges == 1 {
            // Single adjacent face — simple case
            let (ei, adj_fi) = edges_and_adj[0];
            face_descriptors[si].adjacent_faces.push(adj_fi);
            face_descriptors[si].edge_indices.push(ei);
            // Vertex indices: the endpoints of this edge
            let edge = &recon_edges[ei];
            if edge.vertex_indices[0] != usize::MAX {
                face_descriptors[si].vertex_indices.push(edge.vertex_indices[0]);
            }
            if edge.vertex_indices[1] != usize::MAX
                && edge.vertex_indices[1] != edge.vertex_indices[0]
            {
                face_descriptors[si].vertex_indices.push(edge.vertex_indices[1]);
            }
            continue;
        }

        // Build a map from BRepVertex to edges of this face that touch it
        let mut vert_to_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        for &(ei, _) in edges_and_adj {
            let edge = &recon_edges[ei];
            for &vi in &edge.vertex_indices {
                if vi != usize::MAX {
                    vert_to_edges.entry(vi).or_default().push(ei);
                }
            }
        }

        // Walk the edge chain starting from the first edge
        let mut ordered_edges: Vec<usize> = Vec::with_capacity(num_edges);
        let mut ordered_adj: Vec<usize> = Vec::with_capacity(num_edges);
        let mut ordered_verts: Vec<usize> = Vec::new();
        let mut visited_edges: HashSet<usize> = HashSet::new();

        let (first_ei, first_adj) = edges_and_adj[0];
        ordered_edges.push(first_ei);
        ordered_adj.push(first_adj);
        visited_edges.insert(first_ei);

        // Walk forward from the end of the first edge
        let mut current_ei = first_ei;
        for _ in 0..num_edges {
            let edge = &recon_edges[current_ei];
            // Try both endpoints to find the next unvisited edge
            let mut found_next = false;
            for &vi in &edge.vertex_indices {
                if vi == usize::MAX {
                    continue;
                }
                if let Some(connected) = vert_to_edges.get(&vi) {
                    for &next_ei in connected {
                        if !visited_edges.contains(&next_ei) {
                            ordered_verts.push(vi);
                            visited_edges.insert(next_ei);
                            // Find the adjacent face for next_ei as seen from si
                            let next_adj = if recon_edges[next_ei].face_indices[0] == si {
                                recon_edges[next_ei].face_indices[1]
                            } else {
                                recon_edges[next_ei].face_indices[0]
                            };
                            ordered_edges.push(next_ei);
                            ordered_adj.push(next_adj);
                            current_ei = next_ei;
                            found_next = true;
                            break;
                        }
                    }
                }
                if found_next {
                    break;
                }
            }
            if !found_next {
                break;
            }
        }

        // If there are unvisited edges (disconnected boundaries), just append them
        for &(ei, adj_fi) in edges_and_adj {
            if !visited_edges.contains(&ei) {
                ordered_edges.push(ei);
                ordered_adj.push(adj_fi);
            }
        }

        // Close the loop: add the vertex connecting the last edge back to the first
        if ordered_edges.len() >= 2 {
            let last_edge = &recon_edges[*ordered_edges.last().unwrap()];
            let first_edge = &recon_edges[ordered_edges[0]];
            // Find a shared vertex between last and first edges
            for &vi in &last_edge.vertex_indices {
                if vi != usize::MAX {
                    for &fvi in &first_edge.vertex_indices {
                        if vi == fvi {
                            ordered_verts.push(vi);
                        }
                    }
                }
            }
        }

        face_descriptors[si].adjacent_faces = ordered_adj;
        face_descriptors[si].edge_indices = ordered_edges;
        face_descriptors[si].vertex_indices = ordered_verts;
    }

    // Print summary
    if !config.quiet {
        eprintln!(
            "Stage 3.1 ({:.3}s): Created {} OCCT surfaces, {} edges, {} vertices",
            t.elapsed().as_secs_f64(),
            face_descriptors.len(),
            recon_edges.len(),
            brep_vertices.len(),
        );
        if config.verbose {
            for (si, fd) in face_descriptors.iter().enumerate() {
                let stype = match &input.selected_surfaces[si] {
                    SelectedSurface::Planar(i) => {
                        let h = &input.planar_hypotheses[*i];
                        format!(
                            "planar n=[{:.4},{:.4},{:.4}] d={:.4}",
                            h.normal[0], h.normal[1], h.normal[2], h.distance,
                        )
                    }
                    SelectedSurface::Conical(i) => {
                        let h = &input.conical_hypotheses[*i];
                        if h.is_cylinder() {
                            format!(
                                "cylindrical r={:.4} apex=[{:.4},{:.4},{:.4}] axis=[{:.4},{:.4},{:.4}]",
                                h.radius_at_apex, h.apex[0], h.apex[1], h.apex[2],
                                h.axis_direction[0], h.axis_direction[1], h.axis_direction[2],
                            )
                        } else {
                            format!(
                                "conical ha={:.2}deg apex=[{:.4},{:.4},{:.4}] axis=[{:.4},{:.4},{:.4}]",
                                h.half_angle.to_degrees(), h.apex[0], h.apex[1], h.apex[2],
                                h.axis_direction[0], h.axis_direction[1], h.axis_direction[2],
                            )
                        }
                    }
                    SelectedSurface::Spherical(i) => {
                        let h = &input.spherical_hypotheses[*i];
                        format!(
                            "spherical r={:.4} center=[{:.4},{:.4},{:.4}]",
                            h.radius, h.center[0], h.center[1], h.center[2],
                        )
                    }
                    SelectedSurface::Toroidal(i) => {
                        let h = &input.toroidal_hypotheses[*i];
                        format!(
                            "toroidal R={:.4} r={:.4} center=[{:.4},{:.4},{:.4}] axis=[{:.4},{:.4},{:.4}]",
                            h.major_radius, h.minor_radius,
                            h.center[0], h.center[1], h.center[2],
                            h.axis_direction[0], h.axis_direction[1], h.axis_direction[2],
                        )
                    }
                    SelectedSurface::BSpline(i) => {
                        let h = &input.bspline_hypotheses[*i];
                        format!(
                            "bspline deg=({},{}) poles=({},{}) vtx_err_max={:.2e}",
                            h.u_degree, h.v_degree, h.nb_u_poles, h.nb_v_poles,
                            h.error_max,
                        )
                    }
                };
                eprintln!(
                    "  Face {}: {} surface, {} adjacent faces, {} edges, {} vertices",
                    si,
                    stype,
                    fd.adjacent_faces.len(),
                    fd.edge_indices.len(),
                    fd.vertex_indices.len(),
                );
            }
            for (ei, edge) in recon_edges.iter().enumerate() {
                eprintln!(
                    "  Edge {}: faces [{}, {}], vertices [{}, {}], {} boundary mesh verts",
                    ei,
                    edge.face_indices[0],
                    edge.face_indices[1],
                    if edge.vertex_indices[0] == usize::MAX {
                        "none".to_string()
                    } else {
                        edge.vertex_indices[0].to_string()
                    },
                    if edge.vertex_indices[1] == usize::MAX {
                        "none".to_string()
                    } else {
                        edge.vertex_indices[1].to_string()
                    },
                    edge.mesh_boundary_vertices.len(),
                );
            }
            for (vi, vert) in brep_vertices.iter().enumerate() {
                eprintln!(
                    "  Vertex {}: point=[{:.4}, {:.4}, {:.4}], {} adjacent faces, {} adjacent edges",
                    vi,
                    vert.point[0],
                    vert.point[1],
                    vert.point[2],
                    vert.adjacent_faces.len(),
                    vert.adjacent_edges.len(),
                );
            }
        }
    }
    // Note: corner vertex refinement and the stage 3.1 --compare check run
    // from stage3() after tangency detection — refinement applies only to
    // sharp corners, and tangent edges use a looser boundary-vertex
    // tolerance.

    Ok(Stage3Output {
        stage2: input,
        face_descriptors,
        edges: recon_edges,
        vertices: brep_vertices,
        faces: Vec::new(), // populated in stage 3.4
        concave_faces: Vec::new(), // populated in stage 3.4
        shells: Vec::new(), // populated in stage 3.5
        solids: Vec::new(), // populated in stage 3.6
    })}

// ---------------------------------------------------------------------------
// Stage 3.2: Detect tangency relationships
// ---------------------------------------------------------------------------

/// Compute the outward-facing surface normal at a 3D point for a given selected surface.
///
/// For planes: the hypothesis normal (constant everywhere).
/// For cylinders: the radial direction from axis to point, outward if convex.
/// For spheres: the direction from center to point, outward if convex.
///
/// Returns None if the point is degenerate (e.g., on the cylinder axis or sphere center).
fn surface_normal_at_point(
    surface: &SelectedSurface,
    stage2: &Stage2Output,
    point: &[f64; 3],
) -> Option<[f64; 3]> {
    match surface {
        SelectedSurface::Planar(idx) => {
            let hyp = &stage2.planar_hypotheses[*idx];
            Some(hyp.normal)
        }
        SelectedSurface::Conical(idx) => {
            let hyp = &stage2.conical_hypotheses[*idx];
            let dp = [
                point[0] - hyp.apex[0],
                point[1] - hyp.apex[1],
                point[2] - hyp.apex[2],
            ];
            let h = dp[0] * hyp.axis_direction[0]
                + dp[1] * hyp.axis_direction[1]
                + dp[2] * hyp.axis_direction[2];
            let radial = [
                dp[0] - h * hyp.axis_direction[0],
                dp[1] - h * hyp.axis_direction[1],
                dp[2] - h * hyp.axis_direction[2],
            ];
            let radial_len = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            if radial_len < 1e-15 {
                return None; // Point is on the axis
            }
            // Outward normal: radial_unit * cos(half_angle) - axis * sin(half_angle)
            // For cylinders (half_angle=0): simplifies to just radial_unit
            let cos_ha = hyp.half_angle.cos();
            let sin_ha = hyp.half_angle.sin();
            let ru = [radial[0] / radial_len, radial[1] / radial_len, radial[2] / radial_len];
            let mut n = [
                ru[0] * cos_ha - hyp.axis_direction[0] * sin_ha,
                ru[1] * cos_ha - hyp.axis_direction[1] * sin_ha,
                ru[2] * cos_ha - hyp.axis_direction[2] * sin_ha,
            ];
            if !hyp.convex {
                n = [-n[0], -n[1], -n[2]];
            }
            Some(n)
        }
        SelectedSurface::Spherical(idx) => {
            let hyp = &stage2.spherical_hypotheses[*idx];
            let dp = [
                point[0] - hyp.center[0],
                point[1] - hyp.center[1],
                point[2] - hyp.center[2],
            ];
            let len = (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]).sqrt();
            if len < 1e-15 {
                return None; // Point is at the center
            }
            let mut n = [dp[0] / len, dp[1] / len, dp[2] / len];
            if !hyp.convex {
                n = [-n[0], -n[1], -n[2]];
            }
            Some(n)
        }
        SelectedSurface::Toroidal(idx) => {
            let hyp = &stage2.toroidal_hypotheses[*idx];
            let dp = [
                point[0] - hyp.center[0],
                point[1] - hyp.center[1],
                point[2] - hyp.center[2],
            ];
            // Project onto the major circle plane
            let ax = hyp.axis_direction;
            let h = dp[0] * ax[0] + dp[1] * ax[1] + dp[2] * ax[2];
            let radial = [
                dp[0] - h * ax[0],
                dp[1] - h * ax[1],
                dp[2] - h * ax[2],
            ];
            let radial_len = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            if radial_len < 1e-15 {
                return None; // Point is on the axis
            }
            // Point on the major circle closest to our point
            let tube_center = [
                hyp.center[0] + hyp.major_radius * radial[0] / radial_len,
                hyp.center[1] + hyp.major_radius * radial[1] / radial_len,
                hyp.center[2] + hyp.major_radius * radial[2] / radial_len,
            ];
            let tube_vec = [
                point[0] - tube_center[0],
                point[1] - tube_center[1],
                point[2] - tube_center[2],
            ];
            let tube_len = (tube_vec[0] * tube_vec[0] + tube_vec[1] * tube_vec[1] + tube_vec[2] * tube_vec[2]).sqrt();
            if tube_len < 1e-15 {
                return None;
            }
            let mut n = [tube_vec[0] / tube_len, tube_vec[1] / tube_len, tube_vec[2] / tube_len];
            if !hyp.convex {
                n = [-n[0], -n[1], -n[2]];
            }
            Some(n)
        }
        SelectedSurface::BSpline(idx) => {
            // The fit orients the surface's natural normal outward
            // (concave is always false), so d1u × d1v needs no sign logic.
            stage2.bspline_hypotheses[*idx].normal_at(point)
        }
    }
}
///
/// For each ReconEdge, samples mesh boundary vertices and compares the surface
/// normals of the two adjacent surfaces. If all sampled normals agree within
/// a small angle threshold (2°), the edge is marked as tangent.
fn detect_tangency(mut output: Stage3Output, config: &Config) -> Stage3Output {
    let t = Instant::now();
    const TANGENCY_ANGLE_DEG: f64 = 2.0;
    let tangency_cos = (TANGENCY_ANGLE_DEG * std::f64::consts::PI / 180.0).cos();
    let min_samples = 3;
    let sample_step = 5; // sample every 5th boundary vertex

    let mut tangent_count = 0;

    for edge in output.edges.iter_mut() {
        let [fi0, fi1] = edge.face_indices;
        let surf0 = &output.stage2.selected_surfaces[output.face_descriptors[fi0].selected_surface_idx];
        let surf1 = &output.stage2.selected_surfaces[output.face_descriptors[fi1].selected_surface_idx];

        let boundary = &edge.mesh_boundary_vertices;
        if boundary.len() < 2 {
            continue;
        }

        // Choose sample indices: every sample_step-th vertex, at least min_samples,
        // always including first and last.
        let step = if boundary.len() <= min_samples * sample_step {
            (boundary.len() / min_samples).max(1)
        } else {
            sample_step
        };
        let mut sample_indices: Vec<usize> = (0..boundary.len()).step_by(step).collect();
        if let Some(&last) = sample_indices.last() {
            if last != boundary.len() - 1 {
                sample_indices.push(boundary.len() - 1);
            }
        }

        let mut all_tangent = true;
        for &si in &sample_indices {
            let vi = boundary[si];
            let v = &output.stage2.mesh.vertices[vi];
            let pt = [v.x, v.y, v.z];

            let n0 = surface_normal_at_point(surf0, &output.stage2, &pt);
            let n1 = surface_normal_at_point(surf1, &output.stage2, &pt);

            match (n0, n1) {
                (Some(a), Some(b)) => {
                    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                    if dot < tangency_cos {
                        all_tangent = false;
                        break;
                    }
                }
                _ => {
                    // Degenerate point — cannot determine tangency
                    all_tangent = false;
                    break;
                }
            }
        }

        if all_tangent {
            edge.tangent = true;
            tangent_count += 1;
        }
    }

    if !config.quiet {
        eprintln!(
            "Stage 3.2 ({:.3}s): Tangency detection: {} of {} edges are tangent",
            t.elapsed().as_secs_f64(),
            tangent_count,
            output.edges.len(),
        );
    }

    output
}

// ---------------------------------------------------------------------------
// Stage 3.3: Compute edge curves via surface-surface intersection
// ---------------------------------------------------------------------------

/// Compute the 3D intersection curve for a single ReconEdge.
///
/// Uses GeomAPI_IntSS to find intersection curves between the two adjacent surfaces.
/// If multiple curves are returned, selects the one closest to the mesh boundary vertices.
/// Trims the selected curve to the vertex endpoint parameters.
fn compute_edge_curve(
    edge: &mut ReconEdge,
    face_descriptors: &[FaceDescriptor],
    vertices: &[BRepVertex],
    mesh: &ConnectedMesh,
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    let [fi0, fi1] = edge.face_indices;
    let surf0 = &face_descriptors[fi0].surface;
    let surf1 = &face_descriptors[fi1].surface;

    // Compute surface-surface intersection. IntSS can throw rather than
    // report failure (observed: "BSpline curve: Knots interval values too
    // close" approximating a bspline-bspline junction on 3DBenchy) — treat
    // a throw like any other IntSS failure so the caller's mesh-boundary-
    // chain fallback supplies the edge curve.
    let int_ss = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        geom_api::IntSS::new_handlegeomsurface2_real(surf0, surf1, config.vertex_tolerance_mm)
    }))
    .map_err(|e| {
        let msg = e
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| e.downcast_ref::<&str>().copied())
            .unwrap_or("unknown panic");
        format!("IntSS threw for edge between faces {fi0} and {fi1}: {msg}")
    })?;

    if !int_ss.is_done() || int_ss.nb_lines() == 0 {
        return Err(format!(
            "IntSS failed for edge between faces {} and {} (is_done={}, nb_lines={})",
            fi0, fi1, int_ss.is_done(),
            if int_ss.is_done() { int_ss.nb_lines() } else { -1 }
        ));
    }

    let nb_lines = int_ss.nb_lines();

    // Select the best intersection curve (closest to mesh boundary vertices)
    let best_line_idx = if nb_lines == 1 {
        1 // 1-indexed; parameter range doesn't matter — we trim later
    } else {
        select_closest_curve(&int_ss, &edge.mesh_boundary_vertices, mesh)
    };

    if best_line_idx == 0 {
        return Err(format!(
            "no valid intersection curve for edge between faces {} and {} \
             (all {} curves have extreme parameters)",
            fi0, fi1, nb_lines,
        ));
    }

    let curve_handle = int_ss.line(best_line_idx);

    // Validate the selected curve: check if it passes near the boundary vertices.
    // IntSS can return degenerate curves for plane-through-cylinder-axis intersections.
    let curve_valid = {
        let mut valid = true;
        for &vi in edge.mesh_boundary_vertices.iter().take(3) {
            let v = &mesh.vertices[vi];
            let pt = gp::Pnt::new_real3(v.x, v.y, v.z);
            let proj = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(&pt, curve_handle);
            if proj.nb_points() == 0 || proj.lower_distance() > 1.0 {
                valid = false;
                break;
            }
        }
        valid
    };

    // If the IntSS curve is degenerate (doesn't pass near boundary vertices),
    // construct a line from vertex positions. IntSS can return degenerate curves
    // when a plane passes through a cylinder's axis.
    let fallback_curve: Option<OwnedPtr<geom::Handle_Geom_Curve>> = if !curve_valid
        && edge.vertex_indices[0] != usize::MAX
        && edge.vertex_indices[1] != usize::MAX
    {
        let v0 = &vertices[edge.vertex_indices[0]];
        let v1 = &vertices[edge.vertex_indices[1]];
        let dx = v1.point[0] - v0.point[0];
        let dy = v1.point[1] - v0.point[1];
        let dz = v1.point[2] - v0.point[2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if len < 1e-15 {
            return Err(format!(
                "degenerate edge between faces {} and {}: zero-length vertex span",
                fi0, fi1
            ));
        }
        let p0 = gp::Pnt::new_real3(v0.point[0], v0.point[1], v0.point[2]);
        let dir = gp::Dir::new_real3(dx / len, dy / len, dz / len);
        let line = geom::Line::new_pnt_dir(&p0, &dir);
        let line_handle = geom::Line::to_handle(line).to_handle_curve();

        if config.verbose {
            eprintln!("    Edge (f{fi0},f{fi1}): IntSS curve degenerate, using vertex-based line");
        }
        Some(line_handle)
    } else {
        None
    };

    let initial_curve_handle = match fallback_curve.as_ref() {
        Some(h) => &**h,
        None => int_ss.line(best_line_idx),
    };

    // For closed-loop edges, IntSS may return partial arcs (e.g., semicircles)
    // when the intersection circle passes through the sphere's UV poles (the
    // intersection plane contains the sphere axis). In UV space, such a circle
    // splits into two disconnected arcs, which IntSS returns separately.
    // Detect this and reconstruct the full circle from sampled curve points.
    // Partiality is decided by the endpoint gap, NOT the parameter span: a
    // closed junction with a B-spline surface comes back as a closed curve
    // over an arbitrary range (observed [0,1]), and treating it as a
    // partial arc would replace the true junction with a bogus circle.
    let full_circle_curve: Option<OwnedPtr<geom::Handle_Geom_Curve>> =
        if edge.vertex_indices[0] == usize::MAX && edge.vertex_indices[1] == usize::MAX {
            let c = initial_curve_handle.get();
            let span = c.last_parameter() - c.first_parameter();
            let p_first = c.value(c.first_parameter());
            let p_last = c.value(c.last_parameter());
            let endpoint_gap = ((p_first.x() - p_last.x()).powi(2)
                + (p_first.y() - p_last.y()).powi(2)
                + (p_first.z() - p_last.z()).powi(2))
            .sqrt();
            if endpoint_gap > 10.0 * config.vertex_tolerance_mm {
                // Partial arc — reconstruct full circle from 3 sampled points
                let fp = c.first_parameter();
                let p0 = c.value(fp);
                let p1 = c.value(fp + span * 0.25);
                let p2 = c.value(fp + span * 0.75);

                match reconstruct_full_circle(&p0, &p1, &p2) {
                    Some(circle_handle) => {
                        if config.verbose {
                            let cc = circle_handle.get();
                            let r = cc.value(0.0);
                            eprintln!(
                                "    Edge (f{fi0},f{fi1}): IntSS returned partial arc \
                                 (span={span:.4}), reconstructed full circle"
                            );
                            let _ = r; // suppress unused
                        }
                        Some(circle_handle)
                    }
                    None => None,
                }
            } else {
                None
            }
        } else {
            None
        };

    let curve_handle: &geom::Handle_Geom_Curve = match full_circle_curve.as_ref() {
        Some(h) => h,
        None => initial_curve_handle,
    };

    // Trim the curve to vertex endpoint parameters
    let curve = curve_handle.get();
    let first_param = curve.first_parameter();
    let last_param = curve.last_parameter();

    let (t_start, t_end) = if edge.vertex_indices[0] == usize::MAX
        && edge.vertex_indices[1] == usize::MAX
    {
        // Closed loop — use full curve
        (first_param, last_param)
    } else {
        // Trim to vertex endpoints
        let t_start = if edge.vertex_indices[0] != usize::MAX {
            let v = &vertices[edge.vertex_indices[0]];
            project_point_on_curve(&v.point, curve_handle)?
        } else {
            // No start vertex — project first mesh boundary vertex
            let vi = edge.mesh_boundary_vertices[0];
            let v = &mesh.vertices[vi];
            project_point_on_curve(&[v.x, v.y, v.z], curve_handle)?
        };

        let t_end = if edge.vertex_indices[1] != usize::MAX {
            let v = &vertices[edge.vertex_indices[1]];
            project_point_on_curve(&v.point, curve_handle)?
        } else {
            // No end vertex — project last mesh boundary vertex
            let vi = *edge.mesh_boundary_vertices.last().unwrap();
            let v = &mesh.vertices[vi];
            project_point_on_curve(&[v.x, v.y, v.z], curve_handle)?
        };

        // For closed curves (circles from cylinder/sphere intersections), check
        // which arc the mesh boundary vertices lie on. IntSS circle curves report
        // is_periodic()=false, so detect them by parameter span ≈ 2π.
        let (t_lo, t_hi) = if t_start < t_end {
            (t_start, t_end)
        } else {
            (t_end, t_start)
        };

        let param_span = last_param - first_param;
        let is_closed_curve = (param_span - 2.0 * std::f64::consts::PI).abs() < 1e-6;
        if is_closed_curve && edge.mesh_boundary_vertices.len() > 2 {
            let period = param_span;
            // Sample mesh boundary vertices to determine which arc they lie on
            let step = (edge.mesh_boundary_vertices.len() / 5).max(1);
            let mut in_direct = 0;
            let mut in_complement = 0;
            for i in (0..edge.mesh_boundary_vertices.len()).step_by(step) {
                let vi = edge.mesh_boundary_vertices[i];
                let v = &mesh.vertices[vi];
                if let Ok(t_m) = project_point_on_curve(&[v.x, v.y, v.z], curve_handle) {
                    if t_m >= t_lo && t_m <= t_hi {
                        in_direct += 1;
                    } else {
                        in_complement += 1;
                    }
                }
            }
            if in_complement > in_direct {
                // Boundary vertices lie on the complementary arc
                (t_hi, t_lo + period)
            } else {
                (t_lo, t_hi)
            }
        } else if is_closed_curve {
            // Periodic curve with ≤2 boundary vertices (only endpoints, no intermediate
            // samples to determine arc direction). Prefer the shorter arc — polygon edges
            // on CAD surfaces always take the short path between vertices.
            let direct_span = t_hi - t_lo;
            if direct_span > param_span / 2.0 {
                (t_hi, t_lo + param_span)
            } else {
                (t_lo, t_hi)
            }
        } else if t_start < t_end {
            (t_start, t_end)
        } else {
            (t_end, t_start)
        }
    };
    // Create trimmed curve. Geom_TrimmedCurve refuses a degenerate (zero-length)
    // span (Standard_ConstructionError: U1 == U2) — observed when both edge
    // endpoints project to the same curve parameter, e.g. two vertices closer
    // together than the curve's parametrization resolves (a near-duplicate
    // vertex pair from a tiny fallback planar facet). Same guard as the
    // chain-repair path's `candidate_trim_parameters`.
    if (t_end - t_start).abs() < 1e-9 {
        return Err(format!(
            "degenerate trim span for edge between faces {fi0} and {fi1}: \
             t_start={t_start:.6e}, t_end={t_end:.6e}"
        ));
    }
    let trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(
        curve_handle, t_start, t_end,
    );
    let trimmed_handle = geom::TrimmedCurve::to_handle(trimmed).to_handle_curve();

    // Check if trimmed curve endpoints are close to vertex positions.
    // If not (gap > 0.3mm), the IntSS trimming was inaccurate (e.g., the vertex
    // is at a junction of 3+ surfaces and the IntSS parametric trim missed it).
    //
    // Instead of falling back to a straight line (which wouldn't lie on curved
    // surfaces), re-project the vertices onto the original untrimmed IntSS curve
    // and construct the correct arc. For sphere-plane intersections, the IntSS
    // curve is a circle, so the arc lies on both surfaces.
    let final_curve = if edge.vertex_indices[0] != usize::MAX
        && edge.vertex_indices[1] != usize::MAX
    {
        let tc = trimmed_handle.get();
        let p_start = tc.value(tc.first_parameter());
        let p_end = tc.value(tc.last_parameter());
        let v0 = &vertices[edge.vertex_indices[0]];
        let v1 = &vertices[edge.vertex_indices[1]];
        let d0_start = ((v0.point[0] - p_start.x()).powi(2) + (v0.point[1] - p_start.y()).powi(2)
            + (v0.point[2] - p_start.z()).powi(2)).sqrt();
        let d1_start = ((v1.point[0] - p_start.x()).powi(2) + (v1.point[1] - p_start.y()).powi(2)
            + (v1.point[2] - p_start.z()).powi(2)).sqrt();
        let d0_end = ((v0.point[0] - p_end.x()).powi(2) + (v0.point[1] - p_end.y()).powi(2)
            + (v0.point[2] - p_end.z()).powi(2)).sqrt();
        let d1_end = ((v1.point[0] - p_end.x()).powi(2) + (v1.point[1] - p_end.y()).powi(2)
            + (v1.point[2] - p_end.z()).powi(2)).sqrt();
        // Check both possible assignments: (v0↔start, v1↔end) or (v0↔end, v1↔start)
        let assignment_a = d0_start.max(d1_end); // v0 at start, v1 at end
        let assignment_b = d0_end.max(d1_start); // v0 at end, v1 at start
        let max_gap = assignment_a.min(assignment_b); // best of two assignments
        if max_gap > 0.3 {
            // Large gap — re-project vertices onto the original IntSS curve.
            // The untrimmed curve (typically a circle for sphere-plane or
            // cylinder-plane intersections) should extend through both vertices.
            let v0_pt = gp::Pnt::new_real3(v0.point[0], v0.point[1], v0.point[2]);
            let v1_pt = gp::Pnt::new_real3(v1.point[0], v1.point[1], v1.point[2]);
            let proj0 = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(&v0_pt, curve_handle);
            let proj1 = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(&v1_pt, curve_handle);
            if proj0.nb_points() > 0 && proj1.nb_points() > 0
                && proj0.lower_distance() < 0.3
                && proj1.lower_distance() < 0.3
            {
                // Both vertices are close to the original curve — re-trim directly.
                let t0_new = proj0.lower_distance_parameter();
                let t1_new = proj1.lower_distance_parameter();

                // Handle closed curves (circles): choose the correct arc
                let param_span = curve.last_parameter() - curve.first_parameter();
                let is_closed = (param_span - 2.0 * std::f64::consts::PI).abs() < 1e-6;
                let (t_lo_new, t_hi_new) = if t0_new < t1_new {
                    (t0_new, t1_new)
                } else {
                    (t1_new, t0_new)
                };
                let (t_start_new, t_end_new) = if is_closed && edge.mesh_boundary_vertices.len() > 2 {
                    let period = param_span;
                    let step = (edge.mesh_boundary_vertices.len() / 5).max(1);
                    let mut in_direct = 0;
                    let mut in_complement = 0;
                    for i in (0..edge.mesh_boundary_vertices.len()).step_by(step) {
                        let vi = edge.mesh_boundary_vertices[i];
                        let v = &mesh.vertices[vi];
                        if let Ok(t_m) = project_point_on_curve(&[v.x, v.y, v.z], curve_handle) {
                            if t_m >= t_lo_new && t_m <= t_hi_new {
                                in_direct += 1;
                            } else {
                                in_complement += 1;
                            }
                        }
                    }
                    if in_complement > in_direct {
                        (t_hi_new, t_lo_new + period)
                    } else {
                        (t_lo_new, t_hi_new)
                    }
                } else if is_closed {
                    let direct_span = t_hi_new - t_lo_new;
                    if direct_span > param_span / 2.0 {
                        (t_hi_new, t_lo_new + param_span)
                    } else {
                        (t_lo_new, t_hi_new)
                    }
                } else {
                    (t_lo_new, t_hi_new)
                };

                let retrimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(
                    curve_handle, t_start_new, t_end_new,
                );
                if config.verbose {
                    eprintln!(
                        "    Edge (f{fi0},f{fi1}): re-projected onto IntSS curve, gap {max_gap:.3}mm -> proj_dist=[{:.3},{:.3}]mm",
                        proj0.lower_distance(), proj1.lower_distance()
                    );
                }
                geom::TrimmedCurve::to_handle(retrimmed).to_handle_curve()
            } else {
                // At least one vertex is far from the IntSS curve. If one of
                // the surfaces is a sphere, construct a great circle arc on the
                // sphere between the two vertices. This ensures the 3D curve
                // lies on the sphere (pcurve consistency) while still connecting
                // the correct vertices.
                let ss0 = &stage2.selected_surfaces[face_descriptors[fi0].selected_surface_idx];
                let ss1 = &stage2.selected_surfaces[face_descriptors[fi1].selected_surface_idx];
                let sphere_hyp = match (ss0, ss1) {
                    (SelectedSurface::Spherical(idx), _) | (_, SelectedSurface::Spherical(idx)) => {
                        Some(&stage2.spherical_hypotheses[*idx])
                    }
                    _ => None,
                };
                if let Some(sph) = sphere_hyp {
                    // Construct great circle arc on the sphere from v0 to v1
                    let center = sph.center;
                    let radius = sph.radius;
                    let cv0 = [
                        v0.point[0] - center[0],
                        v0.point[1] - center[1],
                        v0.point[2] - center[2],
                    ];
                    let cv1 = [
                        v1.point[0] - center[0],
                        v1.point[1] - center[1],
                        v1.point[2] - center[2],
                    ];
                    let cross = [
                        cv0[1] * cv1[2] - cv0[2] * cv1[1],
                        cv0[2] * cv1[0] - cv0[0] * cv1[2],
                        cv0[0] * cv1[1] - cv0[1] * cv1[0],
                    ];
                    let cross_len = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
                    if cross_len > 1e-12 {
                        let normal = [cross[0] / cross_len, cross[1] / cross_len, cross[2] / cross_len];
                        let center_pnt = gp::Pnt::new_real3(center[0], center[1], center[2]);
                        let normal_dir = gp::Dir::new_real3(normal[0], normal[1], normal[2]);
                        let ax2 = gp::Ax2::new_pnt_dir(&center_pnt, &normal_dir);
                        let circle = geom::Circle::new_ax2_real(&ax2, radius);
                        let circle_handle = geom::Circle::to_handle(circle).to_handle_curve();

                        let t0_c = project_point_on_curve(&v0.point, &circle_handle)?;
                        let t1_c = project_point_on_curve(&v1.point, &circle_handle)?;

                        // Choose the shorter arc
                        let (t_lo_c, t_hi_c) = select_arc_parameters(
                            t0_c, t1_c,
                            &edge.mesh_boundary_vertices,
                            &circle_handle,
                            mesh,
                        );

                        let arc_trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(
                            &circle_handle, t_lo_c, t_hi_c,
                        );
                        if config.verbose {
                            eprintln!(
                                "    Edge (f{fi0},f{fi1}): curve-vertex gap {max_gap:.3}mm, using sphere great circle arc"
                            );
                        }
                        geom::TrimmedCurve::to_handle(arc_trimmed).to_handle_curve()
                    } else {
                        // Vertices are collinear with sphere center — degenerate
                        let dx = v1.point[0] - v0.point[0];
                        let dy = v1.point[1] - v0.point[1];
                        let dz = v1.point[2] - v0.point[2];
                        let len = (dx * dx + dy * dy + dz * dz).sqrt();
                        if len > 1e-15 {
                            let p0 = gp::Pnt::new_real3(v0.point[0], v0.point[1], v0.point[2]);
                            let dir = gp::Dir::new_real3(dx / len, dy / len, dz / len);
                            let line = geom::Line::new_pnt_dir(&p0, &dir);
                            let line_handle = geom::Line::to_handle(line).to_handle_curve();
                            let line_trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(
                                &line_handle, 0.0, len,
                            );
                            geom::TrimmedCurve::to_handle(line_trimmed).to_handle_curve()
                        } else {
                            trimmed_handle
                        }
                    }
                } else {
                    // Neither surface is a sphere — use straight line
                    let dx = v1.point[0] - v0.point[0];
                    let dy = v1.point[1] - v0.point[1];
                    let dz = v1.point[2] - v0.point[2];
                    let len = (dx * dx + dy * dy + dz * dz).sqrt();
                    if len > 1e-15 {
                        let p0 = gp::Pnt::new_real3(v0.point[0], v0.point[1], v0.point[2]);
                        let dir = gp::Dir::new_real3(dx / len, dy / len, dz / len);
                        let line = geom::Line::new_pnt_dir(&p0, &dir);
                        let line_handle = geom::Line::to_handle(line).to_handle_curve();
                        let line_trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(
                            &line_handle, 0.0, len,
                        );
                        if config.verbose {
                            eprintln!("    Edge (f{fi0},f{fi1}): curve-vertex gap {max_gap:.3}mm, using vertex line");
                        }
                        geom::TrimmedCurve::to_handle(line_trimmed).to_handle_curve()
                    } else {
                        trimmed_handle
                    }
                }
            }
        } else {
            trimmed_handle
        }
    } else {
        trimmed_handle
    };

    // Ground-truth validation: the mesh boundary chain lies on the true
    // junction within chordal tolerance, so the trimmed curve must pass
    // close to the chain. When it doesn't — the classic complementary-arc
    // trap on a periodic junction, a degenerate IntSS result with an odd
    // parameterization, or a straight-line stand-in for a curved junction —
    // rebuild the curve as the circle through the chain and re-trim with a
    // chain-supported arc. The chain is the one representation available
    // for every junction, so this guard is uniform across surface pairs.
    let final_curve = repair_edge_curve_from_boundary(
        final_curve, edge, vertices, mesh, config, fi0, fi1,
    );

    // Closed-loop edges must carry a PERIODIC curve: stage 3.4's seam
    // planning re-anchors them with MakeEdge over [t, t + period] at an
    // arbitrary seam azimuth. Analytic closed junctions (circles) already
    // are; IntSS junctions against a fitted B-spline surface arrive as
    // closed-but-non-periodic trimmed curves and are converted (exact —
    // probed at ~1e-14).
    let final_curve = if edge.vertex_indices[0] == usize::MAX
        && edge.vertex_indices[1] == usize::MAX
    {
        periodicize_closed_curve(final_curve, config, fi0, fi1)
    } else {
        final_curve
    };

    edge.curve_3d = Some(final_curve);
    Ok(())
}

/// Convert a closed-loop edge's closed-but-non-periodic curve into a
/// periodic Geom_BSplineCurve (same geometry, same parameter range). Keeps
/// the input curve when it is already periodic, is not actually closed, or
/// the conversion fails.
fn periodicize_closed_curve(
    curve: OwnedPtr<geom::Handle_Geom_Curve>,
    config: &Config,
    fi0: usize,
    fi1: usize,
) -> OwnedPtr<geom::Handle_Geom_Curve> {
    {
        let c = curve.get();
        if c.is_periodic() {
            return curve;
        }
        let p = c.value(c.first_parameter());
        let q = c.value(c.last_parameter());
        let gap = ((p.x() - q.x()).powi(2) + (p.y() - q.y()).powi(2) + (p.z() - q.z()).powi(2))
            .sqrt();
        if gap > config.vertex_tolerance_mm {
            return curve; // not closed; nothing to periodicize
        }
    }
    let converted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut bsc = geom_convert::curve_to_b_spline_curve_handlegeomcurve_parameterisationtype(
            &curve,
            convert::ParameterisationType::Tgtthetaover2,
        );
        bsc.get_mut().set_periodic();
        bsc.to_handle_curve()
    }));
    match converted {
        Ok(periodic) => {
            if config.verbose {
                eprintln!(
                    "    Edge (f{fi0},f{fi1}): closed junction curve made periodic \
                     for seam re-anchoring"
                );
            }
            periodic
        }
        Err(_) => {
            if !config.quiet {
                eprintln!(
                    "  Warning: edge (f{fi0},f{fi1}): could not make the closed junction \
                     curve periodic; seam re-anchoring may keep its original start"
                );
            }
            curve
        }
    }
}

/// Maximum distance from a spread of the edge's mesh boundary vertices to
/// the (trimmed) curve.
fn boundary_chain_deviation(
    curve: &geom::Handle_Geom_Curve,
    mesh_boundary_vertices: &[usize],
    mesh: &ConnectedMesh,
) -> f64 {
    let n = mesh_boundary_vertices.len();
    let step = (n / 9).max(1);
    let mut dev = 0.0_f64;
    for i in (0..n).step_by(step).chain([n - 1]) {
        let v = &mesh.vertices[mesh_boundary_vertices[i]];
        let pt = gp::Pnt::new_real3(v.x, v.y, v.z);
        let proj = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(&pt, curve);
        if proj.nb_points() == 0 {
            // No projection: measure against the endpoints instead.
            let c = curve.get();
            let (p, q) = (c.value(c.first_parameter()), c.value(c.last_parameter()));
            let dp = ((v.x - p.x()).powi(2) + (v.y - p.y()).powi(2) + (v.z - p.z()).powi(2)).sqrt();
            let dq = ((v.x - q.x()).powi(2) + (v.y - q.y()).powi(2) + (v.z - q.z()).powi(2)).sqrt();
            dev = dev.max(dp.min(dq));
        } else {
            dev = dev.max(proj.lower_distance());
        }
    }
    dev
}

/// Validate an edge curve against the mesh boundary chain; on failure, try
/// the circle through three spread chain vertices and a B-spline curve fit
/// through the whole chain, keeping whichever candidate deviates least
/// (genuinely freeform junctions — e.g. a fitted patch meeting a plane
/// along an elliptical profile — are neither lines nor circles).
fn repair_edge_curve_from_boundary(
    curve: OwnedPtr<geom::Handle_Geom_Curve>,
    edge: &ReconEdge,
    vertices: &[BRepVertex],
    mesh: &ConnectedMesh,
    config: &Config,
    fi0: usize,
    fi1: usize,
) -> OwnedPtr<geom::Handle_Geom_Curve> {
    if edge.mesh_boundary_vertices.len() < 3 {
        return curve;
    }
    let tol = config.surface_tolerance_mm;
    let dev = boundary_chain_deviation(&curve, &edge.mesh_boundary_vertices, mesh);
    if dev <= tol {
        return curve;
    }

    let chain = &edge.mesh_boundary_vertices;
    let has_endpoints =
        edge.vertex_indices[0] != usize::MAX && edge.vertex_indices[1] != usize::MAX;

    // Trim a full repair curve to the edge's endpoint vertices (chain-
    // supported side for periodic candidates), or keep it whole for
    // closed-loop edges.
    let trim_to_endpoints = |full: &geom::Handle_Geom_Curve| -> Option<OwnedPtr<geom::Handle_Geom_Curve>> {
        if !has_endpoints {
            return None;
        }
        let v0 = &vertices[edge.vertex_indices[0]];
        let v1 = &vertices[edge.vertex_indices[1]];
        match (
            project_point_on_curve(&v0.point, full),
            project_point_on_curve(&v1.point, full),
        ) {
            (Ok(t0), Ok(t1)) => {
                let (t_lo, t_hi) = candidate_trim_parameters(t0, t1, chain, full, mesh)?;
                let arc = geom::TrimmedCurve::new_handlegeomcurve_real2(full, t_lo, t_hi);
                Some(geom::TrimmedCurve::to_handle(arc).to_handle_curve())
            }
            _ => None,
        }
    };

    let mut best_dev = dev;
    let mut best: Option<(OwnedPtr<geom::Handle_Geom_Curve>, &'static str)> = None;

    // Candidate 1: the circle through three spread chain vertices.
    let pnt_at = |i: usize| {
        let v = &mesh.vertices[chain[i]];
        gp::Pnt::new_real3(v.x, v.y, v.z)
    };
    let (p0, p1, p2) = (pnt_at(0), pnt_at(chain.len() / 2), pnt_at(chain.len() - 1));
    if let Some(circle) = reconstruct_full_circle(&p0, &p1, &p2) {
        let candidate = if has_endpoints {
            trim_to_endpoints(&circle)
        } else {
            Some(circle)
        };
        if let Some(c) = candidate {
            let d = boundary_chain_deviation(&c, chain, mesh);
            if d < best_dev {
                best_dev = d;
                best = Some((c, "circle"));
            }
        }
    }

    // Candidate 2: a B-spline curve through the chain vertices — a
    // least-squares fit trimmed to the endpoints for open chains, a
    // periodic C2 interpolation for closed chains (a closed freeform
    // junction, e.g. a cap meeting a fitted band, that IntSS missed). The
    // chain is ground truth for any junction, so an ordered fit through it
    // handles freeform junction curves no analytic construction covers.
    if has_endpoints && chain.len() >= 4 {
        let mut pts = n_collection::Array1_gp_Pnt::new_int2(1, chain.len() as i32);
        for (k, &vi) in chain.iter().enumerate() {
            let v = &mesh.vertices[vi];
            let p = gp::Pnt::new_real3(v.x, v.y, v.z);
            pts.set_value(k as i32 + 1, &p);
        }
        let fit_tol = (config.surface_tolerance_mm / 10.0).max(config.vertex_tolerance_mm * 10.0);
        let fit = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            geom_api::PointsToBSpline::new_arr1_pnt_int2_shape_real(
                &pts,
                3,
                8,
                geom_abs::Shape::C2,
                fit_tol,
            )
        }));
        if let Ok(fit) = fit {
            if fit.is_done() {
                let full = fit.curve().to_handle_curve();
                if let Some(c) = trim_to_endpoints(&full) {
                    let d = boundary_chain_deviation(&c, chain, mesh);
                    if d < best_dev {
                        best_dev = d;
                        best = Some((c, "B-spline"));
                    }
                }
            }
        }
    } else if !has_endpoints {
        if let Some(c) = closed_chain_periodic_bspline(chain, mesh, config) {
            let d = boundary_chain_deviation(&c, chain, mesh);
            if d < best_dev {
                best_dev = d;
                best = Some((c, "periodic B-spline"));
            }
        }
    }

    match best {
        Some((repaired, kind)) => {
            if !config.quiet {
                eprintln!(
                    "  Edge (f{fi0},f{fi1}): curve was {dev:.3} mm from its mesh boundary \
                     chain; repaired with the chain's {kind} ({best_dev:.2e} mm)"
                );
                if best_dev > tol {
                    eprintln!(
                        "  Warning: edge (f{fi0},f{fi1}) repaired curve still {best_dev:.3} mm \
                         from its mesh boundary chain"
                    );
                }
            }
            repaired
        }
        None => {
            if !config.quiet {
                eprintln!(
                    "  Warning: edge (f{fi0},f{fi1}) curve is {dev:.3} mm from its mesh \
                     boundary chain and no chain-based repair improved it — keeping the curve"
                );
            }
            curve
        }
    }
}

/// Periodic C2 B-spline through a closed mesh boundary chain, via
/// GeomAPI_Interpolate with the periodic flag (the closing point is NOT
/// repeated — the wrap chord is implied). The chain is ground truth for any
/// junction, so this handles closed freeform junctions (e.g. a cap meeting
/// a fitted band) that no analytic construction covers.
fn closed_chain_periodic_bspline(
    chain: &[usize],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Option<OwnedPtr<geom::Handle_Geom_Curve>> {
    if chain.len() < 4 {
        return None;
    }
    let mut loop_chain: &[usize] = chain;
    if loop_chain.first() == loop_chain.last() {
        loop_chain = &loop_chain[..loop_chain.len() - 1];
    }
    let mut pts_h = n_collection::HArray1_gp_Pnt::new_int2(1, loop_chain.len() as i32);
    for (k, &vi) in loop_chain.iter().enumerate() {
        let v = &mesh.vertices[vi];
        let p = gp::Pnt::new_real3(v.x, v.y, v.z);
        pts_h.set_value(k as i32 + 1, &p);
    }
    let pts_handle = n_collection::HArray1_gp_Pnt::to_handle(pts_h);
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let mut interp = geom_api::Interpolate::new_harr1_pnt_bool_real(
            &pts_handle,
            true,
            (config.vertex_tolerance_mm * 0.01).max(1e-9),
        );
        interp.perform();
        if interp.is_done() {
            Some(interp.curve().to_handle_curve())
        } else {
            None
        }
    }))
    .ok()
    .flatten()
}

/// Last-resort edge curve when both the IntSS and tangent constructions
/// fail: build it from the mesh boundary chain — the ground truth for any
/// junction. Starts from the vertex-to-vertex line and lets the chain
/// repair replace it with the chain's circle when the junction is curved.
/// Stores the best curve on the edge and returns its chain deviation.
fn edge_curve_from_boundary_chain(
    edge: &mut ReconEdge,
    vertices: &[BRepVertex],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<f64, String> {
    if edge.vertex_indices[0] == usize::MAX || edge.vertex_indices[1] == usize::MAX {
        // Closed-loop edge — no vertices to anchor a line; interpolate the
        // closed chain directly. Already periodic, as stage 3.4's seam
        // re-anchoring requires. A 3-vertex loop is below the periodic
        // interpolation's minimum, but the circle through its three points
        // is exact (and periodic).
        let curve = closed_chain_periodic_bspline(&edge.mesh_boundary_vertices, mesh, config)
            .or_else(|| {
                let mut chain: &[usize] = &edge.mesh_boundary_vertices;
                if chain.first() == chain.last() && chain.len() > 1 {
                    chain = &chain[..chain.len() - 1];
                }
                if chain.len() != 3 {
                    return None;
                }
                let pnt_at = |i: usize| {
                    let v = &mesh.vertices[chain[i]];
                    gp::Pnt::new_real3(v.x, v.y, v.z)
                };
                reconstruct_full_circle(&pnt_at(0), &pnt_at(1), &pnt_at(2))
            })
            .ok_or_else(|| {
                "closed-loop edge has no vertex endpoints and chain interpolation failed"
                    .to_string()
            })?;
        let dev = boundary_chain_deviation(&curve, &edge.mesh_boundary_vertices, mesh);
        edge.curve_3d = Some(curve);
        return Ok(dev);
    }
    let v0 = &vertices[edge.vertex_indices[0]];
    let v1 = &vertices[edge.vertex_indices[1]];
    let d = [
        v1.point[0] - v0.point[0],
        v1.point[1] - v0.point[1],
        v1.point[2] - v0.point[2],
    ];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if len < 1e-15 {
        return Err("edge has zero-length vertex span".to_string());
    }
    let p0 = gp::Pnt::new_real3(v0.point[0], v0.point[1], v0.point[2]);
    let dir = gp::Dir::new_real3(d[0] / len, d[1] / len, d[2] / len);
    let line = geom::Line::new_pnt_dir(&p0, &dir);
    let line_handle = geom::Line::to_handle(line).to_handle_curve();
    let trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(&line_handle, 0.0, len);
    let trimmed_handle = geom::TrimmedCurve::to_handle(trimmed).to_handle_curve();

    let [fi0, fi1] = edge.face_indices;
    let repaired =
        repair_edge_curve_from_boundary(trimmed_handle, edge, vertices, mesh, config, fi0, fi1);
    let dev = boundary_chain_deviation(&repaired, &edge.mesh_boundary_vertices, mesh);
    edge.curve_3d = Some(repaired);
    Ok(dev)
}

/// Select the intersection curve closest to the mesh boundary vertices.
/// Returns the 1-indexed line number.
///
/// Uses curve midpoint evaluation rather than projection because
/// ProjectPointOnCurve can fail on certain IntSS curve representations.
fn select_closest_curve(
    int_ss: &geom_api::IntSS,
    boundary_vertices: &[usize],
    mesh: &ConnectedMesh,
) -> i32 {
    let nb_lines = int_ss.nb_lines();
    let mut best_idx = 0_i32; // 0 = no valid curve found
    let mut best_total_dist = f64::MAX;

    // Compute centroid of boundary vertices
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for &vi in boundary_vertices {
        let v = &mesh.vertices[vi];
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    let n = boundary_vertices.len() as f64;
    cx /= n;
    cy /= n;
    cz /= n;

    for line_idx in 1..=nb_lines {
        let curve_handle = int_ss.line(line_idx);
        let c = curve_handle.get();
        let fp = c.first_parameter();
        let lp = c.last_parameter();

        // For curves with finite parameter ranges, check for extreme values
        // that could cause overflow (e.g., hyperbolas from cone intersections
        // whose Cosh overflows Standard_Real). Lines have infinite ranges but
        // are safe — we'll use projection instead of parameter sampling.
        let use_projection = !fp.is_finite() || !lp.is_finite()
            || fp.abs() > 1e6 || lp.abs() > 1e6;

        if use_projection {
            // Project centroid onto curve to find distance
            let centroid_pt = gp::Pnt::new_real3(cx, cy, cz);
            let proj = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(
                &centroid_pt, curve_handle,
            );
            if proj.nb_points() > 0 {
                let dist = proj.lower_distance();
                if dist < best_total_dist {
                    best_total_dist = dist;
                    best_idx = line_idx;
                }
            }
            // If projection fails, skip this curve
        } else {
            // Sample 5 points along the curve and compute min distance to centroid.
            let mut min_dist = f64::MAX;
            for i in 0..5 {
                let t = fp + (lp - fp) * (i as f64 / 4.0);
                let p = c.value(t);
                let dx = p.x() - cx;
                let dy = p.y() - cy;
                let dz = p.z() - cz;
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                min_dist = min_dist.min(dist);
            }

            if min_dist < best_total_dist {
                best_total_dist = min_dist;
                best_idx = line_idx;
            }
        }
    }

    best_idx
}

/// Project a 3D point onto a curve and return the curve parameter.
/// Reconstruct a full Geom_Circle from 3 points on a circular arc.
///
/// Used when IntSS returns a partial arc (e.g., semicircle) for a closed-loop
/// edge. This happens when the intersection plane contains the sphere axis,
/// causing the circle to pass through the sphere's UV poles. IntSS splits such
/// circles into two arcs in UV space.
///
/// Returns a Handle(Geom_Curve) for a full circle [0, 2π], or None if the
/// points are nearly collinear.
fn reconstruct_full_circle(
    p0: &gp::Pnt,
    p1: &gp::Pnt,
    p2: &gp::Pnt,
) -> Option<OwnedPtr<geom::Handle_Geom_Curve>> {
    // Vectors from p0 to p1 and p2
    let a = [p1.x() - p0.x(), p1.y() - p0.y(), p1.z() - p0.z()];
    let b = [p2.x() - p0.x(), p2.y() - p0.y(), p2.z() - p0.z()];

    // Normal = cross(a, b)
    let n = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let n2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    if !n2.is_finite() || n2 < 1e-30 {
        // Collinear points, or coordinates so extreme the cross product
        // overflowed (e.g. sampling a degenerate IntSS line at its infinite
        // parameter bounds, ~1e100) — n/|n| would then be 0/inf = 0 and
        // gp_Dir panics on the zero-norm vector.
        return None;
    }

    // Circumcenter = p0 + (|b|²*(n×a) + |a|²*(b×n)) / (2*|n|²)
    let a2 = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    let b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    let nxa = [
        n[1] * a[2] - n[2] * a[1],
        n[2] * a[0] - n[0] * a[2],
        n[0] * a[1] - n[1] * a[0],
    ];
    let bxn = [
        b[1] * n[2] - b[2] * n[1],
        b[2] * n[0] - b[0] * n[2],
        b[0] * n[1] - b[1] * n[0],
    ];
    let denom = 2.0 * n2;
    let cx = p0.x() + (b2 * nxa[0] + a2 * bxn[0]) / denom;
    let cy = p0.y() + (b2 * nxa[1] + a2 * bxn[1]) / denom;
    let cz = p0.z() + (b2 * nxa[2] + a2 * bxn[2]) / denom;
    let radius = ((cx - p0.x()).powi(2) + (cy - p0.y()).powi(2) + (cz - p0.z()).powi(2)).sqrt();
    if !(cx.is_finite() && cy.is_finite() && cz.is_finite() && radius.is_finite()) {
        return None;
    }

    let n_len = n2.sqrt();
    let center = gp::Pnt::new_real3(cx, cy, cz);
    let normal = gp::Dir::new_real3(n[0] / n_len, n[1] / n_len, n[2] / n_len);
    let ax2 = gp::Ax2::new_pnt_dir(&center, &normal);
    let circle = geom::Circle::new_ax2_real(&ax2, radius);
    Some(geom::Circle::to_handle(circle).to_handle_curve())
}

fn project_point_on_curve(
    point: &[f64; 3],
    curve: &geom::Handle_Geom_Curve,
) -> Result<f64, String> {
    let pt = gp::Pnt::new_real3(point[0], point[1], point[2]);
    let projector = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(&pt, curve);
    if projector.nb_points() > 0 {
        Ok(projector.lower_distance_parameter())
    } else {
        Err(format!(
            "Failed to project point [{:.6}, {:.6}, {:.6}] onto curve",
            point[0], point[1], point[2]
        ))
    }
}

/// Compute the edge curve for a tangent edge.
///
/// For tangent edges, `GeomAPI_IntSS` often fails or produces degenerate results.
/// Instead, construct the curve analytically based on the surface pair types:
/// - Plane-cylinder tangent: line parallel to cylinder axis at the tangent point
/// - Sphere-cylinder tangent: circle arc on the sphere, in a plane perpendicular
///   to the cylinder axis through the sphere center's projection on the axis
/// - Cylinder-cylinder tangent: line parallel to shared axis direction
/// - Plane-sphere tangent: degenerate (single point) — construct line from vertices
fn compute_tangent_edge_curve(
    edge: &mut ReconEdge,
    face_descriptors: &[FaceDescriptor],
    vertices: &[BRepVertex],
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    let [fi0, fi1] = edge.face_indices;
    let ss0 = &stage2.selected_surfaces[face_descriptors[fi0].selected_surface_idx];
    let ss1 = &stage2.selected_surfaces[face_descriptors[fi1].selected_surface_idx];

    // Determine which surface types we have
    let has_sphere_cylinder = matches!(
        (ss0, ss1),
        (SelectedSurface::Spherical(_), SelectedSurface::Conical(_))
            | (SelectedSurface::Conical(_), SelectedSurface::Spherical(_))
    );

    let is_closed_loop = edge.vertex_indices[0] == usize::MAX || edge.vertex_indices[1] == usize::MAX;

    if has_sphere_cylinder {
        // Sphere-cone tangencies are NOT great circles on the sphere (the
        // tangent circle sits at radius R·cos α, offset along the cone axis),
        // so handle them with the generic axis-perpendicular circle
        // construction from mesh boundary vertices.
        let cone_axis = match (ss0, ss1) {
            (SelectedSurface::Conical(c), _) | (_, SelectedSurface::Conical(c))
                if !stage2.conical_hypotheses[*c].is_cylinder() =>
            {
                let h = &stage2.conical_hypotheses[*c];
                Some((h.apex, h.axis_direction))
            }
            _ => None,
        };
        if is_closed_loop {
            if let Some((axis_origin, axis_dir)) = cone_axis {
                return compute_tangent_edge_curve_axis_circle_closed(
                    edge, &axis_origin, &axis_dir, stage2, config,
                );
            }
            return compute_tangent_edge_curve_sphere_cylinder_closed(
                edge, ss0, ss1, stage2, config,
            );
        }
        // Open arcs: the circle is constructed from the sphere geometry and
        // the two endpoint vertices, which works for both cylinders and
        // near-cylindrical cones (e.g. fillets fitted as slight cones).
        let v0 = &vertices[edge.vertex_indices[0]];
        let v1 = &vertices[edge.vertex_indices[1]];
        return compute_tangent_edge_curve_sphere_cylinder(
            edge, v0, v1, ss0, ss1, stage2, config,
        );
    }

    // Check for torus involvement
    let has_torus = matches!(ss0, SelectedSurface::Toroidal(_))
        || matches!(ss1, SelectedSurface::Toroidal(_));

    if has_torus {
        if is_closed_loop {
            return compute_tangent_edge_curve_torus_closed(
                edge, ss0, ss1, stage2, config,
            );
        }
        let v0 = &vertices[edge.vertex_indices[0]];
        let v1 = &vertices[edge.vertex_indices[1]];
        return compute_tangent_edge_curve_torus_arc(
            edge, v0, v1, ss0, ss1, stage2, config,
        );
    }

    // Non-sphere-cylinder, non-torus tangent edges require vertex endpoints
    if is_closed_loop {
        return Err("tangent edge has no vertex endpoints".to_string());
    }
    let v0 = &vertices[edge.vertex_indices[0]];
    let v1 = &vertices[edge.vertex_indices[1]];

    // For plane-cylinder, cylinder-cylinder, and fallback cases: construct a line
    compute_tangent_edge_curve_line(edge, v0, v1, ss0, ss1, stage2, config)
}

/// Construct a line-based tangent edge curve (plane-cylinder, cylinder-cylinder, fallback).
fn compute_tangent_edge_curve_line(
    edge: &mut ReconEdge,
    v0: &BRepVertex,
    v1: &BRepVertex,
    ss0: &SelectedSurface,
    ss1: &SelectedSurface,
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    // Get the line direction. For plane-CYLINDER (and cylinder-cylinder)
    // tangencies the line is parallel to the cylinder axis. For plane-CONE
    // tangencies the line is a cone GENERATOR — tilted from the axis by the
    // half-angle — and the two endpoint vertices lie exactly on it, so the
    // vertex-to-vertex direction is used (radius_at_apex is 0 for cones, so
    // the analytic cylinder construction would collapse onto the axis).
    let conical_is_cylinder = |idx: &usize| stage2.conical_hypotheses[*idx].is_cylinder();
    let line_dir = match (ss0, ss1) {
        (SelectedSurface::Conical(idx), SelectedSurface::Planar(_))
        | (SelectedSurface::Planar(_), SelectedSurface::Conical(idx))
            if conical_is_cylinder(idx) =>
        {
            let hyp = &stage2.conical_hypotheses[*idx];
            hyp.axis_direction
        }
        (SelectedSurface::Conical(idx), SelectedSurface::Conical(_))
            if conical_is_cylinder(idx) =>
        {
            let hyp = &stage2.conical_hypotheses[*idx];
            hyp.axis_direction
        }
        _ => {
            // Cone generators and general fallback: direction between the
            // vertex endpoints.
            let dx = v1.point[0] - v0.point[0];
            let dy = v1.point[1] - v0.point[1];
            let dz = v1.point[2] - v0.point[2];
            let len = (dx * dx + dy * dy + dz * dz).sqrt();
            if len < 1e-15 {
                return Err("tangent edge has zero-length vertex span".to_string());
            }
            [dx / len, dy / len, dz / len]
        }
    };

    // Compute the tangent point analytically for plane-cylinder pairs.
    let tangent_point = match (ss0, ss1) {
        (SelectedSurface::Conical(cyl_idx), SelectedSurface::Planar(plane_idx))
        | (SelectedSurface::Planar(plane_idx), SelectedSurface::Conical(cyl_idx))
            if conical_is_cylinder(cyl_idx) =>
        {
            let cyl = &stage2.conical_hypotheses[*cyl_idx];
            let plane = &stage2.planar_hypotheses[*plane_idx];
            let a = cyl.axis_direction;
            let n = plane.normal;
            // Component of plane normal perpendicular to cylinder axis
            let dot = n[0] * a[0] + n[1] * a[1] + n[2] * a[2];
            let perp = [n[0] - dot * a[0], n[1] - dot * a[1], n[2] - dot * a[2]];
            let perp_len = (perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]).sqrt();
            if perp_len < 1e-12 {
                return Err("plane normal parallel to cylinder axis".to_string());
            }
            let radial_unit = [perp[0] / perp_len, perp[1] / perp_len, perp[2] / perp_len];
            let sign = if cyl.convex { 1.0 } else { -1.0 };
            let v0_minus_o = [
                v0.point[0] - cyl.apex[0],
                v0.point[1] - cyl.apex[1],
                v0.point[2] - cyl.apex[2],
            ];
            let t_axis = v0_minus_o[0] * a[0] + v0_minus_o[1] * a[1] + v0_minus_o[2] * a[2];
            [
                cyl.apex[0] + t_axis * a[0] + sign * cyl.radius_at_apex * radial_unit[0],
                cyl.apex[1] + t_axis * a[1] + sign * cyl.radius_at_apex * radial_unit[1],
                cyl.apex[2] + t_axis * a[2] + sign * cyl.radius_at_apex * radial_unit[2],
            ]
        }
        _ => v0.point, // fallback: use vertex position
    };

    // Construct Geom_Line passing through the tangent point in the line direction
    let p0 = gp::Pnt::new_real3(tangent_point[0], tangent_point[1], tangent_point[2]);
    let dir = gp::Dir::new_real3(line_dir[0], line_dir[1], line_dir[2]);
    let line = geom::Line::new_pnt_dir(&p0, &dir);
    let line_handle = geom::Line::to_handle(line).to_handle_curve();

    // Trim to vertex endpoints. Line parameters are computed analytically
    // (t = (p − origin)·dir): GeomAPI_ProjectPointOnCurve can return zero
    // extrema for a point lying exactly on the line (e.g. v0 itself when
    // the line was anchored there for cylinder-cylinder tangencies).
    let line_param = |p: &[f64; 3]| -> f64 {
        (p[0] - tangent_point[0]) * line_dir[0]
            + (p[1] - tangent_point[1]) * line_dir[1]
            + (p[2] - tangent_point[2]) * line_dir[2]
    };
    let t_start = line_param(&v0.point);
    let t_end = line_param(&v1.point);
    let (t_lo, t_hi) = if t_start < t_end {
        (t_start, t_end)
    } else {
        (t_end, t_start)
    };

    let trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(&line_handle, t_lo, t_hi);
    let trimmed_handle = geom::TrimmedCurve::to_handle(trimmed).to_handle_curve();

    validate_tangent_curve(edge, &trimmed_handle, stage2, config)?;
    edge.curve_3d = Some(trimmed_handle);
    Ok(())
}


/// Construct a full-circle tangent edge curve for a closed-loop sphere-cylinder tangency.
///
/// This handles the case where a sphere meets a cylinder with no corner vertices
/// (e.g., a pill/capsule shape where hemispherical caps meet a cylinder body).
/// The tangent curve is a full circle perpendicular to the cylinder axis,
/// centered at the sphere center, with the sphere's radius.
fn compute_tangent_edge_curve_sphere_cylinder_closed(
    edge: &mut ReconEdge,
    ss0: &SelectedSurface,
    ss1: &SelectedSurface,
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    // Extract sphere and cylinder hypotheses
    let (sph_idx, cyl_idx) = match (ss0, ss1) {
        (SelectedSurface::Spherical(s), SelectedSurface::Conical(c)) => (*s, *c),
        (SelectedSurface::Conical(c), SelectedSurface::Spherical(s)) => (*s, *c),
        _ => return Err("not a sphere-cylinder pair".to_string()),
    };
    let sph = &stage2.spherical_hypotheses[sph_idx];
    let cyl = &stage2.conical_hypotheses[cyl_idx];
    let center = sph.center;
    let radius = sph.radius;

    // The circle plane normal is the cylinder axis direction.
    // For a tangent closed loop, the circle lies in a plane perpendicular
    // to the cylinder axis at the sphere center.
    let normal = cyl.axis_direction;

    let center_pnt = gp::Pnt::new_real3(center[0], center[1], center[2]);
    let normal_dir = gp::Dir::new_real3(normal[0], normal[1], normal[2]);
    let ax2 = gp::Ax2::new_pnt_dir(&center_pnt, &normal_dir);
    let circle = geom::Circle::new_ax2_real(&ax2, radius);
    let circle_handle = geom::Circle::to_handle(circle).to_handle_curve();

    // For a closed loop, use the full circle parameter range
    let curve = circle_handle.get();
    let fp = curve.first_parameter();
    let lp = curve.last_parameter();

    let trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(&circle_handle, fp, lp);
    let trimmed_handle = geom::TrimmedCurve::to_handle(trimmed).to_handle_curve();

    validate_tangent_curve(edge, &trimmed_handle, stage2, config)?;
    edge.curve_3d = Some(trimmed_handle);
    Ok(())
}


/// Construct a full-circle tangent edge curve perpendicular to a given axis.
///
/// The circle's axial position and radius are averaged from the edge's mesh
/// boundary vertices, which lie exactly on the junction circle. Used for
/// sphere-cone tangencies (and any coaxial tangency where the junction is a
/// circle around a known axis).
fn compute_tangent_edge_curve_axis_circle_closed(
    edge: &mut ReconEdge,
    axis_origin: &[f64; 3],
    axis_dir: &[f64; 3],
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    let n = edge.mesh_boundary_vertices.len() as f64;
    if n < 3.0 {
        return Err("closed tangent loop has too few boundary vertices".to_string());
    }

    let mut h_sum = 0.0;
    let mut r_sum = 0.0;
    for &vi in &edge.mesh_boundary_vertices {
        let v = &stage2.mesh.vertices[vi];
        let dp = [
            v.x - axis_origin[0],
            v.y - axis_origin[1],
            v.z - axis_origin[2],
        ];
        let h = dp[0] * axis_dir[0] + dp[1] * axis_dir[1] + dp[2] * axis_dir[2];
        let radial = [
            dp[0] - h * axis_dir[0],
            dp[1] - h * axis_dir[1],
            dp[2] - h * axis_dir[2],
        ];
        let r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
        h_sum += h;
        r_sum += r;
    }
    let h_avg = h_sum / n;
    let r_avg = r_sum / n;

    let circle_center = [
        axis_origin[0] + h_avg * axis_dir[0],
        axis_origin[1] + h_avg * axis_dir[1],
        axis_origin[2] + h_avg * axis_dir[2],
    ];

    let center_pnt = gp::Pnt::new_real3(circle_center[0], circle_center[1], circle_center[2]);
    let normal_dir = gp::Dir::new_real3(axis_dir[0], axis_dir[1], axis_dir[2]);
    let ax2 = gp::Ax2::new_pnt_dir(&center_pnt, &normal_dir);
    let circle = geom::Circle::new_ax2_real(&ax2, r_avg);
    let circle_handle = geom::Circle::to_handle(circle).to_handle_curve();

    let curve = circle_handle.get();
    let fp = curve.first_parameter();
    let lp = curve.last_parameter();

    let trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(&circle_handle, fp, lp);
    let trimmed_handle = geom::TrimmedCurve::to_handle(trimmed).to_handle_curve();

    validate_tangent_curve(edge, &trimmed_handle, stage2, config)?;
    edge.curve_3d = Some(trimmed_handle);
    Ok(())
}

/// Extract the torus hypothesis index from one of the two surfaces.
fn torus_idx_from_pair(
    ss0: &SelectedSurface,
    ss1: &SelectedSurface,
) -> Result<usize, String> {
    match (ss0, ss1) {
        (SelectedSurface::Toroidal(t), _) => Ok(*t),
        (_, SelectedSurface::Toroidal(t)) => Ok(*t),
        _ => Err("no torus in surface pair".to_string()),
    }
}

/// Construct a full-circle tangent edge curve for a closed-loop torus tangency.
///
/// Tangent curves involving a torus (with any coaxial surface: plane, cylinder, or
/// another torus) are circles perpendicular to the torus axis. The circle center
/// and radius are computed from the mesh boundary vertices.
fn compute_tangent_edge_curve_torus_closed(
    edge: &mut ReconEdge,
    ss0: &SelectedSurface,
    ss1: &SelectedSurface,
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    let tor_idx = torus_idx_from_pair(ss0, ss1)?;
    let tor = &stage2.toroidal_hypotheses[tor_idx];

    if edge.mesh_boundary_vertices.len() < 3 {
        return Err("torus closed loop has too few boundary vertices".to_string());
    }

    // The junction circle is perpendicular to *some* axis:
    // - Coaxial fillet junctions (torus fillet around a cylinder or against a
    //   plane): perpendicular to the TORUS axis.
    // - Tube-continuation junctions (a straight pipe leaving a pipe-elbow
    //   torus): a tube cross-section, perpendicular to the CYLINDER axis.
    // Try each candidate axis and keep the circle that best fits the mesh
    // boundary vertices.
    let mut candidates: Vec<([f64; 3], [f64; 3])> = vec![(tor.center, tor.axis_direction)];
    for ss in [ss0, ss1] {
        if let SelectedSurface::Conical(c) = ss {
            let hyp = &stage2.conical_hypotheses[*c];
            candidates.push((hyp.apex, hyp.axis_direction));
        }
    }

    let mut best: Option<(f64, [f64; 3], [f64; 3], f64)> = None; // (dev, center, axis, radius)
    for (origin, axis) in candidates {
        // Average axial position and radius of the boundary vertices.
        let mut h_sum = 0.0;
        let mut r_sum = 0.0;
        let n = edge.mesh_boundary_vertices.len() as f64;
        for &vi in &edge.mesh_boundary_vertices {
            let v = &stage2.mesh.vertices[vi];
            let dp = [v.x - origin[0], v.y - origin[1], v.z - origin[2]];
            let h = dp[0] * axis[0] + dp[1] * axis[1] + dp[2] * axis[2];
            let radial = [dp[0] - h * axis[0], dp[1] - h * axis[1], dp[2] - h * axis[2]];
            let r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            h_sum += h;
            r_sum += r;
        }
        let h_avg = h_sum / n;
        let r_avg = r_sum / n;
        let center = [
            origin[0] + h_avg * axis[0],
            origin[1] + h_avg * axis[1],
            origin[2] + h_avg * axis[2],
        ];
        // Max deviation of boundary vertices from this circle.
        let mut max_dev = 0.0_f64;
        for &vi in &edge.mesh_boundary_vertices {
            let v = &stage2.mesh.vertices[vi];
            let dp = [v.x - center[0], v.y - center[1], v.z - center[2]];
            let h = dp[0] * axis[0] + dp[1] * axis[1] + dp[2] * axis[2];
            let radial = [dp[0] - h * axis[0], dp[1] - h * axis[1], dp[2] - h * axis[2]];
            let r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            let dev = ((r - r_avg) * (r - r_avg) + h * h).sqrt();
            max_dev = max_dev.max(dev);
        }
        if best.as_ref().is_none_or(|(d, _, _, _)| max_dev < *d) {
            best = Some((max_dev, center, axis, r_avg));
        }
    }
    let (_, circle_center, axis, r_avg) = best.unwrap();

    let center_pnt = gp::Pnt::new_real3(circle_center[0], circle_center[1], circle_center[2]);
    let normal_dir = gp::Dir::new_real3(axis[0], axis[1], axis[2]);
    let ax2 = gp::Ax2::new_pnt_dir(&center_pnt, &normal_dir);
    let circle = geom::Circle::new_ax2_real(&ax2, r_avg);
    let circle_handle = geom::Circle::to_handle(circle).to_handle_curve();

    // Full circle for closed loop
    let curve = circle_handle.get();
    let fp = curve.first_parameter();
    let lp = curve.last_parameter();

    let trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(&circle_handle, fp, lp);
    let trimmed_handle = geom::TrimmedCurve::to_handle(trimmed).to_handle_curve();

    validate_tangent_curve(edge, &trimmed_handle, stage2, config)?;
    edge.curve_3d = Some(trimmed_handle);
    Ok(())
}

/// Construct a circle-arc tangent edge curve for a torus tangency with vertex endpoints.
///
/// The arc lies on a circle perpendicular to the torus axis. The circle center and
/// radius are computed from the torus hypothesis and vertex positions, then the arc
/// direction is determined from mesh boundary vertices.
fn compute_tangent_edge_curve_torus_arc(
    edge: &mut ReconEdge,
    v0: &BRepVertex,
    v1: &BRepVertex,
    ss0: &SelectedSurface,
    ss1: &SelectedSurface,
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    let tor_idx = torus_idx_from_pair(ss0, ss1)?;
    let tor = &stage2.toroidal_hypotheses[tor_idx];

    // The junction circle is perpendicular to *some* axis — the torus axis
    // for coaxial fillet junctions, or the adjacent cylinder/torus axis for
    // tube-continuation junctions (same reasoning as the closed-loop case).
    // A wrong axis choice projects both endpoint vertices to the SAME circle
    // parameter (they differ only along that axis), producing a degenerate
    // zero-span arc. Pick the candidate axis whose perpendicular circle best
    // fits the edge's mesh boundary vertices.
    let mut candidates: Vec<([f64; 3], [f64; 3])> = vec![(tor.center, tor.axis_direction)];
    for ss in [ss0, ss1] {
        match ss {
            SelectedSurface::Conical(c) => {
                let hyp = &stage2.conical_hypotheses[*c];
                candidates.push((hyp.apex, hyp.axis_direction));
            }
            SelectedSurface::Toroidal(t2) if *t2 != tor_idx => {
                let hyp = &stage2.toroidal_hypotheses[*t2];
                candidates.push((hyp.center, hyp.axis_direction));
            }
            _ => {}
        }
    }

    let boundary_pts: Vec<[f64; 3]> = edge
        .mesh_boundary_vertices
        .iter()
        .map(|&vi| {
            let v = &stage2.mesh.vertices[vi];
            [v.x, v.y, v.z]
        })
        .collect();

    let mut best: Option<(f64, [f64; 3], [f64; 3], f64)> = None; // (dev, center, axis, radius)
    for (origin, axis) in candidates {
        // Circle center and radius from the two endpoint vertices.
        let mut h_avg = 0.0;
        let mut r_avg = 0.0;
        for v in [v0, v1] {
            let dp = [
                v.point[0] - origin[0],
                v.point[1] - origin[1],
                v.point[2] - origin[2],
            ];
            let h = dp[0] * axis[0] + dp[1] * axis[1] + dp[2] * axis[2];
            let radial = [dp[0] - h * axis[0], dp[1] - h * axis[1], dp[2] - h * axis[2]];
            let r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            h_avg += h;
            r_avg += r;
        }
        h_avg /= 2.0;
        r_avg /= 2.0;
        if r_avg < 1e-9 {
            continue;
        }
        let cc = [
            origin[0] + h_avg * axis[0],
            origin[1] + h_avg * axis[1],
            origin[2] + h_avg * axis[2],
        ];
        // Max deviation of boundary vertices from this circle.
        let mut max_dev = 0.0_f64;
        for p in &boundary_pts {
            let dp = [p[0] - cc[0], p[1] - cc[1], p[2] - cc[2]];
            let h = dp[0] * axis[0] + dp[1] * axis[1] + dp[2] * axis[2];
            let radial = [dp[0] - h * axis[0], dp[1] - h * axis[1], dp[2] - h * axis[2]];
            let r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            let dev = ((r - r_avg) * (r - r_avg) + h * h).sqrt();
            max_dev = max_dev.max(dev);
        }
        if best.as_ref().is_none_or(|(d, _, _, _)| max_dev < *d) {
            best = Some((max_dev, cc, axis, r_avg));
        }
    }
    let Some((_, circle_center, axis, r_avg)) = best else {
        return Err("no valid axis candidate for torus tangent arc".to_string());
    };

    let center_pnt = gp::Pnt::new_real3(circle_center[0], circle_center[1], circle_center[2]);
    let normal_dir = gp::Dir::new_real3(axis[0], axis[1], axis[2]);
    let ax2 = gp::Ax2::new_pnt_dir(&center_pnt, &normal_dir);
    let circle = geom::Circle::new_ax2_real(&ax2, r_avg);
    let circle_handle = geom::Circle::to_handle(circle).to_handle_curve();

    // Project vertex endpoints onto the circle curve
    let t0 = project_point_on_curve(&v0.point, &circle_handle)?;
    let t1 = project_point_on_curve(&v1.point, &circle_handle)?;

    // Select the correct arc direction from mesh boundary vertices
    let (t_lo, t_hi) = select_arc_parameters(
        t0,
        t1,
        &edge.mesh_boundary_vertices,
        &circle_handle,
        &stage2.mesh,
    );

    let trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(&circle_handle, t_lo, t_hi);
    let trimmed_handle = geom::TrimmedCurve::to_handle(trimmed).to_handle_curve();

    validate_tangent_curve(edge, &trimmed_handle, stage2, config)?;
    edge.curve_3d = Some(trimmed_handle);
    Ok(())
}

/// Construct a circle-arc tangent edge curve for sphere-cylinder tangencies.
///
/// For a sphere tangent to a cylinder (e.g., at fillet corners of a rounded cube),
/// the tangent curve is a great circle arc on the sphere. Rather than using the
/// cylinder's (potentially imprecise) axis parameters, construct the circle from
/// the sphere geometry and vertex positions:
/// - Center = sphere center
/// - Radius = sphere radius
/// - The plane of the circle is determined by the sphere center and both vertices
fn compute_tangent_edge_curve_sphere_cylinder(
    edge: &mut ReconEdge,
    v0: &BRepVertex,
    v1: &BRepVertex,
    ss0: &SelectedSurface,
    ss1: &SelectedSurface,
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    // Extract sphere hypothesis
    let sph_idx = match (ss0, ss1) {
        (SelectedSurface::Spherical(s), SelectedSurface::Conical(_)) => *s,
        (SelectedSurface::Conical(_), SelectedSurface::Spherical(s)) => *s,
        _ => return Err("not a sphere-cylinder pair".to_string()),
    };
    let sph = &stage2.spherical_hypotheses[sph_idx];
    let center = sph.center;
    let radius = sph.radius;

    // Compute the plane normal for the great circle arc.
    // The plane passes through the sphere center and both vertex endpoints.
    // Normal = normalize(cross(v0 - center, v1 - center))
    let cv0 = [
        v0.point[0] - center[0],
        v0.point[1] - center[1],
        v0.point[2] - center[2],
    ];
    let cv1 = [
        v1.point[0] - center[0],
        v1.point[1] - center[1],
        v1.point[2] - center[2],
    ];
    let cross = [
        cv0[1] * cv1[2] - cv0[2] * cv1[1],
        cv0[2] * cv1[0] - cv0[0] * cv1[2],
        cv0[0] * cv1[1] - cv0[1] * cv1[0],
    ];
    let cross_len = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
    if cross_len < 1e-12 {
        // Vertices are collinear with center — degenerate, fall back to line
        return Err("sphere-cylinder tangent: vertices collinear with sphere center".to_string());
    }
    let normal = [cross[0] / cross_len, cross[1] / cross_len, cross[2] / cross_len];

    // Build the circle in 3D
    let center_pnt = gp::Pnt::new_real3(center[0], center[1], center[2]);
    let normal_dir = gp::Dir::new_real3(normal[0], normal[1], normal[2]);
    let ax2 = gp::Ax2::new_pnt_dir(&center_pnt, &normal_dir);
    let circle = geom::Circle::new_ax2_real(&ax2, radius);
    let circle_handle = geom::Circle::to_handle(circle).to_handle_curve();

    // Project vertex endpoints onto the circle curve to get parameters
    let t0 = project_point_on_curve(&v0.point, &circle_handle)?;
    let t1 = project_point_on_curve(&v1.point, &circle_handle)?;

    // For a circle (periodic curve with period 2π), we need to select the correct arc.
    // Sample mesh boundary vertices to determine which arc they lie on.
    let (t_lo, t_hi) = select_arc_parameters(
        t0,
        t1,
        &edge.mesh_boundary_vertices,
        &circle_handle,
        &stage2.mesh,
    );

    let trimmed = geom::TrimmedCurve::new_handlegeomcurve_real2(&circle_handle, t_lo, t_hi);
    let trimmed_handle = geom::TrimmedCurve::to_handle(trimmed).to_handle_curve();

    validate_tangent_curve(edge, &trimmed_handle, stage2, config)?;
    edge.curve_3d = Some(trimmed_handle);
    Ok(())
}

/// Trim parameters for a chain-repair candidate curve. Periodic candidates
/// (the circle) get the chain-supported arc, which may wrap past one period —
/// Geom_TrimmedCurve normalizes that against a periodic basis. Bounded
/// candidates (the B-spline fit) have exactly one arc between the projected
/// parameters, and Geom_TrimmedCurve panics with "parameters out of range"
/// on anything outside their domain, so the modulo-2π arc selection must not
/// run on them: order the parameters and clamp to the domain instead.
/// Returns None for a degenerate (zero-length) span, which Geom_TrimmedCurve
/// also refuses.
fn candidate_trim_parameters(
    t0: f64,
    t1: f64,
    mesh_boundary_vertices: &[usize],
    curve_handle: &geom::Handle_Geom_Curve,
    mesh: &crate::stage1::ConnectedMesh,
) -> Option<(f64, f64)> {
    let (t_lo, t_hi) = if curve_handle.get().is_periodic() {
        select_arc_parameters(t0, t1, mesh_boundary_vertices, curve_handle, mesh)
    } else {
        let c = curve_handle.get();
        let (lo, hi) = if t0 <= t1 { (t0, t1) } else { (t1, t0) };
        (lo.max(c.first_parameter()), hi.min(c.last_parameter()))
    };
    if t_hi - t_lo < 1e-9 {
        return None;
    }
    Some((t_lo, t_hi))
}

/// Select the correct arc on a periodic curve by sampling mesh boundary vertices.
/// Returns (t_lo, t_hi) such that the arc from t_lo to t_hi contains the boundary vertices.
fn select_arc_parameters(
    t0: f64,
    t1: f64,
    mesh_boundary_vertices: &[usize],
    curve_handle: &geom::Handle_Geom_Curve,
    mesh: &crate::stage1::ConnectedMesh,
) -> (f64, f64) {
    let period = std::f64::consts::TAU; // 2π for circle

    // Normalize t0 and t1 to [0, 2π)
    let t0n = ((t0 % period) + period) % period;
    let t1n = ((t1 % period) + period) % period;

    // The two possible arcs are t0n->t1n (forward) and t1n->t0n (reverse via wrap).
    // Sample mesh boundary vertices, project onto curve, and count which arc they support.
    let mut forward_count = 0;
    let mut reverse_count = 0;

    // Forward arc span: from t0n going forward to t1n
    let forward_span = if t1n > t0n {
        t1n - t0n
    } else {
        t1n + period - t0n
    };

    let step = mesh_boundary_vertices.len().max(1);
    let sample_count = mesh_boundary_vertices.len().min(20);
    let sample_step = step.checked_div(sample_count).unwrap_or(1);

    for (i, &vi) in mesh_boundary_vertices.iter().enumerate() {
        if sample_step > 1 && i % sample_step != 0 {
            continue;
        }
        let v = &mesh.vertices[vi];
        let pt = gp::Pnt::new_real3(v.x, v.y, v.z);
        let proj = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(&pt, curve_handle);
        if proj.nb_points() == 0 {
            continue;
        }
        let tp = proj.lower_distance_parameter();
        let tpn = ((tp % period) + period) % period;

        // Check if this point's parameter falls in the forward arc
        let delta = if tpn >= t0n {
            tpn - t0n
        } else {
            tpn + period - t0n
        };
        if delta <= forward_span {
            forward_count += 1;
        } else {
            reverse_count += 1;
        }
    }

    if forward_count >= reverse_count {
        // Forward arc: t0n to t0n + forward_span
        (t0n, t0n + forward_span)
    } else {
        // Reverse arc: t1n to t1n + (2π - forward_span)
        (t1n, t1n + (period - forward_span))
    }
}

/// Validate that mesh boundary vertices lie reasonably close to a tangent curve.
fn validate_tangent_curve(
    edge: &ReconEdge,
    curve_handle: &geom::Handle_Geom_Curve,
    stage2: &Stage2Output,
    config: &Config,
) -> Result<(), String> {
    let mut max_dist = 0.0_f64;
    for &vi in &edge.mesh_boundary_vertices {
        let v = &stage2.mesh.vertices[vi];
        let pt = gp::Pnt::new_real3(v.x, v.y, v.z);
        let proj = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(&pt, curve_handle);
        if proj.nb_points() > 0 {
            max_dist = max_dist.max(proj.lower_distance());
        }
    }
    // For models with imperfect surface assignments, mesh boundary vertices may be
    // slightly off the analytical tangent curve due to tessellation sagitta.
    // Use surface_tolerance as the threshold (typically 0.4mm) rather than
    // vertex_tolerance (1e-5mm) since the sagitta of curved surfaces can be significant.
    if max_dist > config.surface_tolerance_mm {
        return Err(format!(
            "tangent edge curve too far from mesh boundary (max dist {max_dist:.6}mm)"
        ));
    }
    Ok(())
}


/// Compute the outward surface normal at a point for viz camera orientation.
fn viz_surface_normal_at_point(
    fi: usize,
    point: [f64; 3],
    face_descriptors: &[FaceDescriptor],
    stage2: &Stage2Output,
) -> [f32; 3] {
    let ss_idx = face_descriptors[fi].selected_surface_idx;
    match &stage2.selected_surfaces[ss_idx] {
        SelectedSurface::Planar(i) => {
            let n = stage2.planar_hypotheses[*i].normal;
            [n[0] as f32, n[1] as f32, n[2] as f32]
        }
        SelectedSurface::Conical(i) => {
            let hyp = &stage2.conical_hypotheses[*i];
            let dp = [
                point[0] - hyp.apex[0],
                point[1] - hyp.apex[1],
                point[2] - hyp.apex[2],
            ];
            let h = dp[0] * hyp.axis_direction[0]
                + dp[1] * hyp.axis_direction[1]
                + dp[2] * hyp.axis_direction[2];
            let r = [dp[0] - h * hyp.axis_direction[0], dp[1] - h * hyp.axis_direction[1], dp[2] - h * hyp.axis_direction[2]];
            let rlen = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
            if rlen < 1e-10 {
                return [0.0_f32, 0.0, 1.0];
            }
            let cos_ha = hyp.half_angle.cos();
            let sin_ha = hyp.half_angle.sin();
            let ru = [r[0] / rlen, r[1] / rlen, r[2] / rlen];
            let sign: f64 = if hyp.convex { 1.0 } else { -1.0 };
            let n = [
                sign * (ru[0] * cos_ha - hyp.axis_direction[0] * sin_ha),
                sign * (ru[1] * cos_ha - hyp.axis_direction[1] * sin_ha),
                sign * (ru[2] * cos_ha - hyp.axis_direction[2] * sin_ha),
            ];
            [n[0] as f32, n[1] as f32, n[2] as f32]
        }
        SelectedSurface::Spherical(i) => {
            let hyp = &stage2.spherical_hypotheses[*i];
            let v = [
                point[0] - hyp.center[0],
                point[1] - hyp.center[1],
                point[2] - hyp.center[2],
            ];
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if len < 1e-10 {
                return [0.0_f32, 0.0, 1.0];
            }
            let sign: f64 = if hyp.convex { 1.0 } else { -1.0 };
            [(sign * v[0] / len) as f32, (sign * v[1] / len) as f32, (sign * v[2] / len) as f32]
        }
        SelectedSurface::Toroidal(i) => {
            let hyp = &stage2.toroidal_hypotheses[*i];
            let dp = [
                point[0] - hyp.center[0],
                point[1] - hyp.center[1],
                point[2] - hyp.center[2],
            ];
            let ax = hyp.axis_direction;
            let h = dp[0] * ax[0] + dp[1] * ax[1] + dp[2] * ax[2];
            let radial = [dp[0] - h * ax[0], dp[1] - h * ax[1], dp[2] - h * ax[2]];
            let radial_len = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            if radial_len < 1e-10 {
                return [0.0_f32, 0.0, 1.0];
            }
            let tube_center = [
                hyp.center[0] + hyp.major_radius * radial[0] / radial_len,
                hyp.center[1] + hyp.major_radius * radial[1] / radial_len,
                hyp.center[2] + hyp.major_radius * radial[2] / radial_len,
            ];
            let tv = [
                point[0] - tube_center[0],
                point[1] - tube_center[1],
                point[2] - tube_center[2],
            ];
            let tv_len = (tv[0] * tv[0] + tv[1] * tv[1] + tv[2] * tv[2]).sqrt();
            if tv_len < 1e-10 {
                return [0.0_f32, 0.0, 1.0];
            }
            let sign: f64 = if hyp.convex { 1.0 } else { -1.0 };
            [(sign * tv[0] / tv_len) as f32, (sign * tv[1] / tv_len) as f32, (sign * tv[2] / tv_len) as f32]
        }
        SelectedSurface::BSpline(i) => {
            match stage2.bspline_hypotheses[*i].normal_at(&point) {
                Some(n) => [n[0] as f32, n[1] as f32, n[2] as f32],
                None => [0.0_f32, 0.0, 1.0],
            }
        }
    }
}

/// Compute the vertex centroid of a selected surface for viz camera positioning.
fn viz_selected_surface_centroid(
    fi: usize,
    face_descriptors: &[FaceDescriptor],
    stage2: &Stage2Output,
) -> [f64; 3] {
    let ss_idx = face_descriptors[fi].selected_surface_idx;
    let vertex_indices: &[usize] = match &stage2.selected_surfaces[ss_idx] {
        SelectedSurface::Planar(i) => &stage2.planar_hypotheses[*i].vertices,
        SelectedSurface::Conical(i) => &stage2.conical_hypotheses[*i].vertices,
        SelectedSurface::Spherical(i) => &stage2.spherical_hypotheses[*i].vertices,
        SelectedSurface::Toroidal(i) => &stage2.toroidal_hypotheses[*i].vertices,
        SelectedSurface::BSpline(i) => &stage2.bspline_hypotheses[*i].vertices,
    };
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    for &vi in vertex_indices {
        let v = &stage2.mesh.vertices[vi];
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    let n = vertex_indices.len() as f64;
    if n > 0.0 {
        [cx / n, cy / n, cz / n]
    } else {
        [0.0, 0.0, 0.0]
    }
}

/// Compute edge curves for all ReconEdges.
fn compute_edge_curves_all(
    mut output: Stage3Output,
    config: &Config,
    viz: Option<&crate::viz::VizSender>,
) -> Result<Stage3Output, Stage3Error> {
    let t = Instant::now();
    let mut success_count = 0;
    let mut fail_count = 0;
    let mut chain_count = 0;
    let total = output.edges.len();

    for ei in 0..total {
        // Split borrow: extract edge mutably while borrowing face_descriptors immutably
        let (before, rest) = output.edges.split_at_mut(ei);
        let (edge, _after) = rest.split_first_mut().unwrap();
        let _ = before; // suppress unused warning

        if edge.tangent {
            match compute_tangent_edge_curve(
                edge,
                &output.face_descriptors,
                &output.vertices,
                &output.stage2,
                config,
            ) {
                Ok(()) => {
                    success_count += 1;
                    if config.verbose {
                        let curve = edge.curve_3d.as_ref().unwrap();
                        let c = curve.get();
                        let fp = c.first_parameter();
                        let lp = c.last_parameter();
                        let p_start = c.value(fp);
                        let p_end = c.value(lp);
                        eprintln!(
                            "  Edge {ei}: faces [{}, {}], tangent curve [{:.4},{:.4},{:.4}] -> [{:.4},{:.4},{:.4}], t=[{:.4}, {:.4}]",
                            edge.face_indices[0], edge.face_indices[1],
                            p_start.x(), p_start.y(), p_start.z(),
                            p_end.x(), p_end.y(), p_end.z(),
                            fp, lp,
                        );
                    }
                    if let Some(viz_sender) = viz {
                        let fi0 = edge.face_indices[0];
                        let fi1 = edge.face_indices[1];
                        let mut overlay = crate::viz::VizOverlay::new();
                        let face_colors = [
                            [0.3, 0.3, 0.9, 0.35],
                            [0.7, 0.3, 0.7, 0.35],
                        ];
                        for (k, &fi) in [fi0, fi1].iter().enumerate() {
                            let ss_idx = output.face_descriptors[fi].selected_surface_idx;
                            let ss = &output.stage2.selected_surfaces[ss_idx];
                            let mesh_faces: &[usize] = match ss {
                                SelectedSurface::Planar(i) => &output.stage2.planar_hypotheses[*i].faces,
                                SelectedSurface::Conical(i) => &output.stage2.conical_hypotheses[*i].faces,
                                SelectedSurface::Spherical(i) => &output.stage2.spherical_hypotheses[*i].faces,
                                SelectedSurface::Toroidal(i) => &output.stage2.toroidal_hypotheses[*i].faces,
                                SelectedSurface::BSpline(i) => &output.stage2.bspline_hypotheses[*i].faces,
                            };
                            overlay.face_highlights.push(crate::viz::FaceHighlight {
                                face_indices: mesh_faces.to_vec(),
                                color: face_colors[k],
                            });
                        }
                        let curve = edge.curve_3d.as_ref().unwrap();
                        overlay.lines.push(crate::viz::LineOverlay {
                            positions: crate::viz::sample_curve_for_viz(curve, 64),
                            color: [0.0, 1.0, 0.0, 1.0],
                            no_depth_test: true,
                        });
                        overlay.status_text = format!(
                            "Stage 3.3: Edge {ei}/{total}: faces [{fi0}, {fi1}] (tangent)"
                        );
                        let mid_pt = {
                            let c = curve.get();
                            let mid = (c.first_parameter() + c.last_parameter()) / 2.0;
                            let p = c.value(mid);
                            [p.x(), p.y(), p.z()]
                        };
                        overlay.focus_point = Some([mid_pt[0] as f32, mid_pt[1] as f32, mid_pt[2] as f32]);
                        overlay.focus_normal = Some(viz_surface_normal_at_point(
                            fi0, mid_pt, &output.face_descriptors, &output.stage2,
                        ));
                        viz_sender.show_and_wait(overlay);
                    }
                }
                Err(msg) => {
                    // Last resort, as for non-tangent edges: build the curve
                    // from the mesh boundary chain.
                    match edge_curve_from_boundary_chain(
                        edge,
                        &output.vertices,
                        &output.stage2.mesh,
                        config,
                    ) {
                        Ok(dev) if dev <= config.surface_tolerance_mm => {
                            success_count += 1;
                            chain_count += 1;
                            if !config.quiet {
                                eprintln!(
                                    "  Edge {ei}: tangent construction failed ({msg}); \
                                     using mesh-boundary-chain curve (chain dev {dev:.2e} mm)"
                                );
                            }
                        }
                        Ok(dev) => {
                            fail_count += 1;
                            if !config.quiet {
                                eprintln!(
                                    "  Edge {ei}: FAILED tangent edge - {msg}; best \
                                     chain-based curve still {dev:.3} mm from the mesh boundary"
                                );
                            }
                        }
                        Err(chain_msg) => {
                            fail_count += 1;
                            if config.verbose {
                                eprintln!(
                                    "  Edge {ei}: FAILED tangent edge - {msg} ({chain_msg})"
                                );
                            }
                        }
                    }
                }
            }
            continue;
        }

        match compute_edge_curve(
            edge,
            &output.face_descriptors,
            &output.vertices,
            &output.stage2.mesh,
            &output.stage2,
            config,
        ) {
            Ok(()) => {
                success_count += 1;
                if config.verbose {
                    let curve = edge.curve_3d.as_ref().unwrap();
                    let c = curve.get();
                    let fp = c.first_parameter();
                    let lp = c.last_parameter();
                    let p_start = c.value(fp);
                    let p_end = c.value(lp);
                    eprintln!(
                        "  Edge {ei}: faces [{}, {}], curve [{:.4},{:.4},{:.4}] -> [{:.4},{:.4},{:.4}], t=[{:.4}, {:.4}]",
                        edge.face_indices[0], edge.face_indices[1],
                        p_start.x(), p_start.y(), p_start.z(),
                        p_end.x(), p_end.y(), p_end.z(),
                        fp, lp,
                    );
                }
                    if let Some(viz_sender) = viz {
                        let fi0 = edge.face_indices[0];
                        let fi1 = edge.face_indices[1];
                        let mut overlay = crate::viz::VizOverlay::new();
                        let face_colors = [
                            [0.3, 0.3, 0.9, 0.35],
                            [0.7, 0.3, 0.7, 0.35],
                        ];
                        for (k, &fi) in [fi0, fi1].iter().enumerate() {
                            let ss_idx = output.face_descriptors[fi].selected_surface_idx;
                            let ss = &output.stage2.selected_surfaces[ss_idx];
                            let mesh_faces: &[usize] = match ss {
                                SelectedSurface::Planar(i) => &output.stage2.planar_hypotheses[*i].faces,
                                SelectedSurface::Conical(i) => &output.stage2.conical_hypotheses[*i].faces,
                                SelectedSurface::Spherical(i) => &output.stage2.spherical_hypotheses[*i].faces,
                                SelectedSurface::Toroidal(i) => &output.stage2.toroidal_hypotheses[*i].faces,
                                SelectedSurface::BSpline(i) => &output.stage2.bspline_hypotheses[*i].faces,
                            };
                            overlay.face_highlights.push(crate::viz::FaceHighlight {
                                face_indices: mesh_faces.to_vec(),
                                color: face_colors[k],
                            });
                        }
                        let curve = edge.curve_3d.as_ref().unwrap();
                        overlay.lines.push(crate::viz::LineOverlay {
                            positions: crate::viz::sample_curve_for_viz(curve, 64),
                            color: [0.0, 1.0, 0.0, 1.0],
                            no_depth_test: true,
                        });
                        overlay.status_text = format!(
                            "Stage 3.3: Edge {ei}/{total}: faces [{fi0}, {fi1}]"
                        );
                        let mid_pt = {
                            let c = curve.get();
                            let mid = (c.first_parameter() + c.last_parameter()) / 2.0;
                            let p = c.value(mid);
                            [p.x(), p.y(), p.z()]
                        };
                        overlay.focus_point = Some([mid_pt[0] as f32, mid_pt[1] as f32, mid_pt[2] as f32]);
                        overlay.focus_normal = Some(viz_surface_normal_at_point(
                            fi0, mid_pt, &output.face_descriptors, &output.stage2,
                        ));
                        viz_sender.show_and_wait(overlay);
                    }
            }
            Err(msg) => {
                let fi0 = edge.face_indices[0];
                let fi1 = edge.face_indices[1];
                let s0 = &output.stage2.selected_surfaces[output.face_descriptors[fi0].selected_surface_idx];
                let s1 = &output.stage2.selected_surfaces[output.face_descriptors[fi1].selected_surface_idx];
                let sname = |s: &SelectedSurface| match s {
                    SelectedSurface::Planar(_) => "plane",
                    SelectedSurface::Conical(i) => {
                        if output.stage2.conical_hypotheses[*i].is_cylinder() { "cylinder" } else { "cone" }
                    }
                    SelectedSurface::Spherical(_) => "sphere",
                    SelectedSurface::Toroidal(_) => "torus",
                    SelectedSurface::BSpline(_) => "bspline",
                };
                // Last resort: build the curve from the mesh boundary chain.
                match edge_curve_from_boundary_chain(
                    edge,
                    &output.vertices,
                    &output.stage2.mesh,
                    config,
                ) {
                    Ok(dev) if dev <= config.surface_tolerance_mm => {
                        success_count += 1;
                        chain_count += 1;
                        if !config.quiet {
                            eprintln!(
                                "  Edge {ei}: {}-{} construction failed ({msg}); \
                                 using mesh-boundary-chain curve (chain dev {dev:.2e} mm)",
                                sname(s0),
                                sname(s1),
                            );
                        }
                    }
                    Ok(dev) => {
                        fail_count += 1;
                        if !config.quiet {
                            eprintln!(
                                "  Edge {ei}: FAILED ({}-{}) - {msg}; best chain-based \
                                 curve still {dev:.3} mm from the mesh boundary",
                                sname(s0),
                                sname(s1),
                            );
                        }
                    }
                    Err(msg2) => {
                        fail_count += 1;
                        if !config.quiet {
                            eprintln!(
                                "  Edge {ei}: FAILED ({}-{}) - {msg}; {msg2}",
                                sname(s0),
                                sname(s1),
                            );
                        }
                    }
                }
            }
        }
    }

    if !config.quiet {
        eprintln!(
            "Stage 3.3 ({:.3}s): Computed {success_count}/{total} edge curves ({fail_count} failed)",
            t.elapsed().as_secs_f64(),
        );
        if chain_count > 0 {
            eprintln!(
                "  {chain_count} edge curve(s) built from mesh boundary chains \
                 (surface-surface construction failed)"
            );
        }
    }

    if fail_count > 0 && !config.quiet {
        eprintln!(
            "  Warning: {fail_count} edge curves could not be computed (affected faces will use fallback geometry)"
        );
    }

    Ok(output)
}

// ---------------------------------------------------------------------------
// Stage 3.4: Create OCCT faces from surfaces bounded by edge wires
// ---------------------------------------------------------------------------

/// Compute the centroid of a mesh face.
fn compute_mesh_face_centroid(face_idx: usize, mesh: &ConnectedMesh) -> [f64; 3] {
    let face = &mesh.faces[face_idx];
    let vc = face.vertex_count as usize;
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for i in 0..vc {
        let v = &mesh.vertices[face.vertex_indices[i]];
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    let n = vc as f64;
    [cx / n, cy / n, cz / n]
}

/// Group a face's edges into separate wire loops based on vertex connectivity.
///
/// Edges sharing BRepVertices are grouped together. Closed-loop edges
/// (both vertex_indices == usize::MAX) each form their own wire.
fn group_edges_into_wires(fd: &FaceDescriptor, edges: &[ReconEdge]) -> Vec<Vec<usize>> {
    if fd.edge_indices.is_empty() {
        return Vec::new();
    }

    let mut groups: Vec<Vec<usize>> = Vec::new();
    // Map from BRepVertex index to which group it belongs to
    let mut vertex_group: HashMap<usize, usize> = HashMap::new();

    for &ei in &fd.edge_indices {
        let edge = &edges[ei];
        let [v0, v1] = edge.vertex_indices;

        if v0 == usize::MAX && v1 == usize::MAX {
            // Closed loop — its own wire
            groups.push(vec![ei]);
            continue;
        }

        // Find existing group for v0 or v1
        let g0 = if v0 != usize::MAX {
            vertex_group.get(&v0).copied()
        } else {
            None
        };
        let g1 = if v1 != usize::MAX {
            vertex_group.get(&v1).copied()
        } else {
            None
        };

        match (g0, g1) {
            (None, None) => {
                let gi = groups.len();
                groups.push(vec![ei]);
                if v0 != usize::MAX {
                    vertex_group.insert(v0, gi);
                }
                if v1 != usize::MAX {
                    vertex_group.insert(v1, gi);
                }
            }
            (Some(gi), None) => {
                groups[gi].push(ei);
                if v1 != usize::MAX {
                    vertex_group.insert(v1, gi);
                }
            }
            (None, Some(gi)) => {
                groups[gi].push(ei);
                if v0 != usize::MAX {
                    vertex_group.insert(v0, gi);
                }
            }
            (Some(g0i), Some(g1i)) if g0i == g1i => {
                groups[g0i].push(ei);
            }
            (Some(g0i), Some(g1i)) => {
                // Merge two groups
                let (keep, merge) = if g0i < g1i {
                    (g0i, g1i)
                } else {
                    (g1i, g0i)
                };
                groups[keep].push(ei);
                let merged = std::mem::take(&mut groups[merge]);
                groups[keep].extend(merged);
                // Update vertex_group references
                for g in vertex_group.values_mut() {
                    if *g == merge {
                        *g = keep;
                    }
                }
            }
        }
    }

    // Filter empty groups (from merges)
    groups.into_iter().filter(|g| !g.is_empty()).collect()
}


fn compute_face_area_from_face(face: &topo_ds::Face) -> f64 {
    let mut gprops = g_prop::GProps::new();
    b_rep_g_prop::surface_properties_shape_gprops_bool2(
        face.as_shape(),
        &mut gprops,
        false, // UseTriangulation
        false, // SkipShared
    );
    gprops.mass()
}

/// Human-readable description of a BRepCheck_Status value.
fn brep_check_status_name(status: b_rep_check::Status) -> &'static str {
    use b_rep_check::Status::*;
    match status {
        Noerror => "no error",
        Invalidpointoncurve => "invalid point on curve",
        Invalidpointoncurveonsurface => "invalid point on curve-on-surface",
        Invalidpointonsurface => "invalid point on surface",
        No3dcurve => "no 3D curve",
        Multiple3dcurve => "multiple 3D curves",
        Invalid3dcurve => "invalid 3D curve",
        Nocurveonsurface => "no curve on surface (missing pcurve)",
        Invalidcurveonsurface => "invalid curve on surface",
        Invalidcurveonclosedsurface => "invalid curve on closed surface",
        Invalidsamerangeflag => "invalid SameRange flag",
        Invalidsameparameterflag => "invalid SameParameter flag",
        Invaliddegeneratedflag => "invalid degenerated flag",
        Freeedge => "free edge",
        Invalidmulticonnexity => "invalid multi-connexity",
        Invalidrange => "invalid range",
        Emptywire => "empty wire",
        Redundantedge => "redundant edge",
        Selfintersectingwire => "self-intersecting wire",
        Nosurface => "no surface",
        Invalidwire => "invalid wire",
        Redundantwire => "redundant wire",
        Intersectingwires => "intersecting wires",
        Invalidimbricationofwires => "invalid imbrication of wires",
        Emptyshell => "empty shell",
        Redundantface => "redundant face",
        Invalidimbricationofshells => "invalid imbrication of shells",
        Unorientableshape => "unorientable shape",
        Notclosed => "not closed",
        Notconnected => "not connected",
        Subshapenotinshape => "sub-shape not in shape",
        Badorientation => "bad orientation",
        Badorientationofsubshape => "bad orientation of sub-shape",
        Invalidpolygonontriangulation => "invalid polygon on triangulation",
        Invalidtolerancevalue => "invalid tolerance value",
        Enclosedregion => "enclosed region",
        Checkfail => "check failed",
    }
}

/// Validate created faces using BRepCheck_Analyzer, with detailed BRepCheck_Face diagnostics.
fn validate_faces(output: &Stage3Output, config: &Config) -> Result<(), Stage3Error> {
    let mut invalid_count = 0;

    for (fi, face) in output.faces.iter().enumerate() {
        let analyzer = b_rep_check::Analyzer::new_shape(face.as_shape());
        if !analyzer.is_valid() {
            // The detailed diagnostics below call further OCCT checkers that
            // can throw (e.g. Standard_NullObject on edges without geometry).
            // A diagnostics failure must not abort the pipeline — treat it as
            // an ordinary invalid face.
            let diag_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                validate_face_detailed(fi, face, &analyzer, output, config)
            }));
            match diag_result {
                Ok(counted) => {
                    if counted {
                        invalid_count += 1;
                    }
                }
                Err(_) => {
                    invalid_count += 1;
                    if !config.quiet {
                        eprintln!(
                            "  Warning: face {fi} failed BRepCheck: detailed diagnostics threw an OCCT exception"
                        );
                    }
                }
            }
        }
        if config.verbose {
            // Report face area and orientation for diagnostics
            let mut gprops = g_prop::GProps::new();
            b_rep_g_prop::surface_properties_shape_gprops_bool2(
                face.as_shape(),
                &mut gprops,
                false, // UseTriangulation
                false, // SkipShared
            );
            let area = gprops.mass();
            let orient = face.as_shape().orientation();
            eprintln!("  Face {fi}: area = {area:.4} mm², orientation = {orient:?}");
        }
    }

    if !config.quiet && invalid_count > 0 {
        eprintln!(
            "  Warning: {invalid_count}/{} faces failed BRepCheck validation",
            output.faces.len(),
        );
    }

    Ok(())
}

/// Detailed BRepCheck_Face diagnostics for a face that failed the analyzer.
/// Returns true (the face counts as invalid) and prints what failed.
fn validate_face_detailed(
    fi: usize,
    face: &topo_ds::Face,
    analyzer: &b_rep_check::Analyzer,
    output: &Stage3Output,
    config: &Config,
) -> bool {
    {
        {
            // Run detailed BRepCheck_Face diagnostics
            let mut checker = b_rep_check::Face::new_face(face);
            checker.minimum();
            let intersect = checker.intersect_wires(false);
            let classify = checker.classify_wires(false);
            let orient = checker.orientation_of_wires(false);
            let unorientable = checker.is_unorientable();

            // Check sub-shape edges
            let mut edge_issues = Vec::new();
            {
                let mut exp = top_exp::Explorer::new_shape_shapeenum2(
                    face.as_shape(),
                    top_abs::ShapeEnum::Edge,
                    top_abs::ShapeEnum::Shape,
                );
                let mut ei = 0;
                while exp.more() {
                    let sub = exp.value();
                    if !analyzer.is_valid_shape(sub) {
                        let edge = topo_ds::edge(sub);
                        let mut echk = b_rep_check::Edge::new_edge(edge);
                        echk.minimum();
                        let tol = echk.tolerance();
                        edge_issues.push(format!("edge {ei} (tol={tol:.2e})"));
                    }
                    ei += 1;
                    exp.next();
                }
            }

            // Check sub-shape vertices
            let mut vertex_issues = Vec::new();
            {
                let mut exp = top_exp::Explorer::new_shape_shapeenum2(
                    face.as_shape(),
                    top_abs::ShapeEnum::Vertex,
                    top_abs::ShapeEnum::Shape,
                );
                let mut vi = 0;
                while exp.more() {
                    let sub = exp.value();
                    if !analyzer.is_valid_shape(sub) {
                        vertex_issues.push(format!("vertex {vi}"));
                    }
                    vi += 1;
                    exp.next();
                }
            }

            {
                if !config.quiet {
                    let surface = &output.stage2.selected_surfaces
                        [output.face_descriptors[fi].selected_surface_idx];
                    let surface_desc = match surface {
                        SelectedSurface::Planar(_) => "planar",
                        SelectedSurface::Conical(i) => {
                            if output.stage2.conical_hypotheses[*i].is_cylinder() { "cylindrical" } else { "conical" }
                        }
                        SelectedSurface::Spherical(_) => "spherical",
                        SelectedSurface::Toroidal(_) => "toroidal",
                        SelectedSurface::BSpline(_) => "bspline",
                    };
                    let mut issues = Vec::new();
                    if intersect != b_rep_check::Status::Noerror {
                        issues.push(format!(
                            "intersecting wires ({})",
                            brep_check_status_name(intersect)
                        ));
                    }
                    if classify != b_rep_check::Status::Noerror {
                        issues.push(format!(
                            "wire classification ({})",
                            brep_check_status_name(classify)
                        ));
                    }
                    if orient != b_rep_check::Status::Noerror {
                        issues.push(format!(
                            "wire orientation ({})",
                            brep_check_status_name(orient)
                        ));
                    }
                    if unorientable {
                        issues.push("face is unorientable".to_string());
                    }
                    if !edge_issues.is_empty() {
                        issues.push(format!("bad edges: {}", edge_issues.join(", ")));
                    }
                    if !vertex_issues.is_empty() {
                        issues.push(format!("bad vertices: {}", vertex_issues.join(", ")));
                    }

                    // Check wire sub-shapes for additional diagnostics
                    {
                        let mut exp = top_exp::Explorer::new_shape_shapeenum2(
                            face.as_shape(),
                            top_abs::ShapeEnum::Wire,
                            top_abs::ShapeEnum::Shape,
                        );
                        let mut wi = 0;
                        while exp.more() {
                            let sub = exp.value();
                            if !analyzer.is_valid_shape(sub) {
                                issues.push(format!("wire {wi} invalid"));
                            }
                            wi += 1;
                            exp.next();
                        }
                    }
                    // Check the face shape itself
                    if !analyzer.is_valid_shape(face.as_shape()) {
                        issues.push("face shape invalid".to_string());
                    }

                    if issues.is_empty() {
                        issues.push("edge/vertex consistency issue".to_string());
                    }
                    eprintln!(
                        "  Warning: face {} ({}) failed BRepCheck: {}",
                        fi,
                        surface_desc,
                        issues.join("; ")
                    );
                    if config.verbose {
                        debug_dump_face(face);
                    }
                }
                true
            }
        }
    }
}

/// Dump a face's wire/edge structure for debugging: per edge occurrence its
/// orientation, flags, 3D range, pcurve endpoints on this face's surface,
/// and vertices with parameters.
fn debug_dump_face(face: &topo_ds::Face) {
    let surface = b_rep::Tool::surface_face(face);
    let identity_loc = top_loc::Location::new();
    let mut wire_exp = top_exp::Explorer::new_shape_shapeenum2(
        face.as_shape(),
        top_abs::ShapeEnum::Wire,
        top_abs::ShapeEnum::Shape,
    );
    let mut wi = 0;
    while wire_exp.more() {
        let wire = topo_ds::wire(wire_exp.value());
        eprintln!("    wire {wi} (orient {:?}):", wire.as_shape().orientation());
        let mut we = top_exp::Explorer::new_shape_shapeenum2(
            wire.as_shape(),
            top_abs::ShapeEnum::Edge,
            top_abs::ShapeEnum::Shape,
        );
        while we.more() {
            let edge = topo_ds::edge(we.value());
            let orient = edge.as_shape().orientation();
            let degen = b_rep::Tool::degenerated(edge);
            let tol = b_rep::Tool::tolerance_edge(edge);
            let mut f = 0.0;
            let mut l = 0.0;
            b_rep::Tool::range_edge_real2(edge, &mut f, &mut l);
            let mut pf = 0.0;
            let mut pl = 0.0;
            let c2d = b_rep::Tool::curve_on_surface_edge_handlegeomsurface_location_real2_boolptr(
                edge, &surface, &identity_loc, &mut pf, &mut pl, None,
            );
            let pdesc = if c2d.is_null() {
                "NO PCURVE".to_string()
            } else {
                let c = c2d.get();
                let a = c.value(pf);
                let b = c.value(pl);
                format!(
                    "pcurve ({:.4},{:.4})->({:.4},{:.4})",
                    a.x(),
                    a.y(),
                    b.x(),
                    b.y()
                )
            };
            eprintln!(
                "      edge {orient:?} degen={degen} range=[{f:.4},{l:.4}] tol={tol:.2e} {pdesc}"
            );
            let mut vexp = top_exp::Explorer::new_shape_shapeenum2(
                edge.as_shape(),
                top_abs::ShapeEnum::Vertex,
                top_abs::ShapeEnum::Shape,
            );
            while vexp.more() {
                let v = topo_ds::vertex(vexp.value());
                let p = b_rep::Tool::pnt(v);
                let t = b_rep::Tool::parameter_vertex_edge(v, edge);
                eprintln!(
                    "        vertex {:?} t={t:.4} at ({:.4},{:.4},{:.4}) tol={:.2e}",
                    v.as_shape().orientation(),
                    p.x(),
                    p.y(),
                    p.z(),
                    b_rep::Tool::tolerance_vertex(v),
                );
                vexp.next();
            }
            we.next();
        }
        wi += 1;
        wire_exp.next();
    }
}

fn compare_faces_to_step(
    output: &Stage3Output,
    config: &Config,
) -> Result<(), Stage3CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    // Count faces in reference STEP
    let mut step_face_count = 0;
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        compare_shape,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        step_face_count += 1;
        explorer.next();
    }

    let our_face_count = output.faces.len();
    if !config.quiet {
        eprintln!(
            "  Compare 3.4: {} faces created, {} faces in STEP reference",
            our_face_count, step_face_count,
        );
    }

    // For each face, sample a representative point from mesh face centroids
    // and check distance to the reference STEP shape.
    let mut max_dist = 0.0_f64;
    for (fi, fd) in output.face_descriptors.iter().enumerate() {
        let surface = &output.stage2.selected_surfaces[fd.selected_surface_idx];
        let faces = surface_faces(surface, &output.stage2);
        if faces.is_empty() {
            continue;
        }
        // Use the centroid of the first mesh face as a sample point
        let centroid = compute_mesh_face_centroid(faces[0], &output.stage2.mesh);
        let pt = gp::Pnt::new_real3(centroid[0], centroid[1], centroid[2]);
        let d_ref = min_distance_to_shape(&pt, compare_shape);
        max_dist = max_dist.max(d_ref);

        if d_ref > config.surface_tolerance_mm && config.verbose {
            eprintln!(
                "  Face {fi}: centroid distance to STEP: {:.6e} mm (tolerance: {:.6e})",
                d_ref, config.surface_tolerance_mm,
            );
        }
    }

    if !config.quiet {
        eprintln!(
            "  Compare 3.4: max face sample distance to STEP: {:.6e} mm",
            max_dist,
        );
    }

    if max_dist > config.surface_tolerance_mm {
        return Err(Stage3CompareError {
            substage: 4,
            check_type: "face",
            element_index: 0,
            max_distance: max_dist,
            tolerance: config.surface_tolerance_mm,
        });
    }

    Ok(())
}

/// Compare surface orientations of our deduced surfaces against reference STEP faces.
///
/// For each of our reconstructed faces, we:
/// 1. Sample a mesh face centroid as a representative point
/// 2. Compute our deduced surface normal at that point
/// 3. Find the closest STEP face by distance
/// 4. Evaluate the STEP surface normal at the nearest point (using D1 derivatives)
/// 5. Account for STEP face orientation (REVERSED flips normal)
/// 6. Compare the two normals — they should point in the same direction
fn compare_surface_orientations_to_step(
    output: &Stage3Output,
    config: &Config,
) {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    // Collect all STEP faces
    let mut step_faces: Vec<OwnedPtr<topo_ds::Face>> = Vec::new();
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        compare_shape,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        step_faces.push(topo_ds::face_shape(explorer.value()).to_owned());
        explorer.next();
    }

    let mut mismatch_count = 0;
    let mut checked_count = 0;

    for (fi, fd) in output.face_descriptors.iter().enumerate() {
        let surface = &output.stage2.selected_surfaces[fd.selected_surface_idx];
        let faces = surface_faces(surface, &output.stage2);
        if faces.is_empty() {
            continue;
        }

        // Use a centroid from one of the mesh faces as sample point
        let centroid = compute_mesh_face_centroid(faces[0], &output.stage2.mesh);
        let pt = gp::Pnt::new_real3(centroid[0], centroid[1], centroid[2]);

        // Compute our deduced surface normal at this point
        let our_normal = match surface_normal_at_point(surface, &output.stage2, &centroid) {
            Some(n) => n,
            None => continue, // degenerate point (on axis or center)
        };

        // Also get the mesh face normal for comparison
        let mesh_normal = output.stage2.mesh.faces[faces[0]].normal;

        // Find the closest STEP face
        let mut vertex_shape = b_rep_builder_api::MakeVertex::new_pnt(&pt);
        let mut best_face_idx = None;
        let mut best_dist = f64::MAX;
        for (si, step_face) in step_faces.iter().enumerate() {
            let progress = message::ProgressRange::new();
            let dist_calc = b_rep_extrema::DistShapeShape::new_shape2_extflag_extalgo_progressrange(
                vertex_shape.shape(),
                step_face.as_shape(),
                0,
                extrema::ExtAlgo::Grad,
                &progress,
            );
            if dist_calc.is_done() && dist_calc.nb_solution() > 0 {
                let d = dist_calc.value();
                if d < best_dist {
                    best_dist = d;
                    best_face_idx = Some(si);
                }
            }
        }

        let step_face_idx = match best_face_idx {
            Some(idx) if best_dist <= config.surface_tolerance_mm => idx,
            _ => continue, // no close STEP face found
        };

        let step_face = &step_faces[step_face_idx];

        // Get the underlying surface and project our point to get UV parameters
        let step_surface = b_rep::Tool::surface_face(step_face);
        let projector = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
            &pt,
            &step_surface,
            extrema::ExtAlgo::Grad,
        );

        if !projector.is_done() || projector.nb_points() == 0 {
            continue;
        }

        let mut u = 0.0_f64;
        let mut v = 0.0_f64;
        projector.parameters(1, &mut u, &mut v);

        // Evaluate surface derivatives at (u,v) using BRepAdaptor_Surface
        let adaptor = b_rep_adaptor::Surface::new_face(step_face);
        let mut p_eval = gp::Pnt::new_real3(0.0, 0.0, 0.0);
        let mut d1u = gp::Vec::new_real3(0.0, 0.0, 0.0);
        let mut d1v = gp::Vec::new_real3(0.0, 0.0, 0.0);
        adaptor.d1(u, v, &mut p_eval, &mut d1u, &mut d1v);

        // Normal = D1U × D1V
        let nx = d1u.y() * d1v.z() - d1u.z() * d1v.y();
        let ny = d1u.z() * d1v.x() - d1u.x() * d1v.z();
        let nz = d1u.x() * d1v.y() - d1u.y() * d1v.x();
        let nlen = (nx * nx + ny * ny + nz * nz).sqrt();

        if nlen < 1e-15 {
            continue; // degenerate derivatives
        }

        let mut step_normal = [nx / nlen, ny / nlen, nz / nlen];

        // Account for face orientation: REVERSED means the face normal is flipped
        if step_face.as_shape().orientation() == top_abs::Orientation::Reversed {
            step_normal = [-step_normal[0], -step_normal[1], -step_normal[2]];
        }

        checked_count += 1;

        // Compare our deduced normal with the STEP normal
        let dot_deduced = our_normal[0] * step_normal[0]
            + our_normal[1] * step_normal[1]
            + our_normal[2] * step_normal[2];

        // Compare mesh normal with STEP normal
        let dot_mesh = if let Some(mn) = mesh_normal {
            mn[0] * step_normal[0] + mn[1] * step_normal[1] + mn[2] * step_normal[2]
        } else {
            f64::NAN
        };

        let surface_type_name = match surface {
            SelectedSurface::Planar(_) => "planar",
            SelectedSurface::Conical(i) => {
                if output.stage2.conical_hypotheses[*i].is_cylinder() { "cylindrical" } else { "conical" }
            }
            SelectedSurface::Spherical(_) => "spherical",
            SelectedSurface::Toroidal(_) => "toroidal",
            SelectedSurface::BSpline(_) => "bspline",
        };

        let is_convex = match surface {
            SelectedSurface::Conical(idx) => {
                Some(output.stage2.conical_hypotheses[*idx].convex)
            }
            SelectedSurface::Spherical(idx) => {
                Some(output.stage2.spherical_hypotheses[*idx].convex)
            }
            SelectedSurface::Toroidal(idx) => {
                Some(output.stage2.toroidal_hypotheses[*idx].convex)
            }
            _ => None,
        };

        if dot_deduced < 0.0 {
            mismatch_count += 1;
            let convex_str = match is_convex {
                Some(true) => " convex",
                Some(false) => " concave",
                None => "",
            };
            let step_orient = if step_face.as_shape().orientation() == top_abs::Orientation::Reversed {
                "REVERSED"
            } else {
                "FORWARD"
            };
            eprintln!(
                "  Compare 3.4 orientation: face {fi} ({surface_type_name}{convex_str}): \
                 deduced normal DISAGREES with STEP (dot={dot_deduced:.4}, \
                 mesh_dot={dot_mesh:.4}, step_face={step_face_idx} {step_orient}, dist={best_dist:.2e})"
            );
            if config.verbose {
                eprintln!(
                    "    deduced normal: [{:.4}, {:.4}, {:.4}]",
                    our_normal[0], our_normal[1], our_normal[2]
                );
                eprintln!(
                    "    STEP normal:    [{:.4}, {:.4}, {:.4}]",
                    step_normal[0], step_normal[1], step_normal[2]
                );
                if let Some(mn) = mesh_normal {
                    eprintln!(
                        "    mesh normal:    [{:.4}, {:.4}, {:.4}]",
                        mn[0], mn[1], mn[2]
                    );
                }
            }
        } else if config.verbose {
            eprintln!(
                "  Compare 3.4 orientation: face {fi} ({surface_type_name}): \
                 OK (dot={dot_deduced:.4}, mesh_dot={dot_mesh:.4})"
            );
        }
    }

    if !config.quiet {
        eprintln!(
            "  Compare 3.4 orientation: {checked_count} faces checked, \
             {mismatch_count} orientation mismatch(es)"
        );
    }
}

// ===========================================================================
// Stage 3.4: shared-edge topological assembly (DEVELOPMENT_PLAN §3.x)
//
// Every ReconEdge is exactly one TopoDS_Edge, shared by both adjacent faces.
// For each (edge, face) incidence a pcurve is computed from analytic UV
// coordinates and attached to the shared edge. UV poles (sphere pole, cone
// apex) become degenerate edges and periodic wraps become seam edges —
// first-class topology built by one uniform mechanism for every surface
// type. Faces and shells are assembled directly with BRep_Builder from the
// shared edges, so shells are closed by construction and
// BRepBuilderAPI_Sewing is not used for construction (stage 3.5 validates).
// ===========================================================================

/// Parameterization chart of a face's surface. For analytic kinds it
/// mirrors the gp_Ax3 frame of the corresponding Geom surface, so UVs
/// computed here agree exactly with OCCT's parameterization of that surface.
/// For B-spline faces the chart wraps the exact fitted Geom_BSplineSurface
/// (via ShapeAnalysis_Surface projection) — UVs are OCCT's own.
struct SurfaceChart {
    origin: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
    z: [f64; 3],
    kind: ChartKind,
    /// Populated exactly when `kind == ChartKind::BSpline` (kept out of
    /// `ChartKind` so the enum stays `Copy` for the many by-value matches).
    bspline: Option<BSplineChartData>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum ChartKind {
    Plane,
    Cylinder { radius: f64 },
    /// OCCT convention: apex at origin (reference radius 0), V is the slant
    /// distance from the apex along a generator. The apex is the V=0 pole.
    Cone { half_angle: f64 },
    /// V is latitude; poles at V = ±π/2.
    Sphere { radius: f64 },
    /// U around the main axis, V around the tube; doubly periodic.
    Torus { major_radius: f64, minor_radius: f64 },
    /// Fitted freeform patch: no poles, and non-periodic in both directions
    /// EXCEPT the periodic fit of a closed band (chart data `u_periodic`),
    /// which is U-periodic with period 2π like the analytic charts. For
    /// non-periodic patches every seam/winding/degenerate-edge mechanism
    /// self-disables.
    BSpline,
}

/// Chart data for a B-spline face: the exact fitted surface object (shared
/// with `FaceDescriptor.surface`) and a cached ShapeAnalysis_Surface
/// projector. The projector caches internally (hence RefCell), and each
/// query is seeded from the previous solution (`NextValueOfUV`) — the
/// B-spline analog of the analytic charts' branch-continuity unwrap.
struct BSplineChartData {
    surf: OwnedPtr<geom::Handle_Geom_Surface>,
    sas: RefCell<OwnedPtr<shape_analysis::Surface>>,
    seed: std::cell::Cell<Option<(f64, f64)>>,
    /// Parameter bounds of the surface.
    bounds: (f64, f64, f64, f64),
    /// True for the periodic fit of a closed freeform band: U-periodic with
    /// knots rescaled to [0, 2π], so the chart behaves like the analytic
    /// periodic charts (principal U values, wrapped-face machinery).
    u_periodic: bool,
    /// Parameter-space anchor inside this face's own region. The two halves
    /// of a seam-split closed band (the fallback reconstruction) share a
    /// surface whose parameter-extreme boundary curves are the SAME 3D
    /// curve, so a projection landing on it is ambiguous between the two
    /// extremes; the branch nearer the face's own interior wins.
    home: std::cell::Cell<(f64, f64)>,
}

impl BSplineChartData {
    fn point_uv(&self, p: &[f64; 3]) -> (f64, f64) {
        use crate::stage2::bspline::{PROJ_GAP_RECHECK, PROJ_MAX_PRECISION, PROJ_PRECISION};
        let pnt = gp::Pnt::new_real3(p[0], p[1], p[2]);
        let mut sas = self.sas.borrow_mut();
        let mut best = match self.seed.get() {
            Some((su, sv)) => {
                let prev = gp::Pnt2d::new_real2(su, sv);
                let uv = sas.next_value_of_uv(&prev, &pnt, PROJ_PRECISION, PROJ_MAX_PRECISION);
                (uv.x(), uv.y(), sas.gap())
            }
            None => {
                let uv = sas.value_of_uv(&pnt, PROJ_PRECISION);
                (uv.x(), uv.y(), sas.gap())
            }
        };
        // Re-verify suspicious seeded solutions with a fresh global search:
        // NextValueOfUV's Newton step can settle in a nearby local minimum
        // without triggering its internal fallback.
        if self.seed.get().is_some() && best.2 > PROJ_GAP_RECHECK {
            let uv = sas.value_of_uv(&pnt, PROJ_PRECISION);
            let gap = sas.gap();
            if gap < best.2 {
                best = (uv.x(), uv.y(), gap);
            }
        }
        drop(sas);
        let (u, v) = if self.u_periodic {
            // Match the analytic periodic charts' contract: principal U in
            // [0, 2π). ShapeAnalysis returns principal values on a periodic
            // surface already; the reduction just guards the seam itself.
            // Branch continuity along edges is the caller's unwrap, exactly
            // as for cylinders/cones.
            (best.0.rem_euclid(2.0 * std::f64::consts::PI), best.1)
        } else {
            self.resolve_boundary_branch(best.0, best.1)
        };
        self.seed.set(Some((u, v)));
        (u, v)
    }

    /// A result on a parameter boundary whose opposite boundary carries the
    /// same 3D curve (a seam-split closed band, the fallback) is ambiguous:
    /// keep the branch nearer this face's `home` anchor. Coincidence is
    /// verified in 3D, so ordinary (non-closed) patch boundaries are never
    /// touched.
    fn resolve_boundary_branch(&self, u: f64, v: f64) -> (f64, f64) {
        let (u0, u1, v0, v1) = self.bounds;
        let home = self.home.get();
        let mut out = (u, v);
        let eps_u = 1e-3 * (u1 - u0);
        let eps_v = 1e-3 * (v1 - v0);
        if (v - v0).abs() < eps_v || (v1 - v).abs() < eps_v {
            let here = if (v - v0).abs() < eps_v { v0 } else { v1 };
            let there = if here == v0 { v1 } else { v0 };
            let a = self.point_at(u, here);
            let b = self.point_at(u, there);
            let d2 = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);
            if d2 < 1e-12 && (there - home.1).abs() < (here - home.1).abs() {
                out.1 = there;
            }
        }
        if (u - u0).abs() < eps_u || (u1 - u).abs() < eps_u {
            let here = if (u - u0).abs() < eps_u { u0 } else { u1 };
            let there = if here == u0 { u1 } else { u0 };
            let a = self.point_at(here, out.1);
            let b = self.point_at(there, out.1);
            let d2 = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);
            if d2 < 1e-12 && (there - home.0).abs() < (here - home.0).abs() {
                out.0 = there;
            }
        }
        out
    }

    fn point_at(&self, u: f64, v: f64) -> [f64; 3] {
        let mut pnt = gp::Pnt::new_real3(0.0, 0.0, 0.0);
        self.surf.get().d0(u, v, &mut pnt);
        [pnt.x(), pnt.y(), pnt.z()]
    }
}

const TWO_PI: f64 = 2.0 * std::f64::consts::PI;

fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl SurfaceChart {
    fn u_periodic(&self) -> bool {
        match self.kind {
            ChartKind::Plane => false,
            ChartKind::BSpline => self.bspline.as_ref().is_some_and(|d| d.u_periodic),
            _ => true,
        }
    }

    fn v_periodic(&self) -> bool {
        matches!(self.kind, ChartKind::Torus { .. })
    }

    /// Analytic (u, v) of a 3D point on (or near) the surface, plus the
    /// point's distance from the chart's U-degeneracy locus (pole): U is
    /// meaningless when that distance is ~0. u is reduced to [0, 2π) for
    /// periodic charts; v is reduced to [0, 2π) for the torus.
    fn point_uv(&self, p: &[f64; 3]) -> (f64, f64, f64) {
        if let Some(bs) = &self.bspline {
            // Fitted patches have no UV degeneracies: pole radius is ∞,
            // like the plane chart.
            let (u, v) = bs.point_uv(p);
            return (u, v, f64::INFINITY);
        }
        let d = [p[0] - self.origin[0], p[1] - self.origin[1], p[2] - self.origin[2]];
        let dx = dot3(&d, &self.x);
        let dy = dot3(&d, &self.y);
        let dz = dot3(&d, &self.z);
        match self.kind {
            ChartKind::Plane => (dx, dy, f64::INFINITY),
            ChartKind::Cylinder { radius: _ } => {
                let rad = (dx * dx + dy * dy).sqrt();
                (dy.atan2(dx).rem_euclid(TWO_PI), dz, rad)
            }
            ChartKind::Cone { half_angle } => {
                let rad = (dx * dx + dy * dy).sqrt();
                (dy.atan2(dx).rem_euclid(TWO_PI), dz / half_angle.cos(), rad)
            }
            ChartKind::Sphere { radius } => {
                let rad = (dx * dx + dy * dy).sqrt();
                let v = (dz / radius).clamp(-1.0, 1.0).asin();
                (dy.atan2(dx).rem_euclid(TWO_PI), v, rad)
            }
            ChartKind::Torus { major_radius, minor_radius: _ } => {
                let rad = (dx * dx + dy * dy).sqrt();
                let v = dz.atan2(rad - major_radius).rem_euclid(TWO_PI);
                (dy.atan2(dx).rem_euclid(TWO_PI), v, rad)
            }
            ChartKind::BSpline => unreachable!("handled by the early return above"),
        }
    }

    /// Evaluate the surface at (u, v). Used for pcurve deviation checks.
    fn point_at(&self, u: f64, v: f64) -> [f64; 3] {
        if let Some(bs) = &self.bspline {
            return bs.point_at(u, v);
        }
        let (cu, su) = (u.cos(), u.sin());
        let (rad, h) = match self.kind {
            ChartKind::Plane => {
                return [
                    self.origin[0] + u * self.x[0] + v * self.y[0],
                    self.origin[1] + u * self.x[1] + v * self.y[1],
                    self.origin[2] + u * self.x[2] + v * self.y[2],
                ];
            }
            ChartKind::Cylinder { radius } => (radius, v),
            ChartKind::Cone { half_angle } => (v * half_angle.sin(), v * half_angle.cos()),
            ChartKind::Sphere { radius } => (radius * v.cos(), radius * v.sin()),
            ChartKind::Torus { major_radius, minor_radius } => {
                (major_radius + minor_radius * v.cos(), minor_radius * v.sin())
            }
            ChartKind::BSpline => unreachable!("handled by the early return above"),
        };
        [
            self.origin[0] + rad * (cu * self.x[0] + su * self.y[0]) + h * self.z[0],
            self.origin[1] + rad * (cu * self.x[1] + su * self.y[1]) + h * self.z[1],
            self.origin[2] + rad * (cu * self.x[2] + su * self.y[2]) + h * self.z[2],
        ]
    }

    /// The V parameter and 3D position of this chart's pole on the given
    /// side (sign of V), if the surface has one there. Cones have a single
    /// pole (the apex, V=0); spheres have poles at V=±π/2.
    fn pole(&self, positive_v_side: bool) -> Option<(f64, [f64; 3])> {
        match self.kind {
            ChartKind::Cone { .. } => {
                // The apex (V=0) bounds the face from below only: our cone
                // axes point from apex toward the mesh, so faces live at V>0.
                if positive_v_side {
                    None
                } else {
                    Some((0.0, self.origin))
                }
            }
            ChartKind::Sphere { radius } => {
                let sign = if positive_v_side { 1.0 } else { -1.0 };
                let v = sign * std::f64::consts::FRAC_PI_2;
                Some((
                    v,
                    [
                        self.origin[0] + sign * radius * self.z[0],
                        self.origin[1] + sign * radius * self.z[1],
                        self.origin[2] + sign * radius * self.z[2],
                    ],
                ))
            }
            _ => None,
        }
    }
}

/// Build the analytic chart for a face, together with a replacement OCCT
/// surface when the chart requires one (spheres are re-axised).
///
/// The sphere axis is a pure gauge choice — nothing in the B-rep depends on
/// it — so we point it at the mesh region's mean direction from the center.
/// This guarantees the region's boundary stays away from the UV poles: the
/// poles land either strictly inside the region (a polar cap, handled by a
/// degenerate pole edge) or strictly outside it. That one choice replaces
/// all of the old pole-proximity special cases.
fn build_face_chart(
    fd: &FaceDescriptor,
    stage2: &Stage2Output,
) -> (SurfaceChart, Option<OwnedPtr<geom::Handle_Geom_Surface>>) {
    let surface = &stage2.selected_surfaces[fd.selected_surface_idx];
    let ax3_frame = |origin: [f64; 3], dir: [f64; 3]| -> ([f64; 3], [f64; 3], [f64; 3]) {
        let o = gp::Pnt::new_real3(origin[0], origin[1], origin[2]);
        let d = gp::Dir::new_real3(dir[0], dir[1], dir[2]);
        let ax3 = gp::Ax3::new_pnt_dir(&o, &d);
        let x = ax3.x_direction();
        let y = ax3.y_direction();
        let z = ax3.direction();
        (
            [x.x(), x.y(), x.z()],
            [y.x(), y.y(), y.z()],
            [z.x(), z.y(), z.z()],
        )
    };
    match surface {
        SelectedSurface::Planar(idx) => {
            let hyp = &stage2.planar_hypotheses[*idx];
            let origin = [
                hyp.normal[0] * hyp.distance,
                hyp.normal[1] * hyp.distance,
                hyp.normal[2] * hyp.distance,
            ];
            let (x, y, z) = ax3_frame(origin, hyp.normal);
            (SurfaceChart { origin, x, y, z, kind: ChartKind::Plane, bspline: None }, None)
        }
        SelectedSurface::Conical(idx) => {
            let hyp = &stage2.conical_hypotheses[*idx];
            let (x, y, z) = ax3_frame(hyp.apex, hyp.axis_direction);
            let kind = if hyp.is_cylinder() {
                ChartKind::Cylinder { radius: hyp.radius_at_apex }
            } else {
                ChartKind::Cone { half_angle: hyp.half_angle }
            };
            (SurfaceChart { origin: hyp.apex, x, y, z, kind, bspline: None }, None)
        }
        SelectedSurface::Spherical(idx) => {
            let hyp = &stage2.spherical_hypotheses[*idx];
            // Mean direction of the region's mesh face centroids from the
            // sphere center (deterministic: hyp.faces order is stable).
            let mut mean = [0.0_f64; 3];
            for &mfi in &hyp.faces {
                let c = compute_mesh_face_centroid(mfi, &stage2.mesh);
                let d = [c[0] - hyp.center[0], c[1] - hyp.center[1], c[2] - hyp.center[2]];
                let len = dot3(&d, &d).sqrt();
                if len > 1e-12 {
                    mean[0] += d[0] / len;
                    mean[1] += d[1] / len;
                    mean[2] += d[2] / len;
                }
            }
            let len = dot3(&mean, &mean).sqrt();
            let axis = if len > 1e-9 {
                [mean[0] / len, mean[1] / len, mean[2] / len]
            } else {
                [0.0, 0.0, 1.0] // full sphere: any axis works
            };
            let (x, y, z) = ax3_frame(hyp.center, axis);
            let origin = gp::Pnt::new_real3(hyp.center[0], hyp.center[1], hyp.center[2]);
            let dir = gp::Dir::new_real3(axis[0], axis[1], axis[2]);
            let ax3 = gp::Ax3::new_pnt_dir(&origin, &dir);
            let sphere = geom::SphericalSurface::new_ax3_real(&ax3, hyp.radius);
            let handle = geom::SphericalSurface::to_handle(sphere).to_handle_surface();
            (
                SurfaceChart {
                    origin: hyp.center,
                    x,
                    y,
                    z,
                    kind: ChartKind::Sphere { radius: hyp.radius },
                    bspline: None,
                },
                Some(handle),
            )
        }
        SelectedSurface::Toroidal(idx) => {
            let hyp = &stage2.toroidal_hypotheses[*idx];
            let (x, y, z) = ax3_frame(hyp.center, hyp.axis_direction);
            (
                SurfaceChart {
                    origin: hyp.center,
                    x,
                    y,
                    z,
                    kind: ChartKind::Torus {
                        major_radius: hyp.major_radius,
                        minor_radius: hyp.minor_radius,
                    },
                    bspline: None,
                },
                None,
            )
        }
        SelectedSurface::BSpline(idx) => {
            // The chart must wrap the exact fitted Geom_BSplineSurface
            // stored in the hypothesis (and shared by fd.surface) — never a
            // reconstruction. The pole radius is ∞. For open patches both
            // periodic flags are false, so the seam/winding/degenerate-pole
            // machinery self-disables and they take the plain-loop path;
            // the periodic fit of a closed band reports U-periodic (period
            // 2π after the fit-time knot rescale) and takes the same
            // wrapped-face path as cylinders.
            let hyp = &stage2.bspline_hypotheses[*idx];
            let surf = hyp.surface.to_handle_surface();
            let sas = shape_analysis::Surface::new_handlegeomsurface(&surf);
            let (mut u0, mut u1, mut v0, mut v1) = (0.0, 0.0, 0.0, 0.0);
            surf.get().bounds(&mut u0, &mut u1, &mut v0, &mut v1);
            let data = BSplineChartData {
                surf,
                sas: RefCell::new(sas),
                seed: std::cell::Cell::new(None),
                bounds: (u0, u1, v0, v1),
                u_periodic: hyp.u_periodic,
                home: std::cell::Cell::new((0.5 * (u0 + u1), 0.5 * (v0 + v1))),
            };
            // Anchor the chart inside this face's own region: project the
            // mean of the region's face centroids (for a curved region it
            // lies off-surface, and its nearest surface point is safely
            // mid-region). This both seeds the first query and defines the
            // `home` used to resolve coincident-boundary branches on
            // seam-split closed bands.
            let mut mean = [0.0_f64; 3];
            for &mfi in &hyp.faces {
                let c = compute_mesh_face_centroid(mfi, &stage2.mesh);
                mean = [mean[0] + c[0], mean[1] + c[1], mean[2] + c[2]];
            }
            let n = hyp.faces.len().max(1) as f64;
            let home = data.point_uv(&[mean[0] / n, mean[1] / n, mean[2] / n]);
            data.home.set(home);
            data.seed.set(Some(home));
            (
                SurfaceChart {
                    origin: [0.0; 3],
                    x: [1.0, 0.0, 0.0],
                    y: [0.0, 1.0, 0.0],
                    z: [0.0, 0.0, 1.0],
                    kind: ChartKind::BSpline,
                    bspline: Some(data),
                },
                None,
            )
        }
    }
}

/// UV samples of one edge on one face's chart: ascending edge parameter,
/// U (and V for tori) branch-continuous within the edge, and pole-degenerate
/// U values backfilled from the nearest well-defined sample.
struct EdgeUv {
    ts: Vec<f64>,
    us: Vec<f64>,
    vs: Vec<f64>,
    /// Whether the first/last sample sits at a U pole (U meaningless there).
    start_at_pole: bool,
    end_at_pole: bool,
}

impl EdgeUv {
    fn first_uv(&self) -> (f64, f64) {
        (self.us[0], self.vs[0])
    }
    fn last_uv(&self) -> (f64, f64) {
        (*self.us.last().unwrap(), *self.vs.last().unwrap())
    }
}

/// Sample an edge's stored 3D curve over its stored parameter range and map
/// the points through the chart. Returns None for edges without a 3D curve.
fn sample_edge_uv(
    edge: &topo_ds::Edge,
    chart: &SurfaceChart,
    n: usize,
    pole_eps: f64,
) -> Option<EdgeUv> {
    let mut first = 0.0;
    let mut last = 0.0;
    let curve = b_rep::Tool::curve_edge_real2(edge, &mut first, &mut last);
    if curve.is_null() {
        return None;
    }
    let c = curve.get();
    let mut ts = Vec::with_capacity(n);
    let mut raw: Vec<(f64, f64, f64)> = Vec::with_capacity(n); // (u, v, radial)
    for i in 0..n {
        let t = first + (last - first) * (i as f64) / ((n - 1) as f64);
        let pt = c.value(t);
        let (u, v, rad) = chart.point_uv(&[pt.x(), pt.y(), pt.z()]);
        if std::env::var_os("BREPPER_DEBUG_EDGE_UV").is_some() && chart.bspline.is_some() {
            eprintln!(
                "      edge-uv sample {i}: t={t:.4} p=[{:.3},{:.3},{:.3}] -> uv=({u:.4},{v:.4})",
                pt.x(), pt.y(), pt.z(),
            );
        }
        ts.push(t);
        raw.push((u, v, rad));
    }
    // Backfill pole-degenerate U values and unwrap U (and torus V) so the
    // samples are branch-continuous along the edge. Unwrapping applies only
    // to periodic directions — on a plane chart U/V are Cartesian
    // coordinates and reducing their steps mod 2π corrupts them.
    let good = |i: usize| raw[i].2 >= pole_eps;
    let first_good = (0..raw.len()).find(|&i| good(i))?;
    let mut us = vec![0.0; raw.len()];
    let mut vs = vec![0.0; raw.len()];
    let mut prev_u = raw[first_good].0;
    for i in 0..raw.len() {
        let u = if good(i) {
            let mut u = raw[i].0;
            if chart.u_periodic() {
                u -= TWO_PI * ((u - prev_u) / TWO_PI).round();
            }
            prev_u = u;
            u
        } else {
            prev_u
        };
        us[i] = u;
    }
    // Leading pole samples take the first good sample's (unwrapped) U.
    for i in 0..first_good {
        us[i] = us[first_good];
    }
    let mut prev_v = raw[0].1;
    for i in 0..raw.len() {
        let mut v = raw[i].1;
        if chart.v_periodic() {
            v -= TWO_PI * ((v - prev_v) / TWO_PI).round();
        }
        prev_v = v;
        vs[i] = v;
    }
    let start_at_pole = !good(0);
    let end_at_pole = !good(raw.len() - 1);
    Some(EdgeUv { ts, us, vs, start_at_pole, end_at_pole })
}

/// One edge occurrence in a face boundary loop, in traversal order.
#[derive(Clone, Copy, Debug)]
struct LoopItem {
    /// ReconEdge index.
    ei: usize,
    /// True when the loop traverses the edge in increasing parameter.
    forward: bool,
}

/// A face boundary loop with branch-aligned UV data.
struct BoundaryLoop {
    items: Vec<LoopItem>,
    /// Per item: (u_offset, v_offset) — multiples of 2π added to the edge's
    /// canonical samples to make the loop branch-continuous.
    offsets: Vec<(f64, f64)>,
    /// Net winding around the periodic directions (±1 wraps, 0 plain).
    winding_u: i32,
    winding_v: i32,
    /// U jumps across poles between item j and item j+1 (cyclic): the loop
    /// passes through a pole vertex there and needs a degenerate edge.
    pole_gaps: Vec<PoleGap>,
}

#[derive(Clone, Copy, Debug)]
struct PoleGap {
    /// The gap follows this item index (cyclic).
    after_item: usize,
    u_from: f64,
    u_to: f64,
    /// True when the pole is on the positive-V side of the chart.
    positive_v_side: bool,
}

/// Order the edges of one vertex-connectivity group into a traversal cycle.
fn order_loop(
    group: &[usize],
    edge_vertices: &[[usize; 2]],
) -> Result<Vec<LoopItem>, Stage3Error> {
    if group.len() == 1 {
        return Ok(vec![LoopItem { ei: group[0], forward: true }]);
    }
    let mut items = Vec::with_capacity(group.len());
    let mut used = vec![false; group.len()];
    items.push(LoopItem { ei: group[0], forward: true });
    used[0] = true;
    let start_v = edge_vertices[group[0]][0];
    let mut cur_v = edge_vertices[group[0]][1];
    while items.len() < group.len() {
        let mut advanced = false;
        for (gi, &ei) in group.iter().enumerate() {
            if used[gi] {
                continue;
            }
            let [v0, v1] = edge_vertices[ei];
            if v0 == cur_v {
                items.push(LoopItem { ei, forward: true });
                cur_v = v1;
            } else if v1 == cur_v {
                items.push(LoopItem { ei, forward: false });
                cur_v = v0;
            } else {
                continue;
            }
            used[gi] = true;
            advanced = true;
            break;
        }
        if !advanced {
            return Err(Stage3Error::AdjacencyError(format!(
                "boundary loop is not vertex-connected (edges {group:?})"
            )));
        }
    }
    if cur_v != start_v {
        return Err(Stage3Error::AdjacencyError(format!(
            "boundary loop does not close (edges {group:?})"
        )));
    }
    Ok(items)
}

/// Branch-align a walked loop: choose per-item 2π offsets so consecutive
/// item endpoints agree in UV, detect pole gaps, and compute the winding.
fn align_loop(
    items: Vec<LoopItem>,
    samples: &HashMap<usize, EdgeUv>,
    chart: &SurfaceChart,
) -> Result<BoundaryLoop, Stage3Error> {
    let n = items.len();
    let mut offsets = vec![(0.0, 0.0); n];
    let mut pole_gaps: Vec<PoleGap> = Vec::new();

    // Traversal-order endpoints of item i with its offset applied.
    let ends = |item: &LoopItem, off: (f64, f64)| -> ((f64, f64), (f64, f64), bool, bool) {
        let s = &samples[&item.ei];
        let (a, b) = (s.first_uv(), s.last_uv());
        let (start, end, sp, ep) = if item.forward {
            (a, b, s.start_at_pole, s.end_at_pole)
        } else {
            (b, a, s.end_at_pole, s.start_at_pole)
        };
        (
            (start.0 + off.0, start.1 + off.1),
            (end.0 + off.0, end.1 + off.1),
            sp,
            ep,
        )
    };

    for i in 1..n {
        let (_, prev_end, _, prev_end_pole) = ends(&items[i - 1], offsets[i - 1]);
        let (cur_start, _, cur_start_pole, _) = ends(&items[i], (0.0, 0.0));
        let mut du = 0.0;
        let mut dv = 0.0;
        if chart.u_periodic() {
            du = TWO_PI * ((prev_end.0 - cur_start.0) / TWO_PI).round();
        }
        if chart.v_periodic() {
            dv = TWO_PI * ((prev_end.1 - cur_start.1) / TWO_PI).round();
        }
        offsets[i] = (du, dv);
        if prev_end_pole && cur_start_pole {
            // The junction vertex is a pole: U legitimately jumps there.
            let (cs, _, _, _) = ends(&items[i], offsets[i]);
            if (cs.0 - prev_end.0).abs() > 1e-6 {
                pole_gaps.push(PoleGap {
                    after_item: i - 1,
                    u_from: prev_end.0,
                    u_to: cs.0,
                    positive_v_side: prev_end.1 > 0.0,
                });
            }
        }
    }

    // Winding: net change over the full cycle, including the closing junction.
    let (first_start, _, first_start_pole, _) = ends(&items[0], offsets[0]);
    let (_, last_end, _, last_end_pole) = ends(&items[n - 1], offsets[n - 1]);
    let mut wu = 0i32;
    let mut wv = 0i32;
    if chart.u_periodic() {
        if last_end_pole && first_start_pole {
            if (last_end.0 - first_start.0).abs() > 1e-6 {
                pole_gaps.push(PoleGap {
                    after_item: n - 1,
                    u_from: last_end.0,
                    u_to: first_start.0,
                    positive_v_side: last_end.1 > 0.0,
                });
            }
        } else {
            wu = ((last_end.0 - first_start.0) / TWO_PI).round() as i32;
        }
    }
    if chart.v_periodic() {
        wv = ((last_end.1 - first_start.1) / TWO_PI).round() as i32;
    }
    Ok(BoundaryLoop { items, offsets, winding_u: wu, winding_v: wv, pole_gaps })
}

/// Which periodic direction (if any) a face's boundary wraps all the way
/// around, requiring a seam edge.
#[derive(Clone, Copy, Debug, PartialEq)]
enum WrapKind {
    None,
    U,
    V,
}

/// How the seam endpoint on a winding boundary loop is realized.
#[derive(Clone, Copy, Debug)]
enum AnchorSource {
    /// An existing shared BRepVertex (index into Stage3Output.vertices).
    Vertex(usize),
    /// A closed-loop edge re-ranged to start at the seam; the anchor vertex
    /// lives in AssemblyCtx.closed_anchors under this edge id.
    ClosedEdge(usize),
}

/// Per-face seam plan: one anchor per winding boundary loop, in the face's
/// deterministic loop-scan order.
struct WrapPlan {
    anchors: Vec<AnchorSource>,
}

/// Shared assembly state for stage 3.4.
struct AssemblyCtx {
    topo_vertices: Vec<OwnedPtr<topo_ds::Vertex>>,
    topo_edges: Vec<OwnedPtr<topo_ds::Edge>>,
    /// BRepVertex ids [start, end] of each edge in stored-curve direction
    /// (both usize::MAX for closed-loop edges).
    edge_vertices: Vec<[usize; 2]>,
    /// Seam-anchor vertices added to re-ranged closed-loop edges.
    closed_anchors: HashMap<usize, OwnedPtr<topo_ds::Vertex>>,
}

/// Pole-detection radius: below this distance from the chart's degeneracy
/// locus, a sample's U coordinate is treated as meaningless.
fn pole_radius_eps(config: &Config) -> f64 {
    (10.0 * config.vertex_tolerance_mm).max(1e-6)
}

/// Number of samples per edge for UV analysis and pcurve interpolation.
const PCURVE_BASE_SAMPLES: usize = 33;

/// Scan of one face's boundary in UV: chart, per-edge samples, aligned
/// loops, and the wrap classification.
struct FaceScan {
    chart: SurfaceChart,
    surface_override: Option<OwnedPtr<geom::Handle_Geom_Surface>>,
    samples: HashMap<usize, EdgeUv>,
    loops: Vec<BoundaryLoop>,
    wrap: WrapKind,
}

fn scan_face(
    fi: usize,
    output: &Stage3Output,
    ctx: &AssemblyCtx,
    config: &Config,
) -> Result<FaceScan, Stage3Error> {
    let fd = &output.face_descriptors[fi];
    let (chart, surface_override) = build_face_chart(fd, &output.stage2);
    let pole_eps = pole_radius_eps(config);

    let mut samples = HashMap::new();
    for &ei in &fd.edge_indices {
        let uv = sample_edge_uv(&ctx.topo_edges[ei], &chart, PCURVE_BASE_SAMPLES, pole_eps)
            .ok_or_else(|| {
                Stage3Error::AdjacencyError(format!(
                    "face {fi}: edge {ei} has no stored 3D curve to project"
                ))
            })?;
        samples.insert(ei, uv);
    }

    let groups = group_edges_into_wires(fd, &output.edges);
    let mut loops = Vec::with_capacity(groups.len());
    for group in &groups {
        let items = order_loop(group, &ctx.edge_vertices)?;
        loops.push(align_loop(items, &samples, &chart)?);
    }

    let mut wrap = WrapKind::None;
    for lp in &loops {
        if lp.winding_u != 0 && lp.winding_v != 0 {
            return Err(Stage3Error::AdjacencyError(format!(
                "face {fi}: boundary loop winds both U and V (unsupported)"
            )));
        }
        if lp.winding_u != 0 {
            if wrap == WrapKind::V {
                return Err(Stage3Error::AdjacencyError(format!(
                    "face {fi}: boundary loops wind both U and V (unsupported)"
                )));
            }
            wrap = WrapKind::U;
        }
        if lp.winding_v != 0 {
            if wrap == WrapKind::U {
                return Err(Stage3Error::AdjacencyError(format!(
                    "face {fi}: boundary loops wind both U and V (unsupported)"
                )));
            }
            wrap = WrapKind::V;
        }
    }
    let n_winding = loops
        .iter()
        .filter(|l| l.winding_u != 0 || l.winding_v != 0)
        .count();
    if n_winding > 2 {
        return Err(Stage3Error::AdjacencyError(format!(
            "face {fi}: {n_winding} winding boundary loops (at most 2 supported)"
        )));
    }
    Ok(FaceScan { chart, surface_override, samples, loops, wrap })
}

/// The wrap-direction coordinate of a 3D point on a face's chart (U for
/// U-wrapped faces, V for V-wrapped), reduced mod 2π.
fn wrap_coord(chart: &SurfaceChart, kind: WrapKind, p: &[f64; 3]) -> f64 {
    let (u, v, _) = chart.point_uv(p);
    match kind {
        WrapKind::V => v.rem_euclid(TWO_PI),
        _ => u,
    }
}

/// Congruence distance of two angles mod 2π.
fn angle_dist(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(TWO_PI);
    d.min(TWO_PI - d)
}

/// Choose seam anchors for every wrapped face, mutating closed-loop edges:
/// each re-anchored closed edge is rebuilt to start at its anchor point and
/// gains the anchor vertex at both ends. Because the anchor vertex is part
/// of the shared edge, every face using that edge sees the same anchor —
/// consistency across faces is forced by the sharing, not by tolerances.
fn plan_seam_anchors(
    output: &Stage3Output,
    scans: &[Option<FaceScan>],
    ctx: &mut AssemblyCtx,
    config: &Config,
) -> Result<HashMap<usize, WrapPlan>, Stage3Error> {
    // Committed anchor parameter per closed-loop edge.
    let mut anchor_t: HashMap<usize, f64> = HashMap::new();
    let mut plans: HashMap<usize, WrapPlan> = HashMap::new();

    // The anchor point of a winding loop, if already determined.
    let anchor_point = |src: &AnchorSource, anchor_t: &HashMap<usize, f64>| -> Option<[f64; 3]> {
        match src {
            AnchorSource::Vertex(vid) => Some(output.vertices[*vid].point),
            AnchorSource::ClosedEdge(ei) => anchor_t.get(ei).map(|&t| {
                let mut f = 0.0;
                let mut l = 0.0;
                let curve = b_rep::Tool::curve_edge_real2(&ctx.topo_edges[*ei], &mut f, &mut l);
                let p = curve.get().value(t);
                [p.x(), p.y(), p.z()]
            }),
        }
    };

    // A winding loop's natural anchor source: a real vertex when the loop
    // has one (lowest id for determinism), else its single closed edge.
    let loop_anchor_source = |lp: &BoundaryLoop| -> Result<AnchorSource, Stage3Error> {
        let mut best: Option<usize> = None;
        for item in &lp.items {
            for &vid in &ctx.edge_vertices[item.ei] {
                if vid != usize::MAX && best.is_none_or(|b| vid < b) {
                    best = Some(vid);
                }
            }
        }
        if let Some(vid) = best {
            return Ok(AnchorSource::Vertex(vid));
        }
        if lp.items.len() == 1 {
            return Ok(AnchorSource::ClosedEdge(lp.items[0].ei));
        }
        Err(Stage3Error::AdjacencyError(
            "winding boundary loop has multiple edges but no vertices".into(),
        ))
    };

    // Find the parameter on a closed edge where the wrap coordinate equals
    // `target` (mod 2π), by bisection on the branch-continuous samples.
    let solve_anchor_t = |ei: usize,
                          uv: &EdgeUv,
                          chart: &SurfaceChart,
                          kind: WrapKind,
                          target: f64|
     -> Result<f64, Stage3Error> {
        let coords = match kind {
            WrapKind::V => &uv.vs,
            _ => &uv.us,
        };
        // Lift target into the samples' branch range.
        let (lo, hi) = (coords[0].min(*coords.last().unwrap()), coords[0].max(*coords.last().unwrap()));
        let mut tgt = target;
        while tgt < lo {
            tgt += TWO_PI;
        }
        while tgt > hi {
            tgt -= TWO_PI;
        }
        if tgt < lo - 1e-9 {
            return Err(Stage3Error::AdjacencyError(format!(
                "edge {ei}: seam azimuth not reachable on closed loop"
            )));
        }
        let k = (0..coords.len() - 1)
            .find(|&i| (coords[i] - tgt) * (coords[i + 1] - tgt) <= 0.0)
            .ok_or_else(|| {
                Stage3Error::AdjacencyError(format!(
                    "edge {ei}: seam azimuth not bracketed on closed loop"
                ))
            })?;
        let mut f = 0.0;
        let mut l = 0.0;
        let curve = b_rep::Tool::curve_edge_real2(&ctx.topo_edges[ei], &mut f, &mut l);
        let c = curve.get();
        let (mut ta, mut tb) = (uv.ts[k], uv.ts[k + 1]);
        let (mut ca, _cb) = (coords[k], coords[k + 1]);
        for _ in 0..60 {
            let tm = 0.5 * (ta + tb);
            let pt = c.value(tm);
            let raw = wrap_coord(chart, kind, &[pt.x(), pt.y(), pt.z()]);
            // Lift onto the local branch near ca.
            let cm = raw - TWO_PI * ((raw - ca) / TWO_PI).round();
            if (cm - tgt) * (ca - tgt) <= 0.0 {
                tb = tm;
            } else {
                ta = tm;
                ca = cm;
            }
        }
        Ok(0.5 * (ta + tb))
    };

    // Fixpoint: plan faces that already have a determined anchor first, so
    // chains of wrapped faces adopt a single consistent seam azimuth.
    let wrapped: Vec<usize> = scans
        .iter()
        .enumerate()
        .filter(|(_, s)| s.as_ref().is_some_and(|s| s.wrap != WrapKind::None))
        .map(|(fi, _)| fi)
        .collect();
    let mut pending: Vec<usize> = wrapped.clone();
    while !pending.is_empty() {
        let mut planned_one = false;
        // Two rounds: first faces with a determined anchor, then (if none
        // progressed) force the lowest-index pending face with defaults.
        for force in [false, true] {
            let mut i = 0;
            while i < pending.len() {
                let fi = pending[i];
                let scan = scans[fi].as_ref().unwrap();
                let winding: Vec<&BoundaryLoop> = scan
                    .loops
                    .iter()
                    .filter(|l| l.winding_u != 0 || l.winding_v != 0)
                    .collect();
                let sources: Vec<AnchorSource> = winding
                    .iter()
                    .map(|lp| loop_anchor_source(lp))
                    .collect::<Result<_, _>>()?;
                let det: Vec<Option<[f64; 3]>> =
                    sources.iter().map(|s| anchor_point(s, &anchor_t)).collect();
                if !force && det.iter().all(|d| d.is_none()) {
                    i += 1;
                    continue;
                }
                // Seam azimuth from the first determined anchor, or commit
                // the first loop's default anchor.
                let seam_c = if let Some(p) = det.iter().flatten().next() {
                    wrap_coord(&scan.chart, scan.wrap, p)
                } else {
                    match sources[0] {
                        AnchorSource::ClosedEdge(ei) => {
                            let uv = &scan.samples[&ei];
                            anchor_t.insert(ei, uv.ts[0]);
                            match scan.wrap {
                                WrapKind::V => uv.vs[0].rem_euclid(TWO_PI),
                                _ => uv.us[0].rem_euclid(TWO_PI),
                            }
                        }
                        AnchorSource::Vertex(_) => unreachable!("vertex anchors are always determined"),
                    }
                };
                for li in 0..winding.len() {
                    match det[li] {
                        Some(p) => {
                            let c = wrap_coord(&scan.chart, scan.wrap, &p);
                            if angle_dist(c, seam_c) > 0.01 {
                                return Err(Stage3Error::AdjacencyError(format!(
                                    "face {fi}: winding loops have incompatible seam anchors \
                                     ({c:.4} vs {seam_c:.4} rad); edge splitting not implemented"
                                )));
                            }
                        }
                        None => {
                            let AnchorSource::ClosedEdge(ei) = sources[li] else {
                                unreachable!("vertex anchors are always determined")
                            };
                            let t = solve_anchor_t(ei, &scan.samples[&ei], &scan.chart, scan.wrap, seam_c)?;
                            anchor_t.insert(ei, t);
                        }
                    }
                }
                plans.insert(fi, WrapPlan { anchors: sources });
                pending.remove(i);
                planned_one = true;
                if force {
                    break; // re-enter the determined-anchor round
                }
            }
            if planned_one {
                break;
            }
        }
        if !planned_one {
            break;
        }
    }

    // Apply the committed anchors: rebuild each re-anchored closed edge so
    // its range starts at the anchor and the anchor vertex is at both ends.
    let builder = b_rep::Builder::new();
    let mut anchor_edges: Vec<(usize, f64)> = anchor_t.into_iter().collect();
    anchor_edges.sort_by_key(|&(ei, _)| ei);
    for (ei, t) in anchor_edges {
        let mut f = 0.0;
        let mut l = 0.0;
        let curve = b_rep::Tool::curve_edge_real2(&ctx.topo_edges[ei], &mut f, &mut l);
        let period = l - f;
        let t = if !curve.get().is_periodic() && (t - f).abs() > 1e-12 {
            if !config.quiet {
                eprintln!(
                    "  Warning: edge {ei}: closed curve is not periodic; keeping its \
                     existing start as the seam anchor"
                );
            }
            f
        } else {
            t
        };
        let p = curve.get().value(t);
        let mut mv = b_rep_builder_api::MakeVertex::new_pnt(&p);
        let vtx = mv.vertex().to_owned();
        builder.update_vertex_vertex_real(&vtx, config.vertex_tolerance_mm);
        let make_edge = b_rep_builder_api::MakeEdge::new_handlegeomcurve_vertex2_real2(
            &curve, &vtx, &vtx, t, t + period,
        );
        if !make_edge.is_done() {
            return Err(Stage3Error::AdjacencyError(format!(
                "MakeEdge failed while re-anchoring closed edge {ei}: {:?}",
                make_edge.error(),
            )));
        }
        let mut make_edge = make_edge;
        ctx.topo_edges[ei] = make_edge.edge().to_owned();
        ctx.closed_anchors.insert(ei, vtx);
    }

    Ok(plans)
}

/// Circular window over values mod 2π: the arc complementary to the largest
/// gap between consecutive values, flipped to its complement when a majority
/// of `interior` values (mod 2π) fall outside it. Returns (lo, hi) with
/// lo ∈ [0, 2π) and lo ≤ hi ≤ lo + 2π.
fn circular_window(values: &[f64], interior: &[f64]) -> (f64, f64) {
    let mut sorted: Vec<f64> = values.iter().map(|v| v.rem_euclid(TWO_PI)).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    sorted.dedup_by(|a, b| (*a - *b).abs() < 1e-10);
    if sorted.is_empty() {
        return (0.0, TWO_PI);
    }
    if sorted.len() == 1 {
        return (sorted[0], sorted[0]);
    }
    let mut max_gap = 0.0_f64;
    let mut gap_end = 0;
    for i in 0..sorted.len() {
        let next = (i + 1) % sorted.len();
        let gap = if next > i {
            sorted[next] - sorted[i]
        } else {
            sorted[next] + TWO_PI - sorted[i]
        };
        if gap > max_gap {
            max_gap = gap;
            gap_end = next;
        }
    }
    let gap_start = if gap_end == 0 { sorted.len() - 1 } else { gap_end - 1 };
    let lo = sorted[gap_end];
    let mut hi = sorted[gap_start];
    if hi < lo {
        hi += TWO_PI;
    }
    if !interior.is_empty() {
        let tol = 0.1;
        let inside = |x: f64| -> bool {
            let xx = lo + (x - lo).rem_euclid(TWO_PI);
            xx <= hi + tol || xx >= lo + TWO_PI - tol
        };
        let n_in = interior.iter().filter(|&&x| inside(x)).count();
        if 2 * n_in < interior.len() {
            let new_lo = hi.rem_euclid(TWO_PI);
            let mut new_hi = lo.rem_euclid(TWO_PI);
            if new_hi < new_lo {
                new_hi += TWO_PI;
            }
            return (new_lo, new_hi);
        }
    }
    (lo, hi)
}

/// Winding number of a closed polygon around a point (crossing count).
fn winding_number(poly: &[(f64, f64)], p: (f64, f64)) -> i32 {
    let mut w = 0i32;
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        if a.1 <= p.1 {
            if b.1 > p.1 && cross_z(a, b, p) > 0.0 {
                w += 1;
            }
        } else if b.1 <= p.1 && cross_z(a, b, p) < 0.0 {
            w -= 1;
        }
    }
    w
}

fn cross_z(a: (f64, f64), b: (f64, f64), p: (f64, f64)) -> f64 {
    (b.0 - a.0) * (p.1 - a.1) - (p.0 - a.0) * (b.1 - a.1)
}

/// UV positions of a spread of the face's mesh-region centroids: the ground
/// truth for which side of its boundary a face occupies. `u` is None for
/// pole-degenerate points (their U is meaningless).
fn mesh_interior_uvs(
    fd: &FaceDescriptor,
    stage2: &Stage2Output,
    chart: &SurfaceChart,
    pole_eps: f64,
) -> Vec<(Option<f64>, f64)> {
    const MAX_POINTS: usize = 24;
    let mfaces = surface_faces(&stage2.selected_surfaces[fd.selected_surface_idx], stage2);
    let stride = (mfaces.len() / MAX_POINTS).max(1);
    let mut pts = Vec::new();
    for i in (0..mfaces.len()).step_by(stride) {
        let c = compute_mesh_face_centroid(mfaces[i], &stage2.mesh);
        let (u, v, rad) = chart.point_uv(&c);
        let u = if rad < pole_eps { None } else { Some(u) };
        pts.push((u, v));
    }
    pts
}

/// How to build one edge's pcurve on a face's chart.
#[derive(Default, Clone, Copy)]
struct PcurveSpec {
    /// Pins the U branch: resampled arrays are shifted by whole periods to
    /// match the first sample of the planning-time scan.
    basis_first_u: f64,
    /// Branch/window offsets applied after basis pinning.
    u_off: f64,
    v_off: f64,
    /// Pin every sample's U (or V) exactly — for seam edges, which lie
    /// exactly on an isoline at a window border.
    force_u: Option<f64>,
    force_v: Option<f64>,
    /// Snap the first/last sample (in stored-curve direction) to the shared
    /// vertex's analytic UV, so consecutive pcurves in a wire meet exactly
    /// at their common vertex even when the vertex is slightly off this
    /// face's fitted surface.
    snap_start: Option<(f64, f64)>,
    snap_end: Option<(f64, f64)>,
}

/// Interpolate a pcurve through analytically-computed UV samples of an
/// edge's stored curve, densifying until the pcurve reproduces the 3D curve
/// within a target deviation. Returns (pcurve, max observed deviation).
///
/// The pcurve is parameterized exactly like the edge's stored curve, so
/// SameParameter holds by construction (within the returned deviation).
fn build_pcurve(
    edge: &topo_ds::Edge,
    chart: &SurfaceChart,
    spec: &PcurveSpec,
    config: &Config,
) -> Result<(OwnedPtr<geom2d::Handle_Geom2d_Curve>, f64), Stage3Error> {
    let pole_eps = pole_radius_eps(config);
    let dev_target = (config.vertex_tolerance_mm * 0.1).max(1e-7);
    let mut f = 0.0;
    let mut l = 0.0;
    let curve = b_rep::Tool::curve_edge_real2(edge, &mut f, &mut l);
    let c = curve.get();

    let mut best: Option<(OwnedPtr<geom2d::Handle_Geom2d_Curve>, f64)> = None;
    let mut n = PCURVE_BASE_SAMPLES;
    loop {
        let mut uv = sample_edge_uv(edge, chart, n, pole_eps).ok_or_else(|| {
            Stage3Error::AdjacencyError("edge lost its 3D curve during assembly".into())
        })?;
        // Pin the branch to the basis used during face planning.
        let shift = if chart.u_periodic() {
            TWO_PI * ((spec.basis_first_u - uv.us[0]) / TWO_PI).round()
        } else {
            0.0
        };
        for u in uv.us.iter_mut() {
            *u += shift + spec.u_off;
        }
        for v in uv.vs.iter_mut() {
            *v += spec.v_off;
        }
        if let Some(exact_u) = spec.force_u {
            for u in uv.us.iter_mut() {
                *u = exact_u;
            }
        }
        if let Some(exact_v) = spec.force_v {
            for v in uv.vs.iter_mut() {
                *v = exact_v;
            }
        }
        // Snap endpoints to their shared vertex UVs (branch-lifted in
        // periodic directions only).
        let snap = |i: usize, target: (f64, f64), uv: &mut EdgeUv| {
            let tu = if chart.u_periodic() {
                target.0 + TWO_PI * ((uv.us[i] - target.0) / TWO_PI).round()
            } else {
                target.0
            };
            let tv = if chart.v_periodic() {
                target.1 + TWO_PI * ((uv.vs[i] - target.1) / TWO_PI).round()
            } else {
                target.1
            };
            uv.us[i] = tu;
            uv.vs[i] = tv;
        };
        if let Some(t) = spec.snap_start {
            if !uv.start_at_pole {
                snap(0, t, &mut uv);
            }
        }
        if let Some(t) = spec.snap_end {
            if !uv.end_at_pole {
                let last = uv.us.len() - 1;
                snap(last, t, &mut uv);
            }
        }

        // Geom2dAPI_Interpolate throws Standard_ConstructionError on
        // consecutive points closer than its tolerance, and UV samples are
        // not guaranteed distinct: an edge (nearly) perpendicular to this
        // face's fitted surface projects to (nearly) one UV point — observed
        // on 3DBenchy planar-only fits, where corner refinement legitimately
        // offsets both endpoints of a tangent plane-plane edge along the
        // face normal within the fit slop. Interpolate through the samples
        // that are pairwise-distinct in UV; if the whole projection is
        // pointlike, use a degree-1 pcurve (possibly constant) between the
        // endpoint samples instead — the edge's 3D extent then shows up in
        // the returned deviation, which the caller folds into the edge
        // tolerance like any other pcurve deviation.
        const UV_DISTINCT_EPS: f64 = 1e-8; // 100× the interpolation tolerance
        let eps2 = UV_DISTINCT_EPS * UV_DISTINCT_EPS;
        let d2 = |i: usize, j: usize| -> f64 {
            (uv.us[i] - uv.us[j]).powi(2) + (uv.vs[i] - uv.vs[j]).powi(2)
        };
        let last = uv.ts.len() - 1;
        let mut kept: Vec<usize> = vec![0];
        for i in 1..last {
            if d2(i, *kept.last().unwrap()) >= eps2 {
                kept.push(i);
            }
        }
        // The final (possibly vertex-snapped) sample must be kept exactly:
        // drop earlier samples that crowd it.
        while kept.len() > 1 && d2(last, *kept.last().unwrap()) < eps2 {
            kept.pop();
        }
        let pointlike = kept.len() == 1 && d2(last, 0) < eps2;
        if !pointlike {
            kept.push(last);
        }

        let pcurve = if pointlike {
            let mut poles = n_collection::Array1_gp_Pnt2d::new_int2(1, 2);
            poles.set_value(1, &gp::Pnt2d::new_real2(uv.us[0], uv.vs[0]));
            poles.set_value(2, &gp::Pnt2d::new_real2(uv.us[last], uv.vs[last]));
            let mut knots = n_collection::Array1_double::new_int2(1, 2);
            knots.set_value(1, &uv.ts[0]);
            knots.set_value(2, &uv.ts[last]);
            let mut mults = n_collection::Array1_int::new_int2(1, 2);
            mults.set_value(1, &2);
            mults.set_value(2, &2);
            let bs = geom2d::BSplineCurve::new_arr1_pnt2d_arr1_double_arr1_int_int_bool(
                &poles, &knots, &mults, 1, false,
            );
            geom2d::BSplineCurve::to_handle(bs).to_handle_curve()
        } else {
            let np = kept.len() as i32;
            let mut pts_h = n_collection::HArray1_gp_Pnt2d::new_int2(1, np);
            let mut params_h = n_collection::HArray1_double::new_int2(1, np);
            for (slot, &i) in kept.iter().enumerate() {
                let pnt2d = gp::Pnt2d::new_real2(uv.us[i], uv.vs[i]);
                pts_h.set_value((slot + 1) as i32, &pnt2d);
                params_h.set_value((slot + 1) as i32, &uv.ts[i]);
            }
            let pts_handle = n_collection::HArray1_gp_Pnt2d::to_handle(pts_h);
            let params_handle = n_collection::HArray1_double::to_handle(params_h);
            let mut interp = geom2d_api::Interpolate::new_harr1_pnt2d_harr1_double_bool_real(
                &pts_handle, &params_handle, false, 1e-10,
            );
            interp.perform();
            if !interp.is_done() {
                return Err(Stage3Error::AdjacencyError(
                    "pcurve interpolation failed".into(),
                ));
            }
            interp.curve().to_handle_curve()
        };

        // Deviation of the pcurve (mapped through the chart) from the 3D
        // curve, at the interpolation samples and their midpoints. Endpoints
        // matter: OCCT's in-context BRepCheck validates pcurve/curve
        // agreement over the whole range against the edge tolerance, and a
        // vertex slightly off this face's fitted surface shows up exactly at
        // the (possibly vertex-snapped) pcurve ends.
        let pc = pcurve.get();
        let mut dev = 0.0_f64;
        {
            let mut check = |t: f64| {
                let p2 = pc.value(t);
                let p3 = chart.point_at(p2.x(), p2.y());
                let q = c.value(t);
                let d = ((p3[0] - q.x()).powi(2)
                    + (p3[1] - q.y()).powi(2)
                    + (p3[2] - q.z()).powi(2))
                .sqrt();
                dev = dev.max(d);
            };
            for i in 0..uv.ts.len() {
                check(uv.ts[i]);
                if i + 1 < uv.ts.len() {
                    check(0.5 * (uv.ts[i] + uv.ts[i + 1]));
                }
            }
        }
        // A pointlike projection stays pointlike at any sample density —
        // densification cannot separate it.
        if pointlike {
            return Ok((pcurve, dev));
        }
        // Densify until the target is met, but stop when the deviation has
        // plateaued (e.g. the 3D curve is off this fitted surface by its
        // fit residual, which no pcurve can represent away).
        let plateaued = best.as_ref().is_some_and(|(_, bd)| dev > 0.7 * bd);
        let better = best.as_ref().is_none_or(|(_, bd)| dev < *bd);
        if better {
            best = Some((pcurve, dev));
        }
        if dev <= dev_target || n >= 257 || plateaued {
            break;
        }
        n = n * 2 - 1;
    }
    Ok(best.unwrap())
}

/// Create a degenerate edge at a UV pole: no 3D curve, a straight pcurve at
/// constant V from (u_from, v_pole) to (u_to, v_pole), the pole vertex at
/// both ends. This is the single uniform mechanism that closes a face's UV
/// boundary across cone apexes and sphere poles.
fn make_degenerate_edge(
    builder: &b_rep::Builder,
    pole_vertex: &topo_ds::Vertex,
    u_from: f64,
    u_to: f64,
    v_pole: f64,
    surface: &geom::Handle_Geom_Surface,
    config: &Config,
) -> OwnedPtr<topo_ds::Edge> {
    let identity_loc = top_loc::Location::new();
    let mut degen = topo_ds::Edge::new();
    builder.make_edge_edge(&mut degen);
    builder.degenerated(&degen, true);

    let du = u_to - u_from;
    let pcurve_origin = gp::Pnt2d::new_real2(u_from, v_pole);
    let pcurve_dir = gp::Dir2d::new_real2(if du >= 0.0 { 1.0 } else { -1.0 }, 0.0);
    let pcurve_line = geom2d::Line::new_pnt2d_dir2d(&pcurve_origin, &pcurve_dir);
    let pcurve_handle = geom2d::Line::to_handle(pcurve_line).to_handle_curve();
    builder.update_edge_edge_handlegeom2dcurve_handlegeomsurface_location_real(
        &degen,
        &pcurve_handle,
        surface,
        &identity_loc,
        config.vertex_tolerance_mm,
    );
    builder.range_edge_real2_bool(&degen, 0.0, du.abs(), false);

    let v_fwd = pole_vertex.as_shape().oriented(top_abs::Orientation::Forward);
    let v_rev = pole_vertex.as_shape().oriented(top_abs::Orientation::Reversed);
    builder.add(degen.as_shape_mut(), &v_fwd);
    builder.add(degen.as_shape_mut(), &v_rev);
    builder.update_vertex_vertex_real_edge_real(
        pole_vertex,
        0.0,
        &degen,
        config.vertex_tolerance_mm,
    );
    degen
}

/// One element of a face wire, in traversal order.
enum WireElement {
    /// A shared boundary edge, traversed forward or reversed.
    Edge { ei: usize, forward: bool },
    /// The face's seam edge occurrence.
    Seam { forward: bool },
    /// A face-local degenerate pole edge (index into the face's list).
    Degenerate(usize),
}

/// A boundary loop with its final window placement, ready for emission.
struct PlacedLoop {
    lp: BoundaryLoop,
    /// Global (u, v) shift added to every item of the loop.
    shift: (f64, f64),
}

impl PlacedLoop {
    /// Traversal-order polygon of the loop's UV samples.
    fn polygon(&self, samples: &HashMap<usize, EdgeUv>) -> Vec<(f64, f64)> {
        let mut poly = Vec::new();
        for (idx, item) in self.lp.items.iter().enumerate() {
            let s = &samples[&item.ei];
            let (du, dv) = (
                self.lp.offsets[idx].0 + self.shift.0,
                self.lp.offsets[idx].1 + self.shift.1,
            );
            if item.forward {
                for i in 0..s.us.len() {
                    poly.push((s.us[i] + du, s.vs[i] + dv));
                }
            } else {
                for i in (0..s.us.len()).rev() {
                    poly.push((s.us[i] + du, s.vs[i] + dv));
                }
            }
        }
        poly
    }

    fn start_uv(&self, samples: &HashMap<usize, EdgeUv>) -> (f64, f64) {
        let item = &self.lp.items[0];
        let s = &samples[&item.ei];
        let (u, v) = if item.forward { s.first_uv() } else { s.last_uv() };
        (
            u + self.lp.offsets[0].0 + self.shift.0,
            v + self.lp.offsets[0].1 + self.shift.1,
        )
    }

    /// Reverse the loop's traversal direction in place.
    fn reverse(&mut self, samples: &HashMap<usize, EdgeUv>, chart: &SurfaceChart) {
        let mut items: Vec<LoopItem> = self
            .lp
            .items
            .iter()
            .rev()
            .map(|it| LoopItem { ei: it.ei, forward: !it.forward })
            .collect();
        // Preserve the loop's absolute branch: re-align from the (previously
        // last) item's offset.
        let last_off = *self.lp.offsets.last().unwrap();
        let first = items.remove(0);
        let mut realigned = align_loop(
            {
                let mut v = vec![first];
                v.extend(items);
                v
            },
            samples,
            chart,
        )
        .expect("reversing an aligned loop cannot fail");
        let delta = (last_off.0 - realigned.offsets[0].0, last_off.1 - realigned.offsets[0].1);
        for off in realigned.offsets.iter_mut() {
            off.0 += delta.0;
            off.1 += delta.1;
        }
        realigned.winding_u = -self.lp.winding_u;
        realigned.winding_v = -self.lp.winding_v;
        self.lp = realigned;
    }

    /// Rotate the loop so traversal starts at the item whose start vertex is
    /// `vid`.
    fn rotate_to_vertex(
        &mut self,
        vid: usize,
        edge_vertices: &[[usize; 2]],
        samples: &HashMap<usize, EdgeUv>,
        chart: &SurfaceChart,
    ) -> Result<(), Stage3Error> {
        let k = self
            .lp
            .items
            .iter()
            .position(|it| {
                let [v0, v1] = edge_vertices[it.ei];
                (if it.forward { v0 } else { v1 }) == vid
            })
            .ok_or_else(|| {
                Stage3Error::AdjacencyError(format!(
                    "anchor vertex {vid} is not on the winding loop"
                ))
            })?;
        if k == 0 {
            return Ok(());
        }
        let base_off = self.lp.offsets[k];
        let mut items = self.lp.items.clone();
        items.rotate_left(k);
        let mut realigned = align_loop(items, samples, chart)?;
        let delta = (base_off.0 - realigned.offsets[0].0, base_off.1 - realigned.offsets[0].1);
        for off in realigned.offsets.iter_mut() {
            off.0 += delta.0;
            off.1 += delta.1;
        }
        self.lp = realigned;
        Ok(())
    }
}

/// Result of assembling one face.
struct AssembledFace {
    face: OwnedPtr<topo_ds::Face>,
    concave: bool,
}

/// The seam edge of a wrapped face, plus everything needed to emit it.
struct SeamEdge {
    edge: OwnedPtr<topo_ds::Edge>,
    /// Window borders: pcurve positions of the low and high occurrence.
    border_lo: f64,
    border_hi: f64,
}

/// Assemble one face from the shared edges: pcurves, seam and degenerate
/// edges, wires, and the TopoDS_Face itself — all with BRep_Builder, never
/// copying a shared edge.
fn assemble_face(
    fi: usize,
    output: &Stage3Output,
    ctx: &AssemblyCtx,
    plans: &HashMap<usize, WrapPlan>,
    config: &Config,
) -> Result<AssembledFace, Stage3Error> {
    let fd = &output.face_descriptors[fi];
    let builder = b_rep::Builder::new();
    let identity_loc = top_loc::Location::new();
    let pole_eps = pole_radius_eps(config);
    let scan = scan_face(fi, output, ctx, config)?;
    let FaceScan { chart, surface_override, samples, loops, wrap } = scan;
    let surf: &geom::Handle_Geom_Surface = surface_override.as_deref().unwrap_or(&fd.surface);
    let interior = mesh_interior_uvs(fd, &output.stage2, &chart, pole_eps);

    let mut winding: Vec<PlacedLoop> = Vec::new();
    let mut plain: Vec<PlacedLoop> = Vec::new();
    for lp in loops {
        let placed = PlacedLoop { lp, shift: (0.0, 0.0) };
        if placed.lp.winding_u != 0 || placed.lp.winding_v != 0 {
            winding.push(placed);
        } else {
            plain.push(placed);
        }
    }

    // ------------------------------------------------------------------
    // Window placement, seam, and pole handling.
    // ------------------------------------------------------------------
    // The window is the face's UV extent: [u_lo, u_hi] × [v_lo, v_hi].
    // Wrapped faces span a full period in the wrap direction, bordered by
    // the seam; everything else comes from the boundary samples (with mesh
    // interior points disambiguating periodic windows).
    let mut seam: Option<SeamEdge> = None;
    let mut degenerates: Vec<OwnedPtr<topo_ds::Edge>> = Vec::new();
    // Wire emission plans (filled per wrap case below).
    let mut composite: Vec<WireElement> = Vec::new();

    // Resolve one winding-loop anchor: its vertex handle and 3D position.
    let resolve_anchor = |src: &AnchorSource,
                          ctx: &AssemblyCtx|
     -> Result<(OwnedPtr<topo_ds::Vertex>, [f64; 3], usize), Stage3Error> {
        match src {
            AnchorSource::Vertex(vid) => Ok((
                ctx.topo_vertices[*vid].to_owned(),
                output.vertices[*vid].point,
                *vid,
            )),
            AnchorSource::ClosedEdge(ei) => {
                let v = ctx.closed_anchors.get(ei).ok_or_else(|| {
                    Stage3Error::AdjacencyError(format!(
                        "closed edge {ei} was not re-anchored during seam planning"
                    ))
                })?;
                let p = b_rep::Tool::pnt(v);
                Ok((v.to_owned(), [p.x(), p.y(), p.z()], usize::MAX))
            }
        }
    };

    match wrap {
        WrapKind::None => {
            // Lift every loop (and nothing else) into a consistent window.
            if chart.u_periodic() {
                let all_us: Vec<f64> = winding
                    .iter()
                    .chain(plain.iter())
                    .flat_map(|pl| pl.polygon(&samples))
                    .map(|(u, _)| u)
                    .collect();
                let int_us: Vec<f64> = interior.iter().filter_map(|(u, _)| *u).collect();
                let (u_lo, _u_hi) = circular_window(&all_us, &int_us);
                for pl in plain.iter_mut() {
                    let (u0, _) = pl.start_uv(&samples);
                    pl.shift.0 = u_lo + (u0 - u_lo).rem_euclid(TWO_PI) - u0;
                }
            }
            if chart.v_periodic() {
                let all_vs: Vec<f64> = plain
                    .iter()
                    .flat_map(|pl| pl.polygon(&samples))
                    .map(|(_, v)| v)
                    .collect();
                let int_vs: Vec<f64> = interior.iter().map(|(_, v)| v.rem_euclid(TWO_PI)).collect();
                let (v_lo, _v_hi) = circular_window(&all_vs, &int_vs);
                for pl in plain.iter_mut() {
                    let (_, v0) = pl.start_uv(&samples);
                    pl.shift.1 = v_lo + (v0 - v_lo).rem_euclid(TWO_PI) - v0;
                }
            }
        }
        WrapKind::U => {
            let plan = plans.get(&fi).ok_or_else(|| {
                Stage3Error::AdjacencyError(format!("face {fi}: missing seam plan"))
            })?;
            if plan.anchors.len() != winding.len() {
                return Err(Stage3Error::AdjacencyError(format!(
                    "face {fi}: seam plan does not match winding loops"
                )));
            }
            // Resolve anchors and start each winding loop at its anchor.
            let mut anchors = Vec::new();
            for (pl, src) in winding.iter_mut().zip(&plan.anchors) {
                let (vtx, point, vid) = resolve_anchor(src, ctx)?;
                if vid != usize::MAX {
                    pl.rotate_to_vertex(vid, &ctx.edge_vertices, &samples, &chart)?;
                }
                let (_, v, _) = chart.point_uv(&point);
                anchors.push((vtx, point, v));
            }

            // Decide the loops' V placement: which winding loop is the
            // bottom border, and (torus) which V branch the face occupies.
            let (bottom_idx, v_bottom, top): (usize, f64, Option<(usize, f64)>) = if winding.len() == 2 {
                let v0 = anchors[0].2;
                let v1 = anchors[1].2;
                if chart.v_periodic() {
                    let int_vs: Vec<f64> =
                        interior.iter().map(|(_, v)| v.rem_euclid(TWO_PI)).collect();
                    let (w_lo, w_hi) = circular_window(&[v0, v1], &int_vs);
                    // Bottom = the anchor congruent to the window's low end.
                    let (bi, ti) = if angle_dist(v0, w_lo.rem_euclid(TWO_PI))
                        <= angle_dist(v1, w_lo.rem_euclid(TWO_PI))
                    {
                        (0, 1)
                    } else {
                        (1, 0)
                    };
                    (bi, w_lo, Some((ti, w_hi)))
                } else if v0 <= v1 {
                    (0, v0, Some((1, v1)))
                } else {
                    (1, v1, Some((0, v0)))
                }
            } else {
                // One winding loop: the other border is a pole.
                let v_loop = anchors[0].2;
                let mean_int_v = interior.iter().map(|(_, v)| v).sum::<f64>()
                    / interior.len().max(1) as f64;
                let positive_side = mean_int_v > v_loop;
                let (v_pole, pole_pt) = chart.pole(positive_side).ok_or_else(|| {
                    Stage3Error::AdjacencyError(format!(
                        "face {fi}: single winding loop but the surface has no pole on \
                         the region side"
                    ))
                })?;
                let pt = gp::Pnt::new_real3(pole_pt[0], pole_pt[1], pole_pt[2]);
                let mut mv = b_rep_builder_api::MakeVertex::new_pnt(&pt);
                let pole_vtx = mv.vertex().to_owned();
                builder.update_vertex_vertex_real(&pole_vtx, config.vertex_tolerance_mm);
                if positive_side {
                    // Pole above the loop: the loop is the bottom border.
                    anchors.push((pole_vtx, pole_pt, v_pole));
                    (0, v_loop, None) // top handled below via pole marker
                } else {
                    anchors.push((pole_vtx, pole_pt, v_pole));
                    (usize::MAX, v_pole, Some((0, v_loop)))
                }
            };

            // The seam azimuth: the bottom anchor's U (or the loop anchor's
            // U when the bottom border is a pole).
            let loop_anchor_u = |idx: usize, winding: &[PlacedLoop]| -> f64 {
                winding[idx].start_uv(&samples).0
            };
            let (u_seam, v_lo, v_hi, bottom_loop, top_loop, pole_at) = match (bottom_idx, top) {
                (usize::MAX, Some((ti, v_top))) => {
                    // Pole at the bottom (cone apex cap).
                    let us = loop_anchor_u(ti, &winding).rem_euclid(TWO_PI);
                    (us, v_bottom, v_top, None, Some(ti), Some(false))
                }
                (bi, None) => {
                    // Pole at the top.
                    let us = loop_anchor_u(bi, &winding).rem_euclid(TWO_PI);
                    let v_pole = anchors.last().unwrap().2;
                    (us, v_bottom, v_pole, Some(bi), None, Some(true))
                }
                (bi, Some((ti, v_top))) => {
                    let us = loop_anchor_u(bi, &winding).rem_euclid(TWO_PI);
                    (us, v_bottom, v_top, Some(bi), Some(ti), None)
                }
            };
            let u_hi = u_seam + TWO_PI;

            // Place the winding loops: bottom runs +U starting at u_seam,
            // top runs −U starting at u_seam + 2π. Torus faces also lift V
            // into the window.
            let place = |idx: usize, start_u: f64, want_winding: i32, winding: &mut [PlacedLoop], v_ref: f64| {
                let pl = &mut winding[idx];
                if pl.lp.winding_u != want_winding {
                    pl.reverse(&samples, &chart);
                }
                let (u0, v0) = pl.start_uv(&samples);
                pl.shift.0 += start_u - u0;
                if chart.v_periodic() {
                    pl.shift.1 += TWO_PI * ((v_ref - v0) / TWO_PI).round();
                }
            };
            if let Some(bi) = bottom_loop {
                place(bi, u_seam, 1, &mut winding, v_lo);
            }
            if let Some(ti) = top_loop {
                place(ti, u_hi, -1, &mut winding, v_hi);
            }

            // Build the seam edge along the U isoline, from the bottom
            // border's vertex to the top border's vertex.
            let (bot_vtx, bot_pt, _) = match (bottom_loop, &plan.anchors) {
                (Some(bi), anchorsrc) => resolve_anchor(&anchorsrc[bi], ctx)?,
                (None, _) => {
                    let (v, p, _) = (&anchors.last().unwrap().0, anchors.last().unwrap().1, 0);
                    (v.to_owned(), p, usize::MAX)
                }
            };
            let (top_vtx, top_pt, _) = match (top_loop, &plan.anchors) {
                (Some(ti), anchorsrc) => resolve_anchor(&anchorsrc[ti], ctx)?,
                (None, _) => {
                    let (v, p, _) = (&anchors.last().unwrap().0, anchors.last().unwrap().1, 0);
                    (v.to_owned(), p, usize::MAX)
                }
            };
            let iso = surf.get().u_iso(u_seam);
            let seam_edge = make_seam_edge(
                &builder, &iso, &bot_vtx, &top_vtx, v_lo, v_hi, &bot_pt, &top_pt, config,
            )?;
            seam = Some(SeamEdge { edge: seam_edge, border_lo: u_seam, border_hi: u_hi });

            // Degenerate pole edge spanning the full window at the pole.
            match pole_at {
                Some(true) => {
                    // Pole at the top: traversed −U.
                    let pole_vtx = &anchors.last().unwrap().0;
                    degenerates.push(make_degenerate_edge(
                        &builder, pole_vtx, u_hi, u_seam, v_hi, surf, config,
                    ));
                }
                Some(false) => {
                    // Pole at the bottom (apex cap): traversed +U.
                    let pole_vtx = &anchors.last().unwrap().0;
                    degenerates.push(make_degenerate_edge(
                        &builder, pole_vtx, u_seam, u_hi, v_lo, surf, config,
                    ));
                }
                None => {}
            }

            // Emission order (counterclockwise in UV).
            match pole_at {
                None => {
                    let (bi, ti) = (bottom_loop.unwrap(), top_loop.unwrap());
                    for it in &winding[bi].lp.items {
                        composite.push(WireElement::Edge { ei: it.ei, forward: it.forward });
                    }
                    composite.push(WireElement::Seam { forward: true });
                    for it in &winding[ti].lp.items {
                        composite.push(WireElement::Edge { ei: it.ei, forward: it.forward });
                    }
                    composite.push(WireElement::Seam { forward: false });
                }
                Some(true) => {
                    let bi = bottom_loop.unwrap();
                    for it in &winding[bi].lp.items {
                        composite.push(WireElement::Edge { ei: it.ei, forward: it.forward });
                    }
                    composite.push(WireElement::Seam { forward: true });
                    composite.push(WireElement::Degenerate(0));
                    composite.push(WireElement::Seam { forward: false });
                }
                Some(false) => {
                    let ti = top_loop.unwrap();
                    composite.push(WireElement::Degenerate(0));
                    composite.push(WireElement::Seam { forward: true });
                    for it in &winding[ti].lp.items {
                        composite.push(WireElement::Edge { ei: it.ei, forward: it.forward });
                    }
                    composite.push(WireElement::Seam { forward: false });
                }
            }

            // Lift hole loops into the window.
            for pl in plain.iter_mut() {
                let (u0, _) = pl.start_uv(&samples);
                pl.shift.0 = u_seam + (u0 - u_seam).rem_euclid(TWO_PI) - u0;
                if chart.v_periodic() {
                    let (_, v0) = pl.start_uv(&samples);
                    pl.shift.1 = TWO_PI * ((0.5 * (v_lo + v_hi) - v0) / TWO_PI).round();
                }
            }
        }
        WrapKind::V => {
            // Only tori wrap in V (tube-through faces, e.g. pipe elbows).
            let plan = plans.get(&fi).ok_or_else(|| {
                Stage3Error::AdjacencyError(format!("face {fi}: missing seam plan"))
            })?;
            if winding.len() != 2 {
                return Err(Stage3Error::AdjacencyError(format!(
                    "face {fi}: V-wrapped face must have exactly two winding loops \
                     (found {})",
                    winding.len()
                )));
            }
            let mut anchors = Vec::new();
            for (pl, src) in winding.iter_mut().zip(&plan.anchors) {
                let (vtx, point, vid) = resolve_anchor(src, ctx)?;
                if vid != usize::MAX {
                    pl.rotate_to_vertex(vid, &ctx.edge_vertices, &samples, &chart)?;
                }
                let (u, _, _) = chart.point_uv(&point);
                anchors.push((vtx, point, u));
            }
            // U window between the two circles, disambiguated by interior.
            let int_us: Vec<f64> = interior.iter().filter_map(|(u, _)| *u).collect();
            let (w_lo, w_hi) = circular_window(&[anchors[0].2, anchors[1].2], &int_us);
            let (li, ri) = if angle_dist(anchors[0].2, w_lo.rem_euclid(TWO_PI))
                <= angle_dist(anchors[1].2, w_lo.rem_euclid(TWO_PI))
            {
                (0, 1)
            } else {
                (1, 0)
            };
            let (u_lo, u_hi) = (w_lo, w_hi);
            let v_seam = winding[ri].start_uv(&samples).1.rem_euclid(TWO_PI);
            let v_hi = v_seam + TWO_PI;

            // Right loop runs +V from v_seam; left loop runs −V from v_hi.
            {
                let pl = &mut winding[ri];
                if pl.lp.winding_v != 1 {
                    pl.reverse(&samples, &chart);
                }
                let (u0, v0) = pl.start_uv(&samples);
                pl.shift = (u_hi - u0, v_seam - v0);
            }
            {
                let pl = &mut winding[li];
                if pl.lp.winding_v != -1 {
                    pl.reverse(&samples, &chart);
                }
                let (u0, v0) = pl.start_uv(&samples);
                pl.shift = (u_lo - u0, v_hi - v0);
            }

            if config.verbose {
                eprintln!(
                    "    V-wrap: li={li} ri={ri} u=[{u_lo:.4},{u_hi:.4}] v_seam={v_seam:.6} \
                     right(shift={:?} wv={}) left(shift={:?} wv={})",
                    winding[ri].shift,
                    winding[ri].lp.winding_v,
                    winding[li].shift,
                    winding[li].lp.winding_v,
                );
            }
            let iso = surf.get().v_iso(v_seam);
            let seam_edge = make_seam_edge(
                &builder,
                &iso,
                &anchors[li].0,
                &anchors[ri].0,
                u_lo,
                u_hi,
                &anchors[li].1,
                &anchors[ri].1,
                config,
            )?;
            seam = Some(SeamEdge { edge: seam_edge, border_lo: v_seam, border_hi: v_hi });

            composite.push(WireElement::Seam { forward: true });
            for it in &winding[ri].lp.items {
                composite.push(WireElement::Edge { ei: it.ei, forward: it.forward });
            }
            composite.push(WireElement::Seam { forward: false });
            for it in &winding[li].lp.items {
                composite.push(WireElement::Edge { ei: it.ei, forward: it.forward });
            }

            for pl in plain.iter_mut() {
                let (u0, v0) = pl.start_uv(&samples);
                pl.shift = (
                    u_lo + (u0 - u_lo).rem_euclid(TWO_PI) - u0,
                    v_seam + (v0 - v_seam).rem_euclid(TWO_PI) - v0,
                );
            }
        }
    }

    // ------------------------------------------------------------------
    // Degenerate edges for pole gaps inside plain loops (a boundary loop
    // passing through a pole vertex, e.g. a pie slice ending at the apex).
    // ------------------------------------------------------------------
    // Emission element lists per plain loop, with gap degenerates spliced in.
    let mut plain_elements: Vec<Vec<WireElement>> = Vec::new();
    {
        // Orient plain loops: the outer loop winds around the mesh interior
        // (counterclockwise, winding +1); holes contain no interior points
        // and run clockwise. On wrapped faces the composite wire is the
        // outer boundary, so every plain loop is a hole.
        let int_pts: Vec<(f64, f64)> = interior
            .iter()
            .filter_map(|(u, v)| u.map(|u| (u, *v)))
            .collect();
        let loop_stats: Vec<(i32, usize)> = plain
            .iter()
            .map(|pl| {
                let poly = pl.polygon(&samples);
                let mut w_sum = 0i32;
                let mut n_hit = 0usize;
                for p in &int_pts {
                    // Lift the interior point into the loop's branch.
                    let mut q = *p;
                    if chart.u_periodic() {
                        q.0 = poly[0].0
                            + (q.0 - poly[0].0 + std::f64::consts::PI).rem_euclid(TWO_PI)
                            - std::f64::consts::PI;
                    }
                    if chart.v_periodic() {
                        q.1 = poly[0].1
                            + (q.1 - poly[0].1 + std::f64::consts::PI).rem_euclid(TWO_PI)
                            - std::f64::consts::PI;
                    }
                    let w = winding_number(&poly, q);
                    if w != 0 {
                        w_sum += w;
                        n_hit += 1;
                    }
                }
                (w_sum, n_hit)
            })
            .collect();
        let outer_idx: Option<usize> = if wrap == WrapKind::None {
            let by_hits = (0..plain.len()).max_by_key(|&i| loop_stats[i].1);
            match by_hits {
                Some(i) if loop_stats[i].1 > 0 => Some(i),
                _ => {
                    // Interior sampling failed; fall back to the largest loop.
                    if !config.quiet && plain.len() > 1 {
                        eprintln!(
                            "  Warning: face {fi}: no mesh interior point falls inside any \
                             boundary loop; taking the largest loop as the outer wire"
                        );
                    }
                    (0..plain.len()).max_by(|&a, &b| {
                        let aa = polygon_signed_area(&plain[a].polygon(&samples)).abs();
                        let bb = polygon_signed_area(&plain[b].polygon(&samples)).abs();
                        aa.partial_cmp(&bb).unwrap()
                    })
                }
            }
        } else {
            None
        };
        for (i, pl) in plain.iter_mut().enumerate() {
            let area = polygon_signed_area(&pl.polygon(&samples));
            let flip = if Some(i) == outer_idx {
                if loop_stats[i].1 > 0 {
                    loop_stats[i].0 < 0
                } else {
                    area < 0.0
                }
            } else {
                area > 0.0
            };
            if flip {
                pl.reverse(&samples, &chart);
            }
        }

        for pl in &plain {
            let mut elements = Vec::new();
            for (idx, it) in pl.lp.items.iter().enumerate() {
                elements.push(WireElement::Edge { ei: it.ei, forward: it.forward });
                for gap in &pl.lp.pole_gaps {
                    if gap.after_item == idx {
                        let [v0, v1] = ctx.edge_vertices[it.ei];
                        let vid = if it.forward { v1 } else { v0 };
                        if vid == usize::MAX {
                            return Err(Stage3Error::AdjacencyError(format!(
                                "face {fi}: pole gap after closed-loop edge {}",
                                it.ei
                            )));
                        }
                        let (v_pole, _) =
                            chart.pole(gap.positive_v_side).ok_or_else(|| {
                                Stage3Error::AdjacencyError(format!(
                                    "face {fi}: pole gap on a surface without a pole"
                                ))
                            })?;
                        degenerates.push(make_degenerate_edge(
                            &builder,
                            &ctx.topo_vertices[vid],
                            gap.u_from + pl.shift.0,
                            gap.u_to + pl.shift.0,
                            v_pole,
                            surf,
                            config,
                        ));
                        elements.push(WireElement::Degenerate(degenerates.len() - 1));
                    }
                }
            }
            plain_elements.push(elements);
        }
    }

    // ------------------------------------------------------------------
    // Pcurves: interpolate and attach one per (edge, this surface), plus
    // the two seam pcurves and each degenerate edge's line (already set).
    // ------------------------------------------------------------------
    let all_loops = winding.iter().chain(plain.iter());
    for pl in all_loops {
        for (idx, it) in pl.lp.items.iter().enumerate() {
            let vertex_uv = |vid: usize| -> Option<(f64, f64)> {
                if vid == usize::MAX {
                    return None;
                }
                let (u, v, rad) = chart.point_uv(&output.vertices[vid].point);
                (rad >= pole_eps).then_some((u, v))
            };
            let [vid_start, vid_end] = ctx.edge_vertices[it.ei];
            let spec = PcurveSpec {
                basis_first_u: samples[&it.ei].us[0],
                u_off: pl.lp.offsets[idx].0 + pl.shift.0,
                v_off: pl.lp.offsets[idx].1 + pl.shift.1,
                snap_start: vertex_uv(vid_start),
                snap_end: vertex_uv(vid_end),
                ..Default::default()
            };
            let (pcurve, dev) = build_pcurve(&ctx.topo_edges[it.ei], &chart, &spec, config)?;
            let tol = (dev * 1.01).max(config.vertex_tolerance_mm);
            builder.update_edge_edge_handlegeom2dcurve_handlegeomsurface_location_real(
                &ctx.topo_edges[it.ei],
                &pcurve,
                surf,
                &identity_loc,
                tol,
            );
            if config.verbose {
                let s = &samples[&it.ei];
                eprintln!(
                    "    Edge {}: pcurve attached (dev {dev:.2e} mm, uv ({:.4},{:.4})->({:.4},{:.4}))",
                    it.ei,
                    s.us[0] + spec.u_off,
                    s.vs[0] + spec.v_off,
                    s.us.last().unwrap() + spec.u_off,
                    s.vs.last().unwrap() + spec.v_off,
                );
            }
        }
    }
    if let Some(se) = &seam {
        let basis = sample_edge_uv(&se.edge, &chart, PCURVE_BASE_SAMPLES, pole_eps)
            .ok_or_else(|| Stage3Error::AdjacencyError("seam edge has no curve".into()))?;
        let (is_u_wrap, b0) = match wrap {
            WrapKind::U => (true, basis.us[0]),
            _ => (false, basis.us[0]),
        };
        // The seam is exactly an isoline, so the pcurve's wrap coordinate is
        // pinned to each window border exactly.
        let build_border = |border: f64| -> Result<(OwnedPtr<geom2d::Handle_Geom2d_Curve>, f64), Stage3Error> {
            let spec = PcurveSpec {
                basis_first_u: b0,
                force_u: is_u_wrap.then_some(border),
                force_v: (!is_u_wrap).then_some(border),
                ..Default::default()
            };
            build_pcurve(&se.edge, &chart, &spec, config)
        };
        // PC1 is the pcurve of the FORWARD occurrence: the high border for
        // U wraps (right side, traversed upward) and the low border for V
        // wraps (bottom side, traversed +U).
        let (pc1, dev1) = build_border(if is_u_wrap { se.border_hi } else { se.border_lo })?;
        let (pc2, dev2) = build_border(if is_u_wrap { se.border_lo } else { se.border_hi })?;
        if config.verbose {
            let mut f = 0.0;
            let mut l = 0.0;
            b_rep::Tool::range_edge_real2(&se.edge, &mut f, &mut l);
            let a1 = pc1.get().value(f);
            let b1 = pc1.get().value(l);
            let a2 = pc2.get().value(f);
            let b2 = pc2.get().value(l);
            eprintln!(
                "    Seam: borders=({:.4},{:.4}) basis.vs[0]={:.6} PC1 ({:.4},{:.4})->({:.4},{:.4}) \
                 PC2 ({:.4},{:.4})->({:.4},{:.4})",
                se.border_lo, se.border_hi, basis.vs[0],
                a1.x(), a1.y(), b1.x(), b1.y(),
                a2.x(), a2.y(), b2.x(), b2.y(),
            );
        }
        let tol = (dev1.max(dev2) * 1.01).max(config.vertex_tolerance_mm);
        builder.update_edge_edge_handlegeom2dcurve2_handlegeomsurface_location_real(
            &se.edge, &pc1, &pc2, surf, &identity_loc, tol,
        );
    }

    // ------------------------------------------------------------------
    // Wires and the face.
    // ------------------------------------------------------------------
    let emit_wire = |elements: &[WireElement]| -> OwnedPtr<topo_ds::Wire> {
        let mut wire = topo_ds::Wire::new();
        builder.make_wire(&mut wire);
        for el in elements {
            match el {
                WireElement::Edge { ei, forward } => {
                    let orient = if *forward {
                        top_abs::Orientation::Forward
                    } else {
                        top_abs::Orientation::Reversed
                    };
                    let oriented = ctx.topo_edges[*ei].as_shape().oriented(orient);
                    builder.add(wire.as_shape_mut(), &oriented);
                }
                WireElement::Seam { forward } => {
                    let se = seam.as_ref().expect("seam element requires a seam edge");
                    let orient = if *forward {
                        top_abs::Orientation::Forward
                    } else {
                        top_abs::Orientation::Reversed
                    };
                    let oriented = se.edge.as_shape().oriented(orient);
                    builder.add(wire.as_shape_mut(), &oriented);
                }
                WireElement::Degenerate(di) => {
                    let oriented = degenerates[*di]
                        .as_shape()
                        .oriented(top_abs::Orientation::Forward);
                    builder.add(wire.as_shape_mut(), &oriented);
                }
            }
        }
        wire.as_shape_mut().closed_bool(true);
        wire
    };

    let mut wires: Vec<OwnedPtr<topo_ds::Wire>> = Vec::new();
    if !composite.is_empty() {
        wires.push(emit_wire(&composite));
        for elements in &plain_elements {
            wires.push(emit_wire(elements));
        }
    } else {
        // Outer first (holes were oriented clockwise above; the outer loop
        // is the one that winds around the mesh interior).
        let mut order: Vec<usize> = (0..plain_elements.len()).collect();
        let areas: Vec<f64> = plain
            .iter()
            .map(|pl| polygon_signed_area(&pl.polygon(&samples)))
            .collect();
        order.sort_by(|&a, &b| areas[b].partial_cmp(&areas[a]).unwrap());
        for &i in &order {
            wires.push(emit_wire(&plain_elements[i]));
        }
    }

    let mut face = topo_ds::Face::new();
    builder.make_face_face_handlegeomsurface_real(&mut face, surf, config.vertex_tolerance_mm);
    for wire in &wires {
        builder.add(face.as_shape_mut(), wire.as_shape());
    }
    builder.natural_restriction(&face, false);

    let surface = &output.stage2.selected_surfaces[fd.selected_surface_idx];
    let concave = match surface {
        SelectedSurface::Conical(idx) => !output.stage2.conical_hypotheses[*idx].convex,
        SelectedSurface::Spherical(idx) => !output.stage2.spherical_hypotheses[*idx].convex,
        SelectedSurface::Toroidal(idx) => !output.stage2.toroidal_hypotheses[*idx].convex,
        SelectedSurface::Planar(_) => false,
        // B-spline fits orient the natural normal outward at fit time.
        SelectedSurface::BSpline(_) => false,
    };

    if config.verbose {
        let area = compute_face_area_from_face(&face);
        let n_deg = degenerates.len();
        let has_seam = seam.is_some();
        eprintln!(
            "  Face {fi}: assembled from {} shared edges ({} wires, seam={has_seam}, \
             {n_deg} degenerate), area={area:.4} mm²",
            fd.edge_indices.len(),
            wires.len(),
        );
    }

    Ok(AssembledFace { face, concave })
}

fn polygon_signed_area(poly: &[(f64, f64)]) -> f64 {
    let mut a = 0.0;
    for i in 0..poly.len() {
        let p = poly[i];
        let q = poly[(i + 1) % poly.len()];
        a += p.0 * q.1 - q.0 * p.1;
    }
    0.5 * a
}

/// Build a seam edge along an isoline of the surface, with the given border
/// vertices at parameters (p_lo, p_hi). The vertices' tolerances are bumped
/// to cover their distance from the isoline (the shared boundary curves are
/// only within fit tolerance of this face's surface).
#[allow(clippy::too_many_arguments)]
fn make_seam_edge(
    builder: &b_rep::Builder,
    iso: &geom::Handle_Geom_Curve,
    v_lo: &topo_ds::Vertex,
    v_hi: &topo_ds::Vertex,
    p_lo: f64,
    p_hi: f64,
    pt_lo: &[f64; 3],
    pt_hi: &[f64; 3],
    config: &Config,
) -> Result<OwnedPtr<topo_ds::Edge>, Stage3Error> {
    let c = iso.get();
    for (p, pt, vtx) in [(p_lo, pt_lo, v_lo), (p_hi, pt_hi, v_hi)] {
        let q = c.value(p);
        let d = ((q.x() - pt[0]).powi(2) + (q.y() - pt[1]).powi(2) + (q.z() - pt[2]).powi(2))
            .sqrt();
        let tol = (d * 1.01).max(config.vertex_tolerance_mm);
        builder.update_vertex_vertex_real(vtx, tol);
    }
    let make_edge =
        b_rep_builder_api::MakeEdge::new_handlegeomcurve_vertex2_real2(iso, v_lo, v_hi, p_lo, p_hi);
    if !make_edge.is_done() {
        return Err(Stage3Error::AdjacencyError(format!(
            "MakeEdge failed for seam edge: {:?}",
            make_edge.error(),
        )));
    }
    let mut make_edge = make_edge;
    Ok(make_edge.edge().to_owned())
}

/// Create OCCT faces for all face descriptors (stage 3.4 main entry).
///
/// Shared-edge topological assembly: one TopoDS_Vertex per BRepVertex, one
/// TopoDS_Edge per ReconEdge, pcurves attached to the shared edges for both
/// adjacent faces, seam and degenerate edges as first-class topology, and
/// faces built directly with BRep_Builder (see DEVELOPMENT_PLAN §3.x).
fn create_occt_faces_all(
    mut output: Stage3Output,
    config: &Config,
    viz: Option<&crate::viz::VizSender>,
) -> Result<Stage3Output, Stage3Error> {
    let t = Instant::now();
    let num_faces = output.face_descriptors.len();

    // 1. Create shared TopoDS_Vertex for each BRepVertex.
    // Set tolerance to vertex_tolerance_mm so that mesh vertex positions
    // (which may deviate slightly from fitted analytical surfaces) are accepted
    // by MakeEdge when projecting vertices onto edge curves.
    let builder = b_rep::Builder::new();
    let topo_vertices: Vec<OwnedPtr<topo_ds::Vertex>> = output
        .vertices
        .iter()
        .map(|v| {
            let pt = gp::Pnt::new_real3(v.point[0], v.point[1], v.point[2]);
            let mut mv = b_rep_builder_api::MakeVertex::new_pnt(&pt);
            let vtx = mv.vertex().to_owned();
            builder.update_vertex_vertex_real(&vtx, config.vertex_tolerance_mm);
            vtx
        })
        .collect();

    // 2. Create one shared TopoDS_Edge per ReconEdge.
    let mut topo_edges: Vec<OwnedPtr<topo_ds::Edge>> = Vec::with_capacity(output.edges.len());
    // BRepVertex ids at each edge's curve start/end (stored-curve direction).
    let mut edge_vertices: Vec<[usize; 2]> = Vec::with_capacity(output.edges.len());
    for (ei, edge) in output.edges.iter().enumerate() {
        let curve = match edge.curve_3d.as_ref() {
            Some(c) => c,
            None => {
                // No intersection curve — create a straight-line edge between endpoints
                if edge.vertex_indices[0] != usize::MAX && edge.vertex_indices[1] != usize::MAX {
                    let mut make_edge = b_rep_builder_api::MakeEdge::new_vertex2(
                        &topo_vertices[edge.vertex_indices[0]],
                        &topo_vertices[edge.vertex_indices[1]],
                    );
                    if !make_edge.is_done() {
                        return Err(Stage3Error::AdjacencyError(format!(
                            "Edge {ei}: straight-line fallback failed ({:?}), faces [{}, {}]",
                            make_edge.error(),
                            edge.face_indices[0],
                            edge.face_indices[1]
                        )));
                    }
                    if !config.quiet {
                        eprintln!(
                            "  Edge {ei}: using straight-line fallback (no intersection curve)"
                        );
                    }
                    topo_edges.push(make_edge.edge().to_owned());
                    edge_vertices.push([edge.vertex_indices[0], edge.vertex_indices[1]]);
                } else {
                    return Err(Stage3Error::AdjacencyError(format!(
                        "Edge {ei}: no curve and no vertex endpoints"
                    )));
                }
                continue;
            }
        };
        let (mut make_edge, verts) = if edge.vertex_indices[0] != usize::MAX
            && edge.vertex_indices[1] != usize::MAX
        {
            // Both endpoints have shared vertices — create edge with explicit vertex topology.
            // For periodic curves (circles), MakeEdge strips the TrimmedCurve and uses
            // AdjustPeriodic which forces p1 < p2, always taking the forward arc from V1
            // to V2. We must ensure V1 corresponds to the curve's first_parameter and V2
            // to last_parameter, otherwise MakeEdge will select the complementary arc.
            let c = curve.get();
            let p_start = c.value(c.first_parameter());
            let v0_pt = &output.vertices[edge.vertex_indices[0]].point;
            let v1_pt = &output.vertices[edge.vertex_indices[1]].point;
            let d0_start = (v0_pt[0] - p_start.x()).powi(2)
                + (v0_pt[1] - p_start.y()).powi(2)
                + (v0_pt[2] - p_start.z()).powi(2);
            let d1_start = (v1_pt[0] - p_start.x()).powi(2)
                + (v1_pt[1] - p_start.y()).powi(2)
                + (v1_pt[2] - p_start.z()).powi(2);
            let (vi0, vi1) = if d0_start <= d1_start {
                (edge.vertex_indices[0], edge.vertex_indices[1])
            } else {
                (edge.vertex_indices[1], edge.vertex_indices[0])
            };
            // Use explicit parameter values to avoid OCCT's vertex-to-curve projection,
            // which can fail when vertex tolerance exceeds default precision.
            // vi0 is always the vertex closest to the curve start, vi1 closest to end,
            // so parameters are always (first_parameter, last_parameter).
            let fp = c.first_parameter();
            let lp = c.last_parameter();

            // Update vertex tolerances to accommodate curve-vertex distance.
            // Imprecisely-fitted surfaces produce intersection curves that may not
            // pass exactly through mesh vertex positions.
            let p_fp = c.value(fp);
            let p_lp = c.value(lp);
            let v0_pt = &output.vertices[vi0].point;
            let v1_pt = &output.vertices[vi1].point;
            let d0 = ((v0_pt[0] - p_fp.x()).powi(2)
                + (v0_pt[1] - p_fp.y()).powi(2)
                + (v0_pt[2] - p_fp.z()).powi(2))
            .sqrt();
            let d1 = ((v1_pt[0] - p_lp.x()).powi(2)
                + (v1_pt[1] - p_lp.y()).powi(2)
                + (v1_pt[2] - p_lp.z()).powi(2))
            .sqrt();
            let min_tol = config.vertex_tolerance_mm;
            let tol0 = d0.max(min_tol) * 1.01; // 1% margin
            let tol1 = d1.max(min_tol) * 1.01;
            builder.update_vertex_vertex_real(&topo_vertices[vi0], tol0);
            builder.update_vertex_vertex_real(&topo_vertices[vi1], tol1);

            (
                b_rep_builder_api::MakeEdge::new_handlegeomcurve_vertex2_real2(
                    curve,
                    &topo_vertices[vi0],
                    &topo_vertices[vi1],
                    fp,
                    lp,
                ),
                [vi0, vi1],
            )
        } else {
            // Closed-loop edge (no vertex endpoints) — create without vertices
            (
                b_rep_builder_api::MakeEdge::new_handlegeomcurve(curve),
                [usize::MAX, usize::MAX],
            )
        };
        if !make_edge.is_done() {
            // Degenerate curve geometry (e.g. a micro-edge whose trimmed
            // curve endpoints coincide within Precision::Confusion — MakeEdge
            // then treats the curve as closed and refuses two distinct
            // vertices). Fall back to the straight vertex-to-vertex edge,
            // as for edges with no intersection curve.
            let err = make_edge.error();
            if verts[0] != usize::MAX && verts[1] != usize::MAX {
                let mut line_edge = b_rep_builder_api::MakeEdge::new_vertex2(
                    &topo_vertices[verts[0]],
                    &topo_vertices[verts[1]],
                );
                if line_edge.is_done() {
                    if !config.quiet {
                        eprintln!(
                            "  Edge {ei}: curve-based MakeEdge failed ({err:?}); \
                             using straight-line fallback"
                        );
                    }
                    topo_edges.push(line_edge.edge().to_owned());
                    edge_vertices.push(verts);
                    continue;
                }
            }
            let vtx_detail = if verts[0] != usize::MAX && verts[1] != usize::MAX {
                let p0 = &output.vertices[verts[0]].point;
                let p1 = &output.vertices[verts[1]].point;
                let dist = ((p0[0] - p1[0]).powi(2)
                    + (p0[1] - p1[1]).powi(2)
                    + (p0[2] - p1[2]).powi(2))
                .sqrt();
                format!(
                    "; straight-line fallback also failed, vertices {} \
                     [{:.6},{:.6},{:.6}] and {} [{:.6},{:.6},{:.6}] are {:.3e} mm apart",
                    verts[0], p0[0], p0[1], p0[2],
                    verts[1], p1[0], p1[1], p1[2],
                    dist,
                )
            } else {
                String::new()
            };
            return Err(Stage3Error::AdjacencyError(format!(
                "MakeEdge failed for edge {ei}: {err:?}{vtx_detail}",
            )));
        }
        topo_edges.push(make_edge.edge().to_owned());
        edge_vertices.push(verts);
    }

    let mut ctx = AssemblyCtx {
        topo_vertices,
        topo_edges,
        edge_vertices,
        closed_anchors: HashMap::new(),
    };

    // 3. Scan every face's boundary in UV, then choose seam anchors for the
    // wrapped faces (re-anchoring closed-loop edges so seam endpoints are
    // shared vertices of the shared edges).
    let mut scans: Vec<Option<FaceScan>> = Vec::with_capacity(num_faces);
    for fi in 0..num_faces {
        if output.face_descriptors[fi].edge_indices.is_empty() {
            scans.push(None);
        } else {
            scans.push(Some(scan_face(fi, &output, &ctx, config)?));
        }
    }
    let plans = plan_seam_anchors(&output, &scans, &mut ctx, config)?;
    drop(scans); // re-anchoring invalidated closed-edge samples

    // 4. Assemble each face from the shared edges.
    let mut faces: Vec<OwnedPtr<topo_ds::Face>> = Vec::with_capacity(num_faces);
    let mut concave_faces: Vec<bool> = Vec::with_capacity(num_faces);

    for fi in 0..num_faces {
        let fd = &output.face_descriptors[fi];
        let surface_type = &output.stage2.selected_surfaces[fd.selected_surface_idx];
        let stype = match surface_type {
            SelectedSurface::Planar(_) => "planar",
            SelectedSurface::Conical(i) => {
                if output.stage2.conical_hypotheses[*i].is_cylinder() {
                    "cylindrical"
                } else {
                    "conical"
                }
            }
            SelectedSurface::Spherical(_) => "spherical",
            SelectedSurface::Toroidal(_) => "toroidal",
            SelectedSurface::BSpline(_) => "bspline",
        };

        let (face, is_concave) = if fd.edge_indices.is_empty() {
            // Complete closed surface with no boundary (full sphere or
            // torus): natural bounds. OCCT creates the seam and degenerate
            // edges itself; there are no shared edges to preserve.
            let make_face = b_rep_builder_api::MakeFace::new_handlegeomsurface_real(
                &fd.surface,
                config.vertex_tolerance_mm,
            );
            if !make_face.is_done() {
                return Err(Stage3Error::AdjacencyError(format!(
                    "MakeFace failed for complete closed face {fi}: {:?}",
                    make_face.error(),
                )));
            }
            let concave = match surface_type {
                SelectedSurface::Spherical(idx) => {
                    !output.stage2.spherical_hypotheses[*idx].convex
                }
                SelectedSurface::Toroidal(idx) => {
                    !output.stage2.toroidal_hypotheses[*idx].convex
                }
                _ => false,
            };
            (make_face.face().to_owned(), concave)
        } else {
            let assembled = assemble_face(fi, &output, &ctx, &plans, config)?;
            (assembled.face, assembled.concave)
        };
        concave_faces.push(is_concave);

        if config.verbose {
            eprintln!(
                "  Face {fi}: {stype} — created successfully ({} edges)",
                fd.edge_indices.len(),
            );
        }

        if let Some(viz_sender) = viz {
            let mut overlay = crate::viz::VizOverlay::new();
            // Tessellate the OCCT face and show it translucent
            overlay.shape_meshes.push(crate::viz::tessellate_shape(
                face.as_shape(), [0.3, 0.8, 0.3, 0.4], [0.0, 1.0, 0.0, 1.0], 0.1, 0.5,
            ));
            // Show edge curves in bright green
            for &edge_idx in &output.face_descriptors[fi].edge_indices {
                if let Some(curve) = output.edges[edge_idx].curve_3d.as_ref() {
                    overlay.lines.push(crate::viz::LineOverlay {
                        positions: crate::viz::sample_curve_for_viz(curve, 64),
                        color: [0.0, 1.0, 0.0, 1.0],
                        no_depth_test: true,
                    });
                }
            }
            overlay.status_text = format!(
                "Stage 3.4: Face {fi}/{num_faces}: {stype} ({} edges)",
                output.face_descriptors[fi].edge_indices.len(),
            );
            let centroid = viz_selected_surface_centroid(fi, &output.face_descriptors, &output.stage2);
            overlay.focus_point = Some([centroid[0] as f32, centroid[1] as f32, centroid[2] as f32]);
            overlay.focus_normal = Some(viz_surface_normal_at_point(
                fi, centroid, &output.face_descriptors, &output.stage2,
            ));
            viz_sender.show_and_wait(overlay);
        }

        faces.push(face);
    }

    if !config.quiet {
        eprintln!(
            "Stage 3.4 ({:.3}s): Created {}/{} OCCT faces",
            t.elapsed().as_secs_f64(),
            faces.len(),
            num_faces,
        );
    }

    output.faces = faces;
    output.concave_faces = concave_faces;

    // 5. Validate faces
    validate_faces(&output, config)?;

    // 6. Compare against reference STEP if --compare
    if config.compare_shape.is_some() {
        compare_faces_to_step(&output, config)?;
        compare_surface_orientations_to_step(&output, config);
    }

    Ok(output)
}

// ---------------------------------------------------------------------------
// Stage 3.5: Construct shells
// ---------------------------------------------------------------------------

/// Assemble OCCT faces into shells directly from the shared-edge topology.
///
/// Faces already share their TopoDS_Edge objects (stage 3.4), so a shell is
/// closed by construction: group faces into connected components over the
/// ReconEdge adjacency and add each face to its component's shell with the
/// orientation given by its hypothesis convexity. BRepBuilderAPI_Sewing is
/// no longer used; ShapeAnalysis_Shell remains as validation.
fn construct_shells(
    mut output: Stage3Output,
    config: &Config,
    viz: Option<&crate::viz::VizSender>,
) -> Result<Stage3Output, Stage3Error> {
    let t = Instant::now();
    let num_faces = output.faces.len();
    if num_faces == 0 {
        return Err(Stage3Error::AdjacencyError("no faces to assemble into shells".into()));
    }

    // Connected components of faces over the shared-edge adjacency.
    let mut component: Vec<usize> = vec![usize::MAX; num_faces];
    let mut n_components = 0;
    for start in 0..num_faces {
        if component[start] != usize::MAX {
            continue;
        }
        let comp = n_components;
        n_components += 1;
        let mut queue = vec![start];
        component[start] = comp;
        while let Some(fi) = queue.pop() {
            for &ei in &output.face_descriptors[fi].edge_indices {
                for &nfi in &output.edges[ei].face_indices {
                    if component[nfi] == usize::MAX {
                        component[nfi] = comp;
                        queue.push(nfi);
                    }
                }
            }
        }
    }

    // Build one shell per component. A face whose natural surface normal
    // points into the material (concave) is added Reversed, so the shell is
    // consistently oriented outward by construction (mesh normals point
    // outward by the stage-1 validation).
    let builder = b_rep::Builder::new();
    let mut shells: Vec<OwnedPtr<topo_ds::Shell>> = Vec::new();
    for comp in 0..n_components {
        let mut shell = topo_ds::Shell::new();
        builder.make_shell(&mut shell);
        let mut face_count = 0;
        for (fi, &c) in component.iter().enumerate() {
            if c != comp {
                continue;
            }
            let orient = if output.concave_faces[fi] {
                top_abs::Orientation::Reversed
            } else {
                top_abs::Orientation::Forward
            };
            let oriented = output.faces[fi].as_shape().oriented(orient);
            builder.add(shell.as_shape_mut(), &oriented);
            face_count += 1;
        }
        shell.as_shape_mut().closed_bool(true);
        if config.verbose {
            eprintln!("  Shell {comp}: assembled from {face_count} faces");
        }
        shells.push(shell);
    }

    // Validate shells using ShapeAnalysis_Shell
    // Note: check_oriented_shells returns true if bad edges are FOUND,
    // false if the shell is correctly oriented.
    let mut total_faces_in_shells = 0;
    for (si, shell) in shells.iter().enumerate() {
        let mut sa = shape_analysis::Shell::new();
        let has_bad_edges = sa.check_oriented_shells(shell.as_shape(), true, false);

        // Count faces in this shell
        let mut face_count = 0;
        let mut face_explorer = top_exp::Explorer::new_shape_shapeenum2(
            shell.as_shape(),
            top_abs::ShapeEnum::Face,
            top_abs::ShapeEnum::Shape,
        );
        while face_explorer.more() {
            face_count += 1;
            face_explorer.next();
        }
        total_faces_in_shells += face_count;

        if config.verbose {
            let has_free = sa.has_free_edges();
            eprintln!(
                "  Shell {si}: {face_count} faces, bad_edges={has_bad_edges}, free_edges={has_free}"
            );
        }

        if has_bad_edges {
            eprintln!("  Warning: shell {si} has orientation inconsistencies");
        }
    }

    if let Some(viz_sender) = viz {
        for (si, shell) in shells.iter().enumerate() {
            let mut overlay = crate::viz::VizOverlay::new();
            overlay.shape_meshes.push(crate::viz::tessellate_shape(
                shell.as_shape(), [0.3, 0.8, 0.3, 0.5], [0.0, 1.0, 0.0, 1.0], 0.1, 0.5,
            ));
            overlay.status_text = format!(
                "Stage 3.5: Shell {si}/{} ({total_faces_in_shells} faces)",
                shells.len(),
            );
            viz_sender.show_and_wait(overlay);
        }
    }


    if !config.quiet {
        eprintln!(
            "Stage 3.5 ({:.3}s): Constructed {} shell(s) from {} faces ({} faces in shells)",
            t.elapsed().as_secs_f64(),
            shells.len(),
            num_faces,
            total_faces_in_shells,
        );
    }

    // Compare against reference STEP if --compare
    if config.compare_shape.is_some() {
        compare_shells_to_step(&shells, config)?;
    }

    output.shells = shells;
    Ok(output)
}

/// Compare stage 3.5 shells against reference STEP shape.
fn compare_shells_to_step(
    shells: &[OwnedPtr<topo_ds::Shell>],
    config: &Config,
) -> Result<(), Stage3CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    // Count shells in reference STEP
    let mut step_shell_count = 0;
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        compare_shape,
        top_abs::ShapeEnum::Shell,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        step_shell_count += 1;
        explorer.next();
    }

    if !config.quiet {
        eprintln!(
            "  Compare 3.5: {} shell(s) created, {} shell(s) in STEP reference",
            shells.len(),
            step_shell_count,
        );
    }

    // Check that each shell is closed (no free edges)
    for (si, shell) in shells.iter().enumerate() {
        let mut sa = shape_analysis::Shell::new();
        sa.check_oriented_shells(shell.as_shape(), true, false);
        if sa.has_free_edges() {
            if !config.quiet {
                eprintln!("  Compare 3.5: shell {si} has free (unclosed) edges");
            }
            return Err(Stage3CompareError {
                substage: 5,
                check_type: "shell",
                element_index: si,
                max_distance: f64::INFINITY,
                tolerance: config.vertex_tolerance_mm,
            });
        }
    }

    // Verify orientation (check_oriented_shells returns true if bad edges found)
    for (si, shell) in shells.iter().enumerate() {
        let mut sa = shape_analysis::Shell::new();
        let has_bad_edges = sa.check_oriented_shells(shell.as_shape(), false, false);
        if has_bad_edges && !config.quiet {
            eprintln!("  Compare 3.5: shell {si} has orientation inconsistencies");
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Stage 3.6: Construct solids from shells
// ---------------------------------------------------------------------------

fn construct_solids(
    mut output: Stage3Output,
    config: &Config,
    viz: Option<&crate::viz::VizSender>,
) -> Result<Stage3Output, Stage3Error> {
    let t = Instant::now();
    if output.shells.is_empty() {
        return Err(Stage3Error::AdjacencyError("no shells to make solids from".into()));
    }

    let mut solids: Vec<OwnedPtr<topo_ds::Solid>> = Vec::new();

    for (si, shell) in output.shells.iter().enumerate() {
        // Use ShapeFix_Solid::SolidFromShell which handles orientation automatically
        let mut fixer = shape_fix::Solid::new();
        let mut solid = fixer.solid_from_shell(shell);

        // Apply ShapeFix_Shape to comprehensively fix the solid:
        // fixes faces, wires, edges, shells, and solid orientation.
        let mut shape_fixer = shape_fix::Shape::new_shape(solid.as_shape());
        shape_fixer.set_precision(config.vertex_tolerance_mm);
        let progress = message::ProgressRange::new();
        shape_fixer.perform(&progress);
        let fixed_shape = shape_fixer.shape();

        // Extract the solid from the fixed shape
        let solid_exp = top_exp::Explorer::new_shape_shapeenum2(
            &fixed_shape, top_abs::ShapeEnum::Solid, top_abs::ShapeEnum::Shape,
        );
        if solid_exp.more() {
            solid = topo_ds::solid(solid_exp.value()).to_owned();
        }

        // Fix SameParameter consistency: ensure pcurves and 3D curves agree,
        // and update edge/vertex tolerances to accommodate any gaps.
        // Use forced=false to only recompute edges that are not already flagged
        // as SameParameter, and use the model's vertex tolerance rather than
        // an aggressive tolerance that can corrupt pcurves on closely-spaced edges.
        b_rep_lib::same_parameter_shape_real_bool(solid.as_shape(), config.vertex_tolerance_mm, false);
        b_rep_lib::update_tolerances_shape_bool(solid.as_shape(), true);

        // Orient the solid's faces so all normals point outward
        let oriented = b_rep_lib::orient_closed_solid(&mut solid);
        if config.verbose {
            eprintln!("  Solid {si}: orient_closed_solid returned {oriented}");
        }

        // Validate with BRepCheck_Analyzer
        let analyzer = b_rep_check::Analyzer::new_shape_bool(solid.as_shape(), true);
        if !analyzer.is_valid() {
            eprintln!("  Warning: solid {si} failed BRepCheck validation");
            if config.verbosity >= 2 {
                for shape_type in &[
                    top_abs::ShapeEnum::Face,
                    top_abs::ShapeEnum::Wire,
                    top_abs::ShapeEnum::Edge,
                    top_abs::ShapeEnum::Shell,
                ] {
                    let mut exp = top_exp::Explorer::new_shape_shapeenum2(
                        solid.as_shape(), *shape_type, top_abs::ShapeEnum::Shape,
                    );
                    let mut idx = 0;
                    while exp.more() {
                        if !analyzer.is_valid_shape(exp.value()) {
                            eprintln!("    BRepCheck fail: {shape_type:?} {idx}");
                        }
                        idx += 1;
                        exp.next();
                    }
                }
            }
        }

        // Compute volume for reporting
        let mut gprops = g_prop::GProps::new();
        b_rep_g_prop::volume_properties_shape_gprops_bool3(
            solid.as_shape(),
            &mut gprops,
            true,  // OnlyClosed
            false, // SkipShared
            false, // UseTriangulation
        );
        let volume = gprops.mass();

        if config.verbose {
            // Also compute with OnlyClosed=false for comparison
            let mut gprops2 = g_prop::GProps::new();
            b_rep_g_prop::volume_properties_shape_gprops_bool3(
                solid.as_shape(),
                &mut gprops2,
                false, // OnlyClosed=false
                false,
                false,
            );
            let volume2 = gprops2.mass();
            eprintln!("  Solid {si}: volume = {volume:.6} mm\u{00b3} (open={volume2:.6})");
        }

        solids.push(solid);
    }

    // Merge internal-void solids into their containers. A shell whose faces
    // enclose a cavity gets turned into an ordinary solid by SolidFromShell
    // above; if its bounding box is contained in a larger solid's box it is
    // actually a void of that solid, and the STEP output must be ONE solid
    // whose volume excludes the cavity - not two nested solids.
    if solids.len() > 1 {
        solids = merge_void_solids(solids, config);
    }

    if let Some(viz_sender) = viz {
        for (si, solid) in solids.iter().enumerate() {
            let mut overlay = crate::viz::VizOverlay::new();
            overlay.shape_meshes.push(crate::viz::tessellate_shape(
                solid.as_shape(), [0.3, 0.8, 0.3, 0.5], [0.0, 1.0, 0.0, 1.0], 0.1, 0.5,
            ));
            overlay.status_text = format!(
                "Stage 3.6: Solid {si}/{}",
                solids.len(),
            );
            viz_sender.show_and_wait(overlay);
        }
    }


    if !config.quiet {
        eprintln!(
            "Stage 3.6 ({:.3}s): Constructed {} solid(s) from {} shell(s)",
            t.elapsed().as_secs_f64(),
            solids.len(),
            output.shells.len(),
        );
    }

    // Compare against reference STEP if --compare
    if config.compare_shape.is_some() {
        compare_solids_to_step(&solids, config)?;
    }

    output.solids = solids;
    Ok(output)
}

/// Volume of a solid via BRepGProp.
fn solid_volume(solid: &topo_ds::Solid) -> f64 {
    let mut gprops = g_prop::GProps::new();
    b_rep_g_prop::volume_properties_shape_gprops_bool3(
        solid.as_shape(), &mut gprops, true, false, false,
    );
    gprops.mass()
}

/// Merge void solids (bounding box contained in a larger solid) into their
/// containers by adding their shells to the container solid. Shell
/// orientation is verified by the resulting volume: container minus voids.
fn merge_void_solids(
    solids: Vec<OwnedPtr<topo_ds::Solid>>,
    config: &Config,
) -> Vec<OwnedPtr<topo_ds::Solid>> {
    let n = solids.len();
    let mut volumes = Vec::with_capacity(n);
    let mut boxes: Vec<([f64; 3], [f64; 3])> = Vec::with_capacity(n);
    for solid in &solids {
        volumes.push(solid_volume(solid));
        let mut bb = bnd::Box::new();
        b_rep_bnd_lib::add_shape_box_bool(solid.as_shape(), &mut bb, false);
        let cmin = bb.corner_min();
        let cmax = bb.corner_max();
        boxes.push((
            [cmin.x(), cmin.y(), cmin.z()],
            [cmax.x(), cmax.y(), cmax.z()],
        ));
    }

    // For each solid, find the smallest strictly-larger solid whose bounding
    // box contains it: that is its container.
    let mut container: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let mut best: Option<usize> = None;
        for j in 0..n {
            if i == j || volumes[j] <= volumes[i] {
                continue;
            }
            let contains = (0..3).all(|k| {
                boxes[j].0[k] <= boxes[i].0[k] && boxes[i].1[k] <= boxes[j].1[k]
            });
            if contains && best.is_none_or(|b| volumes[j] < volumes[b]) {
                best = Some(j);
            }
        }
        container[i] = best;
    }

    if container.iter().all(|c| c.is_none()) {
        return solids;
    }

    let mut result: Vec<OwnedPtr<topo_ds::Solid>> = Vec::new();
    for j in 0..n {
        if container[j].is_some() {
            continue; // voids are merged into their containers
        }
        let void_indices: Vec<usize> = (0..n).filter(|&i| container[i] == Some(j)).collect();
        if void_indices.is_empty() {
            result.push(solids[j].to_owned());
            continue;
        }

        let expected: f64 = volumes[j] - void_indices.iter().map(|&i| volumes[i]).sum::<f64>();
        let builder = topo_ds::Builder::new();

        // Try both orientations for the void shells; keep the one whose
        // volume matches container minus voids.
        let mut merged: Option<OwnedPtr<topo_ds::Solid>> = None;
        for reverse in [true, false] {
            let mut candidate = solids[j].to_owned();
            for &i in &void_indices {
                let mut shell_exp = top_exp::Explorer::new_shape_shapeenum2(
                    solids[i].as_shape(), top_abs::ShapeEnum::Shell, top_abs::ShapeEnum::Shape,
                );
                while shell_exp.more() {
                    let shell_shape = shell_exp.value();
                    if reverse {
                        builder.add(candidate.as_shape_mut(), &shell_shape.reversed());
                    } else {
                        builder.add(candidate.as_shape_mut(), shell_shape);
                    }
                    shell_exp.next();
                }
            }
            let vol = solid_volume(&candidate);
            if (vol - expected).abs() <= expected.abs() * 1e-6 + 1e-9 {
                merged = Some(candidate);
                break;
            }
        }

        match merged {
            Some(m) => {
                if !config.quiet {
                    eprintln!(
                        "  Solid {j}: merged {} internal void shell(s), volume {:.6} mm3",
                        void_indices.len(),
                        solid_volume(&m),
                    );
                }
                result.push(m);
            }
            None => {
                eprintln!(
                    "  Warning: could not merge void shell(s) into solid {j} \
(volume mismatch); keeping solids separate"
                );
                result.push(solids[j].to_owned());
                for &i in &void_indices {
                    result.push(solids[i].to_owned());
                }
            }
        }
    }
    result
}

/// Compare stage 3.6 solids against reference STEP shape.
fn compare_solids_to_step(
    solids: &[OwnedPtr<topo_ds::Solid>],
    config: &Config,
) -> Result<(), Stage3CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    // Count solids in reference STEP
    let mut step_solids: Vec<OwnedPtr<topo_ds::Solid>> = Vec::new();
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        compare_shape,
        top_abs::ShapeEnum::Solid,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        let s = topo_ds::solid(explorer.value());
        step_solids.push(s.to_owned());
        explorer.next();
    }

    if !config.quiet {
        eprintln!(
            "  Compare 3.6: {} solid(s) created, {} solid(s) in STEP reference",
            solids.len(),
            step_solids.len(),
        );
    }

    // Compute volumes of constructed and reference solids
    let mut our_volumes: Vec<f64> = Vec::new();
    for solid in solids.iter() {
        let mut gprops = g_prop::GProps::new();
        b_rep_g_prop::volume_properties_shape_gprops_bool3(
            solid.as_shape(),
            &mut gprops,
            true, false, false,
        );
        our_volumes.push(gprops.mass().abs());
    }

    let mut step_volumes: Vec<f64> = Vec::new();
    for solid in step_solids.iter() {
        let mut gprops = g_prop::GProps::new();
        b_rep_g_prop::volume_properties_shape_gprops_bool3(
            solid.as_shape(),
            &mut gprops,
            true, false, false,
        );
        step_volumes.push(gprops.mass().abs());
    }

    // Match solids by volume (sort both lists and compare pairwise)
    let mut our_sorted: Vec<(usize, f64)> = our_volumes.iter().copied().enumerate().collect();
    let mut step_sorted: Vec<(usize, f64)> = step_volumes.iter().copied().enumerate().collect();
    our_sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    step_sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // Check volume agreement for each pair
    let n_pairs = our_sorted.len().min(step_sorted.len());
    for i in 0..n_pairs {
        let (our_idx, our_vol) = our_sorted[i];
        let (_step_idx, step_vol) = step_sorted[i];
        let rel_diff = if step_vol > 0.0 {
            (our_vol - step_vol).abs() / step_vol
        } else {
            0.0
        };

        if !config.quiet {
            eprintln!(
                "  Compare 3.6: solid {our_idx} volume {our_vol:.6} mm\u{00b3} vs STEP {step_vol:.6} mm\u{00b3} (rel diff {rel_diff:.2e})",
            );
        }

        if rel_diff > 0.01 {
            eprintln!(
                "  Compare 3.6: ERROR solid {our_idx} volume differs by {rel_diff:.2e} (>1%)",
            );
            return Err(Stage3CompareError {
                substage: 6,
                check_type: "solid volume",
                element_index: our_idx,
                max_distance: rel_diff,
                tolerance: 0.01,
            });
        }
    }

    // Compute BRepExtrema distance between our solids and STEP shape
    for (si, solid) in solids.iter().enumerate() {
        let progress = message::ProgressRange::new();
        let dist_calc = b_rep_extrema::DistShapeShape::new_shape2_extflag_extalgo_progressrange(
            solid.as_shape(),
            compare_shape,
            0,
            extrema::ExtAlgo::Grad,
            &progress,
        );
        if dist_calc.is_done() && dist_calc.nb_solution() > 0 {
            let dist = dist_calc.value();
            if !config.quiet {
                eprintln!("  Compare 3.6: solid {si} distance to STEP = {dist:.6e} mm");
            }
            if dist > config.surface_tolerance_mm {
                return Err(Stage3CompareError {
                    substage: 6,
                    check_type: "solid",
                    element_index: si,
                    max_distance: dist,
                    tolerance: config.surface_tolerance_mm,
                });
            }
        }
    }

    Ok(())
}

// Stage 3 entry point
// ---------------------------------------------------------------------------

/// Run stage 3: reconstruct B-Rep surfaces, edges, and topology from fitted surfaces.
pub fn stage3(config: &Config, input: Stage2Output, viz: Option<&crate::viz::VizSender>) -> Result<Stage3Output, Stage3Error> {
    // Stage 3.1: Create OCCT surface objects and build adjacency graph
    let output = build_surfaces_and_adjacency(input, config)?;

    // Stage 3.2: Detect tangency relationships along edges. This runs before
    // the stage 3.1 --compare check because boundary vertices on TANGENT
    // edges are checked against the chordal surface tolerance: when the
    // tessellation lacks a precise vertex ring on a tangent junction, the
    // region boundary can sit one row off the junction circle — still within
    // chordal tolerance of both surfaces and harmless to reconstruction
    // (tangent junction curves are built from the fitted hypotheses, not
    // from these vertices).
    let mut output = detect_tangency(output, config);

    // Refine sharp corner vertex positions onto the intersection of their
    // adjacent fitted surfaces (mesh corner vertices can sit a fraction of a
    // facet off the CAD corner, e.g. where a fillet terminates).
    {
        let Stage3Output { stage2, face_descriptors, edges, vertices, .. } = &mut output;
        refine_corner_vertices(vertices, edges, face_descriptors, stage2, config);
    }

    if config.compare_shape.is_some() {
        compare_adjacency_to_step(&output.edges, &output.vertices, &output.stage2, config)?;
    }

    if !config.stage.at_least(3, 2) {
        return Ok(output);
    }


    if !config.stage.at_least(3, 3) {
        return Ok(output);
    }

    // Stage 3.3: Compute edge curves via surface-surface intersection
    let viz_33 = if config.viz_active(3, 3) { viz } else { None };
    let output = compute_edge_curves_all(output, config, viz_33)?;

    if !config.stage.at_least(3, 4) {
        return Ok(output);
    }

    // Stage 3.4: Create OCCT faces
    let viz_34 = if config.viz_active(3, 4) { viz } else { None };
    let output = create_occt_faces_all(output, config, viz_34)?;

    if !config.stage.at_least(3, 5) {
        return Ok(output);
    }

    // Stage 3.5: Construct shells
    let viz_35 = if config.viz_active(3, 5) { viz } else { None };
    let output = construct_shells(output, config, viz_35)?;

    if !config.stage.at_least(3, 6) {
        return Ok(output);
    }

    // Stage 3.6: Construct solids
    let viz_36 = if config.viz_active(3, 6) { viz } else { None };
    let output = construct_solids(output, config, viz_36)?;

    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression for the 3DBenchy stage-3 crash: the chain repair's B-spline
    /// candidate is bounded and non-periodic, but its trim parameters went
    /// through the modulo-2π arc selection, which can return parameters past
    /// the curve's domain — Geom_TrimmedCurve then panics with "parameters
    /// out of range". With the vertex projections in reverse parameter order
    /// and no chain samples to overrule the forward default, the old code
    /// returned t_hi ≈ t1 + 2π on a domain spanning ~1.
    #[test]
    fn bounded_bspline_trim_parameters_stay_in_domain() {
        let points = [
            [0.0, 0.0, 0.0],
            [1.0, 0.4, 0.0],
            [2.0, 0.5, 0.0],
            [3.0, 0.4, 0.0],
            [4.0, 0.0, 0.0],
        ];
        let mut pts = n_collection::Array1_gp_Pnt::new_int2(1, points.len() as i32);
        for (i, p) in points.iter().enumerate() {
            let pnt = gp::Pnt::new_real3(p[0], p[1], p[2]);
            pts.set_value(i as i32 + 1, &pnt);
        }
        let fit = geom_api::PointsToBSpline::new_arr1_pnt_int2_shape_real(
            &pts,
            3,
            8,
            geom_abs::Shape::C2,
            1e-6,
        );
        assert!(fit.is_done());
        let curve = fit.curve().to_handle_curve();
        let (fp, lp) = {
            let c = curve.get();
            assert!(!c.is_periodic());
            (c.first_parameter(), c.last_parameter())
        };

        let mesh = ConnectedMesh::default();
        let t0 = fp + 0.9 * (lp - fp);
        let t1 = fp + 0.1 * (lp - fp);
        let (t_lo, t_hi) = candidate_trim_parameters(t0, t1, &[], &curve, &mesh)
            .expect("distinct projections give a non-degenerate span");
        assert!(
            t_lo >= fp - 1e-12 && t_hi <= lp + 1e-12 && t_lo < t_hi,
            "trim [{t_lo}, {t_hi}] escapes domain [{fp}, {lp}]"
        );

        // The construction that panicked before the fix.
        let arc = geom::TrimmedCurve::new_handlegeomcurve_real2(&curve, t_lo, t_hi);
        let tc = geom::TrimmedCurve::to_handle(arc).to_handle_curve();
        let c = tc.get();
        assert!((c.first_parameter() - t_lo).abs() < 1e-9);
        assert!((c.last_parameter() - t_hi).abs() < 1e-9);
    }

    /// The periodic branch keeps the chain-supported arc selection: a circle
    /// arc crossing the parameter seam comes back with t_hi past 2π, which
    /// Geom_TrimmedCurve normalizes against the periodic basis.
    #[test]
    fn periodic_circle_seam_crossing_arc_kept() {
        let center = gp::Pnt::new_real3(0.0, 0.0, 0.0);
        let normal = gp::Dir::new_real3(0.0, 0.0, 1.0);
        let ax2 = gp::Ax2::new_pnt_dir(&center, &normal);
        let circle = geom::Circle::new_ax2_real(&ax2, 10.0);
        let curve = geom::Circle::to_handle(circle).to_handle_curve();

        // Chain vertices on the short arc crossing the seam, generated from
        // the curve itself so their projected parameters are frame-independent.
        let mut mesh = ConnectedMesh::default();
        for t in [5.5, 5.9, 0.1, 0.5] {
            let p = curve.get().value(t);
            mesh.vertices.push(MeshVertex::from_xyz(p.x(), p.y(), p.z()));
        }
        let chain: Vec<usize> = (0..mesh.vertices.len()).collect();

        let (t_lo, t_hi) = candidate_trim_parameters(5.3, 0.7, &chain, &curve, &mesh)
            .expect("distinct projections give a non-degenerate span");
        let tau = std::f64::consts::TAU;
        assert!(
            (t_lo - 5.3).abs() < 1e-9 && (t_hi - (0.7 + tau)).abs() < 1e-9,
            "expected the seam-crossing arc [5.3, 0.7+2π], got [{t_lo}, {t_hi}]"
        );

        let arc = geom::TrimmedCurve::new_handlegeomcurve_real2(&curve, t_lo, t_hi);
        let tc = geom::TrimmedCurve::to_handle(arc).to_handle_curve();
        let c = tc.get();
        let mid = c.value((c.first_parameter() + c.last_parameter()) / 2.0);
        let expected = curve.get().value((5.3 + 0.7 + tau) / 2.0);
        let dist = ((mid.x() - expected.x()).powi(2)
            + (mid.y() - expected.y()).powi(2)
            + (mid.z() - expected.z()).powi(2))
        .sqrt();
        assert!(dist < 1e-9, "arc midpoint off by {dist}");
    }

    /// Regression for the second 3DBenchy stage-3 crash: sampling a
    /// degenerate IntSS line at its infinite parameter bounds hands
    /// reconstruct_full_circle points with ~1e100 coordinates. The cross
    /// product's squared norm overflows to infinity, slipping past the
    /// collinearity guard, and n/|n| = finite/inf = 0 made gp_Dir panic
    /// with "input vector has zero norm". Such points must yield None.
    #[test]
    fn full_circle_from_overflowing_coordinates_rejected() {
        let p0 = gp::Pnt::new_real3(1e100, 0.0, 0.0);
        let p1 = gp::Pnt::new_real3(0.0, 1e100, 0.0);
        let p2 = gp::Pnt::new_real3(0.0, 0.0, 1e100);
        assert!(reconstruct_full_circle(&p0, &p1, &p2).is_none());

        // Nearly-collinear samples from a line at the same scale (the actual
        // 3DBenchy shape of the failure).
        let q0 = gp::Pnt::new_real3(-1e100, -1e100, 0.0);
        let q1 = gp::Pnt::new_real3(-5e99, -5e99, 1.0);
        let q2 = gp::Pnt::new_real3(5e99, 5e99, -1.0);
        assert!(reconstruct_full_circle(&q0, &q1, &q2).is_none());

        // Sane coordinates still reconstruct a circle.
        let c0 = gp::Pnt::new_real3(10.0, 0.0, 5.0);
        let c1 = gp::Pnt::new_real3(0.0, 10.0, 5.0);
        let c2 = gp::Pnt::new_real3(-10.0, 0.0, 5.0);
        let circle = reconstruct_full_circle(&c0, &c1, &c2).expect("valid circle");
        let c = circle.get();
        assert!(c.is_periodic());
    }

    /// Coincident projections on a bounded candidate are rejected rather than
    /// passed to Geom_TrimmedCurve (which refuses U1 == U2).
    #[test]
    fn bounded_degenerate_span_rejected() {
        let p0 = gp::Pnt::new_real3(0.0, 0.0, 0.0);
        let dir = gp::Dir::new_real3(1.0, 0.0, 0.0);
        let line = geom::Line::new_pnt_dir(&p0, &dir);
        let curve = geom::Line::to_handle(line).to_handle_curve();
        let mesh = ConnectedMesh::default();
        assert!(candidate_trim_parameters(2.0, 2.0, &[], &curve, &mesh).is_none());
    }
}
