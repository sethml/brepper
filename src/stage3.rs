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
    b_rep, b_rep_adaptor, b_rep_builder_api, b_rep_check, b_rep_extrema, b_rep_g_prop, b_rep_lib,
    b_rep_tools, extrema,
    geom, geom2d, geom2d_api, geom_api, g_prop, gp, message, shape_analysis,
    shape_fix,
    t_col_std, t_colgp, top_abs, top_exp, top_loc, topo_ds, OwnedPtr,
};
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
    pub surface: OwnedPtr<geom::HandleGeomSurface>,
    // TODO: OCCT face object (TopoDS_Face) — populated in stage 3.4
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
    pub curve_3d: Option<OwnedPtr<geom::HandleGeomCurve>>,
    // TODO: pcurves on each face (Geom2d_Curve) — populated in stage 3.3
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
    // TODO: UV coordinates on each adjacent face's surface — populated in stage 3.3
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
    /// Each MakeFace owns the resulting TopoDS_Face; access via `.face()`.
    pub make_faces: Vec<OwnedPtr<b_rep_builder_api::MakeFace>>,
    /// Indices of faces that are concave (require reversed orientation).
    pub concave_faces: Vec<bool>,
    /// OCCT shell(s) created in stage 3.5 by sewing faces.
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
            .field("make_faces", &self.make_faces.len())
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
    NotImplemented(String),
    AdjacencyError(String),
    EdgeCurveError(String),
    Compare(Stage3CompareError),
}

impl Display for Stage3Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage3Error::NotImplemented(stage) => write!(f, "stage {stage} not yet implemented"),
            Stage3Error::AdjacencyError(msg) => write!(f, "stage 3.1 adjacency error: {msg}"),
            Stage3Error::EdgeCurveError(msg) => write!(f, "stage 3.3 edge curve error: {msg}"),
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
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage3CompareError> {
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

        for &si in &sample_indices {
            let vi = boundary[si];
            let v = &mesh.vertices[vi];
            let pt = gp::Pnt::new_real3(v.x, v.y, v.z);
            let d = min_distance_to_shape(&pt, edge_compound.as_shape());
            max_dist = max_dist.max(d);
        }

        edge_overall_max = edge_overall_max.max(max_dist);
        if max_dist > config.vertex_tolerance_mm {
            return Err(Stage3CompareError {
                substage: 1,
                check_type: "edge",
                element_index: ei,
                max_distance: max_dist,
                tolerance: config.vertex_tolerance_mm,
            });
        }
    }

    // Check BRepVertices against STEP vertices
    let mut vert_overall_max = 0.0_f64;
    for (vi, vert) in vertices.iter().enumerate() {
        let pt = gp::Pnt::new_real3(vert.point[0], vert.point[1], vert.point[2]);
        let d = min_distance_to_shape(&pt, vertex_compound.as_shape());

        vert_overall_max = vert_overall_max.max(d);
        if d > config.vertex_tolerance_mm {
            return Err(Stage3CompareError {
                substage: 1,
                check_type: "vertex",
                element_index: vi,
                max_distance: d,
                tolerance: config.vertex_tolerance_mm,
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
) -> OwnedPtr<geom::HandleGeomSurface> {
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
        SelectedSurface::Cylindrical(idx) => {
            let hyp = &output.cylindrical_hypotheses[*idx];
            let origin = gp::Pnt::new_real3(
                hyp.axis_origin[0],
                hyp.axis_origin[1],
                hyp.axis_origin[2],
            );
            let dir = gp::Dir::new_real3(
                hyp.axis_direction[0],
                hyp.axis_direction[1],
                hyp.axis_direction[2],
            );
            let ax3 = gp::Ax3::new_pnt_dir(&origin, &dir);
            let cyl = geom::CylindricalSurface::new_ax3_real(&ax3, hyp.radius);
            geom::CylindricalSurface::to_handle(cyl).to_handle_surface()
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
        SelectedSurface::Conical(idx) => {
            let hyp = &output.conical_hypotheses[*idx];
            let origin = gp::Pnt::new_real3(hyp.apex[0], hyp.apex[1], hyp.apex[2]);
            let dir = gp::Dir::new_real3(
                hyp.axis_direction[0],
                hyp.axis_direction[1],
                hyp.axis_direction[2],
            );
            let ax3 = gp::Ax3::new_pnt_dir(&origin, &dir);
            // ConicalSurface::new_ax3_real2(ax3, semi_angle, ref_radius)
            // With origin at apex, ref_radius = 0.0
            let cone = geom::ConicalSurface::new_ax3_real2(&ax3, hyp.half_angle, 0.0);
            geom::ConicalSurface::to_handle(cone).to_handle_surface()
        }
        SelectedSurface::Toroidal(_idx) => {
            todo!("ToroidalSurface OCCT binding not yet available")
        }
    }
}

/// Get the faces belonging to a selected surface.
fn surface_faces<'a>(surface: &SelectedSurface, output: &'a Stage2Output) -> &'a [usize] {
    match surface {
        SelectedSurface::Planar(idx) => &output.planar_hypotheses[*idx].faces,
        SelectedSurface::Cylindrical(idx) => &output.cylindrical_hypotheses[*idx].faces,
        SelectedSurface::Spherical(idx) => &output.spherical_hypotheses[*idx].faces,
        SelectedSurface::Conical(idx) => &output.conical_hypotheses[*idx].faces,
        SelectedSurface::Toroidal(idx) => &output.toroidal_hypotheses[*idx].faces,
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
                    SelectedSurface::Planar(_) => "planar",
                    SelectedSurface::Cylindrical(_) => "cylindrical",
                    SelectedSurface::Spherical(_) => "spherical",
                    SelectedSurface::Conical(_) => "conical",
                    SelectedSurface::Toroidal(_) => "toroidal",
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
    // Compare against STEP edges/vertices if --compare was specified
    if config.compare_shape.is_some() {
        compare_adjacency_to_step(&recon_edges, &brep_vertices, &input.mesh, config)?;
    }

    Ok(Stage3Output {
        stage2: input,
        face_descriptors,
        edges: recon_edges,
        vertices: brep_vertices,
        make_faces: Vec::new(), // populated in stage 3.4
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
        SelectedSurface::Cylindrical(idx) => {
            let hyp = &stage2.cylindrical_hypotheses[*idx];
            // Project point onto the axis to find the closest axis point
            let ax = hyp.axis_direction;
            let ao = hyp.axis_origin;
            let dp = [
                point[0] - ao[0],
                point[1] - ao[1],
                point[2] - ao[2],
            ];
            let t = dp[0] * ax[0] + dp[1] * ax[1] + dp[2] * ax[2];
            let radial = [
                dp[0] - t * ax[0],
                dp[1] - t * ax[1],
                dp[2] - t * ax[2],
            ];
            let len = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            if len < 1e-15 {
                return None; // Point is on the axis
            }
            let mut n = [radial[0] / len, radial[1] / len, radial[2] / len];
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
            // Outward normal on a cone: radial_unit * cos(half_angle) - axis * sin(half_angle)
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

    // Compute surface-surface intersection
    let int_ss = geom_api::IntSS::new_handlegeomsurface2_real(
        surf0, surf1, config.vertex_tolerance_mm,
    );

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
        1 // 1-indexed
    } else {
        select_closest_curve(&int_ss, &edge.mesh_boundary_vertices, mesh)
    };

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
    let fallback_curve: Option<OwnedPtr<geom::HandleGeomCurve>> = if !curve_valid
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
    let full_circle_curve: Option<OwnedPtr<geom::HandleGeomCurve>> =
        if edge.vertex_indices[0] == usize::MAX && edge.vertex_indices[1] == usize::MAX {
            let c = initial_curve_handle.get();
            let span = c.last_parameter() - c.first_parameter();
            if (span - 2.0 * std::f64::consts::PI).abs() > 0.1 {
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

    let curve_handle: &geom::HandleGeomCurve = match full_circle_curve.as_ref() {
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
    // Create trimmed curve
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

    edge.curve_3d = Some(final_curve);
    Ok(())
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
    let mut best_idx = 1_i32;
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

        // Sample 5 points along the curve and compute min distance to centroid
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
) -> Option<OwnedPtr<geom::HandleGeomCurve>> {
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
    if n2 < 1e-30 {
        return None; // points are collinear
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

    let n_len = n2.sqrt();
    let center = gp::Pnt::new_real3(cx, cy, cz);
    let normal = gp::Dir::new_real3(n[0] / n_len, n[1] / n_len, n[2] / n_len);
    let ax2 = gp::Ax2::new_pnt_dir(&center, &normal);
    let circle = geom::Circle::new_ax2_real(&ax2, radius);
    Some(geom::Circle::to_handle(circle).to_handle_curve())
}

fn project_point_on_curve(
    point: &[f64; 3],
    curve: &geom::HandleGeomCurve,
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
        (SelectedSurface::Spherical(_), SelectedSurface::Cylindrical(_))
            | (SelectedSurface::Cylindrical(_), SelectedSurface::Spherical(_))
    );

    let is_closed_loop = edge.vertex_indices[0] == usize::MAX || edge.vertex_indices[1] == usize::MAX;

    if has_sphere_cylinder {
        if is_closed_loop {
            return compute_tangent_edge_curve_sphere_cylinder_closed(
                edge, ss0, ss1, stage2, config,
            );
        }
        let v0 = &vertices[edge.vertex_indices[0]];
        let v1 = &vertices[edge.vertex_indices[1]];
        return compute_tangent_edge_curve_sphere_cylinder(
            edge, v0, v1, ss0, ss1, stage2, config,
        );
    }

    // Non-sphere-cylinder tangent edges require vertex endpoints
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
    // Get the line direction from the cylinder axis (if one of the surfaces is cylindrical)
    let line_dir = match (ss0, ss1) {
        (SelectedSurface::Cylindrical(idx), SelectedSurface::Planar(_))
        | (SelectedSurface::Planar(_), SelectedSurface::Cylindrical(idx)) => {
            let hyp = &stage2.cylindrical_hypotheses[*idx];
            hyp.axis_direction
        }
        (SelectedSurface::Cylindrical(idx), SelectedSurface::Cylindrical(_)) => {
            let hyp = &stage2.cylindrical_hypotheses[*idx];
            hyp.axis_direction
        }
        _ => {
            // Fallback: use direction from vertex to vertex
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
        (SelectedSurface::Cylindrical(cyl_idx), SelectedSurface::Planar(plane_idx))
        | (SelectedSurface::Planar(plane_idx), SelectedSurface::Cylindrical(cyl_idx)) => {
            let cyl = &stage2.cylindrical_hypotheses[*cyl_idx];
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
                v0.point[0] - cyl.axis_origin[0],
                v0.point[1] - cyl.axis_origin[1],
                v0.point[2] - cyl.axis_origin[2],
            ];
            let t_axis = v0_minus_o[0] * a[0] + v0_minus_o[1] * a[1] + v0_minus_o[2] * a[2];
            [
                cyl.axis_origin[0] + t_axis * a[0] + sign * cyl.radius * radial_unit[0],
                cyl.axis_origin[1] + t_axis * a[1] + sign * cyl.radius * radial_unit[1],
                cyl.axis_origin[2] + t_axis * a[2] + sign * cyl.radius * radial_unit[2],
            ]
        }
        _ => v0.point, // fallback: use vertex position
    };

    // Construct Geom_Line passing through the tangent point in the line direction
    let p0 = gp::Pnt::new_real3(tangent_point[0], tangent_point[1], tangent_point[2]);
    let dir = gp::Dir::new_real3(line_dir[0], line_dir[1], line_dir[2]);
    let line = geom::Line::new_pnt_dir(&p0, &dir);
    let line_handle = geom::Line::to_handle(line).to_handle_curve();

    // Trim to vertex endpoints
    let t_start = project_point_on_curve(&v0.point, &line_handle)?;
    let t_end = project_point_on_curve(&v1.point, &line_handle)?;
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
        (SelectedSurface::Spherical(s), SelectedSurface::Cylindrical(c)) => (*s, *c),
        (SelectedSurface::Cylindrical(c), SelectedSurface::Spherical(s)) => (*s, *c),
        _ => return Err("not a sphere-cylinder pair".to_string()),
    };
    let sph = &stage2.spherical_hypotheses[sph_idx];
    let cyl = &stage2.cylindrical_hypotheses[cyl_idx];
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
        (SelectedSurface::Spherical(s), SelectedSurface::Cylindrical(_)) => *s,
        (SelectedSurface::Cylindrical(_), SelectedSurface::Spherical(s)) => *s,
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

/// Select the correct arc on a periodic curve by sampling mesh boundary vertices.
/// Returns (t_lo, t_hi) such that the arc from t_lo to t_hi contains the boundary vertices.
fn select_arc_parameters(
    t0: f64,
    t1: f64,
    mesh_boundary_vertices: &[usize],
    curve_handle: &geom::HandleGeomCurve,
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
    let sample_step = if sample_count > 0 {
        step / sample_count
    } else {
        1
    };

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
    curve_handle: &geom::HandleGeomCurve,
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
        SelectedSurface::Cylindrical(i) => {
            let hyp = &stage2.cylindrical_hypotheses[*i];
            let ao = hyp.axis_origin;
            let ad = hyp.axis_direction;
            let v = [point[0] - ao[0], point[1] - ao[1], point[2] - ao[2]];
            let proj = v[0] * ad[0] + v[1] * ad[1] + v[2] * ad[2];
            let r = [v[0] - proj * ad[0], v[1] - proj * ad[1], v[2] - proj * ad[2]];
            let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
            if len < 1e-10 {
                return [0.0_f32, 0.0, 1.0];
            }
            let sign: f64 = if hyp.convex { 1.0 } else { -1.0 };
            [(sign * r[0] / len) as f32, (sign * r[1] / len) as f32, (sign * r[2] / len) as f32]
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
        SelectedSurface::Cylindrical(i) => &stage2.cylindrical_hypotheses[*i].vertices,
        SelectedSurface::Spherical(i) => &stage2.spherical_hypotheses[*i].vertices,
        SelectedSurface::Conical(i) => &stage2.conical_hypotheses[*i].vertices,
        SelectedSurface::Toroidal(i) => &stage2.toroidal_hypotheses[*i].vertices,
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
                                SelectedSurface::Cylindrical(i) => &output.stage2.cylindrical_hypotheses[*i].faces,
                                SelectedSurface::Spherical(i) => &output.stage2.spherical_hypotheses[*i].faces,
                                SelectedSurface::Conical(i) => &output.stage2.conical_hypotheses[*i].faces,
                                SelectedSurface::Toroidal(i) => &output.stage2.toroidal_hypotheses[*i].faces,
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
                    if config.verbose {
                        eprintln!("  Edge {ei}: FAILED tangent edge - {msg}");
                    }
                    fail_count += 1;
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
                                SelectedSurface::Cylindrical(i) => &output.stage2.cylindrical_hypotheses[*i].faces,
                                SelectedSurface::Spherical(i) => &output.stage2.spherical_hypotheses[*i].faces,
                                SelectedSurface::Conical(i) => &output.stage2.conical_hypotheses[*i].faces,
                                SelectedSurface::Toroidal(i) => &output.stage2.toroidal_hypotheses[*i].faces,
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
                if config.verbose {
                    eprintln!("  Edge {ei}: FAILED - {msg}");
                }
                fail_count += 1;
            }
        }
    }

    if !config.quiet {
        eprintln!(
            "Stage 3.3 ({:.3}s): Computed {success_count}/{total} edge curves ({fail_count} failed)",
            t.elapsed().as_secs_f64(),
        );
    }

    if fail_count > 0 {
        return Err(Stage3Error::EdgeCurveError(format!(
            "{fail_count} edge curves failed to compute"
        )));
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

/// Create an OCCT face for a planar surface, building wires from edge curves.
fn create_planar_face(
    fi: usize,
    output: &Stage3Output,
    topo_edges: &mut [OwnedPtr<b_rep_builder_api::MakeEdge>],
    config: &Config,
) -> Result<OwnedPtr<b_rep_builder_api::MakeFace>, Stage3Error> {
    let fd = &output.face_descriptors[fi];

    // Group edges into wire loops
    let wire_groups = group_edges_into_wires(fd, &output.edges);

    if wire_groups.is_empty() {
        return Err(Stage3Error::AdjacencyError(format!(
            "planar face {fi} has no edges"
        )));
    }

    // Build MakeWire for each group
    let mut wires: Vec<OwnedPtr<b_rep_builder_api::MakeWire>> = Vec::new();
    for (gi, group) in wire_groups.iter().enumerate() {
        let first_edge = topo_edges[group[0]].edge();
        let mut make_wire = b_rep_builder_api::MakeWire::new_edge(first_edge);
        for &ei in &group[1..] {
            make_wire.add_edge(topo_edges[ei].edge());
        }
        if !make_wire.is_done() {
            if config.verbose {
                eprintln!(
                    "  Face {fi} wire {gi}: MakeWire error {:?} ({} edges)",
                    make_wire.error(),
                    group.len(),
                );
            }
            return Err(Stage3Error::AdjacencyError(format!(
                "MakeWire failed for planar face {fi} wire group {gi} ({} edges): {:?}",
                group.len(),
                make_wire.error(),
            )));
        }
        wires.push(make_wire);
    }

    // Identify the outer wire: the one with the most edges (heuristic)
    let outer_idx = wires
        .iter()
        .enumerate()
        .max_by_key(|(i, _)| wire_groups[*i].len())
        .map(|(i, _)| i)
        .unwrap();

    // Create face with the outer wire
    let mut make_face = b_rep_builder_api::MakeFace::new_handlegeomsurface_wire(
        &fd.surface,
        wires[outer_idx].wire(),
    );

    if !make_face.is_done() {
        return Err(Stage3Error::AdjacencyError(format!(
            "MakeFace failed for planar face {fi}: {:?}",
            make_face.error(),
        )));
    }

    // Add inner wires (holes)
    let has_holes = wires.len() > 1;
    for (i, w) in wires.iter_mut().enumerate() {
        if i != outer_idx {
            make_face.add(w.wire());
        }
    }

    // Fix wire orientation: OCCT requires inner wires to be oriented opposite
    // to the outer wire. ShapeFix_Face::fix_orientation() handles this correctly.
    if has_holes {
        let mut fixer = shape_fix::Face::new_face(make_face.face());
        fixer.fix_orientation();
        make_face = b_rep_builder_api::MakeFace::new_face(&fixer.face());
    }

    // Apply ShapeFix_Face to fix any pcurve/edge consistency issues
    {
        let mut fixer = shape_fix::Face::new_face(make_face.face());
        fixer.set_precision(1.0);
        fixer.perform();
        make_face = b_rep_builder_api::MakeFace::new_face(&fixer.face());
    }

    Ok(make_face)
}

/// Compute UV-parameter bounds for a periodic face from its boundary edges.
///
/// Projects edge midpoints and vertices onto the surface to find u and v parameters.
/// For full-revolution faces (all closed loops), u spans [0, 2π].
/// For partial-revolution faces, u spans the arc of the boundary edges.
/// For a hemisphere with one boundary edge, v extends to the pole.
fn compute_uv_bounds_from_edges(
    fi: usize,
    fd: &FaceDescriptor,
    output: &Stage3Output,
    surface: &SelectedSurface,
    config: &Config,
) -> Result<(f64, f64, f64, f64), Stage3Error> {
    let mut u_values: Vec<f64> = Vec::new();
    let mut v_values: Vec<f64> = Vec::new();

    // Project edge midpoints and endpoints onto the surface
    for &ei in &fd.edge_indices {
        let edge = &output.edges[ei];
        let curve = edge.curve_3d.as_ref().unwrap();
        let c = curve.get();

        // Project midpoint
        let mid_t = (c.first_parameter() + c.last_parameter()) / 2.0;
        let mid_pt = c.value(mid_t);
        let proj = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
            &mid_pt,
            &fd.surface,
            extrema::ExtAlgo::Grad,
        );
        if proj.nb_points() > 0 {
            let mut u = 0.0;
            let mut v = 0.0;
            proj.lower_distance_parameters(&mut u, &mut v);
            if config.verbose {
                eprintln!("    edge {ei} midpoint ({:.4},{:.4},{:.4}) -> u={u:.6}, v={v:.6}", mid_pt.x(), mid_pt.y(), mid_pt.z());
            }
            // For spherical surfaces, U is undefined at poles (V=±π/2).
            // Skip U values near poles to avoid corrupting the circular gap algorithm.
            let at_sphere_pole = matches!(surface, SelectedSurface::Spherical(_))
                && v.abs() > std::f64::consts::FRAC_PI_2 - 0.01;
            if !at_sphere_pole {
                u_values.push(u);
            }
            v_values.push(v);
        }

        // Project endpoints (for open edges with vertices)
        for &vi in &edge.vertex_indices {
            if vi != usize::MAX {
                let vtx = &output.vertices[vi];
                let pt = gp::Pnt::new_real3(vtx.point[0], vtx.point[1], vtx.point[2]);
                let proj = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                    &pt,
                    &fd.surface,
                    extrema::ExtAlgo::Grad,
                );
                if proj.nb_points() > 0 {
                    let mut u = 0.0;
                    let mut v = 0.0;
                    proj.lower_distance_parameters(&mut u, &mut v);
                    if config.verbose {
                        eprintln!("    edge {ei} vertex {vi} ({:.4},{:.4},{:.4}) -> u={u:.6}, v={v:.6}", vtx.point[0], vtx.point[1], vtx.point[2]);
                    }
                    // Skip pole U values (same rationale as midpoints above)
                    let at_sphere_pole = matches!(surface, SelectedSurface::Spherical(_))
                        && v.abs() > std::f64::consts::FRAC_PI_2 - 0.01;
                    if !at_sphere_pole {
                        u_values.push(u);
                    }
                    v_values.push(v);
                }
            }
        }
    }

    if v_values.is_empty() {
        return Err(Stage3Error::AdjacencyError(format!(
            "could not project any edge points onto surface for face {fi}"
        )));
    }

    // Compute v bounds from all projected points
    let mut vmin = v_values[0];
    let mut vmax = v_values[0];
    for &v in &v_values[1..] {
        if v < vmin { vmin = v; }
        if v > vmax { vmax = v; }
    }

    // Compute u bounds
    let all_closed_loops = fd.edge_indices.iter().all(|&ei| {
        let e = &output.edges[ei];
        e.vertex_indices[0] == usize::MAX && e.vertex_indices[1] == usize::MAX
    });

    let (umin, umax) = if all_closed_loops {
        // Full revolution: u spans [0, 2π]
        (0.0, 2.0 * std::f64::consts::PI)
    } else {
        // Partial revolution: use circular gap algorithm to handle periodicity.
        // Sort u values, find the largest gap (empty arc), and the face's u range
        // is the complement of that gap.
        let two_pi = 2.0 * std::f64::consts::PI;
        let mut sorted_u: Vec<f64> = u_values.iter().map(|&u| u.rem_euclid(two_pi)).collect();
        sorted_u.sort_by(|a, b| a.partial_cmp(b).unwrap());
        sorted_u.dedup_by(|a, b| (*a - *b).abs() < 1e-10);

        if sorted_u.len() == 1 {
            // Single u value — shouldn't happen for partial revolution, but handle it
            (sorted_u[0], sorted_u[0])
        } else {
            // Find the largest gap between consecutive sorted u values
            let mut max_gap = 0.0_f64;
            let mut gap_end_idx = 0;
            for i in 0..sorted_u.len() {
                let next_i = (i + 1) % sorted_u.len();
                let gap = if next_i > i {
                    sorted_u[next_i] - sorted_u[i]
                } else {
                    sorted_u[next_i] + two_pi - sorted_u[i]
                };
                if gap > max_gap {
                    max_gap = gap;
                    gap_end_idx = next_i;
                }
            }
            // The face spans from sorted_u[gap_end_idx] to sorted_u[gap_start_idx]
            let gap_start_idx = if gap_end_idx == 0 { sorted_u.len() - 1 } else { gap_end_idx - 1 };
            let u_min = sorted_u[gap_end_idx];
            let mut u_max = sorted_u[gap_start_idx];
            if u_max < u_min {
                u_max += two_pi;
            }
            (u_min, u_max)
        }
    };

    // For spherical surfaces with a single boundary edge, the face extends to a pole.
    // Determine which pole by checking where the mesh face centroids are.
    if matches!(surface, SelectedSurface::Spherical(_)) && fd.edge_indices.len() == 1 {
        let hyp_idx = match surface {
            SelectedSurface::Spherical(i) => *i,
            _ => unreachable!(),
        };
        let hyp = &output.stage2.spherical_hypotheses[hyp_idx];

        // Average v of mesh face centroids
        let mut v_sum = 0.0;
        let mut v_count = 0;
        for &mfi in hyp.faces.iter().take(10) {
            let face = &output.stage2.mesh.faces[mfi];
            let v0 = &output.stage2.mesh.vertices[face.vertex_indices[0]];
            let v1 = &output.stage2.mesh.vertices[face.vertex_indices[1]];
            let v2 = &output.stage2.mesh.vertices[face.vertex_indices[2]];
            let cx = (v0.x + v1.x + v2.x) / 3.0;
            let cy = (v0.y + v1.y + v2.y) / 3.0;
            let cz = (v0.z + v1.z + v2.z) / 3.0;
            let pt = gp::Pnt::new_real3(cx, cy, cz);
            let proj = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                &pt,
                &fd.surface,
                extrema::ExtAlgo::Grad,
            );
            if proj.nb_points() > 0 {
                let mut u = 0.0;
                let mut v = 0.0;
                proj.lower_distance_parameters(&mut u, &mut v);
                v_sum += v;
                v_count += 1;
            }
        }

        if v_count > 0 {
            let v_avg = v_sum / v_count as f64;
            let v_edge = vmin; // single edge, so vmin == vmax
            if v_avg > v_edge {
                vmax = std::f64::consts::FRAC_PI_2;
            } else {
                vmin = -std::f64::consts::FRAC_PI_2;
            }
        }
    }

    // For conical surfaces with a single boundary edge, the face extends to the
    // apex (a degenerate point at V=0). Determine direction from mesh centroids.
    if matches!(surface, SelectedSurface::Conical(_)) && fd.edge_indices.len() == 1 {
        let hyp_idx = match surface {
            SelectedSurface::Conical(i) => *i,
            _ => unreachable!(),
        };
        let hyp = &output.stage2.conical_hypotheses[hyp_idx];

        // Average v of mesh face centroids to determine which side of the edge the apex is on
        let mut v_sum = 0.0;
        let mut v_count = 0;
        for &mfi in hyp.faces.iter().take(10) {
            let face = &output.stage2.mesh.faces[mfi];
            let v0 = &output.stage2.mesh.vertices[face.vertex_indices[0]];
            let v1 = &output.stage2.mesh.vertices[face.vertex_indices[1]];
            let v2 = &output.stage2.mesh.vertices[face.vertex_indices[2]];
            let cx = (v0.x + v1.x + v2.x) / 3.0;
            let cy = (v0.y + v1.y + v2.y) / 3.0;
            let cz = (v0.z + v1.z + v2.z) / 3.0;
            let pt = gp::Pnt::new_real3(cx, cy, cz);
            let proj = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                &pt,
                &fd.surface,
                extrema::ExtAlgo::Grad,
            );
            if proj.nb_points() > 0 {
                let mut u = 0.0;
                let mut v = 0.0;
                proj.lower_distance_parameters(&mut u, &mut v);
                v_sum += v;
                v_count += 1;
            }
        }

        if v_count > 0 {
            let v_avg = v_sum / v_count as f64;
            let v_edge = vmin; // single edge, so vmin == vmax
            if v_avg < v_edge {
                // Face extends from apex (V=0) toward the base edge
                vmin = 0.0;
            } else {
                // Face extends from the base edge toward the apex (V=0)
                vmax = 0.0;
            }
        }
    }

    if config.verbose {
        eprintln!("  Face {fi}: computed UV bounds u=[{umin:.6}, {umax:.6}], v=[{vmin:.6}, {vmax:.6}] from {} edge(s)", fd.edge_indices.len());
    }

    Ok((umin, umax, vmin, vmax))
}
/// Compute the area of a face via BRepGProp.
fn compute_face_area(make_face: &b_rep_builder_api::MakeFace) -> f64 {
    let face = make_face.face();
    compute_face_area_from_face(face)
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

/// Create a wire-based OCCT face for a periodic surface with pre-set pcurves.
///
/// The `face_surface` parameter is the surface to use for pcurve computation
/// and face construction. This may differ from `fd.surface` when sphere
/// reorientation is needed to avoid pole singularities.
#[allow(clippy::too_many_arguments)]
fn create_wire_based_periodic_face(
    fi: usize,
    face_surface: &geom::HandleGeomSurface,
    umin: f64,
    umax: f64,
    output: &Stage3Output,
    topo_edges: &mut [OwnedPtr<b_rep_builder_api::MakeEdge>],
    config: &Config,
    verbose_prefix: &str,
    is_sphere: bool,
) -> Result<OwnedPtr<b_rep_builder_api::MakeFace>, Stage3Error> {
    let fd = &output.face_descriptors[fi];
    let two_pi = 2.0 * std::f64::consts::PI;
    let identity_loc = top_loc::Location::new();
    let builder = b_rep::Builder::new();
    let n_samples = 9;
    let u_center = (umin + umax) / 2.0;
    let pole_v_threshold = std::f64::consts::FRAC_PI_2 - 0.1;
    // Threshold for detecting a vertex AT the exact pole (vs. just near it)
    let pole_v_exact = std::f64::consts::FRAC_PI_2 - 0.01;

    for &ei in &fd.edge_indices {
        let edge = &output.edges[ei];
        let curve = edge.curve_3d.as_ref().unwrap();
        let c = curve.get();

        let t_start = c.first_parameter();
        let t_end = c.last_parameter();
        let span = t_end - t_start;
        if span.abs() < 1e-15 {
            continue;
        }

        let mut uv_points: Vec<(f64, f64, f64)> = Vec::new();
        for i in 0..n_samples {
            let t = t_start + span * (i as f64) / ((n_samples - 1) as f64);
            let pt = c.value(t);
            let proj = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                &pt, face_surface, extrema::ExtAlgo::Grad,
            );
            if proj.nb_points() == 0 {
                continue;
            }
            let mut u = 0.0;
            let mut v = 0.0;
            proj.lower_distance_parameters(&mut u, &mut v);
            uv_points.push((t, u, v));
        }
        if uv_points.len() < 2 {
            continue;
        }

        // Fix U at sphere poles (arbitrary projection U).
        // At the pole (|V| ≈ π/2), U is degenerate — use nearest non-pole sample's U.
        {
            let last = uv_points.len() - 1;
            if uv_points[0].2.abs() > pole_v_threshold {
                if let Some(j) = (1..uv_points.len())
                    .find(|&j| uv_points[j].2.abs() <= pole_v_threshold)
                {
                    uv_points[0].1 = uv_points[j].1;
                }
            }
            if last > 0 && uv_points[last].2.abs() > pole_v_threshold {
                if let Some(j) = (0..last)
                    .rev()
                    .find(|&j| uv_points[j].2.abs() <= pole_v_threshold)
                {
                    uv_points[last].1 = uv_points[j].1;
                }
            }
        }

        // Unwrap U for continuity
        while uv_points[0].1 < u_center - std::f64::consts::PI { uv_points[0].1 += two_pi; }
        while uv_points[0].1 > u_center + std::f64::consts::PI { uv_points[0].1 -= two_pi; }
        for i in 1..uv_points.len() {
            let prev_u = uv_points[i - 1].1;
            while uv_points[i].1 < prev_u - std::f64::consts::PI { uv_points[i].1 += two_pi; }
            while uv_points[i].1 > prev_u + std::f64::consts::PI { uv_points[i].1 -= two_pi; }
        }

        let n = uv_points.len() as i32;
        let mut pts_h = t_colgp::HArray1OfPnt2d::new_int2(1, n);
        let mut params_h = t_col_std::HArray1OfReal::new_int2(1, n);
        for (idx, &(t, u, v)) in uv_points.iter().enumerate() {
            let pnt2d = gp::Pnt2d::new_real2(u, v);
            pts_h.change_array1().set_value((idx + 1) as i32, &pnt2d);
            params_h.change_array1().set_value((idx + 1) as i32, &t);
        }

        let pts_handle = t_colgp::HArray1OfPnt2d::to_handle(pts_h);
        let params_handle = t_col_std::HArray1OfReal::to_handle(params_h);
        let mut interp = geom2d_api::Interpolate::new_handletcolgpharray1ofpnt2d_handletcolstdharray1ofreal_bool_real(
            &pts_handle, &params_handle, false, 1e-10,
        );
        interp.perform();
        if !interp.is_done() {
            let (u1, v1) = (uv_points[0].1, uv_points[0].2);
            let (u2, v2) = (uv_points.last().unwrap().1, uv_points.last().unwrap().2);
            let du = u2 - u1;
            let dv = v2 - v1;
            let len = (du * du + dv * dv).sqrt();
            if len < 1e-15 { continue; }
            let dir_x = du / len;
            let dir_y = dv / len;
            let ox = u1 - uv_points[0].0 * dir_x;
            let oy = v1 - uv_points[0].0 * dir_y;
            let pcurve_origin = gp::Pnt2d::new_real2(ox, oy);
            let pcurve_dir = gp::Dir2d::new_real2(dir_x, dir_y);
            let pcurve_line = geom2d::Line::new_pnt2d_dir2d(&pcurve_origin, &pcurve_dir);
            let pcurve_handle = geom2d::Line::to_handle(pcurve_line).to_handle_curve();
            builder.update_edge_edge_handlegeom2dcurve_handlegeomsurface_location_real(
                topo_edges[ei].edge(),
                &pcurve_handle,
                face_surface,
                &identity_loc,
                config.vertex_tolerance_mm,
            );
            if config.verbose {
                eprintln!("    Edge {ei}: set linear pcurve (fallback) ({u1:.4},{v1:.4}) -> ({u2:.4},{v2:.4})");
            }
            continue;
        }

        let pcurve_handle = interp.curve().to_handle_curve();
        builder.update_edge_edge_handlegeom2dcurve_handlegeomsurface_location_real(
            topo_edges[ei].edge(),
            &pcurve_handle,
            face_surface,
            &identity_loc,
            config.vertex_tolerance_mm,
        );

        if config.verbose {
            let (u1, v1) = (uv_points[0].1, uv_points[0].2);
            let (u2, v2) = (uv_points.last().unwrap().1, uv_points.last().unwrap().2);
            eprintln!("    Edge {ei}: set sampled pcurve ({u1:.4},{v1:.4}) -> ({u2:.4},{v2:.4}), {n} pts");
        }
    }

    // Build wire from edges
    let wire_groups = group_edges_into_wires(fd, &output.edges);
    if wire_groups.is_empty() {
        return Err(Stage3Error::AdjacencyError(format!(
            "periodic face {fi} has no edges for wire construction"
        )));
    }

    let mut wires: Vec<OwnedPtr<topo_ds::Wire>> = Vec::new();
    for (gi, group) in wire_groups.iter().enumerate() {
        let first_edge = topo_edges[group[0]].edge();
        let mut make_wire = b_rep_builder_api::MakeWire::new_edge(first_edge);
        for &ei in &group[1..] {
            make_wire.add_edge(topo_edges[ei].edge());
        }
        if !make_wire.is_done() {
            if config.verbose {
                eprintln!(
                    "  Face {fi} wire {gi}: MakeWire error {:?} ({} edges)",
                    make_wire.error(),
                    group.len(),
                );
            }
            return Err(Stage3Error::AdjacencyError(format!(
                "MakeWire failed for periodic face {fi} wire group {gi} ({} edges): {:?}",
                group.len(),
                make_wire.error(),
            )));
        }

        // Apply ShapeFix_Wire to fix pcurve and 3D gaps at shared vertices
        let mut wire_fixer = shape_fix::Wire::new();
        wire_fixer.load_wire(make_wire.wire());
        wire_fixer.set_surface_handlegeomsurface(face_surface);
        wire_fixer.set_precision(1.0);
        *wire_fixer.fix_reorder_mode() = 1;
        *wire_fixer.fix_connected_mode() = 1;
        *wire_fixer.fix_gaps2d_mode() = 1;
        *wire_fixer.fix_gaps3d_mode() = 1;
        *wire_fixer.fix_edge_curves_mode() = 1;
        *wire_fixer.fix_add_p_curve_mode() = 1;
        *wire_fixer.fix_same_parameter_mode() = 1;
        *wire_fixer.fix_shifted_mode() = 1;
        *wire_fixer.fix_seam_mode() = 1;
        *wire_fixer.fix_lacking_mode() = 1;
        *wire_fixer.closed_wire_mode() = true;
        wire_fixer.perform();
        let fixed_wire = wire_fixer.wire();

        // For sphere faces, insert degenerate edges at pole vertices.
        // At a pole (V=±π/2), U is degenerate — two edges meeting at the pole
        // have legitimately different U values. A degenerate edge bridges this
        // gap: a single 3D point but a UV line segment at constant V.
        if is_sphere {
            // Collect pcurve endpoint info for each edge in wire order.
            // For each edge: (start_u, start_v, end_u, end_v)
            #[allow(dead_code)]
            struct EdgeUV { start_u: f64, start_v: f64, end_u: f64, end_v: f64 }
            let mut edge_uvs: Vec<EdgeUV> = Vec::new();
            {
                let mut we = b_rep_tools::WireExplorer::new_wire(&fixed_wire);
                while we.more() {
                    let edge = we.current();
                    let mut first = 0.0;
                    let mut last = 0.0;
                    let c2d = b_rep::Tool::curve_on_surface_edge_handlegeomsurface_location_real2_boolptr(
                        edge, face_surface, &identity_loc, &mut first, &mut last, None,
                    );
                    if !c2d.is_null() {
                        let c = c2d.get();
                        let orient = edge.as_shape().orientation();
                        let (sp, ep) = if orient == top_abs::Orientation::Reversed {
                            (c.value(last), c.value(first))
                        } else {
                            (c.value(first), c.value(last))
                        };
                        if config.verbose {
                            eprintln!("    Wire edge: ({:.4},{:.4}) -> ({:.4},{:.4}) orient={:?}",
                                sp.x(), sp.y(), ep.x(), ep.y(), orient);
                        }
                        edge_uvs.push(EdgeUV {
                            start_u: sp.x(), start_v: sp.y(),
                            end_u: ep.x(), end_v: ep.y(),
                        });
                    } else {
                        if config.verbose {
                            eprintln!("    Wire edge: pcurve NOT FOUND on surface");
                        }
                        edge_uvs.push(EdgeUV {
                            start_u: u_center, start_v: 0.0,
                            end_u: u_center, end_v: 0.0,
                        });
                    }
                    we.next();
                }
            }

            // Find edges whose end is at the pole (gap to next edge's start)
            let n_edges = edge_uvs.len();
            let mut pole_insertions: Vec<(usize, f64, f64, f64)> = Vec::new(); // (after_edge_idx, u_prev, u_next, pole_v)
            for i in 0..n_edges {
                let next = (i + 1) % n_edges;
                let end_v = edge_uvs[i].end_v;
                if end_v.abs() > pole_v_exact {
                    let u_prev = edge_uvs[i].end_u;
                    let u_next = edge_uvs[next].start_u;
                    if (u_prev - u_next).abs() > 0.01 {
                        pole_insertions.push((i, u_prev, u_next, end_v));
                    }
                }
            }

            if !pole_insertions.is_empty() {
                // Rebuild wire with degenerate edges inserted
                let mut new_wire = topo_ds::Wire::new();
                builder.make_wire(&mut new_wire);

                let mut edge_idx = 0usize;
                let mut we = b_rep_tools::WireExplorer::new_wire(&fixed_wire);
                while we.more() {
                    let edge = we.current();
                    // Add the regular edge
                    builder.add(new_wire.as_shape_mut(), edge.as_shape());

                    // Check if a degenerate edge needs to be inserted after this edge
                    if let Some(&(_, u_prev, u_next, pole_v)) = pole_insertions.iter().find(|p| p.0 == edge_idx) {
                        let u_lo = u_prev.min(u_next);
                        let u_hi = u_prev.max(u_next);

                        // Get the pole vertex from this edge's end
                        let pole_vertex = top_exp::last_vertex(edge, true);

                        // Create degenerate edge
                        let mut degen = topo_ds::Edge::new();
                        builder.make_edge_edge(&mut degen);
                        builder.degenerated(&degen, true);

                        // Pcurve: line at constant V from u_prev to u_next
                        let pcurve_origin = gp::Pnt2d::new_real2(u_prev, pole_v);
                        let dir_u = u_next - u_prev;
                        let pcurve_dir = gp::Dir2d::new_real2(if dir_u >= 0.0 { 1.0 } else { -1.0 }, 0.0);
                        let pcurve_line = geom2d::Line::new_pnt2d_dir2d(&pcurve_origin, &pcurve_dir);
                        let pcurve_handle = geom2d::Line::to_handle(pcurve_line).to_handle_curve();

                        builder.update_edge_edge_handlegeom2dcurve_handlegeomsurface_location_real(
                            &degen, &pcurve_handle, face_surface, &identity_loc,
                            config.vertex_tolerance_mm,
                        );

                        // Range: parameterized as distance along U from u_prev
                        builder.range_edge_real2_bool(&degen, 0.0, dir_u.abs(), false);

                        // Add vertices (same vertex at both ends)
                        let v_fwd = pole_vertex.as_shape().oriented(top_abs::Orientation::Forward);
                        let v_rev = pole_vertex.as_shape().oriented(top_abs::Orientation::Reversed);
                        builder.add(degen.as_shape_mut(), &v_fwd);
                        builder.add(degen.as_shape_mut(), &v_rev);
                        builder.update_vertex_vertex_real_edge_real(
                            &pole_vertex, 0.0, &degen, config.vertex_tolerance_mm,
                        );

                        builder.add(new_wire.as_shape_mut(), degen.as_shape());

                        if config.verbose {
                            eprintln!("    Inserted degenerate edge at pole: u=[{u_lo:.4}, {u_hi:.4}], v={pole_v:.4}");
                        }
                    }

                    edge_idx += 1;
                    we.next();
                }

                wires.push(new_wire);
                continue; // Skip the normal wire push below
            }
        }

        wires.push(fixed_wire);
    }

    let outer_idx = wires
        .iter()
        .enumerate()
        .max_by_key(|(i, _)| wire_groups[*i].len())
        .map(|(i, _)| i)
        .unwrap();

    let mut mf = b_rep_builder_api::MakeFace::new_handlegeomsurface_wire(
        face_surface,
        &wires[outer_idx],
    );

    if !mf.is_done() {
        return Err(Stage3Error::AdjacencyError(format!(
            "MakeFace failed for periodic face {fi}: {:?}",
            mf.error(),
        )));
    }

    let has_holes = wires.len() > 1;
    for (i, w) in wires.iter().enumerate() {
        if i != outer_idx {
            mf.add(w);
        }
    }

    if has_holes {
        let mut fixer = shape_fix::Face::new_face(mf.face());
        fixer.fix_orientation();
        mf = b_rep_builder_api::MakeFace::new_face(&fixer.face());
    }

    // Apply ShapeFix_Face to fix pcurve gaps at shared vertices.
    // The pcurves from edge sampling may have UV gaps where 3D edge curves
    // don't share endpoints exactly. ShapeFix_Face fixes wires/edges/pcurves.
    {
        let mut fixer = shape_fix::Face::new_face(mf.face());
        fixer.set_precision(1.0); // cover gaps up to ~0.6mm
        // Enable periodic degenerated fix (adds degenerate edges at sphere poles)
        *fixer.fix_periodic_degenerated_mode() = 1;
        // Enable natural bound addition for single-wire faces (off by default)
        *fixer.fix_add_natural_bound_mode() = 1;
        fixer.perform();

        // Explicitly call fix_periodic_degenerated and fix_add_natural_bound
        let fpd = fixer.fix_periodic_degenerated();
        let fanb = fixer.fix_add_natural_bound();
        if config.verbose {
            eprintln!("  Face {fi}: fix_periodic_degenerated={fpd}, fix_add_natural_bound={fanb}");
        }

        mf = b_rep_builder_api::MakeFace::new_face(&fixer.face());
    }

    if config.verbose {
        let area = compute_face_area(&mf);
        eprintln!("  Face {fi}: {verbose_prefix} wire construction with pcurves, u=[{umin:.4}, {umax:.4}], area={area:.4} mm\u{b2}");
    }
    Ok(mf)
}

/// Create an OCCT face for a periodic surface (cylinder or sphere).
///
/// For full spheres with no edges, uses natural surface bounds.
/// For full-revolution surfaces (all boundary edges are closed loops),
/// uses UV-bounds construction which automatically creates seam edges.
/// For partial-revolution surfaces, uses wire-based construction with
/// pre-computed pcurves to ensure the correct arc is selected.
fn create_periodic_face(
    fi: usize,
    output: &Stage3Output,
    topo_edges: &mut [OwnedPtr<b_rep_builder_api::MakeEdge>],
    config: &Config,
) -> Result<(OwnedPtr<b_rep_builder_api::MakeFace>, bool), Stage3Error> {
    let fd = &output.face_descriptors[fi];
    let surface = &output.stage2.selected_surfaces[fd.selected_surface_idx];

    // Full sphere with no edges — use natural bounds
    if matches!(surface, SelectedSurface::Spherical(_)) && fd.edge_indices.is_empty() {
        let make_face = b_rep_builder_api::MakeFace::new_handlegeomsurface_real(
            &fd.surface,
            config.vertex_tolerance_mm,
        );
        if !make_face.is_done() {
            return Err(Stage3Error::AdjacencyError(format!(
                "MakeFace failed for full sphere face {fi}: {:?}",
                make_face.error(),
            )));
        }
        return Ok((make_face, false)); // Full spheres are always convex
    }

    // Compute UV bounds for all periodic faces
    let (umin, umax, vmin, vmax) = compute_uv_bounds_from_edges(fi, fd, output, surface, config)?;

    // Check if all boundary edges are closed loops (full revolution)
    let all_closed_loops = fd.edge_indices.iter().all(|&ei| {
        let e = &output.edges[ei];
        e.vertex_indices[0] == usize::MAX && e.vertex_indices[1] == usize::MAX
    });

    // For spherical faces with all-closed-loop edges, the boundary circles
    // may pass through the sphere's UV singularities (poles). This happens when
    // the intersection plane contains the sphere axis (e.g., a hemisphere cut by
    // a meridional plane, or a pill/capsule shape). UV-bounds on the original
    // sphere creates wrong face regions. Detect this and reorient the sphere.
    let has_closed_loop_sphere = matches!(surface, SelectedSurface::Spherical(_))
        && all_closed_loops
        && !fd.edge_indices.is_empty();

    let make_face = if has_closed_loop_sphere {
        // Sphere with all-closed-loop edges.
        // The boundary circles may pass through or near the sphere's UV
        // singularities (poles). When within 45°, UV-bounds on the original
        // sphere creates wrong face regions (boundaries are meridians, not
        // circles at constant V).
        //
        // For tangent edges (cylinder-sphere): reorient to align with cylinder axis.
        // For non-tangent edges (plane-sphere): reorient to align with the
        // boundary circle's normal (= plane normal for sphere-plane intersections).
        let sph_idx = match surface {
            SelectedSurface::Spherical(i) => *i,
            _ => unreachable!(),
        };
        let sph = &output.stage2.spherical_hypotheses[sph_idx];
        let sphere_z = [0.0_f64, 0.0, 1.0]; // sphere axis is always Z-up
        let cos_45 = std::f64::consts::FRAC_PI_4.cos(); // cos(45°) ≈ 0.707

        // Check if any boundary edge circle comes within 45° of a pole
        let mut near_pole = false;
        for &ei in &fd.edge_indices {
            let edge = &output.edges[ei];
            let curve_handle = edge.curve_3d.as_ref().unwrap();
            for pole_sign in &[1.0_f64, -1.0] {
                let pole = gp::Pnt::new_real3(
                    sph.center[0] + sph.radius * pole_sign * sphere_z[0],
                    sph.center[1] + sph.radius * pole_sign * sphere_z[1],
                    sph.center[2] + sph.radius * pole_sign * sphere_z[2],
                );
                let proj = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(
                    &pole, curve_handle,
                );
                if proj.nb_points() > 0 {
                    let dist = proj.lower_distance();
                    let cos_angle = 1.0 - dist * dist / (2.0 * sph.radius * sph.radius);
                    if cos_angle > cos_45 {
                        near_pole = true;
                    }
                }
            }
        }

        if !near_pole {
            // Circle doesn't pass near poles — UV-bounds with full U range works
            let mf = b_rep_builder_api::MakeFace::new_handlegeomsurface_real5(
                &fd.surface, umin, umax, vmin, vmax, config.vertex_tolerance_mm,
            );
            if !mf.is_done() {
                return Err(Stage3Error::AdjacencyError(format!(
                    "MakeFace (UV bounds) failed for sphere face {fi}: {:?}",
                    mf.error(),
                )));
            }
            if config.verbose {
                let area = compute_face_area(&mf);
                eprintln!("  Face {fi}: sphere UV-bounds (no pole issue), area={area:.4} mm\u{b2}");
            }
            mf
        } else {
            // Circle passes within 45° of a pole. Reorient the sphere surface
            // so its Z-axis aligns with the boundary circle normal. This makes
            // the boundary circles into iso-V curves, allowing UV-bounds
            // construction to work correctly.
            //
            // Determine the reorientation axis:
            // - For tangent edges: use the adjacent cylinder's axis direction
            // - For non-tangent edges: use the boundary circle's normal
            //   (derived from the edge curve, which is a Geom_Circle)
            let reorient_axis = {
                let mut axis: Option<[f64; 3]> = None;

                // Try cylinder axis from tangent edges first
                for &ei in &fd.edge_indices {
                    let edge = &output.edges[ei];
                    if !edge.tangent { continue; }
                    for &adj_fi in &edge.face_indices {
                        if adj_fi == fi { continue; }
                        let adj_fd = &output.face_descriptors[adj_fi];
                        let adj_ss = &output.stage2.selected_surfaces[adj_fd.selected_surface_idx];
                        if let SelectedSurface::Cylindrical(idx) = adj_ss {
                            axis = Some(output.stage2.cylindrical_hypotheses[*idx].axis_direction);
                            break;
                        }
                    }
                    if axis.is_some() { break; }
                }

                // Fall back to boundary circle normal (for non-tangent edges)
                if axis.is_none() {
                    for &ei in &fd.edge_indices {
                        let edge = &output.edges[ei];
                        let curve = edge.curve_3d.as_ref().unwrap().get();
                        let fp = curve.first_parameter();
                        let lp = curve.last_parameter();
                        let span = lp - fp;
                        // Sample 3 points on the curve to determine circle normal
                        let p0 = curve.value(fp);
                        let p1 = curve.value(fp + span * 0.25);
                        let p2 = curve.value(fp + span * 0.75);
                        let a = [p1.x() - p0.x(), p1.y() - p0.y(), p1.z() - p0.z()];
                        let b = [p2.x() - p0.x(), p2.y() - p0.y(), p2.z() - p0.z()];
                        let n = [
                            a[1] * b[2] - a[2] * b[1],
                            a[2] * b[0] - a[0] * b[2],
                            a[0] * b[1] - a[1] * b[0],
                        ];
                        let n_len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                        if n_len > 1e-10 {
                            axis = Some([n[0] / n_len, n[1] / n_len, n[2] / n_len]);
                            break;
                        }
                    }
                }

                axis.ok_or_else(|| Stage3Error::AdjacencyError(format!(
                    "sphere face {fi} near pole: could not determine reorientation axis"
                )))?
            };

            // Recreate the sphere surface with its axis aligned to the reorientation axis
            let origin = gp::Pnt::new_real3(sph.center[0], sph.center[1], sph.center[2]);
            let dir = gp::Dir::new_real3(reorient_axis[0], reorient_axis[1], reorient_axis[2]);
            let ax3 = gp::Ax3::new_pnt_dir(&origin, &dir);
            let oriented_sphere = geom::SphericalSurface::new_ax3_real(&ax3, sph.radius);
            let oriented_surface = geom::SphericalSurface::to_handle(oriented_sphere).to_handle_surface();

            // Determine which hemisphere the face occupies by projecting a mesh
            // centroid onto the oriented surface and checking the V sign.
            let sample_face = &output.stage2.mesh.faces[sph.faces[0]];
            let sv0 = &output.stage2.mesh.vertices[sample_face.vertex_indices[0]];
            let sv1 = &output.stage2.mesh.vertices[sample_face.vertex_indices[1]];
            let sv2 = &output.stage2.mesh.vertices[sample_face.vertex_indices[2]];
            let sc = gp::Pnt::new_real3(
                (sv0.x + sv1.x + sv2.x) / 3.0,
                (sv0.y + sv1.y + sv2.y) / 3.0,
                (sv0.z + sv1.z + sv2.z) / 3.0,
            );
            let sample_proj = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                &sc, &oriented_surface, extrema::ExtAlgo::Grad,
            );
            let v_side = if sample_proj.nb_points() > 0 {
                let mut _u = 0.0;
                let mut v = 0.0;
                sample_proj.lower_distance_parameters(&mut _u, &mut v);
                v
            } else {
                1.0 // default to positive hemisphere
            };

            let half_pi = std::f64::consts::FRAC_PI_2;
            let two_pi = 2.0 * std::f64::consts::PI;
            let (v_lo, v_hi) = if v_side > 0.0 {
                (0.0, half_pi)
            } else {
                (-half_pi, 0.0)
            };

            // Use UV-bounds construction on the oriented surface
            let mf = b_rep_builder_api::MakeFace::new_handlegeomsurface_real5(
                &oriented_surface,
                0.0,
                two_pi,
                v_lo,
                v_hi,
                config.vertex_tolerance_mm,
            );
            if !mf.is_done() {
                return Err(Stage3Error::AdjacencyError(format!(
                    "MakeFace (UV bounds) failed for oriented sphere face {fi}: {:?}",
                    mf.error(),
                )));
            }

            if config.verbose {
                let area = compute_face_area(&mf);
                eprintln!("  Face {fi}: oriented sphere UV-bounds, u=[0.0000, {two_pi:.4}], v=[{v_lo:.4}, {v_hi:.4}], area={area:.4} mm\u{b2}");
            }
            mf
        }
    } else if matches!(surface, SelectedSurface::Spherical(_))
        && !fd.edge_indices.is_empty()
        && !all_closed_loops
        && (vmin < -(std::f64::consts::FRAC_PI_2 - 0.01)
            || vmax > std::f64::consts::FRAC_PI_2 - 0.01)
    {
        // Spherical face with a boundary vertex at or near a UV pole.
        // Use wire-based construction with pole-aware pcurves: at the pole,
        // U is degenerate, so we derive U from the edge midpoint projection.
        create_wire_based_periodic_face(
            fi, &fd.surface, umin, umax,
            output, topo_edges, config, "sphere pole", true,
        )?
    } else if all_closed_loops && !fd.edge_indices.is_empty() {
        // Full-revolution periodic face: use UV-bounds construction.
        // This automatically creates seam edges needed for proper pcurves.
        let mf = b_rep_builder_api::MakeFace::new_handlegeomsurface_real5(
            &fd.surface,
            umin,
            umax,
            vmin,
            vmax,
            config.vertex_tolerance_mm,
        );
        if !mf.is_done() {
            return Err(Stage3Error::AdjacencyError(format!(
                "MakeFace (UV bounds) failed for periodic face {fi}: {:?}",
                mf.error(),
            )));
        }
        if config.verbose {
            let area = compute_face_area(&mf);
            eprintln!("  Face {fi}: UV-bounds construction, u=[{umin:.4}, {umax:.4}], v=[{vmin:.4}, {vmax:.4}], area={area:.4} mm\u{b2}");
        }
        mf
    } else {
        // Partial-revolution: use wire-based construction with pre-set pcurves.
        create_wire_based_periodic_face(
            fi, &fd.surface, umin, umax,
            output, topo_edges, config, "wire", false,
        )?
    };

    let is_concave = match surface {
        SelectedSurface::Cylindrical(idx) => !output.stage2.cylindrical_hypotheses[*idx].convex,
        SelectedSurface::Spherical(idx) => !output.stage2.spherical_hypotheses[*idx].convex,
        SelectedSurface::Conical(idx) => !output.stage2.conical_hypotheses[*idx].convex,
        SelectedSurface::Toroidal(idx) => !output.stage2.toroidal_hypotheses[*idx].convex,
        _ => false,
    };

    Ok((make_face, is_concave))
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

    for (fi, mf) in output.make_faces.iter().enumerate() {
        let face = mf.face();
        let analyzer = b_rep_check::Analyzer::new_shape(face.as_shape());
        if !analyzer.is_valid() {
            // Run detailed BRepCheck_Face diagnostics
            let mut checker = b_rep_check::Face::new_face(face);
            checker.minimum();
            let intersect = checker.intersect_wires(false);
            let classify = checker.classify_wires(false);
            let orient = checker.orientation_of_wires(false);
            let unorientable = checker.is_unorientable();

            // Check sub-shape edges
            let mut edge_issues = Vec::new();
            let mut only_near_zero_tol_edges = true;
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
                        if tol > config.vertex_tolerance_mm {
                            only_near_zero_tol_edges = false;
                        }
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

            let has_face_issues = intersect != b_rep_check::Status::Noerror
                || classify != b_rep_check::Status::Noerror
                || orient != b_rep_check::Status::Noerror
                || unorientable;

            // If the only failures are edges with near-zero tolerance,
            // skip the warning. This occurs on UV-bounds sphere faces near
            // poles and is resolved by sewing.
            if !has_face_issues
                && vertex_issues.is_empty()
                && !edge_issues.is_empty()
                && only_near_zero_tol_edges
            {
                // Harmless: edge BRepCheck failure with ~zero tolerance on
                // sphere pole face, resolved by sewing
            } else {
                invalid_count += 1;
                if !config.quiet {
                    let surface = &output.stage2.selected_surfaces
                        [output.face_descriptors[fi].selected_surface_idx];
                    let surface_desc = match surface {
                        SelectedSurface::Planar(_) => "planar",
                        SelectedSurface::Cylindrical(_) => "cylindrical",
                        SelectedSurface::Spherical(_) => "spherical",
                        SelectedSurface::Conical(_) => "conical",
                        SelectedSurface::Toroidal(_) => "toroidal",
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
            output.make_faces.len(),
        );
    }

    Ok(())
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

    let our_face_count = output.make_faces.len();
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
            SelectedSurface::Cylindrical(_) => "cylindrical",
            SelectedSurface::Spherical(_) => "spherical",
            SelectedSurface::Conical(_) => "conical",
            SelectedSurface::Toroidal(_) => "toroidal",
        };

        let is_convex = match surface {
            SelectedSurface::Cylindrical(idx) => {
                Some(output.stage2.cylindrical_hypotheses[*idx].convex)
            }
            SelectedSurface::Spherical(idx) => {
                Some(output.stage2.spherical_hypotheses[*idx].convex)
            }
            SelectedSurface::Conical(idx) => {
                Some(output.stage2.conical_hypotheses[*idx].convex)
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

/// Create OCCT faces for all face descriptors (stage 3.4 main entry).
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

    // 2. Create TopoDS_Edge for each ReconEdge, using shared vertices when available
    let mut topo_edges: Vec<OwnedPtr<b_rep_builder_api::MakeEdge>> = Vec::with_capacity(output.edges.len());
    for (ei, edge) in output.edges.iter().enumerate() {
        let curve = edge.curve_3d.as_ref().unwrap_or_else(|| {
            panic!("edge {ei} has no 3D curve \u{2014} stage 3.3 must run first")
        });
        let make_edge = if edge.vertex_indices[0] != usize::MAX
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

            b_rep_builder_api::MakeEdge::new_handlegeomcurve_vertex2_real2(
                curve,
                &topo_vertices[vi0],
                &topo_vertices[vi1],
                fp,
                lp,
            )
        } else {
            // Closed-loop edge (no vertex endpoints) — create without vertices
            b_rep_builder_api::MakeEdge::new_handlegeomcurve(curve)
        };
        if !make_edge.is_done() {
            return Err(Stage3Error::AdjacencyError(format!(
                "MakeEdge failed for edge {ei}: {:?}",
                make_edge.error(),
            )));
        }
        topo_edges.push(make_edge);
    }

    // 2. Create OCCT face for each face descriptor
    let mut make_faces: Vec<OwnedPtr<b_rep_builder_api::MakeFace>> = Vec::with_capacity(num_faces);
    let mut concave_faces: Vec<bool> = vec![false; num_faces];

    for (fi, concave_flag) in concave_faces.iter_mut().enumerate() {
        let surface_type = &output.stage2.selected_surfaces[output.face_descriptors[fi].selected_surface_idx];
        let is_periodic = matches!(
            surface_type,
            SelectedSurface::Cylindrical(_) | SelectedSurface::Spherical(_) | SelectedSurface::Conical(_)
            | SelectedSurface::Toroidal(_)
        );

        let (mut make_face, is_concave) = if is_periodic {
            create_periodic_face(fi, &output, &mut topo_edges, config)?
        } else {
            (create_planar_face(fi, &output, &mut topo_edges, config)?, false)
        };
        *concave_flag = is_concave;

        if config.verbose {
            let stype = match surface_type {
                SelectedSurface::Planar(_) => "planar",
                SelectedSurface::Cylindrical(_) => "cylindrical",
                SelectedSurface::Spherical(_) => "spherical",
                SelectedSurface::Conical(_) => "conical",
                SelectedSurface::Toroidal(_) => "toroidal",
            };
            eprintln!(
                "  Face {fi}: {stype} — created successfully ({} edges)",
                output.face_descriptors[fi].edge_indices.len(),
            );
        }

        if let Some(viz_sender) = viz {
            let mut overlay = crate::viz::VizOverlay::new();
            // Tessellate the OCCT face and show it translucent
            let face_shape = make_face.shape();
            overlay.shape_meshes.push(crate::viz::tessellate_shape(
                face_shape, [0.3, 0.8, 0.3, 0.4], [0.0, 1.0, 0.0, 1.0], 0.1, 0.5,
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
            let stype = match surface_type {
                SelectedSurface::Planar(_) => "planar",
                SelectedSurface::Cylindrical(_) => "cylindrical",
                SelectedSurface::Spherical(_) => "spherical",
                SelectedSurface::Conical(_) => "conical",
                SelectedSurface::Toroidal(_) => "toroidal",
            };
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


        make_faces.push(make_face);
    }

    if !config.quiet {
        eprintln!(
            "Stage 3.4 ({:.3}s): Created {}/{} OCCT faces",
            t.elapsed().as_secs_f64(),
            make_faces.len(),
            num_faces,
        );
    }

    output.make_faces = make_faces;
    output.concave_faces = concave_faces;

    // 3. Validate faces
    validate_faces(&output, config)?;

    // 4. Compare against reference STEP if --compare
    if config.compare_shape.is_some() {
        compare_faces_to_step(&output, config)?;
        compare_surface_orientations_to_step(&output, config);
    }

    Ok(output)
}


// ---------------------------------------------------------------------------
// Stage 3.5: Construct shells
// ---------------------------------------------------------------------------

/// Stitch OCCT faces into shells using BRepBuilderAPI_Sewing.
///
/// Takes all faces created in stage 3.4 and passes them to BRepBuilderAPI_Sewing,
/// which merges shared edges and produces one or more TopoDS_Shell objects.
fn construct_shells(
    mut output: Stage3Output,
    config: &Config,
    viz: Option<&crate::viz::VizSender>,
) -> Result<Stage3Output, Stage3Error> {
    let t = Instant::now();
    let num_faces = output.make_faces.len();
    if num_faces == 0 {
        return Err(Stage3Error::AdjacencyError("no faces to sew into shells".into()));
    }

    // Create sewing operator
    let sewing_tol = config.vertex_tolerance_mm;
    let mut sewing = b_rep_builder_api::Sewing::new_real(sewing_tol);

    // Add all faces to the sewing operator.
    for mf in output.make_faces.iter() {
        sewing.add(mf.face().as_shape());
    }

    // Perform sewing
    let progress = message::ProgressRange::new();
    sewing.perform(&progress);

    let sewn = sewing.sewed_shape();
    let sewn_type = sewn.shape_type();

    if config.verbose {
        eprintln!("  Sewing result shape type: {:?}", sewn_type);
        let n_free = sewing.nb_free_edges();
        let n_multi = sewing.nb_multiple_edges();
        let n_contig = sewing.nb_contigous_edges();
        eprintln!("  Sewing stats: {n_free} free edges, {n_multi} multiple edges, {n_contig} contiguous edges");
        for i in 1..=n_free {
            let free_edge = sewing.free_edge(i);
            let mut edge_explorer = top_exp::Explorer::new_shape_shapeenum2(
                free_edge.as_shape(), top_abs::ShapeEnum::Vertex, top_abs::ShapeEnum::Shape,
            );
            let mut pts: Vec<String> = Vec::new();
            while edge_explorer.more() {
                let vtx = topo_ds::vertex(edge_explorer.value());
                let pt = b_rep::Tool::pnt(vtx);
                pts.push(format!("({:.4},{:.4},{:.4})", pt.x(), pt.y(), pt.z()));
                edge_explorer.next();
            }
            eprintln!("    Free edge {i}: vertices {}", pts.join(" -> "));
        }
    }

    // Extract shells from the sewing result.
    // The result can be a Shell (single), a Solid, or a Compound containing shells.
    let mut shells: Vec<OwnedPtr<topo_ds::Shell>> = Vec::new();

    match sewn_type {
        top_abs::ShapeEnum::Shell => {
            // Single shell result — copy it
            let shell = topo_ds::shell_shape(sewn);
            shells.push(shell.to_owned());
        }
        top_abs::ShapeEnum::Face => {
            // Single face — wrap it in a shell (e.g., a full sphere with no edges)
            let mut shell = topo_ds::Shell::new();
            let builder = topo_ds::Builder::new();
            builder.make_shell(&mut shell);
            builder.add(shell.as_shape_mut(), sewn);
            shells.push(shell);
        }
        top_abs::ShapeEnum::Solid => {
            // Sewing produced a solid directly — extract shells from it
            let mut shell_explorer = top_exp::Explorer::new_shape_shapeenum2(
                sewn,
                top_abs::ShapeEnum::Shell,
                top_abs::ShapeEnum::Shape,
            );
            while shell_explorer.more() {
                let shell = topo_ds::shell_shape(shell_explorer.value());
                shells.push(shell.to_owned());
                shell_explorer.next();
            }
        }
        top_abs::ShapeEnum::Compound | top_abs::ShapeEnum::Compsolid => {
            // Multiple shells — iterate to find them; also check for loose faces
            let mut shell_explorer = top_exp::Explorer::new_shape_shapeenum2(
                sewn,
                top_abs::ShapeEnum::Shell,
                top_abs::ShapeEnum::Shape,
            );
            while shell_explorer.more() {
                let shell = topo_ds::shell_shape(shell_explorer.value());
                shells.push(shell.to_owned());
                shell_explorer.next();
            }
        }
        _ => {
            return Err(Stage3Error::AdjacencyError(format!(
                "unexpected sewing result shape type: {:?}",
                sewn_type
            )));
        }
    }

    if shells.is_empty() {
        return Err(Stage3Error::AdjacencyError(
            "sewing produced no shells".into(),
        ));
    }

    // Apply ShapeFix_Shell to fix pcurves and face orientations.
    // ShapeFix_Shell::Perform() calls ShapeFix_Face::Perform() on each face,
    // re-adding pcurves that sewing may have discarded when merging edges.
    let progress2 = message::ProgressRange::new();
    for shell in shells.iter_mut() {
        let mut fixer = shape_fix::Shell::new_shell(shell);
        fixer.perform(&progress2);
        *shell = fixer.shell();
    }

    // Fix face orientations in each shell.
    // ShapeFix_Shell::FixFaceOrientation uses BFS to propagate consistent
    // face orientations through shared edges.
    for shell in shells.iter_mut() {
        let mut fixer = shape_fix::Shell::new_shell(shell);
        fixer.fix_face_orientation(shell, true, false);
        *shell = fixer.shell();
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
            eprintln!("  Warning: shell {si} has orientation inconsistencies after fixing");
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
                "  Compare 3.6: WARNING solid {our_idx} volume differs by {rel_diff:.2e} (>1%)",
            );
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

    if !config.stage.at_least(3, 2) {
        return Ok(output);
    }

    // Stage 3.2: Detect tangency relationships along edges
    let output = detect_tangency(output, config);


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
