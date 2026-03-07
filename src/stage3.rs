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
    b_rep_builder_api, b_rep_check, geom, geom_api, gp, top_abs, top_exp, topo_ds, OwnedPtr,
};
use std::collections::{BTreeSet, HashMap, HashSet};
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
    // TODO: OCCT shell and solid objects — populated in stages 3.5/3.6
}

impl std::fmt::Debug for Stage3Output {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Stage3Output")
            .field("face_descriptors", &self.face_descriptors.len())
            .field("edges", &self.edges.len())
            .field("vertices", &self.vertices.len())
            .field("make_faces", &self.make_faces.len())
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
    Compare(Stage3CompareError),
}

impl Display for Stage3Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage3Error::NotImplemented(stage) => write!(f, "stage {stage} not yet implemented"),
            Stage3Error::AdjacencyError(msg) => write!(f, "stage 3.1 adjacency error: {msg}"),
            Stage3Error::Compare(e) => write!(f, "stage 3.1 compare: {e}"),
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
    }
}

/// Get the faces belonging to a selected surface.
fn surface_faces<'a>(surface: &SelectedSurface, output: &'a Stage2Output) -> &'a [usize] {
    match surface {
        SelectedSurface::Planar(idx) => &output.planar_hypotheses[*idx].faces,
        SelectedSurface::Cylindrical(idx) => &output.cylindrical_hypotheses[*idx].faces,
        SelectedSurface::Spherical(idx) => &output.spherical_hypotheses[*idx].faces,
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
            "Stage 3.1: Created {} OCCT surfaces, {} edges, {} vertices",
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
    }
}

/// Detect tangency relationships along edges.
///
/// For each ReconEdge, samples mesh boundary vertices and compares the surface
/// normals of the two adjacent surfaces. If all sampled normals agree within
/// a small angle threshold (2°), the edge is marked as tangent.
fn detect_tangency(mut output: Stage3Output, config: &Config) -> Stage3Output {
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
            "Stage 3.2: Tangency detection: {} of {} edges are tangent",
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

        // Ensure t_start < t_end (for non-periodic curves)
        if t_start < t_end {
            (t_start, t_end)
        } else if curve.is_periodic() {
            // For periodic curves, the order matters; keep as-is
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

    edge.curve_3d = Some(trimmed_handle);
    Ok(())
}

/// Select the intersection curve closest to the mesh boundary vertices.
/// Returns the 1-indexed line number.
fn select_closest_curve(
    int_ss: &geom_api::IntSS,
    boundary_vertices: &[usize],
    mesh: &ConnectedMesh,
) -> i32 {
    let nb_lines = int_ss.nb_lines();
    let mut best_idx = 1_i32;
    let mut best_total_dist = f64::MAX;

    // Sample boundary midpoints (edge centroids) for proximity test
    let sample_count = boundary_vertices.len().min(10);
    let step = if boundary_vertices.len() <= sample_count {
        1
    } else {
        boundary_vertices.len() / sample_count
    };
    let sample_indices: Vec<usize> = (0..boundary_vertices.len()).step_by(step).collect();

    for line_idx in 1..=nb_lines {
        let curve_handle = int_ss.line(line_idx);
        let mut total_dist = 0.0;

        for &si in &sample_indices {
            let vi = boundary_vertices[si];
            let v = &mesh.vertices[vi];
            let pt = gp::Pnt::new_real3(v.x, v.y, v.z);
            let projector = geom_api::ProjectPointOnCurve::new_pnt_handlegeomcurve(
                &pt, curve_handle,
            );
            if projector.nb_points() > 0 {
                total_dist += projector.lower_distance();
            } else {
                total_dist += f64::MAX / 2.0; // Penalize curves that can't project
            }
        }

        if total_dist < best_total_dist {
            best_total_dist = total_dist;
            best_idx = line_idx;
        }
    }

    best_idx
}

/// Project a 3D point onto a curve and return the curve parameter.
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

/// Compute edge curves for all ReconEdges.
fn compute_edge_curves_all(
    mut output: Stage3Output,
    config: &Config,
) -> Result<Stage3Output, Stage3Error> {
    let mut success_count = 0;
    let mut fail_count = 0;
    let total = output.edges.len();

    for ei in 0..total {
        // Split borrow: extract edge mutably while borrowing face_descriptors immutably
        let (before, rest) = output.edges.split_at_mut(ei);
        let (edge, _after) = rest.split_first_mut().unwrap();
        let _ = before; // suppress unused warning

        if edge.tangent {
            // TODO: Tangent edges need special handling (stage 3.3 tangent case)
            if config.verbose {
                eprintln!("  Edge {ei}: skipping tangent edge (not yet implemented)");
            }
            fail_count += 1;
            continue;
        }

        match compute_edge_curve(
            edge,
            &output.face_descriptors,
            &output.vertices,
            &output.stage2.mesh,
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
            "Stage 3.3: Computed {success_count}/{total} edge curves ({fail_count} failed)",
        );
    }

    if fail_count > 0 {
        return Err(Stage3Error::AdjacencyError(format!(
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
    for (i, w) in wires.iter_mut().enumerate() {
        if i != outer_idx {
            make_face.add(w.wire());
        }
    }

    Ok(make_face)
}

/// Create an OCCT face for a periodic surface (cylinder or sphere) using UV parameter bounds.
fn create_periodic_face(
    fi: usize,
    output: &Stage3Output,
    config: &Config,
) -> Result<OwnedPtr<b_rep_builder_api::MakeFace>, Stage3Error> {
    let fd = &output.face_descriptors[fi];
    let surface = &output.stage2.selected_surfaces[fd.selected_surface_idx];

    match surface {
        SelectedSurface::Cylindrical(idx) => {
            let hyp = &output.stage2.cylindrical_hypotheses[*idx];
            let ax = hyp.axis_direction;
            let ao = hyp.axis_origin;

            // Compute V bounds from edge boundary vertices
            let mut v_min = f64::MAX;
            let mut v_max = f64::MIN;

            for &ei in &fd.edge_indices {
                let edge = &output.edges[ei];
                for &vi in &edge.mesh_boundary_vertices {
                    let v = &output.stage2.mesh.vertices[vi];
                    let dp = [v.x - ao[0], v.y - ao[1], v.z - ao[2]];
                    let v_param = dp[0] * ax[0] + dp[1] * ax[1] + dp[2] * ax[2];
                    v_min = v_min.min(v_param);
                    v_max = v_max.max(v_param);
                }
            }

            if v_min >= v_max {
                return Err(Stage3Error::AdjacencyError(format!(
                    "cylinder face {fi} has no valid V parameter range"
                )));
            }

            let u_min = 0.0;
            let u_max = 2.0 * std::f64::consts::PI;

            let make_face = b_rep_builder_api::MakeFace::new_handlegeomsurface_real5(
                &fd.surface,
                u_min,
                u_max,
                v_min,
                v_max,
                config.vertex_tolerance_mm,
            );

            if !make_face.is_done() {
                return Err(Stage3Error::AdjacencyError(format!(
                    "MakeFace failed for cylindrical face {fi}: {:?}",
                    make_face.error(),
                )));
            }

            Ok(make_face)
        }
        SelectedSurface::Spherical(idx) => {
            let hyp = &output.stage2.spherical_hypotheses[*idx];

            if fd.edge_indices.is_empty() {
                // Full sphere — use natural bounds
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

                return Ok(make_face);
            }

            // Compute V bounds from edge boundary vertices and mesh face centroids
            // V parameter on sphere: V = asin((P - center) · z_dir / radius)
            // Our Ax3 has z_dir = [0, 0, 1]
            let center = hyp.center;
            let mut v_min = f64::MAX;
            let mut v_max = f64::MIN;

            // Sample from edge boundary vertices
            for &ei in &fd.edge_indices {
                let edge = &output.edges[ei];
                for &vi in &edge.mesh_boundary_vertices {
                    let v = &output.stage2.mesh.vertices[vi];
                    let dp = [v.x - center[0], v.y - center[1], v.z - center[2]];
                    let len = (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]).sqrt();
                    if len < 1e-15 {
                        continue;
                    }
                    let v_param = (dp[2] / len).clamp(-1.0, 1.0).asin();
                    v_min = v_min.min(v_param);
                    v_max = v_max.max(v_param);
                }
            }

            // Extend V range using mesh face centroids to detect pole coverage
            let faces = surface_faces(surface, &output.stage2);
            for &mfi in faces {
                let c = compute_mesh_face_centroid(mfi, &output.stage2.mesh);
                let dp = [c[0] - center[0], c[1] - center[1], c[2] - center[2]];
                let len = (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]).sqrt();
                if len < 1e-15 {
                    continue;
                }
                let v_param = (dp[2] / len).clamp(-1.0, 1.0).asin();
                v_min = v_min.min(v_param);
                v_max = v_max.max(v_param);
            }

            // Snap to poles if close (within ~5.7°)
            let pole_snap = 0.1; // radians
            if v_max > std::f64::consts::FRAC_PI_2 - pole_snap {
                v_max = std::f64::consts::FRAC_PI_2;
            }
            if v_min < -std::f64::consts::FRAC_PI_2 + pole_snap {
                v_min = -std::f64::consts::FRAC_PI_2;
            }

            let u_min = 0.0;
            let u_max = 2.0 * std::f64::consts::PI;

            let make_face = b_rep_builder_api::MakeFace::new_handlegeomsurface_real5(
                &fd.surface,
                u_min,
                u_max,
                v_min,
                v_max,
                config.vertex_tolerance_mm,
            );

            if !make_face.is_done() {
                return Err(Stage3Error::AdjacencyError(format!(
                    "MakeFace failed for spherical face {fi}: {:?}",
                    make_face.error(),
                )));
            }

            Ok(make_face)
        }
        _ => unreachable!("create_periodic_face called for non-periodic surface"),
    }
}

/// Validate created faces using BRepCheck_Analyzer.
fn validate_faces(output: &Stage3Output, config: &Config) -> Result<(), Stage3Error> {
    let mut invalid_count = 0;

    for (fi, mf) in output.make_faces.iter().enumerate() {
        let face = mf.face();
        let analyzer = b_rep_check::Analyzer::new_shape(face.as_shape());
        if !analyzer.is_valid() {
            invalid_count += 1;
            if config.verbose {
                eprintln!("  Face {fi}: BRepCheck_Analyzer reports invalid");
            }
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

/// Compare stage 3.4 faces against reference STEP shape.
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
            check_type: "face",
            element_index: 0,
            max_distance: max_dist,
            tolerance: config.surface_tolerance_mm,
        });
    }

    Ok(())
}

/// Create OCCT faces for all face descriptors (stage 3.4 main entry).
fn create_occt_faces_all(
    mut output: Stage3Output,
    config: &Config,
) -> Result<Stage3Output, Stage3Error> {
    let num_faces = output.face_descriptors.len();

    // 1. Create TopoDS_Edge for each ReconEdge
    let mut topo_edges: Vec<OwnedPtr<b_rep_builder_api::MakeEdge>> = Vec::with_capacity(output.edges.len());
    for (ei, edge) in output.edges.iter().enumerate() {
        let curve = edge.curve_3d.as_ref().unwrap_or_else(|| {
            panic!("edge {ei} has no 3D curve — stage 3.3 must run first")
        });
        let make_edge = b_rep_builder_api::MakeEdge::new_handlegeomcurve(curve);
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

    for fi in 0..num_faces {
        let surface_type = &output.stage2.selected_surfaces[output.face_descriptors[fi].selected_surface_idx];
        let is_periodic = matches!(
            surface_type,
            SelectedSurface::Cylindrical(_) | SelectedSurface::Spherical(_)
        );

        let make_face = if is_periodic {
            create_periodic_face(fi, &output, config)?
        } else {
            create_planar_face(fi, &output, &mut topo_edges, config)?
        };

        if config.verbose {
            let stype = match surface_type {
                SelectedSurface::Planar(_) => "planar",
                SelectedSurface::Cylindrical(_) => "cylindrical",
                SelectedSurface::Spherical(_) => "spherical",
            };
            eprintln!(
                "  Face {fi}: {stype} — created successfully ({} edges)",
                output.face_descriptors[fi].edge_indices.len(),
            );
        }

        make_faces.push(make_face);
    }

    if !config.quiet {
        eprintln!(
            "Stage 3.4: Created {}/{} OCCT faces",
            make_faces.len(),
            num_faces,
        );
    }

    output.make_faces = make_faces;

    // 3. Validate faces
    validate_faces(&output, config)?;

    // 4. Compare against reference STEP if --compare
    if config.compare_shape.is_some() {
        compare_faces_to_step(&output, config)?;
    }

    Ok(output)
}

// Stage 3 entry point
// ---------------------------------------------------------------------------

/// Run stage 3: reconstruct B-Rep surfaces, edges, and topology from fitted surfaces.
pub fn stage3(config: &Config, input: Stage2Output) -> Result<Stage3Output, Stage3Error> {
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
    let output = compute_edge_curves_all(output, config)?;

    if !config.stage.at_least(3, 4) {
        return Ok(output);
    }

    // Stage 3.4: Create OCCT faces
    let output = create_occt_faces_all(output, config)?;

    if !config.stage.at_least(3, 5) {
        return Ok(output);
    }

    // Stage 3.5: Construct shells
    // TODO: BRepBuilderAPI_Sewing to stitch faces
    if !config.quiet {
        eprintln!("Stage 3.5: Construct shells (not yet implemented)");
    }

    if !config.stage.at_least(3, 6) {
        return Err(Stage3Error::NotImplemented("3.5".into()));
    }

    // Stage 3.6: Construct solids
    // TODO: BRepBuilderAPI_MakeSolid from shells, classify voids
    if !config.quiet {
        eprintln!("Stage 3.6: Construct solids (not yet implemented)");
    }

    Err(Stage3Error::NotImplemented("3".into()))
}
