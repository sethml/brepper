//! Stage 1: Mesh Input & Preprocessing
//!
//! 1.1: Read STL file into ConnectedMesh with welded vertices.
//! 1.2: Validate mesh topology, compute normals and neighbors.

use crate::config::Config;
use opencascade_sys::{b_rep_builder_api, b_rep_extrema, extrema, gp, message, rw_stl};
use std::collections::{HashMap, HashSet, VecDeque};
use std::time::Instant;
use std::error::Error;
use std::fmt::{Display, Formatter};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

pub const UNDEDUCED_PLANAR_HYPOTHESIS: i32 = -2;
pub const NO_HYPOTHESIS: i32 = -1;
pub const UNDEDUCED_CYLINDRICAL_HYPOTHESIS: i32 = -2;
pub const UNDEDUCED_SPHERICAL_HYPOTHESIS: i32 = -2;
pub const UNDEDUCED_CONICAL_HYPOTHESIS: i32 = -2;
pub const UNDEDUCED_TOROIDAL_HYPOTHESIS: i32 = -2;

// ---------------------------------------------------------------------------
// Data structures — Stage 1 output
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshVertex {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl MeshVertex {
    pub fn from_xyz(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

#[derive(Debug, Clone)]
pub struct MeshFace {
    pub vertex_count: u8,
    pub vertex_indices: [usize; 4],
    pub neighbors: [i32; 4],
    pub normal: Option<[f64; 3]>,
    pub planar_hypothesis: i32,
    pub cylindrical_hypothesis: i32,
    pub spherical_hypothesis: i32,
    pub conical_hypothesis: i32,
    pub toroidal_hypothesis: i32,
}

#[derive(Debug, Clone, Default)]
pub struct MeshValidationStats {
    pub mesh_faces: usize,
    pub mesh_vertices: usize,
    pub mesh_edges_open: usize,
    pub mesh_edges_non_manifold: usize,
    pub mesh_edges_inconsistent_orientation: usize,
    pub mesh_faces_degenerate: usize,
    pub connected_shells: usize,
    pub solids: usize,
    pub voids_within_solids: usize,
    pub mesh_quads: usize,
}

/// The output of Stage 1: a connected mesh with welded vertices, computed normals,
/// edge neighbors, and validation statistics.
#[derive(Debug, Clone, Default)]
pub struct ConnectedMesh {
    pub vertices: Vec<MeshVertex>,
    pub faces: Vec<MeshFace>,
    pub stats: MeshValidationStats,
}

// ---------------------------------------------------------------------------
// STL reading (Stage 1.1)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct VertexWeldOptions {
    pub tolerance: f64,
}

impl Default for VertexWeldOptions {
    fn default() -> Self {
        Self { tolerance: 1.0e-9 }
    }
}

#[derive(Debug)]
pub enum MeshReadError {
    InvalidTolerance(f64),
    FileNotFound(String),
    EmptyMesh,
}

impl Display for MeshReadError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            MeshReadError::InvalidTolerance(tol) => {
                write!(f, "invalid weld tolerance: {tol}, expected > 0")
            }
            MeshReadError::FileNotFound(path) => {
                write!(f, "STL file not found: {path}")
            }
            MeshReadError::EmptyMesh => write!(f, "STL file did not contain any triangles"),
        }
    }
}

impl Error for MeshReadError {}

#[derive(Debug)]
pub enum MeshValidationError {
    DegenerateFaces { face_indices: Vec<usize> },
    NonManifoldEdges { edge_count: usize },
    InconsistentOrientation { edge_count: usize },
    InvertedNormals { shell_index: usize, signed_volume: f64 },
}

impl Display for MeshValidationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            MeshValidationError::DegenerateFaces { face_indices } => write!(
                f,
                "mesh contains {} degenerate face(s): {:?}",
                face_indices.len(),
                face_indices
            ),
            MeshValidationError::NonManifoldEdges { edge_count } => write!(
                f,
                "mesh contains {edge_count} non-manifold edge(s) with >2 incident faces"
            ),
            MeshValidationError::InconsistentOrientation { edge_count } => {
                write!(f, "mesh contains {edge_count} edge(s) with flipped orientation")
            }
            MeshValidationError::InvertedNormals { shell_index, signed_volume } => {
                write!(f, "shell {shell_index} has inverted face normals (signed volume {signed_volume:.6e}; expected positive)")
            }
        }
    }
}

impl Error for MeshValidationError {}

#[derive(Debug)]
pub enum Stage1Error {
    Read(MeshReadError),
    Validation(MeshValidationError),
    Compare(CompareError),
}

impl Display for Stage1Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage1Error::Read(e) => write!(f, "stage 1.1 (read STL): {e}"),
            Stage1Error::Validation(e) => write!(f, "stage 1.2 (validation): {e}"),
            Stage1Error::Compare(e) => write!(f, "stage 1.2 (compare): {e}"),
        }
    }
}

impl Error for Stage1Error {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Stage1Error::Read(e) => Some(e),
            Stage1Error::Validation(e) => Some(e),
            Stage1Error::Compare(e) => Some(e),
        }
    }
}

impl From<MeshReadError> for Stage1Error {
    fn from(e: MeshReadError) -> Self {
        Stage1Error::Read(e)
    }
}

impl From<MeshValidationError> for Stage1Error {
    fn from(e: MeshValidationError) -> Self {
        Stage1Error::Validation(e)
    }
}

impl From<CompareError> for Stage1Error {
    fn from(e: CompareError) -> Self {
        Stage1Error::Compare(e)
    }
}

// ---------------------------------------------------------------------------
// Vertex welding internals
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct CellKey {
    x: i64,
    y: i64,
    z: i64,
}

impl CellKey {
    fn from_vertex(v: &MeshVertex, tolerance: f64) -> Self {
        Self {
            x: (v.x / tolerance).round() as i64,
            y: (v.y / tolerance).round() as i64,
            z: (v.z / tolerance).round() as i64,
        }
    }
}

#[derive(Debug, Default)]
struct VertexWelder {
    vertices: Vec<MeshVertex>,
    buckets: HashMap<CellKey, Vec<usize>>,
}

impl VertexWelder {
    fn get_or_insert(&mut self, vertex: MeshVertex, tolerance: f64) -> usize {
        let base = CellKey::from_vertex(&vertex, tolerance);
        let tol2 = tolerance * tolerance;

        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let key = CellKey {
                        x: base.x + dx,
                        y: base.y + dy,
                        z: base.z + dz,
                    };
                    if let Some(indices) = self.buckets.get(&key) {
                        for idx in indices {
                            let candidate = self.vertices[*idx];
                            let ddx = candidate.x - vertex.x;
                            let ddy = candidate.y - vertex.y;
                            let ddz = candidate.z - vertex.z;
                            if (ddx * ddx) + (ddy * ddy) + (ddz * ddz) <= tol2 {
                                return *idx;
                            }
                        }
                    }
                }
            }
        }

        let idx = self.vertices.len();
        self.vertices.push(vertex);
        self.buckets.entry(base).or_default().push(idx);
        idx
    }
}

pub fn read_connected_mesh_from_stl(
    stl_path: &str,
    options: VertexWeldOptions,
) -> Result<ConnectedMesh, MeshReadError> {
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        return Err(MeshReadError::InvalidTolerance(options.tolerance));
    }

    if !std::path::Path::new(stl_path).exists() {
        return Err(MeshReadError::FileNotFound(stl_path.to_string()));
    }

    let progress = message::ProgressRange::new();
    let tri_handle = rw_stl::read_file_charptr_progressrange_2(stl_path, &progress);
    let tri = tri_handle.get();
    let num_nodes = tri.nb_nodes();
    let num_triangles = tri.nb_triangles();

    if num_nodes <= 0 || num_triangles <= 0 {
        return Err(MeshReadError::EmptyMesh);
    }

    let mut welder = VertexWelder::default();
    let mut node_to_welded_idx = vec![usize::MAX; (num_nodes as usize) + 1];
    let mut faces = Vec::with_capacity(num_triangles as usize);

    for tri_idx in 1..=num_triangles {
        let triangle = tri.triangle(tri_idx);

        let mut face_vertex_indices = [0_usize; 4];
        for (slot, fvi) in face_vertex_indices.iter_mut().take(3).enumerate() {
            let node_index = triangle.value((slot + 1) as i32);
            let node_slot = node_index as usize;

            let welded_idx = if node_to_welded_idx[node_slot] != usize::MAX {
                node_to_welded_idx[node_slot]
            } else {
                let node = tri.node(node_index);
                let v = MeshVertex::from_xyz(node.x(), node.y(), node.z());
                let idx = welder.get_or_insert(v, options.tolerance);
                node_to_welded_idx[node_slot] = idx;
                idx
            };

            *fvi = welded_idx;
        }

        faces.push(MeshFace {
            vertex_count: 3,
            vertex_indices: face_vertex_indices,
            neighbors: [-1, -1, -1, -1],
            normal: None,
            planar_hypothesis: UNDEDUCED_PLANAR_HYPOTHESIS,
            cylindrical_hypothesis: UNDEDUCED_CYLINDRICAL_HYPOTHESIS,
            spherical_hypothesis: UNDEDUCED_SPHERICAL_HYPOTHESIS,
            conical_hypothesis: UNDEDUCED_CONICAL_HYPOTHESIS,
            toroidal_hypothesis: UNDEDUCED_TOROIDAL_HYPOTHESIS,
        });
    }

    Ok(ConnectedMesh {
        vertices: welder.vertices,
        faces,
        stats: MeshValidationStats::default(),
    })
}

// ---------------------------------------------------------------------------
// Mesh validation (Stage 1.2)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct EdgeKey {
    lo: usize,
    hi: usize,
}

impl EdgeKey {
    fn new(a: usize, b: usize) -> Self {
        if a <= b {
            Self { lo: a, hi: b }
        } else {
            Self { lo: b, hi: a }
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct EdgeUse {
    face_idx: usize,
    edge_idx: usize,
    start: usize,
    end: usize,
}

impl ConnectedMesh {
    pub fn validate_and_populate_topology(&mut self) -> Result<(), MeshValidationError> {
        let mut edge_uses: HashMap<EdgeKey, Vec<EdgeUse>> = HashMap::new();
        let mut degenerate_faces = Vec::new();

        self.stats.mesh_faces = self.faces.len();
        self.stats.mesh_vertices = self.vertices.len();
        self.stats.mesh_edges_open = 0;
        self.stats.mesh_edges_non_manifold = 0;
        self.stats.mesh_edges_inconsistent_orientation = 0;
        self.stats.mesh_faces_degenerate = 0;
        self.stats.connected_shells = 0;
        self.stats.solids = 0;
        self.stats.voids_within_solids = 0;

        for face_idx in 0..self.faces.len() {
            let face = &mut self.faces[face_idx];
            face.neighbors = [-1, -1, -1, -1];

            let vertex_count = face.vertex_count as usize;
            if !(3..=4).contains(&vertex_count) {
                degenerate_faces.push(face_idx);
                continue;
            }

            let normal = compute_face_normal(face, &self.vertices);
            if let Some(dir) = normal {
                face.normal = Some(dir);
            } else {
                face.normal = None;
                degenerate_faces.push(face_idx);
            }

            for edge_idx in 0..vertex_count {
                let start = face.vertex_indices[edge_idx];
                let end = face.vertex_indices[(edge_idx + 1) % vertex_count];
                edge_uses
                    .entry(EdgeKey::new(start, end))
                    .or_default()
                    .push(EdgeUse {
                        face_idx,
                        edge_idx,
                        start,
                        end,
                    });
            }
        }

        let mut adjacency: Vec<HashSet<usize>> = vec![HashSet::new(); self.faces.len()];
        for uses in edge_uses.values() {
            if uses.len() == 1 {
                self.stats.mesh_edges_open += 1;
                continue;
            }

            if uses.len() > 2 {
                self.stats.mesh_edges_non_manifold += 1;
            }

            for i in 0..uses.len() {
                for j in (i + 1)..uses.len() {
                    let a = uses[i].face_idx;
                    let b = uses[j].face_idx;
                    adjacency[a].insert(b);
                    adjacency[b].insert(a);
                }
            }

            if uses.len() == 2 {
                let first = uses[0];
                let second = uses[1];

                self.faces[first.face_idx].neighbors[first.edge_idx] = second.face_idx as i32;
                self.faces[second.face_idx].neighbors[second.edge_idx] = first.face_idx as i32;

                let opposite_orientation = first.start == second.end && first.end == second.start;
                if !opposite_orientation {
                    self.stats.mesh_edges_inconsistent_orientation += 1;
                }
            }
        }

        self.stats.mesh_faces_degenerate = degenerate_faces.len();

        let mut face_component = vec![usize::MAX; self.faces.len()];
        let mut component_count = 0_usize;
        for face_idx in 0..self.faces.len() {
            if face_component[face_idx] != usize::MAX {
                continue;
            }

            let mut queue = VecDeque::new();
            queue.push_back(face_idx);
            face_component[face_idx] = component_count;

            while let Some(cur) = queue.pop_front() {
                for next in &adjacency[cur] {
                    if face_component[*next] == usize::MAX {
                        face_component[*next] = component_count;
                        queue.push_back(*next);
                    }
                }
            }

            component_count += 1;
        }
        self.stats.connected_shells = component_count;

        let mut component_closed = vec![true; component_count];
        for uses in edge_uses.values() {
            let mut comps = HashSet::new();
            for use_ in uses {
                comps.insert(face_component[use_.face_idx]);
            }

            for comp in comps {
                if uses.len() != 2 {
                    component_closed[comp] = false;
                }
            }
        }
        self.stats.solids = component_closed.iter().filter(|closed| **closed).count();

        // Self-intersection checks are intentionally deferred to a later stage where we have
        // richer geometric predicates and can avoid expensive O(N^2) triangle checks here.

        // Check face orientation via signed volume for each closed shell.
        // A positive signed volume means face normals point outward (correct).
        // A negative signed volume means face normals point inward (inverted).
        for (comp, is_closed) in component_closed.iter().enumerate() {
            if !is_closed {
                continue; // Only check closed shells
            }
            let mut signed_volume = 0.0_f64;
            for (fi, &comp_id) in face_component.iter().enumerate() {
                if comp_id != comp {
                    continue;
                }
                let face = &self.faces[fi];
                let vc = face.vertex_count as usize;
                // Signed volume contribution: sum of signed tetrahedra with origin.
                // For triangles: V += (v0 · (v1 × v2)) / 6
                // For quads: split into two triangles (v0,v1,v2) and (v0,v2,v3)
                let v0 = &self.vertices[face.vertex_indices[0]];
                let v1 = &self.vertices[face.vertex_indices[1]];
                let v2 = &self.vertices[face.vertex_indices[2]];
                // v1 × v2
                let cross_x = v1.y * v2.z - v1.z * v2.y;
                let cross_y = v1.z * v2.x - v1.x * v2.z;
                let cross_z = v1.x * v2.y - v1.y * v2.x;
                signed_volume += v0.x * cross_x + v0.y * cross_y + v0.z * cross_z;
                if vc == 4 {
                    let v3 = &self.vertices[face.vertex_indices[3]];
                    let cross_x = v2.y * v3.z - v2.z * v3.y;
                    let cross_y = v2.z * v3.x - v2.x * v3.z;
                    let cross_z = v2.x * v3.y - v2.y * v3.x;
                    signed_volume += v0.x * cross_x + v0.y * cross_y + v0.z * cross_z;
                }
            }
            signed_volume /= 6.0;
            if signed_volume < 0.0 {
                return Err(MeshValidationError::InvertedNormals {
                    shell_index: comp,
                    signed_volume,
                });
            }
        }

        if !degenerate_faces.is_empty() {
            return Err(MeshValidationError::DegenerateFaces {
                face_indices: degenerate_faces,
            });
        }

        if self.stats.mesh_edges_non_manifold > 0 {
            return Err(MeshValidationError::NonManifoldEdges {
                edge_count: self.stats.mesh_edges_non_manifold,
            });
        }

        if self.stats.mesh_edges_inconsistent_orientation > 0 {
            return Err(MeshValidationError::InconsistentOrientation {
                edge_count: self.stats.mesh_edges_inconsistent_orientation,
            });
        }

        Ok(())
    }
}

fn compute_face_normal(face: &MeshFace, vertices: &[MeshVertex]) -> Option<[f64; 3]> {
    let vertex_count = face.vertex_count as usize;
    if !(3..=4).contains(&vertex_count) {
        return None;
    }

    // Newell normal handles both triangles and quads with consistent winding.
    let mut nx = 0.0_f64;
    let mut ny = 0.0_f64;
    let mut nz = 0.0_f64;

    for i in 0..vertex_count {
        let curr = vertices[face.vertex_indices[i]];
        let next = vertices[face.vertex_indices[(i + 1) % vertex_count]];

        nx += (curr.y - next.y) * (curr.z + next.z);
        ny += (curr.z - next.z) * (curr.x + next.x);
        nz += (curr.x - next.x) * (curr.y + next.y);
    }

    let len2 = nx * nx + ny * ny + nz * nz;
    if len2 <= 1.0e-24 {
        return None;
    }

    let inv_len = len2.sqrt().recip();
    Some([nx * inv_len, ny * inv_len, nz * inv_len])
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Coplanar triangle fusion (Stage 1.3)
// ---------------------------------------------------------------------------

/// Check whether the quad `[v0, v1, v2, v3]` is strictly convex, given the face normal.
fn is_convex_quad(vertices: &[MeshVertex], quad_verts: &[usize; 4], normal: &[f64; 3]) -> bool {
    for i in 0..4 {
        let v0 = &vertices[quad_verts[i]];
        let v1 = &vertices[quad_verts[(i + 1) % 4]];
        let v2 = &vertices[quad_verts[(i + 2) % 4]];

        let e1 = [v1.x - v0.x, v1.y - v0.y, v1.z - v0.z];
        let e2 = [v2.x - v1.x, v2.y - v1.y, v2.z - v1.z];

        // Cross product e1 × e2
        let cx = e1[1] * e2[2] - e1[2] * e2[1];
        let cy = e1[2] * e2[0] - e1[0] * e2[2];
        let cz = e1[0] * e2[1] - e1[1] * e2[0];

        // For convex quad, all cross products should agree with the normal direction
        let dot = cx * normal[0] + cy * normal[1] + cz * normal[2];
        if dot <= 0.0 {
            return false;
        }
    }
    true
}

/// Stage 1.3: Fuse adjacent coplanar triangle pairs into quads.
///
/// For each manifold edge shared by two triangles, if both triangles are coplanar
/// (all 4 vertices within `coplanar_tol` of the shared plane) and their merged shape
/// forms a convex quad, merge them into a single quad face. Each triangle can
/// participate in at most one merge (greedy, first valid merge wins).
pub fn fuse_coplanar_triangles(mesh: &mut ConnectedMesh, vertex_tolerance: f64) {
    let coplanar_tol = vertex_tolerance;
    let num_faces = mesh.faces.len();
    let mut merged = vec![false; num_faces];

    // Phase 1: Identify merge pairs.
    // Each entry: (face1_idx, edge1_idx, face2_idx, edge2_idx)
    let mut merge_pairs: Vec<(usize, usize, usize, usize)> = Vec::new();

    for fi in 0..num_faces {
        if merged[fi] {
            continue;
        }
        if mesh.faces[fi].vertex_count != 3 {
            continue;
        }
        let norm = match mesh.faces[fi].normal {
            Some(n) => n,
            None => continue,
        };

        for ei in 0..3_usize {
            let ni = mesh.faces[fi].neighbors[ei];
            if ni < 0 {
                continue;
            }
            let ni = ni as usize;
            // Process each edge only once and skip already-merged faces
            if ni <= fi || merged[ni] {
                continue;
            }
            if mesh.faces[ni].vertex_count != 3 {
                continue;
            }

            // Find which edge in ni points back to fi
            let mut eni = usize::MAX;
            for k in 0..3 {
                if mesh.faces[ni].neighbors[k] == fi as i32 {
                    eni = k;
                    break;
                }
            }
            if eni == usize::MAX {
                continue;
            }

            // Identify the 4 unique vertices
            let s0 = mesh.faces[fi].vertex_indices[ei];
            let s1 = mesh.faces[fi].vertex_indices[(ei + 1) % 3];
            let a = mesh.faces[fi].vertex_indices[(ei + 2) % 3];
            let b = mesh.faces[ni].vertex_indices[(eni + 2) % 3];

            // Coplanarity check: vertex b must be within tolerance of face fi's plane.
            // (Vertices s0, s1, a are on the plane by construction.)
            let v_s0 = &mesh.vertices[s0];
            let plane_d = norm[0] * v_s0.x + norm[1] * v_s0.y + norm[2] * v_s0.z;
            let v_b = &mesh.vertices[b];
            let dist_b = norm[0] * v_b.x + norm[1] * v_b.y + norm[2] * v_b.z - plane_d;
            if dist_b.abs() > coplanar_tol {
                continue;
            }

            // Convexity check: quad [s1, a, s0, b] must be strictly convex.
            let quad_verts = [s1, a, s0, b];
            if !is_convex_quad(&mesh.vertices, &quad_verts, &norm) {
                continue;
            }

            // Valid merge — mark both faces
            merged[fi] = true;
            merged[ni] = true;
            merge_pairs.push((fi, ei, ni, eni));
            break; // fi can only merge once
        }
    }

    if merge_pairs.is_empty() {
        return;
    }

    // Phase 2: Apply merges — transform face1 into a quad, mark face2 as deleted.
    for &(fi, ei, ni, eni) in &merge_pairs {
        let s0 = mesh.faces[fi].vertex_indices[ei];
        let s1 = mesh.faces[fi].vertex_indices[(ei + 1) % 3];
        let a = mesh.faces[fi].vertex_indices[(ei + 2) % 3];
        let b = mesh.faces[ni].vertex_indices[(eni + 2) % 3];

        // Quad vertex order [s1, a, s0, b] follows the right-hand rule,
        // consistent with face1's winding.
        let quad_indices = [s1, a, s0, b];

        // Quad neighbors:
        //   Edge 0 (s1→a): was face1's edge (e1+1)%3
        //   Edge 1 (a→s0): was face1's edge (e1+2)%3
        //   Edge 2 (s0→b): was face2's edge (e2+1)%3
        //   Edge 3 (b→s1): was face2's edge (e2+2)%3
        let n0 = mesh.faces[fi].neighbors[(ei + 1) % 3];
        let n1 = mesh.faces[fi].neighbors[(ei + 2) % 3];
        let n2 = mesh.faces[ni].neighbors[(eni + 1) % 3];
        let n3 = mesh.faces[ni].neighbors[(eni + 2) % 3];

        // Transform fi into the quad
        mesh.faces[fi].vertex_count = 4;
        mesh.faces[fi].vertex_indices = quad_indices;
        mesh.faces[fi].neighbors = [n0, n1, n2, n3];
        mesh.faces[fi].normal = compute_face_normal(&mesh.faces[fi], &mesh.vertices);

        // Update ni's former neighbors to point to fi instead of ni
        for n_ref in [n2, n3] {
            if n_ref >= 0 {
                let n_idx = n_ref as usize;
                for k in 0..mesh.faces[n_idx].vertex_count as usize {
                    if mesh.faces[n_idx].neighbors[k] == ni as i32 {
                        mesh.faces[n_idx].neighbors[k] = fi as i32;
                    }
                }
            }
        }

        // Mark ni as deleted
        mesh.faces[ni].vertex_count = 0;
    }

    // Phase 3: Compact — remove deleted faces and remap neighbor indices.
    let old_len = mesh.faces.len();
    let mut new_index: Vec<i32> = vec![-1; old_len];
    let mut new_faces: Vec<MeshFace> = Vec::with_capacity(old_len - merge_pairs.len());

    for (old_i, face) in mesh.faces.iter().enumerate() {
        if face.vertex_count == 0 {
            continue;
        }
        new_index[old_i] = new_faces.len() as i32;
        new_faces.push(face.clone());
    }

    // Remap neighbor references
    for face in &mut new_faces {
        for k in 0..face.vertex_count as usize {
            if face.neighbors[k] >= 0 {
                let old_ni = face.neighbors[k] as usize;
                let mapped = new_index[old_ni];
                debug_assert!(mapped >= 0, "neighbor points to deleted face");
                face.neighbors[k] = mapped;
            }
        }
    }

    let quads_formed = merge_pairs.len();
    mesh.faces = new_faces;
    mesh.stats.mesh_faces = mesh.faces.len();
    mesh.stats.mesh_quads = quads_formed;
}

// Compare mesh against STEP surfaces
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct CompareError {
    pub vertex_max_distance: f64,
    pub vertex_max_index: usize,
    pub vertex_tolerance: f64,
    pub centroid_max_distance: f64,
    pub centroid_max_face_index: usize,
    pub centroid_tolerance: f64,
}

impl Display for CompareError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut parts = Vec::new();
        if self.vertex_max_distance > self.vertex_tolerance {
            parts.push(format!(
                "vertex {} is {:.6e} mm from nearest STEP surface (tolerance: {:.6e} mm)",
                self.vertex_max_index, self.vertex_max_distance, self.vertex_tolerance
            ));
        }
        if self.centroid_max_distance > self.centroid_tolerance {
            parts.push(format!(
                "face {} centroid is {:.6e} mm from nearest STEP surface (tolerance: {:.6e} mm)",
                self.centroid_max_face_index, self.centroid_max_distance, self.centroid_tolerance
            ));
        }
        write!(f, "comparison failed: {}", parts.join("; "))
    }
}

impl Error for CompareError {}

/// Compute the minimum distance from a point to the bounded compare shape
/// using BRepExtrema_DistShapeShape, which respects face trimming boundaries.
pub fn min_distance_to_shape(
    pt: &gp::Pnt,
    compare_shape: &opencascade_sys::topo_ds::Shape,
) -> f64 {
    let mut vertex = b_rep_builder_api::MakeVertex::new_pnt(pt);
    let vtx_shape = vertex.vertex();
    let progress = message::ProgressRange::new();
    let dist_calc = b_rep_extrema::DistShapeShape::new_shape2_extflag_extalgo_progressrange(
        vtx_shape.as_shape(),
        compare_shape,
        0, // F (unused/obsolete)
        extrema::ExtAlgo::Grad,
        &progress,
    );
    if dist_calc.is_done() && dist_calc.nb_solution() > 0 {
        dist_calc.value()
    } else {
        f64::MAX
    }
}

/// Check mesh vertices and face centroids against comparison STEP shape.
fn compare_mesh_to_step(mesh: &ConnectedMesh, config: &Config) -> Result<(), CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    // Check all vertices
    let mut vtx_max_dist = 0.0_f64;
    let mut vtx_max_idx = 0_usize;
    for (i, v) in mesh.vertices.iter().enumerate() {
        let pt = gp::Pnt::new_real3(v.x, v.y, v.z);
        let d = min_distance_to_shape(&pt, compare_shape);
        if d > vtx_max_dist {
            vtx_max_dist = d;
            vtx_max_idx = i;
        }
    }

    // Check all face centroids
    let mut ctr_max_dist = 0.0_f64;
    let mut ctr_max_face = 0_usize;
    for (fi, face) in mesh.faces.iter().enumerate() {
        let n = face.vertex_count as usize;
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        let mut cz = 0.0_f64;
        for vi in 0..n {
            let v = &mesh.vertices[face.vertex_indices[vi]];
            cx += v.x;
            cy += v.y;
            cz += v.z;
        }
        let inv_n = 1.0 / n as f64;
        let pt = gp::Pnt::new_real3(cx * inv_n, cy * inv_n, cz * inv_n);
        let d = min_distance_to_shape(&pt, compare_shape);
        if d > ctr_max_dist {
            ctr_max_dist = d;
            ctr_max_face = fi;
        }
    }

    if !config.quiet {
        eprintln!(
            "  Compare: vertex max dist {:.6e} mm, centroid max dist {:.6e} mm",
            vtx_max_dist, ctr_max_dist
        );
    }

    let vtx_fail = vtx_max_dist > config.vertex_tolerance_mm;
    let ctr_fail = ctr_max_dist > config.surface_tolerance_mm;

    if vtx_fail || ctr_fail {
        return Err(CompareError {
            vertex_max_distance: vtx_max_dist,
            vertex_max_index: vtx_max_idx,
            vertex_tolerance: config.vertex_tolerance_mm,
            centroid_max_distance: ctr_max_dist,
            centroid_max_face_index: ctr_max_face,
            centroid_tolerance: config.surface_tolerance_mm,
        });
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Stage 1 entry point
// ---------------------------------------------------------------------------
/// Run stage 1: read STL (1.1), validate mesh (1.2), and fuse coplanar triangles (1.3).
pub fn stage1(config: &Config) -> Result<ConnectedMesh, Stage1Error> {
    // Stage 1.1: Read STL
    let t = Instant::now();
    let mut mesh = read_connected_mesh_from_stl(
        &config.input_stl,
        VertexWeldOptions {
            tolerance: config.vertex_tolerance_mm.min(1e-9_f64.max(config.vertex_tolerance_mm * 0.01)),
        },
    )?;

    if !config.quiet {
        eprintln!(
            "Stage 1.1 ({:.3}s): Read {} vertices, {} faces from STL",
            t.elapsed().as_secs_f64(),
            mesh.vertices.len(),
            mesh.faces.len()
        );
    }

    if !config.stage.at_least(1, 2) {
        return Ok(mesh);
    }

    // Stage 1.2: Validate and populate topology
    let t = Instant::now();
    mesh.validate_and_populate_topology()?;

    if !config.quiet {
        let s = &mesh.stats;
        eprintln!("Stage 1.2 ({:.3}s): Mesh validation passed", t.elapsed().as_secs_f64());
        eprintln!(
            "  {} faces, {} vertices, {} shells, {} solids",
            s.mesh_faces, s.mesh_vertices, s.connected_shells, s.solids
        );
        if s.mesh_edges_open > 0 {
            eprintln!("  {} open edges", s.mesh_edges_open);
        }
    }

    if config.stage.at_least(1, 3) {
        // Stage 1.3: Fuse coplanar triangle pairs into quads
        let t = Instant::now();
        fuse_coplanar_triangles(&mut mesh, config.vertex_tolerance_mm);

        if !config.quiet {
            eprintln!(
                "Stage 1.3 ({:.3}s): Fused {} triangle pairs into quads ({} faces remaining)",
                t.elapsed().as_secs_f64(),
                mesh.stats.mesh_quads, mesh.stats.mesh_faces
            );
        }
    }

    // Compare against STEP file if --compare was specified.
    if config.compare_shape.is_some() {
        compare_mesh_to_step(&mesh, config)?;
    }

    Ok(mesh)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_triangle_face(a: usize, b: usize, c: usize) -> MeshFace {
        MeshFace {
            vertex_count: 3,
            vertex_indices: [a, b, c, 0],
            neighbors: [-1, -1, -1, -1],
            normal: None,
            planar_hypothesis: UNDEDUCED_PLANAR_HYPOTHESIS,
            cylindrical_hypothesis: UNDEDUCED_CYLINDRICAL_HYPOTHESIS,
            spherical_hypothesis: UNDEDUCED_SPHERICAL_HYPOTHESIS,
            conical_hypothesis: UNDEDUCED_CONICAL_HYPOTHESIS,
            toroidal_hypothesis: UNDEDUCED_TOROIDAL_HYPOTHESIS,
        }
    }

    #[test]
    fn reads_cube_mesh_and_welds_vertices() {
        let stl_path = format!("{}/tests/manual/cube.stl", env!("CARGO_MANIFEST_DIR"));
        let mesh = read_connected_mesh_from_stl(&stl_path, VertexWeldOptions { tolerance: 1.0e-9 })
            .expect("cube mesh should load");

        assert_eq!(mesh.faces.len(), 12);
        assert_eq!(mesh.vertices.len(), 8);

        for face in &mesh.faces {
            assert_eq!(face.vertex_count, 3);
            assert_eq!(face.neighbors, [-1, -1, -1, -1]);
            assert!(face.normal.is_none());
            assert_eq!(face.planar_hypothesis, UNDEDUCED_PLANAR_HYPOTHESIS);
            assert_eq!(face.cylindrical_hypothesis, UNDEDUCED_CYLINDRICAL_HYPOTHESIS);
            assert_eq!(face.spherical_hypothesis, UNDEDUCED_SPHERICAL_HYPOTHESIS);
        }
    }

    #[test]
    fn validates_cube_mesh_topology() {
        let stl_path = format!("{}/tests/manual/cube.stl", env!("CARGO_MANIFEST_DIR"));
        let mut mesh =
            read_connected_mesh_from_stl(&stl_path, VertexWeldOptions { tolerance: 1.0e-9 })
                .expect("cube mesh should load");

        mesh.validate_and_populate_topology()
            .expect("cube should be a valid closed manifold");

        assert_eq!(mesh.stats.mesh_faces, 12);
        assert_eq!(mesh.stats.mesh_vertices, 8);
        assert_eq!(mesh.stats.mesh_edges_open, 0);
        assert_eq!(mesh.stats.mesh_edges_non_manifold, 0);
        assert_eq!(mesh.stats.mesh_edges_inconsistent_orientation, 0);
        assert_eq!(mesh.stats.mesh_faces_degenerate, 0);
        assert_eq!(mesh.stats.connected_shells, 1);
        assert_eq!(mesh.stats.solids, 1);
        assert_eq!(mesh.stats.voids_within_solids, 0);

        for face in &mesh.faces {
            assert!(face.normal.is_some());
            for edge_idx in 0..(face.vertex_count as usize) {
                assert!(face.neighbors[edge_idx] >= 0);
            }
        }
    }

    #[test]
    fn detects_inverted_normals() {
        // Tetrahedron vertices: 0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(0,0,1)
        // Correct outward winding: (0,2,1), (0,1,3), (0,3,2), (1,2,3)
        // Inverted (swap two vertices per face to flip normals inward):
        let mut mesh = ConnectedMesh {
            vertices: vec![
                MeshVertex::from_xyz(0.0, 0.0, 0.0),
                MeshVertex::from_xyz(1.0, 0.0, 0.0),
                MeshVertex::from_xyz(0.0, 1.0, 0.0),
                MeshVertex::from_xyz(0.0, 0.0, 1.0),
            ],
            faces: vec![
                make_triangle_face(0, 1, 2), // inverted bottom
                make_triangle_face(0, 3, 1), // inverted front
                make_triangle_face(1, 3, 2), // inverted hypotenuse
                make_triangle_face(0, 2, 3), // inverted left
            ],
            ..ConnectedMesh::default()
        };

        let err = mesh
            .validate_and_populate_topology()
            .expect_err("inverted tetrahedron should fail validation");

        assert!(matches!(
            err,
            MeshValidationError::InvertedNormals { shell_index: 0, .. }
        ));
    }

    #[test]
    fn detects_non_manifold_edge() {
        let mut mesh = ConnectedMesh {
            vertices: vec![
                MeshVertex::from_xyz(0.0, 0.0, 0.0),
                MeshVertex::from_xyz(1.0, 0.0, 0.0),
                MeshVertex::from_xyz(0.0, 1.0, 0.0),
                MeshVertex::from_xyz(0.0, 0.0, 1.0),
                MeshVertex::from_xyz(0.0, -1.0, 0.0),
            ],
            faces: vec![
                make_triangle_face(0, 1, 2),
                make_triangle_face(1, 0, 3),
                make_triangle_face(0, 1, 4),
            ],
            ..ConnectedMesh::default()
        };

        let err = mesh
            .validate_and_populate_topology()
            .expect_err("edge (0,1) has 3 incident triangles");

        assert!(matches!(
            err,
            MeshValidationError::NonManifoldEdges { edge_count: 1 }
        ));
        assert_eq!(mesh.stats.mesh_edges_non_manifold, 1);
    }
}
