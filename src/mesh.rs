use opencascade_sys::{message, rw_stl};
use std::collections::{HashMap, HashSet, VecDeque};
use std::error::Error;
use std::fmt::{Display, Formatter};

pub const UNDEDUCED_PLANAR_HYPOTHESIS: i32 = -2;
pub const NO_HYPOTHESIS: i32 = -1;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeshVertex {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl MeshVertex {
    fn from_xyz(x: f64, y: f64, z: f64) -> Self {
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
}

#[derive(Debug, Clone, Default)]
pub struct PlanarHypothesis;

#[derive(Debug, Clone, Default)]
pub struct CylindricalHypothesis;

#[derive(Debug, Clone, Default)]
pub struct SphericalHypothesis;

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
}

#[derive(Debug, Clone, Default)]
pub struct ConnectedMesh {
    pub vertices: Vec<MeshVertex>,
    pub faces: Vec<MeshFace>,
    pub planar_hypotheses: Vec<PlanarHypothesis>,
    pub cylindrical_hypotheses: Vec<CylindricalHypothesis>,
    pub spherical_hypotheses: Vec<SphericalHypothesis>,
    pub stats: MeshValidationStats,
}

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
    EmptyMesh,
}

impl Display for MeshReadError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            MeshReadError::InvalidTolerance(tol) => {
                write!(f, "invalid weld tolerance: {tol}, expected > 0")
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
        }
    }
}

impl Error for MeshValidationError {}

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
        for slot in 0..3 {
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

            face_vertex_indices[slot] = welded_idx;
        }

        faces.push(MeshFace {
            vertex_count: 3,
            vertex_indices: face_vertex_indices,
            neighbors: [-1, -1, -1, -1],
            normal: None,
            planar_hypothesis: UNDEDUCED_PLANAR_HYPOTHESIS,
            cylindrical_hypothesis: NO_HYPOTHESIS,
            spherical_hypothesis: NO_HYPOTHESIS,
        });
    }

    Ok(ConnectedMesh {
        vertices: welder.vertices,
        faces,
        planar_hypotheses: Vec::new(),
        cylindrical_hypotheses: Vec::new(),
        spherical_hypotheses: Vec::new(),
        stats: MeshValidationStats::default(),
    })
}

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
            if vertex_count < 3 || vertex_count > 4 {
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
        self.stats.solids = component_closed.into_iter().filter(|closed| *closed).count();

        // Self-intersection checks are intentionally deferred to a later stage where we have
        // richer geometric predicates and can avoid expensive O(N^2) triangle checks here.

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
    if vertex_count < 3 || vertex_count > 4 {
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

#[cfg(test)]
mod tests {
    use super::{
        read_connected_mesh_from_stl, ConnectedMesh, MeshFace, MeshValidationError, MeshVertex,
        VertexWeldOptions, NO_HYPOTHESIS, UNDEDUCED_PLANAR_HYPOTHESIS,
    };

    fn make_triangle_face(a: usize, b: usize, c: usize) -> MeshFace {
        MeshFace {
            vertex_count: 3,
            vertex_indices: [a, b, c, 0],
            neighbors: [-1, -1, -1, -1],
            normal: None,
            planar_hypothesis: UNDEDUCED_PLANAR_HYPOTHESIS,
            cylindrical_hypothesis: NO_HYPOTHESIS,
            spherical_hypothesis: NO_HYPOTHESIS,
        }
    }

    #[test]
    fn reads_cube_mesh_and_welds_vertices() {
        let stl_path = format!("{}/tests/manual/cube.stl", env!("CARGO_MANIFEST_DIR"));
        let mesh = read_connected_mesh_from_stl(&stl_path, VertexWeldOptions { tolerance: 1.0e-9 })
            .expect("cube mesh should load");

        assert_eq!(mesh.faces.len(), 12);
        assert_eq!(mesh.vertices.len(), 8);
        assert!(mesh.planar_hypotheses.is_empty());
        assert!(mesh.cylindrical_hypotheses.is_empty());
        assert!(mesh.spherical_hypotheses.is_empty());

        for face in &mesh.faces {
            assert_eq!(face.vertex_count, 3);
            assert_eq!(face.neighbors, [-1, -1, -1, -1]);
            assert!(face.normal.is_none());
            assert_eq!(face.planar_hypothesis, UNDEDUCED_PLANAR_HYPOTHESIS);
            assert_eq!(face.cylindrical_hypothesis, NO_HYPOTHESIS);
            assert_eq!(face.spherical_hypothesis, NO_HYPOTHESIS);
        }
    }

    #[test]
    fn validates_cube_mesh_topology() {
        let stl_path = format!("{}/tests/manual/cube.stl", env!("CARGO_MANIFEST_DIR"));
        let mut mesh = read_connected_mesh_from_stl(&stl_path, VertexWeldOptions { tolerance: 1.0e-9 })
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