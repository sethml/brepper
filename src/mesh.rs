use opencascade_sys::{message, rw_stl};
use std::collections::HashMap;
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

#[cfg(test)]
mod tests {
    use super::{read_connected_mesh_from_stl, VertexWeldOptions, NO_HYPOTHESIS, UNDEDUCED_PLANAR_HYPOTHESIS};

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
}