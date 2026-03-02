//! Stage 2: Surface Fitting
//!
//! 2.1: Deduce planar hypotheses
//! 2.2: Deduce cylindrical hypotheses
//! 2.3: Deduce spherical hypotheses
//! 2.4: Deduce ruled surface hypotheses
//! 2.5: Deduce NURBS hypotheses
//! 2.6: Select surfaces for reconstruction

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_PLANAR_HYPOTHESIS};
use opencascade_sys::gp;
use std::collections::{HashSet, VecDeque};
use std::error::Error;
use std::fmt::{Display, Formatter};

// ---------------------------------------------------------------------------
// Hypothesis data structures
// ---------------------------------------------------------------------------

/// A planar surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct PlanarHypothesis {
    /// Unit normal vector pointing outward from the shell/solid.
    pub normal: [f64; 3],
    /// Signed distance from origin to plane along the normal.
    pub distance: f64,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this plane.
    pub vertices: Vec<usize>,
    /// Maximum (positive) distance from any vertex to the plane.
    pub error_max: f64,
    /// Minimum (most negative) distance from any vertex to the plane.
    pub error_min: f64,
    /// Sum of absolute vertex-to-plane distances.
    pub error_abs_sum: f64,
}

/// A cylindrical surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct CylindricalHypothesis {
    /// A point on the cylinder axis.
    pub axis_origin: [f64; 3],
    /// Unit direction vector along the cylinder axis.
    pub axis_direction: [f64; 3],
    /// Radius of the cylinder (always positive).
    pub radius: f64,
    /// Whether the surface normal points away from the axis (convex=true) or toward it (concave=false).
    pub convex: bool,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this cylinder.
    pub vertices: Vec<usize>,
    /// Maximum absolute distance from any vertex to the cylinder surface.
    pub error_max: f64,
    /// Sum of absolute vertex-to-surface distances.
    pub error_abs_sum: f64,
}

/// A spherical surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct SphericalHypothesis {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere (always positive).
    pub radius: f64,
    /// Whether the surface normal points away from center (convex) or toward it (concave).
    pub convex: bool,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this sphere.
    pub vertices: Vec<usize>,
    /// Maximum absolute distance from any vertex to the sphere surface.
    pub error_max: f64,
    /// Sum of absolute vertex-to-surface distances.
    pub error_abs_sum: f64,
}

// TODO: Stage 2.4 - Ruled surface hypothesis
// TODO: Stage 2.5 - NURBS hypothesis

/// Identifies a selected surface hypothesis by type and index.
#[derive(Debug, Clone, Copy)]
pub enum SelectedSurface {
    Planar(usize),
    Cylindrical(usize),
    Spherical(usize),
    // TODO: RuledSurface(usize),
    // TODO: Nurbs(usize),
}

// ---------------------------------------------------------------------------
// Stage 2 output
// ---------------------------------------------------------------------------

/// The output of Stage 2: the mesh with all hypotheses populated and surfaces selected.
#[derive(Debug)]
pub struct Stage2Output {
    /// The mesh from stage 1 with per-face hypothesis indices populated.
    pub mesh: ConnectedMesh,
    /// All planar hypotheses deduced in stage 2.1.
    pub planar_hypotheses: Vec<PlanarHypothesis>,
    /// All cylindrical hypotheses deduced in stage 2.2.
    pub cylindrical_hypotheses: Vec<CylindricalHypothesis>,
    /// All spherical hypotheses deduced in stage 2.3.
    pub spherical_hypotheses: Vec<SphericalHypothesis>,
    /// Surfaces selected in stage 2.6 for reconstruction. Each face should be
    /// covered by exactly one selected surface.
    pub selected_surfaces: Vec<SelectedSurface>,
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Stage2Error {
    NotImplemented(String),
    Compare(Stage2CompareError),
}

impl Display for Stage2Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage2Error::NotImplemented(stage) => write!(f, "stage {stage} not yet implemented"),
            Stage2Error::Compare(e) => write!(f, "stage 2 compare: {e}"),
        }
    }
}

impl Error for Stage2Error {}

impl From<Stage2CompareError> for Stage2Error {
    fn from(e: Stage2CompareError) -> Self {
        Stage2Error::Compare(e)
    }
}

#[derive(Debug)]
pub struct Stage2CompareError {
    pub hypothesis_index: usize,
    pub max_distance: f64,
    pub tolerance: f64,
}

impl Display for Stage2CompareError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "planar hypothesis {} projected centroid is {:.6e} mm from nearest STEP surface (tolerance: {:.6e} mm)",
            self.hypothesis_index, self.max_distance, self.tolerance
        )
    }
}

impl Error for Stage2CompareError {}

// ---------------------------------------------------------------------------
// Stage 2.1: Planar hypothesis deduction
// ---------------------------------------------------------------------------

/// Signed distance from a vertex to a plane defined by (normal, distance).
fn vertex_to_plane_distance(v: &MeshVertex, normal: &[f64; 3], distance: f64) -> f64 {
    normal[0] * v.x + normal[1] * v.y + normal[2] * v.z - distance
}

/// Compute the area of a triangular mesh face.
fn face_area(face: &MeshFace, vertices: &[MeshVertex]) -> f64 {
    let v0 = &vertices[face.vertex_indices[0]];
    let v1 = &vertices[face.vertex_indices[1]];
    let v2 = &vertices[face.vertex_indices[2]];
    let ax = v1.x - v0.x;
    let ay = v1.y - v0.y;
    let az = v1.z - v0.z;
    let bx = v2.x - v0.x;
    let by = v2.y - v0.y;
    let bz = v2.z - v0.z;
    let cx = ay * bz - az * by;
    let cy = az * bx - ax * bz;
    let cz = ax * by - ay * bx;
    0.5 * (cx * cx + cy * cy + cz * cz).sqrt()
}

/// Fit a plane to a set of faces using area-weighted normal averaging and
/// vertex centroid for the distance.
fn fit_plane(
    face_indices: &[usize],
    vertex_set: &HashSet<usize>,
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> ([f64; 3], f64) {
    let mut sum_nx = 0.0_f64;
    let mut sum_ny = 0.0_f64;
    let mut sum_nz = 0.0_f64;

    for &fi in face_indices {
        let face = &faces[fi];
        let n = face.normal.unwrap();
        let area = face_area(face, vertices);
        sum_nx += n[0] * area;
        sum_ny += n[1] * area;
        sum_nz += n[2] * area;
    }

    let len = (sum_nx * sum_nx + sum_ny * sum_ny + sum_nz * sum_nz).sqrt();
    let normal = [sum_nx / len, sum_ny / len, sum_nz / len];

    let mut sum_d = 0.0_f64;
    for &vi in vertex_set {
        let v = &vertices[vi];
        sum_d += normal[0] * v.x + normal[1] * v.y + normal[2] * v.z;
    }
    let distance = sum_d / vertex_set.len() as f64;

    (normal, distance)
}

/// Check that all vertices in a set are within tolerance of a plane.
fn all_vertices_within_tolerance(
    vertex_set: &HashSet<usize>,
    normal: &[f64; 3],
    distance: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        if vertex_to_plane_distance(&vertices[vi], normal, distance).abs() > tolerance {
            return false;
        }
    }
    true
}

/// Deduce planar hypotheses from the mesh using BFS region growing.
///
/// For each unassigned face, creates a new planar hypothesis seeded from that
/// face's plane, then grows it via BFS to neighboring faces that are coplanar
/// (normal alignment and vertex-to-plane distance within tolerance).
fn deduce_planar_hypotheses(
    mesh: &mut ConnectedMesh,
    vertex_tol: f64,
) -> Vec<PlanarHypothesis> {
    let num_faces = mesh.faces.len();
    let mut hypotheses: Vec<PlanarHypothesis> = Vec::new();

    for face in &mut mesh.faces {
        face.planar_hypothesis = UNDEDUCED_PLANAR_HYPOTHESIS;
    }

    for fi in 0..num_faces {
        if mesh.faces[fi].planar_hypothesis != UNDEDUCED_PLANAR_HYPOTHESIS {
            continue;
        }

        let hi = hypotheses.len() as i32;
        let seed_normal = mesh.faces[fi].normal.unwrap();

        // Compute initial plane from seed face: normal from face, distance
        // as average projection of seed face vertices.
        let vc = mesh.faces[fi].vertex_count as usize;
        let mut sum_d = 0.0_f64;
        let mut vertex_set = HashSet::new();
        for vi_idx in 0..vc {
            let vi = mesh.faces[fi].vertex_indices[vi_idx];
            vertex_set.insert(vi);
            let v = &mesh.vertices[vi];
            sum_d += seed_normal[0] * v.x + seed_normal[1] * v.y + seed_normal[2] * v.z;
        }
        let mut current_normal = seed_normal;
        let mut current_distance = sum_d / vc as f64;
        let mut face_list: Vec<usize> = vec![fi];

        mesh.faces[fi].planar_hypothesis = hi;

        // BFS expansion
        let mut queue = VecDeque::new();
        queue.push_back(fi);

        while let Some(current_fi) = queue.pop_front() {
            let vertex_count = mesh.faces[current_fi].vertex_count as usize;
            let neighbors = mesh.faces[current_fi].neighbors;

            for edge_idx in 0..vertex_count {
                let ni = neighbors[edge_idx];
                if ni < 0 {
                    continue;
                }
                let ni = ni as usize;

                if mesh.faces[ni].planar_hypothesis != UNDEDUCED_PLANAR_HYPOTHESIS {
                    continue;
                }
                // Vertex distance check
                let nvc = mesh.faces[ni].vertex_count as usize;
                let nvi = mesh.faces[ni].vertex_indices;
                let mut all_ok = true;
                let mut any_far = false;
                for vi_idx in 0..nvc {
                    let d = vertex_to_plane_distance(
                        &mesh.vertices[nvi[vi_idx]],
                        &current_normal,
                        current_distance,
                    );
                    let abs_d = d.abs();
                    if abs_d > vertex_tol {
                        all_ok = false;
                        if abs_d > 2.0 * vertex_tol {
                            any_far = true;
                            break;
                        }
                    }
                }

                // Skip if any vertex is too far for re-fitting to help
                if any_far {
                    continue;
                }

                if !all_ok {
                    // Try re-fitting with current + new face's vertices
                    let mut trial_vertices = vertex_set.clone();
                    for vi_idx in 0..nvc {
                        trial_vertices.insert(nvi[vi_idx]);
                    }
                    let mut trial_faces = face_list.clone();
                    trial_faces.push(ni);

                    let (new_normal, new_distance) = fit_plane(
                        &trial_faces,
                        &trial_vertices,
                        &mesh.faces,
                        &mesh.vertices,
                    );

                    // Check ALL vertices (existing + new) against the re-fitted plane
                    if !all_vertices_within_tolerance(
                        &trial_vertices,
                        &new_normal,
                        new_distance,
                        vertex_tol,
                        &mesh.vertices,
                    ) {
                        continue;
                    }

                    // Accept re-fit
                    current_normal = new_normal;
                    current_distance = new_distance;
                }

                // Accept this face into the hypothesis
                mesh.faces[ni].planar_hypothesis = hi;
                face_list.push(ni);
                for vi_idx in 0..nvc {
                    vertex_set.insert(nvi[vi_idx]);
                }
                queue.push_back(ni);
            }
        }

        // Final re-fit using all collected faces and vertices
        let (final_normal, final_distance) =
            fit_plane(&face_list, &vertex_set, &mesh.faces, &mesh.vertices);

        // Compute error metrics
        let mut error_max = f64::NEG_INFINITY;
        let mut error_min = f64::INFINITY;
        let mut error_abs_sum = 0.0_f64;
        for &vi in &vertex_set {
            let d = vertex_to_plane_distance(&mesh.vertices[vi], &final_normal, final_distance);
            error_max = error_max.max(d);
            error_min = error_min.min(d);
            error_abs_sum += d.abs();
        }

        hypotheses.push(PlanarHypothesis {
            normal: final_normal,
            distance: final_distance,
            faces: face_list,
            vertices: vertex_set.into_iter().collect(),
            error_max,
            error_min,
            error_abs_sum,
        });
    }

    hypotheses
}

// ---------------------------------------------------------------------------
// Stage 2.1: Compare check
// ---------------------------------------------------------------------------

/// Validate fitted planar hypotheses against a reference STEP shape.
///
/// For each hypothesis, projects face centroids onto the fitted plane and
/// checks that those projected points are within surface_tolerance of the
/// reference STEP surface.
fn compare_planar_hypotheses(
    hypotheses: &[PlanarHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let face = &mesh.faces[fi];
            let n = face.vertex_count as usize;
            let mut cx = 0.0_f64;
            let mut cy = 0.0_f64;
            let mut cz = 0.0_f64;
            for vi_idx in 0..n {
                let v = &mesh.vertices[face.vertex_indices[vi_idx]];
                cx += v.x;
                cy += v.y;
                cz += v.z;
            }
            let inv_n = 1.0 / n as f64;
            cx *= inv_n;
            cy *= inv_n;
            cz *= inv_n;

            // Project centroid onto the fitted plane
            let dist_to_plane = hyp.normal[0] * cx + hyp.normal[1] * cy + hyp.normal[2] * cz
                - hyp.distance;
            let px = cx - dist_to_plane * hyp.normal[0];
            let py = cy - dist_to_plane * hyp.normal[1];
            let pz = cz - dist_to_plane * hyp.normal[2];

            let pt = gp::Pnt::new_real3(px, py, pz);
            let d = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(d);
        }

        if max_dist > config.surface_tolerance_mm {
            return Err(Stage2CompareError {
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance: config.surface_tolerance_mm,
            });
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Stage 2 entry point
// ---------------------------------------------------------------------------

/// Run stage 2: fit surface hypotheses to mesh faces and select surfaces.
pub fn stage2(config: &Config, mut mesh: ConnectedMesh) -> Result<Stage2Output, Stage2Error> {
    // Stage 2.1: Deduce planar hypotheses
    let planar_hypotheses = deduce_planar_hypotheses(&mut mesh, config.vertex_tolerance_mm);

    if !config.quiet {
        let multi_face_count = planar_hypotheses.iter().filter(|h| h.faces.len() > 1).count();
        let single_face_count = planar_hypotheses.len() - multi_face_count;
        let covered_faces: usize = planar_hypotheses.iter().map(|h| h.faces.len()).sum();
        eprintln!(
            "Stage 2.1: Deduced {} planar hypotheses ({} multi-face, {} single-face) covering {} of {} mesh faces",
            planar_hypotheses.len(),
            multi_face_count,
            single_face_count,
            covered_faces,
            mesh.faces.len(),
        );
        if config.verbose {
            for (i, h) in planar_hypotheses.iter().enumerate() {
                if h.faces.len() > 1 {
                    eprintln!(
                        "  Plane {}: {} faces, {} vertices, normal=[{:.4}, {:.4}, {:.4}], d={:.4}, err_max={:.2e}, err_min={:.2e}",
                        i, h.faces.len(), h.vertices.len(),
                        h.normal[0], h.normal[1], h.normal[2], h.distance,
                        h.error_max, h.error_min,
                    );
                }
            }
        }
    }

    // Compare against STEP file if --compare was specified
    if config.compare_shape.is_some() {
        compare_planar_hypotheses(&planar_hypotheses, &mesh, config)?;
        if !config.quiet {
            eprintln!("  Compare: all planar hypothesis centroids within tolerance");
        }
    }

    if !config.stage.at_least(2, 2) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.2: Deduce cylindrical hypotheses
    // TODO: implement cylinder fitting
    if !config.quiet {
        eprintln!("Stage 2.2: Deduce cylindrical hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 3) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.3: Deduce spherical hypotheses
    // TODO: implement sphere fitting
    if !config.quiet {
        eprintln!("Stage 2.3: Deduce spherical hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 4) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.4: Deduce ruled surface hypotheses
    // TODO: optional - detect extruded curve surfaces
    if !config.quiet {
        eprintln!("Stage 2.4: Deduce ruled surface hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 5) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.5: Deduce NURBS hypotheses
    // TODO: fit NURBS to remaining ungrouped faces
    if !config.quiet {
        eprintln!("Stage 2.5: Deduce NURBS hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 6) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.6: Select surfaces for reconstruction
    // TODO: greedy selection of best-fitting hypotheses covering all faces
    if !config.quiet {
        eprintln!("Stage 2.6: Select surfaces (not yet implemented)");
    }

    Ok(Stage2Output {
        mesh,
        planar_hypotheses,
        cylindrical_hypotheses: Vec::new(),
        spherical_hypotheses: Vec::new(),
        selected_surfaces: Vec::new(),
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stage1::{self, MeshFace, MeshVertex, VertexWeldOptions, NO_HYPOTHESIS};

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

    /// Build a simple mesh from vertices and triangle indices, then validate.
    fn build_mesh(verts: Vec<MeshVertex>, tris: Vec<[usize; 3]>) -> ConnectedMesh {
        let faces = tris.iter().map(|t| make_triangle_face(t[0], t[1], t[2])).collect();
        let mut mesh = ConnectedMesh {
            vertices: verts,
            faces,
            stats: Default::default(),
        };
        mesh.validate_and_populate_topology().expect("mesh should be valid");
        mesh
    }

    #[test]
    fn cube_gets_six_planar_hypotheses() {
        let stl_path = format!("{}/tests/manual/cube.stl", env!("CARGO_MANIFEST_DIR"));
        let mut mesh =
            stage1::read_connected_mesh_from_stl(&stl_path, VertexWeldOptions { tolerance: 1e-9 })
                .expect("should load");
        mesh.validate_and_populate_topology().expect("should validate");

        let hypotheses = deduce_planar_hypotheses(&mut mesh, 1e-5);

        // A cube has 6 faces, each with 2 triangles
        assert_eq!(mesh.faces.len(), 12);

        // Should produce exactly 6 planar hypotheses (one per cube face)
        assert_eq!(
            hypotheses.len(),
            6,
            "expected 6 planar hypotheses for cube, got {}: {:?}",
            hypotheses.len(),
            hypotheses.iter().map(|h| h.faces.len()).collect::<Vec<_>>()
        );

        // Each hypothesis should have exactly 2 faces
        for (i, h) in hypotheses.iter().enumerate() {
            assert_eq!(
                h.faces.len(),
                2,
                "hypothesis {i} should have 2 faces, got {}",
                h.faces.len()
            );
        }

        // Total covered faces should be 12
        let total: usize = hypotheses.iter().map(|h| h.faces.len()).sum();
        assert_eq!(total, 12);

        // All errors should be essentially zero for a perfect cube
        for h in &hypotheses {
            assert!(
                h.error_max.abs() < 1e-12,
                "error_max should be ~0, got {}",
                h.error_max
            );
            assert!(
                h.error_min.abs() < 1e-12,
                "error_min should be ~0, got {}",
                h.error_min
            );
        }
    }

    #[test]
    fn two_perpendicular_triangles_get_separate_hypotheses() {
        // Two triangles sharing an edge but on perpendicular planes
        let verts = vec![
            MeshVertex::from_xyz(0.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0), // triangle 0: on XY plane
            MeshVertex::from_xyz(0.0, 0.0, 1.0), // triangle 1: on XZ plane
        ];
        let mesh = build_mesh(verts, vec![[0, 1, 2], [0, 3, 1]]);
        let mut mesh = mesh;
        let hypotheses = deduce_planar_hypotheses(&mut mesh, 1e-5);

        assert_eq!(
            hypotheses.len(),
            2,
            "perpendicular triangles should get 2 hypotheses"
        );
    }

    #[test]
    fn coplanar_triangles_get_one_hypothesis() {
        // Two coplanar triangles sharing an edge on the XY plane
        let verts = vec![
            MeshVertex::from_xyz(0.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 1.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0),
        ];
        let mesh = build_mesh(verts, vec![[0, 1, 2], [0, 2, 3]]);
        let mut mesh = mesh;
        let hypotheses = deduce_planar_hypotheses(&mut mesh, 1e-5);

        assert_eq!(
            hypotheses.len(),
            1,
            "two coplanar triangles should get 1 hypothesis, got {}",
            hypotheses.len()
        );
        assert_eq!(hypotheses[0].faces.len(), 2);
        assert_eq!(hypotheses[0].vertices.len(), 4);
    }
}