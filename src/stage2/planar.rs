//! Stage 2.1: Planar hypothesis deduction and comparison.

use std::collections::HashSet;


use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_PLANAR_HYPOTHESIS};

use super::{
    face_area,
    PlanarHypothesis, Stage2CompareError,
};

// ---------------------------------------------------------------------------
// Planar fitting helpers
// ---------------------------------------------------------------------------

/// Fit a plane to a set of faces using area-weighted normal averaging and
/// vertex centroid for the distance.
pub(super) fn fit_plane(
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

// ---------------------------------------------------------------------------
// Stage 2.1: Compare check
// ---------------------------------------------------------------------------

/// Validate fitted planar hypotheses against a reference STEP shape.
///
/// For each hypothesis, projects face centroids onto the fitted plane and
/// checks that those projected points are within surface_tolerance of the
/// reference STEP surface.
pub(super) fn compare_planar_hypotheses(
    hypotheses: &[PlanarHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let projected_points: Vec<[f64; 3]> = hyp.faces.iter().map(|&fi| {
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
            [
                cx - dist_to_plane * hyp.normal[0],
                cy - dist_to_plane * hyp.normal[1],
                cz - dist_to_plane * hyp.normal[2],
            ]
        }).collect();

        let (max_dist, _) = stage1::max_distance_to_shape(&projected_points, compare_shape);

        if max_dist > config.surface_tolerance_mm {
            return Err(Stage2CompareError {
                hypothesis_type: "planar",
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance: config.surface_tolerance_mm,
            });
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// SurfaceType trait implementation for generic BFS
// ---------------------------------------------------------------------------

use super::bfs::{BfsConfig, Seed, SurfaceType};

/// Fitted plane parameters.
#[derive(Clone)]
pub(super) struct PlaneFit {
    pub normal: [f64; 3],
    pub distance: f64,
}

/// Marker struct for the planar surface type.
pub(super) struct PlanarSurface;

impl SurfaceType for PlanarSurface {
    type Fit = PlaneFit;
    type Hypothesis = PlanarHypothesis;

    fn name() -> &'static str { "planar" }

    fn hypothesis_field(face: &MeshFace) -> i32 { face.planar_hypothesis }

    fn set_hypothesis_field(face: &mut MeshFace, value: i32) {
        face.planar_hypothesis = value;
    }

    fn undeduced_value() -> i32 { UNDEDUCED_PLANAR_HYPOTHESIS }

    fn generate_seeds(mesh: &ConnectedMesh, fi: usize, _config: &BfsConfig) -> Vec<Seed> {
        // Single-face seed, like C++ PlanarTraits
        if mesh.faces[fi].normal.is_none() {
            return Vec::new();
        }
        vec![Seed { faces: vec![fi] }]
    }

    fn fit_surface(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        vertex_list: &[usize],
        _verbosity: u8,
    ) -> Option<Self::Fit> {
        let vertex_set: HashSet<usize> = vertex_list.iter().copied().collect();
        let (normal, distance) = fit_plane(face_list, &vertex_set, &mesh.faces, &mesh.vertices);
        Some(PlaneFit { normal, distance })
    }

    fn vertex_distance(fit: &Self::Fit, pos: &[f64; 3]) -> f64 {
        (fit.normal[0] * pos[0] + fit.normal[1] * pos[1] + fit.normal[2] * pos[2] - fit.distance)
            .abs()
    }

    fn acceptance_test(
        _mesh: &ConnectedMesh,
        face_list: &[usize],
        _fit: &Self::Fit,
        _config: &BfsConfig,
    ) -> bool {
        // Require at least 2 faces for BFS-generated planar hypotheses.
        // Single-face fallbacks are synthesized in stage 2.6.
        face_list.len() >= 2
    }

    fn make_hypothesis(
        mesh: &ConnectedMesh,
        faces: Vec<usize>,
        vertices: Vec<usize>,
        fit: &Self::Fit,
    ) -> PlanarHypothesis {
        let mut error_max = f64::NEG_INFINITY;
        let mut error_min = f64::INFINITY;
        let mut error_abs_sum = 0.0_f64;
        for &vi in &vertices {
            let v = &mesh.vertices[vi];
            let d = fit.normal[0] * v.x + fit.normal[1] * v.y + fit.normal[2] * v.z
                - fit.distance;
            error_max = error_max.max(d);
            error_min = error_min.min(d);
            error_abs_sum += d.abs();
        }
        PlanarHypothesis {
            normal: fit.normal,
            distance: fit.distance,
            faces,
            vertices,
            error_max,
            error_min,
            error_abs_sum,
        }
    }
}
