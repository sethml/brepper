//! Stage 2.1: Planar hypothesis deduction and comparison.

use std::collections::HashSet;


use crate::stage1::{ConnectedMesh, MeshFace, MeshVertex};

use super::{
    face_area,
    PlanarHypothesis,
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

// ---------------------------------------------------------------------------
// SurfaceType trait implementation for generic BFS
// ---------------------------------------------------------------------------

use super::surface::SurfaceType;

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
