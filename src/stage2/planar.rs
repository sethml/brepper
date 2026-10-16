//! Stage 2.1: Planar hypothesis deduction and comparison.

use std::collections::{HashMap, HashSet};


use crate::stage1::{ConnectedMesh, MeshFace, MeshVertex};

use super::{
    face_area,
    PlanarHypothesis,
};

// ---------------------------------------------------------------------------
// Planar fitting helpers
// ---------------------------------------------------------------------------

/// Fit a plane to a set of faces: area-weighted normal averaging, and an
/// area-weighted vertex centroid for the distance.
///
/// The area weighting on the offset matters as much as it does on the normal.
/// An unweighted vertex mean gives every vertex an equal vote, which is a
/// discretization artifact — it over-weights whatever part of the surface the
/// tessellator happened to sample densely. The continuous least-squares plane
/// minimizes ∫(n·x − d)² dA, and weighting each vertex by the area it stands
/// for is that integral's discretization.
///
/// It also makes the fit robust to a small-area leak, which is what exposed
/// this: the end cap of `stepped_block_toroidal_fillets` is a handful of large
/// flat quads ringed by the many tiny quads of its perimeter fillet's first
/// row, 9.5e-4 mm off the plane. Those rows carry ~1% of the area but about
/// half the vertices, so the unweighted mean placed the plane at the *midpoint*
/// of the two levels — every one of its 95 faces 4.7e-4 mm out, with no
/// outlier group for anything downstream to spot.
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

    // Each vertex carries the area of its incident faces within the support,
    // shared equally among that face's corners.
    let mut weight: HashMap<usize, f64> = HashMap::with_capacity(vertex_set.len());
    for &fi in face_indices {
        let face = &faces[fi];
        let vc = face.vertex_count as usize;
        let share = face_area(face, vertices) / vc as f64;
        for &vi in &face.vertex_indices[..vc] {
            if vertex_set.contains(&vi) {
                *weight.entry(vi).or_insert(0.0) += share;
            }
        }
    }
    // Sum in sorted vertex order, never `HashSet` order: the set's iteration
    // order is randomized per process, so summing in it made the fitted
    // distance differ in its last bits from run to run — and with it every
    // decision downstream that compares a residual against a tolerance. Found
    // 2026-10-16 when stage 2 output, fingerprinted at full precision, differed
    // between two runs of the same binary on 5 of 7 models.
    let mut ordered: Vec<usize> = vertex_set.iter().copied().collect();
    ordered.sort_unstable();
    let mut sum_d = 0.0_f64;
    let mut sum_w = 0.0_f64;
    for &vi in &ordered {
        let v = &vertices[vi];
        let w = weight.get(&vi).copied().unwrap_or(0.0);
        sum_d += w * (normal[0] * v.x + normal[1] * v.y + normal[2] * v.z);
        sum_w += w;
    }
    let distance = if sum_w > 0.0 {
        sum_d / sum_w
    } else {
        // Degenerate support (zero total area) — fall back to the plain mean.
        ordered
            .iter()
            .map(|&vi| {
                let v = &vertices[vi];
                normal[0] * v.x + normal[1] * v.y + normal[2] * v.z
            })
            .sum::<f64>()
            / vertex_set.len().max(1) as f64
    };

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
        // The area-weighted normal of a closed or symmetric region (a full
        // cylindrical hole wall) cancels to ~zero and normalizes to NaN. A
        // non-finite plane would then VACUOUSLY validate — every `dist > vtol`
        // test is false for NaN — and be committed over the cylinder the
        // region actually is. Reject it so the ladder falls through.
        if !normal.iter().all(|c| c.is_finite()) || !distance.is_finite() {
            return None;
        }
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
