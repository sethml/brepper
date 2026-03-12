//! Stage 2.4: Conical hypothesis deduction and comparison.

use std::collections::{HashSet, VecDeque};

use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{Const, Dyn, OMatrix, OVector, Owned};
use opencascade_sys::gp;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_CONICAL_HYPOTHESIS, NO_HYPOTHESIS};
use crate::viz::{self, VizAction, VizSender};

use super::{
    bounding_box_diagonal, build_normal_covariance, cross3, dot3, face_area, face_centroid,
    normalize3, perpendicular_basis, smallest_eigenvector_3x3, viz_custom, viz_face_centroid,
    viz_face_normal, ConicalHypothesis, Stage2CompareError,
    MIN_CONE_FACES, REFIT_SKIP_MULTIPLIER,
};

/// Minimum half-angle (degrees) — below this the cone is effectively a cylinder.
const MIN_HALF_ANGLE_DEG: f64 = 2.0;
/// Maximum half-angle (degrees) — above this the cone is effectively a plane.
const MAX_HALF_ANGLE_DEG: f64 = 85.0;
/// Maximum apex distance as a multiple of the mesh bounding-box diagonal.
const MAX_APEX_DISTANCE_FACTOR: f64 = 10.0;

// ---------------------------------------------------------------------------
// Conical fitting helpers
// ---------------------------------------------------------------------------

/// Compute the signed distance from a vertex to a cone surface.
/// Positive if outside the cone, negative if inside.
///
/// For a cone with apex A, axis direction â, half-angle θ:
///   h = (v - A) · â           (height above apex along axis)
///   r = |(v - A) - h·â|       (radial distance from axis)
///   d = r·cos(θ) - h·sin(θ)   (signed distance: 0 on surface)
pub(super) fn vertex_to_cone_distance(
    v: &MeshVertex,
    apex: &[f64; 3],
    axis_direction: &[f64; 3],
    half_angle: f64,
) -> f64 {
    let dx = v.x - apex[0];
    let dy = v.y - apex[1];
    let dz = v.z - apex[2];
    let h = dx * axis_direction[0] + dy * axis_direction[1] + dz * axis_direction[2];
    let rx = dx - h * axis_direction[0];
    let ry = dy - h * axis_direction[1];
    let rz = dz - h * axis_direction[2];
    let r = (rx * rx + ry * ry + rz * rz).sqrt();
    let (sin_a, cos_a) = half_angle.sin_cos();
    r * cos_a - h * sin_a
}

/// Given an axis direction and a set of vertices, project them to (h, r) coordinates
/// relative to the centroid, then fit a linear profile r = m·h + b.
/// Returns (apex, half_angle) or None if degenerate.
///
/// The cone axis goes through the centroid. We compute:
///   h_i = (v_i - centroid) · â
///   r_i = |(v_i - centroid) - h_i·â|
/// Then linear regression r = m·h + b gives:
///   half_angle = atan(|m|)
///   apex = centroid - (b/m)·â  (if m != 0)
///
/// If the linear regression is degenerate (e.g., all vertices at same height),
/// we fall back to a normal-based approach using face normals to determine the
/// half-angle and vertex positions to locate the apex.
fn fit_cone_profile(
    axis_dir: &[f64; 3],
    vertex_set: &HashSet<usize>,
    vertices: &[MeshVertex],
    face_indices: &[usize],
    faces: &[MeshFace],
) -> Option<([f64; 3], f64)> {
    let n = vertex_set.len();
    if n < 3 {
        return None;
    }

    // Compute centroid of all vertices
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    for &vi in vertex_set {
        let v = &vertices[vi];
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    let inv_n = 1.0 / n as f64;
    cx *= inv_n;
    cy *= inv_n;
    cz *= inv_n;

    // Project vertices to (h, r) coordinates
    let mut sum_h = 0.0_f64;
    let mut sum_r = 0.0_f64;
    let mut sum_hh = 0.0_f64;
    let mut sum_hr = 0.0_f64;
    let mut h_min = f64::MAX;
    let mut h_max = f64::MIN;

    for &vi in vertex_set {
        let v = &vertices[vi];
        let dx = v.x - cx;
        let dy = v.y - cy;
        let dz = v.z - cz;
        let h = dx * axis_dir[0] + dy * axis_dir[1] + dz * axis_dir[2];
        let rx = dx - h * axis_dir[0];
        let ry = dy - h * axis_dir[1];
        let rz = dz - h * axis_dir[2];
        let r = (rx * rx + ry * ry + rz * rz).sqrt();
        sum_h += h;
        sum_r += r;
        sum_hh += h * h;
        sum_hr += h * r;
        h_min = h_min.min(h);
        h_max = h_max.max(h);
    }

    let nf = n as f64;
    let denom = nf * sum_hh - sum_h * sum_h;

    // Try linear regression first (works when vertices span different heights)
    if denom.abs() > 1e-30 {
        let m = (nf * sum_hr - sum_h * sum_r) / denom;
        let b = (sum_r - m * sum_h) / nf;

        // Reject if nearly cylindrical (|m| < 0.01) or nearly flat (half-angle > 85°)
        if m.abs() >= 0.01 {
            let half_angle = m.abs().atan();
            if half_angle <= 85.0_f64.to_radians() {
                // Apex position: distance from centroid along axis = -b/m
                let apex_t = -b / m;
                let apex = [
                    cx + apex_t * axis_dir[0],
                    cy + apex_t * axis_dir[1],
                    cz + apex_t * axis_dir[2],
                ];
                return Some((apex, half_angle));
            }
        }
    }

    // Fallback: determine half-angle from face normals.
    // For a cone with half-angle θ, the normal at any surface point satisfies:
    //   |n · axis| = sin(θ)
    // We average this over all seed faces (area-weighted) to get a robust estimate.
    let mut weighted_sin = 0.0_f64;
    let mut total_area = 0.0_f64;
    for &fi in face_indices {
        if let Some(normal) = faces[fi].normal {
            let area = face_area(&faces[fi], vertices);
            let sin_theta = dot3(&normal, axis_dir).abs();
            weighted_sin += sin_theta * area;
            total_area += area;
        }
    }
    if total_area < 1e-30 {
        return None;
    }
    let avg_sin = weighted_sin / total_area;
    // sin(θ) must be in (0, 1) for a valid cone (not cylinder or plane)
    if !(0.01..=0.999).contains(&avg_sin) {
        return None;
    }
    let half_angle = avg_sin.asin();

    // Determine apex: for each vertex, the apex is at distance r/tan(θ)
    // below the vertex along the axis. We average over all vertices.
    // Use axis-projected centroid and average radius to find apex.
    let mean_h = sum_h / nf;
    let mean_r = sum_r / nf;
    // On a cone: r = (h - h_apex) * tan(θ), so h_apex = h - r/tan(θ)
    // The sign depends on which side of the apex we're on. Try both orientations.
    let tan_theta = half_angle.tan();
    if tan_theta.abs() < 1e-10 {
        return None;
    }
    let apex_h_candidate = mean_h - mean_r / tan_theta;

    // Apex in world coordinates: centroid + apex_h * axis
    let apex = [
        cx + apex_h_candidate * axis_dir[0],
        cy + apex_h_candidate * axis_dir[1],
        cz + apex_h_candidate * axis_dir[2],
    ];

    Some((apex, half_angle))
}

// ---------------------------------------------------------------------------
// LM refinement for conical fitting
// ---------------------------------------------------------------------------

/// LM problem for cone fitting.
/// Parameters: [alpha, beta, apex_x, apex_y, apex_z, half_angle]
/// alpha, beta: tilt of axis from initial direction in perpendicular directions
struct ConeLMProblem {
    points: Vec<[f64; 3]>,
    a0: [f64; 3],
    u0: [f64; 3],
    w0: [f64; 3],
    params: OVector<f64, Const<6>>,
}

impl ConeLMProblem {
    fn new(
        points: Vec<[f64; 3]>,
        initial_axis: [f64; 3],
        initial_apex: [f64; 3],
        initial_half_angle: f64,
    ) -> Self {
        let (u0, w0) = perpendicular_basis(&initial_axis);
        let params = OVector::<f64, Const<6>>::new(
            0.0, 0.0,
            initial_apex[0], initial_apex[1], initial_apex[2],
            initial_half_angle,
        );
        Self { points, a0: initial_axis, u0, w0, params }
    }

    fn axis_dir_from_params(&self, params: &OVector<f64, Const<6>>) -> [f64; 3] {
        let (alpha, beta) = (params[0], params[1]);
        let mut v = [
            self.a0[0] + alpha * self.u0[0] + beta * self.w0[0],
            self.a0[1] + alpha * self.u0[1] + beta * self.w0[1],
            self.a0[2] + alpha * self.u0[2] + beta * self.w0[2],
        ];
        normalize3(&mut v);
        v
    }

    fn compute_residuals_for(&self, params: &OVector<f64, Const<6>>) -> OVector<f64, Dyn> {
        let a = self.axis_dir_from_params(params);
        let apex = [params[2], params[3], params[4]];
        let half_angle = params[5];
        let (sin_a, cos_a) = half_angle.sin_cos();
        let n = self.points.len();
        let mut residuals = OVector::<f64, Dyn>::zeros_generic(Dyn(n), Const::<1>);
        for (i, p) in self.points.iter().enumerate() {
            let d = [p[0] - apex[0], p[1] - apex[1], p[2] - apex[2]];
            let h = dot3(&d, &a);
            let radial = [d[0] - h * a[0], d[1] - h * a[1], d[2] - h * a[2]];
            let r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            residuals[i] = r * cos_a - h * sin_a;
        }
        residuals
    }
}

impl LeastSquaresProblem<f64, Dyn, Const<6>> for ConeLMProblem {
    type ParameterStorage = Owned<f64, Const<6>>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Const<6>>;

    fn set_params(&mut self, x: &OVector<f64, Const<6>>) {
        self.params.copy_from(x);
    }

    fn params(&self) -> OVector<f64, Const<6>> {
        self.params
    }

    fn residuals(&self) -> Option<OVector<f64, Dyn>> {
        Some(self.compute_residuals_for(&self.params))
    }

    fn jacobian(&self) -> Option<OMatrix<f64, Dyn, Const<6>>> {
        let n = self.points.len();
        let eps = 1e-8;
        let mut jac = OMatrix::<f64, Dyn, Const<6>>::zeros_generic(
            Dyn(n), Const::<6>,
        );
        for j in 0..6 {
            let mut pp = self.params;
            pp[j] += eps;
            let rp = self.compute_residuals_for(&pp);
            let mut pm = self.params;
            pm[j] -= eps;
            let rm = self.compute_residuals_for(&pm);
            for i in 0..n {
                jac[(i, j)] = (rp[i] - rm[i]) / (2.0 * eps);
            }
        }
        Some(jac)
    }
}

/// Compute the cone axis direction from face normals and/or vertex positions.
///
/// Strategy 1 (many faces): covariance smallest-eigenvector.
/// Strategy 2 (few faces with shared vertex): apex→base-centroid direction.
/// Strategy 3 (few faces, no shared vertex): cross product of normal differences.
fn compute_cone_axis(
    face_indices: &[usize],
    vertex_set: &HashSet<usize>,
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<[f64; 3]> {
    // With many faces, covariance is robust.
    if face_indices.len() >= 8 {
        let cov = build_normal_covariance(face_indices, faces, vertices);
        let mut axis = smallest_eigenvector_3x3(&cov);
        let len = normalize3(&mut axis);
        if len > 1e-15 {
            return Some(axis);
        }
    }

    // Strategy 2: check if faces share a common vertex (likely the apex).
    // Count how many faces each vertex appears in.
    let mut vertex_face_count = std::collections::HashMap::<usize, usize>::new();
    for &fi in face_indices {
        let vc = faces[fi].vertex_count as usize;
        for &vi in &faces[fi].vertex_indices[..vc] {
            *vertex_face_count.entry(vi).or_insert(0) += 1;
        }
    }
    // If a vertex appears in all (or most) faces, treat it as the potential apex
    let threshold = face_indices.len().max(2) - 1;
    let mut apex_candidates: Vec<usize> = vertex_face_count.iter()
        .filter(|(_, &c)| c >= threshold)
        .map(|(&vi, _)| vi)
        .collect();
    apex_candidates.sort(); // deterministic

    for &apex_vi in &apex_candidates {
        // Compute centroid of NON-apex vertices
        let other_verts: Vec<usize> = vertex_set.iter()
            .filter(|&&vi| vi != apex_vi)
            .copied()
            .collect();
        if other_verts.is_empty() { continue; }
        let mut bx = 0.0_f64;
        let mut by = 0.0_f64;
        let mut bz = 0.0_f64;
        for &vi in &other_verts {
            bx += vertices[vi].x;
            by += vertices[vi].y;
            bz += vertices[vi].z;
        }
        let n = other_verts.len() as f64;
        bx /= n; by /= n; bz /= n;
        // Axis from apex to base centroid
        let mut axis = [
            bx - vertices[apex_vi].x,
            by - vertices[apex_vi].y,
            bz - vertices[apex_vi].z,
        ];
        let len = normalize3(&mut axis);
        if len > 1e-10 {
            return Some(axis);
        }
    }

    // Strategy 3: cross product of normal differences.
    let normals: Vec<[f64; 3]> = face_indices.iter()
        .filter_map(|&fi| faces[fi].normal)
        .collect();
    if normals.len() < 2 {
        return None;
    }
    // Find the two most different normals
    let mut best_diff_sq = 0.0_f64;
    let mut best_i = 0;
    let mut best_j = 1;
    for i in 0..normals.len() {
        for j in (i+1)..normals.len() {
            let d = [normals[i][0]-normals[j][0], normals[i][1]-normals[j][1], normals[i][2]-normals[j][2]];
            let dsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            if dsq > best_diff_sq {
                best_diff_sq = dsq;
                best_i = i;
                best_j = j;
            }
        }
    }
    if best_diff_sq < 1e-10 {
        return None;
    }
    let d1 = [normals[best_i][0]-normals[best_j][0],
               normals[best_i][1]-normals[best_j][1],
               normals[best_i][2]-normals[best_j][2]];

    // Try cross(d1, d2) for a third normal
    for k in 0..normals.len() {
        if k == best_i || k == best_j { continue; }
        let d2 = [normals[best_i][0]-normals[k][0],
                   normals[best_i][1]-normals[k][1],
                   normals[best_i][2]-normals[k][2]];
        let mut c = cross3(&d1, &d2);
        let len = normalize3(&mut c);
        if len > 1e-10 {
            return Some(c);
        }
    }

    // Only 2 distinct normals: axis perpendicular to d1, in mean normal plane
    let avg = [normals[best_i][0]+normals[best_j][0],
               normals[best_i][1]+normals[best_j][1],
               normals[best_i][2]+normals[best_j][2]];
    let mut axis = cross3(&d1, &avg);
    let len = normalize3(&mut axis);
    if len > 1e-15 {
        Some(axis)
    } else {
        None
    }
}

/// Fit cone parameters (axis direction, apex, half-angle) from a set of faces.
/// Returns (apex, axis_direction, half_angle).
fn fit_cone(
    face_indices: &[usize],
    vertex_set: &HashSet<usize>,
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64)> {
    if vertex_set.len() < 4 {
        return None;
    }

    // Determine initial axis direction using multiple strategies.
    let mut axis_dir = compute_cone_axis(face_indices, vertex_set, faces, vertices)?;

    // Fit cone profile to get initial apex and half-angle.
    let (mut apex, half_angle) = fit_cone_profile(&axis_dir, vertex_set, vertices, face_indices, faces)?;

    // Ensure axis points from apex toward the bulk of vertices (positive h side).
    // vertex_to_cone_distance assumes h = (v - apex) · axis > 0 for points on the cone.
    let mut h_sum = 0.0_f64;
    for &vi in vertex_set.iter() {
        let v = &vertices[vi];
        h_sum += (v.x - apex[0]) * axis_dir[0]
               + (v.y - apex[1]) * axis_dir[1]
               + (v.z - apex[2]) * axis_dir[2];
    }
    if h_sum < 0.0 {
        axis_dir = [-axis_dir[0], -axis_dir[1], -axis_dir[2]];
        // Re-fit profile with corrected axis (apex position may adjust)
        if let Some((new_apex, new_ha)) = fit_cone_profile(&axis_dir, vertex_set, vertices, face_indices, faces) {
            apex = new_apex;
            // half_angle should be the same
            let _ = new_ha;
        }
    }

    // Only use LM refinement when we have enough vertices for the 6-param problem.
    // With < 10 vertices, the LM is underconstrained and produces bad results.
    if vertex_set.len() >= 10 {
        let points: Vec<[f64; 3]> = vertex_set.iter()
            .map(|&vi| [vertices[vi].x, vertices[vi].y, vertices[vi].z])
            .collect();
        let problem = ConeLMProblem::new(points, axis_dir, apex, half_angle);
        let (result, _report) = LevenbergMarquardt::new().minimize(problem);
        let refined_dir = result.axis_dir_from_params(&result.params);
        let refined_apex = [result.params[2], result.params[3], result.params[4]];
        let refined_half_angle = result.params[5];

        if refined_half_angle > 0.005 && refined_half_angle < 85.0_f64.to_radians() {
            return Some((refined_apex, refined_dir, refined_half_angle));
        }
    }

    Some((apex, axis_dir, half_angle))
}

/// Check if all vertices in a set are within tolerance of a cone surface.
fn all_vertices_within_cone_tolerance(
    vertex_set: &HashSet<usize>,
    apex: &[f64; 3],
    axis_direction: &[f64; 3],
    half_angle: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        let d = vertex_to_cone_distance(&vertices[vi], apex, axis_direction, half_angle);
        if d.abs() > tolerance {
            return false;
        }
    }
    true
}

/// Determine convexity: does the face normal point away from the axis (convex) or toward it (concave)?
fn determine_convexity(
    face: &MeshFace,
    vertices: &[MeshVertex],
    apex: &[f64; 3],
    axis_direction: &[f64; 3],
) -> bool {
    let centroid = face_centroid(face, vertices);
    let d = [
        centroid[0] - apex[0],
        centroid[1] - apex[1],
        centroid[2] - apex[2],
    ];
    let t = dot3(&d, axis_direction);
    let radial = [
        d[0] - t * axis_direction[0],
        d[1] - t * axis_direction[1],
        d[2] - t * axis_direction[2],
    ];
    let n = face.normal.unwrap();
    dot3(&n, &radial) > 0.0
}

/// Validate that faces in a conical hypothesis have sufficient angular
/// coverage around the cone — same algorithm as cylindrical.
fn angular_coverage_valid(
    face_list: &[usize],
    apex: &[f64; 3],
    axis_direction: &[f64; 3],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> bool {
    if face_list.len() < 3 {
        return false;
    }

    let (u, w) = perpendicular_basis(axis_direction);

    let mut thetas: Vec<f64> = Vec::with_capacity(face_list.len());
    for &fi in face_list {
        let c = face_centroid(&faces[fi], vertices);
        let d = [c[0] - apex[0], c[1] - apex[1], c[2] - apex[2]];
        let t = dot3(&d, axis_direction);
        let radial = [
            d[0] - t * axis_direction[0],
            d[1] - t * axis_direction[1],
            d[2] - t * axis_direction[2],
        ];
        let theta = f64::atan2(dot3(&radial, &w), dot3(&radial, &u));
        thetas.push(theta);
    }

    thetas.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = thetas.len();

    let mut gaps: Vec<f64> = Vec::with_capacity(n);
    for i in 1..n {
        gaps.push(thetas[i] - thetas[i - 1]);
    }
    gaps.push(thetas[0] + 2.0 * std::f64::consts::PI - thetas[n - 1]);

    let mut largest = 0.0_f64;
    let mut second_largest = 0.0_f64;
    for &g in &gaps {
        if g > largest {
            second_largest = largest;
            largest = g;
        } else if g > second_largest {
            second_largest = g;
        }
    }

    let span = 2.0 * std::f64::consts::PI - largest;
    if span <= 0.0 {
        return false;
    }

    second_largest <= span / 3.0
}

// ---------------------------------------------------------------------------
// Trial BFS for conical hypothesis evaluation
// ---------------------------------------------------------------------------

struct ConeTrialResult {
    faces: Vec<usize>,
    vertices: HashSet<usize>,
    apex: [f64; 3],
    axis_direction: [f64; 3],
    half_angle: f64,
    convex: bool,
    error_max: f64,
    centroid_error_max: f64,
    error_abs_sum: f64,
    _total_area: f64,
}

/// Run a trial BFS for a conical hypothesis starting from seed faces.
/// Returns None if the trial fails validation.
#[allow(clippy::too_many_arguments)]
fn run_cone_trial_bfs(
    seed_faces: &[usize],
    mesh: &ConnectedMesh,
    _vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
    verbosity: u8,
    viz: Option<&VizSender>,
    viz_quit: &std::cell::Cell<bool>,
) -> Option<ConeTrialResult> {
    // Collect seed vertices
    let mut vertex_set: HashSet<usize> = HashSet::new();
    for &sfi in seed_faces {
        let vc = mesh.faces[sfi].vertex_count as usize;
        for &vi in &mesh.faces[sfi].vertex_indices[..vc] {
            vertex_set.insert(vi);
        }
    }

    let mut face_list: Vec<usize> = seed_faces.to_vec();

    // Fit cone to seed faces
    let (mut current_apex, mut current_dir, mut current_half_angle) = match
        fit_cone(&face_list, &vertex_set, &mesh.faces, &mesh.vertices)
    {
        Some(v) => v,
        None => {
            if verbosity >= 3 {
                let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
                eprintln!("  [BFS-cone] seed=({}) fit_cone returned None", seed_str.join(","));
            }
            return None;
        }
    };

    if verbosity >= 3 {
        eprintln!("  [BFS-cone] fit result: ha={:.4}° apex=[{:.4},{:.4},{:.4}] dir=[{:.4},{:.4},{:.4}]",
            current_half_angle.to_degrees(),
            current_apex[0], current_apex[1], current_apex[2],
            current_dir[0], current_dir[1], current_dir[2]);
    }

    // Verify seed: all vertices within tolerance
    // Use a relaxed tolerance for seed validation because the initial axis
    // estimate from 3 faces is rough. BFS expansion and re-fit will refine.
    let seed_tol = surface_tol * 5.0;
    if !all_vertices_within_cone_tolerance(
        &vertex_set, &current_apex, &current_dir, current_half_angle,
        seed_tol, &mesh.vertices,
    ) {
        if verbosity >= 3 {
            let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
            let worst_dist = vertex_set.iter().map(|&vi| {
                vertex_to_cone_distance(&mesh.vertices[vi], &current_apex, &current_dir, current_half_angle).abs()
            }).fold(0.0_f64, f64::max);
            eprintln!("  [BFS-cone] seed=({}) vertex tolerance failed: worst_dist={:.2e} > tol={:.2e}",
                seed_str.join(","), worst_dist, seed_tol);
        }
        return None;
    }

    // Determine convexity from first seed face
    let convex = determine_convexity(
        &mesh.faces[seed_faces[0]], &mesh.vertices, &current_apex, &current_dir,
    );

    // NOTE: Skip convexity consistency check on seed faces because the initial
    // axis estimate may be inaccurate (especially with only 3 faces). Convexity
    // will be re-evaluated after BFS expansion and final re-fit.

    // Verify all seed face centroids within relaxed tolerance
    for &sfi in seed_faces {
        let centroid = face_centroid(&mesh.faces[sfi], &mesh.vertices);
        let cen_dist = vertex_to_cone_distance(
            &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
            &current_apex, &current_dir, current_half_angle,
        ).abs();
        if cen_dist > seed_tol {
            if verbosity >= 3 {
                let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
                eprintln!("  [BFS-cone] seed=({}) centroid tolerance failed: face {} cen_dist={:.2e}",
                    seed_str.join(","), sfi, cen_dist);
            }
            return None;
        }
    }

    // Level-3 trace
    if verbosity >= 3 {
        let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
        eprintln!(
            "[BFS-cone] Trial seed=({}) → half_angle={:.4}°, apex=[{:.4},{:.4},{:.4}], dir=[{:.4},{:.4},{:.4}], {}",
            seed_str.join(","),
            current_half_angle.to_degrees(),
            current_apex[0], current_apex[1], current_apex[2],
            current_dir[0], current_dir[1], current_dir[2],
            if convex { "convex" } else { "concave" },
        );
    }

    // Viz: show seed faces
    let mut skip_viz = false;
    {
        let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
        let mut highlights = vec![
            viz::FaceHighlight { face_indices: vec![seed_faces[0]], color: [1.0, 0.6, 0.0, 1.0] },
        ];
        if seed_faces.len() > 1 {
            highlights.push(viz::FaceHighlight { face_indices: seed_faces[1..].to_vec(), color: [1.0, 0.8, 0.4, 1.0] });
        }
        if let Some(action) = viz_custom(
            viz, highlights,
            Vec::new(),
            &format!("BFS-cone: seed=({}) θ={:.2}° {} [space=step, shift+space=skip]",
                seed_str.join(","), current_half_angle.to_degrees(),
                if convex { "convex" } else { "concave" }),
            Vec::new(), Vec::new(),
            Some(viz_face_centroid(seed_faces[0], mesh)),
            viz_face_normal(seed_faces[0], mesh),
        ) {
            match action {
                VizAction::Quit => { viz_quit.set(true); return None; }
                VizAction::NextSeed => { skip_viz = true; }
                VizAction::NextStep => {}
            }
        }
    }

    // Track claimed faces in this trial
    let mut trial_claimed: HashSet<usize> = HashSet::new();
    for &sfi in seed_faces {
        trial_claimed.insert(sfi);
    }

    let mut queue = VecDeque::new();
    for &sfi in seed_faces {
        queue.push_back(sfi);
    }

    // BFS expansion
    while let Some(current_fi) = queue.pop_front() {
        let vc = mesh.faces[current_fi].vertex_count as usize;
        let neighbors = mesh.faces[current_fi].neighbors;

        for &cni in &neighbors[..vc] {
            if cni < 0 {
                continue;
            }
            let cni = cni as usize;

            // Skip if already committed to a conical hypothesis or claimed by this trial
            if mesh.faces[cni].conical_hypothesis != UNDEDUCED_CONICAL_HYPOTHESIS {
                continue;
            }
            if trial_claimed.contains(&cni) {
                continue;
            }

            // Skip convexity check during BFS expansion. The initial axis
            // estimate may be inaccurate, making convexity detection unreliable.
            // Convexity is determined after the final re-fit.

            // Angular tolerance check — use 2× angular_tol for cones because
            // tessellated cones have inter-strip dihedral angles that can exceed
            // the default tolerance (e.g., 18° for 20 circumferential divisions
            // vs. 17.5° default). The vertex-to-cone distance check below is the
            // primary discriminator for cone membership.
            let cone_angular_tol = angular_tol * 2.0;
            if let Some(n_cni) = mesh.faces[cni].normal {
                let cni_vc2 = mesh.faces[cni].vertex_count as usize;
                let cni_neighbors = mesh.faces[cni].neighbors;
                let mut angular_reject = false;
                for &adj in &cni_neighbors[..cni_vc2] {
                    if adj < 0 { continue; }
                    let adj = adj as usize;
                    if !trial_claimed.contains(&adj) { continue; }
                    if let Some(n_adj) = mesh.faces[adj].normal {
                        let cos_a = dot3(&n_adj, &n_cni).clamp(-1.0, 1.0);
                        let angle = cos_a.acos();
                        if angle > cone_angular_tol {
                            angular_reject = true;
                            break;
                        }
                    }
                }
                if angular_reject {
                    if verbosity >= 3 {
                        eprintln!("  [BFS-cone] from fi={} try cni={}: angular reject", current_fi, cni);
                    }
                    continue;
                }
            }

            // Vertex distance check
            let cni_vc = mesh.faces[cni].vertex_count as usize;
            let cni_vi = mesh.faces[cni].vertex_indices;
            let mut all_ok = true;
            let mut any_far = false;
            let mut vtx_err_max = 0.0_f64;
            for &vi in &cni_vi[..cni_vc] {
                let d = vertex_to_cone_distance(
                    &mesh.vertices[vi],
                    &current_apex,
                    &current_dir,
                    current_half_angle,
                ).abs();
                vtx_err_max = vtx_err_max.max(d);
                if d > surface_tol {
                    all_ok = false;
                    if d > REFIT_SKIP_MULTIPLIER * surface_tol {
                        any_far = true;
                        break;
                    }
                }
            }

            if any_far {
                if verbosity >= 3 {
                    eprintln!("  [BFS-cone] from fi={} try cni={}: too far vtx_err={:.2e} → REJECT",
                        current_fi, cni, vtx_err_max);
                }
                continue;
            }

            // Centroid distance check
            let centroid = face_centroid(&mesh.faces[cni], &mesh.vertices);
            let centroid_dist = vertex_to_cone_distance(
                &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                &current_apex, &current_dir, current_half_angle,
            ).abs();
            if centroid_dist > REFIT_SKIP_MULTIPLIER * surface_tol {
                if verbosity >= 3 {
                    eprintln!("  [BFS-cone] from fi={} try cni={}: centroid far {:.2e} → REJECT",
                        current_fi, cni, centroid_dist);
                }
                continue;
            }
            let needs_refit = !all_ok || centroid_dist > surface_tol;

            if needs_refit {
                let mut trial_vertices = vertex_set.clone();
                for &vi in &cni_vi[..cni_vc] {
                    trial_vertices.insert(vi);
                }
                let mut trial_faces = face_list.clone();
                trial_faces.push(cni);

                let refit = fit_cone(
                    &trial_faces, &trial_vertices, &mesh.faces, &mesh.vertices,
                );
                let (new_apex, new_dir, new_ha) = match refit {
                    Some(params) => params,
                    None => {
                        if verbosity >= 3 {
                            eprintln!("  [BFS-cone] from fi={} try cni={}: refit failed → REJECT",
                                current_fi, cni);
                        }
                        continue;
                    }
                };

                if !all_vertices_within_cone_tolerance(
                    &trial_vertices, &new_apex, &new_dir, new_ha,
                    surface_tol, &mesh.vertices,
                ) {
                    if verbosity >= 3 {
                        eprintln!("  [BFS-cone] from fi={} try cni={}: refit vertex tol failed → REJECT",
                            current_fi, cni);
                    }
                    continue;
                }

                // Check all face centroids within surface tolerance after re-fit
                let mut centroids_ok = true;
                for &f in &trial_faces {
                    let c = face_centroid(&mesh.faces[f], &mesh.vertices);
                    let d = vertex_to_cone_distance(
                        &MeshVertex::from_xyz(c[0], c[1], c[2]),
                        &new_apex, &new_dir, new_ha,
                    ).abs();
                    if d > surface_tol {
                        centroids_ok = false;
                        break;
                    }
                }
                if !centroids_ok {
                    if verbosity >= 3 {
                        eprintln!("  [BFS-cone] from fi={} try cni={}: refit centroid check failed → REJECT",
                            current_fi, cni);
                    }
                    continue;
                }

                if verbosity >= 3 {
                    eprintln!("  [BFS-cone] from fi={} try cni={}: refit ok θ={:.4}° → ACCEPT[refit]",
                        current_fi, cni, new_ha.to_degrees());
                }

                current_apex = new_apex;
                current_dir = new_dir;
                current_half_angle = new_ha;
            } else if verbosity >= 3 {
                eprintln!("  [BFS-cone] from fi={} try cni={}: vtx_err={:.2e} → ACCEPT",
                    current_fi, cni, vtx_err_max);
            }

            // Accept this face
            trial_claimed.insert(cni);
            face_list.push(cni);
            for &vi in &cni_vi[..cni_vc] {
                vertex_set.insert(vi);
            }
            queue.push_back(cni);

            // Viz: show accepted face
            if !skip_viz {
                let accepted_nonseed: Vec<usize> = face_list.iter()
                    .filter(|f| !seed_faces.contains(f) && **f != cni)
                    .copied().collect();
                let mut highlights = vec![
                    viz::FaceHighlight { face_indices: vec![seed_faces[0]], color: [1.0, 0.6, 0.0, 1.0] },
                ];
                if seed_faces.len() > 1 {
                    highlights.push(viz::FaceHighlight { face_indices: seed_faces[1..].to_vec(), color: [1.0, 0.8, 0.4, 1.0] });
                }
                if !accepted_nonseed.is_empty() {
                    highlights.push(viz::FaceHighlight { face_indices: accepted_nonseed, color: [0.2, 0.4, 1.0, 1.0] });
                }
                highlights.push(viz::FaceHighlight { face_indices: vec![cni], color: [0.1, 0.2, 0.7, 1.0] });
                // Gray background for committed conical faces
                let bg_faces: Vec<usize> = (0..mesh.faces.len())
                    .filter(|f| mesh.faces[*f].conical_hypothesis >= 0 && !face_list.contains(f))
                    .collect();
                if !bg_faces.is_empty() {
                    highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                }
                if let Some(action) = viz_custom(
                    viz, highlights,
                    Vec::new(),
                    &format!("BFS-cone: accepted fi={cni} ({} faces) θ={:.2}° [space=step, shift+space=skip]",
                        face_list.len(), current_half_angle.to_degrees()),
                    Vec::new(), Vec::new(),
                    Some(viz_face_centroid(cni, mesh)),
                    viz_face_normal(cni, mesh),
                ) {
                    match action {
                        VizAction::Quit => { viz_quit.set(true); return None; }
                        VizAction::NextSeed => { skip_viz = true; }
                        VizAction::NextStep => {}
                    }
                }
            }
        }
    }

    // Final re-fit
    if let Some((final_apex, final_dir, final_ha)) =
        fit_cone(&face_list, &vertex_set, &mesh.faces, &mesh.vertices)
    {
        current_apex = final_apex;
        current_dir = final_dir;
        current_half_angle = final_ha;
    }

    // Compute error metrics
    let mut error_max = 0.0_f64;
    let mut error_abs_sum = 0.0_f64;
    for &vi in &vertex_set {
        let d = vertex_to_cone_distance(
            &mesh.vertices[vi], &current_apex, &current_dir, current_half_angle,
        ).abs();
        error_max = error_max.max(d);
        error_abs_sum += d;
    }

    if verbosity >= 3 {
        eprintln!("  [BFS-cone] post-refit: {} faces, error_max={:.2e}, ha={:.4}° dir=[{:.4},{:.4},{:.4}]",
            face_list.len(), error_max, current_half_angle.to_degrees(),
            current_dir[0], current_dir[1], current_dir[2]);
    }

    // Minimum face count
    if face_list.len() < MIN_CONE_FACES {
        if verbosity >= 3 {
            eprintln!("  [BFS-cone] rejected: {} faces < MIN_CONE_FACES={}", face_list.len(), MIN_CONE_FACES);
        }
        return None;
    }

    // Centroid check and compute centroid_error_max
    let mut centroid_error_max = 0.0_f64;
    for &f in &face_list {
        let c = face_centroid(&mesh.faces[f], &mesh.vertices);
        let d = vertex_to_cone_distance(
            &MeshVertex::from_xyz(c[0], c[1], c[2]),
            &current_apex, &current_dir, current_half_angle,
        ).abs();
        centroid_error_max = centroid_error_max.max(d);
        if d > surface_tol {
            if verbosity >= 3 {
                eprintln!("  [BFS-cone] rejected: face {} centroid_err={:.2e} > surface_tol={:.2e}", f, d, surface_tol);
            }
            return None;
        }
    }

    let total_area: f64 = face_list.iter().map(|&f| face_area(&mesh.faces[f], &mesh.vertices)).sum();

    // Determine convexity from the final fit using majority vote
    let mut convex_votes = 0_i32;
    for &f in &face_list {
        if determine_convexity(&mesh.faces[f], &mesh.vertices, &current_apex, &current_dir) {
            convex_votes += 1;
        } else {
            convex_votes -= 1;
        }
    }
    let convex = convex_votes > 0;

    Some(ConeTrialResult {
        faces: face_list,
        vertices: vertex_set,
        apex: current_apex,
        axis_direction: current_dir,
        half_angle: current_half_angle,
        convex,
        error_max,
        centroid_error_max,
        error_abs_sum,
        _total_area: total_area,
    })
}

// ---------------------------------------------------------------------------
// Stage 2.4: Conical hypothesis deduction
// ---------------------------------------------------------------------------

/// Deduce conical hypotheses from the mesh using apex-vertex seeding.
///
/// Instead of triple-seed (which fails for cones because 3 adjacent apex faces
/// don't give enough geometric spread), we identify potential apex vertices —
/// vertices shared by many non-coplanar faces — and use all faces incident to
/// each apex vertex as the seed set. This gives robust axis estimation.
pub(super) fn deduce_conical_hypotheses(
    mesh: &mut ConnectedMesh,
    vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
    verbosity: u8,
    viz: Option<&VizSender>,
) -> (Vec<ConicalHypothesis>, bool) {
    let num_faces = mesh.faces.len();
    let num_vertices = mesh.vertices.len();
    let mut hypotheses: Vec<ConicalHypothesis> = Vec::new();
    let viz_quit = std::cell::Cell::new(false);

    // Compute bounding box diagonal and center for apex distance sanity check
    let bb_diag = bounding_box_diagonal(&mesh.vertices);
    let bb_center = {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for v in &mesh.vertices {
            min[0] = min[0].min(v.x); max[0] = max[0].max(v.x);
            min[1] = min[1].min(v.y); max[1] = max[1].max(v.y);
            min[2] = min[2].min(v.z); max[2] = max[2].max(v.z);
        }
        [(min[0]+max[0])/2.0, (min[1]+max[1])/2.0, (min[2]+max[2])/2.0]
    };

    // Build vertex → face adjacency
    let mut vertex_faces: Vec<Vec<usize>> = vec![Vec::new(); num_vertices];
    for fi in 0..num_faces {
        let vc = mesh.faces[fi].vertex_count as usize;
        for &vi in &mesh.faces[fi].vertex_indices[..vc] {
            vertex_faces[vi].push(fi);
        }
    }

    // Collect apex candidates: vertices with many incident faces that have
    // non-coplanar normals (indicating a cone tip, not a flat fan).
    // Sort by face count descending so we try the best candidates first.
    let mut apex_candidates: Vec<(usize, usize)> = Vec::new(); // (vertex_index, face_count)
    for (vi, vi_faces) in vertex_faces.iter().enumerate() {
        let incident: Vec<usize> = vi_faces.iter()
            .copied()
            .filter(|&fi| {
                mesh.faces[fi].conical_hypothesis == UNDEDUCED_CONICAL_HYPOTHESIS
                    && mesh.faces[fi].normal.is_some()
            })
            .collect();
        if incident.len() < MIN_CONE_FACES {
            continue;
        }
        // Check that normals are non-coplanar: compute normal covariance and
        // check that the smallest eigenvalue is nonzero (normals span 3D).
        let cov = build_normal_covariance(&incident, &mesh.faces, &mesh.vertices);
        let axis = smallest_eigenvector_3x3(&cov);
        // cov is packed symmetric [m00, m01, m02, m11, m12, m22]
        let eigenval = cov[0] * axis[0] * axis[0]
            + cov[3] * axis[1] * axis[1]
            + cov[5] * axis[2] * axis[2]
            + 2.0 * cov[1] * axis[0] * axis[1]
            + 2.0 * cov[2] * axis[0] * axis[2]
            + 2.0 * cov[4] * axis[1] * axis[2];
        // For a cone apex, the smallest eigenvalue should be significantly smaller
        // than the others (normals spread in a plane perpendicular to axis).
        // For a flat region, all eigenvalues are ~0 or normals are parallel.
        let trace = cov[0] + cov[3] + cov[5];
        if trace < 1e-10 || eigenval / trace > 0.3 {
            continue; // Normals too parallel or too isotropic
        }
        apex_candidates.push((vi, incident.len()));
    }
    apex_candidates.sort_by(|a, b| b.1.cmp(&a.1));

    if verbosity >= 3 {
        eprintln!("  [cone] {} apex candidates (from {} vertices)",
            apex_candidates.len(), num_vertices);
    }

    for &(apex_vi, _) in &apex_candidates {
        // Collect all undeduced faces incident to this apex vertex
        let seed_faces: Vec<usize> = vertex_faces[apex_vi].iter()
            .copied()
            .filter(|&fi| mesh.faces[fi].conical_hypothesis == UNDEDUCED_CONICAL_HYPOTHESIS
                && mesh.faces[fi].normal.is_some())
            .collect();
        if seed_faces.len() < MIN_CONE_FACES {
            continue;
        }

        if verbosity >= 3 {
            eprintln!("  [cone] trying apex vertex {} with {} incident faces",
                apex_vi, seed_faces.len());
        }

        let trial = run_cone_trial_bfs(
            &seed_faces, mesh,
            vertex_tol, surface_tol, angular_tol,
            verbosity,
            viz, &viz_quit,
        );

        if viz_quit.get() {
            return (hypotheses, true);
        }

        // Check angular coverage
        let trial = trial.and_then(|t| {
            if angular_coverage_valid(
                &t.faces, &t.apex, &t.axis_direction,
                &mesh.faces, &mesh.vertices,
            ) {
                Some(t)
            } else {
                if verbosity >= 3 {
                    eprintln!("  [cone] apex vertex {} angular coverage failed: {} faces",
                        apex_vi, t.faces.len());
                }
                None
            }
        });

        // Reject degenerate/false-positive cones
        let trial = trial.and_then(|t| {
            let ha_deg = t.half_angle.to_degrees();
            if !(MIN_HALF_ANGLE_DEG..=MAX_HALF_ANGLE_DEG).contains(&ha_deg) {
                if verbosity >= 3 {
                    eprintln!("  [cone] apex vertex {} rejected: half_angle={:.2}° outside [{},{}]",
                        apex_vi, ha_deg, MIN_HALF_ANGLE_DEG, MAX_HALF_ANGLE_DEG);
                }
                return None;
            }
            if t.error_max > surface_tol {
                if verbosity >= 3 {
                    eprintln!("  [cone] apex vertex {} rejected: error_max={:.2e} > surface_tol={:.2e}",
                        apex_vi, t.error_max, surface_tol);
                }
                return None;
            }
            // Reject cones whose apex is absurdly far from the mesh
            let dx = t.apex[0] - bb_center[0];
            let dy = t.apex[1] - bb_center[1];
            let dz = t.apex[2] - bb_center[2];
            let apex_dist = (dx*dx + dy*dy + dz*dz).sqrt();
            let max_apex_dist = bb_diag * MAX_APEX_DISTANCE_FACTOR;
            if apex_dist > max_apex_dist {
                if verbosity >= 3 {
                    eprintln!("  [cone] apex vertex {} rejected: apex_dist={:.1} > max={:.1} ({}x bb_diag)",
                        apex_vi, apex_dist, max_apex_dist, MAX_APEX_DISTANCE_FACTOR);
                }
                return None;
            }
            // Normal-axis consistency: on a true cone, every face normal makes
            // angle (90° - half_angle) with the axis. Compute the standard
            // deviation of acos(|n·axis|) and reject if too large.
            let mut angle_sum = 0.0_f64;
            let mut angle_sq_sum = 0.0_f64;
            let mut n_count = 0usize;
            for &fi in &t.faces {
                if let Some(n) = mesh.faces[fi].normal {
                    let cos_na = dot3(&n, &t.axis_direction).clamp(-1.0, 1.0);
                    let angle = cos_na.abs().acos(); // angle between normal and axis
                    angle_sum += angle;
                    angle_sq_sum += angle * angle;
                    n_count += 1;
                }
            }
            if n_count >= 2 {
                let mean = angle_sum / n_count as f64;
                let variance = angle_sq_sum / n_count as f64 - mean * mean;
                let std_dev = variance.max(0.0).sqrt();
                // For a true cone, std_dev should be ~0. For a sphere patch
                // approximated as a cone, normals vary and std_dev is larger.
                // Use angular_tol as the threshold (17.5° = 0.305 rad).
                let max_std_dev = angular_tol * 0.5;
                if std_dev > max_std_dev {
                    if verbosity >= 3 {
                        eprintln!("  [cone] apex vertex {} rejected: normal-axis std_dev={:.4} rad ({:.2}°) > max={:.4}",
                            apex_vi, std_dev, std_dev.to_degrees(), max_std_dev);
                    }
                    return None;
                }
            }
            Some(t)
        });

        if let Some(candidate) = trial {
            // Viz: post-BFS pause — show accepted hypothesis
            if !viz_quit.get() {
                if let Some(viz) = viz {
                    if verbosity >= 2 {
                        eprintln!("  [BFS-cone] ACCEPTED: {} faces, θ={:.4}°, err_max={:.2e}",
                            candidate.faces.len(), candidate.half_angle.to_degrees(), candidate.error_max);
                    }
                    let mut highlights = vec![
                        viz::FaceHighlight { face_indices: candidate.faces.clone(), color: [0.0, 0.8, 0.0, 1.0] },
                    ];
                    let bg_faces: Vec<usize> = (0..mesh.faces.len())
                        .filter(|f| mesh.faces[*f].conical_hypothesis >= 0 && !candidate.faces.contains(f))
                        .collect();
                    if !bg_faces.is_empty() {
                        highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                    }
                    if let Some(action) = viz_custom(
                        Some(viz), highlights,
                        Vec::new(),
                        &format!("BFS-cone result: ACCEPTED {} faces θ={:.2}° [space=next]",
                            candidate.faces.len(), candidate.half_angle.to_degrees()),
                        Vec::new(), Vec::new(),
                        Some(viz_face_centroid(candidate.faces[0], mesh)),
                        viz_face_normal(candidate.faces[0], mesh),
                    ) {
                        match action {
                            VizAction::Quit => { return (hypotheses, true); }
                            VizAction::NextSeed | VizAction::NextStep => {}
                        }
                    }
                }
            }

            // Commit: assign all faces to this hypothesis
            let hi = hypotheses.len() as i32;
            for &f in &candidate.faces {
                mesh.faces[f].conical_hypothesis = hi;
            }

            hypotheses.push(ConicalHypothesis {
                apex: candidate.apex,
                axis_direction: candidate.axis_direction,
                half_angle: candidate.half_angle,
                convex: candidate.convex,
                faces: candidate.faces,
                vertices: candidate.vertices.into_iter().collect(),
                error_max: candidate.error_max,
                centroid_error_max: candidate.centroid_error_max,
                error_abs_sum: candidate.error_abs_sum,
            });
        }
    }

    // Set remaining UNDEDUCED faces to NO_HYPOTHESIS
    for fi in 0..num_faces {
        if mesh.faces[fi].conical_hypothesis == UNDEDUCED_CONICAL_HYPOTHESIS {
            mesh.faces[fi].conical_hypothesis = NO_HYPOTHESIS;
        }
    }

    (hypotheses, false)
}

// ---------------------------------------------------------------------------
// Stage 2.4: Compare check
// ---------------------------------------------------------------------------

/// Validate fitted conical hypotheses against a reference STEP shape.
pub(super) fn compare_conical_hypotheses(
    hypotheses: &[ConicalHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);

            // Project centroid onto cone surface
            let d = [
                centroid[0] - hyp.apex[0],
                centroid[1] - hyp.apex[1],
                centroid[2] - hyp.apex[2],
            ];
            let h = dot3(&d, &hyp.axis_direction);
            let radial = [
                d[0] - h * hyp.axis_direction[0],
                d[1] - h * hyp.axis_direction[1],
                d[2] - h * hyp.axis_direction[2],
            ];
            let radial_dist = (radial[0] * radial[0] + radial[1] * radial[1]
                + radial[2] * radial[2]).sqrt();

            let projected = if radial_dist > 1e-15 {
                // Project onto cone surface: at height h, the cone radius is h * tan(half_angle)
                let cone_r = h * hyp.half_angle.tan();
                let scale = cone_r / radial_dist;
                [
                    hyp.apex[0] + h * hyp.axis_direction[0] + radial[0] * scale,
                    hyp.apex[1] + h * hyp.axis_direction[1] + radial[1] * scale,
                    hyp.apex[2] + h * hyp.axis_direction[2] + radial[2] * scale,
                ]
            } else {
                centroid
            };

            let pt = gp::Pnt::new_real3(projected[0], projected[1], projected[2]);
            let dist = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(dist);
        }

        let tolerance = config.surface_tolerance_mm;

        if max_dist > tolerance {
            return Err(Stage2CompareError {
                hypothesis_type: "conical",
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance,
            });
        }
    }

    Ok(())
}
