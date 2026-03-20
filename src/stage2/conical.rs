//! Stage 2.4: Conical hypothesis deduction and comparison.

use std::collections::HashSet;

use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{Const, Dyn, OMatrix, OVector, Owned};
use opencascade_sys::gp;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_CONICAL_HYPOTHESIS};
use crate::viz;

use super::{
    build_normal_covariance, centered_cylinder_overlay, cross3, dot3, face_area, face_centroid,
    normalize3, perpendicular_basis, smallest_eigenvector_3x3, ConicalHypothesis,
    Stage2CompareError, MIN_CROSS_THRESHOLD, MIN_CONE_FACES, MIN_CYLINDER_FACES,
};
use super::cylindrical::{vertex_to_cylinder_distance, fit_cylinder};

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
        let mut jac = OMatrix::<f64, Dyn, Const<6>>::zeros_generic(
            Dyn(n), Const::<6>,
        );
        let a = self.axis_dir_from_params(&self.params);
        let apex = [self.params[2], self.params[3], self.params[4]];
        let half_angle = self.params[5];
        let (sa, ca) = half_angle.sin_cos();

        // L = |a0 + α*u0 + β*w0| (unnormalized axis length)
        let alpha = self.params[0];
        let beta = self.params[1];
        let vx = self.a0[0] + alpha * self.u0[0] + beta * self.w0[0];
        let vy = self.a0[1] + alpha * self.u0[1] + beta * self.w0[1];
        let vz = self.a0[2] + alpha * self.u0[2] + beta * self.w0[2];
        let big_l = (vx * vx + vy * vy + vz * vz).sqrt().max(1e-15);

        for (i, p) in self.points.iter().enumerate() {
            let d = [p[0] - apex[0], p[1] - apex[1], p[2] - apex[2]];
            let h = dot3(&d, &a);
            let lat = [d[0] - h * a[0], d[1] - h * a[1], d[2] - h * a[2]];
            let r = (lat[0] * lat[0] + lat[1] * lat[1] + lat[2] * lat[2]).sqrt().max(1e-15);
            let lat_u = dot3(&lat, &self.u0);
            let lat_w = dot3(&lat, &self.w0);
            let axis_sens = ca * h / r + sa;

            // ∂f/∂α
            jac[(i, 0)] = -lat_u * axis_sens / big_l;
            // ∂f/∂β
            jac[(i, 1)] = -lat_w * axis_sens / big_l;
            // ∂f/∂apex_x, ∂f/∂apex_y, ∂f/∂apex_z
            jac[(i, 2)] = -ca * lat[0] / r + sa * a[0];
            jac[(i, 3)] = -ca * lat[1] / r + sa * a[1];
            jac[(i, 4)] = -ca * lat[2] / r + sa * a[2];
            // ∂f/∂half_angle
            jac[(i, 5)] = -r * sa - h * ca;
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

            let projected = if hyp.is_cylinder() {
                // For cylinders: project onto cylinder surface at radius_at_apex
                if radial_dist > 1e-15 {
                    let scale = hyp.radius_at_apex / radial_dist;
                    [
                        hyp.apex[0] + h * hyp.axis_direction[0] + radial[0] * scale,
                        hyp.apex[1] + h * hyp.axis_direction[1] + radial[1] * scale,
                        hyp.apex[2] + h * hyp.axis_direction[2] + radial[2] * scale,
                    ]
                } else {
                    centroid
                }
            } else if radial_dist > 1e-15 {
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

// ---------------------------------------------------------------------------
// SurfaceType trait implementation for generic BFS
// ---------------------------------------------------------------------------

use super::bfs::{BfsConfig, Seed, SurfaceType};

/// Fitted cone/cylinder parameters.
#[derive(Clone)]
pub(super) struct ConeFit {
    pub apex: [f64; 3],
    pub axis_direction: [f64; 3],
    pub half_angle: f64,
    pub radius: f64,  // cylinder radius when half_angle == 0
    pub convex: Option<bool>,
}

/// Marker struct for the conical/cylindrical surface type.
/// Produces ConicalHypothesis for both cones (half_angle > 0) and
/// cylinders (half_angle == 0).
pub(super) struct ConicalSurface;

/// Check if two faces are pairwise qualified for conical/cylindrical seeding.
/// Both must have normals, the normals must be sufficiently different, and the
/// dihedral angle must be within tolerance.
fn conical_pairwise_qualified(
    mesh: &ConnectedMesh, fi: usize, ni: usize, angular_tol: f64,
) -> bool {
    if mesh.faces[ni].conical_hypothesis != UNDEDUCED_CONICAL_HYPOTHESIS {
        return false;
    }
    let fi_normal = match mesh.faces[fi].normal {
        Some(n) => n,
        None => return false,
    };
    let ni_normal = match mesh.faces[ni].normal {
        Some(n) => n,
        None => return false,
    };
    let cross = cross3(&fi_normal, &ni_normal);
    let cross_mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
    if cross_mag < MIN_CROSS_THRESHOLD {
        return false;
    }
    let cos_angle = dot3(&fi_normal, &ni_normal).clamp(-1.0, 1.0);
    cos_angle.acos() <= angular_tol
}

impl SurfaceType for ConicalSurface {
    type Fit = ConeFit;
    type Hypothesis = ConicalHypothesis;

    fn name() -> &'static str { "conical" }

    fn hypothesis_field(face: &MeshFace) -> i32 { face.conical_hypothesis }

    fn set_hypothesis_field(face: &mut MeshFace, value: i32) {
        face.conical_hypothesis = value;
    }

    fn undeduced_value() -> i32 { UNDEDUCED_CONICAL_HYPOTHESIS }

    fn generate_seeds(mesh: &ConnectedMesh, fi: usize, config: &BfsConfig) -> Vec<Seed> {
        if mesh.faces[fi].normal.is_none() {
            return Vec::new();
        }
        if mesh.faces[fi].conical_hypothesis != UNDEDUCED_CONICAL_HYPOTHESIS {
            return Vec::new();
        }

        let vc = mesh.faces[fi].vertex_count as usize;
        let angular_tol = config.angular_tol;

        // Collect qualified neighbors of fi
        let qualified_neighbors: Vec<usize> = mesh.faces[fi].neighbors[..vc]
            .iter()
            .filter(|&&n| n >= 0)
            .map(|&n| n as usize)
            .filter(|&ni| conical_pairwise_qualified(mesh, fi, ni, angular_tol))
            .collect();

        // Generate 3-face seeds: fi + n1 + n2 where n1 and n2 are qualified
        // neighbors of fi but not neighbors of each other
        let mut seeds = Vec::new();
        for (i, &n1) in qualified_neighbors.iter().enumerate() {
            for &n2 in &qualified_neighbors[(i + 1)..] {
                // n2 must not be a neighbor of n1
                if super::bfs::is_neighbor(mesh, n1, n2) { continue; }
                seeds.push(Seed { faces: vec![fi, n1, n2] });
            }
        }
        seeds
    }

    fn fit_surface(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        vertex_list: &[usize],
        _verbosity: u8,
    ) -> Option<Self::Fit> {
        let vertex_set: HashSet<usize> = vertex_list.iter().copied().collect();

        // Try cylinder first, then cone, pick the one with lower max vertex error
        let cyl_result = fit_cylinder(face_list, &vertex_set, &mesh.faces, &mesh.vertices);
        let cone_result = fit_cone(face_list, &vertex_set, &mesh.faces, &mesh.vertices);

        let cyl_fit = cyl_result.map(|(origin, dir, radius)| {
            ConeFit { apex: origin, axis_direction: dir, half_angle: 0.0, radius, convex: None }
        });

        let cone_fit = cone_result.and_then(|(apex, dir, ha)| {
            let ha_deg = ha.to_degrees();
            if !(MIN_HALF_ANGLE_DEG..=MAX_HALF_ANGLE_DEG).contains(&ha_deg) {
                return None;
            }
            // Reject apex too far from mesh
            let bb_diag = super::bounding_box_diagonal(&mesh.vertices);
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
            let apex_dist = ((apex[0] - bb_center[0]).powi(2)
                + (apex[1] - bb_center[1]).powi(2)
                + (apex[2] - bb_center[2]).powi(2)).sqrt();
            if apex_dist > bb_diag * MAX_APEX_DISTANCE_FACTOR {
                return None;
            }
            Some(ConeFit { apex, axis_direction: dir, half_angle: ha, radius: 0.0, convex: None })
        });

        // Compute max vertex error for each successful fit
        let cyl_max = cyl_fit.as_ref().map(|fit| {
            vertex_list.iter().map(|&vi| {
                let v = &mesh.vertices[vi];
                vertex_to_cylinder_distance(v, &fit.apex, &fit.axis_direction, fit.radius).abs()
            }).fold(0.0_f64, f64::max)
        });

        let cone_max = cone_fit.as_ref().map(|fit| {
            vertex_list.iter().map(|&vi| {
                let v = &mesh.vertices[vi];
                vertex_to_cone_distance(v, &fit.apex, &fit.axis_direction, fit.half_angle).abs()
            }).fold(0.0_f64, f64::max)
        });

        // Pick the better fit
        match (cyl_fit, cone_fit) {
            (Some(cf), Some(kf)) => {
                if cyl_max.unwrap() <= cone_max.unwrap() {
                    Some(cf)
                } else {
                    Some(kf)
                }
            }
            (Some(cf), None) => Some(cf),
            (None, Some(kf)) => Some(kf),
            (None, None) => None,
        }
    }

    fn vertex_distance(fit: &Self::Fit, pos: &[f64; 3]) -> f64 {
        let v = MeshVertex::from_xyz(pos[0], pos[1], pos[2]);
        if fit.half_angle < 1e-4 {
            vertex_to_cylinder_distance(&v, &fit.apex, &fit.axis_direction, fit.radius).abs()
        } else {
            vertex_to_cone_distance(&v, &fit.apex, &fit.axis_direction, fit.half_angle).abs()
        }
    }

    fn convexity_check(fit: &Self::Fit, mesh: &ConnectedMesh, fi: usize) -> bool {
        if let Some(convex) = fit.convex {
            let fi_convex = determine_convexity(
                &mesh.faces[fi], &mesh.vertices, &fit.apex, &fit.axis_direction,
            );
            fi_convex == convex
        } else {
            true
        }
    }

    fn determine_convexity(mesh: &ConnectedMesh, face_list: &[usize], fit: &mut Self::Fit) {
        if let Some(&first) = face_list.first() {
            fit.convex = Some(determine_convexity(
                &mesh.faces[first], &mesh.vertices, &fit.apex, &fit.axis_direction,
            ));
        }
    }

    fn acceptance_test(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        fit: &Self::Fit,
        _config: &BfsConfig,
    ) -> bool {
        let min_faces = if fit.half_angle < 1e-4 { MIN_CYLINDER_FACES } else { MIN_CONE_FACES };
        if face_list.len() < min_faces {
            return false;
        }
        angular_coverage_valid(
            face_list, &fit.apex, &fit.axis_direction,
            &mesh.faces, &mesh.vertices,
        )
    }

    fn viz_overlays(
        fit: &Self::Fit,
        mesh: &ConnectedMesh,
        face_list: &[usize],
    ) -> (Vec<viz::CylinderOverlay>, Vec<viz::SphereOverlay>) {
        let cyl = centered_cylinder_overlay(
            fit.apex, fit.axis_direction, fit.radius,
            face_list, mesh,
            [0.2, 0.4, 1.0, 0.3],
        );
        (vec![cyl], Vec::new())
    }

    fn viz_fit_description(fit: &Self::Fit) -> String {
        if fit.half_angle.abs() < 1e-6 {
            format!("r={:.4}", fit.radius)
        } else {
            format!("r={:.4} half_angle={:.2}°", fit.radius, fit.half_angle.to_degrees())
        }
    }

    fn make_hypothesis(
        mesh: &ConnectedMesh,
        faces: Vec<usize>,
        vertices: Vec<usize>,
        fit: &Self::Fit,
    ) -> ConicalHypothesis {
        let mut error_max = 0.0_f64;
        let mut centroid_error_max = 0.0_f64;
        let mut error_abs_sum = 0.0_f64;
        for &vi in &vertices {
            let d = Self::vertex_distance(fit, &[mesh.vertices[vi].x, mesh.vertices[vi].y, mesh.vertices[vi].z]);
            error_max = error_max.max(d);
            error_abs_sum += d;
        }
        for &fi in &faces {
            let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
            let d = Self::vertex_distance(fit, &c);
            centroid_error_max = centroid_error_max.max(d);
        }
        ConicalHypothesis {
            apex: fit.apex,
            axis_direction: fit.axis_direction,
            half_angle: fit.half_angle,
            radius_at_apex: fit.radius,
            convex: fit.convex.unwrap_or(true),
            faces,
            vertices,
            error_max,
            centroid_error_max,
            error_abs_sum,
        }
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Generate points on a cone surface.
    fn cone_points(
        apex: [f64; 3],
        axis_dir: [f64; 3],
        half_angle: f64,
        n_rings: usize,
        n_pts_per_ring: usize,
        h_min: f64,
        h_max: f64,
    ) -> Vec<MeshVertex> {
        let d = {
            let len = (axis_dir[0].powi(2) + axis_dir[1].powi(2) + axis_dir[2].powi(2)).sqrt();
            [axis_dir[0] / len, axis_dir[1] / len, axis_dir[2] / len]
        };
        let (u, w) = super::super::perpendicular_basis(&d);
        let tan_a = half_angle.tan();
        let mut pts = Vec::new();
        for ring in 0..n_rings {
            let t = if n_rings == 1 { 0.5 } else { ring as f64 / (n_rings - 1) as f64 };
            let h = h_min + t * (h_max - h_min);
            let r = h * tan_a;
            for j in 0..n_pts_per_ring {
                let theta = 2.0 * PI * j as f64 / n_pts_per_ring as f64;
                let (sin_t, cos_t) = theta.sin_cos();
                pts.push(MeshVertex::from_xyz(
                    apex[0] + h * d[0] + r * (cos_t * u[0] + sin_t * w[0]),
                    apex[1] + h * d[1] + r * (cos_t * u[1] + sin_t * w[1]),
                    apex[2] + h * d[2] + r * (cos_t * u[2] + sin_t * w[2]),
                ));
            }
        }
        pts
    }

    /// Generate face normals for a cone (radially outward, tilted by half_angle).
    fn cone_normals(axis_dir: [f64; 3], half_angle: f64, n_faces: usize) -> Vec<[f64; 3]> {
        let d = {
            let len = (axis_dir[0].powi(2) + axis_dir[1].powi(2) + axis_dir[2].powi(2)).sqrt();
            [axis_dir[0] / len, axis_dir[1] / len, axis_dir[2] / len]
        };
        let (u, w) = super::super::perpendicular_basis(&d);
        let (sin_a, cos_a) = half_angle.sin_cos();
        (0..n_faces)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n_faces as f64;
                let (sin_t, cos_t) = theta.sin_cos();
                let rad = [cos_t * u[0] + sin_t * w[0],
                           cos_t * u[1] + sin_t * w[1],
                           cos_t * u[2] + sin_t * w[2]];
                // Outward normal on cone: cos(α) * radial - sin(α) * axis
                [cos_a * rad[0] - sin_a * d[0],
                 cos_a * rad[1] - sin_a * d[1],
                 cos_a * rad[2] - sin_a * d[2]]
            })
            .collect()
    }

    /// Build a synthetic cone mesh.
    fn build_cone_mesh(
        apex: [f64; 3],
        axis_dir: [f64; 3],
        half_angle: f64,
        n_rings: usize,
        n_pts_per_ring: usize,
        h_min: f64,
        h_max: f64,
    ) -> (Vec<MeshFace>, Vec<MeshVertex>) {
        let vertices = cone_points(apex, axis_dir, half_angle, n_rings, n_pts_per_ring, h_min, h_max);
        let normals = cone_normals(axis_dir, half_angle, n_rings * n_pts_per_ring);
        let mut faces = Vec::new();
        for ring in 0..n_rings.saturating_sub(1) {
            for j in 0..n_pts_per_ring {
                let j_next = (j + 1) % n_pts_per_ring;
                let i0 = ring * n_pts_per_ring + j;
                let i1 = ring * n_pts_per_ring + j_next;
                let i2 = (ring + 1) * n_pts_per_ring + j;
                let fidx = faces.len();
                faces.push(MeshFace {
                    vertex_count: 3,
                    vertex_indices: [i0, i1, i2, 0],
                    neighbors: [-1, -1, -1, -1],
                    normal: Some(normals[fidx % normals.len()]),
                    planar_hypothesis: -1,
                    spherical_hypothesis: -1,
                    conical_hypothesis: -2,
                    toroidal_hypothesis: -2,
                });
            }
        }
        (faces, vertices)
    }

    // ---- Distance tests ----

    #[test]
    fn distance_to_cone_on_surface() {
        // 45° cone: at (1,0,1), r=1, h=1, so on surface.
        let v = MeshVertex::from_xyz(1.0, 0.0, 1.0);
        let d = vertex_to_cone_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], PI / 4.0);
        assert!(d.abs() < 1e-5, "on-surface distance should be 0, got {d}");
    }

    #[test]
    fn distance_to_cone_on_axis() {
        // 45° cone, point on axis at (0,0,5): r=0, h=5
        // d = 0*cos(45) - 5*sin(45) = -5/√2
        let v = MeshVertex::from_xyz(0.0, 0.0, 5.0);
        let d = vertex_to_cone_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], PI / 4.0);
        let expected = -5.0 / 2.0_f64.sqrt();
        assert!((d - expected).abs() < 1e-4, "on-axis distance should be {expected}, got {d}");
    }

    #[test]
    fn distance_to_cone_frustum_on_surface() {
        // Frustum: apex at (0,0,60), axis (0,0,-1), half_angle = atan(1/6)
        // Point (10,0,0): h from apex = 60 along (0,0,-1), r = 10
        let ha = (1.0_f64 / 6.0).atan();
        let v = MeshVertex::from_xyz(10.0, 0.0, 0.0);
        let d = vertex_to_cone_distance(&v, &[0.0, 0.0, 60.0], &[0.0, 0.0, -1.0], ha);
        assert!(d.abs() < 1e-4, "frustum on-surface distance should be 0, got {d}");
    }

    #[test]
    fn distance_to_cone_off_surface() {
        // Same frustum, point (11,0,0): r=11, h=60
        // d = 11*cos(a) - 60*sin(a) where a=atan(1/6)
        // cos(atan(1/6)) = 6/√37, sin(atan(1/6)) = 1/√37
        // d = 66/√37 - 60/√37 = 6/√37 ≈ 0.986
        let ha = (1.0_f64 / 6.0).atan();
        let v = MeshVertex::from_xyz(11.0, 0.0, 0.0);
        let d = vertex_to_cone_distance(&v, &[0.0, 0.0, 60.0], &[0.0, 0.0, -1.0], ha);
        let expected = 6.0 / 37.0_f64.sqrt();
        assert!((d - expected).abs() < 1e-4, "off-surface distance should be {expected}, got {d}");
    }

    // ---- Fit tests ----

    #[test]
    fn cone_fit_z_axis() {
        // Frustum: r_bottom=10 at z=0, r_top=5 at z=30.
        // Apex at (0,0,60), axis (0,0,-1), half_angle = atan(1/6).
        let ha = (1.0_f64 / 6.0).atan();
        let apex = [0.0, 0.0, 60.0];
        let axis = [0.0, 0.0, -1.0];
        let (faces, vertices) = build_cone_mesh(apex, axis, ha, 8, 20, 30.0, 60.0);
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: HashSet<usize> = (0..vertices.len()).collect();

        let result = fit_cone(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "cone fit should succeed");
        let (fit_apex, fit_dir, fit_ha) = result.unwrap();

        assert!((fit_ha - ha).abs() < 0.01, "half_angle should be ~{:.3}, got {fit_ha:.3}", ha);
        let axis_z = fit_dir[2].abs();
        assert!((axis_z - 1.0).abs() < 0.02, "axis should be along Z, got z={axis_z}");
        assert!(fit_apex[0].abs() < 0.5, "apex.x should be near 0, got {}", fit_apex[0]);
        assert!(fit_apex[1].abs() < 0.5, "apex.y should be near 0, got {}", fit_apex[1]);
        assert!((fit_apex[2] - 60.0).abs() < 1.0, "apex.z should be near 60, got {}", fit_apex[2]);
    }

    #[test]
    fn cone_fit_conical_hole() {
        // Cone: bottom r=6 at z=0, top r=4 at z=20.
        // half_angle = atan(0.1), apex at (0,0,60), axis (0,0,-1).
        let ha = 0.1_f64.atan();
        let apex = [0.0, 0.0, 60.0];
        let axis = [0.0, 0.0, -1.0];
        // h range: from apex distance 40 (z=20, r=4) to distance 60 (z=0, r=6)
        let (faces, vertices) = build_cone_mesh(apex, axis, ha, 6, 16, 40.0, 60.0);
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: HashSet<usize> = (0..vertices.len()).collect();

        let result = fit_cone(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "conical hole fit should succeed");
        let (fit_apex, fit_dir, fit_ha) = result.unwrap();

        assert!((fit_ha - ha).abs() < 0.01, "half_angle should be ~{:.4}, got {fit_ha:.4}", ha);
        let axis_z = fit_dir[2].abs();
        assert!((axis_z - 1.0).abs() < 0.02, "axis should be along Z, got z={axis_z}");
        assert!(fit_apex[0].abs() < 0.5, "apex.x should be near 0, got {}", fit_apex[0]);
        assert!(fit_apex[1].abs() < 0.5, "apex.y should be near 0, got {}", fit_apex[1]);
        assert!((fit_apex[2] - 60.0).abs() < 1.5, "apex.z should be near 60, got {}", fit_apex[2]);
    }

    #[test]
    fn cone_fit_rejects_degenerate() {
        let faces = vec![MeshFace {
            vertex_count: 3,
            vertex_indices: [0, 1, 2, 0],
            neighbors: [-1, -1, -1, -1],
            normal: Some([0.0, 0.0, 1.0]),
            planar_hypothesis: -1,
            spherical_hypothesis: -1,
            conical_hypothesis: -2,
            toroidal_hypothesis: -2,
        }];
        let vertices = vec![
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0),
        ];
        let vertex_set: HashSet<usize> = (0..2).collect();
        assert!(fit_cone(&[0], &vertex_set, &faces, &vertices).is_none(),
                "cone fit should reject degenerate input");
    }
}

