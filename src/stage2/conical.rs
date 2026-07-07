//! Stage 2.4: Conical hypothesis deduction and comparison.

use std::collections::HashSet;

use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{Const, Dyn, OMatrix, OVector, Owned};

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_CONICAL_HYPOTHESIS};
use crate::viz;

use super::{
    build_normal_covariance, centered_cone_overlay, centered_cylinder_overlay, cross3, dot3, face_centroid,
    largest_eigenvector_3x3, normalize3, perpendicular_basis, smallest_eigenvector_3x3,
    ConicalHypothesis, Stage2CompareError, MIN_CROSS_THRESHOLD, MIN_CONE_FACES,
    MIN_CYLINDER_FACES,
};
use super::cylindrical::{vertex_to_cylinder_distance, fit_cylinder};

/// Minimum half-angle (degrees) — below this the cone is effectively a cylinder.
const MIN_HALF_ANGLE_DEG: f64 = 2.0;
/// Maximum half-angle (degrees) — above this the cone is effectively a plane.
const MAX_HALF_ANGLE_DEG: f64 = 85.0;
/// Maximum apex distance as a multiple of the mesh bounding-box diagonal.
const MAX_APEX_DISTANCE_FACTOR: f64 = 10.0;

// ---------------------------------------------------------------------------
// Conical fitting: 4-step algorithm
//
// Step 1: Estimate apex from least-squares line intersection of normal lines.
// Step 2: Estimate axis as eigenvector of smallest eigenvalue of normal covariance.
// Step 3: Estimate half-angle from median of atan(r_i/h_i).
// Step 4: LM refinement of all parameters.
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

// ---------------------------------------------------------------------------
// Step 1: Apex estimation via least-squares line intersection
// ---------------------------------------------------------------------------

/// Estimate the cone apex as the point minimizing sum of squared distances
/// to lines L_i(t) = p_i + t * n_i (one line per face centroid + normal).
///
/// The optimal point minimizes  Σ ||(I - n_i n_i^T)(a - p_i)||^2
/// Leading to normal equations: (Σ (I - n_i n_i^T)) a = Σ (I - n_i n_i^T) p_i
#[cfg(test)]
fn estimate_apex(
    face_indices: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<[f64; 3]> {
    // Check if there is a vertex shared by ALL faces — it's the apex.
    // This gives an exact result for pointed cones and avoids an ill-conditioned system.
    let n_faces = face_indices.len();
    let mut vertex_face_count = std::collections::HashMap::<usize, usize>::new();
    for &fi in face_indices {
        let vc = faces[fi].vertex_count as usize;
        for &vi in &faces[fi].vertex_indices[..vc] {
            *vertex_face_count.entry(vi).or_insert(0) += 1;
        }
    }
    if n_faces >= 3 {
        for (&vi, &count) in &vertex_face_count {
            if count == n_faces {
                return Some([vertices[vi].x, vertices[vi].y, vertices[vi].z]);
            }
        }
    }

    // General case: least-squares intersection of normal lines.
    // Build 3x3 system: A a = b  where A = Σ (I - n n^T), b = Σ (I - n n^T) p
    let mut a_mat = [0.0_f64; 9]; // row-major 3x3
    let mut b_vec = [0.0_f64; 3];

    for &fi in face_indices {
        let n = match faces[fi].normal {
            Some(n) => n,
            None => continue,
        };
        let c = face_centroid(&faces[fi], vertices);

        // (I - n n^T) is the projection matrix perpendicular to n
        for row in 0..3 {
            for col in 0..3 {
                let ident = if row == col { 1.0 } else { 0.0 };
                let proj = ident - n[row] * n[col];
                a_mat[row * 3 + col] += proj;
                b_vec[row] += proj * c[col];
            }
        }
    }

    // Solve 3x3 system via Gaussian elimination with partial pivoting
    solve_3x3(&a_mat, &b_vec)
}

/// Solve a 3x3 linear system Ax = b. Returns None if singular.
#[cfg(test)]
fn solve_3x3(a: &[f64; 9], b: &[f64; 3]) -> Option<[f64; 3]> {
    let mut aug = [[0.0_f64; 4]; 3];
    for i in 0..3 {
        for j in 0..3 {
            aug[i][j] = a[i * 3 + j];
        }
        aug[i][3] = b[i];
    }

    for col in 0..3 {
        // Partial pivoting
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..3 {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            return None;
        }
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        for row in (col + 1)..3 {
            let factor = aug[row][col] / pivot;
            for j in col..4 {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    // Back-substitution
    let mut x = [0.0_f64; 3];
    for i in (0..3).rev() {
        x[i] = aug[i][3];
        for j in (i + 1)..3 {
            x[i] -= aug[i][j] * x[j];
        }
        if aug[i][i].abs() < 1e-30 {
            return None;
        }
        x[i] /= aug[i][i];
    }
    Some(x)
}

// ---------------------------------------------------------------------------
// Step 2: Axis estimation from normal covariance
// ---------------------------------------------------------------------------

/// Estimate axis direction as the eigenvector with smallest eigenvalue of
/// the normal covariance matrix. On a cone, normals have constant dot product
/// with the axis, so the variance of n·u is minimized along the true axis.
/// Works best with many faces spanning a full azimuthal range.
fn estimate_axis_from_normals(
    face_indices: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<[f64; 3]> {
    let cov = build_normal_covariance(face_indices, faces, vertices);
    let mut axis = smallest_eigenvector_3x3(&cov);
    let len = normalize3(&mut axis);
    if len > 1e-15 {
        Some(axis)
    } else {
        None
    }
}

/// Multi-strategy axis estimation for truncated cones.
///
/// Strategy 1 (≥8 faces): covariance smallest-eigenvector.
/// Strategy 2 (few faces, shared vertex): apex→base-centroid direction.
/// Strategy 3 (few faces, no shared vertex): cross product of normal diffs.
fn compute_truncated_cone_axis(
    face_indices: &[usize],
    vertex_set: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<[f64; 3]> {
    // With many faces, covariance is robust.
    if face_indices.len() >= 8 {
        return estimate_axis_from_normals(face_indices, faces, vertices);
    }

    // Strategy 2: vertex shared by most faces → axis from that vertex to centroid.
    let mut vertex_face_count = std::collections::HashMap::<usize, usize>::new();
    for &fi in face_indices {
        let vc = faces[fi].vertex_count as usize;
        for &vi in &faces[fi].vertex_indices[..vc] {
            *vertex_face_count.entry(vi).or_insert(0) += 1;
        }
    }
    let threshold = face_indices.len().max(2) - 1;
    let mut apex_candidates: Vec<usize> = vertex_face_count.iter()
        .filter(|(_, &c)| c >= threshold)
        .map(|(&vi, _)| vi)
        .collect();
    apex_candidates.sort(); // deterministic

    for &apex_vi in &apex_candidates {
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
    estimate_axis_from_normal_diffs(face_indices, faces)
}

/// Estimate axis from cross product of normal differences.
///
/// On a cone, n_i · axis = sin(ha) = constant for all face normals.
/// Therefore (n_i - n_j) · axis = 0 — all normal differences are perpendicular
/// to the axis. The cross product of two independent differences gives the
/// exact axis direction. Works with as few as 3 faces.
fn estimate_axis_from_normal_diffs(
    face_indices: &[usize],
    faces: &[MeshFace],
) -> Option<[f64; 3]> {
    let normals: Vec<[f64; 3]> = face_indices.iter()
        .filter_map(|&fi| faces[fi].normal)
        .collect();
    if normals.len() < 3 {
        return None;
    }

    // Compute differences from first normal
    let diffs: Vec<[f64; 3]> = normals[1..].iter()
        .map(|n| [n[0] - normals[0][0], n[1] - normals[0][1], n[2] - normals[0][2]])
        .collect();

    // Find the pair of differences whose cross product is largest (most independent)
    let mut best_cross_sq = 0.0_f64;
    let mut best_axis = [0.0_f64; 3];
    for i in 0..diffs.len() {
        for j in (i + 1)..diffs.len() {
            let c = cross3(&diffs[i], &diffs[j]);
            let csq = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
            if csq > best_cross_sq {
                best_cross_sq = csq;
                best_axis = c;
            }
        }
    }

    let len = normalize3(&mut best_axis);
    if len > 1e-15 {
        Some(best_axis)
    } else {
        None
    }
}

/// Estimate axis from vertex directions relative to a known apex.
/// The axis is the largest eigenvector of Σ d̂_i⊗d̂_i where d̂_i = (v_i - apex)/|v_i - apex|.
/// This works because all vertex directions cluster around the axis for narrow cones.
fn estimate_axis_from_apex(
    apex: &[f64; 3],
    vertex_set: &[usize],
    vertices: &[MeshVertex],
) -> Option<[f64; 3]> {
    let mut cov = [0.0_f64; 6]; // upper triangle: [xx, xy, xz, yy, yz, zz]
    let mut count = 0;
    for &vi in vertex_set {
        let v = &vertices[vi];
        let mut d = [v.x - apex[0], v.y - apex[1], v.z - apex[2]];
        let len = normalize3(&mut d);
        if len < 1e-10 {
            continue; // skip the apex vertex itself
        }
        cov[0] += d[0] * d[0];
        cov[1] += d[0] * d[1];
        cov[2] += d[0] * d[2];
        cov[3] += d[1] * d[1];
        cov[4] += d[1] * d[2];
        cov[5] += d[2] * d[2];
        count += 1;
    }
    if count < 2 {
        return None;
    }
    let mut axis = largest_eigenvector_3x3(&cov);
    let len = normalize3(&mut axis);
    if len > 1e-15 {
        Some(axis)
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Step 2+3 alternative: profile-based apex and half-angle estimation
// ---------------------------------------------------------------------------

/// Estimate apex and half-angle by projecting vertices to (h, r) along the axis
/// and fitting a linear profile r = m*h + b.
///
/// For a cone: r_i = tan(θ) * (h_i - h_apex), so the linear regression gives:
/// - slope m = tan(θ), i.e. half_angle = atan(|m|)
/// - apex offset along axis = -b/m from the centroid
///
/// This approach is more robust than normal-line intersection for small
/// half-angles, since it doesn't require nearly-parallel normal lines to
/// intersect reliably.
///
/// Also orients the axis so that h increases away from the apex.
fn estimate_apex_and_angle_from_profile(
    axis: &mut [f64; 3],
    vertex_set: &[usize],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], f64)> {
    let n = vertex_set.len();
    if n < 3 {
        return None;
    }

    // Compute centroid of vertex set
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

    // Project vertices into (h, r) space relative to centroid
    let mut h_vals = Vec::with_capacity(n);
    let mut r_vals = Vec::with_capacity(n);
    for &vi in vertex_set {
        let v = &vertices[vi];
        let dx = v.x - cx;
        let dy = v.y - cy;
        let dz = v.z - cz;
        let h = dx * axis[0] + dy * axis[1] + dz * axis[2];
        let rx = dx - h * axis[0];
        let ry = dy - h * axis[1];
        let rz = dz - h * axis[2];
        let r = (rx * rx + ry * ry + rz * rz).sqrt();
        h_vals.push(h);
        r_vals.push(r);
    }

    // Linear regression: r = m*h + b
    let h_mean: f64 = h_vals.iter().sum::<f64>() * inv_n;
    let r_mean: f64 = r_vals.iter().sum::<f64>() * inv_n;
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for i in 0..n {
        let dh = h_vals[i] - h_mean;
        num += dh * (r_vals[i] - r_mean);
        den += dh * dh;
    }
    if den < 1e-30 {
        return None;
    }
    let m = num / den;
    let b = r_mean - m * h_mean;

    // half_angle = atan(|m|)
    let half_angle = m.abs().atan();
    if half_angle < 0.01 || half_angle > 85.0_f64.to_radians() {
        return None;
    }

    // Apex = centroid + apex_offset * axis, where apex_offset = -b/m
    if m.abs() < 1e-15 {
        return None;
    }
    let apex_offset = -b / m;
    let apex = [
        cx + apex_offset * axis[0],
        cy + apex_offset * axis[1],
        cz + apex_offset * axis[2],
    ];

    // Orient axis so it points away from the apex (increasing h from apex)
    // If m > 0, r increases with h, meaning h increases away from apex → axis is correct
    // If m < 0, r decreases with h, meaning apex is on the h+ side → flip axis
    if m < 0.0 {
        axis[0] = -axis[0];
        axis[1] = -axis[1];
        axis[2] = -axis[2];
    }

    Some((apex, half_angle))
}

// ---------------------------------------------------------------------------
// Step 3: Half-angle estimation (for pointed cone path)
// ---------------------------------------------------------------------------

/// Estimate half-angle from apex and axis using median of atan(r_i / h_i).
/// Also orients the axis so that h_i > 0 for the majority of points.
fn estimate_half_angle(
    apex: &[f64; 3],
    axis: &mut [f64; 3],
    vertex_set: &[usize],
    vertices: &[MeshVertex],
) -> Option<f64> {
    // Orient axis so most points have positive h (away from apex)
    let mut h_sum = 0.0_f64;
    for &vi in vertex_set {
        let v = &vertices[vi];
        h_sum += (v.x - apex[0]) * axis[0]
               + (v.y - apex[1]) * axis[1]
               + (v.z - apex[2]) * axis[2];
    }
    if h_sum < 0.0 {
        axis[0] = -axis[0];
        axis[1] = -axis[1];
        axis[2] = -axis[2];
    }

    // Compute r_i / h_i for each vertex
    let mut ratios: Vec<f64> = Vec::with_capacity(vertex_set.len());
    for &vi in vertex_set {
        let v = &vertices[vi];
        let dx = v.x - apex[0];
        let dy = v.y - apex[1];
        let dz = v.z - apex[2];
        let h = dx * axis[0] + dy * axis[1] + dz * axis[2];
        let rx = dx - h * axis[0];
        let ry = dy - h * axis[1];
        let rz = dz - h * axis[2];
        let r = (rx * rx + ry * ry + rz * rz).sqrt();
        // Skip points very close to the apex (h ≈ 0) — they don't constrain θ
        if h.abs() > 1e-10 {
            ratios.push(r / h);
        }
    }

    if ratios.is_empty() {
        return None;
    }

    // Median for robustness
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if ratios.len() % 2 == 1 {
        ratios[ratios.len() / 2]
    } else {
        (ratios[ratios.len() / 2 - 1] + ratios[ratios.len() / 2]) / 2.0
    };

    if median <= 0.0 {
        return None;
    }

    Some(median.atan())
}

// ---------------------------------------------------------------------------
// Step 4: LM refinement
// ---------------------------------------------------------------------------

/// LM problem for cone fitting with fixed apex.
/// Parameters: [alpha, beta, half_angle]
/// Only optimizes axis direction and half-angle; apex is fixed.
/// Used when the apex is known exactly (e.g., shared vertex in pointed cone).
struct ConeFixedApexLMProblem {
    points: Vec<[f64; 3]>,
    apex: [f64; 3],
    a0: [f64; 3],
    u0: [f64; 3],
    w0: [f64; 3],
    params: OVector<f64, Const<3>>,
}

impl ConeFixedApexLMProblem {
    fn new(
        points: Vec<[f64; 3]>,
        apex: [f64; 3],
        initial_axis: [f64; 3],
        initial_half_angle: f64,
    ) -> Self {
        let (u0, w0) = perpendicular_basis(&initial_axis);
        let params = OVector::<f64, Const<3>>::new(0.0, 0.0, initial_half_angle);
        Self { points, apex, a0: initial_axis, u0, w0, params }
    }

    fn axis_dir_from_params(&self, params: &OVector<f64, Const<3>>) -> [f64; 3] {
        let (alpha, beta) = (params[0], params[1]);
        let mut v = [
            self.a0[0] + alpha * self.u0[0] + beta * self.w0[0],
            self.a0[1] + alpha * self.u0[1] + beta * self.w0[1],
            self.a0[2] + alpha * self.u0[2] + beta * self.w0[2],
        ];
        normalize3(&mut v);
        v
    }
}

impl LeastSquaresProblem<f64, Dyn, Const<3>> for ConeFixedApexLMProblem {
    type ParameterStorage = Owned<f64, Const<3>>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Const<3>>;

    fn set_params(&mut self, x: &OVector<f64, Const<3>>) {
        self.params.copy_from(x);
    }

    fn params(&self) -> OVector<f64, Const<3>> {
        self.params
    }

    fn residuals(&self) -> Option<OVector<f64, Dyn>> {
        let a = self.axis_dir_from_params(&self.params);
        let half_angle = self.params[2];
        let (sin_a, cos_a) = half_angle.sin_cos();
        let n = self.points.len();
        let mut residuals = OVector::<f64, Dyn>::zeros_generic(Dyn(n), Const::<1>);
        for (i, p) in self.points.iter().enumerate() {
            let d = [p[0] - self.apex[0], p[1] - self.apex[1], p[2] - self.apex[2]];
            let h = dot3(&d, &a);
            let radial = [d[0] - h * a[0], d[1] - h * a[1], d[2] - h * a[2]];
            let r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
            residuals[i] = r * cos_a - h * sin_a;
        }
        Some(residuals)
    }

    fn jacobian(&self) -> Option<OMatrix<f64, Dyn, Const<3>>> {
        let n = self.points.len();
        let mut jac = OMatrix::<f64, Dyn, Const<3>>::zeros_generic(Dyn(n), Const::<3>);
        let a = self.axis_dir_from_params(&self.params);
        let half_angle = self.params[2];
        let (sa, ca) = half_angle.sin_cos();

        let alpha = self.params[0];
        let beta = self.params[1];
        let vx = self.a0[0] + alpha * self.u0[0] + beta * self.w0[0];
        let vy = self.a0[1] + alpha * self.u0[1] + beta * self.w0[1];
        let vz = self.a0[2] + alpha * self.u0[2] + beta * self.w0[2];
        let big_l = (vx * vx + vy * vy + vz * vz).sqrt().max(1e-15);

        for (i, p) in self.points.iter().enumerate() {
            let d = [p[0] - self.apex[0], p[1] - self.apex[1], p[2] - self.apex[2]];
            let h = dot3(&d, &a);
            let lat = [d[0] - h * a[0], d[1] - h * a[1], d[2] - h * a[2]];
            let r = (lat[0] * lat[0] + lat[1] * lat[1] + lat[2] * lat[2]).sqrt().max(1e-15);
            let lat_u = dot3(&lat, &self.u0);
            let lat_w = dot3(&lat, &self.w0);
            let axis_sens = ca * h / r + sa;

            jac[(i, 0)] = -lat_u * axis_sens / big_l; // ∂f/∂α
            jac[(i, 1)] = -lat_w * axis_sens / big_l; // ∂f/∂β
            jac[(i, 2)] = -r * sa - h * ca;           // ∂f/∂half_angle
        }
        Some(jac)
    }
}

/// LM problem for cone fitting (full 6-parameter).
/// Parameters: [alpha, beta, apex_x, apex_y, apex_z, half_angle]
/// alpha, beta: tilt of axis from initial direction in perpendicular directions.
/// Residual: r_i * cos(θ) - h_i * sin(θ)  (signed distance to cone surface).
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

/// Fit cone parameters (axis direction, apex, half-angle) from a set of faces.
/// Returns (apex, axis_direction, half_angle).
///
/// Two estimation paths depending on cone topology:
///
/// **Pointed cone** (a vertex shared by all faces = apex):
/// 1. Apex = shared vertex (exact).
/// 2. Axis = cross product of normal differences (exact for ideal cones).
/// 3. Half-angle from median of atan(r_i / h_i).
/// 4. Fixed-apex LM refinement (3 params: axis tilt + half_angle).
///
/// **Truncated cone** (no shared vertex):
/// 1. Axis = smallest eigenvector of normal covariance.
/// 2. Apex from normal-line intersection.
/// 3. Half-angle from median of atan(r_i / h_i).
/// 4. Full LM refinement (6 params).
fn fit_cone(
    face_indices: &[usize],
    vertex_set: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64)> {
    if vertex_set.len() < 4 {
        return None;
    }
    let max_err = |fit: &([f64; 3], [f64; 3], f64)| -> f64 {
        vertex_set
            .iter()
            .map(|&vi| vertex_to_cone_distance(&vertices[vi], &fit.0, &fit.1, fit.2).abs())
            .fold(0.0_f64, f64::max)
    };

    let primary = fit_cone_primary(face_indices, vertex_set, faces, vertices);
    let mut best: Option<(f64, ([f64; 3], [f64; 3], f64))> = None;
    if let Some(fit) = primary {
        best = Some((max_err(&fit), fit));
    }

    // Normal-line axis rescue for narrow azimuthal arcs: the primary axis
    // estimates (normal covariance, normal differences) degenerate when the
    // seed's normals are nearly parallel (e.g. side cone patches of a
    // transition piece spanning only a few tessellation columns). Cone
    // surface normals all pass through the axis, so closest-approach
    // midpoints of normal-line pairs trace the axis directly.
    const RESCUE_ERR_THRESHOLD: f64 = 1e-4; // mm
    if best.as_ref().is_none_or(|(e, _)| *e > RESCUE_ERR_THRESHOLD) {
        if let Some(fit) = fit_cone_normal_line_rescue(face_indices, vertex_set, faces, vertices) {
            let err = max_err(&fit);
            if best.as_ref().is_none_or(|(e, _)| err < *e) {
                best = Some((err, fit));
            }
        }
    }

    best.map(|(_, f)| f)
}

/// Closest-approach midpoint between two 3D lines (p1 + t·d1, p2 + s·d2).
/// Returns None for near-parallel lines or when the gap between the closest
/// points is large relative to their distance from the line origins.
fn line_closest_midpoint(
    p1: &[f64; 3],
    d1: &[f64; 3],
    p2: &[f64; 3],
    d2: &[f64; 3],
) -> Option<[f64; 3]> {
    let r = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let a = dot3(d1, d1);
    let b = dot3(d1, d2);
    let c = dot3(d2, d2);
    let d = dot3(d1, &r);
    let e = dot3(d2, &r);
    let denom = a * c - b * b;
    if denom.abs() < 1e-6 {
        return None; // near-parallel
    }
    let t = (d * c - e * b) / denom;
    let sv = (d * b - e * a) / denom;
    let q1 = [p1[0] + t * d1[0], p1[1] + t * d1[1], p1[2] + t * d1[2]];
    let q2 = [p2[0] + sv * d2[0], p2[1] + sv * d2[1], p2[2] + sv * d2[2]];
    let gap = ((q1[0] - q2[0]).powi(2) + (q1[1] - q2[1]).powi(2) + (q1[2] - q2[2]).powi(2)).sqrt();
    let scale = t.abs().max(sv.abs()).max(1e-9);
    if gap > scale * 0.5 {
        return None; // lines don't come close: not a surface of revolution
    }
    Some([
        (q1[0] + q2[0]) / 2.0,
        (q1[1] + q2[1]) / 2.0,
        (q1[2] + q2[2]) / 2.0,
    ])
}

/// Fit a cone by estimating the axis from normal-line closest-approach
/// midpoints, then the profile line in (h, r) coordinates, then LM.
fn fit_cone_normal_line_rescue(
    face_indices: &[usize],
    vertex_set: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64)> {
    // Axis points from normal-line pairs.
    let mut axis_pts: Vec<[f64; 3]> = Vec::new();
    let n_faces = face_indices.len();
    for i in 0..n_faces {
        let fi = face_indices[i];
        let Some(n1) = faces[fi].normal else { continue };
        let c1 = face_centroid(&faces[fi], vertices);
        for &fj in face_indices.iter().skip(i + 1) {
            let Some(n2) = faces[fj].normal else { continue };
            let cr = cross3(&n1, &n2);
            let cr_len = (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt();
            if cr_len < 0.01 {
                continue; // near-parallel normals give no axis information
            }
            let c2 = face_centroid(&faces[fj], vertices);
            if let Some(m) = line_closest_midpoint(&c1, &n1, &c2, &n2) {
                axis_pts.push(m);
            }
            if axis_pts.len() >= 64 {
                break;
            }
        }
        if axis_pts.len() >= 64 {
            break;
        }
    }
    if axis_pts.len() < 2 {
        return None;
    }

    // 3D line through the axis points: centroid + dominant scatter direction.
    let n = axis_pts.len() as f64;
    let mut origin = [0.0_f64; 3];
    for p in &axis_pts {
        origin[0] += p[0];
        origin[1] += p[1];
        origin[2] += p[2];
    }
    origin[0] /= n;
    origin[1] /= n;
    origin[2] /= n;
    let mut cov = [0.0_f64; 6];
    for p in &axis_pts {
        let dx = p[0] - origin[0];
        let dy = p[1] - origin[1];
        let dz = p[2] - origin[2];
        cov[0] += dx * dx;
        cov[1] += dx * dy;
        cov[2] += dx * dz;
        cov[3] += dy * dy;
        cov[4] += dy * dz;
        cov[5] += dz * dz;
    }
    let mut axis = largest_eigenvector_3x3(&cov);
    if normalize3(&mut axis) < 1e-12 {
        return None;
    }

    // Profile line in (h, r): r = slope·h + b  →  half_angle = atan(|slope|).
    let mut sum_h = 0.0_f64;
    let mut sum_r = 0.0_f64;
    let mut sum_hh = 0.0_f64;
    let mut sum_hr = 0.0_f64;
    let np = vertex_set.len() as f64;
    let mut hr_pts: Vec<(f64, f64)> = Vec::with_capacity(vertex_set.len());
    for &vi in vertex_set {
        let v = &vertices[vi];
        let d = [v.x - origin[0], v.y - origin[1], v.z - origin[2]];
        let h = dot3(&d, &axis);
        let radial = [d[0] - h * axis[0], d[1] - h * axis[1], d[2] - h * axis[2]];
        let r = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
        sum_h += h;
        sum_r += r;
        sum_hh += h * h;
        sum_hr += h * r;
        hr_pts.push((h, r));
    }
    let denom = np * sum_hh - sum_h * sum_h;
    if denom.abs() < 1e-12 {
        return None;
    }
    let mut slope = (np * sum_hr - sum_h * sum_r) / denom;
    if slope.abs() < 1e-6 {
        return None; // cylinder, not a cone
    }
    // Validity condition of this rescue: the (h, r) profile must actually be
    // LINEAR. Doubly-curved patches (torus corner blends) have an arc-shaped
    // profile that a line fits poorly; without this gate the LM polish can
    // still land within tolerance for tiny 2-ring subsets and commit spurious
    // cones that block torus seeding.
    {
        let b0 = (sum_r - slope * sum_h) / np;
        let mut ss = 0.0_f64;
        for &(h, r) in &hr_pts {
            let d = r - (slope * h + b0);
            ss += d * d;
        }
        let rms = (ss / np).sqrt();
        let r_mean = (sum_r / np).abs().max(1e-9);
        if rms > r_mean * 0.02 {
            return None;
        }
    }
    // Orient the axis so h increases from apex toward the base (r grows):
    // under h' = -h the slope negates and the h sums flip sign.
    let mut h_mean = sum_h / np;
    if slope < 0.0 {
        axis = [-axis[0], -axis[1], -axis[2]];
        slope = -slope;
        h_mean = -h_mean;
    }
    let r_mean = sum_r / np;
    let intercept = r_mean - slope * h_mean;

    let half_angle = slope.atan();
    if !(2.0_f64.to_radians()..=85.0_f64.to_radians()).contains(&half_angle) {
        return None;
    }
    let apex_h = -intercept / slope;
    let apex = [
        origin[0] + apex_h * axis[0],
        origin[1] + apex_h * axis[1],
        origin[2] + apex_h * axis[2],
    ];

    // LM polish with recentered coordinates.
    let n_pts = vertex_set.len() as f64;
    let mut centroid = [0.0_f64; 3];
    for &vi in vertex_set {
        centroid[0] += vertices[vi].x;
        centroid[1] += vertices[vi].y;
        centroid[2] += vertices[vi].z;
    }
    centroid[0] /= n_pts;
    centroid[1] /= n_pts;
    centroid[2] /= n_pts;
    let centered_points: Vec<[f64; 3]> = vertex_set
        .iter()
        .map(|&vi| [
            vertices[vi].x - centroid[0],
            vertices[vi].y - centroid[1],
            vertices[vi].z - centroid[2],
        ])
        .collect();
    let centered_apex = [
        apex[0] - centroid[0],
        apex[1] - centroid[1],
        apex[2] - centroid[2],
    ];
    let problem = ConeLMProblem::new(centered_points, axis, centered_apex, half_angle);
    let (result, _report) = LevenbergMarquardt::new().with_patience(150).minimize(problem);
    let refined_dir = result.axis_dir_from_params(&result.params);
    let refined_apex = [
        result.params[2] + centroid[0],
        result.params[3] + centroid[1],
        result.params[4] + centroid[2],
    ];
    let refined_ha = result.params[5];
    if refined_ha > 0.01 && refined_ha < 85.0_f64.to_radians() {
        Some((refined_apex, refined_dir, refined_ha))
    } else {
        Some((apex, axis, half_angle))
    }
}

fn fit_cone_primary(
    face_indices: &[usize],
    vertex_set: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64)> {
    if vertex_set.len() < 4 {
        return None;
    }

    // Check for pointed cone: a vertex present in ALL seed faces.
    let shared_apex = find_shared_apex_vertex(face_indices, faces, vertices);

    // Recenter vertices around their centroid for numerical stability.
    let n_pts = vertex_set.len() as f64;
    let mut centroid = [0.0_f64; 3];
    for &vi in vertex_set.iter() {
        centroid[0] += vertices[vi].x;
        centroid[1] += vertices[vi].y;
        centroid[2] += vertices[vi].z;
    }
    centroid[0] /= n_pts;
    centroid[1] /= n_pts;
    centroid[2] /= n_pts;

    let centered_points: Vec<[f64; 3]> = vertex_set.iter()
        .map(|&vi| [
            vertices[vi].x - centroid[0],
            vertices[vi].y - centroid[1],
            vertices[vi].z - centroid[2],
        ])
        .collect();

    if let Some(apex) = shared_apex {
        // --- Pointed cone path ---
        // For fixed-apex LM, recenter both points and apex.
        let centered_apex = [
            apex[0] - centroid[0],
            apex[1] - centroid[1],
            apex[2] - centroid[2],
        ];
        let mut axis = estimate_axis_from_normal_diffs(face_indices, faces)
            .or_else(|| estimate_axis_from_apex(&apex, vertex_set, vertices))?;
        let half_angle = estimate_half_angle(&apex, &mut axis, vertex_set, vertices)?;

        let problem = ConeFixedApexLMProblem::new(centered_points, centered_apex, axis, half_angle);
        let (result, _report) = LevenbergMarquardt::new().with_patience(150).minimize(problem);
        let refined_dir = result.axis_dir_from_params(&result.params);
        let refined_ha = result.params[2];
        if refined_ha > 0.01 && refined_ha < 85.0_f64.to_radians() {
            Some((apex, refined_dir, refined_ha))
        } else {
            Some((apex, axis, half_angle))
        }
    } else {
        // --- Truncated cone path ---
        // Multi-strategy axis: covariance for ≥8 faces, vertex-based for seeds
        let mut axis = compute_truncated_cone_axis(
            face_indices, vertex_set, faces, vertices)?;

        // Step 2+3: Apex and half-angle from linear regression of vertex projections.
        let (apex, half_angle) = estimate_apex_and_angle_from_profile(
            &mut axis, vertex_set, vertices)?;

        // Full LM refinement with recentered coordinates
        let centered_apex = [
            apex[0] - centroid[0],
            apex[1] - centroid[1],
            apex[2] - centroid[2],
        ];
        let problem = ConeLMProblem::new(centered_points, axis, centered_apex, half_angle);
        let (result, _report) = LevenbergMarquardt::new().with_patience(150).minimize(problem);
        let refined_dir = result.axis_dir_from_params(&result.params);
        let refined_apex = [
            result.params[2] + centroid[0],
            result.params[3] + centroid[1],
            result.params[4] + centroid[2],
        ];
        let refined_ha = result.params[5];
        let (final_apex, final_dir, final_ha) =
            if refined_ha > 0.01 && refined_ha < 85.0_f64.to_radians() {
                (refined_apex, refined_dir, refined_ha)
            } else {
                (apex, axis, half_angle)
            };

        /* // DEBUG: normal consistency variance
        let dots: Vec<f64> = face_indices.iter()
            .filter_map(|&fi| faces[fi].normal.map(|n| {
                n[0] * final_dir[0] + n[1] * final_dir[1] + n[2] * final_dir[2]
            }))
            .collect();
        if dots.len() >= 2 {
            let mean = dots.iter().sum::<f64>() / dots.len() as f64;
            let var = dots.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / dots.len() as f64;
            eprintln!("[DEBUG fit_cone truncated] n_faces={} half_angle={:.2}° n·axis var={:.6} std={:.4}",
                face_indices.len(), final_ha.to_degrees(), var, var.sqrt());
        } */

        Some((final_apex, final_dir, final_ha))
    }
}

/// Find a vertex that appears in all faces (the apex of a pointed cone).
fn find_shared_apex_vertex(
    face_indices: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<[f64; 3]> {
    let n_faces = face_indices.len();
    if n_faces < 3 {
        return None;
    }
    let mut vertex_face_count = std::collections::HashMap::<usize, usize>::new();
    for &fi in face_indices {
        let vc = faces[fi].vertex_count as usize;
        for &vi in &faces[fi].vertex_indices[..vc] {
            *vertex_face_count.entry(vi).or_insert(0) += 1;
        }
    }
    for (&vi, &count) in &vertex_face_count {
        if count == n_faces {
            let apex = [vertices[vi].x, vertices[vi].y, vertices[vi].z];
            // Verify that non-apex vertices are at roughly equal distances from the apex,
            // as expected for a cone base. Without this check, any 3 triangles sharing
            // a vertex (common in regular triangle meshes) falsely trigger the pointed-cone path.
            let mut dists: Vec<f64> = Vec::new();
            for &fi in face_indices {
                let vc = faces[fi].vertex_count as usize;
                for &vj in &faces[fi].vertex_indices[..vc] {
                    if vj != vi {
                        let dx = vertices[vj].x - apex[0];
                        let dy = vertices[vj].y - apex[1];
                        let dz = vertices[vj].z - apex[2];
                        dists.push((dx * dx + dy * dy + dz * dz).sqrt());
                    }
                }
            }
            if dists.len() < 2 {
                return Some(apex);
            }
            let mean_dist = dists.iter().sum::<f64>() / dists.len() as f64;
            if mean_dist < 1e-12 {
                return Some(apex);
            }
            // Coefficient of variation: std_dev / mean. For a cone base, distances
            // should be nearly equal (CV < 0.15). For random mesh topology, they vary widely.
            let var = dists.iter().map(|d| (d - mean_dist).powi(2)).sum::<f64>() / dists.len() as f64;
            let cv = var.sqrt() / mean_dist;
            if cv < 0.15 {
                return Some(apex);
            }
        }
    }
    None
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
        let projected_points: Vec<[f64; 3]> = hyp.faces.iter().map(|&fi| {
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

            if hyp.is_cylinder() {
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
            }
        }).collect();

        let (max_dist, _) = stage1::max_distance_to_shape(&projected_points, compare_shape);
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

    /// Vertices should lie precisely on the conical surface, so
    /// the seed vertex tolerance matches the base vertex tolerance.
    fn seed_vtol_multiplier() -> f64 { 1.0 }

    /// Tessellated cones with N circumferential divisions have inter-strip
    /// dihedral angles ≈ 360°/N, which can exceed the default angular
    /// tolerance (e.g. 20° for N=18 vs the 17.5° default). Vertex-to-surface
    /// distance is the primary discriminator for cones/cylinders, so the
    /// angular check is relaxed 2×.
    fn angular_tol_multiplier() -> f64 { 2.0 }

    /// Two passes: pass 1 excludes junction faces from fan seeds; pass 2
    /// retries remaining undeduced faces without the filter (narrow cone
    /// patches between tangent junctions can be all-transition faces).
    fn seeding_passes() -> usize { 2 }

    fn generate_seeds(mesh: &ConnectedMesh, fi: usize, config: &BfsConfig, relaxed: bool) -> Vec<Seed> {
        if mesh.faces[fi].normal.is_none() {
            return Vec::new();
        }
        if mesh.faces[fi].conical_hypothesis != UNDEDUCED_CONICAL_HYPOTHESIS {
            return Vec::new();
        }
        // Skip faces already committed to a spherical hypothesis (stage 2.3
        // runs before this stage). On a tessellated sphere, every band of
        // quad rows lies exactly on some cone/cylinder, so without this check
        // a fine sphere gets tiled with hundreds of bogus band hypotheses —
        // enormously slow and all discarded by selection anyway. Faces near
        // sphere junctions can still be *absorbed* during BFS expansion;
        // only seeding is restricted.
        if mesh.faces[fi].spherical_hypothesis >= 0 {
            return Vec::new();
        }
        // Likewise skip faces committed to a toroidal hypothesis (stage 2.5
        // also runs before this stage): any adjacent pair of quad rows on a
        // torus tube lies exactly on some cone, so fillets would otherwise
        // be tiled with 2-ring cone bands. BFS expansion may still absorb
        // torus-committed faces at junctions.
        if mesh.faces[fi].toroidal_hypothesis >= 0 {
            return Vec::new();
        }
        // Junction faces (curvature transitions) produce contaminated fans;
        // pass 1 skips them as anchors and fan members. Pass 2 (relaxed)
        // lifts the filter for surfaces made entirely of transition faces.
        if config.is_junction(fi, relaxed) {
            return Vec::new();
        }

        let vc = mesh.faces[fi].vertex_count as usize;
        let angular_tol = config.angular_tol;

        // Collect qualified neighbors of fi
        let qualified_neighbors: Vec<usize> = mesh.faces[fi].neighbors[..vc]
            .iter()
            .filter(|&&n| n >= 0)
            .map(|&n| n as usize)
            .filter(|&ni| !config.is_junction(ni, relaxed))
            .filter(|&ni| conical_pairwise_qualified(mesh, fi, ni, angular_tol))
            .collect();

        // Generate 3-face seeds: fi + f1 + f2 where f1 and f2 are qualified
        // neighbors of fi but not neighbors of each other.
        // If the seed has only 5 unique vertices (all triangles), optionally
        // add f3 adjacent to f2 (not sharing an edge with fi or f1) to get
        // at least 6 vertices for reliable cone fitting.
        let mut seeds = Vec::new();
        for (i, &f1) in qualified_neighbors.iter().enumerate() {
            for &f2 in &qualified_neighbors[(i + 1)..] {
                // f2 must not be a neighbor of f1
                if super::bfs::is_neighbor(mesh, f1, f2) { continue; }

                // Count unique vertices across fi, f1, f2
                let mut vset: HashSet<usize> = HashSet::new();
                for &face_idx in &[fi, f1, f2] {
                    let fvc = mesh.faces[face_idx].vertex_count as usize;
                    for &vi in &mesh.faces[face_idx].vertex_indices[..fvc] {
                        vset.insert(vi);
                    }
                }

                if vset.len() >= 6 {
                    // Enough vertices for a reliable cone fit
                    seeds.push(Seed { faces: vec![fi, f1, f2] });
                } else {
                    // Need an f3 adjacent to f2, not sharing an edge with fi, f1, or f2
                    let f2_vc = mesh.faces[f2].vertex_count as usize;
                    let mut found_f3 = false;
                    for &f3n in &mesh.faces[f2].neighbors[..f2_vc] {
                        if f3n < 0 { continue; }
                        let f3 = f3n as usize;
                        if f3 == fi || f3 == f1 || f3 == f2 { continue; }
                        if super::bfs::is_neighbor(mesh, f3, fi) { continue; }
                        if super::bfs::is_neighbor(mesh, f3, f1) { continue; }
                        if mesh.faces[f3].conical_hypothesis != UNDEDUCED_CONICAL_HYPOTHESIS {
                            continue;
                        }
                        if mesh.faces[f3].normal.is_none() { continue; }
                        // Count vertices with f3 added
                        let mut vset4 = vset.clone();
                        let f3_vc = mesh.faces[f3].vertex_count as usize;
                        for &vi in &mesh.faces[f3].vertex_indices[..f3_vc] {
                            vset4.insert(vi);
                        }
                        if vset4.len() >= 6 {
                            seeds.push(Seed { faces: vec![fi, f1, f2, f3] });
                            found_f3 = true;
                            break;
                        }
                    }
                    if !found_f3 {
                        // Still emit the 3-face seed even with only 5 vertices;
                        // fit_surface may succeed for cylinders (4 DOF)
                        seeds.push(Seed { faces: vec![fi, f1, f2] });
                    }
                }
            }
        }
        seeds
    }

    fn fit_surface(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        vertex_list: &[usize],
        verbosity: u8,
    ) -> Option<Self::Fit> {
        // Sorted, deduplicated vertex list: fitting sums must iterate in a
        // deterministic order (float addition is non-associative; HashSet
        // order varies per process and produced run-to-run fit differences).
        let mut vertex_set: Vec<usize> = vertex_list.to_vec();
        vertex_set.sort_unstable();
        vertex_set.dedup();

        // Try the cheap 4-DOF cylinder fit first. If it is essentially exact
        // (max vertex error at CAD-mesh precision), skip the much more
        // expensive cone LM entirely: a cone could at best tie, and ties
        // select the cylinder anyway. This roughly halves fit cost on
        // cylinder-dominated meshes.
        let cyl_result = fit_cylinder(face_list, &vertex_set, &mesh.faces, &mesh.vertices);
        let cyl_fit = cyl_result.map(|(origin, dir, radius)| {
            ConeFit { apex: origin, axis_direction: dir, half_angle: 0.0, radius, convex: None }
        });

        // Compute max vertex error and max surface (centroid) error for cylinder fit
        let (cyl_vtx_max, cyl_sfc_max) = if let Some(ref fit) = cyl_fit {
            let vtx = vertex_list.iter().map(|&vi| {
                let v = &mesh.vertices[vi];
                vertex_to_cylinder_distance(v, &fit.apex, &fit.axis_direction, fit.radius).abs()
            }).fold(0.0_f64, f64::max);
            let sfc = face_list.iter().map(|&fi| {
                let c = super::face_centroid(&mesh.faces[fi], &mesh.vertices);
                let v = MeshVertex::from_xyz(c[0], c[1], c[2]);
                vertex_to_cylinder_distance(&v, &fit.apex, &fit.axis_direction, fit.radius).abs()
            }).fold(0.0_f64, f64::max);
            (Some(vtx), Some(sfc))
        } else {
            (None, None)
        };
        let _ = cyl_sfc_max; // used only for logging below

        const EXACT_CYLINDER_VTX_ERR: f64 = 1e-6; // mm; CAD-mesh vertex precision
        if let (Some(ref fit), Some(vtx)) = (&cyl_fit, cyl_vtx_max) {
            if vtx <= EXACT_CYLINDER_VTX_ERR {
                if verbosity >= 4 {
                    eprintln!(
                        "      cylinder fit exact (radius={:.6e}, max_vtx_err={:.6e}); skipping cone fit",
                        fit.radius, vtx
                    );
                }
                return cyl_fit;
            }
        }

        // Try cone fit
        let cone_result = fit_cone(face_list, &vertex_set, &mesh.faces, &mesh.vertices);
        let cone_fit = cone_result.and_then(|(apex, dir, ha)| {
            let ha_deg = ha.to_degrees();
            if !(MIN_HALF_ANGLE_DEG..=MAX_HALF_ANGLE_DEG).contains(&ha_deg) {
            if verbosity >= 4 {
                    eprintln!(
                        "      cone fit rejected: half_angle={:.2}\u{00b0} outside [{:.0}\u{00b0}, {:.0}\u{00b0}]",
                        ha_deg, MIN_HALF_ANGLE_DEG, MAX_HALF_ANGLE_DEG
                    );
                }
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
                if verbosity >= 4 {
                    eprintln!(
                        "      cone fit rejected: apex dist {:.4e} > {:.1} * bb_diag {:.4e}",
                        apex_dist, MAX_APEX_DISTANCE_FACTOR, bb_diag
                    );
                }
                return None;
            }
            Some(ConeFit { apex, axis_direction: dir, half_angle: ha, radius: 0.0, convex: None })
        });

        if cone_result.is_none() && verbosity >= 4 {
            eprintln!("      cone fit failed: fit_cone returned None");
        }

        // Compute max vertex error and max surface (centroid) error for cone fit
        let (cone_vtx_max, cone_sfc_max) = if let Some(ref fit) = cone_fit {
            let vtx = vertex_list.iter().map(|&vi| {
                let v = &mesh.vertices[vi];
                vertex_to_cone_distance(v, &fit.apex, &fit.axis_direction, fit.half_angle).abs()
            }).fold(0.0_f64, f64::max);
            let sfc = face_list.iter().map(|&fi| {
                let c = super::face_centroid(&mesh.faces[fi], &mesh.vertices);
                let v = MeshVertex::from_xyz(c[0], c[1], c[2]);
                vertex_to_cone_distance(&v, &fit.apex, &fit.axis_direction, fit.half_angle).abs()
            }).fold(0.0_f64, f64::max);
            if verbosity >= 4 {
                eprintln!(
                    "      cone fit: half_angle={:.2}\u{00b0}, max_vtx_err={:.6e}, max_sfc_err={:.6e}",
                    fit.half_angle.to_degrees(), vtx, sfc
                );
            }
            (Some(vtx), Some(sfc))
        } else {
            (None, None)
        };

        if cyl_result.is_none() && verbosity >= 4 && cone_fit.is_none() {
            eprintln!("      cylinder fit failed: fit_cylinder returned None");
        }
        if verbosity >= 4 && cone_fit.is_none() {
            if let (Some(ref fit), Some(vtx), Some(sfc)) = (&cyl_fit, cyl_vtx_max, cyl_sfc_max) {
                eprintln!(
                    "      cylinder fit: radius={:.6e}, max_vtx_err={:.6e}, max_sfc_err={:.6e}",
                    fit.radius, vtx, sfc
                );
            }
        }
        let _ = cone_sfc_max; // used only for logging above

        // Pick the better fit by max vertex error
        match (cyl_fit, cone_fit) {
            (Some(cf), Some(kf)) => {
                let cv = cyl_vtx_max.unwrap();
                let kv = cone_vtx_max.unwrap();
                if cv <= kv {
                    if verbosity >= 4 {
                        eprintln!(
                            "      selected cylinder (vtx_err {:.6e} <= cone {:.6e})",
                            cv, kv
                        );
                    }
                    Some(cf)
                } else {
                    if verbosity >= 4 {
                        eprintln!(
                            "      selected cone (vtx_err {:.6e} < cylinder {:.6e})",
                            kv, cv
                        );
                    }
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
        // Normal-axis consistency: on a cone, every facet normal satisfies
        // n·axis = sin(half_angle) exactly — tessellation coarseness rotates
        // facet normals *around* the axis, preserving the dot product. Faces
        // on spheres and tori violate this as soon as BFS crosses tube rows,
        // so checking here keeps doomed trials small instead of letting them
        // grow across a whole fillet and only failing the post-BFS variance
        // test (which made stage 2.2 quadratic in fillet size).
        if let Some(n) = mesh.faces[fi].normal {
            let n_dot_axis = n[0] * fit.axis_direction[0]
                + n[1] * fit.axis_direction[1]
                + n[2] * fit.axis_direction[2];
            let expected = fit.half_angle.sin();
            if (n_dot_axis.abs() - expected.abs()).abs() > 0.05 {
                return false;
            }
        }
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
        // Reject cones with very large half-angle (nearly a plane). Real CAD
        // cones up to 85° exist (cone_tangent_planes has an 84.25° cone), and
        // the historical worry — near-tangent bogus cones on sphere polar
        // caps — is now prevented upstream: spherical detection runs before
        // conical and sphere-committed faces don't seed cones.
        if fit.half_angle > 85.0_f64.to_radians() {
            return false;
        }
        // Normal-axis consistency: on a cone, n·axis = sin(half_angle) is
        // constant (0 for cylinders). Reject hypotheses where the normal dot
        // products with the axis vary too much — the faces lie on a sphere or
        // torus, not a cone/cylinder. Without this check, stage 2.2 tiles
        // spheres and torus fillets with hundreds of small bogus conical
        // hypotheses. Area-based selection (2.6) would discard them, but they
        // starve torus seeding (which excludes conical-committed faces) and
        // waste large amounts of BFS+LM time on fine meshes.
        let dots: Vec<f64> = face_list.iter()
            .filter_map(|&fi| mesh.faces[fi].normal.map(|n| {
                n[0] * fit.axis_direction[0]
                    + n[1] * fit.axis_direction[1]
                    + n[2] * fit.axis_direction[2]
            }))
            .collect();
        if dots.len() >= 3 {
            let mean = dots.iter().sum::<f64>() / dots.len() as f64;
            let var = dots.iter().map(|d| (d - mean).powi(2)).sum::<f64>()
                / dots.len() as f64;
            if var > 0.0005 {
                return false;
            }
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
    ) -> (Vec<viz::CylinderOverlay>, Vec<viz::ConeOverlay>, Vec<viz::SphereOverlay>, Vec<viz::TorusOverlay>) {
        if fit.half_angle.abs() < 1e-6 {
            // Cylinder
            let cyl = centered_cylinder_overlay(
                fit.apex, fit.axis_direction, fit.radius,
                face_list, mesh,
                [0.2, 0.4, 1.0, 0.3],
            );
            (vec![cyl], Vec::new(), Vec::new(), Vec::new())
        } else {
            // True cone
            let cone = centered_cone_overlay(
                fit.apex, fit.axis_direction, fit.half_angle,
                face_list, mesh,
                [0.2, 0.4, 1.0, 0.3],
            );
            (Vec::new(), vec![cone], Vec::new(), Vec::new())
        }
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
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

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
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

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
        let vertex_set: Vec<usize> = (0..2).collect();
        assert!(fit_cone(&[0], &vertex_set, &faces, &vertices).is_none(),
                "cone fit should reject degenerate input");
    }

    /// Build a pointed cone mesh (comes to a point at the apex).
    /// The cone has `n_slices` triangular faces, radius `base_radius` at height
    /// `height` from the apex, axis along (0,0,1), apex at origin.
    fn build_pointed_cone_mesh(
        n_slices: usize,
        base_radius: f64,
        height: f64,
    ) -> (Vec<MeshFace>, Vec<MeshVertex>) {
        // Vertex 0 is the apex at origin.
        // Vertices 1..=n_slices are base ring.
        let mut vertices = Vec::with_capacity(n_slices + 1);
        vertices.push(MeshVertex::from_xyz(0.0, 0.0, 0.0)); // apex

        for i in 0..n_slices {
            let theta = 2.0 * PI * i as f64 / n_slices as f64;
            let x = base_radius * theta.cos();
            let y = base_radius * theta.sin();
            vertices.push(MeshVertex::from_xyz(x, y, height));
        }

        let mut faces = Vec::with_capacity(n_slices);
        let _half_angle: f64 = (base_radius / height).atan();

        for i in 0..n_slices {
            let v1 = 1 + i;
            let v2 = 1 + (i + 1) % n_slices;

            // Compute outward normal: average of the two edge normals.
            // For a triangular face (apex, v1, v2), the outward normal is:
            //   n = normalize((v1 + v2) / 2 projected outward)
            // More precisely, the face normal from cross product.
            let a = [vertices[v1].x, vertices[v1].y, vertices[v1].z];
            let b = [vertices[v2].x, vertices[v2].y, vertices[v2].z];
            let edge1 = [a[0], a[1], a[2]]; // apex to v1
            let edge2 = [b[0], b[1], b[2]]; // apex to v2
            let mut normal = cross3(&edge1, &edge2);
            normalize3(&mut normal);
            // Ensure normal points outward (radially away from axis).
            // Use radial direction of base midpoint, ignoring axial component.
            let mid_base = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, 0.0];
            if dot3(&normal, &mid_base) < 0.0 {
                normal = [-normal[0], -normal[1], -normal[2]];
            }

            let prev = if i == 0 { n_slices - 1 } else { i - 1 };
            let next = (i + 1) % n_slices;

            faces.push(MeshFace {
                vertex_count: 3,
                vertex_indices: [0, v1, v2, 0],
                neighbors: [prev as i32, next as i32, -1, -1],
                normal: Some(normal),
                planar_hypothesis: -1,
                spherical_hypothesis: -1,
                conical_hypothesis: -2,
                toroidal_hypothesis: -2,
            });
        }

        (faces, vertices)
    }

    #[test]
    fn cone_fit_pointed_cone() {
        // Pointed cone: 10 triangular faces, base radius 5, height 15.
        // Apex at origin, axis (0,0,1), half_angle = atan(5/15) = atan(1/3).
        let n_slices = 10;
        let base_radius: f64 = 5.0;
        let height: f64 = 15.0;
        let expected_ha = (base_radius / height).atan();

        let (faces, vertices) = build_pointed_cone_mesh(n_slices, base_radius, height);
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_cone(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "pointed cone fit should succeed");
        let (fit_apex, fit_dir, fit_ha) = result.unwrap();

        // Half-angle check
        assert!(
            (fit_ha - expected_ha).abs() < 0.05,
            "half_angle should be ~{:.4} rad, got {:.4}",
            expected_ha, fit_ha
        );
        // Axis should be along Z
        let axis_z = fit_dir[2].abs();
        assert!(
            (axis_z - 1.0).abs() < 0.05,
            "axis should be along Z, got dir=({:.3},{:.3},{:.3})",
            fit_dir[0], fit_dir[1], fit_dir[2]
        );
        // Apex should be near origin
        let apex_dist = (fit_apex[0] * fit_apex[0]
            + fit_apex[1] * fit_apex[1]
            + fit_apex[2] * fit_apex[2]).sqrt();
        assert!(
            apex_dist < 1.0,
            "apex should be near origin, got ({:.3},{:.3},{:.3}), dist={:.3}",
            fit_apex[0], fit_apex[1], fit_apex[2], apex_dist
        );
    }

    #[test]
    fn debug_axis_estimation() {
        let n_slices = 10;
        let base_radius: f64 = 5.0;
        let height: f64 = 15.0;
        let (faces, vertices) = build_pointed_cone_mesh(n_slices, base_radius, height);

        // Print normals for first few faces
        for i in 0..3 {
            let n = faces[i].normal.unwrap();
            eprintln!("face {i} normal: ({:.4}, {:.4}, {:.4})", n[0], n[1], n[2]);
        }

        // Test normal-diff axis for faces 0,1,2
        let fi = vec![0, 1, 2];
        let axis = estimate_axis_from_normal_diffs(&fi, &faces);
        eprintln!("axis from normal diffs (0,1,2): {:?}", axis);

        // Also test apex estimation
        let apex = estimate_apex(&fi, &faces, &vertices);
        eprintln!("apex (0,1,2): {:?}", apex);

        // Also test half-angle
        if let (Some(a), Some(mut ax)) = (apex, axis) {
            let mut vset: Vec<usize> = Vec::new();
            for &f in &fi {
                let vc = faces[f].vertex_count as usize;
                vset.extend_from_slice(&faces[f].vertex_indices[..vc]);
            }
            vset.sort_unstable();
            vset.dedup();
            let ha = estimate_half_angle(&a, &mut ax, &vset, &vertices);
            eprintln!("half_angle (0,1,2): {:?} (expected {:.4})", ha, (base_radius/height).atan());
            eprintln!("axis after orient: ({:.4}, {:.4}, {:.4})", ax[0], ax[1], ax[2]);
        }
    }

    #[test]
    fn cone_fit_pointed_cone_subset() {
        // Fit any 3 or 4 consecutive faces of the pointed cone;
        // the fit should still match the reference cone.
        let n_slices = 10;
        let base_radius: f64 = 5.0;
        let height: f64 = 15.0;
        let expected_ha = (base_radius / height).atan();

        let (faces, vertices) = build_pointed_cone_mesh(n_slices, base_radius, height);

        for window_size in [3, 4] {
            for start in 0..n_slices {
                let face_indices: Vec<usize> = (0..window_size)
                    .map(|i| (start + i) % n_slices)
                    .collect();
                let mut vertex_set: Vec<usize> = Vec::new();
                for &fi in &face_indices {
                    let vc = faces[fi].vertex_count as usize;
                    vertex_set.extend_from_slice(&faces[fi].vertex_indices[..vc]);
                }
                vertex_set.sort_unstable();
                vertex_set.dedup();

                let result = fit_cone(&face_indices, &vertex_set, &faces, &vertices);
                assert!(
                    result.is_some(),
                    "pointed cone subset fit failed for window_size={window_size} start={start}"
                );
                let (fit_apex, fit_dir, fit_ha) = result.unwrap();

                assert!(
                    (fit_ha - expected_ha).abs() < 0.1,
                    "subset({window_size},{start}): half_angle {:.4} vs expected {:.4}",
                    fit_ha, expected_ha,
                );
                let axis_z = fit_dir[2].abs();
                assert!(
                    (axis_z - 1.0).abs() < 0.1,
                    "subset({window_size},{start}): axis z={axis_z:.3}",
                );
                let apex_dist = (fit_apex[0].powi(2) + fit_apex[1].powi(2) + fit_apex[2].powi(2)).sqrt();
                assert!(
                    apex_dist < 2.0,
                    "subset({window_size},{start}): apex dist={apex_dist:.3}",
                );
            }
        }
    }
}

