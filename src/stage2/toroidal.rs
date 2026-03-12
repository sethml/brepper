//! Stage 2.5: Toroidal hypothesis deduction and comparison.
//!
//! Uses the medial axis / tube center method:
//! 1. Estimate minor radius from normal-line intersections of adjacent face pairs
//! 2. Compute tube centers k_i = p_i ± r*n_i
//! 3. PCA of tube centers → axis direction
//! 4. 2D algebraic circle fit → center and major radius R
//! 5. LM refinement of 7 parameters [Cx, Cy, Cz, alpha, beta, R, r]
//! 6. BFS expansion with vertex-to-torus distance validation

use std::collections::{HashSet, VecDeque};

use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{Const, Dyn, OMatrix, OVector, Owned};
use opencascade_sys::gp;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_TOROIDAL_HYPOTHESIS, NO_HYPOTHESIS};
use crate::viz::{self, VizAction, VizSender};

use super::{
    build_normal_covariance, dot3, face_centroid,
    normalize3, perpendicular_basis, smallest_eigenvector_3x3, viz_custom, viz_face_centroid,
    viz_face_normal, ToroidalHypothesis, Stage2CompareError,
    MIN_TORUS_FACES, REFIT_SKIP_MULTIPLIER,
};

// ---------------------------------------------------------------------------
// Torus distance
// ---------------------------------------------------------------------------

/// Compute the signed distance from a vertex to a torus surface.
///
/// For a torus with center C, axis â, major radius R, minor radius r:
///   axial  = (v - C) · â
///   radial = |(v - C) - axial·â|
///   d = sqrt((radial - R)² + axial²) - r
///
/// Positive if outside the tube, negative if inside.
pub(super) fn vertex_to_torus_distance(
    v: &MeshVertex,
    center: &[f64; 3],
    axis_direction: &[f64; 3],
    major_radius: f64,
    minor_radius: f64,
) -> f64 {
    let dx = v.x - center[0];
    let dy = v.y - center[1];
    let dz = v.z - center[2];
    let axial = dx * axis_direction[0] + dy * axis_direction[1] + dz * axis_direction[2];
    let rx = dx - axial * axis_direction[0];
    let ry = dy - axial * axis_direction[1];
    let rz = dz - axial * axis_direction[2];
    let radial = (rx * rx + ry * ry + rz * rz).sqrt();
    ((radial - major_radius).powi(2) + axial * axial).sqrt() - minor_radius
}

// ---------------------------------------------------------------------------
// Normal-line intersection for minor radius estimation
// ---------------------------------------------------------------------------

/// For two surface points with normals, compute the closest-approach parameter
/// along the normal lines. Returns (t, s, distance) where t is parameter on
/// line P1+t*n1, s is parameter on P2+s*n2, and distance is the closest-approach
/// distance between the two lines.
fn normal_line_intersection(
    p1: &[f64; 3], n1: &[f64; 3],
    p2: &[f64; 3], n2: &[f64; 3],
) -> Option<(f64, f64, f64)> {
    let a = dot3(n1, n2);
    let denom = 1.0 - a * a;
    if denom.abs() < 1e-12 {
        return None; // Nearly parallel normals
    }
    let dp = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
    let c = dot3(n1, &dp);
    let d = dot3(n2, &dp);
    let t = (a * d - c) / denom;
    let s = (d - a * c) / denom;

    // Closest approach point on each line
    let q1 = [p1[0] + t * n1[0], p1[1] + t * n1[1], p1[2] + t * n1[2]];
    let q2 = [p2[0] + s * n2[0], p2[1] + s * n2[1], p2[2] + s * n2[2]];
    let dist = ((q1[0] - q2[0]).powi(2) + (q1[1] - q2[1]).powi(2) + (q1[2] - q2[2]).powi(2)).sqrt();
    Some((t, s, dist))
}

// ---------------------------------------------------------------------------
// 3D circle fitting from tube centers
// ---------------------------------------------------------------------------

/// Fit a 3D circle to a set of points using PCA + 2D algebraic circle fit.
/// Returns (center, axis_direction, radius) or None if degenerate.
fn fit_circle_3d(points: &[[f64; 3]]) -> Option<([f64; 3], [f64; 3], f64)> {
    let n = points.len();
    if n < 3 {
        return None;
    }

    // Compute centroid
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    for p in points {
        cx += p[0];
        cy += p[1];
        cz += p[2];
    }
    let inv_n = 1.0 / n as f64;
    cx *= inv_n;
    cy *= inv_n;
    cz *= inv_n;

    // Build covariance matrix of positions (for PCA)
    let mut cov = [0.0_f64; 6]; // upper triangle [m00, m01, m02, m11, m12, m22]
    for p in points {
        let dx = p[0] - cx;
        let dy = p[1] - cy;
        let dz = p[2] - cz;
        cov[0] += dx * dx;
        cov[1] += dx * dy;
        cov[2] += dx * dz;
        cov[3] += dy * dy;
        cov[4] += dy * dz;
        cov[5] += dz * dz;
    }

    // The circle normal is the eigenvector of the smallest eigenvalue
    // (the points should lie in a plane, with minimal variance in the normal direction)
    let axis = smallest_eigenvector_3x3(&cov);

    // Project points onto the plane perpendicular to axis, through centroid
    let (u, w) = perpendicular_basis(&axis);

    // 2D coordinates in the plane
    let coords_2d: Vec<(f64, f64)> = points.iter().map(|p| {
        let dx = p[0] - cx;
        let dy = p[1] - cy;
        let dz = p[2] - cz;
        (dot3(&[dx, dy, dz], &u), dot3(&[dx, dy, dz], &w))
    }).collect();

    // Algebraic circle fit in 2D: minimize Σ (xi² + yi² + D*xi + E*yi + F)²
    // Linear system: [Σxi² Σxiyi Σxi] [D]   [-Σxi(xi²+yi²)]
    //                [Σxiyi Σyi² Σyi] [E] = [-Σyi(xi²+yi²)]
    //                [Σxi   Σyi  n  ] [F]   [-Σ(xi²+yi²)  ]
    let mut sx = 0.0_f64;
    let mut sy = 0.0_f64;
    let mut sxx = 0.0_f64;
    let mut syy = 0.0_f64;
    let mut sxy = 0.0_f64;
    let mut sxr = 0.0_f64;
    let mut syr = 0.0_f64;
    let mut sr = 0.0_f64;

    for &(x, y) in &coords_2d {
        let r2 = x * x + y * y;
        sx += x;
        sy += y;
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
        sxr += x * r2;
        syr += y * r2;
        sr += r2;
    }

    let nf = n as f64;
    // Solve 3x3 system using Cramer's rule
    let a11 = sxx; let a12 = sxy; let a13 = sx;
    let a21 = sxy; let a22 = syy; let a23 = sy;
    let a31 = sx;  let a32 = sy;  let a33 = nf;
    let b1 = -sxr; let b2 = -syr; let b3 = -sr;

    let det = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);
    if det.abs() < 1e-30 {
        return None;
    }
    let inv_det = 1.0 / det;

    let d = inv_det * (b1 * (a22 * a33 - a23 * a32) - a12 * (b2 * a33 - a23 * b3) + a13 * (b2 * a32 - a22 * b3));
    let e = inv_det * (a11 * (b2 * a33 - a23 * b3) - b1 * (a21 * a33 - a23 * a31) + a13 * (a21 * b3 - b2 * a31));
    let f = inv_det * (a11 * (a22 * b3 - b2 * a32) - a12 * (a21 * b3 - b2 * a31) + b1 * (a21 * a32 - a22 * a31));

    // Circle center in 2D: (-D/2, -E/2), radius = sqrt(D²/4 + E²/4 - F)
    let center_2d_x = -d / 2.0;
    let center_2d_y = -e / 2.0;
    let r_sq = d * d / 4.0 + e * e / 4.0 - f;
    if r_sq <= 0.0 {
        return None;
    }
    let radius = r_sq.sqrt();

    // Convert center back to 3D
    let center_3d = [
        cx + center_2d_x * u[0] + center_2d_y * w[0],
        cy + center_2d_x * u[1] + center_2d_y * w[1],
        cz + center_2d_x * u[2] + center_2d_y * w[2],
    ];

    Some((center_3d, axis, radius))
}

// ---------------------------------------------------------------------------
// LM refinement for torus fitting
// ---------------------------------------------------------------------------

/// LM problem for torus fitting.
/// Parameters: [alpha, beta, cx, cy, cz, R, r]
///   alpha, beta: tilt of axis from initial direction
///   cx, cy, cz: center position
///   R: major radius
///   r: minor radius
struct TorusLMProblem {
    points: Vec<[f64; 3]>,
    a0: [f64; 3],
    u0: [f64; 3],
    w0: [f64; 3],
    params: OVector<f64, Const<7>>,
}

impl TorusLMProblem {
    fn new(
        points: Vec<[f64; 3]>,
        initial_axis: [f64; 3],
        initial_center: [f64; 3],
        initial_major_r: f64,
        initial_minor_r: f64,
    ) -> Self {
        let (u0, w0) = perpendicular_basis(&initial_axis);
        let params = OVector::<f64, Const<7>>::from_column_slice(&[
            0.0, 0.0,
            initial_center[0], initial_center[1], initial_center[2],
            initial_major_r, initial_minor_r,
        ]);
        Self { points, a0: initial_axis, u0, w0, params }
    }

    fn axis_dir_from_params(&self, params: &OVector<f64, Const<7>>) -> [f64; 3] {
        let (alpha, beta) = (params[0], params[1]);
        let mut v = [
            self.a0[0] + alpha * self.u0[0] + beta * self.w0[0],
            self.a0[1] + alpha * self.u0[1] + beta * self.w0[1],
            self.a0[2] + alpha * self.u0[2] + beta * self.w0[2],
        ];
        normalize3(&mut v);
        v
    }

    fn compute_residuals_for(&self, params: &OVector<f64, Const<7>>) -> OVector<f64, Dyn> {
        let a = self.axis_dir_from_params(params);
        let center = [params[2], params[3], params[4]];
        let major_r = params[5];
        let minor_r = params[6];
        let n = self.points.len();
        let mut residuals = OVector::<f64, Dyn>::zeros_generic(Dyn(n), Const::<1>);
        for (i, p) in self.points.iter().enumerate() {
            let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
            let axial = dot3(&d, &a);
            let rx = d[0] - axial * a[0];
            let ry = d[1] - axial * a[1];
            let rz = d[2] - axial * a[2];
            let radial = (rx * rx + ry * ry + rz * rz).sqrt();
            residuals[i] = ((radial - major_r).powi(2) + axial * axial).sqrt() - minor_r;
        }
        residuals
    }
}

impl LeastSquaresProblem<f64, Dyn, Const<7>> for TorusLMProblem {
    type ParameterStorage = Owned<f64, Const<7>>;
    type ResidualStorage = Owned<f64, Dyn>;
    type JacobianStorage = Owned<f64, Dyn, Const<7>>;

    fn set_params(&mut self, x: &OVector<f64, Const<7>>) {
        self.params.copy_from(x);
    }

    fn params(&self) -> OVector<f64, Const<7>> {
        self.params
    }

    fn residuals(&self) -> Option<OVector<f64, Dyn>> {
        Some(self.compute_residuals_for(&self.params))
    }

    fn jacobian(&self) -> Option<OMatrix<f64, Dyn, Const<7>>> {
        let n = self.points.len();
        let eps = 1e-8;
        let mut jac = OMatrix::<f64, Dyn, Const<7>>::zeros_generic(Dyn(n), Const::<7>);
        for j in 0..7 {
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

// ---------------------------------------------------------------------------
// Torus fitting from faces
// ---------------------------------------------------------------------------

/// Fit a torus to a set of faces using the medial axis / tube center method.
///
/// 1. Estimate minor radius r from normal-line intersections of adjacent face pairs.
/// 2. Compute tube centers k_i = p_i - sign*r*n_i for each vertex.
/// 3. Fit a 3D circle to tube centers → center, axis, major radius R.
/// 4. LM refinement.
///
/// Returns (center, axis_direction, major_radius, minor_radius) or None.
fn fit_torus(
    face_indices: &[usize],
    vertex_set: &HashSet<usize>,
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64, f64)> {
    if vertex_set.len() < 7 || face_indices.len() < 3 {
        return None;
    }

    // Step 1: Estimate minor radius from normal-line intersections of adjacent face pairs.
    let mut r_estimates: Vec<f64> = Vec::new();
    let face_set: HashSet<usize> = face_indices.iter().copied().collect();

    for &fi in face_indices {
        let n1 = match faces[fi].normal {
            Some(n) => n,
            None => continue,
        };
        let c1 = face_centroid(&faces[fi], vertices);
        let vc = faces[fi].vertex_count as usize;

        for &ni in &faces[fi].neighbors[..vc] {
            if ni < 0 { continue; }
            let ni = ni as usize;
            if !face_set.contains(&ni) { continue; }
            if ni <= fi { continue; } // avoid duplicates

            let n2 = match faces[ni].normal {
                Some(n) => n,
                None => continue,
            };
            let c2 = face_centroid(&faces[ni], vertices);

            // Normal-line intersection
            if let Some((t, s, dist)) = normal_line_intersection(&c1, &n1, &c2, &n2) {
                // For a torus, both t and s should have the same sign (both pointing
                // toward the tube center) and similar magnitude ≈ minor_radius.
                // The closest-approach distance should be small relative to |t|.
                let t_abs = t.abs();
                let s_abs = s.abs();
                if t_abs < 1e-10 || s_abs < 1e-10 { continue; }

                // t and s should have the same sign for consistent curvature
                if t * s < 0.0 { continue; }

                // The ratio |t|/|s| should be close to 1 for a torus
                let ratio = t_abs / s_abs;
                if !(0.3..=3.0).contains(&ratio) { continue; }

                // The closest-approach distance should be small relative to the radius
                let avg_r = (t_abs + s_abs) / 2.0;
                if dist > avg_r * 0.5 { continue; }

                r_estimates.push(avg_r);
            }
        }
    }

    if r_estimates.is_empty() {
        // Fallback: try SOR axis + profile fitting approach
        return fit_torus_sor(face_indices, vertex_set, faces, vertices);
    }

    // Robust minor radius: median of estimates
    r_estimates.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let minor_r = r_estimates[r_estimates.len() / 2];

    if minor_r < 1e-10 {
        return None;
    }

    // Step 2: Determine tube center sign from face normals.
    // For convex torus (outer surface), normals point outward from tube center,
    // so tube_center = p - r*n. For concave (inner surface), normals point
    // toward tube center, so tube_center = p + r*n.
    // We try both signs, and use the one that produces a better circle fit.
    let best = [1.0_f64, -1.0].iter().filter_map(|&sign| {
        // Compute tube centers
        let tube_centers: Vec<[f64; 3]> = vertex_set.iter().filter_map(|&vi| {
            // Find a face incident to this vertex to get its normal
            let face = face_indices.iter().find(|&&fi| {
                let vc = faces[fi].vertex_count as usize;
                faces[fi].vertex_indices[..vc].contains(&vi)
            })?;
            let n = faces[*face].normal?;
            let v = &vertices[vi];
            Some([
                v.x + sign * minor_r * n[0],
                v.y + sign * minor_r * n[1],
                v.z + sign * minor_r * n[2],
            ])
        }).collect();

        if tube_centers.len() < 7 { return None; }

        // Step 3: Fit 3D circle to tube centers
        let (center, axis, major_r) = fit_circle_3d(&tube_centers)?;
        if major_r < minor_r * 0.1 {
            return None; // Degenerate: major radius too small
        }

        // Compute RMS residual of tube centers from the fitted circle
        let rms = {
            let mut sum_sq = 0.0_f64;
            for tc in &tube_centers {
                let d = [tc[0] - center[0], tc[1] - center[1], tc[2] - center[2]];
                let axial = dot3(&d, &axis);
                let radial_vec = [d[0] - axial * axis[0], d[1] - axial * axis[1], d[2] - axial * axis[2]];
                let radial = (radial_vec[0].powi(2) + radial_vec[1].powi(2) + radial_vec[2].powi(2)).sqrt();
                let err = ((radial - major_r).powi(2) + axial * axial).sqrt();
                sum_sq += err * err;
            }
            (sum_sq / tube_centers.len() as f64).sqrt()
        };

        Some((center, axis, major_r, rms, sign))
    }).min_by(|a, b| a.3.partial_cmp(&b.3).unwrap())?;

    let (center, axis, major_r, rms, _sign) = best;

    // Quality gate: reject if tube centers don't form a clean circle.
    // This catches mixed convex+concave seeds where tube centers scatter
    // between two different major circles. Use tight threshold since valid
    // torus fits have tube centers very close to the major circle.
    if rms > minor_r * 0.1 {
        return None;
    }

    // Step 4: LM refinement
    let points: Vec<[f64; 3]> = vertex_set.iter()
        .map(|&vi| [vertices[vi].x, vertices[vi].y, vertices[vi].z])
        .collect();

    if points.len() >= 10 {
        let problem = TorusLMProblem::new(points, axis, center, major_r, minor_r);
        let (result, _report) = LevenbergMarquardt::new().minimize(problem);
        let refined_axis = result.axis_dir_from_params(&result.params);
        let refined_center = [result.params[2], result.params[3], result.params[4]];
        let refined_major_r = result.params[5];
        let refined_minor_r = result.params[6];

        if refined_major_r > 0.0 && refined_minor_r > 0.0 {
            return Some((refined_center, refined_axis, refined_major_r, refined_minor_r));
        }
    }

    Some((center, axis, major_r, minor_r))
}

/// Fallback: fit torus using surface-of-revolution axis estimation + (h,r) profile fitting.
/// This works when the normal-line intersection approach fails (e.g., too few face pairs).
fn fit_torus_sor(
    face_indices: &[usize],
    vertex_set: &HashSet<usize>,
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64, f64)> {
    if vertex_set.len() < 7 || face_indices.len() < 3 {
        return None;
    }

    // Estimate axis from normal covariance (smallest eigenvector)
    let cov = build_normal_covariance(face_indices, faces, vertices);
    let mut axis = smallest_eigenvector_3x3(&cov);
    let len = normalize3(&mut axis);
    if len < 1e-15 {
        return None;
    }

    // Compute centroid of vertices
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    let n = vertex_set.len();
    for &vi in vertex_set {
        cx += vertices[vi].x;
        cy += vertices[vi].y;
        cz += vertices[vi].z;
    }
    let inv_n = 1.0 / n as f64;
    cx *= inv_n;
    cy *= inv_n;
    cz *= inv_n;

    // Compute (h, r) coordinates relative to centroid
    let mut hrs: Vec<(f64, f64)> = Vec::with_capacity(n);
    for &vi in vertex_set {
        let dx = vertices[vi].x - cx;
        let dy = vertices[vi].y - cy;
        let dz = vertices[vi].z - cz;
        let h = dot3(&[dx, dy, dz], &axis);
        let rx = dx - h * axis[0];
        let ry = dy - h * axis[1];
        let rz = dz - h * axis[2];
        let r = (rx * rx + ry * ry + rz * rz).sqrt();
        hrs.push((h, r));
    }

    // Fit torus profile: (r - R)² + h² = r_minor²
    // Expand: r² - 2Rr + R² + h² = r_minor²
    // Let s = r² + h², then: s - 2Rr + R² - r_minor² = 0
    // Linear in unknowns: s = 2Rr + (r_minor² - R²) = 2R*r + C
    // So: s_i = 2R * r_i + C  →  linear regression of s vs r.
    let mut sum_r = 0.0_f64;
    let mut sum_s = 0.0_f64;
    let mut sum_rr = 0.0_f64;
    let mut sum_rs = 0.0_f64;
    let nf = n as f64;

    for &(h, r) in &hrs {
        let s = r * r + h * h;
        sum_r += r;
        sum_s += s;
        sum_rr += r * r;
        sum_rs += r * s;
    }

    let denom = nf * sum_rr - sum_r * sum_r;
    if denom.abs() < 1e-30 {
        return None;
    }

    let slope = (nf * sum_rs - sum_r * sum_s) / denom; // slope = 2R
    let intercept = (sum_s - slope * sum_r) / nf; // intercept = r_minor² - R²

    let major_r = slope / 2.0;
    if major_r <= 0.0 {
        return None;
    }

    let r_minor_sq = intercept + major_r * major_r;
    if r_minor_sq <= 0.0 {
        return None;
    }
    let minor_r = r_minor_sq.sqrt();

    // Center is at the vertex centroid projected onto the axis plane
    // Actually center is where the axis passes through the major circle plane.
    // For the centroid-based coordinate system, center is at centroid + correction.
    // The mean radial distance should equal the major radius.
    let mean_r = sum_r / nf;
    // The center offset from centroid is zero in the plane, but we need to check
    // if the centroid is already at the axis.
    let center = [cx, cy, cz];

    // LM refinement
    let points: Vec<[f64; 3]> = vertex_set.iter()
        .map(|&vi| [vertices[vi].x, vertices[vi].y, vertices[vi].z])
        .collect();

    if points.len() >= 10 && mean_r > 1e-10 {
        let problem = TorusLMProblem::new(points, axis, center, major_r, minor_r);
        let (result, _report) = LevenbergMarquardt::new().minimize(problem);
        let refined_axis = result.axis_dir_from_params(&result.params);
        let refined_center = [result.params[2], result.params[3], result.params[4]];
        let refined_major_r = result.params[5];
        let refined_minor_r = result.params[6];

        if refined_major_r > 0.0 && refined_minor_r > 0.0 {
            return Some((refined_center, refined_axis, refined_major_r, refined_minor_r));
        }
    }

    Some((center, axis, major_r, minor_r))
}

// ---------------------------------------------------------------------------
// Validation helpers
// ---------------------------------------------------------------------------

/// Check if all vertices in a set are within tolerance of a torus surface.
fn all_vertices_within_torus_tolerance(
    vertex_set: &HashSet<usize>,
    center: &[f64; 3],
    axis_direction: &[f64; 3],
    major_radius: f64,
    minor_radius: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        let d = vertex_to_torus_distance(&vertices[vi], center, axis_direction, major_radius, minor_radius);
        if d.abs() > tolerance {
            return false;
        }
    }
    true
}

/// Determine convexity: does the face normal point away from the tube center (convex)
/// or toward it (concave)?
fn determine_convexity(
    face: &MeshFace,
    vertices: &[MeshVertex],
    center: &[f64; 3],
    axis_direction: &[f64; 3],
    major_radius: f64,
) -> bool {
    let centroid = face_centroid(face, vertices);
    let d = [
        centroid[0] - center[0],
        centroid[1] - center[1],
        centroid[2] - center[2],
    ];
    let axial = dot3(&d, axis_direction);
    let radial = [
        d[0] - axial * axis_direction[0],
        d[1] - axial * axis_direction[1],
        d[2] - axial * axis_direction[2],
    ];
    let radial_dist = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();
    // Tube center for this face lies on the major circle at distance R in the radial direction
    let tube_radial = if radial_dist > 1e-15 {
        [
            major_radius * radial[0] / radial_dist,
            major_radius * radial[1] / radial_dist,
            major_radius * radial[2] / radial_dist,
        ]
    } else {
        [0.0; 3]
    };
    // Vector from tube center to face centroid
    let to_face = [
        d[0] - tube_radial[0],
        d[1] - tube_radial[1],
        d[2] - tube_radial[2],
    ];
    let n = face.normal.unwrap();
    // Convex: normal points away from tube center (dot > 0)
    dot3(&n, &to_face) > 0.0
}

/// Validate angular coverage around the torus (same gap-based check as cylindrical).
fn angular_coverage_valid(
    face_list: &[usize],
    center: &[f64; 3],
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
        let d = [c[0] - center[0], c[1] - center[1], c[2] - center[2]];
        let radial = [
            d[0] - dot3(&d, axis_direction) * axis_direction[0],
            d[1] - dot3(&d, axis_direction) * axis_direction[1],
            d[2] - dot3(&d, axis_direction) * axis_direction[2],
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
// Trial BFS for toroidal hypothesis evaluation
// ---------------------------------------------------------------------------

struct TorusTrialResult {
    faces: Vec<usize>,
    vertices: HashSet<usize>,
    center: [f64; 3],
    axis_direction: [f64; 3],
    major_radius: f64,
    minor_radius: f64,
    convex: bool,
    error_max: f64,
    centroid_error_max: f64,
    error_abs_sum: f64,
}

/// Run a trial BFS for a toroidal hypothesis starting from seed faces.
#[allow(clippy::too_many_arguments)]
fn run_torus_trial_bfs(
    seed_faces: &[usize],
    mesh: &ConnectedMesh,
    surface_tol: f64,
    angular_tol: f64,
    verbosity: u8,
    viz: Option<&VizSender>,
    viz_quit: &std::cell::Cell<bool>,
) -> Option<TorusTrialResult> {
    // Collect seed vertices
    let mut vertex_set: HashSet<usize> = HashSet::new();
    for &sfi in seed_faces {
        let vc = mesh.faces[sfi].vertex_count as usize;
        for &vi in &mesh.faces[sfi].vertex_indices[..vc] {
            vertex_set.insert(vi);
        }
    }

    let mut face_list: Vec<usize> = seed_faces.to_vec();

    // Fit torus to seed faces
    let (mut current_center, mut current_dir, mut current_major_r, mut current_minor_r) = match
        fit_torus(&face_list, &vertex_set, &mesh.faces, &mesh.vertices)
    {
        Some(v) => v,
        None => {
            if verbosity >= 3 {
                eprintln!("  [BFS-torus] seed fit_torus returned None");
            }
            return None;
        }
    };

    if verbosity >= 3 {
        eprintln!("  [BFS-torus] fit result: R={:.4}, r={:.4}, center=[{:.4},{:.4},{:.4}], dir=[{:.4},{:.4},{:.4}]",
            current_major_r, current_minor_r,
            current_center[0], current_center[1], current_center[2],
            current_dir[0], current_dir[1], current_dir[2]);
    }

    // Verify seed: all vertices within tolerance (relaxed for seed)
    let seed_tol = surface_tol * 5.0;
    if !all_vertices_within_torus_tolerance(
        &vertex_set, &current_center, &current_dir, current_major_r, current_minor_r,
        seed_tol, &mesh.vertices,
    ) {
        if verbosity >= 3 {
            let worst_dist = vertex_set.iter().map(|&vi| {
                vertex_to_torus_distance(&mesh.vertices[vi], &current_center, &current_dir, current_major_r, current_minor_r).abs()
            }).fold(0.0_f64, f64::max);
            eprintln!("  [BFS-torus] seed vertex tolerance failed: worst={:.2e} > tol={:.2e}",
                worst_dist, seed_tol);
        }
        return None;
    }

    // Verify seed face centroids
    for &sfi in seed_faces {
        let centroid = face_centroid(&mesh.faces[sfi], &mesh.vertices);
        let cen_dist = vertex_to_torus_distance(
            &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
            &current_center, &current_dir, current_major_r, current_minor_r,
        ).abs();
        if cen_dist > seed_tol {
            if verbosity >= 3 {
                eprintln!("  [BFS-torus] seed centroid tolerance failed: face {} dist={:.2e}", sfi, cen_dist);
            }
            return None;
        }
    }

    // Determine convexity from first seed face
    let convex = determine_convexity(
        &mesh.faces[seed_faces[0]], &mesh.vertices,
        &current_center, &current_dir, current_major_r,
    );

    if verbosity >= 3 {
        eprintln!("  [BFS-torus] seed accept: {} faces, {}", seed_faces.len(),
            if convex { "convex" } else { "concave" });
    }

    // Viz: show seed faces
    let mut skip_viz = false;
    {
        let highlights = vec![
            viz::FaceHighlight { face_indices: seed_faces.to_vec(), color: [1.0, 0.6, 0.0, 1.0] },
        ];
        let _ = &highlights;
        if let Some(action) = viz_custom(
            viz, vec![
                viz::FaceHighlight { face_indices: seed_faces.to_vec(), color: [1.0, 0.6, 0.0, 1.0] },
            ],
            Vec::new(),
            &format!("BFS-torus: seed ({} faces) R={:.3} r={:.3} {} [space=step]",
                seed_faces.len(), current_major_r, current_minor_r,
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

    // Track claimed faces
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
            if cni < 0 { continue; }
            let cni = cni as usize;

            if mesh.faces[cni].toroidal_hypothesis != UNDEDUCED_TOROIDAL_HYPOTHESIS {
                continue;
            }
            if trial_claimed.contains(&cni) {
                continue;
            }

            // Angular tolerance check (use 2× like cones)
            let torus_angular_tol = angular_tol * 2.0;
            if let Some(n_cni) = mesh.faces[cni].normal {
                let cni_vc = mesh.faces[cni].vertex_count as usize;
                let cni_neighbors = mesh.faces[cni].neighbors;
                let mut angular_reject = false;
                for &adj in &cni_neighbors[..cni_vc] {
                    if adj < 0 { continue; }
                    let adj = adj as usize;
                    if !trial_claimed.contains(&adj) { continue; }
                    if let Some(n_adj) = mesh.faces[adj].normal {
                        let cos_a = dot3(&n_adj, &n_cni).clamp(-1.0, 1.0);
                        let angle = cos_a.acos();
                        if angle > torus_angular_tol {
                            angular_reject = true;
                            break;
                        }
                    }
                }
                if angular_reject {
                    continue;
                }
            }

            // Convexity consistency check
            let cni_convex = determine_convexity(
                &mesh.faces[cni], &mesh.vertices,
                &current_center, &current_dir, current_major_r,
            );
            if cni_convex != convex {
                continue;
            }

            // Vertex distance check
            let cni_vc = mesh.faces[cni].vertex_count as usize;
            let cni_vi = mesh.faces[cni].vertex_indices;
            let mut all_ok = true;
            let mut any_far = false;
            for &vi in &cni_vi[..cni_vc] {
                let d = vertex_to_torus_distance(
                    &mesh.vertices[vi],
                    &current_center, &current_dir, current_major_r, current_minor_r,
                ).abs();
                if d > surface_tol {
                    all_ok = false;
                    if d > REFIT_SKIP_MULTIPLIER * surface_tol {
                        any_far = true;
                        break;
                    }
                }
            }

            if any_far {
                continue;
            }

            // Centroid distance check
            let centroid = face_centroid(&mesh.faces[cni], &mesh.vertices);
            let centroid_dist = vertex_to_torus_distance(
                &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                &current_center, &current_dir, current_major_r, current_minor_r,
            ).abs();
            if centroid_dist > REFIT_SKIP_MULTIPLIER * surface_tol {
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

                let refit = fit_torus(
                    &trial_faces, &trial_vertices, &mesh.faces, &mesh.vertices,
                );
                let (new_center, new_dir, new_major_r, new_minor_r) = match refit {
                    Some(params) => params,
                    None => continue,
                };

                if !all_vertices_within_torus_tolerance(
                    &trial_vertices, &new_center, &new_dir, new_major_r, new_minor_r,
                    surface_tol, &mesh.vertices,
                ) {
                    continue;
                }

                // Check all centroids after refit
                let mut centroids_ok = true;
                for &f in &trial_faces {
                    let c = face_centroid(&mesh.faces[f], &mesh.vertices);
                    let d = vertex_to_torus_distance(
                        &MeshVertex::from_xyz(c[0], c[1], c[2]),
                        &new_center, &new_dir, new_major_r, new_minor_r,
                    ).abs();
                    if d > surface_tol {
                        centroids_ok = false;
                        break;
                    }
                }
                if !centroids_ok {
                    continue;
                }

                current_center = new_center;
                current_dir = new_dir;
                current_major_r = new_major_r;
                current_minor_r = new_minor_r;
            }

            // Accept this face
            trial_claimed.insert(cni);
            face_list.push(cni);
            for &vi in &cni_vi[..cni_vc] {
                vertex_set.insert(vi);
            }
            queue.push_back(cni);

            // Viz
            if !skip_viz {
                let accepted_nonseed: Vec<usize> = face_list.iter()
                    .filter(|f| !seed_faces.contains(f) && **f != cni)
                    .copied().collect();
                let mut highlights = vec![
                    viz::FaceHighlight { face_indices: seed_faces.to_vec(), color: [1.0, 0.6, 0.0, 1.0] },
                ];
                if !accepted_nonseed.is_empty() {
                    highlights.push(viz::FaceHighlight { face_indices: accepted_nonseed, color: [0.2, 0.4, 1.0, 1.0] });
                }
                highlights.push(viz::FaceHighlight { face_indices: vec![cni], color: [0.1, 0.2, 0.7, 1.0] });
                let bg_faces: Vec<usize> = (0..mesh.faces.len())
                    .filter(|f| mesh.faces[*f].toroidal_hypothesis >= 0 && !face_list.contains(f))
                    .collect();
                if !bg_faces.is_empty() {
                    highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                }
                if let Some(action) = viz_custom(
                    viz, highlights,
                    Vec::new(),
                    &format!("BFS-torus: accepted fi={cni} ({} faces) R={:.3} r={:.3}",
                        face_list.len(), current_major_r, current_minor_r),
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
    if let Some((final_center, final_dir, final_major_r, final_minor_r)) =
        fit_torus(&face_list, &vertex_set, &mesh.faces, &mesh.vertices)
    {
        current_center = final_center;
        current_dir = final_dir;
        current_major_r = final_major_r;
        current_minor_r = final_minor_r;
    }

    // Compute error metrics
    let mut error_max = 0.0_f64;
    let mut error_abs_sum = 0.0_f64;
    for &vi in &vertex_set {
        let d = vertex_to_torus_distance(
            &mesh.vertices[vi], &current_center, &current_dir, current_major_r, current_minor_r,
        ).abs();
        error_max = error_max.max(d);
        error_abs_sum += d;
    }

    // Minimum face count
    if face_list.len() < MIN_TORUS_FACES {
        if verbosity >= 3 {
            eprintln!("  [BFS-torus] rejected: {} faces < MIN_TORUS_FACES", face_list.len());
        }
        return None;
    }

    // Centroid check and compute centroid_error_max
    let mut centroid_error_max = 0.0_f64;
    for &f in &face_list {
        let c = face_centroid(&mesh.faces[f], &mesh.vertices);
        let d = vertex_to_torus_distance(
            &MeshVertex::from_xyz(c[0], c[1], c[2]),
            &current_center, &current_dir, current_major_r, current_minor_r,
        ).abs();
        centroid_error_max = centroid_error_max.max(d);
        if d > surface_tol {
            if verbosity >= 3 {
                eprintln!("  [BFS-torus] rejected: face {} centroid_err={:.2e} > tol={:.2e}", f, d, surface_tol);
            }
            return None;
        }
    }

    // Determine convexity from final fit using majority vote
    let mut convex_votes = 0_i32;
    for &f in &face_list {
        if determine_convexity(&mesh.faces[f], &mesh.vertices, &current_center, &current_dir, current_major_r) {
            convex_votes += 1;
        } else {
            convex_votes -= 1;
        }
    }
    let convex = convex_votes > 0;

    Some(TorusTrialResult {
        faces: face_list,
        vertices: vertex_set,
        center: current_center,
        axis_direction: current_dir,
        major_radius: current_major_r,
        minor_radius: current_minor_r,
        convex,
        error_max,
        centroid_error_max,
        error_abs_sum,
    })
}

// ---------------------------------------------------------------------------
// Stage 2.5: Toroidal hypothesis deduction
// ---------------------------------------------------------------------------

/// Deduce toroidal hypotheses from the mesh.
///
/// Seeding strategy: vertex-neighborhood seeding. For each vertex with ≥3
/// incident unclaimed faces, use those faces as a local seed for BFS trial.
/// This naturally separates torus regions from cylinder regions because
/// fit_torus fails on non-torus patches (cylinder faces produce tube centers
/// that are collinear, not circular, causing the circle fit to fail).
///
/// Cylinder hypotheses are deliberately ignored — torus fitting operates
/// independently of cylinder fitting.
pub(super) fn deduce_toroidal_hypotheses(
    mesh: &mut ConnectedMesh,
    planar_hypotheses: &[super::PlanarHypothesis],
    vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
    verbosity: u8,
    viz: Option<&VizSender>,
) -> (Vec<ToroidalHypothesis>, bool) {
    let num_faces = mesh.faces.len();
    let num_vertices = mesh.vertices.len();
    let mut hypotheses: Vec<ToroidalHypothesis> = Vec::new();
    let viz_quit = std::cell::Cell::new(false);
    let _ = vertex_tol;

    // Identify faces committed to multi-face planar, spherical, or conical
    // hypotheses (ineligible for torus fitting).
    // NOTE: Cylindrical hypotheses are deliberately NOT excluded — torus
    // fitting ignores cylinder fits entirely.
    let mut committed = vec![false; num_faces];
    for (fi, face) in mesh.faces.iter().enumerate() {
        // Multi-face planar
        if face.planar_hypothesis >= 0 {
            let pi = face.planar_hypothesis as usize;
            if pi < planar_hypotheses.len() && planar_hypotheses[pi].faces.len() > 1 {
                committed[fi] = true;
            }
        }
        // Spherical
        if face.spherical_hypothesis >= 0 {
            committed[fi] = true;
        }
        // Conical
        if face.conical_hypothesis >= 0 {
            committed[fi] = true;
        }
    }

    // Build vertex → face adjacency for uncommitted faces with normals
    let mut vertex_faces: Vec<Vec<usize>> = vec![Vec::new(); num_vertices];
    for (fi, face) in mesh.faces.iter().enumerate() {
        if committed[fi] || face.normal.is_none() {
            continue;
        }
        let vc = face.vertex_count as usize;
        for &vi in &face.vertex_indices[..vc] {
            vertex_faces[vi].push(fi);
        }
    }

    // Collect candidate seed vertices, sorted by incident face count (descending).
    // Vertices with more incident faces produce better initial fits.
    let mut seed_vertices: Vec<(usize, usize)> = vertex_faces
        .iter()
        .enumerate()
        .filter(|(_, faces)| faces.len() >= 3)
        .map(|(vi, faces)| (vi, faces.len()))
        .collect();
    seed_vertices.sort_by(|a, b| b.1.cmp(&a.1));

    if verbosity >= 2 {
        let uncommitted_faces: usize = committed.iter().filter(|&&c| !c).count();
        eprintln!(
            "  [torus] {} uncommitted faces, {} candidate seed vertices",
            uncommitted_faces,
            seed_vertices.len(),
        );
    }

    for (vi, _face_count) in &seed_vertices {
        if viz_quit.get() {
            break;
        }

        // Collect still-undeduced faces incident to this vertex
        let seed_faces: Vec<usize> = vertex_faces[*vi]
            .iter()
            .filter(|&&fi| mesh.faces[fi].toroidal_hypothesis == UNDEDUCED_TOROIDAL_HYPOTHESIS)
            .copied()
            .collect();
        if seed_faces.len() < 3 {
            continue;
        }

        if verbosity >= 3 {
            eprintln!("  [torus] trying vertex {} seed with {} faces", vi, seed_faces.len());
        }

        // Run BFS trial from this vertex's neighborhood
        let trial = run_torus_trial_bfs(
            &seed_faces,
            mesh,
            surface_tol,
            angular_tol,
            verbosity,
            viz,
            &viz_quit,
        );

        if viz_quit.get() {
            return (hypotheses, true);
        }

        // Validate angular coverage
        let trial = trial.and_then(|t| {
            if angular_coverage_valid(
                &t.faces,
                &t.center,
                &t.axis_direction,
                &mesh.faces,
                &mesh.vertices,
            ) {
                Some(t)
            } else {
                if verbosity >= 3 {
                    eprintln!("  [torus] angular coverage failed: {} faces", t.faces.len());
                }
                None
            }
        });

        // Validate error and radii
        let trial = trial.and_then(|t| {
            if t.error_max > surface_tol {
                if verbosity >= 3 {
                    eprintln!(
                        "  [torus] rejected: error_max={:.2e} > tol={:.2e}",
                        t.error_max, surface_tol
                    );
                }
                return None;
            }
            if t.minor_radius <= 0.0 || t.major_radius <= 0.0 {
                if verbosity >= 3 {
                    eprintln!(
                        "  [torus] rejected: degenerate radii R={:.4} r={:.4}",
                        t.major_radius, t.minor_radius
                    );
                }
                return None;
            }
            // Reject if minor_radius > major_radius (self-intersecting torus)
            if t.minor_radius > t.major_radius {
                if verbosity >= 3 {
                    eprintln!(
                        "  [torus] rejected: r={:.4} > R={:.4}",
                        t.minor_radius, t.major_radius
                    );
                }
                return None;
            }
            Some(t)
        });

        if let Some(candidate) = trial {
            if verbosity >= 2 {
                eprintln!(
                    "  [torus] ACCEPTED: {} faces, R={:.4}, r={:.4}, {}, err_max={:.2e}",
                    candidate.faces.len(),
                    candidate.major_radius,
                    candidate.minor_radius,
                    if candidate.convex { "convex" } else { "concave" },
                    candidate.error_max
                );
            }

            // Viz: show accepted hypothesis
            if !viz_quit.get() {
                if let Some(action) = viz_custom(
                    viz,
                    vec![viz::FaceHighlight {
                        face_indices: candidate.faces.clone(),
                        color: [0.0, 0.8, 0.0, 1.0],
                    }],
                    Vec::new(),
                    &format!(
                        "BFS-torus ACCEPTED: {} faces R={:.3} r={:.3} {} [space=next]",
                        candidate.faces.len(),
                        candidate.major_radius,
                        candidate.minor_radius,
                        if candidate.convex { "convex" } else { "concave" }
                    ),
                    Vec::new(),
                    Vec::new(),
                    Some(viz_face_centroid(candidate.faces[0], mesh)),
                    viz_face_normal(candidate.faces[0], mesh),
                ) {
                    match action {
                        VizAction::Quit => {
                            return (hypotheses, true);
                        }
                        VizAction::NextSeed | VizAction::NextStep => {}
                    }
                }
            }

            // Commit
            let hi = hypotheses.len() as i32;
            for &f in &candidate.faces {
                mesh.faces[f].toroidal_hypothesis = hi;
            }

            hypotheses.push(ToroidalHypothesis {
                center: candidate.center,
                axis_direction: candidate.axis_direction,
                major_radius: candidate.major_radius,
                minor_radius: candidate.minor_radius,
                convex: candidate.convex,
                faces: candidate.faces,
                vertices: candidate.vertices.into_iter().collect(),
                error_max: candidate.error_max,
                centroid_error_max: candidate.centroid_error_max,
                error_abs_sum: candidate.error_abs_sum,
            });
        }
    }

    // Merge hypotheses that describe the same physical torus
    // (vertex-neighborhood seeding can fragment a single torus into many small patches)
    loop {
        let mut merged_any = false;
        let n = hypotheses.len();
        'outer: for i in 0..n {
            if hypotheses[i].faces.is_empty() {
                continue;
            }
            for j in (i + 1)..n {
                if hypotheses[j].faces.is_empty() {
                    continue;
                }
                if !torus_params_compatible(&hypotheses[i], &hypotheses[j]) {
                    continue;
                }

                if verbosity >= 3 {
                    eprintln!(
                        "  [torus] merging hypothesis {} ({} faces) into {} ({} faces)",
                        j, hypotheses[j].faces.len(), i, hypotheses[i].faces.len()
                    );
                }

                // Keep parameters from the larger hypothesis (better fit)
                let (keep, absorb) = if hypotheses[i].faces.len() >= hypotheses[j].faces.len() {
                    (i, j)
                } else {
                    (j, i)
                };
                let absorb_faces = hypotheses[absorb].faces.clone();
                let absorb_verts = hypotheses[absorb].vertices.clone();
                hypotheses[keep].faces.extend(absorb_faces.iter());
                hypotheses[keep].vertices.extend(absorb_verts.iter());
                // Deduplicate vertices
                hypotheses[keep].vertices.sort();
                hypotheses[keep].vertices.dedup();

                // Recompute error metrics for merged hypothesis
                let mut err_max = 0.0_f64;
                let mut err_sum = 0.0_f64;
                for &vi in &hypotheses[keep].vertices {
                    let d = vertex_to_torus_distance(
                        &mesh.vertices[vi],
                        &hypotheses[keep].center,
                        &hypotheses[keep].axis_direction,
                        hypotheses[keep].major_radius,
                        hypotheses[keep].minor_radius,
                    )
                    .abs();
                    err_max = err_max.max(d);
                    err_sum += d;
                }
                let mut cen_err_max = 0.0_f64;
                for &f in &hypotheses[keep].faces {
                    let c = face_centroid(&mesh.faces[f], &mesh.vertices);
                    let d = vertex_to_torus_distance(
                        &MeshVertex::from_xyz(c[0], c[1], c[2]),
                        &hypotheses[keep].center,
                        &hypotheses[keep].axis_direction,
                        hypotheses[keep].major_radius,
                        hypotheses[keep].minor_radius,
                    )
                    .abs();
                    cen_err_max = cen_err_max.max(d);
                }
                hypotheses[keep].error_max = err_max;
                hypotheses[keep].error_abs_sum = err_sum;
                hypotheses[keep].centroid_error_max = cen_err_max;

                // Clear the absorbed hypothesis
                hypotheses[absorb].faces.clear();
                hypotheses[absorb].vertices.clear();

                merged_any = true;
                break 'outer;
            }
        }
        if !merged_any {
            break;
        }
    }

    // Remove empty (merged-away) hypotheses and renumber
    let old_len = hypotheses.len();
    let mut idx_map: Vec<i32> = vec![-1; old_len];
    let mut new_hypotheses: Vec<ToroidalHypothesis> = Vec::new();
    for (old_idx, h) in hypotheses.into_iter().enumerate() {
        if !h.faces.is_empty() {
            idx_map[old_idx] = new_hypotheses.len() as i32;
            new_hypotheses.push(h);
        }
    }
    for fi in 0..num_faces {
        let old = mesh.faces[fi].toroidal_hypothesis;
        if old >= 0 && (old as usize) < old_len {
            mesh.faces[fi].toroidal_hypothesis = idx_map[old as usize];
        }
    }
    hypotheses = new_hypotheses;

    if verbosity >= 2 && hypotheses.len() < old_len {
        eprintln!(
            "  [torus] merged {} → {} hypotheses",
            old_len,
            hypotheses.len()
        );
    }

    // Set remaining UNDEDUCED faces to NO_HYPOTHESIS
    for fi in 0..num_faces {
        if mesh.faces[fi].toroidal_hypothesis == UNDEDUCED_TOROIDAL_HYPOTHESIS {
            mesh.faces[fi].toroidal_hypothesis = NO_HYPOTHESIS;
        }
    }

    (hypotheses, false)
}

/// Check if two torus hypotheses describe the same physical torus.
fn torus_params_compatible(a: &ToroidalHypothesis, b: &ToroidalHypothesis) -> bool {
    // Radii must match within 1% relative (or 0.01 absolute for small radii)
    let r_tol = ((a.major_radius + b.major_radius) * 0.005).max(0.01);
    let r_small_tol = ((a.minor_radius + b.minor_radius) * 0.005).max(0.01);
    if (a.major_radius - b.major_radius).abs() > r_tol {
        return false;
    }
    if (a.minor_radius - b.minor_radius).abs() > r_small_tol {
        return false;
    }
    // Centers must be close (within 10% of minor radius)
    let dc = [
        a.center[0] - b.center[0],
        a.center[1] - b.center[1],
        a.center[2] - b.center[2],
    ];
    let center_dist = (dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2]).sqrt();
    if center_dist > a.minor_radius * 0.1 {
        return false;
    }
    // Axes must be nearly parallel (or anti-parallel)
    let cos = dot3(&a.axis_direction, &b.axis_direction).abs();
    if cos < 0.999 {
        return false;
    }
    true
}

// ---------------------------------------------------------------------------
// Stage 2.5: Compare check
// ---------------------------------------------------------------------------

/// Validate fitted toroidal hypotheses against a reference STEP shape.
pub(super) fn compare_toroidal_hypotheses(
    hypotheses: &[ToroidalHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);

            // Project centroid onto torus surface
            let d = [
                centroid[0] - hyp.center[0],
                centroid[1] - hyp.center[1],
                centroid[2] - hyp.center[2],
            ];
            let axial = dot3(&d, &hyp.axis_direction);
            let radial = [
                d[0] - axial * hyp.axis_direction[0],
                d[1] - axial * hyp.axis_direction[1],
                d[2] - axial * hyp.axis_direction[2],
            ];
            let radial_dist = (radial[0] * radial[0] + radial[1] * radial[1]
                + radial[2] * radial[2]).sqrt();

            let projected = if radial_dist > 1e-15 {
                let tube_center = [
                    hyp.center[0] + hyp.major_radius * radial[0] / radial_dist,
                    hyp.center[1] + hyp.major_radius * radial[1] / radial_dist,
                    hyp.center[2] + hyp.major_radius * radial[2] / radial_dist,
                ];
                let tube_vec = [
                    centroid[0] - tube_center[0],
                    centroid[1] - tube_center[1],
                    centroid[2] - tube_center[2],
                ];
                let tube_dist = (tube_vec[0] * tube_vec[0] + tube_vec[1] * tube_vec[1]
                    + tube_vec[2] * tube_vec[2]).sqrt();
                if tube_dist > 1e-15 {
                    let scale = hyp.minor_radius / tube_dist;
                    [
                        tube_center[0] + tube_vec[0] * scale,
                        tube_center[1] + tube_vec[1] * scale,
                        tube_center[2] + tube_vec[2] * scale,
                    ]
                } else {
                    centroid
                }
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
                hypothesis_type: "toroidal",
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance,
            });
        }
    }

    Ok(())
}
