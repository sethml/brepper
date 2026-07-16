//! Stage 2.5: Toroidal hypothesis deduction and comparison.
//!
//! Uses the medial axis / tube center method:
//! 1. Estimate minor radius from normal-line intersections of adjacent face pairs
//! 2. Compute tube centers k_i = p_i ± r*n_i
//! 3. PCA of tube centers → axis direction
//! 4. 2D algebraic circle fit → center and major radius R
//! 5. LM refinement of 7 parameters [Cx, Cy, Cz, alpha, beta, R, r]
//! 6. BFS expansion with vertex-to-torus distance validation

use std::collections::HashSet;

use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{Const, Dyn, OMatrix, OVector, Owned};

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_TOROIDAL_HYPOTHESIS};
use crate::viz;

use super::{
    dot3, face_centroid,
    normalize3, perpendicular_basis, smallest_eigenvector_3x3,
    ToroidalHypothesis, Stage2CompareError,
    MIN_TORUS_FACES,
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
        let mut jac = OMatrix::<f64, Dyn, Const<7>>::zeros_generic(Dyn(n), Const::<7>);
        let a = self.axis_dir_from_params(&self.params);
        let center = [self.params[2], self.params[3], self.params[4]];
        let major_r = self.params[5];
        let _minor_r = self.params[6];

        // L = |a0 + α*u0 + β*w0| (unnormalized axis length)
        let alpha = self.params[0];
        let beta = self.params[1];
        let vx = self.a0[0] + alpha * self.u0[0] + beta * self.w0[0];
        let vy = self.a0[1] + alpha * self.u0[1] + beta * self.w0[1];
        let vz = self.a0[2] + alpha * self.u0[2] + beta * self.w0[2];
        let big_l = (vx * vx + vy * vy + vz * vz).sqrt().max(1e-15);

        for (i, p) in self.points.iter().enumerate() {
            let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
            let h_i = dot3(&d, &a);
            let lat = [d[0] - h_i * a[0], d[1] - h_i * a[1], d[2] - h_i * a[2]];
            let rho_i = (lat[0] * lat[0] + lat[1] * lat[1] + lat[2] * lat[2]).sqrt();
            let dr = rho_i - major_r;
            let tube_dist = (dr * dr + h_i * h_i).sqrt();

            if tube_dist < 1e-12 {
                jac.row_mut(i).fill(0.0);
                jac[(i, 6)] = -1.0; // ∂f/∂r = -1 always
                continue;
            }

            let s_r = dr / tube_dist;
            let s_h = h_i / tube_dist;

            let rho_hat = if rho_i > 1e-12 {
                [lat[0] / rho_i, lat[1] / rho_i, lat[2] / rho_i]
            } else {
                [0.0; 3]
            };

            // ∂f/∂center_k = -(s_r * rho_hat[k] + s_h * a[k])
            jac[(i, 2)] = -(s_r * rho_hat[0] + s_h * a[0]);
            jac[(i, 3)] = -(s_r * rho_hat[1] + s_h * a[1]);
            jac[(i, 4)] = -(s_r * rho_hat[2] + s_h * a[2]);

            // ∂f/∂R = -s_r
            jac[(i, 5)] = -s_r;

            // ∂f/∂r = -1
            jac[(i, 6)] = -1.0;

            // ∂f/∂α, ∂f/∂β: axis direction perturbations
            // ∂h/∂α = lat·u0/L, ∂ρ/∂α = -h_i*lat·u0/(ρ*L)
            // ∂f/∂α = s_r * ∂ρ/∂α + s_h * ∂h/∂α
            let lat_u = dot3(&lat, &self.u0);
            let lat_w = dot3(&lat, &self.w0);
            let (drho_du, drho_dw) = if rho_i > 1e-12 {
                (-h_i * lat_u / (rho_i * big_l), -h_i * lat_w / (rho_i * big_l))
            } else {
                (0.0, 0.0)
            };
            let dh_du = lat_u / big_l;
            let dh_dw = lat_w / big_l;
            jac[(i, 0)] = s_r * drho_du + s_h * dh_du;
            jac[(i, 1)] = s_r * drho_dw + s_h * dh_dw;
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
/// A torus fit result: (center, axis_direction, major_radius, minor_radius).
type TorusFitTuple = ([f64; 3], [f64; 3], f64, f64);

/// Algebraic (Kasa) circle fit to 2D points. Coordinates are mean-centered
/// internally for conditioning. Returns ((cx, cy), radius).
fn fit_circle_2d(coords: &[(f64, f64)]) -> Option<((f64, f64), f64)> {
    let n = coords.len();
    if n < 3 {
        return None;
    }
    let nf = n as f64;
    let mx = coords.iter().map(|c| c.0).sum::<f64>() / nf;
    let my = coords.iter().map(|c| c.1).sum::<f64>() / nf;

    let mut sx = 0.0_f64;
    let mut sy = 0.0_f64;
    let mut sxx = 0.0_f64;
    let mut syy = 0.0_f64;
    let mut sxy = 0.0_f64;
    let mut sxr = 0.0_f64;
    let mut syr = 0.0_f64;
    let mut sr = 0.0_f64;
    for &(x0, y0) in coords {
        let x = x0 - mx;
        let y = y0 - my;
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

    let a11 = sxx; let a12 = sxy; let a13 = sx;
    let a21 = sxy; let a22 = syy; let a23 = sy;
    let a31 = sx;  let a32 = sy;  let a33 = nf;
    let b1 = -sxr; let b2 = -syr; let b3 = -sr;

    let det = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31)
        + a13 * (a21 * a32 - a22 * a31);
    if det.abs() < 1e-30 {
        return None;
    }
    let inv_det = 1.0 / det;
    let d = inv_det * (b1 * (a22 * a33 - a23 * a32) - a12 * (b2 * a33 - a23 * b3) + a13 * (b2 * a32 - a22 * b3));
    let e = inv_det * (a11 * (b2 * a33 - a23 * b3) - b1 * (a21 * a33 - a23 * a31) + a13 * (a21 * b3 - b2 * a31));
    let f = inv_det * (a11 * (a22 * b3 - b2 * a32) - a12 * (a21 * b3 - b2 * a31) + b1 * (a21 * a32 - a22 * a31));

    let r_sq = d * d / 4.0 + e * e / 4.0 - f;
    if r_sq <= 0.0 {
        return None;
    }
    Some(((mx - d / 2.0, my - e / 2.0), r_sq.sqrt()))
}

/// Third fitting method: multi-axis rescue for small torus patches.
///
/// Small patches (e.g. corner vertex blends spanning ~90° of a small major
/// circle) defeat both the tube-center method (too few tube-direction normal
/// pairs) and position-PCA axis estimation (a compact patch is a blob, not a
/// ring, and its centroid is far off the axis). For each candidate axis
/// direction: locate the axis position with a 2D circle fit of the vertices
/// projected along the axis, fit the tube profile circle in (h, r)
/// coordinates about that axis, and polish with LM. Returns the best result
/// by max vertex distance.
fn fit_torus_axes_rescue(
    face_indices: &[usize],
    vertex_set: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64, f64)> {
    if vertex_set.len() < 7 {
        return None;
    }

    // Candidate axis directions: coordinate axes plus the normal-covariance
    // and position-PCA smallest eigenvectors.
    let mut candidates: Vec<[f64; 3]> =
        vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    {
        let mut ncov = [0.0_f64; 6];
        for &fi in face_indices {
            if let Some(n) = faces[fi].normal {
                ncov[0] += n[0] * n[0];
                ncov[1] += n[0] * n[1];
                ncov[2] += n[0] * n[2];
                ncov[3] += n[1] * n[1];
                ncov[4] += n[1] * n[2];
                ncov[5] += n[2] * n[2];
            }
        }
        let mut a = smallest_eigenvector_3x3(&ncov);
        if normalize3(&mut a) > 1e-12 {
            candidates.push(a);
        }
    }
    {
        let n = vertex_set.len() as f64;
        let mut c = [0.0_f64; 3];
        for &vi in vertex_set {
            c[0] += vertices[vi].x;
            c[1] += vertices[vi].y;
            c[2] += vertices[vi].z;
        }
        c[0] /= n; c[1] /= n; c[2] /= n;
        let mut pcov = [0.0_f64; 6];
        for &vi in vertex_set {
            let dx = vertices[vi].x - c[0];
            let dy = vertices[vi].y - c[1];
            let dz = vertices[vi].z - c[2];
            pcov[0] += dx * dx;
            pcov[1] += dx * dy;
            pcov[2] += dx * dz;
            pcov[3] += dy * dy;
            pcov[4] += dy * dz;
            pcov[5] += dz * dz;
        }
        let mut a = smallest_eigenvector_3x3(&pcov);
        if normalize3(&mut a) > 1e-12 {
            candidates.push(a);
        }
    }

    let points: Vec<[f64; 3]> = vertex_set
        .iter()
        .map(|&vi| [vertices[vi].x, vertices[vi].y, vertices[vi].z])
        .collect();

    let mut best: Option<(f64, TorusFitTuple)> = None;
    for axis in candidates {
        let (u, w) = perpendicular_basis(&axis);
        // Axis position from a 2D circle fit of the projected vertices.
        let coords: Vec<(f64, f64)> = points
            .iter()
            .map(|p| (dot3(p, &u), dot3(p, &w)))
            .collect();
        let Some(((a2, b2), _)) = fit_circle_2d(&coords) else { continue };
        let origin = [
            a2 * u[0] + b2 * w[0],
            a2 * u[1] + b2 * w[1],
            a2 * u[2] + b2 * w[2],
        ];
        // Tube profile circle in (h, r) coordinates about the axis.
        let hr: Vec<(f64, f64)> = points
            .iter()
            .map(|p| {
                let d = [p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]];
                let h = dot3(&d, &axis);
                let rx = d[0] - h * axis[0];
                let ry = d[1] - h * axis[1];
                let rz = d[2] - h * axis[2];
                (h, (rx * rx + ry * ry + rz * rz).sqrt())
            })
            .collect();
        let Some(((h0, major_r), minor_r)) = fit_circle_2d(&hr) else { continue };
        if major_r <= 1e-6 || minor_r <= 1e-6 || minor_r >= major_r {
            continue;
        }
        // Cheap gate before the LM: if the (h,r) profile points sit far off
        // the fitted profile circle, this axis candidate is hopeless (mixed
        // or non-toroidal seed) and LM will not save it. This keeps the
        // rescue nearly free for the many junk seeds it is offered.
        let profile_rms = {
            let mut ss = 0.0_f64;
            for &(h, r) in &hr {
                let d = ((h - h0) * (h - h0) + (r - major_r) * (r - major_r)).sqrt() - minor_r;
                ss += d * d;
            }
            (ss / hr.len() as f64).sqrt()
        };
        if profile_rms > minor_r * 0.2 {
            continue;
        }
        let center = [
            origin[0] + h0 * axis[0],
            origin[1] + h0 * axis[1],
            origin[2] + h0 * axis[2],
        ];

        // LM polish (tightly bounded patience: this runs for several
        // candidates per seed and again per pruning round; the Kasa init is
        // close for real tori, so converging fits converge quickly).
        let problem = TorusLMProblem::new(points.clone(), axis, center, major_r, minor_r);
        let (result, _report) = LevenbergMarquardt::new().with_patience(75).minimize(problem);
        let refined_axis = result.axis_dir_from_params(&result.params);
        let refined_center = [result.params[2], result.params[3], result.params[4]];
        let refined_major = result.params[5];
        let refined_minor = result.params[6];
        let candidate = if refined_major > 0.0
            && refined_minor > 0.0
            && refined_minor < refined_major
        {
            (refined_center, refined_axis, refined_major, refined_minor)
        } else {
            (center, axis, major_r, minor_r)
        };
        let err = vertex_set
            .iter()
            .map(|&vi| {
                vertex_to_torus_distance(
                    &vertices[vi], &candidate.0, &candidate.1, candidate.2, candidate.3,
                )
                .abs()
            })
            .fold(0.0_f64, f64::max);
        if best.as_ref().is_none_or(|(e, _)| err < *e) {
            best = Some((err, candidate));
        }
        // First excellent candidate wins — no need to try the rest.
        if err < 1e-5 {
            break;
        }
    }
    best.map(|(_, c)| c)
}

fn fit_torus(
    face_indices: &[usize],
    vertex_set: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64, f64)> {
    let max_err = |fit: &([f64; 3], [f64; 3], f64, f64)| -> f64 {
        vertex_set.iter().map(|&vi| {
            vertex_to_torus_distance(&vertices[vi], &fit.0, &fit.1, fit.2, fit.3).abs()
        }).fold(0.0_f64, f64::max)
    };
    // Structural validity must be checked BEFORE picking the lowest-error
    // candidate: on small corner patches the tube-center method can return
    // a self-intersecting fit (minor > major) with a locally lower vertex
    // error, which would shadow a valid rescue fit and then be rejected
    // downstream.
    let valid = |fit: &([f64; 3], [f64; 3], f64, f64)| -> bool {
        fit.2 > 0.0 && fit.3 > 0.0 && fit.3 <= fit.2
    };

    // Try the tube-center and SOR methods and keep the better fit.
    let inner = fit_torus_inner(face_indices, vertex_set, faces, vertices);
    let sor = fit_torus_sor(face_indices, vertex_set, faces, vertices);
    let mut best: Option<(f64, TorusFitTuple)> = None;
    for cand in [inner, sor].into_iter().flatten() {
        if !valid(&cand) {
            continue;
        }
        let err = max_err(&cand);
        if best.as_ref().is_none_or(|(e, _)| err < *e) {
            best = Some((err, cand));
        }
    }

    // Multi-axis rescue for small patches, only when the geometric methods
    // fail to reach CAD-mesh precision (bounds the extra LM cost to the
    // rare hard cases).
    const RESCUE_ERR_THRESHOLD: f64 = 1e-4; // mm
    if best.as_ref().is_none_or(|(e, _)| *e > RESCUE_ERR_THRESHOLD) {
        if let Some(cand) = fit_torus_axes_rescue(face_indices, vertex_set, faces, vertices) {
            if valid(&cand) {
                let err = max_err(&cand);
                if best.as_ref().is_none_or(|(e, _)| err < *e) {
                    best = Some((err, cand));
                }
            }
        }
    }

    best.map(|(_, c)| c)
}

fn fit_torus_inner(
    face_indices: &[usize],
    vertex_set: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64, f64)> {
    if vertex_set.len() < 7 || face_indices.len() < 3 {
        return None;
    }

    // Sort vertices for deterministic iteration order (avoids floating-point
    // accumulation-order differences across runs due to HashSet randomization).
    let sorted_verts: Vec<usize> = {
        let mut v: Vec<usize> = vertex_set.to_vec();
        v.sort_unstable();
        v
    };

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
        return None;
    }

    // Robust minor radius: 25th percentile of estimates.
    // On a torus, adjacent-face normal-line intersections measure a mix of:
    //  - tube curvature (≈ r, the minor radius) — from edges along the tube
    //  - azimuthal curvature (≈ R ± r) — from edges along the major circle
    // The tube-curvature estimates are the SMALLEST values, so we use a low
    // percentile instead of median to avoid contamination from azimuthal estimates.
    r_estimates.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let minor_r = r_estimates[r_estimates.len() / 4];

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
        let tube_centers: Vec<[f64; 3]> = sorted_verts.iter().filter_map(|&vi| {
            // Average normal from all incident seed faces for a better
            // surface-normal estimate at the vertex.
            let mut avg_n = [0.0_f64; 3];
            let mut count = 0u32;
            for &fi in face_indices {
                let vc = faces[fi].vertex_count as usize;
                if !faces[fi].vertex_indices[..vc].contains(&vi) { continue; }
                if let Some(n) = faces[fi].normal {
                    avg_n[0] += n[0];
                    avg_n[1] += n[1];
                    avg_n[2] += n[2];
                    count += 1;
                }
            }
            if count == 0 { return None; }
            let len = normalize3(&mut avg_n);
            if len < 1e-15 { return None; }
            let v = &vertices[vi];
            Some([
                v.x + sign * minor_r * avg_n[0],
                v.y + sign * minor_r * avg_n[1],
                v.z + sign * minor_r * avg_n[2],
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
    }).min_by(|a, b| a.3.partial_cmp(&b.3).unwrap());

    let best = match best {
        Some(b) => b,
        None => {
            return None;
        }
    };

    let (center, axis, major_r, rms, _sign) = best;

    // Quality gate: reject if tube centers don't form a clean circle.
    // This catches mixed convex+concave seeds where tube centers scatter
    // between two different major circles. Use tight threshold since valid
    // torus fits have tube centers very close to the major circle.
    if rms > minor_r * 0.1 {
        return None;
    }

    // Step 4: LM refinement
    let points: Vec<[f64; 3]> = sorted_verts.iter()
        .map(|&vi| [vertices[vi].x, vertices[vi].y, vertices[vi].z])
        .collect();

    if points.len() >= 7 {
        let problem = TorusLMProblem::new(points, axis, center, major_r, minor_r);
        let (result, _report) = LevenbergMarquardt::new().with_patience(150).minimize(problem);
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
    vertex_set: &[usize],
    _faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64, f64)> {
    if vertex_set.len() < 7 || face_indices.len() < 3 {
        return None;
    }

    // Sort vertices for deterministic iteration order.
    let sorted_verts: Vec<usize> = {
        let mut v: Vec<usize> = vertex_set.to_vec();
        v.sort_unstable();
        v
    };

    // Estimate axis from PCA of vertex positions (smallest eigenvector).
    // For a torus, axis variance ~ r²/2 while perpendicular variance ~
    // (R² + r²/2)/2 >> r²/2 for any non-self-intersecting torus (R > r).
    // Normal covariance doesn't work well here since torus normals have
    // significant axis components (unlike cylinders/cones).
    let n = sorted_verts.len();
    let inv_n = 1.0 / n as f64;
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    for &vi in &sorted_verts {
        cx += vertices[vi].x;
        cy += vertices[vi].y;
        cz += vertices[vi].z;
    }
    cx *= inv_n;
    cy *= inv_n;
    cz *= inv_n;

    // Build position covariance matrix (symmetric, packed as [xx, xy, xz, yy, yz, zz])
    let mut pos_cov = [0.0_f64; 6];
    for &vi in &sorted_verts {
        let dx = vertices[vi].x - cx;
        let dy = vertices[vi].y - cy;
        let dz = vertices[vi].z - cz;
        pos_cov[0] += dx * dx;
        pos_cov[1] += dx * dy;
        pos_cov[2] += dx * dz;
        pos_cov[3] += dy * dy;
        pos_cov[4] += dy * dz;
        pos_cov[5] += dz * dz;
    }
    let mut axis = smallest_eigenvector_3x3(&pos_cov);
    let len = normalize3(&mut axis);
    if len < 1e-15 {
        return None;
    }

    // Compute (h, r) coordinates relative to centroid
    let mut hrs: Vec<(f64, f64)> = Vec::with_capacity(n);
    for &vi in &sorted_verts {
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

    // Center is at the axis lateral position found by Kasa circle fit
    // on the 2D perpendicular projections of vertices.
    let (u_basis, w_basis) = perpendicular_basis(&axis);
    let mut sum_px = 0.0_f64;
    let mut sum_py = 0.0_f64;
    let mut sum_px2 = 0.0_f64;
    let mut sum_py2 = 0.0_f64;
    let mut sum_pxpy = 0.0_f64;
    let mut sum_pxs = 0.0_f64;
    let mut sum_pys = 0.0_f64;
    let mut sum_ss = 0.0_f64;
    for &vi in &sorted_verts {
        let dx = vertices[vi].x - cx;
        let dy = vertices[vi].y - cy;
        let dz = vertices[vi].z - cz;
        let px = dot3(&[dx, dy, dz], &u_basis);
        let py = dot3(&[dx, dy, dz], &w_basis);
        let s_val = px * px + py * py;
        sum_px += px;
        sum_py += py;
        sum_px2 += px * px;
        sum_py2 += py * py;
        sum_pxpy += px * py;
        sum_pxs += px * s_val;
        sum_pys += py * s_val;
        sum_ss += s_val;
    }
    // Kasa circle fit: solve [2px 2py 1] * [cx_2d; cy_2d; c] = px^2+py^2
    // Normal equations: A^T A x = A^T b
    let a00 = 4.0 * sum_px2;
    let a01 = 4.0 * sum_pxpy;
    let a02 = 2.0 * sum_px;
    let a11 = 4.0 * sum_py2;
    let a12 = 2.0 * sum_py;
    let a22 = nf;
    let b0 = 2.0 * sum_pxs;
    let b1 = 2.0 * sum_pys;
    let b2 = sum_ss;
    // Solve 3x3 system via Cramer's rule
    let det = a00 * (a11 * a22 - a12 * a12)
            - a01 * (a01 * a22 - a12 * a02)
            + a02 * (a01 * a12 - a11 * a02);
    let (center_2d_x, center_2d_y) = if det.abs() > 1e-30 {
        let inv_det = 1.0 / det;
        let kx = (b0 * (a11 * a22 - a12 * a12)
                - a01 * (b1 * a22 - a12 * b2)
                + a02 * (b1 * a12 - a11 * b2)) * inv_det;
        let ky = (a00 * (b1 * a22 - a12 * b2)
                - b0 * (a01 * a22 - a12 * a02)
                + a02 * (a01 * b2 - b1 * a02)) * inv_det;
        (kx, ky)
    } else {
        (0.0, 0.0)
    };
    let center = [
        cx + center_2d_x * u_basis[0] + center_2d_y * w_basis[0],
        cy + center_2d_x * u_basis[1] + center_2d_y * w_basis[1],
        cz + center_2d_x * u_basis[2] + center_2d_y * w_basis[2],
    ];

    // LM refinement
    let points: Vec<[f64; 3]> = sorted_verts.iter()
        .map(|&vi| [vertices[vi].x, vertices[vi].y, vertices[vi].z])
        .collect();

    if points.len() >= 7 && major_r > 1e-10 {
        let problem = TorusLMProblem::new(points, axis, center, major_r, minor_r);
        let (result, _report) = LevenbergMarquardt::new().with_patience(150).minimize(problem);
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

/// Validate angular coverage around the torus minor circle (tube cross-section).
/// For each face centroid, compute the angle φ = atan2(axial, radial - R) in
/// the tube cross-section plane. Reject if the angular span is too narrow.
/// This catches false positives where faces lie on a strip of the torus
/// (e.g., just the outer equator) without wrapping around the tube.
fn minor_circle_coverage_valid(
    face_list: &[usize],
    center: &[f64; 3],
    axis_direction: &[f64; 3],
    major_radius: f64,
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> bool {
    if face_list.len() < MIN_TORUS_FACES {
        return false;
    }

    // Measure tube coverage from face *normals*, not centroid positions.
    // On a torus, the surface normal at tube angle φ makes exactly that
    // angle with the radial direction, regardless of tube radius. Flat
    // annuli and cylinder bands lie exactly on degenerate tori (their
    // centroids can span a nonzero positional arc on a large tube), but
    // their normals are all parallel — tube-angle span ≈ 0 — so the
    // normal-based span cleanly rejects them.
    let _ = major_radius;
    let mut minor_angles: Vec<f64> = Vec::with_capacity(face_list.len());
    for &fi in face_list {
        let n = match faces[fi].normal {
            Some(n) => n,
            None => continue,
        };
        let c = face_centroid(&faces[fi], vertices);
        let d = [c[0] - center[0], c[1] - center[1], c[2] - center[2]];
        let axial = dot3(&d, axis_direction);
        let radial_vec = [
            d[0] - axial * axis_direction[0],
            d[1] - axial * axis_direction[1],
            d[2] - axial * axis_direction[2],
        ];
        let radial_dist = (radial_vec[0] * radial_vec[0]
            + radial_vec[1] * radial_vec[1]
            + radial_vec[2] * radial_vec[2]).sqrt();
        if radial_dist < 1e-10 {
            continue;
        }
        let radial_unit = [
            radial_vec[0] / radial_dist,
            radial_vec[1] / radial_dist,
            radial_vec[2] / radial_dist,
        ];
        let phi = f64::atan2(dot3(&n, axis_direction), dot3(&n, &radial_unit));
        minor_angles.push(phi);
    }

    if minor_angles.len() < MIN_TORUS_FACES {
        return false;
    }

    minor_angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = minor_angles.len();

    // Compute largest gap
    let mut largest_gap = minor_angles[0] + 2.0 * std::f64::consts::PI - minor_angles[n - 1];
    for i in 1..n {
        let gap = minor_angles[i] - minor_angles[i - 1];
        if gap > largest_gap {
            largest_gap = gap;
        }
    }

    let span = 2.0 * std::f64::consts::PI - largest_gap;
    // Require meaningful normal coverage around the tube. A real fillet
    // patch spans a significant tube arc (typically 90°); flat annuli and
    // cylinder bands span ~0° (plus tessellation noise). 0.3 rad (~17°)
    // admits patches of 3+ tessellation rows while rejecting degenerates.
    span >= 0.3
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
        let projected_points: Vec<[f64; 3]> = hyp.faces.iter().map(|&fi| {
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

            if radial_dist > 1e-15 {
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
            }
        }).collect();

        let (max_dist, _) = stage1::max_distance_to_shape(&projected_points, compare_shape);
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

// ---------------------------------------------------------------------------
// SurfaceType trait implementation for generic BFS
// ---------------------------------------------------------------------------

use super::bfs::{BfsConfig, Seed, SurfaceType};

/// Fitted torus parameters.
#[derive(Clone)]
pub(super) struct TorusFit {
    pub center: [f64; 3],
    pub axis_direction: [f64; 3],
    pub major_radius: f64,
    pub minor_radius: f64,
    pub convex: Option<bool>,
}

/// Marker struct for the toroidal surface type.
pub(super) struct ToroidalSurface;

/// Build a torus seed for a face committed to a cylindrical hypothesis by
/// gathering its 3-ring edge neighborhood. Returns an empty vec unless the
/// neighborhood spans ≥ 2 distinct conical hypotheses (evidence of major-arc
/// curvature that a single cylinder cannot represent).
fn strip_merge_seed(mesh: &ConnectedMesh, fi: usize, config: &BfsConfig, relaxed: bool) -> Vec<Seed> {
    let mut ring: Vec<usize> = vec![fi];
    let mut seen = vec![false; mesh.faces.len()];
    seen[fi] = true;
    let mut frontier = vec![fi];
    for _hop in 0..3 {
        let mut next = Vec::new();
        for &f in &frontier {
            let vc = mesh.faces[f].vertex_count as usize;
            for &ni in &mesh.faces[f].neighbors[..vc] {
                if ni < 0 {
                    continue;
                }
                let ni = ni as usize;
                if seen[ni] {
                    continue;
                }
                seen[ni] = true;
                let face = &mesh.faces[ni];
                let ni_large_sphere = if config.large_sphere.is_empty() {
                    face.spherical_hypothesis >= 0
                } else {
                    config.large_sphere.get(ni).copied().unwrap_or(false)
                };
                if face.normal.is_some()
                    && face.toroidal_hypothesis == UNDEDUCED_TOROIDAL_HYPOTHESIS
                    && !ni_large_sphere
                    && !config.is_junction(ni, relaxed)
                {
                    ring.push(ni);
                    next.push(ni);
                }
            }
        }
        frontier = next;
    }

    // Require faces from at least two distinct conical hypotheses.
    let mut hyp_ids: Vec<i32> = ring
        .iter()
        .map(|&f| mesh.faces[f].conical_hypothesis)
        .filter(|&h| h >= 0)
        .collect();
    hyp_ids.sort_unstable();
    hyp_ids.dedup();
    if hyp_ids.len() < 2 {
        return Vec::new();
    }

    let mut verts: Vec<usize> = Vec::new();
    for &sfi in &ring {
        let svc = mesh.faces[sfi].vertex_count as usize;
        verts.extend_from_slice(&mesh.faces[sfi].vertex_indices[..svc]);
    }
    verts.sort_unstable();
    verts.dedup();
    if ring.len() < MIN_TORUS_FACES || verts.len() < 9 {
        return Vec::new();
    }

    ring.sort_unstable();
    vec![Seed { faces: ring }]
}

impl SurfaceType for ToroidalSurface {
    type Fit = TorusFit;
    type Hypothesis = ToroidalHypothesis;

    fn name() -> &'static str { "toroidal" }

    fn hypothesis_field(face: &MeshFace) -> i32 { face.toroidal_hypothesis }

    fn set_hypothesis_field(face: &mut MeshFace, value: i32) {
        face.toroidal_hypothesis = value;
    }

    fn undeduced_value() -> i32 { UNDEDUCED_TOROIDAL_HYPOTHESIS }

    fn face_eligible(face: &MeshFace) -> bool {
        face.toroidal_hypothesis == UNDEDUCED_TOROIDAL_HYPOTHESIS
    }

    /// Torus quads twist (they are doubly curved), so stage 1.3 coplanar
    /// fusion often fails on fillet patches, leaving raw triangles whose
    /// normals differ by up to ~2× the tessellation angular step. Vertex
    /// distance is the primary discriminator, so relax the dihedral check
    /// like the conical stage does.
    fn angular_tol_multiplier() -> f64 { 2.0 }

    /// Corner vertex-blend patches (where two fillets meet) have low face
    /// indices and seed before the adjacent larger fillet tori commit, so
    /// their 1-ring seeds mix faces from both surfaces and fail. A second
    /// pass retries them once everything else is committed and excluded.
    fn seeding_passes() -> usize { 2 }

    /// Torus vertex-neighborhood seeds legitimately mix faces from two
    /// adjacent fillet patches; trimming relative outliers and refitting
    /// isolates the majority patch.
    fn seed_outlier_trimming() -> bool { true }

    fn generate_seeds(mesh: &ConnectedMesh, fi: usize, config: &BfsConfig, relaxed: bool) -> Vec<Seed> {
        if mesh.faces[fi].normal.is_none() {
            return Vec::new();
        }
        // Junction faces (curvature transitions) produce contaminated seeds;
        // pass 1 skips them as anchors and gathering members.
        if config.is_junction(fi, relaxed) {
            return Vec::new();
        }
        // Exclude faces on LARGE spheres only: small spurious spheres (a
        // corner vertex-blend patch locally fits a sphere of the tube
        // radius) must not block torus seeding.
        let on_large_sphere = |f: usize| -> bool {
            if config.large_sphere.is_empty() {
                mesh.faces[f].spherical_hypothesis >= 0
            } else {
                config.large_sphere.get(f).copied().unwrap_or(false)
            }
        };
        if on_large_sphere(fi) {
            return Vec::new();
        }
        // (Flat facets in seeds — e.g. top-face quads pulled into corner
        // vertex-blend seeds — were excluded by a multiface_planar mask
        // until 2026-07-15; the junction mask plus seed outlier trimming
        // proved to cover that failure, so the mask was retired.)
        // Faces committed to a cylinder may still lie on a torus: a coarsely
        // tessellated torus (e.g. a pipe elbow) is tiled by narrow cylinder
        // strips in stage 2.2, each spanning too little of the major arc for
        // a 1-ring torus fit. Seed from a 3-ring neighborhood instead, which
        // spans several strips and provides major-arc curvature — but only
        // when the neighborhood contains faces from ≥ 2 distinct conical
        // hypotheses (a plain cylinder region is a single hypothesis and is
        // skipped cheaply).
        if mesh.faces[fi].conical_hypothesis >= 0 {
            return strip_merge_seed(mesh, fi, config, relaxed);
        }

        // Gather the seed by hop-expansion from fi over edge neighbors,
        // excluding committed/planar/junction faces, growing the radius (up
        // to 3 hops) until the seed is large enough for a solid 7-DOF fit.
        // Small fillets (e.g. r=0.5 tubes only 2-3 rows wide) have little or
        // no interior left after junction exclusion across the tube, but
        // extend far along the fillet axis — hop growth follows the strip.
        let eligible = |fj: usize| -> bool {
            let face = &mesh.faces[fj];
            face.normal.is_some()
                && face.toroidal_hypothesis == UNDEDUCED_TOROIDAL_HYPOTHESIS
                && face.conical_hypothesis < 0
                && !on_large_sphere(fj)
                && !config.is_junction(fj, relaxed)
        };
        let mut seed_set: HashSet<usize> = HashSet::new();
        seed_set.insert(fi);
        let mut frontier: Vec<usize> = vec![fi];
        const SEED_TARGET_FACES: usize = 16;
        for _hop in 0..3 {
            if seed_set.len() >= SEED_TARGET_FACES {
                break;
            }
            let mut next: Vec<usize> = Vec::new();
            for &f in &frontier {
                let vc = mesh.faces[f].vertex_count as usize;
                for &ni in &mesh.faces[f].neighbors[..vc] {
                    if ni >= 0 {
                        let ni = ni as usize;
                        if !seed_set.contains(&ni) && eligible(ni) {
                            seed_set.insert(ni);
                            next.push(ni);
                        }
                    }
                }
            }
            if next.is_empty() {
                break;
            }
            frontier = next;
        }

        // Need enough faces and vertices for 7-DOF fit.
        // Sorted for deterministic face iteration order in fitting.
        let mut faces: Vec<usize> = seed_set.into_iter().collect();
        faces.sort_unstable();
        let mut verts: HashSet<usize> = HashSet::new();
        for &sfi in &faces {
            let svc = mesh.faces[sfi].vertex_count as usize;
            for &vi in &mesh.faces[sfi].vertex_indices[..svc] {
                verts.insert(vi);
            }
        }
        if faces.len() < MIN_TORUS_FACES || verts.len() < 9 {
            return Vec::new();
        }

        vec![Seed { faces }]
    }

    fn fit_surface(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        vertex_list: &[usize],
        verbosity: u8,
    ) -> Option<Self::Fit> {
        // Sorted, deduplicated for deterministic fitting sums.
        let mut vertex_set: Vec<usize> = vertex_list.to_vec();
        vertex_set.sort_unstable();
        vertex_set.dedup();
        let fit = fit_torus(face_list, &vertex_set, &mesh.faces, &mesh.vertices);
        if fit.is_none() && verbosity >= 5 {
            for &vi in &vertex_set {
                let v = &mesh.vertices[vi];
                eprintln!("        torus fit input vtx {vi}: ({:.4}, {:.4}, {:.4})", v.x, v.y, v.z);
            }
        }
        let (center, axis, major_r, minor_r) = fit?;
        if minor_r <= 0.0 || major_r <= 0.0 || minor_r > major_r {
            return None;
        }
        Some(TorusFit {
            center,
            axis_direction: axis,
            major_radius: major_r,
            minor_radius: minor_r,
            convex: None,
        })
    }

    fn vertex_distance(fit: &Self::Fit, pos: &[f64; 3]) -> f64 {
        let v = MeshVertex::from_xyz(pos[0], pos[1], pos[2]);
        vertex_to_torus_distance(&v, &fit.center, &fit.axis_direction, fit.major_radius, fit.minor_radius).abs()
    }

    fn convexity_check(fit: &Self::Fit, mesh: &ConnectedMesh, fi: usize) -> bool {
        if let Some(convex) = fit.convex {
            let fi_convex = determine_convexity(
                &mesh.faces[fi], &mesh.vertices,
                &fit.center, &fit.axis_direction, fit.major_radius,
            );
            fi_convex == convex
        } else {
            true
        }
    }

    fn determine_convexity(mesh: &ConnectedMesh, face_list: &[usize], fit: &mut Self::Fit) {
        if let Some(&first) = face_list.first() {
            fit.convex = Some(determine_convexity(
                &mesh.faces[first], &mesh.vertices,
                &fit.center, &fit.axis_direction, fit.major_radius,
            ));
        }
    }

    fn acceptance_test(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        fit: &Self::Fit,
        _config: &BfsConfig,
    ) -> bool {
        if face_list.len() < MIN_TORUS_FACES {
            return false;
        }
        // Reject degenerate giant-radius tori: a straight-edge fillet
        // (cylinder) is the R→∞ limit of a torus, and a short segment of one
        // can be fitted within tolerance by a torus with an absurdly distant
        // axis (observed R ≈ 1500–3600 mm on ~100 mm models). A legitimate
        // torus's major radius cannot exceed the mesh extent. Analogous to
        // the cone apex-distance bound.
        let bb_diag = super::bounding_box_diagonal(&mesh.vertices);
        if fit.major_radius > bb_diag {
            return false;
        }
        // Major circle angular coverage (azimuthal)
        if !angular_coverage_valid(
            face_list, &fit.center, &fit.axis_direction,
            &mesh.faces, &mesh.vertices,
        ) {
            return false;
        }
        // Minor circle angular coverage (tube cross-section)
        minor_circle_coverage_valid(
            face_list, &fit.center, &fit.axis_direction, fit.major_radius,
            &mesh.faces, &mesh.vertices,
        )
    }

    fn viz_overlays(
        fit: &Self::Fit,
        _mesh: &ConnectedMesh,
        _face_list: &[usize],
    ) -> (Vec<viz::CylinderOverlay>, Vec<viz::ConeOverlay>, Vec<viz::SphereOverlay>, Vec<viz::TorusOverlay>) {
        let torus = viz::TorusOverlay {
            center: fit.center,
            direction: fit.axis_direction,
            major_radius: fit.major_radius,
            minor_radius: fit.minor_radius,
            color: [0.2, 0.4, 1.0, 0.3],
        };
        (Vec::new(), Vec::new(), Vec::new(), vec![torus])
    }

    fn viz_fit_description(fit: &Self::Fit) -> String {
        format!("R={:.4} r={:.4}", fit.major_radius, fit.minor_radius)
    }

    fn make_hypothesis(
        mesh: &ConnectedMesh,
        faces: Vec<usize>,
        vertices: Vec<usize>,
        fit: &Self::Fit,
    ) -> ToroidalHypothesis {
        let mut error_max = 0.0_f64;
        let mut centroid_error_max = 0.0_f64;
        let mut error_abs_sum = 0.0_f64;
        for &vi in &vertices {
            let v = &mesh.vertices[vi];
            let d = vertex_to_torus_distance(
                v, &fit.center, &fit.axis_direction, fit.major_radius, fit.minor_radius,
            ).abs();
            error_max = error_max.max(d);
            error_abs_sum += d;
        }
        for &fi in &faces {
            let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
            let v = MeshVertex::from_xyz(c[0], c[1], c[2]);
            let d = vertex_to_torus_distance(
                &v, &fit.center, &fit.axis_direction, fit.major_radius, fit.minor_radius,
            ).abs();
            centroid_error_max = centroid_error_max.max(d);
        }
        ToroidalHypothesis {
            center: fit.center,
            axis_direction: fit.axis_direction,
            major_radius: fit.major_radius,
            minor_radius: fit.minor_radius,
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

    /// Generate points on a torus surface.
    fn torus_points(
        center: [f64; 3],
        axis_dir: [f64; 3],
        major_r: f64,
        minor_r: f64,
        n_major: usize,
        n_minor: usize,
        theta_min: f64,
        theta_max: f64,
        phi_min: f64,
        phi_max: f64,
    ) -> Vec<MeshVertex> {
        let d = {
            let len = (axis_dir[0].powi(2) + axis_dir[1].powi(2) + axis_dir[2].powi(2)).sqrt();
            [axis_dir[0] / len, axis_dir[1] / len, axis_dir[2] / len]
        };
        let (u, w) = super::super::perpendicular_basis(&d);
        let mut pts = Vec::new();
        for i in 0..n_major {
            let theta = theta_min + (theta_max - theta_min) * i as f64 / n_major as f64;
            let (sin_t, cos_t) = theta.sin_cos();
            let rad = [cos_t * u[0] + sin_t * w[0],
                       cos_t * u[1] + sin_t * w[1],
                       cos_t * u[2] + sin_t * w[2]];
            let tc = [center[0] + major_r * rad[0],
                      center[1] + major_r * rad[1],
                      center[2] + major_r * rad[2]];
            for j in 0..n_minor {
                let phi = phi_min + (phi_max - phi_min) * j as f64 / n_minor as f64;
                let (sin_p, cos_p) = phi.sin_cos();
                pts.push(MeshVertex::from_xyz(
                    tc[0] + minor_r * (cos_p * rad[0] + sin_p * d[0]),
                    tc[1] + minor_r * (cos_p * rad[1] + sin_p * d[1]),
                    tc[2] + minor_r * (cos_p * rad[2] + sin_p * d[2]),
                ));
            }
        }
        pts
    }

    /// Generate normals for a torus.
    fn torus_normals(
        axis_dir: [f64; 3],
        n_major: usize,
        n_minor: usize,
        theta_min: f64,
        theta_max: f64,
        phi_min: f64,
        phi_max: f64,
    ) -> Vec<[f64; 3]> {
        let d = {
            let len = (axis_dir[0].powi(2) + axis_dir[1].powi(2) + axis_dir[2].powi(2)).sqrt();
            [axis_dir[0] / len, axis_dir[1] / len, axis_dir[2] / len]
        };
        let (u, w) = super::super::perpendicular_basis(&d);
        let mut normals = Vec::new();
        for i in 0..n_major {
            let theta = theta_min + (theta_max - theta_min) * i as f64 / n_major as f64;
            let (sin_t, cos_t) = theta.sin_cos();
            let rad = [cos_t * u[0] + sin_t * w[0],
                       cos_t * u[1] + sin_t * w[1],
                       cos_t * u[2] + sin_t * w[2]];
            for j in 0..n_minor {
                let phi = phi_min + (phi_max - phi_min) * j as f64 / n_minor as f64;
                let (sin_p, cos_p) = phi.sin_cos();
                normals.push([
                    cos_p * rad[0] + sin_p * d[0],
                    cos_p * rad[1] + sin_p * d[1],
                    cos_p * rad[2] + sin_p * d[2],
                ]);
            }
        }
        normals
    }

    /// Build a synthetic torus mesh.
    fn build_torus_mesh(
        center: [f64; 3],
        axis_dir: [f64; 3],
        major_r: f64,
        minor_r: f64,
        n_major: usize,
        n_minor: usize,
        theta_min: f64,
        theta_max: f64,
        phi_min: f64,
        phi_max: f64,
    ) -> (Vec<MeshFace>, Vec<MeshVertex>) {
        let vertices = torus_points(center, axis_dir, major_r, minor_r, n_major, n_minor,
                                    theta_min, theta_max, phi_min, phi_max);
        let normals = torus_normals(axis_dir, n_major, n_minor,
                                   theta_min, theta_max, phi_min, phi_max);
        let mut faces = Vec::new();
        for i in 0..n_major {
            let i_next = (i + 1) % n_major;
            for j in 0..n_minor {
                let j_next = (j + 1) % n_minor;
                let i0 = i * n_minor + j;
                let i1 = i * n_minor + j_next;
                let i2 = i_next * n_minor + j;
                let fidx = faces.len();
                faces.push(MeshFace {
                    vertex_count: 3,
                    vertex_indices: [i0, i1, i2, 0],
                    neighbors: [-1, -1, -1, -1],
                    normal: Some(normals[fidx % normals.len()]),
                    planar_hypothesis: -1,
                    spherical_hypothesis: -1,
                    conical_hypothesis: -1,
                    toroidal_hypothesis: -2,
                });
            }
        }
        (faces, vertices)
    }

    // ---- Distance tests ----

    #[test]
    fn distance_to_torus_outer_equator() {
        let r_major = 10.0;
        let r_minor = 3.0;
        let v = MeshVertex::from_xyz(r_major + r_minor, 0.0, 0.0);
        let d = vertex_to_torus_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], r_major, r_minor);
        assert!(d.abs() < 1e-10, "outer equator distance should be 0, got {d}");
    }

    #[test]
    fn distance_to_torus_inner_equator() {
        let r_major = 10.0;
        let r_minor = 3.0;
        let v = MeshVertex::from_xyz(r_major - r_minor, 0.0, 0.0);
        let d = vertex_to_torus_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], r_major, r_minor);
        assert!(d.abs() < 1e-10, "inner equator distance should be 0, got {d}");
    }

    #[test]
    fn distance_to_torus_top() {
        let r_major = 10.0;
        let r_minor = 3.0;
        let v = MeshVertex::from_xyz(r_major, 0.0, r_minor);
        let d = vertex_to_torus_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], r_major, r_minor);
        assert!(d.abs() < 1e-10, "top distance should be 0, got {d}");
    }

    #[test]
    fn distance_to_torus_tube_center() {
        let r_major = 10.0;
        let r_minor = 3.0;
        let v = MeshVertex::from_xyz(r_major, 0.0, 0.0);
        let d = vertex_to_torus_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], r_major, r_minor);
        assert!((d - (-r_minor)).abs() < 1e-10, "tube center distance should be -{r_minor}, got {d}");
    }

    #[test]
    fn distance_to_torus_at_center() {
        let r_major = 10.0;
        let r_minor = 3.0;
        let v = MeshVertex::from_xyz(0.0, 0.0, 0.0);
        let d = vertex_to_torus_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], r_major, r_minor);
        assert!((d - (r_major - r_minor)).abs() < 1e-10, "center distance should be R-r, got {d}");
    }

    #[test]
    fn distance_to_torus_outside() {
        let r_major = 10.0;
        let r_minor = 3.0;
        let v = MeshVertex::from_xyz(r_major + r_minor + 1.0, 0.0, 0.0);
        let d = vertex_to_torus_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], r_major, r_minor);
        assert!((d - 1.0).abs() < 1e-10, "1 unit outside distance should be 1, got {d}");
    }

    // ---- Fit tests ----

    #[test]
    fn torus_fit_z_axis() {
        let (faces, vertices) = build_torus_mesh(
            [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 7.0, 3.0,
            24, 12, 0.0, 2.0 * PI, 0.0, 2.0 * PI,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_torus(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "torus fit should succeed");
        let (center, dir, major_r, minor_r) = result.unwrap();

        let axis_z = dir[2].abs();
        assert!((axis_z - 1.0).abs() < 0.01, "axis should be along Z, got z={axis_z}");
        assert!((major_r - 7.0).abs() < 0.05, "major_radius should be 7, got {major_r}");
        assert!((minor_r - 3.0).abs() < 0.05, "minor_radius should be 3, got {minor_r}");
        assert!(center[0].abs() < 0.05, "center.x should be 0, got {}", center[0]);
        assert!(center[1].abs() < 0.05, "center.y should be 0, got {}", center[1]);
        assert!(center[2].abs() < 0.05, "center.z should be 0, got {}", center[2]);
    }

    #[test]
    fn torus_fit_offset() {
        let (faces, vertices) = build_torus_mesh(
            [5.0, 10.0, -3.0], [0.0, 0.0, 1.0], 9.0, 2.0,
            24, 12, 0.0, 2.0 * PI, 0.0, 2.0 * PI,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_torus(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "offset torus fit should succeed");
        let (center, dir, major_r, minor_r) = result.unwrap();

        let axis_z = dir[2].abs();
        assert!((axis_z - 1.0).abs() < 0.01, "axis should be along Z, got z={axis_z}");
        assert!((major_r - 9.0).abs() < 0.05, "major_radius should be 9, got {major_r}");
        assert!((minor_r - 2.0).abs() < 0.05, "minor_radius should be 2, got {minor_r}");
        assert!((center[0] - 5.0).abs() < 0.1, "center.x should be 5, got {}", center[0]);
        assert!((center[1] - 10.0).abs() < 0.1, "center.y should be 10, got {}", center[1]);
        assert!((center[2] - (-3.0)).abs() < 0.1, "center.z should be -3, got {}", center[2]);
    }

    #[test]
    fn torus_fit_corner_blend_patch() {
        // Corner vertex-blend patch as found on stepped_block_toroidal_fillets:
        // a small torus (R=2.5, r=0.5) patch spanning only 90° of the major
        // circle and 90° of the tube. Both the tube-center method and
        // position-PCA axis estimation fail on such compact doubly-curved
        // patches; the multi-axis (h,r) rescue must recover it.
        let (faces, vertices) = build_torus_mesh(
            [-66.75, -64.5, -18.2], [0.0, 1.0, 0.0], 2.5, 0.5,
            8, 5, 0.0, PI / 2.0, 0.0, PI / 2.0,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_torus(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "corner blend patch fit should succeed");
        let (center, dir, major_r, minor_r) = result.unwrap();
        assert!(dir[1].abs() > 0.999, "axis should be along Y, got {dir:?}");
        assert!((major_r - 2.5).abs() < 0.02, "major_radius should be 2.5, got {major_r}");
        assert!((minor_r - 0.5).abs() < 0.02, "minor_radius should be 0.5, got {minor_r}");
        assert!((center[0] - (-66.75)).abs() < 0.05, "center.x, got {}", center[0]);
        assert!((center[1] - (-64.5)).abs() < 0.05, "center.y, got {}", center[1]);
        assert!((center[2] - (-18.2)).abs() < 0.05, "center.z, got {}", center[2]);
    }

    #[test]
    fn torus_fit_tilted() {
        let (faces, vertices) = build_torus_mesh(
            [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 8.0, 1.5,
            24, 12, 0.0, 2.0 * PI, 0.0, 2.0 * PI,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_torus(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "tilted torus fit should succeed");
        let (_center, dir, major_r, minor_r) = result.unwrap();

        let s3 = 1.0 / 3.0_f64.sqrt();
        let dot = (dir[0] * s3 + dir[1] * s3 + dir[2] * s3).abs();
        assert!((dot - 1.0).abs() < 0.01, "axis should be along (1,1,1), dot={dot}");
        assert!((major_r - 8.0).abs() < 0.05, "major_radius should be 8, got {major_r}");
        assert!((minor_r - 1.5).abs() < 0.05, "minor_radius should be 1.5, got {minor_r}");
    }

    #[test]
    fn torus_fit_quarter_tube() {
        // Quarter-tube cross-section (fillet patch).
        let (faces, vertices) = build_torus_mesh(
            [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 8.5, 1.5,
            24, 8, 0.0, 2.0 * PI, 0.0, PI / 2.0,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_torus(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "quarter-tube torus fit should succeed");
        let (_center, dir, major_r, minor_r) = result.unwrap();

        assert!((major_r - 8.5).abs() < 0.1, "major_radius should be 8.5, got {major_r}");
        assert!((minor_r - 1.5).abs() < 0.1, "minor_radius should be 1.5, got {minor_r}");
        let axis_z = dir[2].abs();
        assert!((axis_z - 1.0).abs() < 0.01, "axis should be along Z, got z={axis_z}");
    }

    #[test]
    fn torus_fit_half_major() {
        // Only half the major circle.
        let (faces, vertices) = build_torus_mesh(
            [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 7.0, 3.0,
            12, 12, 0.0, PI, 0.0, 2.0 * PI,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_torus(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "half-major torus fit should succeed");
        let (_center, _dir, major_r, minor_r) = result.unwrap();

        assert!((major_r - 7.0).abs() < 0.1, "major_radius should be 7, got {major_r}");
        assert!((minor_r - 3.0).abs() < 0.1, "minor_radius should be 3, got {minor_r}");
    }

    #[test]
    fn torus_fit_thin_fillet() {
        // Very thin fillet: large R/r ratio.
        let (faces, vertices) = build_torus_mesh(
            [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 50.0, 0.5,
            36, 12, 0.0, 2.0 * PI, 0.0, 2.0 * PI,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_torus(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "thin fillet torus fit should succeed");
        let (_center, _dir, major_r, minor_r) = result.unwrap();

        assert!((major_r - 50.0).abs() < 0.1, "major_radius should be 50, got {major_r}");
        assert!((minor_r - 0.5).abs() < 0.05, "minor_radius should be 0.5, got {minor_r}");
    }

    #[test]
    fn torus_fit_rejects_degenerate() {
        let faces = vec![MeshFace {
            vertex_count: 3,
            vertex_indices: [0, 1, 2, 0],
            neighbors: [-1, -1, -1, -1],
            normal: Some([0.0, 0.0, 1.0]),
            planar_hypothesis: -1,
            spherical_hypothesis: -1,
            conical_hypothesis: -1,
            toroidal_hypothesis: -2,
        }];
        let vertices = vec![
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0),
            MeshVertex::from_xyz(0.0, 0.0, 1.0),
            MeshVertex::from_xyz(1.0, 1.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 1.0),
            MeshVertex::from_xyz(1.0, 0.0, 1.0),
        ];
        let vertex_set: Vec<usize> = (0..6).collect();
        assert!(fit_torus(&[0], &vertex_set, &faces, &vertices).is_none(),
                "torus fit should reject degenerate input");
    }
}

