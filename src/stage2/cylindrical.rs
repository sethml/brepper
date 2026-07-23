//! Stage 2 cylindrical fitting utilities.
//! Core fitting functions used by ConicalSurface (in conical.rs) for cylinder fitting.


use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{Const, Dyn, OMatrix, OVector, Owned};

use crate::stage1::{MeshFace, MeshVertex};

use super::{
    build_normal_covariance, cross3,
    normalize3, perpendicular_basis, smallest_eigenvector_3x3,
};

// ---------------------------------------------------------------------------
// Cylindrical fitting helpers
// ---------------------------------------------------------------------------

/// Compute the distance from a vertex to a cylinder surface.
/// Returns the signed distance: positive if outside the cylinder, negative if inside.
pub(super) fn vertex_to_cylinder_distance(
    v: &MeshVertex,
    axis_origin: &[f64; 3],
    axis_direction: &[f64; 3],
    radius: f64,
) -> f64 {
    // Vector from axis origin to vertex
    let dx = v.x - axis_origin[0];
    let dy = v.y - axis_origin[1];
    let dz = v.z - axis_origin[2];
    // Project onto axis
    let t = dx * axis_direction[0] + dy * axis_direction[1] + dz * axis_direction[2];
    // Radial vector (perpendicular to axis)
    let rx = dx - t * axis_direction[0];
    let ry = dy - t * axis_direction[1];
    let rz = dz - t * axis_direction[2];
    let radial_dist = (rx * rx + ry * ry + rz * rz).sqrt();
    radial_dist - radius
}

/// Given an axis direction, project vertices onto the perpendicular plane
/// and fit a 2D circle. Returns (axis_origin_3d, radius) or None if degenerate.
fn fit_circle_for_axis(
    axis_dir: &[f64; 3],
    vertex_set: &[usize],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], f64)> {
    let (u, w) = perpendicular_basis(axis_dir);
    let verts: Vec<usize> = vertex_set.to_vec();
    let n = verts.len();
    if n < 3 {
        return None;
    }

    let mut xs = Vec::with_capacity(n);
    let mut ys = Vec::with_capacity(n);
    for &vi in &verts {
        let v = &vertices[vi];
        xs.push(v.x * u[0] + v.y * u[1] + v.z * u[2]);
        ys.push(v.x * w[0] + v.y * w[1] + v.z * w[2]);
    }

    // Center the 2D data to improve numerical conditioning.
    let nf = n as f64;
    let mean_x = xs.iter().sum::<f64>() / nf;
    let mean_y = ys.iter().sum::<f64>() / nf;
    for i in 0..n {
        xs[i] -= mean_x;
        ys[i] -= mean_y;
    }

    // Algebraic circle fit: x² + y² + Dx + Ey + F = 0
    let mut sx = 0.0_f64;
    let mut sy = 0.0_f64;
    let mut sx2 = 0.0_f64;
    let mut sy2 = 0.0_f64;
    let mut sxy = 0.0_f64;
    let mut sx3 = 0.0_f64;
    let mut sy3 = 0.0_f64;
    let mut sx2y = 0.0_f64;
    let mut sxy2 = 0.0_f64;

    for i in 0..n {
        let x = xs[i];
        let y = ys[i];
        sx += x;
        sy += y;
        sx2 += x * x;
        sy2 += y * y;
        sxy += x * y;
        sx3 += x * x * x;
        sy3 += y * y * y;
        sx2y += x * x * y;
        sxy2 += x * y * y;
    }

    let rhs0 = -(sx3 + sxy2);
    let rhs1 = -(sx2y + sy3);
    let rhs2 = -(sx2 + sy2);

    let a = [[sx2, sxy, sx], [sxy, sy2, sy], [sx, sy, nf]];

    let det_a = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
              - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
              + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

    if det_a.abs() < 1e-30 {
        return None;
    }

    let inv_det = 1.0 / det_a;

    let d_val = inv_det * (rhs0 * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                         - a[0][1] * (rhs1 * a[2][2] - a[1][2] * rhs2)
                         + a[0][2] * (rhs1 * a[2][1] - a[1][1] * rhs2));

    let e_val = inv_det * (a[0][0] * (rhs1 * a[2][2] - a[1][2] * rhs2)
                         - rhs0 * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                         + a[0][2] * (a[1][0] * rhs2 - rhs1 * a[2][0]));

    let f_val = inv_det * (a[0][0] * (a[1][1] * rhs2 - rhs1 * a[2][1])
                         - a[0][1] * (a[1][0] * rhs2 - rhs1 * a[2][0])
                         + rhs0 * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));

    let center_x = -d_val / 2.0;
    let center_y = -e_val / 2.0;
    let r_sq = center_x * center_x + center_y * center_y - f_val;
    if r_sq <= 0.0 {
        return None;
    }
    let radius = r_sq.sqrt();

    let abs_center_x = center_x + mean_x;
    let abs_center_y = center_y + mean_y;
    let axis_origin = [
        abs_center_x * u[0] + abs_center_y * w[0],
        abs_center_x * u[1] + abs_center_y * w[1],
        abs_center_x * u[2] + abs_center_y * w[2],
    ];

    Some((axis_origin, radius))
}

/// LM problem for cylinder fitting.
/// Parameters: [alpha, beta, qx, qy, qz, radius]
/// alpha, beta: tilt from initial axis in perpendicular directions
/// qx, qy, qz: axis point
/// radius: cylinder radius
struct CylinderLMProblem {
    points: Vec<[f64; 3]>,
    a0: [f64; 3],
    u0: [f64; 3],
    w0: [f64; 3],
    params: OVector<f64, Const<6>>,
}

impl CylinderLMProblem {
    fn new(
        points: Vec<[f64; 3]>,
        initial_axis: [f64; 3],
        initial_origin: [f64; 3],
        initial_radius: f64,
    ) -> Self {
        let (u0, w0) = perpendicular_basis(&initial_axis);
        let params = OVector::<f64, Const<6>>::new(
            0.0, 0.0,
            initial_origin[0], initial_origin[1], initial_origin[2],
            initial_radius,
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
        let q = [params[2], params[3], params[4]];
        let r = params[5];
        let n = self.points.len();
        let mut residuals = OVector::<f64, Dyn>::zeros_generic(Dyn(n), Const::<1>);
        for (i, p) in self.points.iter().enumerate() {
            let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
            let c = cross3(&d, &a);
            let h = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
            residuals[i] = h - r;
        }
        residuals
    }
}

impl LeastSquaresProblem<f64, Dyn, Const<6>> for CylinderLMProblem {
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

/// Fit cylinder parameters (axis direction, axis origin, radius) from a set of faces.
/// Tries two axis estimation strategies — normal-covariance (good for wide arcs) and
/// vertex-PCA (good for narrow arcs / drum-laced triangles) — picks the better one,
/// then refines the axis direction by minimizing vertex-to-cylinder SSE.
/// Returns (axis_origin, axis_direction, radius).
pub(super) fn fit_cylinder(
    face_indices: &[usize],
    vertex_set: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64)> {
    if vertex_set.len() < 3 {
        return None;
    }

    // Axis from smallest eigenvector of area-weighted normal covariance.
    let cov = build_normal_covariance(face_indices, faces, vertices);
    let mut axis_dir = smallest_eigenvector_3x3(&cov);
    let len = normalize3(&mut axis_dir);
    if len < 1e-15 {
        return None;
    }

    // Initial circle fit to get starting origin and radius.
    let (axis_origin, radius) = fit_circle_for_axis(&axis_dir, vertex_set, vertices)?;

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

    let points: Vec<[f64; 3]> = vertex_set.iter()
        .map(|&vi| [
            vertices[vi].x - centroid[0],
            vertices[vi].y - centroid[1],
            vertices[vi].z - centroid[2],
        ])
        .collect();
    let centered_origin = [
        axis_origin[0] - centroid[0],
        axis_origin[1] - centroid[1],
        axis_origin[2] - centroid[2],
    ];

    // Refine with Levenberg-Marquardt optimization on radial residuals.
    let problem = CylinderLMProblem::new(points, axis_dir, centered_origin, radius);
    let (result, _report) = LevenbergMarquardt::new().with_patience(150).minimize(problem);
    let refined_dir = result.axis_dir_from_params(&result.params);
    let refined_origin = [
        result.params[2] + centroid[0],
        result.params[3] + centroid[1],
        result.params[4] + centroid[2],
    ];
    let refined_radius = result.params[5];

    if refined_radius > 0.0 {
        Some((refined_origin, refined_dir, refined_radius))
    } else {
        Some((axis_origin, axis_dir, radius))
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Generate points on a cylinder: n_rings rings of n_pts_per_ring points.
    fn cylinder_points(
        axis_point: [f64; 3],
        axis_dir: [f64; 3],
        radius: f64,
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
        let mut pts = Vec::new();
        for ring in 0..n_rings {
            let t = if n_rings == 1 { 0.5 } else { ring as f64 / (n_rings - 1) as f64 };
            let h = h_min + t * (h_max - h_min);
            for j in 0..n_pts_per_ring {
                let theta = 2.0 * PI * j as f64 / n_pts_per_ring as f64;
                let (sin_t, cos_t) = theta.sin_cos();
                pts.push(MeshVertex::from_xyz(
                    axis_point[0] + h * d[0] + radius * (cos_t * u[0] + sin_t * w[0]),
                    axis_point[1] + h * d[1] + radius * (cos_t * u[1] + sin_t * w[1]),
                    axis_point[2] + h * d[2] + radius * (cos_t * u[2] + sin_t * w[2]),
                ));
            }
        }
        pts
    }

    /// Generate face normals for a cylinder (radially outward).
    fn cylinder_normals(axis_dir: [f64; 3], n_faces: usize) -> Vec<[f64; 3]> {
        let d = {
            let len = (axis_dir[0].powi(2) + axis_dir[1].powi(2) + axis_dir[2].powi(2)).sqrt();
            [axis_dir[0] / len, axis_dir[1] / len, axis_dir[2] / len]
        };
        let (u, w) = super::super::perpendicular_basis(&d);
        (0..n_faces)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n_faces as f64;
                let (sin_t, cos_t) = theta.sin_cos();
                [cos_t * u[0] + sin_t * w[0],
                 cos_t * u[1] + sin_t * w[1],
                 cos_t * u[2] + sin_t * w[2]]
            })
            .collect()
    }

    /// Build a synthetic cylinder mesh suitable for fit_cylinder.
    /// Returns (faces, vertices) where faces have normals and vertex_indices set.
    fn build_cylinder_mesh(
        axis_point: [f64; 3],
        axis_dir: [f64; 3],
        radius: f64,
        n_rings: usize,
        n_pts_per_ring: usize,
        h_min: f64,
        h_max: f64,
    ) -> (Vec<MeshFace>, Vec<MeshVertex>) {
        let vertices = cylinder_points(axis_point, axis_dir, radius, n_rings, n_pts_per_ring, h_min, h_max);
        let normals = cylinder_normals(axis_dir, n_rings * n_pts_per_ring);
        // Create one face per vertex (using adjacent vertices as a triangle)
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
                });
            }
        }
        (faces, vertices)
    }

    // ---- Distance tests ----

    #[test]
    fn distance_to_cylinder_on_surface() {
        let v = MeshVertex::from_xyz(10.0, 0.0, 5.0);
        let d = vertex_to_cylinder_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], 10.0);
        assert!(d.abs() < 1e-10, "on-surface distance should be 0, got {d}");
    }

    #[test]
    fn distance_to_cylinder_inside() {
        let v = MeshVertex::from_xyz(5.0, 0.0, 5.0);
        let d = vertex_to_cylinder_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], 10.0);
        assert!((d - (-5.0)).abs() < 1e-10, "inside distance should be -5, got {d}");
    }

    #[test]
    fn distance_to_cylinder_outside() {
        let v = MeshVertex::from_xyz(15.0, 0.0, 5.0);
        let d = vertex_to_cylinder_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], 10.0);
        assert!((d - 5.0).abs() < 1e-10, "outside distance should be 5, got {d}");
    }

    #[test]
    fn distance_to_cylinder_diagonal() {
        let s = 10.0 / 2.0_f64.sqrt();
        let v = MeshVertex::from_xyz(s, s, 0.0);
        let d = vertex_to_cylinder_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], 10.0);
        assert!(d.abs() < 1e-3, "diagonal on-surface distance should be ~0, got {d}");
    }

    #[test]
    fn distance_to_cylinder_on_axis() {
        let v = MeshVertex::from_xyz(0.0, 0.0, 5.0);
        let d = vertex_to_cylinder_distance(&v, &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], 10.0);
        assert!((d - (-10.0)).abs() < 1e-10, "on-axis distance should be -10, got {d}");
    }

    // ---- Fit tests ----

    #[test]
    fn cylinder_fit_z_axis() {
        let (faces, vertices) = build_cylinder_mesh(
            [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 10.0, 6, 16, 0.0, 30.0,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_cylinder(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "cylinder fit should succeed");
        let (origin, dir, radius) = result.unwrap();

        assert!((radius - 10.0).abs() < 0.05, "radius should be 10, got {radius}");
        let axis_z = dir[2].abs();
        assert!((axis_z - 1.0).abs() < 0.01, "axis should be along Z, got z={axis_z}");
        assert!(origin[0].abs() < 0.05, "origin.x should be near 0, got {}", origin[0]);
        assert!(origin[1].abs() < 0.05, "origin.y should be near 0, got {}", origin[1]);
    }

    #[test]
    fn cylinder_fit_offset() {
        let (faces, vertices) = build_cylinder_mesh(
            [5.0, 10.0, -3.0], [1.0, 0.0, 0.0], 7.0, 8, 20, -15.0, 15.0,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_cylinder(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "offset cylinder fit should succeed");
        let (origin, dir, radius) = result.unwrap();

        let axis_x = dir[0].abs();
        assert!((axis_x - 1.0).abs() < 0.01, "axis should be along X, got x={axis_x}");
        assert!((radius - 7.0).abs() < 0.05, "radius should be 7, got {radius}");
        assert!((origin[1] - 10.0).abs() < 0.1, "origin.y should be 10, got {}", origin[1]);
        assert!((origin[2] - (-3.0)).abs() < 0.1, "origin.z should be -3, got {}", origin[2]);
    }

    #[test]
    fn cylinder_fit_tilted() {
        let (faces, vertices) = build_cylinder_mesh(
            [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 5.0, 10, 24, -20.0, 20.0,
        );
        let face_indices: Vec<usize> = (0..faces.len()).collect();
        let vertex_set: Vec<usize> = (0..vertices.len()).collect();

        let result = fit_cylinder(&face_indices, &vertex_set, &faces, &vertices);
        assert!(result.is_some(), "tilted cylinder fit should succeed");
        let (_origin, dir, radius) = result.unwrap();

        let s3 = 1.0 / 3.0_f64.sqrt();
        let dot = (dir[0] * s3 + dir[1] * s3 + dir[2] * s3).abs();
        assert!((dot - 1.0).abs() < 0.01, "axis should be along (1,1,1), dot={dot}");
        assert!((radius - 5.0).abs() < 0.05, "radius should be 5, got {radius}");
    }

    #[test]
    fn cylinder_fit_rejects_degenerate() {
        // Too few vertices (< 3)
        let faces = vec![MeshFace {
            vertex_count: 3,
            vertex_indices: [0, 1, 0, 0],
            neighbors: [-1, -1, -1, -1],
            normal: Some([0.0, 0.0, 1.0]),
        }];
        let vertices = vec![
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0),
        ];
        let vertex_set: Vec<usize> = (0..2).collect();
        assert!(fit_cylinder(&[0], &vertex_set, &faces, &vertices).is_none(),
                "cylinder fit should reject 2 vertices");
    }
}

