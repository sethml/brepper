//! Stage 2.3: Spherical hypothesis deduction and comparison.

use std::collections::HashSet;

use opencascade_sys::gp;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_SPHERICAL_HYPOTHESIS};
use crate::viz;

use super::{
    dot3, eigenvalues_3x3_symmetric, face_area, face_centroid, normalize3,
    SphericalHypothesis, Stage2CompareError,
    MAX_SPHERE_RADIUS_FACTOR, MIN_SPHERE_EIGENVALUE_RATIO, MIN_SPHERE_FACES,
};

// ---------------------------------------------------------------------------
// Spherical fitting helpers
// ---------------------------------------------------------------------------

/// Compute the bounding box diagonal of a mesh.
pub(super) fn bounding_box_diagonal(vertices: &[MeshVertex]) -> f64 {
    if vertices.is_empty() {
        return 0.0;
    }
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for v in vertices {
        min[0] = min[0].min(v.x);
        min[1] = min[1].min(v.y);
        min[2] = min[2].min(v.z);
        max[0] = max[0].max(v.x);
        max[1] = max[1].max(v.y);
        max[2] = max[2].max(v.z);
    }
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Compute the distance from a vertex to a sphere surface.
/// Returns signed distance: positive if outside the sphere, negative if inside.
pub(super) fn vertex_to_sphere_distance(
    v: &MeshVertex,
    center: &[f64; 3],
    radius: f64,
) -> f64 {
    let dx = v.x - center[0];
    let dy = v.y - center[1];
    let dz = v.z - center[2];
    (dx * dx + dy * dy + dz * dz).sqrt() - radius
}

/// Fit sphere parameters (center, radius) from a set of vertices.
/// Uses algebraic least-squares: |v - c|² = r², expanded linearly.
/// Returns (center, radius) or None if degenerate.
#[allow(clippy::needless_range_loop)]
fn fit_sphere(
    vertex_set: &HashSet<usize>,
    vertices: &[MeshVertex],
) -> Option<([f64; 3], f64)> {
    let n = vertex_set.len();
    if n < 4 {
        return None;
    }

    // Solve: 2*vx*cx + 2*vy*cy + 2*vz*cz + k = vx² + vy² + vz²
    // where k = r² - |c|²
    // Normal equations: A^T A x = A^T b
    // A = [[2*vx, 2*vy, 2*vz, 1], ...], b = [vx²+vy²+vz², ...]
    let mut ata = [[0.0_f64; 4]; 4];
    let mut atb = [0.0_f64; 4];

    for &vi in vertex_set {
        let v = &vertices[vi];
        let row = [2.0 * v.x, 2.0 * v.y, 2.0 * v.z, 1.0];
        let b_i = v.x * v.x + v.y * v.y + v.z * v.z;

        for i in 0..4 {
            for j in 0..4 {
                ata[i][j] += row[i] * row[j];
            }
            atb[i] += row[i] * b_i;
        }
    }

    // Solve 4x4 system using Gaussian elimination with partial pivoting
    let mut aug = [[0.0_f64; 5]; 4];
    for i in 0..4 {
        for j in 0..4 {
            aug[i][j] = ata[i][j];
        }
        aug[i][4] = atb[i];
    }

    for col in 0..4 {
        // Find pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..4 {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            return None; // Singular
        }
        if max_row != col {
            aug.swap(col, max_row);
        }
        let pivot = aug[col][col];
        for j in col..5 {
            aug[col][j] /= pivot;
        }
        for row in 0..4 {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            for j in col..5 {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    let cx = aug[0][4];
    let cy = aug[1][4];
    let cz = aug[2][4];
    let k = aug[3][4];

    let r_sq = k + cx * cx + cy * cy + cz * cz;
    if r_sq <= 0.0 {
        return None;
    }

    let mut center = [cx, cy, cz];
    let mut radius = r_sq.sqrt();

    // Gauss-Newton refinement minimizing vertex-to-sphere distances.
    // Parameters: (cx, cy, cz, r).
    // f_i = ||p_i - c|| - r
    // Jacobian: ∂f/∂c = -(p_i - c) / ||p_i - c||, ∂f/∂r = -1
    // Solves J^T J delta = -J^T f via normal equations at each iteration.
    if n >= 6 {
        let verts: Vec<usize> = vertex_set.iter().copied().collect();
        for _iter in 0..20 {
            // Build normal equations: (J^T J) delta = -(J^T f)
            let mut jtj = [[0.0_f64; 4]; 4];
            let mut jtf = [0.0_f64; 4];

            for &vi in &verts {
                let v = &vertices[vi];
                let diff = [v.x - center[0], v.y - center[1], v.z - center[2]];
                let dist = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
                let fi = dist - radius;

                let ji = if dist < 1e-12 {
                    [0.0, 0.0, 0.0, -1.0]
                } else {
                    [-diff[0] / dist, -diff[1] / dist, -diff[2] / dist, -1.0]
                };

                for r in 0..4 {
                    for c in 0..4 {
                        jtj[r][c] += ji[r] * ji[c];
                    }
                    jtf[r] += ji[r] * fi;
                }
            }

            // Solve 4x4 system: jtj * delta = -jtf
            let mut aug_gn = [[0.0_f64; 5]; 4];
            for i in 0..4 {
                for j in 0..4 {
                    aug_gn[i][j] = jtj[i][j];
                }
                aug_gn[i][4] = -jtf[i];
            }

            let mut singular = false;
            for col in 0..4 {
                let mut max_row = col;
                let mut max_val = aug_gn[col][col].abs();
                for row in (col + 1)..4 {
                    if aug_gn[row][col].abs() > max_val {
                        max_val = aug_gn[row][col].abs();
                        max_row = row;
                    }
                }
                if max_val < 1e-30 {
                    singular = true;
                    break;
                }
                if max_row != col {
                    aug_gn.swap(col, max_row);
                }
                let pivot = aug_gn[col][col];
                for j in col..5 {
                    aug_gn[col][j] /= pivot;
                }
                for row in 0..4 {
                    if row == col {
                        continue;
                    }
                    let factor = aug_gn[row][col];
                    for j in col..5 {
                        aug_gn[row][j] -= factor * aug_gn[col][j];
                    }
                }
            }
            if singular {
                break;
            }

            let delta = [aug_gn[0][4], aug_gn[1][4], aug_gn[2][4], aug_gn[3][4]];
            center[0] += delta[0];
            center[1] += delta[1];
            center[2] += delta[2];
            radius += delta[3];

            let delta_norm = (delta[0] * delta[0] + delta[1] * delta[1]
                + delta[2] * delta[2] + delta[3] * delta[3])
                .sqrt();
            if delta_norm < 1e-14 {
                break;
            }
        }
    }

    if radius <= 0.0 {
        return None;
    }

    Some((center, radius))
}

/// Determine convexity for a sphere: does the face normal point away from the center?
fn determine_sphere_convexity(
    face: &MeshFace,
    vertices: &[MeshVertex],
    center: &[f64; 3],
) -> bool {
    let centroid = face_centroid(face, vertices);
    let radial = [
        centroid[0] - center[0],
        centroid[1] - center[1],
        centroid[2] - center[2],
    ];
    let n = face.normal.unwrap();
    dot3(&n, &radial) > 0.0
}

/// Validate solid-angle coverage of a spherical hypothesis.
/// Returns true if the centroid-to-center direction vectors span 3D
/// (not a 1D strip like a cylinder fillet).
fn solid_angle_coverage_valid(
    face_list: &[usize],
    center: &[f64; 3],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> bool {
    if face_list.len() < MIN_SPHERE_FACES {
        return false;
    }

    // Compute area-weighted covariance of unit direction vectors from center to centroid.
    let mut cov = [0.0_f64; 6]; // upper triangle: [c00, c01, c02, c11, c12, c22]
    let mut total_weight = 0.0_f64;
    for &fi in face_list {
        let c = face_centroid(&faces[fi], vertices);
        let mut d = [
            c[0] - center[0],
            c[1] - center[1],
            c[2] - center[2],
        ];
        let len = normalize3(&mut d);
        if len < 1e-15 {
            continue;
        }
        let w = face_area(&faces[fi], vertices);
        total_weight += w;
        cov[0] += w * d[0] * d[0];
        cov[1] += w * d[0] * d[1];
        cov[2] += w * d[0] * d[2];
        cov[3] += w * d[1] * d[1];
        cov[4] += w * d[1] * d[2];
        cov[5] += w * d[2] * d[2];
    }
    if total_weight < 1e-30 {
        return false;
    }
    // Normalize
    for v in &mut cov {
        *v /= total_weight;
    }

    // Compute eigenvalues of the 3x3 symmetric matrix.
    let eigenvalues = eigenvalues_3x3_symmetric(&cov);
    let lambda_max = eigenvalues[0].max(eigenvalues[1]).max(eigenvalues[2]);
    let lambda_min = eigenvalues[0].min(eigenvalues[1]).min(eigenvalues[2]);

    if lambda_max < 1e-30 {
        return false;
    }
    lambda_min / lambda_max >= MIN_SPHERE_EIGENVALUE_RATIO
}

// ---------------------------------------------------------------------------
// Stage 2.3: Compare check
// ---------------------------------------------------------------------------

/// Validate fitted spherical hypotheses against a reference STEP shape.
pub(super) fn compare_spherical_hypotheses(
    hypotheses: &[SphericalHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);

            // Project centroid onto sphere surface (nearest point)
            let d = [
                centroid[0] - hyp.center[0],
                centroid[1] - hyp.center[1],
                centroid[2] - hyp.center[2],
            ];
            let dist_to_center = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();

            let projected = if dist_to_center > 1e-15 {
                let scale = hyp.radius / dist_to_center;
                [
                    hyp.center[0] + d[0] * scale,
                    hyp.center[1] + d[1] * scale,
                    hyp.center[2] + d[2] * scale,
                ]
            } else {
                // Centroid is at the center — unlikely but handle gracefully
                centroid
            };

            let pt = gp::Pnt::new_real3(projected[0], projected[1], projected[2]);
            let dist = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(dist);
        }

        if max_dist > config.surface_tolerance_mm {
            return Err(Stage2CompareError {
                hypothesis_type: "spherical",
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

/// Fitted sphere parameters.
#[derive(Clone)]
pub(super) struct SphereFit {
    pub center: [f64; 3],
    pub radius: f64,
    pub convex: Option<bool>,
}

/// Marker struct for the spherical surface type.
pub(super) struct SphericalSurface;

impl SurfaceType for SphericalSurface {
    type Fit = SphereFit;
    type Hypothesis = SphericalHypothesis;

    fn name() -> &'static str { "spherical" }

    fn hypothesis_field(face: &MeshFace) -> i32 { face.spherical_hypothesis }

    fn set_hypothesis_field(face: &mut MeshFace, value: i32) {
        face.spherical_hypothesis = value;
    }

    fn undeduced_value() -> i32 { UNDEDUCED_SPHERICAL_HYPOTHESIS }

    fn face_eligible(face: &MeshFace) -> bool {
        face.spherical_hypothesis == UNDEDUCED_SPHERICAL_HYPOTHESIS
    }

    fn generate_seeds(mesh: &ConnectedMesh, fi: usize, _config: &BfsConfig) -> Vec<Seed> {
        if mesh.faces[fi].normal.is_none() {
            return Vec::new();
        }
        // 1-ring: collect all faces sharing any vertex with fi
        let vc = mesh.faces[fi].vertex_count as usize;
        let mut seed_set: HashSet<usize> = HashSet::new();
        seed_set.insert(fi);
        for &vi in &mesh.faces[fi].vertex_indices[..vc] {
            // Find all faces incident on this vertex
            for (fj, face) in mesh.faces.iter().enumerate() {
                let fj_vc = face.vertex_count as usize;
                for &fj_vi in &face.vertex_indices[..fj_vc] {
                    if fj_vi == vi && face.normal.is_some() {
                        if face.spherical_hypothesis == UNDEDUCED_SPHERICAL_HYPOTHESIS {
                            seed_set.insert(fj);
                        }
                        break;
                    }
                }
            }
        }
        if seed_set.len() < 3 {
            return Vec::new();
        }
        let faces: Vec<usize> = seed_set.into_iter().collect();
        vec![Seed { faces }]
    }

    fn fit_surface(
        mesh: &ConnectedMesh,
        _face_list: &[usize],
        vertex_list: &[usize],
        _verbosity: u8,
    ) -> Option<Self::Fit> {
        let vertex_set: HashSet<usize> = vertex_list.iter().copied().collect();
        let (center, radius) = fit_sphere(&vertex_set, &mesh.vertices)?;
        // Reject obviously wrong radii
        let bb_diag = bounding_box_diagonal(&mesh.vertices);
        if radius > bb_diag * MAX_SPHERE_RADIUS_FACTOR {
            return None;
        }
        Some(SphereFit { center, radius, convex: None })
    }

    fn vertex_distance(fit: &Self::Fit, pos: &[f64; 3]) -> f64 {
        let v = MeshVertex::from_xyz(pos[0], pos[1], pos[2]);
        vertex_to_sphere_distance(&v, &fit.center, fit.radius).abs()
    }

    fn convexity_check(fit: &Self::Fit, mesh: &ConnectedMesh, fi: usize) -> bool {
        if let Some(convex) = fit.convex {
            let fi_convex = determine_sphere_convexity(
                &mesh.faces[fi], &mesh.vertices, &fit.center,
            );
            fi_convex == convex
        } else {
            true
        }
    }

    fn determine_convexity(mesh: &ConnectedMesh, face_list: &[usize], fit: &mut Self::Fit) {
        if let Some(&first) = face_list.first() {
            fit.convex = Some(determine_sphere_convexity(
                &mesh.faces[first], &mesh.vertices, &fit.center,
            ));
        }
    }

    fn acceptance_test(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        fit: &Self::Fit,
        _config: &BfsConfig,
    ) -> bool {
        if face_list.len() < MIN_SPHERE_FACES {
            return false;
        }
        solid_angle_coverage_valid(
            face_list, &fit.center, &mesh.faces, &mesh.vertices,
        )
    }

    fn viz_overlays(
        fit: &Self::Fit,
        _mesh: &ConnectedMesh,
        _face_list: &[usize],
    ) -> (Vec<viz::CylinderOverlay>, Vec<viz::SphereOverlay>) {
        let sphere = viz::SphereOverlay {
            center: fit.center,
            radius: fit.radius,
            color: [0.2, 0.4, 1.0, 0.15],
        };
        (Vec::new(), vec![sphere])
    }

    fn viz_fit_description(fit: &Self::Fit) -> String {
        format!("r={:.4}", fit.radius)
    }

    fn make_hypothesis(
        mesh: &ConnectedMesh,
        faces: Vec<usize>,
        vertices: Vec<usize>,
        fit: &Self::Fit,
    ) -> SphericalHypothesis {
        let mut error_max = 0.0_f64;
        let mut centroid_error_max = 0.0_f64;
        let mut error_abs_sum = 0.0_f64;
        for &vi in &vertices {
            let v = &mesh.vertices[vi];
            let d = vertex_to_sphere_distance(v, &fit.center, fit.radius).abs();
            error_max = error_max.max(d);
            error_abs_sum += d;
        }
        for &fi in &faces {
            let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
            let v = MeshVertex::from_xyz(c[0], c[1], c[2]);
            let d = vertex_to_sphere_distance(&v, &fit.center, fit.radius).abs();
            centroid_error_max = centroid_error_max.max(d);
        }
        SphericalHypothesis {
            center: fit.center,
            radius: fit.radius,
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

    /// Generate points on a sphere.
    fn sphere_points(
        center: [f64; 3],
        radius: f64,
        n_rings: usize,
        n_pts_per_ring: usize,
        elev_min: f64,
        elev_max: f64,
    ) -> Vec<MeshVertex> {
        let mut pts = Vec::new();
        for ring in 0..n_rings {
            let t = if n_rings == 1 { 0.5 } else { ring as f64 / (n_rings - 1) as f64 };
            let elev = elev_min + t * (elev_max - elev_min);
            let (sin_e, cos_e) = elev.sin_cos();
            for j in 0..n_pts_per_ring {
                let theta = 2.0 * PI * j as f64 / n_pts_per_ring as f64;
                let (sin_t, cos_t) = theta.sin_cos();
                pts.push(MeshVertex::from_xyz(
                    center[0] + radius * cos_e * cos_t,
                    center[1] + radius * cos_e * sin_t,
                    center[2] + radius * sin_e,
                ));
            }
        }
        pts
    }

    // ---- Distance tests ----

    #[test]
    fn distance_to_sphere_on_surface() {
        let v = MeshVertex::from_xyz(10.0, 0.0, 0.0);
        let d = vertex_to_sphere_distance(&v, &[0.0, 0.0, 0.0], 10.0);
        assert!(d.abs() < 1e-10, "on-surface distance should be 0, got {d}");
    }

    #[test]
    fn distance_to_sphere_inside() {
        let v = MeshVertex::from_xyz(5.0, 0.0, 0.0);
        let d = vertex_to_sphere_distance(&v, &[0.0, 0.0, 0.0], 10.0);
        assert!((d - (-5.0)).abs() < 1e-10, "inside distance should be -5, got {d}");
    }

    #[test]
    fn distance_to_sphere_outside() {
        let v = MeshVertex::from_xyz(15.0, 0.0, 0.0);
        let d = vertex_to_sphere_distance(&v, &[0.0, 0.0, 0.0], 10.0);
        assert!((d - 5.0).abs() < 1e-10, "outside distance should be 5, got {d}");
    }

    #[test]
    fn distance_to_sphere_at_center() {
        let v = MeshVertex::from_xyz(0.0, 0.0, 0.0);
        let d = vertex_to_sphere_distance(&v, &[0.0, 0.0, 0.0], 10.0);
        assert!((d - (-10.0)).abs() < 1e-10, "at-center distance should be -radius, got {d}");
    }

    #[test]
    fn distance_to_sphere_arbitrary() {
        let v = MeshVertex::from_xyz(3.0, 4.0, 0.0);
        let d = vertex_to_sphere_distance(&v, &[0.0, 0.0, 0.0], 10.0);
        assert!((d - (-5.0)).abs() < 1e-10, "distance should be -5, got {d}");
    }

    #[test]
    fn distance_to_sphere_offset_center() {
        let v = MeshVertex::from_xyz(1.0, 2.0, 10.0);
        let d = vertex_to_sphere_distance(&v, &[1.0, 2.0, 3.0], 7.0);
        assert!(d.abs() < 1e-10, "on-surface distance should be 0, got {d}");
    }

    // ---- Fit tests ----

    #[test]
    fn sphere_fit_origin() {
        let vertices = sphere_points([0.0, 0.0, 0.0], 12.5, 8, 16, -PI / 2.0 + 0.1, PI / 2.0 - 0.1);
        let vertex_set: HashSet<usize> = (0..vertices.len()).collect();

        let result = fit_sphere(&vertex_set, &vertices);
        assert!(result.is_some(), "sphere fit should succeed");
        let (center, radius) = result.unwrap();

        assert!((radius - 12.5).abs() < 0.01, "radius should be 12.5, got {radius}");
        assert!(center[0].abs() < 0.01, "center.x should be 0, got {}", center[0]);
        assert!(center[1].abs() < 0.01, "center.y should be 0, got {}", center[1]);
        assert!(center[2].abs() < 0.01, "center.z should be 0, got {}", center[2]);
    }

    #[test]
    fn sphere_fit_offset() {
        let vertices = sphere_points([5.0, -3.0, 7.0], 10.0, 10, 20, -PI / 2.0 + 0.1, PI / 2.0 - 0.1);
        let vertex_set: HashSet<usize> = (0..vertices.len()).collect();

        let result = fit_sphere(&vertex_set, &vertices);
        assert!(result.is_some(), "offset sphere fit should succeed");
        let (center, radius) = result.unwrap();

        assert!((radius - 10.0).abs() < 0.01, "radius should be 10, got {radius}");
        assert!((center[0] - 5.0).abs() < 0.01, "center.x should be 5, got {}", center[0]);
        assert!((center[1] - (-3.0)).abs() < 0.01, "center.y should be -3, got {}", center[1]);
        assert!((center[2] - 7.0).abs() < 0.01, "center.z should be 7, got {}", center[2]);
    }

    #[test]
    fn sphere_fit_hemisphere() {
        let vertices = sphere_points([0.0, 0.0, 0.0], 10.0, 6, 16, 0.05, PI / 2.0 - 0.05);
        let vertex_set: HashSet<usize> = (0..vertices.len()).collect();

        let result = fit_sphere(&vertex_set, &vertices);
        assert!(result.is_some(), "hemisphere fit should succeed");
        let (center, radius) = result.unwrap();

        assert!((radius - 10.0).abs() < 0.1, "radius should be 10, got {radius}");
        assert!(center[0].abs() < 0.1, "center.x should be 0, got {}", center[0]);
        assert!(center[1].abs() < 0.1, "center.y should be 0, got {}", center[1]);
        assert!(center[2].abs() < 0.1, "center.z should be 0, got {}", center[2]);
    }

    #[test]
    fn sphere_fit_small_cap() {
        // Only a small cap near the north pole: elevation from 60° to 90°.
        let vertices = sphere_points([0.0, 0.0, 0.0], 25.0, 5, 12, PI / 3.0, PI / 2.0 - 0.05);
        let vertex_set: HashSet<usize> = (0..vertices.len()).collect();

        let result = fit_sphere(&vertex_set, &vertices);
        assert!(result.is_some(), "small cap fit should succeed");
        let (center, radius) = result.unwrap();

        assert!((radius - 25.0).abs() < 0.5, "radius should be 25, got {radius}");
        assert!(center[0].abs() < 0.5, "center.x should be 0, got {}", center[0]);
        assert!(center[1].abs() < 0.5, "center.y should be 0, got {}", center[1]);
        assert!(center[2].abs() < 0.5, "center.z should be 0, got {}", center[2]);
    }

    #[test]
    fn sphere_fit_rejects_degenerate() {
        let vertices = vec![
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0),
            MeshVertex::from_xyz(0.0, 0.0, 1.0),
        ];
        let vertex_set: HashSet<usize> = (0..3).collect();
        assert!(fit_sphere(&vertex_set, &vertices).is_none(),
                "sphere fit should reject 3 vertices");
    }

    #[test]
    fn sphere_fit_minimum_points() {
        // 4 points on a unit sphere (regular tetrahedron vertices)
        let s = 1.0_f64 / 3.0_f64.sqrt();
        let vertices = vec![
            MeshVertex::from_xyz(s, s, s),
            MeshVertex::from_xyz(s, -s, -s),
            MeshVertex::from_xyz(-s, s, -s),
            MeshVertex::from_xyz(-s, -s, s),
        ];
        let vertex_set: HashSet<usize> = (0..4).collect();

        let result = fit_sphere(&vertex_set, &vertices);
        assert!(result.is_some(), "4-point sphere fit should succeed");
        let (center, radius) = result.unwrap();

        assert!((radius - 1.0).abs() < 1e-10, "radius should be 1, got {radius}");
        assert!(center[0].abs() < 1e-10, "center.x should be 0, got {}", center[0]);
        assert!(center[1].abs() < 1e-10, "center.y should be 0, got {}", center[1]);
        assert!(center[2].abs() < 1e-10, "center.z should be 0, got {}", center[2]);
    }
}

