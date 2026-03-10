//! Stage 2.2: Cylindrical hypothesis deduction and comparison.

use std::collections::{HashSet, VecDeque};

use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use nalgebra::{Const, Dyn, OMatrix, OVector, Owned};
use opencascade_sys::gp;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_CYLINDRICAL_HYPOTHESIS, NO_HYPOTHESIS};
use crate::viz::{self, VizAction, VizSender};

use super::{
    build_normal_covariance, centered_cylinder_overlay, cross3, dot3, face_area, face_centroid,
    normalize3, perpendicular_basis, smallest_eigenvector_3x3, viz_custom, viz_face_centroid,
    viz_face_normal, CylindricalHypothesis, Stage2CompareError, MIN_CROSS_THRESHOLD,
    MIN_CYLINDER_FACES, REFIT_SKIP_MULTIPLIER,
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
    vertex_set: &HashSet<usize>,
    vertices: &[MeshVertex],
) -> Option<([f64; 3], f64)> {
    let (u, w) = perpendicular_basis(axis_dir);
    let verts: Vec<usize> = vertex_set.iter().copied().collect();
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
fn fit_cylinder(
    face_indices: &[usize],
    vertex_set: &HashSet<usize>,
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

    // Refine with Levenberg-Marquardt optimization on radial residuals.
    let points: Vec<[f64; 3]> = vertex_set.iter()
        .map(|&vi| [vertices[vi].x, vertices[vi].y, vertices[vi].z])
        .collect();
    let problem = CylinderLMProblem::new(points, axis_dir, axis_origin, radius);
    let (result, _report) = LevenbergMarquardt::new().minimize(problem);
    let refined_dir = result.axis_dir_from_params(&result.params);
    let refined_origin = [result.params[2], result.params[3], result.params[4]];
    let refined_radius = result.params[5];

    if refined_radius > 0.0 {
        Some((refined_origin, refined_dir, refined_radius))
    } else {
        Some((axis_origin, axis_dir, radius))
    }
}

/// Check if all vertices in a set are within tolerance of a cylinder surface.
fn all_vertices_within_cylinder_tolerance(
    vertex_set: &HashSet<usize>,
    axis_origin: &[f64; 3],
    axis_direction: &[f64; 3],
    radius: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        let d = vertex_to_cylinder_distance(&vertices[vi], axis_origin, axis_direction, radius);
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
    axis_origin: &[f64; 3],
    axis_direction: &[f64; 3],
) -> bool {
    let centroid = face_centroid(face, vertices);
    // Vector from axis to centroid (radial component)
    let d = [
        centroid[0] - axis_origin[0],
        centroid[1] - axis_origin[1],
        centroid[2] - axis_origin[2],
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

/// Validate that faces in a cylindrical hypothesis have sufficient angular
/// coverage around the cylinder — i.e., they are distributed around the
/// circumference, not clustered on one side.
///
/// Requires at least 3 distinct angular clusters. The check is:
/// - Compute the angular coordinate θ of each face centroid around the axis.
/// - Sort θ values and compute all N inter-face gaps plus the wraparound gap.
/// - The largest gap defines the "empty arc"; span = 2π − largest_gap.
/// - The second-largest gap must be ≤ span / 3.
fn angular_coverage_valid(
    face_list: &[usize],
    axis_origin: &[f64; 3],
    axis_direction: &[f64; 3],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> bool {
    if face_list.len() < 3 {
        return false;
    }

    let (u, w) = perpendicular_basis(axis_direction);

    // Compute angular coordinate for each face centroid
    let mut thetas: Vec<f64> = Vec::with_capacity(face_list.len());
    for &fi in face_list {
        let c = face_centroid(&faces[fi], vertices);
        let d = [
            c[0] - axis_origin[0],
            c[1] - axis_origin[1],
            c[2] - axis_origin[2],
        ];
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

    // Compute all gaps (including wraparound)
    let mut gaps: Vec<f64> = Vec::with_capacity(n);
    for i in 1..n {
        gaps.push(thetas[i] - thetas[i - 1]);
    }
    // Wraparound gap
    gaps.push(thetas[0] + 2.0 * std::f64::consts::PI - thetas[n - 1]);

    // Find largest and second-largest gaps
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

/// Result of a trial BFS for cylindrical hypothesis evaluation.
struct CylinderTrialResult {
    faces: Vec<usize>,
    vertices: HashSet<usize>,
    axis_origin: [f64; 3],
    axis_direction: [f64; 3],
    radius: f64,
    convex: bool,
    error_max: f64,
    centroid_error_max: f64,
    error_abs_sum: f64,
    total_area: f64,
}

/// Run a trial BFS for a cylindrical hypothesis starting from seed faces,
/// using temporary data structures (not mutating mesh face assignments).
/// Returns `None` if the trial fails validation (min faces, centroid check).
/// Sets `viz_quit` if the user quits from viz.
#[allow(clippy::too_many_arguments)]
fn run_cylinder_trial_bfs(
    seed_faces: &[usize],
    mesh: &ConnectedMesh,
    vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
    verbosity: u8,
    viz: Option<&VizSender>,
    viz_quit: &std::cell::Cell<bool>,
) -> Option<CylinderTrialResult> {
    // Collect seed vertices from all seed faces
    let mut vertex_set: HashSet<usize> = HashSet::new();
    for &sfi in seed_faces {
        let vc = mesh.faces[sfi].vertex_count as usize;
        for &vi in &mesh.faces[sfi].vertex_indices[..vc] {
            vertex_set.insert(vi);
        }
    }

    let mut face_list: Vec<usize> = seed_faces.to_vec();

    // Fit cylinder to seed faces
    let (mut current_origin, mut current_dir, mut current_radius) = match
        fit_cylinder(&face_list, &vertex_set, &mesh.faces, &mesh.vertices)
    {
        Some(v) => v,
        None => {
            if verbosity >= 3 {
                let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
                eprintln!("  [BFS-cyl] seed=({}) fit_cylinder returned None (degenerate geometry)", seed_str.join(","));
            }
            return None;
        }
    };

    // Verify seed: all vertices within tolerance
    if !all_vertices_within_cylinder_tolerance(
        &vertex_set, &current_origin, &current_dir, current_radius,
        vertex_tol, &mesh.vertices,
    ) {
        if verbosity >= 3 {
            let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
            let worst_dist = vertex_set.iter().map(|&vi| {
                vertex_to_cylinder_distance(&mesh.vertices[vi], &current_origin, &current_dir, current_radius).abs()
            }).fold(0.0_f64, f64::max);
            eprintln!("  [BFS-cyl] seed=({}) vertex tolerance failed: worst_dist={:.2e} > tol={:.2e}, r={:.6}",
                seed_str.join(","), worst_dist, vertex_tol, current_radius);
        }
        return None;
    }

    // Determine convexity from first seed face
    let convex = determine_convexity(
        &mesh.faces[seed_faces[0]], &mesh.vertices, &current_origin, &current_dir,
    );

    // Verify convexity consistency across all seed faces
    for &sfi in &seed_faces[1..] {
        if determine_convexity(&mesh.faces[sfi], &mesh.vertices, &current_origin, &current_dir) != convex {
            if verbosity >= 3 {
                let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
                eprintln!("  [BFS-cyl] seed=({}) convexity inconsistency: face {} is {} but face {} is {}",
                    seed_str.join(","), seed_faces[0], if convex { "convex" } else { "concave" },
                    sfi, if convex { "concave" } else { "convex" });
            }
            return None;
        }
    }

    // Verify all seed face centroids within surface tolerance
    for &sfi in seed_faces {
        let centroid = face_centroid(&mesh.faces[sfi], &mesh.vertices);
        let cen_dist = vertex_to_cylinder_distance(
            &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
            &current_origin, &current_dir, current_radius,
        ).abs();
        if cen_dist > surface_tol {
            if verbosity >= 3 {
                let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
                eprintln!("  [BFS-cyl] seed=({}) centroid tolerance failed: face {} cen_dist={:.2e} > tol={:.2e}",
                    seed_str.join(","), sfi, cen_dist, surface_tol);
            }
            return None;
        }
    }

    // Level-3 trace: print seed info
    if verbosity >= 3 {
        let print_seed_face = |label: &str, sfi: usize| {
            let svc = mesh.faces[sfi].vertex_count as usize;
            let vis: Vec<usize> = mesh.faces[sfi].vertex_indices[..svc].to_vec();
            let coords: Vec<String> = vis.iter().map(|&vi| {
                let v = &mesh.vertices[vi];
                format!("vi={}:[{:.4},{:.4},{:.4}]", vi, v.x, v.y, v.z)
            }).collect();
            let centroid = face_centroid(&mesh.faces[sfi], &mesh.vertices);
            let cen_err = vertex_to_cylinder_distance(
                &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                &current_origin, &current_dir, current_radius,
            ).abs();
            eprintln!(
                "  {} fi={}: {}v centroid=[{:.4},{:.4},{:.4}] cen_err={:.2e}  vertices: {}",
                label, sfi, svc,
                centroid[0], centroid[1], centroid[2], cen_err,
                coords.join(" "),
            );
        };
        let seed_str: Vec<String> = seed_faces.iter().map(|f| f.to_string()).collect();
        eprintln!(
            "[BFS-cyl] Trial seed=({}) \u{2192} r={:.6}, origin=[{:.4},{:.4},{:.4}], dir=[{:.4},{:.4},{:.4}], {}",
            seed_str.join(","),
            current_radius,
            current_origin[0], current_origin[1], current_origin[2],
            current_dir[0], current_dir[1], current_dir[2],
            if convex { "convex" } else { "concave" },
        );
        for (i, &sfi) in seed_faces.iter().enumerate() {
            print_seed_face(&format!("seed_{}", i), sfi);
        }
    }

    // Viz: show seed faces (first face=orange, others=light orange)
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
            &format!("BFS-cyl: seed=({}) r={:.4} {} [space=step, shift+space=skip]",
                seed_str.join(","), current_radius, if convex { "convex" } else { "concave" }),
            vec![centered_cylinder_overlay(
                current_origin, current_dir, current_radius,
                seed_faces, mesh,
                [0.2, 0.4, 1.0, 0.3],
            )], Vec::new(),
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

    // Track claimed faces in this trial (using a HashSet, not mesh mutation)
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

            // Skip if already committed to a real hypothesis or claimed by this trial
            if mesh.faces[cni].cylindrical_hypothesis != UNDEDUCED_CYLINDRICAL_HYPOTHESIS {
                if verbosity >= 4 {
                    eprintln!("  [BFS-cyl] from fi={} try cni={}: already assigned (hyp={}) → SKIP", current_fi, cni, mesh.faces[cni].cylindrical_hypothesis);
                }
                continue;
            }
            if trial_claimed.contains(&cni) {
                continue;
            }

            // Convexity check
            let cni_convex = determine_convexity(
                &mesh.faces[cni], &mesh.vertices, &current_origin, &current_dir,
            );
            if cni_convex != convex {
                if verbosity >= 3 {
                    eprintln!(
                        "  [BFS-cyl] from fi={} try cni={}: convexity mismatch ({} vs {}) → REJECT(convexity)",
                        current_fi, cni,
                        if cni_convex { "convex" } else { "concave" },
                        if convex { "convex" } else { "concave" },
                    );
                }
                continue;
            }

            // Angular tolerance: reject if dihedral angle with ANY already-claimed
            // neighbor exceeds the limit
            if let Some(n_cni) = mesh.faces[cni].normal {
                let cni_vc2 = mesh.faces[cni].vertex_count as usize;
                let cni_neighbors = mesh.faces[cni].neighbors;
                let mut angular_reject = false;
                let mut worst_angle_deg = 0.0_f64;
                let mut worst_adj = usize::MAX;
                for &adj in &cni_neighbors[..cni_vc2] {
                    if adj < 0 { continue; }
                    let adj = adj as usize;
                    // Check against both committed (same hypothesis doesn't exist yet)
                    // and trial-claimed faces
                    if !trial_claimed.contains(&adj) { continue; }
                    if let Some(n_adj) = mesh.faces[adj].normal {
                        let cos_a = dot3(&n_adj, &n_cni).clamp(-1.0, 1.0);
                        let angle = cos_a.acos();
                        if angle > worst_angle_deg.to_radians() {
                            worst_angle_deg = angle.to_degrees();
                            worst_adj = adj;
                        }
                        if angle > angular_tol {
                            angular_reject = true;
                            break;
                        }
                    }
                }
                if angular_reject {
                    if verbosity >= 3 {
                        eprintln!(
                            "  [BFS-cyl] from fi={} try cni={}: angular {:.2}° > tol {:.2}° (adj fi={}) → REJECT(angular)",
                            current_fi, cni, worst_angle_deg, angular_tol.to_degrees(), worst_adj,
                        );
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
                let d = vertex_to_cylinder_distance(
                    &mesh.vertices[vi],
                    &current_origin,
                    &current_dir,
                    current_radius,
                ).abs();
                vtx_err_max = vtx_err_max.max(d);
                if d > vertex_tol {
                    all_ok = false;
                    if d > REFIT_SKIP_MULTIPLIER * vertex_tol {
                        any_far = true;
                        break;
                    }
                }
            }

            if any_far {
                if verbosity >= 3 {
                    eprintln!(
                        "  [BFS-cyl] from fi={} try cni={}: vtx_err_max={:.2e} > REFIT_SKIP*tol={:.2e} → REJECT(too far)",
                        current_fi, cni, vtx_err_max, REFIT_SKIP_MULTIPLIER * vertex_tol,
                    );
                }
                continue;
            }

            // Centroid validation: check face centroid distance to cylinder
            let centroid = face_centroid(&mesh.faces[cni], &mesh.vertices);
            let centroid_dist = vertex_to_cylinder_distance(
                &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                &current_origin, &current_dir, current_radius,
            ).abs();
            if centroid_dist > REFIT_SKIP_MULTIPLIER * surface_tol {
                if verbosity >= 3 {
                    eprintln!(
                        "  [BFS-cyl] from fi={} try cni={}: cen_err={:.2e} > REFIT_SKIP*stol={:.2e} → REJECT(centroid far)",
                        current_fi, cni, centroid_dist, REFIT_SKIP_MULTIPLIER * surface_tol,
                    );
                }
                continue;
            }
            let needs_refit = !all_ok || centroid_dist > surface_tol;

            if needs_refit {
                // Try re-fitting with new face included
                let mut trial_vertices = vertex_set.clone();
                for &vi in &cni_vi[..cni_vc] {
                    trial_vertices.insert(vi);
                }
                let mut trial_faces = face_list.clone();
                trial_faces.push(cni);

                let refit = fit_cylinder(
                    &trial_faces, &trial_vertices, &mesh.faces, &mesh.vertices,
                );
                let (new_origin, new_dir, new_radius) = match refit {
                    Some(params) => params,
                    None => {
                        if verbosity >= 3 {
                            eprintln!(
                                "  [BFS-cyl] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit degenerate) → REJECT(refit failed)",
                                current_fi, cni, vtx_err_max, centroid_dist,
                            );
                        }
                        continue;
                    }
                };

                if !all_vertices_within_cylinder_tolerance(
                    &trial_vertices, &new_origin, &new_dir, new_radius,
                    vertex_tol, &mesh.vertices,
                ) {
                    if verbosity >= 3 {
                        eprintln!(
                            "  [BFS-cyl] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit vertex check failed) → REJECT(refit tol)",
                            current_fi, cni, vtx_err_max, centroid_dist,
                        );
                    }
                    continue;
                }

                // Check all face centroids within surface tolerance after re-fit
                let mut centroids_ok = true;
                for &f in &trial_faces {
                    let c = face_centroid(&mesh.faces[f], &mesh.vertices);
                    let d = vertex_to_cylinder_distance(
                        &MeshVertex::from_xyz(c[0], c[1], c[2]),
                        &new_origin, &new_dir, new_radius,
                    ).abs();
                    if d > surface_tol {
                        centroids_ok = false;
                        break;
                    }
                }
                if !centroids_ok {
                    if verbosity >= 3 {
                        eprintln!(
                            "  [BFS-cyl] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit centroid check failed) → REJECT(refit centroid)",
                            current_fi, cni, vtx_err_max, centroid_dist,
                        );
                    }
                    continue;
                }

                if verbosity >= 3 {
                    eprintln!(
                        "  [BFS-cyl] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit ok → r={:.6} origin=[{:.4},{:.4},{:.4}]) → ACCEPT[refit]",
                        current_fi, cni, vtx_err_max, centroid_dist,
                        new_radius, new_origin[0], new_origin[1], new_origin[2],
                    );
                }

                // Accept re-fit
                current_origin = new_origin;
                current_dir = new_dir;
                current_radius = new_radius;
            } else if verbosity >= 3 {
                eprintln!(
                    "  [BFS-cyl] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} → ACCEPT",
                    current_fi, cni, vtx_err_max, centroid_dist,
                );
            }

            // Accept this face
            trial_claimed.insert(cni);
            face_list.push(cni);
            for &vi in &cni_vi[..cni_vc] {
                vertex_set.insert(vi);
            }
            queue.push_back(cni);

            // Viz: show accepted face with custom colors
            if !skip_viz {
                let queue_faces: Vec<usize> = queue.iter().copied().collect();
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
                let edge_highlights = if !queue_faces.is_empty() {
                    vec![viz::EdgeHighlight { face_indices: queue_faces, color: [1.0, 0.0, 0.0, 1.0] }]
                } else {
                    Vec::new()
                };
                // Add gray background for faces with existing cylindrical hypotheses
                let bg_faces: Vec<usize> = (0..mesh.faces.len())
                    .filter(|f| mesh.faces[*f].cylindrical_hypothesis >= 0 && !face_list.contains(f))
                    .collect();
                if !bg_faces.is_empty() {
                    highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                }
                if let Some(action) = viz_custom(
                    viz, highlights,
                    edge_highlights,
                    &format!("BFS-cyl: accepted fi={cni} ({} faces) r={:.4} [space=step, shift+space=skip]", face_list.len(), current_radius),
                    vec![centered_cylinder_overlay(
                        current_origin, current_dir, current_radius,
                        &face_list, mesh,
                        [0.2, 0.4, 1.0, 0.3],
                    )], Vec::new(),
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
        // (Early termination disabled — future work per DEVELOPMENT_PLAN.md)
    }

    // Final re-fit
    if let Some((final_origin, final_dir, final_radius)) =
        fit_cylinder(&face_list, &vertex_set, &mesh.faces, &mesh.vertices)
    {
        current_origin = final_origin;
        current_dir = final_dir;
        current_radius = final_radius;
    }

    // Compute error metrics
    let mut error_max = 0.0_f64;
    let mut error_abs_sum = 0.0_f64;
    for &vi in &vertex_set {
        let d = vertex_to_cylinder_distance(
            &mesh.vertices[vi], &current_origin, &current_dir, current_radius,
        ).abs();
        error_max = error_max.max(d);
        error_abs_sum += d;
    }

    // Validate: minimum face count
    if face_list.len() < MIN_CYLINDER_FACES {
        return None;
    }

    // Validate: centroid check; also compute centroid_error_max
    let mut centroid_error_max = 0.0_f64;
    for &f in &face_list {
        let c = face_centroid(&mesh.faces[f], &mesh.vertices);
        let d = vertex_to_cylinder_distance(
            &MeshVertex::from_xyz(c[0], c[1], c[2]),
            &current_origin, &current_dir, current_radius,
        ).abs();
        centroid_error_max = centroid_error_max.max(d);
        if d > surface_tol {
            return None;
        }
    }

    // Compute total area for area-based comparison
    let total_area: f64 = face_list.iter().map(|&f| face_area(&mesh.faces[f], &mesh.vertices)).sum();

    Some(CylinderTrialResult {
        faces: face_list,
        vertices: vertex_set,
        axis_origin: current_origin,
        axis_direction: current_dir,
        radius: current_radius,
        convex,
        error_max,
        centroid_error_max,
        error_abs_sum,
        total_area,
    })
}

/// Helper: check if two faces are pairwise qualified for cylindrical seeding.
/// Both must be unassigned, neighbors, have sufficient normal difference, and
/// dihedral angle within tolerance.
fn is_pairwise_qualified(
    fi: usize,
    ni: usize,
    mesh: &ConnectedMesh,
    angular_tol: f64,
) -> bool {
    if mesh.faces[ni].cylindrical_hypothesis != UNDEDUCED_CYLINDRICAL_HYPOTHESIS {
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
    let cos_angle = dot3(&fi_normal, &ni_normal).clamp(-1.0, 1.0);
    if cross_mag < MIN_CROSS_THRESHOLD {
        return false;
    }
    cos_angle.acos() <= angular_tol
}

// ---------------------------------------------------------------------------
// Stage 2.2: Cylindrical hypothesis deduction
// ---------------------------------------------------------------------------

/// Deduce cylindrical hypotheses from the mesh using 3-face seed BFS evaluation.
///
/// For each face fi, finds triples (fi, n1, n2) where n1 is a pairwise-qualified
/// neighbor of fi, and n2 is a pairwise-qualified neighbor of n1 that is NOT a
/// neighbor of fi. Runs a trial BFS for each triple, keeps the best trial
/// (most area), and commits.
pub(super) fn deduce_cylindrical_hypotheses(
    mesh: &mut ConnectedMesh,
    vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
    verbosity: u8,
    viz: Option<&VizSender>,
) -> (Vec<CylindricalHypothesis>, bool) {
    let num_faces = mesh.faces.len();
    let mut hypotheses: Vec<CylindricalHypothesis> = Vec::new();
    let viz_quit = std::cell::Cell::new(false);

    // Pre-compute neighbor sets for fast "is fi a neighbor of n2" checks
    let neighbor_sets: Vec<HashSet<usize>> = (0..num_faces).map(|fi| {
        let vc = mesh.faces[fi].vertex_count as usize;
        mesh.faces[fi].neighbors[..vc].iter()
            .filter(|&&n| n >= 0)
            .map(|&n| n as usize)
            .collect()
    }).collect();

    #[allow(clippy::needless_range_loop)] // Need fi as index into both neighbor_sets and mesh.faces
    for fi in 0..num_faces {
        if mesh.faces[fi].cylindrical_hypothesis != UNDEDUCED_CYLINDRICAL_HYPOTHESIS {
            continue;
        }

        if mesh.faces[fi].normal.is_none() {
            mesh.faces[fi].cylindrical_hypothesis = NO_HYPOTHESIS;
            continue;
        }

        // Collect pairwise-qualified neighbors of fi
        let fi_vc = mesh.faces[fi].vertex_count as usize;
        let fi_neighbors: Vec<usize> = mesh.faces[fi].neighbors[..fi_vc]
            .iter()
            .filter(|&&n| n >= 0)
            .map(|&n| n as usize)
            .filter(|&n1| is_pairwise_qualified(fi, n1, mesh, angular_tol))
            .collect();

        // Multi-seed evaluation with 3-face seeds: fi -> n1 -> n2
        let mut best_candidate: Option<CylinderTrialResult> = None;

        for &n1 in &fi_neighbors {
            // Find n2: neighbors of n1 that are pairwise qualified with n1,
            // but NOT neighbors of fi (to ensure angular spread)
            let n1_vc = mesh.faces[n1].vertex_count as usize;
            for &n2_raw in &mesh.faces[n1].neighbors[..n1_vc] {
                if n2_raw < 0 { continue; }
                let n2 = n2_raw as usize;
                if n2 == fi { continue; }
                // n2 must NOT be a neighbor of fi
                if neighbor_sets[fi].contains(&n2) { continue; }
                // n2 must be pairwise qualified with n1
                if !is_pairwise_qualified(n1, n2, mesh, angular_tol) { continue; }

                let trial = run_cylinder_trial_bfs(
                    &[fi, n1, n2], mesh,
                    vertex_tol, surface_tol, angular_tol,
                    verbosity,
                    viz, &viz_quit,
                );

                // Reject trial if angular coverage validation fails
                if let Some(ref trial) = trial {
                    if !angular_coverage_valid(
                        &trial.faces, &trial.axis_origin, &trial.axis_direction,
                        &mesh.faces, &mesh.vertices,
                    ) {
                        if verbosity >= 3 {
                            eprintln!("  [BFS-cyl] trial ({},{},{}) angular coverage failed: {} faces",
                                fi, n1, n2, trial.faces.len());
                        }
                        continue;
                    }
                }

                if viz_quit.get() {
                    return (hypotheses, true);
                }

                if let Some(trial) = trial {
                    if best_candidate.as_ref().is_none_or(|b| trial.total_area > b.total_area) {
                        best_candidate = Some(trial);
                    }
                }
            }
        }

        if let Some(candidate) = best_candidate {

            // Viz: post-BFS pause — show accepted hypothesis
            if !viz_quit.get() {
                if let Some(viz) = viz {
                    if verbosity >= 2 {
                        eprintln!("  [BFS-cyl] ACCEPTED: {} faces, r={:.6}, err_max={:.2e}",
                            candidate.faces.len(), candidate.radius, candidate.error_max);
                    }
                    let mut highlights = vec![
                        viz::FaceHighlight { face_indices: candidate.faces.clone(), color: [0.0, 0.8, 0.0, 1.0] },
                    ];
                    // Gray background for already-committed cylindrical faces
                    let bg_faces: Vec<usize> = (0..mesh.faces.len())
                        .filter(|f| mesh.faces[*f].cylindrical_hypothesis >= 0 && !candidate.faces.contains(f))
                        .collect();
                    if !bg_faces.is_empty() {
                        highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                    }
                    if let Some(action) = viz_custom(
                        Some(viz), highlights,
                        Vec::new(),
                        &format!("BFS-cyl result: ACCEPTED {} faces r={:.4} [space=next]",
                            candidate.faces.len(), candidate.radius),
                        vec![centered_cylinder_overlay(
                            candidate.axis_origin, candidate.axis_direction, candidate.radius,
                            &candidate.faces, mesh,
                            [0.0, 0.8, 0.0, 0.3],
                        )], Vec::new(),
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
                mesh.faces[f].cylindrical_hypothesis = hi;
            }

            hypotheses.push(CylindricalHypothesis {
                axis_origin: candidate.axis_origin,
                axis_direction: candidate.axis_direction,
                radius: candidate.radius,
                convex: candidate.convex,
                faces: candidate.faces,
                vertices: candidate.vertices.into_iter().collect(),
                error_max: candidate.error_max,
                centroid_error_max: candidate.centroid_error_max,
                error_abs_sum: candidate.error_abs_sum,
            });
        } else {
            // Viz: post-BFS pause — no accepted hypothesis for this seed face
            if !viz_quit.get() {
                if let Some(viz) = viz {
                    if verbosity >= 2 {
                        eprintln!("  [BFS-cyl] REJECTED: fi={} — no valid cylinder found", fi);
                    }
                    if let Some(action) = viz_custom(
                        Some(viz),
                        vec![viz::FaceHighlight { face_indices: vec![fi], color: [1.0, 0.0, 0.0, 1.0] }],
                        Vec::new(),
                        &format!("BFS-cyl result: REJECTED fi={fi} [space=next]"),
                        Vec::new(), Vec::new(),
                        Some(viz_face_centroid(fi, mesh)),
                        viz_face_normal(fi, mesh),
                    ) {
                        match action {
                            VizAction::Quit => { return (hypotheses, true); }
                            VizAction::NextSeed | VizAction::NextStep => {}
                        }
                    }
                }
            }
        }
        // If no valid triple found, leave fi as UNDEDUCED (may be absorbed by
        // a later face's BFS). Do NOT set NO_HYPOTHESIS yet.
    }

    // After all faces processed: set any remaining UNDEDUCED faces to NO_HYPOTHESIS
    for fi in 0..num_faces {
        if mesh.faces[fi].cylindrical_hypothesis == UNDEDUCED_CYLINDRICAL_HYPOTHESIS {
            mesh.faces[fi].cylindrical_hypothesis = NO_HYPOTHESIS;
        }
    }

    (hypotheses, false)
}

// ---------------------------------------------------------------------------
// Stage 2.2: Compare check
// ---------------------------------------------------------------------------

/// Validate fitted cylindrical hypotheses against a reference STEP shape.
pub(super) fn compare_cylindrical_hypotheses(
    hypotheses: &[CylindricalHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);

            // Project centroid onto cylinder surface (nearest point)
            let d = [
                centroid[0] - hyp.axis_origin[0],
                centroid[1] - hyp.axis_origin[1],
                centroid[2] - hyp.axis_origin[2],
            ];
            let t = dot3(&d, &hyp.axis_direction);
            let radial = [
                d[0] - t * hyp.axis_direction[0],
                d[1] - t * hyp.axis_direction[1],
                d[2] - t * hyp.axis_direction[2],
            ];
            let radial_dist = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();

            let projected = if radial_dist > 1e-15 {
                let scale = hyp.radius / radial_dist;
                [
                    hyp.axis_origin[0] + t * hyp.axis_direction[0] + radial[0] * scale,
                    hyp.axis_origin[1] + t * hyp.axis_direction[1] + radial[1] * scale,
                    hyp.axis_origin[2] + t * hyp.axis_direction[2] + radial[2] * scale,
                ]
            } else {
                // Centroid is on the axis — unlikely but handle gracefully
                centroid
            };

            let pt = gp::Pnt::new_real3(projected[0], projected[1], projected[2]);
            let dist = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(dist);
        }

        if max_dist > config.surface_tolerance_mm {
            return Err(Stage2CompareError {
                hypothesis_type: "cylindrical",
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance: config.surface_tolerance_mm,
            });
        }
    }

    Ok(())
}
