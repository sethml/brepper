//! Stage 2.3: Spherical hypothesis deduction and comparison.

use std::collections::{HashSet, VecDeque};

use opencascade_sys::gp;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_SPHERICAL_HYPOTHESIS, NO_HYPOTHESIS};
use crate::viz::{self, VizAction, VizSender};

use super::{
    dot3, eigenvalues_3x3_symmetric, face_area, face_centroid, normalize3,
    viz_bfs_seed, viz_bfs_step, viz_custom, viz_face_centroid, viz_face_normal,
    PlanarHypothesis, SphericalHypothesis, Stage2CompareError,
    MIN_SPHERE_EIGENVALUE_RATIO, MIN_SPHERE_FACES,
    REFIT_SKIP_MULTIPLIER,
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

    Some(([cx, cy, cz], r_sq.sqrt()))
}

/// Check if all vertices are within tolerance of a sphere surface.
fn all_vertices_within_sphere_tolerance(
    vertex_set: &HashSet<usize>,
    center: &[f64; 3],
    radius: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        if vertex_to_sphere_distance(&vertices[vi], center, radius).abs() > tolerance {
            return false;
        }
    }
    true
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

/// Build a map from vertex index to the set of face indices incident on that vertex.
fn build_vertex_to_faces_map(mesh: &ConnectedMesh) -> Vec<Vec<usize>> {
    let mut vtf: Vec<Vec<usize>> = vec![Vec::new(); mesh.vertices.len()];
    for (fi, face) in mesh.faces.iter().enumerate() {
        let vc = face.vertex_count as usize;
        for &vi in &face.vertex_indices[..vc] {
            vtf[vi].push(fi);
        }
    }
    vtf
}

// ---------------------------------------------------------------------------
// Stage 2.3: Spherical hypothesis deduction
// ---------------------------------------------------------------------------

/// Deduce spherical hypotheses from the mesh using vertex-neighborhood seeding
/// and BFS region growing.
///
/// For each mesh vertex, collect the surrounding faces and fit a sphere to their
/// vertices. If the fit is good, seed a BFS to grow the hypothesis. After BFS,
/// validate solid-angle coverage to reject fillet-strip growth.
#[allow(clippy::too_many_arguments)]
pub(super) fn deduce_spherical_hypotheses(
    mesh: &mut ConnectedMesh,
    planar_hypotheses: &[PlanarHypothesis],
    vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
    max_sphere_radius: f64,
    verbosity: u8,
    viz: Option<&VizSender>,
) -> (Vec<SphericalHypothesis>, bool) {
    let mut hypotheses: Vec<SphericalHypothesis> = Vec::new();
    let mut user_quit = false;

    // Build vertex-to-faces map for vertex-neighborhood seeding
    let vtf = build_vertex_to_faces_map(mesh);
    for face_indices in &vtf {
        // Collect surrounding faces that are valid seed candidates
        let surrounding: Vec<usize> = face_indices
            .iter()
            .copied()
            .filter(|&fi| {
                // Must be undeduced
                if mesh.faces[fi].spherical_hypothesis != UNDEDUCED_SPHERICAL_HYPOTHESIS {
                    return false;
                }
                // Skip multi-face planar faces (genuinely flat)
                let ph = mesh.faces[fi].planar_hypothesis;
                if ph >= 0 && planar_hypotheses[ph as usize].faces.len() > 1 {
                    return false;
                }
                // Must have valid normal
                mesh.faces[fi].normal.is_some()
            })
            .collect();

        if surrounding.len() < 3 {
            continue;
        }

        // Collect all vertices from surrounding faces
        let mut seed_vertex_set: HashSet<usize> = HashSet::new();
        for &fi in &surrounding {
            let vc = mesh.faces[fi].vertex_count as usize;
            for &v in &mesh.faces[fi].vertex_indices[..vc] {
                seed_vertex_set.insert(v);
            }
        }

        // Fit sphere to seed vertices
        let (center, radius) = match fit_sphere(&seed_vertex_set, &mesh.vertices) {
            Some(params) => params,
            None => continue,
        };

        if radius > max_sphere_radius {
            continue;
        }

        // Verify all seed vertices within tolerance
        if !all_vertices_within_sphere_tolerance(
            &seed_vertex_set, &center, radius, vertex_tol, &mesh.vertices,
        ) {
            continue;
        }

        // Angular tolerance: reject if any pair of adjacent seed faces exceeds limit
        let mut angular_reject = false;
        for &fi in &surrounding {
            let fi_n = match mesh.faces[fi].normal {
                Some(n) => n,
                None => continue,
            };
            let fi_vc = mesh.faces[fi].vertex_count as usize;
            for &ni in &mesh.faces[fi].neighbors[..fi_vc] {
                if ni < 0 { continue; }
                let ni = ni as usize;
                if !surrounding.contains(&ni) { continue; }
                if let Some(ni_n) = mesh.faces[ni].normal {
                    let cos_a = dot3(&fi_n, &ni_n).clamp(-1.0, 1.0);
                    if cos_a.acos() > angular_tol {
                        angular_reject = true;
                        break;
                    }
                }
            }
            if angular_reject { break; }
        }
        if angular_reject {
            continue;
        }

        // Centroid validation for seed faces
        let mut centroid_reject = false;
        for &fi in &surrounding {
            let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
            let d = vertex_to_sphere_distance(
                &MeshVertex::from_xyz(c[0], c[1], c[2]),
                &center, radius,
            ).abs();
            if d > surface_tol {
                centroid_reject = true;
                break;
            }
        }
        if centroid_reject {
            continue;
        }

        // Determine convexity from first seed face
        let convex = determine_sphere_convexity(
            &mesh.faces[surrounding[0]], &mesh.vertices, &center,
        );
        // Check convexity consistency across seed faces
        let mut convex_reject = false;
        for &fi in &surrounding[1..] {
            if determine_sphere_convexity(&mesh.faces[fi], &mesh.vertices, &center) != convex {
                convex_reject = true;
                break;
            }
        }
        if convex_reject {
            continue;
        }

        // Seed is valid — assign seed faces and start BFS
        let hi = hypotheses.len() as i32;
        for &fi in &surrounding {
            mesh.faces[fi].spherical_hypothesis = hi;
        }
        let mut face_list: Vec<usize> = surrounding.clone();
        let mut vertex_set: HashSet<usize> = seed_vertex_set;
        let mut current_center = center;
        let mut current_radius = radius;

        let mut queue: VecDeque<usize> = surrounding.iter().copied().collect();

        // Level-3 trace: print seed info
        if verbosity >= 3 {
            eprintln!(
                "[BFS-sph hi={}] Seed: {} faces, center=[{:.4},{:.4},{:.4}], r={:.6}, {}",
                hi, surrounding.len(),
                center[0], center[1], center[2], radius,
                if convex { "convex" } else { "concave" },
            );
            for &fi in &surrounding {
                let vc = mesh.faces[fi].vertex_count as usize;
                let vis: Vec<usize> = mesh.faces[fi].vertex_indices[..vc].to_vec();
                let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);
                let cen_err = vertex_to_sphere_distance(
                    &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                    &center, radius,
                ).abs();
                let coords: Vec<String> = vis.iter().map(|&vi| {
                    let v = &mesh.vertices[vi];
                    format!("vi={}:[{:.4},{:.4},{:.4}]", vi, v.x, v.y, v.z)
                }).collect();
                eprintln!(
                    "  seed fi={}: {}v centroid=[{:.4},{:.4},{:.4}] cen_err={:.2e}  vertices: {}",
                    fi, vc, centroid[0], centroid[1], centroid[2], cen_err,
                    coords.join(" "),
                );
            }
        }

        // Viz: show seed faces
        let mut skip_viz = false;
        let bg_faces: Vec<usize> = (0..mesh.faces.len())
            .filter(|f| mesh.faces[*f].spherical_hypothesis >= 0)
            .collect();
        if let Some(action) = viz_bfs_seed(
            viz, &surrounding,
            &format!("BFS-sph hi={hi}: {} seed faces, r={:.4} {} [space=step, shift+space=skip]",
                surrounding.len(), current_radius, if convex { "convex" } else { "concave" }),
            Vec::new(),
            vec![viz::SphereOverlay {
                center: current_center,
                radius: current_radius,
                color: [0.2, 0.4, 1.0, 0.3],
            }],
            &bg_faces, mesh,
        ) {
            match action {
                VizAction::Quit => { user_quit = true; return (hypotheses, user_quit); }
                VizAction::NextSeed => { skip_viz = true; }
                VizAction::NextStep => {}
            }
        }


        // BFS expansion
        while let Some(current_fi) = queue.pop_front() {
            let vc = mesh.faces[current_fi].vertex_count as usize;
            let neighbors = mesh.faces[current_fi].neighbors;

            for &cni in &neighbors[..vc] {
                if cni < 0 { continue; }
                let cni = cni as usize;

                if mesh.faces[cni].spherical_hypothesis != UNDEDUCED_SPHERICAL_HYPOTHESIS {
                    if verbosity >= 4 {
                        eprintln!("  [BFS-sph] from fi={} try cni={}: already assigned (hyp={}) → SKIP", current_fi, cni, mesh.faces[cni].spherical_hypothesis);
                    }
                    continue;
                }

                // Convexity check
                let cni_convex = determine_sphere_convexity(
                    &mesh.faces[cni], &mesh.vertices, &current_center,
                );
                if cni_convex != convex {
                    if verbosity >= 3 {
                        eprintln!(
                            "  [BFS-sph] from fi={} try cni={}: convexity mismatch ({} vs {}) → REJECT(convexity)",
                            current_fi, cni,
                            if cni_convex { "convex" } else { "concave" },
                            if convex { "convex" } else { "concave" },
                        );
                    }
                    continue;
                }

                // Angular tolerance: reject if dihedral angle with ANY already-assigned
                // neighbor exceeds the limit (defense-in-depth against creased surfaces)
                if let Some(n_cni) = mesh.faces[cni].normal {
                    let cni_vc2 = mesh.faces[cni].vertex_count as usize;
                    let cni_neighbors = mesh.faces[cni].neighbors;
                    let mut ang_reject = false;
                    let mut worst_angle_deg = 0.0_f64;
                    let mut worst_adj = usize::MAX;
                    for &adj in &cni_neighbors[..cni_vc2] {
                        if adj < 0 { continue; }
                        let adj = adj as usize;
                        if mesh.faces[adj].spherical_hypothesis != hi { continue; }
                        if let Some(n_adj) = mesh.faces[adj].normal {
                            let cos_a = dot3(&n_adj, &n_cni).clamp(-1.0, 1.0);
                            let angle = cos_a.acos();
                            if angle > worst_angle_deg.to_radians() {
                                worst_angle_deg = angle.to_degrees();
                                worst_adj = adj;
                            }
                            if angle > angular_tol {
                                ang_reject = true;
                                break;
                            }
                        }
                    }
                    if ang_reject {
                        if verbosity >= 3 {
                            eprintln!(
                                "  [BFS-sph] from fi={} try cni={}: angular {:.2}° > tol {:.2}° (adj fi={}) → REJECT(angular)",
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
                for &v in &cni_vi[..cni_vc] {
                    let d = vertex_to_sphere_distance(
                        &mesh.vertices[v],
                        &current_center,
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
                            "  [BFS-sph] from fi={} try cni={}: vtx_err_max={:.2e} > REFIT_SKIP*tol={:.2e} → REJECT(too far)",
                            current_fi, cni, vtx_err_max, REFIT_SKIP_MULTIPLIER * vertex_tol,
                        );
                    }
                    continue;
                }

                // Centroid validation (surface_tol check during BFS)
                let centroid = face_centroid(&mesh.faces[cni], &mesh.vertices);
                let centroid_dist = vertex_to_sphere_distance(
                    &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                    &current_center, current_radius,
                ).abs();
                if centroid_dist > REFIT_SKIP_MULTIPLIER * surface_tol {
                    if verbosity >= 3 {
                        eprintln!(
                            "  [BFS-sph] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} > REFIT_SKIP*stol={:.2e} → REJECT(centroid far)",
                            current_fi, cni, vtx_err_max, centroid_dist, REFIT_SKIP_MULTIPLIER * surface_tol,
                        );
                    }
                    continue;
                }
                let needs_refit = !all_ok || centroid_dist > surface_tol;

                if needs_refit {
                    // Try re-fitting with new face included
                    let mut trial_vertices = vertex_set.clone();
                    for &v in &cni_vi[..cni_vc] {
                        trial_vertices.insert(v);
                    }

                    let (new_center, new_radius) = match fit_sphere(
                        &trial_vertices, &mesh.vertices,
                    ) {
                        Some(params) => params,
                        None => {
                            if verbosity >= 3 {
                                eprintln!(
                                    "  [BFS-sph] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit degenerate) → REJECT(refit failed)",
                                    current_fi, cni, vtx_err_max, centroid_dist,
                                );
                            }
                            continue;
                        }
                    };

                    if new_radius > max_sphere_radius {
                        if verbosity >= 3 {
                            eprintln!(
                                "  [BFS-sph] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit r={:.4} > max) → REJECT(radius)",
                                current_fi, cni, vtx_err_max, centroid_dist, new_radius,
                            );
                        }
                        continue;
                    }

                    if !all_vertices_within_sphere_tolerance(
                        &trial_vertices, &new_center, new_radius,
                        vertex_tol, &mesh.vertices,
                    ) {
                        if verbosity >= 3 {
                            eprintln!(
                                "  [BFS-sph] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit vertex check failed) → REJECT(refit tol)",
                                current_fi, cni, vtx_err_max, centroid_dist,
                            );
                        }
                        continue;
                    }

                    // Verify all existing centroids still pass
                    let mut refit_centroid_ok = true;
                    for &f in &face_list {
                        let c = face_centroid(&mesh.faces[f], &mesh.vertices);
                        let d = vertex_to_sphere_distance(
                            &MeshVertex::from_xyz(c[0], c[1], c[2]),
                            &new_center, new_radius,
                        ).abs();
                        if d > surface_tol {
                            refit_centroid_ok = false;
                            break;
                        }
                    }
                    if !refit_centroid_ok {
                        if verbosity >= 3 {
                            eprintln!(
                                "  [BFS-sph] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit existing centroid failed) → REJECT(refit centroid)",
                                current_fi, cni, vtx_err_max, centroid_dist,
                            );
                        }
                        continue;
                    }

                    if verbosity >= 3 {
                        eprintln!(
                            "  [BFS-sph] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} (refit ok → r={:.6} center=[{:.4},{:.4},{:.4}]) → ACCEPT[refit]",
                            current_fi, cni, vtx_err_max, centroid_dist,
                            new_radius, new_center[0], new_center[1], new_center[2],
                        );
                    }

                    // Accept re-fit
                    current_center = new_center;
                    current_radius = new_radius;
                } else if verbosity >= 3 {
                    eprintln!(
                        "  [BFS-sph] from fi={} try cni={}: vtx_err_max={:.2e} cen_err={:.2e} → ACCEPT",
                        current_fi, cni, vtx_err_max, centroid_dist,
                    );
                }

                // Accept this face
                mesh.faces[cni].spherical_hypothesis = hi;
                face_list.push(cni);
                for &v in &cni_vi[..cni_vc] {
                    vertex_set.insert(v);
                }
                queue.push_back(cni);

                // Viz: show accepted face
                if !skip_viz {
                    if let Some(action) = viz_bfs_step(
                        viz, &surrounding, &face_list, cni,
                        &format!("BFS-sph hi={hi}: accepted fi={cni} ({} faces) r={:.4} [space=step, shift+space=skip]", face_list.len(), current_radius),
                        Vec::new(),
                        vec![viz::SphereOverlay {
                            center: current_center,
                            radius: current_radius,
                            color: [0.2, 0.4, 1.0, 0.3],
                        }],
                        &bg_faces, mesh,
                    ) {
                        match action {
                            VizAction::Quit => { user_quit = true; return (hypotheses, user_quit); }
                            VizAction::NextSeed => { skip_viz = true; }
                            VizAction::NextStep => {}
                        }
                    }
                }
            }
        }

        // Final re-fit from all accumulated faces
        if let Some((final_center, final_radius)) =
            fit_sphere(&vertex_set, &mesh.vertices)
        {
            current_center = final_center;
            current_radius = final_radius;
        }

        // Compute error metrics
        let mut error_max = 0.0_f64;
        let mut error_abs_sum = 0.0_f64;
        for &v in &vertex_set {
            let d = vertex_to_sphere_distance(
                &mesh.vertices[v], &current_center, current_radius,
            ).abs();
            error_max = error_max.max(d);
            error_abs_sum += d;
        }

        // Validate: minimum face count, centroid check, radius, and solid-angle coverage
        let min_faces_ok = face_list.len() >= MIN_SPHERE_FACES;
        let radius_ok = current_radius <= max_sphere_radius;
        let coverage_ok = min_faces_ok && solid_angle_coverage_valid(
            &face_list, &current_center, &mesh.faces, &mesh.vertices,
        );

        if !min_faces_ok || !radius_ok || !coverage_ok {
            // Viz: post-BFS pause — rejected
            if !user_quit {
                if let Some(viz) = viz {
                    if verbosity >= 2 {
                        eprintln!("  [BFS-sph] REJECTED hi={hi}: {} faces, r={:.6} (min_faces={} radius={} coverage={})",
                            face_list.len(), current_radius, min_faces_ok, radius_ok, coverage_ok);
                    }
                    let mut highlights = vec![
                        viz::FaceHighlight { face_indices: surrounding.clone(), color: [1.0, 0.0, 0.0, 1.0] },
                    ];
                    if !bg_faces.is_empty() {
                        highlights.push(viz::FaceHighlight { face_indices: bg_faces.clone(), color: [0.5, 0.5, 0.5, 1.0] });
                    }
                    if let Some(action) = viz_custom(
                        Some(viz), highlights,
                        Vec::new(),
                        &format!("BFS-sph result: REJECTED hi={hi} {} faces r={:.4} [space=next]",
                            face_list.len(), current_radius),
                        Vec::new(),
                        vec![viz::SphereOverlay {
                            center: current_center, radius: current_radius,
                            color: [1.0, 0.0, 0.0, 0.3],
                        }],
                        Some(viz_face_centroid(surrounding[0], mesh)),
                        viz_face_normal(surrounding[0], mesh),
                    ) {
                        match action {
                            VizAction::Quit => { user_quit = true; return (hypotheses, user_quit); }
                            VizAction::NextSeed | VizAction::NextStep => {}
                        }
                    }
                }
            }
            // Undo assignments
            for &f in &face_list {
                mesh.faces[f].spherical_hypothesis = UNDEDUCED_SPHERICAL_HYPOTHESIS;
            }
            continue;
        }

        // Compute centroid error metrics (after validation to avoid extra work on rejected hypotheses)
        let mut centroid_error_max = 0.0_f64;
        for &f in &face_list {
            let c = face_centroid(&mesh.faces[f], &mesh.vertices);
            let d = vertex_to_sphere_distance(
                &MeshVertex::from_xyz(c[0], c[1], c[2]), &current_center, current_radius,
            ).abs();
            centroid_error_max = centroid_error_max.max(d);
        }

        // Viz: post-BFS pause — accepted
        if !user_quit {
            if let Some(viz) = viz {
                if verbosity >= 2 {
                    eprintln!("  [BFS-sph] ACCEPTED hi={hi}: {} faces, r={:.6}, err_max={:.2e}",
                        face_list.len(), current_radius, error_max);
                }
                let mut highlights = vec![
                    viz::FaceHighlight { face_indices: face_list.clone(), color: [0.0, 0.8, 0.0, 1.0] },
                ];
                if !bg_faces.is_empty() {
                    highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                }
                if let Some(action) = viz_custom(
                    Some(viz), highlights,
                    Vec::new(),
                    &format!("BFS-sph result: ACCEPTED hi={hi} {} faces r={:.4} [space=next]",
                        face_list.len(), current_radius),
                    Vec::new(),
                    vec![viz::SphereOverlay {
                        center: current_center, radius: current_radius,
                        color: [0.0, 0.8, 0.0, 0.3],
                    }],
                    Some(viz_face_centroid(face_list[0], mesh)),
                    viz_face_normal(face_list[0], mesh),
                ) {
                    match action {
                        VizAction::Quit => { user_quit = true; return (hypotheses, user_quit); }
                        VizAction::NextSeed | VizAction::NextStep => {}
                    }
                }
            }
        }

        hypotheses.push(SphericalHypothesis {
            center: current_center,
            radius: current_radius,
            convex,
            faces: face_list,
            vertices: vertex_set.into_iter().collect(),
            error_max,
            centroid_error_max,
            error_abs_sum,
        });
    }

    // Mark remaining undeduced faces as NO_HYPOTHESIS
    for face in &mut mesh.faces {
        if face.spherical_hypothesis == UNDEDUCED_SPHERICAL_HYPOTHESIS {
            face.spherical_hypothesis = NO_HYPOTHESIS;
        }
    }

    (hypotheses, user_quit)
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
