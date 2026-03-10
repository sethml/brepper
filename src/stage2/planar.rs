//! Stage 2.1: Planar hypothesis deduction and comparison.

use std::collections::{HashSet, VecDeque};

use opencascade_sys::gp;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_PLANAR_HYPOTHESIS};
use crate::viz::{self, VizAction, VizSender};

use super::{
    face_area, viz_bfs_seed, viz_bfs_step, viz_custom, viz_face_centroid,
    viz_face_normal, PlanarHypothesis, Stage2CompareError, REFIT_SKIP_MULTIPLIER,
};

// ---------------------------------------------------------------------------
// Planar fitting helpers
// ---------------------------------------------------------------------------

/// Signed distance from a vertex to a plane defined by (normal, distance).
pub(super) fn vertex_to_plane_distance(v: &MeshVertex, normal: &[f64; 3], distance: f64) -> f64 {
    normal[0] * v.x + normal[1] * v.y + normal[2] * v.z - distance
}

/// Fit a plane to a set of faces using area-weighted normal averaging and
/// vertex centroid for the distance.
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

    let mut sum_d = 0.0_f64;
    for &vi in vertex_set {
        let v = &vertices[vi];
        sum_d += normal[0] * v.x + normal[1] * v.y + normal[2] * v.z;
    }
    let distance = sum_d / vertex_set.len() as f64;

    (normal, distance)
}

/// Check that all vertices in a set are within tolerance of a plane.
fn all_vertices_within_tolerance(
    vertex_set: &HashSet<usize>,
    normal: &[f64; 3],
    distance: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        if vertex_to_plane_distance(&vertices[vi], normal, distance).abs() > tolerance {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// Stage 2.1: Planar hypothesis deduction
// ---------------------------------------------------------------------------

/// Deduce planar hypotheses from the mesh using BFS region growing.
///
/// For each unassigned face, creates a new planar hypothesis seeded from that
/// face's plane, then grows it via BFS to neighboring faces that are coplanar
/// (normal alignment and vertex-to-plane distance within tolerance).
pub(super) fn deduce_planar_hypotheses(
    mesh: &mut ConnectedMesh,
    vertex_tol: f64,
    verbosity: u8,
    viz: Option<&VizSender>,
) -> (Vec<PlanarHypothesis>, bool) {
    let num_faces = mesh.faces.len();
    let mut hypotheses: Vec<PlanarHypothesis> = Vec::new();
    let mut user_quit = false;

    for face in &mut mesh.faces {
        face.planar_hypothesis = UNDEDUCED_PLANAR_HYPOTHESIS;
    }

    for fi in 0..num_faces {
        if mesh.faces[fi].planar_hypothesis != UNDEDUCED_PLANAR_HYPOTHESIS {
            continue;
        }

        let hi = hypotheses.len() as i32;
        let seed_normal = mesh.faces[fi].normal.unwrap();

        // Compute initial plane from seed face: normal from face, distance
        // as average projection of seed face vertices.
        let vc = mesh.faces[fi].vertex_count as usize;
        let mut sum_d = 0.0_f64;
        let mut vertex_set = HashSet::new();
        for vi_idx in 0..vc {
            let vi = mesh.faces[fi].vertex_indices[vi_idx];
            vertex_set.insert(vi);
            let v = &mesh.vertices[vi];
            sum_d += seed_normal[0] * v.x + seed_normal[1] * v.y + seed_normal[2] * v.z;
        }
        let mut current_normal = seed_normal;
        let mut current_distance = sum_d / vc as f64;
        let mut face_list: Vec<usize> = vec![fi];

        mesh.faces[fi].planar_hypothesis = hi;

        // Level-3 trace: print seed face info
        if verbosity >= 3 {
            eprintln!(
                "[BFS-plane hi={}] Seed fi={}: normal=[{:.4},{:.4},{:.4}], d={:.4}",
                hi, fi,
                current_normal[0], current_normal[1], current_normal[2],
                current_distance,
            );
            for vi_idx in 0..vc {
                let vi = mesh.faces[fi].vertex_indices[vi_idx];
                let v = &mesh.vertices[vi];
                eprintln!("  seed vertex vi={}: [{:.4},{:.4},{:.4}]", vi, v.x, v.y, v.z);
            }
        }

        // Collect background faces for viz: faces assigned to multi-face planar hypotheses
        let bg_faces: Vec<usize> = if viz.is_some() {
            (0..num_faces).filter(|&f| {
                let ph = mesh.faces[f].planar_hypothesis;
                ph >= 0 && hypotheses.get(ph as usize).is_some_and(|h| h.faces.len() > 1)
            }).collect()
        } else { Vec::new() };

        // Viz: show seed face
        let mut skip_viz = false;
        if let Some(action) = viz_bfs_seed(
            viz,
            &[fi],
            &format!("BFS-plane hi={hi}: seed fi={fi} [space=step, shift+space=skip]"),
            Vec::new(), Vec::new(), &bg_faces, mesh,
        ) {
            match action {
                VizAction::Quit => { user_quit = true; return (hypotheses, user_quit); }
                VizAction::NextSeed => { skip_viz = true; }
                VizAction::NextStep => {}
            }
        }

        // BFS expansion
        let mut queue = VecDeque::new();
        queue.push_back(fi);

        while let Some(current_fi) = queue.pop_front() {
            let vertex_count = mesh.faces[current_fi].vertex_count as usize;
            let neighbors = mesh.faces[current_fi].neighbors;

            for &ni in &neighbors[..vertex_count] {
                if ni < 0 {
                    continue;
                }
                let ni = ni as usize;

                if mesh.faces[ni].planar_hypothesis != UNDEDUCED_PLANAR_HYPOTHESIS {
                    if verbosity >= 4 {
                        eprintln!("  [BFS-plane] fi={ni}: already assigned hyp={} → skip", mesh.faces[ni].planar_hypothesis);
                    }
                    continue;
                }
                // Vertex distance check
                let nvc = mesh.faces[ni].vertex_count as usize;
                let nvi = mesh.faces[ni].vertex_indices;
                let mut all_ok = true;
                let mut any_far = false;
                let mut vtx_err_max = 0.0_f64;
                for &vi in &nvi[..nvc] {
                    let d = vertex_to_plane_distance(
                        &mesh.vertices[vi],
                        &current_normal,
                        current_distance,
                    );
                    let abs_d = d.abs();
                    vtx_err_max = vtx_err_max.max(abs_d);
                    if abs_d > vertex_tol {
                        all_ok = false;
                        if abs_d > REFIT_SKIP_MULTIPLIER * vertex_tol {
                            any_far = true;
                            break;
                        }
                    }
                }

                // Skip if any vertex is too far for re-fitting to help
                if any_far {
                    if verbosity >= 3 {
                        eprintln!(
                            "  [BFS-plane] from fi={} try ni={}: vtx_err_max={:.2e} > REFIT_SKIP*tol={:.2e} → REJECT(too far)",
                            current_fi, ni, vtx_err_max, REFIT_SKIP_MULTIPLIER * vertex_tol,
                        );
                    }
                    continue;
                }

                if !all_ok {
                    // Try re-fitting with current + new face's vertices
                    let mut trial_vertices = vertex_set.clone();
                    for &vi in &nvi[..nvc] {
                        trial_vertices.insert(vi);
                    }
                    let mut trial_faces = face_list.clone();
                    trial_faces.push(ni);

                    let (new_normal, new_distance) = fit_plane(
                        &trial_faces,
                        &trial_vertices,
                        &mesh.faces,
                        &mesh.vertices,
                    );

                    // Check ALL vertices (existing + new) against the re-fitted plane
                    if !all_vertices_within_tolerance(
                        &trial_vertices,
                        &new_normal,
                        new_distance,
                        vertex_tol,
                        &mesh.vertices,
                    ) {
                        if verbosity >= 3 {
                            eprintln!(
                                "  [BFS-plane] from fi={} try ni={}: vtx_err_max={:.2e} (needs refit) → REJECT(refit failed)",
                                current_fi, ni, vtx_err_max,
                            );
                        }
                        continue;
                    }

                    if verbosity >= 3 {
                        eprintln!(
                            "  [BFS-plane] from fi={} try ni={}: vtx_err_max={:.2e} (refit ok, new_normal=[{:.4},{:.4},{:.4}] d={:.4}) → ACCEPT[refit]",
                            current_fi, ni, vtx_err_max,
                            new_normal[0], new_normal[1], new_normal[2], new_distance,
                        );
                    }

                    // Accept re-fit
                    current_normal = new_normal;
                    current_distance = new_distance;
                } else if verbosity >= 3 {
                    eprintln!(
                        "  [BFS-plane] from fi={} try ni={}: vtx_err_max={:.2e} ≤ tol={:.2e} → ACCEPT",
                        current_fi, ni, vtx_err_max, vertex_tol,
                    );
                }

                // Accept this face into the hypothesis
                mesh.faces[ni].planar_hypothesis = hi;
                face_list.push(ni);
                for &vi in &nvi[..nvc] {
                    vertex_set.insert(vi);
                }
                queue.push_back(ni);

                // Viz: show accepted face
                if !skip_viz {
                    if let Some(action) = viz_bfs_step(
                        viz, &[fi], &face_list, ni,
                        &format!("BFS-plane hi={hi}: accepted fi={ni} ({} faces) [space=step, shift+space=skip]", face_list.len()),
                        Vec::new(), Vec::new(), &bg_faces, mesh,
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

        // Final re-fit using all collected faces and vertices
        let (final_normal, final_distance) =
            fit_plane(&face_list, &vertex_set, &mesh.faces, &mesh.vertices);

        // Compute error metrics
        let mut error_max = f64::NEG_INFINITY;
        let mut error_min = f64::INFINITY;
        let mut error_abs_sum = 0.0_f64;
        for &vi in &vertex_set {
            let d = vertex_to_plane_distance(&mesh.vertices[vi], &final_normal, final_distance);
            error_max = error_max.max(d);
            error_min = error_min.min(d);
            error_abs_sum += d.abs();
        }

        if verbosity >= 2 {
            eprintln!(
                "[BFS-plane] hi={hi}: ACCEPTED ({} faces) err_max={:.2e} err_min={:.2e}",
                face_list.len(), error_max, error_min,
            );
        }

        // Viz: post-BFS pause showing accepted hypothesis in green
        if !skip_viz && face_list.len() > 1 {
            let non_seed: Vec<usize> = face_list.iter()
                .filter(|&&f| f != fi)
                .copied().collect();
            let mut highlights = vec![
                viz::FaceHighlight { face_indices: bg_faces.clone(), color: [0.5, 0.5, 0.5, 1.0] },
                viz::FaceHighlight { face_indices: vec![fi], color: [0.0, 0.8, 0.0, 1.0] },
            ];
            if !non_seed.is_empty() {
                highlights.push(viz::FaceHighlight { face_indices: non_seed, color: [0.0, 0.8, 0.0, 1.0] });
            }
            if let Some(action) = viz_custom(
                viz, highlights,
                Vec::new(),
                &format!("BFS-plane hi={hi}: ACCEPTED ({} faces) [space=next]", face_list.len()),
                Vec::new(), Vec::new(),
                Some(viz_face_centroid(fi, mesh)),
                viz_face_normal(fi, mesh),
            ) {
                if matches!(action, VizAction::Quit) { return (hypotheses, true); }
            }
        }

        hypotheses.push(PlanarHypothesis {
            normal: final_normal,
            distance: final_distance,
            faces: face_list,
            vertices: vertex_set.into_iter().collect(),
            error_max,
            error_min,
            error_abs_sum,
        });
    }

    (hypotheses, user_quit)
}

// ---------------------------------------------------------------------------
// Stage 2.1: Compare check
// ---------------------------------------------------------------------------

/// Validate fitted planar hypotheses against a reference STEP shape.
///
/// For each hypothesis, projects face centroids onto the fitted plane and
/// checks that those projected points are within surface_tolerance of the
/// reference STEP surface.
pub(super) fn compare_planar_hypotheses(
    hypotheses: &[PlanarHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let face = &mesh.faces[fi];
            let n = face.vertex_count as usize;
            let mut cx = 0.0_f64;
            let mut cy = 0.0_f64;
            let mut cz = 0.0_f64;
            for vi_idx in 0..n {
                let v = &mesh.vertices[face.vertex_indices[vi_idx]];
                cx += v.x;
                cy += v.y;
                cz += v.z;
            }
            let inv_n = 1.0 / n as f64;
            cx *= inv_n;
            cy *= inv_n;
            cz *= inv_n;

            // Project centroid onto the fitted plane
            let dist_to_plane = hyp.normal[0] * cx + hyp.normal[1] * cy + hyp.normal[2] * cz
                - hyp.distance;
            let px = cx - dist_to_plane * hyp.normal[0];
            let py = cy - dist_to_plane * hyp.normal[1];
            let pz = cz - dist_to_plane * hyp.normal[2];

            let pt = gp::Pnt::new_real3(px, py, pz);
            let d = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(d);
        }

        if max_dist > config.surface_tolerance_mm {
            return Err(Stage2CompareError {
                hypothesis_type: "planar",
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance: config.surface_tolerance_mm,
            });
        }
    }

    Ok(())
}
