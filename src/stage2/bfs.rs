//! Generic BFS hypothesis generation framework.
//!
//! Provides a trait-based abstraction for surface hypothesis BFS, parameterized
//! by surface type. Each surface type (planar, conical, spherical, toroidal)
//! implements the `SurfaceType` trait, and the shared `generate_hypotheses()`
//! function handles the BFS expansion loop.
//!
//! Ported from the C++ `generate_hypotheses<Traits>()` template in brepper-cpp2.

use std::collections::{HashSet, VecDeque};

use crate::stage1::{ConnectedMesh, MeshFace, NO_HYPOTHESIS};
use crate::viz::VizSender;

use super::{dot3, face_area, face_centroid, REFIT_SKIP_MULTIPLIER};

// ---------------------------------------------------------------------------
// BFS configuration
// ---------------------------------------------------------------------------

/// Configuration for the BFS hypothesis generation.
pub(super) struct BfsConfig {
    pub vertex_tol: f64,
    pub surface_tol: f64,
    pub angular_tol: f64,
    pub verbosity: u8,
}

// ---------------------------------------------------------------------------
// Surface type trait
// ---------------------------------------------------------------------------

/// A seed for BFS expansion: a set of face indices to start from.
pub(super) struct Seed {
    pub faces: Vec<usize>,
}

/// Trait that each surface type implements to parameterize the BFS.
///
/// The generic `generate_hypotheses()` function calls these methods to:
/// - Generate seed face sets for each starting face
/// - Fit a surface to a set of faces/vertices
/// - Compute vertex and centroid distances to the fitted surface
/// - Check convexity consistency
/// - Validate angular coverage and other acceptance criteria
/// - Construct the final hypothesis struct
pub(super) trait SurfaceType {
    /// The fitted surface parameters (e.g., axis, radius, center).
    type Fit: Clone;
    /// The output hypothesis struct.
    type Hypothesis;

    /// Human-readable name for logging (e.g., "planar", "conical", "spherical").
    fn name() -> &'static str;

    /// Access the hypothesis index field on a MeshFace (read).
    fn hypothesis_field(face: &MeshFace) -> i32;

    /// Access the hypothesis index field on a MeshFace (write).
    fn set_hypothesis_field(face: &mut MeshFace, value: i32);

    /// The "undeduced" sentinel value for this surface type's hypothesis field.
    fn undeduced_value() -> i32;

    /// Check whether a face is eligible for BFS expansion (not already claimed
    /// by a prior hypothesis or excluded by other surface types).
    fn face_eligible(face: &MeshFace) -> bool {
        Self::hypothesis_field(face) == Self::undeduced_value()
    }

    /// Generate seed face sets anchored at face `fi`.
    /// Returns an empty vec if `fi` is not a valid seed anchor.
    fn generate_seeds(mesh: &ConnectedMesh, fi: usize, config: &BfsConfig) -> Vec<Seed>;

    /// Fit a surface to the given faces and vertices.
    /// Returns None if fitting fails.
    fn fit_surface(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        vertex_list: &[usize],
        verbosity: u8,
    ) -> Option<Self::Fit>;

    /// Compute the absolute distance from a vertex position to the fitted surface.
    fn vertex_distance(fit: &Self::Fit, pos: &[f64; 3]) -> f64;

    /// Compute the absolute distance from a centroid to the fitted surface.
    /// Default: same as vertex_distance.
    fn centroid_distance(fit: &Self::Fit, centroid: &[f64; 3]) -> f64 {
        Self::vertex_distance(fit, centroid)
    }

    /// Check convexity consistency: does `fi` match the convexity of the fit?
    /// Return true to accept, false to reject.
    /// Default implementation always accepts (e.g., for planar surfaces).
    fn convexity_check(_fit: &Self::Fit, _mesh: &ConnectedMesh, _fi: usize) -> bool {
        true
    }

    /// Determine convexity from a set of faces and update the fit.
    /// Default: no-op (for surface types without convexity).
    fn determine_convexity(_mesh: &ConnectedMesh, _face_list: &[usize], _fit: &mut Self::Fit) {}

    /// Post-BFS acceptance test (e.g., angular coverage, minimum face count).
    /// Return true to accept the hypothesis.
    fn acceptance_test(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        fit: &Self::Fit,
        config: &BfsConfig,
    ) -> bool;


    /// Construct the final hypothesis from the BFS result.
    fn make_hypothesis(
        mesh: &ConnectedMesh,
        faces: Vec<usize>,
        vertices: Vec<usize>,
        fit: &Self::Fit,
    ) -> Self::Hypothesis;
}

// ---------------------------------------------------------------------------
// Dihedral angle computation
// ---------------------------------------------------------------------------

/// Compute the dihedral angle in radians between two adjacent faces.
fn dihedral_angle(mesh: &ConnectedMesh, fi: usize, fj: usize) -> f64 {
    let n1 = match mesh.faces[fi].normal {
        Some(n) => n,
        None => return 0.0,
    };
    let n2 = match mesh.faces[fj].normal {
        Some(n) => n,
        None => return 0.0,
    };
    let d = dot3(&n1, &n2).clamp(-1.0, 1.0);
    d.acos()
}

// ---------------------------------------------------------------------------
// Generic BFS hypothesis generation
// ---------------------------------------------------------------------------

/// Run BFS hypothesis generation for a given surface type.
///
/// Iterates through all faces, generating seeds for each eligible face,
/// running trial BFS expansions, and committing the best hypothesis per
/// starting face. Returns the list of accepted hypotheses and whether
/// viz was quit.
pub(super) fn generate_hypotheses<S>(
    mesh: &mut ConnectedMesh,
    config: &BfsConfig,
    _viz: Option<&VizSender>,
) -> (Vec<S::Hypothesis>, bool)
where
    S: SurfaceType,
    S::Hypothesis: HypothesisFaces,
{
    let num_faces = mesh.faces.len();
    let vtol = config.vertex_tol;
    let stol = config.surface_tol;
    let skip_vtol = vtol * REFIT_SKIP_MULTIPLIER;
    let skip_stol = stol * REFIT_SKIP_MULTIPLIER;
    let angular_tol = config.angular_tol;

    let mut hypotheses: Vec<S::Hypothesis> = Vec::new();

    // Statistics
    let mut seeds_generated = 0_u32;
    let mut seeds_tried = 0_u32;
    let mut fit_failures = 0_u32;
    let mut seed_stol_failures = 0_u32;
    let mut seed_angle_failures = 0_u32;
    let mut acceptance_failures = 0_u32;
    let mut committed = 0_u32;

    // Initialize hypothesis fields to UNDEDUCED
    for fi in 0..num_faces {
        S::set_hypothesis_field(&mut mesh.faces[fi], S::undeduced_value());
    }

    for fi in 0..num_faces {
        if !S::face_eligible(&mesh.faces[fi]) {
            continue;
        }

        let seeds = S::generate_seeds(mesh, fi, config);
        if seeds.is_empty() {
            continue;
        }
        seeds_generated += seeds.len() as u32;

        let mut best_area = -1.0_f64;
        let mut best_hypothesis: Option<S::Hypothesis> = None;

        for seed in &seeds {
            let mut seed_vertex_set: HashSet<usize> = HashSet::new();
            for &sfi in &seed.faces {
                let vc = mesh.faces[sfi].vertex_count as usize;
                for &vi in &mesh.faces[sfi].vertex_indices[..vc] {
                    seed_vertex_set.insert(vi);
                }
            }
            let seed_vertices: Vec<usize> = seed_vertex_set.iter().copied().collect();

            let fit = match S::fit_surface(mesh, &seed.faces, &seed_vertices, config.verbosity) {
                Some(f) => f,
                None => {
                    fit_failures += 1;
                    continue;
                }
            };
            seeds_tried += 1;

            // Validate seed: vertex tolerance
            let mut seed_valid = true;
            let mut seed_max_vtol = 0.0_f64;
            for &vi in &seed_vertices {
                let v = &mesh.vertices[vi];
                let d = S::vertex_distance(&fit, &[v.x, v.y, v.z]);
                if d > seed_max_vtol {
                    seed_max_vtol = d;
                }
                if d > vtol {
                    seed_valid = false;
                    break;
                }
            }
            let bootstrapping = !seed_valid;
            if bootstrapping {
                if config.verbosity >= 3 {
                    eprintln!(
                        "    [{}-BFS] seed vtol fail (bootstrapping): max_vtx={:.6e} (limit {:.6e})",
                        S::name(), seed_max_vtol, vtol
                    );
                }
                seed_valid = true;
            }

            if !bootstrapping {
                let mut seed_max_stol = 0.0_f64;
                for &sfi in &seed.faces {
                    let c = face_centroid(&mesh.faces[sfi], &mesh.vertices);
                    let d = S::centroid_distance(&fit, &c);
                    if d > seed_max_stol {
                        seed_max_stol = d;
                    }
                    if d > stol {
                        seed_valid = false;
                        break;
                    }
                }
                if !seed_valid {
                    if config.verbosity >= 3 {
                        eprintln!(
                            "    [{}-BFS] seed stol fail: max_ctr={:.6e} (limit {:.6e})",
                            S::name(), seed_max_stol, stol
                        );
                    }
                    seed_stol_failures += 1;
                    continue;
                }
            }

            let mut angle_ok = true;
            'angle_check: for i in 0..seed.faces.len() {
                for j in (i + 1)..seed.faces.len() {
                    if is_neighbor(mesh, seed.faces[i], seed.faces[j])
                        && dihedral_angle(mesh, seed.faces[i], seed.faces[j]) > angular_tol
                    {
                        angle_ok = false;
                        break 'angle_check;
                    }
                }
            }
            if !angle_ok {
                seed_angle_failures += 1;
                continue;
            }

            let mut current_fit = fit;
            if !bootstrapping {
                S::determine_convexity(mesh, &seed.faces, &mut current_fit);
            }

            // --- BFS expansion ---
            let mut trial_faces: Vec<usize> = seed.faces.clone();
            let mut trial_face_set: HashSet<usize> = seed.faces.iter().copied().collect();
            let mut trial_vertex_set: HashSet<usize> = seed_vertex_set.clone();
            let mut trial_vertices: Vec<usize> = seed_vertices.clone();

            let mut bfs_queue = VecDeque::new();
            for &sfi in &seed.faces {
                bfs_queue.push_back(sfi);
            }

            while let Some(popped) = bfs_queue.pop_front() {
                let vc = mesh.faces[popped].vertex_count as usize;
                let neighbors = mesh.faces[popped].neighbors;

                for &ni in &neighbors[..vc] {
                    if ni < 0 {
                        continue;
                    }
                    let ni = ni as usize;
                    if !S::face_eligible(&mesh.faces[ni]) {
                        continue;
                    }
                    if trial_face_set.contains(&ni) {
                        continue;
                    }

                    let ni_vc = mesh.faces[ni].vertex_count as usize;
                    let ni_neighbors = mesh.faces[ni].neighbors;
                    let mut angle_reject = false;
                    for &nni in &ni_neighbors[..ni_vc] {
                        if nni >= 0 && trial_face_set.contains(&(nni as usize))
                            && dihedral_angle(mesh, ni, nni as usize) > angular_tol
                        {
                            angle_reject = true;
                            break;
                        }
                    }
                    if angle_reject {
                        continue;
                    }

                    if !bootstrapping && !S::convexity_check(&current_fit, mesh, ni) {
                        continue;
                    }

                    let ni_vi = mesh.faces[ni].vertex_indices;
                    let mut new_verts: Vec<usize> = Vec::new();
                    for &vi in &ni_vi[..ni_vc] {
                        if !trial_vertex_set.contains(&vi) {
                            new_verts.push(vi);
                        }
                    }

                    let accepted;

                    if bootstrapping {
                        accepted = true;
                    } else {
                        let mut max_vtx_dist = 0.0_f64;
                        for &vi in &ni_vi[..ni_vc] {
                            let v = &mesh.vertices[vi];
                            let d = S::vertex_distance(&current_fit, &[v.x, v.y, v.z]);
                            if d > max_vtx_dist {
                                max_vtx_dist = d;
                            }
                        }
                        let c = face_centroid(&mesh.faces[ni], &mesh.vertices);
                        let ctr_dist = S::centroid_distance(&current_fit, &c);

                        if max_vtx_dist <= vtol && ctr_dist <= stol {
                            accepted = true;
                        } else if max_vtx_dist > skip_vtol || ctr_dist > skip_stol {
                            continue;
                        } else {
                            let mut test_faces = trial_faces.clone();
                            test_faces.push(ni);
                            let mut test_verts = trial_vertices.clone();
                            for &vi in &new_verts {
                                test_verts.push(vi);
                            }

                            let refit = match S::fit_surface(mesh, &test_faces, &test_verts, 0) {
                                Some(f) => f,
                                None => { continue; },
                            };

                            let mut refit_ok = true;
                            for &tfi in &test_faces {
                                let tf = &mesh.faces[tfi];
                                let tvc = tf.vertex_count as usize;
                                for &vi in &tf.vertex_indices[..tvc] {
                                    let v = &mesh.vertices[vi];
                                    let d = S::vertex_distance(&refit, &[v.x, v.y, v.z]);
                                    if d > vtol {
                                        refit_ok = false;
                                        break;
                                    }
                                }
                                if !refit_ok {
                                    break;
                                }
                                let tc = face_centroid(&mesh.faces[tfi], &mesh.vertices);
                                let td = S::centroid_distance(&refit, &tc);
                                if td > stol {
                                    refit_ok = false;
                                    break;
                                }
                            }
                            if !refit_ok {
                                continue;
                            }

                            current_fit = refit;
                            accepted = true;
                        }
                    }

                    if accepted {
                        trial_faces.push(ni);
                        trial_face_set.insert(ni);
                        for &vi in &new_verts {
                            trial_vertex_set.insert(vi);
                            trial_vertices.push(vi);
                        }
                        bfs_queue.push_back(ni);
                    }
                }
            }

            // Final re-fit
            if let Some(final_fit) = S::fit_surface(
                mesh,
                &trial_faces,
                &trial_vertices,
                if bootstrapping { config.verbosity } else { 0 },
            ) {
                current_fit = final_fit;
            }

            // For bootstrapped hypotheses: validate and determine convexity
            if bootstrapping {
                let mut valid = true;
                let mut bootstrap_max_vtol = 0.0_f64;
                for &tfi in &trial_faces {
                    let tf = &mesh.faces[tfi];
                    let tvc = tf.vertex_count as usize;
                    for &vi in &tf.vertex_indices[..tvc] {
                        let v = &mesh.vertices[vi];
                        let d = S::vertex_distance(&current_fit, &[v.x, v.y, v.z]);
                        if d > bootstrap_max_vtol {
                            bootstrap_max_vtol = d;
                        }
                        if d > vtol {
                            valid = false;
                            break;
                        }
                    }
                    if !valid {
                        break;
                    }
                    let tc = face_centroid(&mesh.faces[tfi], &mesh.vertices);
                    if S::centroid_distance(&current_fit, &tc) > stol {
                        valid = false;
                        break;
                    }
                }
                if !valid {
                    if config.verbosity >= 3 {
                        eprintln!(
                            "    [{}-BFS] bootstrapped vtol fail: max={:.6e} (limit {:.6e}), {} faces",
                            S::name(), bootstrap_max_vtol, vtol, trial_faces.len()
                        );
                    }
                    acceptance_failures += 1;
                    continue;
                }
                // Determine convexity from the final fit
                S::determine_convexity(mesh, &trial_faces, &mut current_fit);
                if config.verbosity >= 2 {
                    eprintln!(
                        "    [{}-BFS] bootstrapped validated: {} faces, max_vtol={:.6e}",
                        S::name(), trial_faces.len(), bootstrap_max_vtol
                    );
                }
            }

            if !S::acceptance_test(mesh, &trial_faces, &current_fit, config) {
                acceptance_failures += 1;
                continue;
            }

            let total_area: f64 = trial_faces
                .iter()
                .map(|&tfi| face_area(&mesh.faces[tfi], &mesh.vertices))
                .sum();

            if total_area > best_area {
                best_area = total_area;
                let trial_verts: Vec<usize> = trial_vertex_set.into_iter().collect();
                best_hypothesis =
                    Some(S::make_hypothesis(mesh, trial_faces, trial_verts, &current_fit));
            }
        }

        // Commit best hypothesis
        if let Some(hyp) = best_hypothesis {
            committed += 1;
            let hyp_idx = hypotheses.len() as i32;
            for &tfi in hyp.faces() {
                S::set_hypothesis_field(&mut mesh.faces[tfi], hyp_idx);
            }
            if config.verbosity >= 1 {
                eprintln!(
                    "  [{}-BFS] committed hyp {} from fi={}: {} faces",
                    S::name(), hyp_idx, fi, hyp.faces().len()
                );
            }
            hypotheses.push(hyp);
        }
    }

    // Set remaining undeduced faces to NO_HYPOTHESIS
    for fi in 0..num_faces {
        if S::hypothesis_field(&mesh.faces[fi]) == S::undeduced_value() {
            S::set_hypothesis_field(&mut mesh.faces[fi], NO_HYPOTHESIS);
        }
    }

    if config.verbosity >= 2 {
        eprintln!(
            "  [{}-BFS] stats: {} seeds generated, {} fit failures, {} tried",
            S::name(), seeds_generated, fit_failures, seeds_tried
        );
        eprintln!(
            "    stol fail: {}, angle fail: {}, accept fail: {}, committed: {}",
            seed_stol_failures, seed_angle_failures, acceptance_failures, committed
        );
    }

    (hypotheses, false)
}

// ---------------------------------------------------------------------------
// Helper trait for accessing faces from a hypothesis
// ---------------------------------------------------------------------------

/// Trait to access the face list from a hypothesis struct.
pub(super) trait HypothesisFaces {
    fn faces(&self) -> &[usize];
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Check if two faces are neighbors (share an edge).
pub(super) fn is_neighbor(mesh: &ConnectedMesh, fi: usize, fj: usize) -> bool {
    let f = &mesh.faces[fi];
    let vc = f.vertex_count as usize;
    for &n in &f.neighbors[..vc] {
        if n >= 0 && n as usize == fj {
            return true;
        }
    }
    false
}
