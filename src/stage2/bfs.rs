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
use crate::viz::{self, VizAction, VizSender};

use super::{dot3, face_area, face_centroid, viz_face_centroid, viz_face_normal, REFIT_SKIP_MULTIPLIER};

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

    /// Multiplier for seed vertex tolerance.  Defaults to 1.0.
    /// Cone seeds are inherently noisier (mixed-face apex neighborhoods)
    /// so they override this to allow a looser seed threshold.
    fn seed_vtol_multiplier() -> f64 { 1.0 }

    /// Multiplier for the angular tolerance used in dihedral-angle checks.
    /// Defaults to 1.0.  Cones override to 2.0 because tessellated cone tips
    /// have inter-strip dihedral angles that exceed the default tolerance.
    fn angular_tol_multiplier() -> f64 { 1.0 }

    /// Generate visualization overlays (cylinder/cone/sphere/torus) from the current fit.
    /// Default: no overlays (e.g., for planar surfaces).
    fn viz_overlays(
        _fit: &Self::Fit,
        _mesh: &ConnectedMesh,
        _face_list: &[usize],
    ) -> (Vec<viz::CylinderOverlay>, Vec<viz::ConeOverlay>, Vec<viz::SphereOverlay>, Vec<viz::TorusOverlay>) {
        (Vec::new(), Vec::new(), Vec::new(), Vec::new())
    }

    /// Short description of the current fit for viz status text.
    /// Default: just the surface type name.
    fn viz_fit_description(_fit: &Self::Fit) -> String {
        Self::name().to_string()
    }

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
    viz: Option<&VizSender>,
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
    let angular_tol = config.angular_tol * S::angular_tol_multiplier();

    let mut hypotheses: Vec<S::Hypothesis> = Vec::new();

    // Statistics
    let mut seeds_generated = 0_u32;
    let mut seeds_tried = 0_u32;
    let mut fit_failures = 0_u32;
    let mut seed_stol_failures = 0_u32;
    let mut seed_angle_failures = 0_u32;
    let mut acceptance_failures = 0_u32;
    let mut committed = 0_u32;

    // PROFILING: timing instrumentation
    let prof_start = std::time::Instant::now();
    let mut prof_faces_visited = 0_u64;
    let mut prof_total_seed_count = 0_u64;
    let mut prof_total_bfs_faces = 0_u64;
    let mut prof_total_fit_calls = 0_u64;
    let mut prof_total_refit_calls = 0_u64;
    let mut prof_time_seed_gen = std::time::Duration::ZERO;
    let mut prof_time_seed_fit = std::time::Duration::ZERO;
    let mut prof_time_bfs = std::time::Duration::ZERO;
    let mut prof_time_refit = std::time::Duration::ZERO;
    let mut prof_time_final_fit = std::time::Duration::ZERO;
    let mut prof_time_acceptance = std::time::Duration::ZERO;
    let mut prof_max_seeds_per_face = 0_usize;
    let mut prof_max_bfs_faces = 0_usize;
    let mut prof_last_report = std::time::Instant::now();

    // Initialize hypothesis fields to UNDEDUCED
    for fi in 0..num_faces {
        S::set_hypothesis_field(&mut mesh.faces[fi], S::undeduced_value());
    }

    // Viz: if the user clicked a face, skip seeds until we reach one containing that face
    let mut skip_to_face: Option<usize> = None;

    for fi in 0..num_faces {
        if !S::face_eligible(&mesh.faces[fi]) {
            continue;
        }

        prof_faces_visited += 1;
        // PROFILING: periodic progress report every 2 seconds
        if prof_last_report.elapsed().as_secs_f64() > 2.0 {
            eprintln!(
                "  [{}-BFS PROFILE] fi={}/{} ({:.1}%), elapsed={:.1}s, seeds_so_far={}, fits={}, refits={}, committed={}",
                S::name(), fi, num_faces, 100.0 * fi as f64 / num_faces as f64,
                prof_start.elapsed().as_secs_f64(),
                prof_total_seed_count, prof_total_fit_calls, prof_total_refit_calls, committed,
            );
            prof_last_report = std::time::Instant::now();
        }

        let t_seed_gen = std::time::Instant::now();
        let seeds = S::generate_seeds(mesh, fi, config);
        prof_time_seed_gen += t_seed_gen.elapsed();
        if seeds.is_empty() {
            continue;
        }
        prof_total_seed_count += seeds.len() as u64;
        if seeds.len() > prof_max_seeds_per_face {
            prof_max_seeds_per_face = seeds.len();
        }
        seeds_generated += seeds.len() as u32;

        let mut best_area = -1.0_f64;
        let mut best_hypothesis: Option<S::Hypothesis> = None;

        for seed in &seeds {
            if config.verbosity >= 3 {
                let seed_str: Vec<String> = seed.faces.iter().map(|f| f.to_string()).collect();
                eprintln!("    [{}-BFS] trying seed ({})", S::name(), seed_str.join(", "));
            }
            let mut seed_vertex_set: HashSet<usize> = HashSet::new();
            for &sfi in &seed.faces {
                let vc = mesh.faces[sfi].vertex_count as usize;
                for &vi in &mesh.faces[sfi].vertex_indices[..vc] {
                    seed_vertex_set.insert(vi);
                }
            }
            let seed_vertices: Vec<usize> = seed_vertex_set.iter().copied().collect();

            let t_fit = std::time::Instant::now();
            let fit = match S::fit_surface(mesh, &seed.faces, &seed_vertices, config.verbosity) {
                Some(f) => f,
                None => {
                    prof_time_seed_fit += t_fit.elapsed();
                    prof_total_fit_calls += 1;
                    fit_failures += 1;
                    if config.verbosity >= 3 {
                        eprintln!("    [{}-BFS] seed fit failed: no valid fit", S::name());
                    }
                    continue;
                }
            };
            prof_time_seed_fit += t_fit.elapsed();
            prof_total_fit_calls += 1;
            seeds_tried += 1;

            // Validate seed: vertex tolerance (using surface-specific multiplier)
            let seed_vtol = vtol * S::seed_vtol_multiplier();
            let mut seed_valid = true;
            let mut seed_max_vtol = 0.0_f64;
            for &vi in &seed_vertices {
                let v = &mesh.vertices[vi];
                let d = S::vertex_distance(&fit, &[v.x, v.y, v.z]);
                if d > seed_max_vtol {
                    seed_max_vtol = d;
                }
            }
            if seed_max_vtol > seed_vtol {
                seed_valid = false;
            }
            if !seed_valid {
                if config.verbosity >= 3 {
                    eprintln!(
                        "    [{}-BFS] seed rejected: vertex err {:.6e} > limit {:.6e}",
                        S::name(), seed_max_vtol, seed_vtol
                    );
                }
                continue;
            }

            let mut seed_max_stol = 0.0_f64;
            for &sfi in &seed.faces {
                let c = face_centroid(&mesh.faces[sfi], &mesh.vertices);
                let d = S::centroid_distance(&fit, &c);
                if d > seed_max_stol {
                    seed_max_stol = d;
                }
            }
            if seed_max_stol > stol {
                seed_valid = false;
            }
            if !seed_valid {
                if config.verbosity >= 3 {
                    eprintln!(
                        "    [{}-BFS] seed rejected: surface err {:.6e} > limit {:.6e}",
                        S::name(), seed_max_stol, stol
                    );
                }
                seed_stol_failures += 1;
                continue;
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
                if config.verbosity >= 3 {
                    eprintln!("    [{}-BFS] seed rejected: dihedral angle exceeds limit", S::name());
                }
                seed_angle_failures += 1;
                continue;
            }

            let mut current_fit = fit;
            S::determine_convexity(mesh, &seed.faces, &mut current_fit);

            // --- Viz: show seed faces ---
            let mut skip_viz = false;
            if let Some(viz_sender) = viz {
                // Auto-skip if we're seeking a specific face
                if let Some(target) = skip_to_face {
                    if !seed.faces.contains(&target) {
                        skip_viz = true;
                    } else {
                        skip_to_face = None; // reached target, stop skipping
                    }
                }
                if !skip_viz {
                let seed_str: Vec<String> = seed.faces.iter().map(|f| f.to_string()).collect();
                let (cylinders, cones, spheres, tori) = S::viz_overlays(&current_fit, mesh, &seed.faces);
                let mut highlights = vec![
                    viz::FaceHighlight { face_indices: vec![seed.faces[0]], color: [1.0, 0.6, 0.0, 1.0] },
                ];
                if seed.faces.len() > 1 {
                    highlights.push(viz::FaceHighlight { face_indices: seed.faces[1..].to_vec(), color: [1.0, 0.8, 0.4, 1.0] });
                }
                // Gray background for already-committed faces
                let bg_faces: Vec<usize> = (0..mesh.faces.len())
                    .filter(|f| {
                        let h = S::hypothesis_field(&mesh.faces[*f]);
                        h >= 0 && !seed.faces.contains(f)
                    })
                    .collect();
                if !bg_faces.is_empty() {
                    highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                }
                let mut overlay = viz::VizOverlay::new();
                overlay.status_text = format!(
                    "BFS-{}: seed=({}) {} [space=step, shift+space=skip, click=skip-to-face]",
                    S::name(), seed_str.join(","), S::viz_fit_description(&current_fit),
                );
                overlay.face_highlights = highlights;
                overlay.cylinders = cylinders;
                overlay.cones = cones;
                overlay.spheres = spheres;
                overlay.tori = tori;
                overlay.focus_point = Some(viz_face_centroid(seed.faces[0], mesh));
                overlay.focus_normal = viz_face_normal(seed.faces[0], mesh);
                match viz_sender.show_and_wait(overlay) {
                    VizAction::Quit => return (hypotheses, true),
                    VizAction::NextSeed => { skip_viz = true; }
                    VizAction::SkipToSeed(face) => { skip_viz = true; skip_to_face = Some(face); }
                    VizAction::NextStep => {}
                }
                }
            }

            // --- BFS expansion ---
            let t_bfs = std::time::Instant::now();
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

                    if !S::convexity_check(&current_fit, mesh, ni) {
                        continue;
                    }

                    let ni_vi = mesh.faces[ni].vertex_indices;
                    let mut new_verts: Vec<usize> = Vec::new();
                    for &vi in &ni_vi[..ni_vc] {
                        if !trial_vertex_set.contains(&vi) {
                            new_verts.push(vi);
                        }
                    }

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
                        // Accept directly
                    } else if max_vtx_dist > skip_vtol || ctr_dist > skip_stol {
                        continue;
                    } else {
                        let t_refit = std::time::Instant::now();
                        let mut test_faces = trial_faces.clone();
                        test_faces.push(ni);
                        let mut test_verts = trial_vertices.clone();
                        for &vi in &new_verts {
                            test_verts.push(vi);
                        }

                        let mut refit = match S::fit_surface(mesh, &test_faces, &test_verts, 0) {
                            Some(f) => f,
                            None => { prof_time_refit += t_refit.elapsed(); prof_total_refit_calls += 1; continue; },
                        };
                        prof_total_refit_calls += 1;

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
                            prof_time_refit += t_refit.elapsed();
                            continue;
                        }

                        // Preserve convexity across refit
                        S::determine_convexity(mesh, &trial_faces, &mut refit);
                        current_fit = refit;
                        prof_time_refit += t_refit.elapsed();
                    }

                    trial_faces.push(ni);
                    trial_face_set.insert(ni);
                    for &vi in &new_verts {
                        trial_vertex_set.insert(vi);
                        trial_vertices.push(vi);
                    }
                    bfs_queue.push_back(ni);

                    // --- Viz: show accepted face ---
                    if let Some(viz_sender) = viz {
                        if !skip_viz {
                        let accepted_nonseed: Vec<usize> = trial_faces.iter()
                            .filter(|f| !seed.faces.contains(f) && **f != ni)
                            .copied().collect();
                        let (cylinders, cones, spheres, tori) = S::viz_overlays(&current_fit, mesh, &trial_faces);
                        let mut highlights = vec![
                            viz::FaceHighlight { face_indices: vec![seed.faces[0]], color: [1.0, 0.6, 0.0, 1.0] },
                        ];
                        if seed.faces.len() > 1 {
                            highlights.push(viz::FaceHighlight { face_indices: seed.faces[1..].to_vec(), color: [1.0, 0.8, 0.4, 1.0] });
                        }
                        if !accepted_nonseed.is_empty() {
                            highlights.push(viz::FaceHighlight { face_indices: accepted_nonseed, color: [0.2, 0.4, 1.0, 1.0] });
                        }
                        highlights.push(viz::FaceHighlight { face_indices: vec![ni], color: [0.1, 0.2, 0.7, 1.0] });
                        // Gray background for already-committed faces
                        let bg_faces: Vec<usize> = (0..mesh.faces.len())
                            .filter(|f| {
                                let h = S::hypothesis_field(&mesh.faces[*f]);
                                h >= 0 && !trial_face_set.contains(f)
                            })
                            .collect();
                        if !bg_faces.is_empty() {
                            highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                        }
                        let mut overlay = viz::VizOverlay::new();
                        overlay.status_text = format!(
                            "BFS-{}: accepted fi={} ({} faces) {} [space=step, shift+space=skip]",
                            S::name(), ni, trial_faces.len(), S::viz_fit_description(&current_fit),
                        );
                        overlay.face_highlights = highlights;
                        overlay.cylinders = cylinders;
                        overlay.cones = cones;
                        overlay.spheres = spheres;
                        overlay.tori = tori;
                        overlay.focus_point = Some(viz_face_centroid(ni, mesh));
                        overlay.focus_normal = viz_face_normal(ni, mesh);
                        match viz_sender.show_and_wait(overlay) {
                            VizAction::Quit => return (hypotheses, true),
                            VizAction::NextSeed => { skip_viz = true; }
                            VizAction::SkipToSeed(face) => { skip_viz = true; skip_to_face = Some(face); }
                            VizAction::NextStep => {}
                        }
                        }
                    }
                }
            }
            prof_time_bfs += t_bfs.elapsed();
            prof_total_bfs_faces += trial_faces.len() as u64;
            if trial_faces.len() > prof_max_bfs_faces {
                prof_max_bfs_faces = trial_faces.len();
            }

            // Final re-fit
            let t_final = std::time::Instant::now();
            if let Some(final_fit) = S::fit_surface(
                mesh,
                &trial_faces,
                &trial_vertices,
                0,
            ) {
                current_fit = final_fit;
            }

            // (Re-)determine convexity after final refit
            S::determine_convexity(mesh, &trial_faces, &mut current_fit);
            prof_time_final_fit += t_final.elapsed();

            let t_accept = std::time::Instant::now();
            if !S::acceptance_test(mesh, &trial_faces, &current_fit, config) {
                prof_time_acceptance += t_accept.elapsed();
                acceptance_failures += 1;
                continue;
            }
            prof_time_acceptance += t_accept.elapsed();

            let total_area: f64 = trial_faces
                .iter()
                .map(|&tfi| face_area(&mesh.faces[tfi], &mesh.vertices))
                .sum();

            if total_area > best_area {
                best_area = total_area;
                let mut sorted_verts: Vec<usize> = trial_vertex_set.into_iter().collect();
                sorted_verts.sort_unstable();
                best_hypothesis =
                    Some(S::make_hypothesis(mesh, trial_faces, sorted_verts, &current_fit));
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

    // PROFILING: final summary
    let prof_total = prof_start.elapsed();
    eprintln!(
        "  [{}-BFS PROFILE SUMMARY] total={:.3}s, faces_visited={}, total_seeds={}, max_seeds_per_face={}",
        S::name(), prof_total.as_secs_f64(), prof_faces_visited, prof_total_seed_count, prof_max_seeds_per_face,
    );
    eprintln!(
        "    seed_gen={:.3}s, seed_fit={:.3}s ({} calls), bfs={:.3}s (total_bfs_faces={}, max={})",
        prof_time_seed_gen.as_secs_f64(), prof_time_seed_fit.as_secs_f64(), prof_total_fit_calls,
        prof_time_bfs.as_secs_f64(), prof_total_bfs_faces, prof_max_bfs_faces,
    );
    eprintln!(
        "    refit={:.3}s ({} calls), final_fit={:.3}s, acceptance={:.3}s",
        prof_time_refit.as_secs_f64(), prof_total_refit_calls,
        prof_time_final_fit.as_secs_f64(), prof_time_acceptance.as_secs_f64(),
    );
    eprintln!(
        "    avg_seed_fit={:.3}ms, avg_refit={:.3}ms",
        if prof_total_fit_calls > 0 { prof_time_seed_fit.as_secs_f64() * 1000.0 / prof_total_fit_calls as f64 } else { 0.0 },
        if prof_total_refit_calls > 0 { prof_time_refit.as_secs_f64() * 1000.0 / prof_total_refit_calls as f64 } else { 0.0 },
    );

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
