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
    /// Viz detail level: 1=commits, 2=seeds, 3=BFS faces
    pub viz_level: u8,
    /// Per-face flag: true if the face belongs to a LARGE spherical
    /// hypothesis from stage 2.3. Torus seed gathering excludes only these:
    /// small spurious spheres (a corner vertex-blend patch locally fits a
    /// sphere of the tube radius) must not block torus seeding — the real
    /// torus beats them in selection once it exists. Empty means "no
    /// information" (unit tests). NOT redundant with `junction_faces`
    /// (verified 2026-07-15): sphere interiors are uniform-curvature and
    /// never junction-flagged, and torus seeding would tile them with
    /// bogus-but-valid band fits (pill_fine stage 2: 1.2 s → 429 s).
    pub large_sphere: Vec<bool>,
    /// Per-face curvature-difference score (see stage2::compute_seed_scores).
    /// Start faces are processed in ascending score order so interior faces
    /// of uniform regions seed before junction faces — contaminated junction
    /// seeds are usually absorbed by BFS growth before they ever run. Empty
    /// means "use mesh index order" (unit tests).
    pub seed_scores: Vec<f64>,
    /// Per-face flag: true if the face sits at a surface transition (high
    /// RAW curvature-difference score, pre-diffusion). Seed gathering
    /// excludes these in pass 1 so seeds stay on one surface; pass 2
    /// retries remaining undeduced faces without the filter, so small
    /// all-junction surfaces (one-row fillets) still get seeded. Empty
    /// means "no information" (unit tests).
    pub junction_faces: Vec<bool>,
}

impl BfsConfig {
    /// Whether `fi` should be excluded from seed gathering as a junction
    /// face. Always false in relaxed (pass-2) seeding or without field data.
    pub fn is_junction(&self, fi: usize, relaxed: bool) -> bool {
        !relaxed && self.junction_faces.get(fi).copied().unwrap_or(false)
    }
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

    /// Generate seed face sets anchored at face `fi`. `relaxed` is true on
    /// retry passes: junction-face exclusion is lifted so surfaces composed
    /// entirely of transition faces can still seed.
    fn generate_seeds(mesh: &ConnectedMesh, fi: usize, config: &BfsConfig, relaxed: bool)
        -> Vec<Seed>;

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

    /// Number of seeding passes over the mesh. Defaults to 1. Tori override
    /// to 2: corner vertex-blend patches seed before the adjacent larger
    /// fillets commit, so their 1-ring seeds are contaminated; a second pass
    /// over still-undeduced faces retries them with clean exclusions.
    fn seeding_passes() -> usize { 1 }

    /// Whether failed seed fits should retry after trimming relative-outlier
    /// faces. Defaults to false: for surface types with strong acceptance
    /// gates this salvages spurious mini-hypotheses from junk seeds (e.g.
    /// 10-face spheres on cone apex fans that then block cone seeding).
    /// Tori override to true — their vertex-neighborhood seeds legitimately
    /// mix faces from two adjacent fillet patches.
    fn seed_outlier_trimming() -> bool { false }

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

    // Two seeding passes: seeds tried early in pass 1 can be contaminated by
    // neighboring faces that belong to a *different* surface which has not
    // been committed yet (e.g. corner vertex-blend torus patches whose 1-ring
    // pulls in faces of the adjacent, larger fillet torus — the corner faces
    // have lower indices and seed first). Pass 2 revisits only faces still
    // undeduced, whose seed gathering now excludes everything committed in
    // pass 1.
    // Start-face order: ascending curvature-difference score (interior faces
    // of uniform regions first), ties broken by index for determinism.
    let seed_order: Vec<usize> = if config.seed_scores.len() == num_faces {
        let mut order: Vec<usize> = (0..num_faces).collect();
        order.sort_by(|&a, &b| {
            config.seed_scores[a]
                .partial_cmp(&config.seed_scores[b])
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.cmp(&b))
        });
        order
    } else {
        (0..num_faces).collect()
    };

    'passes: for pass in 0..S::seeding_passes() {
        if pass > 0 && config.verbosity >= 3 {
            eprintln!("    {} pass 2: retrying remaining undeduced faces", S::name());
        }
    for &fi in &seed_order {
        if !S::face_eligible(&mesh.faces[fi]) {
            continue;
        }

        prof_faces_visited += 1;
        // Periodic progress report every 2 seconds
        if config.verbosity >= 3 && prof_last_report.elapsed().as_secs_f64() > 2.0 {
            eprintln!(
                "    {} fi={}/{} ({:.1}%), elapsed={:.1}s, seeds={}, fits={}, refits={}, committed={}",
                S::name(), fi, num_faces, 100.0 * fi as f64 / num_faces as f64,
                prof_start.elapsed().as_secs_f64(),
                prof_total_seed_count, prof_total_fit_calls, prof_total_refit_calls, committed,
            );
            prof_last_report = std::time::Instant::now();
        }

        let t_seed_gen = std::time::Instant::now();
        let seeds = S::generate_seeds(mesh, fi, config, pass > 0);
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
        let mut best_seed_faces: Vec<usize> = Vec::new();
        let mut best_fit_desc = String::new();
        let mut best_fit: Option<S::Fit> = None;

        for seed in &seeds {
            if config.verbosity >= 4 {
                let seed_str: Vec<String> = seed.faces.iter().map(|f| f.to_string()).collect();
                eprintln!("      trying seed ({})", seed_str.join(", "));
            }
            let collect_vertices = |faces: &[usize]| -> Vec<usize> {
                let mut verts: Vec<usize> = Vec::new();
                for &sfi in faces {
                    let vc = mesh.faces[sfi].vertex_count as usize;
                    verts.extend_from_slice(&mesh.faces[sfi].vertex_indices[..vc]);
                }
                verts.sort_unstable();
                verts.dedup();
                verts
            };

            // Fit with outlier pruning: seeds gathered from vertex
            // neighborhoods can mix faces from two adjacent surfaces (e.g.
            // two tangent corner-fillet tori, where the leftover faces of
            // both patches interleave). When the fit fails vertex tolerance,
            // drop the faces whose vertices exceed it and refit from the
            // inliers, up to two pruning rounds.
            let seed_vtol = vtol * S::seed_vtol_multiplier();
            let mut seed_faces: Vec<usize> = seed.faces.clone();
            let mut seed_vertices = collect_vertices(&seed_faces);
            let mut fit_result: Option<S::Fit> = None;
            let mut seed_max_vtol = 0.0_f64;
            let max_rounds = if S::seed_outlier_trimming() { 4 } else { 1 };
            for _round in 0..max_rounds {
                let t_fit = std::time::Instant::now();
                let fit = S::fit_surface(mesh, &seed_faces, &seed_vertices, config.verbosity);
                prof_time_seed_fit += t_fit.elapsed();
                prof_total_fit_calls += 1;
                let fit = match fit {
                    Some(f) => f,
                    None => break,
                };
                seed_max_vtol = seed_vertices
                    .iter()
                    .map(|&vi| {
                        let v = &mesh.vertices[vi];
                        S::vertex_distance(&fit, &[v.x, v.y, v.z])
                    })
                    .fold(0.0_f64, f64::max);
                if seed_max_vtol <= seed_vtol {
                    fit_result = Some(fit);
                    break;
                }
                // Trim relative outliers: a contaminated seed produces a
                // compromise fit under which even the true-surface faces sit
                // above the (very tight) vertex tolerance, so absolute
                // thresholding finds no inliers. Instead drop the faces whose
                // error is a large fraction of the worst error and refit; a
                // genuinely bad seed has uniform errors, trims to nothing,
                // and is rejected. Only the final acceptance check above uses
                // the strict tolerance.
                let face_errs: Vec<f64> = seed_faces
                    .iter()
                    .map(|&sfi| {
                        let vc = mesh.faces[sfi].vertex_count as usize;
                        mesh.faces[sfi].vertex_indices[..vc]
                            .iter()
                            .map(|&vi| {
                                let v = &mesh.vertices[vi];
                                S::vertex_distance(&fit, &[v.x, v.y, v.z])
                            })
                            .fold(0.0_f64, f64::max)
                    })
                    .collect();
                let max_face_err = face_errs.iter().copied().fold(0.0_f64, f64::max);
                let keep_thresh = seed_vtol.max(max_face_err * 0.5);
                let inliers: Vec<usize> = seed_faces
                    .iter()
                    .copied()
                    .zip(face_errs.iter())
                    .filter(|&(_, &e)| e <= keep_thresh)
                    .map(|(f, _)| f)
                    .collect();
                if inliers.len() < 3 || inliers.len() == seed_faces.len() {
                    break;
                }
                if config.verbosity >= 4 {
                    eprintln!(
                        "      seed pruned {} outlier face(s), refitting from {}",
                        seed_faces.len() - inliers.len(),
                        inliers.len(),
                    );
                }
                seed_faces = inliers;
                seed_vertices = collect_vertices(&seed_faces);
            }
            seeds_tried += 1;
            let fit = match fit_result {
                Some(f) => f,
                None => {
                    fit_failures += 1;
                    if config.verbosity >= 4 {
                        if seed_max_vtol > 0.0 {
                            eprintln!(
                                "      seed rejected: max vertex err {:.6e} > limit {:.6e}",
                                seed_max_vtol, seed_vtol
                            );
                        } else {
                            eprintln!("      seed fit failed: no valid fit");
                        }
                    }
                    continue;
                }
            };

            let mut seed_max_stol = 0.0_f64;
            for &sfi in &seed_faces {
                let c = face_centroid(&mesh.faces[sfi], &mesh.vertices);
                let d = S::centroid_distance(&fit, &c);
                if d > seed_max_stol {
                    seed_max_stol = d;
                }
            }
            if seed_max_stol > stol {
                if config.verbosity >= 4 {
                    eprintln!(
                        "      seed rejected: max surface err {:.6e} > limit {:.6e}",
                        seed_max_stol, stol
                    );
                }
                seed_stol_failures += 1;
                continue;
            }

            let mut angle_ok = true;
            'angle_check: for i in 0..seed_faces.len() {
                for j in (i + 1)..seed_faces.len() {
                    if is_neighbor(mesh, seed_faces[i], seed_faces[j])
                        && dihedral_angle(mesh, seed_faces[i], seed_faces[j]) > angular_tol
                    {
                        angle_ok = false;
                        break 'angle_check;
                    }
                }
            }
            if !angle_ok {
                if config.verbosity >= 4 {
                    eprintln!("      seed rejected: dihedral angle exceeds limit");
                }
                seed_angle_failures += 1;
                continue;
            }

            let mut current_fit = fit;
            S::determine_convexity(mesh, &seed_faces, &mut current_fit);

            if config.verbosity >= 4 {
                let seed_str: Vec<String> = seed.faces.iter().map(|f| f.to_string()).collect();
                eprintln!(
                    "      seed accepted ({}): {} vertices, max_vtx_err={:.6e}, max_sfc_err={:.6e}, {}",
                    seed_str.join(","), seed_vertices.len(),
                    seed_max_vtol, seed_max_stol,
                    S::viz_fit_description(&current_fit),
                );
            }

            // --- Viz: show seed faces (viz_level >= 2) ---
            let mut skip_viz = false;
            if config.viz_level >= 2 {
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
            }

            // --- BFS expansion ---
            let t_bfs = std::time::Instant::now();
            let mut trial_faces: Vec<usize> = seed_faces.clone();
            let mut trial_face_set: HashSet<usize> = seed_faces.iter().copied().collect();
            let mut trial_vertex_set: HashSet<usize> = seed_vertices.iter().copied().collect();
            let mut trial_vertices: Vec<usize> = seed_vertices.clone();

            let mut bfs_queue = VecDeque::new();
            for &sfi in &seed_faces {
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
                        if config.verbosity >= 5 {
                            eprintln!("        face {}: rejected (dihedral angle exceeds limit)", ni);
                        }
                        continue;
                    }

                    if !S::convexity_check(&current_fit, mesh, ni) {
                        if config.verbosity >= 5 {
                            eprintln!("        face {}: rejected (convexity mismatch)", ni);
                        }
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
                        if config.verbosity >= 5 {
                            eprintln!(
                                "        face {}: accepted (vtx_err={:.6e}, sfc_err={:.6e})",
                                ni, max_vtx_dist, ctr_dist,
                            );
                        }
                    } else if max_vtx_dist > skip_vtol || ctr_dist > skip_stol {
                        if config.verbosity >= 5 {
                            eprintln!(
                                "        face {}: skipped (vtx_err={:.6e}, sfc_err={:.6e}, beyond refit range)",
                                ni, max_vtx_dist, ctr_dist,
                            );
                        }
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
                            None => {
                                if config.verbosity >= 5 {
                                    eprintln!("        face {}: refit failed (no valid fit)", ni);
                                }
                                prof_time_refit += t_refit.elapsed(); prof_total_refit_calls += 1; continue;
                            },
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
                            if config.verbosity >= 5 {
                                eprintln!("        face {}: rejected after refit (tolerance exceeded)", ni);
                            }
                            prof_time_refit += t_refit.elapsed();
                            continue;
                        }

                        // Preserve convexity across refit
                        S::determine_convexity(mesh, &trial_faces, &mut refit);
                        current_fit = refit;
                        prof_time_refit += t_refit.elapsed();
                        if config.verbosity >= 5 {
                            // Compute refit max errors for logging
                            let mut refit_max_vtx = 0.0_f64;
                            let mut refit_max_ctr = 0.0_f64;
                            for &tfi in &test_faces {
                                let tf = &mesh.faces[tfi];
                                let tvc = tf.vertex_count as usize;
                                for &vi in &tf.vertex_indices[..tvc] {
                                    let v = &mesh.vertices[vi];
                                    let d = S::vertex_distance(&current_fit, &[v.x, v.y, v.z]);
                                    refit_max_vtx = refit_max_vtx.max(d);
                                }
                                let tc = face_centroid(&mesh.faces[tfi], &mesh.vertices);
                                let td = S::centroid_distance(&current_fit, &tc);
                                refit_max_ctr = refit_max_ctr.max(td);
                            }
                            eprintln!(
                                "        face {}: accepted after refit (vtx_err={:.6e}, sfc_err={:.6e}, {})",
                                ni, refit_max_vtx, refit_max_ctr,
                                S::viz_fit_description(&current_fit),
                            );
                        }
                    }

                    trial_faces.push(ni);
                    trial_face_set.insert(ni);
                    for &vi in &new_verts {
                        trial_vertex_set.insert(vi);
                        trial_vertices.push(vi);
                    }
                    bfs_queue.push_back(ni);

                    // --- Viz: show accepted face (viz_level >= 3) ---
                    if config.viz_level >= 3 {
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
            }
            prof_time_bfs += t_bfs.elapsed();
            prof_total_bfs_faces += trial_faces.len() as u64;
            if trial_faces.len() > prof_max_bfs_faces {
                prof_max_bfs_faces = trial_faces.len();
            }

            // Final re-fit. Adopt the fresh fit only if it beats the
            // incrementally-maintained fit on the accumulated trial
            // vertices: fitting from scratch can drift to a fit whose
            // error is orders of magnitude above tolerance (observed:
            // tilted r=3.37 "cylinders" with 2.6e-3 mm error committed on
            // corner wall patches, which then poisoned all torus seeding
            // in the region), or converge to a wrong local minimum even
            // when current_fit is exact (observed: the full-ring R=3.5
            // pocket-rim fillet torus on stepped_block MEDIUM refits to
            // 0.87 mm error although every trial vertex is within 1e-5 of
            // current_fit). current_fit already satisfies the vertex
            // tolerance for every accepted face by construction, so it is
            // always a valid fallback.
            let trial_max_err = |fit: &S::Fit| -> f64 {
                trial_vertices
                    .iter()
                    .map(|&vi| {
                        let v = &mesh.vertices[vi];
                        S::vertex_distance(fit, &[v.x, v.y, v.z])
                    })
                    .fold(0.0_f64, f64::max)
            };
            let t_final = std::time::Instant::now();
            let mut final_max_err = trial_max_err(&current_fit);
            if let Some(final_fit) = S::fit_surface(
                mesh,
                &trial_faces,
                &trial_vertices,
                0,
            ) {
                let refit_err = trial_max_err(&final_fit);
                if refit_err < final_max_err {
                    current_fit = final_fit;
                    final_max_err = refit_err;
                }
            }

            // (Re-)determine convexity after final refit
            S::determine_convexity(mesh, &trial_faces, &mut current_fit);
            prof_time_final_fit += t_final.elapsed();
            if final_max_err > vtol * S::seed_vtol_multiplier() {
                if config.verbosity >= 4 {
                    eprintln!(
                        "      trial rejected after final refit: max vertex err {:.6e} > {:.6e}",
                        final_max_err,
                        vtol * S::seed_vtol_multiplier(),
                    );
                }
                acceptance_failures += 1;
                continue;
            }

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
                best_seed_faces = seed_faces.clone();
                best_fit_desc = S::viz_fit_description(&current_fit);
                best_fit = Some(current_fit.clone());
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
            if config.verbosity >= 3 {
                let seed_str: Vec<String> = best_seed_faces.iter().map(|f| f.to_string()).collect();
                eprintln!(
                    "    committed {} hyp {}: {} faces, seed=({}), {}",
                    S::name(), hyp_idx, hyp.faces().len(), seed_str.join(","), best_fit_desc,
                );
            }

            // --- Viz: show committed hypothesis (viz_level >= 1) ---
            if config.viz_level >= 1 {
                if let Some(viz_sender) = viz {
                    let hyp_faces: Vec<usize> = hyp.faces().to_vec();
                    let (cylinders, cones, spheres, tori) = if let Some(ref fit) = best_fit {
                        S::viz_overlays(fit, mesh, &hyp_faces)
                    } else {
                        (Vec::new(), Vec::new(), Vec::new(), Vec::new())
                    };
                    let mut highlights = vec![
                        viz::FaceHighlight { face_indices: hyp_faces.clone(), color: [0.2, 0.8, 0.2, 1.0] },
                    ];
                    // Gray background for all other committed faces
                    let bg_faces: Vec<usize> = (0..mesh.faces.len())
                        .filter(|f| {
                            let h = S::hypothesis_field(&mesh.faces[*f]);
                            h >= 0 && !hyp_faces.contains(f)
                        })
                        .collect();
                    if !bg_faces.is_empty() {
                        highlights.push(viz::FaceHighlight { face_indices: bg_faces, color: [0.5, 0.5, 0.5, 1.0] });
                    }
                    let mut overlay = viz::VizOverlay::new();
                    overlay.status_text = format!(
                        "BFS-{}: committed hyp {} ({} faces) {} [space=next, click=skip-to-face]",
                        S::name(), hyp_idx, hyp.faces().len(), best_fit_desc,
                    );
                    overlay.face_highlights = highlights;
                    overlay.cylinders = cylinders;
                    overlay.cones = cones;
                    overlay.spheres = spheres;
                    overlay.tori = tori;
                    if !hyp_faces.is_empty() {
                        overlay.focus_point = Some(viz_face_centroid(hyp_faces[0], mesh));
                        overlay.focus_normal = viz_face_normal(hyp_faces[0], mesh);
                    }
                    match viz_sender.show_and_wait(overlay) {
                        VizAction::Quit => return (hypotheses, true),
                        VizAction::SkipToSeed(face) => { skip_to_face = Some(face); }
                        _ => {}
                    }
                }
            }

            hypotheses.push(hyp);
        }
    }

        // Skip pass 2 when pass 1 covered every face (common case).
        if (0..num_faces).all(|fi| !S::face_eligible(&mesh.faces[fi])) {
            break 'passes;
        }
    }

    // Set remaining undeduced faces to NO_HYPOTHESIS
    for fi in 0..num_faces {
        if S::hypothesis_field(&mesh.faces[fi]) == S::undeduced_value() {
            S::set_hypothesis_field(&mut mesh.faces[fi], NO_HYPOTHESIS);
        }
    }

    if config.verbosity >= 3 {
        eprintln!(
            "    {} stats: {} seeds generated, {} fit failures, {} tried",
            S::name(), seeds_generated, fit_failures, seeds_tried
        );
        eprintln!(
            "      stol fail: {}, angle fail: {}, accept fail: {}, committed: {}",
            seed_stol_failures, seed_angle_failures, acceptance_failures, committed
        );
    }

    // PROFILING: final summary
    let prof_total = prof_start.elapsed();
    if config.verbosity >= 3 {
    eprintln!(
        "    {} profile: total={:.3}s, faces_visited={}, total_seeds={}, max_seeds_per_face={}",
        S::name(), prof_total.as_secs_f64(), prof_faces_visited, prof_total_seed_count, prof_max_seeds_per_face,
    );
    eprintln!(
        "      seed_gen={:.3}s, seed_fit={:.3}s ({} calls), bfs={:.3}s (total_bfs_faces={}, max={})",
        prof_time_seed_gen.as_secs_f64(), prof_time_seed_fit.as_secs_f64(), prof_total_fit_calls,
        prof_time_bfs.as_secs_f64(), prof_total_bfs_faces, prof_max_bfs_faces,
    );
    eprintln!(
        "      refit={:.3}s ({} calls), final_fit={:.3}s, acceptance={:.3}s",
        prof_time_refit.as_secs_f64(), prof_total_refit_calls,
        prof_time_final_fit.as_secs_f64(), prof_time_acceptance.as_secs_f64(),
    );
    eprintln!(
        "      avg_seed_fit={:.3}ms, avg_refit={:.3}ms",
        if prof_total_fit_calls > 0 { prof_time_seed_fit.as_secs_f64() * 1000.0 / prof_total_fit_calls as f64 } else { 0.0 },
        if prof_total_refit_calls > 0 { prof_time_refit.as_secs_f64() * 1000.0 / prof_total_refit_calls as f64 } else { 0.0 },
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
