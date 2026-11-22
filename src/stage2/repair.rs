//! Stage 2.4: resolution & repair.
//!
//! Turns stage 2.3's per-region hypotheses into a strict, connected,
//! validated partition. All operations are fit-validated and deterministic
//! (fixed iteration order — by surface index, smallest first for merges —
//! bounded rounds).
//!
//! - **Boundary relabel** (±1-row polish): move a region-boundary face to an
//!   adjacent surface that fits its vertices better, refitting every surface a
//!   move touched before the next round — assignment and fit are one problem,
//!   since a surface holding its neighbour's row is pulled off true by exactly
//!   that row. Iterates to a fixpoint, bounded by each face being allowed to
//!   move at most once.
//! - **Boundary snap**: put every *tangent* junction on the mesh-edge chain
//!   the two fitted surfaces' junction geometry picks out (`locus.rs`), rather
//!   than on whichever chain a residual comparison happens to prefer. Runs
//!   after relabel and merge, so the exact criterion has the last word on
//!   whatever survives them; junctions it cannot build a locus for are
//!   reported, not guessed.
//! - **Merge**: absorb a small region S into an adjacent, larger, real
//!   surface T when every vertex of S lies within `--vertex-tolerance` of T's
//!   fitted surface and every centroid within `--surface-tolerance`; refit
//!   the union, revalidate, keep only if valid. Eats junction-circle
//!   singletons and chord-degenerate slivers with a principled criterion —
//!   vertex tolerance, at which "on the same surface" is decisive.
//! - **Connectivity re-split**: relabel/merge can disconnect a surface;
//!   re-split into mesh-connected components (identical fitted geometry on
//!   two non-adjacent regions is legal).
//! - **Invariants** (hard checks before stage 3): total coverage and
//!   per-surface mesh connectivity.
//!
//! Claim resolution (multi-claimed junction rows) is a no-op here until
//! stage 2.3(b) scoped competition exists — the whole-region ladder produces
//! disjoint regions, so no face is claimed twice.

use crate::config::Config;
use crate::stage1::{ConnectedMesh, MeshVertex};

use std::collections::BTreeSet;

use super::classify::{fit_validate_joined, Classification};
use super::locus::{analyze_pairs, PairOutcome, Surf};
use super::surface::SurfaceType;
use super::conical::ConicalSurface;
use super::planar::PlanarSurface;
use super::spherical::SphericalSurface;
use super::toroidal::ToroidalSurface;
use super::{
    face_centroid, vertex_distance_to_selected, ConicalHypothesis, PlanarHypothesis,
    SelectedSurface, SphericalHypothesis, ToroidalHypothesis,
};

/// Largest region (in faces) eligible to be absorbed by a merge. Junction
/// singletons and chord-degenerate slivers are small; a real surface is not
/// dissolved into a neighbor.
const MERGE_MAX_SMALL: usize = 16;

/// Repair-operation counts, for logging.
#[derive(Default)]
pub(super) struct RepairStats {
    pub snapped: usize,
    pub snap_problems: usize,
    pub relabeled: usize,
    pub merged: usize,
    pub split: usize,
}

/// Run the stage-2.4 repair sequence in place, then assert the hard
/// partition invariants.
pub(super) fn repair(
    mesh: &ConnectedMesh,
    cls: &mut Classification,
    face_region: &[usize],
    config: &Config,
) -> RepairStats {
    if config.boundary_join() {
        cls.joined.resize(mesh.faces.len(), false);
    }
    let centroids: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|f| face_centroid(f, &mesh.vertices))
        .collect();

    // Snap runs late, after relabel and merge: relabel is a residual
    // comparison that gets the boundary approximately right everywhere,
    // including the junctions snap does not own, merge dissolves the slivers
    // that would otherwise be reported as junctions with no locus, and snap is
    // the exact criterion that should have the last word on what survives.
    // Re-split still follows it, since a snap can in principle disconnect a
    // surface.
    let relabeled = relabel_boundaries(mesh, cls, config);
    let merged = merge_regions(mesh, cls, &centroids, face_region, config)
        + merge_same_surface(mesh, cls, &centroids, config);
    let (snapped, snap_problems) = snap_boundaries(mesh, cls, config);
    let stats = RepairStats {
        snapped,
        snap_problems,
        relabeled,
        merged,
        split: resplit_disconnected(mesh, cls),
    };
    drop_empty_surfaces(cls);
    assert_invariants(mesh, cls);
    stats
}

// ---------------------------------------------------------------------------
// Surface accessors
// ---------------------------------------------------------------------------

fn surface_faces(cls: &Classification, si: usize) -> &[usize] {
    match cls.selected[si] {
        SelectedSurface::Planar(i) => &cls.planar[i].faces,
        SelectedSurface::Conical(i) => &cls.conical[i].faces,
        SelectedSurface::Spherical(i) => &cls.spherical[i].faces,
        SelectedSurface::Toroidal(i) => &cls.toroidal[i].faces,
        SelectedSurface::BSpline(i) => &cls.bspline[i].faces,
    }
}

fn surface_faces_mut(cls: &mut Classification, si: usize) -> &mut Vec<usize> {
    match cls.selected[si] {
        SelectedSurface::Planar(i) => &mut cls.planar[i].faces,
        SelectedSurface::Conical(i) => &mut cls.conical[i].faces,
        SelectedSurface::Spherical(i) => &mut cls.spherical[i].faces,
        SelectedSurface::Toroidal(i) => &mut cls.toroidal[i].faces,
        SelectedSurface::BSpline(i) => &mut cls.bspline[i].faces,
    }
}

/// A B-spline surface is inert in repair: its fit vouches only for its own
/// data, so relabel/merge/re-split never move faces into or out of it (a
/// moved face would invalidate the fit and its cached projector). All three
/// operations skip B-spline surfaces via this predicate.
fn is_bspline(cls: &Classification, si: usize) -> bool {
    matches!(cls.selected[si], SelectedSurface::BSpline(_))
}

/// Whether an analytic surface's own strict facets determine it — at least the
/// per-type minimum (`CompeteHyp::min_support`: two for a plane, since one
/// facet is exactly a plane). A B-spline is never counted: its fit vouches only
/// for its own data.
fn determined(cls: &Classification, si: usize) -> bool {
    let min = match cls.selected[si] {
        SelectedSurface::Planar(_) => 2,
        SelectedSurface::Conical(_) => super::MIN_CYLINDER_FACES,
        SelectedSurface::Spherical(_) => super::MIN_SPHERE_FACES,
        SelectedSurface::Toroidal(_) => super::MIN_TORUS_FACES,
        SelectedSurface::BSpline(_) => return false,
    };
    surface_faces(cls, si).iter().filter(|&&fi| !cls.is_joined(fi)).count() >= min
}

/// Distance from a point to a selected surface's fitted geometry.
fn point_dist(cls: &Classification, si: usize, p: &[f64; 3]) -> f64 {
    let v = MeshVertex::from_xyz(p[0], p[1], p[2]);
    vertex_distance_to_selected(
        &v,
        &cls.selected[si],
        &cls.planar,
        &cls.conical,
        &cls.spherical,
        &cls.toroidal,
        &cls.bspline,
    )
}

/// Max distance from a face's vertices to a selected surface.
fn face_vertex_dist(mesh: &ConnectedMesh, cls: &Classification, si: usize, fi: usize) -> f64 {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    face.vertex_indices[..vc]
        .iter()
        .map(|&vi| {
            let v = &mesh.vertices[vi];
            point_dist(cls, si, &[v.x, v.y, v.z])
        })
        .fold(0.0_f64, f64::max)
}

/// Build the face → selected-surface map (usize::MAX where unassigned).
fn face_to_surface(mesh: &ConnectedMesh, cls: &Classification) -> Vec<usize> {
    let mut map = vec![usize::MAX; mesh.faces.len()];
    for si in 0..cls.selected.len() {
        for &fi in surface_faces(cls, si) {
            map[fi] = si;
        }
    }
    map
}

/// Sorted, deduplicated vertex indices of a face set.
fn vertices_of(mesh: &ConnectedMesh, faces: &[usize]) -> Vec<usize> {
    let mut v: Vec<usize> = faces
        .iter()
        .flat_map(|&fi| {
            let vc = mesh.faces[fi].vertex_count as usize;
            mesh.faces[fi].vertex_indices[..vc].iter().copied()
        })
        .collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// Refit a surface's geometry to its current face set, in place, keeping the
/// existing fit if the fitter cannot produce one. Unlike `refit_union` this
/// does not require the result to validate: it is used *during* boundary
/// repair, where the support is mid-correction and a validating fit is the
/// goal rather than the precondition.
fn refit_geometry(mesh: &ConnectedMesh, cls: &mut Classification, si: usize) {
    let faces = surface_faces(cls, si).to_vec();
    if faces.is_empty() {
        return;
    }
    let verts = vertices_of(mesh, &faces);
    // The fit comes from the faces that validated strictly; a facet that joined
    // at the boundary bound is a member, not evidence (`Classification::joined`).
    let strict: Vec<usize> = faces.iter().copied().filter(|&fi| !cls.is_joined(fi)).collect();
    let (fit_faces, fit_verts) = if strict.is_empty() || strict.len() == faces.len() {
        (faces.clone(), verts.clone())
    } else {
        let v = vertices_of(mesh, &strict);
        (strict, v)
    };
    match cls.selected[si] {
        SelectedSurface::Planar(i) => {
            if let Some(fit) = super::surface::cached_fit::<PlanarSurface>(mesh, &fit_faces, &fit_verts, 0) {
                cls.planar[i] = PlanarSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        SelectedSurface::Conical(i) => {
            if let Some(mut fit) = super::surface::cached_fit::<ConicalSurface>(mesh, &fit_faces, &fit_verts, 0) {
                ConicalSurface::determine_convexity(mesh, &faces, &mut fit);
                cls.conical[i] = ConicalSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        SelectedSurface::Spherical(i) => {
            if let Some(mut fit) = super::surface::cached_fit::<SphericalSurface>(mesh, &fit_faces, &fit_verts, 0) {
                SphericalSurface::determine_convexity(mesh, &faces, &mut fit);
                cls.spherical[i] = SphericalSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        SelectedSurface::Toroidal(i) => {
            if let Some(mut fit) = super::surface::cached_fit::<ToroidalSurface>(mesh, &fit_faces, &fit_verts, 0) {
                ToroidalSurface::determine_convexity(mesh, &faces, &mut fit);
                cls.toroidal[i] = ToroidalSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        // A B-spline's fit vouches only for its own data and repair never
        // changes its face set.
        SelectedSurface::BSpline(_) => {}
    }
}

// ---------------------------------------------------------------------------
// Boundary snap
// ---------------------------------------------------------------------------

/// Put every shared boundary on the mesh-edge chain the junction geometry
/// picks out, rather than wherever flood and claim resolution happened to
/// leave it.
///
/// For each adjacent pair of committed surfaces, `locus.rs` builds the curve
/// where they meet and turns it into a signed side test; a face on the wrong
/// side, that the other surface also fits, moves. This is a *selection*, not a
/// search: the seam assumption guarantees the correct boundary already exists
/// as a chain of mesh edges, so the only question is which one.
///
/// Assignment and fit are one problem (the same lesson `relabel_boundaries`
/// records), so each round refits every surface its moves touched and rebuilds
/// the loci from the corrected geometry. A face may move at most once per
/// invocation, which bounds the alternation without a round cap: the locus
/// each round is computed from strictly better-supported fits, so a face that
/// wants to move back is evidence the locus was wrong rather than a step to
/// iterate on.
///
/// Pairs the module cannot build a locus for, or whose locus fails validation,
/// are **reported and left alone**. Transversal junctions are not snap's
/// business at all — they are sharp mesh edges, which stage 2.1 already
/// partitions exactly.
fn snap_boundaries(
    mesh: &ConnectedMesh,
    cls: &mut Classification,
    config: &Config,
) -> (usize, usize) {
    let mut already_moved = vec![false; mesh.faces.len()];
    let mut total = 0usize;
    let mut reported: BTreeSet<String> = BTreeSet::new();

    loop {
        let f2s = face_to_surface(mesh, cls);
        let surfs: Vec<Surf<'_>> = cls
            .selected
            .iter()
            .map(|s| Surf::of(s, &cls.planar, &cls.conical, &cls.spherical, &cls.toroidal))
            .collect();
        let pairs = analyze_pairs(mesh, &surfs, &f2s, config.vertex_tolerance_mm);

        let mut round: Vec<(usize, usize, usize)> = Vec::new();
        let mut claimed = vec![false; mesh.faces.len()];
        for p in &pairs {
            if p.is_problem() {
                reported.insert(format!(
                    "  Stage 2.4 snap: surfaces {} ({}) and {} ({}), {} interface faces — {}",
                    p.a, p.kind_a, p.b, p.kind_b, p.interface_faces, p.describe(),
                ));
                continue;
            }
            if let PairOutcome::Snapped { moves, .. } = &p.outcome {
                for &(fi, from, to) in moves {
                    // One move per face per invocation, and one claimant per
                    // face per round: a face at a triple point borders two
                    // junctions, and the first pair in index order takes it.
                    if !already_moved[fi] && !claimed[fi] {
                        claimed[fi] = true;
                        round.push((fi, from, to));
                    }
                }
            }
        }

        if round.is_empty() {
            for line in &reported {
                eprintln!("{line}");
            }
            return (total, reported.len());
        }

        let mut touched: BTreeSet<usize> = BTreeSet::new();
        for &(fi, from, _) in &round {
            surface_faces_mut(cls, from).retain(|&f| f != fi);
            touched.insert(from);
        }
        for &(fi, _, to) in &round {
            surface_faces_mut(cls, to).push(fi);
            touched.insert(to);
            already_moved[fi] = true;
        }
        for si in touched {
            refit_geometry(mesh, cls, si);
        }
        total += round.len();
    }
}

// ---------------------------------------------------------------------------
// Boundary relabel (port of the old refine_boundary_assignments)
// ---------------------------------------------------------------------------

/// Move region-boundary faces to an adjacent surface that fits them at least
/// better, refitting the surfaces each round's moves touched.
///
/// The 2× anti-flapping margin and the 10-round cap this used to carry are
/// gone. Both were the price of a criterion that compares a face against two
/// fits *both* contaminated by that face: it could oscillate, so the margin
/// held it back and the cap stopped it. Refitting inside the loop (added
/// 2026-07-26) removed the contamination, and allowing each face to move at
/// most once bounds the alternation directly, which is a statement about
/// convergence rather than a tuned number. See DEVELOPMENT_PLAN, "Snap shared
/// boundaries to the junction, not to a residual comparison".
fn relabel_boundaries(mesh: &ConnectedMesh, cls: &mut Classification, _config: &Config) -> usize {
    let num_faces = mesh.faces.len();
    let mut total_moved = 0usize;
    // A face moves at most once. This is what makes the loop terminate
    // without a round cap: each round refits from a cleaner support, so a
    // face asking to move back is not converging, it is oscillating between
    // two fits that each still hold some of the other's data.
    let mut already_moved = vec![false; num_faces];

    loop {
        let f2s = face_to_surface(mesh, cls);
        let mut moves: Vec<(usize, usize, usize)> = Vec::new(); // (face, from, to)
        for fi in 0..num_faces {
            let si = f2s[fi];
            // Unassigned faces and B-spline faces never move (a B-spline's
            // fit vouches only for its own data).
            if si == usize::MAX || is_bspline(cls, si) || already_moved[fi] {
                continue;
            }
            let face = &mesh.faces[fi];
            let vc = face.vertex_count as usize;
            let mut candidates: Vec<usize> = face.neighbors[..vc]
                .iter()
                .filter(|&&ni| ni >= 0)
                .map(|&ni| f2s[ni as usize])
                .filter(|&sj| sj != usize::MAX && sj != si && !is_bspline(cls, sj))
                .collect();
            if candidates.is_empty() {
                continue;
            }
            candidates.sort_unstable();
            candidates.dedup();

            let d_own = face_vertex_dist(mesh, cls, si, fi);
            let mut best: Option<(f64, usize)> = None;
            for &sj in &candidates {
                let d = face_vertex_dist(mesh, cls, sj, fi);
                if best.as_ref().is_none_or(|(bd, _)| d < *bd) {
                    best = Some((d, sj));
                }
            }
            if let Some((d_best, sj)) = best {
                if d_best < d_own {
                    moves.push((fi, si, sj));
                }
            }
        }

        if moves.is_empty() {
            break;
        }
        total_moved += moves.len();
        for &(fi, _, _) in &moves {
            already_moved[fi] = true;
        }
        for &(fi, si, _) in &moves {
            surface_faces_mut(cls, si).retain(|&f| f != fi);
        }
        for &(fi, _, sj) in &moves {
            surface_faces_mut(cls, sj).push(fi);
        }

        // Refit every surface a move touched before deciding the next round.
        // Boundary repair and fitting are one problem, not two: a surface
        // holding a row that belongs to its neighbour is pulled off true by
        // exactly that row, and comparing a face against two fits both
        // contaminated by it is what stalls the round. Refitting reopens the
        // gap the 2× margin needs, and each round starts from cleaner
        // geometry — alternating assignment and fit, which is what makes the
        // fixpoint worth iterating to at all.
        let mut touched: Vec<usize> = moves.iter().flat_map(|&(_, si, sj)| [si, sj]).collect();
        touched.sort_unstable();
        touched.dedup();
        for si in touched {
            let mut faces = surface_faces(cls, si).to_vec();
            faces.sort_unstable();
            set_surface_faces(mesh, cls, si, faces);
            refit_geometry(mesh, cls, si);
        }
    }

    if total_moved > 0 {
        rebuild_all(mesh, cls);
    }
    total_moved
}

// ---------------------------------------------------------------------------
// Merge
// ---------------------------------------------------------------------------

/// A refit union hypothesis, tagged by type so the caller can store it back
/// into the matching typed vector.
enum Refit {
    Planar(PlanarHypothesis),
    Conical(ConicalHypothesis),
    Spherical(SphericalHypothesis),
    Toroidal(ToroidalHypothesis),
}

/// Whether two selected surfaces are the same kind for a union refit. The
/// conical family covers cylinder and cone together, as its fitter does.
fn same_kind(a: SelectedSurface, b: SelectedSurface) -> bool {
    matches!(
        (a, b),
        (SelectedSurface::Planar(_), SelectedSurface::Planar(_))
            | (SelectedSurface::Conical(_), SelectedSurface::Conical(_))
            | (SelectedSurface::Spherical(_), SelectedSurface::Spherical(_))
            | (SelectedSurface::Toroidal(_), SelectedSurface::Toroidal(_))
    )
}

/// Merge adjacent surfaces of the same kind whose union validates as ONE
/// surface of that kind.
///
/// `merge_regions` absorbs a small surface whose vertices lie on a larger
/// neighbour's CURRENT fit. That misses pieces of one surface that were each
/// fitted to a partial support: neither fit extrapolates onto the other piece
/// to vertex tolerance, though one surface explains both exactly. Measured on
/// `Open-Task-Light__led_cover_-_top`: a 358 mm² reference cone came out as
/// seven cones of 9-55 facets, every facet within tolerance of the true cone.
///
/// So ask the question directly: refit the union and validate it end to end,
/// exactly as a merge does. Validation at `--vertex-tolerance` is the whole
/// test — two adjacent surfaces that are geometrically different cannot both
/// lie on one fit — so this needs no size limit and no new constant.
/// Smallest pair first, to a fixpoint. B-splines are left alone, as in
/// `merge_regions`.
fn merge_same_surface(
    mesh: &ConnectedMesh,
    cls: &mut Classification,
    centroids: &[[f64; 3]],
    config: &Config,
) -> usize {
    let mut merged = 0usize;
    // Pairs already tried and refused, as (smaller, larger) face counts would
    // not be stable across merges; keyed on the two face sets' first faces and
    // sizes, which change whenever either surface does.
    let mut refused: std::collections::HashSet<(usize, usize, usize, usize)> =
        std::collections::HashSet::new();
    loop {
        let f2s = face_to_surface(mesh, cls);
        let mut pairs: Vec<(usize, usize, usize)> = Vec::new(); // (union size, a, b)
        for si in 0..cls.selected.len() {
            let faces = surface_faces(cls, si);
            if faces.is_empty() || is_bspline(cls, si) {
                continue;
            }
            let mut nbrs: Vec<usize> = Vec::new();
            for &fi in faces {
                let f = &mesh.faces[fi];
                for &ni in &f.neighbors[..f.vertex_count as usize] {
                    if ni >= 0 {
                        let sj = f2s[ni as usize];
                        if sj != usize::MAX && sj > si {
                            nbrs.push(sj);
                        }
                    }
                }
            }
            nbrs.sort_unstable();
            nbrs.dedup();
            for sj in nbrs {
                if !is_bspline(cls, sj) && same_kind(cls.selected[si], cls.selected[sj]) {
                    pairs.push((faces.len() + surface_faces(cls, sj).len(), si, sj));
                }
            }
        }
        pairs.sort_unstable();
        let mut did = false;
        for (_, si, sj) in pairs {
            let (a, b) = (surface_faces(cls, si), surface_faces(cls, sj));
            let key = (a[0], a.len(), b[0], b.len());
            if refused.contains(&key) {
                continue;
            }
            let mut union: Vec<usize> = a.iter().chain(b.iter()).copied().collect();
            union.sort_unstable();
            let verts = vertices_of(mesh, &union);
            // Keep the larger surface's slot; empty the smaller.
            let (keep, drop) = if b.len() >= a.len() { (sj, si) } else { (si, sj) };
            match refit_union(mesh, cls.selected[keep], &union, &verts, centroids, config, &cls.joined) {
                Some(refit) => {
                    store_refit(cls, keep, refit);
                    surface_faces_mut(cls, drop).clear();
                    merged += 1;
                    did = true;
                    break;
                }
                None => {
                    refused.insert(key);
                }
            }
        }
        if !did {
            break;
        }
    }
    merged
}

/// Refit a union of faces with the type of `sel` and validate it end to end
/// (joined facets left out of the fit and held to the boundary bound).
fn refit_union(
    mesh: &ConnectedMesh,
    sel: SelectedSurface,
    faces: &[usize],
    verts: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
    joined: &[bool],
) -> Option<Refit> {
    match sel {
        SelectedSurface::Planar(_) => {
            fit_validate_joined::<PlanarSurface>(mesh, faces, verts, centroids, config, joined)
                .map(Refit::Planar)
        }
        SelectedSurface::Conical(_) => {
            fit_validate_joined::<ConicalSurface>(mesh, faces, verts, centroids, config, joined)
                .map(Refit::Conical)
        }
        SelectedSurface::Spherical(_) => {
            fit_validate_joined::<SphericalSurface>(mesh, faces, verts, centroids, config, joined)
                .map(Refit::Spherical)
        }
        SelectedSurface::Toroidal(_) => {
            fit_validate_joined::<ToroidalSurface>(mesh, faces, verts, centroids, config, joined)
                .map(Refit::Toroidal)
        }
        SelectedSurface::BSpline(_) => None,
    }
}

/// Absorb small regions into adjacent larger real surfaces (fit-validated).
fn merge_regions(
    mesh: &ConnectedMesh,
    cls: &mut Classification,
    centroids: &[[f64; 3]],
    face_region: &[usize],
    config: &Config,
) -> usize {
    let vtol = config.vertex_tolerance_mm;
    let btol = config.boundary_tolerance();
    // Stage-2.2 region of a face; faces of an unknown region (an empty map, as
    // in unit tests) are never across a sharp edge from each other.
    let region = |fi: usize| face_region.get(fi).copied();
    let stol = config.surface_tolerance_mm;
    let mut merged = 0usize;

    loop {
        let f2s = face_to_surface(mesh, cls);
        // Smallest surfaces first (then index) — deterministic.
        let mut order: Vec<usize> = (0..cls.selected.len())
            .filter(|&si| !surface_faces(cls, si).is_empty())
            .collect();
        order.sort_by_key(|&si| (surface_faces(cls, si).len(), si));

        let mut did = false;
        for &si in &order {
            // A B-spline is never a merge source: its fit vouches only for
            // its own data.
            if is_bspline(cls, si) {
                continue;
            }
            let s_faces: Vec<usize> = surface_faces(cls, si).to_vec();
            if s_faces.is_empty() || s_faces.len() > MERGE_MAX_SMALL {
                continue;
            }
            // Adjacent surfaces (real, strictly larger than S).
            let mut neighbors: Vec<usize> = Vec::new();
            for &fi in &s_faces {
                let vc = mesh.faces[fi].vertex_count as usize;
                for &ni in &mesh.faces[fi].neighbors[..vc] {
                    if ni >= 0 {
                        let sj = f2s[ni as usize];
                        if sj != usize::MAX && sj != si {
                            neighbors.push(sj);
                        }
                    }
                }
            }
            neighbors.sort_unstable();
            neighbors.dedup();

            for &ti in &neighbors {
                // A B-spline is never a merge TARGET: a fit vouches only for
                // its data (the 2.3(e) island-dissolution guard is the
                // fit-validated way faces join a freeform patch).
                if is_bspline(cls, ti) {
                    continue;
                }
                let t_len = surface_faces(cls, ti).len();
                if t_len <= 2 || t_len <= s_faces.len() {
                    continue; // T must be a larger, real surface
                }
                // Every S vertex on T within vertex tol, every S centroid
                // within surface tol — "on the same surface" is decisive.
                //
                // With the boundary join on, a facet of S may instead come
                // within `--boundary-tolerance` of T, provided every vertex of
                // it beyond vertex tolerance is on a real edge of S ∪ T — a
                // facet around it is across a sharp edge, or belongs to a third
                // analytic surface of at least two facets — the rule
                // `classify::join_boundary_facets` applies. S is then a
                // sliver of boundary facets that competition fitted on its own
                // before the join could claim them — measured on
                // `VoronUsers__12mm_strip_mini_mount`: cones of 3-8 facets along
                // a tolerant edge of a cylinder — not a surface of its own. Such
                // facets join T: left out of its refit, held to the boundary
                // bound.
                let mut loose: Vec<usize> = Vec::new();
                let verts_ok = s_faces.iter().all(|&fi| {
                    let d = face_vertex_dist(mesh, cls, ti, fi);
                    if d <= vtol {
                        return true;
                    }
                    if !config.boundary_join() || d > btol {
                        return false;
                    }
                    let f = &mesh.faces[fi];
                    let on_edge = f.vertex_indices[..f.vertex_count as usize].iter().all(|&vi| {
                        let v = &mesh.vertices[vi];
                        point_dist(cls, ti, &[v.x, v.y, v.z]) <= vtol
                            || super::vertex_fan(mesh, fi, vi).iter().any(|&g| {
                                let o = f2s[g];
                                (region(g).is_some() && region(g) != region(fi))
                                    || (o != usize::MAX && o != si && o != ti && determined(cls, o))
                            })
                    });
                    if on_edge {
                        loose.push(fi);
                    }
                    on_edge
                });
                if !verts_ok {
                    continue;
                }
                let cents_ok = s_faces.iter().all(|&fi| point_dist(cls, ti, &centroids[fi]) <= stol);
                if !cents_ok {
                    continue;
                }
                // Refit the union and revalidate; keep only if valid.
                let mut union_faces = surface_faces(cls, ti).to_vec();
                union_faces.extend_from_slice(&s_faces);
                union_faces.sort_unstable();
                union_faces.dedup();
                let union_verts = vertices_of(mesh, &union_faces);
                let sel_t = cls.selected[ti];
                let newly: Vec<usize> =
                    loose.iter().copied().filter(|&fi| !cls.is_joined(fi)).collect();
                for &fi in &newly {
                    cls.joined[fi] = true;
                }
                if let Some(refit) =
                    refit_union(mesh, sel_t, &union_faces, &union_verts, centroids, config, &cls.joined)
                {
                    store_refit(cls, ti, refit);
                    surface_faces_mut(cls, si).clear();
                    merged += 1;
                    did = true;
                    break;
                }
                for &fi in &newly {
                    cls.joined[fi] = false;
                }
            }
        }
        if !did {
            break;
        }
    }
    merged
}

/// Store a refit union hypothesis back into the vector for surface `ti`.
fn store_refit(cls: &mut Classification, ti: usize, refit: Refit) {
    match (cls.selected[ti], refit) {
        (SelectedSurface::Planar(i), Refit::Planar(h)) => cls.planar[i] = h,
        (SelectedSurface::Conical(i), Refit::Conical(h)) => cls.conical[i] = h,
        (SelectedSurface::Spherical(i), Refit::Spherical(h)) => cls.spherical[i] = h,
        (SelectedSurface::Toroidal(i), Refit::Toroidal(h)) => cls.toroidal[i] = h,
        _ => unreachable!("refit type must match surface type"),
    }
}

// ---------------------------------------------------------------------------
// Connectivity re-split
// ---------------------------------------------------------------------------

/// Re-split any surface whose faces are not mesh-connected into one surface
/// per connected component (each an independent copy of the fitted geometry).
fn resplit_disconnected(mesh: &ConnectedMesh, cls: &mut Classification) -> usize {
    let mut extra = 0usize;
    let original = cls.selected.len();
    for si in 0..original {
        // A B-spline's faces are the fitted component — mesh-connected by
        // construction — and its geometry is not cloneable; skip it.
        if is_bspline(cls, si) {
            continue;
        }
        let faces = surface_faces(cls, si).to_vec();
        if faces.len() < 2 {
            continue;
        }
        let components = connected_components(mesh, &faces);
        if components.len() < 2 {
            continue;
        }
        // Keep the largest component on the original surface; clone the rest.
        let mut comps = components;
        comps.sort_by_key(|c| std::cmp::Reverse(c.len()));
        for comp in comps.iter().skip(1) {
            clone_surface_with_faces(mesh, cls, si, comp);
            extra += 1;
        }
        set_surface_faces(mesh, cls, si, comps[0].clone());
    }
    extra
}

/// Connected components of a face set over shared mesh edges.
fn connected_components(mesh: &ConnectedMesh, faces: &[usize]) -> Vec<Vec<usize>> {
    let member: std::collections::HashSet<usize> = faces.iter().copied().collect();
    let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut components: Vec<Vec<usize>> = Vec::new();
    for &start in faces {
        if seen.contains(&start) {
            continue;
        }
        let mut comp = Vec::new();
        let mut stack = vec![start];
        seen.insert(start);
        while let Some(fi) = stack.pop() {
            comp.push(fi);
            let vc = mesh.faces[fi].vertex_count as usize;
            for &ni in &mesh.faces[fi].neighbors[..vc] {
                if ni >= 0 {
                    let ni = ni as usize;
                    if member.contains(&ni) && seen.insert(ni) {
                        stack.push(ni);
                    }
                }
            }
        }
        comp.sort_unstable();
        components.push(comp);
    }
    components
}

/// Append a new selected surface cloning `si`'s geometry but owning `faces`.
fn clone_surface_with_faces(
    mesh: &ConnectedMesh,
    cls: &mut Classification,
    si: usize,
    faces: &[usize],
) {
    let verts = vertices_of(mesh, faces);
    match cls.selected[si] {
        SelectedSurface::Planar(i) => {
            let mut h = cls.planar[i].clone();
            h.faces = faces.to_vec();
            h.vertices = verts;
            cls.selected.push(SelectedSurface::Planar(cls.planar.len()));
            cls.planar.push(h);
        }
        SelectedSurface::Conical(i) => {
            let mut h = cls.conical[i].clone();
            h.faces = faces.to_vec();
            h.vertices = verts;
            cls.selected.push(SelectedSurface::Conical(cls.conical.len()));
            cls.conical.push(h);
        }
        SelectedSurface::Spherical(i) => {
            let mut h = cls.spherical[i].clone();
            h.faces = faces.to_vec();
            h.vertices = verts;
            cls.selected.push(SelectedSurface::Spherical(cls.spherical.len()));
            cls.spherical.push(h);
        }
        SelectedSurface::Toroidal(i) => {
            let mut h = cls.toroidal[i].clone();
            h.faces = faces.to_vec();
            h.vertices = verts;
            cls.selected.push(SelectedSurface::Toroidal(cls.toroidal.len()));
            cls.toroidal.push(h);
        }
        // resplit_disconnected skips B-splines (their geometry is not
        // cloneable and a fitted component is already mesh-connected).
        SelectedSurface::BSpline(_) => unreachable!("resplit skips B-spline surfaces"),
    }
}

/// Replace surface `si`'s faces and rebuild its vertex list.
fn set_surface_faces(mesh: &ConnectedMesh, cls: &mut Classification, si: usize, faces: Vec<usize>) {
    let verts = vertices_of(mesh, &faces);
    match cls.selected[si] {
        SelectedSurface::Planar(i) => {
            cls.planar[i].faces = faces;
            cls.planar[i].vertices = verts;
        }
        SelectedSurface::Conical(i) => {
            cls.conical[i].faces = faces;
            cls.conical[i].vertices = verts;
        }
        SelectedSurface::Spherical(i) => {
            cls.spherical[i].faces = faces;
            cls.spherical[i].vertices = verts;
        }
        SelectedSurface::Toroidal(i) => {
            cls.toroidal[i].faces = faces;
            cls.toroidal[i].vertices = verts;
        }
        // rebuild_all re-sorts every surface's own face list; for an inert
        // B-spline this is a no-op re-sort of its fitted faces (repair never
        // gives it a different face set).
        SelectedSurface::BSpline(i) => {
            cls.bspline[i].faces = faces;
            cls.bspline[i].vertices = verts;
        }
    }
}

// ---------------------------------------------------------------------------
// Cleanup + invariants
// ---------------------------------------------------------------------------

/// Sort face lists and rebuild vertex lists for every surface.
fn rebuild_all(mesh: &ConnectedMesh, cls: &mut Classification) {
    for si in 0..cls.selected.len() {
        let mut faces = surface_faces(cls, si).to_vec();
        faces.sort_unstable();
        set_surface_faces(mesh, cls, si, faces);
    }
}

/// Drop selected surfaces that lost all their faces (e.g. merged away).
fn drop_empty_surfaces(cls: &mut Classification) {
    let keep: Vec<bool> = (0..cls.selected.len())
        .map(|si| !surface_faces(cls, si).is_empty())
        .collect();
    let mut si = 0;
    cls.selected.retain(|_| {
        let k = keep[si];
        si += 1;
        k
    });
}

/// Hard partition invariants required before stage 3.
fn assert_invariants(mesh: &ConnectedMesh, cls: &Classification) {
    let n = mesh.faces.len();

    // Total coverage: every face in exactly one selected surface.
    let mut count = vec![0u32; n];
    for si in 0..cls.selected.len() {
        for &fi in surface_faces(cls, si) {
            count[fi] += 1;
        }
    }
    for (fi, &c) in count.iter().enumerate() {
        assert_eq!(c, 1, "stage 2.4 invariant: face {fi} covered {c} times, expected 1");
    }

    // Connectivity: each surface's faces form one mesh-connected component.
    for si in 0..cls.selected.len() {
        let faces = surface_faces(cls, si);
        if faces.len() < 2 {
            continue;
        }
        let comps = connected_components(mesh, faces);
        assert_eq!(
            comps.len(),
            1,
            "stage 2.4 invariant: surface {si} has {} disconnected components",
            comps.len()
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use crate::stage1::{
        ConnectedMesh, MeshFace, MeshVertex,
    };

    fn tri(a: usize, b: usize, c: usize) -> MeshFace {
        MeshFace {
            vertex_count: 3,
            vertex_indices: [a, b, c, 0],
            neighbors: [-1, -1, -1, -1],
            normal: None,
        }
    }

    /// A connected strip of 3 coplanar quads (z=0), each split into two CCW
    /// (+z) triangles — 6 faces, all on the plane z=0.
    fn strip_mesh() -> ConnectedMesh {
        let verts = vec![
            MeshVertex::from_xyz(0.0, 0.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0),
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 1.0, 0.0),
            MeshVertex::from_xyz(2.0, 0.0, 0.0),
            MeshVertex::from_xyz(2.0, 1.0, 0.0),
            MeshVertex::from_xyz(3.0, 0.0, 0.0),
            MeshVertex::from_xyz(3.0, 1.0, 0.0),
        ];
        let faces = vec![
            tri(0, 2, 3), // quad0
            tri(0, 3, 1),
            tri(2, 4, 5), // quad1
            tri(2, 5, 3),
            tri(4, 6, 7), // quad2
            tri(4, 7, 5),
        ];
        let mut mesh = ConnectedMesh {
            vertices: verts,
            faces,
            stats: Default::default(),
        };
        mesh.validate_and_populate_topology().expect("valid mesh");
        mesh
    }

    fn planar(faces: Vec<usize>, verts: Vec<usize>) -> PlanarHypothesis {
        PlanarHypothesis {
            normal: [0.0, 0.0, 1.0],
            distance: 0.0,
            faces,
            vertices: verts,
            error_max: 0.0,
            error_min: 0.0,
            error_abs_sum: 0.0,
        }
    }

    fn test_config() -> config::Config {
        config::parse_config_from(["brepper", "dummy.stl", "-q"]).unwrap()
    }

    /// A strip of 4 quads along x (z=0), 8 triangles, with the two vertices at
    /// x=3 raised `lift` in z: as if the edge between quads 2 and 3 were on an
    /// approximated edge curve, off the plane.
    fn lifted_strip(lift: f64) -> ConnectedMesh {
        let mut verts = Vec::new();
        for x in 0..=4 {
            let z = if x == 3 { lift } else { 0.0 };
            verts.push(MeshVertex::from_xyz(x as f64, 0.0, z));
            verts.push(MeshVertex::from_xyz(x as f64, 1.0, z));
        }
        let mut faces = Vec::new();
        for q in 0..4 {
            let (a, b, c, d) = (2 * q, 2 * q + 1, 2 * q + 2, 2 * q + 3);
            faces.push(tri(a, c, d));
            faces.push(tri(a, d, b));
        }
        let mut mesh = ConnectedMesh { vertices: verts, faces, stats: Default::default() };
        mesh.validate_and_populate_topology().expect("valid mesh");
        mesh
    }

    /// T = quads 0-1, S = quad 2, U = quad 3, all fitted as z=0.
    fn lifted_cls(n_faces: usize) -> Classification {
        Classification {
            planar: vec![
                planar(vec![0, 1, 2, 3], vec![0, 1, 2, 3, 4, 5]),
                planar(vec![4, 5], vec![4, 5, 6, 7]),
                planar(vec![6, 7], vec![6, 7, 8, 9]),
            ],
            conical: Vec::new(),
            spherical: Vec::new(),
            toroidal: Vec::new(),
            bspline: Vec::new(),
            selected: vec![
                SelectedSurface::Planar(0),
                SelectedSurface::Planar(1),
                SelectedSurface::Planar(2),
            ],
            joined: vec![false; n_faces],
        }
    }

    #[test]
    fn boundary_merge_joins_a_sliver_whose_off_vertices_stay_on_an_edge() {
        // S's lifted vertices (5e-4 off T, within the 1e-3 boundary bound) are
        // shared with U, so they stay on an edge of S+T: S merges into T, and
        // its facets are marked joined. U's would then be inside T+U — every
        // facet around them on T or U — so U is refused.
        let mesh = lifted_strip(5e-4);
        let mut cls = lifted_cls(mesh.faces.len());
        let centroids: Vec<[f64; 3]> =
            mesh.faces.iter().map(|f| face_centroid(f, &mesh.vertices)).collect();
        let merged = merge_regions(&mesh, &mut cls, &centroids, &[], &test_config());
        assert_eq!(merged, 1);
        assert_eq!(cls.planar[0].faces.len(), 6, "S joined T");
        assert!(cls.planar[1].faces.is_empty());
        assert_eq!(cls.planar[2].faces, vec![6, 7], "U kept: its lifted vertices would be inside");
        assert!(cls.is_joined(4) && cls.is_joined(5));
        assert!(!cls.is_joined(6) && !cls.is_joined(7));
        // T's fit came from its strict facets only: still z=0.
        assert!(cls.planar[0].distance.abs() < 1e-12);
    }

    #[test]
    fn boundary_merge_is_off_at_the_vertex_tolerance_and_past_the_bound() {
        let off = config::parse_config_from([
            "brepper", "dummy.stl", "-q", "--boundary-tolerance=5e-5",
        ])
        .unwrap();
        let mesh = lifted_strip(5e-4);
        let mut cls = lifted_cls(mesh.faces.len());
        let centroids: Vec<[f64; 3]> =
            mesh.faces.iter().map(|f| face_centroid(f, &mesh.vertices)).collect();
        assert_eq!(merge_regions(&mesh, &mut cls, &centroids, &[], &off), 0);

        // 2e-3 off is past the default 1e-3 bound.
        let mesh = lifted_strip(2e-3);
        let mut cls = lifted_cls(mesh.faces.len());
        let centroids: Vec<[f64; 3]> =
            mesh.faces.iter().map(|f| face_centroid(f, &mesh.vertices)).collect();
        assert_eq!(merge_regions(&mesh, &mut cls, &centroids, &[], &test_config()), 0);
    }

    #[test]
    fn merge_absorbs_coplanar_small_region() {
        let mesh = strip_mesh();
        // T = quads 0+1 (4 faces), S = quad 2 (2 faces), all coplanar.
        let mut cls = Classification {
            planar: vec![
                planar(vec![0, 1, 2, 3], vec![0, 1, 2, 3, 4, 5]),
                planar(vec![4, 5], vec![4, 5, 6, 7]),
            ],
            conical: Vec::new(),
            spherical: Vec::new(),
            toroidal: Vec::new(),
            bspline: Vec::new(),
            selected: vec![SelectedSurface::Planar(0), SelectedSurface::Planar(1)],
            joined: Vec::new(),
        };
        let stats = repair(&mesh, &mut cls, &[], &test_config());
        assert_eq!(stats.merged, 1, "the coplanar small region should merge");
        assert_eq!(cls.selected.len(), 1, "one surface remains after merge");
        let faces = surface_faces(&cls, 0);
        assert_eq!(faces.len(), 6, "merged surface covers all 6 faces");
    }

    /// Two horizontal coplanar quads (z=0) plus one vertical quad (x=2)
    /// sharing their far edge — an L. Faces 0..3 lie on z=0, faces 4..5 on
    /// x=2, all one connected mesh.
    fn l_mesh() -> ConnectedMesh {
        let verts = vec![
            MeshVertex::from_xyz(0.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 1.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0),
            MeshVertex::from_xyz(2.0, 0.0, 0.0),
            MeshVertex::from_xyz(2.0, 1.0, 0.0),
            MeshVertex::from_xyz(2.0, 0.0, 1.0),
            MeshVertex::from_xyz(2.0, 1.0, 1.0),
        ];
        let faces = vec![
            tri(0, 1, 2), // quad0, z=0
            tri(0, 2, 3),
            tri(1, 4, 5), // quad1, z=0
            tri(1, 5, 2),
            tri(5, 4, 6), // quad2, x=2 (vertical); wound so edge 4-5 is
            tri(5, 6, 7), // traversed opposite to quad1 for consistent orientation
        ];
        let mut mesh = ConnectedMesh {
            vertices: verts,
            faces,
            stats: Default::default(),
        };
        mesh.validate_and_populate_topology().expect("valid mesh");
        mesh
    }

    #[test]
    fn merge_leaves_perpendicular_region_alone() {
        let mesh = l_mesh();
        // T = the horizontal plane (4 faces), S = the vertical plane (2
        // faces). S shares an edge with T but its top vertices sit 1 mm off
        // T's plane, so the vertex-tolerance check must reject the merge.
        let mut cls = Classification {
            planar: vec![
                planar(vec![0, 1, 2, 3], vec![0, 1, 2, 3, 4, 5]),
                PlanarHypothesis {
                    normal: [1.0, 0.0, 0.0],
                    distance: 2.0,
                    faces: vec![4, 5],
                    vertices: vec![4, 5, 6, 7],
                    error_max: 0.0,
                    error_min: 0.0,
                    error_abs_sum: 0.0,
                },
            ],
            conical: Vec::new(),
            spherical: Vec::new(),
            toroidal: Vec::new(),
            bspline: Vec::new(),
            selected: vec![SelectedSurface::Planar(0), SelectedSurface::Planar(1)],
            joined: Vec::new(),
        };
        let stats = repair(&mesh, &mut cls, &[], &test_config());
        assert_eq!(stats.merged, 0, "perpendicular planes must not merge");
        assert_eq!(cls.selected.len(), 2, "both surfaces survive");
        assert_invariants(&mesh, &cls);
    }

    #[test]
    fn resplit_breaks_disconnected_surface() {
        let mesh = strip_mesh();
        // One "surface" owning quad0 (faces 0,1) and quad2 (faces 4,5), which
        // are NOT mesh-adjacent (quad1 sits between them). Re-split must break
        // it into two connected surfaces.
        let mut cls = Classification {
            planar: vec![
                planar(vec![0, 1, 4, 5], vec![0, 1, 2, 3, 4, 5, 6, 7]),
                planar(vec![2, 3], vec![2, 3, 4, 5]),
            ],
            conical: Vec::new(),
            spherical: Vec::new(),
            toroidal: Vec::new(),
            bspline: Vec::new(),
            selected: vec![SelectedSurface::Planar(0), SelectedSurface::Planar(1)],
            joined: Vec::new(),
        };
        let extra = resplit_disconnected(&mesh, &mut cls);
        assert_eq!(extra, 1, "the disconnected surface splits into two");
        // After the split every surface is one connected component.
        for si in 0..cls.selected.len() {
            let faces = surface_faces(&cls, si);
            if faces.len() >= 2 {
                assert_eq!(connected_components(&mesh, faces).len(), 1);
            }
        }
    }
}
