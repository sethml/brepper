//! Stage 2.4: resolution & repair.
//!
//! Turns stage 2.3's per-region hypotheses into a strict, connected,
//! validated partition. All operations are fit-validated and deterministic
//! (fixed iteration order — by surface index, smallest first for merges —
//! bounded rounds).
//!
//! - **Boundary relabel** (±1-row polish): move a region-boundary face to an
//!   adjacent surface that fits its vertices at least 2× better, refitting
//!   every surface a move touched before the next round — assignment and fit
//!   are one problem, since a surface holding its neighbour's row is pulled
//!   off true by exactly that row. Iterates to a fixpoint (≤ 10 rounds;
//!   boundaries only move locally).
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

use super::classify::{fit_validate, Classification};
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
    pub relabeled: usize,
    pub merged: usize,
    pub split: usize,
}

/// Run the stage-2.4 repair sequence in place, then assert the hard
/// partition invariants.
pub(super) fn repair(mesh: &ConnectedMesh, cls: &mut Classification, config: &Config) -> RepairStats {
    let centroids: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|f| face_centroid(f, &mesh.vertices))
        .collect();

    let stats = RepairStats {
        relabeled: relabel_boundaries(mesh, cls, config),
        merged: merge_regions(mesh, cls, &centroids, config),
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
    match cls.selected[si] {
        SelectedSurface::Planar(i) => {
            if let Some(fit) = PlanarSurface::fit_surface(mesh, &faces, &verts, 0) {
                cls.planar[i] = PlanarSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        SelectedSurface::Conical(i) => {
            if let Some(mut fit) = ConicalSurface::fit_surface(mesh, &faces, &verts, 0) {
                ConicalSurface::determine_convexity(mesh, &faces, &mut fit);
                cls.conical[i] = ConicalSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        SelectedSurface::Spherical(i) => {
            if let Some(mut fit) = SphericalSurface::fit_surface(mesh, &faces, &verts, 0) {
                SphericalSurface::determine_convexity(mesh, &faces, &mut fit);
                cls.spherical[i] = SphericalSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        SelectedSurface::Toroidal(i) => {
            if let Some(mut fit) = ToroidalSurface::fit_surface(mesh, &faces, &verts, 0) {
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
// Boundary relabel (port of the old refine_boundary_assignments)
// ---------------------------------------------------------------------------

/// Move region-boundary faces to an adjacent surface that fits them at least
/// 2× better, refitting the surfaces each round's moves touched.
fn relabel_boundaries(mesh: &ConnectedMesh, cls: &mut Classification, _config: &Config) -> usize {
    let num_faces = mesh.faces.len();
    let mut total_moved = 0usize;

    for _iter in 0..10 {
        let f2s = face_to_surface(mesh, cls);
        let mut moves: Vec<(usize, usize, usize)> = Vec::new(); // (face, from, to)
        for fi in 0..num_faces {
            let si = f2s[fi];
            // Unassigned faces and B-spline faces never move (a B-spline's
            // fit vouches only for its own data).
            if si == usize::MAX || is_bspline(cls, si) {
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
                if d_best < d_own * 0.5 {
                    moves.push((fi, si, sj));
                }
            }
        }

        if moves.is_empty() {
            break;
        }
        total_moved += moves.len();
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

/// Refit a union of faces with the type of `sel` and validate it end to end.
fn refit_union(
    mesh: &ConnectedMesh,
    sel: SelectedSurface,
    faces: &[usize],
    verts: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
) -> Option<Refit> {
    match sel {
        SelectedSurface::Planar(_) => {
            fit_validate::<PlanarSurface>(mesh, faces, verts, centroids, config).map(Refit::Planar)
        }
        SelectedSurface::Conical(_) => {
            fit_validate::<ConicalSurface>(mesh, faces, verts, centroids, config).map(Refit::Conical)
        }
        SelectedSurface::Spherical(_) => {
            fit_validate::<SphericalSurface>(mesh, faces, verts, centroids, config)
                .map(Refit::Spherical)
        }
        SelectedSurface::Toroidal(_) => {
            fit_validate::<ToroidalSurface>(mesh, faces, verts, centroids, config)
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
    config: &Config,
) -> usize {
    let vtol = config.vertex_tolerance_mm;
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

            let s_verts = vertices_of(mesh, &s_faces);
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
                let verts_ok = s_verts.iter().all(|&vi| {
                    let v = &mesh.vertices[vi];
                    point_dist(cls, ti, &[v.x, v.y, v.z]) <= vtol
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
                if let Some(refit) =
                    refit_union(mesh, sel_t, &union_faces, &union_verts, centroids, config)
                {
                    store_refit(cls, ti, refit);
                    surface_faces_mut(cls, si).clear();
                    merged += 1;
                    did = true;
                    break;
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
        };
        let stats = repair(&mesh, &mut cls, &test_config());
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
        };
        let stats = repair(&mesh, &mut cls, &test_config());
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
