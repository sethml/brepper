//! Stage 2.3(b): scoped fit competition.
//!
//! For a smooth region the whole-region ladder (2.3(a)) cannot explain — a
//! tangent chain like the rounded-cube shell (6 planes, 12 fillet cylinders,
//! 8 corner spheres) or the nosecone's cone+sphere — decompose it by
//! seed-first multi-type competition run INSIDE the region. The sharp
//! boundary cordons the arena, so competition cost is bounded and
//! cross-region accidents are impossible.
//!
//! - **Seeds**: region faces not yet claimed by any committed hypothesis, in
//!   jet-interior-first order (faces whose (κ1,κ2) best match their
//!   neighbors' seed before junction rows). Claimed faces never seed.
//! - **Per seed**: grow each candidate surface type (jet-pruned) from a
//!   small jet-similar seed patch — fit, flood-accept every region face the
//!   fit validates, refit on the grown support, to a fixpoint — and commit
//!   the hypothesis with the **largest validated area** (simplest-DOF
//!   tiebreak). The jet only orders and prunes; vertex-tolerance validation
//!   decides every face.
//! - **Claims**: a flood may grow through faces already claimed by an
//!   earlier hypothesis (junction rows that genuinely lie on two surfaces),
//!   so a face can be claimed by several hypotheses. `resolve_claims`
//!   assigns each contested face to the claimant that fits its vertices best
//!   (2× margin to leave the larger claimant).
//!
//! **Soundness** (the sub-region lemma, under exactness): an accidental fit
//! is supported on a *sub-region*; the true surface fitted from an interior
//! seed floods — vertex tolerance the arbiter — to its entire true extent,
//! and largest-validated-area then selects it over the accidental one.
//!
//! Faces no hypothesis claims are returned as leftovers; the caller degrades
//! them to per-facet planes (freeform fitting is stage 2.3(e), a later step).

use std::collections::HashSet;

use crate::config::Config;
use crate::stage1::ConnectedMesh;

use super::bfs::SurfaceType;
use super::conical::ConicalSurface;
use super::edges::JetField;
use super::planar::PlanarSurface;
use super::regions::Region;
use super::spherical::SphericalSurface;
use super::{
    face_area, ConicalHypothesis, PlanarHypothesis, SphericalHypothesis, MIN_CYLINDER_FACES,
    MIN_SPHERE_FACES,
};

/// Max faces in the jet-similar seed patch used to condition the initial fit.
const SEED_PATCH_CAP: usize = 32;

/// Max BFS depth of the seed patch. Kept shallow so the patch stays local to
/// the seed and does not bleed across a tangent junction between two surfaces
/// of the same curvature class (e.g. a cone meeting a cylinder, where the
/// curvature is continuous across the seam) — a straddling patch fits neither
/// surface. An interior seed's shallow neighbourhood is on one surface, which
/// is all the initial fit needs; the flood then grows it to the true extent.
const MAX_PATCH_DEPTH: usize = 3;

/// Refit-and-reflood rounds while growing one hypothesis.
const MAX_GROW_ROUNDS: usize = 8;

/// Flood-acceptance tolerance as a multiple of `--vertex-tolerance`. A fit
/// grown from partial support has ~1e-4 mm residual and a coarsely-tessellated
/// "flat" face reads a small nonzero curvature (plane residual ~1e-3 mm), both
/// well above the 5e-5 mm default vertex tolerance yet far below the ~0.1 mm a
/// neighbouring surface diverges by one row past a tangent junction — so a
/// modest loosening bootstraps the grow without leaking across junctions.
/// Any residual one-row leak is corrected by claim resolution, and every
/// committed hypothesis is refit on its full support (accurate to ~1e-6 mm).
const GROW_VTOL_FACTOR: f64 = 30.0;

/// A committed analytic hypothesis from competition, tagged by type.
pub(super) enum CompeteHyp {
    Planar(PlanarHypothesis),
    Conical(ConicalHypothesis),
    Spherical(SphericalHypothesis),
}

impl CompeteHyp {
    fn faces(&self) -> &[usize] {
        match self {
            CompeteHyp::Planar(h) => &h.faces,
            CompeteHyp::Conical(h) => &h.faces,
            CompeteHyp::Spherical(h) => &h.faces,
        }
    }
    fn set_faces(&mut self, faces: Vec<usize>, verts: Vec<usize>) {
        match self {
            CompeteHyp::Planar(h) => {
                h.faces = faces;
                h.vertices = verts;
            }
            CompeteHyp::Conical(h) => {
                h.faces = faces;
                h.vertices = verts;
            }
            CompeteHyp::Spherical(h) => {
                h.faces = faces;
                h.vertices = verts;
            }
        }
    }
}

/// Decompose a tangent-chain region. Returns the committed analytic
/// hypotheses (a strict partition after claim resolution) plus the leftover
/// faces no hypothesis explained.
pub(super) fn compete_region(
    mesh: &ConnectedMesh,
    region: &Region,
    jet: &[JetField],
    centroids: &[[f64; 3]],
    config: &Config,
) -> (Vec<CompeteHyp>, Vec<usize>) {
    let n = mesh.faces.len();
    let debug = config.verbosity >= 3 || std::env::var("BREPPER_DEBUG_COMPETE").is_ok();
    let mut in_region = vec![false; n];
    for &fi in &region.faces {
        in_region[fi] = true;
    }

    // Region curvature scale: a face reads "flat" when its principal
    // curvature radius exceeds a couple of region diameters.
    let extent = region_extent(mesh, region);
    let flat_eps = if extent > 0.0 { 1.0 / (2.0 * extent) } else { 1e-6 };

    // Seed order: jet-interior-first (ascending curvature mismatch with
    // region neighbors), ties by face index for determinism.
    let mut seeds = region.faces.clone();
    let score: Vec<f64> = (0..n)
        .map(|fi| {
            if !in_region[fi] {
                return f64::INFINITY;
            }
            neighbor_jet_mismatch(mesh, fi, jet, &in_region)
        })
        .collect();
    seeds.sort_by(|&a, &b| {
        score[a]
            .partial_cmp(&score[b])
            .unwrap()
            .then_with(|| a.cmp(&b))
    });

    let mut claimed = vec![false; n]; // out of the seed pool once claimed
    let mut committed: Vec<CompeteHyp> = Vec::new();
    let mut claim_lists: Vec<Vec<usize>> = vec![Vec::new(); n]; // face -> committed idxs

    for &seed in &seeds {
        if claimed[seed] {
            continue;
        }
        let jf = jet[seed];
        let a1 = jf.k1.abs();
        let a2 = jf.k2.abs();

        // Grow the jet-pruned candidate types; keep the largest validated
        // area (simplest-DOF tiebreak).
        let mut best: Option<(CompeteHyp, f64, u8)> = None; // (hyp, area, dof)
        let consider = |hyp_area: Option<(CompeteHyp, f64)>, dof: u8, best: &mut Option<(CompeteHyp, f64, u8)>| {
            if let Some((h, area)) = hyp_area {
                let better = match best {
                    None => true,
                    Some((_, ba, bd)) => area > *ba + 1e-9 || ((area - *ba).abs() <= 1e-9 && dof < *bd),
                };
                if better {
                    *best = Some((h, area, dof));
                }
            }
        };

        // Plane is always a candidate (cheap; catches misread flats).
        let plane = grow::<PlanarSurface>(mesh, &in_region, seed, jet, centroids, config, 1);
        let mut cone_area = None;
        let mut sphere_area = None;
        if a1 >= flat_eps {
            // Curved: a developable (cylinder/cone) is always plausible.
            let cone = grow::<ConicalSurface>(mesh, &in_region, seed, jet, centroids, config, MIN_CYLINDER_FACES);
            cone_area = cone.as_ref().map(|(_, a, f)| (*a, f.len()));
            consider(cone.map(|(h, a, _)| (CompeteHyp::Conical(h), a)), 2, &mut best);
            // Curved in both directions → also try a sphere.
            if a2 >= 0.3 * a1 {
                let sphere = grow::<SphericalSurface>(mesh, &in_region, seed, jet, centroids, config, MIN_SPHERE_FACES);
                sphere_area = sphere.as_ref().map(|(_, a, f)| (*a, f.len()));
                consider(sphere.map(|(h, a, _)| (CompeteHyp::Spherical(h), a)), 1, &mut best);
            }
        }
        let plane_area = plane.as_ref().map(|(_, a, f)| (*a, f.len()));
        consider(plane.map(|(h, a, _)| (CompeteHyp::Planar(h), a)), 0, &mut best);

        if debug && best.is_none() {
            eprintln!(
                "    compete: seed {seed} jet=({:.3},{:.3}) class {:?} → NO GROWTH (plane {plane_area:?} cone {cone_area:?} sphere {sphere_area:?})",
                jf.k1, jf.k2, curvature_class(jf),
            );
        }

        match best {
            Some((hyp, area, dof)) => {
                if debug {
                    let kind = match hyp {
                        CompeteHyp::Planar(_) => "plane",
                        CompeteHyp::Conical(ref h) => if h.is_cylinder() { "cyl" } else { "cone" },
                        CompeteHyp::Spherical(_) => "sphere",
                    };
                    eprintln!(
                        "    compete: seed {seed} → {kind} ({} faces, area {area:.2}, dof {dof})",
                        hyp.faces().len()
                    );
                }
                let idx = committed.len();
                for &f in hyp.faces() {
                    claimed[f] = true;
                    claim_lists[f].push(idx);
                }
                committed.push(hyp);
            }
            None => {
                // The seed grows nothing analytic; drop it from the pool so
                // the loop terminates. It becomes a leftover.
                claimed[seed] = true;
            }
        }
    }

    resolve_claims(mesh, &region.faces, &mut committed, &claim_lists)
}

/// Resolve multi-claimed faces to their best-fitting claimant (2× margin to
/// leave the larger claimant), rebuild each hypothesis's face/vertex lists,
/// drop hypotheses that lost all faces, and collect unclaimed leftovers.
fn resolve_claims(
    mesh: &ConnectedMesh,
    region_faces: &[usize],
    committed: &mut Vec<CompeteHyp>,
    claim_lists: &[Vec<usize>],
) -> (Vec<CompeteHyp>, Vec<usize>) {
    let hyp_area: Vec<f64> = committed
        .iter()
        .map(|h| h.faces().iter().map(|&fi| face_area(&mesh.faces[fi], &mesh.vertices)).sum())
        .collect();

    let mut assigned: Vec<Vec<usize>> = vec![Vec::new(); committed.len()];
    let mut leftovers: Vec<usize> = Vec::new();

    for &fi in region_faces {
        let claimants = &claim_lists[fi];
        if claimants.is_empty() {
            leftovers.push(fi);
            continue;
        }
        if claimants.len() == 1 {
            assigned[claimants[0]].push(fi);
            continue;
        }
        // Best (smallest max-vertex-distance) claimant; the larger-area
        // claimant wins ties and only loses by a 2× margin (a true junction
        // face fits both sides equally — the margin prevents flapping).
        let dist = |ci: usize| max_vertex_dist(mesh, &committed[ci], fi);
        let mut order: Vec<usize> = claimants.clone();
        order.sort_by(|&a, &b| {
            hyp_area[b]
                .partial_cmp(&hyp_area[a])
                .unwrap()
                .then_with(|| a.cmp(&b))
        });
        let incumbent = order[0]; // largest area
        let mut winner = incumbent;
        let mut winner_d = dist(incumbent);
        for &ci in &order[1..] {
            let d = dist(ci);
            if d * 2.0 < winner_d {
                winner = ci;
                winner_d = d;
            }
        }
        assigned[winner].push(fi);
    }

    let mut out = Vec::new();
    for (ci, mut hyp) in committed.drain(..).enumerate() {
        let faces = std::mem::take(&mut assigned[ci]);
        if faces.is_empty() {
            continue;
        }
        let verts = vertices_of(mesh, &faces);
        hyp.set_faces(faces, verts);
        out.push(hyp);
    }
    (out, leftovers)
}

/// Grow a hypothesis of type `S` from `seed`: condition an initial fit on a
/// jet-similar seed patch, then flood-accept every region face the fit
/// validates and refit on the grown support, to a fixpoint. Returns the
/// hypothesis, its validated area, and its face list.
fn grow<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    seed: usize,
    jet: &[JetField],
    centroids: &[[f64; 3]],
    config: &Config,
    min_faces: usize,
) -> Option<(S::Hypothesis, f64, Vec<usize>)> {
    let vtol = config.vertex_tolerance_mm * GROW_VTOL_FACTOR;
    let stol = config.surface_tolerance_mm;

    let mut accepted = seed_patch(mesh, in_region, seed, jet);
    let mut prev = 0usize;
    for round in 0..MAX_GROW_ROUNDS {
        let verts = vertices_of(mesh, &accepted);
        let mut fit = S::fit_surface(mesh, &accepted, &verts, config.verbosity)?;
        S::determine_convexity(mesh, &accepted, &mut fit);
        let comp = flood::<S>(mesh, in_region, seed, &fit, centroids, vtol, stol);
        if comp.is_empty() {
            return None; // the seed does not lie on its own fit — reject
        }
        if comp.len() == accepted.len() && round > 0 {
            accepted = comp;
            break;
        }
        if comp.len() == prev {
            // Oscillating between two sizes — stop at the larger.
            if comp.len() >= accepted.len() {
                accepted = comp;
            }
            break;
        }
        prev = accepted.len();
        accepted = comp;
    }

    // Final fit + validated flood.
    let verts = vertices_of(mesh, &accepted);
    let mut fit = S::fit_surface(mesh, &accepted, &verts, config.verbosity)?;
    S::determine_convexity(mesh, &accepted, &mut fit);
    let faces = flood::<S>(mesh, in_region, seed, &fit, centroids, vtol, stol);
    if faces.len() < min_faces {
        return None;
    }
    let area: f64 = faces.iter().map(|&fi| face_area(&mesh.faces[fi], &mesh.vertices)).sum();
    let verts = vertices_of(mesh, &faces);
    let hyp = S::make_hypothesis(mesh, faces.clone(), verts, &fit);
    Some((hyp, area, faces))
}

/// Connected region-faces reachable from `seed` that all validate `fit`
/// (vertex tolerance + centroid tolerance + orientation coherence).
fn flood<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    seed: usize,
    fit: &S::Fit,
    centroids: &[[f64; 3]],
    vtol: f64,
    stol: f64,
) -> Vec<usize> {
    if !face_validates::<S>(mesh, seed, fit, centroids, vtol, stol) {
        return Vec::new();
    }
    let mut seen: HashSet<usize> = HashSet::new();
    seen.insert(seed);
    let mut out = Vec::new();
    let mut stack = vec![seed];
    while let Some(fi) = stack.pop() {
        out.push(fi);
        let vc = mesh.faces[fi].vertex_count as usize;
        for &ni in &mesh.faces[fi].neighbors[..vc] {
            if ni < 0 {
                continue;
            }
            let ni = ni as usize;
            if in_region[ni]
                && !seen.contains(&ni)
                && face_validates::<S>(mesh, ni, fit, centroids, vtol, stol)
            {
                seen.insert(ni);
                stack.push(ni);
            }
        }
    }
    out.sort_unstable();
    out
}

/// Does face `fi` validate `fit` (all vertices ≤ vtol, centroid ≤ stol,
/// orientation coherent)?
fn face_validates<S: SurfaceType>(
    mesh: &ConnectedMesh,
    fi: usize,
    fit: &S::Fit,
    centroids: &[[f64; 3]],
    vtol: f64,
    stol: f64,
) -> bool {
    if !S::convexity_check(fit, mesh, fi) {
        return false;
    }
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    for &vi in &face.vertex_indices[..vc] {
        let v = &mesh.vertices[vi];
        if S::vertex_distance(fit, &[v.x, v.y, v.z]) > vtol {
            return false;
        }
    }
    S::centroid_distance(fit, &centroids[fi]) <= stol
}

/// A small jet-similar patch around `seed` (same curvature class, magnitude
/// within 2×), capped, for a well-conditioned initial fit.
fn seed_patch(mesh: &ConnectedMesh, in_region: &[bool], seed: usize, jet: &[JetField]) -> Vec<usize> {
    let seed_class = curvature_class(jet[seed]);
    let seed_mag = jet[seed].k1.abs().max(jet[seed].k2.abs());
    let mut seen: HashSet<usize> = HashSet::new();
    seen.insert(seed);
    let mut out = vec![seed];
    // Depth-limited BFS (frontier by ring).
    let mut frontier = vec![seed];
    for _depth in 0..MAX_PATCH_DEPTH {
        if out.len() >= SEED_PATCH_CAP {
            break;
        }
        let mut next = Vec::new();
        for &fi in &frontier {
            let vc = mesh.faces[fi].vertex_count as usize;
            for &ni in &mesh.faces[fi].neighbors[..vc] {
                if ni < 0 {
                    continue;
                }
                let ni = ni as usize;
                if !in_region[ni] || seen.contains(&ni) {
                    continue;
                }
                if curvature_class(jet[ni]) != seed_class {
                    continue;
                }
                let mag = jet[ni].k1.abs().max(jet[ni].k2.abs());
                if seed_mag > 1e-9 && (mag > 1.5 * seed_mag || mag < 0.667 * seed_mag) {
                    continue;
                }
                seen.insert(ni);
                out.push(ni);
                next.push(ni);
                if out.len() >= SEED_PATCH_CAP {
                    break;
                }
            }
        }
        frontier = next;
    }
    out.sort_unstable();
    out
}

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
enum CurvClass {
    Flat,
    Developable,
    Spherelike,
}

fn curvature_class(jf: JetField) -> CurvClass {
    let a1 = jf.k1.abs();
    let a2 = jf.k2.abs();
    // Absolute floor is generous; only used to cluster the seed patch.
    if a1 < 1e-4 {
        CurvClass::Flat
    } else if a2 >= 0.4 * a1 {
        CurvClass::Spherelike
    } else {
        CurvClass::Developable
    }
}

/// Curvature mismatch of a face vs its region neighbors (low = interior).
fn neighbor_jet_mismatch(mesh: &ConnectedMesh, fi: usize, jet: &[JetField], in_region: &[bool]) -> f64 {
    let jf = jet[fi];
    let vc = mesh.faces[fi].vertex_count as usize;
    let mut worst = 0.0_f64;
    for &ni in &mesh.faces[fi].neighbors[..vc] {
        if ni < 0 {
            continue;
        }
        let ni = ni as usize;
        if !in_region[ni] {
            continue;
        }
        let jn = jet[ni];
        let d = (jf.k1 - jn.k1).abs() + (jf.k2 - jn.k2).abs();
        worst = worst.max(d);
    }
    worst
}

/// Max distance from face `fi`'s vertices to a committed hypothesis.
fn max_vertex_dist(mesh: &ConnectedMesh, hyp: &CompeteHyp, fi: usize) -> f64 {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    face.vertex_indices[..vc]
        .iter()
        .map(|&vi| {
            let v = &mesh.vertices[vi];
            let p = [v.x, v.y, v.z];
            match hyp {
                CompeteHyp::Planar(h) => {
                    (h.normal[0] * p[0] + h.normal[1] * p[1] + h.normal[2] * p[2] - h.distance).abs()
                }
                CompeteHyp::Conical(h) => {
                    if h.is_cylinder() {
                        super::cylindrical::vertex_to_cylinder_distance(
                            v, &h.apex, &h.axis_direction, h.radius_at_apex,
                        )
                        .abs()
                    } else {
                        super::conical::vertex_to_cone_distance(v, &h.apex, &h.axis_direction, h.half_angle)
                            .abs()
                    }
                }
                CompeteHyp::Spherical(h) => {
                    super::spherical::vertex_to_sphere_distance(v, &h.center, h.radius).abs()
                }
            }
        })
        .fold(0.0_f64, f64::max)
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

/// Largest bounding-box dimension of a region's vertices.
fn region_extent(mesh: &ConnectedMesh, region: &Region) -> f64 {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for &vi in &region.vertices {
        let v = &mesh.vertices[vi];
        for (k, c) in [v.x, v.y, v.z].into_iter().enumerate() {
            min[k] = min[k].min(c);
            max[k] = max[k].max(c);
        }
    }
    (max[0] - min[0]).max(max[1] - min[1]).max(max[2] - min[2])
}
