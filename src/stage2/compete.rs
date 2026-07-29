//! Stage 2.3(b): scoped fit competition.
//!
//! For a smooth region the whole-region ladder (2.3(a)) cannot explain — a
//! tangent chain like the rounded-cube shell (6 planes, 12 fillet cylinders,
//! 8 corner spheres) or the nosecone's cone+sphere — decompose it by
//! seed-first multi-type competition run INSIDE the region. The sharp
//! boundary cordons the arena, so competition cost is bounded and
//! cross-region accidents are impossible.
//!
//! - **Seeds**: every region face seeds exactly once, in jet-interior-first
//!   order (faces whose (κ1,κ2) best match their neighbors' seed before
//!   junction rows). Committing a hypothesis does NOT retire the faces it grew
//!   into — growing into a face is a hypothesis, not a conclusion, and a face
//!   swallowed by a wrong surface must still be able to seed the right one
//!   (`resolve_claims` arbitrates the resulting multi-claims). Only faces in
//!   the *interior* of the explained set are skipped, since a fit grown from
//!   the middle of a validated surface merely refits that surface.
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
//! - **Toroidal candidates (stage 2.3(c))**: a torus competes like any other
//!   type — grown directly from curved seeds (fillet tori, bend tubes), and
//!   additionally *promoted* from a cylinder whose flood stalled against
//!   faces whose jet continues the tube (κ1 ≈ the tube's 1/r, |κ2| well
//!   below κ1 — a bent tube's signature). The promotion seeds its fit from
//!   the pure-bend patch just past the stall line: the cylinder fit itself
//!   locates the tangency, so the patch cannot straddle it the way a
//!   jet-class patch can (bend and straight tube share κ1).
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

use super::surface::SurfaceType;
use super::conical::{ConeFit, ConicalSurface};
use super::edges::JetField;
use super::planar::PlanarSurface;
use super::regions::Region;
use super::spherical::{bounding_box_diagonal, SphericalSurface};
use super::toroidal::ToroidalSurface;
use super::{
    face_area, ConicalHypothesis, PlanarHypothesis, SphericalHypothesis, ToroidalHypothesis,
    MIN_CYLINDER_FACES, MIN_SPHERE_FACES, MIN_TORUS_FACES,
};

/// Max faces in the jet-similar seed patch used to condition the initial fit.
const SEED_PATCH_CAP: usize = 32;

/// Max trim rounds when conditioning a seed patch (one worst face per round).
const PATCH_TRIM_ROUNDS: usize = 16;

/// Floor on the patch size while conditioning — below this an analytic fit is
/// under-determined and trimming would chase noise.
const MIN_CONDITION_FACES: usize = 6;

/// Max BFS depth of the seed patch. Kept shallow so the patch stays local to
/// the seed and does not bleed across a tangent junction between two surfaces
/// of the same curvature class (e.g. a cone meeting a cylinder, where the
/// curvature is continuous across the seam) — a straddling patch fits neither
/// surface. An interior seed's shallow neighbourhood is on one surface, which
/// is all the initial fit needs; the flood then grows it to the true extent.
const MAX_PATCH_DEPTH: usize = 3;

/// Seed-patch radius bound, measured in radii of curvature of the seed's own
/// jet. `MAX_PATCH_DEPTH` bounds the patch in *hops*, which stops bounding it
/// in space once the facets are large: on a fillet network swept coarsely along
/// its own length, three hops walk many millimetres around the fillet chain and
/// straddle several of its segments. A surface of radius r is only locally
/// itself over a few r, and a few r of support is all any fit needs — the flood
/// grows the rest. Flats (jet magnitude ~0) are unbounded, as they should be.
const PATCH_RADII: f64 = 4.0;

/// Refit-and-reflood rounds while growing one hypothesis.
const MAX_GROW_ROUNDS: usize = 8;

// `GROW_VTOL_FACTOR` (30× vertex tolerance) and `RESIDUAL_SIGMAS` (5) lived
// here, together with `acceptance_tol`/`residual_scale`. The flood bar is now
// `--vertex-tolerance` itself. See DEVELOPMENT_PLAN, "The flood bar is the
// vertex tolerance".

/// Max BFS depth of the pure-bend fit patch grown past a stalled cylinder's
/// frontier (bent-tube promotion). A few rows of bend arc condition the
/// 7-DOF torus fit; the flood then claims the rest.
const BEND_PATCH_DEPTH: usize = 4;

/// Face cap of the pure-bend fit patch.
const BEND_PATCH_CAP: usize = 64;

/// How far a surface of revolution's support must sweep around its own axis, in
/// degrees, before the radius counts as determined by the data.
///
/// This is an *identifiability* test, and worth naming as such because the
/// distinction is easy to lose: it does not ask which surface fits best, it
/// asks whether the data pins the parameter at all. A fit whose parameter is
/// unpinned is not a slightly-worse explanation, it is a free variable that
/// happened to land somewhere, and no amount of comparing residuals or coverage
/// between candidates can detect it — see the note on complexity pricing in
/// TODO.md for why the obvious model-selection fix does not work here.
///
/// A patch subtending a negligible angle about the axis never observes the
/// revolution, so the radius is free — and at the limit of a free radius each
/// type collapses into the simpler one below it: a torus with R → ∞ is a
/// cylinder, a cylinder with r → ∞ is a plane. Such a fit is not merely
/// useless, it actively wins. Being a super-model of the simpler surface it
/// covers everything that surface covers, then floods a little further on the
/// slack — over the first tessellation row past a tangent junction, whose
/// chordal sag the simpler surface cannot take — and so carries the
/// largest-validated-area contest every time.
///
/// Measured on `stepped_block_toroidal_fillets` fine, where neither population
/// overlaps by anything close to this bar:
///
/// - **Tori.** The 18 real corner blends sweep 84.7°–103.4° and fit R to the
///   reference's exact 2.0000 / 2.5000 / 14.0000. Six impostors astride the
///   r = 0.5 and r = 1 perimeter fillets sweep 3.1°–6.7° and report
///   R = 20, 26, 29, 37, 87, 114 — junk, each an unconstrained parameter's
///   arbitrary landing place, and every one of them inside the model's own
///   121 mm bounding box, which is why bounding the radius never caught them.
/// - **Cylinders.** Every real one sweeps ≥ 64.8°. Four impostors sweep
///   2.8°–2.9°, each a flat wall plus the first tessellation row of the four
///   fillets bordering it (those rows sit 2.9e-3 mm off the wall plane, over
///   the flood tolerance, so the plane that should have won them cannot).
///
/// 20° sits 3× above every impostor and 3× below every real surface. Note that
/// sagitta-versus-tolerance does *not* separate them — the toroidal impostors'
/// spines still deviate 0.016–0.042 mm, well above the flood tolerance. The
/// angle is what matters, not the sag.
///
/// This is the counterpart to `MIN_HALF_ANGLE_DEG` in conical.rs, which rejects
/// a cone too shallow to be told apart from a cylinder.
const MIN_AXIAL_SWEEP_DEG: f64 = 20.0;

/// How many flood-start faces a toroidal candidate tries before giving up
/// (each attempt refits from the same initial patch, so extras only matter
/// when the best start happens not to validate its own fit).
const MAX_TORUS_STARTS: usize = 4;

/// A committed analytic hypothesis from competition, tagged by type.
pub(super) enum CompeteHyp {
    Planar(PlanarHypothesis),
    Conical(ConicalHypothesis),
    Spherical(SphericalHypothesis),
    Toroidal(ToroidalHypothesis),
}

impl CompeteHyp {
    pub(super) fn faces(&self) -> &[usize] {
        match self {
            CompeteHyp::Planar(h) => &h.faces,
            CompeteHyp::Conical(h) => &h.faces,
            CompeteHyp::Spherical(h) => &h.faces,
            CompeteHyp::Toroidal(h) => &h.faces,
        }
    }
    /// Short surface-type name for diagnostics.
    fn kind_name(&self) -> &'static str {
        match self {
            CompeteHyp::Planar(_) => "plane",
            CompeteHyp::Conical(h) => {
                if h.is_cylinder() {
                    "cyl"
                } else {
                    "cone"
                }
            }
            CompeteHyp::Spherical(_) => "sphere",
            CompeteHyp::Toroidal(_) => "torus",
        }
    }
    /// Simplest-first rank for claim resolution (a face both a cylinder and a
    /// torus explain within tolerance is the cylinder's). Matches the ladder's
    /// ascending-DOF order.
    fn dof_rank(&self) -> u8 {
        match self {
            CompeteHyp::Planar(_) => 0,
            CompeteHyp::Spherical(_) => 1,
            CompeteHyp::Conical(_) => 2,
            CompeteHyp::Toroidal(_) => 3,
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
            CompeteHyp::Toroidal(h) => {
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

    // Every region face seeds exactly once — `seeds` is a permutation of them,
    // which is what terminates the loop. Committing a hypothesis does NOT
    // remove the faces it grew into from the seed pool.
    //
    // Growing into a face is a hypothesis, not a conclusion. The committed
    // face set is always a `flood` output, so every member validates — but
    // against a fit whose parameters came from an unvalidated seed
    // neighbourhood, so "validated" does not mean "the right surface". Letting
    // those faces seed too means a wrong hypothesis can be out-competed by a
    // right one grown from a face it had swallowed; `resolve_claims` already
    // exists to arbitrate the multi-claims that produces.
    //
    // It also makes the long-standing claim that seed misordering "costs only
    // time" actually true: previously an early bad hypothesis consumed faces a
    // later, better seed needed.
    let mut committed: Vec<CompeteHyp> = Vec::new();
    let mut claim_lists: Vec<Vec<usize>> = vec![Vec::new(); n]; // face -> committed idxs
    let mut claimed = vec![false; n];

    for &seed in &seeds {
        // Skip only faces in the INTERIOR of what is already explained: every
        // neighbour claimed too. Seeding there cannot produce a competitor —
        // a fit grown from the middle of a validated surface refits that same
        // surface — so it is pure cost. Faces at the margin still seed.
        if claimed[seed] && interior_of_claims(mesh, seed, &claimed) {
            continue;
        }
        let jf = jet[seed];
        // Dominant / secondary principal-curvature MAGNITUDES. The jet orders
        // k1 ≥ k2 by signed value, so for a concave feature (a hole wall,
        // κ = −1/r) the dominant curvature is negative and lands in k2 —
        // taking |k1| as the dominant would read the wall's flat axial
        // direction and treat the hole as planar. Order by magnitude instead.
        let a1 = jf.k1.abs().max(jf.k2.abs());
        let a2 = jf.k1.abs().min(jf.k2.abs());

        // Grow the jet-pruned candidate types; keep the largest validated
        // area (simplest-DOF tiebreak).
        let mut best: Option<(CompeteHyp, f64, u8)> = None; // (hyp, area, dof)
        let consider = |hyp_area: Option<(CompeteHyp, f64)>, dof: u8, best: &mut Option<(CompeteHyp, f64, u8)>| {
            if let Some((h, area)) = hyp_area {
                let sweep = candidate_sweep(mesh, &h);
                if sweep < MIN_AXIAL_SWEEP_DEG.to_radians() {
                    if debug {
                        eprintln!(
                            "    compete: seed {seed} rejected {} sweeping only {:.1}\u{00b0} about its axis — radius undetermined ({} faces, area {area:.2})",
                            h.kind_name(),
                            sweep.to_degrees(),
                            h.faces().len(),
                        );
                    }
                    return;
                }
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
        let mut torus_area = None;
        if a1 >= flat_eps {
            // Curved: a developable (cylinder/cone) is always plausible.
            let cone = grow::<ConicalSurface>(mesh, &in_region, seed, jet, centroids, config, MIN_CYLINDER_FACES);
            cone_area = cone.as_ref().map(|(_, a, f)| (*a, f.len()));

            // Toroidal candidate: a direct grow from this seed (fillet tori,
            // bend tubes), skipped when the seed reads sphere-like — every
            // adjacent quad-row band of a sphere lies near some torus, so a
            // (κ,κ) seed's torus is at best a bogus band that loses to the
            // sphere anyway. A torus seed reads |κ2| < |κ1| (κ2 varies
            // azimuthally around the tube; it only approaches κ1 near the
            // degenerate R≈r hole of a donut — a whole-region ladder case,
            // not a competition case).
            let mut torus = if a2 < 0.8 * a1 {
                let patch = seed_patch(mesh, &in_region, seed, jet);
                torus_candidate(mesh, &in_region, &patch, &[seed], centroids, config)
            } else {
                None
            };
            // Bent-tube promotion: a cylinder that stalled against
            // tube-continuation faces is the developable prefix of a bent
            // tube; grow a torus from the pure-bend patch past the stall.
            if let Some((cyl, _, cyl_faces)) = cone.as_ref() {
                if cyl.is_cylinder() {
                    let promoted =
                        promote_bend(mesh, &in_region, cyl, cyl_faces, jet, centroids, config);
                    if let Some(p) = promoted {
                        if torus.as_ref().is_none_or(|(_, a, _)| p.1 > *a) {
                            torus = Some(p);
                        }
                    }
                }
            }
            torus_area = torus.as_ref().map(|(_, a, f)| (*a, f.len()));

            consider(cone.map(|(h, a, _)| (CompeteHyp::Conical(h), a)), 2, &mut best);
            // Curved in both directions → also try a sphere.
            if a2 >= 0.3 * a1 {
                let sphere = grow::<SphericalSurface>(mesh, &in_region, seed, jet, centroids, config, MIN_SPHERE_FACES);
                sphere_area = sphere.as_ref().map(|(_, a, f)| (*a, f.len()));
                consider(sphere.map(|(h, a, _)| (CompeteHyp::Spherical(h), a)), 1, &mut best);
            }
            consider(torus.map(|(h, a, _)| (CompeteHyp::Toroidal(h), a)), 3, &mut best);
        }
        let plane_area = plane.as_ref().map(|(_, a, f)| (*a, f.len()));
        consider(plane.map(|(h, a, _)| (CompeteHyp::Planar(h), a)), 0, &mut best);

        if debug && best.is_none() {
            eprintln!(
                "    compete: seed {seed} jet=({:.3},{:.3}) class {:?} → NO GROWTH (plane {plane_area:?} cone {cone_area:?} sphere {sphere_area:?} torus {torus_area:?})",
                jf.k1, jf.k2, curvature_class(jf),
            );
        }

        if let Some((hyp, area, dof)) = best {
            if debug {
                eprintln!(
                    "    compete: seed {seed} → {} ({} faces, area {area:.2}, dof {dof}, sweep {:.1}\u{00b0})",
                    hyp.kind_name(),
                    hyp.faces().len(),
                    match &hyp {
                        CompeteHyp::Toroidal(t) => torus_sweep(mesh, t).to_degrees(),
                        CompeteHyp::Conical(c) => conical_sweep(mesh, c).to_degrees(),
                        _ => 0.0,
                    }
                );
            }
            let idx = committed.len();
            for &f in hyp.faces() {
                claim_lists[f].push(idx);
                claimed[f] = true;
            }
            committed.push(hyp);
        }
        // A seed that grows nothing analytic costs only its own turn.
    }

    resolve_claims(
        mesh,
        &region.faces,
        &mut committed,
        &claim_lists,
        config.vertex_tolerance_mm,
        config,
    )
}

/// Resolve multi-claimed faces to their best-fitting claimant (2× margin to
/// leave the larger claimant), rebuild each hypothesis's face/vertex lists,
/// drop hypotheses that lost all faces, and collect unclaimed leftovers.
/// Is every neighbour of `fi` claimed too? Such a face is interior to the
/// explained set and cannot seed anything new.
fn interior_of_claims(mesh: &ConnectedMesh, fi: usize, claimed: &[bool]) -> bool {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    face.neighbors[..vc]
        .iter()
        .all(|&ni| ni < 0 || claimed[ni as usize])
}

/// Grow a cylinder whose axis DIRECTION is supplied by an already-committed
/// neighbour instead of being estimated from an unvalidated seed neighbourhood.
///
/// This is the point of stage 2.3(d2). A free 5-DOF cylinder fit on a leftover
/// arc fails on the **axis**, not the radius: measured on the three fillet
/// cylinders `stepped_block_toroidal` fine loses, fixing only the axis
/// direction recovers r = 0.99998 / 0.99999 / 0.50008 against true 1.0 / 1.0 /
/// 0.5, with residuals five to ten times *under* vertex tolerance. With the
/// axis fixed the fit is a 3-DOF circle in the perpendicular plane, which is
/// well conditioned on a short arc — so growth can start from the seed facet
/// alone and there is no blind patch to straddle a junction.
pub(super) fn grow_with_axis(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    seed: usize,
    axis: [f64; 3],
    centroids: &[[f64; 3]],
    config: &Config,
) -> Option<(ConicalHypothesis, f64, Vec<usize>)> {
    let vtol = config.vertex_tolerance_mm;
    let stol = config.surface_tolerance_mm;
    let mut accepted = vec![seed];
    let mut fit = None;
    for _ in 0..MAX_GROW_ROUNDS {
        let verts = vertices_of(mesh, &accepted);
        let (center, radius) =
            super::cylindrical::fit_circle_for_axis(&axis, &verts, &mesh.vertices)?;
        let mut f = ConeFit {
            apex: center,
            axis_direction: axis,
            half_angle: 0.0,
            radius,
            convex: None,
        };
        ConicalSurface::determine_convexity(mesh, &accepted, &mut f);
        let grown = flood::<ConicalSurface>(mesh, in_region, seed, &f, centroids, vtol, stol);
        if grown.is_empty() {
            return None; // the seed does not lie on its own constrained fit
        }
        let done = grown.len() == accepted.len();
        accepted = grown;
        fit = Some(f);
        if done {
            break;
        }
    }
    let fit = fit?;
    if accepted.len() < MIN_CYLINDER_FACES {
        return None;
    }
    let area: f64 = accepted
        .iter()
        .map(|&fi| face_area(&mesh.faces[fi], &mesh.vertices))
        .sum();
    let verts = vertices_of(mesh, &accepted);
    let hyp = ConicalSurface::make_hypothesis(mesh, accepted.clone(), verts, &fit);
    Some((hyp, area, accepted))
}

fn resolve_claims(
    mesh: &ConnectedMesh,
    region_faces: &[usize],
    committed: &mut Vec<CompeteHyp>,
    claim_lists: &[Vec<usize>],
    vtol: f64,
    config: &Config,
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
        let dist = |ci: usize| max_vertex_dist(mesh, &committed[ci], fi);
        // Among claimants that ALREADY explain the face within vertex
        // tolerance, prefer the simplest surface (lowest DOF), then largest
        // area, then lowest index: a face on both a cylinder and a torus near
        // their tangency is the cylinder's; a face on a real fillet and a
        // spurious thin junction torus (both sub-micron) is the larger
        // fillet's. Comparing absolute best-fit here instead would let a tiny
        // tighter sliver steal it and fragment one surface into a chain.
        let within: Vec<usize> = claimants.iter().copied().filter(|&ci| dist(ci) <= vtol).collect();
        let winner = if !within.is_empty() {
            *within
                .iter()
                .min_by(|&&a, &&b| {
                    committed[a]
                        .dof_rank()
                        .cmp(&committed[b].dof_rank())
                        .then_with(|| hyp_area[b].partial_cmp(&hyp_area[a]).unwrap())
                        .then_with(|| a.cmp(&b))
                })
                .unwrap()
        } else {
            // No claimant fits within tolerance — best-fit, with a 2× margin
            // to leave the larger-area claimant (a true junction face fits
            // both sides equally; the margin prevents flapping).
            let mut order: Vec<usize> = claimants.clone();
            order.sort_by(|&a, &b| {
                hyp_area[b].partial_cmp(&hyp_area[a]).unwrap().then_with(|| a.cmp(&b))
            });
            let mut winner = order[0];
            let mut winner_d = dist(winner);
            for &ci in &order[1..] {
                let d = dist(ci);
                if d * 2.0 < winner_d {
                    winner = ci;
                    winner_d = d;
                }
            }
            winner
        };
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
        // The fit was computed on the support the flood grew, which included
        // every face this hypothesis could reach — junction rows it has now
        // lost to a better claimant among them. Those rows are exactly the ones
        // that pull a fit off true, so a hypothesis that has just been trimmed
        // is carrying the error of faces it no longer owns. Refit on what it
        // actually kept. (Measured on rounded_cube fine, whose r = 2 edge
        // fillets read 1.9976–1.9984 with no vertex closer than 5e-4 mm to
        // their own surface, against a mesh whose vertices are 1.9e-6 mm from
        // the true CAD surfaces; refitting takes them to 2.0000.)
        refit_in_place(mesh, &mut hyp, config);
        out.push(hyp);
    }
    (out, leftovers)
}


/// Refit a committed hypothesis to its current face set, keeping the existing
/// fit if the fitter cannot produce one.
fn refit_in_place(mesh: &ConnectedMesh, hyp: &mut CompeteHyp, config: &Config) {
    let faces = hyp.faces().to_vec();
    let verts = vertices_of(mesh, &faces);
    let v = config.verbosity;
    match hyp {
        CompeteHyp::Planar(h) => {
            if let Some(fit) = PlanarSurface::fit_surface(mesh, &faces, &verts, v) {
                *h = PlanarSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        CompeteHyp::Conical(h) => {
            if let Some(mut fit) = ConicalSurface::fit_surface(mesh, &faces, &verts, v) {
                ConicalSurface::determine_convexity(mesh, &faces, &mut fit);
                *h = ConicalSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        CompeteHyp::Spherical(h) => {
            if let Some(mut fit) = SphericalSurface::fit_surface(mesh, &faces, &verts, v) {
                SphericalSurface::determine_convexity(mesh, &faces, &mut fit);
                *h = SphericalSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
        CompeteHyp::Toroidal(h) => {
            if let Some(mut fit) = ToroidalSurface::fit_surface(mesh, &faces, &verts, v) {
                ToroidalSurface::determine_convexity(mesh, &faces, &mut fit);
                *h = ToroidalSurface::make_hypothesis(mesh, faces, verts, &fit);
            }
        }
    }
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
    let patch = seed_patch(mesh, in_region, seed, jet);
    grow_from::<S>(mesh, in_region, seed, patch, centroids, config, min_faces)
}

/// The grow loop proper, from an explicit initial fit patch and flood start:
/// fit the patch, flood-accept every region face the fit validates, refit on
/// the grown support, to a fixpoint. Returns the hypothesis, its validated
/// area, and its face list.
fn grow_from<S: SurfaceType>(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    start: usize,
    initial: Vec<usize>,
    centroids: &[[f64; 3]],
    config: &Config,
    min_faces: usize,
) -> Option<(S::Hypothesis, f64, Vec<usize>)> {
    let stol = config.surface_tolerance_mm;
    let vtol = config.vertex_tolerance_mm;
    let mut accepted =
        condition_patch::<S>(mesh, start, initial, centroids, config, vtol, stol, min_faces);
    let mut prev = 0usize;
    for round in 0..MAX_GROW_ROUNDS {
        let verts = vertices_of(mesh, &accepted);
        let mut fit = S::fit_surface(mesh, &accepted, &verts, config.verbosity)?;
        S::determine_convexity(mesh, &accepted, &mut fit);
        let comp = flood::<S>(mesh, in_region, start, &fit, centroids, vtol, stol);
        if comp.is_empty() {
            return None; // the start face does not lie on its own fit — reject
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
    let faces = flood::<S>(mesh, in_region, start, &fit, centroids, vtol, stol);
    if faces.len() < min_faces {
        return None;
    }
    let area: f64 = faces.iter().map(|&fi| face_area(&mesh.faces[fi], &mesh.vertices)).sum();
    let verts = vertices_of(mesh, &faces);
    let hyp = S::make_hypothesis(mesh, faces.clone(), verts, &fit);
    Some((hyp, area, faces))
}


/// A committed torus's minor radius may not exceed this fraction of its major
/// radius. Real fillets and bend tubes are slim (r/R ≤ 0.5 across the corpus:
/// full donut 0.3, pipe-elbow 0.5, fillets 0.15–0.2); an r ≈ R "fat" torus is
/// a bogus fit a flat cap or wall band accepts loosely and then wins on area
/// (filleted_pipe: an r/R = 0.82 torus ate an annular end cap).
const MAX_TORUS_RADIUS_RATIO: f64 = 0.7;

/// How far a support sweeps around an axis, in radians: the angular extent of
/// `verts` about the line through `origin` along `axis`, taken as 2π minus the
/// largest angular gap between consecutive vertices, so a patch straddling the
/// angle origin measures correctly.
///
/// This is the observable that says whether a rotational parameter is
/// determined by the data at all. A surface of revolution fitted to a patch
/// that never goes around its axis has its radius free, and a free parameter is
/// worse than useless in competition: the fit becomes a super-model of the
/// simpler surface it should have lost to — covering everything that surface
/// covers, then flooding a little further on the slack — so it takes the
/// largest-validated-area contest every time. See `MIN_TORUS_SWEEP_DEG` and
/// `MIN_CYLINDER_SWEEP_DEG`.
fn axial_sweep(
    mesh: &ConnectedMesh,
    verts: &[usize],
    origin: &[f64; 3],
    axis: &[f64; 3],
) -> f64 {
    let (u, w) = super::perpendicular_basis(axis);
    let mut angles: Vec<f64> = verts
        .iter()
        .filter_map(|&vi| {
            let v = &mesh.vertices[vi];
            let d = [v.x - origin[0], v.y - origin[1], v.z - origin[2]];
            let (x, y) = (
                d[0] * u[0] + d[1] * u[1] + d[2] * u[2],
                d[0] * w[0] + d[1] * w[1] + d[2] * w[2],
            );
            // A vertex on the axis has no meaningful angle.
            (x.hypot(y) > 1e-12).then(|| y.atan2(x))
        })
        .collect();
    if angles.len() < 2 {
        return 0.0;
    }
    angles.sort_by(|p, q| p.partial_cmp(q).unwrap());
    let mut gap = angles[0] + std::f64::consts::TAU - angles[angles.len() - 1];
    for pair in angles.windows(2) {
        gap = gap.max(pair[1] - pair[0]);
    }
    (std::f64::consts::TAU - gap).max(0.0)
}

/// Sweep of a candidate's support about the axis whose radius it determines.
/// Planes and spheres have no axis — a sphere's radius is pinned by any patch —
/// so they report a full turn and are never rejected on this ground.
fn candidate_sweep(mesh: &ConnectedMesh, hyp: &CompeteHyp) -> f64 {
    match hyp {
        CompeteHyp::Toroidal(h) => torus_sweep(mesh, h),
        CompeteHyp::Conical(h) => conical_sweep(mesh, h),
        CompeteHyp::Planar(_) | CompeteHyp::Spherical(_) => std::f64::consts::TAU,
    }
}

/// Sweep of a toroidal hypothesis's support about its own axis.
fn torus_sweep(mesh: &ConnectedMesh, hyp: &ToroidalHypothesis) -> f64 {
    axial_sweep(mesh, &hyp.vertices, &hyp.center, &hyp.axis_direction)
}

/// Sweep of a cylinder's support about its own axis. Meaningless for a true
/// cone, whose apex pins the axis regardless of sweep, so cones report a full
/// turn.
fn conical_sweep(mesh: &ConnectedMesh, hyp: &ConicalHypothesis) -> f64 {
    if !hyp.is_cylinder() {
        return std::f64::consts::TAU;
    }
    axial_sweep(mesh, &hyp.vertices, &hyp.apex, &hyp.axis_direction)
}

/// Grow a toroidal candidate from an initial fit patch, trying flood starts
/// in order until one grows. Rejects degenerate giant-major-radius fits: a
/// short segment of straight cylinder fits a torus with an absurdly distant
/// axis (observed R ≈ 1500–3600 mm on ~100 mm models under the old BFS), and
/// a legitimate torus's major radius cannot exceed the mesh extent — the
/// same guard as the ladder's acceptance test. Also rejects "fat" tori
/// (`r > MAX_TORUS_RADIUS_RATIO·R`). The sweep guard that rejects a torus too
/// straight to have a major radius at all is applied to every candidate type
/// alike, in `compete_region` (see `MIN_AXIAL_SWEEP_DEG`).
fn torus_candidate(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    initial: &[usize],
    starts: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
) -> Option<(ToroidalHypothesis, f64, Vec<usize>)> {
    let bb_diag = bounding_box_diagonal(&mesh.vertices);
    for &start in starts.iter().take(MAX_TORUS_STARTS) {
        if let Some((hyp, area, faces)) = grow_from::<ToroidalSurface>(
            mesh,
            in_region,
            start,
            initial.to_vec(),
            centroids,
            config,
            MIN_TORUS_FACES,
        ) {
            // A grow that succeeded but produced an implausible torus says the
            // *type* is wrong here, so stop — retrying other flood starts from
            // the same fit patch would only rediscover it. Extra starts are
            // for the case where the grow itself failed.
            return (hyp.major_radius <= bb_diag
                && hyp.minor_radius <= MAX_TORUS_RADIUS_RATIO * hyp.major_radius)
            .then_some((hyp, area, faces));
        }
    }
    None
}

/// Stage 2.3(c) bent-tube promotion. A cylinder flood that stalled against
/// region faces whose jet continues the tube — κ1 ≈ the tube's 1/r, |κ2|
/// well below κ1 (the same tube swinging around a bend, not a cap plane's
/// (0,0) or an end sphere's (κ,κ)) — is the developable prefix of a bent
/// tube. Fit a torus to the pure-bend patch just past the stall line
/// (frontier faces plus their off-cylinder tube-jet neighbours) and grow it;
/// the result competes on validated area like any other candidate.
///
/// The cylinder fit itself locates the tangency — frontier faces failed
/// cylinder validation — so unlike a jet-class seed patch (which cannot tell
/// bend from straight tube: both read the same κ1) this patch cannot
/// straddle the junction.
fn promote_bend(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    cyl: &ConicalHypothesis,
    cyl_faces: &[usize],
    jet: &[JetField],
    centroids: &[[f64; 3]],
    config: &Config,
) -> Option<(ToroidalHypothesis, f64, Vec<usize>)> {
    if cyl.radius_at_apex <= 0.0 {
        return None;
    }
    let k_tube = 1.0 / cyl.radius_at_apex;
    let grow_vtol = config.vertex_tolerance_mm;

    // Max vertex distance to the stalled cylinder — the "past the tangency"
    // test (frontier faces failed the flood at exactly this tolerance).
    let off_cyl_dist = |fi: usize| -> f64 {
        let face = &mesh.faces[fi];
        let vc = face.vertex_count as usize;
        face.vertex_indices[..vc]
            .iter()
            .map(|&vi| {
                super::cylindrical::vertex_to_cylinder_distance(
                    &mesh.vertices[vi],
                    &cyl.apex,
                    &cyl.axis_direction,
                    cyl.radius_at_apex,
                )
                .abs()
            })
            .fold(0.0_f64, f64::max)
    };
    // Tube-continuation jet signature. The κ2 bound is generous: a bend's
    // inner (saddle) side reads |κ2| = 1/(R−r), which approaches κ1 = 1/r
    // only for very tight bends; junction rows have contaminated jets and
    // are excluded from the fit patch either way (the flood claims them).
    let tube_jet = |fi: usize| -> bool {
        let a1 = jet[fi].k1.abs();
        a1 >= 0.5 * k_tube && a1 <= 2.0 * k_tube && jet[fi].k2.abs() <= 0.8 * a1
    };

    // Frontier: in-region neighbours of the cylinder support that are past
    // the tangency and continue the tube.
    let mut seen: HashSet<usize> = cyl_faces.iter().copied().collect();
    let mut frontier: Vec<usize> = Vec::new();
    for &fi in cyl_faces {
        let vc = mesh.faces[fi].vertex_count as usize;
        for &ni in &mesh.faces[fi].neighbors[..vc] {
            if ni < 0 {
                continue;
            }
            let ni = ni as usize;
            if !in_region[ni] || !seen.insert(ni) {
                continue;
            }
            if off_cyl_dist(ni) > grow_vtol && tube_jet(ni) {
                frontier.push(ni);
            }
        }
    }
    if frontier.is_empty() {
        return None; // finished a straight pipe (or stalled at a cap/sphere)
    }
    frontier.sort_unstable();

    // Pure-bend fit patch: depth-limited multi-source BFS from the frontier
    // through off-cylinder tube-jet faces. (Neighbours examined and rejected
    // above stay excluded — junction rows with contaminated jets would only
    // pollute the fit.)
    let mut patch: Vec<usize> = frontier.clone();
    let mut ring = frontier.clone();
    for _depth in 0..BEND_PATCH_DEPTH {
        if patch.len() >= BEND_PATCH_CAP {
            break;
        }
        let mut next: Vec<usize> = Vec::new();
        for &fi in &ring {
            let vc = mesh.faces[fi].vertex_count as usize;
            for &ni in &mesh.faces[fi].neighbors[..vc] {
                if ni < 0 {
                    continue;
                }
                let ni = ni as usize;
                if !in_region[ni] || seen.contains(&ni) {
                    continue;
                }
                if off_cyl_dist(ni) <= grow_vtol || !tube_jet(ni) {
                    continue;
                }
                seen.insert(ni);
                patch.push(ni);
                next.push(ni);
                if patch.len() >= BEND_PATCH_CAP {
                    break;
                }
            }
        }
        if next.is_empty() {
            break;
        }
        ring = next;
    }
    patch.sort_unstable();

    // Flood starts: deepest-off-cylinder first — the faces most surely
    // interior to the bend (frontier faces sit on the junction row and may
    // not validate the torus fit themselves).
    let mut starts = patch.clone();
    starts.sort_by(|&a, &b| {
        off_cyl_dist(b)
            .partial_cmp(&off_cyl_dist(a))
            .unwrap()
            .then_with(|| a.cmp(&b))
    });

    torus_candidate(mesh, in_region, &patch, &starts, centroids, config)
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

/// Max distance from any vertex of face `fi` to `fit` — the trim criterion for
/// conditioning the seed patch (`condition_patch`).
fn face_max_vertex_dist<S: SurfaceType>(mesh: &ConnectedMesh, fi: usize, fit: &S::Fit) -> f64 {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    face.vertex_indices[..vc]
        .iter()
        .map(|&vi| {
            let v = &mesh.vertices[vi];
            S::vertex_distance(fit, &[v.x, v.y, v.z])
        })
        .fold(0.0_f64, f64::max)
}

/// Condition the seed patch so the seed lies on the patch's own fit. A
/// depth-limited jet-similar patch can bleed a row or two across a tangent
/// junction between two surfaces of the same curvature class — a short
/// cylinder wall into its fillet, where the transition faces' κ1 ramps from
/// the wall's toward the fillet's — which pulls the fit off far enough that
/// the seed fails its own fit and the grow gives up (leaving the wall to be
/// eaten by an over-grown torus). Trim the worst-fitting patch face (never the
/// seed) and refit until the seed validates, bottoming out at a size floor so
/// the fit stays well-determined. Fit-validated, in the spirit of "validation
/// decides geometry"; a no-op when the patch is already clean.
///
/// When the patch admits no fit of this type at all there is no residual to
/// rank by, so the round instead drops the face geometrically farthest from
/// the seed: a straddling patch is unfittable precisely because its far end
/// sits on another surface, and shrinking toward the seed is the one move that
/// is monotonically local.
#[allow(clippy::too_many_arguments)]
fn condition_patch<S: SurfaceType>(
    mesh: &ConnectedMesh,
    seed: usize,
    mut patch: Vec<usize>,
    centroids: &[[f64; 3]],
    config: &Config,
    vtol: f64,
    stol: f64,
    min_faces: usize,
) -> Vec<usize> {
    let floor = min_faces.max(MIN_CONDITION_FACES);
    for _ in 0..PATCH_TRIM_ROUNDS {
        if patch.len() <= floor {
            break;
        }
        let verts = vertices_of(mesh, &patch);
        let mut fit = match S::fit_surface(mesh, &patch, &verts, config.verbosity) {
            Some(f) => f,
            None => {
                // No fit to rank by — shrink toward the seed instead.
                let far = patch
                    .iter()
                    .copied()
                    .filter(|&f| f != seed)
                    .max_by(|&a, &b| {
                        centroid_gap(centroids, seed, a)
                            .partial_cmp(&centroid_gap(centroids, seed, b))
                            .unwrap()
                    });
                match far {
                    Some(f) => {
                        patch.retain(|&p| p != f);
                        continue;
                    }
                    None => break,
                }
            }
        };
        S::determine_convexity(mesh, &patch, &mut fit);
        if face_validates::<S>(mesh, seed, &fit, centroids, vtol, stol) {
            break;
        }
        let worst = patch
            .iter()
            .copied()
            .filter(|&f| f != seed)
            .max_by(|&a, &b| {
                face_max_vertex_dist::<S>(mesh, a, &fit)
                    .partial_cmp(&face_max_vertex_dist::<S>(mesh, b, &fit))
                    .unwrap()
            });
        match worst {
            Some(w) => patch.retain(|&f| f != w),
            None => break,
        }
    }
    patch
}

/// Distance between two faces' centroids.
fn centroid_gap(centroids: &[[f64; 3]], a: usize, b: usize) -> f64 {
    let (p, q) = (centroids[a], centroids[b]);
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
}

/// A small patch around `seed`, capped and depth-limited, for a well-conditioned
/// initial fit: jet-similar faces (same curvature class, magnitude within 2×,
/// and within `PATCH_RADII` radii of curvature) where the jet is informative.
///
/// Falls back to the plain topological neighbourhood when that leaves too
/// little to fit. The jet is a *prior*, not a gate — validation decides every
/// face — and on a small fillet tessellated into a handful of facets it is
/// simply not reliable enough to gate on: measured on
/// `stepped_block_toroidal_fillets` medium, adjacent facets of one exact
/// r = 0.5 arc read κ₁ = 0.023, 0.162 and 2.2 (radii of 43 mm, 6 mm and
/// 0.45 mm for the same 0.5 mm fillet). The similarity filter then admits 2–4
/// faces, too few to determine a 5-DOF cylinder, so every candidate fails and
/// the arc degrades to one plane per facet. Handing the fit an unfiltered
/// neighbourhood instead costs nothing in soundness: a patch straddling a
/// junction is what `condition_patch` trims and what flood validation rejects.
fn seed_patch(mesh: &ConnectedMesh, in_region: &[bool], seed: usize, jet: &[JetField]) -> Vec<usize> {
    let patch = grow_patch(mesh, in_region, seed, jet, true);
    if patch.len() >= MIN_CONDITION_FACES {
        return patch;
    }
    grow_patch(mesh, in_region, seed, jet, false)
}

/// Depth-limited BFS around `seed`, optionally filtered by jet similarity.
fn grow_patch(
    mesh: &ConnectedMesh,
    in_region: &[bool],
    seed: usize,
    jet: &[JetField],
    use_jet: bool,
) -> Vec<usize> {
    let seed_class = curvature_class(jet[seed]);
    let seed_mag = jet[seed].k1.abs().max(jet[seed].k2.abs());
    let seed_c = super::face_centroid(&mesh.faces[seed], &mesh.vertices);
    let radius_limit = if use_jet && seed_mag > 1e-9 {
        PATCH_RADII / seed_mag
    } else {
        f64::INFINITY
    };
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
                if use_jet {
                    if curvature_class(jet[ni]) != seed_class {
                        continue;
                    }
                    let mag = jet[ni].k1.abs().max(jet[ni].k2.abs());
                    if seed_mag > 1e-9 && (mag > 1.5 * seed_mag || mag < 0.667 * seed_mag) {
                        continue;
                    }
                    let c = super::face_centroid(&mesh.faces[ni], &mesh.vertices);
                    let d = ((c[0] - seed_c[0]).powi(2)
                        + (c[1] - seed_c[1]).powi(2)
                        + (c[2] - seed_c[2]).powi(2))
                    .sqrt();
                    if d > radius_limit {
                        continue;
                    }
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
    // Magnitude-ordered so a concave feature (dominant curvature negative, in
    // k2) is not misread — see the note in `compete_region`.
    let a1 = jf.k1.abs().max(jf.k2.abs());
    let a2 = jf.k1.abs().min(jf.k2.abs());
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
                CompeteHyp::Toroidal(h) => super::toroidal::vertex_to_torus_distance(
                    v,
                    &h.center,
                    &h.axis_direction,
                    h.major_radius,
                    h.minor_radius,
                )
                .abs(),
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
