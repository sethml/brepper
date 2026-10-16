//! Stage 2.3(a): the whole-region DOF ladder.
//!
//! Try to explain each smooth region (stage 2.2) with ONE surface, ascending
//! degrees of freedom; the first fit that validates wins. Simplest-first
//! resolves degeneracy (any plane also lies on an infinite-radius cylinder; a
//! single facet is exactly a plane) without coverage heuristics.
//!
//! | Rung | Surface  | DOF |
//! |------|----------|-----|
//! | 1    | Plane    | 3   |
//! | 2    | Sphere   | 4   |
//! | 3    | Conical  | 5/6 | (cylinder-or-cone; the fitter prefers cylinder)
//! | 4    | Torus    | 7   |
//!
//! The conical rung fits a cylinder first and only falls to a cone when the
//! cylinder does not validate (`ConicalSurface::fit_surface`), so it covers
//! both the 5-DOF cylinder and 6-DOF cone cases in DOF order.
//!
//! **Jet-informed order between the sphere and conical rungs.** Strict
//! simplest-first breaks on a short cylinder/cone tessellated as a few rings:
//! two same-radius rings at heights ±h lie exactly on a sphere of radius
//! √(R²+h²), so the lower-DOF sphere would validate and wrongly win. The jet
//! resolves it — a developable region (median `min/max` principal-curvature
//! ratio ~0) tries the conical rung before the sphere rung; a bi-curved
//! region (ratio ~1) keeps sphere-first. The jet only ORDERS; a genuine
//! sphere it misreads still falls through to the sphere rung, since a
//! cylinder fit cannot validate a sphere's vertices (`region_is_developable`).
//!
//! **Validation decides all geometry** (Design Principle #1): every claimed
//! vertex within `--vertex-tolerance`, every claimed face centroid within
//! `--surface-tolerance`, and orientation coherence (every face normal
//! agrees in sign with the surface normal at its centroid — the per-type
//! `convexity_check`, which also fixes the hypothesis `convex` flag). The
//! fitters apply their own numeric-validity guards (cone half-angle bounds,
//! torus `0 < r ≤ R`, sphere radius cap) and return `None` outside them.
//!
//! A region the ladder cannot explain with one surface is a tangent chain
//! (pill capsule, rounded-cube shell, pipe-elbow tube); it goes to stage
//! 2.3(b)/(c) scoped competition (`compete.rs`). Whatever competition leaves
//! unclaimed **terminally degrades to per-facet planes** — total coverage
//! can never fail, and a single triangle is exactly a plane.

use crate::config::Config;
use crate::stage1::ConnectedMesh;

use super::surface::SurfaceType;
use super::bspline::{self, BSplineHypothesis};
use super::compete::{self, CompeteHyp};
use super::conical::ConicalSurface;
use super::edges::JetField;
use super::planar::{fit_plane, PlanarSurface};
use super::regions::Region;
use super::spherical::SphericalSurface;
use super::toroidal::ToroidalSurface;
use super::{
    face_centroid, ConicalHypothesis, PlanarHypothesis, SelectedSurface, SphericalHypothesis,
    ToroidalHypothesis, MIN_CYLINDER_FACES, MIN_SPHERE_FACES, MIN_TORUS_FACES,
};

use std::collections::HashSet;

/// Regions larger than this many faces are fitted on a deterministic
/// stride subsample (validated against the whole region, refit-augmented on
/// a near-miss). Bounds fit cost on huge regions — which, lacking scoped
/// competition, are the tangent chains that fail the ladder anyway.
const MAX_FIT_FACES: usize = 3000;

/// Bounded refit rounds: after a subsample fit, the worst validation
/// offenders are added to the fit set and the surface is refit, so a large
/// clean surface still validates end to end.
const MAX_REFIT: usize = 3;

/// How many worst-offender faces to fold into the fit set per refit round.
const REFIT_ADD: usize = 32;

/// Give-up threshold for a refit near-miss: if the worst offender exceeds
/// this multiple of `--surface-tolerance`, the region is decisively not this
/// surface and no amount of augmenting will help.
const REFIT_GIVE_UP: f64 = 8.0;

// The former `LEFTOVER_BLOB_MAX` whole-region discard lived here. It is gone;
// `dissolve_fragmented_surfaces` replaces it. See DEVELOPMENT_PLAN,
// "The leftover blob is a boundary problem, not a size problem".

/// The output of stage 2.3(a): committed hypotheses in the typed vectors plus
/// one `SelectedSurface` per committed surface. Regions partition the mesh,
/// so every mesh face lands in exactly one selected surface.
pub(super) struct Classification {
    pub planar: Vec<PlanarHypothesis>,
    pub conical: Vec<ConicalHypothesis>,
    pub spherical: Vec<SphericalHypothesis>,
    pub toroidal: Vec<ToroidalHypothesis>,
    pub bspline: Vec<BSplineHypothesis>,
    pub selected: Vec<SelectedSurface>,
}

impl Classification {
    fn new() -> Self {
        Classification {
            planar: Vec::new(),
            conical: Vec::new(),
            spherical: Vec::new(),
            toroidal: Vec::new(),
            bspline: Vec::new(),
            selected: Vec::new(),
        }
    }

    fn push_planar(&mut self, hyp: PlanarHypothesis) {
        self.selected.push(SelectedSurface::Planar(self.planar.len()));
        self.planar.push(hyp);
    }
    fn push_conical(&mut self, hyp: ConicalHypothesis) {
        self.selected
            .push(SelectedSurface::Conical(self.conical.len()));
        self.conical.push(hyp);
    }
    fn push_spherical(&mut self, hyp: SphericalHypothesis) {
        self.selected
            .push(SelectedSurface::Spherical(self.spherical.len()));
        self.spherical.push(hyp);
    }
    fn push_toroidal(&mut self, hyp: ToroidalHypothesis) {
        self.selected
            .push(SelectedSurface::Toroidal(self.toroidal.len()));
        self.toroidal.push(hyp);
    }
    fn push_bspline(&mut self, hyp: BSplineHypothesis) {
        self.selected
            .push(SelectedSurface::BSpline(self.bspline.len()));
        self.bspline.push(hyp);
    }
    fn push_compete(&mut self, hyp: CompeteHyp) {
        match hyp {
            CompeteHyp::Planar(h) => self.push_planar(h),
            CompeteHyp::Conical(h) => self.push_conical(h),
            CompeteHyp::Spherical(h) => self.push_spherical(h),
            CompeteHyp::Toroidal(h) => self.push_toroidal(h),
        }
    }
}

/// Run the whole-region ladder over every smooth region; regions the ladder
/// cannot explain go to scoped competition (2.3(b)), and whatever competition
/// leaves unclaimed degrades to per-facet planes.
pub(super) fn classify_regions(
    mesh: &ConnectedMesh,
    regions: &[Region],
    jet: &[JetField],
    face_region: &[usize],
    config: &Config,
) -> Classification {
    let centroids: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|f| face_centroid(f, &mesh.vertices))
        .collect();

    let mut cls = Classification::new();
    let mut ladder_wins = [0usize; 4]; // plane, sphere, conical, torus
    let mut competed_regions = 0usize;
    let mut competed_surfaces = 0usize;
    let mut freeform_surfaces = 0usize;
    let mut facet_faces = 0usize;
    for region in regions {
        if classify_one_region(mesh, region, &centroids, jet, config, &mut cls, &mut ladder_wins) {
            continue;
        }
        // The ladder cannot explain this region — scoped competition.
        competed_regions += 1;
        let (mut hyps, leftovers) = compete::compete_region(mesh, region, jet, &centroids, &compete::CompeteConfig::from(config));

        // Stage 2.3(d): run the whole-region ladder on each connected
        // component of the leftovers.
        //
        // Competition seeds a fit from an unvalidated neighbourhood of one
        // face, so a surface whose every seed patch straddles a junction is
        // never fitted at all — measured on `stepped_block_toroidal` fine,
        // three fillet cylinders are lost this way and reach stage 3 as
        // B-splines. But the *component* those faces form is exactly the kind
        // of input the ladder handles: given the right face set the existing
        // conical fitter estimates its axis from the normal covariance and
        // recovers r = 0.99998 / 0.99999 / 0.50008 against true 1.0 / 1.0 /
        // 0.5, at five to ten times under vertex tolerance. Nothing new is
        // needed — only feeding the fitter a face set that is not
        // contaminated by a neighbouring surface.
        //
        // This can only *add* surfaces: it runs on faces nothing else
        // explained, and every rung is vertex-tolerance validated, so a
        // component that is genuinely freeform simply falls through to
        // 2.3(e) as before.
        // Stages 2.3(d) and 2.3(d2), alternated to a fixpoint. Each changes
        // the leftover set the other sees: (d2) splitting a blob can leave a
        // remainder one ladder rung explains, and (d) committing a surface
        // gives (d2) another axis to borrow. Terminating — a round that
        // commits nothing ends the loop, and every round strictly shrinks the
        // leftovers.
        let mut leftovers = leftovers;
        loop {
            let before = leftovers.len();
            leftovers =
                rescue_leftover_components(mesh, &leftovers, &centroids, config, &mut hyps);
            leftovers =
                rescue_with_neighbour_axes(mesh, &leftovers, &centroids, config, &mut hyps);
            if leftovers.len() == before {
                break;
            }
        }

        // Stage 2.3(e): freeform-fit the leftover components competition
        // could not explain analytically, dissolving any small analytic
        // island enclosed by a validating B-spline patch.
        let ff = bspline::fit_leftover_components(mesh, &leftovers, &hyps, face_region, config);

        // A committed analytic surface whose boundary with the per-facet
        // residual is more fragmented than the surface itself reconstructs to
        // an invalid wire; dissolve exactly those, and leave the rest of the
        // region's decomposition alone.
        let mut residual = ff.residual;
        let fragmented = dissolve_fragmented_surfaces(mesh, &hyps, &ff.dissolved, &mut residual);
        let surviving_hyps = hyps.len() - ff.dissolved.len() - fragmented.len();
        if config.verbosity >= 3 || std::env::var("BREPPER_DEBUG_COMPETE").is_ok() {
            eprintln!(
                "    compete: region {} faces → {surviving_hyps} surviving hyps, {} leftovers, freeform residual {} (largest blob {}, {} dissolved as over-fragmented)",
                region.faces.len(),
                leftovers.len(),
                residual.len(),
                largest_component(mesh, &residual),
                fragmented.len(),
            );
        }

        // Commit surviving analytic hypotheses (skipping dissolved islands),
        // the freeform surfaces, and the residual leftovers as per-facet
        // planes (total coverage can never fail; a triangle is exactly a
        // plane).
        for (i, hyp) in hyps.into_iter().enumerate() {
            if !ff.dissolved.contains(&i) && !fragmented.contains(&i) {
                cls.push_compete(hyp);
                competed_surfaces += 1;
            }
        }
        for plane in ff.planes {
            cls.push_planar(plane);
        }
        for hyp in ff.bsplines {
            cls.push_bspline(hyp);
            freeform_surfaces += 1;
        }
        facet_faces += residual.len();
        for fi in residual {
            cls.push_planar(facet_plane(mesh, fi));
        }
    }

    if config.verbosity >= 1 {
        eprintln!(
            "  ladder: {} plane, {} sphere, {} cyl/cone, {} torus; competition: {} region(s) → {} analytic + {} freeform + {} per-facet",
            ladder_wins[0], ladder_wins[1], ladder_wins[2], ladder_wins[3],
            competed_regions, competed_surfaces, freeform_surfaces, facet_faces,
        );
    }
    cls
}

/// Median `min(|κ1|,|κ2|)/max(|κ1|,|κ2|)` at or below which a region is treated
/// as developable (cylinder/cone) for rung ordering: a cylinder/cone reads ~0
/// (one principal curvature vanishes), a sphere reads ~1 (both equal). Well
/// clear of both classes' jet noise. Ordering only — never rejection.
const DEVELOPABLE_MAX_RATIO: f64 = 0.4;

/// Whether the region's jet reads developable (one principal curvature ~0),
/// used only to order the sphere and conical rungs. A short cylinder/cone
/// tessellated as a few rings is vertex-ambiguous with a sphere — two
/// same-radius rings at heights ±h lie exactly on a sphere of radius
/// √(R²+h²) — so simplest-DOF-first would wrongly commit the sphere; the jet
/// breaks the tie by curvature class. Faces too flat to have a curvature
/// signal are skipped; with no signal at all the region is not treated as
/// developable (keeps the plane→sphere→conical default order).
fn region_is_developable(region: &Region, jet: &[JetField]) -> bool {
    let mut ratios: Vec<f64> = region
        .faces
        .iter()
        .filter_map(|&fi| {
            let j = &jet[fi];
            let (a, b) = (j.k1.abs(), j.k2.abs());
            let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
            if hi > 1e-9 {
                Some(lo / hi)
            } else {
                None
            }
        })
        .collect();
    if ratios.is_empty() {
        return false;
    }
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ratios[ratios.len() / 2] <= DEVELOPABLE_MAX_RATIO
}

/// Try the DOF ladder on one region; commit the first surface that validates.
/// Returns true on a committed surface, false if the region fell through
/// (caller applies the per-facet fallback).
fn classify_one_region(
    mesh: &ConnectedMesh,
    region: &Region,
    centroids: &[[f64; 3]],
    jet: &[JetField],
    config: &Config,
    cls: &mut Classification,
    wins: &mut [usize; 4],
) -> bool {
    let faces = &region.faces;
    let verts = &region.vertices;

    if let Some(hyp) = fit_validate::<PlanarSurface>(mesh, faces, verts, centroids, config) {
        cls.push_planar(hyp);
        wins[0] += 1;
        return true;
    }

    let try_sphere = |cls: &mut Classification, wins: &mut [usize; 4]| {
        if let Some(hyp) = fit_validate::<SphericalSurface>(mesh, faces, verts, centroids, config) {
            cls.push_spherical(hyp);
            wins[1] += 1;
            true
        } else {
            false
        }
    };
    let try_conical = |cls: &mut Classification, wins: &mut [usize; 4]| {
        if let Some(hyp) = fit_validate::<ConicalSurface>(mesh, faces, verts, centroids, config) {
            cls.push_conical(hyp);
            wins[2] += 1;
            true
        } else {
            false
        }
    };

    // Jet-informed rung order: a developable region tries the conical
    // (cylinder/cone) rung before the sphere rung, so a short ring-tessellated
    // cylinder is not claimed by an ambiguous equal-radius sphere. The jet
    // only orders — a genuine sphere the jet misread as developable still
    // falls through to the sphere rung (a cylinder fit cannot validate a
    // sphere's vertices).
    if region_is_developable(region, jet) {
        if try_conical(cls, wins) || try_sphere(cls, wins) {
            return true;
        }
    } else if try_sphere(cls, wins) || try_conical(cls, wins) {
        return true;
    }

    if let Some(hyp) = fit_validate::<ToroidalSurface>(mesh, faces, verts, centroids, config) {
        cls.push_toroidal(hyp);
        wins[3] += 1;
        return true;
    }
    false
}

/// Fit surface type `S` to a region and validate it against the whole region.
/// Large regions fit on a stride subsample, refit-augmented with the worst
/// validation offenders on a near-miss. Returns the committed hypothesis, or
/// `None` if `S` cannot explain the region within tolerance.
///
/// Also the stage-2.4 merge revalidation entry point (`repair.rs`): refit a
/// candidate union and confirm it validates end to end before committing.
pub(super) fn fit_validate<S: SurfaceType>(
    mesh: &ConnectedMesh,
    region_faces: &[usize],
    region_verts: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
) -> Option<S::Hypothesis> {
    let vtol = config.vertex_tolerance_mm;
    let stol = config.surface_tolerance_mm;

    // Initial fit set: the whole region, or a deterministic stride subsample
    // when the region is large.
    let mut fit_faces: Vec<usize> = if region_faces.len() > MAX_FIT_FACES {
        let stride = region_faces.len().div_ceil(MAX_FIT_FACES);
        region_faces.iter().copied().step_by(stride).collect()
    } else {
        region_faces.to_vec()
    };

    for _round in 0..=MAX_REFIT {
        let fit_verts = vertices_of(mesh, &fit_faces);
        let mut fit = super::surface::cached_fit::<S>(mesh, &fit_faces, &fit_verts, config.verbosity)?;
        S::determine_convexity(mesh, region_faces, &mut fit);

        // Validate every region face; collect near-miss offenders for a refit.
        let mut offenders: Vec<(f64, usize)> = Vec::new();
        for &fi in region_faces {
            // Orientation coherence: a face whose normal disagrees with the
            // fitted surface's convexity means the region is not one clean
            // surface of type S — no refit rescues that.
            if !S::convexity_check(&fit, mesh, fi) {
                return None;
            }
            let face = &mesh.faces[fi];
            let vc = face.vertex_count as usize;
            let mut excess = 0.0_f64;
            for &vi in &face.vertex_indices[..vc] {
                let v = &mesh.vertices[vi];
                let d = S::vertex_distance(&fit, &[v.x, v.y, v.z]);
                if d > vtol {
                    excess = excess.max(d - vtol);
                }
            }
            let cd = S::centroid_distance(&fit, &centroids[fi]);
            if cd > stol {
                excess = excess.max(cd - stol);
            }
            if excess > 0.0 {
                offenders.push((excess, fi));
            }
        }

        if offenders.is_empty() {
            return Some(S::make_hypothesis(
                mesh,
                region_faces.to_vec(),
                region_verts.to_vec(),
                &fit,
            ));
        }

        // Near-miss handling: if the worst offender is far out, the region is
        // decisively not this surface. Otherwise fold the worst offenders
        // into the fit set and refit (bounded rounds).
        offenders.sort_by(|a, b| {
            b.0.partial_cmp(&a.0)
                .unwrap()
                .then_with(|| a.1.cmp(&b.1))
        });
        if offenders[0].0 > REFIT_GIVE_UP * stol {
            return None;
        }
        let before = fit_faces.len();
        let existing: HashSet<usize> = fit_faces.iter().copied().collect();
        for &(_, fi) in offenders.iter().take(REFIT_ADD) {
            if !existing.contains(&fi) {
                fit_faces.push(fi);
            }
        }
        fit_faces.sort_unstable();
        fit_faces.dedup();
        if fit_faces.len() == before {
            // No new faces to learn from — cannot improve.
            return None;
        }
    }
    None
}

/// Size of the largest mesh-connected component within a face set.
/// Try the whole-region DOF ladder on each connected component of the
/// leftovers, appending any that validate to `hyps` and returning the faces
/// still unexplained.
///
/// Rungs in ascending DOF, first validating fit wins, exactly as
/// `classify_one_region` does — a leftover component *is* a region that
/// competition failed to explain, so it deserves the same treatment before
/// freeform gets it.
fn rescue_leftover_components(
    mesh: &ConnectedMesh,
    leftovers: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
    hyps: &mut Vec<CompeteHyp>,
) -> Vec<usize> {
    let mut still: Vec<usize> = Vec::new();
    for comp in components_of(mesh, leftovers) {
        if comp.len() < MIN_CYLINDER_FACES {
            still.extend(comp);
            continue;
        }
        let verts = vertices_of(mesh, &comp);
        if let Some(h) = fit_validate::<PlanarSurface>(mesh, &comp, &verts, centroids, config) {
            hyps.push(CompeteHyp::Planar(h));
        } else if let Some(h) = (comp.len() >= MIN_SPHERE_FACES)
            .then(|| fit_validate::<SphericalSurface>(mesh, &comp, &verts, centroids, config))
            .flatten()
        {
            hyps.push(CompeteHyp::Spherical(h));
        } else if let Some(h) =
            fit_validate::<ConicalSurface>(mesh, &comp, &verts, centroids, config)
        {
            hyps.push(CompeteHyp::Conical(h));
        } else if let Some(h) = (comp.len() >= MIN_TORUS_FACES)
            .then(|| fit_validate::<ToroidalSurface>(mesh, &comp, &verts, centroids, config))
            .flatten()
        {
            hyps.push(CompeteHyp::Toroidal(h));
        } else {
            still.extend(comp);
        }
    }
    still.sort_unstable();
    still
}

/// Grow cylinders through the remaining leftovers using axis directions
/// borrowed from adjacent committed surfaces, largest validated area first.
fn rescue_with_neighbour_axes(
    mesh: &ConnectedMesh,
    leftovers: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
    hyps: &mut Vec<CompeteHyp>,
) -> Vec<usize> {
    let mut open: HashSet<usize> = leftovers.iter().copied().collect();
    // face -> committed hypothesis index, for finding a component's neighbours
    let mut owner: Vec<usize> = vec![usize::MAX; mesh.faces.len()];
    for (hi, h) in hyps.iter().enumerate() {
        for &fi in h.faces() {
            owner[fi] = hi;
        }
    }

    loop {
        // Candidate axes: the axis directions of every committed surface
        // adjacent to something still open.
        let mut axes: Vec<[f64; 3]> = Vec::new();
        // In sorted face order, not `HashSet` order: near-parallel axes are
        // merged by keeping the FIRST one seen, so the surviving axis — and the
        // order the axes are tried in — would otherwise vary per process.
        let mut open_sorted: Vec<usize> = open.iter().copied().collect();
        open_sorted.sort_unstable();
        for &fi in &open_sorted {
            let vc = mesh.faces[fi].vertex_count as usize;
            for &ni in &mesh.faces[fi].neighbors[..vc] {
                if ni < 0 {
                    continue;
                }
                let hi = owner[ni as usize];
                if hi == usize::MAX {
                    continue;
                }
                let ax = match &hyps[hi] {
                    CompeteHyp::Conical(c) => c.axis_direction,
                    CompeteHyp::Toroidal(t) => t.axis_direction,
                    _ => continue,
                };
                if !axes.iter().any(|a| {
                    let d = a[0] * ax[0] + a[1] * ax[1] + a[2] * ax[2];
                    d.abs() > 0.999_999
                }) {
                    axes.push(ax);
                }
            }
        }
        if axes.is_empty() {
            break;
        }

        // The arena for this round is what remains open.
        let mut in_open = vec![false; mesh.faces.len()];
        for &fi in open.iter() {
            in_open[fi] = true;
        }
        let mut seeds: Vec<usize> = open.iter().copied().collect();
        seeds.sort_unstable();

        let mut best: Option<(ConicalHypothesis, f64, Vec<usize>)> = None;
        for &seed in &seeds {
            for &axis in &axes {
                if let Some((h, area, faces)) =
                    compete::grow_with_axis(mesh, &in_open, seed, axis, centroids, &compete::CompeteConfig::from(config))
                {
                    if best.as_ref().is_none_or(|(_, ba, _)| area > *ba + 1e-12) {
                        best = Some((h, area, faces));
                    }
                }
            }
        }

        match best {
            Some((h, _, faces)) => {
                let hi = hyps.len();
                for &fi in &faces {
                    open.remove(&fi);
                    owner[fi] = hi;
                }
                hyps.push(CompeteHyp::Conical(h));
            }
            None => break,
        }
    }

    let mut still: Vec<usize> = open.into_iter().collect();
    still.sort_unstable();
    still
}

/// Connected components of a face set, over mesh adjacency.
fn components_of(mesh: &ConnectedMesh, faces: &[usize]) -> Vec<Vec<usize>> {
    let member: HashSet<usize> = faces.iter().copied().collect();
    let mut seen: HashSet<usize> = HashSet::new();
    let mut out: Vec<Vec<usize>> = Vec::new();
    let mut starts = faces.to_vec();
    starts.sort_unstable();
    for &start in &starts {
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
        out.push(comp);
    }
    out
}

/// Dissolve into the per-facet residual any committed surface whose boundary
/// with that residual is more fragmented than the surface itself, and return
/// the indices dissolved.
///
/// This replaces the old whole-region discard, which threw away a region's
/// entire decomposition whenever freeform left a connected leftover blob
/// bigger than a fixed size. The premise it protected is real — a surface
/// abutting a jagged per-facet blob reconstructs to an invalid wire — but the
/// remedy was wildly out of proportion to it, and *size was the wrong
/// quantity*. Measured: `cone_tangent_planes` leaves a blob covering 73% of
/// its region and reconstructs with zero BRepCheck failures, while
/// `bspline_blister` leaves one covering 43% and produces sixteen invalid
/// faces — so no threshold on the blob, absolute or relative, separates them.
///
/// What does separate them is *where the invalidity comes from*: a face's wire
/// gains one edge per distinct per-facet neighbour, so the face that goes
/// invalid is the one whose boundary has been shattered into more pieces than
/// the face has substance. Comparing a surface's per-facet neighbour count
/// against its own face count says exactly that, carries no length and no
/// tuned number, and dissolves only the surfaces that are actually in trouble
/// rather than every surface that happens to share a region with them.
///
/// Iterated to a fixpoint, because dissolving one surface enlarges the
/// residual and can shatter its neighbour's boundary in turn. In the limit
/// this reproduces the old whole-region behaviour — but only when every
/// surface in the region individually earns it.
fn dissolve_fragmented_surfaces(
    mesh: &ConnectedMesh,
    hyps: &[CompeteHyp],
    already_dissolved: &HashSet<usize>,
    residual: &mut Vec<usize>,
) -> HashSet<usize> {
    let mut dissolved: HashSet<usize> = HashSet::new();
    let mut in_residual: HashSet<usize> = residual.iter().copied().collect();
    loop {
        let mut round: Vec<usize> = Vec::new();
        for (i, hyp) in hyps.iter().enumerate() {
            if already_dissolved.contains(&i) || dissolved.contains(&i) {
                continue;
            }
            let own: HashSet<usize> = hyp.faces().iter().copied().collect();
            let mut neighbours: HashSet<usize> = HashSet::new();
            for &fi in hyp.faces() {
                let vc = mesh.faces[fi].vertex_count as usize;
                for &ni in &mesh.faces[fi].neighbors[..vc] {
                    if ni >= 0 {
                        let ni = ni as usize;
                        if in_residual.contains(&ni) && !own.contains(&ni) {
                            neighbours.insert(ni);
                        }
                    }
                }
            }
            if neighbours.len() > hyp.faces().len() {
                round.push(i);
            }
        }
        if round.is_empty() {
            return dissolved;
        }
        for i in round {
            dissolved.insert(i);
            for &fi in hyps[i].faces() {
                if in_residual.insert(fi) {
                    residual.push(fi);
                }
            }
        }
    }
}

fn largest_component(mesh: &ConnectedMesh, faces: &[usize]) -> usize {
    let member: HashSet<usize> = faces.iter().copied().collect();
    let mut seen: HashSet<usize> = HashSet::new();
    let mut largest = 0usize;
    for &start in faces {
        if seen.contains(&start) {
            continue;
        }
        let mut size = 0usize;
        let mut stack = vec![start];
        seen.insert(start);
        while let Some(fi) = stack.pop() {
            size += 1;
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
        largest = largest.max(size);
    }
    largest
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

/// A single-face planar hypothesis (the terminal per-facet fallback).
fn facet_plane(mesh: &ConnectedMesh, fi: usize) -> PlanarHypothesis {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    let vertices: Vec<usize> = face.vertex_indices[..vc].to_vec();
    let vertex_set: HashSet<usize> = vertices.iter().copied().collect();
    let (normal, distance) = fit_plane(&[fi], &vertex_set, &mesh.faces, &mesh.vertices);
    let mut error_max = f64::NEG_INFINITY;
    let mut error_min = f64::INFINITY;
    let mut error_abs_sum = 0.0_f64;
    for &vi in &vertices {
        let v = &mesh.vertices[vi];
        let d = normal[0] * v.x + normal[1] * v.y + normal[2] * v.z - distance;
        error_max = error_max.max(d);
        error_min = error_min.min(d);
        error_abs_sum += d.abs();
    }
    PlanarHypothesis {
        normal,
        distance,
        faces: vec![fi],
        vertices,
        error_max,
        error_min,
        error_abs_sum,
    }
}
