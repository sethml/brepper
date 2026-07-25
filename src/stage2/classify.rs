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
    ToroidalHypothesis,
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

/// Largest connected leftover blob competition may leave before its result is
/// discarded in favour of a clean per-facet polyhedron (see the fallback in
/// `classify_regions`). A junction sliver is 1–2 faces; a whole undecomposed
/// sub-surface is many more.
const LEFTOVER_BLOB_MAX: usize = 8;

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
        let (hyps, leftovers) = compete::compete_region(mesh, region, jet, &centroids, config);

        // Stage 2.3(e): freeform-fit the leftover components competition
        // could not explain analytically, dissolving any small analytic
        // island enclosed by a validating B-spline patch.
        let ff = bspline::fit_leftover_components(mesh, &leftovers, &hyps, face_region, config);

        // A committed analytic surface islanded among a large jagged
        // per-facet blob reconstructs to invalid wires where they abut —
        // worse than a clean polyhedron. When an analytic hypothesis
        // survives AND freeform could not clear a big leftover blob, discard
        // competition and per-facet the whole region. (Healthy competition
        // leaves only tiny scattered slivers, and freeform absorbs genuine
        // freeform patches, so this fires only on undecomposable coarse
        // fillets.)
        let surviving_hyps = hyps.len() - ff.dissolved.len();
        if config.verbosity >= 3 || std::env::var("BREPPER_DEBUG_COMPETE").is_ok() {
            eprintln!(
                "    compete: region {} faces → {surviving_hyps} surviving hyps, {} leftovers, freeform residual {} (largest blob {}, discard above {LEFTOVER_BLOB_MAX})",
                region.faces.len(),
                leftovers.len(),
                ff.residual.len(),
                largest_component(mesh, &ff.residual),
            );
        }
        if surviving_hyps > 0 && largest_component(mesh, &ff.residual) > LEFTOVER_BLOB_MAX {
            for &fi in &region.faces {
                cls.push_planar(facet_plane(mesh, fi));
            }
            facet_faces += region.faces.len();
            continue;
        }

        // Commit surviving analytic hypotheses (skipping dissolved islands),
        // the freeform surfaces, and the residual leftovers as per-facet
        // planes (total coverage can never fail; a triangle is exactly a
        // plane).
        for (i, hyp) in hyps.into_iter().enumerate() {
            if !ff.dissolved.contains(&i) {
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
        facet_faces += ff.residual.len();
        for fi in ff.residual {
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
        let mut fit = S::fit_surface(mesh, &fit_faces, &fit_verts, config.verbosity)?;
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
