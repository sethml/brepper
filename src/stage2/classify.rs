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
//! **The sphere and conical rungs are decided together, not in DOF order.**
//! A sphere passes exactly through every vertex of a single band of facets
//! between two coaxial circles — a chamfer, a countersink, a short hole — so
//! simplest-first would give all of those to the sphere. See
//! `sphere_or_conical` for the rule and the measurement behind it. (Until
//! 2026-11-16 the jet ordered the two rungs: a region reading developable,
//! by a tuned ratio `DEVELOPABLE_MAX_RATIO = 0.4`, tried the conical rung
//! first. It read chamfers as not developable; the constant is gone.)
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
    /// Per mesh face: the face joined its surface at `--boundary-tolerance`
    /// ([`join_boundary_facets`]) rather than validating at vertex tolerance.
    /// Such a face is a member for coverage but never evidence for a fit, and
    /// anything that revalidates it holds it to the boundary bound. Empty means
    /// no face joined.
    pub joined: Vec<bool>,
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
            joined: Vec::new(),
        }
    }

    /// Whether face `fi` joined its surface at the boundary bound.
    pub(super) fn is_joined(&self, fi: usize) -> bool {
        self.joined.get(fi).copied().unwrap_or(false)
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
    if config.boundary_join() {
        cls.joined = vec![false; mesh.faces.len()];
    }
    let mut joined_faces = 0usize;
    let mut ladder_wins = [0usize; 4]; // plane, sphere, conical, torus
    let mut competed_regions = 0usize;
    let mut competed_surfaces = 0usize;
    let mut freeform_surfaces = 0usize;
    let mut facet_faces = 0usize;
    for region in regions {
        if classify_one_region(mesh, region, &centroids, config, &mut cls, &mut ladder_wins) {
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
        //
        // Each round first lets boundary facets join the surfaces already
        // found ([`join_boundary_facets`]), so that the rescue below never sees
        // a row of them: a ring of facets between an edge circle and a surface
        // circle is exactly a cone frustum, and the ladder would commit it.
        let mut leftovers = leftovers;
        loop {
            let before = leftovers.len();
            leftovers = join_boundary_facets(
                mesh, &leftovers, &centroids, config, &mut hyps, &mut cls.joined,
            );
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
                joined_faces += hyp.faces().iter().filter(|&&fi| cls.is_joined(fi)).count();
                let hyp = match hyp {
                    CompeteHyp::Toroidal(torus) => {
                        retype_torus(mesh, torus, &centroids, config, &cls.joined)
                    }
                    _ => retype_sphere(mesh, hyp, &centroids, config, &cls.joined),
                };
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
        if config.boundary_join() {
            eprintln!(
                "  boundary join: {joined_faces} facet(s) joined an adjacent surface at --boundary-tolerance {:.1e}",
                config.boundary_tolerance(),
            );
        }
    }
    cls
}

/// Stage 2.3(c2): let a leftover facet at a face boundary join the adjacent
/// surface it lies on, to within `--boundary-tolerance`.
///
/// Some exporters write the vertices of an edge on their approximation of the
/// edge curve rather than on the two surfaces. Measured on six Fusion 360
/// parts: 18-72 vertices per part 1.7e-4 to 5.2e-4 mm off every untrimmed
/// reference surface, and within 7e-6 mm of the trimmed faces. A facet holding
/// one is on no surface at `--vertex-tolerance`, so competition leaves it, and
/// the rescue and freeform stages turned those rows into cone and torus
/// impostors, single-facet planes and small B-spline patches.
///
/// Discovery stays at the strict tolerance: raising `--vertex-tolerance` to
/// 1e-3 instead fixed these parts and broke exact ones (measured 2026-11-21:
/// impostor tori and spheres, a pipe elbow from good to invalid). Here a facet
/// may only join a surface that was already found strictly, and only if it
/// shares an edge with one of its faces, passes the same per-face test a member
/// passes (orientation, centroid within `--surface-tolerance`) with the vertex
/// bound raised to `--boundary-tolerance`, and fits no adjacent surface better.
/// Rounds repeat until none joins, since a facet whose vertices are all on
/// edges (a corner, a narrow strip) may only touch the surface through another
/// joined facet. The fit is not changed: it stays the one made from facets that
/// validate strictly, and `joined` marks the facets so stage 2.4 never refits
/// from them either.
///
/// Only a surface its own strict facets determine can take joined facets: at
/// least [`CompeteHyp::min_support`] of them. Competition's claim resolution
/// can leave a surface with a single facet, and such a fit extrapolates
/// anywhere — measured on `ToolChanger__blower_fan_duct`: a one-facet
/// "cylinder" took 69 facets of a freeform face at 1e-3, and the freeform patch
/// they came from no longer fitted.
fn join_boundary_facets(
    mesh: &ConnectedMesh,
    leftovers: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
    hyps: &mut [CompeteHyp],
    joined: &mut [bool],
) -> Vec<usize> {
    if !config.boundary_join() || leftovers.is_empty() || hyps.is_empty() {
        return leftovers.to_vec();
    }
    let btol = config.boundary_tolerance();
    let stol = config.surface_tolerance_mm;
    let is_joined = |fi: usize| joined.get(fi).copied().unwrap_or(false);
    // Whether a surface's strict facets determine it (see above).
    let determined: Vec<bool> = hyps
        .iter()
        .map(|h| h.faces().iter().filter(|&&fi| !is_joined(fi)).count() >= h.min_support())
        .collect();
    let mut owner = vec![usize::MAX; mesh.faces.len()];
    for (hi, h) in hyps.iter().enumerate() {
        for &fi in h.faces() {
            owner[fi] = hi;
        }
    }
    let fits: Vec<compete::HypFit> = hyps.iter().map(|h| compete::HypFit::of(mesh, h)).collect();
    let mut open = vec![false; mesh.faces.len()];
    for &fi in leftovers {
        open[fi] = true;
    }
    let mut added: Vec<Vec<usize>> = vec![Vec::new(); hyps.len()];
    loop {
        // Every open facet is judged against the surfaces as they stood at the
        // start of the round, so the outcome does not depend on visiting order.
        let mut round: Vec<(usize, usize)> = Vec::new();
        for &fi in leftovers {
            if !open[fi] {
                continue;
            }
            let f = &mesh.faces[fi];
            let mut cands: Vec<usize> = f.neighbors[..f.vertex_count as usize]
                .iter()
                .filter(|&&ni| ni >= 0)
                .map(|&ni| owner[ni as usize])
                .filter(|&hi| hi != usize::MAX && determined[hi])
                .collect();
            cands.sort_unstable();
            cands.dedup();
            let best = cands
                .into_iter()
                .filter_map(|hi| fits[hi].face_fits(mesh, fi, centroids, btol, stol).map(|d| (d, hi)))
                .min_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            if let Some((_, hi)) = best {
                round.push((fi, hi));
            }
        }
        if round.is_empty() {
            break;
        }
        for (fi, hi) in round {
            open[fi] = false;
            owner[fi] = hi;
            added[hi].push(fi);
        }
    }
    // A vertex more than vertex tolerance off the surface must be on a real
    // edge of it — the case the bound exists for, a vertex on an approximated
    // edge curve between two CAD faces. It is on one when a facet around it is
    // across a sharp edge (outside this region: neither a surface's nor a
    // leftover here) or belongs to another surface its own facets determine. A
    // leftover facet does not count. Bordering only leftovers is what a
    // near-tangent blend looks like: measured on `ToolChanger__blower_fan_duct`,
    // counting leftovers let 244 blend facets join analytic surfaces, the
    // freeform patch they belonged to stopped fitting, and per-facet planes
    // went from 767 to 2039. Facets that joined the other surface count: on a
    // tangent edge between two surfaces both sides' edge facets straddle. A
    // joined facet with an off vertex on no real edge goes back to the
    // leftovers, repeated until none does, since sending one back can take
    // away another's evidence.
    let vtol = config.vertex_tolerance_mm;
    loop {
        let mut evict: Vec<(usize, usize)> = Vec::new();
        for (hi, extra) in added.iter().enumerate() {
            for &fi in extra {
                let f = &mesh.faces[fi];
                let loose = f.vertex_indices[..f.vertex_count as usize].iter().any(|&vi| {
                    let v = &mesh.vertices[vi];
                    fits[hi].distance(&[v.x, v.y, v.z]) > vtol
                        && !super::vertex_fan(mesh, fi, vi).iter().any(|&g| {
                            let o = owner[g];
                            (o == usize::MAX && !open[g]) || (o != usize::MAX && o != hi && determined[o])
                        })
                });
                if loose {
                    evict.push((fi, hi));
                }
            }
        }
        if evict.is_empty() {
            break;
        }
        for (fi, hi) in evict {
            open[fi] = true;
            owner[fi] = usize::MAX;
            added[hi].retain(|&g| g != fi);
        }
    }
    if std::env::var("BREPPER_DEBUG_JOIN").is_ok() {
        for (hi, h) in hyps.iter().enumerate() {
            eprintln!(
                "    join: hyp {hi} {} {} faces (+{} joined)",
                h.kind_name(),
                h.faces().len(),
                added[hi].len()
            );
        }
        for &fi in leftovers.iter().filter(|&&fi| open[fi]) {
            let f = &mesh.faces[fi];
            let near: Vec<String> = f.neighbors[..f.vertex_count as usize]
                .iter()
                .filter(|&&ni| ni >= 0 && owner[ni as usize] != usize::MAX)
                .map(|&ni| {
                    let hi = owner[ni as usize];
                    let d = fits[hi].face_fits(mesh, fi, centroids, f64::INFINITY, f64::INFINITY);
                    let ok = fits[hi].face_fits(mesh, fi, centroids, f64::INFINITY, stol).is_some();
                    format!("hyp {hi} {} dmax {:?} centroid+orientation ok {ok}", hyps[hi].kind_name(), d)
                })
                .collect();
            eprintln!("    join: facet {fi} stays open; adjacent: {}", near.join("; "));
        }
    }
    for (hi, extra) in added.iter().enumerate() {
        if extra.is_empty() {
            continue;
        }
        for &fi in extra {
            joined[fi] = true;
        }
        hyps[hi].add_faces(mesh, extra);
    }
    leftovers.iter().copied().filter(|&fi| open[fi]).collect()
}

/// The facets just outside a face set — sharing an edge with it — that a fit
/// of the set also validates, by the same test a member facet must pass.
fn validated_past<S: SurfaceType>(
    mesh: &ConnectedMesh,
    faces: &[usize],
    fit: &S::Fit,
    centroids: &[[f64; 3]],
    config: &Config,
) -> Vec<usize> {
    let inside: std::collections::HashSet<usize> = faces.iter().copied().collect();
    let mut outside: Vec<usize> = Vec::new();
    for &fi in faces {
        let f = &mesh.faces[fi];
        for &nb in &f.neighbors[..f.vertex_count as usize] {
            if nb >= 0 && !inside.contains(&(nb as usize)) {
                outside.push(nb as usize);
            }
        }
    }
    outside.sort_unstable();
    outside.dedup();
    outside.retain(|&fi| {
        super::compete::face_validates::<S>(
            mesh,
            fi,
            fit,
            centroids,
            config.vertex_tolerance_mm,
            config.surface_tolerance_mm,
        )
    });
    outside
}

/// Does the sphere go on where no cylinder or cone can follow?
///
/// A single band is the same mesh whether it was cut from a cone or from a
/// sphere, so the band cannot say which — but its neighbours can. A sphere
/// tessellated coarsely enough that every latitude row is its own sharp-edged
/// region presents each row as an exact cone frustum, ruling edges and all.
/// What gives it away is the NEXT row: the sphere passes through it too, and
/// no cone passes through three circles of a sphere.
///
/// So: take the facets just past `faces` that the sphere validates and that
/// reach a vertex `faces` does not have, and ask
/// whether a sphere still explains `faces` plus those while a cylinder-or-cone
/// does not. Both are asked by REFITTING on the larger set, never by
/// extrapolating the fits of `faces`: a cone fitted to part of a band does not
/// extrapolate around the rest of it to vertex tolerance, where a sphere's
/// four parameters do, and comparing extrapolations read every partial band as
/// a sphere (measured on `Open-Task-Light__led_cover_-_top`).
fn sphere_continues(
    mesh: &ConnectedMesh,
    faces: &[usize],
    sphere_fit: &<SphericalSurface as SurfaceType>::Fit,
    centroids: &[[f64; 3]],
    config: &Config,
) -> bool {
    // Only neighbours that bring a vertex the set does not already have. A
    // facet whose vertices are all on the band's own boundary circle — a
    // planar cap, the flat around a hole — "validates" against the sphere
    // through that circle without being evidence of anything, and with it in
    // the union the cylinder is refused for its orientation, not its shape
    // (measured: four short holes on `Open-Task-Light__led_cover_-_top` turned
    // into spheres that way).
    let own: std::collections::HashSet<usize> = vertices_of(mesh, faces).into_iter().collect();
    let mut past = validated_past::<SphericalSurface>(mesh, faces, sphere_fit, centroids, config);
    past.retain(|&fi| {
        let f = &mesh.faces[fi];
        f.vertex_indices[..f.vertex_count as usize].iter().any(|vi| !own.contains(vi))
    });
    if past.is_empty() {
        return false;
    }
    let mut union: Vec<usize> = faces.iter().copied().chain(past).collect();
    union.sort_unstable();
    union.dedup();
    let verts = vertices_of(mesh, &union);
    fit_validate::<SphericalSurface>(mesh, &union, &verts, centroids, config).is_some()
        && fit_validate::<ConicalSurface>(mesh, &union, &verts, centroids, config).is_none()
}

/// What the sphere and conical rungs of the ladder settle on.
enum SphereOrConical {
    Sphere(SphericalHypothesis),
    Conical(ConicalHypothesis),
}

/// The sphere rung and the cylinder-or-cone rung, decided together.
///
/// These two cannot simply be tried in DOF order, because a sphere passes
/// exactly through every vertex of a single band of facets between two
/// coaxial circles — which is how a chamfer, a countersink and a short hole
/// are all tessellated. Measured over the corpus (2026-11-16, 151 models): a
/// sphere validates the true facets of 180 of 223 reference cones and 136 of
/// 858 reference cylinders, and before this rule every chamfer on the wild
/// corpus came out as a sphere.
///
/// So: fit both. If only one validates, that is the answer. If both do, and
/// the sphere carries on into neighbouring facets where no cylinder or cone
/// can follow ([`sphere_continues`]), it is a sphere. Otherwise — the usual
/// case, a band between a plane and a hole — the
/// one that contains more of the facets' edge length wins
/// ([`edges_on_surface`]): a cone or cylinder contains its ruling edges whole,
/// a sphere contains no edge at all. On the same 316 ambiguous faces that
/// picks the reference's type 316 times, with no ties; and no cylinder or cone
/// validates the facets of any of the 143 reference spheres, so the rule is
/// never asked about a real one.
///
/// A tie — no edge on either, as on a cone triangulated across its rulings —
/// goes to the cylinder-or-cone. That is a prior, not a deduction: the band is
/// then genuinely the same mesh either way. It is the prior the corpus
/// supports (every one of the 316 faces that both surfaces validate is a cone
/// or a cylinder), and a sphere that is more than one band never reaches it.
/// The jet's curvature reading used to break this tie and got it wrong on
/// exactly such a cone (`Open-Task-Light__led_cover_-_top`, 167 facets).
fn sphere_or_conical(
    mesh: &ConnectedMesh,
    faces: &[usize],
    verts: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
    try_sphere: bool,
) -> Option<SphereOrConical> {
    let sphere = try_sphere
        .then(|| fit_validate_with_fit::<SphericalSurface>(mesh, faces, verts, centroids, config))
        .flatten();
    let conical = fit_validate_with_fit::<ConicalSurface>(mesh, faces, verts, centroids, config);
    match (sphere, conical) {
        (None, None) => None,
        (Some((s, _)), None) => Some(SphereOrConical::Sphere(s)),
        (None, Some((c, _))) => Some(SphereOrConical::Conical(c)),
        (Some((s, s_fit)), Some((c, c_fit))) => {
            let vtol = config.vertex_tolerance_mm;
            // First: does the sphere carry on into neighbouring facets that no
            // cylinder or cone can also explain? Then it is a sphere, whatever
            // this one band looks like.
            if sphere_continues(mesh, faces, &s_fit, centroids, config) {
                if std::env::var("BREPPER_DEBUG_COMPETE").is_ok() {
                    eprintln!("    sphere_or_conical: {} faces — sphere continues past them", faces.len());
                }
                return Some(SphereOrConical::Sphere(s));
            }
            if std::env::var("BREPPER_DEBUG_COMPETE").is_ok() {
                let es = edges_on_surface::<SphericalSurface>(mesh, faces, &s_fit, vtol);
                let ec = edges_on_surface::<ConicalSurface>(mesh, faces, &c_fit, vtol);
                eprintln!(
                    "    sphere_or_conical: {} faces — edges on sphere {:.3} mm ({}), on conical {:.3} mm ({}), of {:.3} mm",
                    faces.len(), es.0, es.2, ec.0, ec.2, ec.1,
                );
            }
            let on_sphere = edges_on_surface::<SphericalSurface>(mesh, faces, &s_fit, vtol).0;
            let on_conical = edges_on_surface::<ConicalSurface>(mesh, faces, &c_fit, vtol).0;
            if on_conical >= on_sphere {
                Some(SphereOrConical::Conical(c))
            } else {
                Some(SphereOrConical::Sphere(s))
            }
        }
    }
}

/// Ask the ladder's sphere-or-conical question of a sphere that competition
/// committed, on the faces it ended up with.
///
/// Competition's per-seed rule is "largest validated area", and the edge
/// evidence is only its tie-break. That is not enough: from a seed on a cone,
/// the cone candidate is fitted to a small patch and often grows less far than
/// the sphere candidate does (the cone is the weak fitter at small support), so
/// the sphere wins on area and never reaches the tie-break. Measured on
/// `Open-Task-Light__led_cover_-_top`: with the tie-break alone, three pieces of
/// two reference cones (94-104 mm^2 each) still came out as spheres.
///
/// Once the faces are settled the question can be asked properly: fit a
/// cylinder-or-cone to exactly those faces and let `sphere_or_conical` decide.
/// A sphere that no cylinder or cone validates is returned untouched.
///
/// Facets that joined the sphere at `--boundary-tolerance` are left out of the
/// question — they would fail any strict refit — and a cylinder-or-cone verdict
/// is taken only if it also holds them at that bound.
fn retype_sphere(
    mesh: &ConnectedMesh,
    hyp: CompeteHyp,
    centroids: &[[f64; 3]],
    config: &Config,
    joined: &[bool],
) -> CompeteHyp {
    let CompeteHyp::Spherical(sphere) = hyp else {
        return hyp;
    };
    let is_joined = |fi: usize| joined.get(fi).copied().unwrap_or(false);
    let strict: Vec<usize> = sphere.faces.iter().copied().filter(|&fi| !is_joined(fi)).collect();
    if strict.len() < sphere.faces.len() {
        let verts = vertices_of(mesh, &strict);
        let Some(SphereOrConical::Conical(mut c)) =
            sphere_or_conical(mesh, &strict, &verts, centroids, config, true)
        else {
            return CompeteHyp::Spherical(sphere);
        };
        let fit = compete::HypFit::of(mesh, &CompeteHyp::Conical(c.clone()));
        let (btol, stol) = (config.boundary_tolerance(), config.surface_tolerance_mm);
        let holds = sphere
            .faces
            .iter()
            .filter(|&&fi| is_joined(fi))
            .all(|&fi| fit.face_fits(mesh, fi, centroids, btol, stol).is_some());
        if !holds {
            return CompeteHyp::Spherical(sphere);
        }
        c.faces = sphere.faces;
        c.vertices = sphere.vertices;
        return CompeteHyp::Conical(c);
    }
    let verdict =
        sphere_or_conical(mesh, &sphere.faces, &sphere.vertices, centroids, config, true);
    if std::env::var("BREPPER_DEBUG_COMPETE").is_ok() {
        let conical = fit_validate::<ConicalSurface>(
            mesh, &sphere.faces, &sphere.vertices, centroids, config,
        );
        eprintln!(
            "    retype: sphere r={:.3} over {} faces — conical validates: {}, verdict: {}",
            sphere.radius,
            sphere.faces.len(),
            conical.is_some(),
            match &verdict {
                Some(SphereOrConical::Conical(_)) => "conical",
                Some(SphereOrConical::Sphere(_)) => "sphere",
                None => "neither refits",
            },
        );
    }
    match verdict {
        Some(SphereOrConical::Conical(c)) => CompeteHyp::Conical(c),
        _ => CompeteHyp::Spherical(sphere),
    }
}

/// Refit a torus's settled support as a cylinder or cone. A torus can pass
/// through both vertex rings of a chamfer. The chamfer's ruling edges supply
/// evidence between those vertices, using the same test as competition.
fn retype_torus(
    mesh: &ConnectedMesh,
    torus: ToroidalHypothesis,
    centroids: &[[f64; 3]],
    config: &Config,
    joined: &[bool],
) -> CompeteHyp {
    let Some(conical) = fit_validate_joined::<ConicalSurface>(
        mesh, &torus.faces, &torus.vertices, centroids, config, joined,
    ) else {
        return CompeteHyp::Toroidal(torus);
    };
    let strict: Vec<_> = torus.faces.iter().copied()
        .filter(|&fi| !joined.get(fi).copied().unwrap_or(false))
        .collect();
    let conical = CompeteHyp::Conical(conical);
    let torus = CompeteHyp::Toroidal(torus);
    let vtol = config.vertex_tolerance_mm;
    let (conical_length, _) = compete::edges_on_hyp(mesh, &conical, &strict, vtol);
    let (torus_length, _) = compete::edges_on_hyp(mesh, &torus, &strict, vtol);
    if config.verbosity >= 3 {
        eprintln!(
            "    retype: torus over {} faces — edges on conical {conical_length:.3} mm, on torus {torus_length:.3} mm",
            torus.faces().len(),
        );
    }
    if conical_length > torus_length { conical } else { torus }
}

/// Try the DOF ladder on one region; commit the first surface that validates.
/// Returns true on a committed surface, false if the region fell through
/// (caller applies the per-facet fallback).
fn classify_one_region(
    mesh: &ConnectedMesh,
    region: &Region,
    centroids: &[[f64; 3]],
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

    // The sphere and conical rungs together: where both validate, the mesh's
    // own edges decide (`sphere_or_conical`).
    match sphere_or_conical(mesh, faces, verts, centroids, config, true) {
        Some(SphereOrConical::Sphere(hyp)) => {
            cls.push_spherical(hyp);
            wins[1] += 1;
            return true;
        }
        Some(SphereOrConical::Conical(hyp)) => {
            cls.push_conical(hyp);
            wins[2] += 1;
            return true;
        }
        None => {}
    }

    if let Some(hyp) = fit_validate::<ToroidalSurface>(mesh, faces, verts, centroids, config) {
        cls.push_toroidal(hyp);
        wins[3] += 1;
        return true;
    }
    false
}

/// How much of a face set's mesh EDGES lies on a fitted surface: the total
/// length of the edges whose midpoint is within `vtol` of it, the total length
/// of all the set's edges, and the number on the surface.
///
/// This is evidence a vertex test cannot give. Every vertex of a validated
/// support is on the surface by definition, for the true surface and for an
/// impostor alike. But a tessellator that meshes a cylinder or a cone lays
/// facet edges along its rulings, and a ruling is a straight line IN the
/// surface: the whole edge lies on it, midpoint included, to the same
/// precision as its ends. No straight edge lies in a sphere; the midpoint of a
/// chord of length `L` is `L^2 / 8R` inside it. So where a sphere and a cone
/// both pass through every vertex of a band of quads — which they always can,
/// for a single band between two circles — the cone contains the band's
/// ruling edges and the sphere contains none of its edges.
///
/// The bar is the vertex tolerance itself: an edge is on the surface exactly
/// when its midpoint validates as a vertex would. No new constant.
pub(super) fn edges_on_surface<S: SurfaceType>(
    mesh: &ConnectedMesh,
    faces: &[usize],
    fit: &S::Fit,
    vtol: f64,
) -> (f64, f64, usize) {
    let mut edges: Vec<(usize, usize)> = Vec::with_capacity(faces.len() * 3);
    for &fi in faces {
        let f = &mesh.faces[fi];
        let n = f.vertex_count as usize;
        for k in 0..n {
            let (a, b) = (f.vertex_indices[k], f.vertex_indices[(k + 1) % n]);
            edges.push((a.min(b), a.max(b)));
        }
    }
    edges.sort_unstable();
    edges.dedup();
    let (mut on, mut total, mut count) = (0.0_f64, 0.0_f64, 0usize);
    for (a, b) in edges {
        let (p, q) = (&mesh.vertices[a], &mesh.vertices[b]);
        let len = ((p.x - q.x).powi(2) + (p.y - q.y).powi(2) + (p.z - q.z).powi(2)).sqrt();
        total += len;
        let mid = [(p.x + q.x) / 2.0, (p.y + q.y) / 2.0, (p.z + q.z) / 2.0];
        if S::vertex_distance(fit, &mid) <= vtol {
            on += len;
            count += 1;
        }
    }
    (on, total, count)
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
    fit_validate_with_fit::<S>(mesh, region_faces, region_verts, centroids, config).map(|r| r.0)
}

/// [`fit_validate`], also returning the fit the hypothesis was made from, for
/// callers that go on to measure against the surface ([`edges_on_surface`]).
pub(super) fn fit_validate_with_fit<S: SurfaceType>(
    mesh: &ConnectedMesh,
    region_faces: &[usize],
    region_verts: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
) -> Option<(S::Hypothesis, S::Fit)> {
    fit_validate_impl::<S>(mesh, region_faces, region_verts, centroids, config, &[])
}

/// [`fit_validate`] for a face set that may hold facets that joined their
/// surface at `--boundary-tolerance` (`joined`, per mesh face): those are left
/// out of the fit and validated at the boundary bound; every other face is held
/// to vertex tolerance as usual. The stage-2.4 refit entry point.
pub(super) fn fit_validate_joined<S: SurfaceType>(
    mesh: &ConnectedMesh,
    region_faces: &[usize],
    region_verts: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
    joined: &[bool],
) -> Option<S::Hypothesis> {
    fit_validate_impl::<S>(mesh, region_faces, region_verts, centroids, config, joined)
        .map(|r| r.0)
}

fn fit_validate_impl<S: SurfaceType>(
    mesh: &ConnectedMesh,
    region_faces: &[usize],
    region_verts: &[usize],
    centroids: &[[f64; 3]],
    config: &Config,
    joined: &[bool],
) -> Option<(S::Hypothesis, S::Fit)> {
    let stol = config.surface_tolerance_mm;
    let is_joined = |fi: usize| joined.get(fi).copied().unwrap_or(false);
    let vtol_of = |fi: usize| {
        if is_joined(fi) {
            config.boundary_tolerance()
        } else {
            config.vertex_tolerance_mm
        }
    };

    // Initial fit set: the whole region (less any joined facets, unless that
    // leaves nothing), or a deterministic stride subsample when it is large.
    let strict: Vec<usize> = region_faces.iter().copied().filter(|&fi| !is_joined(fi)).collect();
    let support: &[usize] = if strict.is_empty() { region_faces } else { &strict };
    let mut fit_faces: Vec<usize> = if support.len() > MAX_FIT_FACES {
        let stride = support.len().div_ceil(MAX_FIT_FACES);
        support.iter().copied().step_by(stride).collect()
    } else {
        support.to_vec()
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
            let vtol = vtol_of(fi);
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
            let hyp =
                S::make_hypothesis(mesh, region_faces.to_vec(), region_verts.to_vec(), &fit);
            return Some((hyp, fit));
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
        // A joined facet is never evidence for the fit, even as an offender.
        for &(_, fi) in offenders.iter().filter(|o| !is_joined(o.1)).take(REFIT_ADD) {
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
        } else if let Some(h) = sphere_or_conical(
            mesh,
            &comp,
            &verts,
            centroids,
            config,
            comp.len() >= MIN_SPHERE_FACES,
        ) {
            hyps.push(match h {
                SphereOrConical::Sphere(h) => CompeteHyp::Spherical(h),
                SphereOrConical::Conical(h) => CompeteHyp::Conical(h),
            });
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
