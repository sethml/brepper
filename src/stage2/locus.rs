//! Where two committed surfaces meet, as a side test (DEVELOPMENT_PLAN §2.4).
//!
//! The seam assumption says a B-rep face boundary is always a chain of mesh
//! edges, tangent junctions included. So once two surfaces are committed, the
//! *correct* boundary between them exists in the mesh already and the job is
//! to select it, not to search for it. This module supplies the selection
//! criterion: a scalar field `s(p)` that is zero on the junction, negative on
//! one surface's side and positive on the other. `snap_boundaries`
//! (`repair.rs`) then assigns each interface face by the sign of `s`.
//!
//! # Why not the intersection curve
//!
//! At a *tangent* junction the two surfaces graze rather than cross, and their
//! intersection is catastrophically ill-conditioned: a plane offset by δ from
//! the torus it should graze meets it in a circle displaced by √(2rδ), so a
//! 4.7e-4 mm plane-fit error became the recorded 4.4e-2 mm edge error. That
//! amplification is a property of the intersection, not of the junction.
//!
//! The locus of *minimum distance* has no such amplification. It is the set
//! of points where the two surfaces' normal fields agree, and it depends on
//! the surfaces' directions and radii but **not** on the offset that is
//! amplified: displacing the plane along its own normal leaves the tangency
//! circle exactly where it was. Only a direction error moves it, and that
//! moves it linearly. Positions depart quadratically from a tangency while
//! normals turn linearly, which is the whole reason to key the boundary off
//! the latter.
//!
//! # The meridian reduction
//!
//! Nearly every tangent junction in CAD joins two surfaces of revolution
//! about a *common* axis (a fillet torus against the plane it dies into,
//! against the coaxial cylinder it continues, against a coaxial cone; a
//! sphere against the cylinder it caps). Reduce such a pair to the meridian
//! half-plane (ρ, h) of the shared axis, where each surface becomes a line
//! (plane ⊥ axis, cylinder, cone) or a circle (sphere, torus), and two
//! coaxial surfaces are tangent exactly when their meridian curves are. The
//! junction is then the circle swept by a single meridian point q, and the
//! side coordinate is displacement along the curves' common tangent t̂ there:
//!
//! ```text
//! s(p) = ([ρ(p), h(p)] − q) · t̂
//! ```
//!
//! This subsumes the cases that would otherwise be written out one by one:
//! for a plane ⊥ the torus axis the tangency sits at tube angle ±90°, where
//! ρ = R, and `s` reduces to ρ − R (the tangency circle is the major-radius
//! circle); for a coaxial cylinder it sits at tube angle 0 or 180°, where the
//! axial coordinate is the non-degenerate one, and `s` reduces to h − h_c.
//! Note the complementarity — at each junction one of ρ, h is stationary and
//! useless and the other varies linearly — which is why the common tangent,
//! rather than either coordinate, is the thing to project onto.
//!
//! The one common tangency that is *not* axisymmetric is a plane meeting a
//! cylinder along a ruling (a flat milled on a shaft), handled separately by
//! `Locus::Planar`. Everything else — cone-plane ruling tangency, non-coaxial
//! pairs, point tangencies — is deliberately unhandled and reported loudly by
//! the caller rather than guessed at.

use crate::stage1::{ConnectedMesh, MeshVertex};

use super::{
    ConicalHypothesis, PlanarHypothesis, SelectedSurface, SphericalHypothesis, ToroidalHypothesis,
};

/// Two fitted surfaces count as sharing an axis when their axis directions are
/// parallel to within this many radians and their axis lines are separated by
/// no more than this fraction of the pair's characteristic radius. A *validity*
/// precondition on the meridian reduction — a pair that fails it has no
/// meridian plane, so no locus can be built — not a decision threshold: real
/// coaxial features agree to fit accuracy (~1e-6 relative), and features that
/// are not coaxial disagree by orders of magnitude more.
const COAXIAL_TOL: f64 = 1e-3;

/// Surface normals this close (in radians) at every sampled point of a shared
/// boundary make the junction tangent. Same criterion and same value as stage
/// 3.2's tangency detection, which decides whether stage 3 may relax an edge
/// to the chordal tolerance; the two stages must agree about which junctions
/// are tangent, so they share one definition.
pub(crate) const TANGENCY_TOL_RAD: f64 = 2.0 * std::f64::consts::PI / 180.0;

// ---------------------------------------------------------------------------
// A view of one committed surface's fitted geometry
// ---------------------------------------------------------------------------

/// The fitted geometry of one committed surface, borrowed from the typed
/// hypothesis vectors. B-splines carry no closed-form axis or normal here and
/// are never part of a locus.
#[derive(Clone, Copy)]
pub(crate) enum Surf<'a> {
    Planar(&'a PlanarHypothesis),
    Conical(&'a ConicalHypothesis),
    Spherical(&'a SphericalHypothesis),
    Toroidal(&'a ToroidalHypothesis),
    BSpline,
}

impl<'a> Surf<'a> {
    pub(crate) fn of(
        sel: &SelectedSurface,
        planar: &'a [PlanarHypothesis],
        conical: &'a [ConicalHypothesis],
        spherical: &'a [SphericalHypothesis],
        toroidal: &'a [ToroidalHypothesis],
    ) -> Surf<'a> {
        match *sel {
            SelectedSurface::Planar(i) => Surf::Planar(&planar[i]),
            SelectedSurface::Conical(i) => Surf::Conical(&conical[i]),
            SelectedSurface::Spherical(i) => Surf::Spherical(&spherical[i]),
            SelectedSurface::Toroidal(i) => Surf::Toroidal(&toroidal[i]),
            SelectedSurface::BSpline(_) => Surf::BSpline,
        }
    }

    /// Short type name, for diagnostics.
    pub(crate) fn kind(&self) -> &'static str {
        match self {
            Surf::Planar(_) => "plane",
            Surf::Conical(h) if h.is_cylinder() => "cyl",
            Surf::Conical(_) => "cone",
            Surf::Spherical(_) => "sphere",
            Surf::Toroidal(_) => "torus",
            Surf::BSpline => "bspline",
        }
    }

    /// The surface's axis of revolution as (point, unit direction), for the
    /// types that have exactly one. A plane is a surface of revolution about
    /// any axis parallel to its normal and a sphere about any axis through its
    /// centre, so neither pins one down; both take the partner's.
    fn axis(&self) -> Option<([f64; 3], [f64; 3])> {
        match self {
            Surf::Conical(h) => Some((h.apex, h.axis_direction)),
            Surf::Toroidal(h) => Some((h.center, h.axis_direction)),
            _ => None,
        }
    }

    /// Outward unit normal of the fitted surface at `p`, or None where the
    /// normal is undefined (on a cone's axis) or unavailable (B-spline).
    pub(crate) fn normal_at(&self, p: &[f64; 3]) -> Option<[f64; 3]> {
        match self {
            Surf::Planar(h) => Some(h.normal),
            Surf::Conical(h) => {
                let (rho, radial, _) = axis_frame(p, &h.apex, &h.axis_direction);
                if rho < 1e-15 {
                    return None;
                }
                let (ca, sa) = (h.half_angle.cos(), h.half_angle.sin());
                let n = [
                    radial[0] * ca - h.axis_direction[0] * sa,
                    radial[1] * ca - h.axis_direction[1] * sa,
                    radial[2] * ca - h.axis_direction[2] * sa,
                ];
                Some(if h.convex { n } else { neg(n) })
            }
            Surf::Spherical(h) => {
                let d = sub(*p, h.center);
                let l = norm(d);
                if l < 1e-15 {
                    return None;
                }
                let n = scale(d, 1.0 / l);
                Some(if h.convex { n } else { neg(n) })
            }
            Surf::Toroidal(h) => {
                let (rho, radial, height) = axis_frame(p, &h.center, &h.axis_direction);
                if rho < 1e-15 {
                    return None;
                }
                // Vector from the tube centreline to p, in the meridian frame.
                let (dr, dh) = (rho - h.major_radius, height);
                let l = (dr * dr + dh * dh).sqrt();
                if l < 1e-15 {
                    return None;
                }
                let n = [
                    (radial[0] * dr + h.axis_direction[0] * dh) / l,
                    (radial[1] * dr + h.axis_direction[1] * dh) / l,
                    (radial[2] * dr + h.axis_direction[2] * dh) / l,
                ];
                Some(if h.convex { n } else { neg(n) })
            }
            Surf::BSpline => None,
        }
    }

    /// Absolute distance from `p` to the fitted surface. B-splines have no
    /// closed form here and never reach a locus.
    pub(crate) fn distance_at(&self, p: &[f64; 3]) -> f64 {
        match self {
            Surf::Planar(h) => (dot(*p, h.normal) - h.distance).abs(),
            Surf::Conical(h) => {
                let (rho, _, height) = axis_frame(p, &h.apex, &h.axis_direction);
                if h.is_cylinder() {
                    (rho - h.radius_at_apex).abs()
                } else {
                    (rho * h.half_angle.cos() - height * h.half_angle.sin()).abs()
                }
            }
            Surf::Spherical(h) => (norm(sub(*p, h.center)) - h.radius).abs(),
            Surf::Toroidal(h) => {
                let (rho, _, height) = axis_frame(p, &h.center, &h.axis_direction);
                let dr = rho - h.major_radius;
                (((dr * dr + height * height).sqrt()) - h.minor_radius).abs()
            }
            Surf::BSpline => f64::INFINITY,
        }
    }

    /// A radius characterising the surface's size, used to scale the
    /// coaxiality test. Zero when the type has no radius (a plane).
    fn characteristic_radius(&self) -> f64 {
        match self {
            Surf::Planar(_) | Surf::BSpline => 0.0,
            Surf::Conical(h) => h.radius_at_apex.abs(),
            Surf::Spherical(h) => h.radius,
            Surf::Toroidal(h) => h.major_radius + h.minor_radius,
        }
    }
}

// ---------------------------------------------------------------------------
// The locus
// ---------------------------------------------------------------------------

/// A signed side test for the junction between two surfaces: negative on the
/// first surface's side, positive on the second's, zero on the junction.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Locus {
    /// Axisymmetric junction: a circle about `axis`. `s` is the displacement
    /// in the meridian half-plane from the tangency point `q` along the two
    /// meridian curves' common tangent `t`.
    Meridian {
        origin: [f64; 3],
        axis: [f64; 3],
        q: [f64; 2],
        t: [f64; 2],
    },
    /// A straight junction: a line where a plane grazes a cylinder. `s` is the
    /// signed distance to the plane through `p` with normal `n`, where `n`
    /// lies in the common tangent plane perpendicular to the line.
    Planar { p: [f64; 3], n: [f64; 3] },
}

impl Locus {
    /// The signed side of `p`: negative belongs to the first surface of the
    /// pair the locus was built from, positive to the second.
    pub(crate) fn side(&self, p: &[f64; 3]) -> f64 {
        match self {
            Locus::Meridian { origin, axis, q, t } => {
                let (rho, _, h) = axis_frame(p, origin, axis);
                (rho - q[0]) * t[0] + (h - q[1]) * t[1]
            }
            Locus::Planar { p: p0, n } => dot(sub(*p, *p0), *n),
        }
    }

    /// Short description, for diagnostics.
    pub(crate) fn kind(&self) -> &'static str {
        match self {
            Locus::Meridian { .. } => "circle",
            Locus::Planar { .. } => "line",
        }
    }

    /// Flip the sign convention, so that `side` is negative on the other
    /// surface. Used once the caller has established which support lies on
    /// which side.
    pub(crate) fn flipped(&self) -> Locus {
        match *self {
            Locus::Meridian { origin, axis, q, t } => Locus::Meridian {
                origin,
                axis,
                q,
                t: [-t[0], -t[1]],
            },
            Locus::Planar { p, n } => Locus::Planar { p, n: neg(n) },
        }
    }
}

/// Why no locus could be built for a pair, for the loud report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LocusFailure {
    /// One side is a B-spline: no closed-form geometry to intersect.
    Freeform,
    /// Two planes, or a plane and a sphere. Not an unimplemented junction —
    /// no junction *exists*. Two tangent planes are the same plane (so any
    /// boundary between them is merge's business, not snap's, and a smooth
    /// surface left as per-facet planes produces thousands of such pairs), and
    /// a plane grazes a sphere at a single point rather than along a curve.
    NoJunctionExists,
    /// A plane grazing a cone along a ruling. Constructible in principle;
    /// deliberately not implemented until a corpus model needs it.
    ConePlaneRuling,
    /// The two surfaces do not share an axis, so there is no meridian plane.
    NotCoaxial,
    /// Both meridian curves are straight and so cannot be tangent unless they
    /// coincide — the pair is one surface split in two.
    MeridianLinesParallel,
    /// The tangency point landed on the axis: the junction degenerates to a
    /// single point rather than a circle.
    DegenerateOnAxis,
    /// The construction produced a non-finite result.
    NonFinite,
}

impl LocusFailure {
    pub(crate) fn as_str(&self) -> &'static str {
        match self {
            LocusFailure::Freeform => "one side is a B-spline",
            LocusFailure::NoJunctionExists => "plane/plane or plane/sphere: no junction curve exists",
            LocusFailure::ConePlaneRuling => "cone-plane ruling tangency is not implemented",
            LocusFailure::NotCoaxial => "surfaces do not share an axis",
            LocusFailure::MeridianLinesParallel => "meridian curves are parallel lines",
            LocusFailure::DegenerateOnAxis => "tangency degenerates to a point on the axis",
            LocusFailure::NonFinite => "construction produced a non-finite result",
        }
    }
}

/// Build the junction locus for a tangent pair `(a, b)`, oriented so that
/// `side` is negative on `a`'s side before the caller's own orientation check.
///
/// `sample` is a point on the shared boundary, used only to pick between the
/// two branches a construction can offer (which side of a cone's apex, which
/// of a cylinder's two tangent rulings) — never to place the locus.
pub(crate) fn build_locus(
    a: Surf<'_>,
    b: Surf<'_>,
    sample: &[f64; 3],
) -> Result<Locus, LocusFailure> {
    if matches!(a, Surf::BSpline) || matches!(b, Surf::BSpline) {
        return Err(LocusFailure::Freeform);
    }

    // A plane grazing a cylinder or cone along a ruling is not axisymmetric:
    // recognised by the plane's normal being perpendicular to the axis.
    if let Some((plane, cone, flip)) = as_plane_conical(a, b) {
        // Along a ruling the conical surface's normal equals ±the plane's, and
        // a cone's normal makes the angle (90° − half-angle) with its axis —
        // so the test is |n·a| = sin(half_angle), which is 0 only in the
        // cylinder case. Testing for the cylinder condition alone reported
        // every cone-plane ruling tangency as "not coaxial".
        let axial = dot(plane.normal, cone.axis_direction).abs();
        if (axial - cone.half_angle.sin()).abs() < COAXIAL_TOL {
            if !cone.is_cylinder() {
                return Err(LocusFailure::ConePlaneRuling);
            }
            let locus = cylinder_ruling_locus(plane, cone)?;
            return Ok(if flip { locus.flipped() } else { locus });
        }
    }

    // A fillet running along a straight edge is a cylinder; where the edge
    // turns, the same tube becomes a torus. The two share the *tube*, not an
    // axis of revolution — their axes are perpendicular — so the meridian
    // reduction does not apply and the junction is a cross-section of the
    // tube rather than a circle about either axis.
    if let Some((torus, cyl, flip)) = as_torus_conical(a, b) {
        if cyl.is_cylinder()
            && dot(torus.axis_direction, cyl.axis_direction).abs() < COAXIAL_TOL
        {
            let locus = tube_junction_locus(torus, cyl)?;
            return Ok(if flip { locus.flipped() } else { locus });
        }
    }

    let (origin, axis) = shared_axis(a, b, sample)?;
    let ma = meridian(a, &origin, &axis, sample)?;
    let mb = meridian(b, &origin, &axis, sample)?;
    let (q, t) = meridian_tangency(&ma, &mb)?;
    if q[0] <= 0.0 {
        return Err(LocusFailure::DegenerateOnAxis);
    }
    if !q.iter().chain(t.iter()).all(|v| v.is_finite()) {
        return Err(LocusFailure::NonFinite);
    }
    Ok(Locus::Meridian { origin, axis, q, t })
}

/// The pair as (plane, conical, flipped) when exactly one side is planar and
/// the other conical.
fn as_plane_conical<'a>(
    a: Surf<'a>,
    b: Surf<'a>,
) -> Option<(&'a PlanarHypothesis, &'a ConicalHypothesis, bool)> {
    match (a, b) {
        (Surf::Planar(p), Surf::Conical(c)) => Some((p, c, false)),
        (Surf::Conical(c), Surf::Planar(p)) => Some((p, c, true)),
        _ => None,
    }
}

/// The pair as (torus, conical, flipped) when it is one of each.
fn as_torus_conical<'a>(
    a: Surf<'a>,
    b: Surf<'a>,
) -> Option<(&'a ToroidalHypothesis, &'a ConicalHypothesis, bool)> {
    match (a, b) {
        (Surf::Toroidal(t), Surf::Conical(c)) => Some((t, c, false)),
        (Surf::Conical(c), Surf::Toroidal(t)) => Some((t, c, true)),
        _ => None,
    }
}

/// Where a straight tube (cylinder) continues into a bent one (torus): the
/// cross-section of the tube at the point where its spine stops being
/// straight, as a `Planar` locus with `side` negative on the torus's side.
///
/// C¹ continuity of the tube forces the cylinder's axis to be *tangent* to the
/// torus's spine circle, and a line tangent to a circle touches it at the foot
/// of the perpendicular from the centre — so the junction is pinned by the
/// torus centre and the cylinder axis alone, both of which their supports
/// determine well. Neither fit says where its own surface *ends*, which is why
/// the junction has to be read off the pair rather than either surface.
fn tube_junction_locus(
    torus: &ToroidalHypothesis,
    cyl: &ConicalHypothesis,
) -> Result<Locus, LocusFailure> {
    let axis = unit(cyl.axis_direction);
    let d = sub(torus.center, cyl.apex);
    let foot = add(cyl.apex, scale(axis, dot(d, axis)));
    if !foot.iter().all(|v| v.is_finite()) {
        return Err(LocusFailure::NonFinite);
    }
    Ok(Locus::Planar { p: foot, n: axis })
}

/// The tangency line where a plane grazes a cylinder, as a `Planar` locus with
/// `side` negative on the plane's side.
fn cylinder_ruling_locus(
    plane: &PlanarHypothesis,
    cyl: &ConicalHypothesis,
) -> Result<Locus, LocusFailure> {
    let r = cyl.radius_at_apex;
    if r <= 0.0 || !r.is_finite() {
        return Err(LocusFailure::NonFinite);
    }
    // The cylinder touches the plane where its outward radial direction is
    // ±the plane normal; the offset decides which of the two rulings.
    let gap = plane.distance - dot(cyl.apex, plane.normal);
    let sigma = if gap >= 0.0 { 1.0 } else { -1.0 };
    let p0 = add(cyl.apex, scale(plane.normal, sigma * r));
    // In the shared tangent plane, the direction across the ruling.
    let m = cross(cyl.axis_direction, plane.normal);
    let l = norm(m);
    if l < COAXIAL_TOL {
        return Err(LocusFailure::NonFinite);
    }
    Ok(Locus::Planar {
        p: p0,
        n: scale(m, 1.0 / l),
    })
}

/// The axis both surfaces revolve about, as (point, unit direction).
fn shared_axis(
    a: Surf<'_>,
    b: Surf<'_>,
    sample: &[f64; 3],
) -> Result<([f64; 3], [f64; 3]), LocusFailure> {
    let scale_r = a.characteristic_radius().max(b.characteristic_radius());
    match (a.axis(), b.axis()) {
        (Some((oa, da)), Some((ob, db))) => {
            // Both pin an axis: the lines must coincide.
            if norm(cross(da, db)) > COAXIAL_TOL {
                return Err(LocusFailure::NotCoaxial);
            }
            let (_, _, _) = (oa, ob, db);
            let off = point_axis_distance(&ob, &oa, &da);
            if off > COAXIAL_TOL * scale_r.max(1.0) {
                return Err(LocusFailure::NotCoaxial);
            }
            Ok((oa, unit(da)))
        }
        (Some((o, d)), None) => partner_on_axis(b, &o, &d, scale_r).map(|_| (o, unit(d))),
        (None, Some((o, d))) => partner_on_axis(a, &o, &d, scale_r).map(|_| (o, unit(d))),
        (None, None) => match (a, b) {
            // Two spheres: the axis is their line of centres.
            (Surf::Spherical(sa), Surf::Spherical(sb)) => {
                let d = sub(sb.center, sa.center);
                let l = norm(d);
                if l < COAXIAL_TOL * scale_r.max(1.0) {
                    return Err(LocusFailure::DegenerateOnAxis);
                }
                Ok((sa.center, scale(d, 1.0 / l)))
            }
            // A sphere and a plane graze at a single point; two planes have no
            // junction at all unless they are the same plane.
            _ => {
                let _ = sample;
                Err(LocusFailure::NoJunctionExists)
            }
        },
    }
}

/// Check that the axis-less partner is compatible with the axis the other
/// surface pins: a plane must be perpendicular to it, a sphere centred on it.
fn partner_on_axis(
    s: Surf<'_>,
    o: &[f64; 3],
    d: &[f64; 3],
    scale_r: f64,
) -> Result<(), LocusFailure> {
    match s {
        Surf::Planar(h) => {
            if norm(cross(h.normal, unit(*d))) > COAXIAL_TOL {
                Err(LocusFailure::NotCoaxial)
            } else {
                Ok(())
            }
        }
        Surf::Spherical(h) => {
            if point_axis_distance(&h.center, o, d) > COAXIAL_TOL * scale_r.max(1.0) {
                Err(LocusFailure::NotCoaxial)
            } else {
                Ok(())
            }
        }
        _ => Err(LocusFailure::NotCoaxial),
    }
}

// ---------------------------------------------------------------------------
// Meridian curves
// ---------------------------------------------------------------------------

/// A surface of revolution's profile in the (ρ, h) half-plane of the shared
/// axis: a line for planes/cylinders/cones, a circle for spheres/tori.
#[derive(Debug, Clone, Copy)]
enum Meridian {
    Line { p: [f64; 2], d: [f64; 2] },
    Circle { c: [f64; 2], r: f64 },
}

/// The surface's meridian curve in the frame of the shared axis.
fn meridian(
    s: Surf<'_>,
    origin: &[f64; 3],
    axis: &[f64; 3],
    sample: &[f64; 3],
) -> Result<Meridian, LocusFailure> {
    let (rho_s, _, h_s) = axis_frame(sample, origin, axis);
    match s {
        Surf::Planar(hyp) => {
            // Perpendicular to the axis: a horizontal line at the plane's own
            // height. (A plane parallel to the axis never reaches here — that
            // is the ruling case, handled before the meridian reduction.)
            let on_plane = scale(hyp.normal, hyp.distance);
            let (_, _, h0) = axis_frame(&on_plane, origin, axis);
            Ok(Meridian::Line {
                p: [0.0, h0],
                d: [1.0, 0.0],
            })
        }
        Surf::Conical(hyp) if hyp.is_cylinder() => Ok(Meridian::Line {
            p: [hyp.radius_at_apex, 0.0],
            d: [0.0, 1.0],
        }),
        Surf::Conical(hyp) => {
            let (_, _, h_apex) = axis_frame(&hyp.apex, origin, axis);
            // The cone is two rays from the apex; the support picks the branch.
            let branch = if h_s >= h_apex { 1.0 } else { -1.0 };
            let (sa, ca) = (hyp.half_angle.sin(), hyp.half_angle.cos());
            Ok(Meridian::Line {
                p: [0.0, h_apex],
                d: [sa, branch * ca],
            })
        }
        Surf::Spherical(hyp) => {
            let (rho_c, _, h_c) = axis_frame(&hyp.center, origin, axis);
            let _ = rho_s;
            Ok(Meridian::Circle {
                c: [rho_c, h_c],
                r: hyp.radius,
            })
        }
        Surf::Toroidal(hyp) => {
            let (_, _, h_c) = axis_frame(&hyp.center, origin, axis);
            Ok(Meridian::Circle {
                c: [hyp.major_radius, h_c],
                r: hyp.minor_radius,
            })
        }
        Surf::BSpline => Err(LocusFailure::Freeform),
    }
}

/// The point where two meridian curves touch, and their common unit tangent
/// there.
///
/// Taken as the point of closest approach, which *is* the tangency point when
/// the curves are tangent and degrades gracefully as the fits drift apart.
/// Two straight meridians cannot be tangent without coinciding, so that case
/// is rejected rather than solved.
fn meridian_tangency(a: &Meridian, b: &Meridian) -> Result<([f64; 2], [f64; 2]), LocusFailure> {
    match (a, b) {
        (Meridian::Line { .. }, Meridian::Line { .. }) => Err(LocusFailure::MeridianLinesParallel),
        (Meridian::Line { p, d }, Meridian::Circle { c, .. })
        | (Meridian::Circle { c, .. }, Meridian::Line { p, d }) => {
            // Closest approach is the foot of the perpendicular from the
            // circle's centre; the common tangent is the line's direction.
            let n = [-d[1], d[0]];
            let off = (c[0] - p[0]) * n[0] + (c[1] - p[1]) * n[1];
            let q = [c[0] - off * n[0], c[1] - off * n[1]];
            Ok((q, *d))
        }
        (Meridian::Circle { c: c1, r: r1 }, Meridian::Circle { c: c2, r: r2 }) => {
            let u = [c2[0] - c1[0], c2[1] - c1[1]];
            let l = (u[0] * u[0] + u[1] * u[1]).sqrt();
            if l < 1e-12 {
                // Concentric: no isolated tangency.
                return Err(LocusFailure::DegenerateOnAxis);
            }
            let u = [u[0] / l, u[1] / l];
            // Tangency lies on the line of centres; internal and external
            // tangency put it on opposite sides of the first circle.
            let cand = [
                [c1[0] + r1 * u[0], c1[1] + r1 * u[1]],
                [c1[0] - r1 * u[0], c1[1] - r1 * u[1]],
            ];
            let gap = |q: &[f64; 2]| {
                (((q[0] - c2[0]).powi(2) + (q[1] - c2[1]).powi(2)).sqrt() - r2).abs()
            };
            let q = if gap(&cand[0]) <= gap(&cand[1]) {
                cand[0]
            } else {
                cand[1]
            };
            Ok((q, [-u[1], u[0]]))
        }
    }
}

// ---------------------------------------------------------------------------
// Tangency
// ---------------------------------------------------------------------------

/// Whether two committed surfaces meet tangentially, judged by their fitted
/// normals at every sampled point of their shared boundary.
///
/// The same test stage 3.2 applies to reconstructed edges; a junction the two
/// stages disagreed about would be one where stage 2 placed the boundary by a
/// criterion stage 3 then declines to honour.
pub(crate) fn is_tangent_pair(
    a: Surf<'_>,
    b: Surf<'_>,
    boundary_points: &[[f64; 3]],
) -> bool {
    if boundary_points.is_empty() {
        return false;
    }
    let cos_tol = TANGENCY_TOL_RAD.cos();
    boundary_points.iter().all(|p| match (a.normal_at(p), b.normal_at(p)) {
        (Some(na), Some(nb)) => dot(na, nb) >= cos_tol,
        _ => false,
    })
}

// ---------------------------------------------------------------------------
// Pair analysis: what the locus says about the boundary the mesh currently has
// ---------------------------------------------------------------------------

/// What building a locus for one adjacent pair of surfaces concluded.
#[derive(Debug, Clone)]
pub(crate) enum PairOutcome {
    /// The surfaces cross rather than graze. Out of scope: a transversal
    /// junction is a sharp mesh edge, which stage 2.1 already partitions
    /// exactly, so there is no boundary to correct.
    Transversal,
    /// Tangent, with a locus that passed validation. `moves` lists the faces
    /// it would reassign, as (face, from, to).
    Snapped {
        /// Oriented so `side` is negative on the lower-indexed surface.
        locus: Locus,
        /// The junction's own length scale (median interface facet).
        scale: f64,
        /// How far the locus sits from the current boundary chain.
        worst_on_chain: f64,
        moves: Vec<(usize, usize, usize)>,
    },
    /// Tangent, but no locus could be constructed for this pair of types.
    NoLocus(LocusFailure),
    /// Tangent, locus constructed, but it failed a validation check.
    Rejected(String),
}

/// One adjacent pair of committed surfaces and what the locus concluded.
#[derive(Debug, Clone)]
pub(crate) struct PairAnalysis {
    pub a: usize,
    pub b: usize,
    pub kind_a: &'static str,
    pub kind_b: &'static str,
    /// Faces of either surface that touch the other across a mesh edge.
    pub interface_faces: usize,
    pub outcome: PairOutcome,
}

impl PairAnalysis {
    /// True when snap took responsibility for the pair and something is wrong
    /// — the cases the caller reports loudly.
    ///
    /// A pair whose surfaces simply have no junction between them is not one
    /// of those: it is out of scope in the same way a transversal junction is.
    /// This matters in practice, not just in principle — a smooth surface left
    /// as per-facet planes makes every adjacent facet pair a tangent
    /// plane/plane pair, which would bury a real report under thousands of
    /// lines on exactly the models most likely to have one.
    pub(crate) fn is_problem(&self) -> bool {
        match &self.outcome {
            PairOutcome::NoLocus(LocusFailure::NoJunctionExists) => false,
            PairOutcome::NoLocus(_) | PairOutcome::Rejected(_) => true,
            PairOutcome::Transversal | PairOutcome::Snapped { .. } => false,
        }
    }

    pub(crate) fn describe(&self) -> String {
        match &self.outcome {
            PairOutcome::Transversal => "transversal (out of scope)".to_string(),
            PairOutcome::Snapped {
                locus,
                scale,
                worst_on_chain,
                moves,
            } => format!(
                "tangent, {} locus {:.2e} from the chain (facet {:.2e}), {} face(s) to move",
                locus.kind(),
                worst_on_chain,
                scale,
                moves.len(),
            ),
            PairOutcome::NoLocus(f) => format!("tangent, NO LOCUS: {}", f.as_str()),
            PairOutcome::Rejected(why) => format!("tangent, locus REJECTED: {why}"),
        }
    }
}

/// Analyse every adjacent pair of committed surfaces.
///
/// `surfs` and `f2s` index the same selected-surface space; faces not assigned
/// to any surface carry `usize::MAX`. `vtol` is `--vertex-tolerance`, used for
/// the soundness precondition on a move: a face is only ever handed to a
/// surface that actually fits it.
pub(crate) fn analyze_pairs(
    mesh: &ConnectedMesh,
    surfs: &[Surf<'_>],
    f2s: &[usize],
    vtol: f64,
) -> Vec<PairAnalysis> {
    let mut counts = vec![0usize; surfs.len()];
    for &si in f2s {
        if si != usize::MAX {
            counts[si] += 1;
        }
    }
    let mut out = Vec::new();
    for (a, b, interface) in interfaces(mesh, f2s) {
        out.push(analyze_pair(mesh, surfs, f2s, &counts, a, b, &interface, vtol));
    }
    out
}

/// Per-face detail for one pair's interface: why each face does or does not
/// move. Diagnostic only.
pub(crate) fn explain_pair(
    mesh: &ConnectedMesh,
    surfs: &[Surf<'_>],
    f2s: &[usize],
    vtol: f64,
    a: usize,
    b: usize,
) -> Vec<String> {
    let mut counts = vec![0usize; surfs.len()];
    for &si in f2s {
        if si != usize::MAX {
            counts[si] += 1;
        }
    }
    let Some((_, _, interface)) = interfaces(mesh, f2s)
        .into_iter()
        .find(|&(x, y, _)| (x, y) == (a.min(b), a.max(b)))
    else {
        return vec![format!("surfaces {a} and {b} do not share a boundary")];
    };
    let (a, b) = (a.min(b), a.max(b));
    let analysis = analyze_pair(mesh, surfs, f2s, &counts, a, b, &interface, vtol);
    let mut out = vec![
        format!(
            "pair {a} ({}) | {b} ({}): {}",
            surfs[a].kind(),
            surfs[b].kind(),
            analysis.describe()
        ),
        format!(
            "  {} faces on {a}'s side, {} on {b}'s; {} chain vertices",
            interface.faces_a.len(),
            interface.faces_b.len(),
            interface.vertices.len(),
        ),
    ];
    let PairOutcome::Snapped { locus, moves, .. } = &analysis.outcome else {
        return out;
    };
    let moving: std::collections::BTreeSet<usize> = moves.iter().map(|m| m.0).collect();
    out.push(
        "  face      owner  diameter    s(centroid)   s(min)     s(max)  fits:a b  move".to_string(),
    );
    for &fi in interface.faces_a.iter().chain(interface.faces_b.iter()) {
        let face = &mesh.faces[fi];
        let vc = face.vertex_count as usize;
        let sv: Vec<f64> = face.vertex_indices[..vc]
            .iter()
            .map(|&vi| locus.side(&point_of(&mesh.vertices[vi])))
            .collect();
        out.push(format!(
            "  {:>7}  {:>5}  {:>9.3e}  {:>11.3e}  {:>9.2e}  {:>9.2e}   {} {}   {}",
            fi,
            f2s[fi],
            face_diameter(mesh, fi),
            face_side(mesh, fi, locus),
            sv.iter().cloned().fold(f64::INFINITY, f64::min),
            sv.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            if face_fits(mesh, surfs[a], fi, vtol) { "Y" } else { "." },
            if face_fits(mesh, surfs[b], fi, vtol) { "Y" } else { "." },
            if moving.contains(&fi) { "YES" } else { "" },
        ));
    }
    out
}

/// The faces on each side of every surface-to-surface interface, keyed by the
/// ordered pair (a < b) of selected-surface indices.
fn interfaces(mesh: &ConnectedMesh, f2s: &[usize]) -> Vec<(usize, usize, Interface)> {
    use std::collections::BTreeMap;
    let mut map: BTreeMap<(usize, usize), Interface> = BTreeMap::new();
    for fi in 0..mesh.faces.len() {
        let si = f2s[fi];
        if si == usize::MAX {
            continue;
        }
        let face = &mesh.faces[fi];
        let vc = face.vertex_count as usize;
        for k in 0..vc {
            let ni = face.neighbors[k];
            if ni < 0 {
                continue;
            }
            let nj = ni as usize;
            let sj = f2s[nj];
            if sj == usize::MAX || sj == si {
                continue;
            }
            let e = map.entry((si.min(sj), si.max(sj))).or_default();
            if si < sj {
                e.faces_a.push(fi);
            } else {
                e.faces_b.push(fi);
            }
            // The two vertices of the shared edge lie on the boundary chain.
            e.vertices.push(face.vertex_indices[k]);
            e.vertices.push(face.vertex_indices[(k + 1) % vc]);
        }
    }
    map.into_iter()
        .map(|((a, b), mut i)| {
            for v in [&mut i.faces_a, &mut i.faces_b, &mut i.vertices] {
                v.sort_unstable();
                v.dedup();
            }
            (a, b, i)
        })
        .collect()
}

#[derive(Default)]
struct Interface {
    /// Faces of the lower-indexed surface that touch the other.
    faces_a: Vec<usize>,
    /// Faces of the higher-indexed surface that touch the other.
    faces_b: Vec<usize>,
    /// Vertices of the shared mesh-edge chain.
    vertices: Vec<usize>,
}

#[allow(clippy::too_many_arguments)]
fn analyze_pair(
    mesh: &ConnectedMesh,
    surfs: &[Surf<'_>],
    f2s: &[usize],
    counts: &[usize],
    a: usize,
    b: usize,
    interface: &Interface,
    vtol: f64,
) -> PairAnalysis {
    let (sa, sb) = (surfs[a], surfs[b]);
    let mut result = PairAnalysis {
        a,
        b,
        kind_a: sa.kind(),
        kind_b: sb.kind(),
        interface_faces: interface.faces_a.len() + interface.faces_b.len(),
        outcome: PairOutcome::Transversal,
    };

    let pts: Vec<[f64; 3]> = interface
        .vertices
        .iter()
        .map(|&vi| point_of(&mesh.vertices[vi]))
        .collect();
    if !is_tangent_pair(sa, sb, &pts) {
        return result;
    }

    // Build the locus about the interface's own centroid, so that a
    // construction with two branches picks the one the boundary is on.
    let locus = match build_locus(sa, sb, &mean_point(&pts)) {
        Ok(l) => l,
        Err(f) => {
            result.outcome = PairOutcome::NoLocus(f);
            return result;
        }
    };

    // The junction's own length scale: the *median* facet across the
    // interface, not the largest. The two sides are often tessellated very
    // differently — a whole flat wall can be one fused quad against a hundred
    // facets of fillet — and it is the finely tessellated side that carries
    // any misplaced row, so the fine scale is the one that matters. (Median
    // for the same reason `residual_scale` uses one: a few outsized faces on
    // one side must not set the scale for the junction.)
    let scale = median_face_diameter(mesh, &interface.faces_a, &interface.faces_b);

    // The seam assumption says the true boundary is a mesh-edge chain at most
    // about a row from the one competition produced, so a locus further than
    // that from the current chain is not a correction to it — it is describing
    // some other junction. The bar is the mesh's own local scale.
    let worst = pts.iter().map(|p| locus.side(p).abs()).fold(0.0_f64, f64::max);
    if worst > scale {
        result.outcome = PairOutcome::Rejected(format!(
            "locus is {worst:.3e} from the boundary, more than the {scale:.3e} facet scale"
        ));
        return result;
    }

    // Orient the locus so `side` is negative on a's side. The two surfaces'
    // *supports* cannot settle this between them: a locus circle is complete
    // but the junction is only the arc where the surfaces actually meet, so a
    // large planar face wanders to the far side of its own tangency circle
    // while still being unambiguously the plane's. Read the orientation off
    // the finely tessellated side alone, whose extent is comparable to the
    // junction, and only from its faces the *other* surface cannot explain —
    // those are unambiguously its own.
    let locus = match orient(mesh, surfs, a, b, interface, &locus, vtol) {
        Some(l) => l,
        None => {
            result.outcome =
                PairOutcome::Rejected("cannot orient: every nearby face fits both".to_string());
            return result;
        }
    };

    let moves = movable_faces(mesh, surfs, f2s, a, b, interface, &locus, vtol);
    // A boundary correction never dissolves a surface; if the locus says it
    // should, the locus is wrong about which junction this is.
    for si in [a, b] {
        let leaving = moves.iter().filter(|m| m.1 == si).count();
        if leaving >= counts[si] {
            result.outcome = PairOutcome::Rejected(format!(
                "would empty surface {si} ({leaving} of {} faces)",
                counts[si]
            ));
            return result;
        }
    }

    result.outcome = PairOutcome::Snapped {
        locus,
        scale,
        worst_on_chain: worst,
        moves,
    };
    result
}

/// Orient the locus so that `side` is negative on surface `a`.
///
/// Uses the interface faces of whichever side is more finely tessellated, and
/// among those only the ones the opposite surface does not fit — a face both
/// surfaces explain is exactly the ambiguous kind whose side is in question.
fn orient(
    mesh: &ConnectedMesh,
    surfs: &[Surf<'_>],
    a: usize,
    b: usize,
    interface: &Interface,
    locus: &Locus,
    vtol: f64,
) -> Option<Locus> {
    let med_a = median_face_diameter(mesh, &interface.faces_a, &[]);
    let med_b = median_face_diameter(mesh, &[], &interface.faces_b);
    // Try the finer side first; fall back to the coarser one if every face
    // there is ambiguous.
    let order: [(usize, &Vec<usize>, usize); 2] = if med_a <= med_b {
        [
            (a, &interface.faces_a, b),
            (b, &interface.faces_b, a),
        ]
    } else {
        [
            (b, &interface.faces_b, a),
            (a, &interface.faces_a, b),
        ]
    };
    for (own, faces, other) in order {
        let net: f64 = faces
            .iter()
            .filter(|&&fi| !face_fits(mesh, surfs[other], fi, vtol))
            .map(|&fi| face_side(mesh, fi, locus).signum())
            .sum();
        if net == 0.0 {
            continue;
        }
        // `net > 0` means this surface sits on the locus's positive side.
        let own_is_positive = net > 0.0;
        let a_is_positive = if own == a { own_is_positive } else { !own_is_positive };
        return Some(if a_is_positive { locus.flipped() } else { locus.clone() });
    }
    None
}

/// The faces the locus reassigns, as (face, from, to).
///
/// A face moves only when three things hold at once: the surface it would
/// move to actually fits it (soundness — a boundary correction must never
/// hand a face to a surface that does not explain it); the locus puts the
/// whole face on that surface's side; and it is reachable from the interface
/// through other moving faces (a boundary correction moves a band attached to
/// the boundary, not a detached patch).
#[allow(clippy::too_many_arguments)]
fn movable_faces(
    mesh: &ConnectedMesh,
    surfs: &[Surf<'_>],
    f2s: &[usize],
    a: usize,
    b: usize,
    interface: &Interface,
    locus: &Locus,
    vtol: f64,
) -> Vec<(usize, usize, usize)> {
    let wants_move = |fi: usize, from: usize, to: usize| -> bool {
        if !face_fits(mesh, surfs[to], fi, vtol) {
            return false;
        }
        // `a` owns the negative side, so a face leaving `a` must lie wholly on
        // the positive side, and vice versa.
        let want_positive = from == a;
        let face = &mesh.faces[fi];
        let vc = face.vertex_count as usize;
        for &vi in &face.vertex_indices[..vc] {
            let s = locus.side(&point_of(&mesh.vertices[vi]));
            if (want_positive && s < 0.0) || (!want_positive && s > 0.0) {
                return false;
            }
        }
        // The corners alone are not enough: a face can have *every* vertex on
        // the junction and still lie squarely on one side of it — a disc
        // triangulated with all its vertices on the rim is exactly that, and
        // those vertices sit on both surfaces, so they pass the fit test too.
        // The interior is what decides, so require the centroid to be on the
        // target's side outright.
        let s = face_side(mesh, fi, locus);
        if want_positive {
            s > 0.0
        } else {
            s < 0.0
        }
    };

    // Grow outward from the interface through faces that want to move, so a
    // band more than one row deep is corrected in full while a detached patch
    // elsewhere on either surface is never touched.
    let mut moves: Vec<(usize, usize, usize)> = Vec::new();
    let mut queue: Vec<usize> = interface
        .faces_a
        .iter()
        .chain(interface.faces_b.iter())
        .copied()
        .collect();
    let mut seen = vec![false; mesh.faces.len()];
    while let Some(fi) = queue.pop() {
        if seen[fi] {
            continue;
        }
        seen[fi] = true;
        let from = f2s[fi];
        if from != a && from != b {
            continue;
        }
        let to = if from == a { b } else { a };
        if !wants_move(fi, from, to) {
            continue;
        }
        moves.push((fi, from, to));
        let face = &mesh.faces[fi];
        for k in 0..face.vertex_count as usize {
            let ni = face.neighbors[k];
            if ni >= 0 && !seen[ni as usize] {
                queue.push(ni as usize);
            }
        }
    }
    moves.sort_unstable();
    moves
}

/// Whether every vertex of `fi` lies within `vtol` of `s`.
fn face_fits(mesh: &ConnectedMesh, s: Surf<'_>, fi: usize, vtol: f64) -> bool {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    face.vertex_indices[..vc]
        .iter()
        .all(|&vi| s.distance_at(&point_of(&mesh.vertices[vi])) <= vtol)
}

/// Median of the faces' diameters over both sides of an interface.
fn median_face_diameter(mesh: &ConnectedMesh, fa: &[usize], fb: &[usize]) -> f64 {
    let mut d: Vec<f64> = fa
        .iter()
        .chain(fb.iter())
        .map(|&fi| face_diameter(mesh, fi))
        .collect();
    if d.is_empty() {
        return 0.0;
    }
    d.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    d[d.len() / 2]
}

/// Largest distance between any two of a face's vertices — the mesh's own
/// local length scale at that face.
fn face_diameter(mesh: &ConnectedMesh, fi: usize) -> f64 {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    let mut d: f64 = 0.0;
    for i in 0..vc {
        for j in (i + 1)..vc {
            let p = point_of(&mesh.vertices[face.vertex_indices[i]]);
            let q = point_of(&mesh.vertices[face.vertex_indices[j]]);
            d = d.max(norm(sub(p, q)));
        }
    }
    d
}

fn mean_point(pts: &[[f64; 3]]) -> [f64; 3] {
    let n = pts.len().max(1) as f64;
    let mut s = [0.0; 3];
    for p in pts {
        s = add(s, *p);
    }
    scale(s, 1.0 / n)
}

// ---------------------------------------------------------------------------
// Small vector helpers (kept local: stage 2 has no shared vector module)
// ---------------------------------------------------------------------------

/// Decompose `p` in the frame of the axis `(origin, dir)`, returning the
/// radial distance ρ, the unit radial direction, and the axial coordinate h.
fn axis_frame(p: &[f64; 3], origin: &[f64; 3], dir: &[f64; 3]) -> (f64, [f64; 3], f64) {
    let d = sub(*p, *origin);
    let h = dot(d, *dir);
    let radial = [
        d[0] - h * dir[0],
        d[1] - h * dir[1],
        d[2] - h * dir[2],
    ];
    let rho = norm(radial);
    let unit_radial = if rho > 1e-15 {
        scale(radial, 1.0 / rho)
    } else {
        [0.0, 0.0, 0.0]
    };
    (rho, unit_radial, h)
}

/// Perpendicular distance from `p` to the axis line through `o` along `d`.
fn point_axis_distance(p: &[f64; 3], o: &[f64; 3], d: &[f64; 3]) -> f64 {
    axis_frame(p, o, &unit(*d)).0
}

pub(crate) fn point_of(v: &MeshVertex) -> [f64; 3] {
    [v.x, v.y, v.z]
}

/// A face's centroid.
///
/// The locus is evaluated here rather than averaged over the vertices: `side`
/// is not affine in position (the meridian radius is not), and the difference
/// is exactly what matters for a face whose vertices all sit *on* the junction.
fn face_centroid(mesh: &ConnectedMesh, fi: usize) -> [f64; 3] {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    let mut c = [0.0; 3];
    for &vi in &face.vertex_indices[..vc] {
        c = add(c, point_of(&mesh.vertices[vi]));
    }
    scale(c, 1.0 / vc as f64)
}

/// Which side of the locus a face's interior lies on.
fn face_side(mesh: &ConnectedMesh, fi: usize, locus: &Locus) -> f64 {
    locus.side(&face_centroid(mesh, fi))
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn neg(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = norm(a);
    if l > 1e-15 {
        scale(a, 1.0 / l)
    } else {
        a
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane(normal: [f64; 3], distance: f64) -> PlanarHypothesis {
        PlanarHypothesis {
            normal,
            distance,
            faces: Vec::new(),
            vertices: Vec::new(),
            error_max: 0.0,
            error_min: 0.0,
            error_abs_sum: 0.0,
        }
    }

    fn torus(center: [f64; 3], axis: [f64; 3], r_major: f64, r_minor: f64) -> ToroidalHypothesis {
        ToroidalHypothesis {
            center,
            axis_direction: axis,
            major_radius: r_major,
            minor_radius: r_minor,
            convex: true,
            faces: Vec::new(),
            vertices: Vec::new(),
            error_max: 0.0,
            centroid_error_max: 0.0,
            error_abs_sum: 0.0,
        }
    }

    fn cylinder(apex: [f64; 3], axis: [f64; 3], radius: f64) -> ConicalHypothesis {
        ConicalHypothesis {
            apex,
            axis_direction: axis,
            half_angle: 0.0,
            radius_at_apex: radius,
            convex: true,
            faces: Vec::new(),
            vertices: Vec::new(),
            error_max: 0.0,
            centroid_error_max: 0.0,
            error_abs_sum: 0.0,
        }
    }

    /// A fillet torus dying into the floor it is tangent to: the junction is
    /// the major-radius circle, so the side test must reduce to ρ − R.
    #[test]
    fn plane_torus_locus_is_the_major_radius_circle() {
        let p = plane([0.0, 0.0, 1.0], 0.0);
        let t = torus([0.0, 0.0, 1.0], [0.0, 0.0, 1.0], 5.0, 1.0);
        let locus = build_locus(Surf::Planar(&p), Surf::Toroidal(&t), &[5.0, 0.0, 0.0]).unwrap();
        // Zero on the tangency circle, and the magnitude is the radial offset.
        assert!(locus.side(&[5.0, 0.0, 0.0]).abs() < 1e-12);
        assert!(locus.side(&[0.0, 5.0, 0.0]).abs() < 1e-12);
        assert!((locus.side(&[6.0, 0.0, 0.0]).abs() - 1.0).abs() < 1e-12);
        // The two supports land on opposite sides.
        let inside = locus.side(&[4.0, 0.0, 0.2]);
        let outside = locus.side(&[7.0, 0.0, 0.0]);
        assert!(inside * outside < 0.0);
    }

    /// Displacing the plane along its own normal must not move the locus:
    /// this is the amplification the minimum-distance locus exists to avoid.
    #[test]
    fn plane_torus_locus_is_immune_to_plane_offset_error() {
        let t = torus([0.0, 0.0, 1.0], [0.0, 0.0, 1.0], 5.0, 1.0);
        let exact = build_locus(
            Surf::Planar(&plane([0.0, 0.0, 1.0], 0.0)),
            Surf::Toroidal(&t),
            &[5.0, 0.0, 0.0],
        )
        .unwrap();
        // 4.7e-4 mm of plane error is what produced the recorded 4.4e-2 mm
        // edge error through the intersection; √(2rδ) is here ≈ 3.1e-2.
        let off = build_locus(
            Surf::Planar(&plane([0.0, 0.0, 1.0], 4.7e-4)),
            Surf::Toroidal(&t),
            &[5.0, 0.0, 0.0],
        )
        .unwrap();
        let probe = [5.3, 0.0, 0.0];
        assert!((exact.side(&probe) - off.side(&probe)).abs() < 1e-12);
    }

    /// A fillet torus continuing into the coaxial cylinder it is tangent to:
    /// the junction is at the tube's centre height, where the radial
    /// coordinate is stationary and the axial one is the informative one.
    #[test]
    fn cylinder_torus_locus_is_the_tube_centre_circle() {
        let c = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 6.0);
        let t = torus([0.0, 0.0, 3.0], [0.0, 0.0, 1.0], 5.0, 1.0);
        let locus =
            build_locus(Surf::Conical(&c), Surf::Toroidal(&t), &[6.0, 0.0, 3.0]).unwrap();
        assert!(locus.side(&[6.0, 0.0, 3.0]).abs() < 1e-12);
        // Moving along the axis is what crosses this junction.
        assert!((locus.side(&[6.0, 0.0, 4.0]).abs() - 1.0).abs() < 1e-12);
        assert!(locus.side(&[6.0, 0.0, 2.0]) * locus.side(&[6.0, 0.0, 4.0]) < 0.0);
    }

    /// A flat milled on a shaft: the junction is a ruling, and the side test
    /// runs across it in the shared tangent plane.
    #[test]
    fn plane_cylinder_ruling_locus() {
        let c = cylinder([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2.0);
        let p = plane([1.0, 0.0, 0.0], 2.0);
        let locus = build_locus(Surf::Planar(&p), Surf::Conical(&c), &[2.0, 0.0, 0.0]).unwrap();
        assert_eq!(locus.kind(), "line");
        // Zero all along the ruling, and sign-separating across it.
        assert!(locus.side(&[2.0, 0.0, 0.0]).abs() < 1e-12);
        assert!(locus.side(&[2.0, 0.0, 9.0]).abs() < 1e-12);
        assert!(locus.side(&[2.0, 1.0, 0.0]) * locus.side(&[2.0, -1.0, 0.0]) < 0.0);
    }

    #[test]
    fn non_coaxial_and_freeform_pairs_are_rejected() {
        let t = torus([0.0, 0.0, 1.0], [0.0, 0.0, 1.0], 5.0, 1.0);
        let tilted = torus([0.0, 0.0, 1.0], [0.0, 0.7, 0.7], 5.0, 1.0);
        assert_eq!(
            build_locus(Surf::Toroidal(&t), Surf::Toroidal(&tilted), &[5.0, 0.0, 0.0]),
            Err(LocusFailure::NotCoaxial)
        );
        assert_eq!(
            build_locus(Surf::Toroidal(&t), Surf::BSpline, &[5.0, 0.0, 0.0]),
            Err(LocusFailure::Freeform)
        );
        let p1 = plane([0.0, 0.0, 1.0], 0.0);
        let p2 = plane([0.0, 0.0, 1.0], 3.0);
        assert_eq!(
            build_locus(Surf::Planar(&p1), Surf::Planar(&p2), &[0.0, 0.0, 0.0]),
            Err(LocusFailure::NoJunctionExists)
        );
    }

    #[test]
    fn tangency_matches_normal_agreement() {
        let t = torus([0.0, 0.0, 1.0], [0.0, 0.0, 1.0], 5.0, 1.0);
        let floor = plane([0.0, 0.0, -1.0], 0.0);
        // At the tangency circle the torus normal points straight down.
        assert!(is_tangent_pair(
            Surf::Planar(&floor),
            Surf::Toroidal(&t),
            &[[5.0, 0.0, 0.0], [0.0, 5.0, 0.0]],
        ));
        // A quarter turn away it points straight out, and the pair is not.
        assert!(!is_tangent_pair(
            Surf::Planar(&floor),
            Surf::Toroidal(&t),
            &[[6.0, 0.0, 1.0]],
        ));
    }
}
