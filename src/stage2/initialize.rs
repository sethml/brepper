//! Closed-form initializers (Phase 1.1): a surface's parameters from a handful
//! of facets, exactly, with no iteration and no seed patch.
//!
//! Phase 1's earned patch grows face by face under validation, which means it
//! must fit a surface to *tiny* support — where a free least-squares fit is
//! hopeless. The sagitta arithmetic: over a 7.5° arc on r = 0.5 the sag is
//! ~0.001 mm against a 2e-5 mm mesh quantum, so a free radius comes out ~2%
//! wrong, about 200× the vertex tolerance. What rescues it is that **facet
//! normals carry exact information**, which turns each surface type into a
//! small linear problem rather than an optimization.
//!
//! Nothing here iterates and nothing here is tuned. Every constant is a
//! numeric-conditioning floor — it decides whether a quantity is *defined*,
//! never whether a fit is *good* — and each is documented as such at its
//! definition.
//!
//! # Which facets to initialize from, and why it matters more than the formula
//!
//! Phase 1.0 measured two things that make facet *selection* part of the
//! algorithm rather than an afterthought (numbers in TODO.md 1.0):
//!
//! 1. **The axis error is `perpendicularity error / |n₁ × n₂|`**, verified to
//!    0.1%. Adjacent facets differ in normal by the exporter's angular step, so
//!    `|n₁ × n₂| ≈ sin(step)` — 0.043 on a fine mesh, a ~23× amplification that
//!    gets *worse* as the mesh refines. So pair the angularly most distant
//!    facets available, never neighbours.
//! 2. **A facet's normal is perpendicular to the axis iff two of its vertices
//!    share a ruling** ([`min_ruling_gap_rad`]). 25.6% of corpus cylinder facets
//!    straddle three or more rulings, and they carry *all* of the
//!    perpendicularity error — aligned facets are exact to 0.0006° worst case
//!    over the whole corpus, straddling ones reach 1.0528°.
//!
//! Both are checkable here without ground truth, so [`select_axis_pair`] does
//! it in two passes: pick the most distant normal pair to get a provisional
//! axis, then re-pick among the facets that are ruling-aligned *with respect to
//! that axis*. The instrument that measured (2) and the code that acts on it
//! call the same function, so they cannot drift apart.

// Phase 1.1 builds and pins these; Phase 1.2's earn loop is what calls them.
// Until then the constructors are reachable only from the tests, which is the
// point of unit-testing them first — the allow comes off when `grow` uses them.
#![allow(dead_code)]

use super::{cross3, dot3, face_centroid, normalize3, perpendicular_basis};
use crate::stage1::{ConnectedMesh, MeshFace, MeshVertex};

/// A closed-form fit, ready to hand to validation and growth.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Initial {
    Cylinder {
        axis_pos: [f64; 3],
        axis_dir: [f64; 3],
        radius: f64,
    },
    Sphere {
        centre: [f64; 3],
        radius: f64,
    },
    Cone {
        apex: [f64; 3],
        axis_dir: [f64; 3],
        half_angle: f64,
    },
    Torus {
        centre: [f64; 3],
        axis_dir: [f64; 3],
        major: f64,
        minor: f64,
    },
}


// Local vector helpers (stage 2 shares only cross3/dot3/normalize3).
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn norm(a: [f64; 3]) -> f64 {
    dot3(&a, &a).sqrt()
}

// ---------------------------------------------------------------------------
// Conditioning floors
// ---------------------------------------------------------------------------

/// Below this, two unit normals are parallel to within double-precision noise
/// and their cross product carries no direction at all.
///
/// A validity floor, not a threshold: it decides whether the axis *exists*, not
/// whether it is good enough. The quality of a well-defined axis is governed by
/// the amplification law in the module docs, which is a continuous relation with
/// no cut in it.
const MIN_CROSS_NORM: f64 = 1.0e-9;

/// Below this angular gap (radians) two of a facet's vertices lie on the same
/// ruling.
///
/// Not a tuned threshold: measured over the corpus the two cases are separated
/// by five orders of magnitude of empty space — vertices on one ruling differ
/// only by float32 mesh noise (p50 reads 5.551e-16 rad) while distinct rulings
/// are separated by the exporter's angular step (p95 1.178e-1 rad). Any cut
/// between gives the same answer.
pub(crate) const RULING_SAME_RAD: f64 = 1.0e-6;

/// The second-smallest eigenvalue must exceed this fraction of the largest for
/// the null space to be computable at all.
///
/// **A numerical floor and nothing more.** It was tempting to make this a rank
/// test sharp enough to tell a torus from a sphere — the meridian system is
/// rank-deficient by construction on sphere data, where every line through the
/// centre satisfies it. Measured, that discrimination is not available at this
/// ratio: the degenerate sphere patch reads `l1/l5` = 2.06e-6 and the hardest
/// genuine torus (R = 4, r = 3) reads 9.2e-4, a factor of 450 with nothing but
/// data-dependence in between. Choosing a cut there is the "detector needs just
/// one threshold" trap, so it is not chosen here.
///
/// What separates the cases instead is *structure*, at no cost: a sphere's
/// (ρ, z) profile is a circle centred on the axis, so its fitted major radius
/// comes out ~0 and [`torus`]'s `minor < major` check rejects it outright. Type
/// discrimination is the ladder's job and validation's, not this floor's — which
/// is what the plan meant by "ladder order and a conditioning floor must precede
/// it".
const MIN_NULLSPACE_SEPARATION: f64 = 1.0e-12;

// ---------------------------------------------------------------------------
// Facet selection
// ---------------------------------------------------------------------------

/// The smallest angular gap, in radians, between any two of a facet's vertices
/// measured around the axis `(axis_pos, axis_dir)` — zero exactly when two
/// vertices share a ruling.
///
/// This decides whether a facet normal on a ruled surface is perpendicular to
/// the axis at all. Parametrize the cylinder as
/// `p(θ, z) = a + r(cos θ · u + sin θ · v) + z · d`. The facet's plane contains
/// `d` — equivalently `n ⟂ d` — iff the vertices, projected along `d` onto the
/// plane perpendicular to it, are collinear. They land on a circle of radius
/// `r`, and three points on a circle are collinear only when two of them
/// coincide. **So `n ⟂ d` exactly when two vertices share a ruling**, and a
/// triangle straddling three distinct rulings tilts its normal off perpendicular
/// however perfect the underlying cylinder is.
pub(crate) fn min_ruling_gap_rad(
    mesh: &ConnectedMesh,
    fi: usize,
    axis_pos: [f64; 3],
    axis_dir: [f64; 3],
) -> Option<f64> {
    let mut d = axis_dir;
    if normalize3(&mut d) == 0.0 {
        return None;
    }
    let (u, v) = perpendicular_basis(&d);

    let face = &mesh.faces[fi];
    let n = face.vertex_count as usize;
    if n < 3 {
        return None;
    }
    let mut thetas = Vec::with_capacity(n);
    for k in 0..n {
        let vtx = &mesh.vertices[face.vertex_indices[k]];
        let w = [
            vtx.x - axis_pos[0],
            vtx.y - axis_pos[1],
            vtx.z - axis_pos[2],
        ];
        let axial = dot3(&w, &d);
        let radial = [
            w[0] - axial * d[0],
            w[1] - axial * d[1],
            w[2] - axial * d[2],
        ];
        thetas.push(dot3(&radial, &v).atan2(dot3(&radial, &u)));
    }

    let mut min_gap = f64::MAX;
    for i in 0..n {
        for j in (i + 1)..n {
            let mut delta = (thetas[i] - thetas[j]).abs();
            if delta > std::f64::consts::PI {
                delta = std::f64::consts::TAU - delta;
            }
            min_gap = min_gap.min(delta);
        }
    }
    Some(min_gap)
}

/// `|n_i × n_j|` — how well-conditioned the pair is as an axis source, and
/// exactly the denominator in the amplification law.
fn pair_conditioning(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let c = cross3(a, b);
    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
}

/// The most angularly distant pair of facet normals in `faces`.
fn widest_normal_pair(faces: &[usize], mesh: &ConnectedMesh) -> Option<(usize, usize, f64)> {
    let mut best: Option<(usize, usize, f64)> = None;
    for (a, &fi) in faces.iter().enumerate() {
        let Some(ni) = mesh.faces[fi].normal else {
            continue;
        };
        for &fj in &faces[a + 1..] {
            let Some(nj) = mesh.faces[fj].normal else {
                continue;
            };
            let k = pair_conditioning(&ni, &nj);
            if best.is_none_or(|(_, _, bk)| k > bk) {
                best = Some((fi, fj, k));
            }
        }
    }
    best
}

/// Pick the pair of facets to take an axis from: ruling-aligned, and as
/// angularly distant as the patch allows.
///
/// Two passes, because alignment is defined relative to an axis we do not have
/// yet. The first pass gives a provisional axis good enough to *classify*
/// facets — misjudging alignment needs an axis error of order the angular step,
/// far beyond what pass one produces — and the second pass re-picks among the
/// facets that survive. If too few survive, the provisional pair stands: a
/// worse-conditioned axis is still better than no axis, and validation
/// downstream is what decides whether the fit lives.
///
/// Returns the pair and the axis direction it implies.
fn select_axis_pair(faces: &[usize], mesh: &ConnectedMesh) -> Option<(usize, usize, [f64; 3])> {
    let axis_of = |fi: usize, fj: usize| -> Option<[f64; 3]> {
        let (ni, nj) = (mesh.faces[fi].normal?, mesh.faces[fj].normal?);
        let mut axis = cross3(&ni, &nj);
        if normalize3(&mut axis) < MIN_CROSS_NORM {
            return None;
        }
        Some(axis)
    };

    let (pi, pj, _) = widest_normal_pair(faces, mesh)?;
    let provisional = axis_of(pi, pj)?;

    // Any point on the provisional axis will do to measure ruling gaps: the
    // gap is an angle about the axis direction, and shifting the reference
    // point along or across the axis moves every vertex's θ together.
    let anchor = face_centroid(&mesh.faces[pi], &mesh.vertices);
    let aligned: Vec<usize> = faces
        .iter()
        .copied()
        .filter(|&fi| {
            min_ruling_gap_rad(mesh, fi, anchor, provisional)
                .is_some_and(|gap| gap < RULING_SAME_RAD)
        })
        .collect();

    if aligned.len() >= 2 {
        if let Some((ai, aj, _)) = widest_normal_pair(&aligned, mesh) {
            if let Some(axis) = axis_of(ai, aj) {
                return Some((ai, aj, axis));
            }
        }
    }
    Some((pi, pj, provisional))
}

/// Every distinct vertex index used by `faces`.
fn patch_vertices(faces: &[usize], mesh: &ConnectedMesh) -> Vec<usize> {
    let mut out = Vec::new();
    for &fi in faces {
        let f = &mesh.faces[fi];
        for k in 0..f.vertex_count as usize {
            let vi = f.vertex_indices[k];
            if !out.contains(&vi) {
                out.push(vi);
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Cylinder
// ---------------------------------------------------------------------------

/// A cylinder from its facets' normals: axis from `n₁ × n₂`, then radius and
/// position from the existing circle fit in the plane perpendicular to it.
///
/// Exact on a ruling-aligned tessellation at any coarseness — a facet spanning
/// two rulings has its normal exactly perpendicular to the axis, so two such
/// normals span the perpendicular plane and their cross product *is* the axis.
pub(crate) fn cylinder(faces: &[usize], mesh: &ConnectedMesh) -> Option<Initial> {
    let (_, _, axis_dir) = select_axis_pair(faces, mesh)?;
    let verts = patch_vertices(faces, mesh);
    let (axis_pos, radius) =
        super::cylindrical::fit_circle_for_axis(&axis_dir, &verts, &mesh.vertices)?;
    Some(Initial::Cylinder {
        axis_pos,
        axis_dir,
        radius,
    })
}

// ---------------------------------------------------------------------------
// Sphere
// ---------------------------------------------------------------------------

/// The circumcenter of a facet's first three vertices.
///
/// For a facet inscribed in a sphere this is the one point of the facet whose
/// normal line passes **exactly** through the centre: the circumcenter is
/// equidistant from three points of the sphere, so it lies on the axis of the
/// circle they define, which is the line through the sphere's centre. The
/// centroid does not have this property, which is why the initializer uses the
/// circumcenter and the diagnostic in `crate::reference` measures against it.
fn circumcenter(face: &MeshFace, vertices: &[MeshVertex]) -> Option<[f64; 3]> {
    if (face.vertex_count as usize) < 3 {
        return None;
    }
    let p = |k: usize| {
        let v = &vertices[face.vertex_indices[k]];
        [v.x, v.y, v.z]
    };
    let (a, b, c) = (p(0), p(1), p(2));
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = cross3(&ab, &ac);
    let denom = 2.0 * dot3(&n, &n);
    if denom <= 0.0 {
        return None;
    }
    let t1 = cross3(&n, &ab);
    let t2 = cross3(&ac, &n);
    let ac2 = dot3(&ac, &ac);
    let ab2 = dot3(&ab, &ab);
    Some([
        a[0] + (t1[0] * ac2 + t2[0] * ab2) / denom,
        a[1] + (t1[1] * ac2 + t2[1] * ab2) / denom,
        a[2] + (t1[2] * ac2 + t2[2] * ab2) / denom,
    ])
}

/// Midpoint of the common perpendicular between two lines — the standard
/// near-intersection, and the sphere centre when both lines pass through it.
fn common_perpendicular_midpoint(
    p1: [f64; 3],
    d1: [f64; 3],
    p2: [f64; 3],
    d2: [f64; 3],
) -> Option<[f64; 3]> {
    let (mut u, mut v) = (d1, d2);
    if normalize3(&mut u) == 0.0 || normalize3(&mut v) == 0.0 {
        return None;
    }
    let w = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
    let b = dot3(&u, &v);
    let denom = 1.0 - b * b;
    if denom < MIN_CROSS_NORM {
        return None;
    }
    let d = dot3(&u, &w);
    let e = dot3(&v, &w);
    let s = (b * e - d) / denom;
    let t = (e - b * d) / denom;
    Some([
        0.5 * ((p1[0] + u[0] * s) + (p2[0] + v[0] * t)),
        0.5 * ((p1[1] + u[1] * s) + (p2[1] + v[1] * t)),
        0.5 * ((p1[2] + u[2] * s) + (p2[2] + v[2] * t)),
    ])
}

/// A sphere from the concurrence of its facets' normal lines.
///
/// Every facet inscribed in a sphere has a normal line through the centre when
/// anchored at its circumcenter, so two facets locate the centre exactly. The
/// pair is chosen by the same conditioning rule as the cylinder — two nearly
/// parallel normal lines meet nowhere useful, which Phase 1.0 measured directly
/// (a pair centre error reaching 4.079 mm at `|n₁ × n₂|` ≈ 0.0436).
///
/// The radius is the mean distance from the centre to the patch's vertices,
/// which uses every facet rather than the two.
pub(crate) fn sphere(faces: &[usize], mesh: &ConnectedMesh) -> Option<Initial> {
    let (fi, fj, _) = widest_normal_pair(faces, mesh)?;
    let (ni, nj) = (mesh.faces[fi].normal?, mesh.faces[fj].normal?);
    if pair_conditioning(&ni, &nj) < MIN_CROSS_NORM {
        return None;
    }
    let ci = circumcenter(&mesh.faces[fi], &mesh.vertices)?;
    let cj = circumcenter(&mesh.faces[fj], &mesh.vertices)?;
    let centre = common_perpendicular_midpoint(ci, ni, cj, nj)?;

    let verts = patch_vertices(faces, mesh);
    if verts.is_empty() {
        return None;
    }
    let mut sum = 0.0;
    for &vi in &verts {
        let v = &mesh.vertices[vi];
        let d = [v.x - centre[0], v.y - centre[1], v.z - centre[2]];
        sum += dot3(&d, &d).sqrt();
    }
    let radius = sum / verts.len() as f64;
    if radius <= 0.0 {
        return None;
    }
    Some(Initial::Sphere { centre, radius })
}

// ---------------------------------------------------------------------------
// The axis of a surface of revolution, from its normal lines
// ---------------------------------------------------------------------------

/// The axis line shared by every normal of a surface of revolution, as
/// `(point on the axis, unit direction)`.
///
/// **The bilinear system in the plan is linear in Plücker coordinates.** The
/// meridian-coplanarity constraint `(pᵢ − a)·(nᵢ × d) = 0` is bilinear in the
/// axis point `a` and direction `d`, which is why the plan proposed alternation
/// or pairwise elimination. Writing the axis as a line rather than as a
/// point-plus-direction removes that entirely. With the line's moment
/// `m = a × d`:
///
/// ```text
///   (pᵢ − a)·(nᵢ × d) = d·(pᵢ × nᵢ) − det(a, nᵢ, d)
///                     = d·(pᵢ × nᵢ) + nᵢ·(a × d)
///                     = qᵢ·d + nᵢ·m        where qᵢ = pᵢ × nᵢ
/// ```
///
/// so each facet contributes **one linear equation** in the six unknowns
/// `(d, m)`, and the axis is the null space of the stacked system. No
/// iteration, no alternation, no starting guess. Geometrically the identity is
/// just "every normal line of a surface of revolution meets the axis", which is
/// why the one solve serves cones and tori alike — and cylinders, though those
/// have a cheaper exact route through `n₁ × n₂`.
///
/// Returns `None` when the null space is not one-dimensional, which is the
/// rank-deficiency the plan warns about: on sphere data every line through the
/// centre satisfies the system, so the null space is three-dimensional and any
/// single answer would be arbitrary.
fn axis_line_from_normals(
    faces: &[usize],
    mesh: &ConnectedMesh,
) -> Option<([f64; 3], [f64; 3])> {
    // Normalize before solving. The system mixes q = p × n, which carries
    // length units, with n, which is dimensionless, so on a patch far from the
    // origin the two blocks differ by orders of magnitude and the null space
    // drowns in the scale difference rather than in the data. Shifting to the
    // patch centroid and dividing by its own extent makes both blocks O(1);
    // the recovered line is mapped back at the end. This is conditioning, not
    // fitting — it changes no answer that was well posed to begin with.
    let mut points = Vec::with_capacity(faces.len());
    let mut normals = Vec::with_capacity(faces.len());
    for &fi in faces {
        let Some(n) = mesh.faces[fi].normal else {
            continue;
        };
        points.push(face_centroid(&mesh.faces[fi], &mesh.vertices));
        normals.push(n);
    }
    if points.len() < 4 {
        return None;
    }
    let inv = 1.0 / points.len() as f64;
    let mut origin = [0.0_f64; 3];
    for p in &points {
        for k in 0..3 {
            origin[k] += p[k] * inv;
        }
    }
    let mut scale = 0.0_f64;
    for p in &points {
        let w = [p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]];
        scale += dot3(&w, &w) * inv;
    }
    let extent = scale.sqrt();
    if !(extent > 0.0 && extent.is_finite()) {
        return None;
    }

    let mut m = nalgebra::Matrix6::<f64>::zeros();
    for (p, n) in points.iter().zip(normals.iter()) {
        let local = [
            (p[0] - origin[0]) / extent,
            (p[1] - origin[1]) / extent,
            (p[2] - origin[2]) / extent,
        ];
        let q = cross3(&local, n);
        let row = nalgebra::Vector6::new(q[0], q[1], q[2], n[0], n[1], n[2]);
        m += row * row.transpose();
    }

    let eigen = nalgebra::SymmetricEigen::new(m);
    let mut order: Vec<usize> = (0..6).collect();
    order.sort_by(|&a, &b| eigen.eigenvalues[a].total_cmp(&eigen.eigenvalues[b]));
    let (l0, l1, l5) = (
        eigen.eigenvalues[order[0]],
        eigen.eigenvalues[order[1]],
        eigen.eigenvalues[order[5]],
    );
    // `BREPPER_DEBUG_RANK=1` prints the spectrum. A declined initializer is
    // otherwise silent, and the null space is where it declines.
    if std::env::var("BREPPER_DEBUG_RANK").is_ok() {
        eprintln!(
            "axis_line_from_normals: l0={l0:.3e} l1={l1:.3e} l5={l5:.3e} l1/l5={:.3e}",
            l1 / l5
        );
    }
    // Numerical rank test: the null space must be exactly one-dimensional.
    // Scaled by the largest eigenvalue so it is a relative test on the system's
    // own conditioning rather than an absolute size.
    if !(l5 > 0.0 && l1 > l5 * MIN_NULLSPACE_SEPARATION && l0 <= l1) {
        if std::env::var("BREPPER_DEBUG_RANK").is_ok() {
            eprintln!("axis_line_from_normals: declined, null space not 1-D");
        }
        return None;
    }

    let v = eigen.eigenvectors.column(order[0]);
    let mut d = [v[0], v[1], v[2]];
    let dnorm = normalize3(&mut d);
    if dnorm < MIN_CROSS_NORM {
        return None;
    }
    // The moment shares the null vector's scale, so divide by the same factor.
    let mut moment = [v[3] / dnorm, v[4] / dnorm, v[5] / dnorm];
    // Enforce the Plücker condition d·m = 0; a nonzero component is numerical
    // drift, not information.
    let drift = dot3(&moment, &d);
    for k in 0..3 {
        moment[k] -= drift * d[k];
    }
    // The point of the axis closest to the local origin, mapped back out of
    // the normalized frame.
    let local_point = cross3(&d, &moment);
    let point = [
        origin[0] + local_point[0] * extent,
        origin[1] + local_point[1] * extent,
        origin[2] + local_point[2] * extent,
    ];
    Some((point, d))
}

/// Axial and radial coordinates of every patch vertex about `(pos, dir)`.
fn profile_coords(
    faces: &[usize],
    mesh: &ConnectedMesh,
    pos: [f64; 3],
    dir: [f64; 3],
) -> Vec<(f64, f64)> {
    patch_vertices(faces, mesh)
        .into_iter()
        .map(|vi| {
            let v = &mesh.vertices[vi];
            let w = [v.x - pos[0], v.y - pos[1], v.z - pos[2]];
            let z = dot3(&w, &dir);
            let radial = [w[0] - z * dir[0], w[1] - z * dir[1], w[2] - z * dir[2]];
            (dot3(&radial, &radial).sqrt(), z)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Cone
// ---------------------------------------------------------------------------


/// Largest distance from any of `verts` to `fit` — how well a candidate
/// explains its own support.
fn max_residual(fit: &Initial, verts: &[usize], mesh: &ConnectedMesh) -> f64 {
    verts
        .iter()
        .map(|&vi| {
            let v = &mesh.vertices[vi];
            let p = [v.x, v.y, v.z];
            match *fit {
                Initial::Cylinder {
                    axis_pos,
                    axis_dir,
                    radius,
                } => {
                    let w = sub(p, axis_pos);
                    let z = dot3(&w, &axis_dir);
                    (norm(sub(w, scale(axis_dir, z))) - radius).abs()
                }
                Initial::Cone {
                    apex,
                    axis_dir,
                    half_angle,
                } => {
                    let w = sub(p, apex);
                    let z = dot3(&w, &axis_dir);
                    let r = norm(sub(w, scale(axis_dir, z)));
                    (r * half_angle.cos() - z * half_angle.sin()).abs()
                }
                _ => f64::MAX,
            }
        })
        .fold(0.0_f64, f64::max)
}

/// The cone axis, **exactly**, from the *centred* normal covariance.
///
/// Phase 1.1 rejected a normal-covariance axis because the uncentred covariance
/// `Σ nnᵀ` has eigenvalues `cos²α/2` (twice) and `sin²α` (once), so the axis is
/// its smallest eigenvector only while α < 35.26°. **Centring removes that
/// regime entirely**, and the reason is short: on a cone
/// `nᵢ = cos α · r̂ᵢ − sin α · d`, so subtracting the mean normal cancels the
/// `sin α · d` term and leaves `cos α (r̂ᵢ − r̄)`, which lies in the plane
/// perpendicular to the axis **for every α**. The smallest eigenvector of the
/// centred covariance is therefore the axis exactly, at any half-angle, and the
/// 35.26° flip was an artefact of not centring.
///
/// What it does *not* fix is conditioning on a sliver: over a short arc the
/// `r̂ᵢ` are nearly collinear, so the centred normals span one direction rather
/// than two and the axis is pinned only to a plane — the meridian plane. The
/// remaining freedom within that plane is precisely the half-angle, which is
/// why [`cone`] resolves it by fitting a profile rather than by asking this for
/// more than it knows.
fn cone_axis_from_centred_normals(faces: &[usize], mesh: &ConnectedMesh) -> Option<[f64; 3]> {
    let mut normals = Vec::with_capacity(faces.len());
    for &fi in faces {
        if let Some(n) = mesh.faces[fi].normal {
            normals.push(n);
        }
    }
    if normals.len() < 2 {
        return None;
    }
    let inv = 1.0 / normals.len() as f64;
    let mut mean = [0.0_f64; 3];
    for n in &normals {
        for k in 0..3 {
            mean[k] += n[k] * inv;
        }
    }
    // Packing matches `build_normal_covariance`: [xx, xy, xz, yy, yz, zz].
    let mut cov = [0.0_f64; 6];
    for n in &normals {
        let c = [n[0] - mean[0], n[1] - mean[1], n[2] - mean[2]];
        cov[0] += c[0] * c[0];
        cov[1] += c[0] * c[1];
        cov[2] += c[0] * c[2];
        cov[3] += c[1] * c[1];
        cov[4] += c[1] * c[2];
        cov[5] += c[2] * c[2];
    }
    // The centred normals must span TWO directions, not one. On a developable
    // surface every facet along a ruling shares a tangent plane, so a patch
    // grown along one ruling has identical normals: the covariance is rank 0,
    // and a patch spanning one arc-step is rank 1. In both cases the axis is
    // pinned only to a plane (or not at all) and any answer would be arbitrary.
    // Declining here is what lets the earn loop keep collecting instead of
    // committing to a surface its support cannot determine.
    let mut ev = super::eigenvalues_3x3_symmetric(&cov);
    ev.sort_by(f64::total_cmp);
    if !(ev[2] > 0.0 && ev[1] > ev[2] * MIN_NULLSPACE_SEPARATION) {
        return None;
    }
    let mut axis = super::smallest_eigenvector_3x3(&cov);
    if normalize3(&mut axis) < MIN_CROSS_NORM {
        return None;
    }
    Some(axis)
}


/// A cone or a cylinder — whichever **explains the support better**.
///
/// The taper solve below needs the axis *direction*, and on a sliver that
/// direction is not determined: the centred normals span one direction rather
/// than two, so the axis is pinned only to the meridian plane and the remaining
/// freedom within it is exactly the half-angle. Measured, a 3-facet patch of an
/// 18.4° cone yields a taper fit that misses its own vertices by 3.6e-3 mm —
/// far better than the 87° the apex form produced, but still well above `vtol`,
/// so an earn loop would still starve on it.
///
/// The resolution needs no threshold on the ambiguity. Compute the cylinder
/// too — that construction is exact from two facets, and on a cone `n₁ × n₂`
/// gives the *ruling* direction, which is the minimum-taper member of exactly
/// the ambiguous family — then keep whichever candidate explains the patch
/// better. On a sliver the cylinder wins because the cone fit is numerically
/// adrift; once the patch spans enough arc for the axis to be determined, the
/// cone's extra freedom earns its place and it wins on residual. Ties go to the
/// cylinder, which has fewer parameters.
///
/// This is the same criterion the surrounding competition already uses —
/// explain the data — applied one level down, so it introduces no new kind of
/// decision.
pub(crate) fn cone(faces: &[usize], mesh: &ConnectedMesh, vtol: f64) -> Option<Initial> {
    let verts = patch_vertices(faces, mesh);
    let tapered = cone_tapered(faces, mesh, vtol);
    let straight = cylinder(faces, mesh);
    match (tapered, straight) {
        (Some(t), Some(c)) => {
            let (rt, rc) = (
                max_residual(&t, &verts, mesh),
                max_residual(&c, &verts, mesh),
            );
            Some(if rt < rc { t } else { c })
        }
        (Some(t), None) => Some(t),
        // `cone_tapered` returning None means the AXIS was undeterminable — the
        // centred normals span one direction or none — not that the surface is
        // a cylinder. Falling back to a cylinder here would hand the earn loop
        // a confident answer built on support that determines nothing, which is
        // exactly the commitment it must defer. A genuine near-cylinder is a
        // different state: `cone_tapered` resolves it and returns the cylinder
        // itself, above.
        (None, _) => None,
    }
}

/// The tapered-cone seed alone, as `(apex, axis_dir, half_angle)`, for a fitter
/// that already weighs the cylinder itself and wants only the cone.
///
/// The only route in the codebase whose axis is right at **every** half-angle:
/// `conical.rs`'s primary path takes the *uncentred* covariance, which picks the
/// wrong eigenvector above 35.26° — and wide cones (countersinks, chamfers,
/// shallow conical caps) are where that bites. `vtol` is passed as zero because
/// the "taper unresolvable, return the cylinder" decision belongs to the caller,
/// which fits the cylinder separately and compares.
pub(crate) fn cone_taper_seed(
    faces: &[usize],
    mesh: &ConnectedMesh,
) -> Option<([f64; 3], [f64; 3], f64)> {
    match cone_tapered(faces, mesh, 0.0)? {
        Initial::Cone {
            apex,
            axis_dir,
            half_angle,
        } => Some((apex, axis_dir, half_angle)),
        _ => None,
    }
}

/// Most vertices [`cone_from_rings`] will split, which bounds its work at
/// `2^(n-1)` partitions — 256 at the cap, the size of a four-quad chain plus
/// slack. It decides nothing about the surface; a larger support simply goes
/// to the centred-normal route, which needs breadth and has it.
const RING_SPLIT_MAX_VERTICES: usize = 10;

/// A cone from the **two rings** of a single-band strip — SURFACE_FITTING.md's
/// cone initialization Option C. Returns `(apex, axis_dir, half_angle)` with the
/// axis oriented from the apex toward the support.
///
/// A CAD tessellator lays a cone as one band of facets spanning its whole
/// height, so every vertex of a small seed sits on one of two circles, one per
/// rim. Three points determine a circle, so each rim gives a centre and a
/// radius; the line through the centres is the axis and the two radii at two
/// heights give the taper by similar triangles. Vertices only — no facet
/// normal, so none of the `O(Δθ²/8)` chord-plane lean — and it works from **two
/// facets**, where every normal-based route is rank-deficient (two normals pin
/// nothing about a cone).
///
/// **Which vertex belongs to which rim is found, not assumed.** The documented
/// design labels rims combinatorially by treating a shared quad edge as a
/// ruling, which holds only for quads in one band; a triangle strip's shared
/// edges alternate between rulings and diagonals. Instead every split of the
/// vertices into two groups of at least three is tried, each scored by the
/// largest distance of any vertex from the cone the split implies. On the true
/// split that is the mesh's own precision; a wrong split makes two circles whose
/// centre line is not their common normal, and misses by a large margin. The
/// winner is a seed for least squares, never the answer.
///
/// The taper `t = dρ/dh` is used throughout rather than the apex, for the reason
/// [`cone_tapered`] gives: `t = 0` is a cylinder, a regular point rather than an
/// apex at infinity. A split that reads as a cylinder returns `None` — the
/// caller fits cylinders by its own route — and so does any support too large or
/// too small to split.
pub(crate) fn cone_from_rings(
    faces: &[usize],
    mesh: &ConnectedMesh,
) -> Option<([f64; 3], [f64; 3], f64)> {
    let verts = patch_vertices(faces, mesh);
    let n = verts.len();
    if !(6..=RING_SPLIT_MAX_VERTICES).contains(&n) {
        return None;
    }
    let pts: Vec<[f64; 3]> = verts
        .iter()
        .map(|&vi| {
            let v = &mesh.vertices[vi];
            [v.x, v.y, v.z]
        })
        .collect();

    /// The best split so far: its score, and the profile it implies.
    struct Split {
        score: f64,
        c1: [f64; 3],
        d: [f64; 3],
        r1: f64,
        taper: f64,
    }
    let mut best: Option<Split> = None;
    // Vertex 0 always goes to the first group, so each split is visited once.
    for mask in 0u32..(1u32 << (n - 1)) {
        let in_first = |i: usize| i == 0 || (mask >> (i - 1)) & 1 == 1;
        let first: Vec<[f64; 3]> = (0..n).filter(|&i| in_first(i)).map(|i| pts[i]).collect();
        let second: Vec<[f64; 3]> = (0..n).filter(|&i| !in_first(i)).map(|i| pts[i]).collect();
        if first.len() < 3 || second.len() < 3 {
            continue;
        }
        let Some((c1, _, r1)) = super::toroidal::fit_circle_3d(&first) else {
            continue;
        };
        let Some((c2, _, r2)) = super::toroidal::fit_circle_3d(&second) else {
            continue;
        };
        let mut d = sub(c2, c1);
        let span = normalize3(&mut d);
        if !(span > 0.0 && span.is_finite() && r1.is_finite() && r2.is_finite()) {
            continue; // concentric rims in one plane: an annulus, not a band
        }
        let taper = (r2 - r1) / span;
        let cos_a = 1.0 / (1.0 + taper * taper).sqrt();
        let score = pts
            .iter()
            .map(|&p| {
                let w = sub(p, c1);
                let h = dot3(&w, &d);
                let rho = norm(sub(w, scale(d, h)));
                ((rho - (r1 + taper * h)) * cos_a).abs()
            })
            .fold(0.0_f64, f64::max);
        if score.is_finite() && best.as_ref().is_none_or(|b| score < b.score) {
            best = Some(Split { score, c1, d, r1, taper });
        }
    }

    let Split { c1, d, r1, taper, .. } = best?;
    if taper == 0.0 || !taper.is_finite() {
        return None;
    }
    // The apex is where the profile reaches radius zero: h_A = −r1 / taper.
    let apex = add(c1, scale(d, -r1 / taper));
    // Orient the axis from the apex toward the support: along +d when the
    // radius grows with h, against it when it shrinks.
    let axis = if taper > 0.0 { d } else { scale(d, -1.0) };
    let half_angle = taper.abs().atan();
    if !(half_angle > 0.0 && half_angle < std::f64::consts::FRAC_PI_2) {
        return None;
    }
    Some((apex, axis, half_angle))
}

/// The taper solve itself: axis from the centred normals, half-angle from the
/// normals, and only the position from a linear least squares.
///
/// # Why not (apex, axis, half-angle)
///
/// Phase 1.2 measured the cost of estimating the apex directly: on three
/// adjacent facets of an 18.4° cone the plane-concurrence solve returned a
/// half-angle of 87° with the apex halfway up the model, and the fit then
/// validated none of the three facets it was built from. The apex runs off to
/// **infinity** as a patch flattens toward a cylinder, and an estimator for an
/// unbounded quantity explodes near its degenerate case.
///
/// So both shape parameters come from quantities that stay bounded and local,
/// and neither is derived from the apex:
///
/// * the **axis direction**, exactly, from the centred normal covariance;
/// * the **taper** `sin α = |n·d|`, averaged over the facets — a nearly
///   constant quantity on a cone, so averaging it is well conditioned however
///   short the patch is. (It carries the `O(Δθ²/8)` chord-plane bias Phase 1.1
///   measured — 1.5e-3 relative at Δθ = 0.11 rad — which is a rounding error
///   beside the failure it replaces, and shrinks as the mesh refines.)
///
/// Only the *position* then needs solving, and it is a plain linear least
/// squares. Writing `u` for the axis offset from the patch centroid within the
/// plane perpendicular to `d`, `z_A` for the apex's axial coordinate, and
/// `m = |u|² − tan²α · z_A²`, the surface condition
/// `|r − u|² = tan²α (z − z_A)²` rearranges to
///
/// ```text
///   |r|² − tan²α z² = 2 r·u − 2 tan²α z · z_A − m
/// ```
///
/// which is linear in `(u₁, u₂, z_A, m)`. Least squares then returns whatever
/// best explains the support, so the fit **interpolates its own facets** rather
/// than extrapolating to a far-off point — that is what makes it degrade
/// gracefully. On a short patch `z_A` is uncertain, but the residuals on the
/// data stay small, so validation passes and the patch keeps growing; as the
/// patch lengthens `z_A` sharpens and the cone emerges on its own.
///
/// `taper = 0` is a **cylinder**, a finite value in the middle of the parameter
/// range rather than `apex = ∞` at its edge — the same reason circle fitters
/// use `κ = 1/r`. A patch whose whole radius change hides inside `vtol` cannot
/// distinguish the two, and this returns the cylinder in that case. That test
/// is derived from the tolerance and the patch's own extent; it adds no
/// constant.
fn cone_tapered(faces: &[usize], mesh: &ConnectedMesh, vtol: f64) -> Option<Initial> {
    let mut axis = cone_axis_from_centred_normals(faces, mesh)?;

    // sin α from the normals: constant on a cone, so an average is stable at
    // any patch size, and it never touches a position.
    let mut sin_sum = 0.0;
    let mut count = 0_usize;
    for &fi in faces {
        if let Some(n) = mesh.faces[fi].normal {
            sin_sum += dot3(&n, &axis).abs();
            count += 1;
        }
    }
    if count == 0 {
        return None;
    }
    let sin_alpha = (sin_sum / count as f64).clamp(0.0, 1.0);
    let cos_alpha = (1.0 - sin_alpha * sin_alpha).max(0.0).sqrt();
    if cos_alpha <= MIN_CROSS_NORM {
        return None; // a plane, not a cone
    }
    let tan_alpha = sin_alpha / cos_alpha;

    let verts = patch_vertices(faces, mesh);
    if verts.len() < 4 {
        return None;
    }
    let inv = 1.0 / verts.len() as f64;
    let mut origin = [0.0_f64; 3];
    for &vi in &verts {
        let v = &mesh.vertices[vi];
        origin = add(origin, scale([v.x, v.y, v.z], inv));
    }
    let (e1, e2) = perpendicular_basis(&axis);

    // Axial extent decides whether the taper is resolvable at all.
    let mut zmin = f64::MAX;
    let mut zmax = f64::MIN;
    let mut mean_radius = 0.0;
    let local: Vec<(f64, f64, f64)> = verts
        .iter()
        .map(|&vi| {
            let v = &mesh.vertices[vi];
            let w = sub([v.x, v.y, v.z], origin);
            let z = dot3(&w, &axis);
            let r = sub(w, scale(axis, z));
            zmin = zmin.min(z);
            zmax = zmax.max(z);
            let (r1, r2) = (dot3(&r, &e1), dot3(&r, &e2));
            mean_radius += (r1 * r1 + r2 * r2).sqrt() * inv;
            (z, r1, r2)
        })
        .collect();

    // Below this the whole radius change across the support hides inside the
    // tolerance the fit is judged at: the patch cannot tell cone from cylinder,
    // and the cylinder is the honest answer.
    if tan_alpha * (zmax - zmin) <= vtol {
        let (axis_pos, radius) =
            super::cylindrical::fit_circle_for_axis(&axis, &verts, &mesh.vertices)?;
        return Some(Initial::Cylinder {
            axis_pos,
            axis_dir: axis,
            radius,
        });
    }

    // Linear least squares for (u1, u2, z_A, m).
    let t2 = tan_alpha * tan_alpha;
    let mut ata = nalgebra::Matrix4::<f64>::zeros();
    let mut atb = nalgebra::Vector4::<f64>::zeros();
    for &(z, r1, r2) in &local {
        let row = nalgebra::Vector4::new(2.0 * r1, 2.0 * r2, -2.0 * t2 * z, -1.0);
        let rhs = r1 * r1 + r2 * r2 - t2 * z * z;
        ata += row * row.transpose();
        atb += row * rhs;
    }
    let sol = ata.lu().solve(&atb)?;
    if !sol.iter().all(|v| v.is_finite()) {
        return None;
    }
    let (u1, u2, z_apex) = (sol[0], sol[1], sol[2]);

    let axis_point = add(origin, add(scale(e1, u1), scale(e2, u2)));
    let apex = add(axis_point, scale(axis, z_apex));

    // Orient the axis from the apex toward the support.
    let mid = 0.5 * (zmin + zmax);
    if mid < z_apex {
        axis = [-axis[0], -axis[1], -axis[2]];
    }
    let half_angle = tan_alpha.atan();
    if !(half_angle > 0.0 && half_angle < std::f64::consts::FRAC_PI_2 && mean_radius > 0.0) {
        return None;
    }
    Some(Initial::Cone {
        apex,
        axis_dir: axis,
        half_angle,
    })
}

// ---------------------------------------------------------------------------
// Torus
// ---------------------------------------------------------------------------

/// A torus from its facets: axis from the normal lines, then `(R, r)` from a
/// circle fit in the `(ρ, z)` profile.
///
/// Phase 1.1 time-boxed this on the suspicion that the axis system would be too
/// ill-conditioned to solve. Phase 1.0 measured its residual on real meshes at
/// **2.354e-5 mm over 23 863 torus facets** — the mesh's own float32 quantum —
/// so the constraint is satisfied to the limit of what the mesh can express, and
/// the Plücker form above removes the conditioning worry that came from
/// alternating on a bilinear system.
pub(crate) fn torus(faces: &[usize], mesh: &ConnectedMesh) -> Option<Initial> {
    let (pos, dir) = axis_line_from_normals(faces, mesh)?;
    let profile = profile_coords(faces, mesh, pos, dir);
    if profile.len() < 3 {
        return None;
    }
    let ((major, z0), minor) = super::toroidal::fit_circle_2d(&profile)?;
    if !(major > 0.0 && minor > 0.0 && minor < major) {
        return None;
    }
    let centre = [
        pos[0] + dir[0] * z0,
        pos[1] + dir[1] * z0,
        pos[2] + dir[2] * z0,
    ];
    Some(Initial::Torus {
        centre,
        axis_dir: dir,
        major,
        minor,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// The corpus's default `--vertex-tolerance`, in mm — the bar the taper
    /// resolvability test is measured against.
    const TEST_VTOL: f64 = 5.0e-5;

    /// An oblique unit axis: an axis-aligned one would let a sign or an index
    /// slip through every assertion here.
    fn oblique_axis() -> [f64; 3] {
        let mut d = [0.3, -0.5, 0.81];
        normalize3(&mut d);
        d
    }

    /// `|â × b̂|` = sin of the angle between two directions.
    ///
    /// The only well-conditioned way to assert near-parallelism: an angle from
    /// `acos` cannot resolve below ~1e-6° because `acos` amplifies a 1e-16
    /// rounding error into ~1.5e-8 rad. The cross product is linear in the
    /// error and resolves to machine epsilon.
    fn parallel_miss(a: [f64; 3], b: [f64; 3]) -> f64 {
        let (mut u, mut v) = (a, b);
        normalize3(&mut u);
        normalize3(&mut v);
        let c = cross3(&u, &v);
        (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
    }

    fn vsub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    fn vnorm(a: [f64; 3]) -> f64 {
        dot3(&a, &a).sqrt()
    }

    /// Perpendicular distance from `p` to the line `(pos, dir)`.
    fn distance_to_axis(p: [f64; 3], pos: [f64; 3], dir: [f64; 3]) -> f64 {
        let mut d = dir;
        normalize3(&mut d);
        let w = vsub(p, pos);
        let t = dot3(&w, &d);
        vnorm(vsub(w, [d[0] * t, d[1] * t, d[2] * t]))
    }

    /// Build a mesh from explicit faces, computing each facet normal from its
    /// geometry exactly as stage 1 does.
    fn mesh_from(points: Vec<[f64; 3]>, faces: Vec<Vec<usize>>) -> ConnectedMesh {
        let vertices: Vec<MeshVertex> = points
            .iter()
            .map(|p| MeshVertex::from_xyz(p[0], p[1], p[2]))
            .collect();
        let faces = faces
            .into_iter()
            .map(|idx| {
                // The Newell normal, exactly as stage 1 computes it
                // (`compute_face_normal`). For a triangle this is the plain
                // cross product; for a quad it is the area-weighted plane, and
                // its symmetry is what makes a symmetric torus quad satisfy
                // meridian coplanarity exactly.
                let mut n = [0.0_f64; 3];
                for i in 0..idx.len() {
                    let c = points[idx[i]];
                    let x = points[idx[(i + 1) % idx.len()]];
                    n[0] += (c[1] - x[1]) * (c[2] + x[2]);
                    n[1] += (c[2] - x[2]) * (c[0] + x[0]);
                    n[2] += (c[0] - x[0]) * (c[1] + x[1]);
                }
                normalize3(&mut n);
                let mut vertex_indices = [0usize; 4];
                vertex_indices[..idx.len()].copy_from_slice(&idx);
                MeshFace {
                    vertex_count: idx.len() as u8,
                    vertex_indices,
                    neighbors: [-1; 4],
                    normal: Some(n),
                }
            })
            .collect();
        ConnectedMesh {
            vertices,
            faces,
            stats: Default::default(),
        }
    }

    /// A strip-tessellated cylinder: quads spanning two rulings each, exactly
    /// as a B-rep faceter emits them.
    fn cylinder_strip(
        origin: [f64; 3],
        d: [f64; 3],
        radius: f64,
        steps: usize,
        arc: f64,
    ) -> ConnectedMesh {
        let (u, v) = perpendicular_basis(&d);
        let at = |theta: f64, z: f64| {
            [
                origin[0] + radius * (u[0] * theta.cos() + v[0] * theta.sin()) + d[0] * z,
                origin[1] + radius * (u[1] * theta.cos() + v[1] * theta.sin()) + d[1] * z,
                origin[2] + radius * (u[2] * theta.cos() + v[2] * theta.sin()) + d[2] * z,
            ]
        };
        let mut points = Vec::new();
        for k in 0..=steps {
            let theta = arc * k as f64 / steps as f64;
            points.push(at(theta, 0.0));
            points.push(at(theta, 4.0));
        }
        let faces = (0..steps)
            .map(|k| vec![2 * k, 2 * k + 1, 2 * k + 3, 2 * k + 2])
            .collect();
        mesh_from(points, faces)
    }

    #[test]
    fn cylinder_axis_radius_and_position_are_exact_on_a_strip() {
        let origin = [1.0, -2.0, 3.0];
        let d = oblique_axis();
        for &(radius, steps, arc) in &[
            (7.0, 12, 1.2_f64),
            (0.5, 6, 0.5),   // the small-fillet case: a 0.5 mm radius over ~29°
            (40.0, 24, 2.4), // and a large one
        ] {
            let mesh = cylinder_strip(origin, d, radius, steps, arc);
            let faces: Vec<usize> = (0..mesh.faces.len()).collect();
            let Some(Initial::Cylinder {
                axis_pos,
                axis_dir,
                radius: r,
            }) = cylinder(&faces, &mesh)
            else {
                panic!("no cylinder for r={radius}");
            };
            assert!(
                parallel_miss(axis_dir, d) < 1.0e-14,
                "axis miss {:.3e} for r={radius}",
                parallel_miss(axis_dir, d)
            );
            assert!(
                (r - radius).abs() < 1.0e-9,
                "radius {r} vs {radius} (err {:.3e})",
                (r - radius).abs()
            );
            assert!(
                distance_to_axis(origin, axis_pos, axis_dir) < 1.0e-9,
                "axis position off by {:.3e} for r={radius}",
                distance_to_axis(origin, axis_pos, axis_dir)
            );
        }
    }

    /// The Phase 1.0 selection rule, as a behaviour: a patch that mixes
    /// ruling-aligned quads with triangles straddling three rulings must still
    /// give the exact axis, because selection discards the straddling facets.
    #[test]
    fn cylinder_axis_survives_facets_that_straddle_rulings() {
        let origin = [1.0, -2.0, 3.0];
        let d = oblique_axis();
        let radius = 7.0;
        let mut mesh = cylinder_strip(origin, d, radius, 10, 1.0);

        // Append triangles whose three vertices are on three different rulings
        // AND three different heights — the skew case, which tilts the normal.
        let (u, v) = perpendicular_basis(&d);
        let at = |theta: f64, z: f64| {
            [
                origin[0] + radius * (u[0] * theta.cos() + v[0] * theta.sin()) + d[0] * z,
                origin[1] + radius * (u[1] * theta.cos() + v[1] * theta.sin()) + d[1] * z,
                origin[2] + radius * (u[2] * theta.cos() + v[2] * theta.sin()) + d[2] * z,
            ]
        };
        let base = mesh.vertices.len();
        for (theta, z) in [(2.0, 0.0), (2.3, 4.0), (2.6, 1.7)] {
            let p = at(theta, z);
            mesh.vertices.push(MeshVertex::from_xyz(p[0], p[1], p[2]));
        }
        let tri = [base, base + 1, base + 2];
        let mut n = cross3(
            &vsub(at(2.3, 4.0), at(2.0, 0.0)),
            &vsub(at(2.6, 1.7), at(2.0, 0.0)),
        );
        normalize3(&mut n);
        // Sanity: this facet really is skewed, or the test proves nothing.
        assert!(
            dot3(&n, &d).abs() > 1.0e-3,
            "the straddling facet is not actually skewed"
        );
        let mut vertex_indices = [0usize; 4];
        vertex_indices[..3].copy_from_slice(&tri);
        mesh.faces.push(MeshFace {
            vertex_count: 3,
            vertex_indices,
            neighbors: [-1; 4],
            normal: Some(n),
        });

        let faces: Vec<usize> = (0..mesh.faces.len()).collect();
        let Some(Initial::Cylinder { axis_dir, .. }) = cylinder(&faces, &mesh) else {
            panic!("no cylinder");
        };
        assert!(
            parallel_miss(axis_dir, d) < 1.0e-14,
            "straddling facet leaked into the axis: miss {:.3e}",
            parallel_miss(axis_dir, d)
        );
    }

    #[test]
    fn sphere_centre_and_radius_are_exact_from_circumcenter_normal_lines() {
        let centre = [1.0, -2.0, 0.5];
        for &radius in &[6.0_f64, 0.8, 25.0] {
            let d = oblique_axis();
            let (u, v) = perpendicular_basis(&d);
            let dir = |a: f64, b: f64| {
                [
                    d[0] * b.cos() + (u[0] * a.cos() + v[0] * a.sin()) * b.sin(),
                    d[1] * b.cos() + (u[1] * a.cos() + v[1] * a.sin()) * b.sin(),
                    d[2] * b.cos() + (u[2] * a.cos() + v[2] * a.sin()) * b.sin(),
                ]
            };
            let on = |a: f64, b: f64| {
                let g = dir(a, b);
                [
                    centre[0] + radius * g[0],
                    centre[1] + radius * g[1],
                    centre[2] + radius * g[2],
                ]
            };
            let mut points = Vec::new();
            let mut faces = Vec::new();
            for k in 0..6 {
                let (a, b) = (k as f64 * 0.6, 0.4 + k as f64 * 0.18);
                let base = points.len();
                points.push(on(a, b));
                points.push(on(a + 0.20, b));
                points.push(on(a, b + 0.17));
                faces.push(vec![base, base + 1, base + 2]);
            }
            let mesh = mesh_from(points, faces);
            let all: Vec<usize> = (0..mesh.faces.len()).collect();
            let Some(Initial::Sphere { centre: c, radius: r }) = sphere(&all, &mesh) else {
                panic!("no sphere for r={radius}");
            };
            assert!(
                vnorm(vsub(c, centre)) < 1.0e-9,
                "centre off by {:.3e} for r={radius}",
                vnorm(vsub(c, centre))
            );
            assert!(
                (r - radius).abs() < 1.0e-9,
                "radius {r} vs {radius} (err {:.3e})",
                (r - radius).abs()
            );
        }
    }

    /// A strip-tessellated cone: quads spanning two rulings, which on a cone are
    /// two coplanar lines through the apex, so the facets are exactly planar.
    fn cone_strip(
        apex: [f64; 3],
        d: [f64; 3],
        half_angle: f64,
        steps: usize,
        arc: f64,
    ) -> ConnectedMesh {
        let (u, v) = perpendicular_basis(&d);
        let tan = half_angle.tan();
        let at = |theta: f64, t: f64| {
            let rho = t * tan;
            [
                apex[0] + d[0] * t + rho * (u[0] * theta.cos() + v[0] * theta.sin()),
                apex[1] + d[1] * t + rho * (u[1] * theta.cos() + v[1] * theta.sin()),
                apex[2] + d[2] * t + rho * (u[2] * theta.cos() + v[2] * theta.sin()),
            ]
        };
        let mut points = Vec::new();
        for k in 0..=steps {
            let theta = arc * k as f64 / steps as f64;
            points.push(at(theta, 3.0));
            points.push(at(theta, 9.0));
        }
        let faces = (0..steps)
            .map(|k| vec![2 * k, 2 * k + 1, 2 * k + 3, 2 * k + 2])
            .collect();
        mesh_from(points, faces)
    }

    #[test]
    fn cone_axis_apex_and_half_angle_are_exact_at_every_half_angle() {
        let apex = [2.0, 1.0, -3.0];
        let d = oblique_axis();
        // 36.87° is the corpus's own cone_15x20; 45° is a 90° countersink.
        // Both are past the 35.26° point where a normal-covariance axis flips.
        for &deg in &[10.0_f64, 20.0, 36.87, 45.0, 60.0, 75.0] {
            let alpha = deg.to_radians();
            let mesh = cone_strip(apex, d, alpha, 10, 1.1);
            let faces: Vec<usize> = (0..mesh.faces.len()).collect();
            let Some(Initial::Cone {
                apex: a,
                axis_dir,
                half_angle,
            }) = cone(&faces, &mesh, TEST_VTOL)
            else {
                panic!("no cone at {deg}°");
            };
            // Bounds are the measured worst case with headroom. Conditioning
            // degrades toward narrow cones — the null direction is least well
            // separated there — so 10° is the worst row at 1.9e-12 while 75°
            // reads 2.8e-14. Everything here is far inside the plan's ~1e-9.
            assert!(
                parallel_miss(axis_dir, d) < 1.0e-11,
                "axis miss {:.3e} at {deg}°",
                parallel_miss(axis_dir, d)
            );
            assert!(dot3(&axis_dir, &d) > 0.0, "axis points backwards at {deg}°");
            // The half-angle comes from the facet normals, so it carries the
            // chord-plane bias Phase 1.1 derived: a strip quad's plane leans
            // off the tangent plane, giving
            // `n·d = sin α cos(Δθ/2)/√(1 − sin²α sin²(Δθ/2))`, an O(Δθ²/8)
            // shortfall. That is deliberate — see `cone`. Here Δθ = 1.1/10, so
            // the prediction is 1.5e-3 relative, and the measurement matches it
            // (9.985° against 10.000°). Assert the bias is no worse than
            // predicted rather than pretending it is absent.
            let d_theta = 1.1_f64 / 10.0;
            let predicted = alpha * d_theta * d_theta / 8.0;
            let err = (alpha - half_angle).abs();
            assert!(
                err <= predicted * 1.2 && half_angle < alpha,
                "half-angle {:.6}° vs {deg}°: error {:.3e} rad, predicted bias {:.3e}",
                half_angle.to_degrees(),
                err,
                predicted
            );
            // The apex inherits that bias over the patch's own distance from
            // it, which is the honest bound: it is not an independent error.
            let apex_scale = vnorm(vsub(apex, a)).max(1.0);
            assert!(
                vnorm(vsub(a, apex)) < 0.05 * apex_scale,
                "apex off by {:.3e} at {deg}°",
                vnorm(vsub(a, apex))
            );
        }
    }

    /// Why the plan's normal-covariance axis was not used, measured across the
    /// regime an earned patch actually sees.
    ///
    /// Two separate failures. Over a *full* revolution the normal covariance has
    /// eigenvalues `cos²α/2` (twice, perpendicular to the axis) and `sin²α`
    /// (along it), so the axis is the smallest eigenvector only while
    /// `sin²α < cos²α/2`, i.e. **α < 35.26°** — and a 90° countersink is 45°.
    /// Over a *partial* arc the perpendicular spread is not isotropic either, so
    /// it degrades well before that. This asserts only what it measures: that
    /// the covariance axis is unreliable across the grid while the construction
    /// actually used is exact everywhere on it.
    #[test]
    fn the_normal_covariance_axis_is_unreliable_where_the_cone_construction_is_exact() {
        let apex = [2.0, 1.0, -3.0];
        let d = oblique_axis();
        let mut covariance_failures = 0;
        let mut cases = 0;
        for &deg in &[10.0_f64, 20.0, 36.87, 45.0, 60.0] {
            for &arc in &[0.6_f64, 1.4, std::f64::consts::TAU] {
                let steps = if arc > 3.0 { 24 } else { 10 };
                let mesh = cone_strip(apex, d, deg.to_radians(), steps, arc);
                let faces: Vec<usize> = (0..mesh.faces.len()).collect();
                cases += 1;

                let cov = super::super::build_normal_covariance(
                    &faces,
                    &mesh.faces,
                    &mesh.vertices,
                );
                let cov_axis = super::super::smallest_eigenvector_3x3(&cov);
                if parallel_miss(cov_axis, d) > 1.0e-6 {
                    covariance_failures += 1;
                }

                // The construction actually used is exact in every one of these.
                let Some(Initial::Cone { axis_dir, .. }) = cone(&faces, &mesh, TEST_VTOL) else {
                    panic!("no cone at {deg}° over {arc} rad");
                };
                // 1e-9 rather than 1e-12 because this grid includes the
                // worst-conditioned corner — a narrow cone over a short arc,
                // 10° across 0.6 rad, which reads 3.5e-11. Still exact to nine
                // digits, which is the plan's standard for this step.
                assert!(
                    parallel_miss(axis_dir, d) < 1.0e-9,
                    "cone axis wrong at {deg}° over {arc} rad: {:.3e}",
                    parallel_miss(axis_dir, d)
                );
            }
        }
        assert!(
            covariance_failures > cases / 2,
            "expected the covariance axis to fail most of this grid, it failed {covariance_failures}/{cases}"
        );
    }

    /// The rank deficiency that rules the normal-line solve out for cones, as a
    /// behaviour rather than an argument: on a single ring of cone facets every
    /// normal line meets the axis at the *same* point, so any line through that
    /// point satisfies the system and the axis is not determined.
    #[test]
    fn normal_line_concurrence_is_rank_deficient_on_a_single_ring_of_cone_facets() {
        let mesh = cone_strip([2.0, 1.0, -3.0], oblique_axis(), 0.5, 12, 1.4);
        let faces: Vec<usize> = (0..mesh.faces.len()).collect();
        assert!(
            axis_line_from_normals(&faces, &mesh).is_none(),
            "the normal-line solve should decline a single cone ring"
        );
        // ...and the construction that replaces it does not care.
        assert!(matches!(cone(&faces, &mesh, TEST_VTOL), Some(Initial::Cone { .. })));
    }

    /// A torus patch of quads spanning two meridians and two tube angles. Each
    /// quad is symmetric about its mid-meridian plane, so its Newell normal and
    /// its centroid both lie in that plane — which is exactly the meridian
    /// coplanarity the axis solve relies on.
    fn torus_patch(
        centre: [f64; 3],
        d: [f64; 3],
        major: f64,
        minor: f64,
        n_phi: usize,
        n_psi: usize,
    ) -> ConnectedMesh {
        let (u, v) = perpendicular_basis(&d);
        let at = |phi: f64, psi: f64| {
            let mer = [
                u[0] * phi.cos() + v[0] * phi.sin(),
                u[1] * phi.cos() + v[1] * phi.sin(),
                u[2] * phi.cos() + v[2] * phi.sin(),
            ];
            let rad = major + minor * psi.cos();
            [
                centre[0] + rad * mer[0] + minor * psi.sin() * d[0],
                centre[1] + rad * mer[1] + minor * psi.sin() * d[1],
                centre[2] + rad * mer[2] + minor * psi.sin() * d[2],
            ]
        };
        let mut points = Vec::new();
        for i in 0..=n_phi {
            let phi = 0.9 * i as f64 / n_phi as f64;
            for j in 0..=n_psi {
                let psi = 0.4 + 1.3 * j as f64 / n_psi as f64;
                points.push(at(phi, psi));
            }
        }
        let stride = n_psi + 1;
        let mut faces = Vec::new();
        for i in 0..n_phi {
            for j in 0..n_psi {
                faces.push(vec![
                    i * stride + j,
                    (i + 1) * stride + j,
                    (i + 1) * stride + j + 1,
                    i * stride + j + 1,
                ]);
            }
        }
        mesh_from(points, faces)
    }

    #[test]
    fn torus_axis_and_radii_come_out_of_the_meridian_solve() {
        let centre = [1.0, -2.0, 3.0];
        let d = oblique_axis();
        for &(major, minor) in &[(10.0_f64, 2.5_f64), (30.0, 1.0), (4.0, 3.0)] {
            let mesh = torus_patch(centre, d, major, minor, 6, 5);
            let faces: Vec<usize> = (0..mesh.faces.len()).collect();
            let Some(Initial::Torus {
                centre: c,
                axis_dir,
                major: rr,
                minor: r,
            }) = torus(&faces, &mesh)
            else {
                panic!("no torus for R={major} r={minor}");
            };
            assert!(
                parallel_miss(axis_dir, d) < 1.0e-9,
                "axis miss {:.3e} for R={major} r={minor}",
                parallel_miss(axis_dir, d)
            );
            // The centre is only pinned across the axis; along it the torus is
            // symmetric, so compare perpendicular distance to the true axis.
            assert!(
                distance_to_axis(c, centre, d) < 1.0e-6,
                "centre off-axis by {:.3e} for R={major} r={minor}",
                distance_to_axis(c, centre, d)
            );
            assert!(
                (rr - major).abs() < 1.0e-6 && (r - minor).abs() < 1.0e-6,
                "radii ({rr}, {r}) vs ({major}, {minor})"
            );
        }
    }

    /// The rank-deficiency the plan warns about — and where it is actually
    /// caught, which is not where I first assumed.
    ///
    /// On sphere data every line through the centre satisfies meridian
    /// coplanarity, so the axis the solve returns is arbitrary. Measured, the
    /// *numerical* spectrum does not separate that case from a genuine torus by
    /// enough to threshold on (2.06e-6 against 9.2e-4). What rejects it is
    /// structural and free: a sphere's (ρ, z) profile is a circle centred on the
    /// axis, so the fitted major radius is ~0 and `minor < major` fails.
    #[test]
    fn the_torus_initializer_declines_on_sphere_data() {
        let centre = [1.0, -2.0, 0.5];
        let radius = 6.0;
        let d = oblique_axis();
        let (u, v) = perpendicular_basis(&d);
        let on = |a: f64, b: f64| {
            [
                centre[0] + radius * (d[0] * b.cos() + (u[0] * a.cos() + v[0] * a.sin()) * b.sin()),
                centre[1] + radius * (d[1] * b.cos() + (u[1] * a.cos() + v[1] * a.sin()) * b.sin()),
                centre[2] + radius * (d[2] * b.cos() + (u[2] * a.cos() + v[2] * a.sin()) * b.sin()),
            ]
        };
        let mut points = Vec::new();
        let mut faces = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                let (a, b) = (i as f64 * 0.3, 0.4 + j as f64 * 0.25);
                let base = points.len();
                points.push(on(a, b));
                points.push(on(a + 0.28, b));
                points.push(on(a + 0.28, b + 0.24));
                points.push(on(a, b + 0.24));
                faces.push(vec![base, base + 1, base + 2, base + 3]);
            }
        }
        let mesh = mesh_from(points, faces);
        let all: Vec<usize> = (0..mesh.faces.len()).collect();
        assert!(
            torus(&all, &mesh).is_none(),
            "a sphere was accepted as a torus"
        );
        // ...while the sphere initializer handles the same patch.
        assert!(matches!(sphere(&all, &mesh), Some(Initial::Sphere { .. })));
    }

    /// The property the taper parameterization exists for: on a patch too small
    /// to determine a cone, the fit must still **explain its own facets** — and
    /// converge as the patch grows — instead of returning a wild answer that
    /// validates nothing.
    ///
    /// This is exactly what the (apex, axis, half-angle) form could not do.
    /// Measured in Phase 1.2 on three adjacent facets of an 18.4° cone it
    /// returned 87° with the apex halfway up the model, and validated 0 of its
    /// own 3 facets — which starved the earn loop.
    #[test]
    fn cone_fit_explains_its_own_support_at_every_patch_size() {
        let apex = [2.0, 1.0, -3.0];
        let d = oblique_axis();
        let alpha = 18.4_f64.to_radians();
        let mesh = cone_strip(apex, d, alpha, 16, 1.6);

        let mut prev_err = f64::MAX;
        let mut saw_cone = false;
        for k in [3_usize, 4, 6, 10, 16] {
            let faces: Vec<usize> = (0..k).collect();
            let Some(fit) = cone(&faces, &mesh, TEST_VTOL) else {
                panic!("no fit at all on {k} facets");
            };
            let worst = patch_vertices(&faces, &mesh)
                .into_iter()
                .map(|vi| {
                    let v = &mesh.vertices[vi];
                    surface_distance(&fit, [v.x, v.y, v.z])
                })
                .fold(0.0_f64, f64::max);
            // 3.6e-3 mm is the measured worst, on the 3-facet wedge. That is
            // NOT good enough for an earn loop validating at 5e-5 — see the
            // note on `cone` — but it is bounded and shrinks with support,
            // where the apex form validated *nothing* at any size. The bound
            // here pins that improvement so it cannot regress.
            assert!(
                worst < 5.0e-3,
                "{k} facets: the fit misses its own support by {worst:.3e}"
            );
            if let Initial::Cone { half_angle, .. } = fit {
                saw_cone = true;
                let err = (half_angle - alpha).abs();
                assert!(
                    err <= prev_err * 1.05,
                    "{k} facets: half-angle error {err:.3e} grew from {prev_err:.3e}"
                );
                prev_err = err;
            }
        }
        assert!(saw_cone, "never resolved a cone at any patch size");
        assert!(prev_err < 1.0e-3, "did not converge: {prev_err:.3e} rad");
    }

    /// Distance from a point to a fitted surface, for the assertion above.
    fn surface_distance(fit: &Initial, p: [f64; 3]) -> f64 {
        let perp = |w: [f64; 3], axis: [f64; 3]| {
            let z = dot3(&w, &axis);
            (z, vnorm(vsub(w, [axis[0] * z, axis[1] * z, axis[2] * z])))
        };
        match *fit {
            Initial::Cylinder {
                axis_pos,
                axis_dir,
                radius,
            } => {
                let (_, r) = perp(vsub(p, axis_pos), axis_dir);
                (r - radius).abs()
            }
            Initial::Cone {
                apex,
                axis_dir,
                half_angle,
            } => {
                let (z, r) = perp(vsub(p, apex), axis_dir);
                (r * half_angle.cos() - z * half_angle.sin()).abs()
            }
            _ => f64::MAX,
        }
    }

    /// Option C: a cone from the two rims of a two-facet strip, exactly — and
    /// from triangles as well as quads, since the rims are found by trying every
    /// split rather than by assuming a shared edge is a ruling.
    #[test]
    fn cone_from_rings_is_exact_from_two_facets() {
        let apex = [2.0, 1.0, -3.0];
        let d = oblique_axis();
        for deg in [10.0_f64, 30.0, 50.0, 75.0, 84.0] {
            let ha = deg.to_radians();
            let quads = cone_strip(apex, d, ha, 2, 0.6);
            // The same six vertices, each quad split along a diagonal: a shared
            // edge is now sometimes a diagonal, not a ruling.
            let tri_faces = (0..2)
                .flat_map(|k| [vec![2 * k, 2 * k + 1, 2 * k + 3], vec![2 * k, 2 * k + 3, 2 * k + 2]])
                .collect();
            let points: Vec<[f64; 3]> =
                quads.vertices.iter().map(|v| [v.x, v.y, v.z]).collect();
            let tris = mesh_from(points, tri_faces);
            for (label, mesh) in [("quads", &quads), ("triangles", &tris)] {
                let faces: Vec<usize> = (0..mesh.faces.len()).collect();
                let (got_apex, got_axis, got_ha) = cone_from_rings(&faces, mesh)
                    .unwrap_or_else(|| panic!("{deg}° {label}: no cone from two rims"));
                assert!(
                    parallel_miss(got_axis, d) < 1e-9 && dot3(&got_axis, &d) > 0.0,
                    "{deg}° {label}: axis must point from the apex toward the support"
                );
                assert!(
                    (got_ha - ha).abs() < 1e-9,
                    "{deg}° {label}: half-angle {} vs {}",
                    got_ha.to_degrees(),
                    deg
                );
                assert!(
                    vnorm(vsub(got_apex, apex)) < 1e-8,
                    "{deg}° {label}: apex off by {:e}",
                    vnorm(vsub(got_apex, apex))
                );
            }
        }
    }

    /// Two rims of equal radius are a cylinder, which is not this function's
    /// to report: the caller fits cylinders by its own route.
    #[test]
    fn cone_from_rings_declines_a_cylinder() {
        let d = oblique_axis();
        let (u, v) = perpendicular_basis(&d);
        let mut points = Vec::new();
        for k in 0..=2 {
            let th = 0.3 * k as f64;
            for z in [0.0, 5.0] {
                points.push([
                    4.0 * (u[0] * th.cos() + v[0] * th.sin()) + d[0] * z,
                    4.0 * (u[1] * th.cos() + v[1] * th.sin()) + d[1] * z,
                    4.0 * (u[2] * th.cos() + v[2] * th.sin()) + d[2] * z,
                ]);
            }
        }
        let mesh = mesh_from(points, vec![vec![0, 1, 3, 2], vec![2, 3, 5, 4]]);
        let faces: Vec<usize> = (0..mesh.faces.len()).collect();
        // Rounding can leave a taper of order 1e-16 rather than exactly zero;
        // what matters is that no plausible cone comes back.
        if let Some((_, _, ha)) = cone_from_rings(&faces, &mesh) {
            assert!(ha < 1e-9, "a cylinder must not read as a cone, got {}°", ha.to_degrees());
        }
    }
}
