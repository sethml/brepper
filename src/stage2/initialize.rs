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

/// Axis direction (**unoriented**) and half-angle from directions that lie
/// exactly on the surface.
///
/// Every point `p` of a cone satisfies `(p − A)/|p − A| · d = cos α` — one
/// linear equation per point in the four unknowns `(d, cos α)` — so the answer
/// is the null space of the stacked `[ûᵢ | 1]`. It uses **vertices, not facet
/// normals**, and that is the whole point:
///
/// * Vertices lie exactly on the surface at any coarseness, so this is exact.
/// * A facet normal does *not* give the half-angle exactly. A strip quad
///   spanning two rulings is planar — two rulings meet at the apex, so they
///   span a plane — but that chord plane leans away from the tangent plane, and
///   working the geometry through gives
///   `n·d = sin α · cos(Δθ/2) / √(1 − sin²α sin²(Δθ/2))`, an `O(Δθ²/8)` bias.
///   Measured on a 10-step, 1.1 rad strip at α = 10°: Δθ = 0.11 rad predicts
///   1.5e-3 relative, and the half-angle came out 9.985° against 10.000° — the
///   same 1.5e-3. Exactly the kind of "nearly right" an initializer must not
///   inherit.
///
/// The sign is deliberately not resolved here; [`cone`] orients the axis
/// geometrically from the apex, so it does not become a property of the mesh's
/// winding.
fn cone_direction_and_angle(
    apex: [f64; 3],
    faces: &[usize],
    mesh: &ConnectedMesh,
) -> Option<([f64; 3], f64)> {
    let mut m = nalgebra::Matrix4::<f64>::zeros();
    let mut rows = 0_usize;
    for vi in patch_vertices(faces, mesh) {
        let v = &mesh.vertices[vi];
        let mut u = [v.x - apex[0], v.y - apex[1], v.z - apex[2]];
        if normalize3(&mut u) < MIN_CROSS_NORM {
            continue; // a vertex at the apex says nothing about the axis
        }
        let row = nalgebra::Vector4::new(u[0], u[1], u[2], 1.0);
        m += row * row.transpose();
        rows += 1;
    }
    if rows < 3 {
        return None;
    }
    let eigen = nalgebra::SymmetricEigen::new(m);
    let mut order: Vec<usize> = (0..4).collect();
    order.sort_by(|&a, &b| eigen.eigenvalues[a].total_cmp(&eigen.eigenvalues[b]));
    let (l1, l3) = (eigen.eigenvalues[order[1]], eigen.eigenvalues[order[3]]);
    if !(l3 > 0.0 && l1 > l3 * MIN_NULLSPACE_SEPARATION) {
        return None; // the null space is not one-dimensional
    }

    let v = eigen.eigenvectors.column(order[0]);
    let mut d = [v[0], v[1], v[2]];
    let dnorm = normalize3(&mut d);
    if dnorm < MIN_CROSS_NORM {
        return None;
    }
    let cos_alpha = (v[3] / dnorm).abs();
    if !(cos_alpha > 0.0 && cos_alpha < 1.0) {
        return None; // a plane (α = 90°) or a degenerate solve
    }
    Some((d, cos_alpha.acos()))
}

/// A cone from its facets: the apex from where their planes concur, then the
/// direction and half-angle from the vertex directions about it.
///
/// **Neither of the two routes the plan and the module's own torus solve use
/// would work here**, and the reasons are worth recording:
///
/// * The *normal-covariance* axis the plan named is not the axis in general.
///   A cone's normals lie on a cone of the same half-angle about the axis, and
///   over a full revolution the covariance has eigenvalues `cos²α/2` twice and
///   `sin²α` once, so the axis is the smallest eigenvector only while
///   **α < 35.26°** — a 90° countersink is 45°, and the corpus's `cone_15x20`
///   is 36.87°. Over the partial arc an earned patch actually has, the
///   perpendicular spread is not even isotropic and it degrades further. The
///   test measures both regimes.
/// * The *normal-line concurrence* solve used for the torus is
///   **rank-deficient on a single ring of cone facets**, which is exactly what
///   a real cone tessellation is: a cone is ruled, so one row of quads spans its
///   whole height, every facet centroid sits at the same axial position, and
///   their normal lines then all meet the axis at the *same* point. Any line
///   through that point satisfies the system. (On a torus the meeting point
///   slides with the tube angle, which is why the same solve is well posed
///   there.)
///
/// What does work is exact and needs no arc at all: **a cone's tangent plane
/// along a ruling contains that ruling, and every ruling passes through the
/// apex, so every facet plane contains the apex** — including the chord planes
/// of strip quads, since two rulings meet at the apex. That is one linear
/// equation per facet, `nᵢ·A = nᵢ·pᵢ`, and three independent normals pin `A`
/// exactly, whatever the tessellation.
pub(crate) fn cone(faces: &[usize], mesh: &ConnectedMesh) -> Option<Initial> {
    // The apex is the common point of the facet planes.
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    let mut b = nalgebra::Vector3::<f64>::zeros();
    for &fi in faces {
        let Some(n) = mesh.faces[fi].normal else {
            continue;
        };
        let p = face_centroid(&mesh.faces[fi], &mesh.vertices);
        let row = nalgebra::Vector3::new(n[0], n[1], n[2]);
        m += row * row.transpose();
        b += row * dot3(&n, &p);
    }
    let eigen = nalgebra::SymmetricEigen::new(m);
    let mut evs: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
    evs.sort_by(f64::total_cmp);
    // The planes must span all three directions, or their intersection is a
    // line or a plane rather than a point.
    if !(evs[2] > 0.0 && evs[0] > evs[2] * MIN_NULLSPACE_SEPARATION) {
        return None;
    }
    let apex_v = m.lu().solve(&b)?;
    let apex = [apex_v[0], apex_v[1], apex_v[2]];

    let (mut axis_dir, half_angle) = cone_direction_and_angle(apex, faces, mesh)?;

    // Orient the axis from the apex toward the support — geometric, so it does
    // not inherit the mesh's normal convention.
    let mut toward = [0.0_f64; 3];
    for &fi in faces {
        let p = face_centroid(&mesh.faces[fi], &mesh.vertices);
        for k in 0..3 {
            toward[k] += p[k] - apex[k];
        }
    }
    if dot3(&toward, &axis_dir) < 0.0 {
        axis_dir = [-axis_dir[0], -axis_dir[1], -axis_dir[2]];
    }

    Some(Initial::Cone {
        apex,
        axis_dir,
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
            }) = cone(&faces, &mesh)
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
            assert!(
                (half_angle - alpha).to_degrees().abs() < 1.0e-9,
                "half-angle {:.9} vs {:.9} at {deg}°",
                half_angle.to_degrees(),
                deg
            );
            assert!(
                vnorm(vsub(a, apex)) < 1.0e-11,
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
                let Some(Initial::Cone { axis_dir, .. }) = cone(&faces, &mesh) else {
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
        assert!(matches!(cone(&faces, &mesh), Some(Initial::Cone { .. })));
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
}
