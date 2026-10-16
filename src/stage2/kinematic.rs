//! The kinematic rank ladder: one linear system whose **null-space dimension**
//! names the surface type, so that "the simplest surface that fits the support
//! so far" is a rank statement rather than a residual statement.
//!
//! See SURFACE_FITTING.md, "Nested fitting: a rank ladder instead of a type
//! ladder", for the survey this implements and for why the alternatives (general
//! quadric, quartic, Dupin cyclide) were declined.
//!
//! # The system
//!
//! Fit the instantaneous similarity velocity field
//!
//! ```text
//! v(p) = γ p + c × p + c̄
//! ```
//!
//! and require it tangent to the surface, `v(p)·n(p) = 0`. Each oriented facet
//! gives one homogeneous linear row in the seven unknowns `(γ, c, c̄)`:
//!
//! ```text
//! γ (p·n) + c·(p × n) + c̄·n = 0
//! ```
//!
//! The set of similarity fields tangent to a surface is a property of the
//! surface, and its dimension is the ladder:
//!
//! | surface | null dim | rank | facets to *confirm* it |
//! |---|---|---|---|
//! | plane | 4 | 3 | 4 |
//! | sphere | 3 | 4 | 5 |
//! | cylinder | 2 | 5 | 6 |
//! | cone | 2 | 5 | 6 |
//! | torus / surface of revolution | 1 | 6 | 7 |
//! | freeform | 0 | 7 | never — only "nothing simpler was confirmed" |
//!
//! Confirmation costs one facet more than the rank because `k` rows leave
//! `7 − k` eigenvalues zero whatever the surface is; see [`Spectrum::null_dim`].
//!
//! Under-determination therefore behaves correctly with no special case: `k`
//! rows can only pin rank `k`, so a patch that has not yet earned its type reads
//! as something simpler. Cylinder and cone share a dimension and are separated by
//! the null space's *content* — a cone admits a dilation (`γ ≠ 0`) from its apex
//! and a cylinder does not.
//!
//! # What this module is and is not
//!
//! It is a **classifier and an initializer**, never a validator. Every row is
//! built from a facet normal, which carries the `O(Δθ²/8)` chord-plane lean, so a
//! parameter read out of here is a seed for the vertex-based least squares and
//! not an answer — the lesson Phase 1.2 paid for. Nothing here decides
//! membership.
//!
//! There is deliberately **no threshold that converts a spectrum into a type**.
//! [`Spectrum::null_dim`] reads a dimension and [`Spectrum::candidates`] returns
//! every type that dimension admits, ordered; validation decides which one lives.
//! See `initialize::MIN_NULLSPACE_SEPARATION` for the measurement that says no
//! such threshold exists, and `BLOCK_PRESENT` here for a second one that had to
//! be learned the same way.
//!
//! # Measured status (2026-08-06, family K)
//!
//! Over 138 models and 1360 reference faces: plane 94.9%, sphere **100.0%** (at
//! the circumcenter), torus 81.5% by dimension, cone 66.3%, and **cylinder
//! 4.5%** — because a cylinder is tessellated as a single axial band of quads and
//! one band has a plane's rank at any circumferential extent
//! ([`tests::a_single_axial_band_cannot_pin_a_cylinder_at_any_width`]). So this
//! is a **complement** to `initialize::cylinder`, not a replacement for it; the
//! full result and what it changes are in SURFACE_FITTING.md.

// Phase K is instrumented before it is wired in: today the only callers are the
// unit tests and the family-K probe in `invariants.rs`.
#![allow(dead_code)]

use nalgebra::{SMatrix, SVector};

use super::{cross3, dot3, face_centroid, normalize3};
use crate::stage1::ConnectedMesh;

type Mat7 = SMatrix<f64, 7, 7>;
type Vec7 = SVector<f64, 7>;

/// Where on a facet the tangency constraint is evaluated.
///
/// The row wants a point *on the surface* paired with the *true* normal, and a
/// facet supplies neither exactly. Which point makes the two errors cancel
/// depends on the surface, which is why this is a parameter rather than a
/// constant:
///
/// - [`Eval::Centroid`] is exact for a ruling-aligned quad on a cylinder or
///   cone. The chord of a circular arc is parallel to the tangent at the arc's
///   *midpoint*, and the centroid sits at that midpoint angle, so the rotation
///   generator annihilates the row exactly. A vertex would be `O(Δθ)` wrong.
/// - [`Eval::Circumcenter`] is exact for a facet inscribed in a sphere, where
///   the normal line through the circumcenter passes through the centre — the
///   construction `initialize::sphere` already relies on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Eval {
    Centroid,
    Circumcenter,
}

/// The frame a system was normalized in, needed to map a recovered field back.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Frame {
    pub origin: [f64; 3],
    pub extent: f64,
}

/// The accumulated `7×7` Gram matrix `Σ rowᵢ rowᵢᵀ` and its frame.
///
/// Built incrementally — one rank-1 update per facet — which is the property
/// that lets a promotion check run at every face admission.
#[derive(Debug, Clone)]
pub(crate) struct System {
    m: Mat7,
    frame: Frame,
    rows: usize,
}

/// The circumcenter of a facet's vertices, in the facet's own plane.
///
/// Least squares over the facet's vertices, so it is the exact circumcenter for
/// a triangle and the best-fit circle centre for a fused quad (whose vertices
/// are concyclic when they lie on a sphere, since a plane meets a sphere in a
/// circle).
fn facet_circumcenter(mesh: &ConnectedMesh, fi: usize) -> Option<[f64; 3]> {
    let face = &mesh.faces[fi];
    let n = face.vertex_count as usize;
    if n < 3 {
        return None;
    }
    let pts: Vec<[f64; 3]> = (0..n)
        .map(|k| {
            let v = &mesh.vertices[face.vertex_indices[k]];
            [v.x, v.y, v.z]
        })
        .collect();
    let mut centre = [0.0_f64; 3];
    for p in &pts {
        for k in 0..3 {
            centre[k] += p[k] / n as f64;
        }
    }
    let mut normal = face.normal?;
    if normalize3(&mut normal) == 0.0 {
        return None;
    }
    // An orthonormal basis of the facet plane, from its longest edge so the
    // frame is well conditioned on a sliver.
    let mut e1 = [
        pts[1][0] - pts[0][0],
        pts[1][1] - pts[0][1],
        pts[1][2] - pts[0][2],
    ];
    if normalize3(&mut e1) == 0.0 {
        return None;
    }
    let e2 = cross3(&normal, &e1);

    // |p - c|² equal for all p: subtract pairwise to linearize in c.
    // 2(pᵢ - p₀)·c = |pᵢ|² - |p₀|², in plane coordinates about the mean.
    let local: Vec<(f64, f64)> = pts
        .iter()
        .map(|p| {
            let w = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
            (dot3(&w, &e1), dot3(&w, &e2))
        })
        .collect();
    let mut ata = nalgebra::Matrix2::<f64>::zeros();
    let mut atb = nalgebra::Vector2::<f64>::zeros();
    let (x0, y0) = local[0];
    for &(x, y) in &local[1..] {
        let row = nalgebra::Vector2::new(2.0 * (x - x0), 2.0 * (y - y0));
        let rhs = (x * x + y * y) - (x0 * x0 + y0 * y0);
        ata += row * row.transpose();
        atb += row * rhs;
    }
    let sol = ata.lu().solve(&atb)?;
    if !sol.iter().all(|v| v.is_finite()) {
        return None;
    }
    Some([
        centre[0] + e1[0] * sol[0] + e2[0] * sol[1],
        centre[1] + e1[1] * sol[0] + e2[1] * sol[1],
        centre[2] + e1[2] * sol[0] + e2[2] * sol[1],
    ])
}

/// The evaluation point and unit normal a facet contributes.
fn facet_sample(mesh: &ConnectedMesh, fi: usize, eval: Eval) -> Option<([f64; 3], [f64; 3])> {
    let mut n = mesh.faces[fi].normal?;
    if normalize3(&mut n) == 0.0 {
        return None;
    }
    let p = match eval {
        Eval::Centroid => face_centroid(&mesh.faces[fi], &mesh.vertices),
        Eval::Circumcenter => facet_circumcenter(mesh, fi)?,
    };
    Some((p, n))
}

impl System {
    /// Build the system for `faces`, normalized to their own centroid and
    /// extent.
    ///
    /// **Normalization is provably rank-preserving.** Under `p = o + L p'` the
    /// field `γp + c×p + c̄` becomes `L(γp' + c×p') + (γo + c×o + c̄)`, so in the
    /// primed coordinates the parameters are `(γ, c, (γo + c×o + c̄)/L)` — an
    /// invertible linear map of the 7-parameter space onto itself, leaving unit
    /// normals untouched. So this changes conditioning and never an answer,
    /// which matters because the `γ` and `c` blocks carry length units and the
    /// `c̄` block does not: on a patch far from the origin the null space would
    /// otherwise drown in the scale mismatch rather than in the data. This is
    /// the same defect, and the same fix, as `initialize::axis_line_from_normals`.
    pub(crate) fn build(faces: &[usize], mesh: &ConnectedMesh, eval: Eval) -> Option<System> {
        let samples: Vec<([f64; 3], [f64; 3])> = faces
            .iter()
            .filter_map(|&fi| facet_sample(mesh, fi, eval))
            .collect();
        if samples.is_empty() {
            return None;
        }
        let inv = 1.0 / samples.len() as f64;
        let mut origin = [0.0_f64; 3];
        for (p, _) in &samples {
            for k in 0..3 {
                origin[k] += p[k] * inv;
            }
        }
        let mut var = 0.0_f64;
        for (p, _) in &samples {
            let w = [p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]];
            var += dot3(&w, &w) * inv;
        }
        // A single facet, or a patch whose points coincide, has no extent of its
        // own to normalize by. Fall back to unit scale: the system is still
        // correct, just conditioned by the caller's units.
        let extent = if var > 0.0 { var.sqrt() } else { 1.0 };
        if !extent.is_finite() {
            return None;
        }
        let frame = Frame { origin, extent };

        let mut m = Mat7::zeros();
        let mut rows = 0usize;
        for (p, n) in &samples {
            let row = Self::row(*p, *n, frame);
            m += row * row.transpose();
            rows += 1;
        }
        Some(System { m, frame, rows })
    }

    /// One facet's row `[p·n, (p × n), n]` in the normalized frame.
    fn row(p: [f64; 3], n: [f64; 3], frame: Frame) -> Vec7 {
        let local = [
            (p[0] - frame.origin[0]) / frame.extent,
            (p[1] - frame.origin[1]) / frame.extent,
            (p[2] - frame.origin[2]) / frame.extent,
        ];
        let q = cross3(&local, &n);
        Vec7::from_row_slice(&[dot3(&local, &n), q[0], q[1], q[2], n[0], n[1], n[2]])
    }

    pub(crate) fn rows(&self) -> usize {
        self.rows
    }

    pub(crate) fn frame(&self) -> Frame {
        self.frame
    }

    /// Add one more facet — the incremental path, for a patch that grows.
    pub(crate) fn push(&mut self, fi: usize, mesh: &ConnectedMesh, eval: Eval) {
        if let Some((p, n)) = facet_sample(mesh, fi, eval) {
            let row = Self::row(p, n, self.frame);
            self.m += row * row.transpose();
            self.rows += 1;
        }
    }

    pub(crate) fn spectrum(&self) -> Spectrum {
        let eigen = nalgebra::SymmetricEigen::new(self.m);
        let mut order: Vec<usize> = (0..7).collect();
        order.sort_by(|&a, &b| eigen.eigenvalues[a].total_cmp(&eigen.eigenvalues[b]));
        let mut values = [0.0_f64; 7];
        let mut vectors = Mat7::zeros();
        for (slot, &src) in order.iter().enumerate() {
            values[slot] = eigen.eigenvalues[src].max(0.0);
            vectors.set_column(slot, &eigen.eigenvectors.column(src));
        }
        Spectrum {
            values,
            vectors,
            rows: self.rows,
            frame: self.frame,
        }
    }

    /// `‖M v‖ / λ_max` for a unit `v` — how well a *claimed* tangent field is
    /// annihilated by the measured rows, relative to the system's own scale.
    ///
    /// This is the family-K analogue of family F: fed the reference surface's own
    /// generators it asks whether the signal is there to be read, before asking
    /// whether the read succeeds.
    pub(crate) fn generator_residual(&self, v: Vec7) -> f64 {
        // `<= 0.0` rather than a negated `>`: these must also reject NaN, and
        // an explicit `is_finite` says so where `!(x > 0.0)` only implies it.
        let scale = self.m.trace();
        if !scale.is_finite() || scale <= 0.0 {
            return f64::NAN;
        }
        let mut u = v;
        let norm = u.norm();
        if !norm.is_finite() || norm <= 0.0 {
            return f64::NAN;
        }
        u /= norm;
        (self.m * u).norm() / scale
    }
}

/// The eigen-decomposition of a [`System`], ascending.
#[derive(Debug, Clone)]
pub(crate) struct Spectrum {
    /// Eigenvalues, ascending. The null space is at the low end.
    pub values: [f64; 7],
    /// Matching eigenvectors, column `i` for `values[i]`.
    pub vectors: Mat7,
    pub rows: usize,
    pub frame: Frame,
}

/// What a null space of a given dimension and content names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum KinematicType {
    /// Fewer independent rows than the ladder needs — nothing is pinned yet.
    Underdetermined,
    Plane,
    Sphere,
    Cylinder,
    Cone,
    /// A rotation survives: torus, or any other surface of revolution.
    Revolution,
    /// A translation survives and no rotation does.
    Extrusion,
    /// One field survives carrying a dilation — a spiral or conical helix.
    Spiral,
    Freeform,
}

impl KinematicType {
    pub(crate) fn name(self) -> &'static str {
        match self {
            KinematicType::Underdetermined => "under",
            KinematicType::Plane => "plane",
            KinematicType::Sphere => "sphere",
            KinematicType::Cylinder => "cylinder",
            KinematicType::Cone => "cone",
            KinematicType::Revolution => "revolution",
            KinematicType::Extrusion => "extrusion",
            KinematicType::Spiral => "spiral",
            KinematicType::Freeform => "freeform",
        }
    }
}

/// Below this share of a null vector's norm, a block is treated as absent when
/// **ordering** candidates.
///
/// It was originally documented as a validity floor on the grounds that "a
/// cylinder's null space has `γ` identically zero by construction, a cone's has
/// `γ` of order one — so any cut well below 1 gives the same answer". **That is
/// measured false**, and by a wide margin: on the corpus the rank reads a true
/// torus's dimension correctly 79–81% of the time while this floor labels 70% of
/// those readings `Spiral`, because a real facet normal carries `O(Δθ²/8)` error
/// and the `γ` block picks it up at ~1e-3 — a thousand times over the cut. The
/// clean structural separation exists on exact data only, which is exactly where
/// the unit tests live and is why they did not catch it.
///
/// So it no longer *decides* anything. [`Spectrum::candidates`] returns every
/// type a dimension admits, ordered by this quantity, and validation picks —
/// which is the same rule the rest of stage 2 follows and the same rule
/// `initialize::MIN_NULLSPACE_SEPARATION` records from the other direction.
const BLOCK_PRESENT: f64 = 1.0e-6;

impl Spectrum {
    /// The null dimension **that row starvation alone does not already explain**,
    /// with the ratio gap that produced it.
    ///
    /// The subtlety this exists for, found by measurement rather than by
    /// derivation: with `k` rows the Gram matrix has rank at most `k`, so
    /// `7 − k` of its eigenvalues are exactly zero **for a reason that has
    /// nothing to do with the surface**. A naive widest-gap read therefore just
    /// reports `7 − rows` and looks impressively confident doing it — measured
    /// on `simple_sphere` and `simple_cone`, it walked 4, 3, 2, 1, 0 as the
    /// support went 3, 4, 5, 6, 8, matching the true type exactly once each, by
    /// coincidence, at the size where `7 − k` happened to equal it.
    ///
    /// So only positions **above** the trivial floor are considered. The
    /// consequence is a support cost of one row more than the type's rank:
    ///
    /// | type | rank | rows to *confirm* it |
    /// |---|---|---|
    /// | plane | 3 | 4 |
    /// | sphere | 4 | 5 |
    /// | cylinder, cone | 5 | 6 |
    /// | torus / revolution | 6 | 7 |
    ///
    /// with freeform never confirmable — only ever "nothing simpler was
    /// confirmed", which is the correct epistemic status for it.
    ///
    /// **The returned ratio is the whole point and must not be thresholded
    /// here.** On exact synthetic data it is ~1e13; on a real mesh it is set by
    /// facet-normal error and can be as low as ~10. Report it, order candidates
    /// by it, and let validation decide — see
    /// `initialize::MIN_NULLSPACE_SEPARATION` for the measurement that killed
    /// the idea of a cut.
    pub(crate) fn null_dim(&self) -> (usize, f64) {
        let top = self.values[6];
        if !top.is_finite() || top <= 0.0 {
            return (7, f64::INFINITY);
        }
        // Floor-clamp so an exact zero gives a large finite ratio rather than
        // an infinity that cannot be compared.
        let floor = top * 1.0e-18;
        let mut best = (0usize, 0.0_f64);
        for k in 1..7 {
            let lo = self.values[k - 1].max(floor);
            let hi = self.values[k].max(floor);
            let ratio = hi / lo;
            if ratio > best.1 {
                best = (k, ratio);
            }
        }
        // The gap is read over the whole spectrum and *then* compared against
        // the trivial floor, rather than searched for above it: restricting the
        // search would force a wrong answer whenever the true dimension sits at
        // or below the floor, instead of admitting that nothing is pinned yet.
        best
    }

    /// Every type a null space of dimension `dim` admits, **best-ordered first**.
    ///
    /// The dimension is the measurement; the content only *orders* the types that
    /// dimension permits. Returning a list rather than a winner is not hedging —
    /// it is what the corpus forced (see [`BLOCK_PRESENT`]): a real facet normal
    /// puts ~1e-3 into a block that is exactly zero in theory, so a content test
    /// sharp enough to name one type gets 70% of true tori wrong while the
    /// dimension underneath is right. Downstream, the caller fits each candidate
    /// and lets vertex validation choose, exactly as stage 2.3's ladder does.
    pub(crate) fn candidates(&self, dim: usize) -> Vec<KinematicType> {
        // `rows` rows leave `7 - rows` eigenvalues zero whatever the surface is,
        // so a dimension at or below that floor is row starvation and names
        // nothing.
        if dim + self.rows <= 7 {
            return vec![KinematicType::Underdetermined];
        }
        // The largest dilation component available anywhere in the null space.
        // Test the whole subspace, not one basis vector: the eigenvectors are an
        // arbitrary basis of it.
        let gamma = (0..dim.min(7))
            .map(|k| self.vectors[(0, k)].abs())
            .fold(0.0_f64, f64::max);
        match dim {
            0 => vec![KinematicType::Freeform],
            1 => {
                let v = self.vectors.column(0);
                let rot = (v[1] * v[1] + v[2] * v[2] + v[3] * v[3]).sqrt();
                if rot <= BLOCK_PRESENT {
                    vec![KinematicType::Extrusion, KinematicType::Revolution]
                } else if gamma > rot {
                    vec![KinematicType::Spiral, KinematicType::Revolution]
                } else {
                    vec![KinematicType::Revolution, KinematicType::Spiral]
                }
            }
            // A cone's 2-space contains the dilation from its apex; a cylinder's
            // is rotation + translation along the axis, with γ zero in theory.
            2 if gamma > BLOCK_PRESENT => vec![KinematicType::Cone, KinematicType::Cylinder],
            2 => vec![KinematicType::Cylinder, KinematicType::Cone],
            3 => vec![KinematicType::Sphere],
            4 => vec![KinematicType::Plane],
            _ => vec![KinematicType::Underdetermined],
        }
    }

    /// The single best-ordered candidate — for reporting and for the tests.
    /// Prefer [`candidates`](Self::candidates) anywhere a decision follows.
    pub(crate) fn classify(&self, dim: usize) -> KinematicType {
        self.candidates(dim)[0]
    }

    /// The type the spectrum reads on its own.
    pub(crate) fn read(&self) -> (KinematicType, usize, f64) {
        let (dim, ratio) = self.null_dim();
        (self.classify(dim), dim, ratio)
    }
}

// ---------------------------------------------------------------------------
// The reference's own generators, for family K's satisfiability question
// ---------------------------------------------------------------------------

/// A tangent similarity field in the *world* frame, mapped into a system's
/// normalized frame so it can be tested against the measured rows.
///
/// The map is the one derived in [`System::build`]: `(γ, c, c̄)` in world
/// coordinates becomes `(γ, c, (γ o + c × o + c̄) / L)` in the frame.
pub(crate) fn generator_in_frame(gamma: f64, c: [f64; 3], cbar: [f64; 3], frame: Frame) -> Vec7 {
    let co = cross3(&c, &frame.origin);
    let t = [
        (gamma * frame.origin[0] + co[0] + cbar[0]) / frame.extent,
        (gamma * frame.origin[1] + co[1] + cbar[1]) / frame.extent,
        (gamma * frame.origin[2] + co[2] + cbar[2]) / frame.extent,
    ];
    Vec7::from_row_slice(&[gamma, c[0], c[1], c[2], t[0], t[1], t[2]])
}

/// Rotation about the line through `pos` with direction `dir`: `v(p) = d × (p − a)`.
pub(crate) fn rotation_about(pos: [f64; 3], dir: [f64; 3], frame: Frame) -> Vec7 {
    let mut d = dir;
    normalize3(&mut d);
    // d × (p − a) = d × p + (a × d), so c = d and c̄ = a × d.
    let cbar = cross3(&pos, &d);
    generator_in_frame(0.0, d, cbar, frame)
}

/// Translation along `dir`.
pub(crate) fn translation_along(dir: [f64; 3], frame: Frame) -> Vec7 {
    let mut d = dir;
    normalize3(&mut d);
    generator_in_frame(0.0, [0.0; 3], d, frame)
}

/// Dilation about `centre`: `v(p) = p − a`.
pub(crate) fn dilation_about(centre: [f64; 3], frame: Frame) -> Vec7 {
    generator_in_frame(1.0, [0.0; 3], [-centre[0], -centre[1], -centre[2]], frame)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stage1::{MeshFace, MeshVertex};

    /// An oblique unit axis: an axis-aligned one would let a sign or an index
    /// slip through every assertion here.
    fn oblique_axis() -> [f64; 3] {
        let mut d = [0.3, -0.5, 0.81];
        normalize3(&mut d);
        d
    }

    /// An origin well away from the coordinate origin, so the normalization is
    /// actually exercised — an untranslated patch would pass even with the
    /// scaling bug the frame exists to prevent.
    fn far_origin() -> [f64; 3] {
        [120.0, -45.0, 260.0]
    }

    fn mesh_from(points: Vec<[f64; 3]>, faces: Vec<Vec<usize>>) -> ConnectedMesh {
        let vertices: Vec<MeshVertex> = points
            .iter()
            .map(|p| MeshVertex::from_xyz(p[0], p[1], p[2]))
            .collect();
        let faces = faces
            .into_iter()
            .map(|idx| {
                // The Newell normal, exactly as stage 1 computes it.
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

    fn perp_basis(d: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        super::super::perpendicular_basis(&d)
    }

    /// A ruling-aligned quad strip: `rings` circles of `divs` points about the
    /// axis, at radius `radius(z)`.
    fn revolution_strip(
        origin: [f64; 3],
        axis: [f64; 3],
        zs: &[f64],
        divs: usize,
        radius: impl Fn(f64) -> f64,
    ) -> ConnectedMesh {
        let (u, v) = perp_basis(axis);
        let mut points = Vec::new();
        for &z in zs {
            let r = radius(z);
            for k in 0..divs {
                let t = std::f64::consts::TAU * k as f64 / divs as f64;
                let (c, s) = (t.cos(), t.sin());
                points.push([
                    origin[0] + axis[0] * z + r * (c * u[0] + s * v[0]),
                    origin[1] + axis[1] * z + r * (c * u[1] + s * v[1]),
                    origin[2] + axis[2] * z + r * (c * u[2] + s * v[2]),
                ]);
            }
        }
        let mut faces = Vec::new();
        for ring in 0..zs.len() - 1 {
            for k in 0..divs - 1 {
                let a = ring * divs + k;
                let b = ring * divs + k + 1;
                faces.push(vec![a, b, b + divs, a + divs]);
            }
        }
        mesh_from(points, faces)
    }

    fn all(mesh: &ConnectedMesh) -> Vec<usize> {
        (0..mesh.faces.len()).collect()
    }

    fn read(mesh: &ConnectedMesh) -> (KinematicType, usize, f64) {
        System::build(&all(mesh), mesh, Eval::Centroid)
            .expect("system")
            .spectrum()
            .read()
    }

    #[test]
    fn a_plane_reads_null_dimension_four() {
        // Six triangles of a planar patch, deliberately irregular so no
        // symmetry can hide an error.
        let d = oblique_axis();
        let (u, v) = perp_basis(d);
        let o = far_origin();
        let pt = |s: f64, t: f64| {
            [
                o[0] + s * u[0] + t * v[0],
                o[1] + s * u[1] + t * v[1],
                o[2] + s * u[2] + t * v[2],
            ]
        };
        let points = vec![
            pt(0.0, 0.0),
            pt(3.0, 0.4),
            pt(1.1, 2.7),
            pt(-2.2, 1.9),
            pt(-1.7, -2.4),
            pt(2.3, -3.1),
        ];
        let mesh = mesh_from(
            points,
            vec![
                vec![0, 1, 2],
                vec![0, 2, 3],
                vec![0, 3, 4],
                vec![0, 4, 5],
                vec![0, 5, 1],
            ],
        );
        let (kind, dim, _) = read(&mesh);
        assert_eq!(dim, 4, "a plane's tangent similarity fields are 4-dimensional");
        assert_eq!(kind, KinematicType::Plane);
    }

    #[test]
    fn a_sphere_reads_null_dimension_three_at_the_circumcenter() {
        // A sphere patch as a quad grid in spherical coordinates. Its facets are
        // NOT equilateral, so centroid and circumcenter genuinely differ — which
        // is the property under test.
        let centre = far_origin();
        let radius = 7.0;
        let mut points = Vec::new();
        let (nu, nv) = (5usize, 6usize);
        for i in 0..nu {
            let phi = 0.6 + 0.9 * i as f64 / (nu - 1) as f64;
            for j in 0..nv {
                let th = 0.3 + 1.4 * j as f64 / (nv - 1) as f64;
                points.push([
                    centre[0] + radius * phi.sin() * th.cos(),
                    centre[1] + radius * phi.sin() * th.sin(),
                    centre[2] + radius * phi.cos(),
                ]);
            }
        }
        let mut faces = Vec::new();
        for i in 0..nu - 1 {
            for j in 0..nv - 1 {
                let a = i * nv + j;
                // Triangles: a quad on a sphere is not planar, so a quad facet
                // here would have a meaningless normal.
                faces.push(vec![a, a + 1, a + nv]);
                faces.push(vec![a + 1, a + nv + 1, a + nv]);
            }
        }
        let mesh = mesh_from(points, faces);

        let sys = System::build(&all(&mesh), &mesh, Eval::Circumcenter).expect("system");
        let (kind, dim, _) = sys.spectrum().read();
        assert_eq!(dim, 3, "a sphere's tangent fields are its three rotations");
        assert_eq!(kind, KinematicType::Sphere);

        // And the circumcenter is exact where the centroid is not: the true
        // generators are annihilated to machine precision by one and not the
        // other. This is the module's evaluation-point claim, in a test.
        let g = rotation_about(centre, oblique_axis(), sys.frame());
        let exact = sys.generator_residual(g);
        let centroid_sys = System::build(&all(&mesh), &mesh, Eval::Centroid).expect("system");
        let sloppy =
            centroid_sys.generator_residual(rotation_about(centre, oblique_axis(), centroid_sys.frame()));
        assert!(
            exact < 1.0e-14,
            "circumcenter rows annihilate a true rotation exactly, got {exact:e}"
        );
        assert!(
            sloppy > exact * 100.0,
            "the centroid should be measurably worse on a sphere: {sloppy:e} vs {exact:e}"
        );
    }

    #[test]
    fn a_cylinder_reads_null_dimension_two_with_no_dilation() {
        let axis = oblique_axis();
        let mesh = revolution_strip(far_origin(), axis, &[-3.0, 0.0, 2.5, 6.0], 16, |_| 4.0);
        let (kind, dim, _) = read(&mesh);
        assert_eq!(dim, 2, "rotation about the axis plus translation along it");
        assert_eq!(
            kind,
            KinematicType::Cylinder,
            "a cylinder's null space has no dilation in it"
        );
    }

    #[test]
    fn a_cone_reads_null_dimension_two_with_a_dilation() {
        let axis = oblique_axis();
        let apex = far_origin();
        // r = z tan α, so the apex is at z = 0 and the patch is a real cone.
        let tan_a = (18.0_f64).to_radians().tan();
        let mesh = revolution_strip(apex, axis, &[4.0, 7.0, 11.0, 16.0], 16, |z| z * tan_a);
        let (kind, dim, _) = read(&mesh);
        assert_eq!(dim, 2, "rotation about the axis plus dilation from the apex");
        assert_eq!(
            kind,
            KinematicType::Cone,
            "the dilation from the apex is what separates a cone from a cylinder"
        );
    }

    #[test]
    fn a_torus_reads_null_dimension_one() {
        let axis = oblique_axis();
        let centre = far_origin();
        let (u, v) = perp_basis(axis);
        let (major, minor) = (9.0, 2.0);
        let (nt, np) = (12usize, 9usize);
        let mut points = Vec::new();
        for i in 0..nt {
            let t = std::f64::consts::TAU * i as f64 / nt as f64;
            for j in 0..np {
                // A partial tube sweep, which is what a fillet actually is.
                let s = -0.9 + 1.8 * j as f64 / (np - 1) as f64;
                let rho = major + minor * s.cos();
                let z = minor * s.sin();
                points.push([
                    centre[0] + rho * (t.cos() * u[0] + t.sin() * v[0]) + z * axis[0],
                    centre[1] + rho * (t.cos() * u[1] + t.sin() * v[1]) + z * axis[1],
                    centre[2] + rho * (t.cos() * u[2] + t.sin() * v[2]) + z * axis[2],
                ]);
            }
        }
        let mut faces = Vec::new();
        for i in 0..nt - 1 {
            for j in 0..np - 1 {
                let a = i * np + j;
                faces.push(vec![a, a + np, a + np + 1, a + 1]);
            }
        }
        let mesh = mesh_from(points, faces);
        let (kind, dim, _) = read(&mesh);
        assert_eq!(dim, 1, "only rotation about the torus axis survives");
        assert_eq!(kind, KinematicType::Revolution);
    }

    #[test]
    fn the_reference_generators_are_annihilated_on_a_ruled_strip() {
        // The claim that makes the centroid the right evaluation point on a
        // developable: the chord of an arc is parallel to the tangent at the
        // arc's MIDPOINT, and a strip quad's centroid sits at that angle.
        let axis = oblique_axis();
        let pos = far_origin();
        let mesh = revolution_strip(pos, axis, &[-3.0, 0.0, 2.5, 6.0], 16, |_| 4.0);
        let sys = System::build(&all(&mesh), &mesh, Eval::Centroid).expect("system");
        let rot = sys.generator_residual(rotation_about(pos, axis, sys.frame()));
        let tra = sys.generator_residual(translation_along(axis, sys.frame()));
        assert!(rot < 1.0e-14, "rotation about the true axis, got {rot:e}");
        assert!(tra < 1.0e-14, "translation along the true axis, got {tra:e}");

        // And a cone's apex dilation, on the same construction.
        let apex = far_origin();
        let tan_a = (18.0_f64).to_radians().tan();
        let cone = revolution_strip(apex, axis, &[4.0, 7.0, 11.0, 16.0], 16, |z| z * tan_a);
        let csys = System::build(&all(&cone), &cone, Eval::Centroid).expect("system");
        let dil = csys.generator_residual(dilation_about(apex, csys.frame()));
        assert!(
            dil < 1.0e-14,
            "every facet plane of a cone contains the apex, got {dil:e}"
        );
    }

    #[test]
    fn support_below_the_rank_the_type_needs_reads_simpler() {
        // The property the whole ladder rests on: too few rows cannot pin a
        // type, and what comes back is a SIMPLER surface rather than a wrong
        // one. Two facets of a cylinder read as a plane's dimension or more.
        let axis = oblique_axis();
        let mesh = revolution_strip(far_origin(), axis, &[-3.0, 0.0, 2.5, 6.0], 16, |_| 4.0);
        let two = vec![0usize, 1usize];
        let sys = System::build(&two, &mesh, Eval::Centroid).expect("system");
        let (_, dim, _) = sys.spectrum().read();
        assert!(
            dim >= 5,
            "two rows can pin at most rank 2, so the null space is 5-D or more; got {dim}"
        );
        assert_eq!(
            sys.spectrum().classify(dim),
            KinematicType::Underdetermined,
            "and that must be reported as under-determination, not as a type"
        );

        // A cylinder has rank 5, so five rows reach its null dimension — and
        // cannot *confirm* it, because row starvation alone would produce the
        // same two zeros. Confirmation costs one row more than the rank. Faces
        // are laid out `ring * (divs - 1) + k`, so 0,1,2 are one band and
        // 15,16,17 the next, which is also the axial spread the band test below
        // shows is mandatory.
        let five = vec![0usize, 1, 2, 15, 16];
        let sys5 = System::build(&five, &mesh, Eval::Centroid).expect("system");
        assert_eq!(
            sys5.spectrum().read().0,
            KinematicType::Underdetermined,
            "five rows cannot distinguish a cylinder's null space from row starvation"
        );

        let six = vec![0usize, 1, 2, 15, 16, 17];
        let sys6 = System::build(&six, &mesh, Eval::Centroid).expect("system");
        let (kind6, dim6, _) = sys6.spectrum().read();
        assert_eq!(dim6, 2, "six independent rows confirm the cylinder");
        assert_eq!(kind6, KinematicType::Cylinder);
    }

    #[test]
    fn a_single_axial_band_cannot_pin_a_cylinder_at_any_width() {
        // The ladder's sharpest limitation, and it is structural rather than
        // numerical. In a band of ruling-aligned quads every centroid sits at
        // the same axial position `z_c` and the same radial distance `r_c`, so
        //
        //     p·n = r_c                      (the same constant for every facet)
        //     p × n = z_c (d × n)            (a fixed multiple of the normal)
        //
        // and every row is an affine function of `n` alone. The row space is
        // spanned by {constant, n_x, n_y} — **rank 3, whatever the arc** — so a
        // band reads as a plane no matter how far around the cylinder it goes.
        // The band carries no axial information at all: it is "a short arc
        // cannot determine a surface of revolution" applied to the axial
        // direction instead of the circumferential one.
        //
        // This is why the ladder needs a two-dimensional patch and not a chain,
        // the same conclusion the measured seed rules reached for the torus.
        let axis = oblique_axis();
        let mesh = revolution_strip(far_origin(), axis, &[-3.0, 0.0, 2.5, 6.0], 16, |_| 4.0);
        for width in [5usize, 12] {
            let band: Vec<usize> = (0..width).collect();
            let sp = System::build(&band, &mesh, Eval::Centroid)
                .expect("system")
                .spectrum();
            let (kind, dim, _) = sp.read();
            assert_eq!(
                (kind, dim),
                (KinematicType::Plane, 4),
                "a {width}-facet single band still reads plane"
            );
        }
    }

    #[test]
    fn normalization_does_not_change_the_reading() {
        // The rank-preservation argument in `System::build`, as a test: the same
        // cylinder a metre from the origin reads identically.
        let axis = oblique_axis();
        let near = revolution_strip([0.0, 0.0, 0.0], axis, &[-3.0, 0.0, 2.5, 6.0], 16, |_| 4.0);
        let far = revolution_strip(
            [1000.0, -2000.0, 4000.0],
            axis,
            &[-3.0, 0.0, 2.5, 6.0],
            16,
            |_| 4.0,
        );
        let (kn, dn, _) = read(&near);
        let (kf, df, _) = read(&far);
        assert_eq!((kn, dn), (kf, df), "a translated patch must read the same");
    }
}
