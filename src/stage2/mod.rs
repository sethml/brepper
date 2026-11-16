//! Stage 2: Segmentation & Classification (partition-then-classify)
//!
//! Sub-stage numbering (DEVELOPMENT_PLAN):
//!
//! 2.1: Sharp-edge classification (fit-free) + osculating-jet prior
//!      (`edges.rs`).
//! 2.2: Region extraction — connected components over interior edges
//!      (`regions.rs`).
//! 2.3: Per-region surface discovery (`classify.rs` driver) — whole-region
//!      DOF ladder (2.3(a)) then scoped fit competition (2.3(b), `compete.rs`)
//!      with bent-tube/torus promotion (2.3(c)) for the tangent chains the
//!      ladder cannot explain, and B-spline fitting of the leftover components
//!      (2.3(e), `bspline.rs`). Whatever remains degrades to per-facet planes.
//! 2.4: Resolution & repair — boundary snap / relabel / merge / re-split +
//!      invariants (`repair.rs`, with the junction geometry in `locus.rs`).
//! 2.5: Selected-surface output.
//!
//! The analytic fitters (`planar.rs`, `conical.rs`, `spherical.rs`,
//! `toroidal.rs`) share the `surface::SurfaceType` trait, driven generically
//! by the ladder and competition. The old BFS hypothesis-generation pipeline
//! that predated this rewrite was removed (Rewrite Roadmap step 7).

pub mod edges;
pub(crate) mod initialize;
pub mod invariants;
pub(crate) mod kinematic;
pub mod regions;

mod surface;
pub(crate) mod bspline;
mod classify;
mod compete;
mod conical;
mod cylindrical;
mod locus;
mod planar;
mod repair;
mod spherical;
mod toroidal;

pub use self::bspline::BSplineHypothesis;
pub use self::edges::{EdgeClass, EdgeClassification, JetField};
pub use self::regions::Region;

use crate::config::Config;
use crate::stage1::{ConnectedMesh, MeshFace, MeshVertex};
use crate::viz::{self, VizAction};
use std::time::Instant;
use std::error::Error;
use std::fmt::{Display, Formatter};

// Re-export used by conical.rs via `super::bounding_box_diagonal`.
use self::spherical::bounding_box_diagonal;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Minimum number of mesh faces required to accept a cylindrical hypothesis.
/// Any real CAD tessellation of a cylinder produces at least 3 facets around
/// the circumference; this rejects spurious 2-face fits (e.g. adjacent torus
/// facets that locally approximate a cylinder).
pub(crate) const MIN_CYLINDER_FACES: usize = 3;

/// Minimum number of mesh faces required to accept a spherical hypothesis.
/// Sphere fitting has 4 degrees of freedom (cx, cy, cz, r), so at least 4
/// non-degenerate faces are needed for a well-determined fit. This also rejects
/// spurious fits from small patches that are consistent with many surface types.
pub(crate) const MIN_SPHERE_FACES: usize = 4;

/// Minimum number of mesh faces required to accept a toroidal hypothesis.
/// A torus has 7 degrees of freedom (center xyz, axis direction 2, R, r),
/// so at least 7 faces are needed for a well-determined fit.
pub(crate) const MIN_TORUS_FACES: usize = 7;

/// Maximum sphere radius as a multiple of the mesh bounding-box diagonal.
/// With solid-angle coverage validation and surface-tolerance validation during
/// BFS, this no longer needs to be tight — those checks prevent pathological
/// growth. A large value serves as a numerical guardrail against degenerate fits.
pub(crate) const MAX_SPHERE_RADIUS_FACTOR: f64 = 1000.0;

// ---------------------------------------------------------------------------
// Hypothesis data structures
// ---------------------------------------------------------------------------

/// A planar surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct PlanarHypothesis {
    /// Unit normal vector pointing outward from the shell/solid.
    pub normal: [f64; 3],
    /// Signed distance from origin to plane along the normal.
    pub distance: f64,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this plane.
    pub vertices: Vec<usize>,
    /// Maximum (positive) distance from any vertex to the plane.
    pub error_max: f64,
    /// Minimum (most negative) distance from any vertex to the plane.
    pub error_min: f64,
    /// Sum of absolute vertex-to-plane distances.
    pub error_abs_sum: f64,
}

/// A cylindrical or conical surface hypothesis fitted to a set of mesh faces.
/// Cylinders are represented as cones with half_angle == 0.
#[derive(Debug, Clone)]
pub struct ConicalHypothesis {
    /// The apex of the cone, or a point on the axis for cylinders.
    pub apex: [f64; 3],
    /// Unit direction vector along the cone/cylinder axis.
    pub axis_direction: [f64; 3],
    /// Half-angle of the cone in radians (0 for cylinders).
    pub half_angle: f64,
    /// Cylinder radius (when half_angle == 0). 0 for true cones.
    pub radius_at_apex: f64,
    /// Whether the surface normal points away from the axis (convex=true) or toward it (concave=false).
    pub convex: bool,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this surface.
    pub vertices: Vec<usize>,
    /// Maximum absolute distance from any vertex to the surface.
    pub error_max: f64,
    /// Maximum absolute distance from any face centroid to the surface.
    pub centroid_error_max: f64,
    /// Sum of absolute vertex-to-surface distances.
    pub error_abs_sum: f64,
}

impl ConicalHypothesis {
    /// Returns true if this hypothesis represents a cylinder (half_angle ≈ 0).
    pub fn is_cylinder(&self) -> bool {
        self.half_angle < 1e-4
    }
}

/// A spherical surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct SphericalHypothesis {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere (always positive).
    pub radius: f64,
    /// Whether the surface normal points away from center (convex) or toward it (concave).
    pub convex: bool,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this sphere.
    pub vertices: Vec<usize>,
    /// Maximum absolute distance from any vertex to the sphere surface.
    pub error_max: f64,
    /// Maximum absolute distance from any face centroid to the sphere surface.
    pub centroid_error_max: f64,
    /// Sum of absolute vertex-to-surface distances.
    pub error_abs_sum: f64,
}

/// A toroidal surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct ToroidalHypothesis {
    /// Center of the torus (on the axis, in the plane of the major circle).
    pub center: [f64; 3],
    /// Unit direction vector along the torus rotation axis.
    pub axis_direction: [f64; 3],
    /// Distance from center to the tube centerline (major radius R).
    pub major_radius: f64,
    /// Radius of the tube cross-section (minor radius r).
    pub minor_radius: f64,
    /// Whether face normals point away from the tube center (convex=true, outer surface)
    /// or toward the tube center (concave=false, inner surface).
    pub convex: bool,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this torus.
    pub vertices: Vec<usize>,
    /// Maximum absolute distance from any vertex to the torus surface.
    pub error_max: f64,
    /// Maximum absolute distance from any face centroid to the torus surface.
    pub centroid_error_max: f64,
    /// Sum of absolute vertex-to-surface distances.
    pub error_abs_sum: f64,
}

/// Identifies a selected surface hypothesis by type and index.
/// `Conical` indexes into `conical_hypotheses` and covers both cylinders (half_angle≈0) and true cones.
#[derive(Debug, Clone, Copy)]
pub enum SelectedSurface {
    Planar(usize),
    /// Index into conical_hypotheses. Covers both cylinders (half_angle≈0) and true cones.
    Conical(usize),
    Spherical(usize),
    Toroidal(usize),
    /// Index into bspline_hypotheses (stage 2.7 NURBS fallback).
    BSpline(usize),
}


// ---------------------------------------------------------------------------
// Stage 2 output
// ---------------------------------------------------------------------------

/// The output of Stage 2: the mesh with all hypotheses populated and surfaces selected.
///
/// After stage 2.6 (surface selection), the hypothesis vectors contain only the
/// hypotheses that were actually selected. Unselected hypotheses are pruned.
/// Stage 2.1/2.2 results: the fit-free partition of the mesh
/// (DEVELOPMENT_PLAN §2.1–2.2).
#[derive(Debug, Default)]
pub struct Segmentation {
    /// Per-edge sharp/interior classes (stage 2.1(a)).
    pub edge_classes: EdgeClassification,
    /// Per-face osculating-jet prior (stage 2.1(b)), consumed by 2.3 for
    /// effort allocation. Empty when stage 2.1(b) was not run.
    pub jet: Vec<JetField>,
    /// Region id per mesh face (stage 2.2).
    pub face_region: Vec<usize>,
    /// Regions, ordered by smallest member face index (stage 2.2).
    pub regions: Vec<Region>,
}

#[derive(Debug)]
pub struct Stage2Output {
    /// The mesh from stage 1 with per-face hypothesis indices populated.
    pub mesh: ConnectedMesh,
    /// Selected planar hypotheses (after stage 2.6 pruning).
    pub planar_hypotheses: Vec<PlanarHypothesis>,
    /// Selected conical/cylindrical hypotheses (after stage 2.6 pruning).
    /// Cylinders are represented as cones with half_angle == 0.
    pub conical_hypotheses: Vec<ConicalHypothesis>,
    /// Selected spherical hypotheses (after stage 2.6 pruning).
    pub spherical_hypotheses: Vec<SphericalHypothesis>,
    /// Selected toroidal hypotheses (after stage 2.6 pruning).
    pub toroidal_hypotheses: Vec<ToroidalHypothesis>,
    /// B-spline hypotheses fitted by the stage 2.7 NURBS fallback. These own
    /// the actual OCCT surfaces (stage 3 reuses the same objects).
    pub bspline_hypotheses: Vec<BSplineHypothesis>,
    /// Surfaces selected in stage 2.6 for reconstruction. Each face is
    /// covered by exactly one selected surface. Indices point into the
    /// (pruned) hypothesis vectors above.
    pub selected_surfaces: Vec<SelectedSurface>,
    /// The fit-free partition from stages 2.1/2.2 (partition-then-classify
    /// rewrite). Empty (default) in output produced by the nerfed BFS
    /// pipeline.
    pub segmentation: Segmentation,
}

impl Stage2Output {
    /// The mesh faces owned by each selected surface, in `selected_surfaces`
    /// order — our partition of the mesh.
    pub fn selected_face_sets(&self) -> Vec<Vec<usize>> {
        self.selected_surfaces
            .iter()
            .map(|s| match s {
                SelectedSurface::Planar(i) => self.planar_hypotheses[*i].faces.clone(),
                SelectedSurface::Conical(i) => self.conical_hypotheses[*i].faces.clone(),
                SelectedSurface::Spherical(i) => self.spherical_hypotheses[*i].faces.clone(),
                SelectedSurface::Toroidal(i) => self.toroidal_hypotheses[*i].faces.clone(),
                SelectedSurface::BSpline(i) => self.bspline_hypotheses[*i].faces.clone(),
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Stage2Error {
    Compare(Stage2CompareError),
}

impl Display for Stage2Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage2Error::Compare(e) => write!(f, "stage 2 compare: {e}"),
        }
    }
}

impl Error for Stage2Error {}

impl From<Stage2CompareError> for Stage2Error {
    fn from(e: Stage2CompareError) -> Self {
        Stage2Error::Compare(e)
    }
}

#[derive(Debug)]
pub struct Stage2CompareError {
    pub hypothesis_type: &'static str,
    pub hypothesis_index: usize,
    pub max_distance: f64,
    pub tolerance: f64,
}

impl Display for Stage2CompareError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} hypothesis {} projected centroid is {:.6e} mm from nearest STEP surface (tolerance: {:.6e} mm)",
            self.hypothesis_type, self.hypothesis_index, self.max_distance, self.tolerance
        )
    }
}

impl Error for Stage2CompareError {}

// ---------------------------------------------------------------------------
// Shared math utilities
// ---------------------------------------------------------------------------

/// Compute the area of a triangular mesh face.
pub(crate) fn face_area(face: &MeshFace, vertices: &[MeshVertex]) -> f64 {
    let v0 = &vertices[face.vertex_indices[0]];
    let v1 = &vertices[face.vertex_indices[1]];
    let v2 = &vertices[face.vertex_indices[2]];
    let ax = v1.x - v0.x;
    let ay = v1.y - v0.y;
    let az = v1.z - v0.z;
    let bx = v2.x - v0.x;
    let by = v2.y - v0.y;
    let bz = v2.z - v0.z;
    let cx = ay * bz - az * by;
    let cy = az * bx - ax * bz;
    let cz = ax * by - ay * bx;
    let mut area = 0.5 * (cx * cx + cy * cy + cz * cz).sqrt();
    if face.vertex_count == 4 {
        let v3 = &vertices[face.vertex_indices[3]];
        let dx = v3.x - v0.x;
        let dy = v3.y - v0.y;
        let dz = v3.z - v0.z;
        let ex = by * dz - bz * dy;
        let ey = bz * dx - bx * dz;
        let ez = bx * dy - by * dx;
        area += 0.5 * (ex * ex + ey * ey + ez * ez).sqrt();
    }
    area
}

/// Compute the face centroid for a mesh face.
pub(crate) fn face_centroid(face: &MeshFace, vertices: &[MeshVertex]) -> [f64; 3] {
    let n = face.vertex_count as usize;
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    for vi_idx in 0..n {
        let v = &vertices[face.vertex_indices[vi_idx]];
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    let inv_n = 1.0 / n as f64;
    [cx * inv_n, cy * inv_n, cz * inv_n]
}

/// Normalize a 3D vector in place, return its length.
pub(crate) fn normalize3(v: &mut [f64; 3]) -> f64 {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 0.0 {
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
    }
    len
}

/// Cross product of two 3D vectors.
pub(crate) fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two 3D vectors.
pub(crate) fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute eigenvalues of a 3x3 symmetric matrix stored as upper triangle
/// [m00, m01, m02, m11, m12, m22]. Returns eigenvalues in unspecified order.
pub(crate) fn eigenvalues_3x3_symmetric(m: &[f64; 6]) -> [f64; 3] {
    let a00 = m[0]; let a01 = m[1]; let a02 = m[2];
    let a11 = m[3]; let a12 = m[4]; let a22 = m[5];

    let c2 = a00 + a11 + a22;
    let c1 = a00 * a11 - a01 * a01 + a00 * a22 - a02 * a02 + a11 * a22 - a12 * a12;
    let c0 = a00 * (a11 * a22 - a12 * a12)
            - a01 * (a01 * a22 - a12 * a02)
            + a02 * (a01 * a12 - a11 * a02);

    let p = c1 - c2 * c2 / 3.0;
    let q = -2.0 * c2 * c2 * c2 / 27.0 + c1 * c2 / 3.0 - c0;

    if p.abs() < 1e-30 {
        let ev = c2 / 3.0;
        return [ev, ev, ev];
    }

    let neg_p_3 = (-p / 3.0).max(0.0);
    let r = neg_p_3.sqrt();
    let cos_arg = (-q / (2.0 * neg_p_3 * r)).clamp(-1.0, 1.0);
    let theta = cos_arg.acos() / 3.0;
    let shift = c2 / 3.0;
    let two_pi_3 = 2.0 * std::f64::consts::PI / 3.0;
    [2.0 * r * theta.cos() + shift,
     2.0 * r * (theta - two_pi_3).cos() + shift,
     2.0 * r * (theta - 2.0 * two_pi_3).cos() + shift]
}

/// Compute the eigenvector of a 3x3 symmetric matrix for a given eigenvalue.
/// Uses cross products of (M - λI) rows to find the null space.
fn eigenvector_for_eigenvalue(m: &[f64; 6], lambda: f64) -> [f64; 3] {
    let b00 = m[0] - lambda;
    let b11 = m[3] - lambda;
    let b22 = m[5] - lambda;

    let row0 = [b00, m[1], m[2]];
    let row1 = [m[1], b11, m[4]];
    let row2 = [m[2], m[4], b22];

    let candidates = [
        cross3(&row0, &row1),
        cross3(&row0, &row2),
        cross3(&row1, &row2),
    ];

    let mut best_idx = 0;
    let mut best_len_sq = 0.0;
    for (i, c) in candidates.iter().enumerate() {
        let len_sq = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
        if len_sq > best_len_sq {
            best_len_sq = len_sq;
            best_idx = i;
        }
    }

    let mut ev = candidates[best_idx];
    let len = normalize3(&mut ev);
    if len < 1e-15 {
        return [1.0, 0.0, 0.0];
    }
    ev
}

/// Returns the eigenvector corresponding to the smallest eigenvalue.
pub(crate) fn smallest_eigenvector_3x3(m: &[f64; 6]) -> [f64; 3] {
    let eigenvalues = eigenvalues_3x3_symmetric(m);
    let (min_idx, _) = eigenvalues.iter().enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();
    eigenvector_for_eigenvalue(m, eigenvalues[min_idx])
}

/// Given a 3x3 symmetric matrix, return the eigenvector for the largest eigenvalue.
pub(crate) fn largest_eigenvector_3x3(m: &[f64; 6]) -> [f64; 3] {
    let eigenvalues = eigenvalues_3x3_symmetric(m);
    let (max_idx, _) = eigenvalues.iter().enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();
    eigenvector_for_eigenvalue(m, eigenvalues[max_idx])
}

/// Build the area-weighted normal covariance matrix M = Σ wᵢ nᵢ nᵢᵀ.
/// Returns the upper triangle [m00, m01, m02, m11, m12, m22].
pub(crate) fn build_normal_covariance(
    face_indices: &[usize],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> [f64; 6] {
    let mut m = [0.0_f64; 6];
    for &fi in face_indices {
        let face = &faces[fi];
        let n = face.normal.unwrap();
        let w = face_area(face, vertices);
        m[0] += w * n[0] * n[0];
        m[1] += w * n[0] * n[1];
        m[2] += w * n[0] * n[2];
        m[3] += w * n[1] * n[1];
        m[4] += w * n[1] * n[2];
        m[5] += w * n[2] * n[2];
    }
    m
}

/// Build an orthonormal basis (u, w) perpendicular to a unit direction `d`.
pub(crate) fn perpendicular_basis(d: &[f64; 3]) -> ([f64; 3], [f64; 3]) {
    let (ax, ay, az) = (d[0].abs(), d[1].abs(), d[2].abs());
    let perp = if ax <= ay && ax <= az {
        [1.0, 0.0, 0.0]
    } else if ay <= az {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let mut u = cross3(d, &perp);
    normalize3(&mut u);
    let w = cross3(d, &u);
    (u, w)
}

// ---------------------------------------------------------------------------
// BFS visualization helpers
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Stage 2.6: Select surfaces for reconstruction
// ---------------------------------------------------------------------------

/// Greedy area-based surface selection.
///
/// Selects hypotheses greedily by total mesh face area (largest first).
/// This naturally resolves conflicts: real surface hypotheses cover many faces
/// with large total area, while bogus hypotheses cover minimal area and lose.
///
/// After selection, updates hypothesis face/vertex lists to reflect only the
/// faces actually assigned to each selected surface.
pub(super) fn vertex_distance_to_selected(
    v: &MeshVertex,
    surface: &SelectedSurface,
    planar_hypotheses: &[PlanarHypothesis],
    conical_hypotheses: &[ConicalHypothesis],
    spherical_hypotheses: &[SphericalHypothesis],
    toroidal_hypotheses: &[ToroidalHypothesis],
    bspline_hypotheses: &[BSplineHypothesis],
) -> f64 {
    match surface {
        SelectedSurface::BSpline(i) => {
            bspline_hypotheses[*i].distance(&[v.x, v.y, v.z])
        }
        SelectedSurface::Planar(i) => {
            let h = &planar_hypotheses[*i];
            (h.normal[0] * v.x + h.normal[1] * v.y + h.normal[2] * v.z - h.distance).abs()
        }
        SelectedSurface::Conical(i) => {
            let h = &conical_hypotheses[*i];
            if h.is_cylinder() {
                cylindrical::vertex_to_cylinder_distance(
                    v, &h.apex, &h.axis_direction, h.radius_at_apex,
                ).abs()
            } else {
                conical::vertex_to_cone_distance(
                    v, &h.apex, &h.axis_direction, h.half_angle,
                ).abs()
            }
        }
        SelectedSurface::Spherical(i) => {
            let h = &spherical_hypotheses[*i];
            spherical::vertex_to_sphere_distance(v, &h.center, h.radius).abs()
        }
        SelectedSurface::Toroidal(i) => {
            let h = &toroidal_hypotheses[*i];
            toroidal::vertex_to_torus_distance(
                v, &h.center, &h.axis_direction, h.major_radius, h.minor_radius,
            ).abs()
        }
    }
}

/// Diagnostic: what the stage-2.4 boundary locus (`locus.rs`) concludes about
/// every adjacent pair of committed surfaces — whether the pair is tangent,
/// whether a junction locus could be built for it, and how many faces the
/// locus places on the other surface's side. One line per pair, sorted with
/// the problems first. Scratch-tool support for `dump_regions locus`; it reads
/// stage-2 output and changes nothing.
pub fn locus_report(out: &Stage2Output, config: &Config) -> Vec<String> {
    let (surfs, f2s) = locus_inputs(out);
    let mut pairs = locus::analyze_pairs(&out.mesh, &surfs, &f2s, config.vertex_tolerance_mm);
    pairs.sort_by_key(|p| (!p.is_problem(), p.a, p.b));
    pairs
        .iter()
        .map(|p| {
            format!(
                "{:>4} {:<7} | {:>4} {:<7} | {:>4} interface faces | {}",
                p.a,
                p.kind_a,
                p.b,
                p.kind_b,
                p.interface_faces,
                p.describe(),
            )
        })
        .collect()
}

/// How many adjacent surface pairs stage 2.4's boundary snap could not handle,
/// and how many faces are still on the wrong side of a junction it could.
///
/// Both should be zero on a well-decomposed model: the second is the property
/// snap exists to establish, and the first is the set of junction geometries
/// `locus.rs` does not implement. Exposed so tests can assert on it.
pub fn locus_summary(out: &Stage2Output, config: &Config) -> (usize, usize) {
    let (surfs, f2s) = locus_inputs(out);
    let pairs = locus::analyze_pairs(&out.mesh, &surfs, &f2s, config.vertex_tolerance_mm);
    let problems = pairs.iter().filter(|p| p.is_problem()).count();
    let wrong_side = pairs
        .iter()
        .map(|p| match &p.outcome {
            locus::PairOutcome::Snapped { moves, .. } => moves.len(),
            _ => 0,
        })
        .sum();
    (problems, wrong_side)
}

/// Diagnostic: per-face detail for one pair's shared boundary — which faces
/// the locus would move and, for the ones it would not, which test stopped it.
pub fn locus_pair_detail(out: &Stage2Output, config: &Config, a: usize, b: usize) -> Vec<String> {
    let (surfs, f2s) = locus_inputs(out);
    locus::explain_pair(&out.mesh, &surfs, &f2s, config.vertex_tolerance_mm, a, b)
}

/// The surface views and face→surface map the locus analysis works over.
fn locus_inputs<'a>(out: &'a Stage2Output) -> (Vec<locus::Surf<'a>>, Vec<usize>) {
    let surfs: Vec<locus::Surf<'_>> = out
        .selected_surfaces
        .iter()
        .map(|s| {
            locus::Surf::of(
                s,
                &out.planar_hypotheses,
                &out.conical_hypotheses,
                &out.spherical_hypotheses,
                &out.toroidal_hypotheses,
            )
        })
        .collect();

    let mut f2s = vec![usize::MAX; out.mesh.faces.len()];
    for (si, sel) in out.selected_surfaces.iter().enumerate() {
        for &fi in selected_surface_faces(out, sel) {
            f2s[fi] = si;
        }
    }

    (surfs, f2s)
}

/// The mesh faces a selected surface owns.
fn selected_surface_faces<'a>(out: &'a Stage2Output, sel: &SelectedSurface) -> &'a [usize] {
    match *sel {
        SelectedSurface::Planar(i) => &out.planar_hypotheses[i].faces,
        SelectedSurface::Conical(i) => &out.conical_hypotheses[i].faces,
        SelectedSurface::Spherical(i) => &out.spherical_hypotheses[i].faces,
        SelectedSurface::Toroidal(i) => &out.toroidal_hypotheses[i].faces,
        SelectedSurface::BSpline(i) => &out.bspline_hypotheses[i].faces,
    }
}

/// Refine surface assignments at region boundaries (DEVELOPMENT_PLAN §2.6
/// "Surface refitting"). BFS greedily absorbs border faces that fit within
/// tolerance but belong more naturally to the adjacent surface — e.g. a face
/// one row past a fillet junction absorbed by the plane. The misassignment
/// shifts the mesh region boundary by a row, which shifts reconstructed edge
/// boundary vertices and corner vertices off the true feature curves.
///
/// For each face on a region boundary, compare the max vertex distance to
/// its own surface against each adjacent selected surface, and move the face
/// when an adjacent surface fits at least 2× better. The margin prevents
/// oscillation for true junction faces (whose vertices lie on both surfaces
/// equally). Iterates to a fixpoint (boundaries only ever move locally).
/// Hypothesis geometry is NOT re-fitted: the fits were validated during BFS;
/// only the region map (face membership) changes.
#[allow(clippy::too_many_arguments)]
pub fn stage2(
    config: &Config,
    mesh: ConnectedMesh,
    viz: Option<&crate::viz::VizSender>,
) -> Result<Stage2Output, Stage2Error> {
    // --- Stage 2.1(a): sharp-edge classification (fit-free) ---
    let t = Instant::now();
    let edge_classes = edges::classify_edges(&mesh, config.angular_tolerance_rad);
    if config.verbosity >= 1 {
        let total = edge_classes.sharp_count + edge_classes.interior_count;
        eprintln!(
            "Stage 2.1(a) ({:.3}s): Sharp-edge classification: {} manifold edges — {} sharp, {} interior",
            t.elapsed().as_secs_f64(),
            total,
            edge_classes.sharp_count,
            edge_classes.interior_count,
        );
    }
    // --- Stage 2.1(b): osculating-jet curvature prior ---
    let t = Instant::now();
    let jet = edges::compute_jet_field(&mesh, &edge_classes);
    if config.verbosity >= 1 {
        eprintln!(
            "Stage 2.1(b) ({:.3}s): Jet curvature prior: {} faces summarized",
            t.elapsed().as_secs_f64(),
            jet.len(),
        );
    }
    let mut segmentation = Segmentation {
        edge_classes,
        jet,
        face_region: Vec::new(),
        regions: Vec::new(),
    };
    if config.viz_active(2, 1) {
        if let Some(vs) = viz {
            let overlay = segmentation_edges_overlay(&mesh, &segmentation.edge_classes);
            if matches!(vs.show_and_wait(overlay), VizAction::Quit) {
                return Ok(segmentation_only_output(mesh, segmentation));
            }
        }
    }
    if !config.stage.at_least(2, 2) {
        return Ok(segmentation_only_output(mesh, segmentation));
    }

    // --- Stage 2.2: region extraction ---
    let t = Instant::now();
    let (face_region, region_list) =
        regions::extract_regions(&mesh, &segmentation.edge_classes);
    if config.verbosity >= 1 {
        let single = region_list.iter().filter(|r| r.faces.len() == 1).count();
        let largest = region_list.iter().map(|r| r.faces.len()).max().unwrap_or(0);
        eprintln!(
            "Stage 2.2 ({:.3}s): Region extraction: {} regions over {} faces ({} single-face, largest {} faces)",
            t.elapsed().as_secs_f64(),
            region_list.len(),
            mesh.faces.len(),
            single,
            largest,
        );
        if config.verbosity >= 2 {
            for (i, r) in region_list.iter().enumerate() {
                if r.faces.len() > 1 {
                    eprintln!(
                        "  Region {}: {} faces, {} vertices (first face {})",
                        i,
                        r.faces.len(),
                        r.vertices.len(),
                        r.faces[0],
                    );
                }
            }
        }
    }
    segmentation.face_region = face_region;
    segmentation.regions = region_list;
    if config.viz_active(2, 2) {
        if let Some(vs) = viz {
            let overlay = segmentation_regions_overlay(&mesh, &segmentation);
            if matches!(vs.show_and_wait(overlay), VizAction::Quit) {
                return Ok(segmentation_only_output(mesh, segmentation));
            }
        }
    }

    // Invariant scoring (--invariant-report). Diagnostics only, and placed here
    // because this is the first point where BOTH halves of what it needs exist:
    // the jet and the smooth regions (2.2) and the untouched mesh. Running it
    // before 2.3 keeps it usable on the models whose decomposition fails, which
    // are the ones a fragility measurement most wants to see.
    if config.invariant_report {
        match config.compare_shape.as_ref() {
            Some(shape) => {
                let model = crate::reference::ReferenceModel::load(shape);
                let labels = crate::reference::label_faces(&mesh, &model);
                invariants::report(
                    config,
                    &mesh,
                    &model,
                    &labels,
                    &segmentation.jet,
                    &segmentation.edge_classes,
                    &segmentation.face_region,
                    &segmentation.regions,
                    viz,
                );
            }
            None => eprintln!("  Invariants: --invariant-report needs --compare; skipped"),
        }
    }

    if !config.stage.at_least(2, 3) {
        return Ok(segmentation_only_output(mesh, segmentation));
    }

    // --- Stage 2.3: per-region surface discovery (ladder + competition) ---
    let t = Instant::now();
    // Fit memoization (surface::with_fit_cache) — bit-identical, see there.
    // Off at verbosity >= 4, where the fitters print per-fit diagnostics that a
    // cache hit would silently skip.
    let memoize = config.verbosity < 4;
    let mut cls = surface::with_fit_cache(memoize, || {
        classify::classify_regions(
            &mesh,
            &segmentation.regions,
            &segmentation.jet,
            &segmentation.face_region,
            config,
        )
    });
    if config.verbosity >= 1 {
        eprintln!(
            "Stage 2.3 ({:.3}s): Region classification: {} surfaces ({} planar, {} cyl/cone, {} spherical, {} toroidal)",
            t.elapsed().as_secs_f64(),
            cls.selected.len(),
            cls.planar.len(),
            cls.conical.len(),
            cls.spherical.len(),
            cls.toroidal.len(),
        );
    }

    // --- Stage 2.4: resolution & repair (snap, relabel, merge, re-split) ---
    if config.stage.at_least(2, 4) {
        let t = Instant::now();
        let stats = surface::with_fit_cache(memoize, || repair::repair(&mesh, &mut cls, config));
        if config.verbosity >= 1 {
            eprintln!(
                "Stage 2.4 ({:.3}s): Repair: {} snapped ({} pair(s) unhandled), {} relabeled, {} merged, {} re-split → {} surfaces",
                t.elapsed().as_secs_f64(),
                stats.snapped,
                stats.snap_problems,
                stats.relabeled,
                stats.merged,
                stats.split,
                cls.selected.len(),
            );
        }
    }

    // Per-type breakdown of the final selection. Printed because surface
    // COUNTS hide regressions — a three-torus loss once hid behind an
    // unchanged total of 71 — and because obtaining this line used to require
    // a second full stage-2 run via `dump_regions surfaces`.
    if config.verbosity >= 1 {
        let (mut planes, mut cyls, mut cones, mut spheres, mut tori, mut bsplines) =
            (0, 0, 0, 0, 0, 0);
        for sel in &cls.selected {
            match sel {
                SelectedSurface::Planar(_) => planes += 1,
                SelectedSurface::Conical(i) => {
                    if cls.conical[*i].is_cylinder() {
                        cyls += 1;
                    } else {
                        cones += 1;
                    }
                }
                SelectedSurface::Spherical(_) => spheres += 1,
                SelectedSurface::Toroidal(_) => tori += 1,
                SelectedSurface::BSpline(_) => bsplines += 1,
            }
        }
        eprintln!(
            "Stage 2 types: {planes} plane, {cyls} cyl, {cones} cone, {spheres} sphere, \
             {tori} torus, {bsplines} bspline ({} total)",
            cls.selected.len(),
        );
    }

    // --- Stage 2.5: selected surfaces ---
    let output = Stage2Output {
        mesh,
        planar_hypotheses: cls.planar,
        conical_hypotheses: cls.conical,
        spherical_hypotheses: cls.spherical,
        toroidal_hypotheses: cls.toroidal,
        bspline_hypotheses: cls.bspline,
        selected_surfaces: cls.selected,
        segmentation,
    };

    // Segmentation scoring against the reference (--reference-report).
    // Diagnostics only. This is the half that needs our fits; the mesh-side half
    // runs in stage 1. It re-labels rather than carrying labels forward, because
    // `Config` is immutable through the pipeline — the same reason the Phase 0.6
    // quantum clamp is applied in the driver.
    if config.reference_report {
        if let Some(shape) = config.compare_shape.as_ref() {
            let model = crate::reference::ReferenceModel::load(shape);
            let labels = crate::reference::label_faces(&output.mesh, &model);
            crate::reference::score_segmentation(
                &output.mesh,
                &model,
                &labels,
                &output.selected_face_sets(),
            )
            .print();

            // Type agreement: the share of labelled mesh area on which our
            // surface is the same KIND of surface as the reference's. Coverage
            // and purity above cannot see a cone recognized as a sphere (one
            // surface for one face, both score 1); this can. And unlike a
            // comparison of face counts it does not care that a kernel splits
            // a cylinder into two half faces.
            let (mut agree, mut total) = (0.0_f64, 0.0_f64);
            for (surface, faces) in output.selected_surfaces.iter().zip(output.selected_face_sets()) {
                let ours = match surface {
                    SelectedSurface::Planar(_) => "plane",
                    SelectedSurface::Conical(i) if output.conical_hypotheses[*i].is_cylinder() => {
                        "cylinder"
                    }
                    SelectedSurface::Conical(_) => "cone",
                    SelectedSurface::Spherical(_) => "sphere",
                    SelectedSurface::Toroidal(_) => "torus",
                    SelectedSurface::BSpline(_) => "freeform",
                };
                for fi in faces {
                    let Some(ri) = labels[fi].face() else { continue };
                    let theirs = match &model.surfaces[ri] {
                        crate::reference::RefSurface::Other { .. } => "freeform",
                        analytic => analytic.type_name(),
                    };
                    let f = &output.mesh.faces[fi];
                    let p = |k: usize| {
                        let v = &output.mesh.vertices[f.vertex_indices[k]];
                        [v.x, v.y, v.z]
                    };
                    let mut area = 0.0;
                    for k in 1..(f.vertex_count as usize).saturating_sub(1) {
                        let (a, b, c) = (p(0), p(k), p(k + 1));
                        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                        let x = [
                            u[1] * w[2] - u[2] * w[1],
                            u[2] * w[0] - u[0] * w[2],
                            u[0] * w[1] - u[1] * w[0],
                        ];
                        area += 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
                    }
                    total += area;
                    if ours == theirs {
                        agree += area;
                    }
                }
            }
            if total > 0.0 {
                eprintln!(
                    "    type agreement {:.4} (share of labelled area whose surface type matches the reference's)",
                    agree / total
                );
            }
        }
    }

    Ok(output)
}

/// Stage2Output for the current rewrite frontier: segmentation only, no
/// hypotheses or selected surfaces (stages 3/4 are unreachable until the
/// classification ladder lands).
fn segmentation_only_output(mesh: ConnectedMesh, segmentation: Segmentation) -> Stage2Output {
    Stage2Output {
        mesh,
        planar_hypotheses: Vec::new(),
        conical_hypotheses: Vec::new(),
        spherical_hypotheses: Vec::new(),
        toroidal_hypotheses: Vec::new(),
        bspline_hypotheses: Vec::new(),
        selected_surfaces: Vec::new(),
        segmentation,
    }
}

// ---------------------------------------------------------------------------
// Segmentation viz (--viz-stages=2.1 / 2.2)
// ---------------------------------------------------------------------------

const SHARP_EDGE_COLOR: [f32; 4] = [0.9, 0.1, 0.1, 1.0];

/// Boundary-edge line segments for viz, one segment per undirected sharp
/// edge, offset slightly along the adjacent faces' average normal so they
/// sit visibly on top of the rendered mesh instead of z-fighting it.
/// Returns positions as [p0, p1, p0, p1, ...] pairs.
fn boundary_line_positions(
    mesh: &ConnectedMesh,
    ec: &EdgeClassification,
) -> Vec<[f32; 3]> {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for v in &mesh.vertices {
        for (k, c) in [v.x, v.y, v.z].into_iter().enumerate() {
            min[k] = min[k].min(c);
            max[k] = max[k].max(c);
        }
    }
    let diag = ((max[0] - min[0]).powi(2)
        + (max[1] - min[1]).powi(2)
        + (max[2] - min[2]).powi(2))
    .sqrt();
    let offset = if diag.is_finite() { diag * 5e-4 } else { 0.0 };

    let mut sharp = Vec::new();
    for fi in 0..mesh.faces.len() {
        let vc = mesh.faces[fi].vertex_count as usize;
        for e in 0..vc {
            let ni = mesh.faces[fi].neighbors[e];
            if ni < fi as i32 {
                continue; // open edge, or emitted from the other side
            }
            if ec.class[fi][e] != EdgeClass::Sharp {
                continue;
            }
            let n1 = mesh.faces[fi].normal.unwrap_or([0.0; 3]);
            let n2 = mesh.faces[ni as usize].normal.unwrap_or(n1);
            let mut nn = [n1[0] + n2[0], n1[1] + n2[1], n1[2] + n2[2]];
            let len = (nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]).sqrt();
            if len > 1e-12 {
                for c in &mut nn {
                    *c /= len;
                }
            } else {
                nn = [0.0; 3];
            }
            let va = &mesh.vertices[mesh.faces[fi].vertex_indices[e]];
            let vb = &mesh.vertices[mesh.faces[fi].vertex_indices[(e + 1) % vc]];
            let pa = [
                (va.x + nn[0] * offset) as f32,
                (va.y + nn[1] * offset) as f32,
                (va.z + nn[2] * offset) as f32,
            ];
            let pb = [
                (vb.x + nn[0] * offset) as f32,
                (vb.y + nn[1] * offset) as f32,
                (vb.z + nn[2] * offset) as f32,
            ];
            sharp.push(pa);
            sharp.push(pb);
        }
    }
    sharp
}

fn push_boundary_lines(overlay: &mut viz::VizOverlay, mesh: &ConnectedMesh, ec: &EdgeClassification) {
    let sharp = boundary_line_positions(mesh, ec);
    if !sharp.is_empty() {
        overlay.lines.push(viz::LineOverlay {
            positions: sharp,
            color: SHARP_EDGE_COLOR,
            no_depth_test: false,
        });
    }
}

fn segmentation_edges_overlay(mesh: &ConnectedMesh, ec: &EdgeClassification) -> viz::VizOverlay {
    let mut overlay = viz::VizOverlay::new();
    overlay.status_text = format!(
        "Stage 2.1: sharp-edge classification — {} sharp (red), {} interior [space=continue, q=quit]",
        ec.sharp_count, ec.interior_count,
    );
    push_boundary_lines(&mut overlay, mesh, ec);
    overlay
}

/// Deterministic distinct color for region `i`: golden-ratio hue stepping.
fn region_color(i: usize) -> [f32; 4] {
    let h = (i as f64 * 0.618_033_988_749_895).fract() * 6.0;
    let sector = (h.floor() as i32).rem_euclid(6);
    let f = h - h.floor();
    let (s, v) = (0.55, 0.9);
    let p = v * (1.0 - s);
    let q = v * (1.0 - s * f);
    let t = v * (1.0 - s * (1.0 - f));
    let (r, g, b) = match sector {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        _ => (v, p, q),
    };
    [r as f32, g as f32, b as f32, 1.0]
}

fn segmentation_regions_overlay(mesh: &ConnectedMesh, seg: &Segmentation) -> viz::VizOverlay {
    let mut overlay = viz::VizOverlay::new();
    let single = seg.regions.iter().filter(|r| r.faces.len() == 1).count();
    let largest = seg.regions.iter().map(|r| r.faces.len()).max().unwrap_or(0);
    overlay.status_text = format!(
        "Stage 2.2: {} regions ({} single-face, largest {} faces) — sharp boundaries in red [space=continue, q=quit]",
        seg.regions.len(),
        single,
        largest,
    );
    overlay.face_highlights = seg
        .regions
        .iter()
        .enumerate()
        .map(|(i, r)| viz::FaceHighlight {
            face_indices: r.faces.clone(),
            color: region_color(i),
        })
        .collect();
    push_boundary_lines(&mut overlay, mesh, &seg.edge_classes);
    overlay
}

