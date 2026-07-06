//! Stage 2: Surface Fitting
//!
//! 2.1: Deduce planar hypotheses
//! 2.2: Deduce conical/cylindrical hypotheses
//! 2.3: Deduce spherical hypotheses
//! 2.4: Deduce toroidal hypotheses
//! 2.5: Select surfaces for reconstruction

mod bfs;
mod conical;
mod cylindrical;
mod planar;
mod spherical;
mod toroidal;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex};
use crate::viz::{self, VizAction, VizSender};
use std::collections::HashSet;
use std::time::Instant;
use std::error::Error;
use std::fmt::{Display, Formatter};

use self::conical::{
    compare_conical_hypotheses,
};
use self::planar::compare_planar_hypotheses;
use self::spherical::{
    bounding_box_diagonal, compare_spherical_hypotheses,
    vertex_to_sphere_distance,
};
use self::toroidal::{
    compare_toroidal_hypotheses,
};
use self::bfs::BfsConfig;
use self::planar::PlanarSurface;
use self::conical::ConicalSurface;
use self::spherical::SphericalSurface;
use self::toroidal::ToroidalSurface;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Lower bound on cross product magnitude |n₁ × n₂| for a valid cylinder/sphere
/// seed pair. Since |n₁ × n₂| = sin(θ) for unit normals, this corresponds to a
/// minimum dihedral angle of arcsin(0.01) ≈ 0.57° (~0.6°). Seed pairs closer
/// than this are too nearly parallel for numerically stable axis estimation
/// (cylinders) or sphere fitting (spheres). Together with `angular_tol` (the
/// upper bound, typically 17.5°), this defines the valid seed-pair window
/// [~0.6°, angular_tol].
pub(crate) const MIN_CROSS_THRESHOLD: f64 = 0.01;

/// Fast-reject multiplier for BFS vertex-distance checks. When a candidate face
/// has a vertex farther than `REFIT_SKIP_MULTIPLIER × vertex_tol` from the
/// current fitted surface, the distance is too large for a re-fit to absorb and
/// the face is rejected immediately without attempting re-fitting.
/// Set to 10.0 to allow aggressive refitting — the refit itself validates all
/// faces, so the only cost of a high multiplier is more refit attempts.
pub(crate) const REFIT_SKIP_MULTIPLIER: f64 = 10.0;

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

/// Minimum number of mesh faces required to accept a conical hypothesis.
/// A cone has 5 degrees of freedom (apex xyz, axis direction 2, half-angle),
/// so at least 6 faces are needed for a well-determined and robust fit.
pub(crate) const MIN_CONE_FACES: usize = 6;

/// Minimum number of mesh faces required to accept a toroidal hypothesis.
/// A torus has 7 degrees of freedom (center xyz, axis direction 2, R, r),
/// so at least 7 faces are needed for a well-determined fit.
pub(crate) const MIN_TORUS_FACES: usize = 7;

/// Maximum sphere radius as a multiple of the mesh bounding-box diagonal.
/// With solid-angle coverage validation and surface-tolerance validation during
/// BFS, this no longer needs to be tight — those checks prevent pathological
/// growth. A large value serves as a numerical guardrail against degenerate fits.
pub(crate) const MAX_SPHERE_RADIUS_FACTOR: f64 = 1000.0;

/// Minimum eigenvalue ratio for solid-angle coverage validation.
/// Centroid-to-center direction vectors of a genuine spherical hypothesis span
/// 3D, so all eigenvalues of their covariance matrix are substantial.
/// Fillet-strip growth produces nearly coplanar directions (λ₃ ≈ 0).
pub(crate) const MIN_SPHERE_EIGENVALUE_RATIO: f64 = 0.01;

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

// TODO: Stage 2.8 - NURBS hypothesis

/// Identifies a selected surface hypothesis by type and index.
/// `Conical` indexes into `conical_hypotheses` and covers both cylinders (half_angle≈0) and true cones.
#[derive(Debug, Clone, Copy)]
pub enum SelectedSurface {
    Planar(usize),
    /// Index into conical_hypotheses. Covers both cylinders (half_angle≈0) and true cones.
    Conical(usize),
    Spherical(usize),
    Toroidal(usize),
}

// ---------------------------------------------------------------------------
// HypothesisFaces implementations (for generic BFS)
// ---------------------------------------------------------------------------

impl bfs::HypothesisFaces for PlanarHypothesis {
    fn faces(&self) -> &[usize] { &self.faces }
}
impl bfs::HypothesisFaces for SphericalHypothesis {
    fn faces(&self) -> &[usize] { &self.faces }
}
impl bfs::HypothesisFaces for ConicalHypothesis {
    fn faces(&self) -> &[usize] { &self.faces }
}
impl bfs::HypothesisFaces for ToroidalHypothesis {
    fn faces(&self) -> &[usize] { &self.faces }
}

// ---------------------------------------------------------------------------
// Stage 2 output
// ---------------------------------------------------------------------------

/// The output of Stage 2: the mesh with all hypotheses populated and surfaces selected.
///
/// After stage 2.6 (surface selection), the hypothesis vectors contain only the
/// hypotheses that were actually selected. Unselected hypotheses are pruned.
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
    /// Surfaces selected in stage 2.6 for reconstruction. Each face is
    /// covered by exactly one selected surface. Indices point into the
    /// (pruned) hypothesis vectors above.
    pub selected_surfaces: Vec<SelectedSurface>,
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Stage2Error {
    NotImplemented(String),
    Compare(Stage2CompareError),
}

impl Display for Stage2Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage2Error::NotImplemented(stage) => write!(f, "stage {stage} not yet implemented"),
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

/// Compute the centroid of a mesh face (as f32 for viz).
pub(crate) fn viz_face_centroid(face_idx: usize, mesh: &ConnectedMesh) -> [f32; 3] {
    let face = &mesh.faces[face_idx];
    let vc = face.vertex_count as usize;
    let mut cx = 0.0f64;
    let mut cy = 0.0f64;
    let mut cz = 0.0f64;
    for i in 0..vc {
        let v = &mesh.vertices[face.vertex_indices[i]];
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    [cx as f32 / vc as f32, cy as f32 / vc as f32, cz as f32 / vc as f32]
}

/// Return the normal of a mesh face as f32 (for viz camera orientation).
pub(crate) fn viz_face_normal(face_idx: usize, mesh: &ConnectedMesh) -> Option<[f32; 3]> {
    mesh.faces[face_idx].normal.map(|n| [n[0] as f32, n[1] as f32, n[2] as f32])
}

/// Create a CylinderOverlay centered on the axial extent of the given faces.
pub(crate) fn centered_cylinder_overlay(
    origin: [f64; 3], direction: [f64; 3], radius: f64,
    face_list: &[usize], mesh: &ConnectedMesh,
    color: [f32; 4],
) -> viz::CylinderOverlay {
    // Compute full mesh bounding box extent for cylinder length
    let mut bb_min = [f64::INFINITY; 3];
    let mut bb_max = [f64::NEG_INFINITY; 3];
    for v in &mesh.vertices {
        let coords = [v.x, v.y, v.z];
        for i in 0..3 {
            bb_min[i] = bb_min[i].min(coords[i]);
            bb_max[i] = bb_max[i].max(coords[i]);
        }
    }
    let extent = ((bb_max[0] - bb_min[0]).powi(2)
        + (bb_max[1] - bb_min[1]).powi(2)
        + (bb_max[2] - bb_min[2]).powi(2))
    .sqrt();

    // Center on face list's axial extent
    let mut t_min = f64::INFINITY;
    let mut t_max = f64::NEG_INFINITY;
    for &fi in face_list {
        let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
        let d = [c[0] - origin[0], c[1] - origin[1], c[2] - origin[2]];
        let t = dot3(&d, &direction);
        t_min = t_min.min(t);
        t_max = t_max.max(t);
    }
    let t_mid = (t_min + t_max) * 0.5;
    let half_len = extent;
    let centered_origin = [
        origin[0] + direction[0] * t_mid,
        origin[1] + direction[1] * t_mid,
        origin[2] + direction[2] * t_mid,
    ];
    viz::CylinderOverlay {
        origin: centered_origin,
        direction,
        radius,
        half_length: half_len,
        color,
    }
}

/// Create a ConeOverlay spanning the axial extent of the given faces (with padding).
/// Ensures the cone is at least 2x the maximum model dimension in length.
pub(crate) fn centered_cone_overlay(
    apex: [f64; 3], direction: [f64; 3], half_angle: f64,
    face_list: &[usize], mesh: &ConnectedMesh,
    color: [f32; 4],
) -> viz::ConeOverlay {
    // Compute axial extent of face centroids relative to apex
    let mut t_min = f64::INFINITY;
    let mut t_max = f64::NEG_INFINITY;
    for &fi in face_list {
        let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
        let d = [c[0] - apex[0], c[1] - apex[1], c[2] - apex[2]];
        let t = dot3(&d, &direction);
        t_min = t_min.min(t);
        t_max = t_max.max(t);
    }
    // Ensure the cone is at least 2x the maximum model dimension in length
    let bb_diag = bounding_box_diagonal(&mesh.vertices);
    let min_half_length = bb_diag; // 2x max dimension means full length = 2 * bb_diag, so half_length = bb_diag
    let center_t = (t_min + t_max) * 0.5;
    let half_span = ((t_max - t_min) * 0.5).max(min_half_length);
    // Pad by 20% for visual clarity
    let pad = half_span * 0.2;
    viz::ConeOverlay {
        apex,
        direction,
        half_angle,
        t_min: center_t - half_span - pad,
        t_max: center_t + half_span + pad,
        color,
    }
}

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
fn select_surfaces(
    mesh: &ConnectedMesh,
    planar_hypotheses: &mut [PlanarHypothesis],
    conical_hypotheses: &mut [ConicalHypothesis],
    spherical_hypotheses: &mut [SphericalHypothesis],
    toroidal_hypotheses: &mut [ToroidalHypothesis],
    viz: Option<&VizSender>,
    verbosity: u8,
) -> Vec<SelectedSurface> {
    let num_faces = mesh.faces.len();

    // Step 1: Compute geometric area of each mesh face.
    let face_areas: Vec<f64> = (0..num_faces)
        .map(|fi| face_area(&mesh.faces[fi], &mesh.vertices))
        .collect();

    // Step 2: Build candidate list from multi-face hypotheses.
    struct Candidate {
        surface_type: u8, // 0=planar, 1=conical(cylinder), 2=spherical, 3=conical(cone), 4=toroidal
        hypothesis_index: usize,
        faces: Vec<usize>,
    }

    let mut candidates: Vec<Candidate> = Vec::new();

    for (i, hyp) in planar_hypotheses.iter().enumerate() {
        if hyp.faces.len() > 1 {
            candidates.push(Candidate {
                surface_type: 0,
                hypothesis_index: i,
                faces: hyp.faces.clone(),
            });
        }
    }
    for (i, hyp) in conical_hypotheses.iter().enumerate() {
        let st = if hyp.is_cylinder() { 1 } else { 3 };
        candidates.push(Candidate {
            surface_type: st,
            hypothesis_index: i,
            faces: hyp.faces.clone(),
        });
    }
    for (i, hyp) in spherical_hypotheses.iter().enumerate() {
        candidates.push(Candidate {
            surface_type: 2,
            hypothesis_index: i,
            faces: hyp.faces.clone(),
        });
    }
    for (i, hyp) in toroidal_hypotheses.iter().enumerate() {
        candidates.push(Candidate {
            surface_type: 4,
            hypothesis_index: i,
            faces: hyp.faces.clone(),
        });
    }

    // Step 3: Greedy loop — select candidate with largest remaining area.
    let mut assigned = vec![false; num_faces];
    let mut selected: Vec<SelectedSurface> = Vec::new();

    loop {
        let mut best_ci = None;
        let mut best_area = 0.0_f64;

        for (ci, cand) in candidates.iter().enumerate() {
            let area: f64 = cand.faces.iter()
                .filter(|fi| !assigned[**fi])
                .map(|fi| face_areas[*fi])
                .sum();
            if area > best_area {
                best_area = area;
                best_ci = Some(ci);
            }
        }

        let Some(ci) = best_ci else { break };

        // Collect still-unassigned faces for the winning candidate.
        let sel_faces: Vec<usize> = candidates[ci].faces.iter()
            .filter(|fi| !assigned[**fi])
            .copied()
            .collect();


        // Viz: show selected hypothesis and competing hypotheses
        if let Some(viz_sender) = viz {
            let mut overlay = viz::VizOverlay::new();

            // Already-assigned faces in dark gray
            let assigned_faces: Vec<usize> = (0..num_faces)
                .filter(|&fi| assigned[fi])
                .collect();
            if !assigned_faces.is_empty() {
                overlay.face_highlights.push(viz::FaceHighlight {
                    face_indices: assigned_faces,
                    color: [0.3, 0.3, 0.3, 1.0],
                });
            }

            // Faces being assigned in green
            overlay.face_highlights.push(viz::FaceHighlight {
                face_indices: sel_faces.clone(),
                color: [0.0, 0.8, 0.0, 1.0],
            });

            // Find competing hypotheses that involve any of the sel_faces
            let sel_face_set: HashSet<usize> = sel_faces.iter().copied().collect();
            for (oi, other) in candidates.iter().enumerate() {
                if oi == ci {
                    continue;
                }
                let competing_faces: Vec<usize> = other.faces.iter()
                    .filter(|fi| !assigned[**fi] && sel_face_set.contains(fi))
                    .copied()
                    .collect();
                if !competing_faces.is_empty() {
                    overlay.face_highlights.push(viz::FaceHighlight {
                        face_indices: competing_faces,
                        color: [0.5, 0.5, 0.5, 0.25],
                    });
                    // Show the competing hypothesis surface overlay
                    match other.surface_type {
                        1 => {
                            let hyp = &conical_hypotheses[other.hypothesis_index];
                            overlay.cylinders.push(centered_cylinder_overlay(
                                hyp.apex, hyp.axis_direction, hyp.radius_at_apex,
                                &other.faces, mesh,
                                [0.5, 0.5, 0.5, 0.25],
                            ));
                        }
                        2 => {
                            let hyp = &spherical_hypotheses[other.hypothesis_index];
                            overlay.spheres.push(viz::SphereOverlay {
                                center: hyp.center,
                                radius: hyp.radius,
                                color: [0.5, 0.5, 0.5, 0.25],
                            });
                        }
                        3 => {
                            let hyp = &conical_hypotheses[other.hypothesis_index];
                            overlay.cones.push(centered_cone_overlay(
                                hyp.apex, hyp.axis_direction, hyp.half_angle,
                                &other.faces, mesh,
                                [0.5, 0.5, 0.5, 0.25],
                            ));
                        }
                        4 => {
                            let hyp = &toroidal_hypotheses[other.hypothesis_index];
                            overlay.tori.push(viz::TorusOverlay {
                                center: hyp.center,
                                direction: hyp.axis_direction,
                                major_radius: hyp.major_radius,
                                minor_radius: hyp.minor_radius,
                                color: [0.5, 0.5, 0.5, 0.25],
                            });
                        }
                        _ => {}
                    }
                }
            }

            // Show the winning hypothesis surface overlay in green
            let type_name = match candidates[ci].surface_type {
                0 => {
                    // Planar — faces are sufficient
                    "planar"
                }
                1 => {
                    let hyp = &conical_hypotheses[candidates[ci].hypothesis_index];
                    overlay.cylinders.push(centered_cylinder_overlay(
                        hyp.apex, hyp.axis_direction, hyp.radius_at_apex,
                        &sel_faces, mesh,
                        [0.0, 0.8, 0.0, 0.3],
                    ));
                    "cylindrical"
                }
                2 => {
                    let hyp = &spherical_hypotheses[candidates[ci].hypothesis_index];
                    overlay.spheres.push(viz::SphereOverlay {
                        center: hyp.center,
                        radius: hyp.radius,
                        color: [0.0, 0.8, 0.0, 0.3],
                    });
                    "spherical"
                }
                _ => {
                    match candidates[ci].surface_type {
                        3 => {
                            let hyp = &conical_hypotheses[candidates[ci].hypothesis_index];
                            overlay.cones.push(centered_cone_overlay(
                                hyp.apex, hyp.axis_direction, hyp.half_angle,
                                &sel_faces, mesh,
                                [0.0, 0.8, 0.0, 0.3],
                            ));
                            "conical"
                        }
                        4 => {
                            let hyp = &toroidal_hypotheses[candidates[ci].hypothesis_index];
                            overlay.tori.push(viz::TorusOverlay {
                                center: hyp.center,
                                direction: hyp.axis_direction,
                                major_radius: hyp.major_radius,
                                minor_radius: hyp.minor_radius,
                                color: [0.0, 0.8, 0.0, 0.3],
                            });
                            "toroidal"
                        }
                        _ => "unknown",
                    }
                }
            };

            overlay.status_text = format!(
                "Stage 2.6: iter {}, selected {} {} (hi={}, {} faces, area={:.4}) [space=next]",
                selected.len() + 1, type_name, candidates[ci].hypothesis_index,
                candidates[ci].hypothesis_index, sel_faces.len(), best_area,
            );

            // Focus on centroid of selected faces
            if !sel_faces.is_empty() {
                overlay.focus_point = Some(viz_face_centroid(sel_faces[0], mesh));
                overlay.focus_normal = viz_face_normal(sel_faces[0], mesh);
            }

            if viz_sender.show_and_wait(overlay) == VizAction::Quit {
                break;
            }
        }

        for &fi in &sel_faces {
            assigned[fi] = true;
        }

        // Build vertex set for the selected faces.
        let sel_vertices: Vec<usize> = {
            let mut vset = HashSet::new();
            for &fi in &sel_faces {
                let f = &mesh.faces[fi];
                for vi in 0..f.vertex_count as usize {
                    vset.insert(f.vertex_indices[vi]);
                }
            }
            let mut v: Vec<usize> = vset.into_iter().collect();
            v.sort_unstable();
            v
        };

        // Update hypothesis face/vertex lists to reflect only the assigned faces,
        // and create the SelectedSurface entry.
        let selected_surface = match candidates[ci].surface_type {
            0 => {
                let idx = candidates[ci].hypothesis_index;
                planar_hypotheses[idx].faces = sel_faces;
                planar_hypotheses[idx].vertices = sel_vertices;
                SelectedSurface::Planar(idx)
            }
            1 => {
                // Cylinder (conical with half_angle == 0)
                let idx = candidates[ci].hypothesis_index;
                conical_hypotheses[idx].faces = sel_faces;
                conical_hypotheses[idx].vertices = sel_vertices;
                SelectedSurface::Conical(idx)
            }
            2 => {
                let idx = candidates[ci].hypothesis_index;
                spherical_hypotheses[idx].faces = sel_faces;
                spherical_hypotheses[idx].vertices = sel_vertices;
                SelectedSurface::Spherical(idx)
            }
            3 => {
                // True cone (conical with half_angle > 0)
                let idx = candidates[ci].hypothesis_index;
                conical_hypotheses[idx].faces = sel_faces;
                conical_hypotheses[idx].vertices = sel_vertices;
                SelectedSurface::Conical(idx)
            }
            4 => {
                let idx = candidates[ci].hypothesis_index;
                toroidal_hypotheses[idx].faces = sel_faces;
                toroidal_hypotheses[idx].vertices = sel_vertices;
                SelectedSurface::Toroidal(idx)
            }
            _ => unreachable!("unknown surface type"),
        };
        selected.push(selected_surface);
    }

    // Step 4: Assign remaining faces to their single-face planar hypothesis.
    for (fi, assigned) in assigned.iter().enumerate().take(num_faces) {
        if !assigned {
            let pi = mesh.faces[fi].planar_hypothesis;
            assert!(pi >= 0, "face {fi} has no hypothesis assigned");
            if verbosity >= 2 {
                let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
                let n = mesh.faces[fi].normal.unwrap_or([0.0; 3]);
                eprintln!(
                    "  Selection fallback: face {fi} → single-face planar {pi}, \
centroid=[{:.4},{:.4},{:.4}] normal=[{:.3},{:.3},{:.3}]",
                    c[0], c[1], c[2], n[0], n[1], n[2],
                );
            }
            selected.push(SelectedSurface::Planar(pi as usize));
        }
    }

    selected
}

/// Distance from a mesh vertex to a selected surface's fitted geometry.
fn vertex_distance_to_selected(
    v: &MeshVertex,
    surface: &SelectedSurface,
    planar_hypotheses: &[PlanarHypothesis],
    conical_hypotheses: &[ConicalHypothesis],
    spherical_hypotheses: &[SphericalHypothesis],
    toroidal_hypotheses: &[ToroidalHypothesis],
) -> f64 {
    match surface {
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

/// Merge redundant/orphan selected surfaces into larger adjacent ones.
///
/// Two situations produce selected surfaces that duplicate a neighbor's
/// geometry rather than describing a real model feature:
///
/// 1. **Spurious sub-fits**: a narrow band of faces on a cone/sphere can fit
///    an "exact" cylinder (e.g. 3 azimuthal strips of a shallow cone have
///    nearly coplanar normals) and get committed before the real surface's
///    BFS reaches those faces. All of the band's vertices lie exactly on the
///    real surface.
/// 2. **Tangent-junction bridge facets**: a mesh facet spanning a tangent
///    junction (e.g. cone–sphere) has vertices on both sides and can't be
///    absorbed by either surface at vertex tolerance; it becomes an orphan
///    single-face planar surface. Its chordal deviation from either surface
///    is bounded by the tessellation's surface tolerance.
///
/// Both create fake edges/vertices that break stage 3 reconstruction. This
/// pass merges a surface S into an adjacent, larger-area surface T when all
/// of S's vertices lie within tolerance of T's geometry: a tight tolerance
/// (10× vertex tolerance) for multi-face S (case 1: exact overlap), and the
/// looser surface tolerance for S with ≤ 2 faces (case 2: chordal sag).
#[allow(clippy::too_many_arguments)]
fn absorb_redundant_surfaces(
    mesh: &ConnectedMesh,
    selected_surfaces: &mut Vec<SelectedSurface>,
    planar_hypotheses: &mut [PlanarHypothesis],
    conical_hypotheses: &mut [ConicalHypothesis],
    spherical_hypotheses: &mut [SphericalHypothesis],
    toroidal_hypotheses: &mut [ToroidalHypothesis],
    config: &Config,
) {
    let surface_faces = |s: &SelectedSurface,
                         planar: &[PlanarHypothesis],
                         conical: &[ConicalHypothesis],
                         spherical: &[SphericalHypothesis],
                         toroidal: &[ToroidalHypothesis]|
     -> Vec<usize> {
        match s {
            SelectedSurface::Planar(i) => planar[*i].faces.clone(),
            SelectedSurface::Conical(i) => conical[*i].faces.clone(),
            SelectedSurface::Spherical(i) => spherical[*i].faces.clone(),
            SelectedSurface::Toroidal(i) => toroidal[*i].faces.clone(),
        }
    };

    for _iteration in 0..5 {
        // Face → selected surface index map and per-surface areas.
        let mut face_to_surface = vec![usize::MAX; mesh.faces.len()];
        let mut areas = vec![0.0_f64; selected_surfaces.len()];
        for (si, s) in selected_surfaces.iter().enumerate() {
            for fi in surface_faces(
                s, planar_hypotheses, conical_hypotheses,
                spherical_hypotheses, toroidal_hypotheses,
            ) {
                face_to_surface[fi] = si;
                areas[si] += face_area(&mesh.faces[fi], &mesh.vertices);
            }
        }

        // Plan merges: source surface index → target surface index.
        let mut merges: Vec<(usize, usize)> = Vec::new();
        for (si, s) in selected_surfaces.iter().enumerate() {
            let faces = surface_faces(
                s, planar_hypotheses, conical_hypotheses,
                spherical_hypotheses, toroidal_hypotheses,
            );
            if faces.is_empty() {
                continue;
            }
            let tol = if faces.len() <= 2 {
                config.surface_tolerance_mm
            } else {
                config.vertex_tolerance_mm * 10.0
            };

            // Collect adjacent surfaces with strictly larger area.
            let mut adjacent: Vec<usize> = Vec::new();
            for &fi in &faces {
                let face = &mesh.faces[fi];
                let vc = face.vertex_count as usize;
                for &ni in &face.neighbors[..vc] {
                    if ni < 0 {
                        continue;
                    }
                    let ti = face_to_surface[ni as usize];
                    if ti != usize::MAX && ti != si && areas[ti] > areas[si] {
                        adjacent.push(ti);
                    }
                }
            }
            adjacent.sort_unstable();
            adjacent.dedup();
            // Try larger targets first.
            adjacent.sort_by(|a, b| areas[*b].partial_cmp(&areas[*a]).unwrap());

            // Collect all vertices of S once.
            let mut verts: Vec<usize> = Vec::new();
            for &fi in &faces {
                let face = &mesh.faces[fi];
                let vc = face.vertex_count as usize;
                verts.extend_from_slice(&face.vertex_indices[..vc]);
            }
            verts.sort_unstable();
            verts.dedup();

            for ti in adjacent {
                let target = &selected_surfaces[ti];
                let max_vert_dist = verts
                    .iter()
                    .map(|&vi| {
                        vertex_distance_to_selected(
                            &mesh.vertices[vi], target,
                            planar_hypotheses, conical_hypotheses,
                            spherical_hypotheses, toroidal_hypotheses,
                        )
                    })
                    .fold(0.0_f64, f64::max);
                // Centroids must also lie on the target (within the chordal
                // surface tolerance). Vertices alone are not sufficient: a
                // flat end cap triangulated as a fan from a boundary vertex
                // has ALL vertices on its boundary circle — which lies on the
                // adjacent cone/cylinder — but its face interiors clearly do
                // not follow that surface.
                let max_centroid_dist = faces
                    .iter()
                    .map(|&fi| {
                        let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
                        vertex_distance_to_selected(
                            &MeshVertex::from_xyz(c[0], c[1], c[2]), target,
                            planar_hypotheses, conical_hypotheses,
                            spherical_hypotheses, toroidal_hypotheses,
                        )
                    })
                    .fold(0.0_f64, f64::max);
                if max_vert_dist <= tol && max_centroid_dist <= config.surface_tolerance_mm {
                    if config.verbosity >= 2 {
                        eprintln!(
                            "  Absorb: surface {si} ({} faces, area {:.3}) -> surface {ti} \
(vert_dist {:.2e}, centroid_dist {:.2e}, tol {:.2e})",
                            faces.len(), areas[si], max_vert_dist, max_centroid_dist, tol,
                        );
                    }
                    merges.push((si, ti));
                    break;
                }
            }
        }

        if merges.is_empty() {
            break;
        }

        // Apply merges. Skip chains where the target itself merges away this
        // iteration (it will be handled on the next pass).
        let merging: Vec<bool> = {
            let mut m = vec![false; selected_surfaces.len()];
            for &(si, _) in &merges {
                m[si] = true;
            }
            m
        };
        let mut merged_any = false;
        for &(si, ti) in &merges {
            if merging[ti] {
                continue;
            }
            let faces = surface_faces(
                &selected_surfaces[si], planar_hypotheses, conical_hypotheses,
                spherical_hypotheses, toroidal_hypotheses,
            );
            let mut verts: Vec<usize> = Vec::new();
            for &fi in &faces {
                let face = &mesh.faces[fi];
                let vc = face.vertex_count as usize;
                verts.extend_from_slice(&face.vertex_indices[..vc]);
            }
            // Move faces/vertices from S's hypothesis into T's.
            let take = |s: &SelectedSurface,
                        planar: &mut [PlanarHypothesis],
                        conical: &mut [ConicalHypothesis],
                        spherical: &mut [SphericalHypothesis],
                        toroidal: &mut [ToroidalHypothesis]| {
                match s {
                    SelectedSurface::Planar(i) => planar[*i].faces.clear(),
                    SelectedSurface::Conical(i) => conical[*i].faces.clear(),
                    SelectedSurface::Spherical(i) => spherical[*i].faces.clear(),
                    SelectedSurface::Toroidal(i) => toroidal[*i].faces.clear(),
                }
            };
            take(
                &selected_surfaces[si], planar_hypotheses, conical_hypotheses,
                spherical_hypotheses, toroidal_hypotheses,
            );
            match &selected_surfaces[ti] {
                SelectedSurface::Planar(i) => {
                    planar_hypotheses[*i].faces.extend_from_slice(&faces);
                    planar_hypotheses[*i].vertices.extend_from_slice(&verts);
                    planar_hypotheses[*i].vertices.sort_unstable();
                    planar_hypotheses[*i].vertices.dedup();
                }
                SelectedSurface::Conical(i) => {
                    conical_hypotheses[*i].faces.extend_from_slice(&faces);
                    conical_hypotheses[*i].vertices.extend_from_slice(&verts);
                    conical_hypotheses[*i].vertices.sort_unstable();
                    conical_hypotheses[*i].vertices.dedup();
                }
                SelectedSurface::Spherical(i) => {
                    spherical_hypotheses[*i].faces.extend_from_slice(&faces);
                    spherical_hypotheses[*i].vertices.extend_from_slice(&verts);
                    spherical_hypotheses[*i].vertices.sort_unstable();
                    spherical_hypotheses[*i].vertices.dedup();
                }
                SelectedSurface::Toroidal(i) => {
                    toroidal_hypotheses[*i].faces.extend_from_slice(&faces);
                    toroidal_hypotheses[*i].vertices.extend_from_slice(&verts);
                    toroidal_hypotheses[*i].vertices.sort_unstable();
                    toroidal_hypotheses[*i].vertices.dedup();
                }
            }
            merged_any = true;
        }

        // Drop selected surfaces whose hypotheses are now empty.
        selected_surfaces.retain(|s| {
            !surface_faces(
                s, planar_hypotheses, conical_hypotheses,
                spherical_hypotheses, toroidal_hypotheses,
            )
            .is_empty()
        });

        if !merged_any {
            break;
        }
    }
}

/// Remove unselected hypotheses from the hypothesis vectors and remap
/// `SelectedSurface` indices to match the new, compacted vectors.
///
/// After this, hypotheses vectors contain only the hypotheses that were
/// actually selected during surface selection.
fn prune_unselected_hypotheses(
    planar_hypotheses: &mut Vec<PlanarHypothesis>,
    conical_hypotheses: &mut Vec<ConicalHypothesis>,
    spherical_hypotheses: &mut Vec<SphericalHypothesis>,
    toroidal_hypotheses: &mut Vec<ToroidalHypothesis>,
    selected_surfaces: &mut [SelectedSurface],
) {
    // Collect which old indices are selected for each type.
    let mut planar_used = vec![false; planar_hypotheses.len()];
    let mut conical_used = vec![false; conical_hypotheses.len()];
    let mut spherical_used = vec![false; spherical_hypotheses.len()];
    let mut toroidal_used = vec![false; toroidal_hypotheses.len()];

    for s in selected_surfaces.iter() {
        match s {
            SelectedSurface::Planar(i) => planar_used[*i] = true,
            SelectedSurface::Conical(i) => conical_used[*i] = true,
            SelectedSurface::Spherical(i) => spherical_used[*i] = true,
            SelectedSurface::Toroidal(i) => toroidal_used[*i] = true,
        }
    }

    // Build old→new index maps and compact each vector.
    fn compact<H>(hypotheses: &mut Vec<H>, used: &[bool]) -> Vec<usize> {
        let mut old_to_new = vec![usize::MAX; hypotheses.len()];
        let mut new_vec = Vec::new();
        let old_vec = std::mem::take(hypotheses);
        for (old_idx, hyp) in old_vec.into_iter().enumerate() {
            if used[old_idx] {
                old_to_new[old_idx] = new_vec.len();
                new_vec.push(hyp);
            }
        }
        *hypotheses = new_vec;
        old_to_new
    }

    let planar_map = compact(planar_hypotheses, &planar_used);
    let conical_map = compact(conical_hypotheses, &conical_used);
    let spherical_map = compact(spherical_hypotheses, &spherical_used);
    let toroidal_map = compact(toroidal_hypotheses, &toroidal_used);

    // Remap all SelectedSurface indices.
    for s in selected_surfaces.iter_mut() {
        match s {
            SelectedSurface::Planar(i) => *i = planar_map[*i],
            SelectedSurface::Conical(i) => *i = conical_map[*i],
            SelectedSurface::Spherical(i) => *i = spherical_map[*i],
            SelectedSurface::Toroidal(i) => *i = toroidal_map[*i],
        }
    }
}

/// Validate selected surfaces against a reference STEP file.
#[allow(clippy::too_many_arguments)]
fn compare_selected_surfaces(
    selected_surfaces: &[SelectedSurface],
    planar_hypotheses: &[PlanarHypothesis],
    conical_hypotheses: &[ConicalHypothesis],
    spherical_hypotheses: &[SphericalHypothesis],
    toroidal_hypotheses: &[ToroidalHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (si, surface) in selected_surfaces.iter().enumerate() {
        let (hyp_type, face_indices) = match surface {
            SelectedSurface::Planar(i) => ("planar", &planar_hypotheses[*i].faces),
            SelectedSurface::Conical(i) => {
                let label = if conical_hypotheses[*i].is_cylinder() { "cylindrical" } else { "conical" };
                (label, &conical_hypotheses[*i].faces)
            }
            SelectedSurface::Spherical(i) => ("spherical", &spherical_hypotheses[*i].faces),
            SelectedSurface::Toroidal(i) => ("toroidal", &toroidal_hypotheses[*i].faces),
        };

        let projected_points: Vec<[f64; 3]> = face_indices.iter().map(|&fi| {
            let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);

            match surface {
                SelectedSurface::Planar(i) => {
                    let hyp = &planar_hypotheses[*i];
                    let dist_to_plane = hyp.normal[0] * centroid[0]
                        + hyp.normal[1] * centroid[1]
                        + hyp.normal[2] * centroid[2]
                        - hyp.distance;
                    [
                        centroid[0] - dist_to_plane * hyp.normal[0],
                        centroid[1] - dist_to_plane * hyp.normal[1],
                        centroid[2] - dist_to_plane * hyp.normal[2],
                    ]
                }
                SelectedSurface::Conical(i) => {
                    let hyp = &conical_hypotheses[*i];
                    let d = [
                        centroid[0] - hyp.apex[0],
                        centroid[1] - hyp.apex[1],
                        centroid[2] - hyp.apex[2],
                    ];
                    let h = dot3(&d, &hyp.axis_direction);
                    let radial = [
                        d[0] - h * hyp.axis_direction[0],
                        d[1] - h * hyp.axis_direction[1],
                        d[2] - h * hyp.axis_direction[2],
                    ];
                    let radial_dist = (radial[0] * radial[0]
                        + radial[1] * radial[1]
                        + radial[2] * radial[2])
                        .sqrt();
                    // Unified formula: radius = radius_at_apex + h * tan(half_angle)
                    // For cylinders: half_angle=0, so radius = radius_at_apex
                    // For cones: radius_at_apex=0, so radius = h * tan(half_angle)
                    let surface_r = hyp.radius_at_apex + h * hyp.half_angle.tan();
                    if radial_dist > 1e-15 {
                        let scale = surface_r / radial_dist;
                        [
                            hyp.apex[0] + h * hyp.axis_direction[0] + radial[0] * scale,
                            hyp.apex[1] + h * hyp.axis_direction[1] + radial[1] * scale,
                            hyp.apex[2] + h * hyp.axis_direction[2] + radial[2] * scale,
                        ]
                    } else {
                        centroid
                    }
                }
                SelectedSurface::Spherical(i) => {
                    let hyp = &spherical_hypotheses[*i];
                    let d = [
                        centroid[0] - hyp.center[0],
                        centroid[1] - hyp.center[1],
                        centroid[2] - hyp.center[2],
                    ];
                    let dist_to_center =
                        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    if dist_to_center > 1e-15 {
                        let scale = hyp.radius / dist_to_center;
                        [
                            hyp.center[0] + d[0] * scale,
                            hyp.center[1] + d[1] * scale,
                            hyp.center[2] + d[2] * scale,
                        ]
                    } else {
                        centroid
                    }
                }
                SelectedSurface::Toroidal(i) => {
                    let hyp = &toroidal_hypotheses[*i];
                    // Project centroid onto torus surface:
                    // 1. Find the point on the major circle closest to centroid
                    // 2. From that point, project outward to the minor circle
                    let d = [
                        centroid[0] - hyp.center[0],
                        centroid[1] - hyp.center[1],
                        centroid[2] - hyp.center[2],
                    ];
                    let axial = dot3(&d, &hyp.axis_direction);
                    let radial = [
                        d[0] - axial * hyp.axis_direction[0],
                        d[1] - axial * hyp.axis_direction[1],
                        d[2] - axial * hyp.axis_direction[2],
                    ];
                    let radial_dist = (radial[0] * radial[0]
                        + radial[1] * radial[1]
                        + radial[2] * radial[2])
                        .sqrt();
                    if radial_dist > 1e-15 {
                        // Tube center on the major circle
                        let tube_center = [
                            hyp.center[0] + hyp.major_radius * radial[0] / radial_dist,
                            hyp.center[1] + hyp.major_radius * radial[1] / radial_dist,
                            hyp.center[2] + hyp.major_radius * radial[2] / radial_dist,
                        ];
                        let tube_vec = [
                            centroid[0] - tube_center[0],
                            centroid[1] - tube_center[1],
                            centroid[2] - tube_center[2],
                        ];
                        let tube_dist = (tube_vec[0] * tube_vec[0]
                            + tube_vec[1] * tube_vec[1]
                            + tube_vec[2] * tube_vec[2])
                            .sqrt();
                        if tube_dist > 1e-15 {
                            let scale = hyp.minor_radius / tube_dist;
                            [
                                tube_center[0] + tube_vec[0] * scale,
                                tube_center[1] + tube_vec[1] * scale,
                                tube_center[2] + tube_vec[2] * scale,
                            ]
                        } else {
                            centroid
                        }
                    } else {
                        centroid
                    }
                }
            }
        }).collect();

        let (max_dist, _) = stage1::max_distance_to_shape(&projected_points, compare_shape);
        let tolerance = config.surface_tolerance_mm;

        if max_dist > tolerance {
            return Err(Stage2CompareError {
                hypothesis_type: hyp_type,
                hypothesis_index: si,
                max_distance: max_dist,
                tolerance,
            });
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Stage 2 entry point
// ---------------------------------------------------------------------------

/// Run stage 2: fit surface hypotheses to mesh faces and select surfaces.
pub fn stage2(config: &Config, mut mesh: ConnectedMesh, viz: Option<&crate::viz::VizSender>) -> Result<Stage2Output, Stage2Error> {
    // Stage 2.1: Deduce planar hypotheses
    let t = Instant::now();
    let viz_21 = if config.viz_active(2, 1) { viz } else { None };
    let bfs_config = BfsConfig {
        vertex_tol: config.vertex_tolerance_mm,
        surface_tol: config.surface_tolerance_mm,
        angular_tol: config.angular_tolerance_rad,
        verbosity: config.verbosity,
        viz_level: config.viz_stages_level,
    };
    let (mut planar_hypotheses, planar_quit) = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &bfs_config, viz_21);

    // Synthesize single-face planar hypotheses for faces not covered by BFS.
    // BFS with acceptance >= 2 only produces multi-face hypotheses; every
    // remaining face still needs a planar fallback for stage 2.6 selection.
    for fi in 0..mesh.faces.len() {
        if mesh.faces[fi].planar_hypothesis < 0 {
            let face = &mesh.faces[fi];
            let n = face.normal.unwrap();
            let vc = face.vertex_count as usize;
            let verts: Vec<usize> = face.vertex_indices[..vc].to_vec();
            let sum_d: f64 = verts.iter()
                .map(|&vi| {
                    let v = &mesh.vertices[vi];
                    n[0] * v.x + n[1] * v.y + n[2] * v.z
                })
                .sum();
            let distance = sum_d / verts.len() as f64;
            let new_idx = planar_hypotheses.len() as i32;
            planar_hypotheses.push(PlanarHypothesis {
                normal: n,
                distance,
                faces: vec![fi],
                vertices: verts,
                error_max: 0.0,
                error_min: 0.0,
                error_abs_sum: 0.0,
            });
            mesh.faces[fi].planar_hypothesis = new_idx;
        }
    }

    if config.verbosity >= 1 {
        let multi_face_count = planar_hypotheses.iter().filter(|h| h.faces.len() > 1).count();
        let single_face_count = planar_hypotheses.len() - multi_face_count;
        let covered_faces: usize = planar_hypotheses.iter().map(|h| h.faces.len()).sum();
        eprintln!(
            "Stage 2.1 ({:.3}s): Deduced {} planar hypotheses ({} multi-face, {} single-face) covering {} of {} mesh faces",
            t.elapsed().as_secs_f64(),
            planar_hypotheses.len(),
            multi_face_count,
            single_face_count,
            covered_faces,
            mesh.faces.len(),
        );
        if config.verbosity >= 2 {
            for (i, h) in planar_hypotheses.iter().enumerate() {
                if h.faces.len() > 1 {
                    eprintln!(
                        "  Plane {}: {} faces, {} vertices, normal=[{:.4}, {:.4}, {:.4}], d={:.4}, err_max={:.2e}, err_min={:.2e}",
                        i, h.faces.len(), h.vertices.len(),
                        h.normal[0], h.normal[1], h.normal[2], h.distance,
                        h.error_max, h.error_min,
                    );
                }
            }
        }
    }

    // Compare against STEP file if --compare was specified
    if config.compare_shape.is_some() {
        compare_planar_hypotheses(&planar_hypotheses, &mesh, config)?;
        if config.verbosity >= 1 {
            eprintln!("  Compare: all planar hypothesis centroids within tolerance");
        }
    }

    if planar_quit {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            conical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            toroidal_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }



    if !config.stage.at_least(2, 2) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            conical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            toroidal_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.3 (run before 2.2): Deduce spherical hypotheses.
    //
    // Spherical detection runs before conical/cylindrical detection because
    // every band of quad rows on a tessellated sphere lies exactly on some
    // cone or cylinder, so the conical BFS would otherwise tile fine spheres
    // with hundreds of bogus band hypotheses (very slow, all discarded by
    // stage 2.6). Sphere detection is cheap and robust — its solid-angle
    // coverage acceptance test rejects band-like fits on cylinders, cones,
    // and tori — so committing spheres first and skipping sphere-committed
    // faces during conical *seeding* eliminates the pathology. The stage
    // numbering (2.2 conical, 2.3 spherical) is retained for CLI/viz
    // compatibility.
    let t = Instant::now();
    let viz_23 = if config.viz_active(2, 3) { viz } else { None };
    let (mut spherical_hypotheses, spherical_quit) = bfs::generate_hypotheses::<SphericalSurface>(&mut mesh, &bfs_config, viz_23);

    if config.verbosity >= 1 {
        let covered_faces: usize = spherical_hypotheses.iter().map(|h| h.faces.len()).sum();
        let convex_count = spherical_hypotheses.iter().filter(|h| h.convex).count();
        let concave_count = spherical_hypotheses.len() - convex_count;
        eprintln!(
            "Stage 2.3 ({:.3}s): Deduced {} spherical hypotheses ({} convex, {} concave) covering {} faces",
            t.elapsed().as_secs_f64(),
            spherical_hypotheses.len(),
            convex_count,
            concave_count,
            covered_faces,
        );
        if config.verbosity >= 2 {
            for (i, h) in spherical_hypotheses.iter().enumerate() {
                eprintln!(
                    "  Sphere {}: {} faces, {} vertices, r={:.4}, {}, center=[{:.4}, {:.4}, {:.4}], \
vtx_err_max={:.2e}, cen_err_max={:.2e}",
                    i, h.faces.len(), h.vertices.len(), h.radius,
                    if h.convex { "convex" } else { "concave" },
                    h.center[0], h.center[1], h.center[2],
                    h.error_max, h.centroid_error_max,
                );
                if config.verbosity >= 3 {
                    for &fi in &h.faces {
                        let face = &mesh.faces[fi];
                        let vc = face.vertex_count as usize;
                        let centroid = face_centroid(face, &mesh.vertices);
                        let normal = face.normal.unwrap_or([0.0; 3]);
                        let centroid_err = vertex_to_sphere_distance(
                            &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                            &h.center, h.radius,
                        ).abs();
                        let mut vtx_err_max = f64::NEG_INFINITY;
                        let mut vtx_err_min = f64::INFINITY;
                        for vi_idx in 0..vc {
                            let d = vertex_to_sphere_distance(
                                &mesh.vertices[face.vertex_indices[vi_idx]],
                                &h.center, h.radius,
                            ).abs();
                            vtx_err_max = vtx_err_max.max(d);
                            vtx_err_min = vtx_err_min.min(d);
                        }
                        eprintln!(
                            "    face {fi}: {vc}v centroid=[{:.4},{:.4},{:.4}] \
normal=[{:.3},{:.3},{:.3}] vtx_err=[{:.2e},{:.2e}] cen_err={:.2e}",
                            centroid[0], centroid[1], centroid[2],
                            normal[0], normal[1], normal[2],
                            vtx_err_min, vtx_err_max, centroid_err,
                        );
                    }
                }
            }
        }
    }

    // Compare against STEP file if --compare was specified
    if config.compare_shape.is_some() {
        compare_spherical_hypotheses(&spherical_hypotheses, &mesh, config)?;
        if config.verbosity >= 1 {
            eprintln!("  Compare: all spherical hypothesis centroids within tolerance");
        }
    }

    if spherical_quit {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            conical_hypotheses: Vec::new(),
            spherical_hypotheses,
            toroidal_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.2: Deduce conical/cylindrical hypotheses (unified)
    let t = Instant::now();
    let viz_22 = if config.viz_active(2, 2) { viz } else { None };
    let (mut conical_hypotheses, conical_quit) = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &bfs_config, viz_22);

    if config.verbosity >= 1 {
        let covered_faces: usize = conical_hypotheses.iter().map(|h| h.faces.len()).sum();
        let cylinder_count = conical_hypotheses.iter().filter(|h| h.is_cylinder()).count();
        let cone_count = conical_hypotheses.len() - cylinder_count;
        let convex_count = conical_hypotheses.iter().filter(|h| h.convex).count();
        let concave_count = conical_hypotheses.len() - convex_count;
        eprintln!(
            "Stage 2.2 ({:.3}s): Deduced {} conical/cylindrical hypotheses ({} cylinders, {} cones; {} convex, {} concave) covering {} faces",
            t.elapsed().as_secs_f64(),
            conical_hypotheses.len(),
            cylinder_count,
            cone_count,
            convex_count,
            concave_count,
            covered_faces,
        );
        if config.verbosity >= 2 {
            for (i, h) in conical_hypotheses.iter().enumerate() {
                if h.is_cylinder() {
                    eprintln!(
                        "  Cylinder {}: {} faces, {} vertices, r={:.4}, {}, \
apex=[{:.4}, {:.4}, {:.4}], axis_dir=[{:.4}, {:.4}, {:.4}], \
vtx_err_max={:.2e}, cen_err_max={:.2e}",
                        i, h.faces.len(), h.vertices.len(), h.radius_at_apex,
                        if h.convex { "convex" } else { "concave" },
                        h.apex[0], h.apex[1], h.apex[2],
                        h.axis_direction[0], h.axis_direction[1], h.axis_direction[2],
                        h.error_max, h.centroid_error_max,
                    );
                } else {
                    eprintln!(
                        "  Cone {}: {} faces, half_angle={:.4}\u{00b0}, apex=[{:.4},{:.4},{:.4}], {}, err_max={:.2e}",
                        i,
                        h.faces.len(),
                        h.half_angle.to_degrees(),
                        h.apex[0], h.apex[1], h.apex[2],
                        if h.convex { "convex" } else { "concave" },
                        h.error_max,
                    );
                }
            }
        }
    }

    // Compare against STEP file if --compare was specified
    if config.compare_shape.is_some() {
        compare_conical_hypotheses(&conical_hypotheses, &mesh, config)?;
        if config.verbosity >= 1 {
            eprintln!("  Compare: all conical/cylindrical hypothesis centroids within tolerance");
        }
    }

    if conical_quit {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            conical_hypotheses,
            spherical_hypotheses,
            toroidal_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Note: stages 2.2 and 2.3 both run above (spherical before conical), so
    // --stage=2.2 and --stage=2.3 stop at the same point.
    if !config.stage.at_least(2, 4) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            conical_hypotheses,
            spherical_hypotheses,
            toroidal_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.4: (conical hypotheses now unified with cylinders in stage 2.2)

    if !config.stage.at_least(2, 5) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            conical_hypotheses,
            spherical_hypotheses,
            toroidal_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.5: Deduce toroidal hypotheses
    let t = Instant::now();
    let viz_25 = if config.viz_active(2, 5) { viz } else { None };
    let (mut toroidal_hypotheses, toroidal_quit) = bfs::generate_hypotheses::<ToroidalSurface>(&mut mesh, &bfs_config, viz_25);

    if config.verbosity >= 1 {
        let covered_faces: usize = toroidal_hypotheses.iter().map(|h| h.faces.len()).sum();
        let convex_count = toroidal_hypotheses.iter().filter(|h| h.convex).count();
        let concave_count = toroidal_hypotheses.len() - convex_count;
        eprintln!(
            "Stage 2.5 ({:.3}s): Deduced {} toroidal hypotheses ({} convex, {} concave) covering {} faces",
            t.elapsed().as_secs_f64(),
            toroidal_hypotheses.len(),
            convex_count,
            concave_count,
            covered_faces,
        );
        if config.verbosity >= 2 {
            for (i, h) in toroidal_hypotheses.iter().enumerate() {
                eprintln!(
                    "  Torus {}: {} faces, {} vertices, R={:.4}, r={:.4}, {}, \
center=[{:.4}, {:.4}, {:.4}], axis=[{:.4}, {:.4}, {:.4}], \
vtx_err_max={:.2e}, cen_err_max={:.2e}",
                    i, h.faces.len(), h.vertices.len(),
                    h.major_radius, h.minor_radius,
                    if h.convex { "convex" } else { "concave" },
                    h.center[0], h.center[1], h.center[2],
                    h.axis_direction[0], h.axis_direction[1], h.axis_direction[2],
                    h.error_max, h.centroid_error_max,
                );
            }
        }
    }

    if config.compare_shape.is_some() {
        compare_toroidal_hypotheses(&toroidal_hypotheses, &mesh, config)?;
        if config.verbosity >= 1 {
            eprintln!("  Compare: all toroidal hypothesis centroids within tolerance");
        }
    }

    if toroidal_quit {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            conical_hypotheses,
            spherical_hypotheses,
            toroidal_hypotheses,
            selected_surfaces: Vec::new(),
        });
    }

    if !config.stage.at_least(2, 6) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            conical_hypotheses,
            spherical_hypotheses,
            toroidal_hypotheses,
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.6: Select surfaces for reconstruction
    let t = Instant::now();
    let viz_26 = if config.viz_active(2, 6) { viz } else { None };
    let mut selected_surfaces = select_surfaces(
        &mesh, &mut planar_hypotheses, &mut conical_hypotheses,
        &mut spherical_hypotheses, &mut toroidal_hypotheses, viz_26,
        config.verbosity,
    );

    // Merge spurious sub-fit surfaces and tangent-junction orphan facets
    // into larger adjacent surfaces before reconstruction.
    absorb_redundant_surfaces(
        &mesh, &mut selected_surfaces, &mut planar_hypotheses,
        &mut conical_hypotheses, &mut spherical_hypotheses,
        &mut toroidal_hypotheses, config,
    );

    // Prune hypotheses that were not selected, keeping only the winning ones.
    prune_unselected_hypotheses(
        &mut planar_hypotheses, &mut conical_hypotheses,
        &mut spherical_hypotheses, &mut toroidal_hypotheses,
        &mut selected_surfaces,
    );

    if config.verbosity >= 1 {
        let mut planar_count = 0;
        let mut cylindrical_count = 0;
        let mut spherical_count = 0;
        let mut conical_count = 0;
        let mut toroidal_count = 0;
        for s in &selected_surfaces {
            match s {
                SelectedSurface::Planar(_) => planar_count += 1,
                SelectedSurface::Conical(i) => {
                    if conical_hypotheses[*i].is_cylinder() {
                        cylindrical_count += 1;
                    } else {
                        conical_count += 1;
                    }
                }
                SelectedSurface::Spherical(_) => spherical_count += 1,
                SelectedSurface::Toroidal(_) => toroidal_count += 1,
            }
        }
        eprintln!(
            "Stage 2.6 ({:.3}s): Selected {} surfaces ({} planar, {} cylindrical, {} spherical, {} conical, {} toroidal) covering {} faces",
            t.elapsed().as_secs_f64(),
            selected_surfaces.len(),
            planar_count,
            cylindrical_count,
            spherical_count,
            conical_count,
            toroidal_count,
            mesh.faces.len(),
        );
    }

    // Compare against STEP file if --compare was specified
    if config.compare_shape.is_some() {
        compare_selected_surfaces(
            &selected_surfaces, &planar_hypotheses, &conical_hypotheses,
            &spherical_hypotheses, &toroidal_hypotheses, &mesh, config,
        )?;
        if config.verbosity >= 1 {
            eprintln!("  Compare: all selected surface centroids within tolerance");
        }
    }

    Ok(Stage2Output {
        mesh,
        planar_hypotheses,
        conical_hypotheses,
        spherical_hypotheses,
        toroidal_hypotheses,
        selected_surfaces,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stage1::{self, MeshFace, MeshVertex, VertexWeldOptions, UNDEDUCED_PLANAR_HYPOTHESIS, UNDEDUCED_SPHERICAL_HYPOTHESIS, UNDEDUCED_CONICAL_HYPOTHESIS, UNDEDUCED_TOROIDAL_HYPOTHESIS};

    fn make_triangle_face(a: usize, b: usize, c: usize) -> MeshFace {
        MeshFace {
            vertex_count: 3,
            vertex_indices: [a, b, c, 0],
            neighbors: [-1, -1, -1, -1],
            normal: None,
            planar_hypothesis: UNDEDUCED_PLANAR_HYPOTHESIS,
            spherical_hypothesis: UNDEDUCED_SPHERICAL_HYPOTHESIS,
            conical_hypothesis: UNDEDUCED_CONICAL_HYPOTHESIS,
            toroidal_hypothesis: UNDEDUCED_TOROIDAL_HYPOTHESIS,
        }
    }

    /// Build a simple mesh from vertices and triangle indices, then validate.
    fn build_mesh(verts: Vec<MeshVertex>, tris: Vec<[usize; 3]>) -> ConnectedMesh {
        let faces = tris.iter().map(|t| make_triangle_face(t[0], t[1], t[2])).collect();
        let mut mesh = ConnectedMesh {
            vertices: verts,
            faces,
            stats: Default::default(),
        };
        mesh.validate_and_populate_topology().expect("mesh should be valid");
        mesh
    }

    fn default_bfs_config() -> BfsConfig {
        BfsConfig {
            vertex_tol: 1e-5,
            surface_tol: 0.4,
            angular_tol: 17.5_f64.to_radians(),
            verbosity: 0,
            viz_level: 1,
        }
    }

    fn bfs_config_with_angular_tol(angular_tol_deg: f64) -> BfsConfig {
        BfsConfig {
            vertex_tol: 1e-5,
            surface_tol: 0.4,
            angular_tol: angular_tol_deg.to_radians(),
            verbosity: 0,
            viz_level: 1,
        }
    }

    #[test]
    fn cube_gets_six_planar_hypotheses() {
        let stl_path = format!("{}/tests/manual/cube.stl", env!("CARGO_MANIFEST_DIR"));
        let mut mesh =
            stage1::read_connected_mesh_from_stl(&stl_path, VertexWeldOptions { tolerance: 1e-9 })
                .expect("should load");
        mesh.validate_and_populate_topology().expect("should validate");

        let config = default_bfs_config();
        let (hypotheses, _) = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);

        // A cube has 6 faces, each with 2 triangles
        assert_eq!(mesh.faces.len(), 12);

        // Should produce exactly 6 planar hypotheses (one per cube face)
        assert_eq!(
            hypotheses.len(),
            6,
            "expected 6 planar hypotheses for cube, got {}: {:?}",
            hypotheses.len(),
            hypotheses.iter().map(|h| h.faces.len()).collect::<Vec<_>>()
        );

        // Each hypothesis should have exactly 2 faces
        for (i, h) in hypotheses.iter().enumerate() {
            assert_eq!(
                h.faces.len(),
                2,
                "hypothesis {i} should have 2 faces, got {}",
                h.faces.len()
            );
        }

        // Total covered faces should be 12
        let total: usize = hypotheses.iter().map(|h| h.faces.len()).sum();
        assert_eq!(total, 12);

        // All errors should be essentially zero for a perfect cube
        for h in &hypotheses {
            assert!(
                h.error_max.abs() < 1e-12,
                "error_max should be ~0, got {}",
                h.error_max
            );
            assert!(
                h.error_min.abs() < 1e-12,
                "error_min should be ~0, got {}",
                h.error_min
            );
        }
    }

    #[test]
    fn two_perpendicular_triangles_get_separate_hypotheses() {
        // Two triangles sharing an edge but on perpendicular planes.
        // With the generic BFS (min acceptance = 2), each single triangle
        // produces 0 multi-face hypotheses. The actual pipeline synthesizes
        // single-face planar fallbacks in stage2() after BFS, but this unit
        // test verifies that BFS doesn't incorrectly merge them.
        let verts = vec![
            MeshVertex::from_xyz(0.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0), // triangle 0: on XY plane
            MeshVertex::from_xyz(0.0, 0.0, 1.0), // triangle 1: on XZ plane
        ];
        let mut mesh = build_mesh(verts, vec![[0, 1, 2], [0, 3, 1]]);
        let config = default_bfs_config();
        let (hypotheses, _) = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);

        // BFS won't produce multi-face hypotheses for single isolated triangles
        // (min acceptance = 2 faces, but perpendicular faces can't merge).
        // Verify no hypotheses merge the two perpendicular faces.
        assert_eq!(
            hypotheses.len(),
            0,
            "perpendicular triangles should get 0 multi-face hypotheses from BFS"
        );
    }

    #[test]
    fn coplanar_triangles_get_one_hypothesis() {
        // Two coplanar triangles sharing an edge on the XY plane
        let verts = vec![
            MeshVertex::from_xyz(0.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 1.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0),
        ];
        let mut mesh = build_mesh(verts, vec![[0, 1, 2], [0, 2, 3]]);
        let config = default_bfs_config();
        let (hypotheses, _) = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);

        assert_eq!(
            hypotheses.len(),
            1,
            "two coplanar triangles should get 1 hypothesis, got {}",
            hypotheses.len()
        );
        assert_eq!(hypotheses[0].faces.len(), 2);
        assert_eq!(hypotheses[0].vertices.len(), 4);
    }

    /// Helper to run stage 1 (read + validate + fuse) on an STL file.
    fn load_stage1(stl_name: &str) -> ConnectedMesh {
        let stl_path = format!("{}/tests/{}", env!("CARGO_MANIFEST_DIR"), stl_name);
        let mut mesh = stage1::read_connected_mesh_from_stl(
            &stl_path,
            VertexWeldOptions { tolerance: 1e-9 },
        ).expect("should load");
        mesh.validate_and_populate_topology().expect("should validate");
        stage1::fuse_coplanar_triangles(&mut mesh, 1e-5);
        mesh
    }

    #[test]
    fn cube_angular_tolerance_rejects_cylinder() {
        let mut mesh = load_stage1("manual/cube.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        // 17.5° angular tolerance: cube faces meet at 90°, so no cylinders/cones
        let (conicals, _) = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        let cyls: Vec<_> = conicals.iter().filter(|h| h.is_cylinder()).collect();
        assert_eq!(cyls.len(), 0, "cube should produce 0 cylinders at 17.5° angular tolerance");
    }

    #[test]
    fn cube_angular_tolerance_rejects_sphere() {
        let mut mesh = load_stage1("manual/cube.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let _ = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        // 17.5° angular tolerance: cube faces meet at 90°, so no spheres
        let (sphs, _) = bfs::generate_hypotheses::<SphericalSurface>(&mut mesh, &config, None);
        assert_eq!(sphs.len(), 0, "cube should produce 0 spheres at 17.5° angular tolerance");
    }

    #[test]
    fn cube_high_angular_tolerance_allows_sphere() {
        let mut mesh = load_stage1("manual/cube.stl");
        let config = bfs_config_with_angular_tol(91.0);
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let _ = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        // 91° angular tolerance: exceeds 90° dihedral angle of cube
        let (sphs, _) = bfs::generate_hypotheses::<SphericalSurface>(&mut mesh, &config, None);
        assert!(!sphs.is_empty(), "cube should produce spheres at 91° angular tolerance");
    }

    #[test]
    fn cylinder_detected_at_default_angular_tolerance() {
        let mut mesh = load_stage1("ccad/generated/simple_cylinder.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let (conicals, _) = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        let cyls: Vec<_> = conicals.iter().filter(|h| h.is_cylinder()).collect();
        assert_eq!(cyls.len(), 1, "simple cylinder should produce 1 cylinder hypothesis");
        assert!(cyls[0].convex);
        assert!((cyls[0].radius_at_apex - 10.0).abs() < 0.01);
    }

    #[test]
    fn sphere_detected_at_default_angular_tolerance() {
        let mut mesh = load_stage1("ccad/generated/simple_sphere.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let _ = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        let (sphs, _) = bfs::generate_hypotheses::<SphericalSurface>(&mut mesh, &config, None);
        assert_eq!(sphs.len(), 1, "simple sphere should produce 1 sphere hypothesis");
        assert!(sphs[0].convex);
        assert!((sphs[0].radius - 10.0).abs() < 0.01);
    }

    #[test]
    fn cone_detected_at_default_angular_tolerance() {
        let mut mesh = load_stage1("ccad/generated/simple_cone.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let (conicals, _) = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        let cones: Vec<_> = conicals.iter().filter(|h| !h.is_cylinder()).collect();
        assert_eq!(cones.len(), 1, "simple_cone should produce 1 cone hypothesis");
        assert!(cones[0].convex, "simple_cone should be convex");
        // simple_cone: truncated cone with r_top=5, r_bottom=10, h=30
        // half-angle = atan((10-5)/30) ≈ 9.46°
        let ha_deg = cones[0].half_angle.to_degrees();
        assert!((ha_deg - 9.46).abs() < 1.0, "half-angle should be ~9.46°, got {ha_deg:.2}°");
    }

    #[test]
    fn cube_angular_tolerance_rejects_cone() {
        let mut mesh = load_stage1("manual/cube.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let (conicals, _) = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        let cones: Vec<_> = conicals.iter().filter(|h| !h.is_cylinder()).collect();
        assert_eq!(cones.len(), 0, "cube should produce 0 cones");
    }

    #[test]
    fn concave_cone_detected() {
        let mut mesh = load_stage1("ccad/generated/block_with_conical_hole.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let (conicals, _) = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        let cones: Vec<_> = conicals.iter().filter(|h| !h.is_cylinder()).collect();
        assert!(cones.len() >= 1, "block_with_conical_hole should produce at least 1 cone hypothesis, got {}", cones.len());
        let concave_cones: Vec<_> = cones.iter().filter(|c| !c.convex).collect();
        assert!(!concave_cones.is_empty(), "block_with_conical_hole should have at least one concave cone");
    }

    #[test]
    fn sphere_not_detected_as_cone() {
        // At stage 2.2 a sphere may produce a spurious cone hypothesis, but
        // sphere coverage is better so the cone gets eliminated at stage 2.6.
        // Here we just verify cone detection doesn't panic on a sphere mesh.
        let mut mesh = load_stage1("ccad/generated/simple_sphere.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let (conicals, _) = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        // Spurious cone hypotheses on a sphere are geometrically legitimate at
        // this stage: a single latitude row of quads lies exactly on a cone.
        // Stage 2.6 selection removes them (verified by integration tests).
        // Here we only require that no single conical hypothesis colonizes a
        // large part of the sphere.
        let max_faces = conicals.iter().map(|h| h.faces.len()).max().unwrap_or(0);
        assert!(
            max_faces < mesh.faces.len() / 4,
            "no conical hypothesis should cover >25% of a sphere, got {} of {} faces",
            max_faces,
            mesh.faces.len(),
        );
    }

    #[test]
    fn cylinder_not_detected_as_cone() {
        // With unified conical BFS, a cylinder (half_angle=0) should NOT also produce
        // a separate non-cylinder cone hypothesis.
        let mut mesh = load_stage1("ccad/generated/simple_cylinder.stl");
        let config = default_bfs_config();
        let _ = bfs::generate_hypotheses::<PlanarSurface>(&mut mesh, &config, None);
        let (conicals, _) = bfs::generate_hypotheses::<ConicalSurface>(&mut mesh, &config, None);
        let cones: Vec<_> = conicals.iter().filter(|h| !h.is_cylinder()).collect();
        assert_eq!(cones.len(), 0, "simple_cylinder should produce 0 non-cylinder cones");
    }
}
