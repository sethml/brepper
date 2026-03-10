//! Stage 2: Surface Fitting
//!
//! 2.1: Deduce planar hypotheses
//! 2.2: Deduce cylindrical hypotheses
//! 2.3: Deduce spherical hypotheses
//! 2.4: Deduce ruled surface hypotheses
//! 2.5: Deduce NURBS hypotheses
//! 2.6: Select surfaces for reconstruction

mod cylindrical;
mod planar;
mod spherical;

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex};
use crate::viz::{self, VizAction, VizSender};
use opencascade_sys::gp;
use std::collections::HashSet;
use std::error::Error;
use std::fmt::{Display, Formatter};

use self::cylindrical::{
    compare_cylindrical_hypotheses, deduce_cylindrical_hypotheses, vertex_to_cylinder_distance,
};
use self::planar::{compare_planar_hypotheses, deduce_planar_hypotheses};
use self::spherical::{
    bounding_box_diagonal, compare_spherical_hypotheses, deduce_spherical_hypotheses,
    vertex_to_sphere_distance,
};

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
pub(crate) const REFIT_SKIP_MULTIPLIER: f64 = 2.0;

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

/// A cylindrical surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct CylindricalHypothesis {
    /// A point on the cylinder axis.
    pub axis_origin: [f64; 3],
    /// Unit direction vector along the cylinder axis.
    pub axis_direction: [f64; 3],
    /// Radius of the cylinder (always positive).
    pub radius: f64,
    /// Whether the surface normal points away from the axis (convex=true) or toward it (concave=false).
    pub convex: bool,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this cylinder.
    pub vertices: Vec<usize>,
    /// Maximum absolute distance from any vertex to the cylinder surface.
    pub error_max: f64,
    /// Maximum absolute distance from any face centroid to the cylinder surface.
    pub centroid_error_max: f64,
    /// Sum of absolute vertex-to-surface distances.
    pub error_abs_sum: f64,
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

// TODO: Stage 2.4 - Ruled surface hypothesis
// TODO: Stage 2.5 - NURBS hypothesis

/// Identifies a selected surface hypothesis by type and index.
#[derive(Debug, Clone, Copy)]
pub enum SelectedSurface {
    Planar(usize),
    Cylindrical(usize),
    Spherical(usize),
    // TODO: RuledSurface(usize),
    // TODO: Nurbs(usize),
}

// ---------------------------------------------------------------------------
// Stage 2 output
// ---------------------------------------------------------------------------

/// The output of Stage 2: the mesh with all hypotheses populated and surfaces selected.
#[derive(Debug)]
pub struct Stage2Output {
    /// The mesh from stage 1 with per-face hypothesis indices populated.
    pub mesh: ConnectedMesh,
    /// All planar hypotheses deduced in stage 2.1.
    pub planar_hypotheses: Vec<PlanarHypothesis>,
    /// All cylindrical hypotheses deduced in stage 2.2.
    pub cylindrical_hypotheses: Vec<CylindricalHypothesis>,
    /// All spherical hypotheses deduced in stage 2.3.
    pub spherical_hypotheses: Vec<SphericalHypothesis>,
    /// Surfaces selected in stage 2.6 for reconstruction. Each face should be
    /// covered by exactly one selected surface.
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

/// Send a BFS step visualization overlay and return the user's action.
#[allow(clippy::too_many_arguments)]
pub(crate) fn viz_bfs_step(
    viz: Option<&VizSender>,
    seed_faces: &[usize],
    hyp_faces: &[usize],
    explored_face: usize,
    status: &str,
    cylinders: Vec<viz::CylinderOverlay>,
    spheres: Vec<viz::SphereOverlay>,
    background_faces: &[usize],
    mesh: &ConnectedMesh,
) -> Option<VizAction> {
    let viz_sender = viz?;
    let mut overlay = viz::VizOverlay::new();
    overlay.status_text = status.to_string();
    // Background (already-hypothesized) faces in gray
    if !background_faces.is_empty() {
        overlay.face_highlights.push(viz::FaceHighlight {
            face_indices: background_faces.to_vec(),
            color: [0.5, 0.5, 0.5, 1.0],
        });
    }
    // Seed faces in green
    if !seed_faces.is_empty() {
        overlay.face_highlights.push(viz::FaceHighlight {
            face_indices: seed_faces.to_vec(),
            color: [0.0, 0.8, 0.0, 1.0],
        });
    }
    // Hypothesis faces (accepted so far) in blue
    if !hyp_faces.is_empty() {
        overlay.face_highlights.push(viz::FaceHighlight {
            face_indices: hyp_faces.to_vec(),
            color: [0.2, 0.4, 1.0, 1.0],
        });
    }
    // Face being explored in yellow
    overlay.face_highlights.push(viz::FaceHighlight {
        face_indices: vec![explored_face],
        color: [1.0, 0.9, 0.0, 1.0],
    });
    overlay.cylinders = cylinders;
    overlay.spheres = spheres;
    overlay.focus_point = Some(viz_face_centroid(explored_face, mesh));
    overlay.focus_normal = viz_face_normal(explored_face, mesh);
    Some(viz_sender.show_and_wait(overlay))
}

/// Send a BFS seed visualization (before expansion starts).
pub(crate) fn viz_bfs_seed(
    viz: Option<&VizSender>,
    seed_faces: &[usize],
    status: &str,
    cylinders: Vec<viz::CylinderOverlay>,
    spheres: Vec<viz::SphereOverlay>,
    background_faces: &[usize],
    mesh: &ConnectedMesh,
) -> Option<VizAction> {
    let viz_sender = viz?;
    let mut overlay = viz::VizOverlay::new();
    overlay.status_text = status.to_string();
    // Background faces in gray
    if !background_faces.is_empty() {
        overlay.face_highlights.push(viz::FaceHighlight {
            face_indices: background_faces.to_vec(),
            color: [0.5, 0.5, 0.5, 1.0],
        });
    }
    overlay.face_highlights.push(viz::FaceHighlight {
        face_indices: seed_faces.to_vec(),
        color: [0.0, 0.8, 0.0, 1.0],
    });
    overlay.cylinders = cylinders;
    overlay.spheres = spheres;
    // Focus on first seed face
    if !seed_faces.is_empty() {
        overlay.focus_point = Some(viz_face_centroid(seed_faces[0], mesh));
        overlay.focus_normal = viz_face_normal(seed_faces[0], mesh);
    }
    Some(viz_sender.show_and_wait(overlay))
}

/// Send a custom viz overlay with full control over highlights.
#[allow(clippy::too_many_arguments)]
pub(crate) fn viz_custom(
    viz: Option<&VizSender>,
    highlights: Vec<viz::FaceHighlight>,
    edge_highlights: Vec<viz::EdgeHighlight>,
    status: &str,
    cylinders: Vec<viz::CylinderOverlay>,
    spheres: Vec<viz::SphereOverlay>,
    focus_point: Option<[f32; 3]>,
    focus_normal: Option<[f32; 3]>,
) -> Option<VizAction> {
    let viz_sender = viz?;
    let mut overlay = viz::VizOverlay::new();
    overlay.status_text = status.to_string();
    overlay.face_highlights = highlights;
    overlay.edge_highlights = edge_highlights;
    overlay.cylinders = cylinders;
    overlay.spheres = spheres;
    overlay.focus_point = focus_point;
    overlay.focus_normal = focus_normal;
    Some(viz_sender.show_and_wait(overlay))
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
    cylindrical_hypotheses: &mut [CylindricalHypothesis],
    spherical_hypotheses: &mut [SphericalHypothesis],
    viz: Option<&VizSender>,
) -> Vec<SelectedSurface> {
    let num_faces = mesh.faces.len();

    // Step 1: Compute geometric area of each mesh face.
    let face_areas: Vec<f64> = (0..num_faces)
        .map(|fi| face_area(&mesh.faces[fi], &mesh.vertices))
        .collect();

    // Step 2: Build candidate list from multi-face hypotheses.
    struct Candidate {
        surface_type: u8, // 0=planar, 1=cylindrical, 2=spherical
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
    for (i, hyp) in cylindrical_hypotheses.iter().enumerate() {
        candidates.push(Candidate {
            surface_type: 1,
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
                            let hyp = &cylindrical_hypotheses[other.hypothesis_index];
                            overlay.cylinders.push(centered_cylinder_overlay(
                                hyp.axis_origin, hyp.axis_direction, hyp.radius,
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
                    let hyp = &cylindrical_hypotheses[candidates[ci].hypothesis_index];
                    overlay.cylinders.push(centered_cylinder_overlay(
                        hyp.axis_origin, hyp.axis_direction, hyp.radius,
                        &sel_faces, mesh,
                        [0.0, 0.8, 0.0, 0.3],
                    ));
                    "cylindrical"
                }
                _ => {
                    let hyp = &spherical_hypotheses[candidates[ci].hypothesis_index];
                    overlay.spheres.push(viz::SphereOverlay {
                        center: hyp.center,
                        radius: hyp.radius,
                        color: [0.0, 0.8, 0.0, 0.3],
                    });
                    "spherical"
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
                let idx = candidates[ci].hypothesis_index;
                cylindrical_hypotheses[idx].faces = sel_faces;
                cylindrical_hypotheses[idx].vertices = sel_vertices;
                SelectedSurface::Cylindrical(idx)
            }
            _ => {
                let idx = candidates[ci].hypothesis_index;
                spherical_hypotheses[idx].faces = sel_faces;
                spherical_hypotheses[idx].vertices = sel_vertices;
                SelectedSurface::Spherical(idx)
            }
        };
        selected.push(selected_surface);
    }

    // Step 4: Assign remaining faces to their single-face planar hypothesis.
    for (fi, assigned) in assigned.iter().enumerate().take(num_faces) {
        if !assigned {
            let pi = mesh.faces[fi].planar_hypothesis;
            assert!(pi >= 0, "face {fi} has no hypothesis assigned");
            selected.push(SelectedSurface::Planar(pi as usize));
        }
    }

    selected
}

/// Validate selected surfaces against a reference STEP file.
fn compare_selected_surfaces(
    selected_surfaces: &[SelectedSurface],
    planar_hypotheses: &[PlanarHypothesis],
    cylindrical_hypotheses: &[CylindricalHypothesis],
    spherical_hypotheses: &[SphericalHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (si, surface) in selected_surfaces.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        let (hyp_type, face_indices) = match surface {
            SelectedSurface::Planar(i) => ("planar", &planar_hypotheses[*i].faces),
            SelectedSurface::Cylindrical(i) => ("cylindrical", &cylindrical_hypotheses[*i].faces),
            SelectedSurface::Spherical(i) => ("spherical", &spherical_hypotheses[*i].faces),
        };

        for &fi in face_indices {
            let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);

            let projected = match surface {
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
                SelectedSurface::Cylindrical(i) => {
                    let hyp = &cylindrical_hypotheses[*i];
                    let d = [
                        centroid[0] - hyp.axis_origin[0],
                        centroid[1] - hyp.axis_origin[1],
                        centroid[2] - hyp.axis_origin[2],
                    ];
                    let t = dot3(&d, &hyp.axis_direction);
                    let radial = [
                        d[0] - t * hyp.axis_direction[0],
                        d[1] - t * hyp.axis_direction[1],
                        d[2] - t * hyp.axis_direction[2],
                    ];
                    let radial_dist = (radial[0] * radial[0]
                        + radial[1] * radial[1]
                        + radial[2] * radial[2])
                        .sqrt();
                    if radial_dist > 1e-15 {
                        let scale = hyp.radius / radial_dist;
                        [
                            hyp.axis_origin[0] + t * hyp.axis_direction[0] + radial[0] * scale,
                            hyp.axis_origin[1] + t * hyp.axis_direction[1] + radial[1] * scale,
                            hyp.axis_origin[2] + t * hyp.axis_direction[2] + radial[2] * scale,
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
            };

            let pt = gp::Pnt::new_real3(projected[0], projected[1], projected[2]);
            let dist = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(dist);
        }

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
    let viz_21 = if config.viz_active(2, 1) { viz } else { None };
    let (mut planar_hypotheses, planar_quit) = deduce_planar_hypotheses(&mut mesh, config.vertex_tolerance_mm, config.verbosity, viz_21);

    if !config.quiet {
        let multi_face_count = planar_hypotheses.iter().filter(|h| h.faces.len() > 1).count();
        let single_face_count = planar_hypotheses.len() - multi_face_count;
        let covered_faces: usize = planar_hypotheses.iter().map(|h| h.faces.len()).sum();
        eprintln!(
            "Stage 2.1: Deduced {} planar hypotheses ({} multi-face, {} single-face) covering {} of {} mesh faces",
            planar_hypotheses.len(),
            multi_face_count,
            single_face_count,
            covered_faces,
            mesh.faces.len(),
        );
        if config.verbose {
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
        if !config.quiet {
            eprintln!("  Compare: all planar hypothesis centroids within tolerance");
        }
    }

    if planar_quit {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }



    if !config.stage.at_least(2, 2) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses: Vec::new(),
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.2: Deduce cylindrical hypotheses
    let viz_22 = if config.viz_active(2, 2) { viz } else { None };
    let (mut cylindrical_hypotheses, cylindrical_quit) = deduce_cylindrical_hypotheses(
        &mut mesh, config.vertex_tolerance_mm,
        config.surface_tolerance_mm, config.angular_tolerance_rad,
        config.verbosity,
        viz_22,
    );

    if !config.quiet {
        let covered_faces: usize = cylindrical_hypotheses.iter().map(|h| h.faces.len()).sum();
        let convex_count = cylindrical_hypotheses.iter().filter(|h| h.convex).count();
        let concave_count = cylindrical_hypotheses.len() - convex_count;
        eprintln!(
            "Stage 2.2: Deduced {} cylindrical hypotheses ({} convex, {} concave) covering {} faces",
            cylindrical_hypotheses.len(),
            convex_count,
            concave_count,
            covered_faces,
        );
        if config.verbose {
            for (i, h) in cylindrical_hypotheses.iter().enumerate() {
                eprintln!(
                    "  Cylinder {}: {} faces, {} vertices, r={:.4}, {}, \
axis_origin=[{:.4}, {:.4}, {:.4}], axis_dir=[{:.4}, {:.4}, {:.4}], \
vtx_err_max={:.2e}, cen_err_max={:.2e}",
                    i, h.faces.len(), h.vertices.len(), h.radius,
                    if h.convex { "convex" } else { "concave" },
                    h.axis_origin[0], h.axis_origin[1], h.axis_origin[2],
                    h.axis_direction[0], h.axis_direction[1], h.axis_direction[2],
                    h.error_max, h.centroid_error_max,
                );
                if config.verbosity >= 2 {
                    for &fi in &h.faces {
                        let face = &mesh.faces[fi];
                        let vc = face.vertex_count as usize;
                        let centroid = face_centroid(face, &mesh.vertices);
                        let normal = face.normal.unwrap_or([0.0; 3]);
                        let centroid_err = vertex_to_cylinder_distance(
                            &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                            &h.axis_origin, &h.axis_direction, h.radius,
                        ).abs();
                        let mut vtx_err_max = f64::NEG_INFINITY;
                        let mut vtx_err_min = f64::INFINITY;
                        for vi_idx in 0..vc {
                            let d = vertex_to_cylinder_distance(
                                &mesh.vertices[face.vertex_indices[vi_idx]],
                                &h.axis_origin, &h.axis_direction, h.radius,
                            ).abs();
                            vtx_err_max = vtx_err_max.max(d);
                            vtx_err_min = vtx_err_min.min(d);
                        }
                        eprintln!(
                            "    face {fi}: {vc}v centroid=[{:.4},{:.4},{:.4}] \
Normal=[{:.3},{:.3},{:.3}] vtx_err=[{:.2e},{:.2e}] cen_err={:.2e}",
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
        compare_cylindrical_hypotheses(&cylindrical_hypotheses, &mesh, config)?;
        if !config.quiet {
            eprintln!("  Compare: all cylindrical hypothesis centroids within tolerance");
        }
    }

    if cylindrical_quit {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses,
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }



    if !config.stage.at_least(2, 3) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses,
            spherical_hypotheses: Vec::new(),
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.3: Deduce spherical hypotheses
    let bb_diag = bounding_box_diagonal(&mesh.vertices);
    let max_sphere_radius = bb_diag * MAX_SPHERE_RADIUS_FACTOR;
    let viz_23 = if config.viz_active(2, 3) { viz } else { None };
    let (mut spherical_hypotheses, spherical_quit) = deduce_spherical_hypotheses(
        &mut mesh, &planar_hypotheses, config.vertex_tolerance_mm,
        config.surface_tolerance_mm, config.angular_tolerance_rad, max_sphere_radius,
        config.verbosity,
        viz_23,
    );

    if !config.quiet {
        let covered_faces: usize = spherical_hypotheses.iter().map(|h| h.faces.len()).sum();
        let convex_count = spherical_hypotheses.iter().filter(|h| h.convex).count();
        let concave_count = spherical_hypotheses.len() - convex_count;
        eprintln!(
            "Stage 2.3: Deduced {} spherical hypotheses ({} convex, {} concave) covering {} faces",
            spherical_hypotheses.len(),
            convex_count,
            concave_count,
            covered_faces,
        );
        if config.verbose {
            for (i, h) in spherical_hypotheses.iter().enumerate() {
                eprintln!(
                    "  Sphere {}: {} faces, {} vertices, r={:.4}, {}, center=[{:.4}, {:.4}, {:.4}], \
vtx_err_max={:.2e}, cen_err_max={:.2e}",
                    i, h.faces.len(), h.vertices.len(), h.radius,
                    if h.convex { "convex" } else { "concave" },
                    h.center[0], h.center[1], h.center[2],
                    h.error_max, h.centroid_error_max,
                );
                if config.verbosity >= 2 {
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
        if !config.quiet {
            eprintln!("  Compare: all spherical hypothesis centroids within tolerance");
        }
    }

    if spherical_quit {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses,
            spherical_hypotheses,
            selected_surfaces: Vec::new(),
        });
    }



    if !config.stage.at_least(2, 4) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses,
            spherical_hypotheses,
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.4: Deduce ruled surface hypotheses
    // TODO: optional - detect extruded curve surfaces
    if !config.quiet {
        eprintln!("Stage 2.4: Deduce ruled surface hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 5) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses,
            spherical_hypotheses,
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.5: Deduce NURBS hypotheses
    // TODO: fit NURBS to remaining ungrouped faces
    if !config.quiet {
        eprintln!("Stage 2.5: Deduce NURBS hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 6) {
        return Ok(Stage2Output {
            mesh,
            planar_hypotheses,
            cylindrical_hypotheses,
            spherical_hypotheses,
            selected_surfaces: Vec::new(),
        });
    }

    // Stage 2.6: Select surfaces for reconstruction
    let viz_26 = if config.viz_active(2, 6) { viz } else { None };
    let selected_surfaces = select_surfaces(
        &mesh, &mut planar_hypotheses, &mut cylindrical_hypotheses,
        &mut spherical_hypotheses, viz_26,
    );

    if !config.quiet {
        let mut planar_count = 0;
        let mut cylindrical_count = 0;
        let mut spherical_count = 0;
        for s in &selected_surfaces {
            match s {
                SelectedSurface::Planar(_) => planar_count += 1,
                SelectedSurface::Cylindrical(_) => cylindrical_count += 1,
                SelectedSurface::Spherical(_) => spherical_count += 1,
            }
        }
        eprintln!(
            "Stage 2.6: Selected {} surfaces ({} planar, {} cylindrical, {} spherical) covering {} faces",
            selected_surfaces.len(),
            planar_count,
            cylindrical_count,
            spherical_count,
            mesh.faces.len(),
        );
    }

    // Compare against STEP file if --compare was specified
    if config.compare_shape.is_some() {
        compare_selected_surfaces(
            &selected_surfaces, &planar_hypotheses, &cylindrical_hypotheses,
            &spherical_hypotheses, &mesh, config,
        )?;
        if !config.quiet {
            eprintln!("  Compare: all selected surface centroids within tolerance");
        }
    }

    Ok(Stage2Output {
        mesh,
        planar_hypotheses,
        cylindrical_hypotheses,
        spherical_hypotheses,
        selected_surfaces,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stage1::{self, MeshFace, MeshVertex, VertexWeldOptions, UNDEDUCED_PLANAR_HYPOTHESIS, UNDEDUCED_CYLINDRICAL_HYPOTHESIS, UNDEDUCED_SPHERICAL_HYPOTHESIS};

    fn make_triangle_face(a: usize, b: usize, c: usize) -> MeshFace {
        MeshFace {
            vertex_count: 3,
            vertex_indices: [a, b, c, 0],
            neighbors: [-1, -1, -1, -1],
            normal: None,
            planar_hypothesis: UNDEDUCED_PLANAR_HYPOTHESIS,
            cylindrical_hypothesis: UNDEDUCED_CYLINDRICAL_HYPOTHESIS,
            spherical_hypothesis: UNDEDUCED_SPHERICAL_HYPOTHESIS,
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

    #[test]
    fn cube_gets_six_planar_hypotheses() {
        let stl_path = format!("{}/tests/manual/cube.stl", env!("CARGO_MANIFEST_DIR"));
        let mut mesh =
            stage1::read_connected_mesh_from_stl(&stl_path, VertexWeldOptions { tolerance: 1e-9 })
                .expect("should load");
        mesh.validate_and_populate_topology().expect("should validate");

        let (hypotheses, _) = deduce_planar_hypotheses(&mut mesh, 1e-5, 0, None);

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
        // Two triangles sharing an edge but on perpendicular planes
        let verts = vec![
            MeshVertex::from_xyz(0.0, 0.0, 0.0),
            MeshVertex::from_xyz(1.0, 0.0, 0.0),
            MeshVertex::from_xyz(0.0, 1.0, 0.0), // triangle 0: on XY plane
            MeshVertex::from_xyz(0.0, 0.0, 1.0), // triangle 1: on XZ plane
        ];
        let mesh = build_mesh(verts, vec![[0, 1, 2], [0, 3, 1]]);
        let mut mesh = mesh;
        let (hypotheses, _) = deduce_planar_hypotheses(&mut mesh, 1e-5, 0, None);

        assert_eq!(
            hypotheses.len(),
            2,
            "perpendicular triangles should get 2 hypotheses"
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
        let mesh = build_mesh(verts, vec![[0, 1, 2], [0, 2, 3]]);
        let mut mesh = mesh;
        let (hypotheses, _) = deduce_planar_hypotheses(&mut mesh, 1e-5, 0, None);

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
        let (_planar, _) = deduce_planar_hypotheses(&mut mesh, 1e-5, 0, None);
        // 17.5° angular tolerance: cube faces meet at 90°, so no cylinders
        let (cyls, _) = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 17.5_f64.to_radians(), 0, None
        );
        assert_eq!(cyls.len(), 0, "cube should produce 0 cylinders at 17.5° angular tolerance");
    }

    #[test]
    fn cube_angular_tolerance_rejects_sphere() {
        let mut mesh = load_stage1("manual/cube.stl");
        let (planar, _) = deduce_planar_hypotheses(&mut mesh, 1e-5, 0, None);
        let _ = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 17.5_f64.to_radians(), 0, None
        );
        // 17.5° angular tolerance: cube faces meet at 90°, so no spheres
        let (sphs, _) = deduce_spherical_hypotheses(
            &mut mesh, &planar, 1e-5, 0.4, 17.5_f64.to_radians(), 100.0, 0, None,
        );
        assert_eq!(sphs.len(), 0, "cube should produce 0 spheres at 17.5° angular tolerance");
    }

    #[test]
    fn cube_high_angular_tolerance_allows_sphere() {
        let mut mesh = load_stage1("manual/cube.stl");
        let (planar, _) = deduce_planar_hypotheses(&mut mesh, 1e-5, 0, None);
        let _ = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 91.0_f64.to_radians(), 0, None
        );
        // 91° angular tolerance: exceeds 90° dihedral angle of cube
        let (sphs, _) = deduce_spherical_hypotheses(
            &mut mesh, &planar, 1e-5, 0.4, 91.0_f64.to_radians(), 100.0, 0, None,
        );
        assert!(!sphs.is_empty(), "cube should produce spheres at 91° angular tolerance");
    }

    #[test]
    fn cylinder_detected_at_default_angular_tolerance() {
        let mut mesh = load_stage1("ccad/generated/simple_cylinder.stl");
        let (_planar, _) = deduce_planar_hypotheses(&mut mesh, 1e-5, 0, None);
        let (cyls, _) = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 17.5_f64.to_radians(), 0, None
        );
        assert_eq!(cyls.len(), 1, "simple cylinder should produce 1 cylinder hypothesis");
        assert!(cyls[0].convex);
        assert!((cyls[0].radius - 10.0).abs() < 0.01);
    }

    #[test]
    fn sphere_detected_at_default_angular_tolerance() {
        let mut mesh = load_stage1("ccad/generated/simple_sphere.stl");
        let (planar, _) = deduce_planar_hypotheses(&mut mesh, 1e-5, 0, None);
        let _ = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 17.5_f64.to_radians(), 0, None
        );
        let (sphs, _) = deduce_spherical_hypotheses(
            &mut mesh, &planar, 1e-5, 0.4, 17.5_f64.to_radians(), 1000.0, 0, None,
        );
        assert_eq!(sphs.len(), 1, "simple sphere should produce 1 sphere hypothesis");
        assert!(sphs[0].convex);
        assert!((sphs[0].radius - 10.0).abs() < 0.01);
    }
}
