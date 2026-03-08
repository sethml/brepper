//! Stage 2: Surface Fitting
//!
//! 2.1: Deduce planar hypotheses
//! 2.2: Deduce cylindrical hypotheses
//! 2.3: Deduce spherical hypotheses
//! 2.4: Deduce ruled surface hypotheses
//! 2.5: Deduce NURBS hypotheses
//! 2.6: Select surfaces for reconstruction

use crate::config::Config;
use crate::stage1::{self, ConnectedMesh, MeshFace, MeshVertex, UNDEDUCED_PLANAR_HYPOTHESIS, UNDEDUCED_CYLINDRICAL_HYPOTHESIS, UNDEDUCED_SPHERICAL_HYPOTHESIS, NO_HYPOTHESIS};
use opencascade_sys::gp;
use std::collections::{HashSet, VecDeque};
use std::error::Error;
use std::fmt::{Display, Formatter};


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
const MIN_CROSS_THRESHOLD: f64 = 0.01;

/// Fast-reject multiplier for BFS vertex-distance checks. When a candidate face
/// has a vertex farther than `REFIT_SKIP_MULTIPLIER × vertex_tol` from the
/// current fitted surface, the distance is too large for a re-fit to absorb and
/// the face is rejected immediately without attempting re-fitting.
const REFIT_SKIP_MULTIPLIER: f64 = 2.0;

/// Minimum number of mesh faces required to accept a cylindrical hypothesis.
/// Any real CAD tessellation of a cylinder produces at least 3 facets around
/// the circumference; this rejects spurious 2-face fits (e.g. adjacent torus
/// facets that locally approximate a cylinder).
const MIN_CYLINDER_FACES: usize = 3;

/// Minimum number of mesh faces required to accept a spherical hypothesis.
/// Sphere fitting has 4 degrees of freedom (cx, cy, cz, r), so at least 4
/// non-degenerate faces are needed for a well-determined fit. This also rejects
/// spurious fits from small patches that are consistent with many surface types.
const MIN_SPHERE_FACES: usize = 4;

/// Maximum sphere radius as a multiple of the mesh bounding-box diagonal.
/// With solid-angle coverage validation and surface-tolerance validation during
/// BFS, this no longer needs to be tight — those checks prevent pathological
/// growth. A large value serves as a numerical guardrail against degenerate fits.
const MAX_SPHERE_RADIUS_FACTOR: f64 = 1000.0;

/// Minimum eigenvalue ratio for solid-angle coverage validation.
/// Centroid-to-center direction vectors of a genuine spherical hypothesis span
/// 3D, so all eigenvalues of their covariance matrix are substantial.
/// Fillet-strip growth produces nearly coplanar directions (λ₃ ≈ 0).
const MIN_SPHERE_EIGENVALUE_RATIO: f64 = 0.01;

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
// Stage 2.1: Planar hypothesis deduction
// ---------------------------------------------------------------------------

/// Signed distance from a vertex to a plane defined by (normal, distance).
fn vertex_to_plane_distance(v: &MeshVertex, normal: &[f64; 3], distance: f64) -> f64 {
    normal[0] * v.x + normal[1] * v.y + normal[2] * v.z - distance
}

/// Compute the area of a triangular mesh face.
fn face_area(face: &MeshFace, vertices: &[MeshVertex]) -> f64 {
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
    0.5 * (cx * cx + cy * cy + cz * cz).sqrt()
}

/// Fit a plane to a set of faces using area-weighted normal averaging and
/// vertex centroid for the distance.
fn fit_plane(
    face_indices: &[usize],
    vertex_set: &HashSet<usize>,
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> ([f64; 3], f64) {
    let mut sum_nx = 0.0_f64;
    let mut sum_ny = 0.0_f64;
    let mut sum_nz = 0.0_f64;

    for &fi in face_indices {
        let face = &faces[fi];
        let n = face.normal.unwrap();
        let area = face_area(face, vertices);
        sum_nx += n[0] * area;
        sum_ny += n[1] * area;
        sum_nz += n[2] * area;
    }

    let len = (sum_nx * sum_nx + sum_ny * sum_ny + sum_nz * sum_nz).sqrt();
    let normal = [sum_nx / len, sum_ny / len, sum_nz / len];

    let mut sum_d = 0.0_f64;
    for &vi in vertex_set {
        let v = &vertices[vi];
        sum_d += normal[0] * v.x + normal[1] * v.y + normal[2] * v.z;
    }
    let distance = sum_d / vertex_set.len() as f64;

    (normal, distance)
}

/// Check that all vertices in a set are within tolerance of a plane.
fn all_vertices_within_tolerance(
    vertex_set: &HashSet<usize>,
    normal: &[f64; 3],
    distance: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        if vertex_to_plane_distance(&vertices[vi], normal, distance).abs() > tolerance {
            return false;
        }
    }
    true
}

/// Deduce planar hypotheses from the mesh using BFS region growing.
///
/// For each unassigned face, creates a new planar hypothesis seeded from that
/// face's plane, then grows it via BFS to neighboring faces that are coplanar
/// (normal alignment and vertex-to-plane distance within tolerance).
fn deduce_planar_hypotheses(
    mesh: &mut ConnectedMesh,
    vertex_tol: f64,
) -> Vec<PlanarHypothesis> {
    let num_faces = mesh.faces.len();
    let mut hypotheses: Vec<PlanarHypothesis> = Vec::new();

    for face in &mut mesh.faces {
        face.planar_hypothesis = UNDEDUCED_PLANAR_HYPOTHESIS;
    }

    for fi in 0..num_faces {
        if mesh.faces[fi].planar_hypothesis != UNDEDUCED_PLANAR_HYPOTHESIS {
            continue;
        }

        let hi = hypotheses.len() as i32;
        let seed_normal = mesh.faces[fi].normal.unwrap();

        // Compute initial plane from seed face: normal from face, distance
        // as average projection of seed face vertices.
        let vc = mesh.faces[fi].vertex_count as usize;
        let mut sum_d = 0.0_f64;
        let mut vertex_set = HashSet::new();
        for vi_idx in 0..vc {
            let vi = mesh.faces[fi].vertex_indices[vi_idx];
            vertex_set.insert(vi);
            let v = &mesh.vertices[vi];
            sum_d += seed_normal[0] * v.x + seed_normal[1] * v.y + seed_normal[2] * v.z;
        }
        let mut current_normal = seed_normal;
        let mut current_distance = sum_d / vc as f64;
        let mut face_list: Vec<usize> = vec![fi];

        mesh.faces[fi].planar_hypothesis = hi;

        // BFS expansion
        let mut queue = VecDeque::new();
        queue.push_back(fi);

        while let Some(current_fi) = queue.pop_front() {
            let vertex_count = mesh.faces[current_fi].vertex_count as usize;
            let neighbors = mesh.faces[current_fi].neighbors;

            for &ni in &neighbors[..vertex_count] {
                if ni < 0 {
                    continue;
                }
                let ni = ni as usize;

                if mesh.faces[ni].planar_hypothesis != UNDEDUCED_PLANAR_HYPOTHESIS {
                    continue;
                }
                // Vertex distance check
                let nvc = mesh.faces[ni].vertex_count as usize;
                let nvi = mesh.faces[ni].vertex_indices;
                let mut all_ok = true;
                let mut any_far = false;
                for &vi in &nvi[..nvc] {
                    let d = vertex_to_plane_distance(
                        &mesh.vertices[vi],
                        &current_normal,
                        current_distance,
                    );
                    let abs_d = d.abs();
                    if abs_d > vertex_tol {
                        all_ok = false;
                        if abs_d > REFIT_SKIP_MULTIPLIER * vertex_tol {
                            any_far = true;
                            break;
                        }
                    }
                }

                // Skip if any vertex is too far for re-fitting to help
                if any_far {
                    continue;
                }

                if !all_ok {
                    // Try re-fitting with current + new face's vertices
                    let mut trial_vertices = vertex_set.clone();
                    for &vi in &nvi[..nvc] {
                        trial_vertices.insert(vi);
                    }
                    let mut trial_faces = face_list.clone();
                    trial_faces.push(ni);

                    let (new_normal, new_distance) = fit_plane(
                        &trial_faces,
                        &trial_vertices,
                        &mesh.faces,
                        &mesh.vertices,
                    );

                    // Check ALL vertices (existing + new) against the re-fitted plane
                    if !all_vertices_within_tolerance(
                        &trial_vertices,
                        &new_normal,
                        new_distance,
                        vertex_tol,
                        &mesh.vertices,
                    ) {
                        continue;
                    }

                    // Accept re-fit
                    current_normal = new_normal;
                    current_distance = new_distance;
                }

                // Accept this face into the hypothesis
                mesh.faces[ni].planar_hypothesis = hi;
                face_list.push(ni);
                for &vi in &nvi[..nvc] {
                    vertex_set.insert(vi);
                }
                queue.push_back(ni);
            }
        }

        // Final re-fit using all collected faces and vertices
        let (final_normal, final_distance) =
            fit_plane(&face_list, &vertex_set, &mesh.faces, &mesh.vertices);

        // Compute error metrics
        let mut error_max = f64::NEG_INFINITY;
        let mut error_min = f64::INFINITY;
        let mut error_abs_sum = 0.0_f64;
        for &vi in &vertex_set {
            let d = vertex_to_plane_distance(&mesh.vertices[vi], &final_normal, final_distance);
            error_max = error_max.max(d);
            error_min = error_min.min(d);
            error_abs_sum += d.abs();
        }

        hypotheses.push(PlanarHypothesis {
            normal: final_normal,
            distance: final_distance,
            faces: face_list,
            vertices: vertex_set.into_iter().collect(),
            error_max,
            error_min,
            error_abs_sum,
        });
    }

    hypotheses
}

// ---------------------------------------------------------------------------
// Stage 2.1: Compare check
// ---------------------------------------------------------------------------

/// Validate fitted planar hypotheses against a reference STEP shape.
///
/// For each hypothesis, projects face centroids onto the fitted plane and
/// checks that those projected points are within surface_tolerance of the
/// reference STEP surface.
fn compare_planar_hypotheses(
    hypotheses: &[PlanarHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let face = &mesh.faces[fi];
            let n = face.vertex_count as usize;
            let mut cx = 0.0_f64;
            let mut cy = 0.0_f64;
            let mut cz = 0.0_f64;
            for vi_idx in 0..n {
                let v = &mesh.vertices[face.vertex_indices[vi_idx]];
                cx += v.x;
                cy += v.y;
                cz += v.z;
            }
            let inv_n = 1.0 / n as f64;
            cx *= inv_n;
            cy *= inv_n;
            cz *= inv_n;

            // Project centroid onto the fitted plane
            let dist_to_plane = hyp.normal[0] * cx + hyp.normal[1] * cy + hyp.normal[2] * cz
                - hyp.distance;
            let px = cx - dist_to_plane * hyp.normal[0];
            let py = cy - dist_to_plane * hyp.normal[1];
            let pz = cz - dist_to_plane * hyp.normal[2];

            let pt = gp::Pnt::new_real3(px, py, pz);
            let d = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(d);
        }

        if max_dist > config.surface_tolerance_mm {
            return Err(Stage2CompareError {
                hypothesis_type: "planar",
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance: config.surface_tolerance_mm,
            });
        }
    }

    Ok(())
}
// ---------------------------------------------------------------------------
// Stage 2.2: Cylindrical hypothesis deduction
// ---------------------------------------------------------------------------

/// Compute the face centroid for a mesh face.
fn face_centroid(face: &MeshFace, vertices: &[MeshVertex]) -> [f64; 3] {
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

/// Compute the distance from a vertex to a cylinder surface.
/// Returns the signed distance: positive if outside the cylinder, negative if inside.
fn vertex_to_cylinder_distance(
    v: &MeshVertex,
    axis_origin: &[f64; 3],
    axis_direction: &[f64; 3],
    radius: f64,
) -> f64 {
    // Vector from axis origin to vertex
    let dx = v.x - axis_origin[0];
    let dy = v.y - axis_origin[1];
    let dz = v.z - axis_origin[2];
    // Project onto axis
    let t = dx * axis_direction[0] + dy * axis_direction[1] + dz * axis_direction[2];
    // Radial vector (perpendicular to axis)
    let rx = dx - t * axis_direction[0];
    let ry = dy - t * axis_direction[1];
    let rz = dz - t * axis_direction[2];
    let radial_dist = (rx * rx + ry * ry + rz * rz).sqrt();
    radial_dist - radius
}

/// Normalize a 3D vector in place, return its length.
fn normalize3(v: &mut [f64; 3]) -> f64 {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 0.0 {
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
    }
    len
}

/// Cross product of two 3D vectors.
fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two 3D vectors.
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the smallest eigenvector of a 3x3 symmetric matrix.
/// The matrix is stored as [m00, m01, m02, m11, m12, m22] (upper triangle).
/// Returns the eigenvector corresponding to the smallest eigenvalue.
fn smallest_eigenvector_3x3(m: &[f64; 6]) -> [f64; 3] {
    // m = [[m[0], m[1], m[2]],
    //      [m[1], m[3], m[4]],
    //      [m[2], m[4], m[5]]]
    // Use the characteristic polynomial approach for 3x3 symmetric matrices.
    // Characteristic equation: λ³ - c2*λ² + c1*λ - c0 = 0
    // where c2 = trace, c1 = sum of 2x2 minors, c0 = determinant
    let a00 = m[0]; let a01 = m[1]; let a02 = m[2];
    let a11 = m[3]; let a12 = m[4]; let a22 = m[5];

    let c2 = a00 + a11 + a22;
    let c1 = a00 * a11 - a01 * a01 + a00 * a22 - a02 * a02 + a11 * a22 - a12 * a12;
    let c0 = a00 * (a11 * a22 - a12 * a12)
            - a01 * (a01 * a22 - a12 * a02)
            + a02 * (a01 * a12 - a11 * a02);

    // Depressed cubic: t³ + pt + q = 0 where λ = t + c2/3
    let p = c1 - c2 * c2 / 3.0;
    let q = -2.0 * c2 * c2 * c2 / 27.0 + c1 * c2 / 3.0 - c0;

    let eigenvalues = if p.abs() < 1e-30 {
        // All eigenvalues equal
        let ev = c2 / 3.0;
        [ev, ev, ev]
    } else {
        // p < 0 for real symmetric matrices with distinct eigenvalues
        let neg_p_3 = (-p / 3.0).max(0.0);
        let r = neg_p_3.sqrt();
        let cos_arg = (-q / (2.0 * neg_p_3 * r)).clamp(-1.0, 1.0);
        let theta = cos_arg.acos() / 3.0;
        let shift = c2 / 3.0;
        let two_pi_3 = 2.0 * std::f64::consts::PI / 3.0;
        let ev0 = 2.0 * r * theta.cos() + shift;
        let ev1 = 2.0 * r * (theta - two_pi_3).cos() + shift;
        let ev2 = 2.0 * r * (theta - 2.0 * two_pi_3).cos() + shift;
        [ev0, ev1, ev2]
    };

    // Find the smallest eigenvalue
    let (min_idx, _) = eigenvalues.iter().enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();

    // Compute eigenvector for the smallest eigenvalue using inverse iteration:
    // (M - λI)v = 0, find null space
    // Shift slightly to avoid exact singularity issues
    let lambda = eigenvalues[min_idx];
    let b00 = a00 - lambda;
    let b11 = a11 - lambda;
    let b22 = a22 - lambda;

    // Try cross products of rows to find the null space
    let row0 = [b00, a01, a02];
    let row1 = [a01, b11, a12];
    let row2 = [a02, a12, b22];

    let candidates = [
        cross3(&row0, &row1),
        cross3(&row0, &row2),
        cross3(&row1, &row2),
    ];

    // Pick the cross product with largest magnitude
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
        // Fallback: eigenvalues are degenerate, return arbitrary unit vector
        return [1.0, 0.0, 0.0];
    }
    ev
}

/// Build the area-weighted normal covariance matrix M = Σ wᵢ nᵢ nᵢᵀ.
/// Returns the upper triangle [m00, m01, m02, m11, m12, m22].
fn build_normal_covariance(
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

/// Fit cylinder parameters (axis direction, axis origin, radius) from a set of faces.
/// Returns (axis_origin, axis_direction, radius).
fn fit_cylinder(
    face_indices: &[usize],
    vertex_set: &HashSet<usize>,
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> Option<([f64; 3], [f64; 3], f64)> {
    // Step 1: Axis direction from smallest eigenvector of normal covariance.
    let cov = build_normal_covariance(face_indices, faces, vertices);
    let mut axis_dir = smallest_eigenvector_3x3(&cov);
    let len = normalize3(&mut axis_dir);
    if len < 1e-15 {
        return None;
    }

    // Step 2: Compute orthogonal basis perpendicular to axis.
    let (u, w) = perpendicular_basis(&axis_dir);

    // Step 3: Project vertices onto 2D plane perpendicular to axis and fit circle.
    let verts: Vec<usize> = vertex_set.iter().copied().collect();
    let n = verts.len();
    if n < 3 {
        return None;
    }

    // Collect 2D coordinates
    let mut xs = Vec::with_capacity(n);
    let mut ys = Vec::with_capacity(n);
    for &vi in &verts {
        let v = &vertices[vi];
        let vx = v.x;
        let vy = v.y;
        let vz = v.z;
        xs.push(vx * u[0] + vy * u[1] + vz * u[2]);
        ys.push(vx * w[0] + vy * w[1] + vz * w[2]);
    }

    // Center the 2D data to improve numerical conditioning.
    // Without centering, the normal equations matrix becomes ill-conditioned
    // when absolute coordinates are large relative to the arc span (e.g.,
    // 6 vertices spanning 2.5° of a r=2mm cylinder at position (15, 12.5)).
    let nf = n as f64;
    let mean_x = xs.iter().sum::<f64>() / nf;
    let mean_y = ys.iter().sum::<f64>() / nf;
    for i in 0..n {
        xs[i] -= mean_x;
        ys[i] -= mean_y;
    }

    // Algebraic circle fit: x² + y² + Dx + Ey + F = 0
    // Least squares: minimize Σ(x² + y² + Dx + Ey + F)²
    // Normal equations: [Σx²  Σxy  Σx ] [D]   [-Σx(x²+y²)]
    //                   [Σxy  Σy²  Σy ] [E] = [-Σy(x²+y²)]
    //                   [Σx   Σy   n  ] [F]   [-Σ(x²+y²) ]
    let mut sx = 0.0_f64;
    let mut sy = 0.0_f64;
    let mut sx2 = 0.0_f64;
    let mut sy2 = 0.0_f64;
    let mut sxy = 0.0_f64;
    let mut sx3 = 0.0_f64;
    let mut sy3 = 0.0_f64;
    let mut sx2y = 0.0_f64;
    let mut sxy2 = 0.0_f64;

    for i in 0..n {
        let x = xs[i];
        let y = ys[i];
        sx += x;
        sy += y;
        sx2 += x * x;
        sy2 += y * y;
        sxy += x * y;
        sx3 += x * x * x;
        sy3 += y * y * y;
        sx2y += x * x * y;
        sxy2 += x * y * y;
    }

    // RHS
    let rhs0 = -(sx3 + sxy2);
    let rhs1 = -(sx2y + sy3);
    let rhs2 = -(sx2 + sy2);

    // 3x3 system: solve using Cramer's rule
    let a = [[sx2, sxy, sx], [sxy, sy2, sy], [sx, sy, nf]];

    let det_a = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
              - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
              + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);

    if det_a.abs() < 1e-30 {
        return None; // Degenerate — vertices are collinear in projection
    }

    let inv_det = 1.0 / det_a;

    let d_val = inv_det * (rhs0 * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                         - a[0][1] * (rhs1 * a[2][2] - a[1][2] * rhs2)
                         + a[0][2] * (rhs1 * a[2][1] - a[1][1] * rhs2));

    let e_val = inv_det * (a[0][0] * (rhs1 * a[2][2] - a[1][2] * rhs2)
                         - rhs0 * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                         + a[0][2] * (a[1][0] * rhs2 - rhs1 * a[2][0]));

    let f_val = inv_det * (a[0][0] * (a[1][1] * rhs2 - rhs1 * a[2][1])
                         - a[0][1] * (a[1][0] * rhs2 - rhs1 * a[2][0])
                         + rhs0 * (a[1][0] * a[2][1] - a[1][1] * a[2][0]));

    let center_x = -d_val / 2.0;
    let center_y = -e_val / 2.0;
    let r_sq = center_x * center_x + center_y * center_y - f_val;
    if r_sq <= 0.0 {
        return None;
    }
    let radius = r_sq.sqrt();

    // Convert 2D center back to 3D (undo centering offset)
    let abs_center_x = center_x + mean_x;
    let abs_center_y = center_y + mean_y;
    let axis_origin = [
        abs_center_x * u[0] + abs_center_y * w[0],
        abs_center_x * u[1] + abs_center_y * w[1],
        abs_center_x * u[2] + abs_center_y * w[2],
    ];

    Some((axis_origin, axis_dir, radius))
}

/// Check if all vertices in a set are within tolerance of a cylinder surface.
fn all_vertices_within_cylinder_tolerance(
    vertex_set: &HashSet<usize>,
    axis_origin: &[f64; 3],
    axis_direction: &[f64; 3],
    radius: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        let d = vertex_to_cylinder_distance(&vertices[vi], axis_origin, axis_direction, radius);
        if d.abs() > tolerance {
            return false;
        }
    }
    true
}

/// Determine convexity: does the face normal point away from the axis (convex) or toward it (concave)?
fn determine_convexity(
    face: &MeshFace,
    vertices: &[MeshVertex],
    axis_origin: &[f64; 3],
    axis_direction: &[f64; 3],
) -> bool {
    let centroid = face_centroid(face, vertices);
    // Vector from axis to centroid (radial component)
    let d = [
        centroid[0] - axis_origin[0],
        centroid[1] - axis_origin[1],
        centroid[2] - axis_origin[2],
    ];
    let t = dot3(&d, axis_direction);
    let radial = [
        d[0] - t * axis_direction[0],
        d[1] - t * axis_direction[1],
        d[2] - t * axis_direction[2],
    ];
    let n = face.normal.unwrap();
    dot3(&n, &radial) > 0.0
}

/// Build an orthonormal basis (u, w) perpendicular to a unit direction `d`.
fn perpendicular_basis(d: &[f64; 3]) -> ([f64; 3], [f64; 3]) {
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

/// Validate that faces in a cylindrical hypothesis have sufficient angular
/// coverage around the cylinder — i.e., they are distributed around the
/// circumference, not clustered on one side.
///
/// Requires at least 3 distinct angular clusters. The check is:
/// - Compute the angular coordinate θ of each face centroid around the axis.
/// - Sort θ values and compute all N inter-face gaps plus the wraparound gap.
/// - The largest gap defines the "empty arc"; span = 2π − largest_gap.
/// - The second-largest gap must be ≤ span / 3.
fn angular_coverage_valid(
    face_list: &[usize],
    axis_origin: &[f64; 3],
    axis_direction: &[f64; 3],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> bool {
    if face_list.len() < 3 {
        return false;
    }

    let (u, w) = perpendicular_basis(axis_direction);

    // Compute angular coordinate for each face centroid
    let mut thetas: Vec<f64> = Vec::with_capacity(face_list.len());
    for &fi in face_list {
        let c = face_centroid(&faces[fi], vertices);
        let d = [
            c[0] - axis_origin[0],
            c[1] - axis_origin[1],
            c[2] - axis_origin[2],
        ];
        let t = dot3(&d, axis_direction);
        let radial = [
            d[0] - t * axis_direction[0],
            d[1] - t * axis_direction[1],
            d[2] - t * axis_direction[2],
        ];
        let theta = f64::atan2(dot3(&radial, &w), dot3(&radial, &u));
        thetas.push(theta);
    }

    thetas.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = thetas.len();

    // Compute all gaps (including wraparound)
    let mut gaps: Vec<f64> = Vec::with_capacity(n);
    for i in 1..n {
        gaps.push(thetas[i] - thetas[i - 1]);
    }
    // Wraparound gap
    gaps.push(thetas[0] + 2.0 * std::f64::consts::PI - thetas[n - 1]);

    // Find largest and second-largest gaps
    let mut largest = 0.0_f64;
    let mut second_largest = 0.0_f64;
    for &g in &gaps {
        if g > largest {
            second_largest = largest;
            largest = g;
        } else if g > second_largest {
            second_largest = g;
        }
    }

    let span = 2.0 * std::f64::consts::PI - largest;
    if span <= 0.0 {
        return false;
    }

    second_largest <= span / 3.0
}

/// Result of a trial BFS for cylindrical hypothesis evaluation.
struct CylinderTrialResult {
    faces: Vec<usize>,
    vertices: HashSet<usize>,
    axis_origin: [f64; 3],
    axis_direction: [f64; 3],
    radius: f64,
    convex: bool,
    error_max: f64,
    error_abs_sum: f64,
}

/// Run a trial BFS for a cylindrical hypothesis starting from a seed pair,
/// using temporary data structures (not mutating mesh face assignments).
/// Returns `None` if the trial fails validation (min faces, centroid check).
#[allow(clippy::too_many_arguments)]
fn run_cylinder_trial_bfs(
    seed_fi: usize,
    seed_ni: usize,
    mesh: &ConnectedMesh,
    vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
    best_candidate: Option<&CylinderTrialResult>,
) -> Option<CylinderTrialResult> {
    let fi_vc = mesh.faces[seed_fi].vertex_count as usize;
    let ni_vc = mesh.faces[seed_ni].vertex_count as usize;

    // Collect seed vertices
    let mut vertex_set: HashSet<usize> = mesh.faces[seed_fi].vertex_indices[..fi_vc]
        .iter().copied().collect();
    for &vi in &mesh.faces[seed_ni].vertex_indices[..ni_vc] {
        vertex_set.insert(vi);
    }

    let mut face_list = vec![seed_fi, seed_ni];

    // Fit cylinder to seed pair
    let (mut current_origin, mut current_dir, mut current_radius) =
        fit_cylinder(&face_list, &vertex_set, &mesh.faces, &mesh.vertices)?;

    // Verify seed: all vertices within tolerance
    if !all_vertices_within_cylinder_tolerance(
        &vertex_set, &current_origin, &current_dir, current_radius,
        vertex_tol, &mesh.vertices,
    ) {
        return None;
    }

    // Determine convexity from seed face
    let convex = determine_convexity(
        &mesh.faces[seed_fi], &mesh.vertices, &current_origin, &current_dir,
    );

    // Track claimed faces in this trial (using a HashSet, not mesh mutation)
    let mut trial_claimed: HashSet<usize> = HashSet::new();
    trial_claimed.insert(seed_fi);
    trial_claimed.insert(seed_ni);

    let mut queue = VecDeque::new();
    queue.push_back(seed_fi);
    queue.push_back(seed_ni);

    // BFS expansion
    while let Some(current_fi) = queue.pop_front() {
        let vc = mesh.faces[current_fi].vertex_count as usize;
        let neighbors = mesh.faces[current_fi].neighbors;

        for &cni in &neighbors[..vc] {
            if cni < 0 {
                continue;
            }
            let cni = cni as usize;

            // Skip if already committed to a real hypothesis or claimed by this trial
            if mesh.faces[cni].cylindrical_hypothesis != UNDEDUCED_CYLINDRICAL_HYPOTHESIS {
                continue;
            }
            if trial_claimed.contains(&cni) {
                continue;
            }

            // Convexity check
            let cni_convex = determine_convexity(
                &mesh.faces[cni], &mesh.vertices, &current_origin, &current_dir,
            );
            if cni_convex != convex {
                continue;
            }

            // Angular tolerance: reject if dihedral angle with ANY already-claimed
            // neighbor exceeds the limit
            if let Some(n_cni) = mesh.faces[cni].normal {
                let cni_vc2 = mesh.faces[cni].vertex_count as usize;
                let cni_neighbors = mesh.faces[cni].neighbors;
                let mut angular_reject = false;
                for &adj in &cni_neighbors[..cni_vc2] {
                    if adj < 0 { continue; }
                    let adj = adj as usize;
                    // Check against both committed (same hypothesis doesn't exist yet)
                    // and trial-claimed faces
                    if !trial_claimed.contains(&adj) { continue; }
                    if let Some(n_adj) = mesh.faces[adj].normal {
                        let cos_a = dot3(&n_adj, &n_cni).clamp(-1.0, 1.0);
                        if cos_a.acos() > angular_tol {
                            angular_reject = true;
                            break;
                        }
                    }
                }
                if angular_reject { continue; }
            }

            // Vertex distance check
            let cni_vc = mesh.faces[cni].vertex_count as usize;
            let cni_vi = mesh.faces[cni].vertex_indices;
            let mut all_ok = true;
            let mut any_far = false;
            for &vi in &cni_vi[..cni_vc] {
                let d = vertex_to_cylinder_distance(
                    &mesh.vertices[vi],
                    &current_origin,
                    &current_dir,
                    current_radius,
                ).abs();
                if d > vertex_tol {
                    all_ok = false;
                    if d > REFIT_SKIP_MULTIPLIER * vertex_tol {
                        any_far = true;
                        break;
                    }
                }
            }

            if any_far {
                continue;
            }

            // Centroid validation: check face centroid distance to cylinder
            let centroid = face_centroid(&mesh.faces[cni], &mesh.vertices);
            let centroid_dist = vertex_to_cylinder_distance(
                &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                &current_origin, &current_dir, current_radius,
            ).abs();
            if centroid_dist > REFIT_SKIP_MULTIPLIER * surface_tol {
                continue;
            }
            let needs_refit = !all_ok || centroid_dist > surface_tol;

            if needs_refit {
                // Try re-fitting with new face included
                let mut trial_vertices = vertex_set.clone();
                for &vi in &cni_vi[..cni_vc] {
                    trial_vertices.insert(vi);
                }
                let mut trial_faces = face_list.clone();
                trial_faces.push(cni);

                let refit = fit_cylinder(
                    &trial_faces, &trial_vertices, &mesh.faces, &mesh.vertices,
                );
                let (new_origin, new_dir, new_radius) = match refit {
                    Some(params) => params,
                    None => continue,
                };

                if !all_vertices_within_cylinder_tolerance(
                    &trial_vertices, &new_origin, &new_dir, new_radius,
                    vertex_tol, &mesh.vertices,
                ) {
                    continue;
                }

                // Check all face centroids within surface tolerance after re-fit
                let mut centroids_ok = true;
                for &f in &trial_faces {
                    let c = face_centroid(&mesh.faces[f], &mesh.vertices);
                    let d = vertex_to_cylinder_distance(
                        &MeshVertex::from_xyz(c[0], c[1], c[2]),
                        &new_origin, &new_dir, new_radius,
                    ).abs();
                    if d > surface_tol {
                        centroids_ok = false;
                        break;
                    }
                }
                if !centroids_ok {
                    continue;
                }

                // Accept re-fit
                current_origin = new_origin;
                current_dir = new_dir;
                current_radius = new_radius;
            }

            // Accept this face
            trial_claimed.insert(cni);
            face_list.push(cni);
            for &vi in &cni_vi[..cni_vc] {
                vertex_set.insert(vi);
            }
            queue.push_back(cni);
        }

        // Early termination: if we have 3+ faces and are rediscovering the
        // same cylinder as the current best candidate, abandon.
        if face_list.len() >= 3 {
            if let Some(best) = best_candidate {
                let axis_dot = dot3(&current_dir, &best.axis_direction).abs();
                if axis_dot > 1.0 - 1e-6
                    && (current_radius - best.radius).abs() < vertex_tol
                {
                    // Check axis-to-axis distance
                    let d = [
                        current_origin[0] - best.axis_origin[0],
                        current_origin[1] - best.axis_origin[1],
                        current_origin[2] - best.axis_origin[2],
                    ];
                    let t = dot3(&d, &current_dir);
                    let perp_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - t * t;
                    let perp_dist = if perp_sq > 0.0 { perp_sq.sqrt() } else { 0.0 };

                    if perp_dist < vertex_tol
                        && face_list.iter().all(|f| best.faces.contains(f))
                    {
                        return None; // Redundant trial
                    }
                }
            }
        }
    }

    // Final re-fit
    if let Some((final_origin, final_dir, final_radius)) =
        fit_cylinder(&face_list, &vertex_set, &mesh.faces, &mesh.vertices)
    {
        current_origin = final_origin;
        current_dir = final_dir;
        current_radius = final_radius;
    }

    // Compute error metrics
    let mut error_max = 0.0_f64;
    let mut error_abs_sum = 0.0_f64;
    for &vi in &vertex_set {
        let d = vertex_to_cylinder_distance(
            &mesh.vertices[vi], &current_origin, &current_dir, current_radius,
        ).abs();
        error_max = error_max.max(d);
        error_abs_sum += d;
    }

    // Validate: minimum face count
    if face_list.len() < MIN_CYLINDER_FACES {
        return None;
    }

    // Validate: centroid check
    for &f in &face_list {
        let c = face_centroid(&mesh.faces[f], &mesh.vertices);
        let d = vertex_to_cylinder_distance(
            &MeshVertex::from_xyz(c[0], c[1], c[2]),
            &current_origin, &current_dir, current_radius,
        ).abs();
        if d > surface_tol {
            return None;
        }
    }

    Some(CylinderTrialResult {
        faces: face_list,
        vertices: vertex_set,
        axis_origin: current_origin,
        axis_direction: current_dir,
        radius: current_radius,
        convex,
        error_max,
        error_abs_sum,
    })
}

/// Deduce cylindrical hypotheses from the mesh using multi-seed BFS evaluation.
///
/// For each face, collects all valid seed partners, runs a trial BFS for each,
/// keeps the best trial (most faces), validates with angular coverage, and commits.
fn deduce_cylindrical_hypotheses(
    mesh: &mut ConnectedMesh,
    vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
) -> Vec<CylindricalHypothesis> {
    let num_faces = mesh.faces.len();
    let mut hypotheses: Vec<CylindricalHypothesis> = Vec::new();

    for fi in 0..num_faces {
        if mesh.faces[fi].cylindrical_hypothesis != UNDEDUCED_CYLINDRICAL_HYPOTHESIS {
            continue;
        }


        let fi_normal = match mesh.faces[fi].normal {
            Some(n) => n,
            None => {
                mesh.faces[fi].cylindrical_hypothesis = NO_HYPOTHESIS;
                continue;
            }
        };

        // Collect ALL valid seed partners
        let fi_vc = mesh.faces[fi].vertex_count as usize;
        let fi_neighbors = mesh.faces[fi].neighbors;
        let mut seed_partners: Vec<usize> = Vec::new();

        for &ni in &fi_neighbors[..fi_vc] {
            if ni < 0 {
                continue;
            }
            let ni = ni as usize;

            if mesh.faces[ni].cylindrical_hypothesis != UNDEDUCED_CYLINDRICAL_HYPOTHESIS {
                continue;
            }


            let ni_normal = match mesh.faces[ni].normal {
                Some(n) => n,
                None => continue,
            };

            // Check cross product magnitude for sufficient angular difference
            let cross = cross3(&fi_normal, &ni_normal);
            let cross_mag = (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
            if cross_mag < MIN_CROSS_THRESHOLD {
                continue;
            }

            // Angular tolerance: reject seed pairs with dihedral angle > limit
            let cos_angle = dot3(&fi_normal, &ni_normal).clamp(-1.0, 1.0);
            if cos_angle.acos() > angular_tol {
                continue;
            }

            seed_partners.push(ni);
        }

        if seed_partners.is_empty() {
            mesh.faces[fi].cylindrical_hypothesis = NO_HYPOTHESIS;
            continue;
        }

        // Multi-seed evaluation: try each seed partner, keep best
        let mut best_candidate: Option<CylinderTrialResult> = None;

        for &ni in &seed_partners {
            let trial = run_cylinder_trial_bfs(
                fi, ni, mesh,
                vertex_tol, surface_tol, angular_tol,
                best_candidate.as_ref(),
            );

            if let Some(trial) = trial {
                if best_candidate.as_ref().is_none_or(|b| trial.faces.len() > b.faces.len()) {
                    best_candidate = Some(trial);
                }
            }
        }

        if let Some(candidate) = best_candidate {
            // Angular coverage validation
            if !angular_coverage_valid(
                &candidate.faces, &candidate.axis_origin, &candidate.axis_direction,
                &mesh.faces, &mesh.vertices,
            ) {
                mesh.faces[fi].cylindrical_hypothesis = NO_HYPOTHESIS;
                continue;
            }

            // Commit: assign all faces to this hypothesis
            let hi = hypotheses.len() as i32;
            for &f in &candidate.faces {
                mesh.faces[f].cylindrical_hypothesis = hi;
            }

            hypotheses.push(CylindricalHypothesis {
                axis_origin: candidate.axis_origin,
                axis_direction: candidate.axis_direction,
                radius: candidate.radius,
                convex: candidate.convex,
                faces: candidate.faces,
                vertices: candidate.vertices.into_iter().collect(),
                error_max: candidate.error_max,
                error_abs_sum: candidate.error_abs_sum,
            });
        } else {
            mesh.faces[fi].cylindrical_hypothesis = NO_HYPOTHESIS;
        }
    }

    hypotheses
}

/// Validate fitted cylindrical hypotheses against a reference STEP shape.
fn compare_cylindrical_hypotheses(
    hypotheses: &[CylindricalHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);

            // Project centroid onto cylinder surface (nearest point)
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
            let radial_dist = (radial[0] * radial[0] + radial[1] * radial[1] + radial[2] * radial[2]).sqrt();

            let projected = if radial_dist > 1e-15 {
                let scale = hyp.radius / radial_dist;
                [
                    hyp.axis_origin[0] + t * hyp.axis_direction[0] + radial[0] * scale,
                    hyp.axis_origin[1] + t * hyp.axis_direction[1] + radial[1] * scale,
                    hyp.axis_origin[2] + t * hyp.axis_direction[2] + radial[2] * scale,
                ]
            } else {
                // Centroid is on the axis — unlikely but handle gracefully
                centroid
            };

            let pt = gp::Pnt::new_real3(projected[0], projected[1], projected[2]);
            let dist = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(dist);
        }

        if max_dist > config.surface_tolerance_mm {
            return Err(Stage2CompareError {
                hypothesis_type: "cylindrical",
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance: config.surface_tolerance_mm,
            });
        }
    }

    Ok(())
}


// ---------------------------------------------------------------------------
// Stage 2.3: Spherical hypothesis deduction
// ---------------------------------------------------------------------------

/// Compute the bounding box diagonal of a mesh.
fn bounding_box_diagonal(vertices: &[MeshVertex]) -> f64 {
    if vertices.is_empty() {
        return 0.0;
    }
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for v in vertices {
        min[0] = min[0].min(v.x);
        min[1] = min[1].min(v.y);
        min[2] = min[2].min(v.z);
        max[0] = max[0].max(v.x);
        max[1] = max[1].max(v.y);
        max[2] = max[2].max(v.z);
    }
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Compute the distance from a vertex to a sphere surface.
/// Returns signed distance: positive if outside the sphere, negative if inside.
fn vertex_to_sphere_distance(
    v: &MeshVertex,
    center: &[f64; 3],
    radius: f64,
) -> f64 {
    let dx = v.x - center[0];
    let dy = v.y - center[1];
    let dz = v.z - center[2];
    (dx * dx + dy * dy + dz * dz).sqrt() - radius
}

/// Fit sphere parameters (center, radius) from a set of vertices.
/// Uses algebraic least-squares: |v - c|² = r², expanded linearly.
/// Returns (center, radius) or None if degenerate.
#[allow(clippy::needless_range_loop)]
fn fit_sphere(
    vertex_set: &HashSet<usize>,
    vertices: &[MeshVertex],
) -> Option<([f64; 3], f64)> {
    let n = vertex_set.len();
    if n < 4 {
        return None;
    }

    // Solve: 2*vx*cx + 2*vy*cy + 2*vz*cz + k = vx² + vy² + vz²
    // where k = r² - |c|²
    // Normal equations: A^T A x = A^T b
    // A = [[2*vx, 2*vy, 2*vz, 1], ...], b = [vx²+vy²+vz², ...]
    let mut ata = [[0.0_f64; 4]; 4];
    let mut atb = [0.0_f64; 4];

    for &vi in vertex_set {
        let v = &vertices[vi];
        let row = [2.0 * v.x, 2.0 * v.y, 2.0 * v.z, 1.0];
        let b_i = v.x * v.x + v.y * v.y + v.z * v.z;

        for i in 0..4 {
            for j in 0..4 {
                ata[i][j] += row[i] * row[j];
            }
            atb[i] += row[i] * b_i;
        }
    }

    // Solve 4x4 system using Gaussian elimination with partial pivoting
    let mut aug = [[0.0_f64; 5]; 4];
    for i in 0..4 {
        for j in 0..4 {
            aug[i][j] = ata[i][j];
        }
        aug[i][4] = atb[i];
    }

    for col in 0..4 {
        // Find pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..4 {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            return None; // Singular
        }
        if max_row != col {
            aug.swap(col, max_row);
        }
        let pivot = aug[col][col];
        for j in col..5 {
            aug[col][j] /= pivot;
        }
        for row in 0..4 {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            for j in col..5 {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    let cx = aug[0][4];
    let cy = aug[1][4];
    let cz = aug[2][4];
    let k = aug[3][4];

    let r_sq = k + cx * cx + cy * cy + cz * cz;
    if r_sq <= 0.0 {
        return None;
    }

    Some(([cx, cy, cz], r_sq.sqrt()))
}

/// Check if all vertices are within tolerance of a sphere surface.
fn all_vertices_within_sphere_tolerance(
    vertex_set: &HashSet<usize>,
    center: &[f64; 3],
    radius: f64,
    tolerance: f64,
    vertices: &[MeshVertex],
) -> bool {
    for &vi in vertex_set {
        if vertex_to_sphere_distance(&vertices[vi], center, radius).abs() > tolerance {
            return false;
        }
    }
    true
}

/// Determine convexity for a sphere: does the face normal point away from the center?
fn determine_sphere_convexity(
    face: &MeshFace,
    vertices: &[MeshVertex],
    center: &[f64; 3],
) -> bool {
    let centroid = face_centroid(face, vertices);
    let radial = [
        centroid[0] - center[0],
        centroid[1] - center[1],
        centroid[2] - center[2],
    ];
    let n = face.normal.unwrap();
    dot3(&n, &radial) > 0.0
}

/// Deduce spherical hypotheses from the mesh using BFS region growing.
/// Validate solid-angle coverage of a spherical hypothesis.
/// Returns true if the centroid-to-center direction vectors span 3D
/// (not a 1D strip like a cylinder fillet).
fn solid_angle_coverage_valid(
    face_list: &[usize],
    center: &[f64; 3],
    faces: &[MeshFace],
    vertices: &[MeshVertex],
) -> bool {
    if face_list.len() < MIN_SPHERE_FACES {
        return false;
    }

    // Compute area-weighted covariance of unit direction vectors from center to centroid.
    let mut cov = [0.0_f64; 6]; // upper triangle: [c00, c01, c02, c11, c12, c22]
    let mut total_weight = 0.0_f64;
    for &fi in face_list {
        let c = face_centroid(&faces[fi], vertices);
        let mut d = [
            c[0] - center[0],
            c[1] - center[1],
            c[2] - center[2],
        ];
        let len = normalize3(&mut d);
        if len < 1e-15 {
            continue;
        }
        let w = face_area(&faces[fi], vertices);
        total_weight += w;
        cov[0] += w * d[0] * d[0];
        cov[1] += w * d[0] * d[1];
        cov[2] += w * d[0] * d[2];
        cov[3] += w * d[1] * d[1];
        cov[4] += w * d[1] * d[2];
        cov[5] += w * d[2] * d[2];
    }
    if total_weight < 1e-30 {
        return false;
    }
    // Normalize
    for v in &mut cov {
        *v /= total_weight;
    }

    // Compute eigenvalues of the 3x3 symmetric matrix.
    // Use the cubic formula for eigenvalues of a symmetric matrix.
    let a = cov[0];
    let b = cov[3];
    let c = cov[5];
    let d = cov[1];
    let e = cov[4];
    let f = cov[2];

    let p1 = d * d + f * f + e * e;
    if p1 < 1e-30 {
        // Already diagonal
        let mut eigs = [a, b, c];
        eigs.sort_by(|x, y| y.partial_cmp(x).unwrap());
        if eigs[0] < 1e-30 {
            return false;
        }
        return eigs[2] / eigs[0] >= MIN_SPHERE_EIGENVALUE_RATIO;
    }

    let q = (a + b + c) / 3.0;
    let p2 = (a - q) * (a - q) + (b - q) * (b - q) + (c - q) * (c - q) + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();

    // B = (1/p) * (A - q*I)
    let b00 = (a - q) / p;
    let b11 = (b - q) / p;
    let b22 = (c - q) / p;
    let b01 = d / p;
    let b02 = f / p;
    let b12 = e / p;

    let det_b = b00 * (b11 * b22 - b12 * b12)
              - b01 * (b01 * b22 - b12 * b02)
              + b02 * (b01 * b12 - b11 * b02);
    let r = det_b / 2.0;

    let phi = if r <= -1.0 {
        std::f64::consts::PI / 3.0
    } else if r >= 1.0 {
        0.0
    } else {
        r.acos() / 3.0
    };

    let eig1 = q + 2.0 * p * phi.cos();
    let eig3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).cos();
    let eig2 = 3.0 * q - eig1 - eig3;

    let lambda_max = eig1.max(eig2).max(eig3);
    let lambda_min = eig1.min(eig2).min(eig3);

    if lambda_max < 1e-30 {
        return false;
    }
    lambda_min / lambda_max >= MIN_SPHERE_EIGENVALUE_RATIO
}

/// Build a map from vertex index to the set of face indices incident on that vertex.
fn build_vertex_to_faces_map(mesh: &ConnectedMesh) -> Vec<Vec<usize>> {
    let mut vtf: Vec<Vec<usize>> = vec![Vec::new(); mesh.vertices.len()];
    for (fi, face) in mesh.faces.iter().enumerate() {
        let vc = face.vertex_count as usize;
        for &vi in &face.vertex_indices[..vc] {
            vtf[vi].push(fi);
        }
    }
    vtf
}

/// Deduce spherical hypotheses from the mesh using vertex-neighborhood seeding
/// and BFS region growing.
///
/// For each mesh vertex, collect the surrounding faces and fit a sphere to their
/// vertices. If the fit is good, seed a BFS to grow the hypothesis. After BFS,
/// validate solid-angle coverage to reject fillet-strip growth.
fn deduce_spherical_hypotheses(
    mesh: &mut ConnectedMesh,
    planar_hypotheses: &[PlanarHypothesis],
    vertex_tol: f64,
    surface_tol: f64,
    angular_tol: f64,
    max_sphere_radius: f64,
) -> Vec<SphericalHypothesis> {
    let mut hypotheses: Vec<SphericalHypothesis> = Vec::new();

    // Build vertex-to-faces map for vertex-neighborhood seeding
    let vtf = build_vertex_to_faces_map(mesh);
    for face_indices in &vtf {
        // Collect surrounding faces that are valid seed candidates
        let surrounding: Vec<usize> = face_indices
            .iter()
            .copied()
            .filter(|&fi| {
                // Must be undeduced
                if mesh.faces[fi].spherical_hypothesis != UNDEDUCED_SPHERICAL_HYPOTHESIS {
                    return false;
                }
                // Skip multi-face planar faces (genuinely flat)
                let ph = mesh.faces[fi].planar_hypothesis;
                if ph >= 0 && planar_hypotheses[ph as usize].faces.len() > 1 {
                    return false;
                }
                // Must have valid normal
                mesh.faces[fi].normal.is_some()
            })
            .collect();

        if surrounding.len() < 3 {
            continue;
        }

        // Collect all vertices from surrounding faces
        let mut seed_vertex_set: HashSet<usize> = HashSet::new();
        for &fi in &surrounding {
            let vc = mesh.faces[fi].vertex_count as usize;
            for &v in &mesh.faces[fi].vertex_indices[..vc] {
                seed_vertex_set.insert(v);
            }
        }

        // Fit sphere to seed vertices
        let (center, radius) = match fit_sphere(&seed_vertex_set, &mesh.vertices) {
            Some(params) => params,
            None => continue,
        };

        if radius > max_sphere_radius {
            continue;
        }

        // Verify all seed vertices within tolerance
        if !all_vertices_within_sphere_tolerance(
            &seed_vertex_set, &center, radius, vertex_tol, &mesh.vertices,
        ) {
            continue;
        }

        // Angular tolerance: reject if any pair of adjacent seed faces exceeds limit
        let mut angular_reject = false;
        for &fi in &surrounding {
            let fi_n = match mesh.faces[fi].normal {
                Some(n) => n,
                None => continue,
            };
            let fi_vc = mesh.faces[fi].vertex_count as usize;
            for &ni in &mesh.faces[fi].neighbors[..fi_vc] {
                if ni < 0 { continue; }
                let ni = ni as usize;
                if !surrounding.contains(&ni) { continue; }
                if let Some(ni_n) = mesh.faces[ni].normal {
                    let cos_a = dot3(&fi_n, &ni_n).clamp(-1.0, 1.0);
                    if cos_a.acos() > angular_tol {
                        angular_reject = true;
                        break;
                    }
                }
            }
            if angular_reject { break; }
        }
        if angular_reject {
            continue;
        }

        // Centroid validation for seed faces
        let mut centroid_reject = false;
        for &fi in &surrounding {
            let c = face_centroid(&mesh.faces[fi], &mesh.vertices);
            let d = vertex_to_sphere_distance(
                &MeshVertex::from_xyz(c[0], c[1], c[2]),
                &center, radius,
            ).abs();
            if d > surface_tol {
                centroid_reject = true;
                break;
            }
        }
        if centroid_reject {
            continue;
        }

        // Determine convexity from first seed face
        let convex = determine_sphere_convexity(
            &mesh.faces[surrounding[0]], &mesh.vertices, &center,
        );
        // Check convexity consistency across seed faces
        let mut convex_reject = false;
        for &fi in &surrounding[1..] {
            if determine_sphere_convexity(&mesh.faces[fi], &mesh.vertices, &center) != convex {
                convex_reject = true;
                break;
            }
        }
        if convex_reject {
            continue;
        }

        // Seed is valid — assign seed faces and start BFS
        let hi = hypotheses.len() as i32;
        for &fi in &surrounding {
            mesh.faces[fi].spherical_hypothesis = hi;
        }
        let mut face_list: Vec<usize> = surrounding.clone();
        let mut vertex_set: HashSet<usize> = seed_vertex_set;
        let mut current_center = center;
        let mut current_radius = radius;

        let mut queue: VecDeque<usize> = surrounding.iter().copied().collect();

        // BFS expansion
        while let Some(current_fi) = queue.pop_front() {
            let vc = mesh.faces[current_fi].vertex_count as usize;
            let neighbors = mesh.faces[current_fi].neighbors;

            for &cni in &neighbors[..vc] {
                if cni < 0 { continue; }
                let cni = cni as usize;

                if mesh.faces[cni].spherical_hypothesis != UNDEDUCED_SPHERICAL_HYPOTHESIS {
                    continue;
                }

                // Convexity check
                let cni_convex = determine_sphere_convexity(
                    &mesh.faces[cni], &mesh.vertices, &current_center,
                );
                if cni_convex != convex {
                    continue;
                }

                // Angular tolerance: reject if dihedral angle with ANY already-assigned
                // neighbor exceeds the limit (defense-in-depth against creased surfaces)
                if let Some(n_cni) = mesh.faces[cni].normal {
                    let cni_vc2 = mesh.faces[cni].vertex_count as usize;
                    let cni_neighbors = mesh.faces[cni].neighbors;
                    let mut ang_reject = false;
                    for &adj in &cni_neighbors[..cni_vc2] {
                        if adj < 0 { continue; }
                        let adj = adj as usize;
                        if mesh.faces[adj].spherical_hypothesis != hi { continue; }
                        if let Some(n_adj) = mesh.faces[adj].normal {
                            let cos_a = dot3(&n_adj, &n_cni).clamp(-1.0, 1.0);
                            if cos_a.acos() > angular_tol {
                                ang_reject = true;
                                break;
                            }
                        }
                    }
                    if ang_reject { continue; }
                }

                // Vertex distance check
                let cni_vc = mesh.faces[cni].vertex_count as usize;
                let cni_vi = mesh.faces[cni].vertex_indices;
                let mut all_ok = true;
                let mut any_far = false;
                for &v in &cni_vi[..cni_vc] {
                    let d = vertex_to_sphere_distance(
                        &mesh.vertices[v],
                        &current_center,
                        current_radius,
                    ).abs();
                    if d > vertex_tol {
                        all_ok = false;
                        if d > REFIT_SKIP_MULTIPLIER * vertex_tol {
                            any_far = true;
                            break;
                        }
                    }
                }
                if any_far { continue; }

                // Centroid validation (surface_tol check during BFS)
                let centroid = face_centroid(&mesh.faces[cni], &mesh.vertices);
                let centroid_dist = vertex_to_sphere_distance(
                    &MeshVertex::from_xyz(centroid[0], centroid[1], centroid[2]),
                    &current_center, current_radius,
                ).abs();
                if centroid_dist > REFIT_SKIP_MULTIPLIER * surface_tol {
                    continue;
                }
                let needs_refit = !all_ok || centroid_dist > surface_tol;

                if needs_refit {
                    // Try re-fitting with new face included
                    let mut trial_vertices = vertex_set.clone();
                    for &v in &cni_vi[..cni_vc] {
                        trial_vertices.insert(v);
                    }

                    let (new_center, new_radius) = match fit_sphere(
                        &trial_vertices, &mesh.vertices,
                    ) {
                        Some(params) => params,
                        None => continue,
                    };

                    if new_radius > max_sphere_radius {
                        continue;
                    }

                    if !all_vertices_within_sphere_tolerance(
                        &trial_vertices, &new_center, new_radius,
                        vertex_tol, &mesh.vertices,
                    ) {
                        continue;
                    }

                    // Verify all existing centroids still pass
                    let mut refit_centroid_ok = true;
                    for &f in &face_list {
                        let c = face_centroid(&mesh.faces[f], &mesh.vertices);
                        let d = vertex_to_sphere_distance(
                            &MeshVertex::from_xyz(c[0], c[1], c[2]),
                            &new_center, new_radius,
                        ).abs();
                        if d > surface_tol {
                            refit_centroid_ok = false;
                            break;
                        }
                    }
                    if !refit_centroid_ok {
                        continue;
                    }

                    // Accept re-fit
                    current_center = new_center;
                    current_radius = new_radius;
                }

                // Accept this face
                mesh.faces[cni].spherical_hypothesis = hi;
                face_list.push(cni);
                for &v in &cni_vi[..cni_vc] {
                    vertex_set.insert(v);
                }
                queue.push_back(cni);
            }
        }

        // Final re-fit from all accumulated faces
        if let Some((final_center, final_radius)) =
            fit_sphere(&vertex_set, &mesh.vertices)
        {
            current_center = final_center;
            current_radius = final_radius;
        }

        // Compute error metrics
        let mut error_max = 0.0_f64;
        let mut error_abs_sum = 0.0_f64;
        for &v in &vertex_set {
            let d = vertex_to_sphere_distance(
                &mesh.vertices[v], &current_center, current_radius,
            ).abs();
            error_max = error_max.max(d);
            error_abs_sum += d;
        }

        // Validate: minimum face count, centroid check, radius, and solid-angle coverage
        let min_faces_ok = face_list.len() >= MIN_SPHERE_FACES;
        let radius_ok = current_radius <= max_sphere_radius;
        let coverage_ok = min_faces_ok && solid_angle_coverage_valid(
            &face_list, &current_center, &mesh.faces, &mesh.vertices,
        );

        if !min_faces_ok || !radius_ok || !coverage_ok {
            // Undo assignments
            for &f in &face_list {
                mesh.faces[f].spherical_hypothesis = UNDEDUCED_SPHERICAL_HYPOTHESIS;
            }
            continue;
        }

        hypotheses.push(SphericalHypothesis {
            center: current_center,
            radius: current_radius,
            convex,
            faces: face_list,
            vertices: vertex_set.into_iter().collect(),
            error_max,
            error_abs_sum,
        });
    }

    // Mark remaining undeduced faces as NO_HYPOTHESIS
    for face in &mut mesh.faces {
        if face.spherical_hypothesis == UNDEDUCED_SPHERICAL_HYPOTHESIS {
            face.spherical_hypothesis = NO_HYPOTHESIS;
        }
    }

    hypotheses
}

/// Validate fitted spherical hypotheses against a reference STEP shape.
fn compare_spherical_hypotheses(
    hypotheses: &[SphericalHypothesis],
    mesh: &ConnectedMesh,
    config: &Config,
) -> Result<(), Stage2CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    for (hi, hyp) in hypotheses.iter().enumerate() {
        let mut max_dist = 0.0_f64;

        for &fi in &hyp.faces {
            let centroid = face_centroid(&mesh.faces[fi], &mesh.vertices);

            // Project centroid onto sphere surface (nearest point)
            let d = [
                centroid[0] - hyp.center[0],
                centroid[1] - hyp.center[1],
                centroid[2] - hyp.center[2],
            ];
            let dist_to_center = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();

            let projected = if dist_to_center > 1e-15 {
                let scale = hyp.radius / dist_to_center;
                [
                    hyp.center[0] + d[0] * scale,
                    hyp.center[1] + d[1] * scale,
                    hyp.center[2] + d[2] * scale,
                ]
            } else {
                // Centroid is at the center — unlikely but handle gracefully
                centroid
            };

            let pt = gp::Pnt::new_real3(projected[0], projected[1], projected[2]);
            let dist = stage1::min_distance_to_shape(&pt, compare_shape);
            max_dist = max_dist.max(dist);
        }

        if max_dist > config.surface_tolerance_mm {
            return Err(Stage2CompareError {
                hypothesis_type: "spherical",
                hypothesis_index: hi,
                max_distance: max_dist,
                tolerance: config.surface_tolerance_mm,
            });
        }
    }

    Ok(())
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
) -> Vec<SelectedSurface> {
    let num_faces = mesh.faces.len();

    // Step 1: Compute geometric area of each mesh face.
    let face_areas: Vec<f64> = (0..num_faces)
        .map(|fi| face_area(&mesh.faces[fi], &mesh.vertices))
        .collect();

    // Step 2: Build candidate list from multi-face hypotheses.
    // Each candidate tracks its hypothesis type/index and original face set.
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
/// For each selected surface, project face centroids onto the fitted surface
/// and check that the projected points are within surface_tolerance of the STEP file.
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
pub fn stage2(config: &Config, mut mesh: ConnectedMesh) -> Result<Stage2Output, Stage2Error> {
    // Stage 2.1: Deduce planar hypotheses
    let mut planar_hypotheses = deduce_planar_hypotheses(&mut mesh, config.vertex_tolerance_mm);

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
    let mut cylindrical_hypotheses = deduce_cylindrical_hypotheses(
        &mut mesh, config.vertex_tolerance_mm,
        config.surface_tolerance_mm, config.angular_tolerance_rad,
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
                    "  Cylinder {}: {} faces, {} vertices, r={:.4}, {}, axis=[{:.4}, {:.4}, {:.4}], err_max={:.2e}",
                    i, h.faces.len(), h.vertices.len(), h.radius,
                    if h.convex { "convex" } else { "concave" },
                    h.axis_direction[0], h.axis_direction[1], h.axis_direction[2],
                    h.error_max,
                );
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
    let mut spherical_hypotheses = deduce_spherical_hypotheses(
        &mut mesh, &planar_hypotheses, config.vertex_tolerance_mm,
        config.surface_tolerance_mm, config.angular_tolerance_rad, max_sphere_radius,
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
                    "  Sphere {}: {} faces, {} vertices, r={:.4}, {}, center=[{:.4}, {:.4}, {:.4}], err_max={:.2e}",
                    i, h.faces.len(), h.vertices.len(), h.radius,
                    if h.convex { "convex" } else { "concave" },
                    h.center[0], h.center[1], h.center[2],
                    h.error_max,
                );
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
    let selected_surfaces = select_surfaces(
        &mesh, &mut planar_hypotheses, &mut cylindrical_hypotheses,
        &mut spherical_hypotheses,
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
    use crate::stage1::{self, MeshFace, MeshVertex, VertexWeldOptions, UNDEDUCED_CYLINDRICAL_HYPOTHESIS, UNDEDUCED_SPHERICAL_HYPOTHESIS};

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

        let hypotheses = deduce_planar_hypotheses(&mut mesh, 1e-5);

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
        let hypotheses = deduce_planar_hypotheses(&mut mesh, 1e-5);

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
        let hypotheses = deduce_planar_hypotheses(&mut mesh, 1e-5);

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
        let planar = deduce_planar_hypotheses(&mut mesh, 1e-5);
        // 17.5° angular tolerance: cube faces meet at 90°, so no cylinders
        let cyls = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 17.5_f64.to_radians(),
        );
        assert_eq!(cyls.len(), 0, "cube should produce 0 cylinders at 17.5° angular tolerance");
    }

    #[test]
    fn cube_angular_tolerance_rejects_sphere() {
        let mut mesh = load_stage1("manual/cube.stl");
        let planar = deduce_planar_hypotheses(&mut mesh, 1e-5);
        let _ = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 17.5_f64.to_radians(),
        );
        // 17.5° angular tolerance: cube faces meet at 90°, so no spheres
        let sphs = deduce_spherical_hypotheses(
            &mut mesh, &planar, 1e-5, 0.4, 17.5_f64.to_radians(), 100.0,
        );
        assert_eq!(sphs.len(), 0, "cube should produce 0 spheres at 17.5° angular tolerance");
    }

    #[test]
    fn cube_high_angular_tolerance_allows_sphere() {
        let mut mesh = load_stage1("manual/cube.stl");
        let planar = deduce_planar_hypotheses(&mut mesh, 1e-5);
        let _ = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 91.0_f64.to_radians(),
        );
        // 91° angular tolerance: exceeds 90° dihedral angle of cube
        let sphs = deduce_spherical_hypotheses(
            &mut mesh, &planar, 1e-5, 0.4, 91.0_f64.to_radians(), 100.0,
        );
        assert!(!sphs.is_empty(), "cube should produce spheres at 91° angular tolerance");
    }

    #[test]
    fn cylinder_detected_at_default_angular_tolerance() {
        let mut mesh = load_stage1("ccad/generated/simple_cylinder.stl");
        let planar = deduce_planar_hypotheses(&mut mesh, 1e-5);
        let cyls = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 17.5_f64.to_radians(),
        );
        assert_eq!(cyls.len(), 1, "simple cylinder should produce 1 cylinder hypothesis");
        assert!(cyls[0].convex);
        assert!((cyls[0].radius - 10.0).abs() < 0.01);
    }

    #[test]
    fn sphere_detected_at_default_angular_tolerance() {
        let mut mesh = load_stage1("ccad/generated/simple_sphere.stl");
        let planar = deduce_planar_hypotheses(&mut mesh, 1e-5);
        let _ = deduce_cylindrical_hypotheses(
            &mut mesh, 1e-5, 0.4, 17.5_f64.to_radians(),
        );
        let sphs = deduce_spherical_hypotheses(
            &mut mesh, &planar, 1e-5, 0.4, 17.5_f64.to_radians(), 1000.0,
        );
        assert_eq!(sphs.len(), 1, "simple sphere should produce 1 sphere hypothesis");
        assert!(sphs[0].convex);
        assert!((sphs[0].radius - 10.0).abs() < 0.01);
    }
}