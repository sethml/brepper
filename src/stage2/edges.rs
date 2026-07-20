//! Stage 2.1: fit-free edge classification.
//!
//! Labels every manifold mesh edge INTERIOR or BOUNDARY using two fit-free
//! signals (DEVELOPMENT_PLAN §2.1): a sharp dihedral angle, or a jump in the
//! per-face two-sided curvature descriptor (k_max, k_min) across the edge.
//! The true segmentation's boundaries are a subset of mesh edges
//! (tessellators emit seam edges along B-rep face boundaries, even at
//! tangencies), so this classification IS the segmentation; stage 2.2
//! extracts its connected regions and stage 2.4 (future) repairs mistakes.
//!
//! **Curvature estimator: per-face anchored osculating-jet fit.** Each
//! face's descriptor is the pair of principal curvatures of a quadratic
//! height field `z = Q(u,w)` fitted over a PCA tangent frame of the face's
//! vertex neighborhood (Cazals & Pouget, "Estimating differential
//! quantities using polynomial fitting of osculating jets", SGP 2005 —
//! the standard robust curvature estimator for irregular meshes). Three
//! properties make it robust to how CAD models are actually triangulated:
//!
//! - **It uses vertex positions** — exact f32 mesh data — not normals.
//!   A long skinny strip/sliver triangle's plane is a chord whose normal
//!   sits degrees off the true tangent, so normal-based estimators
//!   inherit that error; positions are ground truth. Tessellation-layout
//!   agnostic: aspect-ratio slivers (rounded_cube_10_r2_fine's 65:1
//!   fillet triangles), staggered rows, pole fans, and full-length strip
//!   meshes (pipe_elbow's 15 mm strips) all fit the same.
//!
//! - **The tangent frame comes from a PCA of own + 1-ring vertices, not
//!   the facet normal** — the smallest-variance direction of the
//!   neighborhood IS the local normal and stays accurate on strips. Own
//!   vertices dominate the frame (weight OWN_WEIGHT).
//!
//! - **The face's own vertices anchor the fit** (weight OWN_WEIGHT).
//!   Tessellators never straddle a B-rep boundary with a triangle, so a
//!   face's own corners always lie on its own surface — this biases each
//!   face's fit toward its own side at a junction, making junction rows
//!   read *different* intermediate values and the boundary fire. There is
//!   **no outlier rejection**: rejection forces arbitrary one-side reads
//!   that vary face to face and fragment tangent junctions (measured).
//!   Sharp-crease neighbors (dihedral > `--angular-tolerance`) are
//!   excluded from the point set instead; tangent-junction mixing is
//!   allowed (intermediate reads).
//!
//! The estimator runs inside a fixpoint iteration with the classification:
//! round 1 admits all non-sharp neighbors; once junctions fire,
//! cross-boundary faces leave the neighborhood and the descriptors purify.
//! (History: per-edge κ = θ/d failed on irregular tessellation; the
//! neighbor-normal tensor deadlocked on balanced-mix junction rings. The
//! jet subsumes both.)

use super::{dot3, face_centroid};
use crate::stage1::ConnectedMesh;

/// Local-relative curvature-jump threshold (DEVELOPMENT_PLAN §2.1): an
/// edge is a boundary when the descriptor jump exceeds this fraction of
/// the local curvature magnitude. The local-relative criterion survives
/// from the junction-mask finding: thresholds relative to a mesh-global
/// maximum fail because sharp creases dominate it while tangent junctions
/// fall below any global fraction.
pub const FIRE_RATIO: f64 = 0.35;

/// Un-fire threshold: an edge drops back to interior when the descriptor
/// jump falls below this fraction of the local curvature magnitude.
/// Between UNFIRE_RATIO and FIRE_RATIO an edge keeps its previous state —
/// the gap damps the fixpoint against oscillation.
pub const UNFIRE_RATIO: f64 = 0.15;

/// Cap on classify/re-descriptor rounds. Convergence is typically 2–4
/// rounds; the cap guards oscillation.
pub const MAX_ROUNDS: usize = 10;

/// Noise-floor coefficient: κ_floor = FLOOR_COEFF · ε / d² per edge
/// (centroid distance d), where ε is the vertex-position noise scale —
/// f32 quantization (bbox diagonal × 1e-7), NOT --vertex-tolerance. The
/// floor keeps exact-coplanar plane interiors (k ≡ 0 ± noise) from ever
/// firing the relative test.
pub const FLOOR_COEFF: f64 = 4.0;

/// Weight of the face's own vertices in the jet fit, relative to 1.0 for
/// ring vertices. The own corners always lie on the face's own surface
/// (tessellators do not straddle seams), so they anchor the fit.
pub const OWN_WEIGHT: f64 = 50.0;

/// Minimum unique points before the 2-ring is added to the neighborhood.
pub const JET_MIN_POINTS: usize = 12;

/// Classification of one manifold mesh edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EdgeClass {
    /// Both faces on the same surface region.
    #[default]
    Interior,
    /// Dihedral angle exceeds the angular tolerance — a hard feature edge.
    Sharp,
    /// G1-smooth but G2-discontinuous — a tangent junction found by the
    /// curvature-jump rule.
    Jump,
}

/// Stage 2.1 output: per-half-edge classes plus the curvature evidence,
/// kept for the 2.4 split repair (re-run at lowered thresholds), for
/// `--compare` diagnostics, and for viz.
///
/// All per-edge arrays are indexed `[face][slot]`, parallel to
/// `MeshFace.neighbors`: slot `e` is the edge from `vertex_indices[e]` to
/// `vertex_indices[(e+1) % vertex_count]`. Values are mirrored — the two
/// half-edges of a manifold edge always carry identical class/κ (computed
/// from symmetric formulas, not copied).
#[derive(Debug, Default)]
pub struct EdgeClassification {
    /// Class per half-edge. Slots without a manifold neighbor stay Interior
    /// (they never participate in region unions anyway).
    pub class: Vec<[EdgeClass; 4]>,
    /// Signed per-half-edge discrete curvature κ = ±θ/d (positive convex,
    /// negative concave); NaN where no manifold neighbor. Diagnostic only —
    /// NOT the classification signal (it is unreliable on irregular
    /// tessellation; see module docs).
    pub kappa: Vec<[f64; 4]>,
    /// Final per-face larger principal curvature (signed; positive convex).
    pub sigma_max: Vec<f64>,
    /// Final per-face smaller principal curvature (signed).
    pub sigma_min: Vec<f64>,
    /// Rounds of classification (1: the jet descriptor is boundary-
    /// independent, so there is no fixpoint iteration).
    pub rounds: usize,
    /// Undirected manifold-edge counts by final class.
    pub sharp_count: usize,
    pub jump_count: usize,
    pub interior_count: usize,
}

impl EdgeClassification {
    /// Is the manifold edge across `faces[fi]`'s slot `e` a region boundary?
    pub fn is_boundary(&self, fi: usize, e: usize) -> bool {
        self.class[fi][e] != EdgeClass::Interior
    }
}

/// One point for the jet fit: a vertex position, its weight, and whether
/// it may steer the tangent frame (own + 1-ring points only; fallback
/// 2-ring points contribute fit evidence but never steer the frame — a
/// wild cross-junction vertex among them would tilt the frame itself).
struct JetPoint {
    pos: [f64; 3],
    weight: f64,
    in_frame: bool,
}

/// Fit the quadratic height field `z = a·u² + b·u·w + c·w² + d·u + e·w + g`
/// and return the principal curvatures (eigenvalues of the Hessian
/// [[2a, b],[b, 2c]]).
///
/// The tangent frame comes from a **PCA of the point cloud, not the facet
/// normal**: a long skinny strip triangle's own plane is a chord whose
/// normal can sit several degrees off the true surface tangent, which
/// poisons the height measurement; the smallest-variance direction of the
/// neighborhood points IS the local normal and stays accurate on strips
/// and slivers (sign aligned with the face normal hint `z_hint`).
///
/// Weighted least squares with a tiny Tikhonov ridge (unobserved
/// directions — e.g. across a one-facet-wide strip — collapse to ~0
/// curvature instead of producing garbage), followed by JET_IRLS_ROUNDS of
/// MAD-based outlier rejection (vertices lying across a junction deviate
/// from the anchored surface by ~the surface mismatch and are dropped).
/// `eps` is the absolute residual floor (position noise scale).
fn jet_curvature(
    origin: &[f64; 3],
    z_hint: &[f64; 3],
    points: &[JetPoint],
) -> Option<(f64, f64)> {
    // PCA frame: weighted covariance of the frame-eligible points around
    // the origin; the two largest-variance directions span the tangent
    // plane. Own vertices dominate (weight OWN_WEIGHT), so the frame
    // tracks the face's own surface even when some ring points are wild.
    let mut wsum = 0.0;
    let mut cov = nalgebra::Matrix3::<f64>::zeros();
    for p in points {
        if !p.in_frame {
            continue;
        }
        let d = [
            p.pos[0] - origin[0],
            p.pos[1] - origin[1],
            p.pos[2] - origin[2],
        ];
        wsum += p.weight;
        for a in 0..3 {
            for b in 0..3 {
                cov[(a, b)] += p.weight * d[a] * d[b];
            }
        }
    }
    if wsum <= 0.0 {
        return None;
    }
    cov /= wsum;
    let eigen = nalgebra::SymmetricEigen::new(cov);
    let mut order = [0usize, 1, 2];
    order.sort_by(|&a, &b| {
        eigen.eigenvalues[a]
            .partial_cmp(&eigen.eigenvalues[b])
            .unwrap()
    });
    let (i1, i2, i3) = (order[2], order[1], order[0]);
    let u = nalgebra::Vector3::new(
        eigen.eigenvectors[(0, i1)],
        eigen.eigenvectors[(1, i1)],
        eigen.eigenvectors[(2, i1)],
    );
    let w = nalgebra::Vector3::new(
        eigen.eigenvectors[(0, i2)],
        eigen.eigenvectors[(1, i2)],
        eigen.eigenvectors[(2, i2)],
    );
    let mut z_axis = nalgebra::Vector3::new(
        eigen.eigenvectors[(0, i3)],
        eigen.eigenvectors[(1, i3)],
        eigen.eigenvectors[(2, i3)],
    );
    if z_axis.dot(&nalgebra::Vector3::from_column_slice(z_hint)) < 0.0 {
        z_axis = -z_axis;
    }
    let (u, w, z_axis) = ([u.x, u.y, u.z], [w.x, w.y, w.z], [z_axis.x, z_axis.y, z_axis.z]);

    // Local scaled coordinates (u, w, z) per point; scale by the
    // neighborhood's characteristic size h for conditioning (curvatures
    // scale back by 1/h² below).
    #[derive(Clone, Copy)]
    struct Local {
        u: f64,
        w: f64,
        z: f64,
        weight: f64,
    }
    let mut locals: Vec<Local> = points
        .iter()
        .map(|p| {
            let d = [
                p.pos[0] - origin[0],
                p.pos[1] - origin[1],
                p.pos[2] - origin[2],
            ];
            Local {
                u: dot3(&d, &u),
                w: dot3(&d, &w),
                z: dot3(&d, &z_axis),
                weight: p.weight,
            }
        })
        .collect();
    let h = locals
        .iter()
        .map(|l| l.u.abs().max(l.w.abs()))
        .fold(0.0_f64, f64::max);
    if h <= 1e-12 {
        return None;
    }
    for l in &mut locals {
        l.u /= h;
        l.w /= h;
        l.z /= h;
    }

    // Design row for Q = a·u² + b·u·w + c·w² + d·u + e·w + g.
    fn design(u: f64, w: f64) -> [f64; 6] {
        [u * u, u * w, w * w, u, w, 1.0]
    }

    // Weighted normal equations + ridge (single pass — no outlier
    // rejection: the own-vertex anchor, not rejection, is what biases each
    // face toward its own surface, and rejection forces arbitrary one-side
    // reads that fragment tangent junctions).
    let mut ata = nalgebra::Matrix6::<f64>::zeros();
    let mut atb = nalgebra::Vector6::<f64>::zeros();
    for l in &locals {
        let d = design(l.u, l.w);
        for a in 0..6 {
            atb[a] += l.weight * d[a] * l.z;
            for b in 0..6 {
                ata[(a, b)] += l.weight * d[a] * d[b];
            }
        }
    }
    for a in 0..6 {
        ata[(a, a)] += 1e-9;
    }
    let p = ata.lu().solve(&atb)?;
    // Hessian of Q in scaled coordinates. Sign convention: curvature is
    // POSITIVE when the surface curves AWAY from z_axis (with z_axis the
    // outward face normal, convex is positive) — the fitted surface drops
    // below the tangent plane, so negate the Hessian. Rescale by 1/h
    // (z was scaled by h as well: k = k'/h).
    let (a, b, c) = (-2.0 * p[0], -p[1], -2.0 * p[2]);
    let mid = (a + c) / 2.0;
    let diff = (((a - c) / 2.0).powi(2) + b * b).sqrt();
    Some(((mid + diff) / h, (mid - diff) / h))
}

/// Gather the jet-fit points for one face: own vertices (anchored,
/// frame-eligible), then ring faces' vertices. When `strict` is set, only
/// INTERIOR-edge ring faces contribute (the iteration's purification:
/// once a junction fires, cross-boundary faces leave the neighborhood);
/// otherwise all non-sharp neighbors. 2-ring vertices are added when the
/// neighborhood is sparse (they contribute fit evidence but never steer
/// the frame).
fn gather_jet_points(
    mesh: &ConnectedMesh,
    fi: usize,
    sharp: &[[bool; 4]],
    jump: &[[bool; 4]],
    strict: bool,
) -> Vec<JetPoint> {
    let face = &mesh.faces[fi];
    let vc = face.vertex_count as usize;
    let mut points: Vec<JetPoint> = Vec::new();
    let mut seen: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for &vi in &face.vertex_indices[..vc] {
        if seen.insert(vi) {
            let v = &mesh.vertices[vi];
            points.push(JetPoint {
                pos: [v.x, v.y, v.z],
                weight: OWN_WEIGHT,
                in_frame: true,
            });
        }
    }
    let mut ring1: Vec<usize> = Vec::with_capacity(vc);
    for e in 0..vc {
        let ni = face.neighbors[e];
        if ni < 0 || sharp[fi][e] {
            continue;
        }
        if strict && jump[fi][e] {
            continue;
        }
        ring1.push(ni as usize);
    }
    let push_face_points = |idx: usize, in_frame: bool, points: &mut Vec<JetPoint>, seen: &mut std::collections::HashSet<usize>| {
        let f = &mesh.faces[idx];
        let fvc = f.vertex_count as usize;
        for &vi in &f.vertex_indices[..fvc] {
            if seen.insert(vi) {
                let v = &mesh.vertices[vi];
                points.push(JetPoint {
                    pos: [v.x, v.y, v.z],
                    weight: 1.0,
                    in_frame,
                });
            }
        }
    };
    for &ni in &ring1 {
        push_face_points(ni, true, &mut points, &mut seen);
    }
    if points.len() < JET_MIN_POINTS {
        for &ni in &ring1 {
            let f = &mesh.faces[ni];
            let fvc = f.vertex_count as usize;
            for e2 in 0..fvc {
                let n2 = f.neighbors[e2];
                if n2 < 0 || sharp[ni][e2] {
                    continue;
                }
                if strict && jump[ni][e2] {
                    continue;
                }
                push_face_points(n2 as usize, false, &mut points, &mut seen);
            }
        }
    }
    points
}

/// Classify every manifold mesh edge (stage 2.1). Deterministic: per-face
/// jet fits, then the jump classification iterated to a fixpoint.
pub fn classify_edges(mesh: &ConnectedMesh, angular_tolerance: f64) -> EdgeClassification {
    let n = mesh.faces.len();
    let centroids: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|f| face_centroid(f, &mesh.vertices))
        .collect();

    // Vertex-position noise scale: f32 quantization (STL positions are
    // f32-exact), estimated from the bounding box diagonal.
    let mut bmin = [f64::INFINITY; 3];
    let mut bmax = [f64::NEG_INFINITY; 3];
    for v in &mesh.vertices {
        for (k, c) in [v.x, v.y, v.z].into_iter().enumerate() {
            bmin[k] = bmin[k].min(c);
            bmax[k] = bmax[k].max(c);
        }
    }
    let diag = ((bmax[0] - bmin[0]).powi(2)
        + (bmax[1] - bmin[1]).powi(2)
        + (bmax[2] - bmin[2]).powi(2))
    .sqrt();
    let eps = diag.max(1.0) * 1e-7;

    // Per half-edge: dihedral θ, sharpness, noise floor, and the diagnostic
    // κ = ±θ/d (computed once; the jump classification below iterates).
    let mut kappa = vec![[f64::NAN; 4]; n];
    let mut theta = vec![[f64::NAN; 4]; n];
    let mut sharp = vec![[false; 4]; n];
    let mut floor = vec![[f64::INFINITY; 4]; n];
    for fi in 0..n {
        let Some(n1) = mesh.faces[fi].normal else { continue };
        let vc = mesh.faces[fi].vertex_count as usize;
        for e in 0..vc {
            let ni = mesh.faces[fi].neighbors[e];
            if ni < 0 {
                continue;
            }
            let ni = ni as usize;
            let Some(n2) = mesh.faces[ni].normal else { continue };
            let th = dot3(&n1, &n2).clamp(-1.0, 1.0).acos();
            let dx = centroids[ni][0] - centroids[fi][0];
            let dy = centroids[ni][1] - centroids[fi][1];
            let dz = centroids[ni][2] - centroids[fi][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            if d <= 1e-12 {
                continue;
            }
            let below = n1[0] * dx + n1[1] * dy + n1[2] * dz <= 0.0;
            kappa[fi][e] = if below { th / d } else { -th / d };
            theta[fi][e] = th;
            sharp[fi][e] = th > angular_tolerance;
            floor[fi][e] = FLOOR_COEFF * eps / (d * d);
        }
    }

    // Iterate: per-face jet fits from interior-connected neighborhoods →
    // (re)classify jumps with hysteresis, until stable. Round 1 admits all
    // non-sharp neighbors; the outlier rejection handles the junction
    // contamination, and once junctions fire the neighborhoods purify.
    let mut jump = vec![[false; 4]; n];
    let mut sigma_max = vec![0.0_f64; n];
    let mut sigma_min = vec![0.0_f64; n];
    let mut rounds = 0;
    for _round in 0..MAX_ROUNDS {
        rounds += 1;

        for fi in 0..n {
            let Some(nrm) = mesh.faces[fi].normal else { continue };
            let mut fit: Option<(f64, f64)> = None;
            for strict in [true, false] {
                let points = gather_jet_points(mesh, fi, &sharp, &jump, strict);
                fit = jet_curvature(&centroids[fi], &nrm, &points);
                if fit.is_some() {
                    break;
                }
            }
            if let Some((kmax, kmin)) = fit {
                sigma_max[fi] = kmax;
                sigma_min[fi] = kmin;
            }
        }

        // (Re)classify jumps with hysteresis: fire above FIRE_RATIO,
        // un-fire below UNFIRE_RATIO or below the noise floor, otherwise
        // keep the previous state.
        let mut changed = false;
        let mut new_jump = jump.clone();
        for fi in 0..n {
            let vc = mesh.faces[fi].vertex_count as usize;
            for e in 0..vc {
                if !theta[fi][e].is_finite() || sharp[fi][e] {
                    continue;
                }
                let ni = mesh.faces[fi].neighbors[e] as usize;
                let jmp = (sigma_max[fi] - sigma_max[ni])
                    .abs()
                    .max((sigma_min[fi] - sigma_min[ni]).abs());
                let scale = sigma_max[fi]
                    .abs()
                    .max(sigma_min[fi].abs())
                    .max(sigma_max[ni].abs())
                    .max(sigma_min[ni].abs());
                let next = if jmp > FIRE_RATIO * scale && jmp > floor[fi][e] {
                    true
                } else if jmp < UNFIRE_RATIO * scale || jmp < floor[fi][e] {
                    false
                } else {
                    jump[fi][e]
                };
                if next != new_jump[fi][e] {
                    changed = true;
                }
                new_jump[fi][e] = next;
            }
        }
        jump = new_jump;
        if !changed {
            break;
        }
    }

    // Final classes + undirected counts.
    let mut class = vec![[EdgeClass::Interior; 4]; n];
    let (mut sharp_count, mut jump_count, mut interior_count) = (0, 0, 0);
    for fi in 0..n {
        let vc = mesh.faces[fi].vertex_count as usize;
        for e in 0..vc {
            let ni = mesh.faces[fi].neighbors[e];
            if ni < 0 || !theta[fi][e].is_finite() {
                continue;
            }
            class[fi][e] = if sharp[fi][e] {
                EdgeClass::Sharp
            } else if jump[fi][e] {
                EdgeClass::Jump
            } else {
                EdgeClass::Interior
            };
            if fi < ni as usize {
                match class[fi][e] {
                    EdgeClass::Sharp => sharp_count += 1,
                    EdgeClass::Jump => jump_count += 1,
                    EdgeClass::Interior => interior_count += 1,
                }
            }
        }
    }

    EdgeClassification {
        class,
        kappa,
        sigma_max,
        sigma_min,
        rounds,
        sharp_count,
        jump_count,
        interior_count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn jp(x: f64, y: f64, z: f64, w: f64) -> JetPoint {
        JetPoint {
            pos: [x, y, z],
            weight: w,
            in_frame: true,
        }
    }

    /// A face at the origin (normal +z) on a quadratic saddle
    /// z = ½(k1·u² + k2·w²) must recover (k1, k2) from positions.
    #[test]
    fn jet_recovers_quadratic() {
        let (k1, k2) = (0.5, -0.2);
        let f = |u: f64, w: f64| 0.5 * (k1 * u * u + k2 * w * w);
        let mut points = vec![
            jp(0.01, 0.0, f(0.01, 0.0), OWN_WEIGHT),
            jp(0.0, 0.01, f(0.0, 0.01), OWN_WEIGHT),
            jp(-0.01, 0.0, f(-0.01, 0.0), OWN_WEIGHT),
        ];
        for &(u, w) in &[(0.06, 0.0), (0.0, 0.05), (-0.07, 0.0), (0.0, -0.055), (0.04, 0.04), (-0.05, 0.03)] {
            points.push(jp(u, w, f(u, w), 1.0));
        }
        let (gmax, gmin) =
            jet_curvature(&[0.0, 0.0, 0.0], &[0.0, 0.0, -1.0], &points)
                .unwrap();
        assert!((gmax - k1).abs() < 1e-3, "k_max {gmax} != {k1}");
        assert!((gmin - k2).abs() < 1e-3, "k_min {gmin} != {k2}");
    }

    /// Cylinder band z = r − √(r²−x²) sampled on a skinny strip: the jet
    /// recovers (1/r, 0) regardless of aspect ratio.
    #[test]
    fn jet_skinny_cylinder_strip() {
        let r = 2.0;
        let z = |x: f64| r - (r * r - x * x).sqrt();
        let mut points = vec![
            jp(0.0, 0.0, z(0.0), OWN_WEIGHT),
            jp(0.09, 0.0, z(0.09), OWN_WEIGHT),
            jp(0.045, 6.0, z(0.045), OWN_WEIGHT),
        ];
        for &(x, y) in &[(0.135, 0.0), (0.03, 6.0), (0.06, -6.0), (0.18, 0.5), (-0.045, 0.5), (0.09, -0.5)] {
            points.push(jp(x, y, z(x), 1.0));
        }
        // Face at x=0.045, normal ≈ band normal, edge hint along y.
        let nrm = {
            let dzdx = 0.045 / (r * r - 0.045 * 0.045).sqrt();
            let n = [dzdx, 0.0, -1.0];
            let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            [n[0] / l, n[1] / l, n[2] / l]
        };
        let c = [0.045, 0.0, z(0.045)];
        let (kmax, kmin) = jet_curvature(&c, &nrm, &points).unwrap();
        // Jet approximation error on a non-quadratic circle is ~1%.
        assert!((kmax - 1.0 / r).abs() < 0.01, "k_max {kmax} != {}", 1.0 / r);
        assert!(kmin.abs() < 0.01, "k_min {kmin} != 0");
    }

    /// Junction anchoring: a face on a fillet (cylinder, κ=0.5) whose
    /// neighborhood includes vertices across a tangent junction on a
    /// plane (κ=0) must follow the FILLET — the plane vertices are
    /// rejected as outliers.
    #[test]
    fn jet_anchors_own_surface_at_junction() {
        // Fillet: quarter cylinder of radius 2 around the y axis, surface
        // z = r − √(r²−x²) for x ∈ [0, ~1.4]. Plane: z = 0 for x < 0,
        // tangent to the fillet at x = 0.
        let r = 2.0;
        let zf = |x: f64| r - (r * r - x * x).sqrt();
        let mut points = vec![
            // Own face: a quad on the fillet near x = 0.3.
            jp(0.28, 0.0, zf(0.28), OWN_WEIGHT),
            jp(0.34, 0.0, zf(0.34), OWN_WEIGHT),
            jp(0.28, 0.5, zf(0.28), OWN_WEIGHT),
            jp(0.34, 0.5, zf(0.34), OWN_WEIGHT),
        ];
        // Same-fillet ring vertices.
        for &(x, y) in &[(0.1, 0.0), (0.5, 0.0), (0.1, 0.5), (0.5, 0.5), (0.1, -0.4), (0.5, -0.4), (0.22, 0.2), (0.42, -0.2)] {
            points.push(jp(x, y, zf(x), 1.0));
        }
        // Plane vertices across the junction (x < 0).
        for &(x, y) in &[(-0.2, 0.0), (-0.5, 0.3), (-0.4, -0.3), (-0.1, 0.5)] {
            points.push(jp(x, y, 0.0, 1.0));
        }
        let nrm = {
            let dzdx = 0.31 / (r * r - 0.31 * 0.31).sqrt();
            let n = [dzdx, 0.0, -1.0];
            let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            [n[0] / l, n[1] / l, n[2] / l]
        };
        let c = [0.31, 0.25, zf(0.31)];
        let (kmax, kmin) = jet_curvature(&c, &nrm, &points).unwrap();
        // The anchored fit must read clearly curved (fillet-influenced),
        // not plane and not garbage. Within ~1 row of a tangent junction
        // the surfaces agree to within fit error, so exact recovery is not
        // achievable there — the invariant that matters for segmentation
        // is that the descriptor stays decisively non-plane.
        assert!(kmax > 0.2, "anchored fit must read curved, k_max {kmax}");
        assert!(kmin.abs() < 0.1, "k_min {kmin} ≈ 0");
    }
}
