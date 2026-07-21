//! Stage 2.1: sharp-edge classification + the osculating-jet curvature prior.
//!
//! Two products, both deterministic, both a single cheap pass
//! (DEVELOPMENT_PLAN §2.1). The classification fixpoint is gone: nothing
//! iterates on classification anymore.
//!
//! **(a) Sharp-edge classification — the only *decision* this stage makes.**
//! A manifold edge is BOUNDARY when its dihedral angle exceeds
//! `--angular-tolerance`; everything else is INTERIOR. This rule is
//! measured-exact on every sharp-edged analytic model. Tangent junctions
//! (G1-smooth, G2-discontinuous) are NOT classified here — the fit-free
//! curvature-jump rule that used to do that was retired (it cannot separate
//! genuine tangent junctions from the smooth curvature gradients of cones);
//! they are discovered by stage 2.3's fit competition instead.
//!
//! **(b) The jet curvature prior — data for 2.3, never a decision.** Each
//! face carries the differential summary of an **anchored osculating jet**
//! (Cazals & Pouget, SGP 2005 — the standard robust curvature estimator for
//! irregular meshes): a quadratic height field `z = Q(u,w)` fitted to vertex
//! *positions* over a PCA tangent frame of the face's own + 1-ring vertex
//! neighborhood. Three properties make it robust to how CAD models are
//! actually triangulated:
//!
//! - **It uses vertex positions** — exact f32 mesh data — not normals.
//!   A long skinny strip/sliver triangle's plane is a chord whose normal
//!   sits degrees off the true tangent, so normal-based estimators inherit
//!   that error; positions are ground truth. Tessellation-layout agnostic:
//!   aspect-ratio slivers, staggered rows, pole fans, and full-length strip
//!   meshes all fit the same.
//!
//! - **The tangent frame comes from a PCA of own + 1-ring vertices, not
//!   the facet normal** — the smallest-variance direction of the
//!   neighborhood IS the local normal and stays accurate on strips. Own
//!   vertices dominate the frame (weight OWN_WEIGHT).
//!
//! - **The face's own vertices anchor the fit** (weight OWN_WEIGHT).
//!   Tessellators never straddle a B-rep boundary with a triangle, so a
//!   face's own corners always lie on its own surface. There is **no outlier
//!   rejection** (it forces arbitrary one-side reads that vary face to
//!   face). Sharp-crease neighbors are excluded from the point set;
//!   tangent-junction mixing is allowed — junction rows read intermediate
//!   values, which is fine because no decision ever keys on them.
//!
//! The **division of labor** (Design Principle #4): the jet allocates
//! effort — seed order, rung order, initial parameters, refit triggers — and
//! never decides membership or boundaries; vertex-tolerance fit validation
//! (stage 2.3) decides all geometry. Every jet error therefore converts to
//! wasted or saved time, never to wrong output.

use super::{dot3, face_centroid};
use crate::stage1::ConnectedMesh;

/// Weight of the face's own vertices in the jet fit, relative to 1.0 for
/// ring vertices. The own corners always lie on the face's own surface
/// (tessellators do not straddle seams), so they anchor the fit.
pub const OWN_WEIGHT: f64 = 50.0;

/// Minimum unique points before the 2-ring is added to the neighborhood.
pub const JET_MIN_POINTS: usize = 12;

/// Noise-floor coefficient for the medial-point cutoff: `|κ1|` must clear
/// `FLOOR_COEFF · ε / h²` (ε the vertex-position noise scale — f32
/// quantization, bbox diagonal × 1e-7 — and h the jet neighborhood size) for
/// the curvature center to be trustworthy. Below it the face reads as flat
/// (radius → ∞) and the medial point is left undefined.
pub const FLOOR_COEFF: f64 = 4.0;

/// Classification of one manifold mesh edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EdgeClass {
    /// Both faces on the same smooth region (dihedral within tolerance).
    #[default]
    Interior,
    /// Dihedral angle exceeds the angular tolerance — a hard feature edge.
    Sharp,
}

/// Stage 2.1(a) output: per-half-edge sharp/interior classes.
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
    /// NOT a classification signal (unreliable on irregular tessellation).
    pub kappa: Vec<[f64; 4]>,
    /// Undirected manifold-edge counts by final class.
    pub sharp_count: usize,
    pub interior_count: usize,
}

impl EdgeClassification {
    /// Is the manifold edge across `faces[fi]`'s slot `e` a region boundary?
    pub fn is_boundary(&self, fi: usize, e: usize) -> bool {
        self.class[fi][e] != EdgeClass::Interior
    }
}

/// The per-face osculating-jet prior (stage 2.1(b)): the differential
/// summary consumed by stage 2.3 for effort allocation (seed order, rung
/// order, fit initialization). Never a membership or boundary decision.
#[derive(Debug, Clone, Copy)]
pub struct JetField {
    /// Larger signed principal curvature (positive convex). Few-% accurate
    /// on region interiors; contaminated on junction rows (by design — the
    /// anchor mixes). 0 where the fit was degenerate.
    pub k1: f64,
    /// Smaller signed principal curvature.
    pub k2: f64,
    /// Principal direction for `k1`, lifted to 3D in the PCA tangent frame
    /// (unit; the κ1 direction wraps a tube). The κ≈0 ruling direction is
    /// `normal × dir1`.
    pub dir1: [f64; 3],
    /// The PCA normal (unit, aligned with the facet's outward normal).
    pub normal: [f64; 3],
    /// The medial point (curvature center for the k1 direction),
    /// `centroid − n̂/k1`, defined when `|k1|` clears the noise floor. A
    /// curvature-center / tube-center estimate feeding torus initialization
    /// and bent-tube promotion.
    pub medial: Option<[f64; 3]>,
}

impl Default for JetField {
    fn default() -> Self {
        JetField {
            k1: 0.0,
            k2: 0.0,
            dir1: [0.0; 3],
            normal: [0.0; 3],
            medial: None,
        }
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

/// The full differential result of one jet fit, in absolute units.
struct JetResult {
    /// Larger signed principal curvature (positive convex).
    k1: f64,
    /// Smaller signed principal curvature.
    k2: f64,
    /// Unit 3D principal direction for k1.
    dir1: [f64; 3],
    /// Unit 3D PCA normal (aligned with the facet normal hint).
    normal: [f64; 3],
    /// Characteristic neighborhood size h (for the medial noise floor).
    h: f64,
}

/// Eigen-decomposition of the symmetric 2×2 matrix `[[a, b], [b, c]]`.
/// Returns `(l1, e1, l2)` with `l1 >= l2` and `e1` the unit 2D eigenvector
/// of the larger eigenvalue.
fn eig_sym_2x2(a: f64, b: f64, c: f64) -> (f64, [f64; 2], f64) {
    let mid = (a + c) / 2.0;
    let diff = (((a - c) / 2.0).powi(2) + b * b).sqrt();
    let l1 = mid + diff;
    let l2 = mid - diff;
    // Eigenvector for l1: solve (a-l1)x + b y = 0 and b x + (c-l1) y = 0.
    // Use whichever row is better conditioned.
    let mut v = [b, l1 - a];
    if v[0] * v[0] + v[1] * v[1] < 1e-30 {
        v = [l1 - c, b];
    }
    let n = (v[0] * v[0] + v[1] * v[1]).sqrt();
    let e1 = if n > 1e-15 { [v[0] / n, v[1] / n] } else { [1.0, 0.0] };
    (l1, e1, l2)
}

/// Fit the quadratic height field `z = a·u² + b·u·w + c·w² + d·u + e·w + g`
/// over a PCA tangent frame and return the full differential summary
/// (principal curvatures, the k1 principal direction in 3D, the PCA normal,
/// and the neighborhood scale h).
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
/// curvature instead of producing garbage). Single pass, no outlier
/// rejection: the own-vertex anchor, not rejection, is what biases each face
/// toward its own surface.
fn jet_curvature(
    origin: &[f64; 3],
    z_hint: &[f64; 3],
    points: &[JetPoint],
) -> Option<JetResult> {
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

    // Weighted normal equations + ridge.
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

    // Hessian of Q in scaled coordinates, sign-flipped so curvature is
    // POSITIVE when the surface curves AWAY from z_axis (z_axis the outward
    // face normal → convex positive; the fitted surface drops below the
    // tangent plane, so negate). The 2×2 shape matrix is [[m_a, m_b],[m_b,
    // m_c]] in the (u, w) tangent basis.
    let (m_a, m_b, m_c) = (-2.0 * p[0], -p[1], -2.0 * p[2]);
    let (l1, e1, l2) = eig_sym_2x2(m_a, m_b, m_c);
    // Rescale by 1/h (z was scaled by h as well: k = k'/h).
    let k1 = l1 / h;
    let k2 = l2 / h;
    // Lift the 2D k1 eigenvector into 3D via the tangent basis (u, w).
    let mut dir1 = [
        e1[0] * u[0] + e1[1] * w[0],
        e1[0] * u[1] + e1[1] * w[1],
        e1[0] * u[2] + e1[1] * w[2],
    ];
    let dn = (dir1[0] * dir1[0] + dir1[1] * dir1[1] + dir1[2] * dir1[2]).sqrt();
    if dn > 1e-15 {
        for c in &mut dir1 {
            *c /= dn;
        }
    }
    Some(JetResult {
        k1,
        k2,
        dir1,
        normal: z_axis,
        h,
    })
}

/// Gather the jet-fit points for one face: own vertices (anchored,
/// frame-eligible), then non-sharp ring faces' vertices. Sharp-crease
/// neighbors are excluded (their vertices lie across a hard feature edge, on
/// a different surface). 2-ring vertices are added when the neighborhood is
/// sparse (they contribute fit evidence but never steer the frame).
fn gather_jet_points(
    mesh: &ConnectedMesh,
    fi: usize,
    ec: &EdgeClassification,
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
        if ni < 0 || ec.class[fi][e] == EdgeClass::Sharp {
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
                if n2 < 0 || ec.class[ni][e2] == EdgeClass::Sharp {
                    continue;
                }
                push_face_points(n2 as usize, false, &mut points, &mut seen);
            }
        }
    }
    points
}

/// Classify every manifold mesh edge sharp/interior (stage 2.1(a)).
/// Deterministic single pass: a purely local dihedral-angle test.
pub fn classify_edges(mesh: &ConnectedMesh, angular_tolerance: f64) -> EdgeClassification {
    let n = mesh.faces.len();
    let centroids: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|f| face_centroid(f, &mesh.vertices))
        .collect();

    let mut class = vec![[EdgeClass::Interior; 4]; n];
    let mut kappa = vec![[f64::NAN; 4]; n];
    let (mut sharp_count, mut interior_count) = (0, 0);
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
            if th > angular_tolerance {
                class[fi][e] = EdgeClass::Sharp;
            }
            if fi < ni {
                match class[fi][e] {
                    EdgeClass::Sharp => sharp_count += 1,
                    EdgeClass::Interior => interior_count += 1,
                }
            }
        }
    }

    EdgeClassification {
        class,
        kappa,
        sharp_count,
        interior_count,
    }
}

/// Compute the per-face osculating-jet prior (stage 2.1(b)). Deterministic
/// single pass: one weighted jet fit per face over its own + non-sharp
/// 1-ring neighborhood. The sharp classification cordons each fit to its own
/// smooth region's side of any hard feature edge.
pub fn compute_jet_field(mesh: &ConnectedMesh, ec: &EdgeClassification) -> Vec<JetField> {
    let n = mesh.faces.len();
    let centroids: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|f| face_centroid(f, &mesh.vertices))
        .collect();

    // Vertex-position noise scale: f32 quantization (STL positions are
    // f32-exact), estimated from the bounding box diagonal. Feeds the
    // medial-point noise floor.
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

    let mut jet = vec![JetField::default(); n];
    for fi in 0..n {
        let Some(nrm) = mesh.faces[fi].normal else {
            jet[fi] = JetField::default();
            continue;
        };
        let points = gather_jet_points(mesh, fi, ec);
        let Some(res) = jet_curvature(&centroids[fi], &nrm, &points) else {
            jet[fi] = JetField {
                normal: nrm,
                ..JetField::default()
            };
            continue;
        };
        // Medial point = curvature center for the k1 direction. Defined only
        // when |k1| clears the noise floor κ_floor = FLOOR_COEFF·ε/h²
        // (below it the surface reads flat and the center runs to infinity).
        let floor = FLOOR_COEFF * eps / (res.h * res.h);
        let medial = if res.k1.abs() > floor {
            let r = 1.0 / res.k1;
            Some([
                centroids[fi][0] - res.normal[0] * r,
                centroids[fi][1] - res.normal[1] * r,
                centroids[fi][2] - res.normal[2] * r,
            ])
        } else {
            None
        };
        jet[fi] = JetField {
            k1: res.k1,
            k2: res.k2,
            dir1: res.dir1,
            normal: res.normal,
            medial,
        };
    }
    jet
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
        let res =
            jet_curvature(&[0.0, 0.0, 0.0], &[0.0, 0.0, -1.0], &points)
                .unwrap();
        assert!((res.k1 - k1).abs() < 1e-3, "k_max {} != {k1}", res.k1);
        assert!((res.k2 - k2).abs() < 1e-3, "k_min {} != {k2}", res.k2);
    }

    /// Cylinder band z = r − √(r²−x²) sampled on a skinny strip: the jet
    /// recovers (1/r, 0) regardless of aspect ratio, and the k1 principal
    /// direction points across the band (along x, not the strip's long y).
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
        let res = jet_curvature(&c, &nrm, &points).unwrap();
        // Jet approximation error on a non-quadratic circle is ~1%.
        assert!((res.k1 - 1.0 / r).abs() < 0.01, "k_max {} != {}", res.k1, 1.0 / r);
        assert!(res.k2.abs() < 0.01, "k_min {} != 0", res.k2);
        // The wrapping direction is (roughly) the x axis; the y component is
        // small because the band is straight along y.
        assert!(res.dir1[0].abs() > 0.9, "dir1 {:?} should point across the band", res.dir1);
        assert!(res.dir1[1].abs() < 0.2, "dir1 {:?} should not run along the strip", res.dir1);
    }

    /// Junction anchoring: a face on a fillet (cylinder, κ=0.5) whose
    /// neighborhood includes vertices across a tangent junction on a
    /// plane (κ=0) must follow the FILLET — the anchor biases the fit.
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
        let res = jet_curvature(&c, &nrm, &points).unwrap();
        // The anchored fit must read clearly curved (fillet-influenced),
        // not plane and not garbage. Within ~1 row of a tangent junction
        // the surfaces agree to within fit error, so exact recovery is not
        // achievable there — the invariant that matters is that the
        // descriptor stays decisively non-plane.
        assert!(res.k1 > 0.2, "anchored fit must read curved, k_max {}", res.k1);
        assert!(res.k2.abs() < 0.1, "k_min {} ≈ 0", res.k2);
    }
}
