//! Stage 2.1: fit-free edge classification.
//!
//! Labels every manifold mesh edge INTERIOR or BOUNDARY using two fit-free
//! signals (DEVELOPMENT_PLAN §2.1): a sharp dihedral angle, or a jump in the
//! per-face two-sided discrete-curvature descriptor (σ_max, σ_min) across
//! the edge. The true segmentation's boundaries are a subset of mesh edges
//! (tessellators emit seam edges along B-rep face boundaries, even at
//! tangencies), so this classification IS the segmentation; stage 2.2
//! extracts its connected regions and stage 2.4 (future) repairs mistakes.
//!
//! The curvature-jump rule runs iteratively: descriptors are computed over a
//! face's NON-boundary edges only, then edges are (re)classified, until the
//! classification is stable. The iteration exists to kill "ghost"
//! boundaries: the plane face bordering a tangent fillet inherits the
//! junction edge's transition curvature into its descriptor, which on round
//! one makes it differ from its plane neighbors one row in. Once the
//! junction edge is classified boundary it leaves the descriptor, the
//! plane face's descriptor returns to ≈0, and the ghost un-fires (the
//! junction edge itself keeps firing — the fillet side's descriptor is
//! supported by its interior cross-tube edges).

use super::{dot3, face_centroid};
use crate::stage1::ConnectedMesh;

/// Local-relative curvature-jump threshold J (DEVELOPMENT_PLAN §2.1): an
/// edge is a boundary when the descriptor jump exceeds this fraction of the
/// local curvature magnitude. Carried over from the junction-mask finding:
/// thresholds relative to a mesh-global maximum fail because sharp creases
/// dominate it while tangent junctions fall below any global fraction.
pub const JUMP_RATIO: f64 = 0.5;

/// Noise-floor coefficient: κ_floor = FLOOR_COEFF · vertex_tolerance / d²
/// per edge (centroid distance d). Vertex position noise ε perturbs the
/// dihedral by ~ε/d and the discrete curvature by ~ε/d²; the floor keeps
/// exact-coplanar plane interiors (κ ≡ 0 ± noise) from ever firing the
/// relative test. Calibration constant — revisit against the corpus.
pub const FLOOR_COEFF: f64 = 4.0;

/// Cap on classify/re-descriptor rounds. Convergence is typically 2–3
/// rounds (fire junctions, ghosts un-fire); the cap guards oscillation.
pub const MAX_ROUNDS: usize = 10;

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
    /// Signed discrete normal curvature κ = ±θ/d per half-edge (positive
    /// convex, negative concave); NaN where no manifold neighbor.
    pub kappa: Vec<[f64; 4]>,
    /// Final per-face curvature descriptor: max signed κ over the face's
    /// non-boundary manifold edges (0.0 when it has none).
    pub sigma_max: Vec<f64>,
    /// Final per-face descriptor minimum (see `sigma_max`).
    pub sigma_min: Vec<f64>,
    /// Rounds until the jump classification stabilized (or MAX_ROUNDS).
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

/// Classify every manifold mesh edge (stage 2.1). Deterministic, O(E) per
/// round with a small bounded round count.
pub fn classify_edges(
    mesh: &ConnectedMesh,
    vertex_tolerance: f64,
    angular_tolerance: f64,
) -> EdgeClassification {
    let n = mesh.faces.len();
    let centroids: Vec<[f64; 3]> = mesh
        .faces
        .iter()
        .map(|f| face_centroid(f, &mesh.vertices))
        .collect();

    // Per half-edge: signed κ, sharpness, and the noise floor. Each side is
    // computed independently from symmetric formulas, so mirrored slots get
    // identical values without explicit pairing.
    let mut kappa = vec![[f64::NAN; 4]; n];
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
            let theta = dot3(&n1, &n2).clamp(-1.0, 1.0).acos();
            let dx = centroids[ni][0] - centroids[fi][0];
            let dy = centroids[ni][1] - centroids[fi][1];
            let dz = centroids[ni][2] - centroids[fi][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            if d <= 1e-12 {
                continue;
            }
            // Sign: convex (+) when the neighbor's centroid lies below this
            // face's plane. Symmetric — for a convex edge each centroid is
            // below the other face's plane.
            let below = n1[0] * dx + n1[1] * dy + n1[2] * dz <= 0.0;
            kappa[fi][e] = if below { theta / d } else { -theta / d };
            sharp[fi][e] = theta > angular_tolerance;
            floor[fi][e] = FLOOR_COEFF * vertex_tolerance / (d * d);
        }
    }

    // Iterate: descriptors over non-boundary edges → (re)classify jumps.
    let mut jump = vec![[false; 4]; n];
    let mut sigma_max = vec![0.0_f64; n];
    let mut sigma_min = vec![0.0_f64; n];
    let mut rounds = 0;
    for _round in 0..MAX_ROUNDS {
        rounds += 1;
        for fi in 0..n {
            let vc = mesh.faces[fi].vertex_count as usize;
            let mut smax = f64::NEG_INFINITY;
            let mut smin = f64::INFINITY;
            for e in 0..vc {
                let k = kappa[fi][e];
                if k.is_finite() && !sharp[fi][e] && !jump[fi][e] {
                    smax = smax.max(k);
                    smin = smin.min(k);
                }
            }
            if smax == f64::NEG_INFINITY {
                smax = 0.0;
                smin = 0.0;
            }
            sigma_max[fi] = smax;
            sigma_min[fi] = smin;
        }

        let mut changed = false;
        let mut new_jump = vec![[false; 4]; n];
        for fi in 0..n {
            let vc = mesh.faces[fi].vertex_count as usize;
            for e in 0..vc {
                if !kappa[fi][e].is_finite() || sharp[fi][e] {
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
                let fire = jmp > JUMP_RATIO * scale && jmp > floor[fi][e];
                new_jump[fi][e] = fire;
                if fire != jump[fi][e] {
                    changed = true;
                }
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
            if ni < 0 || !kappa[fi][e].is_finite() {
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
