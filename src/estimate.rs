//! Tolerance estimators (Phase 0.7): read the exporter's settings back out of
//! the mesh it produced.
//!
//! Two of brepper's three tolerances are guesses. `--surface-tolerance` defaults
//! to 0.4 mm against a measured true sag around 1e-3 mm, and
//! `--angular-tolerance` stands in for an exporter setting discarded when the
//! STL was written. Both leave evidence in the dihedral angles, because a CAD
//! tessellator subdivides until a chordal bound and an angular bound are both
//! met.
//!
//! **Diagnostics only.** Nothing downstream reads these. Phase 3 decides whether
//! to flip the defaults, on corpus measurements rather than on this printout.
//!
//! # What is reported, and what is deliberately not
//!
//! The **sagitta** estimate is fit-free: for a chord of length `h` across an arc
//! turning by `θ`, the radius is `R ≈ h/θ` and the sagitta is
//! `R(1 − cos(θ/2)) ≈ Rθ²/8 = hθ/8`. So `sag ≈ dihedral × facet size / 8`, with
//! the dihedral from adjacent facet normals and the facet size from the distance
//! between their centroids.
//!
//! Both quantities are reported as QUANTILES OVER EVERY INTERIOR EDGE, with no
//! attempt to separate smooth edges from sharp ones. That separation was tried
//! and does not survive contact with the corpus:
//!
//! * Otsu's method splits a bimodal distribution, but this one is trimodal —
//!   after quad fusion a flat region still contributes ~0° edges, a curved
//!   surface contributes a cluster at the angular setting, and real edges sit
//!   near 90°. On `rounded_cube_10_r2_fine` Otsu split at 1.75°, calling 15534
//!   of 21586 edges sharp on a model that HAS NO SHARP EDGES.
//! * Widest-empty-band fails too: on `stepped_block_toroidal_fillets_fine` the
//!   true band is 3°–5° with only 442 of 76443 edges above it, but the sparse
//!   sharp tail from 5° to 103° contains wider empty runs.
//! * Every remaining formulation needed one arbitrary ratio, which is the
//!   "detector needs just one threshold" case TODO.md rule 5 says to stop on.
//!
//! Quantiles sidestep it — but only on a curvature-dominated mesh, and Phase 0.8
//! showed that is a property of our corpus rather than of CAD models. On the
//! synthetic corpus the dihedral distribution is a SINGLE SPIKE at the
//! exporter's angular setting (rounded_cube fine reads p50 2.49° / p95 2.50°;
//! sphere_25 fine 2.50° / 2.50°), because those models are spheres, rounded
//! cubes and fillets by construction. Real print models are mostly flat panels
//! with sharp rims: 3DBenchy reads p50 0.79° / p95 89.21°, OrcaCube 1.13° /
//! 90.00°, xyz-cube 3.00° / 90.00° — **the p95 is a model rim on all three**,
//! and the sag p95 with it (2.0e-2 to 6.2e-1 mm, implausible on inspection).
//!
//! No fixed quantile fixes this: which one lands on the curved mode depends on
//! each model's flat:curved:sharp proportions, so picking one is the "detector
//! needs just one threshold" case again. **So this stays a diagnostic and
//! nothing derives a tolerance from it.** What makes that safe is that the
//! dihedral quantiles are printed beside the sag and announce the failure
//! directly — on a polyhedron like `chamfered_cube_10_c1`, or on any of the
//! three real models above, the median or p95 dihedral is itself a model edge
//! and no sag estimate from that mesh means anything.
//!
//! Measured accuracy against the true sag that stage 1 `--compare` reports, over
//! 16 SYNTHETIC corpus models: p95 lands within 5× on 13. The two failure modes
//! are distinct, and the dihedral quantiles printed beside the sag announce
//! which one a given mesh is in:
//!
//! * **2× low on doubly-curved geometry** — sphere 0.50×, hemisphere 0.43×,
//!   rounded_cube 0.50×, tori 0.49×. `hθ/8` is the sagitta of a curve bending
//!   ONE way; a sphere bends two ways and both curvatures displace the facet
//!   centre. Adding the second term would be an exact correction — and is
//!   deliberately NOT done, because it would sharpen a number that is already
//!   measuring the wrong edges on every real model tested above.
//! * **Wildly high where sharp edges are a large fraction of all edges**, since
//!   the p95 is then a model rim rather than a curvature step: `simple_cylinder`
//!   has a dihedral p50 of 90° (a short cylinder is mostly cap-to-wall rim) and
//!   reads 62×; `plate_with_hole` coarse reads 151×.
//!
//! So read the dihedral quantiles first. 2.49/2.50/2.50 means every edge is a
//! curvature step and the sag estimate is meaningful; 12.41/90.00/90.00 means it
//! is not.

use crate::stage1::ConnectedMesh;

#[derive(Debug, Clone)]
pub struct ToleranceEstimates {
    /// Dihedral angle between adjacent facets, degrees. The p95 approximates the
    /// exporter's angular setting on a model whose edges are mostly smooth;
    /// a large p50 means the mesh is mostly creases, so neither estimate here
    /// describes a tessellated surface.
    pub dihedral_p50_deg: f64,
    pub dihedral_p95_deg: f64,
    pub dihedral_max_deg: f64,
    /// Sagitta, `dihedral × facet size / 8`, in mm. Read p95.
    pub sag_p50: f64,
    pub sag_p95: f64,
    pub sag_max: f64,
    pub interior_edges: usize,
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    sorted[(((sorted.len() - 1) as f64) * q).round() as usize]
}

fn centroid(mesh: &ConnectedMesh, face_idx: usize) -> [f64; 3] {
    let face = &mesh.faces[face_idx];
    let n = face.vertex_count as usize;
    let mut c = [0.0_f64; 3];
    for k in 0..n {
        let v = &mesh.vertices[face.vertex_indices[k]];
        c[0] += v.x;
        c[1] += v.y;
        c[2] += v.z;
    }
    let inv = 1.0 / n as f64;
    [c[0] * inv, c[1] * inv, c[2] * inv]
}

/// Estimate the exporter's angular and chordal settings from the mesh alone.
pub fn estimate(mesh: &ConnectedMesh) -> ToleranceEstimates {
    let mut dihedrals: Vec<f64> = Vec::new();
    let mut sags: Vec<f64> = Vec::new();
    for (fi, face) in mesh.faces.iter().enumerate() {
        let Some(n1) = face.normal else { continue };
        for k in 0..face.vertex_count as usize {
            let nb = face.neighbors[k];
            if nb < 0 {
                continue;
            }
            let fj = nb as usize;
            if fj <= fi {
                continue; // each undirected edge once
            }
            let Some(n2) = mesh.faces[fj].normal else { continue };
            let dot = (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]).clamp(-1.0, 1.0);
            let dihedral = dot.acos();
            let (c1, c2) = (centroid(mesh, fi), centroid(mesh, fj));
            let gap = ((c1[0] - c2[0]).powi(2) + (c1[1] - c2[1]).powi(2) + (c1[2] - c2[2]).powi(2))
                .sqrt();
            dihedrals.push(dihedral.to_degrees());
            sags.push(dihedral * gap / 8.0);
        }
    }
    dihedrals.sort_by(f64::total_cmp);
    sags.sort_by(f64::total_cmp);
    ToleranceEstimates {
        dihedral_p50_deg: quantile(&dihedrals, 0.50),
        dihedral_p95_deg: quantile(&dihedrals, 0.95),
        dihedral_max_deg: quantile(&dihedrals, 1.0),
        sag_p50: quantile(&sags, 0.50),
        sag_p95: quantile(&sags, 0.95),
        sag_max: quantile(&sags, 1.0),
        interior_edges: dihedrals.len(),
    }
}
