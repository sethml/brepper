# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

(completed and removed: stage-3 shared-edge assembly rework (DEVELOPMENT_PLAN §3.x)
including cone_tangent_planes apex-cap face construction — volume now matches to
~1e-7 relative. See DONE.md.)

(completed and removed: stepped_block_toroidal_fillets MEDIUM — all 18 tori + 37
cylinders now selected and the model passes stage 4.1 (dedicated
onshape_stepped_block_toroidal_medium_stage41_compare test). Three fixes: (a) the BFS
final re-fit now keeps the incrementally-maintained fit when a fresh from-scratch LM
converges worse — recovered the full-ring R=3.5 pocket-rim fillet torus that refit to
0.87 mm; (b) the cylinder angular-coverage test admits coarse tessellation gaps up to
the dihedral limit (medium fillets step ~24° per row, four rows across a 90° tube arc,
which the span/3 rule rejected); (c) the default --angular-tolerance was raised from
17.5° to 25° so seed-pair qualification accepts those 24° fillet rows. The coverage
allowance is gated twice so it can't admit bogus bands: to cylinders (half_angle≈0), so
narrow corner sphere-blend patches can't win selection as 6-face cones (else
rounded_cube regresses); and to hypotheses already spanning ≥45°, so thin cylinder
bands of a smooth cone stay rejected (else nosecone_truncated over-segments). See
LEARNINGS "Medium-resolution tessellation".)

(items completed and removed: multi-shell voids — hollow_cube passes stage 4.1 with
one 7000 mm³ solid; full-torus seam coverage — full_torus (tests/occt/, R=10 r=3)
passes stage 4.1 via natural-bounds face construction; partial-cone detection and
plane-cone tangent generators — all 3 cone_tangent_planes cones detected, 11/11 edge
curves, stage 3.4 face construction issue remains, see item 1.)

(completed and removed: curvature-guided seeding stages 1+2 — per-face field with
ascending seed ordering AND pass-1 junction exclusion from gathering, plus the
toroidal-before-conical stage reorder they enabled. See DEVELOPMENT_PLAN
"Curvature-guided seed ordering". Candidate follow-up: the multiface_planar /
large_sphere masks may now be redundant with the junction mask and could be retired
after verification.)

(completed and removed: NURBS fallback (DEVELOPMENT_PLAN §2.7) — SHIPPED 2026-07-10.
Freeform regions are gathered from leftover single-face planar selections (plus
small accidental analytic fits) across smooth edges, Floater-mean-value
parameterized, resampled at mesh-vertex-aligned stations, fitted with
`GeomAPI_PointsToBSplineSurface` (interpolation first), and reconstructed through
stage 3 via a B-spline chart wrapping the exact fitted surface. All 8
tests/nurbs models pass full `--compare` (volumes 2.7e-6…5.2e-3 relative;
elliptical_cylinder via the seam-split, see below); see DEVELOPMENT_PLAN §2.7
for the as-built algorithm and its measured design rules.)

# Deferred

## NURBS follow-ups
- Non-compromised closed-band reconstruction: `elliptical_cylinder` passes via the
  seam-split (two B-spline faces on one fitted surface, two tagged artificial edges
  exempted from the stage 3.1 compare). CAD kernels and the reference write ONE
  periodic face with a seam edge instead. DEVELOPMENT_PLAN §2.7 "Path to the
  non-compromised design" lists the steps (periodic Interpolate, U-knot rescale to
  2π so the existing wrapped-face machinery applies, periodic chart flag, seam
  anchors, periodic chain repair); the acceptance test is removing the compare
  exemption.
- Periodic chain repair: `repair_edge_curve_from_boundary`'s B-spline arm covers open
  chains only; a closed freeform junction that IntSS misses would need a periodic
  interpolation (`GeomAPI_Interpolate` with the periodic flag).

## Other
- Support for OBJ, PLY input formats.
- IGES export option.

# Dropped (motivating problems solved by other means)

- **Gaussian map pre-classification**: was a performance idea from when stage 2.2
  wasted ~45 s tiling spheres with bogus cylinder bands. Fixed instead by running
  spherical detection first, per-face normal-axis consistency, and early rejection.
  No current model's runtime justifies it.
- **Unified SOR framework refactor**: its concrete motivations — pipe_elbow torus
  recovery and cylinder/cone/sphere disambiguation — were fixed by strip-merge torus
  seeding, stage ordering, and the absorption pass. A minimal (h,r) profile check may
  still be added where needed, as a targeted tool rather than a framework (the
  normal-line axis rescue for narrow cone arcs and the (h,r) profile fits in the
  torus multi-axis rescue are exactly such tools).
