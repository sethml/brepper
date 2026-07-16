# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

(No active items — everything below is deferred or dropped. Completed work is
recorded in DONE.md.)

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
- Verify whether the multiface_planar / large_sphere seeding masks are now redundant
  with the pass-1 junction mask (from curvature-guided seeding, see DONE.md) and
  retire them if so.
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
