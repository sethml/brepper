# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

(No active items — everything below is deferred or dropped. Completed work is
recorded in DONE.md.)

# Deferred

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
