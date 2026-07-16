# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

(No active items — everything below is deferred or dropped. Completed work is
recorded in DONE.md.)

# Deferred

## Other
- Retire `strip_merge_seed` (src/stage2/toroidal.rs) if it is truly dead: its
  trigger is `conical_hypothesis >= 0` on the start face, but since the
  toroidal-before-conical reorder no face has a conical hypothesis when the
  toroidal stage runs — pipe_elbow now recovers its torus via the direct
  seeding path. Verify with an execution probe across the suite, then delete
  (~70 lines). Discovered during the 2026-07-15 seeding-mask verification.
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
