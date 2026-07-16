# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

(Completed work is recorded in DONE.md.)

# Active

- Fix the stage 3 crash on 3DBenchy: OCCT panic
  `Standard_ConstructionError: Geom_TrimmedCurve::parameters out of range`
  after ~2027 edges are built. Repro:
  `./target/release/brepper tests/downloaded/3DBenchy.stl -o /tmp/benchy.step`
  (a few minutes; stage 1 mesh repair and all of stage 2 complete fine; run
  with -v to watch edge progress). Hints from the 2026-07-16 run: the panic
  is raised in opencascade-sys's generated geom.rs from a Geom_TrimmedCurve
  constructor, i.e. stage 3 trims an edge curve to endpoint parameters
  outside the curve's admissible domain. Immediately before the crash the -v
  log shows (a) edges with large negative parameter ranges (e.g.
  `t=[-21.8688, -21.4696]`), which are fine for infinite analytic curves but
  out of range for a bounded B-spline, and (b) a boundary-edge repair that
  swaps in a different curve: `Edge (f289,f4913): curve was 0.402 mm from
  its mesh boundary chain; repaired with the chain's B-spline (4.37e-2 mm)`.
  Prime suspects: trimming a repaired/replacement curve with parameters that
  were computed by projecting onto the *original* curve, or failing to
  normalize t0/t1 into a periodic curve's period. Find the stage3 call sites
  constructing trimmed curves, re-project endpoint parameters onto the curve
  actually being trimmed (clamping/normalizing to its domain), and extract a
  small regression repro if possible.

# Deferred (until a model needs it)

- Boundary crack stitching (mesh repair): close hairline cracks by pairing
  nearly-coincident open boundary edges and re-welding their endpoints at a
  progressively larger tolerance (bounded by --vertex-tolerance) — what
  netfabb "stitch triangles" and CGAL `stitch_borders` do. Common in CAD
  exports where facets along a shared edge disagree at ~1e-6. No current
  test model exhibits this; add when one does. Implementation slot:
  `ConnectedMesh::repair()` in src/stage1.rs, between the T-junction pass
  (`repair_degenerate_faces`) and `repair_coincident_faces`.

# Deferred

## Other
- Support for OBJ, PLY input formats.
- IGES export option.
- `EdgeHighlight` struct and `VizOverlay.edge_highlights` field (src/viz.rs)
  are never populated by any producer (unlike the sibling overlay types
  `FaceHighlight`, `CylinderOverlay`, `ConeOverlay`, `SphereOverlay`,
  `TorusOverlay`, `LineOverlay`, `ShapeMeshOverlay`, all of which have real
  producers in stage2/stage3). The consuming render path
  (`run_viz_window`/render loop) always iterates an empty vec. Investigate
  whether this was meant to highlight boundary edges and never got wired up,
  or is pure vestige — then either wire it up or delete it (struct + field +
  init + consume + render, ~20 lines across viz.rs). Discovered during the
  2026-07-15 dead-code audit.
- `error_abs_sum` is accumulated but never read on `PlanarHypothesis`,
  `ConicalHypothesis`, `SphericalHypothesis`, `ToroidalHypothesis`
  (src/stage2/mod.rs) and on `BSplineHypothesis` (src/stage2/bspline.rs,
  where it's also excluded from the manual `Debug` impl). Unlike sibling
  fields (`centroid_error_max`, `error_max`, `error_min`), which are read in
  verbosity diagnostics or test assertions, this field has no reader
  anywhere. Either wire it into a diagnostic (e.g. mean-absolute-error
  alongside the existing max-error prints) or delete the field and its
  accumulation code across all five structs. Discovered during the
  2026-07-15 dead-code audit.
- `--debug` CLI flag (`Config.debug`, src/config.rs) is parsed but nothing
  branches on it — the flag's doc comment promises "debug output and
  intermediate files" but no code path currently delivers that. Either
  implement the promised behavior or remove the flag. Discovered during the
  2026-07-15 dead-code audit.
- `--step-units` CLI flag (`Config.step_units`, src/config.rs) is parsed but
  never consulted — the STEP writer doesn't scale output by it. Either wire
  it into `write_step_file` or remove the flag. Discovered during the
  2026-07-15 dead-code audit.
- `impl Display for Units` (src/config.rs) is never invoked (only the derived
  `Debug` is used, e.g. in `Config`'s custom `Debug` impl). Delete it unless
  something is meant to print units in a user-facing way. Discovered during
  the 2026-07-15 dead-code audit.

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
