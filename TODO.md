# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

(Completed work is recorded in DONE.md.)

# Active

- 3DBenchy stage 3.4 abort: `face 5841: 3 winding boundary loops (at most 2
  supported)` (scan_face, src/stage3.rs). Repro:
  `./target/release/brepper tests/downloaded/3DBenchy.stl -o /tmp/benchy.step`
  (~25 min; stages 1-3.3 now complete cleanly — 68924/68924 edge curves, 0
  failed — after the 2026-07-16 crash-fix session recorded in DONE.md).
  Diagnosis so far (2026-07-16 exploration, no fix yet): face 5841 is a
  fitted cylinder r=0.2048mm, axis [0.8356,-0.1769,0.5200] through
  [-3.0714,17.9168,11.3308]. Its ordinary boundary edges cluster near
  [10.7,15.1,19.7] (axial station +16.4), but its two closed-loop edges —
  25777 (f5841,f9597) and 25778 (f5841,f9598) — sit at stations +28.5 and
  -14.7, all three ON the cylinder (radial 0.19-0.205mm), 43mm apart along
  one line. Faces 9597/9598 are single-FACET islands (1 edge, 1 neighbor,
  0 vertices — selection step-4 single-face-planar fallbacks) surrounded
  by region 5841. Their edge curves are full 2π IntSS circles (param spans
  6.2832) that passed chain validation vacuously: the feature's whole
  diameter (0.41mm) is about equal to surface_tolerance (0.4mm), so ANY
  small chain near a station is "within tolerance" of the full ring. Two
  candidate root causes, possibly combined: (a) genuine disconnection —
  BFS regions are connected at deduction, but stage 2.6 greedy selection
  keeps only a winner's still-unassigned faces (select_surfaces,
  src/stage2/mod.rs ~line 700) and overwrites hyp.faces with that remnant,
  with NO connectivity recheck (refine_boundary_assignments can likewise
  sever 1-face-wide strips; absorb preserves adjacency per merge but
  propagates existing disconnection); (b) sub-tolerance feature scale —
  fit-based growth and chain validation are both untethered when the
  feature is smaller than surface_tolerance, so contractible island
  junctions got full-circle curves that classify as winding in UV. Next
  step: add a debug dump that BFSes each selected hypothesis's final face
  list to report disconnected components (distinguishes (a) from (b)),
  then either split hypotheses into mesh-connected components before face
  construction and/or gate fits below a minimum feature size. Extending
  the seam planner to N winding loops is the wrong fix (≥3 winding loops
  cannot bound a connected region). Note stage 3.4+ is uncharted for this
  model — expect further issues behind this one. Also
  from that session, quality follow-ups worth revisiting: 48 micro-edges
  whose trimmed curves collapse to a point fall back to straight vertex-to-
  vertex edges (MakeEdge `DifferentPointsOnClosedCurve`), and 684 of 68924
  edge curves come from mesh-boundary-chain fits because surface-surface
  intersection failed.

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
