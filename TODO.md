# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

(Completed work is recorded in DONE.md.)

# Active

- **stage-2 rewrite (partition-then-classify) is COMPLETE.** All 7 roadmap
  steps landed: the old BFS pipeline is deleted, the analytic fitters are
  driven by the whole-region ladder + scoped competition + freeform rung via
  the trimmed `surface::SurfaceType` trait, and the suite is green. See the
  **Rewrite Roadmap** in DEVELOPMENT_PLAN.md and the **Stage-2 Retrospective**
  for rationale. Remaining classification-quality gaps (disabled test
  blockers) are listed in the item below; the obsolete old-BFS-pipeline TODO
  items (absorb pass, `--skip-stage`, planar-only mode) were pruned 2026-07-23.

- **Re-enable the models disabled in step 7a** (each `#[ignore]`/`ignore =`
  in tests/stage2_integration.rs, stage3_integration.rs, stage4_integration.rs
  carries its blocker string). These are PRE-EXISTING new-pipeline
  classification-quality gaps (steps 2–5), not step-6 regressions:
  - **hole → sphere: FIXED 2026-07-23** (jet-informed rung ordering in
    classify.rs). A short cylinder tessellated as a few rings is
    vertex-ambiguous with a sphere (two same-radius rings at heights ±h lie
    exactly on a sphere of radius √(R²+h²)), so simplest-DOF-first wrongly
    committed the sphere. `region_is_developable` now uses the jet's
    `min/max` principal-curvature ratio (~0 for cylinder/cone, ~1 for sphere)
    to try the conical rung before the sphere rung for developable regions;
    the jet only orders, so a genuine sphere the jet misreads still falls
    through. Fusion `plate_with_hole` (9 tests) re-enabled and passing
    (`--compare` volume rel diff 1.1e-9, was 10.8%); genuine spheres
    unregressed. Two related tests stay disabled for DIFFERENT gaps: (a)
    **l_bracket** now classes its hole as a cylinder but only 1 of the
    reference's 2 cylinders is found (a second cylindrical feature classifies
    as planar/split); (b) **onshape plate_with_hole (coarse)** faceted its
    hole into planar strips (coarse facet angles read as sharp) so no
    cylindrical region forms — a coarse-mesh angular-tolerance issue.
  - **fillet over-segmentation** (filleted_cylinder/hole_block/pipe): a
    fillet chain decomposes into extra tori + a freeform patch and the
    straight cylinder/cone body is absorbed rather than kept as one primitive
    — compete.rs growth/promotion.
  - **cone/cylinder + concave-fillet** (cone_cylinder, cone_tangent_planes):
    decomposition differs from the old BFS.
  - **coarse-mesh leftover blob** (part_rounded_cube, stepped_block_toroidal
    coarse/medium): the undecomposed fillet blob exceeds `LEFTOVER_BLOB_MAX`
    so the whole region degrades to per-facet planes (documented in
    classify.rs). Freeform (step 6) does not rescue it — the blob is not a
    clean freeform patch.
  - **ccad reference mismatches** (ccad cone_cylinder/filleted_* stage-4
    `--compare`): the generated reference STEP disagrees with the
    reconstruction at ~mm scale; environmental / reference-quality.

- **FIXED (2026-07-22, OCCT 8.0.0_p1 upgrade): partial-torus STEP round-trip
  volume loss (`pipe_elbow_10_fine`).** Under OCCT 7.9.3 the written STEP for
  the fine pipe-elbow re-read 14% low (3087 vs 3590 mm³) because OCCT's
  pcurve computation for the partial (bend-segment) torus degenerated on
  both the write and read paths — a known OCCT bug (the 8.0.0 notes list
  "partial torus creation handles inverted V range"). The OCCT 8.0 port
  (commit 4b20904) fixes it with no brepper geometry change: round-trip rel
  diff is now 7.6e-9, and `pipe_elbow_10_fine --compare` is a tracked
  regression guard (commit 71a2735). The `write.surfacecurve.mode=0`
  workaround did NOT help and was abandoned. Full diagnosis is in the
  commit history / memory if the failure ever recurs.

- Add full-pipeline tests for nosecone_truncated and hill_with_hole — both
  fail full-pipeline `--compare` volume gates with default flags but have no
  (or only stage-1) test coverage, so the suite cannot catch regressions on
  them (re-measure on the new pipeline first — these figures predate the
  rewrite):
  - `tests/onshape/nosecone_truncated_medium.stl` vs
    `tests/onshape/nosecone_truncated.step`: stage 4.1 volume differed by
    ~1.4e-2 relative (gate 1e-2). No test file references this model at all.
  - `tests/nurbs/hill_with_hole.stl` vs `tests/nurbs/hill_with_hole.step`:
    stage 3.6 volume differed by ~1.1e-2 relative (gate 1e-2). Only tested at
    stage 1 (tests/stage1_integration.rs).
  Diagnose each miss on the current pipeline, fix if tractable, then add
  stage-4.1 compare tests for both (or `#[ignore]` with a blocker comment,
  per the migration policy).

- **3DBenchy full-pipeline is UNCHARTED on the new pipeline.** All prior
  3DBenchy diagnoses (the face-5841 "3 winding boundary loops" abort, the
  planar-only `--skip-stage` mode, the absorb-pass BRepCheck failures) were
  OLD-BFS-pipeline-specific and are obsolete: stage-2.6 greedy selection,
  the absorb pass, and `--skip-stage` are all deleted, and the new stage-2.4
  repair enforces per-surface mesh-connectivity as a hard invariant
  (`resplit_disconnected` + `assert_invariants`) — which was a leading root
  cause of the winding-loop abort (a selected hypothesis with disconnected
  faces). Re-run `brepper tests/downloaded/3DBenchy.stl -o out.step` and
  re-diagnose from scratch when the large-real-model stress case becomes a
  priority. (The three shared stage-3 crash fixes from 2026-07-16 —
  degenerate trim span, corner-refinement separation guard, pointlike-pcurve
  dedupe in stage3.rs — are committed and remain valid.)

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
