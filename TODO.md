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
  classification-quality gaps (steps 2–5), not step-6 regressions.
  **Status 2026-07-25: 17 disabled → 4.** stage3_integration is fully
  re-enabled; what remains is cone_tangent_planes + cone_15x20 (stage2) and
  stepped_block_toroidal medium/fine (stage4). The stepped_block blocker
  changed character on 2026-07-25 — it is no longer a per-facet collapse but a
  decomposition-exactness gap; see its entry below.
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
    unregressed. **l_bracket FIXED 2026-07-24**: its second hole is a full
    cylindrical band whose area-weighted plane normal cancels to NaN;
    `PlanarSurface::fit_surface` returned that NaN plane and the ladder
    VACUOUSLY validated it (`dist > vtol` is false for NaN), committing a
    degenerate plane over the cylinder. Now the planar fit rejects a
    non-finite plane so the ladder falls through to the cylinder rung; both
    l_bracket holes are cylinders and it passes `--compare`. Still disabled:
    **onshape plate_with_hole (coarse)** faceted its hole into planar strips
    (coarse facet angles read as sharp) so no cylindrical region forms — a
    coarse-mesh angular-tolerance issue.
  - **fillet over-segmentation: FIXED 2026-07-24** (compete.rs). All three
    filleted models (filleted_cylinder/hole_block/pipe, 13 tests) re-enabled
    and passing `--compare`. Four competition fixes: (1) **magnitude-ordered
    principal curvatures** — a concave hole wall's dominant curvature is
    negative and lands in `k2`, so gating "is this curved?" on `|k1|` read the
    wall as flat and skipped every curved candidate (the deep hole fell to
    freeform); now `a1=max(|k1|,|k2|)`, `a2=min(...)`. (2) **seed-patch
    conditioning** (`condition_patch`): a jet-similar patch bleeds a row or
    two across a short cylinder→fillet tangent junction, pulling the fit off
    so the seed fails its own fit; trim fit-outlier patch faces until the seed
    validates. (3) **DOF-aware claim resolution**: a face several surfaces
    explain within vertex tolerance goes to the simplest (cylinder over torus;
    larger torus over a spurious thin junction ring) — absolute best-fit let a
    sub-micron sliver steal a fillet's rows and fragment it. (4) **fat-torus
    rejection** (`r ≤ 0.7·R`): an r≈R torus loosely fits a flat cap and wins
    on area (filleted_pipe).
  - **cone/cylinder** (cone_cylinder): FIXED 2026-07-24 as a bonus of the
    fillet-competition fixes (now 2 planes + 1 cylinder + 1 cone, passing).
    Still open: **cone_tangent_planes** — one cone of a tangent-cone chain
    fits as a freeform patch plus a degenerate huge-radius cone; and
    **cone_15x20** whose reference STEP is faceted (59 planes, not analytic),
    so the surface-coverage test cannot match.
  - **degenerate huge-radius analytic fits on coarse meshes: FIXED
    2026-07-24** (`radius_plausible` in compete.rs). A large-radius surface is
    locally almost flat, so a nearly-planar patch fits one within vertex
    tolerance while leaving its radius essentially unconstrained — and because
    the bogus fit floods further than the plane it displaces, it wins the area
    contest. On `part_rounded_cube` coarse (a 10 mm cube rounded on its 4
    vertical edges: 4 cylinders + 6 planes) each flat side wall grew an
    r ≈ 32.6 cylinder — nearly twice the model's 17.3 mm bounding-box diagonal
    — over the correct plane, reconstructing 0.39 mm off. A surface belonging
    to the model must fit inside it, so a radius past the bounding-box
    diagonal is degenerate by construction; this is the guard
    `torus_candidate` already applied to a torus's major radius, extended to
    cylinders/cones (via `conical_local_radius`, since a true cone's
    `radius_at_apex` is 0) and spheres. part_rounded_cube coarse now
    decomposes to the reference 4 cylinders + 6 planes and passes `--compare`
    (volume rel diff 2.4e-7); all 6 of its tests re-enabled.
    **This is a validation fix, not a tolerance fix** — relaxing
    `--vertex-tolerance` makes the error *grow* (0.39 mm at 1e-1, 1.79 mm at
    5e-1), because a looser bar admits still more degenerate fits.
    (Correction: earlier notes described these patches as "3-face corner
    blends" that should be "the r=2 corner sphere". The reference STEP has
    **no spheres at all** — only the 4 vertical edges are rounded — and the
    mis-fit patches are the flat side walls.)
  - **stepped_block_toroidal medium/fine: the per-facet collapse is FIXED
    2026-07-25; the remaining gap is decomposition *exactness*.** The whole
    model tangent-connects into ONE region (2615 faces medium, 26122 fine), so
    every surface must come out of competition. It used to leave a 30-face
    leftover blob, which tripped the `LEFTOVER_BLOB_MAX` whole-region discard
    and degraded all 26122 faces to per-facet planes. Three fixes (below) took
    it to 74 analytic surfaces on fine (reference: 71 = 37 cylinders + 18 tori
    + 16 planes) and 124 on medium; both now build real B-rep topology and
    reach stage 3.
    **The surface decomposition is now essentially exact (2026-07-26).** Fine
    commits **37 cylinders and 18 tori — the reference's counts exactly** — plus
    17 planes against its 16, for 72 surfaces against 71. The extra one is a
    single 0.253 mm² junction sliver left as a per-facet plane. This came from
    one more guard, `MIN_AXIAL_SWEEP_DEG` in compete.rs: ten surfaces were
    degenerate super-models of the surface they displaced (six near-cylindrical
    tori over perimeter fillets, four near-planar cylinders over wall strips).
    A patch that sweeps only a few degrees about its own axis never observes the
    revolution, so the radius is free — and a free radius is worse than useless
    in competition, because the richer model covers everything the simpler one
    covers and then floods one row further on the slack, taking the
    largest-area contest every time. The two populations do not overlap: real
    tori sweep 84.7°–103.4° (fitting R to the reference's exact 2.0000 / 2.5000
    / 14.0000), the impostors 3.1°–6.7° (reporting R = 20, 26, 29, 37, 87, 114);
    real cylinders sweep ≥ 64.8°, the impostors 2.8°–2.9°.

    **Still open — and now precisely diagnosed.** Both `--compare` tests still
    fail at the stage 3.1 reference-edge check (fine: edge 27 at 4.4e-2 mm;
    medium: 2.1e-1 mm). The chain, measured end to end on fine:
      1. Face 11548 is the first tessellation row of the y = -164 end cap's
         1 mm perimeter fillet. Its far vertices sit 9.5e-4 mm off the end-cap
         plane — *under* the 1.5e-3 mm flood tolerance — so the plane's grow
         swallows the row.
      2. That drags the least-squares plane to y = -163.99953, 4.7e-4 mm off
         true. The fillet cylinder is contaminated the same way (r = 1.0009).
      3. Stage 2.4 `relabel_boundaries` should hand the row back — the cylinder
         explains it better — but only by 3.2e-4 against the plane's 4.7e-4,
         inside the 2× margin that exists to stop relabel flapping. Both fits
         being contaminated by the *same* row is what closes the gap; the
         distances would separate if the loop refitted between rounds instead
         of only at the end.
      4. Edge 27 is the plane↔torus junction, which is tangent in the
         reference. Stage 3.2 does not classify it tangent, because the
         boundary sits one row into the fillet where the torus normal has
         already turned ~2.5°, so the compare applies the strict
         `--vertex-tolerance` bar instead of the chordal one it reserves for
         tangent edges.
      5. At that bar the check is unachievable in principle: a plane offset δ
         moves its intersection with a torus it should graze by √(2rδ), so
         landing the edge within 5e-5 mm would need the plane fitted to
         ~1e-9 mm.
    So the fix is upstream of the compare — place the region boundary exactly
    on the tangent circle — and the leads are (3) refit inside the relabel
    loop, and a growth tolerance that scales with facet size rather than the
    fixed `GROW_VTOL_FACTOR` × vertex tolerance (one row past a tangent
    junction departs by ~h²/2r, which shrinks quadratically as the mesh
    refines, so a fixed absolute bar is always too loose at fine tessellation).
    A tightening pass at the strict vertex tolerance was tried and does NOT
    work: 115 of 116 grows lose their own start face, because the fits
    themselves are only ~1e-4 mm accurate on full support.

    Fixes that got it here, each a real bug with a measured mechanism:
    (a) **degenerate cone on a thin support** (`axis_normals_consistent`,
    conical.rs). An r = 0.5 fillet band one quad tall between two corner tori
    leaves the axis barely constrained by the vertex residual, so the
    least-squares cone tips its axis *into the band's own plane* and beats the
    true cylinder: a bogus 12.2° cone at 3.8e-7 mm against the correct
    r = 0.5003 cylinder at 3.9e-6 mm, its axis nearly anti-parallel to the
    facet normals. `convexity_check` then rejected every face of that fit, so
    the seed grew nothing and the whole band fell through to the leftover
    blob. The facet normals are unambiguous exactly where the positions are
    degenerate (on a cone `|n·axis| = sin(half_angle)` identically), so the
    fit is now accepted only if its own support's *median* facet agrees.
    (b) **seed patch bounded in space, not just in hops** (`PATCH_RADII`,
    compete.rs). `MAX_PATCH_DEPTH` = 3 hops walks ~18 mm along a fillet chain
    swept with 6 mm facets, straddling several segments; the patch is now also
    capped at 4 radii of curvature of the seed's own jet (flats unbounded).
    Ablation: without it the medium mesh collapses back to 2615 per-facet
    planes. `condition_patch` also shrinks toward the seed when the patch
    admits no fit at all, instead of giving up with nothing to rank by.
    (c) **freeform fitted a disconnected face set** (bspline.rs). When the
    small analytic islands are *not* absorbed, stage 2.3(e) was fitting the
    leftover faces alone — which the islands may have been the only thing
    connecting. `build_region` then produced a region whose boundary loops
    bound different pieces, and `cut_region_seam` failed with `closed band: no
    path between boundary loops` (visible on stepped_block, part_rounded_cube
    and bspline_wave). Each leftover piece is now fitted separately; and when
    the islands split the leftovers into pieces that are each patch-sized, the
    islands are the accident and get absorbed. The freeform domain also now
    admits only *enclosed* islands — one bordering a surface competition
    explained at scale is a real small feature, not an accident adrift in
    freeform — which is what stops a densely-filleted model pulling its whole
    decomposition into one nonsensical freeform domain.
  - **stale blockers pruned 2026-07-24**: 7 disabled tests were re-verified
    and already passed (fixed incidentally by the hole→sphere, fillet
    competition, and NaN-plane fixes) — onshape `plate_with_hole` coarse
    (params_match/stage22/stage23), `rounded_cube_tangent_edges`, and the 3
    ccad stage-4 `--compare` tests (cone_cylinder, filleted_cylinder,
    filleted_pipe), whose "reference geometry mismatch" no longer reproduces.
    Re-enabled. Lesson: re-run `-- --ignored` before assuming a blocker is
    still live.

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
