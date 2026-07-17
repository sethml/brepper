# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

(Completed work is recorded in DONE.md.)

# Active

- Add full-pipeline tests for nosecone_truncated and hill_with_hole — both
  fail full-pipeline `--compare` volume gates with default flags but have no
  (or only stage-1) test coverage, so the suite cannot catch regressions on
  them (found 2026-07-16 while measuring the absorb pass; both failures are
  absorb-independent):
  - `tests/onshape/nosecone_truncated_medium.stl` vs
    `tests/onshape/nosecone_truncated.step`: stage 4.1 volume differs by
    1.38e-2 relative (gate 1e-2). No test file references this model at all.
  - `tests/nurbs/hill_with_hole.stl` vs `tests/nurbs/hill_with_hole.step`:
    stage 3.6 volume differs by 1.11e-2 relative (gate 1e-2). Only tested at
    stage 1 (tests/stage1_integration.rs). A 2026-07-15 measurement recorded
    ~4e-3 (inside the gate), so this may be a regression from the
    intervening stage-3 changes (commits 5f8d03a, eb7dd65: corner-refinement
    separation guard, pcurve dedupe/pointlike handling, 2.8/2.9 gating) —
    worth a bisect before adding the test.
  Diagnose each miss, fix if tractable, then add stage-4.1 compare tests for
  both (or document in the test files why a model is excluded).

- Absorb (2.8) loose rule is a deletion/tightening candidate. Its motivating
  "tangent-junction bridge facet" story was tested 2026-07-17 and is FALSE
  for the current corpus — CAD tessellators emit seam edges along B-rep face
  boundaries even where surfaces are tangent (nosecone cone–sphere tangency:
  0 straddling STL triangles, 104 vertices exactly on the tangency ring;
  probe method: classify STL vertices against the tangency plane computed
  from the fitted cone/sphere). Every beneficial merge in the corpus is the
  TIGHT rule (nosecone_medium: bogus 6-face cylinder band on the sphere;
  nosecone_truncated: 2 bogus sub-fits + a 31-face split fragment of the
  real cone). The loose rule's only current firings (stepped_block_toroidal,
  2 merges) absorb REAL narrow planar faces — walls tessellated as a single
  rectangle → one fused quad, below the ≥2-face planar hypothesis minimum —
  replacing them with neighbor geometry 0.12–0.27mm off; it also caused the
  planar-only 3DBenchy geometry errors (794 BRepCheck failures). Before
  deleting, fix the underlying gap: let a single fused-quad face form a
  legitimate planar hypothesis (pre-fusion it was two coplanar triangles,
  which met the minimum) so real single-rectangle faces stop arriving as
  orphans. Acceptance: suite green, stepped_block keeps passing --compare
  with the strips as real planar faces (it already passes with
  --skip-stage=2.8: 71 surfaces vs 69).

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

- 3DBenchy planar-only mode (`--skip-stage=2.2,2.3,2.5,2.7`) now runs
  END-TO-END (~5 min, 704MB STEP, 1 solid, volume 15393.7mm³ ≈ mesh volume;
  all 196393 mesh faces covered by 41324 planar faces, 41324/41324 OCCT faces
  created), after three fixes landed 2026-07-16 (uncommitted; suite 494 green
  + clippy clean after each): (1) `compute_edge_curve` (src/stage3.rs ~1577)
  guarded against a degenerate trim span — `Geom_TrimmedCurve` threw
  `U1 == U2` when both edge endpoints projected to the same curve parameter;
  now Errs into the mesh-boundary-chain fallback, matching
  `candidate_trim_parameters`. (2) `refine_corner_vertices` (src/stage3.rs
  ~712) rejects a refined corner position landing within
  max(vertex_tolerance, 1e-6mm) of an edge-adjacent corner's current
  position (vertex keeps the raw mesh position) — cyclic projection had
  pulled BOTH endpoints of edge 5715 onto the same near-parallel-plane
  intersection fixpoint (4.25e-14mm apart), collapsing the edge so MakeEdge
  failed `DifferentPointsOnClosedCurve` and even its straight-line fallback
  failed. NOT a mesh-weld gap: weld is grid+27-neighbor tolerance-based and
  repair creates no vertices, so distinct mesh vertices are always ≥ weld
  tolerance apart; only refinement can create sub-tolerance pairs. (3)
  `build_pcurve` (src/stage3.rs ~5454) dedupes UV interpolation samples and
  emits a degree-1 (possibly constant) 2-pole B-spline pcurve when an edge's
  projection onto a face's chart is pointlike — edge 25813 runs purely along
  face 7145's plane normal (corner refinement legitimately offset both
  endpoints ~0.02mm along x within the 2×surface_tolerance contract), so all
  UV samples coincided and `Geom2dAPI_Interpolate`'s ctor threw
  `Standard_ConstructionError` (it rejects consecutive points closer than
  its 1e-10 tolerance); the edge's 3D extent goes into the returned
  deviation → edge tolerance. NOTE these three fixes harden shared stage-3
  paths, so they may also help the normal curved-fit pipeline's stage 3.4+
  (untested here; its known blocker — face 5841, 3 winding loops — is the
  separate item above).

- **Planar-only geometry errors ROOT-CAUSED to the absorb pass (2026-07-16
  follow-up session; absorb/refine now numbered + skippable).** With absorb
  on, 794/41324 faces + the solid failed BRepCheck; the absorb pass's loose
  ≤2-face-orphan rule merges facets into neighbors within
  `--surface-tolerance` (0.4mm) chordal slop, and on benchy it merged ~128k
  of 169,890 facet surfaces — each absorbed facet then claims to lie on a
  plane up to 0.4mm away (corner refinement compounds this, moving vertices
  up to 2× surface tolerance onto sloppy-plane intersections). Change:
  absorb and boundary refinement are now CLI sub-stages **2.8** and **2.9**
  (`Stage::last_minor(2)` = 9; `--stage=2.7` now stops before them; the
  `--stage=2.6`/`2.7` integration-test helpers migrated to `--stage=2`,
  preserving their post-absorb assertions; suite 494 green, clippy clean).
  Exact planar-only mode is `--skip-stage=2.2,2.3,2.5,2.7,2.8`: BRepCheck
  face failures 794 → **4** (0.002%), the SOLID passes BRepCheck, and
  volume = 15550.53 mm³ = the official 3DBenchy 15.55 cm³ (absorb-on had
  lost 1% to slop: 15393.7). Cost: 169,890 faces, 2.3GB STEP, ~10 min.
  Stage 2.9 is harmless in this mode (196 tolerance-level moves;
  corner-refinement guard fired 0 times). Residual: the 4 remaining
  BRepCheck failures are NOT absorb-related — face 40 (planar, 13 edges,
  TWO wires, 1.32mm²) fails "invalid imbrication of wires / unorientable"
  (the multi-loop / disconnected-region family, cf. face 5841), and faces
  19767/160616/163493 are micro-slivers (0.0004–0.0017 mm², 4 edges) with
  "wire 0 invalid" — likely micro-edge / tolerance-stretched wire
  construction.

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
