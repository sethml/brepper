# Brepper TODO — the Robustness Rework (plan set 2026-07-29)

This file is the **step-by-step execution plan**. Work the phases in order;
within a phase, work the numbered steps in order unless a step is marked
**[parallel-ok]**. The rationale, design detail, corpus strategy, and the
choices deliberately not made live in DEVELOPMENT_PLAN.md under **"The
Robustness Rework (2026-07-29)"** — read that section AND "Rewriting stage
2.3(b): the blind patch is the root cause" before starting any phase.
Completed work is recorded in DONE.md. When you finish a step, annotate it
here in place (`**DONE <date>**` plus a one-line result with the key numbers),
matching the project's existing convention.

**Goal.** Maximum robustness and minimum fragility over all plausible
CAD-derived meshes — not score on the current corpus. Fragility is
operationalized: the stage-2 tuned-constant count (**34** at plan time) and
the **cliff count** (constants whose ±2× perturbation collapses any model).
Report both in every commit message that moves either.

---

## Execution rules — apply to every commit

1. **Invariants after every commit**: every mesh face in exactly one surface
   (`assert_invariants`); committed supports are flood/validation outputs at
   `vtol`/`stol`, never hand-adjusted; build clean; `cargo clippy` clean;
   suite green (`cargo test --release --no-fail-fast`, ~2.5 min) or newly
   `#[ignore]`d only with a precise blocker string.
2. **Measurement protocol**: save a `tmp/probe.sh` baseline before each phase;
   diff **per-type breakdowns** (plane/cyl/cone/sphere/torus/bspline), never
   surface counts alone (a 3-torus regression once hid behind an unchanged
   total of 71). The full suite is the decider. Once Phase 0 lands, the sweep
   harness joins the protocol.
3. **No new tolerance-like or threshold-like constant without user sign-off.**
   Numeric-conditioning floors (determinant floors, minimum |n₁×n₂| before
   normalizing) are not thresholds — say so in their doc comments when adding.
4. **Nothing is deleted without ablation evidence**: disable → probe → suite
   (→ harness) → delete, recording in DEVELOPMENT_PLAN what it was for and
   what superseded it.
5. **Stop and ask the user**, with measurements written down, when:
   - a step's listed assumption check fails;
   - results become **ordering-dependent** (ordering changes output);
   - you are about to implement a freeform-candidacy or model-selection rule
     (five have failed; none lands unmeasured);
   - a detector seems to need "just one threshold";
   - a guard proposed for deletion still fires anywhere;
   - you want to start anything in the **Reserve** section.
6. **Dead ends — do not re-derive** (details in DEVELOPMENT_PLAN):
   curvature/dihedral tangent-seam *detection* (dead twice; a tangent junction
   is locally the smoothest place on the mesh); acceptance-tolerance scaling
   with mesh quantum (chaotic, no plateau); complexity pricing among nested
   analytic types; clustering/snapping fitted radii; per-seed borrowed-axis
   commits without global selection (fragments surfaces, measured twice);
   face-size-ratio boundary tests (anti-correlated); count/fraction-based
   freeform candidacy rules (four dead); retuning flood/conditioning floors
   (non-monotonic everywhere).

**Diagnostics toolbox**: `tmp/probe.sh` (14-model per-type breakdown);
`dump_regions <model> jetdump` (per-face k1/k2/area/adjacency/normal);
`dump_regions <model> locus [a b]` (junction loci); `BREPPER_DEBUG_COMPETE=1`
(per-seed compete log); `src/bin/step_measure_tolerance.rs`.

---

## Phase 0 — Make robustness measurable

No behavioral changes in this phase (except the tiny 0.6). Build the
instrument first; you cannot claim generality you cannot measure.

- **0.1 Coverage lock.** **DONE 2026-07-29** — audited all 65 STLs under
  `tests/` (excluding `bad/`); 19 had a reference STEP and no coverage past
  stage 1 or 2. All 19 pass end-to-end `--compare` today, so all 19 are now
  non-ignored stage-4 tests (suite 513/9 → 533/9). `nosecone_truncated`
  re-measured at **1.08e-7** relative volume, not the recorded 1.4e-2 — that
  number was stale. Two data defects recorded in comments, not chased:
  `cone_15x20.step` is a faceted reference (59 planes; its volume is the
  inscribed polyhedron's, 1.96e-3 under true πr²h/3, which brepper hits
  exactly) and `cone_tangent_planes_medium` passes at 5.2e-3 against the 1%
  gate where every other model is ≤1.7e-7 — the Phase-5c canary, now pinned.
  Orphan `chamfered_cube_10_c1_medium.ascii.stl` pinned by an ASCII-vs-binary
  stage-1 equivalence test instead of a duplicate pipeline run.
- **0.2 OCCT re-mesh utility.** Extend `src/bin/gen_test_models.rs` (it
  already drives `BRepMesh` via `write_tessellated_stl(shape, dir, name,
  deflection, angle_deg)`) with a mode that loads an existing reference STEP
  and emits STL variants at several (deflection, angle) settings, plus
  **off-origin variants**: translate the shape ~1 m, write both the STL and a
  translated reference STEP so `--compare` still works. Wire a handful into
  the suite. This adds constrained-Delaunay provenance and makes the float32
  quantum testable.
- **0.3 Sweep harness.** A tracked script (`tools/sweep.sh`; check
  `git ls-files tmp/probe.sh` first — if probe.sh is untracked, move both into
  `tools/` so the methodology is versioned) that runs every corpus + generated
  model through the full pipeline with `--compare` and emits one
  machine-readable row per model: per-type surface breakdown, stage 3.1 edge
  distance, stage 4.1 volume error, runtime, exit status. Supports saving a
  baseline and diffing against it. Extend `gen_test_models` with seeded
  randomized variants (radii, fillet sizes, tangent cone chains, tori across
  r/R, densities) so the model population is not memorizable.
  **Assumption check**: the harness reproduces today's known-good results on
  the existing models; if not, the harness is wrong — fix it before anything
  else.
- **0.4 Fragility meter.** A script that, one constant at a time (list them
  via the grep in the header of DEVELOPMENT_PLAN's rework section), perturbs
  ±2× or disables, rebuilds, runs the sweep, and reports which models changed.
  Record the baseline **cliff count**. This makes "minimize fragility"
  a number.
- **0.5 Real-pair corpus. [parallel-ok, user-involved]**
  - (a) **Fusion re-export**: ask the user to batch-import reference STEPs
    into Fusion 360 and export STLs at low/medium/high; drop under
    `tests/fusion_reexport/`.
  - (b) **Wild pairs**: harvest 10–20 models published with BOTH a
    STEP (or F3D the user can convert) and an STL — Thingiverse API (token
    from user), Printables, GrabCAD (browser session if no API). Keep files
    in `tests/external/` (gitignored); commit a manifest (URL, author,
    license, sha256). Do not commit model files unless the license permits
    redistribution.
  - (c) **Pair-verification gate**: for each pair, run stage 1 `--compare`;
    the pair counts as ground truth only if mesh vertices lie on the STEP
    surfaces within a few quanta. Otherwise quarantine as crash-test-only
    (different revision / decimated / repaired). **Record the verification
    rate** — it is the first field measurement of the provenance premise in
    the wild. **Stop and ask** if most wild pairs fail it.
- **0.6 Quantum clamp + report** (small, safe): compute the mesh's float32
  quantum (`max|coordinate| × 1.19e-7`) in stage 1; if `--vertex-tolerance`
  is below it, clamp and warn; print the quantum in the stage-1 summary.
  Honor the recorded DO-NOTs: no acceptance-tolerance scaling (measured
  chaotic); the quantum for a *distance* is direction-weighted, not |p|∞;
  least squares averages quantization down ~√N.
- **0.7 Tolerance estimators as diagnostics** (no behavior change): implement
  the fit-free sag estimator (`sag ≈ dihedral × facet size / 8`, computable in
  stage 2.1) and the dihedral-distribution angular-gap estimator
  (Otsu-style). Print both next to the stage-1 summary. Cross-check the sag
  estimator against the true sag stage-1 `--compare` already reports
  (8.9e-4 mm stepped fine, 9.5e-4 rounded_cube fine; Onshape angular settings
  are 12.5°/6.25°/2.5° for coarse/medium/fine). Wire a **provenance-violation
  warning** for estimator disagreement or a missing dihedral gap.
  **Stop and ask** if the sag estimator misses the reference truth by >5× on
  any corpus model.
- **0.8 Baseline runs of real meshes**: run 3DBenchy, OrcaCube, and
  xyz-cube-test end to end once each; record the diagnoses as harness rows.
  Fix nothing here; just know the starting point.

**Phase acceptance**: one-command harness with saved baselines; cliff count
recorded; coverage holes closed; estimator cross-check table produced.

---

## Phase 1 — Earned patches (the core rework)

Replace the blind seed patch with a patch **grown face-by-face under
validation**, promoting the surface type as the accumulated support allows.
Full design + math in DEVELOPMENT_PLAN ("The Robustness Rework", Phase 1
design). Commit semantics deliberately unchanged: earned (type, params,
patch) hands off to the existing `grow_from` flood/refit path, and every
commit gate (`MIN_AXIAL_SWEEP_DEG`, torus ratio, `MIN_*_FACES`,
`resolve_claims`) keeps its exact meaning.

- **1.1 Closed-form initializers, unit-tested first** (synthetic strips/tori;
  the corpus's exactness precedents — r = 0.99998-style — make ~1e-9
  assertions appropriate):
  - Cylinder: axis = n₁ × n₂ from two admitted faces with non-parallel
    normals (exact on strip/quad tessellation), then `fit_circle_for_axis`
    (cylindrical.rs — exists, measured well-conditioned on 40–65° arcs).
    Guard with a numeric conditioning floor on |n₁×n₂|, documented as
    validity-not-threshold.
  - Sphere: normal lines anchored at facet **circumcenters** (exact through
    the centre for inscribed facets); centre = common-perpendicular midpoint
    of two lines; radius = mean distance.
  - Cone: existing normal-covariance axis over the earned patch
    (`smallest_eigenvector_3x3`) + `axis_normals_consistent`.
  - Torus (**time-boxed — skip without blocking the phase if conditioning
    fails**): axis from the meridian-coplanarity system
    `(pᵢ − a)·(nᵢ × d) = 0` over ≥4 faces (bilinear; solve by alternation or
    pairwise elimination; rank-deficient on sphere/cylinder data by
    construction — ladder order and a conditioning floor must precede it);
    then (R, r) via the existing `fit_circle_2d` profile fit (toroidal.rs).
- **1.2 The earn loop.** Replace `seed_patch` + `condition_patch` inside
  `grow` (compete.rs): BFS from the seed; admit a frontier face iff it
  validates the current fit at vtol (`face_validates`); on failure attempt
  **promotion** in ladder order (plane → cylinder/cone or sphere → torus),
  refitting accepted+candidate and requiring **all already-accepted faces**
  to validate under the promoted fit; if nothing validates, the face is
  across a seam — exclude it and **record it as a frontier face** (Phase 2
  consumes these). **Branch on ambiguity** (e.g. two same-height cylinder
  faces validating as a sphere): grow both; false branches should self-prune
  within ~2 faces. Hand each surviving (type, patch, params) to the unchanged
  `grow_from`. Add counters (behind verbosity/env): NO-GROWTH count, branches
  spawned / died-within-2, promotions attempted/accepted, frontier faces.
- **1.3 Drift-triggered cylinder→torus promotion** (time-boxed). A bend is
  locally a valid cylinder, so admission never fails; the signal is
  **sign-consistent precession** of the fitted axis across successive refits
  as support accumulates — not a magnitude test (if a magnitude constant
  seems needed, stop and ask). If it works, `promote_bend` + `BEND_*` become
  Phase-4 deletion candidates; if not, keep them and move on.
- **1.4 Validate.**
  - *Initializer accuracy per provenance*: measure cylinder-axis error vs
    references across Parasolid strips AND the Phase-0 OCCT Delaunay meshes.
    If median error exceeds ~1° on a provenance, degrade that provenance to
    the covariance axis over the earned patch (still clean, still not blind).
    **Stop and ask** if even that underperforms the blind patch anywhere.
  - *Branch lifetime*: if branches persist past ~2 faces, validation is less
    sharp than believed — stop and ask before the cost model changes.
  - *NO-GROWTH count* (baseline 323 on stepped_block fine) should collapse
    toward zero.
  - *Acceptance*: probe per-type breakdowns equal-or-better corpus-wide;
    stepped_block fine ≥ today's 75-surface state; **stepped_block medium and
    part_rounded_cube coarse not worse** (the identifiability stress cases —
    the one place the blind patch did irreplaceable work); suite green;
    harness ≥ baseline.
- **1.5 Delete after ablation** (rule 4): `condition_patch`, `seed_patch`'s
  jet filtering, `SEED_PATCH_CAP`, `MAX_PATCH_DEPTH`, `PATCH_RADII`,
  `PATCH_TRIM_ROUNDS`, `MIN_CONDITION_FACES`; likely `MAX_TORUS_STARTS`.
  Report the constant count.

**Explicitly NOT fixed by this phase — do not chase here**: stepped fine's 4
surplus surfaces, medium's remaining gap, hill_with_hole's junk torus
(Phase 5c), runtime (Phase 2).

---

## Phase 2 — Distance-guided exploration and the runtime bound

- **2.1 Profile first.** Instrument the stepped_block fine run (103 s
  baseline): time in earn/patch construction vs LM fits vs floods vs
  `resolve_claims`. Report the split to the user. **Only build 2.2–2.4 if the
  seed sweep dominates**; otherwise apply the cheapest fix the profile names.
- **2.2 Distance field + farthest-first seeding.** Per-face geodesic distance
  to the region boundary (multi-source Dijkstra over face adjacency,
  **centroid-gap weights in mm, not hops** — hop counts lie on graded
  meshes). Order seeds farthest-first, replacing `neighbor_jet_mismatch`.
- **2.3 Duplicate abort.** Skip a seed whose reachable ball lies entirely
  inside one committed surface (it can only refit that surface) — the
  topological generalization of `interior_of_claims`; no numeric tolerance.
- **2.4 Frontier feedback** (only if the profile still warrants): add Phase
  1's excluded frontier faces as new distance sources incrementally — the
  exact tangent-junction map, discovered by fitting, never detected.
- **Acceptance**: results bit-identical to Phase 1's output, or every diff
  explained (**ordering-dependence = stop and surface it**); stepped fine
  runtime well under 103 s; record the number.

---

## Phase 3 — Flip the tolerance regime to measured values

Deferred until seeding is stable so the two changes are not debugged
concurrently (the earn loop admits on *vertex* tolerance, so Phase 1 did not
depend on this).

- **3.1** Add the fit-based sag estimator (per-facet centroid distance to its
  own committed fit on unambiguous regions — possible now that commits are
  trustworthy) and cross-check both estimators against reference truth across
  the full Phase-0 corpus including OCCT/Fusion/external pairs. **Stop and
  ask** on >5× disagreement anywhere.
- **3.2** Flip the `--surface-tolerance` default from 0.4 mm to the estimate
  (high quantile + safety margin — **never a median**; exporters also bind
  angular/min-width constraints, so legitimate facets sit above the median).
  Audit every dependent before flipping: `REFIT_GIVE_UP` (8×stol), merge's
  centroid gate, `face_validates`. Full harness sweep before/after; every
  regression explained or the flip reverted. The flag remains a user
  override; effective value = estimator unless overridden.
- **3.3** Angular tolerance from the dihedral gap where a clean gap exists;
  fall back to the current default otherwise, and *report which path was
  taken*. **Stop and ask** if no clean gap on more than 2 corpus models.
- **3.4** Document in README: recentring before export improves precision
  (quantum grows with |p|); decimated/repaired meshes violate the provenance
  premise and are detected, not silently mangled.

---

## Phase 4 — Ablation harvest and consolidation

- **4.1** Re-baseline the sweep harness and fragility meter under the new
  seeding + tolerance regime.
- **4.2** One candidate at a time, the standard protocol (disable → probe →
  suite → harness → delete-or-keep), in this order: `promote_bend` + `BEND_*`
  (if 1.3 landed); `rescue_with_neighbour_axes` (d2) and
  `rescue_leftover_components` (d) — earned patches may end the "never handed
  a clean face set" problem they exist for; `MAX_TORUS_STARTS`;
  `MAX_GROW_ROUNDS`; `MERGE_MAX_SMALL`; `SMALL_SELECTION_ABSORB_FACES`;
  `MIN_CYLINDER/SPHERE/TORUS_FACES`; then the degeneracy guards
  (`MIN_AXIAL_SWEEP_DEG`, `MAX_TORUS_RADIUS_RATIO`) — **expect these to still
  fire; firing means keep, not tune**. Finally ask whether the jet is needed
  beyond `region_is_developable` rung ordering.
- **4.3** `snap_boundaries` keep/delete: re-run its ablation and bring the
  numbers to the user (it is inert today but is the per-face assignment
  criterion a future VSA needs, and `dump_regions locus` shares the module —
  this is the user's call, flagged since 2026-07-28).
- **4.4** Radius-cloud diagnostic: fitted radii on stepped_block should
  tighten (was 2.9996/2.9997/3.0004, 1.0009/1.0018) now that supports are
  clean. If clouds persist, the supports still disagree — investigate that;
  never cluster, never snap to round numbers (recorded decisions).
- **4.5** Report the headline pair (constant count from 34; cliff count) in
  the phase-closing commit.

---

## Phase 5 — The research problems (each behind a measurement gate)

- **5a. Freeform fit conditioning. [parallel-ok — independent of Phases 1–4]**
  `bspline_blister`'s 926-face component fails the §2.7 pole-bound check
  ("control points 2749.5 mm outside the data bounding box"), saddle's three
  components likewise; nothing recovers and the regions fall to per-facet.
  The check is correct — fix the conditioning (Floater parameterization /
  knot placement / domain), do NOT loosen the check. Acceptance: both models
  decompose and pass their `--compare` blockers; then re-examine their
  ignored invalid-face tests (4 and 3 micro-sliver faces).
- **5b. Quantum as the unit of fit quality.** Report committed-fit residuals
  in quanta as a diagnostic (genuine ≤ ~1 q across a 7× quantum range;
  spurious ~80 q). Any use as a commit criterion is a user decision made on
  corpus measurements — bring numbers, not a patch.
- **5c. Analytic-vs-freeform model selection** (hill_with_hole's junk torus;
  the strip-vs-B-spline residue). **Hard gate**: before implementing
  anything, measure the relative-area separation across the whole corpus +
  harness (the junk torus is 0.095% of its patch; `cone_tangent_planes`'
  real cones are large fractions — and that model plus the ellipticals are
  what killed four prior rules; re-check them specifically). Bring the
  distribution to the user; only then choose a mechanism. Also add a test
  pinning `cone_tangent_planes`' currently *accidental* pass (its freeform
  rung fails as ill-conditioned; if 5a makes freeform stronger, this model is
  the canary). Afterwards: hill_with_hole's downstream stage-4.1 volume miss
  (6.96e-2 when the decomposition is right) is a separate freeform-accuracy/
  export bug — diagnose once candidacy is fixed.

---

## Phase 6 — Real-world hardening

- **6.1 Graceful degradation.** When the Phase-0 provenance detector fires
  (decimated/remeshed/faceted-by-construction input): warn precisely, run
  best-effort (freeform + per-facet), report what the data could support,
  never crash, never emit silent garbage. Exercise it on the quarantined
  external models and 3DBenchy.
- **6.2 Runtime scaling.** Measure on 100k–1M-face meshes from the external
  corpus; set and document a budget; profile-driven fixes only.
- **6.3 OpenSCAD / faceted-by-construction policy.** Generate two corpus
  models (low/high `$fn`), then implement the policy as an explicit, reported
  decision (flag or auto-detect-and-announce) per DEVELOPMENT_PLAN "Mesh
  provenance" — the prism-vs-cylinder ambiguity is genuine and cannot be
  measured away.
- **6.4** Promote backlog items as models demand (OBJ/PLY input, boundary
  crack stitching, IGES export).

---

## Reserve — partition-level architecture (user sign-off required)

If, **as measured by the harness**, Phases 1–4 leave seeding quality or sweep
cost unacceptable, the remaining lever is replacing the per-face sweep with a
partition: recursive split-and-merge (earned-patch fits at nodes,
union-revalidate merges via `fit_validate`, identifiability guards becoming
"defer, don't reject" until merge time) and/or VSA-style global alternation.
Both are deliberately reserved: they rewrite known-working machinery
(`resolve_claims`, relabel, merge) and their known hard parts (cut lines vs
seams; inverted identifiability) are documented in DEVELOPMENT_PLAN
("Considered and not adopted"). **Do not start without bringing Phase 1–4
measurements to the user.**

---

# Backlog (deferred until a model needs it)

- Boundary crack stitching (mesh repair): pair nearly-coincident open
  boundary edges and re-weld at progressive tolerance (≤ vtol), as netfabb
  "stitch triangles" / CGAL `stitch_borders` do. Slot:
  `ConnectedMesh::repair()` in src/stage1.rs, between the T-junction pass and
  `repair_coincident_faces`. No current model exhibits it.
- OBJ, PLY input formats; IGES export option.
- Dead-code audit leftovers (2026-07-15): `EdgeHighlight`/`edge_highlights`
  in src/viz.rs (never populated — wire up or delete); `error_abs_sum`
  accumulated but never read on five hypothesis structs (wire into a
  diagnostic or delete); `--debug` flag parsed but unused; `--step-units`
  parsed but never consulted by the STEP writer; `impl Display for Units`
  never invoked.

# Dropped (motivating problems solved by other means — kept as a record)

- **Gaussian map pre-classification**: performance idea from when stage 2.2
  wasted ~45 s tiling spheres with bogus cylinder bands; fixed by spherical
  detection order, per-face normal-axis consistency, early rejection.
- **Unified SOR framework refactor**: its motivations (pipe-elbow torus
  recovery, cyl/cone/sphere disambiguation) were fixed by targeted tools
  (strip-merge torus seeding, stage ordering, normal-line rescues); add a
  minimal (h,r) profile check where needed, never a framework.
- **Self-calibrating acceptance tolerance** (built 2026-07-27, removed
  2026-07-28): masked a seed-conditioning defect; the flood bar is
  `--vertex-tolerance`. See DEVELOPMENT_PLAN "The flood bar is the vertex
  tolerance".
- **LEFTOVER_BLOB_MAX** whole-region discard (removed 2026-07-28): replaced
  by `dissolve_fragmented_surfaces`; size is the wrong quantity, measured.
