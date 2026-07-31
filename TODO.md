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
2. **Measurement protocol**: save a `tools/sweep.sh --baseline` before each
   phase (and `tools/probe.sh` for a fast look);
   diff **per-type breakdowns** (plane/cyl/cone/sphere/torus/bspline), never
   surface counts alone (a 3-torus regression once hid behind an unchanged
   total of 71 — and measured again in Phase 0.3, the `--compare` gates pass a
   9-surface decomposition of a 26-surface part, so counts are the only thing
   that catches over-segmentation). The full suite is the decider.
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

**Diagnostics toolbox**: `tools/sweep.sh` (109-model harness, baselines and
diffs); `tools/fragility.sh` (per-constant ±2× perturbation, cliff count);
`tools/probe.sh` (14-model per-type breakdown);
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
- **0.2 OCCT re-mesh utility.** **DONE 2026-07-29** — `gen-test-models remesh
  [<out_dir>] [<step>...]` re-tessellates any reference STEP at
  coarse/medium/fine (deflection 0.5/0.1/0.02 mm, 30°/15°/6°) plus a
  ~1 m off-origin variant with its own translated reference STEP. 8 models ×
  4 variants committed under `tests/remesh/` (5.4 MB, byte-deterministic), all
  32 wired into the suite: **558 passed / 16 ignored**, stage-4 target
  107 s → 117 s. **25 pass, 7 fail, and the failures are the finding**:
  `stepped_block_toroidal` fails at *every* OCCT density — 183 surfaces coarse,
  88 medium, **97 fine at 35596 triangles** — where the Parasolid fine mesh
  reaches 75 against a reference 71. All 18 tori are recovered every time; the
  planes fragment. (Read at the time as the price of the strip assumption;
  **0.5a's Fusion re-export refuted that** — all three tessellators fail and
  triangle count does not order the results. The shared quantity is narrow flat
  regions. See 0.5a.) Also: coarse `sphere_25` (6 planes + 2
  spheres), coarse `rounded_cube` (34 vs 26, with 4 spurious tori), coarse
  `cone_tangent_planes` (10 vs 7, "cone-plane ruling tangency is not
  implemented"). Off-origin quantum measured 8.8e-5–1.02e-4 mm vs the 5e-5
  default vtol, yet max vertex distance is only 2.7e-5–4.2e-5 ≈ **0.42×
  quantum** — so the naive quantum is not the bound on a distance; 0.6 derives
  the one that is, and the default vtol has ~14% margin at a metre out. Phase 5c
  note: on the Delaunay mesh `cone_tangent_planes` resolves as 4 planes +
  1 cone + **2 B-splines**, a different decomposition from the strip mesh's
  4 planes + 2 cones + 1 sphere — both inside the 1% volume gate.
- **0.3 Sweep harness.** **DONE 2026-07-29** — `tools/sweep.sh` (probe.sh moved
  to `tools/probe.sh`; both now tracked). One TSV row per model: status, ref,
  total, per-type breakdown, stage 3.1 edge distance, stage 4.1 volume error,
  runtime, first error. `--baseline` / `--diff` / `--diff-only` / `--list`;
  parallel with `--jobs` (runtime column only comparable at `--jobs 1`, and a
  diff refuses to compare runtimes across job counts). Model/reference pairing
  is rule-based with the re-mesh sources **recorded** by the generator in
  `tests/remesh/sources.tsv` — searching `tests/` for a matching stem is
  ambiguous (two directories hold a `plate_with_hole_100x50.step`) and
  filesystem-order dependent. Needed one pipeline change: `Stage 2 types:` now
  prints the per-type breakdown at `-v`, because getting it previously meant a
  second full stage-2 run via `dump_regions`.
  `gen-test-models random` adds 5 seeded families × 3 (filleted boxes, cone
  chains, tori across r/R **including 0.84 > the 0.7 ratio cap**, tubes,
  ball-on-cylinder) under `tests/random/` with all parameters in `params.tsv`.
  **Assumption check PASSED**: all 14 `probe.sh` models match the harness
  exactly, and the harness's failure set over suite models is exactly the
  suite's 12 ignored stage-4 tests.
  **Baseline: 109 models, 92 ok, 16 fail, 1 timeout** (`tmp/sweep-baseline.tsv`).
  Four new findings, none chased — all are Phase 1/4/5 material:
  - **Filleted boxes collapse.** All three `rand_filleted_box` give
    **0 planes / 0 cylinders / 8 spheres / 1 B-spline** against a reference
    6 planes + 12 cylinders + 8 spheres. Fillet-to-face ratios are 0.03–0.10
    where `rounded_cube_10_r2` (the corpus's only filleted box, which passes
    exactly) is 0.20 — so the corpus tests a ratio 2–7× larger than typical,
    and at typical ratios the whole flat structure is swallowed by one
    freeform patch.
  - **A fat torus is shredded.** `rand_torus_1` (r/R = 0.84, deliberately over
    `MAX_TORUS_RADIUS_RATIO = 0.7`) yields **40 surfaces: 17 cones + 20
    spheres + 2 tori + 1 B-spline** for a single-face reference. The guard
    aimed at spurious fits rejects a legitimate solid; ratios 0.36 and 0.10
    are perfect. Phase 4.2 should read this before deciding that guard's fate.
  - **The compare gates do not catch over-segmentation.** `rand_filleted_box_0`
    (9 surfaces vs 26, no planes at all) PASSES at 6.7e-3 relative volume, and
    `rand_cone_chain_1` (a cone became a B-spline) passes at 2.0e-3. Stage
    3.1's edge compare is one-sided — our edges to the nearest reference edge —
    so a decomposition with *fewer* edges passes trivially: all these rows read
    `edge_mm = 0.0`. **Per-type counts are the only thing that catches it**,
    which is exactly why the protocol demands them.
  - **A near-tangent sphere/cylinder junction fragments.**
    `rand_ball_on_cylinder_0` (ball 13.98 on post 12.30) splits the sphere into
    3 and adds a plane: 6 surfaces vs 3.
  Also fixed a generator flaw the first sweep exposed: two of three
  ball-on-cylinder draws had the ball CONTAINING the post, so the reference was
  a plain sphere and the row measured nothing.
- **0.4 Fragility meter.** **DONE 2026-07-29** — `tools/fragility.sh` rewrites
  one constant at a time as `(original expression) × 2` and `÷ 2` (so the form
  of the expression does not matter), rebuilds, sweeps a 45-model fragility set,
  and diffs. Restores the source from a byte copy on any exit path. Reports
  *changed* (sensitivity) and *cliff* (a passing model that now fails, or a
  surface count that at least doubled — the per-facet-explosion signature).
  **Constant count is 32, not 34**: the grep counts a `usize::MAX` sentinel
  (`CUT_OUTSIDE`) and a viz colour (`SHARP_EDGE_COLOR`), both skipped with a
  recorded reason.
  **BASELINE CLIFF COUNT: 9 of 32.** Full table in `tmp/fragility-full.tsv`;
  the headline is that the cliffs cluster on the blind patch, which is the
  rework's central claim now measured rather than argued:
  | constant | ÷2 chg/cliff | ×2 chg/cliff |
  |---|---|---|
  | `MAX_PATCH_DEPTH` | 12 / **9** | 7 / **3** |
  | `MIN_CONDITION_FACES` | 8 / **5** | 10 / **7** |
  | `PATCH_TRIM_ROUNDS` | 7 / **5** | 6 / **3** |
  | `SEED_PATCH_CAP` | 6 / **3** | 0 / 0 |
  | `NORMAL_AXIS_TOL` | 6 / **5** | 3 / 0 |
  | `MAX_TORUS_RADIUS_RATIO` | 5 / **4** | 4 / 0 |
  | `MIN_HALF_ANGLE_DEG` | 0 / 0 | 3 / **2** |
  | `PROJ_GAP_RECHECK` | 0 / 0 | 2 / **2** |
  | `MIN_BSPLINE_FACES` | 3 / **1** | 1 / 0 |
  Four of the nine (`MAX_PATCH_DEPTH`, `MIN_CONDITION_FACES`,
  `PATCH_TRIM_ROUNDS`, `SEED_PATCH_CAP`) are blind-patch constants **already on
  Phase 1.5's deletion list**, and they account for most of the collapses. The
  collapsing models are overwhelmingly the fillet family — `rounded_cube` at all
  three provenances, `filleted_pipe`/`_cylinder`/`_hole_block`, all three
  `rand_filleted_box` — matching 0.3's finding that fillet-to-face ratio is
  where the corpus was thin.
  **`NORMAL_AXIS_TOL` (conical.rs, 0.05) is the unexpected one**: halving it
  collapses 5 models (all three `rand_cone_chain`, `block_with_conical_hole`,
  `remesh/cone_tangent_planes_medium`) and it is on no deletion list. Cone-axis
  work should treat it as load-bearing and fragile.
  **16 of 32 constants are completely insensitive to ±2×**, including
  `REFIT_GIVE_UP`, `MIN_AXIAL_SWEEP_DEG`, `MAX_TORUS_STARTS`, both `BEND_*`,
  `COAXIAL_TOL`, `TANGENCY_TOL_RAD`, `MAX_FIT_FACES`, `MAX_REFIT`, `REFIT_ADD`,
  `DEVELOPABLE_MAX_RATIO`, `MAX_GROW_ROUNDS`, `MAX_SPHERE_RADIUS_FACTOR`.
  **Read that as "its value is not load-bearing", NOT as a deletion list** — a
  guard that fires far from its boundary is exactly what a guard should look
  like; `MIN_AXIAL_SWEEP_DEG` is measured essential and ±2×-insensitive at once.
  Phase 4.2's ablations still have to disable, not perturb.
  **Cadence**: run the full matrix at *phase boundaries*, not per commit; per
  change use `--only` on the constants the change touches (minutes). The meter
  now parallelizes *across* perturbations — waves of `--wave` stashed binaries
  swept concurrently at `--jobs` each — measured **61 s per perturbation vs
  105 s**, so the full matrix is **~39 min, down from ~110** (utilization 39% →
  63%). 5×2 is optimal at 10 cores: a sweep's wall is
  `max(101 s longest model, 315 s serial / jobs)`, so J=4 would idle the wave on
  one model. Phase 2.1's runtime work on `stepped_block_medium` is what lifts
  the ceiling. Validated by reproducing the serial run's counts *and*
  collapsed-model lists on 4 constants / 8 perturbations.
- **0.5 Real-pair corpus. [parallel-ok, user-involved]**
  - (a) **Fusion re-export**: **DONE 2026-07-29** — `tools/fusion/BrepperReexport.py`
    (a Fusion Python script, one click, installed by symlink) imported 13
    reference STEPs and re-exported 31 STLs; user ran it on Fusion 2704.1.36.
    **Verification rate 15/15** (0.5c gate applied): every mesh's vertices lie on
    the ORIGINAL reference surfaces at 4.1e-7–7.9e-6 mm, so Fusion's STEP import
    neither healed nor reparameterized the geometry. Committed as
    `tests/fusion_reexport/` (1.8 MB, 15 distinct meshes) and wired into the
    suite: **suite 559/16 → 570/20**, stage-4 target 117 s → 119 s.
    **THE KEY EXPERIMENT CAME BACK NEGATIVE — AND MORE USEFUL FOR IT.** The
    framing was "if Fusion also fails `stepped_block_toroidal`, the problem is
    anything-but-strips; if it passes, it is Delaunay-specific". Neither:
    | provenance | triangles | surfaces | plane/cyl/cone/sph/torus/bspl |
    |---|---|---|---|
    | Parasolid fine | 50962 | **75** | 18/38/0/0/18/1 |
    | Parasolid medium | 4576 | 84 | 33/31/0/0/19/1 |
    | OCCT fine | 35596 | 97 | 41/36/0/0/18/2 |
    | OCCT medium | 6464 | 88 | 34/33/0/0/18/3 |
    | OCCT coarse | 1948 | 183 | 149/16/0/0/18/0 |
    | Fusion fine (0.02 mm / 6°) | 17598 | **82** | 27/36/0/0/18/1 |
    | Fusion medium (0.1 mm / 15°) | 2948 | 129 | 82/29/0/0/18/0 |
    | Fusion coarse (0.5 mm / 30°) | 866 | — | rejected at stage 1: 0.51 mm sag > 0.4 mm |

    (OCCT and Fusion rows are at IDENTICAL requested deviations. Fusion yields
    about half the triangles for the same request, so a nominal tolerance is not
    comparable across tessellators — and at fine it reaches 82 against OCCT's 97
    with half the triangles.)
    All three tessellators fail (best is Parasolid fine at 75 vs the reference
    71), so it is not Delaunay-specific. **And it is not triangle count either**:
    OCCT fine has 7.8× Parasolid medium's triangles and does worse (97 vs 84);
    two Delaunay-family meshes at ~6.5K triangles land at 88 and 109 while
    Parasolid at 4.6K lands at 84. **The surplus is almost entirely PLANES** —
    18 at best, 232 at worst — while the 18 tori survive every provenance above
    ~3.7K triangles. Narrow flat regions between fillets are the failure, which
    is exactly what the blind patch straddles and what Phase 1's earned patch
    targets. Stop attributing this to Delaunay.
    Two more findings: **`cone_tangent_planes` agrees with Parasolid** (4 planes
    + 2 cones + 1 sphere), so OCCT's B-spline reading is provenance-specific, not
    the model's nature — Phase 5c should read it that way; and **`nosecone` is a
    NEW failure**, 10 surfaces vs 3, from ~3000 mm² of its flat base being
    absorbed into the cone *plus* the remainder fragmenting into six planes with
    byte-identical parameters (connectivity, not fitting — the
    'absorbing bridge facets' family in LEARNINGS.md, now a hard failure).
    **DONE 2026-07-30 after six runs.** `tests/fusion_reexport/`, 13 models × 31
    exports, 1.8 MB, all **31/31 verified** by the 0.5c gate (mesh vertices on
    the original reference surfaces, 1.8e-7–7.3e-6 mm). Wired into the suite:
    15 pass, 16 ignored with numbers. Densities are **numeric and identical to
    what `tests/remesh` asks of OCCT** — 0.5/0.1/0.02 mm surface deviation at
    30/15/6° — so Fusion-vs-OCCT is a controlled comparison, not a comparison of
    two vendors' preset names.
    **THE HEADLINE, at matched requested tolerances**: Fusion produces roughly
    HALF the triangles OCCT does for the same requested deviation (2948 vs 6464
    at medium; 17598 vs 35596 at fine) — so nominal tolerance settings are not
    comparable across tessellators even when numerically identical. And at fine,
    Fusion's 82 surfaces beat OCCT's 97 with half the triangles, while Parasolid
    fine still leads at 75 with 3× more.
    **`rounded_cube_10_r2` reproduces the random-corpus collapse on a CORPUS
    model.** At Fusion medium (684 triangles) it gives **9 surfaces — 0 planes,
    0 cylinders, 8 spheres, 1 B-spline** — failing volume by 1.94e-2, which is
    the exact signature the seeded `rand_filleted_box` models showed in 0.3
    (0/0/8/1 against a reference 6/12/8). `part_rounded_cube` medium collapses
    the same way to 3 surfaces. Fine is exact at 26. This is no longer a
    property of generated geometry: the corpus's canonical fillet model does it
    too, at an ordinary export density.
    **Four models are rejected at STAGE 1, not by the decomposition**: a mesh
    exported at 0.5 mm chordal deviation has facet-centroid sag of 0.51–0.83 mm,
    over the 0.4 mm default `--surface-tolerance`. That is the tolerance regime
    refusing an ordinary coarse export, and it is a **Phase 3 datum** — the
    default rejects any mesh coarser than ~0.4 mm deviation. It only became
    visible once the centroid check stopped measuring against a solid.
    **Six runs to get here, and every failure was the same mistake**: trusting a
    stated interface over evidence already in hand. Presets that silently did
    nothing; `unitType` defaulting to the Design's units (metres) so 31 files
    came out 1000× too small; and `normalDeviation` documented as radians when
    Fusion's own preset readback in my probe log said 30/15/10 — degrees.
    `STLExportOptions` inherits ambient state, so every property is now pinned.
    The structural fix is a **canary**: the script exports `sphere_25` alone,
    checks scale and strict monotonicity, and aborts in seconds before touching
    `tests/` — sharing `configure_export()` with the batch, since a pre-flight
    that configures differently proves nothing. `BrepperProbe.py` and
    `BrepperBisect.py` are kept as the record of what was ruled out.
  - (b) **Wild pairs**: **HARVESTED AND MEASURED 2026-07-30 — STOP AND ASK, per
    this step's own criterion (most wild pairs fail).** `tools/harvest_github.py`
    takes STEP+STL pairs from open-source hardware repos on GitHub rather than a
    model-sharing site, because those repos commit both files from one author and
    one revision under an explicit, machine-readable licence. 18 pairs from **6
    unrelated projects** (Voron printer mods, E3D ToolChanger, OpenBot robotics,
    White-Knight, NucDeck, Open-Task-Light), 16 MB, **all redistributable**
    (GPL-3.0 / MIT / CERN-OHL-S-2.0). Models gitignored; `manifest.tsv` records
    URL, repo, author, licence, path, size, sha256 and the pairing rule. The
    Thingiverse token is no longer needed — `gh` supplies auth.
    **RAW VERIFICATION RATE: 4 of 18 (22%)** — but that number mostly measured
    COORDINATE FRAMES. Alignment was added (user-approved) and the corpus
    resolved; the raw rate is kept above because the gap between it and the final
    one is the finding.

    **Alignment is discrete and never fitted.** `gen-test-models align` derives a
    transform from bounding-box centres plus the first of the **24 proper
    rotations of a cube** whose *unsorted* extents match, writing
    `<stem>_aligned.step`. It never sees the distance it is helping to pass, so
    it cannot make a wrong pair look right — the failure mode the question
    flagged. 10 of 18 were re-framed (8 translation-only, 2 with a rotation).

    **A FIRST PASS OF THIS GOT IT WRONG AND DELETED REAL GROUND TRUTH. Read this
    before trusting any bounding box.** `shape_bounds` called
    `BRepBndLib::Add(shape, box, useTriangulation = true)`, and that flag means
    "use a triangulation IF THE FACE HAS ONE". A shape straight from STEP has
    none, so OCCT fell back to bounds enclosing the surface POLES — the control
    hull, not the surface. On `OrcaCube v2.step`, a 30 mm cube, the padded box
    reads **34.8 x 56.0 x 30.9**: one axis 4.8 mm too large. Alignment compared
    that against an *exact* mesh box, so both the size test and the
    centre-derived translation were wrong, and 7 pairs were deleted as "different
    objects" on the strength of it. **`VoronUsers__10x11...` was one of them and
    is genuine ground truth: it verifies at 4.49e-5 mm once rotated.** All 7 were
    restored from the manifest's URLs and sha256 — which is the reason to record
    those. Two fixes: mesh the shape before bounding it, and do not apply a
    translation below a noise floor (see finding 3). Nothing is deleted now; the
    corpus is gitignored, so a failing pair costs nothing and is evidence.

    **FINAL: 7 of 18 verify (39%); 11 of 18 correspond (61%).** The denominator
    is 18, not 11: deleting the failures and then reporting a rate over what
    survived is how the first pass talked itself into "6 of 11".
    | outcome | n | what it is |
    |---|---|---|
    | VERIFIED at `vtol` 5e-5 | 7 | **4 at or below one float32 quantum** |
    | corresponds, stray-vertex tail | 3 | ToolChanger ducts: p95 6.9e-6 **below** their 9.3e-6 f32 bound, max 4.3e-4. `NucDeck__back_cover`: p95 1.3e-4, max 3.5e-1 |
    | corresponds, alignment residual | 1 | `holder_arm`: p95 4.8e-4 ≈ max — distributed, 10% of its own sag |
    | does not correspond | 5 | belt clips ×2 at a uniform **exactly 0.2 mm** (an offset/clearance variant); White-Knight ×3 at p95 8–24 mm |
    | mesh invalid | 2 | `OpenBot__100deg_joint` ±mirror: 2 non-manifold edges, fails stage 1 before compare |
    **The premise holds where the pairs are real.** On files we did not export,
    the median vertex sits on the CAD surface at the float32 quantum. Findings:
    1. **About a third of published STEP+STL pairs are not usable ground truth**
       (5 different objects + 2 invalid meshes of 18), and **10 more of the 18
       needed re-framing first**. Coordinate frame is the single largest cause,
       and a wild-pair gate must align before it can measure anything. Proven,
       not inferred: translating `NucDeck__+_volume_button` by (−0.5, 0, −6.5) mm
       — round numbers, a deliberate placement — drops its vertex distance from
       **6.52 mm to 4.58e-7 mm**.
    2. **A max-based vertex gate is the wrong shape for wild meshes.** Both
       ToolChanger ducts are exact to float32 at p95 and fail only on a thin
       tail. Stage 1's compare now prints **p50 / p95 / max** so a local mesh
       defect can be told from a systematic offset; Phase 3 should decide whether
       the gate itself becomes a quantile. See §3.2.
    3. **Bounding-box alignment cannot beat the mesh's own sag, so below that it
       must not act.** A tessellated mesh's bbox is inset from the true solid
       asymmetrically, so a centre-derived translation carries a sag-scale error.
       Applying it where none was needed makes a good pair *worse*: the
       ToolChanger ducts went from p95 6.9e-6 to 4.4e-4 mm under a 4.4e-4 mm
       "correction". `align` therefore ignores translations below a noise floor.
       The floor is not tuned — the corpus's shifts are bimodal with a **4000×
       empty gap** (noise at 8e-8, 4.4e-4, 4.9e-4, 1.6e-3 mm; real displacements
       at 6.5, 24.2, 101.9, 112.1 mm), so any cut inside it gives the same
       answer, and the skipped magnitude is printed so a future corpus landing
       *in* the gap is visible rather than silently decided. Above the floor it
       recovers exact round-number placements perfectly (NucDeck, 4.58e-7 mm) and
       nothing finer — what remains is the sag-scale residual on `holder_arm`.
       Sub-sag alignment would require fitting, which is why it is not done.
    4. **Text precision is not a factor either way** — the SolidWorks STEPs carry
       20 significant digits, the ST-Developer ones 15. Checked because the
       residuals were near ASCII-STL resolution; they are not that.
    5. **Real-world STLs are not always valid**, and **some STEPs are assemblies**
       where the STL is one part — the gate handles the latter correctly
       (`Open-Task-Light__counterweight` verifies at 1.9e-7 mm while covering a
       fraction of the STEP's extent).
    6. **Verifying and being worth testing are different properties, and the
       first 18 were selected for neither.** The harvester took the first three
       pairs per repo alphabetically, so what passed was what had least to get
       wrong: of the models that verified, one was a 52-triangle part with **no
       curved surface at all**, one was tessellated at 18° (5 facets per 90° of
       arc), and the best was 17. The three densest meshes in the corpus were all
       in quarantine. Two numbers now grade every mesh, in `manifest.tsv` as
       columns — `facets_per_90` (how many facets span a quarter turn, which
       decides whether a cylinder or torus is recoverable) and `smooth_pct` (how
       much of the mesh is on curved surface at all). **Never quote a pass rate
       over the coarse rows.**
    7. **The vertex verdict is taken on p95, not max** (`verify_pairs.sh`), with
       both reported plus their ratio. A wild mesh can be exact to the file
       format across almost all vertices and carry a thin tail of strays from
       repair or decimation; the max throws the whole pair away for those few.
       The ToolChanger ducts — the densest meshes here at 64–71 facets/90° — have
       p95 6.9e-6 mm, **below their own 9.3e-6 float32 bound**, with a max of
       4.3e-4. That took the corpus from 7 verified to **9 of 18**, and those two
       were the *only* GOOD-density models in the verified set. The ratio is a
       description, not a verdict: ≈1 means a uniform offset (wrong frame, wrong
       revision, an offset variant), in the tens or hundreds means a few bad
       vertices on an otherwise exact mesh. `holder_arm` is correctly NOT
       rescued — ratio 1, and its bbox centres match its reference to 4.1e-12 mm,
       so the uniform 4.8e-4 offset is a mesh that is not a tessellation of that
       STEP (decimated, optimized, or a near-identical revision).
    8. **Corpus expanded on density 2026-07-30**: `harvest_github.py` now grades
       each STL and downloads the STL *first*, so a reject costs one file instead
       of two. 266 pairs were available across the six usable repos; selecting
       for ≥20 facets/90° and ≥10% smooth took the corpus from 18 pairs to **53**,
       and the grade mix from 4 GOOD / 18 to **39 GOOD, 7 MODERATE, 6 WEAK, 1
       TRIVIAL**. The rejection log is itself a finding: **29 of 30 OpenBot
       candidates were rejected**, that whole project exporting at 5–10
       facets/90°, while every ToolChanger and Voron candidate passed. Mesh
       density is a property of the publisher, not of the part.
       **Result: 26 of 53 verify (49%), and 19 of those are GOOD-density —
       against 2 before.** That is the number that matters; the rate barely
       moved, the usable dense ground truth went up 10×.
       | grade | verified | quarantined | rate |
       |---|---|---|---|
       | GOOD | 19 | 20 | 49% |
       | MODERATE | 4 | 3 | 57% |
       | WEAK | 2 | 4 | 33% |
       **Verification rate is independent of mesh grade**, so pair correctness
       and mesh density are independent properties — density selection buys
       better tests, not a better pass rate. The 27 quarantined split into:
       | cause | n | |
       |---|---|---|
       | gross mismatch, p95 > 1 mm | 11 | 4 of them are **an STL print-plate against a single-part STEP** (`10x11_chain_top`: STL 41.8×20.2×10.0 vs STEP 20.0×20.2×20.1) — a new category, correctly refused by the aligner with extents reported |
       | near miss, p95 < 1e-2 mm | 8 | **5 are NucDeck clustering at 1.3–1.7e-4 mm**, ~16× their f32 bound — systematic to that publisher, so likely an export setting rather than five independent revision drifts |
       | large offset | 4 | includes both belt clips at a uniform **exactly 0.2 mm** |
       | mesh invalid (non-manifold) | 4 | 2 are GOOD-density Voron ziptie mounts — real crash-test material Phase 6.1 has none of |
    9. **`shape_bounds` meshes at a deflection that scales with the part.** A
       fixed 1e-4 mm was wasted precision on a 10 mm button and did not finish a
       4.8 MB, 200 mm enclosure — it stalled the corpus run. It now sizes the
       part from the (over-estimating, so safe-direction) pole box and meshes at
       `diag × 5e-6` clamped to [1e-4, 5e-3], an order of magnitude tighter than
       the noise floor it feeds. Whole 53-model corpus in 246 s, with transforms
       byte-identical to the slow version on every shared model.
  - (c) **Pair-verification gate**: **DONE 2026-07-29** — `tools/verify_pairs.sh`,
    versioned rather than an ad-hoc loop, because the ad-hoc loop **reported
    15/15 on a corpus that was entirely unusable**. Its hole: stage 1's vertex
    compare returns distance **exactly 0** for a vertex strictly *inside* the
    reference solid, so a mesh at the wrong scale or position passes trivially —
    measured on 3 metres-scale exports, where a 25 mm sphere spanning 0.025 mm
    reported 0.0 while a correct mesh of the same reference reports 1.1e-6. The
    gate now classifies **VERIFIED** (passes with a nonzero distance),
    **CONTAINED** (passes with exactly zero — suspect scale or placement), and
    **QUARANTINE** (fails), and reports the f32 position bound, which reveals
    coordinate scale directly (5.7e-10 mm for the metres sphere against 7.2e-7 mm
    for the correct one). Verified against the known-bad corpus: 0/31, catching
    all 28 failures *and* the 3 that the old loop passed. Apply to (b)'s wild
    pairs, where the rate is the actual measurement. **Stop and ask** if most
    wild pairs fail it.
- **0.6 Quantum clamp + report** **DONE 2026-07-29** — stage 1 reports
  `float32 precision: position bound B mm (quantum Q mm), vertex tolerance T`
  and warns when T < B; the CLI driver then raises the tolerance to B
  (`stage1::effective_vertex_tolerance_mm`; applied in the driver, not inside
  stage 1, because `Config` is immutable through the pipeline and ~40 sites
  read `vertex_tolerance_mm` — threading a derived value belongs with Phase
  3.2's audit).
  The bound is **not** `max|coord| × 1.19e-7`, honoring the recorded DO-NOT:
  it is `max over vertices of ‖(ulpₓ, ulp_y, ulp_z)‖ / 2`, the largest
  projection of the per-coordinate rounding box onto any direction. Measured on
  the off-origin corpus: bound 4.32e-5 mm where the true max vertex-to-surface
  distance is 3.93e-5 (bound holds, 0.91 tight) and the naive quantum says
  6.10e-5 — so the naive figure is 1.4× too loose here. Both are reported;
  only the bound is used.
  **The clamp is inert on the entire corpus, including at 1 m** (4.32e-5 <
  5e-5 default). It first bites past ~1.02 m, where the f32 exponent steps and
  the bound roughly doubles; pinned by a synthetic 5 m mesh in
  `float32_position_bound_scales_with_distance_from_origin`. Sweep diff after
  landing: no model changed.
- **0.7 Tolerance estimators as diagnostics.** **IMPLEMENTED AND MEASURED
  2026-07-30 — BLOCKED ON THE USER, per this step's own stop-and-ask.**
  `src/estimate.rs`, printed in the stage-1 summary, no behaviour change:
  dihedral p50/p95/max and the fit-free sagitta `dihedral × facet size / 8` as
  p50/p95/max over every interior edge.
  **The angular estimate works.** On the fused-quad mesh the dihedral p95 reads
  the exporter's setting directly: **2.50° on `rounded_cube_10_r2_fine`, against
  Onshape's documented fine setting of 2.5°**, and 2.52° on stepped_block fine.
  A polyhedron is self-evident too — `chamfered_cube_10_c1`'s p50 is 45°.
  **The sag estimate misses by >5× on 3 of 16 corpus models**, which is exactly
  the "stop and ask" this step defines. Cross-check against the true sag stage-1
  `--compare` reports (all ratios est_p95 / truth):
  | within 5× (13) | ratio | | over 5× (3) | ratio |
  |---|---|---|---|---|
  | stepped_block fine | 0.93× | | **stepped_block medium** | **22.6×** |
  | pipe_elbow medium | 1.26× | | **plate_with_hole coarse** | **151×** |
  | rounded_cube coarse | 1.61× | | **remesh sphere medium** | **0.20×** |
  | fusion rounded_cube fine | 1.78× | | | |
  | pill coarse | 1.96× | | | |
  | remesh rounded_cube medium | 2.11× | | | |
  | rounded_cube medium | 2.53× | | | |
  | nosecone medium | 4.89× | | | |
  | hemisphere / filleted_pipe | 0.43× / 0.49× | | | |
  | rounded_cube fine, sphere fine, pill fine | **0.50× each** | | | |
  **Two effects, cleanly separated once the dihedral quantiles are read
  correctly** (an earlier note here said `plate_with_hole coarse` had a p95 of
  0.00° — that was a broken shell parse on my side; it is 90.00°, which tells
  the opposite story):
  1. **Where smooth edges dominate, the estimate is systematically 2× LOW on
     doubly-curved geometry.** `sphere_25_fine` 0.50×, `hemisphere` 0.43×,
     `simple_sphere` 0.42×, `rounded_cube_fine` 0.50× (its worst sag is on the
     spherical corners), `filleted_pipe` 0.49× (tori). The reason is geometric,
     not a coding slip: `hθ/8` is the sagitta of a curve bending ONE way, and a
     sphere bends two ways at once, so both principal curvatures push the facet
     centre off the surface. Adding the second term is an exact fix, not a tune.
  2. **Where sharp edges are a large fraction of all edges, the p95 IS a sharp
     edge and the estimate is meaningless.** `simple_cylinder` has a dihedral
     **p50 of 90°** — the median edge of a coarse cylinder is a cap-to-wall rim,
     not a curvature step — giving 62×; `cylinder_10x30_medium` p95 90° gives
     404×; `plate_with_hole coarse` p95 90° gives 151×. A flat plate or a short
     cylinder is mostly rim.
  **The dihedral quantiles already announce which case a mesh is in**, at no
  extra cost: 2.49/2.50/2.50 on rounded_cube fine says every edge is a curvature
  step, while 12.41/90.00/90.00 on the plate says the p95 is a rim and no sag
  statistic from it means anything.
  **ANSWERED BY 0.8, and the answer is neither option.** The two effects above
  are real, but 0.8's real print models expose a third that makes the quantile
  itself the problem. Measured side by side:
  | corpus | model | dih p50 | dih p95 | sag p50 | sag p95 |
  |---|---|---|---|---|---|
  | synthetic | rounded_cube fine | 2.49° | 2.50° | 4.755e-4 | 4.764e-4 |
  | synthetic | sphere_25 fine | 2.50° | 2.50° | 2.972e-3 | 2.974e-3 |
  | synthetic | stepped_block_toroidal fine | 2.43° | 2.52° | 2.382e-4 | 3.561e-3 |
  | **real** | 3DBenchy | **0.79°** | **89.21°** | 1.805e-4 | 2.032e-2 |
  | **real** | OrcaCube v2 | **1.13°** | **90.00°** | 4.337e-4 | 2.030e-1 |
  | **real** | xyz-cube-test | **3.00°** | **90.00°** | 1.713e-3 | 6.233e-1 |
  On the synthetic corpus the dihedral distribution is a SINGLE SPIKE at the
  exporter's angular setting — p50 and p95 are the same number, because those
  models are curvature-dominated by construction (spheres, rounded cubes,
  fillets). That is the only reason a fixed quantile works there. Real print
  models are mostly flat panels with sharp rims and a modest curved fraction, so
  the distribution is genuinely trimodal and **p95 lands on a 90° rim in all
  three**. No fixed quantile lands on the curved mode, because which quantile
  would depends on each model's flat:curved:sharp proportions — that is the
  "detector needs just one threshold" trap in rule 5, reached from a new
  direction. Chasing the 0.50× factor would improve a number that is already
  measuring the wrong edges on exactly the meshes brepper is meant to handle.
  **Disposition**: keep the estimator as a diagnostic, unchanged. It is honest —
  it prints the dihedral quantiles that announce its own failure, and they fire
  correctly on all three real models. Phase 3.2 must NOT read the sag p95 as an
  automatic tolerance source. The **angular** half remains sound where the mesh
  is curvature-dominated (it recovers 2.5° to two digits) and useless elsewhere,
  which is the same condition. Deriving a per-face tolerance from local geometry
  (§3.2) is the route that does not depend on a global quantile at all.
- **0.8 Baseline runs of real meshes** **DONE 2026-07-30 — nothing fixed, as
  specified.** `tools/sweep.sh` already covers `tests/downloaded/` and already
  resolves both STEP references, so these are harness rows; the runs below add
  what a 300 s harness timeout cannot show.
  | model | tri | outcome |
  |---|---|---|
  | `xyz-cube-test-20mm` | 1 128 | **ok** — 58 surfaces, 1 solid, 2.3 s |
  | `OrcaCube v2` | 44 606 | **PANIC after 40 min** (2407 s wall) |
  | `3DBenchy` | 225 706 | **PANIC after 8 h 21 min** (30 089 s wall, 28 592 s user) |
  1. **Both real models crash in stage 2.4's `assert_invariants`, and between
     them they violate BOTH invariants it checks.** These are crashes, not
     timeouts; the harness rows said `fail` and `timeout`, which hid that.
     | model | assertion | meaning |
     |---|---|---|
     | OrcaCube | `repair.rs:741` — face 1460 covered 0 times | a face ends up in **no** selected surface |
     | 3DBenchy | `repair.rs:751` — surface 5 has 7 disconnected components | a surface's faces are **not mesh-connected** |
     The assertion is the symptom; the defect is upstream, in how surfaces claim
     faces. OrcaCube reaches it after 2263 tangency-snap attempts, nearly all
     reporting `NO LOCUS: surfaces do not share an axis`. **This is Phase 1's
     first concrete target** — the whole corpus of synthetic models never
     produces either violation.
  1b. **Runtime is superlinear**: 5.06× the triangles costs 12.5× the time
     (44 606 → 225 706 tri, 2407 → 30 089 s), about O(n^1.6). A Phase 6.2
     datapoint, and the reason the 300 s harness timeout hides these as
     `timeout` rather than reporting the crash.
  2. **Real meshes need repair before they can be read.** 3DBenchy: dropped 2
     duplicate-vertex faces, removed 11 slivers, split 11 at T-junctions.
     OrcaCube: removed 6 coincident faces. Both still report a dihedral **max of
     180°** afterwards, so folded facet pairs survive repair.
  3. **Multi-body is normal in the wild** — OrcaCube reads as 3 shells / 3
     solids. The synthetic corpus is almost entirely single-solid.
  4. **OrcaCube's `--compare` failure is a corpus fact, not a brepper defect,
     and the new distribution readout is what settles it**: vertex p50 3.83e-4 /
     p95 5.69e-3 / max 1.54e-2 mm against an f32 bound of 1.98e-6. The *median*
     vertex is 193× the bound off the surface — systematic, not a tail — and the
     p50 nearly equals the mesh's own estimated sag p50 (4.34e-4), i.e. the
     vertices do not lie on that STEP's surface; the facets straddle it. Its
     extents match to a 1e-3 mm translation, so this is not a frame difference.
     Treat that pair as not-ground-truth. **Finding the padded-bbox bug in 0.5b
     came out of chasing this**, since the first bounds reading claimed the two
     differed by 4.8 mm.

**Phase acceptance**: one-command harness with saved baselines; cliff count
recorded; coverage holes closed; estimator cross-check table produced.

**BASELINE REFRESHED 2026-07-30 at commit 09534ab** — 158 models (was 109; the
external, fusion_reexport, remesh and random corpora landed during Phase 0):
**107 ok, 50 fail, 1 timeout**. The old baseline predated the compare fixes, when
three gates read zero for every input, so it is not comparable row-for-row.
Exactly **three** models regressed, and all three are the instrument starting to
work rather than new defects — each is a true deviation that exceeds
`--surface-tolerance` 0.4 mm and used to be measured as 0.0:
| model | gate | now reports |
|---|---|---|
| `remesh/sphere_25_occt_coarse` | stage 1.2 centroid | 6.675e-1 mm |
| `random/rand_torus_1` | stage 1.2 centroid | 4.128e-1 mm |
| `random/rand_filleted_box_0` | stage 3.6 boundary | 4.299e-1 mm |
The first is predicted to 0.1% by theory — see §3.2's 4/3 derivation — so it is
the *tolerance* that is wrong, not the mesh. All three are Phase 3.2's business.

---

## Phase 1 — Earned patches (the core rework)

- **1.0 Reference labeling — ground truth per mesh face.** **DONE 2026-07-30.**
  Added after Phase 0 closed, because 1.4 requires "measure cylinder-axis error
  vs references across Parasolid strips AND the Phase-0 OCCT Delaunay meshes"
  and there was no instrument for it, and because the over-segmentation hole
  recorded in DEVELOPMENT_PLAN ("belongs in its own step, not smuggled into
  Phase 0") could now have that step. `src/reference.rs`, `--reference-report`
  (requires `--compare`), diagnostics only — nothing downstream reads it.
  **It costs no extra geometric work**: stage 1's compare already runs
  `DistShapeShape` from every face centroid to the reference face compound, and
  that call already decides which reference face won —
  `min_distance_to_shape` was discarding it. `SupportOnShape2` hands it back.
  Mesh-side scoring runs at stage 1 (so it still works where stage 2 crashes or
  times out); segmentation scoring runs at the end of stage 2.
  **Measured over 79 meshes in three provenances (157 483 facets):**
  1. **The closed-form axis error is `perpendicularity error / |n₁ × n₂|`** —
     verified to 0.1% on `pipe_elbow_10` Parasolid medium (max 0.0280° over p50
     0.0271 predicts 1.033°; measured 1.0342°, constant across all 348 pairs).
     This one relation explains every cylinder outlier in the sweep.
  2. **The Delaunay caveat did not reproduce, so 1.1 needs no provenance
     fallback.** OCCT re-meshes read perpendicularity p50 = p95 = **0.0000°**
     and axis error **0.0000°** at every density on 5 of 6 models: OCCT's
     Delaunay is *constrained* by the pre-discretized face boundaries, so a
     cylinder's triangles still span two rulings. Fusion agrees. The single
     genuine non-perpendicular case corpus-wide is `l_bracket_simple` Fusion
     medium (perpendicularity p50 0.2740°) on a small hole with 96 facets.
  3. **Adjacent facets are the WORST pairing, and finer meshes are worse.**
     `|n₁ × n₂| ≈ sin(angular step)` — p50 0.0432 on stepped_block fine, which
     is sin(2.48°) and recovers Onshape's 2.5° setting — so amplification is
     ~23× at fine. **1.2 should pair the seed with the angularly most distant
     admitted facet**: at 90° separation the factor is 1 and the axis error
     falls to the perpendicularity error itself, under 0.032° on every corpus
     mesh. (A 90° max reading is the degenerate two-facets-on-one-ruling pair —
     no information, not a wrong answer.)
  4. **The torus initializer should NOT be skipped.** Its meridian-coplanarity
     residual — the constraint 1.1 time-boxes — reads p50 **2.354e-5 mm** over
     23 863 torus facets on stepped_block fine, i.e. at the mesh quantum.
  5. **The seam assumption holds**, and this is its first measurement (0.5b
     measured only that *vertices* lie on the surface). 51 of 157 483 facets
     disagree across their own interior — **0.032%** — and the only recurring
     cluster is exactly 12 facets on `stepped_block_toroidal` under OCCT
     coarse/medium and Fusion medium/fine, absent under Parasolid and OCCT fine.
     Consistent with coarse facets spanning that model's known narrow flat
     regions; not chased. The rate is an *estimate, not a bound* — it
     over-reports near concave junctions and under-reports a straddle confined
     to a facet's rim; what makes it usable is that both error terms are
     bounded by a number this small.
  6. **Segmentation scoring separates the two failure directions**, which the
     existing gates cannot see at all: `rand_filleted_box_0` (9 surfaces
     against a reference 26, passes every gate at 6.7e-3 relative volume)
     scores ref coverage **1.0000** and surface purity **0.2201** — correctly
     diagnosed as merging, not splitting. A correct model scores 1.0000/1.0000.
  **Extended 2026-07-30 to freeform reference faces** once `GeomLProp_SLProps`
  was bound (the normal report previously skipped them, covering only 32–66% of
  a NURBS model's area). Over the 8-model NURBS corpus every freeform facet
  evaluated and none landed on a pole or seam. The result that matters is
  **ruled vs doubly-curved, not freeform vs analytic**: `bspline_wave`'s
  extrusion faces read p50 0.0507° at 16.7 mm²/facet while `bspline_terrain`'s
  B-spline faces read 0.9396° at 0.53 mm²/facet — the 30× coarser surface is
  19× more exact, because a facet spanning two rulings of an extrusion has an
  exactly perpendicular normal for the same reason it does on a cylinder. So the
  closed-form initializers' exactness is a property of ruledness, and Phase 5's
  freeform work should not expect facet normals to be reliable on
  doubly-curved patches.
  Suite 585 passed / 0 failed / 32 ignored; 11 new unit tests assert the
  initializer identities to machine precision. One precision fact worth
  keeping: **an angle assertion near 0° cannot be tightened past ~1e-6°**
  because `acos` amplifies a 1e-16 rounding error into ~1.5e-8 rad; assert on
  `|â × b̂|` instead, which is linear in the error (measured — an exact strip
  cylinder failed a 1e-9° assertion at 1 sample in 16).

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
  - Torus (**1.0 measured the constraint's residual at 2.354e-5 mm over 23 863
    facets — at the mesh quantum — so this is no longer expected to fail;
    keep the time-box, but the prior is now that it works**): axis from the
    meridian-coplanarity system
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
  **From 1.0: re-initialize from the angularly MOST DISTANT admitted facet, not
  the newest neighbour.** The axis error is `perpendicularity / |n₁ × n₂|`, and
  adjacent facets put `|n₁ × n₂|` at sin(angular step) — 0.043 on a fine mesh,
  a ~23× amplification, worse the finer the mesh. This costs nothing: the patch
  already knows its own facets' normals.
- **1.3 Drift-triggered cylinder→torus promotion** (time-boxed). A bend is
  locally a valid cylinder, so admission never fails; the signal is
  **sign-consistent precession** of the fitted axis across successive refits
  as support accumulates — not a magnitude test (if a magnitude constant
  seems needed, stop and ask). If it works, `promote_bend` + `BEND_*` become
  Phase-4 deletion candidates; if not, keep them and move on.
- **1.4 Validate.**
  - *Initializer accuracy per provenance*: **the instrument is 1.0's
    `--reference-report`, and the baseline it already produced is above.**
    Re-run it on the earned patches and compare against those numbers. The
    provenance fallback this step was written to trigger is **not needed** —
    OCCT Delaunay reads 0.0000° — so what to watch instead is the
    *conditioning* of the pairs the earn loop actually chooses: the error is
    `perpendicularity / |n₁ × n₂|`, so a loop that pairs neighbours inherits a
    ~23× amplification that has nothing to do with provenance. **Stop and ask**
    if the earned initializer underperforms the blind patch anywhere.
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
  **Consider a PER-FACE tolerance instead of a global one** (user's proposal,
  2026-07-30). Measured in DEVELOPMENT_PLAN "What a tessellator actually
  bounds": achieved sag is `min(chordal bound, R·(1−cos θ/2))`, so it varies by
  an order of magnitude *within one model* — 0.08 of the requested deflection on
  a 2 mm fillet against 1.0 on a 12.5 mm sphere. Any single global value is
  therefore simultaneously too loose on small features and too tight on large
  ones. The local sag is already computable per face from the mesh alone
  (`estimate.rs`: `θ·h/8` from the facet's own dihedral and size), so
  `tol(face) = max(f32 position bound, margin × local sag)` needs no new
  measurement machinery — 0.6 supplies the floor, 0.7 the estimate.
  **Note the quantity is sag, not angle**: a 30° turn across a 0.05 mm facet is
  a tiny sag and across a 5 mm facet a large one, so the tolerance must scale
  with `θ × h`, not with θ.
  **This is NOT the 2026-07-27 self-calibrating tolerance that was removed.**
  That one scaled the bar by *each fit's own median residual* — self-referential,
  so a degenerate fit on a small patch loosened its own bar, and the recorded
  post-mortem is "the statistic measured the wrong thing". A mesh-derived local
  sag is independent of the fit being judged, which is the whole difference.
  Still needs the full Phase 3 treatment before landing: audit every dependent
  (`REFIT_GIVE_UP` is 8×stol, merge's centroid gate, `face_validates`), and
  check against the fragility meter, since four of the nine cliffs are
  seed-conditioning constants that a changed bar could move.
  **The dependent audit is DONE** (2026-07-30): all 29 uses read and classified
  in DEVELOPMENT_PLAN "Phase 3.2 dependent audit". Six purposes, three different
  length scales, and `h/sag = 8/θ` — so one constant provably cannot serve the
  two dominant uses. One use looks like a live bug rather than a tuning question:
  the §2.7 freeform fits are held to `stol/8` when their input is mesh vertices
  accurate to ~1e-6 mm. Full harness sweep before/after; every
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
