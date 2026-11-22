# Next steps for brepper (written 2026-11-23)

## Where things stand

- Branch `experiment-fable`, HEAD "Decide cone versus sphere by edge evidence,
  merge same-surface pieces, add --boundary-tolerance". The working tree was
  clean at that commit; this file is new and uncommitted.
- The working plan is the top of `TODO.md` (rules, yardstick, baseline, steps
  1-5). `DEVELOPMENT_PLAN.md` and the "Superseded" half of `TODO.md` are
  history: accurate about the pipeline, but argue for a replaced plan.
- Yardstick: `tools/scoreboard.py`. It judges output against the INPUT MESH
  (self-check), and uses a reference STEP only for face-type labels. The
  baseline, `tests/scoreboard-baseline.tsv`, was run on 2026-11-22 with 4 jobs
  and a 120 s limit:

  | | wild (26) | synthetic (137) |
  |---|---|---|
  | good | 15 | 110 |
  | unfaithful | 0 | 12 |
  | invalid | 2 | 6 |
  | fail | 4 | 7 |
  | timeout | 5 | 2 |

- Suite: 121 unit tests and 378 integration tests pass, 32 ignored
  (`time cargo test --release --no-fail-fast`, about 3.5 min). Clippy is clean.

## Read before working

1. `AGENTS.md` — all of it.
   - The shell is zsh.
   - Prefix builds and tests with `time`.
   - Temporary code goes in `tmp/`.
   - Every commit needs clippy and the build clean, plus updated
     `DEVELOPMENT_PLAN.md` and `README.md`.
   - Commit messages follow the format it gives: model, provider and year,
     and the user prompts.
   - If an OCCT API you need isn't bound in opencascade-sys, complain and
     STOP.
2. `TODO.md`, from the top through "Steps".
3. `README.md`, the sections "Self-check", "Boundary tolerance" and
   "Scoreboard" (which covers the probes).

## Rules the user cares about

- **No new tolerance-like or threshold-like constant without the user's
  sign-off** (TODO rule 4). Reusing existing constants is fine; say which ones.
  The last sign-off was `--boundary-tolerance` (default 1e-3 mm).
- **State measurements as facts and inferences as inferences** (TODO rule 3).
  The user challenges claims that go beyond what was measured, so name the
  test that would separate the alternatives.
- **The user prefers simple remedies.** For example, "the user can raise
  `--vertex-tolerance`" beats building machinery for an unproven problem.
- **Ask before a global change of policy.** Raising `--vertex-tolerance`
  globally was measured on 2026-11-21 and rejected: it creates impostor
  surfaces on exact meshes.
- **Commit only when asked.**

## Next steps, in suggested order

### 1. TODO step 4: freeform cost and runtime bounds (the timeouts)

The five wild timeouts:

- Freeform-heavy parts: `downloaded/3DBenchy` (one 139k-facet freeform
  region), `ToolChanger__blower_fan_duct_(m)`,
  `Open-Task-Light__penultimate_link`, `…drag_chain_termination_(base)`.
- `ToolChanger__docking_station_(m)`, which is different. Stage 3.3 hangs
  inside OCCT: `GeomAPI_IntSS` → `IntPatch_WLine::ComputeVertexParameters`,
  reached from `stage3::compute_edge_curve`. It was sampled with
  `sample <pid>` after 10 min of single-core CPU. Before the boundary-tolerance
  commit this part failed at stage 3.1 (in 47 s); it now gets past 3.1.

Also timing out: synthetic `nurbs/bspline_hill` and `nurbs/bspline_saddle`.

What step 4 asks for:
- A stage-2 time budget per region: whatever isn't classified within the
  budget stays faceted.
- Budgeted seeding, and an early "this region is freeform" exit.
- Faceted output by default for freeform regions, with B-spline optional.
- A bound on stage-3 OCCT intersections.

The intersection bound is the most self-contained piece. Options:
- Run `IntSS` under a time or iteration limit, if OCCT exposes one through the
  bindings; check `../opencascade-rs/crates/opencascade-sys/generated/` first.
- Detect the pair types (quadric–quadric ALine cases) that hang.
- Fall back to a mesh-derived edge.

`blower_fan_duct` (not `_(m)`) is borderline. It took 119 s of the 120 s
limit and then failed at stage 3.1; a faster stage 2 would turn its timeout
into a stage-3 failure.

### 2. Stage-3.1 adjacency failures

There are 4 wild failures and 7 synthetic. All are "unsupported winding-loop /
pole / seam layout" errors:

- `led_cover_drag_chain_termination` (face 20) and `blower_fan_holder_arm_(m)`
  (face 26): "single winding loop but the surface has no pole on the region
  side".
- `blower_fan_holder_arm` (face 27): "V-wrapped face must have exactly two
  winding loops".
- `blower_fan_duct` (face 124).
- Synthetic: several stepped-block, rounded-cube and cone-tangent coarse
  models, and `random/rand_torus_1`.

TODO step 2 records that every wild stage-3 hard failure, measured on
2026-11-15, sat on a surface stage 2 typed wrong. Re-check that against the
current output. Trace a stage-3 face index back to the stage-2 decision with:

```
./target/release/brepper <stl> --compare <ref.step> --compare-labels-only --reference-report --stage 2 -vv
```

This lists every surface whose type the reference disagrees with, by the
surface index stage 3 uses, and every reference face that was split. If the
failing faces are correctly typed, a per-face fallback in stage 3 (the
dropped step 2) becomes worth building.

### 3. Remaining invalid outputs

- `ToolChanger__docking_station` declares a tolerance of 7.7 mm.
- `White-Knight__1-5i_case_front_v4` declares 4.7e5 mm.

Both fail `BRepCheck`.

### 4. Face-type mismatches on parts that are now `good`

- **The three VoronUsers strip mounts:** each has 3 B-splines where the
  reference has cylinders, and more planes than the reference (mini mount 16
  against 8).
  - Measured: two 2.8-2.9 mm² counterbore cylinders come out as one B-spline
    each.
  - Inference: every facet there straddles a tolerant edge, so no surface is
    found at the strict tolerance and the boundary join has nothing to join.
    It is the same probe-based inference recorded in TODO.
  - Possible fix: when every facet of a small component straddles, discover it
    at the boundary bound. That is a NEW use of a tolerance, so it needs the
    user's sign-off. The 2026-11-21 global-1e-3 run showed what loose discovery
    does to such counterbores: they came out as tori.
  - The 18-facet B-spline on the mini mount's 743 mm² cylinder is
    unexplained.
  - The extra planes are not investigated. Coverage is 0.996 and purity 1.0,
    so they are probably small fragments; check this.
- **`Open-Task-Light__drag_chain_link`:** 13 cones where the reference has 4,
  and 1 B-spline (type agreement 0.975).
- The synthetic models with low type agreement are listed in the baseline TSV
  (column `type_agree`). Among them: `fusion_reexport/*_coarse`,
  `onshape/cone_15x20_medium`, `part_rounded_cube_10_r2_fusion_medium`.

### 5. TODO step 5: shrink the docs

Archive `DEVELOPMENT_PLAN` / `LEARNINGS` / `SURFACE_FITTING` and the
superseded half of `TODO.md`. Keep a short design document and a dead-ends
list.

## How the boundary tolerance works (so you don't undo it by accident)

- `--boundary-tolerance` (`config.rs`): read it through
  `Config::boundary_tolerance()`, which is never below the vertex tolerance.
  `Config::boundary_join()` says whether it is looser.
- **Stage 2.3, `classify::join_boundary_facets`.** It runs in each competed
  region, inside the rescue loop and BEFORE `rescue_leftover_components`; a
  ring of straddling facets is exactly a cone frustum, and the rescue would
  commit it.
  - A leftover facet sharing an edge with a surface joins the best-fitting
    one, if it passes `compete::face_validates` at the boundary bound
    (orientation, centroid within `--surface-tolerance`).
  - Rounds repeat until no facet joins.
- **Stage 2.4, `repair::merge_regions`.** A small surface (at most
  `MERGE_MAX_SMALL` = 16 facets) whose facets all lie within the bound of a
  larger neighbour merges into it. Competition commits 3-8 facet cones along
  tolerant edges before the join can run.
- **Both places share two guards. Each was added after a measured
  regression:**
  - Only a surface its own strict facets determine can take joined facets: at
    least `CompeteHyp::min_support` (2 for a plane, then `MIN_CYLINDER_FACES`,
    `MIN_SPHERE_FACES`, `MIN_TORUS_FACES`). Without it, a one-facet
    competition "cylinder" took 69 freeform facets on `blower_fan_duct`.
  - A vertex beyond the vertex tolerance must be on a real edge: some facet
    around it (`stage2::vertex_fan`) is across a sharp edge (another stage-2.2
    region) or belongs to another determined surface. Leftover facets do NOT
    count. Without this, `blower_fan_duct`'s single-facet planes went from 767
    to 2039.
- **Joined facets** are flagged in `Classification::joined`. Fits never use
  them (`repair::refit_geometry`, `classify::fit_validate_joined`), and
  revalidation holds them to the boundary bound. `retype_sphere` leaves them
  out of the question.
- **The self-check** (`selfcheck.rs`) allows the boundary bound only at
  vertices shared by facets of two different surfaces.
- **Turning it off:** `--boundary-tolerance=5e-5`, equal to the default
  vertex tolerance.

## Practical information

- **Build:** `time cargo build --release` (the first OCCT build takes about
  8 min; after that, seconds). Binary: `target/release/brepper <stl> -o
  <out.step> -v`.
- **Scoreboard:**
  - `time python3 tools/scoreboard.py --diff --out tmp/sb.tsv` takes about
    1 hour awake for 163 models with 4 jobs.
  - `--filter REGEX` narrows the run.
  - `--diff-only A.tsv B.tsv` compares two saved runs.
  - `--baseline` saves a run as the baseline.
  - Run it in the background, and don't run the test suite at the same time:
    the extra load pushes borderline models past the 120 s limit.
- **Machine sleep:** the user's Mac sleeps when the lid is closed. On the
  night of 2026-11-21 that stretched runs from minutes to hours; `pmset -g
  log` shows "Clamshell Sleep" and 45 s DarkWakes every 15 min.
  - The user just installed Amphetamine Enhancer; whether closed-lid mode is
    active is unverified. If wall-clock time far exceeds CPU time, check
    `pmset -g log` before suspecting the code.
  - The scoreboard now times on the monotonic clock, which stops during
    sleep, so its `secs` and timeouts stay meaningful.
- **Probes** (stderr, or stdout where noted):
  - `--reference-report -vv`: type disagreements and split reference faces.
  - `BREPPER_TYPE_PROBE=1`: sphere and cylinder-or-cone validation and edge
    evidence for each reference face (`TYPEPROBE` rows, on stdout).
  - `BREPPER_FACET_PROBE=1`: vertices off every incident reference surface,
    straddling facets, and fragments; `=2` prints each straddling facet.
  - `BREPPER_DEBUG_JOIN=1`: each region's surfaces with their joined counts,
    and every facet the join left out with its distance to each adjacent
    surface.
  - `BREPPER_DEBUG_COMPETE=1`: per-region competition, leftovers and freeform
    summary, and the sphere-or-conical verdicts.
- **Stuck process:** `sample <pid> 3 -file tmp/x.txt` shows where it is stuck.
- **Corpus:**
  - `tests/external` holds 24 verified single-body wild pairs.
    `tests/external/README.md` covers provenance, which exporter wrote each
    one (SolidWorks 2021 or Fusion 360), alignment, and the stage-1.2 boundary
    rule.
  - `tests/downloaded` has 2 meshes.
  - `tests/external/rejected.tsv` lists the 29 removed pairs and why.
  - Synthetic models are under `tests/{ccad,onshape,fusion,fusion_reexport,remesh,random,nurbs,occt,manual}`.
- **Gaps to watch:**
  - The scoreboard counts "types match" over the good parts only.
  - Coverage and purity count a kernel's two half-cylinder faces as two
    labels; type agreement is the number to steer by.
