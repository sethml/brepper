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

- **NEXT MAJOR WORK — exploit the seam assumption. Read this first.** This is
  the recommended plan for the next several sessions, in order (A → B → C).
  Background reading before starting: **Mesh provenance and the seam
  assumption** and **Identifiability is not model selection** in
  DEVELOPMENT_PLAN.md, plus the ranked list in the item below (which records
  two ideas already tried and refuted — do not re-derive them).

  **The framing.** The assumption is that a B-rep face boundary is always a
  chain of mesh edges, tangent junctions included. There are two ways to use
  it, with very different risk: as a **constraint** ("the answer's boundaries
  must be mesh-edge chains" — provenance-independent, true for every faceter)
  or as a **detector** ("here is where the seams are" — provenance-dependent,
  weak on OCCT's Delaunay interiors, information-free on mesh-native input).
  **Stage 2 currently uses it as neither.** Stage 3 uses it fully, and is the
  proof of concept: §3.x "shared-edge topological assembly instead of geometric
  sewing" *deleted* a heuristic-laden alternative outright.

  **The prize, quantified.** Almost every tuned parameter in stage 2 is the
  price of not knowing the tangent seams. If the partition were known, each
  region goes to the whole-region DOF ladder (`classify_one_region`, which
  already exists and whose few constants are cost-control — subsampling, refit
  rounds — not decision thresholds). The evidence is already in hand: **sharp-
  edged models pass the ladder with 0 competition regions and no tuning at
  all.** Everything below exists only because tangent seams are invisible:
  `SEED_PATCH_CAP`, `MAX_PATCH_DEPTH`, `PATCH_RADII`, `PATCH_TRIM_ROUNDS`,
  `MIN_CONDITION_FACES`, `MAX_GROW_ROUNDS`, `GROW_VTOL_FACTOR`,
  `RESIDUAL_SIGMAS`, `MIN_AXIAL_SWEEP_DEG`, `MAX_TORUS_RADIUS_RATIO`,
  `BEND_PATCH_DEPTH`, `BEND_PATCH_CAP`, `MAX_TORUS_STARTS`,
  ~~`LEFTOVER_BLOB_MAX`~~ (deleted 2026-07-28, item 4), `MERGE_MAX_SMALL`,
  `SMALL_SELECTION_ABSORB_FACES`, ~~relabel's 2× margin and 10-round cap~~
  (retired by A, 2026-07-28), and the freeform island rules — roughly **two-thirds of stage 2's behavioural
  thresholds**.

  ---

  **A. Snap shared region boundaries to mesh-edge chains. — DONE 2026-07-28.**
  Implemented as `src/stage2/locus.rs` + `snap_boundaries` (repair.rs), with a
  `dump_regions locus [a b]` diagnostic. Full write-up in DEVELOPMENT_PLAN,
  "Snap shared boundaries to the junction, not to a residual comparison".
  Read that before touching any of it; the summary and the caveats:
  - **Retired relabel's 2× anti-flapping margin and its 10-round cap**, exactly
    as predicted here. Both were the price of comparing a face against two fits
    each contaminated by that face; in-loop refitting removed the
    contamination, and "each face moves at most once" bounds the loop as a
    convergence statement rather than a number. Suite 509 passing / 3 ignored,
    corpus decomposition byte-identical, clippy clean.
  - **`MERGE_MAX_SMALL` was NOT retired** — untested; it is a merge criterion,
    not a boundary one, and nothing this session touched merge.
  - **Snap itself measures as inert on the corpus**: with the margin and cap
    gone, unconstrained relabel reaches the same boundary on its own, and the
    suite passes with snap ablated. It is a guarantee rather than a fix today.
    Decide deliberately whether to keep it — it is the per-face assignment
    criterion item C needs, and the diagnostic that verifies boundaries are
    right (zero wrong-side faces corpus-wide) depends on the same module.
  - **The validation target did not move, and the reason is worth knowing.**
    `stepped_block_toroidal` fine still reports a 4.36e-2 mm `Compare 3.1` edge distance
    and medium still fails at 1.56e-1. Neither is a wrong-side-face problem:
    fine's worst edge is a torus↔cylinder tube junction whose boundary sits on
    the locus to ~1e-6 mm, and medium's is a torus↔B-spline edge, i.e. the
    known decomposition gap. Stage 3.1's compare now names the worst edge by
    its two surfaces and tangency even when it passes, so this is one run to
    re-check rather than an investigation.
  - **Still unhandled, reported loudly if hit:** cone-plane ruling tangency;
    non-coaxial pairs other than the cylinder→torus tube junction (notably
    torus↔torus, which the corpus does not contain); point tangencies. Nothing
    in the corpus fires today, and `boundaries_sit_on_their_junctions`
    (stage2_competition.rs) pins that at zero across seven models along with
    the zero-wrong-side-faces invariant.

  Original plan, kept for the reasoning:

  Every remaining `--compare` failure has been a boundary sitting one
  tessellation row off, and nothing in the pipeline reasons about where a
  boundary *should* be — it falls out of flood plus claim resolution. The
  assumption guarantees the correct boundary **exists** as a mesh-edge chain,
  so this is selection from a guaranteed-nonempty set, not a search.

  - Once two adjacent surfaces are committed, compute the curve where they
    meet, then reassign the boundary faces by which side of it they fall on,
    choosing the mesh-edge chain that best follows it.
  - **The hard case is exactly the one that matters.** At a *tangent* junction
    the surfaces graze rather than cross, so `GeomAPI_IntSS` is the wrong tool
    — it is numerically fragile there and that fragility is what produced the
    original 4.4e-2 mm error (a plane offset δ moves its intersection with a
    torus it should graze by √(2rδ), so tiny fit error is hugely amplified).
    For tangency, find the locus of minimum distance between the two fitted
    surfaces instead — for a plane grazing a torus that is the tangency circle,
    computable in closed form from the fitted parameters. Consider handling
    tangent and transversal junctions by different code paths; stage 3.2
    already classifies edges as tangent, and stage 3's compare already exempts
    tangent edges to the chordal tolerance for this reason.
  - Note OCCT gotchas already recorded in DONE/DEVELOPMENT_PLAN: `IntSS`
    construction can *throw* (knots too close) rather than report `!IsDone`.
  - **Retires:** relabel's 2× anti-flapping margin and 10-round cap, probably
    `MERGE_MAX_SMALL`.
  - **Validation:** `stepped_block_toroidal_fillets_medium` is the target (a
    1.6e-1 mm reference-edge compare, `#[ignore]`d in stage4_integration);
    fine must stay passing. Use `tmp/probe.sh` first, full suite to decide.

  ---

  **Next after A: B.** Note one input B now has that it did not before — stage
  3.1's `--compare` names its worst edge on success, so the chordal-tolerance
  question ("how far off may a tangent boundary be?") can be asked per junction
  rather than per model.

  ---

  **B-prime. Compute the mesh's precision floor and use it as the UNIT, not as
  a knob.** *Slots in ahead of B; the two share a premise and this half is
  exact rather than estimated.* Full write-up in DEVELOPMENT_PLAN, "Mesh
  precision is a float32 quantum, not a constant".

  The corpus's STLs are binary, so vertices are float32 and the floor on vertex
  accuracy is `max|coordinate| × 1.19e-7` — 2.76e-6 mm for rounded_cube,
  **1.96e-5 mm for stepped_block_toroidal**, matching their measured stage-1
  vertex distances. So `--vertex-tolerance` at 5e-5 demands 18× the quantum of
  one model and 2.6× of another.

  - **DO NOT make the acceptance tolerance scale with it — measured wrong.**
    Sweeping `--vertex-tolerance` on stepped_block_toroidal fine gives 138 / 87
    / 148 / 161 / 111 / 76 surfaces at 5e-5 / 1e-4 / 2e-4 / 4e-4 / 8e-4 /
    1.5e-3. Chaotic, no plateau. A scale-adaptive bar lands that model at
    ~3.5e-4, which is worse than today.
  - **DO use it as the unit for fit quality.** Committed analytic surfaces fit
    their support to under ~1 quantum across a 7× range of quanta (medians
    0.31–0.47, maxima 0.6–1.2), while the spurious fits `hill_with_hole` used
    to commit sat at 80×. That is the scale-free commit criterion the
    post-refit tightening idea needs, and it is what "a genuine fit explains
    its support to the precision the data has" actually means.
  - **DO clamp and report.** A `--vertex-tolerance` below the quantum is
    unsatisfiable — the data does not carry that precision. Clamp, warn, and
    print the quantum in the stage-1 summary next to the existing vertex and
    centroid compare numbers.
  - **Not a flag.** The quantum is a measurable property of the input file, not
    a user preference; a flag would only let someone misdeclare it.
    `--vertex-tolerance` stays the user's *request*, the quantum is what the
    data can *support*, and the effective value is the larger of the two.
  - Note the quantum depends on distance from the origin, so the same part
    authored 200 mm off-origin carries ~10× worse data than one centred. Worth
    saying in the README: recentring before export improves precision, and
    nothing downstream can recover it afterwards.

  - **Add an off-origin corpus model so the quantum question is testable.**
    `max(--vertex-tolerance, 2·quantum(p))` — a *per-point* floor, tight near
    the origin and loose far from it — is physically well-motivated but
    **inert on today's corpus**: only `cone_truncated_medium` binds, and
    marginally (2q = 5.8e-5 against a 5e-5 default). On
    `stepped_block_toroidal` fine the whole model spans |p|∞ 64…164, so 2q
    spans 1.5e-5…3.9e-5 and never reaches the default. Generate the same part
    at the origin and translated ~1 m out (`gen_test_models` can do it): there
    2q ≈ 2.4e-4, five times the default, so the floor would actually bind and
    could be *validated* rather than argued about. It would also test the
    "recentre before export" advice directly.
    Two things to get right when it is built: the quantum for a *distance* is
    direction-weighted, `Σᵢ|n̂ᵢ||pᵢ|ε`, not `|p|∞ ε` — on stepped_block the
    y = −164 end-cap plane really does eat the full 1.96e-5 while a fillet
    cylinder measuring radially in xz does not, and a blanket `|p|∞` would
    loosen the latter ~20× more than the data justifies. And least squares
    averages quantization down by ~√N, which is why fits land at 0.3–0.5 q
    while individual vertices need ~2 q of headroom.

  ---

  **B. Estimate `--surface-tolerance` and `--angular-tolerance` from the mesh
  instead of defaulting them.** Both are stand-ins for settings the *exporter*
  used and then discarded; neither is a property of the fitting problem.
  `--vertex-tolerance` is different — it is a genuine precision request from
  the user and should stay a plain flag.

  *Why it is measurable without a reference STEP.* The mesh is its own
  reference, because **vertices lie on the true surface (~1e-6 mm) while
  centroids sit inside it by the chordal sag**. For a facet whose vertices lie
  on a surface of radius R, with circumradius ρ perpendicular to the axis:
  `depth = R − √(R² − ρ²) ≈ ρ²/(2R)`. The tessellator sized the facet so this
  stays under its deflection setting ε, so the sag distribution over the mesh
  *is* a sample of ε.

  *How far off the default is:* stage 1's `--compare` already reports the true
  sag — **8.9e-4 mm** on stepped_block_toroidal fine, **9.5e-4 mm** on rounded_cube
  fine — against a default of **0.4 mm**, i.e. ~400× too loose. (0.4 was a hand
  edit in March 2026 with no recorded rationale; see the tolerance-history
  notes.) Earlier measurements: worst centroid sag 0.026 mm Onshape-coarse,
  0.038 mm Fusion-low.

  Three estimators, increasingly independent:
  1. **Fit-based, bootstrapped from the unambiguous regions.** Sharp-bounded
     regions are solved by the ladder with no tolerance choices beyond vertex
     validation. For each of their facets, distance from centroid to its own
     committed fit is a direct sag sample. High quantile → ε. Available before
     competition runs.
  2. **Fit-free, from the mesh's second-order structure.** For adjacent facets
     on a cylinder of radius R spanning dihedral θ with width h = Rθ:
     `sag ≈ (h/2)²/(2R) = h²/(8R) = h·θ/8`. So **sag ≈ dihedral × facet size /
     8**, per facet, no fits and no jet inversion — computable in stage 2.1.
  3. **Cross-check.** (1) and (2) should agree. Disagreement means the
     provenance premise is violated — a decimated, repaired or remeshed STL —
     which is itself worth detecting and reporting.

  *Prove it before trusting it:* the corpus has reference STEPs and stage 1
  `--compare` reports the true centroid max distance, so run the reference-free
  estimator across the corpus and check it reproduces those numbers. The
  reference is the **validation instrument, never the production input**.

  *The angular tolerance falls out of the same reasoning.* Its own doc comment
  admits it is a stand-in ("Medium-resolution CAD exports step small-radius
  fillets up to ~24° between rows… the default must admit those"), and Onshape
  exposes the setting directly (12.5° / 6.25° / 2.5° for coarse/medium/fine).
  Recover it as a **clustering** rather than a threshold: the dihedral
  distribution should be bimodal — tessellation steps near zero, real sharp
  edges well separated — so find the gap (Otsu-style). Plot it before
  believing it; on stepped_block_toroidal fine only 384 of 51603 edges are called sharp,
  so the populations are very lopsided. **Bonus:** the absence of a clean gap,
  or a spiky constant-dihedral distribution, is exactly the OpenSCAD
  faceted-by-construction detector described in the item below.

  *Risks to measure first:*
  - ε is not truly global — exporters also cap angular deviation and minimum
    facet width, so different features are bound by different constraints. Take
    a high quantile or max, **never a median**, or legitimate faces get
    rejected.
  - Tightening 0.4 → ~1e-3 is a 400× change that moves behaviour everywhere.
    `REFIT_GIVE_UP` is defined as 8×stol; merge's centroid gate and every
    `face_validates` use it. Do the sweep harness (item 6 below) first.
  - Coarse meshes genuinely have larger sag — the corpus's coarse variants are
    the real test, not the fine ones.

  ---

  **C. Reformulate competition as alternating partition/fit.** *The big one —
  where the real simplification lives. Do it last, and only after the sweep
  harness, because it is a rewrite of stage 2.3(b)–2.4 and the before/after
  must be measurable on more than the handful of models one inspects by hand.*

  The assumption says face-level assignment is the right granularity, which is
  exactly the setting for Lloyd-style iteration — this is **Variational Shape
  Approximation** (Cohen-Steiner et al.): alternate (assign faces to proxies)
  and (refit proxies) to a fixpoint. Half of it already exists and is known to
  work here: `relabel_boundaries` + `refit_geometry` converged this session and
  took stepped_block_toroidal fine from 73 surfaces to exactly the reference's 71.

  - **Subsumes:** seed ordering (`neighbor_jet_mismatch`), `seed_patch`,
    `condition_patch`, `grow`/`grow_from`/`flood`, `resolve_claims`,
    `relabel_boundaries`, `merge_regions`, `resplit_disconnected`.
  - **Needs:** an error metric (vertex distance, self-calibrated as in
    `acceptance_tol`); a proxy-count rule (add a proxy where regional error is
    worst, drop one whose faces are all absorbed); and the existing DOF ladder
    to choose each proxy's *type*.
  - **Keep:** the DOF ladder itself, `axis_normals_consistent`, and the
    identifiability guards until measured inert (see below).
  - Beware the known failure mode of greedy selection that VSA is meant to fix:
    a richer model stealing a neighbour's support. Global iteration should
    dissolve it — that is the main hypothesis to test.

  ---

  **After A/B/C: measure how much of the rest can be deleted.** Do not guess —
  this session found three of five guards already inert. Protocol, per
  candidate, one at a time:
  1. Disable it (set the constant to a no-op value, or short-circuit the
     predicate) and rebuild.
  2. `./tmp/probe.sh` — per-model plane/cyl/cone/sphere/torus/bspline
     breakdown over 14 models; diff against a saved baseline. **Surface counts
     alone are not enough: a 3-torus regression once hid behind an unchanged
     total of 71.**
  3. Full suite (`cargo test --release --no-fail-fast`, ~2.5 min) — decisive.
  4. Inert on both ⇒ delete it, and record in DEVELOPMENT_PLAN what it was for
     and what superseded it, so the knowledge is not lost with the code.

  Re-test in roughly this order, most-likely-dead first:
  - `MIN_AXIAL_SWEEP_DEG`, `MAX_TORUS_RADIUS_RATIO` — degeneracy guards that
    exist because competition builds wrong supports; correct supports should
    make them unnecessary.
  - `LEFTOVER_BLOB_MAX` — its rationale is that an analytic surface islanded in
    a jagged per-facet blob makes invalid wires; if boundaries are clean chains
    after A, that rationale may simply evaporate.
  - `GROW_VTOL_FACTOR` + `RESIDUAL_SIGMAS` — both belong to flood, which C
    replaces.
  - The whole `PATCH_*` / `SEED_PATCH_CAP` / `MIN_CONDITION_FACES` family —
    gone with seeds.
  - `MERGE_MAX_SMALL`, `SMALL_SELECTION_ABSORB_FACES` and the freeform island
    rules (`enclosed`, `islands_split_a_patch`).
  - `MIN_CYLINDER_FACES` / `MIN_SPHERE_FACES` / `MIN_TORUS_FACES`.
  - Finally, ask whether the **jet** is still needed at all beyond the
    ladder's rung ordering (`region_is_developable`). It was demoted to a prior
    in the 2026-07-20 pivot; C may remove its last real job.

  Track the constant count as a headline number: stage 2 carried **35** at the
  time of writing (`grep -h "^const \|^pub(crate) const \|^pub(super) const "
  src/stage2/*.rs | wc -l`). Report it in the commit each time it moves.

  *Deliberately not recommended:* seam **detection** proper (finding tangent
  seams directly). Curvature-based detection was measured dead in the
  2026-07-20 pivot — no threshold separates real junctions from real curvature
  gradients, the classes overlap. Structural signals (grid-direction
  discontinuity, density change) were never tried but are provenance-dependent
  and near-absent on OCCT's Delaunay interiors, so a detector would have to
  fall back to competition anyway and would not retire any parameter.

- **The NURBS corpus: 5 of 8 models pass, 3 fail, and the failures were
  invisible.** Surveyed 2026-07-28 against the reference STEPs that were
  already in the repo. **`bspline_hill`, `bspline_terrain`, `bspline_wave`,
  `elliptical_arch`, `elliptical_cylinder` all pass full `--compare` at
  0.000000e0 mm.** Failing: **`hill_with_hole`** (6.73 mm at stage 3.1),
  **`bspline_blister`** (8.04 mm), **`bspline_saddle`** (1.41 mm).

  **The correlation is the finding: every NURBS model with a test past stage 1
  passes, and every one that fails had stage-1-only coverage.** Five of the
  eight were untested past stage 1. These are not hard cases that resisted
  effort, they are unwatched ones — and every one has a reference STEP sitting
  next to it, so the tests could have been written at any time. Two of them
  (`blister`, `saddle`) had silently collapsed to all-per-facet polyhedra;
  that half is now pinned. **Adding `--compare` coverage for the whole NURBS
  corpus is the first thing to do here, before any fix.**

  Two distinct bugs behind the three failures:

  **(a) `hill_with_hole` — competition carves the freeform hill into 40
  spurious analytic patches** (47 surfaces against a reference 7: 5 planes +
  1 B-spline + 1 cylinder). The correct 7 are all present and correct — the
  hole cylinder is r=5.0000 fitted to 1.8e-6 — with 40 cones/spheres of 4–23
  faces and 2–8 mm² area littered over the 1541 mm² hill. `bspline_hill` is
  the *same hill without the hole* and decomposes perfectly, so the pair is a
  controlled reproducer; the difference is that hill_with_hole's mesh is
  coarser (1740 faces for the hill vs 3200), and analytic surfaces fit
  freeform patches happily at coarse density.

  The mechanism that should catch this is the stage-2.3(e) island absorption,
  whose candidacy rule (`SMALL_SELECTION_ABSORB_FACES = 16` plus enclosure) is
  too weak two ways: the accidents run to 23 faces, and because the same test
  decides which surfaces count as "real" for enclosure, the accidents just over
  the cap mark *themselves* real and disqualify their smaller neighbours from
  absorption too — so a whole cluster survives.

  **Four candidacy formulations were tried. Three are measured wrong — do not
  re-derive them:**
  - **Fit quality** (`error_max > --vertex-tolerance`). Tempting, because on
    hill_with_hole the accidents sit 200–800× worse than the genuine surfaces
    (4.0e-4…1.4e-3 against 1.8e-6), with vertex tolerance in the gap. **Wrong
    in general: an accidental fit is frequently an *excellent* fit** — a run of
    coplanar facets along an inflection really is planar, a band of a singly
    curved patch really is cylindrical. `bspline_hill`: 6 surfaces → **2719**.
  - **A fraction of the region.** Fixed both hills; broke `elliptical_arch`
    and `elliptical_cylinder`, because a fraction is *stricter* than a count on
    small regions and looser on large ones.
  - **Smaller than the region's total leftover count.** Fixed hill_with_hole
    (exactly 7 surfaces) and kept the ellipticals; broke `stepped_block`
    medium, 81 → 84 surfaces with a spurious 19th torus — scattered slivers
    along real seams are not a freeform patch.
  - **Smaller than the largest *connected* leftover patch.** Probe clean, all
    eight NURBS models at their best; but it dissolves `cone_tangent_planes`'
    real 16-face cones into a B-spline and re-breaks the test just re-enabled.

  **Why no such rule works, and what it actually needs.** A B-spline can fit
  anything, a cone included — so "the freeform patch validates over these
  faces" is *not* evidence that the region is freeform, and the absorption
  guard has no way to tell a real analytic surface adrift in leftovers from an
  accident. The window the corpus demands is 13.5%–18.6% of the patch, which is
  overfitting, not a criterion. **This needs genuine model selection between a
  B-spline and an analytic surface on the same support.** Note the relation to
  item 2: complexity pricing was refuted there for *nested analytic* types,
  where the DOF difference is 1–2; here it is a B-spline's dozens against a
  cone's six, so the argument that killed it there does not obviously apply and
  it is worth re-deriving for this case specifically.

  **Also learned: `cone_tangent_planes` has been passing for an accidental
  reason.** Its region's freeform fit fails as ill-conditioned, and that is the
  only thing leaving its cones alone. Give the freeform rung a
  better-conditioned domain and it swallows a real cone. Treat that test as
  fragile.

  **And a second, downstream problem on hill_with_hole**: with the decomposition
  corrected to exactly 7 surfaces it clears stage 3.1 and then fails stage 4.1
  on **volume, 6.96e-2 relative**. So candidacy is necessary but not sufficient
  there — expect a freeform-accuracy or export issue behind it.

  **(b) `bspline_blister` and `bspline_saddle` — the freeform fit itself
  fails.** Their large components are rejected by the §2.7 pole-bound check
  (blister: one 926-face component, "control points 2749.5 mm outside the data
  bounding box"; saddle: three components). The check is doing its job — the
  fit really is ill-conditioned — but nothing recovers, so the region falls to
  per-facet. This is where to start: it is a fitting problem with a clear
  signal, not a model-selection argument.

- **REWRITE stage 2.3(b) seeding: kill the blind patch. — NEXT MAJOR WORK on
  stage 2.** Design and all supporting measurements are in DEVELOPMENT_PLAN,
  "Rewriting stage 2.3(b): the blind patch is the root cause". Read it first;
  it records four measured dead ends that must not be re-derived.

  **Why.** `grow_from` derives its initial fit from an *unvalidated* 32-face BFS
  neighbourhood. The patch itself is discarded (the committed face set is always
  a `flood` output, so every member is validated) — but its contamination lives
  on in the fit *parameters*, which produces both observed failures: 323 seeds
  returning NO GROWTH, and plausible-but-wrong fits whose members all genuinely
  validate. Validation is sound; what it validates against is wrong.

  **Step 0 DONE 2026-07-28 — leftover-component ladder rescue (2.3(d)).**
  `rescue_leftover_components` (classify.rs) runs the whole-region DOF ladder on
  each connected component of competition's leftovers, before freeform. No new
  math: the existing conical fitter already estimates its axis from the normal
  covariance (`smallest_eigenvector_3x3`), and given an *uncontaminated* face
  set it recovers r = 0.99998 / 0.99999 / 0.50008 against true 1.0 / 1.0 / 0.5.
  Competition was simply never handing it a clean face set.
  Result on `stepped_block_toroidal` fine: **138 → 108 surfaces**, cylinders
  34 → 36 (reference 37), tori 21 → 19 (18), B-splines 3 → 1, per-facet planes
  80 → 52. Twelve of fourteen probe models bit-identical; suite 508/0/6.
  It can only add surfaces — it runs on faces nothing else explained and every
  rung is vertex-tolerance validated.

  **Step 1 DONE 2026-07-28 — visited ≠ claimed.** Committing a hypothesis no
  longer retires the faces it grew into; every region face seeds exactly once.
  `stepped_block_toroidal` fine **108 → 95 surfaces** (planes 52 → 39), other
  13 probe models bit-identical, suite 508/0/6.
  **Cost, and it is real: that model goes 42 s → 103 s (~2.5×).** Reseeding
  every face unbounded was far worse (>10 min, abandoned), so faces in the
  *interior* of the explained set — every neighbour claimed — are skipped:
  a fit grown from the middle of a validated surface only refits that surface,
  so it is pure cost. Faces at the margin still seed, which is where a
  competitor can actually arise. If the runtime becomes a problem, that
  interior test is the knob, not the principle.

  **Build order for the rest, each step measurable on its own:**
  1. ~~**Visited ≠ claimed.**~~ DONE, see above. A seed that commits nothing costs one entry from the
     untried list; the faces it touched stay seedable. Makes seed misordering
     cost only time (as the plan has always claimed but never delivered) and
     makes an unreliable bootstrap survivable — a 36-face surface gets 36
     attempts. Small, self-contained, and testable alone.
  2. **Identifiability as the commit gate, never residual.** Measured: the three
     spurious tori fit at 0.52/0.72/0.87 quanta against genuine surfaces at mean
     0.39, max 0.56 — overlapping, because a 7-DOF torus on 7 faces is exactly
     determined. Use sweep/conditioning; keep the quantum as the *unit* for
     already-pinned fits.
  3. ~~**Junction-constrained expansion**~~ **DONE 2026-07-28.**
     `grow_with_axis` (compete.rs) + `rescue_with_neighbour_axes` (2.3(d2),
     classify.rs): grow a cylinder whose axis DIRECTION is borrowed from an
     adjacent committed surface, largest validated area first, iterating until
     nothing more grows. With the axis fixed the fit is a 3-DOF circle in the
     perpendicular plane — well conditioned on a short arc — so growth starts
     from the seed facet alone and there is no blind patch to straddle.
     `stepped_block_toroidal` fine **95 → 75 surfaces**: planes 39 → 18
     (reference 16), cylinders 36 → 38 (37), **tori 19 → 18, exactly the
     reference**, B-splines 1. Other 13 probe models bit-identical, suite
     508/0/6.
     Was: **Junction-constrained expansion** from committed neighbours (`locus.rs`).
     Measured on the three cylinders this model loses: fixing the axis alone
     recovers r = 0.99998 / 0.99999 / 0.50008 at 5–10× under vertex tolerance.
     The free fit fails on the axis, not the radius.
  4. **Fixpoint over seed passes**, so pass 2 can use what pass 1 committed.
  5. Only then delete `seed_patch`/`condition_patch` and the `PATCH_*`/
     `SEED_PATCH_CAP`/`MIN_CONDITION_FACES` family.

  **The one real risk**, and it is narrower than "no anchor" (seeding stays
  universal, so with nothing committed the design degrades to today's behaviour,
  which already finds 34/37 cylinders and 18/18 tori): the **coarse-mesh curved
  bootstrap**. A plane is exact from one facet; on fine meshes the jet is good
  (4 of 6 sampled strips read 2.005 against a true 2.0); but the κ₁ = 0.023 /
  0.162 / 2.2 spread over one r=0.5 arc was measured on stepped_block **medium**.
  That is where the blind patch is irreplaceable today. **Validate any
  replacement against stepped_block_toroidal medium and part_rounded_cube coarse
  before deleting anything.**

  **Also re-enables**: the three tests disabled 2026-07-28 (two in
  stage2_competition.rs, one in stage4_integration.rs) all carry this as their
  blocker.

- **stepped_block_toroidal: the seeding is scale-blind — and face size is NOT
  the fix.** Diagnosis in DEVELOPMENT_PLAN, "Why stepped_block_toroidal is hard".
  This is the model every stage-2 threshold sits on a cliff for; it is one
  tangent-connected region of 26122 faces, so all 71 reference surfaces come out
  of competition.

  **The mechanism is established.** 89% of faces are fillet facets covering
  224 mm²; the planes are 2.5% of faces covering 13783 mm². `seed_patch` grows
  a depth-3 BFS neighbourhood capped at `SEED_PATCH_CAP = 32` from surfaces that
  often have fewer than 32 faces (a 3232 mm² plane is **12 faces**; 72 of 75
  connected flat components are under 32). Measured: of 200 sampled flat seeds,
  **zero** produced a patch that stayed on the flat — median 50% fillet, worst
  96%. `PATCH_RADII` is the one spatial bound and is explicitly disabled for
  flats. So the constants are face counts standing in for lengths, which is why
  every sweep of them is non-monotonic.

  **MEASURED DEAD END — do not implement a face-size test.** The obvious fix
  ("a facet much smaller than its neighbour is across a boundary") is not merely
  unreliable, it is **anti-correlated**. Neighbour size ratio (√area):
  *across* a flat|curved junction the median is **1.26** (p90 15); *within* the
  flats it is p75 **18.7**, p90 **47.5**. The discontinuity is inside the
  planes, not at their edges — a 4380 mm² plane spans 0.048…54.36 mm (1100×
  within one surface): one huge interior quad ringed by fine slivers, and of the
  170 faces touching a fillet the median size is 0.118 mm, *finer* than the
  0.52 mm fillet facets they abut. This is constrained-Delaunay grading, exactly
  as predicted by **Mesh provenance** — the plane's triangulation must match the
  fillet's edge discretization along the shared edge. A size-ratio detector
  would fire inside planes and stay silent at junctions. The weaker version
  (keep the seed patch *scale-homogeneous*) dies with it: a 0.118 mm plane
  sliver and a 0.52 mm fillet facet are "comparable".

  **Partially explored, unresolved:** classify by *exact face normals* instead
  of the jet. Motivation: the jet is unreliable here (adjacent jet-"flat" faces
  have a median dihedral of 1.24°, so the jet's flat set is not the planes),
  while facet normals are exact and `curved|curved` dihedrals are a tight 2.44°
  (p90 2.51) against the exporter's angular tolerance. But connected components
  over near-zero-dihedral edges give **342 components above 10 mm² where the
  reference has 16 planes**, and the count is insensitive to the threshold from
  0.001° to 0.2°. The top ~16 sum to 14301 mm² against an expected 13783, so the
  planes are in there with a 326-component tail that is not yet explained.
  **Next step: identify that tail.** Classify mesh faces against the reference
  STEP's 16 planes to get ground truth, then ask what the extra components are —
  large-radius fillets whose dihedral steps are below the threshold
  (`curved|curved` p10 is 0.70°, so some curved edges are much flatter than
  2.5°) is the leading hypothesis. If it is, the separation is between "constant
  dihedral" and "zero dihedral", not between flat and curved, and wants a
  different statistic.

  Diagnostic: `dump_regions <model> jetdump` prints per-face k1/k2/area,
  adjacency and normal — all the measurements above came from it.

- **Faceted-by-construction sources (OpenSCAD, Tinkercad, Blender): a flag or
  auto-detection.** Not a bug today, and not being worried about yet — recorded
  so the eventual decision is informed. These have no B-rep at any stage, so
  there are no tangent junctions to recover, only discretization steps, and the
  seam assumption that carries the whole pipeline holds trivially but conveys
  no information (every edge is a facet boundary). Full findings, with the
  fragment-count arithmetic and the thresholds, are under **Mesh provenance and
  the seam assumption** in DEVELOPMENT_PLAN.md. The headlines:
  - `--angular-tolerance` (default 25°) silently decides whether an OpenSCAD
    revolution comes out as an analytic surface or as the prism it literally
    is: `$fn ≥ 15` reads smooth, `$fn ≤ 14` reads sharp. On OpenSCAD's own
    defaults that is radius ≳ 4.5 mm. A `$fn = 6` hex nut is handled correctly
    by accident of the same rule.
  - The ambiguity is *genuine* — a 15-gon prism and a coarse cylinder are the
    same mesh — so this is a policy choice, not something to be measured. That
    argues for an explicit flag over silent auto-detection, or at least for
    auto-detection that reports what it decided.
  - Expected failure mode where a fit *is* attempted: large radii at default
    resolution, where the sagitta `r(1 − cos(180/n))` breaks
    `--surface-tolerance` (0.55 mm at r = 100, n = 30) and the surface falls to
    per-facet planes. Vertex validation never fails, because OpenSCAD inscribes
    its polygons and the vertices sit exactly on the ideal surface.
  - Detection signals if we go that way: a constant dihedral of exactly 360/n
    all the way around a revolution (B-rep tessellation steps vary, and density
    changes at seams), and vertices exactly on the circumscribed circle.
  - Worth a corpus model or two before deciding anything.

- **Make the stage-2 decomposition robust rather than tuned.** Stage 2 carries
  ~25 behavioural thresholds, and several sit one step from catastrophe: on
  `stepped_block_toroidal` medium, four *unrelated* knobs — `GROW_VTOL_FACTOR`
  ×3.3, `MIN_CONDITION_FACES` 6→4, `PATCH_TRIM_ROUNDS` 16→32, and capping the
  fallback seed patch — each flip the model from 83 surfaces to a full collapse
  to per-facet planes. When unrelated parameters all sit on a cliff the values
  are not the problem, the structure is. (This is the same failure mode the
  2026-07-17 pivot was called for; it re-grew on the new architecture within a
  fortnight, so it is worth watching for deliberately.) Ranked candidates, each
  aimed at *replacing* guards rather than adding one:
  0. **Self-calibrating acceptance tolerance — REMOVED 2026-07-28.** Deleted
     outright (`acceptance_tol`, `residual_scale`, `GROW_VTOL_FACTOR`,
     `RESIDUAL_SIGMAS`); the flood bar is `--vertex-tolerance`. Stage-2
     constant count 36 → 34. Full write-up in DEVELOPMENT_PLAN, "The flood bar
     is the vertex tolerance". Headlines:
     - **It was masking a seed-conditioning defect, not doing its stated job.**
       Its response is non-monotonic (stepped_block_toroidal fine: 138 / 90 / 134 / 113
       / 71 surfaces at 1× / 2× / 5× / 12× / 30×), which is a lever, not a
       bound. And `residual_scale` measured interpolation error over the fit's
       own support and spent it on extrapolation to faces outside it.
     - **The 2026-07-21 evidence for it was stale**: "everything returned NO
       GROWTH" without it is no longer true — 12 of 14 probe models are
       bit-identical at a bar of exactly `--vertex-tolerance`.
     - **Won**: `bspline_blister` 968 → 7 surfaces (reference 6),
       `hill_with_hole` 47 → 8 (reference 7). The spurious analytic patches
       were what had prevented their freeform fits from taking the region.
     - **Cost, now blocking 3 disabled tests**: stepped_block_toroidal fine 71 → 138,
       via 323 no-growth seeds — including flat seeds failing to grow a plane.
       `condition_patch` gives up at `MIN_CONDITION_FACES` without checking the
       seed validates, and stepped_block's flats are narrow strips between
       fillets that a 32-face BFS patch always straddles.
     - **DO NOT retune the floor — measured dead end.** Conical floor 2 gives
       stepped_block 75 with cylinders/tori exactly right (37/18) but hands
       filleted_pipe a spurious 18-face cone at 6e-7; floor 6 is the reverse;
       floor 4 breaks both. A per-type `SurfaceType::CONDITION_FLOOR` does not
       help — the conflict is *within* the conical type. **Fix `seed_patch` so
       the patch cannot straddle, rather than repairing it afterwards.**

  1. **Self-calibrating acceptance tolerance — SUPERSEDED (was DONE 2026-07-27,
     removed 2026-07-28; see item 0)**
     (`acceptance_tol`, compete.rs). See the DEVELOPMENT_PLAN section
     "Calibrate acceptance on the fit, not on a constant". Note what it did and
     did not achieve: the operative bar is now scale-free and insensitive to its
     parameter over a 3× range, but it did not remove a constant —
     `GROW_VTOL_FACTOR` is demoted from operative value to a load-bearing
     safety ceiling and `RESIDUAL_SIGMAS` is added in its place.
  2. **Price model complexity instead of maximizing area — TRIED 2026-07-28,
     PREMISE DOES NOT HOLD. Do not retry as stated.** The idea was that
     `radius_plausible`, `MAX_TORUS_RADIUS_RATIO`, `MIN_AXIAL_SWEEP_DEG` and
     `MIN/MAX_HALF_ANGLE_DEG` are five guards fighting one bias — that the
     surface types nest, so a richer model can always match a simpler one's
     coverage and add a little more — and that an MDL-style score charging for
     parameters would subsume all five. Two measurements killed it.
     *First, the arithmetic.* Coverage pricing cannot reach the case it was
     aimed at. The degenerate torus on `stepped_block_toroidal` covers 180
     faces where the cylinder it displaces covers 94: 86 extra faces for 2
     extra degrees of freedom. Any price that rejects it must charge more than
     ~43 faces per DOF, which would also forbid a cylinder from ever beating a
     plane by fewer than 86 faces — wrong everywhere else. The degenerate fit
     genuinely *is* the better description of its own support; MDL is right and
     unhelpful. The fault is not in the scoring, it is that the support is
     stolen from a neighbour, which greedy per-seed selection cannot see.
     *Second, the guards had already stopped being five.* Ablated individually
     against per-model surface-type breakdowns (`tmp/probe.sh`, 14 models) plus
     the full suite: `radius_plausible` is **completely inert** — deleted, along
     with `hyp_radius` and `conical_local_radius`; `MIN/MAX_HALF_ANGLE_DEG` are
     inert on the corpus but kept as *validity* guards, since a cone with a
     near-zero half-angle has its apex at infinity and cannot be exported
     sanely; `MAX_TORUS_RADIUS_RATIO` still fires (medium 18→19 tori,
     filleted_pipe 4→5); `MIN_AXIAL_SWEEP_DEG` is firmly load-bearing (without
     it stepped_block_toroidal fine goes 16/37/18 planes/cylinders/tori → 12/38/21).
     Most of the pile had already been dissolved by the tolerance and refit
     work rather than by anything aimed at it.
     **What replaced the framing** is the identifiability/validity/selection
     split now in DEVELOPMENT_PLAN.md. The lesson for the next attempt: check
     whether a guard still fires before designing its replacement, and use
     per-type breakdowns — a surface *count* is too coarse a metric and hid a
     3-torus regression here.
     **Methodology note:** `tmp/probe.sh` prints the per-model
     plane/cyl/cone/sphere/torus/bspline breakdown across 14 models and is the
     cheap first check for any stage-2 change; the full suite remains the
     decisive one.

  3. **Refine fitted radii — but NOT by clustering them.** *The problem is
     real; the two obvious answers are both rejected. Open, unassigned.*
     The observation stands: every surface is fitted independently, which
     throws away the strongest prior a CAD model has — that its features
     repeat. `stepped_block_toroidal` has exactly four distinct radii (0.5, 1,
     3, 15) and the pipeline produces a cloud around each (2.9996, 2.9997,
     3.0004, 1.0009, 1.0018).
     - **Clustering the fitted parameters and refitting under a shared value
       is NOT the approach to take** (user's call, 2026-07-28). It treats the
       symptom: the cloud exists because individual fits are being pulled off
       by supports that are still slightly wrong, and averaging a cloud of
       biased estimates gives a well-centred biased estimate. It also needs a
       cluster-width threshold, which is the kind of parameter this whole
       effort is trying to remove, and it silently destroys the case where a
       part legitimately has two nearly-equal radii.
     - **Snapping to round numbers of mm or inches is TABLED**, not rejected —
       a plausible future improvement once radii are accurate enough that
       "extremely close to round" is a safe judgement rather than a guess, and
       it would need a unit inference to go with it. Do not build it now.
     - **What to look for instead: why an individual radius is off at all.**
       The precedent is the 2026-07-26 stale-fit work, which took rounded_cube's
       r=2 fillets from 1.9976–1.9984 to 2.0000 (1000×) purely by refitting on
       corrected supports — no consensus, no snapping, and it improved
       everything downstream rather than just the radii. Treat a spread of
       2.9996/3.0004 as evidence that those two supports still disagree about
       which faces they own, and go find that.
  4. **Make the `LEFTOVER_BLOB_MAX` fallback proportionate — DONE 2026-07-28.**
     The constant is deleted and `dissolve_fragmented_surfaces` (classify.rs)
     replaces it; full write-up in DEVELOPMENT_PLAN, "The leftover blob is a
     boundary problem, not a size problem". Headlines:
     - **Size is the wrong quantity — measured, do not retry a size gate.**
       `cone_tangent_planes` leaves a blob covering 73% of its region and
       reconstructs with zero BRepCheck failures; `bspline_blister` covers 43%
       and produces sixteen invalid faces. Neither an absolute nor a relative
       blob threshold separates them, and nor does "does the decomposition
       explain a minority of the region". Both were tried.
     - **What replaced it**: a face's wire gains one edge per distinct
       per-facet neighbour, so dissolve exactly those surfaces whose per-facet
       neighbour count exceeds their own face count. No length, no tuned
       number, local. Iterated to a fixpoint.
     - **`onshape_cone_tangent_planes_three_cones` re-enabled and passing** —
       its third cone was fitted correctly all along and thrown away with the
       region around it. Stage-2 constant count 37 → 36.
     - **Two corpus models had no coverage past stage 1** and had silently
       collapsed to all-per-facet polyhedra: `bspline_blister` (2137 planar,
       zero analytic) and `bspline_saddle` (1605). They now decompose, pinned
       by `freeform_models_decompose_rather_than_collapsing`. **Look for more
       of this** — a model with only a stage-1 test is a model nobody is
       watching.
     - **Open, tracked by an `#[ignore]`d test**: those two now build 4 and 3
       invalid faces. All are single-mesh-face per-facet planes (plus one
       B-spline) — the micro-sliver family, not the one the discard addressed.
       Blister's root cause is upstream: its 926-face freeform component fails
       the pole-bound check, so nothing can explain it. **Fix the freeform fit
       first.**
  5. **Snap region boundaries to computed intersection curves.** → **SUPERSEDED
     by item A** in "NEXT MAJOR WORK — exploit the seam assumption" above,
     which adds the reason it is sound (the assumption guarantees the answer
     exists as a mesh-edge chain, so this is selection, not search) and the
     warning about tangent junctions, where `IntSS` is the wrong tool.
  6. **Make robustness measurable instead of anecdotal.** *Prerequisite for
     item C above, and advisable before item B.* Guards accumulate because the
     feedback loop is one model at a time, inspected by hand, which is how
     overfitting happens. `gen_test_models` already exists: a sweep over
     randomized radii, feature sizes and tessellation densities, scored
     automatically against the generated references, would say whether a guard
     generalizes or just memorizes the corpus. Note the corpus is Onshape-heavy
     (Parasolid) and so is not representative of the faceter mix in the wild —
     see **Mesh provenance** in DEVELOPMENT_PLAN.md.

- **Re-enable the models disabled in step 7a** (each `#[ignore]`/`ignore =`
  in tests/stage2_integration.rs, stage3_integration.rs, stage4_integration.rs
  carries its blocker string). These are PRE-EXISTING new-pipeline
  classification-quality gaps (steps 2–5), not step-6 regressions.
  **Status 2026-07-26: 17 disabled → 3.** stage3_integration is fully
  re-enabled and stepped_block_toroidal *fine* now passes; what remains is
  cone_tangent_planes + cone_15x20 (stage2) and stepped_block_toroidal
  *medium* (stage4).
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
  - **stepped_block_toroidal FINE PASSES `--compare` (2026-07-26); only the
    medium mesh is still disabled.** The heaviest model in the corpus (51K
    triangles, 26122 fused faces, 71 reference surfaces) tangent-connects into
    ONE region, so every surface has to come out of competition. It used to
    degrade *entirely* to per-facet planes. It now reconstructs to **exactly
    the reference's 71 surfaces — 37 cylinders + 18 tori + 16 planes** — and
    passes the full stage-4 `--compare`. Four fixes, in the order they were
    found; each is a real bug with a measured mechanism, and each unblocked the
    next:
    1. **Degenerate cone on a thin support** (`axis_normals_consistent`,
       conical.rs). An r = 0.5 fillet band one quad tall between two corner
       tori leaves the axis barely constrained by the vertex residual, so the
       least-squares cone tips its axis *into the band's own plane* and beats
       the true cylinder on residual: a bogus 12.2° cone at 3.8e-7 mm against
       the correct r = 0.5003 cylinder at 3.9e-6 mm, its axis nearly
       anti-parallel to the facet normals. `convexity_check` then rejected
       every face of that fit, so the seed grew nothing and the band fell into
       the leftover blob that tripped the whole-region discard. On a cone
       `|n·axis| = sin(half_angle)` holds identically, so the normals pin the
       axis exactly where the positions are degenerate; a fit is now accepted
       only if its own support's *median* facet agrees.
    2. **Seed patch bounded in hops, not in space** (`PATCH_RADII`,
       compete.rs) — three hops walk ~18 mm along a fillet chain swept with
       6 mm facets, straddling several segments. Now also capped at 4 radii of
       curvature of the seed's own jet. Ablation-verified load-bearing.
    3. **Unobserved rotations** (`MIN_AXIAL_SWEEP_DEG`, compete.rs). Ten
       committed surfaces were degenerate super-models of the surface they
       displaced. See the constant's doc for the measurements.
    4. **Two stale-geometry bugs, which is what actually closed it.** Both are
       the same oversight — a surface whose *support* changed kept the fit
       computed from the old support:
       - `resolve_claims` (compete.rs) trims a hypothesis back to the faces it
         won, but the fit still came from everything its flood reached,
         junction rows included — and those rows are exactly what pulls a fit
         off true. Now refitted on what it kept. This is corpus-wide, not
         specific to this model: rounded_cube's r = 2 edge fillets read
         1.9976–1.9984 with no vertex closer than 5e-4 mm to their own
         surface, against a mesh whose vertices are 1.9e-6 mm from the true CAD
         surfaces. They now read 2.0000 at 1.6e-6 mm — a thousandfold. (The
         fitter was innocent: `cylinder_fit_exact_on_partial_arc` pins that
         down.)
       - `relabel_boundaries` (repair.rs) iterated ten rounds but only rebuilt
         face lists between them, never refitting. Boundary repair and fitting
         are one problem: a surface holding a row that belongs to its neighbour
         is pulled off true by exactly that row, and comparing a face against
         two fits both contaminated by it is what stalled the loop inside its
         2× margin. Refitting each round reopens the gap the margin needs.
         Relabel went from 0 moves to 105, and the surface count from 73 to
         exactly 71.
    Unblocking those needed one more fix, in the plane fitter: **`fit_plane`
    area-weights its normal but took a plain vertex mean for the offset**
    (planar.rs). An unweighted mean gives every vertex an equal vote, which
    over-weights whatever the tessellator sampled densely — here the end cap is
    a few large flat quads ringed by the many tiny quads of its perimeter
    fillet's first row, 9.5e-4 mm off the plane and carrying ~1% of the area but
    half the vertices. The plane landed at the *midpoint* of the two levels,
    every one of its 95 faces 4.7e-4 mm out with no outlier group for relabel to
    see. Area-weighting is also just the correct discretization of the
    continuous least-squares plane (∫(n·x − d)² dA).
  - **stepped_block_toroidal MEDIUM is still disabled** — 82 surfaces against
    the reference's 71 (was 96), and stage 3.1's reference-edge compare 1.6e-1
    mm out (was 2.1e-1). The tori are exactly right (18); what is missing is
    ~4 fillet cylinders, and ~14 facets of small fillet arcs survive as
    individual per-facet planes.
    **The improvement came from the seed patch, and it corrected a real
    mis-use of the jet.** `seed_patch` gated its BFS on jet similarity — same
    curvature class, magnitude within 2×. On a small fillet tessellated into a
    handful of facets the jet is simply not reliable enough to gate on:
    adjacent facets of one *exact* r = 0.5 arc read κ₁ = 0.023, 0.162 and 2.2
    (radii of 43 mm, 6 mm and 0.45 mm for the same 0.5 mm fillet). The filter
    then admitted 2–4 faces, too few to determine a 5-DOF cylinder, so every
    candidate failed and the arc degraded to one plane per facet. The jet is a
    *prior, not a gate* — validation decides every face — so the patch now
    falls back to the plain topological neighbourhood when the filtered one is
    too small to fit. Costs nothing in soundness: a patch straddling a junction
    is what `condition_patch` trims and what flood validation rejects.
    **The facet-scaled-growth-tolerance lead is REFUTED.** It was the obvious
    next move — one row past a tangent junction departs by ~h²/2r, which the
    fixed `GROW_VTOL_FACTOR` bar cannot track — but measured on medium the
    response is not monotonic and not a matter of degree: ×1 gives 112
    surfaces, ×3.3 **collapses the model to 2615 per-facet planes**, ×10 gives
    90 with a spurious sphere and a 19th torus. Loosening does not recover
    coarse fillets; it just changes which accidents win.
    **What actually blocks medium now is `LEFTOVER_BLOB_MAX`'s whole-region
    discard, and it is the next thing to attack.** Medium currently sits just
    inside it (largest residual blob 5, threshold 8), and *every* knob tried —
    growth tolerance ×3.3, `MIN_CONDITION_FACES` 6→4, `PATCH_TRIM_ROUNDS`
    16→32, capping the fallback patch at 8/10/14 — pushes the blob over 8 and
    nukes all 2615 faces to per-facet. That is a catastrophic, all-or-nothing
    response to a local defect, and it makes the coarse path untunable: no
    incremental improvement can be measured, because any perturbation lands on
    the cliff. (It is also why several of those knobs *look* load-bearing at
    their current values when they are really just on the safe side of it.)
    The discard earned its place — an analytic surface islanded in a jagged
    per-facet blob reconstructs to invalid wires (rounded_cube coarse: 43
    BRepCheck failures → 0) — so the fix is not to delete it but to make it
    proportionate: scope it to the offending neighbourhood rather than the
    whole region, or gate it on the blob's size *relative* to the region. That
    needs its own session with stage-3/4 validity measured directly, since
    validity is the property it exists to protect.
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
