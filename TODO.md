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
  `LEFTOVER_BLOB_MAX`, `MERGE_MAX_SMALL`, `SMALL_SELECTION_ABSORB_FACES`,
  relabel's 2× margin and 10-round cap, and the freeform island rules —
  roughly **two-thirds of stage 2's behavioural thresholds**.

  ---

  **A. Snap shared region boundaries to mesh-edge chains.** *Start here:
  narrow, provenance-independent, and it removes a heuristic rather than adding
  one.*

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
  sag — **8.9e-4 mm** on stepped_block fine, **9.5e-4 mm** on rounded_cube
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
  believing it; on stepped_block fine only 384 of 51603 edges are called sharp,
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
  took stepped_block fine from 73 surfaces to exactly the reference's 71.

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
  1. **Self-calibrating acceptance tolerance — DONE 2026-07-27**
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
     it stepped_block fine goes 16/37/18 planes/cylinders/tori → 12/38/21).
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

  3. **Global parameter consensus.** *Still open and still worth doing; it is
     independent of A/B/C above and could be slotted in anywhere.* Every
     surface is fitted independently,
     which throws away the strongest prior a CAD model has — that its features
     repeat. `stepped_block_toroidal` has exactly four distinct radii (0.5, 1,
     3, 15) and the pipeline produces a cloud around each (2.9996, 2.9997,
     3.0004, 1.0009, 1.0018). Cluster the fitted parameters, recover the shared
     radii and axes, and refit under those constraints. Should improve accuracy
     *and* remove spurious surfaces that belong to no cluster; it is what
     commercial reverse-engineering tools do as "beautification".
  4. **Make the `LEFTOVER_BLOB_MAX` fallback proportionate.** It is the cliff
     above: a local defect triggers a global catastrophe, which makes coarse
     meshes untunable because no incremental improvement can be measured. Scope
     it to the offending neighbourhood, or gate it on the blob's size relative
     to the region. It earned its place (rounded_cube coarse: 43 BRepCheck
     failures → 0), so measure stage-3/4 *validity* directly when changing it —
     that is the property it protects. Needs stage-3 robustness at
     analytic/per-facet boundaries first.
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
