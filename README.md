# Brepper - STL to STEP Converter

Convert triangulated STL meshes from CAD exports into parametric STEP files with fitted analytic and freeform surfaces.

## Usage

```
brepper - Convert STL mesh to STEP with fitted surfaces

USAGE:
    brepper [OPTIONS] <input.stl> [-o <output.step>]

REQUIRED:
    <input.stl>              Input STL file (binary or ASCII)

GENERAL OPTIONS:
    -o <step>, --output=<step>  Output STEP file
    --stl-units=<units>      Units used by the STL file (mm, cm, m, in, ft, um; default: mm)
    --step-units=<units>     Units to use in exported STEP file (default: mm)
    --compare=<step>         STEP file to compare to at each step
    --reference-report       Score the mesh and our segmentation against the --compare
                             reference: facet-normal exactness per surface type, closed-form
                             initializer error, the seam assumption, and over/under-segmentation.
                             Diagnostics only. Requires --compare.
    --compare-labels-only    Use the --compare reference for labels only: the stage-1.2 check
                             that the mesh lies on it, and --reference-report's scoring. The
                             stage 3/4 gates that fail a run whose edges, faces or volume
                             differ from the reference's are skipped. tools/scoreboard.py
                             runs in this mode.
    --no-self-check          Skip the end-of-run self-check (see "Self-check" below).
    --invariant-report       Score the DEVELOPMENT_PLAN "Invariants" against the --compare
                             reference: whether the fitting requirement is satisfiable, whether
                             an angular tolerance can separate adjacent surfaces, whether a wrong
                             surface type also explains a true support, and whether growth from
                             different seeds converges on the same surface. Runs after stage 2.2.
                             Diagnostics only. Requires --compare.
    --vertex-tolerance=<value>  Fitting tolerance in STL units (default: 5e-5)
    --surface-tolerance=<value> Surface-to-face offset tolerance in STL units (default: 0.4)
    --angular-tolerance=<deg>   Max dihedral angle between adjacent triangles on the same surface (default: 25)
    -v, --verbose[=LEVEL]    Verbosity level: 1 = stage summaries, 2 = per-hypothesis details,
                             3 = commit details/stats/profiles, 4 = per-seed details, 5 = BFS face trace
    -q, --quiet              Suppress non-error output
    --threads=<N>            Worker threads for stage 2's fit competition (default: all
                             cores; 1 = serial). Output is identical for every value.
    --strict-mesh            Fail on mesh defects instead of repairing them (dropping
                             duplicate-vertex and coincident faces, stitching T-junction
                             slivers, unifying winding, reorienting inverted shells)
    --debug                  Enable debug output and intermediate files
    --stage=<stage>          Stop after stage, e.g. 2.2 (default: 4.1). --stage=2 stops after all of stage2
                             Stage 2 sub-stages (partition-then-classify rewrite): 2.1 sharp-edge
                             classification + jet prior, 2.2 smooth-region extraction, 2.3 per-region
                             surface ladder, 2.4 resolution & repair, 2.5 selected surfaces
    --viz-stages=<stages>    Open interactive 3D visualization at specified stages (e.g. 2.1,2.2)
    --viz-stages-level=<N>   Viz detail level: 1 = committed hypotheses (default), 2 = region colors
    --viz-seeds              Interactively visualize how SEED SETS grow, against the --compare
                             ground truth. Stops on each facet for which every conforming seed
                             set fails, and names the validation clause that refused it.
                             Requires --compare. See "Seed-growth visualization" below.
    --viz-seed-show=<what>   Which seed sets to stop on: failing (default), all, covering
    --viz-seed-sizes=<list>  Facet counts the seed generators may use; only the maximum
                             is read, as the bound on chain length (default 2,3,4 — the
                             sizes the DEVELOPMENT_PLAN Invariants section specifies)
    --viz-seed-samples=<N>   Anchor facets sampled per reference face (default 12)
    --viz-seed-max-attempts=<N>  Conforming seed sets tried per anchor facet (default 12)
    --viz-seed-spec-only     Try only the seed patterns the Invariants section names, without
                             the grown-ball fallback — shows where those rules alone fall short
    --viz-seed-coplanar-clause   Enforce the section's "no more than three coplanar" clause
                             (measured: it rejects every seed containing a fused quad)
    --viz-seed-confine-region    Confine seed growth to the stage-2.1 smooth region instead of
                             letting validation be the only limit
```

### Examples

```bash
# Full conversion
./target/debug/brepper input.stl -o output.step

# Just load and validate mesh (stage 1.2)
./target/debug/brepper input.stl --stage 1.2 -v

# Compare mesh against reference STEP at each stage
./target/debug/brepper input.stl --compare reference.step --stage 1 -v

# Score the mesh against ground truth read from the reference (see src/reference.rs)
./target/debug/brepper input.stl --compare reference.step --stage 1 --reference-report -v

# Score the Invariants section's clauses against that ground truth
./target/release/brepper input.stl --compare reference.step --stage 2.2 --invariant-report -v

# Watch which seed sets fail to grow to their surface, and why
./target/release/brepper input.stl --compare reference.step --stage 2.2 --viz-seeds -v
```

### Self-check

Every conversion ends by checking the solid it wrote against the mesh it was
given, with no reference involved, and prints one line:

```
Self-check: status=good solids=1 faces=44 valid=true vtx_max=1.249e-6 vtx_off=0 uncovered=0 cover_worst=0.000e0 vol_rel=-1.574e-4 tol_max=2.856e-3
```

- `status=good` — a valid closed solid (`BRepCheck`), every mesh vertex on it
  (`vtx_off=0`), and output under every mesh facet (`uncovered=0`).
- `status=unfaithful` — a valid solid that is not the mesh: vertices off it,
  facets with nothing under them, or a declared tolerance (`tol_max`) above
  `--surface-tolerance`.
- `status=invalid` — no solid, or `BRepCheck` fails.

Anything but `good` exits non-zero. The STEP file is still written, so it can
be inspected, and the message names the worst offending vertices or facets.
`vol_rel` is the solid's volume against the mesh's and is reported, not gated:
a mesh legitimately encloses slightly less than the solid it came from.
`--no-self-check` skips the check. See `src/selfcheck.rs`.

### Scoreboard

```bash
tools/scoreboard.py            # run the corpus, write tmp/scoreboard.tsv, print a summary
tools/scoreboard.py --diff     # ...and diff against tests/scoreboard-baseline.tsv
tools/scoreboard.py --baseline # ...and save the result as that baseline
```

One row per corpus model: the self-check status, runtime, face types out, and —
where the model has a reference — the reference's face types and three label
scores for stage 2 (coverage, purity, and `type_agree`, the share of mesh area
whose surface is the same kind as the reference's). The verdict on geometry
never involves the reference. This is the yardstick for the plan in TODO.md;
`tools/sweep.sh`, which scores runs by per-stage agreement with a reference
STEP, remains for the synthetic suite.

### Reference check (`--compare`, stage 1.2)

With `--compare`, stage 1.2 requires every mesh vertex to be within
`--vertex-tolerance` of the reference's faces and every facet centroid within
`--surface-tolerance`, and prints the vertex-distance p50 / p95 / max.

One class of vertex is not counted: a boundary vertex that lies on both of the
surfaces meeting at its edge and misses only the reference's *trim*, by no more
than that edge's own tolerance. STEP stores an edge as an explicit curve that
approximates the true intersection, so such a vertex can be exact and still read
1e-4 mm from the trimmed faces. When any are excused, a second `Compare:` line
says how many and how far. To see every over-tolerance vertex classified:

```bash
./target/release/gen-test-models strays input.stl reference.step
```

### Reference scoring (`--reference-report`)

With a `--compare` reference, brepper can label every mesh face with the
reference face it came from — recovered from the distance computation the
compare already performs — and then score three things against that ground
truth. All of it is diagnostics: nothing downstream reads it.

At **stage 1**, properties of the mesh alone (so they are still measurable on
models where stage 2 fails):

- **Facet-normal exactness** per reference surface type: the angle between the
  mesh facet normal and the true surface normal. Analytic faces use a closed
  form; B-spline, Bezier, revolution, extrusion and offset faces are evaluated
  through OCCT at the footpoint, and facets at a pole or seam — where the normal
  does not exist — are excluded rather than scored.
- **Closed-form initializer error**: the cylinder axis from `n₁ × n₂`, the
  sphere centre from circumcenter normal lines, and the torus
  meridian-coplanarity residual, each against the reference's own parameters,
  reported with the conditioning that governs them.
- **The seam assumption**: how often a facet's own interior disagrees about
  which reference face it belongs to.

At the end of **stage 2**, our decomposition against the reference's:

- **`ref coverage`** — 1.0 when no reference face was split across several of
  our surfaces (no over-segmentation).
- **`surface purity`** — 1.0 when none of our surfaces spans several reference
  faces (no merging).

Both are area-weighted and threshold-free. They catch what the distance and
volume gates cannot: a 9-surface decomposition of a 26-surface part passes every
existing gate but scores `surface purity 0.2201`.

### Invariant scoring (`--invariant-report`)

Scores the clauses of DEVELOPMENT_PLAN's **Invariants** section against the same
ground truth. Also diagnostics only, and it runs after **stage 2.2** rather than at
the end of stage 2, so it still works on models whose decomposition crashes or
times out. One TSV row per reference face on stdout, a per-model summary on stderr.

```bash
./target/release/brepper input.stl --compare reference.step --stage 2.2 --invariant-report -v
```

Five families, one per clause:

- **F — is the fitting requirement satisfiable?** With the reference's *own* exact
  parameters, how far off it the facets' vertices, edge midpoints, quad diagonals,
  centroids and interiors actually land. This asks whether a clause can be met at
  all before asking whether the algorithm meets it, and it is what showed that
  facet centroids are not a bound on the whole facet while edge midpoints are.
- **A — can an angular tolerance separate adjacent surfaces?** The largest dihedral
  *inside* a reference face against the smallest *across* its boundary, plus the
  count of edges stage 2.1 calls sharp strictly inside a true surface.
- **I — does a wrong surface type also explain a true support?** All four analytic
  types fitted to each reference face's own facets, with residuals in float32
  quanta.
- **S — seed independence.** Growth run from many seeds inside one true surface:
  how many recover it, how many distinct answers come back, and whether a failure
  leaked past the true boundary or stopped short of it.
- **G — growth independence.** The flood re-run with the neighbour iteration
  reversed, and the refit fixpoint restarted from the whole true face set.

Corpus sweep and aggregation:

```bash
./tools/invariants.sh --out tmp/inv.tsv --jobs 8 && python3 tools/inv_analyze.py tmp/inv.tsv
```

A sixth family is off by default because it multiplies the row count:

- **K — the kinematic rank ladder** (`--kinematic`, or `BREPPER_INV_KINEMATIC=1`).
  Scores the proposal in SURFACE_FITTING.md, "Nested fitting: a rank ladder
  instead of a type ladder", which reads the surface type from the *null-space
  dimension* of the tangent-velocity-field system so that too little support
  reads as a simpler surface rather than a wrong one. Per reference face, over a
  patch grown from several starts at support sizes 3…48: the spectrum, the
  dimension and the type it names, whether that matches the truth, the straddling
  fraction and the number of distinct axial bands in the support, and the
  residual of the reference's *own* tangent generators — the family-F question
  asked of rank rather than of distance — at both candidate evaluation points.

```bash
./tools/invariants.sh --out tmp/kin.tsv --jobs 8 --kinematic && python3 tools/kin_analyze.py tmp/kin.tsv
```

And a seventh, also off by default:

- **J — the junction race** (`--junction`, or `BREPPER_INV_JUNCTION=1`). When a
  patch grown on a surface crosses a tangent junction, does the kinematic rank
  lose a tangent field before the surface's own refit stops validating? Per
  analytic reference face with a smooth junction: unconfined growth from facets
  2, 4, 8 and 16 hops from the junction, recording at each ladder size whether
  the patch has crossed, whether the true type still validates at
  `--vertex-tolerance`, and whether the rank has dropped.

```bash
./tools/invariants.sh --out tmp/junction.tsv --jobs 8 --junction && python3 tools/junction_analyze.py tmp/junction.tsv
```

`tools/invariants.sh` reuses `tools/sweep.sh`'s model/reference pairing, so the two
instruments always cover the same corpus.

Env-gated probes on the same path, for when a single fit misbehaves:
`BREPPER_INV_FITDUMP=1` prints the reference's own parameters beside the fitted
ones for any seed that fails admission, with a verdict on whether the failure is
non-identifiability or non-convergence; `BREPPER_INV_BASIN=1` sweeps the basin of
attraction of grow-with-refit; `BREPPER_INV_PATCHSIZE=1` measures fit accuracy
against support size; `BREPPER_INV_REFDUMP=1` prints every reference face's exact
parameters to stderr; `BREPPER_INV_FITVERBOSE=N` passes verbosity `N` to the
fitters that score each reference face's own support, so a `NaN` in the `res_*`
columns names the guard that refused it.

### Seed-growth visualization (`--viz-seeds`)

Opens the 3D window and walks **anchor facets**, stopping on each one for which
*every* conforming seed set fails to grow to its surface. Requires `--compare`.

```bash
./target/release/brepper tests/ccad/generated/simple_cone.stl   --compare tests/ccad/generated/simple_cone.step --stage 2.2   --viz-seeds --viz-seed-spec-only --viz-seed-samples 12 -v
```

| colour | meaning |
|--------|---------|
| white | the anchor facet |
| green | the rest of the seed set |
| blue | correctly covered |
| red | overshoot — grown past the true boundary |
| orange | true facets that **blocked** further growth |
| yellow | true facets missed for some other reason |

Controls are the same as `--viz-stages`: `Space` next seed set, `Shift+Space` next
reference face, `Q` quit. The status line names the seed shape, how far the seed
sits from its own fit, and the diagnosis — which validation clause refused what,
since `face_validates` is a conjunction of vertex distance, centroid distance and
orientation coherence and naming only the first is misleading.

`BREPPER_VIZ_SEEDS_TEXT=1` prints the same diagnoses with no window, which is what
the corpus catalogues are built from.

See [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) for full documentation.

### Tolerance Settings

#### Surface Tolerance

This is how far the surface of a mesh triangle is allowed to deviate from the surface it is approximating.

Onshape calls this "Chordal tolerance (mm)", and the preset values are:
- Coarse: 0.24 mm
- Medium: 0.12 mm
- Fine: 0.06 mm

Fusion 360 calls this "Surface Deviation", and the preset values seem to vary based on model size:
- Low: 1x
- Medium: 2.5x
- High: 8x

#### Angular Tolerance

This is the maximum angle between adjacent triangle faces that approximate a smooth surface.

Onshape calls this "Angular deviation (deg)", and the present values are:
- Coarse: 12.5 deg
- Medium: 6.25 deg
- Fine: 2.5 deg

Fusion 360 calls this "Normal Deviation", and the preset values are:
- Low: 30.0
- Medium: 15.0
- High: 10.0
In theory "Normal Deviation" should be the maximum angular error of the face relative to the surface normal, but in practice it seems to result in the given maximum angular deviation.

### OCCT-Generated Test Models

Models in `tests/occt/` and `tests/nurbs/` are generated directly via OCCT by
the `gen-test-models` binary (`./scripts/generate_occt_models.sh`) — used for
models that ccad cannot express. `tests/occt/` holds analytic models such as
`hollow_cube` (a solid with a fully enclosed internal void, for multi-shell
solid support). `tests/nurbs/` holds freeform models (exact B-spline patches,
conic/B-spline extrusions) exercising the stage 2.7 NURBS fallback; all 8
pass full `--compare` — see `tests/nurbs/README.md`. The generated STL/STEP pairs are committed.

### Generating Test Models

Some tests use models generated via [CodeCAD](https://codecad.xyz/). The generated STL/STEP files are committed to the repo, so CodeCAD is not required to run tests.

To add or modify these models:
1. Install CodeCAD
2. Edit files in `tests/ccad/parts/` — see `tests/ccad/README.md` for primitive conventions and placement rules
3. Run `./tests/ccad/generate_models.sh`

**Important:** `box(w,d,h)` places a corner at the origin, not centered. Use `center_xyz()` to center. Each part file includes analytical surface area and volume in its header comments.

### Interactive Debugging Visualization

Use `--viz-stages` to open a 3D window and step through BFS hypothesis deduction interactively.
Use `--viz-stages-level` to control how often the visualization pauses:
- Level 1 (default): Pauses at each committed hypothesis
- Level 2: Also pauses at each accepted seed before BFS expansion
- Level 3: Also pauses at each face accepted during BFS expansion

```bash
# Visualize planar and cylindrical hypothesis stages (committed hypotheses only)
brepper input.stl --viz-stages=2.1,2.2

# Visualize with seed-level detail
brepper input.stl --viz-stages=2.2 --viz-stages-level=2

# Visualize all stage 2 hypothesis types with STEP comparison overlay
brepper input.stl --compare reference.step --viz-stages=2.1,2.2,2.3
```

At each BFS step, the window shows the STL mesh with highlighted faces: seed faces in green, previously accepted faces in blue, the just-accepted face in dark blue, and BFS queue frontier edges in red. Cylindrical and spherical hypotheses also display translucent geometry overlays. When `--compare` is provided, the STEP surface is shown at 35% opacity for reference.

| Action | Control |
|--------|:--------|
| Next step | `Space` |
| Skip to next seed/hypothesis | `Shift+Space` |
| Skip to face seed | `Click` on a face |
| Quit visualization | `Q` |

Supported stages: `2.1` (sharp edges + jet prior), `2.2` (smooth regions), `3.3` (edge curves), `3.4` (faces), `3.5` (shells), `3.6` (solids). (The old per-type BFS viz stages are retired with the stage-2 rewrite.)

`--viz-seeds` opens the same window on a different question — how seed sets grow rather than how the pipeline proceeds — and is documented under [Seed-growth visualization](#seed-growth-visualization---viz-seeds) above.


## Viewer

`viewer` is an interactive 3D viewer for STL and STEP files.

```
USAGE:
    viewer [OPTIONS] <file>...

ARGUMENTS:
    <file>...   One or more STL or STEP files to view (detected by extension)

OPTIONS:
    --deflection <value>          Linear deflection for STEP tessellation (default: 0.1)
    --angular-deflection <value>  Angular deflection for STEP tessellation in degrees (default: 28.6)
```

Multiple files can be loaded simultaneously; each is rendered in a distinct color.

### Examples

```bash
# View a single STEP file
./target/release/viewer part.step

# View a STEP and an STL side by side
./target/release/viewer reference.step output.stl

# Finer tessellation for curved surfaces
./target/release/viewer part.step --deflection 0.01 --angular-deflection 5.7
```

### Controls

| Action | Control |
|--------|:--------|
| Rotate (around model center) | Left-click + drag |
| Pan | Right-click + drag |
| Zoom (at cursor) | Scroll wheel |
| Quit | `q` or `Esc` |
| Toggle perspective/orthographic | `p` |
| Toggle wireframe edges | `e` |
| Toggle soft (mesh) edges | `Shift+E` |
| Toggle solid faces | `s` |
| Toggle hidden-edge removal | `h` |

Rendering uses PBR shading with ambient and two directional lights. Sharp feature edges (boundary edges and edges where the dihedral angle exceeds ~45°) are drawn as a solid wireframe overlay. Softer mesh edges (dihedral angle below ~45°) are drawn at 50% opacity and can be toggled independently with `Shift+E`. Zooming targets the point under the mouse cursor.

### Running Tests

```bash
# Run all tests (unit + integration)
cargo test

# Run only integration tests
cargo test --test integration stage1_integration::
```

Integration tests in `tests/integration/stage1_integration.rs` (one test binary, `tests/integration/main.rs`) verify that:
- All STL files in `tests/manual/`, `tests/onshape/`, `tests/ccad/generated/`, and `tests/fusion/` pass stage 1 mesh validation and coplanar triangle fusion
- All STL/STEP pairs pass `--compare` surface distance checks (after fusion)
- Bad files in `tests/bad/` fail with expected error types (degenerate faces, non-manifold edges, inconsistent winding, compare failures)

Integration tests in `tests/integration/stage2_integration.rs` (one test binary, `tests/integration/main.rs`) verify that:
- Planar hypothesis deduction produces correct counts for known planar models (cube=6, wedge=6, chamfered cube=26)
- Cylindrical hypothesis deduction produces correct counts and properties (simple_cylinder=1 convex, block_with_hole=1 concave, pipe=2, two_holes=2 concave)
- Spherical hypothesis deduction produces correct counts and properties (simple_sphere=1 convex r=10, hemisphere=1 convex r=10, spherical_pocket=1 concave r=8, ball_on_cylinder=1 sphere+1 cylinder)
- Surface selection (stage 2.6) produces correct surface type counts for known models (block_with_hole=6 planar+1 cylindrical, pipe=2 planar+2 cylindrical, ball_on_cylinder=1 planar+1 cylindrical+1 spherical, hemisphere=1 planar+1 spherical, simple_sphere=1 spherical, cone_15x20=1 planar+1 conical, simple_cone=2 planar+1 conical)
- All STL/STEP pairs pass `--compare` validation at stages 2.1, 2.2, 2.3, and 2.6
- Deduced cylinder parameters (axis direction, axis position, radius) match STEP cylinder parameters for all cylindrical test models
- Deduced sphere parameters (center, radius) match STEP sphere parameters for all spherical test models
- Angular tolerance enforcement prevents spurious cylindrical/spherical hypotheses on planar-only models (e.g., cube faces meet at 90°, exceeding the 25° default)
Integration tests in `tests/integration/stage3_integration.rs` (one test binary, `tests/integration/main.rs`) verify that:
- Adjacency graph topology (faces, edges, vertices) is correct for known models (cube=6F/12E/8V, cylinder=3F/2E/0V, sphere=1F/0E/0V, etc.)
- Euler's formula V-E+F=2 holds for genus-0 models
- Edge validity, adjacency symmetry, and edge-face consistency are maintained
- All STL/STEP pairs pass `--compare` validation at stage 3.1, verifying that reconstructed edge boundary vertices lie on STEP edges and BRep vertices coincide with STEP vertices
- Tangency detection (stage 3.2) correctly identifies no tangent edges for models composed of planar, cylindrical, and spherical surfaces meeting at angles > 2°
- Edge curve computation (stage 3.3) produces 3D curves for all edges via surface-surface intersection, verified across all test models (cube, cylinder, hemisphere, ball_on_cylinder, block_with_hole, chamfered_cube, spherical_pocket)
- Face creation (stage 3.4) successfully creates OCCT `TopoDS_Face` objects for all test models via shared-edge topological assembly: every edge is one `TopoDS_Edge` shared by both adjacent faces, with per-face pcurves, seam edges, and degenerate pole edges built uniformly for all surface types. All STL/STEP pairs pass `--compare` validation at stage 3.4
- Shell construction (stage 3.5) assembles faces into shells directly from the shared-edge topology (closed by construction; no sewing) and produces correctly oriented shells. All STL/STEP pairs pass `--compare` validation
- Solid construction (stage 3.6) creates solids from shells via `ShapeFix_Solid::SolidFromShell`. Volume comparison against STEP reference and `BRepExtrema_DistShapeShape` distance checks pass for all test models
End-to-end tests in `tests/integration/end_to_end.rs` (one test binary, `tests/integration/main.rs`) verify that:
- STEP output files are written successfully for cube, cylinder, and sphere models
- All STL/STEP pairs pass `--compare` validation at stage 4.1, verifying volume agreement and `BRepExtrema_DistShapeShape` distance between written output and reference STEP
- Missing output path produces the expected error


Bad test STL files can be regenerated with `python3 scripts/generate_bad_tests.py`.

## Funny

https://www.smbc-comics.com/comic/crank

## License

MIT License
