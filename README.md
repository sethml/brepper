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
    --vertex-tolerance=<value>  Fitting tolerance in STL units (default: 1e-5)
    --surface-tolerance=<value> Surface-to-face offset tolerance in STL units (default: 0.4)
    --angular-tolerance=<deg>   Max dihedral angle between adjacent triangles on the same surface (default: 17.5)
    -v, --verbose[=LEVEL]    Verbosity level: 1 (or bare -v) = summaries, -vv or --verbose=2 = per-face details, -vvv or --verbose=3 = full BFS trace
    -q, --quiet              Suppress non-error output
    --debug                  Enable debug output and intermediate files
    --stage=<stage>          Stop after stage, e.g. 2.2 (default: 4.1). --stage=2 stops after all of stage2
    --viz-stages=<stages>    Open interactive 3D visualization at specified stages (e.g. 2.1,2.2,2.3)
```

### Examples

```bash
# Full conversion
./target/debug/brepper input.stl -o output.step

# Just load and validate mesh (stage 1.2)
./target/debug/brepper input.stl --stage 1.2 -v

# Compare mesh against reference STEP at each stage
./target/debug/brepper input.stl --compare reference.step --stage 1 -v
```

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

### Generating Test Models

Some tests use models generated via [CodeCAD](https://codecad.xyz/). The generated STL/STEP files are committed to the repo, so CodeCAD is not required to run tests.

To add or modify these models:
1. Install CodeCAD
2. Edit files in `tests/ccad/parts/` — see `tests/ccad/README.md` for primitive conventions and placement rules
3. Run `./tests/ccad/generate_models.sh`

**Important:** `box(w,d,h)` places a corner at the origin, not centered. Use `center_xyz()` to center. Each part file includes analytical surface area and volume in its header comments.


## Viewer

`viewer` is an interactive 3D viewer for STL and STEP files.

```
USAGE:
    viewer [OPTIONS] <file>...

ARGUMENTS:
    <file>...   One or more STL or STEP files to view (detected by extension)

OPTIONS:
    --deflection <value>          Linear deflection for STEP tessellation (default: 0.1)
    --angular-deflection <value>  Angular deflection for STEP tessellation in radians (default: 0.5)
```

Multiple files can be loaded simultaneously; each is rendered in a distinct color.

### Examples

```bash
# View a single STEP file
./target/release/viewer part.step

# View a STEP and an STL side by side
./target/release/viewer reference.step output.stl

# Finer tessellation for curved surfaces
./target/release/viewer part.step --deflection 0.01 --angular-deflection 0.1
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

### Interactive Debugging Visualization

Use `--viz-stages` to open a 3D window and step through BFS hypothesis deduction interactively:

```bash
# Visualize planar and cylindrical hypothesis stages
brepper input.stl --viz-stages=2.1,2.2

# Visualize all stage 2 hypothesis types with STEP comparison overlay
brepper input.stl --compare reference.step --viz-stages=2.1,2.2,2.3
```

At each BFS step, the window shows the STL mesh with highlighted faces: seed faces in green, accepted hypothesis faces in blue, and the face being evaluated in yellow. Cylindrical and spherical hypotheses also display translucent geometry overlays. When `--compare` is provided, the STEP surface is shown at 35% opacity for reference.

| Action | Control |
|--------|:--------|
| Next BFS step | `Space` |
| Skip to next seed/hypothesis | `Shift+Space` |
| Quit visualization | `Q` |

Supported stages: `2.1` (planar), `2.2` (cylindrical), `2.3` (spherical), `2.6` (surface selection).

### Running Tests

```bash
# Run all tests (unit + integration)
cargo test

# Run only integration tests
cargo test --test stage1_integration
```

Integration tests in `tests/stage1_integration.rs` verify that:
- All STL files in `tests/manual/`, `tests/onshape/`, `tests/ccad/generated/`, and `tests/fusion/` pass stage 1 mesh validation and coplanar triangle fusion
- All STL/STEP pairs pass `--compare` surface distance checks (after fusion)
- Bad files in `tests/bad/` fail with expected error types (degenerate faces, non-manifold edges, inconsistent winding, compare failures)

Integration tests in `tests/stage2_integration.rs` verify that:
- Planar hypothesis deduction produces correct counts for known planar models (cube=6, wedge=6, chamfered cube=26)
- Cylindrical hypothesis deduction produces correct counts and properties (simple_cylinder=1 convex, block_with_hole=1 concave, pipe=2, two_holes=2 concave)
- Spherical hypothesis deduction produces correct counts and properties (simple_sphere=1 convex r=10, hemisphere=1 convex r=10, spherical_pocket=1 concave r=8, ball_on_cylinder=1 sphere+1 cylinder)
- Surface selection (stage 2.6) produces correct surface type counts for known models (block_with_hole=6 planar+1 cylindrical, pipe=2 planar+2 cylindrical, ball_on_cylinder=1 planar+1 cylindrical+1 spherical, hemisphere=1 planar+1 spherical, simple_sphere=1 spherical)
- All STL/STEP pairs pass `--compare` validation at stages 2.1, 2.2, 2.3, and 2.6
- Deduced cylinder parameters (axis direction, axis position, radius) match STEP cylinder parameters for all cylindrical test models
- Deduced sphere parameters (center, radius) match STEP sphere parameters for all spherical test models
- Angular tolerance enforcement prevents spurious cylindrical/spherical hypotheses on planar-only models (e.g., cube faces meet at 90°, exceeding the 17.5° default)
Integration tests in `tests/stage3_integration.rs` verify that:
- Adjacency graph topology (faces, edges, vertices) is correct for known models (cube=6F/12E/8V, cylinder=3F/2E/0V, sphere=1F/0E/0V, etc.)
- Euler's formula V-E+F=2 holds for genus-0 models
- Edge validity, adjacency symmetry, and edge-face consistency are maintained
- All STL/STEP pairs pass `--compare` validation at stage 3.1, verifying that reconstructed edge boundary vertices lie on STEP edges and BRep vertices coincide with STEP vertices
- Tangency detection (stage 3.2) correctly identifies no tangent edges for models composed of planar, cylindrical, and spherical surfaces meeting at angles > 2°
- Edge curve computation (stage 3.3) produces 3D curves for all edges via surface-surface intersection, verified across all test models (cube, cylinder, hemisphere, ball_on_cylinder, block_with_hole, chamfered_cube, spherical_pocket)
- Face creation (stage 3.4) successfully creates OCCT `TopoDS_Face` objects for all test models: planar faces built from wire loops, cylindrical/spherical faces from UV parameter bounds. All STL/STEP pairs pass `--compare` validation at stage 3.4
- Shell construction (stage 3.5) stitches faces via `BRepBuilderAPI_Sewing` and produces correctly oriented shells. All STL/STEP pairs pass `--compare` validation
- Solid construction (stage 3.6) creates solids from shells via `ShapeFix_Solid::SolidFromShell`. Volume comparison against STEP reference and `BRepExtrema_DistShapeShape` distance checks pass for all test models
Integration tests in `tests/stage4_integration.rs` verify that:
- STEP output files are written successfully for cube, cylinder, and sphere models
- All STL/STEP pairs pass `--compare` validation at stage 4.1, verifying volume agreement and `BRepExtrema_DistShapeShape` distance between written output and reference STEP
- Missing output path produces the expected error


Bad test STL files can be regenerated with `python3 scripts/generate_bad_tests.py`.

## License

MIT License