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
    -v, --verbose            Enable verbose output
    -q, --quiet              Suppress non-error output
    --debug                  Enable debug output and intermediate files
    --stage=<stage>          Stop after stage, e.g. 2.2 (default: 4.1). --stage=2 stops after all of stage2
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

### Deviation Settings

Onshape calls this "Chordal tolerance (mm)", and the preset values are:
- Coarse: 0.24 mm
- Medium: 0.12 mm
- Fine: 0.06 mm

Fusion 360 calls this "Surface Deviation", and the preset values seem to vary based on model size:
- Low: 1x
- Medium: 2.5x
- High: 8x

### Generating Test Models

Some tests use models generated via [CodeCAD](https://codecad.xyz/). The generated STL/STEP files are committed to the repo, so CodeCAD is not required to run tests. 

To add or modify these models:
1. Install CodeCAD
2. Edit files in `tests/ccad/parts/`
3. Run `./tests/ccad/generate_models.sh`


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
- All STL/STEP pairs pass `--compare` validation at stage 2.1 (fitted planes are close to reference STEP surfaces)

Bad test STL files can be regenerated with `python3 scripts/generate_bad_tests.py`.

## License

MIT License