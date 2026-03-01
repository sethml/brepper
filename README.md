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
    --compare=<step>         STEP file to compare to at each step
    --tolerance=<meters>     Fitting tolerance - default 1e-6 = 1 micron
    -v, --verbose            Enable verbose output
    -q, --quiet              Suppress non-error output
    --debug                  Enable debug output and intermediate files
    --stage=<stage>          Stop after stage, e.g. 2.2 (default: 4.1). --stage=2 stops after all of stage2
```

### Examples

```bash
# Full conversion
./build/brepper input.stl -o output.step

# Just load and validate mesh (stage 1.2)
./build/brepper input.stl --stage 1.2 -v

# Warn and return with error code if mesh differs from 
./build/brepper input.stl --stage 2 -v
```

See [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) for full documentation.

### Generating Test Models

Some tests use models generated via [CodeCAD](https://codecad.xyz/). The generated STL/STEP files are committed to the repo, so CodeCAD is not required to run tests. 

To add or modify these models:
1. Install CodeCAD
2. Edit files in `tests/ccad/parts/`
3. Run `./tests/ccad/generate_models.sh`

## License

MIT License