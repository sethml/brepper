# Brepper - STL to STEP Converter

Convert triangulated STL meshes from CAD exports into parametric STEP files with fitted analytic and freeform surfaces.

## Quick Start

```bash
# Install dependencies (macOS)
brew install cmake pcl opencascade eigen cli11

# Build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(sysctl -n hw.ncpu)

# Run
./build/brepper input.stl -o output.step
```

## Usage

```
brepper - Convert STL mesh to STEP with fitted surfaces

USAGE:
    brepper [OPTIONS] <input.stl> -o <output.step>

REQUIRED:
    <input.stl>              Input STL file (binary or ASCII)
    -o, --output <file>      Output STEP file

GENERAL OPTIONS:
    -v, --verbose            Enable verbose output
    -q, --quiet              Suppress non-error output
    --debug                  Enable debug output and intermediate files
    --threads <N>            Number of threads (default: auto)
```

See [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) for full documentation.

## Dependencies

- PCL (Point Cloud Library) ≥1.12
- OpenCASCADE (OCCT) ≥7.6  
- Eigen ≥3.4
- CLI11 (header-only)
- OpenMP (required for parallel processing)

On macOS, install OpenMP via: `brew install libomp`

## License

MIT License