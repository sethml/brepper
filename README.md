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
    brepper [OPTIONS] <input.stl> [-o <output.step>]

REQUIRED:
    <input.stl>              Input STL file (binary or ASCII)
    -o, --output <file>      Output STEP file (required unless --stage < 6)

GENERAL OPTIONS:
    -v, --verbose            Enable verbose output
    -q, --quiet              Suppress non-error output
    --debug                  Enable debug output and intermediate files
    --threads <N>            Number of threads (default: auto)
    --stage <1-6>            Stop after stage (default: 6=export)
                             1=load, 2=segment, 3=assign, 4=boundary, 5=brep, 6=export
    --dimensions             Print mesh bounding box dimensions
```

### Examples

```bash
# Full conversion
./build/brepper input.stl -o output.step

# Just load and sample mesh (stage 1)
./build/brepper input.stl --stage 1 -v

# Run through segmentation (stage 2)
./build/brepper input.stl --stage 2 -v
```

See [DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md) for full documentation.

## Testing

```bash
# Build with tests enabled
cmake -B build -DBUILD_TESTS=ON
cmake --build build

# Run tests in parallel (recommended - ~4x faster)
ctest --test-dir build -j8 --output-on-failure

# Run tests sequentially
ctest --test-dir build --output-on-failure
```

### Debug Builds with Sanitizers

For catching memory corruption, undefined behavior, and concurrency issues:

```bash
# Address Sanitizer + Undefined Behavior Sanitizer (recommended)
cmake -B build-sanitize -DBUILD_TESTS=ON -DENABLE_ASAN=ON -DENABLE_UBSAN=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-sanitize && ctest --test-dir build-sanitize -j4 --output-on-failure

# Thread Sanitizer (for data races in OpenMP code)
cmake -B build-tsan -DBUILD_TESTS=ON -DENABLE_TSAN=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-tsan && ctest --test-dir build-tsan -j4 --output-on-failure
```

Note: TSAN cannot be combined with ASAN. Use `-j4` with sanitizers due to increased memory usage.

## Dependencies

- PCL (Point Cloud Library) ≥1.12
- OpenCASCADE (OCCT) ≥7.6  
- Eigen ≥3.4
- CLI11 (header-only)
- OpenMP (required for parallel processing)

On macOS, install OpenMP via: `brew install libomp`

## License

MIT License