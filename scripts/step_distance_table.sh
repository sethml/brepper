#!/bin/bash
# Run step-distance on all STL/STEP pairs under tests/ and print a table.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BINARY="$PROJECT_DIR/target/release/step-distance"

# Build release binary
echo "Building step-distance (release)..." >&2
cargo build --release --bin step-distance --manifest-path "$PROJECT_DIR/Cargo.toml" 2>&1 | tail -1 >&2

# Print header
printf "%-45s %15s %15s %8s %8s\n" "Test Case" "Max Distance" "Avg Distance" "Nodes" "Faces"
printf "%-45s %15s %15s %8s %8s\n" \
    "$(printf '%0.s-' {1..45})" \
    "$(printf '%0.s-' {1..15})" \
    "$(printf '%0.s-' {1..15})" \
    "$(printf '%0.s-' {1..8})" \
    "$(printf '%0.s-' {1..8})"

# Find all STL files and look for matching STEP files
find "$PROJECT_DIR/tests" -name "*.stl" -not -name "*.ascii.stl" | sort | while read -r stl_file; do
    step_file="${stl_file%.stl}.step"
    if [ -f "$step_file" ]; then
        rel_path="${stl_file#"$PROJECT_DIR"/tests/}"
        base="${rel_path%.stl}"

        # Run tool, capture stderr for diagnostics and stdout for max distance
        diag_file=$(mktemp)
        max_dist=$("$BINARY" "$stl_file" "$step_file" 2>"$diag_file") || {
            printf "%-45s %15s\n" "$base" "ERROR"
            rm -f "$diag_file"
            continue
        }

        # Parse diagnostic output
        nodes=$(grep "^STL:" "$diag_file" | sed 's/STL: \([0-9]*\) nodes.*/\1/')
        faces=$(grep "^STEP:" "$diag_file" | sed 's/STEP: \([0-9]*\) faces/\1/')
        avg_dist=$(grep "^Average" "$diag_file" | sed 's/Average distance: //')
        rm -f "$diag_file"

        printf "%-45s %15s %15s %8s %8s\n" "$base" "$max_dist" "$avg_dist" "$nodes" "$faces"
    fi
done
