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
printf "%-40s %15s %15s %15s %15s %8s %8s\n" "Test Case" "Vtx Max Dist" "Vtx Avg Dist" "Ctr Max Dist" "Ctr Avg Dist" "Nodes" "Faces"
printf "%-40s %15s %15s %15s %15s %8s %8s\n" \
    "$(printf '%0.s-' {1..40})" \
    "$(printf '%0.s-' {1..15})" \
    "$(printf '%0.s-' {1..15})" \
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

        # Run tool, capture stderr for diagnostics and stdout for distances
        diag_file=$(mktemp)
        output=$($BINARY "$stl_file" "$step_file" 2>"$diag_file") || {
            printf "%-40s %15s\n" "$base" "ERROR"
            rm -f "$diag_file"
            continue
        }

        # Parse stdout (tab-separated: vtx_max_dist ctr_max_dist)
        vtx_max_dist=$(echo "$output" | cut -f1)
        ctr_max_dist=$(echo "$output" | cut -f2)

        # Parse diagnostic output
        nodes=$(grep "^STL:" "$diag_file" | sed 's/STL: \([0-9]*\) nodes.*/\1/')
        faces=$(grep "^STEP:" "$diag_file" | sed 's/STEP: \([0-9]*\) faces/\1/')
        vtx_avg_dist=$(grep "^Vertex" "$diag_file" | sed 's/.*avg: \([^ ]*\).*/\1/')
        ctr_avg_dist=$(grep "^Centroid" "$diag_file" | sed 's/.*avg: \([^ ]*\).*/\1/')
        rm -f "$diag_file"

        printf "%-40s %15s %15s %15s %15s %8s %8s\n" "$base" "$vtx_max_dist" "$vtx_avg_dist" "$ctr_max_dist" "$ctr_avg_dist" "$nodes" "$faces"
    fi
done
