#!/bin/zsh

# Run step-measure-tolerance on all STL/STEP pairs under tests/ and print a table
# with per-group (suffix) summary rows.
# Distances displayed in mm. OCCT converts STEP units to mm internally.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BINARY="$PROJECT_DIR/target/release/step-measure-tolerance"

# Build release binary
echo "Building step-measure-tolerance (release)..." >&2
cargo build --release --bin step-measure-tolerance --manifest-path "$PROJECT_DIR/Cargo.toml" 2>&1 | tail -1 >&2

ROW_FMT="%-40s %12s %12s %12s %12s %10s %10s %10s %8s %8s %10s %10s %10s %10s\n"

print_header() {
    printf "$ROW_FMT" "Test Case" "Vtx Max mm" "Vtx Avg mm" "Ctr Max mm" "Ctr Avg mm" "Max Dim" "Ang Max \u00b0" "Ang Avg \u00b0" "Nodes" "Faces" "STEP Area" "STEP Vol" "STL Area" "STL Vol"
    print_sep
}

print_sep() {
    printf "$ROW_FMT" \
        "$(printf '%0.s-' {1..40})" \
        "$(printf '%0.s-' {1..12})" \
        "$(printf '%0.s-' {1..12})" \
        "$(printf '%0.s-' {1..12})" \
        "$(printf '%0.s-' {1..12})" \
        "$(printf '%0.s-' {1..10})" \
        "$(printf '%0.s-' {1..10})" \
        "$(printf '%0.s-' {1..10})" \
        "$(printf '%0.s-' {1..8})" \
        "$(printf '%0.s-' {1..8})" \
        "$(printf '%0.s-' {1..10})" \
        "$(printf '%0.s-' {1..10})" \
        "$(printf '%0.s-' {1..10})" \
        "$(printf '%0.s-' {1..10})"
}

# Format vertex distance (mm, 6 decimal places)
fmt_vtx() {
    awk "BEGIN { v = $1 + 0; if (v == 0) printf \"0\"; else printf \"%.6f\", v }"
}

# Format centroid distance (mm, 3 decimal places)
fmt_ctr() {
    awk "BEGIN { v = $1 + 0; if (v == 0) printf \"0\"; else printf \"%.3f\", v }"
}
# Format angular deviation (degrees, 4 decimal places)
fmt_ang() {
    awk "BEGIN { v = $1 + 0; if (v == 0) printf \"0\"; else printf \"%.4f\", v }"
}


# Collect all STL files into an array (avoids subshell from pipe)
stl_files=(${(f)"$(find "$PROJECT_DIR/tests" -name '*.stl' -not -name '*.ascii.stl' | sort)"})

# Track per-group maximums
typeset -A group_vtx_max group_vtx_avg group_ctr_max group_ctr_avg group_ang_max group_ang_avg
prev_dir=""
group_order=()

flush_groups() {
    if [[ ${#group_order[@]} -gt 0 ]]; then
        for g in "${group_order[@]}"; do
            local vm=$(fmt_vtx "${group_vtx_max[$g]}")
            local va=$(fmt_vtx "${group_vtx_avg[$g]}")
            local cm=$(fmt_ctr "${group_ctr_max[$g]}")
            local am=$(fmt_ang "${group_ang_max[$g]}")
            local aa=$(fmt_ang "${group_ang_avg[$g]}")
            printf "$ROW_FMT" "  ** MAX ($g) **" "$vm" "$va" "$cm" "$ca" "" "$am" "$aa" "" "" "" "" "" ""
        done
    fi
    group_vtx_max=()
    group_vtx_avg=()
    group_ctr_max=()
    group_ctr_avg=()
    group_ang_max=()
    group_ang_avg=()
    group_order=()
}

update_group_max() {
    local group="$1" vm="$2" va="$3" cm="$4" ca="$5" am="$6" aa="$7"
    if [[ -z "${group_vtx_max[$group]:-}" ]]; then
        group_vtx_max[$group]="$vm"
        group_vtx_avg[$group]="$va"
        group_ctr_max[$group]="$cm"
        group_ctr_avg[$group]="$ca"
        group_ang_max[$group]="$am"
        group_ang_avg[$group]="$aa"
        group_order+=("$group")
    else
        group_vtx_max[$group]=$(awk "BEGIN { a=${group_vtx_max[$group]}; b=$vm; print (a>b) ? a : b }")
        group_vtx_avg[$group]=$(awk "BEGIN { a=${group_vtx_avg[$group]}; b=$va; print (a>b) ? a : b }")
        group_ctr_max[$group]=$(awk "BEGIN { a=${group_ctr_max[$group]}; b=$cm; print (a>b) ? a : b }")
        group_ctr_avg[$group]=$(awk "BEGIN { a=${group_ctr_avg[$group]}; b=$ca; print (a>b) ? a : b }")
        group_ang_max[$group]=$(awk "BEGIN { a=${group_ang_max[$group]}; b=$am; print (a>b) ? a : b }")
        group_ang_avg[$group]=$(awk "BEGIN { a=${group_ang_avg[$group]}; b=$aa; print (a>b) ? a : b }")
    fi
}

print_header

for stl_file in "${stl_files[@]}"; do
    step_file="${stl_file%.stl}.step"
    [ -f "$step_file" ] || continue

    rel_path="${stl_file#${PROJECT_DIR}/tests/}"
    base="${rel_path%.stl}"
    dir=$(dirname "$base")

    # On directory change, print group summaries
    if [ "$dir" != "$prev_dir" ] && [ -n "$prev_dir" ]; then
        flush_groups
        echo
    fi
    prev_dir="$dir"

    # Extract group suffix (after last _)
    filename=$(basename "$base")
    group="${filename##*_}"

    # Run tool, capture both stdout and stderr
    diag_file=$(mktemp)
    output=$("$BINARY" "$stl_file" "$step_file" 2>"$diag_file") || {
        printf "$ROW_FMT" "$base" "ERROR" "" "" "" "" "" "" "" "" "" "" "" ""

        rm -f "$diag_file"
        continue
    }

    # Parse tab-separated stdout: vtx_max vtx_avg ctr_max ctr_avg max_dim ang_max ang_avg step_area step_vol stl_area stl_vol
    vtx_max_raw=$(echo "$output" | cut -f1)
    vtx_avg_raw=$(echo "$output" | cut -f2)
    ctr_max_raw=$(echo "$output" | cut -f3)
    ctr_avg_raw=$(echo "$output" | cut -f4)
    max_dim=$(echo "$output" | cut -f5)
    ang_max_raw=$(echo "$output" | cut -f6)
    ang_avg_raw=$(echo "$output" | cut -f7)
    step_area_raw=$(echo "$output" | cut -f8)
    step_vol_raw=$(echo "$output" | cut -f9)
    stl_area_raw=$(echo "$output" | cut -f10)
    stl_vol_raw=$(echo "$output" | cut -f11)

    # Parse stderr for node/face counts
    nodes=$(grep "^STL:" "$diag_file" | sed 's/STL: \([0-9]*\) nodes.*/\1/')
    faces=$(grep "^STEP:" "$diag_file" | sed 's/STEP: \([0-9]*\) faces/\1/')
    rm -f "$diag_file"

    # Format distances for display (mm)
    vm=$(fmt_vtx "$vtx_max_raw")
    va=$(fmt_vtx "$vtx_avg_raw")
    cm=$(fmt_ctr "$ctr_max_raw")
    ca=$(fmt_ctr "$ctr_avg_raw")
    am=$(fmt_ang "$ang_max_raw")
    aa=$(fmt_ang "$ang_avg_raw")
    sa=$(fmt_ctr "$step_area_raw")
    sv=$(fmt_ctr "$step_vol_raw")
    la=$(fmt_ctr "$stl_area_raw")
    lv=$(fmt_ctr "$stl_vol_raw")

    printf "$ROW_FMT" "$base" "$vm" "$va" "$cm" "$ca" "$max_dim" "$am" "$aa" "$nodes" "$faces" "$sa" "$sv" "$la" "$lv"

    update_group_max "$group" "$vtx_max_raw" "$vtx_avg_raw" "$ctr_max_raw" "$ctr_avg_raw" "$ang_max_raw" "$ang_avg_raw"

done

# Print final group summaries
flush_groups
