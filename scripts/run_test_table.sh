#!/bin/bash
# Run all STL test files through brepper and collect results
# Usage: ./scripts/run_test_table.sh <brepper_binary>

set -euo pipefail
BREPPER="${1:?Usage: $0 <brepper_binary>}"
cd "$(dirname "$0")/.."

# Map STL file to its corresponding STEP file
find_step() {
    local stl="$1"
    local dir=$(dirname "$stl")
    local base=$(basename "$stl" .stl)
    
    # Strip .ascii suffix if present
    base="${base%.ascii}"
    
    # Try exact match first
    if [ -f "$dir/$base.step" ]; then
        echo "$dir/$base.step"
        return
    fi
    
    # Strip resolution suffixes (_coarse, _medium, _fine, _low, _high)
    local stripped=$(echo "$base" | sed -E 's/_(coarse|medium|fine|low|high)$//')
    if [ -f "$dir/$stripped.step" ]; then
        echo "$dir/$stripped.step"
        return
    fi
    
    echo ""
}

# Special flags for specific tests
get_extra_flags() {
    local stl="$1"
    case "$stl" in
        *rounded_cube_10_r2_coarse*)
            echo "--angular-tolerance 20"
            ;;
        *)
            echo ""
            ;;
    esac
}

echo "STL_FILE|STEP_FILE|STATUS|TIME_S"

for stl in $(find tests -name '*.stl' | sort); do
    step=$(find_step "$stl")
    extra_flags=$(get_extra_flags "$stl")
    
    if [ -z "$step" ]; then
        compare_flag=""
    else
        compare_flag="--compare $step"
    fi
    
    tmpout=$(mktemp)
    start_time=$(python3 -c 'import time; print(time.monotonic())')
    
    if timeout 120 "$BREPPER" "$stl" -o /dev/null $compare_flag $extra_flags > "$tmpout" 2>&1; then
        status="ok"
    else
        exitcode=$?
        if [ $exitcode -eq 124 ]; then
            status="timeout"
        else
            # Extract failure stage from output
            last_line=$(tail -5 "$tmpout" | grep -iE "error|panic|fail|abort|thread" | tail -1 || true)
            if [ -n "$last_line" ]; then
                # Truncate long error messages
                status="FAIL: $(echo "$last_line" | head -c 80)"
            else
                status="FAIL (exit $exitcode)"
            fi
        fi
    fi
    
    end_time=$(python3 -c 'import time; print(time.monotonic())')
    elapsed=$(python3 -c "print(f'{$end_time - $start_time:.2f}')")
    
    echo "$stl|$step|$status|${elapsed}"
    rm -f "$tmpout"
done
