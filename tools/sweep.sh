#!/bin/bash
# The sweep harness (Phase 0.3): every corpus model through the full pipeline
# with --compare, one machine-readable row per model, with baseline save/diff.
#
# This is the instrument the Robustness Rework is measured with. It replaces
# eyeballing tools/probe.sh output on 14 hand-picked models. See TODO.md
# ("Measurement protocol") and DEVELOPMENT_PLAN.md ("The Robustness Rework").
#
# Usage:
#   tools/sweep.sh [--out FILE] [--jobs N] [--filter REGEX] [--timeout SECS]
#   tools/sweep.sh --baseline [FILE]   run, then save the result as the baseline
#   tools/sweep.sh --diff [FILE]       run, then diff against the baseline
#   tools/sweep.sh --diff-only A B     diff two existing result files
#   tools/sweep.sh --list              print the model/reference pairing and exit
#
# Defaults: --out tmp/sweep.tsv, --jobs 8, --timeout 300, baseline
# tmp/sweep-baseline.tsv (tmp/ is gitignored, so baselines stay local).
#
# Output is TSV with a '#'-comment header. Columns:
#   model    path under tests/, without the .stl extension
#   status   ok | fail | timeout
#   ref      step (compared against a reference) | none (crash-test only)
#   n        total surfaces in the final selection
#   plane cyl cone sphere torus bspline    per-type breakdown
#   edge_mm  stage 3.1 max edge-boundary distance to the reference
#   vol_rel  stage 4.1 output-volume relative error
#   secs     wall time
#   note     first error line, or empty
#
# RUNTIME CAVEAT: with --jobs > 1 the secs column is inflated by contention.
# Use --jobs 1 for any runtime measurement (Phase 2). The header records the
# job count so a diff cannot silently compare the two.

set -uo pipefail
cd "$(dirname "$0")/.."

BREPPER=./target/release/brepper
OUT=tmp/sweep.tsv
JOBS=8
FILTER=""
# 300 s is above every corpus model's runtime (worst: bspline_hill ~285 s) and
# below 3DBenchy's, which does not finish — the one row that is meant to read
# `timeout` until Phase 6.2 addresses runtime on real-world meshes.
TIMEOUT=300
BASELINE=tmp/sweep-baseline.tsv
MODE=run

# ---------------------------------------------------------------------------
# Model/reference pairing
# ---------------------------------------------------------------------------
#
# Rules, in order. Anything that falls through runs without --compare and is
# reported as ref=none, which is correct for crash-test-only inputs but is
# also how a new model with a perfectly good reference could go unmeasured —
# so `--list` prints the pairing for review.

find_step() {
    local stl="$1"
    local dir base stripped
    dir=$(dirname "$stl")
    base=$(basename "$stl" .stl)

    # tests/remesh: <model>_offset_occt_<density> pairs with the translated
    # reference written next to it; <model>_occt_<density> with the original
    # reference in another directory, which tests/remesh/sources.tsv records
    # (the generator writes it — do not guess by searching for a matching stem,
    # plate_with_hole_100x50.step exists under two directories).
    if [[ "$base" == *_offset_occt_* ]]; then
        local off="${base%_occt_*}"
        [ -f "$dir/$off.step" ] && { echo "$dir/$off.step"; return; }
    fi
    if [[ "$base" == *_occt_* && -f tests/remesh/sources.tsv ]]; then
        local model="${base%_occt_*}" hit
        hit=$(awk -F'\t' -v m="$model" '$1 == m { print $2; exit }' tests/remesh/sources.tsv)
        [ -n "$hit" ] && { echo "$hit"; return; }
    fi

    base="${base%.ascii}"
    [ -f "$dir/$base.step" ] && { echo "$dir/$base.step"; return; }

    stripped=$(echo "$base" | sed -E 's/_(coarse|medium|fine|low|high)$//')
    [ -f "$dir/$stripped.step" ] && { echo "$dir/$stripped.step"; return; }

    echo ""
}

# Per-model CLI flags. rounded_cube_10_r2_coarse is tessellated past the 17.5°
# default dihedral limit, so it needs the same flag its suite test uses.
extra_flags() {
    case "$1" in
        *rounded_cube_10_r2_coarse*) echo "--angular-tolerance 20" ;;
        *) echo "" ;;
    esac
}

# Every STL under tests/ except:
#   bad/         inputs that must FAIL; the suite pins their exact errors
#   *.ascii.stl  duplicate geometry, pinned by a reader test instead
model_list() {
    find tests -name '*.stl' \
        -not -path 'tests/bad/*' \
        -not -name '*.ascii.stl' \
    | sort \
    | { [ -n "$FILTER" ] && grep -E "$FILTER" || cat; }
}

# ---------------------------------------------------------------------------
# One model
# ---------------------------------------------------------------------------

run_one() {
    local stl="$1"
    local model="${stl#tests/}"; model="${model%.stl}"
    local step flags log outstep start end secs rc
    step=$(find_step "$stl")
    flags=$(extra_flags "$stl")
    log=$(mktemp); outstep="$log.step"

    start=$(date +%s.%N)
    if [ -n "$step" ]; then
        timeout "$TIMEOUT" "$BREPPER" "$stl" -o "$outstep" --compare "$step" -v $flags \
            >"$log" 2>&1
    else
        timeout "$TIMEOUT" "$BREPPER" "$stl" -o "$outstep" -v $flags >"$log" 2>&1
    fi
    rc=$?
    end=$(date +%s.%N)
    secs=$(awk "BEGIN{printf \"%.2f\", $end - $start}")

    local status=ok
    [ $rc -eq 0 ] || status=fail
    [ $rc -eq 124 ] && status=timeout

    # "Stage 2 types: 2 plane, 2 cyl, 0 cone, 0 sphere, 4 torus, 0 bspline (8 total)"
    local types
    types=$(grep -m1 '^Stage 2 types:' "$log" \
        | sed -E 's/^Stage 2 types: ([0-9]+) plane, ([0-9]+) cyl, ([0-9]+) cone, ([0-9]+) sphere, ([0-9]+) torus, ([0-9]+) bspline \(([0-9]+) total\)/\7\t\1\t\2\t\3\t\4\t\5\t\6/')
    [ -n "$types" ] || types=$(printf -- '-\t-\t-\t-\t-\t-\t-')

    local edge vol note
    edge=$(grep -m1 'Compare 3.1: edge boundary max dist' "$log" \
        | sed -E 's/.*max dist ([0-9.e+-]+) mm.*/\1/')
    [ -n "$edge" ] || edge=-
    vol=$(grep -m1 'Compare 4.1: output volume' "$log" \
        | sed -E 's/.*rel diff ([0-9.e+-]+)\).*/\1/')
    [ -n "$vol" ] || vol=-
    note=""
    if [ "$status" != ok ]; then
        note=$(grep -m1 -iE '^Error|panicked at' "$log" | tr '\t' ' ' | cut -c1-160)
        [ "$status" = timeout ] && note="timed out after ${TIMEOUT}s"
    fi

    # $types is itself seven tab-separated fields, so this is one row of the
    # 14 documented columns.
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$model" "$status" "$([ -n "$step" ] && echo step || echo none)" \
        "$types" "$edge" "$vol" "$secs" "$note"
    # Progress on stderr: rows are sorted at the end, so without this a
    # ten-minute sweep looks like a hang.
    printf '  %-8s %6ss  %s\n' "$status" "$secs" "$model" >&2
    rm -f "$log" "$outstep"
}

# ---------------------------------------------------------------------------
# Diff
# ---------------------------------------------------------------------------
#
# Counts and status are compared exactly; edge/vol only when they move by more
# than 1% and above a picometre floor. Runtime is compared ONLY when both runs
# used --jobs 1: under parallelism the contention noise is larger than any
# regression worth reporting (measured: the same model, same binary, 27.8 s vs
# 63.3 s across two jobs=8 runs).

diff_results() {
    local old="$1" new="$2"
    [ -f "$old" ] || { echo "no baseline at $old" >&2; return 2; }
    local oj nj
    oj=$(grep -m1 '^# jobs' "$old" | awk '{print $3}')
    nj=$(grep -m1 '^# jobs' "$new" | awk '{print $3}')
    awk -v OLD="$old" -v serial="$([ "$oj" = 1 ] && [ "$nj" = 1 ] && echo 1 || echo 0)" '
        function num(x) { return (x == "-" ? -1 : x + 0) }
        # A measurement moved if one side is missing or they differ by >1%.
        #
        # NOISE is 1e-9 mm: an exact fit reports an edge distance around 1e-15 mm
        # — cancellation noise in a machine-epsilon quantity — and comparing two
        # of those RELATIVELY flags every run. A picometre is nine orders below
        # the f32 quantum of any real mesh, so below it, zero and 1.8e-15 are
        # the same measurement.
        function moved(a, b,   x, y) {
            x = num(a); y = num(b)
            if ((x < 0) != (y < 0)) return 1
            if (x < 0) return 0
            if (x < 1e-9 && y < 1e-9) return 0
            return (y > x * 1.01 || y < x * 0.99)
        }
        BEGIN {
            FS = OFS = "\t"
            split("model status ref n plane cyl cone sphere torus bspline " \
                  "edge_mm vol_rel secs note", name, " ")
            while ((getline line < OLD) > 0) {
                if (line ~ /^#/) continue
                split(line, f, FS)
                base[f[1]] = line
            }
        }
        /^#/ { next }
        {
            key = $1
            if (!(key in base)) { printf "NEW      %s (%s, n=%s)\n", key, $2, $4; changed++; next }
            split(base[key], b, FS)
            seen[key] = 1
            msg = ""
            # Status and the integer counts are identities: any change matters.
            for (i = 2; i <= 10; i++)
                if (b[i] != $i) msg = msg sprintf(" %s %s->%s", name[i], b[i], $i)
            for (i = 11; i <= 12; i++)
                if (moved(b[i], $i)) msg = msg sprintf(" %s %s->%s", name[i], b[i], $i)
            # Runtime: only meaningful serially, and only large moves.
            if (serial && num(b[13]) > 0.5 && \
                (num($13) > num(b[13]) * 1.5 || num($13) < num(b[13]) / 1.5))
                msg = msg sprintf(" secs %s->%s", b[13], $13)
            if (msg != "") { printf "CHANGED  %s:%s\n", key, msg; changed++ }
        }
        END {
            for (key in base)
                if (!(key in seen)) {
                    split(base[key], b, FS)
                    printf "GONE     %s (was %s)\n", key, b[2]
                    changed++
                }
            if (changed == 0) print "no changes vs baseline"
            else printf "%d model(s) changed\n", changed
        }
    ' "$new"
    if [ "$oj" != 1 ] || [ "$nj" != 1 ]; then
        echo "(runtime not compared: baseline jobs=$oj, this run jobs=$nj; use --jobs 1 for runtime)"
    fi
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

while [ $# -gt 0 ]; do
    case "$1" in
        --out) OUT="$2"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --filter) FILTER="$2"; shift 2 ;;
        --timeout) TIMEOUT="$2"; shift 2 ;;
        --baseline) MODE=baseline; [ $# -ge 2 ] && [[ "$2" != --* ]] && { BASELINE="$2"; shift; }; shift ;;
        --diff) MODE=diff; [ $# -ge 2 ] && [[ "$2" != --* ]] && { BASELINE="$2"; shift; }; shift ;;
        --diff-only) diff_results "$2" "$3"; exit $? ;;
        --list) MODE=list; shift ;;
        --one) shift; run_one "$1"; exit 0 ;;   # internal: one model, for xargs
        -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [ "$MODE" = list ]; then
    while read -r stl; do
        printf '%-58s %s\n' "${stl#tests/}" "$(find_step "$stl" || echo '(no reference)')"
    done < <(model_list)
    exit 0
fi

[ -x "$BREPPER" ] || { echo "build first: cargo build --release" >&2; exit 1; }
mkdir -p "$(dirname "$OUT")"

export -f run_one find_step extra_flags
export BREPPER TIMEOUT

count=$(model_list | wc -l | tr -d ' ')
echo "sweeping $count models, $JOBS at a time" >&2
# Rows land in $OUT.partial as they finish (the final file is sorted, so it
# cannot be written incrementally); it is deleted on success.
model_list | xargs -P "$JOBS" -I{} bash "$0" --one {} > "$OUT.partial"
{
    echo "# brepper sweep $(date '+%Y-%m-%d %H:%M:%S')"
    echo "# commit $(git rev-parse --short HEAD 2>/dev/null || echo unknown)$(git diff --quiet 2>/dev/null || echo '-dirty')"
    echo "# jobs $JOBS   models $count   timeout ${TIMEOUT}s"
    printf '#model\tstatus\tref\tn\tplane\tcyl\tcone\tsphere\ttorus\tbspline\tedge_mm\tvol_rel\tsecs\tnote\n'
    sort "$OUT.partial"
} > "$OUT"
rm -f "$OUT.partial"

echo "wrote $OUT" >&2
awk -F'\t' '!/^#/{ n[$2]++ } END { for (s in n) printf "%s: %d\n", s, n[s] }' "$OUT" >&2

case "$MODE" in
    baseline) cp "$OUT" "$BASELINE"; echo "baseline saved to $BASELINE" >&2 ;;
    diff)     diff_results "$BASELINE" "$OUT" ;;
esac
