#!/bin/bash
# Invariant harness: run `--invariant-report` over every corpus model that has
# a reference STEP, and concatenate the per-reference-face TSV rows.
#
# Reuses tools/sweep.sh's model/reference pairing (`--list`) so the two
# instruments always measure the same corpus.
#
# Usage:
#   tools/invariants.sh [--out FILE] [--jobs N] [--filter REGEX] [--timeout SECS]
#                       [--seeds N] [--no-compete]
#
# Output TSV columns are documented in src/stage2/invariants.rs. Per-model
# stderr summaries go to <out>.log.

set -uo pipefail
cd "$(dirname "$0")/.."

BREPPER=${BREPPER:-./target/release/brepper}
OUT=tmp/invariants.tsv
JOBS=6
FILTER=""
TIMEOUT=900
SEEDS=""
NO_COMPETE=""

while [ $# -gt 0 ]; do
    case "$1" in
        --out) OUT="$2"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --filter) FILTER="$2"; shift 2 ;;
        --timeout) TIMEOUT="$2"; shift 2 ;;
        --seeds) SEEDS="$2"; shift 2 ;;
        --no-compete) NO_COMPETE=1; shift ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "$(dirname "$OUT")"
LOG="$OUT.log"
: > "$LOG"

command -v gtimeout >/dev/null && TO=gtimeout || TO=timeout

run_one() {
    local model="$1" step="$2"
    local stl="tests/$model"
    local name="${model}"
    local env_args=("BREPPER_INVARIANTS=1")
    [ -n "$SEEDS" ] && env_args+=("BREPPER_INV_SEEDS=$SEEDS")
    [ -n "$NO_COMPETE" ] && env_args+=("BREPPER_INV_NO_COMPETE=1")
    local out err
    err=$(mktemp)
    out=$(env "${env_args[@]}" $TO "$TIMEOUT" "$BREPPER" "$stl" --compare "$step" \
              --stage 2.2 --invariant-report -v 2>"$err" )
    local rc=$?
    {
        echo "=== $name (rc=$rc)"
        cat "$err"
    } >> "$LOG"
    rm -f "$err"
    # Rewrite the model column to the full corpus-relative stem so rows from
    # two directories with the same basename stay distinguishable.
    echo "$out" | awk -v m="${model%.stl}" -F'\t' 'BEGIN{OFS="\t"} /^(INV|SEEDRULE)\t/ {$2=m; print}'
}

# sweep.sh --list pads the model column with spaces, and some model paths
# contain single spaces ("OrcaCube v2.stl"), so the separator is 2+ spaces.
pairs=$(./tools/sweep.sh --list 2>/dev/null \
        | sed -E 's/[[:space:]][[:space:]]+/\t/' \
        | awk -F'\t' 'NF==2 && $2 != ""')
[ -n "$FILTER" ] && pairs=$(echo "$pairs" | grep -E "$FILTER")

n=$(echo "$pairs" | grep -c .)
echo "invariants: $n model(s), jobs=$JOBS, timeout=${TIMEOUT}s" >&2

PAIRFILE=$(mktemp)
echo "$pairs" > "$PAIRFILE"

# bash 3.2 (macOS) has no `wait -n`, so hold the job count with a poll.
while IFS=$'\t' read -r model step; do
    [ -z "$model" ] && continue
    [ -z "$step" ] && continue
    while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "$JOBS" ]; do sleep 0.3; done
    run_one "$model" "$step" >> "$OUT" &
done < "$PAIRFILE"
wait
rm -f "$PAIRFILE"

echo "wrote $OUT ($(grep -c '^INV' "$OUT") reference-face rows); stderr in $LOG" >&2
