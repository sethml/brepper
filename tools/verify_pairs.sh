#!/bin/bash
# The pair-verification gate (Phase 0.5c): does a mesh actually lie on the
# reference surfaces it claims as ground truth?
#
# A STEP+STL pair is only ground truth if the mesh was tessellated from that
# exact geometry. A pair can fail for many reasons — a different revision, a
# decimated or repaired mesh, a healed or reparameterized import, a units or
# placement error — and every one of them silently poisons a --compare result.
# Run this on any newly acquired corpus before trusting it.
#
# Usage:
#   tools/verify_pairs.sh <dir> [--jobs N]
#
# <dir> must contain a sources.tsv mapping each STL stem to its reference STEP,
# the same manifest tools/sweep.sh reads.
#
# Two outcomes per file:
#   VERIFIED   the mesh's VERTICES lie on the reference surfaces
#   QUARANTINE they do not. Keep for crash-testing, never as truth.
#
# The gate is defined on vertices, per the plan: "the pair counts as ground truth
# only if mesh vertices lie on the STEP surfaces within a few quanta". Stage 1
# also gates facet CENTROID distance against --surface-tolerance, which measures
# something else entirely — the tessellation's sag — and a legitimately coarse
# mesh fails it while being perfect ground truth. Measured: four Fusion exports
# requested at 0.5 mm chordal deviation have centroid sags of 0.51-0.83 mm, over
# the 0.4 mm default, with vertex distances of 1.8e-7 to 5.4e-6 mm. Those are
# impeccable pairs from a coarse tessellation, not bad provenance.
#
# So the centroid gate is disabled here and the sag REPORTED instead. Note what
# it implies for Phase 3: brepper's default --surface-tolerance rejects any mesh
# exported at a chordal deviation above ~0.4 mm, which is an ordinary setting.
#
# This script briefly carried a third state, CONTAINED, for a pass reporting a
# distance of exactly zero — the signature of a mesh lying wholly inside its
# reference, which is how a corpus written 1000x too small once scored 15/15.
# That was a workaround for a bug in stage 1, which measured against the
# reference SOLID rather than its faces and so answered "is this point inside?".
# Stage 1 now measures against a face compound and reports the true distance
# (12.4875 mm for that corpus), so the heuristic is redundant and is gone. See
# src/compare.rs.
#
# The stage-1 summary's float32 position bound is reported alongside, because it
# reveals coordinate scale directly: the metres-scale sphere read 5.7e-10 mm
# where a correct mesh of the same model reads 7.2e-7 mm.

set -uo pipefail
cd "$(dirname "$0")/.."

BREPPER=./target/release/brepper
JOBS=8

DIR=""
while [ $# -gt 0 ]; do
    case "$1" in
        --jobs) JOBS="$2"; shift 2 ;;
        -h|--help) sed -n '2,44p' "$0"; exit 0 ;;
        *) DIR="${1%/}"; shift ;;
    esac
done

[ -n "$DIR" ] || { echo "usage: tools/verify_pairs.sh <dir> [--jobs N]" >&2; exit 2 ;}
[ -f "$DIR/sources.tsv" ] || { echo "no $DIR/sources.tsv manifest" >&2; exit 2; }
[ -x "$BREPPER" ] || { echo "build first: cargo build --release" >&2; exit 1; }

verify_one() {
    local stl="$1" dir
    dir=$(dirname "$stl")
    local stem ref out rc dist p95 bound
    stem=$(basename "$stl" .stl)
    ref=$(awk -F'\t' -v m="$stem" '$1 == m { print $2; exit }' "$dir/sources.tsv")
    if [ -z "$ref" ]; then
        printf 'NO-REFERENCE\t%s\t-\t-\t-\t-\n' "$stem"
        return
    fi
    # --surface-tolerance 1e6 disables the centroid gate; see the header. The
    # vertex tolerance stays at its default, which is what this gate is about.
    out=$("$BREPPER" "$stl" --stage=1 --compare "$ref" --surface-tolerance 1e6 2>&1)
    rc=$?
    dist=$(printf '%s' "$out" | grep -oE 'vertex dist p50 [0-9.e+-]+ / p95 [0-9.e+-]+ / max [0-9.e+-]+' \
        | head -1 | awk '{print $10}')
    p95=$(printf '%s' "$out" | grep -oE 'vertex dist p50 [0-9.e+-]+ / p95 [0-9.e+-]+' \
        | head -1 | awk '{print $7}')
    sag=$(printf '%s' "$out" | grep -oE 'centroid max dist [0-9.e+-]+' | head -1 | awk '{print $4}')
    bound=$(printf '%s' "$out" | grep -oE 'position bound [0-9.e+-]+' | head -1 | awk '{print $3}')
    if [ $rc -ne 0 ]; then
        local why
        why=$(printf '%s' "$out" | grep -oE 'comparison failed.*' | head -1 | cut -c1-90)
        printf 'QUARANTINE\t%s\t%s\t%s\t%s\t%s\t%s\n' "$stem" "${dist:--}" "${p95:--}" "${sag:--}" "${bound:--}" "$why"
    else
        printf 'VERIFIED\t%s\t%s\t%s\t%s\t%s\n' "$stem" "${dist:--}" "${p95:--}" "${sag:--}" "${bound:--}"
    fi
}
export -f verify_one
export BREPPER

results=$(find "$DIR" -name '*.stl' | sort | xargs -P "$JOBS" -I{} bash -c 'verify_one "$@"' _ {})

printf '%-12s %-46s %-14s %-13s %-12s %-12s %s\n' \
    STATUS MODEL VERTEX_MAX_MM VERTEX_P95_MM SAG_MM F32_BOUND_MM NOTE
printf '%s\n' "$results" | sort | awk -F'\t' \
    '{ printf "%-12s %-46s %-14s %-13s %-12s %-12s %s\n", \
        $1, $2, $3, $4, $5, $6, ($7 == "" ? "" : $7) }'

echo
total=$(printf '%s\n' "$results" | grep -c .)
for state in VERIFIED QUARANTINE NO-REFERENCE; do
    n=$(printf '%s\n' "$results" | grep -c "^$state" || true)
    [ "$n" -gt 0 ] && printf '%s: %d of %d\n' "$state" "$n" "$total"
done
ok=$(printf '%s\n' "$results" | grep -c '^VERIFIED' || true)
printf 'verification rate: %d/%d\n' "$ok" "$total"
[ "$ok" = "$total" ]
