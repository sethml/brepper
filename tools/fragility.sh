#!/bin/bash
# The fragility meter (Phase 0.4): perturb one stage-2 constant at a time and
# report which models change. Turns "minimize fragility" into a number.
#
# For each constant, the line is rewritten as `= (original expression) * 2` and
# `... / 2` (f64 constants use * 2.0 and * 0.5), the binary is rebuilt, and a
# sweep of the fragility model set is diffed against a baseline of the same set.
# The source is always restored, including on interrupt.
#
# Builds are serial — only one perturbation can be in the source at a time — but
# the sweeps are not, and that is where the time goes. Each wave builds --wave
# perturbed binaries, stashes them, then runs their sweeps CONCURRENTLY. Waves
# rather than all-at-once because sixty release binaries are over a gigabyte.
#
# Two numbers come out:
#   CHANGED   models whose result moved at all — sensitivity
#   CLIFF     models that COLLAPSED: a passing model that now fails, or a
#             surface count that at least doubled. This is the cliff count the
#             Robustness Rework tracks (four unrelated knobs once each flipped
#             stepped_block medium from 83 surfaces to a per-facet explosion).
#
# Usage:
#   tools/fragility.sh --baseline                     build the fragility baseline
#   tools/fragility.sh [--only REGEX] [--wave W] [--jobs J] [--out FILE]
#
# --wave W concurrent sweeps, each with --jobs J models in parallel: W x J is the
# core count to aim for. Defaults 5 x 2 on a 10-core machine.
#
# Cost, measured 2026-07-29: 61 s per perturbation at 5x2 against 105 s when each
# perturbation swept alone at --jobs 8, so the 60-perturbation full matrix is
# ~39 min (12 waves x ~157 s plus 8 s per build) against ~110 min. Core
# utilisation went from 39% to 63%.
#
# Why 5x2 and not another split: a sweep's wall time is
# max(longest single model, serial cost / jobs) = max(101 s, 315 s / J), and
# throughput is W/that with W*J fixed at the core count. J=2 and J=1 tie at
# ~0.032 perturbations/s; J=4 drops to 0.020 because the wave then waits on one
# 101 s model. Phase 2.1's runtime work on that model is what lifts the ceiling.
#
# RUN THE FULL MATRIX AT PHASE BOUNDARIES, not per commit; --only makes a
# targeted run minutes. Validated against the serial implementation on
# MAX_PATCH_DEPTH, NORMAL_AXIS_TOL, MIN_AXIAL_SWEEP_DEG and MIN_HALF_ANGLE_DEG:
# all eight perturbations reproduced the same changed/cliff counts and the same
# collapsed-model lists.

set -uo pipefail
cd "$(dirname "$0")/.."

JOBS=2
WAVE=5
OUT=tmp/fragility.tsv
BASELINE=tmp/fragility-baseline.tsv
ONLY=""
MODE=run

# ---------------------------------------------------------------------------
# The fragility model set
# ---------------------------------------------------------------------------
#
# Not the whole corpus: a full sweep costs 2861 s of CPU, and 68 of them is a
# day. This set keeps every model under ~30 s plus the one expensive model that
# matters — stepped_block_toroidal_fillets_medium, where four separate knobs
# have each been observed to cause a total collapse. It spans fillet networks,
# tangent junctions, both tessellators, all six surface types, and the seeded
# random families whose parameters were never tuned against.
#
# Excluded: the four models over 140 s (3DBenchy, bspline_hill/saddle/terrain,
# stepped fine, OrcaCube, pipe_elbow fine, hill_with_hole). If a constant is
# suspected of affecting those specifically, run tools/sweep.sh directly.
FRAG_SET='tests/(onshape/(stepped_block_toroidal_fillets_medium|rounded_cube_10_r2_(coarse|medium)|part_rounded_cube_10_r2_coarse|pipe_elbow_10_medium|pill_fine|chamfered_cube_10_c1_medium|plate_with_hole_100x50_coarse|sphere_25_fine|nosecone_medium|cone_tangent_planes_medium)|ccad/generated/(filleted_pipe|filleted_cylinder|filleted_hole_block|cone_cylinder|ball_on_cylinder|spherical_pocket|stepped_cylinder|two_holes|block_with_conical_hole)|nurbs/(bspline_wave|elliptical_arch|elliptical_cylinder)|occt/full_torus|remesh/(rounded_cube_10_r2_occt_medium|cone_tangent_planes_occt_medium|pill_occt_medium|sphere_25_occt_medium|stepped_block_toroidal_occt_medium|pipe_elbow_10_occt_medium)|random/)'

# Constants that are not tunables. Perturbing them measures nothing and one of
# them overflows.
skip_reason() {
    case "$1" in
        CUT_OUTSIDE) echo "sentinel index (usize::MAX - 1), not a threshold" ;;
        SHARP_EDGE_COLOR) echo "viz RGBA colour, not a threshold" ;;
        *) echo "" ;;
    esac
}

constants() {
    grep -Hn '^const \|^pub(crate) const \|^pub(super) const ' src/stage2/*.rs \
    | sed -E 's/^([^:]+):([0-9]+):(const|pub\(crate\) const|pub\(super\) const) ([A-Z0-9_]+): *([^=]+)= *(.*);$/\1\t\2\t\4\t\5\t\6/' \
    | { [ -n "$ONLY" ] && grep -E "$ONLY" || cat; }
}

# ---------------------------------------------------------------------------
# Perturb / restore
# ---------------------------------------------------------------------------

# Restore from a byte copy taken at perturbation time rather than from git, so
# the file comes back exactly as it was even if it was legitimately dirty.
DIRTY_FILE=""
SAVED_FILE=""
restore() {
    if [ -n "$DIRTY_FILE" ] && [ -f "$SAVED_FILE" ]; then
        cp "$SAVED_FILE" "$DIRTY_FILE"
        rm -f "$SAVED_FILE"
        DIRTY_FILE=""
    fi
}
trap 'restore; echo "restored source and stopped" >&2; exit 130' INT TERM
trap restore EXIT

# Rewrite the constant's value as (original) scaled, so the form of the
# original expression does not matter — TANGENCY_TOL_RAD is
# `2.0 * std::f64::consts::PI / 180.0` and scales exactly the same way.
perturb() {
    local file="$1" line="$2" ty="$3" expr="$4" dir="$5" new
    case "$ty" in
        *usize*|*u32*|*i32*) [ "$dir" = up ] && new="($expr) * 2" || new="($expr) / 2" ;;
        *) [ "$dir" = up ] && new="($expr) * 2.0" || new="($expr) * 0.5" ;;
    esac
    SAVED_FILE=$(mktemp)
    cp "$file" "$SAVED_FILE"
    DIRTY_FILE="$file"
    # The value runs to the end of the line, so replace from the first '=' on.
    python3 - "$file" "$line" "$new" <<'PY'
import sys
path, line, new = sys.argv[1], int(sys.argv[2]), sys.argv[3]
lines = open(path).read().splitlines(keepends=True)
head, _, _ = lines[line - 1].partition("=")
lines[line - 1] = f"{head}= {new};\n"
open(path, "w").writelines(lines)
PY
}

# ---------------------------------------------------------------------------
# Build one perturbed binary and stash it
# ---------------------------------------------------------------------------
#
# Builds must be serial — only one perturbation can be in the source at a time —
# but the SWEEPS do not have to be, and that is where the time goes. Each
# perturbed binary is copied aside so its sweep can run later, concurrently with
# other perturbations' sweeps.

build_perturbed() {
    local file="$1" line="$2" ty="$3" expr="$4" dir="$5" bin="$6"
    perturb "$file" "$line" "$ty" "$expr" "$dir"
    local rc=0
    cargo build --release --bin brepper >/dev/null 2>&1 || rc=1
    [ $rc -eq 0 ] && cp target/release/brepper "$bin"
    restore
    return $rc
}

# ---------------------------------------------------------------------------
# Classify one perturbation's sweep against the baseline
# ---------------------------------------------------------------------------
#
# Only fields 2..10 (status and the per-type counts) are compared, never the
# runtime — which is what makes it safe to run these sweeps at a different
# parallelism from the baseline's. brepper is single-threaded and deterministic,
# so concurrency changes timing and nothing else.

classify() {
    local out="$1"
    awk -F'\t' -v OLD="$BASELINE" '
        BEGIN {
            while ((getline l < OLD) > 0) {
                if (l ~ /^#/) continue
                split(l, f, "\t"); base[f[1]] = l
            }
        }
        /^#/ { next }
        {
            if (!($1 in base)) next
            split(base[$1], b, "\t")
            moved = 0
            for (i = 2; i <= 10; i++) if (b[i] != $i) moved = 1
            if (moved) changed++
            # Collapse: a passing model now fails, or the surface count at
            # least doubled (the per-facet explosion signature).
            if (b[2] == "ok" && $2 != "ok") { collapse++; which = which " " $1 }
            else if (b[4] != "-" && $4 != "-" && $4 + 0 >= 2 * (b[4] + 0) && b[4] + 0 > 0) {
                collapse++; which = which " " $1
            }
        }
        END { printf "%d\t%d\t%s\n", changed + 0, collapse + 0, (which == "" ? "-" : substr(which, 2)) }
    ' "$out"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

while [ $# -gt 0 ]; do
    case "$1" in
        --baseline) MODE=baseline; shift ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --wave) WAVE="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --only) ONLY="$2"; shift 2 ;;
        -h|--help) sed -n '2,34p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

git diff --quiet -- src/stage2/ || {
    echo "src/stage2/ has uncommitted changes; the meter restores files with" >&2
    echo "git checkout and would discard them. Commit or stash first." >&2
    exit 1
}

mkdir -p tmp/frag

if [ "$MODE" = baseline ]; then
    cargo build --release --bin brepper >/dev/null 2>&1 || { echo "build failed" >&2; exit 1; }
    ./tools/sweep.sh --filter "$FRAG_SET" --jobs $((WAVE * JOBS)) --out "$BASELINE" 2>&1 | tail -4
    echo "fragility baseline: $(grep -vc '^#' "$BASELINE") models -> $BASELINE" >&2
    exit 0
fi

[ -f "$BASELINE" ] || { echo "run --baseline first" >&2; exit 1; }

{
    echo "# brepper fragility meter $(date '+%Y-%m-%d %H:%M:%S')"
    echo "# commit $(git rev-parse --short HEAD)   models $(grep -vc '^#' "$BASELINE")   jobs $JOBS"
    printf '#constant\tfile\tvalue\tdown_changed\tdown_cliff\tup_changed\tup_cliff\tcollapsed_models\n'
} > "$OUT"

# Work list (one line per perturbation) and metadata (one line per constant).
# Plain files and indexed arrays only: /bin/bash on macOS is 3.2, with no
# associative arrays and no mapfile.
WORK=$(mktemp)
METAF=$(mktemp)
while IFS=$'\t' read -r file line name ty expr; do
    reason=$(skip_reason "$name")
    if [ -n "$reason" ]; then
        printf '%s\tskip\t%s\t%s\t%s\n' "$name" "${file#src/stage2/}" "$expr" "$reason" >> "$METAF"
        continue
    fi
    printf '%s\tdo\t%s\t%s\t\n' "$name" "${file#src/stage2/}" "$expr" >> "$METAF"
    for dir in down up; do
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$dir" "$file" "$line" "$ty" "$expr" >> "$WORK"
    done
done < <(constants)

WL=()
while IFS= read -r l; do WL[${#WL[@]}]="$l"; done < "$WORK"
count=${#WL[@]}
echo "$count perturbations: builds serial, sweeps $WAVE at a time x $JOBS jobs" >&2

mkdir -p tmp/frag/bin
i=0
while [ $i -lt $count ]; do
    bins=(); names=(); j=$i
    while [ $j -lt $count ] && [ $((j - i)) -lt "$WAVE" ]; do
        IFS=$'\t' read -r pname pdir pfile pline pty pexpr <<< "${WL[$j]}"
        bin="tmp/frag/bin/${pname}_${pdir}"
        if build_perturbed "$pfile" "$pline" "$pty" "$pexpr" "$pdir" "$bin"; then
            bins[${#bins[@]}]="$bin"
            names[${#names[@]}]="${pname}_${pdir}"
        else
            echo "  BUILD FAILED $pname $pdir" >&2
            printf 'build-failed\t-\t-\n' > "tmp/frag/${pname}_${pdir}.class"
        fi
        j=$((j + 1))
    done
    k=0
    while [ $k -lt ${#bins[@]} ]; do
        ( BREPPER="${bins[$k]}" ./tools/sweep.sh --filter "$FRAG_SET" --jobs "$JOBS" \
              --out "tmp/frag/${names[$k]}.tsv" >/dev/null 2>&1
          classify "tmp/frag/${names[$k]}.tsv" > "tmp/frag/${names[$k]}.class" ) &
        k=$((k + 1))
    done
    wait
    [ ${#bins[@]} -gt 0 ] && rm -f "${bins[@]}"
    k=0
    while [ $k -lt ${#names[@]} ]; do
        IFS=$'\t' read -r c cl w < "tmp/frag/${names[$k]}.class"
        printf '  %-40s changed=%-3s cliff=%-3s %s\n' "${names[$k]}" "$c" "$cl" "$w" >&2
        k=$((k + 1))
    done
    i=$j
done
rm -f "$WORK"

# Rebuild the unperturbed binary so the tree is left usable.
cargo build --release --bin brepper >/dev/null 2>&1

# Assemble the table in constant order.
cliffs=0
total=0
while IFS=$'\t' read -r name kind file expr reason; do
    if [ "$kind" = skip ]; then
        printf '%s\t%s\t%s\tskipped\tskipped\tskipped\tskipped\t%s\n' \
            "$name" "$file" "$expr" "$reason" >> "$OUT"
        continue
    fi
    total=$((total + 1))
    row=""; collapsed_all=""
    for dir in down up; do
        IFS=$'\t' read -r c cl w < "tmp/frag/${name}_${dir}.class"
        row="$row	$c	$cl"
        if [ -n "$w" ] && [ "$w" != "-" ]; then
            for m in $w; do collapsed_all="$collapsed_all ${dir}:$m"; done
        fi
    done
    [ -n "$collapsed_all" ] && cliffs=$((cliffs + 1))
    printf '%s\t%s\t%s%s\t%s\n' "$name" "$file" "$expr" "$row" \
        "$([ -n "$collapsed_all" ] && echo "${collapsed_all# }" || echo -)" >> "$OUT"
done < "$METAF"
rm -f "$METAF"

echo "# cliff count: $cliffs of $total constants perturbed" >> "$OUT"
echo "wrote $OUT" >&2
echo "CLIFF COUNT: $cliffs of $total constants perturbed" >&2
