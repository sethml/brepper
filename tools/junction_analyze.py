#!/usr/bin/env python3
"""Aggregate family-J rows (the junction race) from tools/invariants.sh.

    ./tools/invariants.sh --out tmp/junction.tsv --jobs 8 --junction
    python3 tools/junction_analyze.py tmp/junction.tsv

Question: when a patch grown on surface A crosses a tangent junction, does the
kinematic rank lose a tangent field before A's own refit stops validating?

A run is a RACE only when both signals had established A before the crossing
(the true-type fit validated the last pre-crossing patch, and the rank read A's
null dimension there). Otherwise one signal never had anything to lose.
"""

import sys
from collections import defaultdict

COLS = ("tag model ref kind start dist k_cross foreign_kind pre_val pre_rank "
        "k_val foreign_at_val k_rank foreign_at_rank").split()
TYPES = ["plane", "sphere", "cylinder", "cone", "torus"]


def num(s):
    return None if s == "-" else int(s)


def pct(a, b):
    return f"{100.0 * a / b:5.1f}%" if b else "    -"


def median(v):
    if not v:
        return "-"
    v = sorted(v)
    return str(v[len(v) // 2])


def main(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("JUNCTION\t"):
                f = line.rstrip("\n").split("\t")
                if len(f) >= len(COLS):
                    rows.append(dict(zip(COLS, f)))
    if not rows:
        sys.exit(f"no JUNCTION rows in {path}")

    faces = {(r["model"], r["ref"]) for r in rows}
    print(f"family J: {len(rows)} runs from {len(faces)} reference faces, "
          f"{len({r['model'] for r in rows})} models\n")

    print("Eligibility: did each signal establish A before the crossing?")
    print(f"    {'type':9} {'runs':>6} {'crossed':>8} {'val ok':>8} {'rank ok':>8} {'both (race)':>12}")
    for t in TYPES:
        sel = [r for r in rows if r["kind"] == t]
        if not sel:
            continue
        crossed = [r for r in sel if r["k_cross"] != "-"]
        v = sum(r["pre_val"] == "1" for r in crossed)
        k = sum(r["pre_rank"] == "1" for r in crossed)
        b = sum(r["pre_val"] == "1" and r["pre_rank"] == "1" for r in crossed)
        print(f"    {t:9} {len(sel):6d} {len(crossed):8d} {pct(v, len(crossed)):>8} "
              f"{pct(k, len(crossed)):>8} {pct(b, len(crossed)):>12}")
    print()

    print("The race, among runs where both signals had established A:")
    print(f"    {'type':9} {'races':>6} {'rank 1st':>9} {'tie':>7} {'val 1st':>8} "
          f"{'neither':>8}   foreign facets when it fires (p50): rank / val")
    tot = defaultdict(int)
    for t in TYPES:
        sel = [r for r in rows if r["kind"] == t and r["k_cross"] != "-"
               and r["pre_val"] == "1" and r["pre_rank"] == "1"]
        if not sel:
            continue
        rank_first = tie = val_first = neither = 0
        fr, fv = [], []
        for r in sel:
            kr, kv = num(r["k_rank"]), num(r["k_val"])
            if kr is not None:
                fr.append(int(r["foreign_at_rank"]))
            if kv is not None:
                fv.append(int(r["foreign_at_val"]))
            if kr is None and kv is None:
                neither += 1
            elif kv is None or (kr is not None and kr < kv):
                rank_first += 1
            elif kr is None or kv < kr:
                val_first += 1
            else:
                tie += 1
        n = len(sel)
        for key, val in (("n", n), ("rank", rank_first), ("tie", tie),
                         ("val", val_first), ("neither", neither)):
            tot[key] += val
        print(f"    {t:9} {n:6d} {pct(rank_first, n):>9} {pct(tie, n):>7} "
              f"{pct(val_first, n):>8} {pct(neither, n):>8}   {median(fr):>6} / {median(fv)}")
    n = tot["n"]
    if n:
        print(f"    {'all':9} {n:6d} {pct(tot['rank'], n):>9} {pct(tot['tie'], n):>7} "
              f"{pct(tot['val'], n):>8} {pct(tot['neither'], n):>8}")
    print()

    print("Validation alone (any run where the true-type fit established A):")
    print(f"    {'type':9} {'runs':>6} {'fires':>7}   foreign facets when it fires (p50)")
    for t in TYPES:
        sel = [r for r in rows if r["kind"] == t and r["k_cross"] != "-" and r["pre_val"] == "1"]
        if not sel:
            continue
        fired = [int(r["foreign_at_val"]) for r in sel if r["k_val"] != "-"]
        print(f"    {t:9} {len(sel):6d} {pct(len(fired), len(sel)):>7}   {median(fired)}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tmp/junction.tsv")
