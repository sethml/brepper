#!/usr/bin/env python3
"""Aggregate family-K rows (the kinematic rank ladder) from tools/invariants.sh.

    ./tools/invariants.sh --out tmp/kin.tsv --jobs 8 --kinematic
    python3 tools/kin_analyze.py tmp/kin.tsv

Answers the four questions the probe is built around, in order:

  1. Is the signal there?      — the reference's own generators' residual,
                                 per true type, at each evaluation point.
  2. Does the read converge?   — first support size at which the read is right,
                                 and whether it stays right.
  3. Is it believable?         — the ratio gap, split by right/wrong reads.
  4. Does straddling govern?   — gap and correctness against straddle fraction.
"""

import sys
from collections import defaultdict

COLS = (
    "tag model ref kind start support dim_c read_c ok_c gap_c "
    "dim_x read_x ok_x gap_x straddle gen_c gen_x l0 l1 l2 l3 l4 l5 l6"
).split()

TYPES = ["plane", "sphere", "cylinder", "cone", "torus"]


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def pct(a, b):
    return f"{100.0 * a / b:5.1f}%" if b else "    - "


def quantiles(v, qs=(0.5, 0.95)):
    if not v:
        return [float("nan")] * len(qs)
    v = sorted(v)
    return [v[min(len(v) - 1, int(round(q * (len(v) - 1))))] for q in qs]


def main(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("KINEMATIC\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < len(COLS):
                continue
            rows.append(dict(zip(COLS, f)))
    if not rows:
        sys.exit(f"no KINEMATIC rows in {path}")

    models = {r["model"] for r in rows}
    faces = {(r["model"], r["ref"]) for r in rows}
    print(f"family K: {len(rows)} readings, {len(faces)} reference faces, "
          f"{len(models)} models\n")

    # --- Q1: is the signal there? ------------------------------------------
    print("Q1  the reference's OWN generators against the measured rows")
    print("    (relative residual; the null space cannot be cleaner than this)")
    print(f"    {'type':9} {'n':>7}  {'centroid p50':>12} {'p95':>10}"
          f"  {'circumcenter p50':>16} {'p95':>10}   winner")
    for t in TYPES:
        gc = [num(r["gen_c"]) for r in rows if r["kind"] == t]
        gx = [num(r["gen_x"]) for r in rows if r["kind"] == t]
        gc = [v for v in gc if v == v]
        gx = [v for v in gx if v == v]
        if not gc:
            continue
        c50, c95 = quantiles(gc)
        x50, x95 = quantiles(gx)
        win = "circumcenter" if x50 < c50 * 0.5 else (
            "centroid" if c50 < x50 * 0.5 else "tie")
        print(f"    {t:9} {len(gc):7d}  {c50:12.2e} {c95:10.2e}"
              f"  {x50:16.2e} {x95:10.2e}   {win}")
    print()

    # --- Q2: does the read converge, and through what? ---------------------
    for evl, dimk, readk, okk, gapk in (
        ("centroid", "dim_c", "read_c", "ok_c", "gap_c"),
        ("circumcenter", "dim_x", "read_x", "ok_x", "gap_x"),
    ):
        print(f"Q2  correctness by support size — {evl}")
        sizes = sorted({int(r["support"]) for r in rows})
        hdr = "    " + f"{'type':9}" + "".join(f"{s:>7}" for s in sizes)
        print(hdr)
        for t in TYPES:
            cells = []
            for s in sizes:
                sel = [r for r in rows if r["kind"] == t and int(r["support"]) == s]
                ok = sum(1 for r in sel if r[okk] == "1")
                cells.append(f"{100.0*ok/len(sel):6.0f}%" if sel else "     -")
            print(f"    {t:9}" + "".join(cells))
        print()

    # What do the wrong reads say? A wrong read that is SIMPLER than the truth
    # is the designed behaviour; a wrong read that is a different type at the
    # same or greater complexity is the failure the ladder must not have.
    order = {"under": -1, "plane": 0, "sphere": 1, "cylinder": 2, "cone": 2,
             "revolution": 3, "extrusion": 3, "spiral": 3, "freeform": 4}
    truth_order = {"plane": 0, "sphere": 1, "cylinder": 2, "cone": 2, "torus": 3}
    print("Q2b what the wrong reads actually say (centroid / circumcenter)")
    print(f"    {'type':9} {'simpler':>9} {'under':>9} {'WRONG TYPE':>12}"
          f"  {'top wrong reads'}")
    for t in TYPES:
        for evl, readk, okk in (("cen", "read_c", "ok_c"),
                                ("cir", "read_x", "ok_x")):
            sel = [r for r in rows if r["kind"] == t and r[okk] != "1"]
            if not sel:
                continue
            simpler = under = wrong = 0
            hist = defaultdict(int)
            for r in sel:
                g = r[readk]
                hist[g] += 1
                if g == "under":
                    under += 1
                elif order.get(g, 9) < truth_order[t]:
                    simpler += 1
                else:
                    wrong += 1
            top = ", ".join(f"{k}:{v}" for k, v in
                            sorted(hist.items(), key=lambda kv: -kv[1])[:3])
            n = len(sel)
            print(f"    {t+'/'+evl:9} {pct(simpler,n):>9} {pct(under,n):>9}"
                  f" {pct(wrong,n):>12}  {top}")
    print()

    # --- Q3: how believable is a reading? ----------------------------------
    print("Q3  the ratio gap, split by whether the read was right")
    print(f"    {'type':9} {'eval':5} {'right p50':>11} {'wrong p50':>11}"
          f" {'wrong p95':>11}  separable?")
    for t in TYPES:
        for evl, okk, gapk in (("cen", "ok_c", "gap_c"), ("cir", "ok_x", "gap_x")):
            right = [num(r[gapk]) for r in rows
                     if r["kind"] == t and r[okk] == "1"]
            wrong = [num(r[gapk]) for r in rows
                     if r["kind"] == t and r[okk] != "1" and r[readk_of(evl)] != "under"]
            right = [v for v in right if v == v]
            wrong = [v for v in wrong if v == v]
            if not right and not wrong:
                continue
            r50 = quantiles(right)[0] if right else float("nan")
            w50, w95 = quantiles(wrong) if wrong else (float("nan"),) * 2
            sep = "yes" if (right and wrong and r50 > w95) else "NO"
            print(f"    {t:9} {evl:5} {r50:11.2e} {w50:11.2e} {w95:11.2e}  {sep}")
    print()

    # --- Q4: straddling ----------------------------------------------------
    print("Q4  correctness against straddle fraction (ruled types only)")
    print(f"    {'type':9} {'eval':5} {'straddle=0':>11} {'0<s<=0.25':>11}"
          f" {'0.25<s<=0.5':>12} {'s>0.5':>9}")
    bins = [(0.0, 0.0), (0.0, 0.25), (0.25, 0.5), (0.5, 1.01)]
    for t in ("cylinder", "cone", "torus"):
        for evl, okk in (("cen", "ok_c"), ("cir", "ok_x")):
            cells = []
            for i, (lo, hi) in enumerate(bins):
                sel = [r for r in rows if r["kind"] == t
                       and num(r["straddle"]) == num(r["straddle"])
                       and (num(r["straddle"]) == 0.0 if i == 0
                            else lo < num(r["straddle"]) <= hi)]
                ok = sum(1 for r in sel if r[okk] == "1")
                cells.append(f"{100.0*ok/len(sel):9.0f}% " if sel else "        - ")
            print(f"    {t:9} {evl:5}" + "".join(cells))
    print()

    # --- the headline: best achievable read per face -----------------------
    print("Headline  per reference face, does SOME support size read it right?")
    print(f"    {'type':9} {'faces':>7} {'centroid':>10} {'circumcenter':>14}"
          f" {'either':>9}   {'first size that works (p50)':>28}")
    for t in TYPES:
        byface = defaultdict(list)
        for r in rows:
            if r["kind"] == t:
                byface[(r["model"], r["ref"])].append(r)
        if not byface:
            continue
        okc = okx = oke = 0
        firsts = []
        for _, rs in byface.items():
            c = any(r["ok_c"] == "1" for r in rs)
            x = any(r["ok_x"] == "1" for r in rs)
            okc += c
            okx += x
            oke += c or x
            good = [int(r["support"]) for r in rs
                    if r["ok_c"] == "1" or r["ok_x"] == "1"]
            if good:
                firsts.append(min(good))
        n = len(byface)
        f50 = quantiles(firsts)[0] if firsts else float("nan")
        print(f"    {t:9} {n:7d} {pct(okc,n):>10} {pct(okx,n):>14}"
              f" {pct(oke,n):>9}   {f50:28.0f}")


def readk_of(evl):
    return "read_c" if evl == "cen" else "read_x"


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tmp/kin.tsv")
