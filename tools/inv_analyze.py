#!/usr/bin/env python3
"""Aggregate tools/invariants.sh output into the per-invariant verdicts.

Reads one or more TSV files produced by tools/invariants.sh and prints, for
each proposed invariant, how well it holds and which models break it. Nothing
here decides anything — it only summarises the measurement.
"""
import sys
import math
from collections import defaultdict

COLS = ("kind model ref surf facets area vmax emax cmax imax edmax dih_in dih_out "
        "sharp_in res_plane res_sphere res_cyl res_torus seeds grew exact distinct "
        "jacc_min over_max under_max seeds_int grew_int exact_int distinct_int seam_group exact_union ctypes ccorrect chist order_diffs oracle "
        "vtol stol quantum").split()
SR = ("kind model ref surf pairs six_verts all_coplanar four_coplanar "
      "ax_p50 ax_p95 ax_max n_admitted adm_p50 adm_max selected_err").split()


def f(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def load(paths):
    inv, sr = [], []
    for p in paths:
        for line in open(p):
            fs = line.rstrip("\n").split("\t")
            if fs[0] == "INV" and len(fs) >= len(COLS):
                inv.append(dict(zip(COLS, fs)))
            elif fs[0] == "SEEDRULE" and len(fs) >= len(SR):
                sr.append(dict(zip(SR, fs)))
    return inv, sr


def q(vals, p):
    v = sorted(x for x in vals if not math.isnan(x))
    if not v:
        return float("nan")
    return v[min(len(v) - 1, int(round((len(v) - 1) * p)))]


def hdr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def main():
    inv, sr = load(sys.argv[1:])
    models = {r["model"] for r in inv}
    print(f"{len(inv)} reference faces over {len(models)} models; "
          f"{len(sr)} developable faces seed-rule-probed")

    analytic = [r for r in inv if not math.isnan(f(r["vmax"]))]

    # ---------------------------------------------------------------- F ----
    hdr("F. Surface Fitting Requirement — is it satisfiable on GROUND TRUTH?")
    by_t = defaultdict(list)
    for r in analytic:
        by_t[r["surf"]].append(r)
    print(f"{'type':10} {'faces':>6} {'vmax/quantum p50':>17} {'p95':>9} {'max':>10} "
          f"{'over vtol':>10}")
    for t, rs in sorted(by_t.items()):
        rel = [f(r["vmax"]) / f(r["quantum"]) for r in rs if f(r["quantum"]) > 0]
        over = sum(1 for r in rs if f(r["vmax"]) > f(r["vtol"]))
        print(f"{t:10} {len(rs):6d} {q(rel,.5):17.2f} {q(rel,.95):9.2f} "
              f"{max(rel, default=float('nan')):10.1f} {over:10d}")

    print("\n  probe adequacy: whole-facet max vs the section's chosen probe")
    ratios_e, ratios_c, ratios_ed = [], [], []
    for r in analytic:
        fm = max(f(r["imax"]), f(r["cmax"]), f(r["emax"]), f(r["edmax"]))
        if f(r["emax"]) > 0:
            ratios_e.append(fm / f(r["emax"]))
        if f(r["edmax"]) > 0:
            ratios_ed.append(fm / f(r["edmax"]))
        if f(r["cmax"]) > 0:
            ratios_c.append(fm / f(r["cmax"]))
    print(f"    facet_max / edge_mid (BOUNDARY only, as the section says): "
          f"p95 {q(ratios_e,.95):.4f}  max {max(ratios_e, default=float('nan')):.4f}")
    print(f"    facet_max / edge_mid (+ quad DIAGONALS): "
          f"p95 {q(ratios_ed,.95):.4f}  max {max(ratios_ed, default=float('nan')):.4f}"
          f"   (theory: <= 1.3333)")
    print(f"    facet_max / centroid     : p50 {q(ratios_c,.5):.4f}  p95 {q(ratios_c,.95):.4f}  "
          f"max {max(ratios_c, default=float('nan')):.4f}")
    pairs = [(max(f(r['imax']), f(r['cmax']), f(r['emax']), f(r['edmax'])) / f(r['edmax']), r)
             for r in analytic if f(r['edmax']) > 0]
    pairs.sort(key=lambda t: -t[0])
    n_over = sum(1 for v, _ in pairs if v > 4 / 3 + 1e-9)
    print(f"    {n_over} of {len(pairs)} faces exceed the 4/3 bound; worst:")
    for v, r in pairs[:6]:
        print(f"      {r['model']}#{r['ref']} {r['surf']} {v:.4f} "
              f"({r['facets']} facets)")

    print("\n  centroid probe blindness (centroid reads ~0 where the facet does not)")
    blind = [r for r in analytic
             if f(r["emax"]) > 0 and f(r["cmax"]) < 0.25 * f(r["emax"])]
    print(f"    {len(blind)} of {len(analytic)} faces have centroid < 1/4 of the edge-mid max")
    for r in blind[:5]:
        print(f"      {r['model']}#{r['ref']} {r['surf']}: "
              f"centroid {f(r['cmax']):.3e} vs edge {f(r['emax']):.3e}")

    print("\n  is --surface-tolerance 0.4 mm binding?")
    over_stol = [r for r in analytic if f(r["edmax"]) > f(r["stol"])]
    allsag = [max(f(r["edmax"]), f(r["imax"])) for r in analytic]
    print(f"    true worst facet sag over the corpus: {max(allsag, default=0):.4e} mm; "
          f"p50 {q(allsag,.5):.3e}, p95 {q(allsag,.95):.3e}")
    print(f"    {len(over_stol)} reference faces exceed stol on GROUND TRUTH "
          f"(a true surface its own tolerance rejects)")
    for r in over_stol[:8]:
        print(f"      {r['model']}#{r['ref']} {r['surf']}: {f(r['emax']):.4e} > {f(r['stol']):.2e}")

    # ---------------------------------------------------------------- A ----
    hdr("A. 'All facet relative angles within angle tolerance'")
    sep = [r for r in inv if f(r["dih_out"]) >= 0 and f(r["dih_in"]) > 0]
    insep = [r for r in sep if f(r["dih_in"]) >= f(r["dih_out"])]
    print(f"  {len(insep)} of {len(sep)} reference faces have an INTERIOR dihedral >= their "
          f"smallest BOUNDARY dihedral.")
    print("  Where that holds, no angular threshold can both keep the surface together and "
          "cut it from its neighbour.")
    worst = sorted(sep, key=lambda r: -(f(r["dih_in"]) / max(f(r["dih_out"]), 1e-9)))[:8]
    for r in worst:
        print(f"    {r['model']}#{r['ref']} {r['surf']:9} interior {f(r['dih_in']):7.3f}deg  "
              f"boundary {f(r['dih_out']):7.4f}deg  ratio {f(r['dih_in'])/max(f(r['dih_out']),1e-9):8.1f}x")
    split = [r for r in inv if int(r.get("sharp_in", 0) or 0) > 0]
    print(f"\n  the pipeline's OWN classifier (stage 2.1) calls "
          f"{sum(int(r['sharp_in']) for r in split)} edge(s) sharp strictly INSIDE a reference "
          f"face, cutting {len(split)} true surfaces into pieces before any fitting:")
    for r in sorted(split, key=lambda r: -int(r["sharp_in"]))[:12]:
        print(f"    {r['model']}#{r['ref']} {r['surf']:9} {r['sharp_in']:>4} sharp interior edges "
              f"({r['facets']} facets, worst interior dihedral {f(r['dih_in']):.2f}deg)")
    over25 = [r for r in inv if f(r["dih_in"]) > 25.0]
    print(f"\n  {len(over25)} reference faces have an interior dihedral over the 25deg default "
          f"--angular-tolerance (the sharp-edge partition SPLITS these true surfaces):")
    for r in over25[:10]:
        print(f"    {r['model']}#{r['ref']} {r['surf']:9} interior {f(r['dih_in']):.2f}deg "
              f"({r['facets']} facets)")

    # ---------------------------------------------------------------- I ----
    hdr("I. Surface Identification — does a WRONG type also explain the true support?")
    idx = {"plane": "res_plane", "sphere": "res_sphere", "cylinder": "res_cyl",
           "cone": "res_cyl", "torus": "res_torus"}
    names = {"res_plane": "plane", "res_sphere": "sphere", "res_cyl": "cyl/cone",
             "res_torus": "torus"}
    own_fail, typed = [], 0
    impostor = defaultdict(lambda: defaultdict(int))
    quanta = defaultdict(list)
    for r in inv:
        own = idx.get(r["surf"])
        if own is None:
            continue
        typed += 1
        vtol, qm = f(r["vtol"]), f(r["quantum"])
        if not (f(r[own]) <= vtol):
            own_fail.append(r)
        for c in ("res_plane", "res_sphere", "res_cyl", "res_torus"):
            v = f(r[c])
            if c == own or math.isnan(v):
                continue
            if v <= vtol:
                impostor[r["surf"]][names[c]] += 1
                if qm > 0:
                    quanta[(r["surf"], names[c])].append(v / qm)
    print(f"  {len(own_fail)} of {typed} reference faces are NOT explained at vtol by a fit of "
          f"their own type:")
    for r in own_fail[:10]:
        own = idx[r["surf"]]
        print(f"    {r['model']}#{r['ref']} {r['surf']:9} own-type residual {f(r[own]):.3e} "
              f"({f(r[own])/max(f(r['quantum']),1e-30):.0f} q, {r['facets']} facets)")
    print("\n  wrong types that ALSO validate the true support, and how far they sit "
          "(in float32 quanta):")
    for t in sorted(impostor):
        n = sum(1 for r in inv if r["surf"] == t)
        for c, k in sorted(impostor[t].items()):
            qs = quanta[(t, c)]
            print(f"    true {t:9} also fits as {c:9}: {k:4d} of {n:4d} faces "
                  f"(residual p50 {q(qs,.5):8.2f} q, max {max(qs, default=0):8.1f} q)")
    print("\n  own-type residual, for contrast (in quanta):")
    for t in sorted(set(r["surf"] for r in inv if r["surf"] in idx)):
        qs = [f(r[idx[t]]) / f(r["quantum"]) for r in inv
              if r["surf"] == t and f(r["quantum"]) > 0 and not math.isnan(f(r[idx[t]]))]
        print(f"    {t:9}: p50 {q(qs,.5):8.2f} q  p95 {q(qs,.95):8.2f} q  max {max(qs, default=0):9.1f} q")

    # ---------------------------------------------------------------- S ----
    hdr("S. Surface Seed Independence")
    tested = [r for r in inv if int(r["seeds"]) > 0]
    seeds = sum(int(r["seeds"]) for r in tested)
    grew = sum(int(r["grew"]) for r in tested)
    exact = sum(int(r["exact"]) for r in tested)
    stable = [r for r in tested if int(r["distinct"]) <= 1 and int(r["grew"]) == int(r["seeds"])]
    print(f"  {len(tested)} reference faces seed-tested with {seeds} seeds")
    print(f"    seeds that grew anything          : {grew:5d} ({100*grew/max(seeds,1):.1f}%)")
    exact_u = sum(int(r["exact_union"]) for r in tested)
    seamed = [r for r in tested if int(r["seam_group"]) > 1]
    print(f"    seeds that recovered the TRUE set : {exact:5d} ({100*exact/max(seeds,1):.1f}%)")
    print(f"    ... or its whole SEAM GROUP       : {exact_u:5d} ({100*exact_u/max(seeds,1):.1f}%)  "
          f"[{len(seamed)} of {len(tested)} faces are one half of a seam-split surface]")
    best = sum(max(int(r["exact"]), int(r["exact_union"])) for r in tested)
    print(f"    either                            : {best:5d} ({100*best/max(seeds,1):.1f}%)")
    print(f"    faces where EVERY seed grew and all seeds agreed: {len(stable)} of {len(tested)} "
          f"({100*len(stable)/max(len(tested),1):.1f}%)")
    byt = defaultdict(lambda: [0, 0, 0, 0])
    for r in tested:
        b = byt[r["surf"]]
        b[0] += int(r["seeds"]); b[1] += int(r["grew"]); b[2] += int(r["exact"])
        b[3] += 1 if (int(r["distinct"]) <= 1 and int(r["grew"]) == int(r["seeds"])) else 0
    print(f"\n  {'type':10} {'faces':>6} {'seeds':>6} {'grew%':>7} {'exact%':>7} {'stable faces':>13}")
    for t in sorted(byt):
        s, g, e, st = byt[t]
        n = sum(1 for r in tested if r["surf"] == t)
        print(f"  {t:10} {n:6d} {s:6d} {100*g/max(s,1):7.1f} {100*e/max(s,1):7.1f} "
              f"{st:6d}/{n:<6d}")
    si = sum(int(r["seeds_int"]) for r in tested)
    gi = sum(int(r["grew_int"]) for r in tested)
    ei = sum(int(r["exact_int"]) for r in tested)
    with_int = [r for r in tested if int(r["seeds_int"]) > 0]
    stable_int = [r for r in with_int
                  if int(r["distinct_int"]) <= 1 and int(r["grew_int"]) == int(r["seeds_int"])]
    print(f"\n  RESTRICTED TO INTERIOR SEEDS (every mesh neighbour in the same reference face):")
    print(f"    {si} of {seeds} seeds are interior")
    print(f"    grew                        : {gi:5d} ({100*gi/max(si,1):.1f}%)   "
          f"[all seeds: {100*grew/max(seeds,1):.1f}%]")
    print(f"    recovered the TRUE face set : {ei:5d} ({100*ei/max(si,1):.1f}%)   "
          f"[all seeds: {100*exact/max(seeds,1):.1f}%]")
    print(f"    faces where every interior seed grew and agreed: {len(stable_int)} of "
          f"{len(with_int)} ({100*len(stable_int)/max(len(with_int),1):.1f}%)   "
          f"[all seeds: {100*len(stable)/max(len(tested),1):.1f}%]")
    byt2 = defaultdict(lambda: [0, 0, 0, 0, 0])
    for r in with_int:
        b = byt2[r["surf"]]
        b[0] += int(r["seeds_int"]); b[1] += int(r["grew_int"]); b[2] += int(r["exact_int"])
        b[3] += 1
        b[4] += 1 if (int(r["distinct_int"]) <= 1 and int(r["grew_int"]) == int(r["seeds_int"])) else 0
    print(f"\n    {'type':10}{'faces':>7}{'int seeds':>11}{'grew%':>8}{'exact%':>8}{'stable faces':>15}")
    for t in sorted(byt2):
        s2, g2, e2, n2, st2 = byt2[t]
        print(f"    {t:10}{n2:7d}{s2:11d}{100*g2/max(s2,1):8.1f}{100*e2/max(s2,1):8.1f}"
              f"{st2:9d}/{n2:<5d}")
    unstable = sorted(tested, key=lambda r: (f(r["jacc_min"]) if not math.isnan(f(r["jacc_min"])) else 2))
    print("\n  worst seed-dependence (lowest Jaccard of a grown set vs the truth):")
    for r in unstable[:12]:
        if math.isnan(f(r["jacc_min"])):
            continue
        print(f"    {r['model']}#{r['ref']} {r['surf']:9} {r['facets']:>6} facets: "
              f"{r['grew']}/{r['seeds']} grew, {r['exact']} exact, {r['distinct']} distinct, "
              f"jacc {f(r['jacc_min']):.3f}, leak-out {100*f(r['over_max']):.0f}%, "
              f"stop-short {100*f(r['under_max']):.0f}%")
    leaks = [r for r in tested if f(r["over_max"]) > 0.02]
    shorts = [r for r in tested if f(r["under_max"]) > 0.02]
    print(f"\n  {len(leaks)} faces where some seed leaked >2% of its set ACROSS a true boundary")
    print(f"  {len(shorts)} faces where some seed stopped >2% short of the true surface")

    # ---------------------------------------------------------------- G ----
    hdr("G. Growth Independence")
    od = sum(int(r["order_diffs"]) for r in tested)
    orc = sum(1 for r in tested if r["oracle"] == "1")
    print(f"  reversing the flood's neighbour order changed the result on {od} of {seeds} seeds")
    print(f"  starting the refit fixpoint from the WHOLE true face set reproduces the truth on "
          f"{orc} of {len(tested)} faces ({100*orc/max(len(tested),1):.1f}%)")
    disagree = [r for r in tested
                if (r["oracle"] == "1") != (int(r["exact"]) == int(r["seeds"]) and int(r["seeds"]) > 0)]
    print(f"  {len(disagree)} faces where the oracle start and the seed starts disagree "
          f"(the fixpoint is start-dependent):")
    for r in disagree[:10]:
        print(f"    {r['model']}#{r['ref']} {r['surf']:9} oracle={r['oracle']} "
              f"seeds exact {r['exact']}/{r['seeds']}")

    # ------------------------------------------------------- competition ----
    comp = [r for r in inv if r.get("chist", "-") not in ("-", "", None)]
    if comp:
        hdr("C. Competition from a seed inside a known surface (the pipeline's own rule)")
        seeds_c = sum(int(r["seeds"]) for r in comp)
        cc = sum(int(r["ccorrect"]) for r in comp)
        one = sum(1 for r in comp if int(r["ctypes"]) == 1)
        print(f"  {len(comp)} reference faces, {seeds_c} seeds")
        print(f"    seeds whose winner had the reference face's own type: {cc} "
              f"({100*cc/max(seeds_c,1):.1f}%)")
        print(f"    faces where every seed agreed on ONE winning type: {one} of {len(comp)} "
              f"({100*one/max(len(comp),1):.1f}%)")
        byk = defaultdict(lambda: defaultdict(int))
        for r in comp:
            for part in r["chist"].split(","):
                if ":" in part:
                    k, v = part.rsplit(":", 1)
                    byk[r["surf"]][k] += int(v)
        print(f"\n  what the competition actually picks, per true type:")
        for t in sorted(byk):
            tot = sum(byk[t].values())
            got = ", ".join(f"{k} {100*v/tot:.0f}%" for k, v in
                            sorted(byk[t].items(), key=lambda kv: -kv[1]))
            print(f"    true {t:9} (n={tot:5d} winning seeds): {got}")

    # ---------------------------------------------------------------- E ----
    hdr("E. Surface Seed Requirements (the stated cylinder/cone seed rule)")
    if sr:
        tot = sum(int(x["pairs"]) for x in sr)
        six = sum(int(x["six_verts"]) for x in sr)
        allc = sum(int(x["all_coplanar"]) for x in sr)
        four = sum(int(x["four_coplanar"]) for x in sr)
        print(f"  {len(sr)} developable reference faces, {tot} facet chains (<=4 facets) enumerated")
        print(f"    with exactly 6 distinct vertices        : {six} ({100*six/max(tot,1):.1f}%)")
        print(f"    of those, ALL SIX coplanar (degenerate) : {allc} ({100*allc/max(six,1):.2f}%)")
        print(f"    of those, SOME FOUR coplanar            : {four} ({100*four/max(six,1):.2f}%)")
        print("    -> read literally, 'no more than three coplanar' rejects that last group.")
        amax = [f(x["ax_max"]) for x in sr]
        sel = [f(x["selected_err"]) for x in sr]
        print(f"\n  closed-form axis error from a bare 6-vertex chain: "
              f"p50 {q(amax,.5):.4f}deg p95 {q(amax,.95):.4f}deg max {max(amax, default=0):.4f}deg")
        print(f"  same, from the whole face with Phase-1.1 pair SELECTION: "
              f"p50 {q(sel,.5):.4f}deg p95 {q(sel,.95):.4f}deg max {max(sel, default=0):.4f}deg")
        worst = sorted(sr, key=lambda x: -f(x["ax_max"]))[:8]
        print("\n  worst bare-chain axis errors:")
        for x in worst:
            print(f"    {x['model']}#{x['ref']} {x['surf']:9} chain max {f(x['ax_max']):9.4f}deg  "
                  f"selected {f(x['selected_err']):9.4f}deg  ({x['six_verts']} chains)")


if __name__ == "__main__":
    main()
