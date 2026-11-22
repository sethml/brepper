#!/usr/bin/env python3
"""The scoreboard: what brepper produces from each corpus mesh, judged against
that mesh and nothing else.

This is the yardstick for the plan at the top of TODO.md. It replaces
tools/sweep.sh as the thing progress is measured by, for one reason: the sweep
scores a run by whether our edges, faces and volume match a REFERENCE STEP at
every stage, and those gates failed good output and passed bad (a re-framed
reference, a trimmed-face artifact, a cone recognized as a sphere). Here the
verdict on the geometry comes from brepper's own self-check against the input
mesh. A reference, where one exists, contributes labels only.

Each model is run twice:

  1. `brepper <stl> -o <step> -v` with no reference. The status, the runtime
     and the self-check numbers come from this run, so the runtime is the
     runtime a user would see.
  2. If the model has a reference: `brepper <stl> --compare <ref>
     --compare-labels-only --reference-report --stage 2`. This contributes the
     reference's face types and three label scores for our stage-2
     decomposition (coverage, purity, type agreement). Stage 2 is
     deterministic, so these describe the decomposition run 1 went on to build.

status (run 1):
  good        valid closed solid, every mesh vertex on it, every facet covered
  unfaithful  valid solid, but not the mesh (or it declares itself too loose)
  invalid     no solid, or BRepCheck fails
  fail        no output at all (the note says which stage gave up)
  timeout     did not finish in --timeout seconds

Usage:
  tools/scoreboard.py [--filter REGEX] [--jobs N] [--timeout SECS] [--out FILE]
  tools/scoreboard.py --baseline        run, then save as the committed baseline
  tools/scoreboard.py --diff            run, then diff against the baseline
  tools/scoreboard.py --diff-only A B   diff two result files
  tools/scoreboard.py --summary FILE    summarize an existing result file

The model list and the model -> reference pairing are tools/sweep.sh's
(`--list`), so the two instruments cover the same corpus. "wild" means
tests/external and tests/downloaded: meshes exported by someone else.

RUNTIME CAVEAT: with --jobs > 1 the secs column is inflated by contention
(brepper is itself multi-threaded). Use --jobs 1 to measure runtime.
"""

import argparse
import concurrent.futures
import datetime
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BREPPER = os.environ.get("BREPPER", os.path.join(ROOT, "target/release/brepper"))
BASELINE = os.path.join(ROOT, "tests", "scoreboard-baseline.tsv")
OUT_DIR = os.path.join(ROOT, "tmp", "scoreboard")
COLS = ["model", "status", "secs", "tris", "faces", "types", "vtx_max", "vtx_off",
        "uncovered", "vol_rel", "tol_max", "ref", "ref_faces", "ref_types",
        "coverage", "purity", "type_agree", "note"]
TYPE_RE = r"(\d+) plane, (\d+) cyl, (\d+) cone, (\d+) sphere, (\d+) torus, (\d+) (?:bspline|freeform)"


def models():
    """[(relative stl path, reference step or None)] from tools/sweep.sh --list."""
    out = subprocess.run([os.path.join(ROOT, "tools/sweep.sh"), "--list"],
                         capture_output=True, text=True, cwd=ROOT).stdout
    result = []
    for line in out.splitlines():
        # Names contain spaces ("OrcaCube v2.stl"), so split on the extension.
        m = re.match(r"(.*\.stl)\s*(.*)$", line)
        if not m:
            continue
        ref = m.group(2).strip()
        result.append((m.group(1), ref if ref.startswith("tests/") else None))
    return result


def run(cmd, timeout):
    """(return code or None on timeout, combined output, seconds).

    Timed on the monotonic clock, the one subprocess's timeout uses: on macOS
    it stops while the machine sleeps, so a run that spans a sleep reports the
    time it actually ran (wall-clock time reported 895 s for a 4.5 s run)."""
    t0 = time.monotonic()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT, timeout=timeout)
        return p.returncode, p.stdout + p.stderr, time.monotonic() - t0
    except subprocess.TimeoutExpired as e:
        text = (e.stdout or b"") + (e.stderr or b"")
        return None, text.decode(errors="replace") if isinstance(text, bytes) else text, time.monotonic() - t0


def types(text, label):
    m = re.search(label + r": " + TYPE_RE, text)
    return ":".join(m.groups()) if m else "-"


def score(model, ref, timeout):
    stl = os.path.join("tests", model)
    key = model[:-4].replace("/", "__")
    step = os.path.join(OUT_DIR, key + ".step")
    if os.path.exists(step):
        os.remove(step)
    row = {c: "-" for c in COLS}
    row.update(model=model[:-4], ref="none", note="")

    # --- Run 1: no reference. ---
    rc, text, secs = run([BREPPER, stl, "-o", step, "-v"], timeout)
    row["secs"] = f"{secs:.2f}"
    m = re.search(r"Read \d+ vertices, (\d+) faces", text)
    if m:
        row["tris"] = m.group(1)
    row["types"] = types(text, "Stage 2 types")
    check = re.search(r"Self-check: (status=.*)", text)
    if rc is None:
        row["status"] = "timeout"
    elif check:
        kv = dict(p.split("=", 1) for p in check.group(1).split() if "=" in p)
        row["status"] = kv.get("status", "fail")
        for col, k in (("faces", "faces"), ("vtx_max", "vtx_max"), ("vtx_off", "vtx_off"),
                       ("uncovered", "uncovered"), ("vol_rel", "vol_rel"), ("tol_max", "tol_max")):
            row[col] = kv.get(k, "-")
    else:
        row["status"] = "fail"
    if row["status"] in ("fail", "invalid", "unfaithful"):
        err = [l for l in text.splitlines() if l.startswith("Error") or "panicked" in l]
        row["note"] = re.sub(r"\s+", " ", err[0])[:160] if err else ""

    # --- Run 2: labels from the reference. ---
    # Not after a timeout: there is no decomposition to label, and the run
    # would only time out again, at three times the cost.
    if ref and row["status"] == "timeout":
        row["ref"] = "skipped"
    elif ref:
        # Stage 2 only, but with three passes of nearest-reference-face
        # queries on top, so it gets more time than the conversion itself.
        rc, text, _ = run([BREPPER, stl, "--compare", ref, "--compare-labels-only",
                           "--reference-report", "--stage", "2", "-v"], timeout * 3)
        if rc is None:
            row["ref"] = "timeout"
        elif "stage 1.2 (compare)" in text:
            # The mesh is not on this reference: not a pair, so no labels.
            row["ref"] = "not-a-pair"
        elif "Reference types:" not in text:
            # The labelling run died first (the same early failure as run 1).
            row["ref"] = "failed"
        else:
            row["ref"] = "ok"
            row["ref_types"] = types(text, "Reference types")
            m = re.search(r"Reference types: .*\((\d+) total\)", text)
            if m:
                row["ref_faces"] = m.group(1)
            m = re.search(r"ref coverage ([\d.]+) .* surface purity ([\d.]+)", text)
            if m:
                row["coverage"], row["purity"] = m.groups()
            m = re.search(r"type agreement ([\d.]+)", text)
            if m:
                row["type_agree"] = m.group(1)
    return row


def read(path):
    rows = {}
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("#")]
    cols = lines[0].split("\t")
    for line in lines[1:]:
        vals = line.split("\t")
        vals += [""] * (len(cols) - len(vals))
        r = dict(zip(cols, vals))
        rows[r["model"]] = r
    return rows


def is_wild(model):
    return model.startswith(("external/", "downloaded/"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def summarize(rows, out=sys.stdout):
    order = ["good", "unfaithful", "invalid", "fail", "timeout"]
    for name, pick in (("wild", is_wild), ("synthetic", lambda m: not is_wild(m))):
        group = [r for r in rows.values() if pick(r["model"])]
        if not group:
            continue
        counts = {s: sum(1 for r in group if r["status"] == s) for s in order}
        good = [r for r in group if r["status"] == "good"]
        labelled = [r for r in good if num(r["type_agree"]) is not None]
        # 0.9995 = agreement on all but a twentieth of a percent of the area:
        # the label for a facet whose samples straddle two reference faces is
        # dropped, so 1.0000 exactly is not always reachable.
        match = [r for r in labelled if num(r["type_agree"]) >= 0.9995]
        print(f"{name}: {len(group)} models: " +
              ", ".join(f"{counts[s]} {s}" for s in order), file=out)
        if labelled:
            mean = sum(num(r["type_agree"]) for r in labelled) / len(labelled)
            print(f"  of the {len(good)} good, {len(labelled)} have labels: "
                  f"{len(match)} match the reference's face types, "
                  f"mean type agreement {mean:.4f}", file=out)


def diff(old, new, out=sys.stdout):
    changed = 0
    for model in sorted(set(old) | set(new)):
        a, b = old.get(model), new.get(model)
        if a is None or b is None:
            print(f"  {'added' if a is None else 'removed'}: {model}", file=out)
            changed += 1
            continue
        notes = []
        if a["status"] != b["status"]:
            notes.append(f"status {a['status']} -> {b['status']}")
        for col in ("type_agree", "coverage", "purity"):
            x, y = num(a.get(col)), num(b.get(col))
            if x is not None and y is not None and abs(x - y) > 5e-4:
                notes.append(f"{col} {x:.4f} -> {y:.4f}")
        if a.get("types") != b.get("types") and a["status"] == b["status"]:
            notes.append(f"types {a.get('types')} -> {b.get('types')}")
        x, y = num(a.get("secs")), num(b.get("secs"))
        if x and y and max(x, y) > 5 and (y > 2 * x or x > 2 * y):
            notes.append(f"secs {x:.1f} -> {y:.1f}")
        if notes:
            print(f"  {model}: " + "; ".join(notes), file=out)
            changed += 1
    print(f"{changed} model(s) changed", file=out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--filter", default="", help="regex on the model path")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=120, help="seconds per conversion")
    ap.add_argument("--out", default=os.path.join(ROOT, "tmp", "scoreboard.tsv"))
    ap.add_argument("--baseline", action="store_true", help="save the result as the baseline")
    ap.add_argument("--diff", action="store_true", help="diff the result against the baseline")
    ap.add_argument("--diff-only", nargs=2, metavar=("A", "B"))
    ap.add_argument("--summary", metavar="FILE")
    args = ap.parse_args()

    if args.diff_only:
        diff(read(args.diff_only[0]), read(args.diff_only[1]))
        return
    if args.summary:
        summarize(read(args.summary))
        return

    todo = [(m, r) for m, r in models() if re.search(args.filter, m)]
    os.makedirs(OUT_DIR, exist_ok=True)
    rows = {}
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = {pool.submit(score, m, r, args.timeout): m for m, r in todo}
        for fut in concurrent.futures.as_completed(futures):
            row = fut.result()
            rows[row["model"]] = row
            print(f"  {row['status']:10s} {float(row['secs']):7.2f}s  {row['model']}",
                  file=sys.stderr)

    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                            text=True, cwd=ROOT).stdout.strip()
    dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                           capture_output=True, text=True, cwd=ROOT).stdout.strip()
    header = [f"# brepper scoreboard {datetime.datetime.now():%Y-%m-%d %H:%M}",
              f"# commit {commit}{'-dirty' if dirty else ''}",
              f"# jobs {args.jobs}   models {len(rows)}   timeout {args.timeout}s",
              "# types and ref_types are plane:cyl:cone:sphere:torus:freeform face counts"]
    body = ["\t".join(COLS)] + ["\t".join(str(rows[m][c]) for c in COLS) for m in sorted(rows)]
    for path in [args.out] + ([BASELINE] if args.baseline else []):
        with open(path, "w") as fh:
            fh.write("\n".join(header + body) + "\n")
        print(f"wrote {os.path.relpath(path, ROOT)}", file=sys.stderr)

    summarize(rows)
    if args.diff:
        if os.path.exists(BASELINE):
            print(f"against {os.path.relpath(BASELINE, ROOT)}:")
            diff(read(BASELINE), rows)
        else:
            print("no baseline to diff against; run with --baseline first")


if __name__ == "__main__":
    main()
