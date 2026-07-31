#!/usr/bin/env python3
"""Harvest wild STEP+STL pairs from open-source hardware repos (Phase 0.5b).

The point of 0.5b is pairs WE DID NOT EXPORT. Our own corpora verify at 31/31
by construction, so they cannot measure whether the provenance premise — mesh
vertices lie on the true CAD surface — survives contact with the wild. A pair
published by a third party can fail for reasons we never see in-house: the STL
was decimated for printing, the mesh was repaired, or the CAD was updated and
the mesh never re-exported. Measuring how often that happens is the deliverable.

GitHub is used rather than a model-sharing site because open-source hardware
repos commit STEP and STL side by side, from one author and one revision, under
a licence that is explicit and machine-readable. `gh` supplies auth.

Discovery walks each repo's git tree in ONE API call rather than crawling
directories, then pairs by stem. Pairing rule is recorded per pair so a
surprising verification result can be audited:
  same-dir   both files in one directory (jig.step + jig.stl)
  sibling    stem matches across sibling directories (STEP/lock.step +
             STL/lock.stl), which is the other common layout
A stem that matches more than once in a repo is skipped rather than guessed.

Writes to tests/external/ (gitignored) plus a committed manifest recording URL,
repo, author, licence, path, size and sha256 for every file taken.
"""

import argparse
import hashlib
import json
import math
import os
import re
import struct
import subprocess
import sys
import urllib.parse
from collections import defaultdict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(REPO_ROOT, "tests", "external")
MANIFEST = os.path.join(OUT_DIR, "manifest.tsv")

# Unrelated projects, so no single CAD workflow dominates the sample. Found via
# `gh search code --extension=step` across several unrelated query terms.
CANDIDATE_REPOS = [
    "bit-bots/bitbots_3d_print",         # humanoid robotics (university team)
    "VoronDesign/VoronUsers",            # 3D printer mods, hundreds of authors
    "e3donline/ToolChanger",             # printer toolchanger (E3D, commercial)
    "ob-f/OpenBot",                      # robotics (Intel research)
    "NAK3DDesigns/White-Knight",         # printer design
    "OpenAstroTech/OpenAstroTracker",    # astronomy mount
    "Rat-Rig/V-core-3",                  # printer (commercial)
    "NikodemBartnik/IndyMill",           # CNC mill
    "dmcke5/NucDeck",                    # handheld console enclosure
    "stevenbennett/Open-Task-Light",     # lighting
]

# Licences whose terms permit redistributing the files, so a pair could be
# committed rather than only analysed locally. Anything else is fetched for
# local analysis only and flagged in the manifest.
REDISTRIBUTABLE = {
    "MIT", "Apache-2.0", "BSD-2-Clause", "BSD-3-Clause", "Unlicense", "CC0-1.0",
    "GPL-2.0", "GPL-3.0", "LGPL-2.1", "LGPL-3.0", "AGPL-3.0", "MPL-2.0",
    "CC-BY-4.0", "CC-BY-SA-4.0", "CERN-OHL-S-2.0", "CERN-OHL-W-2.0",
}

MAX_FILE_BYTES = 20 * 1024 * 1024


def gh_api(path, jq=None):
    cmd = ["gh", "api", path]
    if jq:
        cmd += ["--jq", jq]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        return None
    return out.stdout


def repo_info(repo):
    raw = gh_api(f"repos/{repo}")
    if raw is None:
        return None
    d = json.loads(raw)
    return {
        "repo": repo,
        "owner": d["owner"]["login"],
        "license": (d.get("license") or {}).get("spdx_id") or "NONE",
        "branch": d["default_branch"],
        "stars": d["stargazers_count"],
    }


def find_pairs(repo, branch):
    """Every (step, stl) pair in the repo, by stem, with the rule that matched."""
    raw = gh_api(f"repos/{repo}/git/trees/{branch}?recursive=1")
    if raw is None:
        return []
    tree = json.loads(raw).get("tree", [])
    steps, stls = defaultdict(list), defaultdict(list)
    for node in tree:
        if node.get("type") != "blob":
            continue
        path = node["path"]
        stem, _, ext = path.rpartition("/")[2].rpartition(".")
        ext = ext.lower()
        key = stem.lower()
        if ext in ("step", "stp"):
            steps[key].append(node)
        elif ext == "stl":
            stls[key].append(node)

    pairs = []
    for stem in sorted(set(steps) & set(stls)):
        # Ambiguous stems are skipped, not guessed: two different parts both
        # called "base" would otherwise be paired across mods.
        if len(steps[stem]) != 1 or len(stls[stem]) != 1:
            continue
        step, stl = steps[stem][0], stls[stem][0]
        sdir, mdir = step["path"].rpartition("/")[0], stl["path"].rpartition("/")[0]
        if sdir == mdir:
            rule = "same-dir"
        elif sdir.rpartition("/")[0] == mdir.rpartition("/")[0]:
            rule = "sibling"
        else:
            continue
        if step.get("size", 0) > MAX_FILE_BYTES or stl.get("size", 0) > MAX_FILE_BYTES:
            continue
        pairs.append({"stem": stem, "step": step, "stl": stl, "rule": rule})
    return pairs


def download(repo, branch, path, dest):
    """Fetch one blob. CAD paths routinely contain spaces, plus signs and other
    URL-significant characters ("Buttons - STEP/+ Volume Button.STEP"), so the
    path must be percent-encoded; an unencoded URL fails and, if the failure is
    swallowed, four repos vanish from the sample without a word."""
    quoted = urllib.parse.quote(path, safe="/")
    url = f"https://raw.githubusercontent.com/{repo}/{branch}/{quoted}"
    out = subprocess.run(
        ["curl", "-sSL", "--fail", "--max-time", "120", "-o", dest, url],
        capture_output=True, text=True)
    if out.returncode != 0:
        print(f"    DOWNLOAD FAILED {path}: curl exit {out.returncode}", file=sys.stderr)
        return False
    if not os.path.exists(dest) or os.path.getsize(dest) == 0:
        print(f"    DOWNLOAD EMPTY {path}", file=sys.stderr)
        return False
    return True


def mesh_adequacy(path):
    """(triangles, smooth %, facets per 90 deg of arc) for an STL.

    A pair can be impeccable ground truth and still be a weak test: verification
    says the vertices lie ON the surface, not that there is enough mesh to
    RECOVER it. The first harvest took the first three pairs per repo
    alphabetically and landed a 52-triangle part with no curved surface at all
    and three tessellated at 18 degrees -- so the models that passed were the
    ones with least to get wrong. These two numbers are the selection criteria
    instead.

    `facets per 90 deg` is the reciprocal of the median smooth dihedral, scaled:
    it says how many facets span a quarter turn, which is what decides whether a
    cylinder or torus is recoverable at all. `smooth %` is the share of interior
    edges under 20 degrees, i.e. how much of the mesh sits on curved surface.
    """
    d = open(path, "rb").read()
    if len(d) < 84:
        return 0, 0.0, float("nan")
    n = struct.unpack("<I", d[80:84])[0]
    if len(d) == 84 + 50 * n:
        tris = [tuple(struct.unpack_from("<3f", d, 84 + 50 * i + 12 + 12 * k)
                      for k in range(3)) for i in range(n)]
    else:
        v = [tuple(map(float, m.groups()))
             for m in re.finditer(rb"vertex\s+(\S+)\s+(\S+)\s+(\S+)", d)]
        tris = [tuple(v[i:i + 3]) for i in range(0, len(v) - 2, 3)]
    if not tris:
        return 0, 0.0, float("nan")

    normals = []
    for a, b, c in tris:
        u = [b[i] - a[i] for i in range(3)]
        w = [c[i] - a[i] for i in range(3)]
        cr = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
              u[0] * w[1] - u[1] * w[0]]
        L = math.sqrt(sum(x * x for x in cr))
        normals.append([x / L for x in cr] if L > 0 else None)

    edges = defaultdict(list)
    for i, t in enumerate(tris):
        for k in range(3):
            key = tuple(sorted([tuple(round(x, 4) for x in t[k]),
                                tuple(round(x, 4) for x in t[(k + 1) % 3])]))
            edges[key].append(i)
    dihedrals = []
    for faces in edges.values():
        if len(faces) != 2:
            continue
        n1, n2 = normals[faces[0]], normals[faces[1]]
        if not n1 or not n2:
            continue
        dot = max(-1.0, min(1.0, sum(n1[i] * n2[i] for i in range(3))))
        dihedrals.append(math.degrees(math.acos(dot)))
    if not dihedrals:
        return len(tris), 0.0, float("nan")
    smooth = sorted(d for d in dihedrals if 0.5 < d < 20.0)
    pct = 100.0 * len(smooth) / len(dihedrals)
    if not smooth:
        return len(tris), pct, float("nan")
    median = smooth[len(smooth) // 2]
    return len(tris), pct, (90.0 / median if median > 0 else float("nan"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-repo", type=int, default=3,
                    help="cap per repo so one project cannot dominate the sample")
    ap.add_argument("--target", type=int, default=20, help="total pairs wanted")
    ap.add_argument("--min-facets-per-90", type=float, default=20.0,
                    help="reject meshes too coarse to constrain a surface fit "
                         "(20 = the GOOD grade in tests/external/manifest.tsv)")
    ap.add_argument("--min-smooth-pct", type=float, default=10.0,
                    help="reject meshes with too little curved surface to test")
    ap.add_argument("--max-candidates", type=int, default=25,
                    help="STLs to try per repo before giving up, bounding "
                         "bandwidth spent on rejects")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    rows, sources, taken = [], [], 0

    for repo in CANDIDATE_REPOS:
        if taken >= args.target:
            break
        info = repo_info(repo)
        if info is None:
            print(f"  {repo}: unreachable", file=sys.stderr)
            continue
        pairs = find_pairs(repo, info["branch"])
        if not pairs:
            print(f"  {repo}: licence={info['license']} — no stem-matched pairs",
                  file=sys.stderr)
            continue
        redist = info["license"] in REDISTRIBUTABLE
        print(f"  {repo}: licence={info['license']} "
              f"({'redistributable' if redist else 'LOCAL ANALYSIS ONLY'}), "
              f"{len(pairs)} pairs, taking {min(args.per_repo, len(pairs))}",
              file=sys.stderr)
        kept = tried = rejected = 0
        for pair in pairs:
            if taken >= args.target or kept >= args.per_repo:
                break
            if tried >= args.max_candidates:
                break
            slug = f"{repo.split('/')[1]}__{pair['stem']}".replace(" ", "_")[:80]
            if args.dry_run:
                taken += 1
                kept += 1
                continue
            if os.path.exists(os.path.join(OUT_DIR, f"{slug}.stl")):
                continue  # already harvested by an earlier run
            tried += 1

            # The STL is fetched and graded FIRST: a reject then costs one file
            # instead of two, and the grade is a property of the mesh alone.
            stl_dest = os.path.join(OUT_DIR, f"{slug}.stl")
            if not download(repo, info["branch"], pair["stl"]["path"], stl_dest):
                continue
            tris, smooth, per90 = mesh_adequacy(stl_dest)
            too_coarse = not (per90 == per90) or per90 < args.min_facets_per_90
            if too_coarse or smooth < args.min_smooth_pct:
                why = ("no curved surface" if not (per90 == per90)
                       else f"{per90:.1f} facets/90deg" if too_coarse
                       else f"{smooth:.1f}% smooth")
                print(f"    reject {slug}: {why}", file=sys.stderr)
                os.remove(stl_dest)
                rejected += 1
                continue

            step_dest = os.path.join(OUT_DIR, f"{slug}.step")
            if not download(repo, info["branch"], pair["step"]["path"], step_dest):
                # Half a pair is worse than none: it would be scored against
                # whatever reference happened to remain.
                os.remove(stl_dest)
                continue

            for kind, node, dest in (("step", pair["step"], step_dest),
                                     ("stl", pair["stl"], stl_dest)):
                with open(dest, "rb") as fh:
                    digest = hashlib.sha256(fh.read()).hexdigest()
                rows.append({
                    "file": f"{slug}.{kind}",
                    "repo": repo,
                    "owner": info["owner"],
                    "license": info["license"],
                    "redistributable": "yes" if redist else "no",
                    "path": node["path"],
                    "url": f"https://github.com/{repo}/blob/{info['branch']}/{node['path']}",
                    "bytes": str(os.path.getsize(dest)),
                    "sha256": digest,
                    "pair_rule": pair["rule"],
                    "tris": str(tris) if kind == "stl" else "",
                    "smooth_pct": f"{smooth:.1f}" if kind == "stl" else "",
                    "facets_per_90": f"{per90:.1f}" if kind == "stl" else "",
                })
            sources.append(f"{slug}\ttests/external/{slug}.step")
            print(f"    keep   {slug}: {tris} tris, {smooth:.1f}% smooth, "
                  f"{per90:.1f} facets/90deg", file=sys.stderr)
            taken += 1
            kept += 1
        if not args.dry_run:
            print(f"    -> kept {kept}, rejected {rejected} of {tried} tried",
                  file=sys.stderr)

    if args.dry_run:
        print(f"dry run: {taken} pairs available", file=sys.stderr)
        return

    cols = ["file", "repo", "owner", "license", "redistributable", "path", "url",
            "bytes", "sha256", "pair_rule", "tris", "smooth_pct", "facets_per_90"]
    # Merge with whatever a previous run recorded, then emit one row per file
    # actually on disk. Rewriting from this run's rows alone would discard the
    # URL and sha256 of every earlier pair -- the only reason a deleted model can
    # be recovered, which is not hypothetical: seven were, once.
    merged = {}
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as fh:
            lines = [l.rstrip("\n") for l in fh if not l.startswith("#")]
        if lines:
            old_cols = lines[0].split("\t")
            for line in lines[1:]:
                vals = line.split("\t")
                if len(vals) < len(old_cols):
                    vals += [""] * (len(old_cols) - len(vals))
                merged[vals[0]] = dict(zip(old_cols, vals))
    for r in rows:
        merged[r["file"]] = r

    final = []
    for name in sorted(merged):
        path = os.path.join(OUT_DIR, name)
        if not os.path.exists(path):
            continue                      # recorded but no longer on disk
        row = {c: merged[name].get(c, "") for c in cols}
        if name.endswith(".stl") and not row["facets_per_90"]:
            tris, smooth, per90 = mesh_adequacy(path)
            row["tris"] = str(tris)
            row["smooth_pct"] = f"{smooth:.1f}"
            row["facets_per_90"] = f"{per90:.1f}" if per90 == per90 else "-"
        final.append(row)

    with open(MANIFEST, "w") as fh:
        fh.write("# Wild STEP+STL pairs harvested from GitHub (Phase 0.5b).\n")
        fh.write("# Model files are gitignored; this manifest is the record, and\n")
        fh.write("# url + sha256 are what make a deleted pair recoverable.\n")
        fh.write("#\n")
        fh.write("# facets_per_90: facets spanning a quarter turn of curvature.\n")
        fh.write("# smooth_pct: share of interior edges under 20 deg.\n")
        fh.write("# Below ~20 facets/90deg a mesh is weak test material however\n")
        fh.write("# well it verifies -- passing then means having little to get\n")
        fh.write("# wrong. Do not quote a pass rate over the coarse rows.\n")
        fh.write("\t".join(cols) + "\n")
        for r in final:
            fh.write("\t".join(r[c] for c in cols) + "\n")

    stems = sorted({r["file"].rsplit(".", 1)[0] for r in final})
    with open(os.path.join(OUT_DIR, "sources.tsv"), "w") as fh:
        fh.write("# STL stem\treference STEP (…_aligned.step = rigidly re-framed)\n")
        for stem in stems:
            aligned = os.path.join(OUT_DIR, f"{stem}_aligned.step")
            ref = f"{stem}_aligned.step" if os.path.exists(aligned) else f"{stem}.step"
            fh.write(f"{stem}\ttests/external/{ref}\n")
    print(f"{taken} new pairs; {len(stems)} total -> {OUT_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
