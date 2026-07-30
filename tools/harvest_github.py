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
import os
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-repo", type=int, default=3,
                    help="cap per repo so one project cannot dominate the sample")
    ap.add_argument("--target", type=int, default=20, help="total pairs wanted")
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
        for pair in pairs[: args.per_repo]:
            if taken >= args.target:
                break
            slug = f"{repo.split('/')[1]}__{pair['stem']}".replace(" ", "_")[:80]
            if args.dry_run:
                taken += 1
                continue
            ok = True
            for kind, node in (("step", pair["step"]), ("stl", pair["stl"])):
                dest = os.path.join(OUT_DIR, f"{slug}.{kind}")
                if not download(repo, info["branch"], node["path"], dest):
                    ok = False
                    break
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
                })
            if ok:
                sources.append(f"{slug}\ttests/external/{slug}.step")
                taken += 1
            else:
                # Half a pair is worse than none: it would be scored against
                # whatever reference happened to remain.
                for kind in ("step", "stl"):
                    stale = os.path.join(OUT_DIR, f"{slug}.{kind}")
                    if os.path.exists(stale):
                        os.remove(stale)
                rows = [r for r in rows if not r["file"].startswith(slug + ".")]

    if args.dry_run:
        print(f"dry run: {taken} pairs available", file=sys.stderr)
        return

    cols = ["file", "repo", "owner", "license", "redistributable", "path", "url",
            "bytes", "sha256", "pair_rule"]
    with open(MANIFEST, "w") as fh:
        fh.write("# Wild STEP+STL pairs harvested from GitHub (Phase 0.5b).\n")
        fh.write("# Model files are gitignored; this manifest is the record.\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in cols) + "\n")
    with open(os.path.join(OUT_DIR, "sources.tsv"), "w") as fh:
        fh.write("# STL stem\treference STEP it is paired with\n")
        fh.write("\n".join(sources) + "\n")
    print(f"{taken} pairs -> {OUT_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
