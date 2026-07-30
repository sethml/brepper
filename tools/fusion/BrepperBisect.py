# Which line in the batch's export block broke the refinement? (Phase 0.5a)
#
# BrepperProbe settled that presets DO work: on sphere_25, Low/Medium/High give
# 168 / 528 / 1224 triangles, identically in a shared document and in a fresh one,
# with the deviation readbacks confirming each assignment landed (Low ->
# surfaceDeviation 0.00398 mm / normalDeviation 30 deg, High -> 0.000501 mm / 10 deg).
#
# So the batch's identical meshes were caused by the batch's own code. The batch
# exported sphere_25 at 9024 triangles — finer than Fusion's own High — at all
# three requested presets. Its export block differed from the probe's working
# block in exactly two ways:
#
#   1. it assigned `stl.filename = out_path` AFTER setting meshRefinement
#      (the probe passes the path to createSTLExportOptions only), and
#   2. it read three deviation properties in one try block, the third of which
#      raised AttributeError because Fusion spells it `maximumEdgeLength`.
#
# This API is known to have assignment side effects — Autodesk documents that
# setting surfaceDeviation forces meshRefinement to Custom — so a property
# assignment resetting the mesh settings is entirely plausible, and a read is
# not. This script bisects the two candidates instead of assuming.
#
# All five cases request LOW, whose correct answer (168) is unmistakable against
# the batch's 9024. The deviation values are logged immediately before execute,
# so a reset shows up directly rather than being inferred from a triangle count.
#
# Run: Fusion → Utilities → Add-Ins → Scripts → BrepperBisect → Run.
# Writes REPO/tmp/fusion_bisect.txt. Nothing touches tests/.

import os
import traceback

import adsk.core
import adsk.fusion

REPO = "/Users/sethml/src/brepper"
STEP = os.path.join(REPO, "tests", "onshape", "sphere_25.step")
OUT_DIR = os.path.join(REPO, "tmp", "fusion_bisect")
LOG = os.path.join(REPO, "tmp", "fusion_bisect.txt")

EXPECTED_LOW = 168      # measured by BrepperProbe on this exact model and build
BATCH_RESULT = 9024     # what the broken batch produced for sphere_25


def describe(stl, units):
    out = []
    try:
        out.append("refinement={}".format(int(stl.meshRefinement)))
    except Exception as exc:
        out.append("refinement=err({})".format(type(exc).__name__))
    for label, getter, is_len in (
        ("surfDev", lambda: stl.surfaceDeviation, True),
        ("normDev", lambda: stl.normalDeviation, False),
        ("maxEdge", lambda: stl.maximumEdgeLength, True),
    ):
        try:
            v = getter()
            if is_len:
                try:
                    v = units.convert(v, "cm", "mm")
                except Exception:
                    v = v * 10.0
                out.append("{}={:.6g}mm".format(label, v))
            else:
                out.append("{}={:.6g}".format(label, v))
        except Exception as exc:
            out.append("{}=err({})".format(label, type(exc).__name__))
    return " ".join(out)


def tri_count(path):
    return max(0, (os.path.getsize(path) - 84) // 50) if os.path.exists(path) else -1


def main():
    app = adsk.core.Application.get()
    ui = app.userInterface
    os.makedirs(OUT_DIR, exist_ok=True)
    log = ["brepper Fusion bisect: which line broke the refinement?",
           "Fusion version: {}".format(app.version),
           "all five cases request LOW; correct = {} triangles, "
           "broken batch gave {}".format(EXPECTED_LOW, BATCH_RESULT), ""]

    options = app.importManager.createSTEPImportOptions(STEP)
    doc = app.importManager.importToNewDocument(options)
    if doc is None:
        doc = app.activeDocument
    design = adsk.fusion.Design.cast(app.activeProduct)
    results = []
    try:
        mgr, units = design.exportManager, design.unitsManager
        low = adsk.fusion.MeshRefinementSettings.MeshRefinementLow

        def case(name, build):
            path = os.path.join(OUT_DIR, "{}.stl".format(name))
            stl = mgr.createSTLExportOptions(design.rootComponent, path)
            build(stl, path)
            log.append("  {:26s} before execute: {}".format(name, describe(stl, units)))
            mgr.execute(stl)
            n = tri_count(path)
            results.append((name, n))
            log.append("  {:26s} -> {} triangles{}".format(
                name, n, "  <-- BROKEN" if n != EXPECTED_LOW else "  ok"))

        # 1. The probe's working sequence, as a control.
        def c1(stl, path):
            stl.meshRefinement = low
            stl.isBinaryFormat = True
        case("1_control", c1)

        # 2. Add only the filename assignment, after the refinement.
        def c2(stl, path):
            stl.meshRefinement = low
            stl.isBinaryFormat = True
            stl.filename = path
        case("2_filename_after", c2)

        # 3. Add only the failing property read (no filename assignment).
        def c3(stl, path):
            stl.meshRefinement = low
            stl.isBinaryFormat = True
            try:
                _ = stl.surfaceDeviation
                _ = stl.normalDeviation
                _ = stl.maxEdgeLength      # AttributeError, as in the batch
            except Exception:
                pass
        case("3_failed_read", c3)

        # 4. Exactly what the batch did: both.
        def c4(stl, path):
            stl.meshRefinement = low
            stl.isBinaryFormat = True
            stl.filename = path
            try:
                _ = stl.surfaceDeviation
                _ = stl.normalDeviation
                _ = stl.maxEdgeLength
            except Exception:
                pass
        case("4_batch_exact", c4)

        # 5. filename first, refinement after — the fix if ordering is the issue.
        def c5(stl, path):
            stl.filename = path
            stl.isBinaryFormat = True
            stl.meshRefinement = low
        case("5_filename_first", c5)
    finally:
        if doc is not None:
            try:
                doc.close(False)
            except Exception:
                pass

    log.append("")
    broken = [n for n, c in results if c != EXPECTED_LOW]
    ok = [n for n, c in results if c == EXPECTED_LOW]
    log.append("RESULT  correct: {}".format(", ".join(ok) or "none"))
    log.append("        broken:  {}".format(", ".join(broken) or "none"))
    log.append("")
    if "2_filename_after" in broken and "3_failed_read" in ok:
        log.append("VERDICT: assigning `filename` AFTER meshRefinement resets the mesh "
                   "settings. Fix: pass the path to createSTLExportOptions only, never "
                   "assign filename afterwards"
                   + (" (case 5 shows assigning it FIRST is safe)."
                      if "5_filename_first" in ok else "."))
    elif "3_failed_read" in broken and "2_filename_after" in ok:
        log.append("VERDICT: the failing property read mutates the options object — "
                   "reading a nonexistent attribute is not side-effect free on this "
                   "binding. Fix: never touch maxEdgeLength; read maximumEdgeLength.")
    elif "4_batch_exact" in broken and not broken[:-1]:
        log.append("VERDICT: only the two together break it; investigate the interaction.")
    elif not broken:
        log.append("VERDICT: none of the five reproduces the batch failure, so the cause "
                   "is elsewhere in the batch — most likely per-model rather than "
                   "per-export. Stop and re-read the batch loop.")
    else:
        log.append("VERDICT: see the per-case table above.")

    with open(LOG, "w") as fh:
        fh.write("\n".join(log) + "\n")
    ui.messageBox("\n".join(log[-7:]) + "\n\nFull log: {}".format(LOG))


def run(context):
    try:
        main()
    except Exception:
        app = adsk.core.Application.get()
        if app and app.userInterface:
            app.userInterface.messageBox("BrepperBisect failed:\n{}".format(
                traceback.format_exc()))
