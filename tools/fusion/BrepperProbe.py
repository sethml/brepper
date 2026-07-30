# Why did Low/Medium/High produce byte-identical STLs? (Phase 0.5a follow-up)
#
# In the 2026-07-29 batch, 8 of 9 multi-preset models produced a BYTE-IDENTICAL
# mesh at all three refinements. Two are legitimate — chamfered_cube (44
# triangles) and stepped_block (52) are entirely planar, so refinement cannot
# change them. The other six have curved faces, sphere_25 among them, and must
# respond. The one model that DID respond, stepped_block_toroidal_fillets
# (2436 / 3702 / 6724 triangles, monotonic in the requested order), was the
# FIRST exported, which points at state rather than geometry.
#
# Autodesk documents that setting meshRefinement "will reset the
# surfaceDeviation, normalDeviation, maximumEdgeLength, and aspectRatio to
# values that correspond to the specified mesh refinement", and that the default
# is MeshRefinementMedium. That makes the deviation values a DIRECT probe of
# whether the assignment took effect, which is the measurement the batch lost:
# its three property reads shared one try block, and Fusion has no
# `maxEdgeLength` (it is `maximumEdgeLength`), so that AttributeError discarded
# surfaceDeviation and normalDeviation with it.
#
# Four experiments on sphere_25, which cannot legitimately be preset-invariant:
#
#   A. PRESET READBACK — set each preset, read meshRefinement and the four
#      deviation values back. If they do not change per preset, the assignment
#      never reached the options object.
#   B. REVERSED ORDER — export high, low, medium (the batch used low, medium,
#      high). If counts differ now, order or first-export state is the cause.
#   C. FRESH DOCUMENT PER EXPORT — re-import for every preset, so no two exports
#      share a document. If this works and B does not, per-document state is it.
#   D. EXPLICIT CUSTOM DEVIATIONS — set surfaceDeviation numerically instead of
#      using a preset. If D varies while A/B/C do not, the preset path is broken
#      and this is the workaround.
#
# Run: Fusion → Utilities → Add-Ins → Scripts → BrepperProbe → Run. Writes
# REPO/tmp/fusion_probe.txt. Nothing goes into tests/, so the committed corpus
# cannot be disturbed.

import os
import traceback

import adsk.core
import adsk.fusion

REPO = "/Users/sethml/src/brepper"
STEP = os.path.join(REPO, "tests", "onshape", "sphere_25.step")
OUT_DIR = os.path.join(REPO, "tmp", "fusion_probe")
LOG = os.path.join(REPO, "tmp", "fusion_probe.txt")

PRESETS_REVERSED = ["high", "low", "medium"]
# Surface deviations in mm, spanning coarse to fine for a 25 mm sphere.
CUSTOM_DEVIATIONS_MM = [0.5, 0.05, 0.005]


def refinement_enum(name):
    s = adsk.fusion.MeshRefinementSettings
    return {"low": s.MeshRefinementLow,
            "medium": s.MeshRefinementMedium,
            "high": s.MeshRefinementHigh,
            "custom": s.MeshRefinementCustom}[name]


def mm(units, value):
    try:
        return units.convert(value, "cm", "mm")
    except Exception:
        return value * 10.0


def describe(stl, units):
    """One try per property: a missing attribute must cost one value, not four."""
    out = []
    try:
        out.append("meshRefinement={}".format(int(stl.meshRefinement)))
    except Exception as exc:
        out.append("meshRefinement=unavailable({})".format(type(exc).__name__))
    for label, getter, is_len in (
        ("surfaceDeviation", lambda: stl.surfaceDeviation, True),
        ("normalDeviation", lambda: stl.normalDeviation, False),
        ("maximumEdgeLength", lambda: stl.maximumEdgeLength, True),
        ("aspectRatio", lambda: stl.aspectRatio, False),
    ):
        try:
            v = getter()
            out.append("{}={:.6g}{}".format(label, mm(units, v) if is_len else v,
                                            "mm" if is_len else ""))
        except Exception as exc:
            out.append("{}=unavailable({})".format(label, type(exc).__name__))
    return " ".join(out)


def tri_count(path):
    if not os.path.exists(path):
        return -1
    return max(0, (os.path.getsize(path) - 84) // 50)


def import_sphere(app):
    options = app.importManager.createSTEPImportOptions(STEP)
    doc = app.importManager.importToNewDocument(options)
    if doc is None:
        doc = app.activeDocument
    design = adsk.fusion.Design.cast(app.activeProduct)
    if design is None:
        raise RuntimeError("imported document is not a Design")
    return doc, design


def main():
    app = adsk.core.Application.get()
    ui = app.userInterface
    os.makedirs(OUT_DIR, exist_ok=True)
    log = ["brepper Fusion refinement probe",
           "Fusion version: {}".format(app.version),
           "model: {}".format(STEP), ""]

    s = adsk.fusion.MeshRefinementSettings
    log.append("enum: low={} medium={} high={} custom={}".format(
        int(s.MeshRefinementLow), int(s.MeshRefinementMedium),
        int(s.MeshRefinementHigh), int(s.MeshRefinementCustom)))
    log.append("")

    results = {}

    # --- A and B: presets in reversed order, one document ---------------------
    log.append("A/B: presets in REVERSED order, sharing one document")
    doc, design = import_sphere(app)
    try:
        mgr, units = design.exportManager, design.unitsManager
        counts = []
        for preset in PRESETS_REVERSED:
            path = os.path.join(OUT_DIR, "shared_{}.stl".format(preset))
            stl = mgr.createSTLExportOptions(design.rootComponent, path)
            stl.meshRefinement = refinement_enum(preset)
            stl.isBinaryFormat = True
            log.append("  requested {:6s} -> {}".format(preset, describe(stl, units)))
            mgr.execute(stl)
            counts.append((preset, tri_count(path)))
            log.append("  {:6s} exported {} triangles".format(preset, counts[-1][1]))
        results["shared_document"] = counts
    finally:
        doc.close(False)
    log.append("")

    # --- C: a fresh import per export ---------------------------------------
    log.append("C: FRESH document per export")
    counts = []
    for preset in PRESETS_REVERSED:
        doc, design = import_sphere(app)
        try:
            mgr, units = design.exportManager, design.unitsManager
            path = os.path.join(OUT_DIR, "fresh_{}.stl".format(preset))
            stl = mgr.createSTLExportOptions(design.rootComponent, path)
            stl.meshRefinement = refinement_enum(preset)
            stl.isBinaryFormat = True
            log.append("  requested {:6s} -> {}".format(preset, describe(stl, units)))
            mgr.execute(stl)
            counts.append((preset, tri_count(path)))
            log.append("  {:6s} exported {} triangles".format(preset, counts[-1][1]))
        finally:
            doc.close(False)
    results["fresh_document"] = counts
    log.append("")

    # --- D: explicit numeric deviations, no preset ---------------------------
    log.append("D: explicit surfaceDeviation, no preset (expect Custom on readback)")
    doc, design = import_sphere(app)
    try:
        mgr, units = design.exportManager, design.unitsManager
        counts = []
        for dev_mm in CUSTOM_DEVIATIONS_MM:
            path = os.path.join(OUT_DIR, "custom_{}.stl".format(dev_mm))
            stl = mgr.createSTLExportOptions(design.rootComponent, path)
            stl.isBinaryFormat = True
            try:
                stl.surfaceDeviation = units.convert(dev_mm, "mm", "cm")
            except Exception:
                stl.surfaceDeviation = dev_mm / 10.0
            log.append("  requested {:g} mm -> {}".format(dev_mm, describe(stl, units)))
            mgr.execute(stl)
            counts.append((dev_mm, tri_count(path)))
            log.append("  {:g} mm exported {} triangles".format(dev_mm, counts[-1][1]))
        results["custom"] = counts
    finally:
        doc.close(False)
    log.append("")

    # --- verdict -------------------------------------------------------------
    def varies(key):
        return len(set(c for _, c in results.get(key, [])) ) > 1

    log.append("RESULT")
    log.append("  shared-document presets vary: {}  {}".format(
        varies("shared_document"), results.get("shared_document")))
    log.append("  fresh-document presets vary: {}  {}".format(
        varies("fresh_document"), results.get("fresh_document")))
    log.append("  explicit deviations vary:    {}  {}".format(
        varies("custom"), results.get("custom")))
    log.append("")
    if varies("shared_document"):
        log.append("VERDICT: presets work in reversed order too, so the batch's identical "
                   "meshes were NOT a broken assignment. Suspect first-export or "
                   "cross-document state; compare the deviation readbacks above.")
    elif varies("fresh_document"):
        log.append("VERDICT: presets work ONLY with a fresh document per export. The batch "
                   "reused one document per model, which is the bug. Fix: re-import per "
                   "export.")
    elif varies("custom"):
        log.append("VERDICT: the preset path does not reach the exporter, but explicit "
                   "surfaceDeviation does. Fix: drive the batch with numeric deviations "
                   "instead of presets -- which is better anyway, since the deviation is "
                   "then a stated input rather than a named preset.")
    else:
        log.append("VERDICT: nothing varies the mesh, not even explicit deviations. The "
                   "export is not honouring options at all on this build; stop and "
                   "reconsider before generating more corpus.")

    with open(LOG, "w") as fh:
        fh.write("\n".join(log) + "\n")
    ui.messageBox("\n".join(log[-8:]) + "\n\nFull log: {}".format(LOG))


def run(context):
    try:
        main()
    except Exception:
        app = adsk.core.Application.get()
        if app and app.userInterface:
            app.userInterface.messageBox("BrepperProbe failed:\n{}".format(
                traceback.format_exc()))
