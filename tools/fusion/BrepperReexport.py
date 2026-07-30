# Batch re-export of brepper's reference STEPs as STL, from inside Fusion.
#
# Phase 0.5(a) of the Robustness Rework: add a THIRD tessellator provenance to a
# corpus that is Parasolid strips (Onshape) plus constrained Delaunay (OCCT).
# The ground truth does not change — every model is a re-tessellation of a
# reference STEP that already lives in tests/ — so any behaviour difference is
# attributable to the mesh alone.
#
# Run it: Fusion → Utilities → Add-Ins → Scripts → BrepperReexport → Run.
# Output goes to REPO/tests/fusion_reexport/. Safe to re-run: it skips STLs
# already on disk, so an interrupted run resumes.
#
# ---------------------------------------------------------------------------
# WHY THIS DRIVES NUMERIC DEVIATIONS AND NOT THE LOW/MEDIUM/HIGH PRESETS
# ---------------------------------------------------------------------------
#
# The first run used presets and produced BYTE-IDENTICAL meshes at Low, Medium
# and High for 8 of 9 models — including sphere_25, which cannot be
# refinement-invariant. Two diagnostic scripts settled what was and was not
# happening:
#
#   * BrepperProbe: presets work in isolation. sphere_25 gives 168 / 528 / 1224
#     triangles for Low / Medium / High, identically in a shared and in a fresh
#     document, with deviation readbacks confirming each assignment landed.
#   * BrepperBisect: every variant of this script's export block also works in
#     isolation — filename assigned after the refinement, the failing property
#     read, both together, filename assigned first — all five giving 168.
#   * Yet the batch produced 9024 triangles for sphere_25 at all three requested
#     presets: finer than Fusion's own High, so not a silent fallback either.
#
# Then a fourth run exposed the actual root cause, which explains both anomalies:
# EVERY export from every one of these scripts was written in METRES — a 25 mm
# sphere spanning 0.025 — because `unitType` "initializes to the default units
# specified for the Design", and an imported STEP's design defaulted to metres.
# `printUtility` likewise defaults to "last UI setting". So an STLExportOptions
# object inherits ambient state, and an export that does not set EVERY property
# explicitly is not reproducible. That is one explanation covering both the
# byte-identical presets and the silent 1000x scale error.
#
# This version therefore (a) sets every property explicitly, unitType included,
# (b) drives numeric deviations rather than presets so nothing can substitute a
# preset, and (c) VERIFIES THE GEOMETRY of each written file against the design's
# own bounding box, not just that the triangle count changed. The earlier
# self-check passed while every file was unusable, because it checked a proxy.
#
# That is an improvement, not a workaround. The deviations become stated inputs
# rather than named presets; they are comparable across models, where Fusion's
# presets scale with bounding-box size so "Low" means different things on a 10 mm
# cube and a 160 mm block; and they are chosen to MIRROR tests/remesh's OCCT
# settings exactly — 0.5 / 0.1 / 0.02 mm at 30 / 15 / 6 degrees. Fusion versus
# OCCT then becomes a controlled comparison at identical requested tolerances
# instead of a comparison of two vendors' preset names.
#
# Both deviations must be set. BrepperProbe measured that sweeping
# surfaceDeviation from 0.5 mm to 0.005 mm changed nothing on a sphere while
# normalDeviation stayed at 15 degrees: on curved primitives the ANGULAR
# tolerance binds, not the chordal one.

import os
import pathlib
import struct
import traceback

import adsk.core
import adsk.fusion

REPO = "/Users/sethml/src/brepper"
OUT_DIR = os.path.join(REPO, "tests", "fusion_reexport")
LOG_PATH = os.path.join(OUT_DIR, "fusion_export_log.txt")

# Skip a model whose STL is already on disk, so an interrupted run resumes.
RESUME = True

# (label, surfaceDeviation mm, normalDeviation DEGREES), identical to the
# deflection/angle pairs tests/remesh uses for OCCT.
#
# DEGREES, and assigned as-is. The API reference says normalDeviation is in
# radians; it is not, and believing it cost a run. BrepperProbe read Fusion's OWN
# presets back as normalDeviation 30 / 15 / 10 for Low / Medium / High, which are
# degrees (30 radians would be five full turns). Converting 30 degrees to 0.5236
# was therefore read as 0.52 DEGREES -- 57x tighter than intended -- and every
# model hit an internal ceiling and produced identical 600k-triangle meshes, a
# 148 MB corpus. Trust the readback over the documentation.
DENSITIES = [
    ("coarse", 0.5, 30.0),
    ("medium", 0.1, 15.0),
    ("fine", 0.02, 6.0),
]
ALL = [d[0] for d in DENSITIES]
COARSE_ONLY = ["coarse"]

# (reference STEP relative to REPO/tests, density labels to export). Chosen from
# the OCCT re-mesh results rather than by taste: that experiment sorted eight
# models into ones where a new tessellator changed the outcome and ones where it
# did not.
MODELS = [
    # --- discriminators: provenance already changed the result here ---
    # Fails at every OCCT density (88-183 surfaces vs a reference 71) where the
    # Parasolid fine mesh reaches 75. The key row in this batch.
    ("onshape/stepped_block_toroidal_fillets", ALL),
    # Fillet network. OCCT coarse gives 34 surfaces vs 26, with 4 spurious tori.
    ("onshape/rounded_cube_10_r2", ALL),
    # OCCT coarse fails; OCCT medium/fine pass with a DIFFERENT decomposition
    # (4 planes + 1 cone + 2 B-splines) than Parasolid's (4 planes + 2 cones +
    # 1 sphere). The first Fusion run agreed with Parasolid.
    ("onshape/cone_tangent_planes", ALL),
    # OCCT coarse yields 6 planes + 2 spheres for a one-face reference.
    ("onshape/sphere_25", ALL),
    # --- geometry the OCCT re-mesh set never covered ---
    ("onshape/chamfered_cube_10_c1", ALL),
    ("onshape/part_rounded_cube_10_r2", ALL),
    # Failed on the first Fusion run at 10 surfaces against a reference 3.
    ("onshape/nosecone", ALL),
    ("onshape/stepped_block", ALL),
    ("onshape/l_bracket_simple", ALL),
    # --- controls: passed at every OCCT density, so only coarse can teach ---
    ("onshape/pipe_elbow_10", COARSE_ONLY),
    ("onshape/pill", COARSE_ONLY),
    ("onshape/cylinder_10x30", COARSE_ONLY),
    # --- calibration against the Fusion export already in tests/fusion/ ---
    ("fusion/plate_with_hole_100x50", COARSE_ONLY),
]


class Log(object):
    """Appends and flushes every line. The first bisect run built its log in
    memory and wrote it at the end, so when the reporting block raised, five
    successful experiments left a zero-byte file behind."""

    def __init__(self, path):
        self.handle = open(path, "w")

    def __call__(self, line=""):
        self.handle.write(line + "\n")
        self.handle.flush()

    def close(self):
        try:
            self.handle.close()
        except Exception:
            pass


def to_internal(units, value_mm):
    try:
        return units.convert(value_mm, "mm", "cm")
    except Exception:
        return value_mm / 10.0


def to_mm(units, value):
    try:
        return units.convert(value, "cm", "mm")
    except Exception:
        return value * 10.0


def readback(stl, units):
    """One try per property. Fusion spells it maximumEdgeLength; the first run
    used maxEdgeLength inside a try shared with two other reads, and that
    AttributeError discarded the values that would have diagnosed this at once."""
    out = []
    try:
        out.append("refinement={}".format(int(stl.meshRefinement)))
    except Exception as exc:
        out.append("refinement=err({})".format(type(exc).__name__))
    for label, getter, is_len in (
        ("surfDev", lambda: stl.surfaceDeviation, True),
        ("normDev", lambda: stl.normalDeviation, False),
        ("maxEdge", lambda: stl.maximumEdgeLength, True),
        ("aspect", lambda: stl.aspectRatio, False),
    ):
        try:
            v = getter()
            out.append("{}={:.6g}{}".format(
                label, to_mm(units, v) if is_len else v, "mm" if is_len else ""))
        except Exception as exc:
            out.append("{}=err({})".format(label, type(exc).__name__))
    return " ".join(out)


def tri_count(path):
    return max(0, (os.path.getsize(path) - 84) // 50) if os.path.exists(path) else -1


def stl_span_mm(path):
    """Largest bounding-box extent of a written binary STL, in whatever units it
    was written. Parsed from the file rather than asked of Fusion, because the
    question is what actually landed on disk."""
    data = pathlib.Path(path).read_bytes()
    if len(data) < 84:
        return -1.0
    n = struct.unpack("<I", data[80:84])[0]
    lo = [1e30] * 3
    hi = [-1e30] * 3
    for t in range(n):
        base = 84 + 50 * t + 12
        for v in range(3):
            for c in range(3):
                off = base + 12 * v + 4 * c
                x = struct.unpack("<f", data[off:off + 4])[0]
                lo[c] = min(lo[c], x)
                hi[c] = max(hi[c], x)
    return max(h - l for l, h in zip(lo, hi))


def design_span_mm(design, units):
    """Largest bounding-box extent of the design itself, in mm."""
    box = design.rootComponent.boundingBox
    lo, hi = box.minPoint, box.maxPoint
    return max(to_mm(units, b - a) for a, b in
               ((lo.x, hi.x), (lo.y, hi.y), (lo.z, hi.z)))


def has_curved_face(design):
    """Whether refinement can legitimately change this model's mesh at all. An
    all-planar body — chamfered_cube at 44 triangles, stepped_block at 52 — is
    correctly density-invariant and must not be reported as a failure."""
    plane = adsk.core.SurfaceTypes.PlaneSurfaceType
    for body in design.rootComponent.bRepBodies:
        for face in body.faces:
            try:
                if face.geometry.surfaceType != plane:
                    return True
            except Exception:
                return True          # unknown geometry: assume it can vary
    return False


def configure_export(mgr, design, out_path, units, surf_mm, norm_deg):
    """Every export option, set explicitly.

    Shared with the canary deliberately: a pre-flight check that configures the
    export differently from the batch proves nothing about the batch.

    Nothing is left to inherit. `unitType` initializes to the Design's default
    units, which for an imported STEP is METRES, and 31 exports were written
    1000x too small before it was pinned. `printUtility` likewise defaults to
    the last UI setting."""
    stl = mgr.createSTLExportOptions(design.rootComponent, out_path)
    stl.unitType = adsk.fusion.DistanceUnits.MillimeterDistanceUnits
    stl.isBinaryFormat = True
    stl.isOneFilePerBody = False
    stl.sendToPrintUtility = False
    # Numeric deviations, never a preset: with no preset involved nothing can
    # silently substitute one. Setting either deviation moves meshRefinement to
    # Custom, which is the intent. normalDeviation is in DEGREES despite the API
    # reference saying radians -- see the DENSITIES comment.
    stl.surfaceDeviation = to_internal(units, surf_mm)
    stl.normalDeviation = norm_deg
    # Pinned so they neither bind nor drift with ambient defaults: the observed
    # ambient maximumEdgeLength was 10 mm, small enough to dominate a 164 mm
    # model and to subdivide flat faces for no benefit.
    stl.maximumEdgeLength = to_internal(units, 1.0e6)
    stl.aspectRatio = 21.5
    return stl


def canary(app, log):
    """Export one curved model at every density into a scratch directory and
    require strictly increasing triangle counts. Returns (ok, reason).

    sphere_25 is the canary because it cannot legitimately be
    density-invariant, and it is small enough to cost seconds."""
    scratch = os.path.join(REPO, "tmp", "fusion_canary")
    os.makedirs(scratch, exist_ok=True)
    step = os.path.join(REPO, "tests", "onshape", "sphere_25.step")
    doc = None
    try:
        options = app.importManager.createSTEPImportOptions(step)
        doc = app.importManager.importToNewDocument(options)
        if doc is None:
            doc = app.activeDocument
        design = adsk.fusion.Design.cast(app.activeProduct)
        mgr, units = design.exportManager, design.unitsManager
        counts = []
        for label, surf_mm, norm_deg in DENSITIES:
            path = os.path.join(scratch, "sphere_25_{}.stl".format(label))
            stl = configure_export(mgr, design, path, units, surf_mm, norm_deg)
            log("canary {:6s}: {}".format(label, readback(stl, units)))
            mgr.execute(stl)
            counts.append(tri_count(path))
            span = stl_span_mm(path)
            expected = design_span_mm(design, units)
            ratio = span / expected if expected > 0 else 0.0
            log("canary {:6s}: {} triangles, span {:.4g} mm (ratio {:.4g})".format(
                label, counts[-1], span, ratio))
            if not 0.97 < ratio < 1.03:
                return False, ("wrong scale: sphere_25 {} spans {:.4g} mm against a "
                               "design span of {:.4g} mm. unitType is not being "
                               "honoured.".format(label, span, expected))
        if not all(b > a for a, b in zip(counts, counts[1:])):
            return False, ("densities do not vary: sphere_25 gives {} triangles for "
                           "coarse/medium/fine, and a sphere cannot be "
                           "density-invariant. The deviation request is not reaching "
                           "the tessellator.".format(counts))
        log("canary OK: {} triangles across coarse/medium/fine".format(counts))
        log()
        return True, ""
    except Exception:
        return False, "canary raised:\n{}".format(traceback.format_exc())
    finally:
        if doc is not None:
            try:
                doc.close(False)
            except Exception:
                pass


def main():
    app = adsk.core.Application.get()
    ui = app.userInterface
    os.makedirs(OUT_DIR, exist_ok=True)
    log = Log(LOG_PATH)

    by_label = dict((d[0], d) for d in DENSITIES)
    total = sum(len(labels) for _, labels in MODELS)
    log("brepper Fusion re-export")
    log("Fusion version: {}".format(app.version))
    log("densities (surfaceDeviation mm / normalDeviation deg): {}".format(
        ", ".join("{}={:g}/{:g}".format(*d) for d in DENSITIES)))
    log("{} models, {} exports".format(len(MODELS), total))
    log()

    # CANARY. Five runs of this script have now produced unusable corpora --
    # metres instead of millimetres, degrees mistaken for radians, presets that
    # silently did nothing -- and each cost a full batch to discover. Prove the
    # parameterization varies the mesh on ONE model first, and abort in seconds
    # rather than after 31 exports and 148 MB.
    ok, why = canary(app, log)
    if not ok:
        log("ABORTED: {}".format(why))
        log.close()
        ui.messageBox("brepper re-export ABORTED before writing the corpus.\n\n{}\n\n"
                      "Nothing in tests/ was touched.\n\nLog: {}".format(why, LOG_PATH))
        return

    progress = ui.createProgressDialog()
    progress.isCancelButtonShown = True
    progress.show("brepper re-export", "%v of %m: %p%", 0, total, 0)

    done = failures = skipped = 0
    unverified = []
    scale_failures = []
    sources = []
    try:
        for rel_step, labels in MODELS:
            stem = rel_step.split("/")[-1]
            step_path = os.path.join(REPO, "tests", rel_step + ".step")
            for label in labels:
                sources.append("{}_fusion_{}\ttests/{}.step".format(stem, label, rel_step))

            wanted = []
            for label in labels:
                out_path = os.path.join(OUT_DIR, "{}_fusion_{}.stl".format(stem, label))
                if RESUME and os.path.exists(out_path):
                    skipped += 1
                    done += 1
                else:
                    wanted.append((label, out_path))
            if not wanted:
                log("{}: already exported, skipped".format(stem))
                progress.progressValue = done
                continue
            if not os.path.exists(step_path):
                log("{}: MISSING reference {}".format(stem, step_path))
                failures += len(wanted)
                done += len(wanted)
                progress.progressValue = done
                continue

            doc = None
            try:
                options = app.importManager.createSTEPImportOptions(step_path)
                doc = app.importManager.importToNewDocument(options)
                if doc is None:
                    doc = app.activeDocument
                design = adsk.fusion.Design.cast(app.activeProduct)
                if design is None:
                    raise RuntimeError("imported document is not a Design")
                mgr, units = design.exportManager, design.unitsManager
                curved = has_curved_face(design)

                counts = {}
                for label, out_path in wanted:
                    if progress.wasCancelled:
                        raise KeyboardInterrupt()
                    _, surf_mm, norm_deg = by_label[label]
                    stl = configure_export(
                        mgr, design, out_path, units, surf_mm, norm_deg)
                    log("{}_{}: requested surfDev={:g}mm normDev={:g}deg -> {}".format(
                        stem, label, surf_mm, norm_deg, readback(stl, units)))
                    mgr.execute(stl)
                    counts[label] = tri_count(out_path)
                    # Verify the GEOMETRY, not a proxy. The previous run's
                    # self-check confirmed triangle counts increased while every
                    # file was written in metres and therefore useless.
                    written = stl_span_mm(out_path)
                    expected = design_span_mm(design, units)
                    ratio = written / expected if expected > 0 else 0.0
                    scale_ok = 0.97 < ratio < 1.03
                    log("{}_{}: {} triangles, span {:.4g} mm vs design {:.4g} mm "
                        "(ratio {:.4g}){}".format(
                            stem, label, counts[label], written, expected, ratio,
                            "" if scale_ok else "  *** WRONG SCALE ***"))
                    if not scale_ok:
                        scale_failures.append("{}_{}".format(stem, label))
                    done += 1
                    progress.progressValue = done

                # Self-verification. Distinct densities are the entire point of
                # the batch, so a curved model whose counts do not strictly
                # increase means the request never reached the tessellator — and
                # that has to be loud, not silently written to disk.
                ordered = [counts[l] for l in ALL if l in counts]
                increasing = all(b > a for a, b in zip(ordered, ordered[1:]))
                if not curved:
                    log("{}: all-planar body, density-invariant by construction "
                        "({} triangles) — expected".format(stem, ordered[0]))
                elif len(ordered) > 1 and not increasing:
                    unverified.append(stem)
                    log("{}: *** VERIFICATION FAILED *** curved body but counts are {}. "
                        "The deviation request did not reach the tessellator; do NOT "
                        "treat these as distinct densities.".format(stem, ordered))
                log()
            except KeyboardInterrupt:
                log("CANCELLED after {} of {} exports".format(done, total))
                break
            except Exception:
                failures += len(wanted)
                done += len(wanted)
                progress.progressValue = done
                log("{}: FAILED\n{}".format(stem, traceback.format_exc()))
            finally:
                if doc is not None:
                    try:
                        doc.close(False)
                    except Exception:
                        pass
    finally:
        progress.hide()

    with open(os.path.join(OUT_DIR, "sources.tsv"), "w") as fh:
        fh.write("# STL stem\treference STEP it was tessellated from\n")
        fh.write("\n".join(sources) + "\n")

    log("exports written: {}, skipped: {}, failed: {}".format(
        done - skipped - failures, skipped, failures))
    summary = "written {}, skipped {}, failed {}".format(
        done - skipped - failures, skipped, failures)
    if scale_failures:
        log("WRONG SCALE on: {}".format(", ".join(scale_failures)))
        summary += "\n\nWRONG SCALE ({} files) — the export did not honour "\
                   "unitType. Do not use these files.".format(len(scale_failures))
    if unverified:
        log("DENSITY VERIFICATION FAILED on: {}".format(", ".join(unverified)))
        summary += "\n\nDENSITY VERIFICATION FAILED on:\n  {}".format(
            "\n  ".join(unverified))
    if not scale_failures and not unverified:
        log("all curved models produced strictly increasing densities at the "
            "correct scale")
        summary += "\n\nAll models verified: correct scale, and strictly "\
                   "increasing densities on every curved body."
    log.close()
    ui.messageBox("brepper re-export finished.\n\n{}\n\nLog: {}".format(summary, LOG_PATH))


def run(context):
    try:
        main()
    except Exception:
        tb = traceback.format_exc()
        try:
            with open(LOG_PATH, "a") as fh:
                fh.write("\nFATAL\n{}\n".format(tb))
        except Exception:
            pass
        app = adsk.core.Application.get()
        if app and app.userInterface:
            app.userInterface.messageBox("BrepperReexport failed:\n{}".format(tb))
