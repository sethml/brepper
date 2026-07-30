# Batch re-export of brepper's reference STEPs as STL, from inside Fusion.
#
# Phase 0.5(a) of the Robustness Rework: add a THIRD tessellator provenance to a
# corpus that is Parasolid strips (Onshape) plus constrained Delaunay (OCCT).
# The ground truth does not change — every model here is a re-tessellation of a
# reference STEP that already lives in tests/ — so any behaviour difference is
# attributable to the mesh alone.
#
# Run it: Fusion → Utilities → Add-Ins → Scripts → BrepperReexport → Run.
# Everything is written under REPO/tests/fusion_reexport/. Safe to re-run: it
# skips STLs that already exist, so an interrupted run resumes.
#
# It also records what no GUI export dialog reports: the numeric surface and
# normal deviation each refinement preset actually uses. tools/sweep.sh's
# Phase 0.7 estimator cross-check has four sag truth values, all inferred from
# --compare output; these are stated inputs.

import os
import traceback

import adsk.core
import adsk.fusion

REPO = "/Users/sethml/src/brepper"
OUT_DIR = os.path.join(REPO, "tests", "fusion_reexport")

# Skip a model whose STL is already on disk, so an interrupted run resumes
# instead of redoing an hour of work.
RESUME = True

# (reference STEP relative to REPO/tests, refinement presets to export)
#
# Chosen from the OCCT re-mesh results rather than by taste: that experiment
# sorted eight models into ones where a new tessellator changed the outcome and
# ones where it did not. The four discriminators get all three presets, the
# geometry OCCT never covered gets all three, and the models that passed at
# every OCCT density get "low" only — coarse is where every failure but
# stepped_block's appeared.
MODELS = [
    # --- discriminators: provenance already changed the result here ---
    # Fails at ALL FOUR OCCT densities (88-183 surfaces vs a reference 71) but
    # passes on the Parasolid fine mesh. The key experiment in this batch: if
    # Fusion also fails, the problem is "anything but strips"; if it passes, it
    # is specific to Delaunay. Those imply different Phase 1 work.
    ("onshape/stepped_block_toroidal_fillets", ["low", "medium", "high"]),
    # Fillet network. OCCT coarse gives 34 surfaces vs 26 with 4 spurious tori,
    # where the Parasolid coarse mesh passes.
    ("onshape/rounded_cube_10_r2", ["low", "medium", "high"]),
    # OCCT coarse fails; OCCT medium/fine pass with a DIFFERENT decomposition
    # (4 planes + 1 cone + 2 B-splines) than Parasolid's (4 planes + 2 cones +
    # 1 sphere), both inside the 1% volume gate. A third vote decides whether
    # that split is provenance-driven.
    ("onshape/cone_tangent_planes", ["low", "medium", "high"]),
    # OCCT coarse yields 6 planes + 2 spheres for a one-face reference.
    ("onshape/sphere_25", ["low", "medium", "high"]),
    # --- geometry the OCCT re-mesh set never covered ---
    ("onshape/chamfered_cube_10_c1", ["low", "medium", "high"]),
    ("onshape/part_rounded_cube_10_r2", ["low", "medium", "high"]),
    ("onshape/nosecone", ["low", "medium", "high"]),
    ("onshape/stepped_block", ["low", "medium", "high"]),
    ("onshape/l_bracket_simple", ["low", "medium", "high"]),
    # --- controls: passed at every OCCT density, so only coarse can teach ---
    ("onshape/pipe_elbow_10", ["low"]),
    ("onshape/pill", ["low"]),
    ("onshape/cylinder_10x30", ["low"]),
    # --- calibration: tests/fusion/ already holds a Fusion export of this at
    # low/medium/high from an earlier session. Re-exporting at low says whether
    # this Fusion version and these settings match what is already committed,
    # so a version difference cannot be mistaken for a pipeline finding.
    ("fusion/plate_with_hole_100x50", ["low"]),
]


def refinement_enum(name):
    settings = adsk.fusion.MeshRefinementSettings
    return {
        "low": settings.MeshRefinementLow,
        "medium": settings.MeshRefinementMedium,
        "high": settings.MeshRefinementHigh,
    }[name]


def to_mm(units_manager, value):
    """Fusion's API speaks internal units (cm) for lengths; the log needs mm."""
    try:
        return units_manager.convert(value, "cm", "mm")
    except Exception:
        return value * 10.0


def main():
    app = adsk.core.Application.get()
    ui = app.userInterface

    os.makedirs(OUT_DIR, exist_ok=True)
    log_path = os.path.join(OUT_DIR, "fusion_export_log.txt")
    log_lines = []
    # Pairing manifest, the same mechanism tests/remesh/sources.tsv uses: the
    # reference lives in another directory, and having the harness guess it by
    # searching for a matching stem is ambiguous and filesystem-order dependent.
    sources = []

    def log(line):
        log_lines.append(line)

    log("brepper Fusion re-export")
    log("Fusion version: {}".format(app.version))
    total = sum(len(presets) for _, presets in MODELS)
    log("{} models, {} exports".format(len(MODELS), total))
    log("")

    progress = ui.createProgressDialog()
    progress.cancelButtonText = "Cancel"
    progress.isBackgroundTranslucent = False
    progress.isCancelButtonShown = True
    progress.show("brepper re-export", "%v of %m: %p%", 0, total, 0)

    done = 0
    failures = 0
    skipped = 0
    for rel_step, presets in MODELS:
        stem = rel_step.split("/")[-1]
        step_path = os.path.join(REPO, "tests", rel_step + ".step")

        wanted = []
        for preset in presets:
            out_path = os.path.join(OUT_DIR, "{}_{}.stl".format(stem, preset))
            sources.append("{}_{}\ttests/{}.step".format(stem, preset, rel_step))
            if RESUME and os.path.exists(out_path):
                skipped += 1
                done += 1
            else:
                wanted.append((preset, out_path))
        if not wanted:
            log("{}: all presets already exported, skipped".format(stem))
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
            # Import into its own document, so nothing touches the user's open
            # designs and the document can be discarded unsaved afterwards.
            doc = app.importManager.importToNewDocument(options)
            if doc is None:
                doc = app.activeDocument

            design = adsk.fusion.Design.cast(app.activeProduct)
            if design is None:
                raise RuntimeError("imported document is not a Design")
            export_mgr = design.exportManager
            units = design.unitsManager

            for preset, out_path in wanted:
                if progress.wasCancelled:
                    raise KeyboardInterrupt()
                stl = export_mgr.createSTLExportOptions(design.rootComponent, out_path)
                stl.meshRefinement = refinement_enum(preset)
                stl.isBinaryFormat = True
                stl.filename = out_path
                # Read the preset's own deviation values back. Whether Fusion
                # populates these for a preset (as opposed to only for Custom)
                # is unverified, so the log records whatever it reports rather
                # than asserting a meaning: a zero means "not reported", not
                # "zero deviation".
                try:
                    surf = to_mm(units, stl.surfaceDeviation)
                    norm = stl.normalDeviation
                    edge = to_mm(units, stl.maxEdgeLength)
                    deviations = (
                        "surfaceDeviation={:.6g} mm normalDeviation={:.6g} "
                        "maxEdgeLength={:.6g} mm".format(surf, norm, edge)
                    )
                except Exception as exc:
                    deviations = "deviations unavailable ({})".format(exc)

                export_mgr.execute(stl)
                size = os.path.getsize(out_path) if os.path.exists(out_path) else 0
                tris = max(0, (size - 84) // 50)
                log("{}_{}: {} triangles, {} bytes, {}".format(
                    stem, preset, tris, size, deviations))
                done += 1
                progress.progressValue = done
        except KeyboardInterrupt:
            log("CANCELLED by user after {} of {} exports".format(done, total))
            break
        except Exception:
            failures += len(wanted)
            done += len(wanted)
            progress.progressValue = done
            log("{}: FAILED\n{}".format(stem, traceback.format_exc()))
        finally:
            # Close without saving: these documents are throwaway imports.
            if doc is not None:
                try:
                    doc.close(False)
                except Exception:
                    pass

    progress.hide()

    log("")
    log("exports written: {}, skipped (already present): {}, failed: {}".format(
        done - skipped - failures, skipped, failures))
    with open(log_path, "w") as fh:
        fh.write("\n".join(log_lines) + "\n")
    with open(os.path.join(OUT_DIR, "sources.tsv"), "w") as fh:
        fh.write("# STL stem\treference STEP it was tessellated from\n")
        fh.write("\n".join(sources) + "\n")

    ui.messageBox(
        "brepper re-export finished.\n\n"
        "written: {}\nskipped: {}\nfailed: {}\n\nLog: {}".format(
            done - skipped - failures, skipped, failures, log_path))


def run(context):
    try:
        main()
    except Exception:
        app = adsk.core.Application.get()
        if app and app.userInterface:
            app.userInterface.messageBox("BrepperReexport failed:\n{}".format(
                traceback.format_exc()))
