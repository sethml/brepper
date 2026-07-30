# Fusion 360 batch re-export (Phase 0.5a)

`BrepperReexport.py` is a Fusion script that imports brepper's reference STEPs
and re-exports them as STL at Low/Medium/High refinement, writing everything to
`tests/fusion_reexport/`.

## Running it

Fusion → **Utilities → Add-Ins → Scripts** → `BrepperReexport` → **Run**.

The script is installed by symlink, so it stays versioned here:

```bash
ln -sf "$PWD/tools/fusion/BrepperReexport.py" "$PWD/tools/fusion/BrepperReexport.manifest" \
  "$HOME/Library/Application Support/Autodesk/Autodesk Fusion 360/API/Scripts/BrepperReexport/"
```

13 models, 31 exports. It shows a progress dialog with a working Cancel, and it
skips any STL already on disk, so an interrupted run resumes where it stopped.
Imported documents are closed without saving — nothing lands in the Fusion hub.

## Why

The corpus was Parasolid strips (Onshape) plus constrained Delaunay (OCCT).
Fusion is a third tessellator over the *same* reference geometry, so any
behaviour difference is attributable to the mesh alone. The model list is not a
matter of taste: the OCCT re-mesh experiment (`tests/remesh/`) already sorted
eight models into ones where a new tessellator changed the outcome and ones
where it did not, so the discriminators get all three presets, geometry OCCT
never covered gets all three, and models that passed at every OCCT density get
Low only — coarse is where every failure except `stepped_block`'s appeared.

The single most valuable row is `stepped_block_toroidal_fillets`: it fails at
all four OCCT densities (88–183 surfaces against a reference 71) yet passes on
the Parasolid fine mesh. If Fusion fails too, the problem is "anything but
strips"; if it passes, it is specific to Delaunay. Those point at different
Phase 1 work.

## Outputs

* `<stem>_{low,medium,high}.stl` — named after the reference's stem.
* `sources.tsv` — the STL-stem → reference-STEP pairing. `tools/sweep.sh` reads
  a `sources.tsv` from a model's own directory in preference to any naming rule,
  because guessing a reference by stem search is ambiguous (two directories hold
  a `plate_with_hole_100x50.step`) and resolves in filesystem order.
* `fusion_export_log.txt` — the Fusion version, plus per-export triangle count
  and the **numeric surface/normal deviation each preset reports**. No export
  dialog states those, and Phase 0.7's sag-estimator cross-check currently has
  four truth values, all inferred from `--compare` output rather than stated as
  inputs. Treat a zero as "not reported", not as "zero deviation" — whether
  Fusion populates these for a preset rather than only for Custom is unverified.

## After the run

Nothing needs checking by hand. Fusion may heal or reparameterize geometry on
STEP import, which would silently invalidate the pair, so each export is run
through stage 1 `--compare` against the **original** reference: anything whose
vertices do not lie on those surfaces within a few quanta is quarantined as
crash-test-only rather than treated as ground truth, and the verification rate
is itself recorded (Phase 0.5c).
