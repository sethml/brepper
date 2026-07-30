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

13 models, 31 exports, named `<model>_fusion_{coarse,medium,fine}.stl`. It shows
a progress dialog with a working Cancel, skips any STL already on disk so an
interrupted run resumes, closes imported documents without saving so nothing
lands in the Fusion hub, and flushes its log line by line so a crash cannot lose
the record.

**It drives numeric deviations, not the Low/Medium/High presets**, at values
mirroring `tests/remesh`'s OCCT settings (0.5 / 0.1 / 0.02 mm surface deviation
at 30 / 15 / 6° normal deviation). And it **self-verifies**: a curved model whose
triangle counts do not strictly increase across densities is reported as
VERIFICATION FAILED rather than written to disk as three densities, with
all-planar bodies exempted by inspecting face geometry. The header comment
records why — the preset mechanism produced byte-identical meshes in the first
run and three diagnostic runs could not localise the cause, while presets and
every variant of the export block work correctly in isolation.

`BrepperProbe.py` and `BrepperBisect.py` are those diagnostics, kept because they
document what was ruled out: presets work in a shared and a fresh document
(168 / 528 / 1224 triangles on sphere_25), and filename-assignment order, the
failing property read, and both together all export correctly.

The root cause turned out to be neither: **`STLExportOptions` inherits ambient
state.** `unitType` initializes to the Design's default units — metres for an
imported STEP — so every export from every one of these scripts was written 1000×
too small until it was pinned. Re-measuring the probe and bisect outputs confirmed
they were in metres too, so the scale error predated the preset question. Pin every
property; verify the artifact, not a proxy. Validate a fresh corpus with
`tools/verify_pairs.sh` before trusting it.

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
