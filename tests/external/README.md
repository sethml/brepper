# tests/external — wild STEP+STL pairs

Harvested from GitHub by `tools/harvest_github.py`. Model files are gitignored;
`manifest.tsv` records url + sha256, `sources.tsv` maps each STL to its reference.

## The corpus now: 32 pairs (pruned 2026-11-11)

53 pairs were harvested. 21 were removed after measuring each mesh against its
own STEP, and are listed in `rejected.tsv` with the reason and the url + sha256
that recover them; `tools/harvest_github.py` skips those stems. Every pair that
remains is within tolerance.

| Verdict | Pairs | Meaning |
|---|---|---|
| exact | 28 | every vertex within 5e-5 mm of the trimmed reference |
| trim | 4 | every vertex beyond it misses only the reference's trim (see below) |

The stage-1.2 `--compare` gate passes **28 of 32**. The other four are exact
meshes that stage 1 rejects as non-manifold before it compares
(`OpenBot__100deg_joint`, `…100deg_joint_mirror`,
`VoronUsers__1515_snapin_ziptie_ldo`, `…ziptie_makerbeam`).

## Which CAD package exported each model (read 2026-11-11)

Evidence is the STEP header (`FILE_NAME` preprocessor / originating system) and
the 80-byte binary STL header. The STEP and the STL of a pair were not
necessarily written by the same tool or in the same session, so where the STL
header carries no signature the STL exporter is an inference from the STEP.

| Repo (pairs kept) | STEP header | STL header | Exporter |
|---|---|---|---|
| NucDeck volume buttons (2) | `SwSTEP 2.0` / `SolidWorks 2021` | `solid <part name>` | SolidWorks 2021 (Parasolid) |
| ToolChanger, all but one (8) | `SwSTEP 2.0` / `SolidWorks 2021` | `solid <part name>` | SolidWorks 2021 (Parasolid) |
| ToolChanger `hemera_bowden_bracket` | Autodesk Translation Framework v9.7.1.1290 | `STLB ATF 9.7.1.1290` | Fusion 360 |
| Open-Task-Light (8) | Autodesk Translation Framework v11.7.0.108 / v11.17.0.187 | `STLB ATF` same versions | Fusion 360 |
| OpenBot joints (3) | `SIEMENS PLM Software NX 9.0` | blank | Siemens NX 9 (Parasolid); STL writer unconfirmed |
| OpenBot `front_light_spacer` | Autodesk Translation Framework v10.10.0.1391 (cloud translation path) | `STLB ATF 10.10.0.1391` | Fusion 360 |
| VoronUsers strip mounts, snap-ins (5) | Autodesk Translation Framework v9.3.0.1241 | `STLB ATF 9.3.0.1241` | Fusion 360 |
| VoronUsers `10x11_befenybay_chain_top_anchor` | `CATIA Version 5-6 Release 2018` | no signature | CATIA V5 STEP; STL writer unknown |
| VoronUsers `2020-wago_221-412_tray-5x2` | originating system blank, ST-DEVELOPER v18.1 | no signature | unknown |
| White-Knight (2) | Autodesk Translation Framework v8.6.0.1094 | `STLB ATF 8.6.0.1094` (the original exports; see below) | Fusion 360 |

Of the removed pairs, the NucDeck and ToolChanger ones are SolidWorks 2021 and
the Open-Task-Light, VoronUsers chain-anchor and White-Knight ones are Fusion
360 STEPs (`rejected.tsv` names each).

`tests/downloaded/`: `3DBenchy` STEP is a FreeCAD / Open CASCADE export (STL
header `solid Shape0`); `OrcaCube v2` is Fusion 360 (ATF v12.4.0.73, both
files); `xyz-cube-test-20mm.stl` is Fusion 360 (`STLB ATF 8.12.0.6`).

So the wild corpus is two kernels: Parasolid (SolidWorks, NX) and Autodesk
ShapeManager (Fusion 360), plus one CATIA part. Nothing here is Onshape,
FreeCAD or Rhino.

## Alignment (re-done 2026-11-11)

About half the pairs are published in two frames: the STL laid out for
printing, the STEP left in design coordinates. `gen-test-models align
tests/external [<stem substring>]` re-frames the reference, writes
`<stem>_aligned.step`, re-points `sources.tsv`, and records every pair in
`alignment.tsv` (p50 / p95 / max over every vertex, as published and after the
fit). How it works is documented on `align_corpus` in
`src/bin/gen_test_models.rs`; in short:

1. A coarse frame chosen by score from derived candidates: identity, the
   axis-aligned rotations whose box extents match, and the principal-axes frames.
2. A rigid least-squares refinement (point-to-surface, six degrees of freedom).
3. Multi-body parts are fitted body by body, paired by volume and then by fit.
4. A pair that already verifies as published is left untouched.

The first version stopped at bounding-box centres and the first extents-matching
rotation. That left five genuine pairs 6e-5..1e-2 mm off (a mesh's bounding box
is inset by its chord error), six White-Knight pairs 10-17 mm off (a box has the
same extents after a 180 degree flip), and the four two-body chain anchors
200-300 mm off (no single motion exists). Several of those pairs have since
been removed for other reasons; `alignment.tsv` covers the 32 that remain.

## The boundary rule ("trim")

A bounded point-to-face distance measures "off the surface" and "outside the
trim" at once. On `ToolChanger__docking_station`, `…docking_station_(m)`,
`…blower_fan_duct` and `…blower_fan_duct_(m)` every over-tolerance vertex (12 to
60 per part, up to 4.3e-4 mm) is a boundary vertex on both of its surfaces to
~1e-6 mm, sitting beside an edge that OCCT itself reads as good to only that
much. Stage 1.2 excuses exactly these (`reference::BoundaryMiss::excused`);
`gen-test-models strays <stl> <step>` prints the per-vertex classification, and
with `BREPPER_STRAYS_DUMP=<path>` writes every vertex's distance instead.

## White-Knight: the STLs were re-saved after export

Upstream commit "Updated STL file orientation" (2019-10-22) rewrote every STL
two days after the STEPs were committed. Against those STLs the fitted rotations
come out at 89.9995 and 120.0461 degrees and 18-97% of vertices sit
1e-4..5e-3 mm off across the whole part.

The original exports (commit `4f0178691b`, "Added STL files.", 2019-10-20) are
direct Fusion 360 exports in the STEP's own frame: no alignment is needed.
Measured against their STEPs:

| Part | Original export | Outcome |
|---|---|---|
| `1-5i_case_front_v4` | exact, max 3.7e-6 mm | kept, as the original export (manifest url pinned to the commit) |
| `1-belt_stepper_roller_endcap_v1` | exact, max 7.8e-6 mm | kept, as the original export (pinned) |
| `1-5i_case_rear_v7` | median 4e-7 mm; 320 of 2213 vertices beyond tolerance, max 6.0e-3 mm | removed |
| `1-berdair_mount_clamp_v3` | median 5e-6 mm; 986 of 3016 beyond, max 4.9e-3 mm | removed |
| `1-direct_drive_hotend_carriage_v10` | median 2e-6 mm; 1757 of 6085 beyond, max 5.3e-3 mm | removed |
| `1-right_rear_idler_part_2_…_v3` | median 1e-5 mm; 1213 of 9888 beyond, max 5.5e-3 mm | removed |

In the four that are not exact, the vertices beyond tolerance are almost all in
face interiors, up to about 5e-3 mm — the same pattern as the removed
`…mount2_for_e3d_v2`, also a direct Fusion export from this author. Re-saving
is therefore not what makes them inexact, and the cause is not established.

## Removed pairs (21, in `rejected.tsv`)

- **STEP and STL are different geometry (7).** Vertices more than 0.01 mm
  off, confined to one feature, while the rest of the part is exact.
  - `Open-Task-Light__belt_clip`, `…carriage_belt_clip`: 530 of 1394 vertices
    (38%), every one exactly 0.200 mm off; volume +11%.
  - `NucDeck__joystick_surround_ring`: 145 vertices (8%) 0.06-3.3 mm off;
    volume +27%. Upstream, the STL was committed six weeks after the STEP.
  - `NucDeck__right_side_grip`: 72 vertices 0.44-0.85 mm off in one
    35 x 11 mm planar patch.
  - `White-Knight__1-direct_drive_extruder_mount_v9`: 615 vertices (8%)
    0.01-0.24 mm off in a 10 x 13 x 13 mm region, in the original export too.
  - `NucDeck__back_cover…` (6 vertices, up to 0.35 mm, one 1 mm spot) and
    `White-Knight__1-5i_case_front_knight_v14` (30 vertices, up to 0.10 mm,
    one spot, in the original export too): a small local difference.
- **Genuine pairs that are not within tolerance (14).** No vertex more than
  0.01 mm off, but 0.6-33% of vertices beyond 5e-5 mm after the boundary rule:
  the NucDeck trigger buttons and left side grip, `ToolChanger__dock_brace` and
  `…dock_brace_(m)`, the four two-body VoronUsers chain anchors,
  `White-Knight__…mount2_for_e3d_v2`, and the four White-Knight parts in the
  table above. In the first eight the misses are a mix
  of vertices on one surface but up to 4e-3 mm off its neighbour and vertices
  in face interiors. The largest miss across all 14 is 6.0e-3 mm, far below
  the 0.1 mm scale of a resized feature. Which side is off was NOT established:
  the measurement is mesh against the STEP as OCCT reads it, so the mesh may
  be exact and the STEP (or its translation) approximate. If such meshes do
  turn out to trouble the converter, `--vertex-tolerance` can be raised.
