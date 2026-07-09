# NURBS Recognizer Test Models

Targets for the DEVELOPMENT_PLAN §2.7 NURBS fallback. Every model contains at
least one surface outside brepper's analytic set (plane / cylinder / sphere /
cone / torus), so **all of these are expected to FAIL full `--compare` until
§2.7 lands** — they are the goal posts for that work, listed in RESULTS.md as
their own group. Stage 1 (mesh validity and STL↔STEP vertex agreement) passes
today and is covered by `stage1_integration.rs`.

Regenerate with `./scripts/generate_occt_models.sh` (see
`src/bin/gen_test_models.rs`).

## Construction principles

- **Exact references, no fitting.** Every STEP is built from exact geometry:
  explicit-pole `Geom_BSplineSurface` lattices, conic / B-spline extrusions,
  and booleans of those. The recognizer will use OCCT fitting APIs
  (`GeomAPI_PointsToBSplineSurface`), so a fitted reference would be circular
  (shared biases cancel) and only approximate. Here the reference *is* ground
  truth.
- **Meshes sampled on the exact surfaces.** The scripted STLs place vertices
  exactly on the reference geometry (measured ≤2e-6 mm, i.e. f32 STL
  quantization; see `step-measure-tolerance`). Ruled/extruded surfaces are
  meshed on a parameter grid because OCCT tessellation would emit 2-ring
  strip meshes with no interior samples along the ruling. `hill_with_hole`
  instead uses OCCT tessellation deliberately, to provide one realistically
  irregular CAD-export-style mesh.
- **Slab models** (`bspline_*`): the exact surface overhangs a box footprint
  (pole lattice padded beyond it); the solid is Cut(box, Prism(surface face))
  so walls and bottom are box planes and the top is the untouched exact
  surface. The generator cross-checks B-rep volume against the closed
  scripted mesh's divergence-theorem volume. (Neither BRepGProp quadrature
  variant alone is trustworthy on freeform faces: fixed-order Gauss is ~0.5%
  low on bspline_hill's bump; the adaptive variant is ~4% high on the
  extrusion solid. The generator accepts agreement with either.)

## Models (increasing recognizer difficulty)

| Model | Freeform surface | Exercises |
|-------|------------------|-----------|
| elliptical_arch | half-ellipse (15×8) linear extrusion, 20 long | Simplest target: one non-periodic singly-curved patch, sharp 90° boundaries, planar neighbors. |
| bspline_hill | bicubic patch, Gaussian-shaped 15×15 pole lattice | Canonical doubly-curved dome on a slab; smooth curvature decay into near-flat corners (plane-absorption risk at loose tolerance). |
| bspline_saddle | bicubic patch, hyperbolic-paraboloid pole net | Negative Gaussian curvature everywhere (only torus inner regions have it today). |
| bspline_wave | cubic B-spline profile extruded along Y | Curvature sign alternates; crest/trough radii ≈3 mm are cylinder-detector bait — the over-segmentation trap. |
| bspline_terrain | bicubic patch: bump + dip + diagonal ripple | Multi-feature single patch; curvature varies in sign and magnitude. |
| bspline_blister | bicubic patch, raised poles only near center | Outer pole rings flat ⇒ margin is EXACTLY planar (B-spline local support), C² tangent join to the bump: segmentation boundary is ambiguous, junction curve is freeform. |
| hill_with_hole | bspline_hill minus off-center vertical cylinder r=5 | Freeform × analytic interaction: inner loop on the B-spline face, cylinder∩B-spline edge curves, irregular tessellated mesh. |
| elliptical_cylinder | full ellipse (12×7) prism, 18 tall | Stretch case: closed periodic freeform wall needs a periodic B-spline face + seam (out of scope for the first §2.7 pass). |

Current failure mode (recorded 2026-07-09): stage 2 leaves the freeform
region as hundreds of single-face planar hypotheses (e.g. 51 surfaces for
elliptical_arch, 935 for bspline_hill) and the stage 3.1 compare hard-errors
with reconstructed edges millimeters from the reference. Runtimes are ~15 s
(arch) to ~70 s (hill) — near the 120 s run_test_table timeout; they will
drop once stage 2.7 collapses those regions into single surfaces.
