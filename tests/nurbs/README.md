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
- **Two mesh styles, both from the exact geometry.** The doubly-curved patch
  models ship scripted regular-grid STLs with vertices sampled exactly on
  the reference surface (measured ≤2e-6 mm, i.e. f32 STL quantization; see
  `step-measure-tolerance`) — the cleanest possible fitter input. The
  ruled/extruded models (`elliptical_arch`, `elliptical_cylinder`,
  `bspline_wave`) arrive as full-length STRIP meshes with no interior
  vertices along the ruling — how CAD exporters actually tessellate
  extrusions. That is deliberate: the NURBS recognizer has to handle real
  exports, so it gets built against strip meshes from day one instead of
  conveniently gridded ones. The conic models' strips come from genuine
  OCCT tessellation; the wave's are scripted (vertices exactly on the
  curve) because BRepMesh internally grid-refines B-spline-based surfaces
  at any deflection/angle setting. `hill_with_hole` is a plain OCCT
  tessellation — one realistically irregular CAD-export-style mesh.
- **Slab models** (`bspline_*`): the exact surface overhangs a box footprint
  (pole lattice padded beyond it); the solid is Cut(box, Prism(surface face))
  so walls and bottom are box planes and the top is the untouched exact
  surface. The generator validates every written STL by parsing it back:
  closed oriented manifold, and divergence-theorem volume against the B-rep.
  (Neither BRepGProp quadrature variant alone is trustworthy on freeform
  faces: fixed-order Gauss is ~0.5% low on bspline_hill's bump and ~4e-5
  off on B-spline extrusions; the adaptive variant measured 4% high on a
  boolean-trimmed extrusion sheet and 0.3% low on the same surface as a
  natural prism face. The generator accepts agreement with either, plus a
  machine-exact profile integral for the wave.)

## Models (increasing recognizer difficulty)

| Model | Freeform surface | Exercises |
|-------|------------------|-----------|
| elliptical_arch | half-ellipse (15×8) linear extrusion, 20 long | Simplest target: one non-periodic singly-curved patch, sharp 90° boundaries, planar neighbors. Strip mesh: 2 vertex rings along the extrusion. |
| bspline_hill | bicubic patch, Gaussian-shaped 15×15 pole lattice | Canonical doubly-curved dome on a slab; smooth curvature decay into near-flat corners (plane-absorption risk at loose tolerance). |
| bspline_saddle | bicubic patch, hyperbolic-paraboloid pole net | Negative Gaussian curvature everywhere (only torus inner regions have it today). |
| bspline_wave | cubic B-spline profile extruded along Y | Curvature sign alternates; crest/trough radii ≈3 mm are cylinder-detector bait — the over-segmentation trap. Strip mesh across Y. |
| bspline_terrain | bicubic patch: bump + dip + diagonal ripple | Multi-feature single patch; curvature varies in sign and magnitude. |
| bspline_blister | bicubic patch, raised poles only near center | Outer pole rings flat ⇒ margin is EXACTLY planar (B-spline local support), C² tangent join to the bump: segmentation boundary is ambiguous, junction curve is freeform. |
| hill_with_hole | bspline_hill minus off-center vertical cylinder r=5 | Freeform × analytic interaction: inner loop on the B-spline face, cylinder∩B-spline edge curves, irregular tessellated mesh. |
| elliptical_cylinder | full ellipse (12×7) prism, 18 tall | Stretch case: closed periodic freeform wall needs a periodic B-spline face + seam (out of scope for the first §2.7 pass). Strip mesh: full-height wall strips. |

Current failure mode (recorded 2026-07-09): stage 2 leaves each freeform
region as single-face planar hypotheses (39 surfaces for the elliptical_arch
strip mesh, 53 for bspline_wave, 935 for the gridded bspline_hill) and the
stage 3.1 compare hard-errors with reconstructed vertices up to millimeters
from the reference. The strip-mesh models fail fast (~1 s); the dense grid
patches take ~35–105 s (bspline_hill ~104 s, near the 120 s run_test_table
timeout). Both will drop once stage 2.7 collapses those regions into single
surfaces.
