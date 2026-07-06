# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

## 1. Fix stepped_block_toroidal_fillets (medium + fine)
The most representative remaining model: a block with a network of interacting toroidal
fillets. PROGRESS: 17 of its 18 tori are now detected (was 16; two-pass toroidal
seeding, seed outlier trimming, corner-blend multi-axis rescue, final-refit
validation), and the fine variant's stage 1-2 runtime dropped ~119 s → ~43 s.
REMAINING: one corner vertex-blend torus (R=2.5, r=0.5 @ (−66.75,−64.5,−18.2)) still
fails — its seeds are contaminated by leftover faces of ~4 undetected corner WALL
CYLINDERS (narrow 90° arcs at the step corners; same partial-arc fitting problem as
the cone_tangent_planes side patches — the normal-line axis rescue now exists for
cones and could be extended to cylinders). Fixing the wall cylinders should
de-contaminate the torus seeds.

## 2. cone_tangent_planes: stage 3.4 partial-cone face construction
DONE at stage 2: all 3 cones are now detected (normal-line axis rescue for narrow
arcs; half-angle cap raised to 85°), and plane-cone tangent generator lines work in
stage 3.3 (11/11 edge curves). REMAINING: the wire-based partial-revolution conical
face construction produces invalid faces (BRepCheck: face shape invalid / wire
invalid) and the solid volume comes out 6617 mm³ vs 14676 mm³ in STEP — one or more
cone faces cover the wrong UV region (complementary arc or wrong pcurves on narrow
cone patches with vertices).

(items completed and removed: multi-shell voids — hollow_cube passes stage 4.1 with
one 7000 mm³ solid; full-torus seam coverage — full_torus (tests/occt/, R=10 r=3)
passes stage 4.1 via natural-bounds face construction; partial-cone detection and
plane-cone tangent generators — all 3 cone_tangent_planes cones detected, 11/11 edge
curves, stage 3.4 face construction issue remains, see item 2.)

# Deferred

## NURBS fallback (the remaining project-defining feature)
Deferred by explicit decision for now. When picked up: this is the gap between "passes
the test zoo" and the README's stated goal — any region that isn't
plane/cylinder/sphere/cone/torus currently degenerates into single-facet planar
surfaces and breaks reconstruction. Plan sketch:
- Test models: ccad cannot produce freeform surfaces (polyline profiles only). Instead,
  script-generate STL height-field meshes from a known smooth function, and build the
  matching reference STEP via OCCT `GeomAPI_PointsToBSplineSurface` from the same grid,
  giving real `--compare` coverage without any external CAD tool.
- Fit B-splines to connected components of leftover single-face planar regions
  (`GeomAPI_PointsToBSplineSurface`), validate against `--surface-tolerance`.
- Integrate into stage 2.6 selection and stage 3 reconstruction (edge curves via IntSS
  with B-spline surfaces, face creation from UV bounds/wires).

## nosecone_truncated_medium
An angled truncation cutting through a tangent sphere–cone junction. Bridging facets
along the cut produce spurious planes and fake corners. Unusual construction; low
value relative to effort.

## Other
- Support for OBJ, PLY input formats.
- IGES export option.

# Dropped (motivating problems solved by other means)

- **Gaussian map pre-classification**: was a performance idea from when stage 2.2
  wasted ~45 s tiling spheres with bogus cylinder bands. Fixed instead by running
  spherical detection first, per-face normal-axis consistency, and early rejection.
  No current model's runtime justifies it.
- **Unified SOR framework refactor**: its concrete motivations — pipe_elbow torus
  recovery and cylinder/cone/sphere disambiguation — were fixed by strip-merge torus
  seeding, stage ordering, and the absorption pass. A minimal (h,r) profile check may
  still be added where needed, as a targeted tool rather than a framework (the
  normal-line axis rescue for narrow cone arcs and the (h,r) profile fits in the
  torus multi-axis rescue are exactly such tools).
