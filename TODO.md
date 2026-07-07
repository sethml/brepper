# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

## 1. stepped_block_toroidal_fillets MEDIUM: one undetected corner-blend torus
The FINE variant now passes the full pipeline (boundary refinement fixed its corner
vertex; all 18 tori + 37 cylinders selected; ~95 s total). MEDIUM still fails: one
corner vertex-blend torus (R=2.5, r=0.5 @ (−66.75,−64.5,−18.2)) is undetected — at
this resolution the blend patch is ~30 faces with almost no interior after junction
masking, so pass-1 seeds don't form and pass-2 (relaxed) seeds are contaminated. A
dedicated small-patch path (seed from the patch's connected component of leftover
faces with outlier trimming) may close it.

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

(completed and removed: curvature-guided seeding stages 1+2 — per-face field with
ascending seed ordering AND pass-1 junction exclusion from gathering, plus the
toroidal-before-conical stage reorder they enabled. See DEVELOPMENT_PLAN
"Curvature-guided seed ordering". Candidate follow-up: the multiface_planar /
large_sphere masks may now be redundant with the junction mask and could be retired
after verification.)

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
