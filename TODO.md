# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

## 1. Fix stepped_block_toroidal_fillets (medium + fine)
The most representative remaining model: a block with a network of interacting toroidal
fillets, i.e. what real CAD parts look like. Two problems:
- **Selection**: 73 surfaces selected vs 72 in STEP (21/17 planar, 36/37 cylindrical,
  16/18 toroidal); stage 3.1 edge compare fails on a sharp edge (~0.36 mm). A few
  fillet junction regions mis-assign.
- **Performance**: the fine variant (51K triangles) takes ~119 s, missing the plan's
  <30 s target for the 10K–100K band. Cost is dominated by repeated failed conical
  BFS trials and strip-merge torus seeds across the fillet network.

## 2. Multi-shell solids (internal voids)
Stage 3.6 turns every closed shell into its own solid, so a part with an internal
cavity produces a wrong STEP (the void surfaces become a second, nested solid instead
of a cavity in the outer one). Manifold input with voids is legal per the project's
assumptions. Plan:
- Add a ccad test model: box minus a fully-interior box (2 shells, 1 solid with void).
- Classify shells by containment (signed volume + point-in-solid test).
- Build the outer solid, then add inner shells as voids (`BRep_Builder::Add` of
  reversed inner shells, or `ShapeFix_Solid` with holes).

## 3. Partial cones + plane-cone tangent lines (cone_tangent_planes_medium)
STEP has 3 cones + 5 planes; only 1 cone is detected. The two side cone patches span a
narrow azimuthal arc and fail the angular coverage acceptance (which requires wide
azimuthal span to reject bogus fits). Real transition pieces have such partial cones.
- Add a partial-cone acceptance path: narrow-span cone patches are acceptable when the
  fit is tight and the patch is bounded by tangent continuations (planes with matching
  slope). A minimal (h,r) profile check (line fit residual vs circle fit residual) may
  serve as the discriminator — implement only what this needs, not a general framework.
- Implement plane-cone tangent line construction in stage 3.3 (documented TODO:
  `Geom_Line` along the cone generator at the tangent azimuth).

## 4. Full-torus seam coverage
No test model exercises a complete donut (needs both U and V seams per the stage 3.4
notes). Add a full-torus model (Onshape export or scripted STL + OCCT-generated
reference STEP) and fix whatever the seam handling needs.

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
  still be added if partial cones (item 3) need it, as a targeted tool rather than a
  framework.
