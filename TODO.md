# Prioritized by value toward the project goal (README: "analytic and freeform surfaces")

## 1. cone_tangent_planes: stage 3.4 conical apex-cap face construction
DONE at stage 2: all 3 cones detected (normal-line axis rescue for narrow arcs;
half-angle cap 85°); stage 3.3 computes 11/11 edge curves incl. plane-cone tangent
generator lines. The model completes and writes a STEP; the only symptom is a stage-3.6
volume WARNING (6617 mm³ vs 14676 mm³) — the surfaces are all correct (max face sample
distance 1.16e-6 mm), only the trimming of two cone faces is wrong.

Sharpened diagnosis (2026-07-07): faces 1 (cone ha 84.25°) and 3 (cone ha 73.04°) are
"pie-slice-with-apex" patches — the cone APEX is one of their boundary vertices, i.e. a
UV U-pole. Both of each face's two boundary edges lie at the SAME iso-V line (same
distance from the apex): e.g. face 1's edge 0 (junction with cone 0) and edge 4
(junction with plane 2, a plane perpendicular to the cone axis tangent to the cone)
both sit at V=11.35, U∈[4.76,7.34]. The generic wire builder (create_wire_based_periodic_
face, is_sphere=false) then collapses the wire: both edges project onto the same iso-V
arc, so the enclosed area is ~0 (face 3 came out 0.0001 mm²; face 1 mis-built to 130 mm²).

Two dead ends already explored:
- Pure UV-bounds "apex cap" (U=[0,2π], V=[0,V_arc], apex as degenerate V-pole, like the
  spherical pole path): gets the AREA right (403 / 163 mm²) but breaks sewing — the cap's
  boundary is OCCT-generated, not the shared edge objects, and it only geometrically
  matches the neighbours along the 148° wedge, leaving the rest free → open shell →
  stage 3.5 hard error. Strictly worse than the current warning-only state.
- Detecting the cap via mesh U-span is unreliable: the apex is a U-pole, so apex-adjacent
  mesh vertices scatter U across the full circle and a 148° wedge falsely reads as a 354°
  full revolution.

The real fix needs a wire-based cap that (a) REUSES the shared edges (edge 0/edge 4 for
face 1) as the outer boundary — unwrapped to the correct arc so they don't stack on the
same U-range — (b) inserts a degenerate edge at the apex U-pole (the sphere pole-insertion
machinery in create_wire_based_periodic_face, currently gated behind is_sphere, generalized
to the cone apex), and (c) picks the correct V-side. Caution: existing code notes ShapeFix
segfaults on cones, so pole/degenerate handling for cones must be built without ShapeFix.
Low value (one unusual model, already emits valid-ish output) — but this is now a concrete,
well-scoped plan rather than a mystery.

(completed and removed: stepped_block_toroidal_fillets MEDIUM — all 18 tori + 37
cylinders now selected and the model passes stage 4.1 (dedicated
onshape_stepped_block_toroidal_medium_stage41_compare test). Three fixes: (a) the BFS
final re-fit now keeps the incrementally-maintained fit when a fresh from-scratch LM
converges worse — recovered the full-ring R=3.5 pocket-rim fillet torus that refit to
0.87 mm; (b) the cylinder angular-coverage test admits coarse tessellation gaps up to
the dihedral limit (medium fillets step ~24° per row, four rows across a 90° tube arc,
which the span/3 rule rejected); (c) the default --angular-tolerance was raised from
17.5° to 25° so seed-pair qualification accepts those 24° fillet rows. The coverage
allowance is gated twice so it can't admit bogus bands: to cylinders (half_angle≈0), so
narrow corner sphere-blend patches can't win selection as 6-face cones (else
rounded_cube regresses); and to hypotheses already spanning ≥45°, so thin cylinder
bands of a smooth cone stay rejected (else nosecone_truncated over-segments). See
LEARNINGS "Medium-resolution tessellation".)

(items completed and removed: multi-shell voids — hollow_cube passes stage 4.1 with
one 7000 mm³ solid; full-torus seam coverage — full_torus (tests/occt/, R=10 r=3)
passes stage 4.1 via natural-bounds face construction; partial-cone detection and
plane-cone tangent generators — all 3 cone_tangent_planes cones detected, 11/11 edge
curves, stage 3.4 face construction issue remains, see item 1.)

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
