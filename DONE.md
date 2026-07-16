# Foundation

## Port surface fit unit tests from C++
Added comprehensive unit tests for all four surface fitting modules (cylindrical, conical, spherical, toroidal). Each module now has distance tests and fit tests using synthetic meshes. 53 tests pass.

## Project setup
Crates, cargo.toml, dependencies.

## Test utility
Read an STL and a STEP, compute maximum distance between STL vertices and STEP surfaces, and print it out. Created a script in scripts/ to apply it to all of the stl/step file pairs under tests/ and print out a table of maximum distances.

## Program skeleton
Implemented main program, flag parsing, and created stub source files for each stage. Created stage output data structures. Took existing stage 1.1 and stage 1.2 implementations in mesh.rs and reformulated them into stage1.rs and appropriate output data structures.

# Stage 1 Mesh Input

## Stage 1.1: Read STL file
Read STL file into `ConnectedMesh`, including welded vertices and per-face placeholder fields for neighbors, normals, and hypotheses.

## Stage 1.2: Mesh validation
Mesh validation pass to compute face normals, edge neighbors, manifold stats, connected shells, and orientation consistency checks.

## Stage 1.3: Coplanar triangle fusion
Merge adjacent coplanar triangle pairs into quads. This simplifies surface fitting by ensuring curved-surface facets (cylinder/sphere quad strips) become single faces rather than 2-face planar hypothesis groups.

## Test STL/STEP file pairs
Implemented tests that all stl/step file pairs in tests/ pass consistency checks and pass when fed to brepper with the --compare flag. Also invented file pairs in tests/bad that fail with the --compare flag and bad cases that fail the mesh validation tests in various ways.

# Stage 2: Surface Fitting

## Stage 2.1: Deduce planar hypotheses
Understood stage 2.1 and existing test models. Created additional challenging test shapes in tests/ccad/. Implemented stage 2.1 with --compare passing for all test shapes composed only of planar surfaces.

## Stage 2.2: Deduce cylindrical hypotheses
Understood stage 2.2, created challenging test shapes, and implemented cylindrical hypothesis deduction.

## Stage 2.3: Deduce spherical hypotheses
Understood stage 2.3, created challenging test shapes, and implemented spherical hypothesis deduction.

## Stage 2.6: Select surfaces for reconstruction
Implemented surface selection using per-face priority rule, then revisited and replaced with greedy area-based selection.

## --angular-tolerance flag
Implemented the --angular-tolerance flag for cylindrical and spherical surface hypothesis generation.

## Revisit stage 2.2: Cylinder improvements
Updated cylindrical hypothesis to match updated algorithm description from commits 2ac6d4f onward. Tested that deduced cylinders on test models match the step cylinders.

## Revisit stage 2.3: Sphere improvements
Updated spherical hypothesis to match updated algorithm description (vertex-based seeding, solid-angle coverage validation, surface-tolerance in BFS, relaxed max_sphere_radius). Fixed the sphere overgrowth problem blocking `onshape_rounded_cube_stage26`.

# Stage 3: Surface Reconstruction

## Stage 3.1: Create OCCT surface objects and build adjacency graph
Created OCCT surface objects and built adjacency graph from mesh connectivity. Implemented --compare check.

## Stage 3.2: Detect tangency relationships
Detect tangency relationships between adjacent surfaces. Computes outward-facing surface normals at sampled boundary points and compares angle; marks edge as tangent if normals agree within 2°.

## Stage 3.3: Compute edge curves
Computed edge curves via surface-surface intersection (`GeomAPI_IntSS`), trim to vertex endpoints via `GeomAPI_ProjectPointOnCurve` + `Geom_TrimmedCurve`. Multi-curve selection picks closest to mesh boundary vertices.

## Stage 3.4: Create OCCT faces
Created OCCT faces from surfaces bounded by edge wires. Handles planar faces with wire-based construction, cylindrical/spherical faces with UV-bounds or wire-based construction, conical and toroidal faces. Validates all faces with `BRepCheck_Analyzer`.

## Stage 3.5: Construct shells
Stitched faces into shells using `BRepBuilderAPI_Sewing`. Applies `ShapeFix_Shell` orientation fixing.

## Stage 3.6: Construct solids
Constructed solids from shells using `ShapeFix_Solid::SolidFromShell`. Validates with `BRepCheck_Analyzer`, computes volume via `BRepGProp::VolumeProperties`.

# Stage 4: Output

## Stage 4.1: Write STEP file
Write solids to STEP file using `STEPControl_Writer` with `StepModelType::Asis`. Builds compound from all solids, sets AP214 schema. Compare re-reads written file and validates volume agreement and `BRepExtrema_DistShapeShape` distance. 15 integration tests.

# Stage 3 Refinements

## Revisit stage 3.2: Tangency validation
Tangency detection correctly identifies 8 tangent edges on part_rounded_cube_10_r2 (cylinder-plane tangencies at fillet boundaries). No changes needed to the 2° threshold.

## Revisit stage 3.3: Tangent edge curves
Implemented tangent edge curve computation. For tangent edges, constructs curves analytically rather than using `GeomAPI_IntSS`. Plane-cylinder tangencies produce lines parallel to the cylinder axis. Fixed `MakeEdge` with adaptive vertex tolerance.

## Sphere-cylinder tangent edge curves
Great circle arcs on the sphere surface. Adaptive per-vertex tolerance in MakeEdge accommodates surface fit imprecision. Unblocked rounded_cube_10_r2_coarse full pipeline.

## Closed-loop sphere-cylinder tangent edges + pole-safe hemispheres
Closed-loop tangent circles constructed as full `Geom_Circle` curves. Sphere surface reoriented when boundary circles pass within 45° of a pole. Unblocked pill_coarse.stl and pill_fine.stl.

## Evaluated stage 2.7 (surface refitting)
Not needed. All test models pass --compare with max face-to-STEP distances of ~3e-6 mm. Un-ignored 7 tests that were previously marked ignored.

# Stage 2 Extensions: Conical Surfaces

## Create ccad test models for cones
Created `simple_cone.lua`, `block_with_conical_hole.lua`, `cone_cylinder.lua`, `nosecone.lua` with STL+STEP pairs.

## ConicalHypothesis data structure
Added `ConicalHypothesis` struct with fields for apex, axis_direction, half_angle, convex, faces, vertices, errors. Extended `MeshFace` and `Stage2Output`.

## Cone fitting algorithm
Axis estimation from normal covariance, (h,r) profile fitting via linear regression. Apex-vertex seeding strategy. Levenberg-Marquardt refinement for 6-parameter optimization. False-positive rejection filters.

## Cone BFS region growing
Vertex-to-cone distance validation and angular coverage check. Multi-axis LM rescue for seeds. BFS expansion with refit-if-close pattern.

## Cone surface selection
Extended stage 2.6 surface selection to include conical candidates in greedy area-based selection.

## Cone/cylinder/sphere/plane discrimination
Half-angle bounds, apex distance sanity check, normal-axis angle standard deviation check. All 325 tests pass with zero false positives.

## Cone unit and integration tests
5 unit tests and 57 integration tests covering hypothesis counts, parameter matching, compare tests, topology tests, and output tests. Total: 382 tests.

## Cone full pipeline tests
All cone models pass `--compare` through stage 4.1. Three bug fixes: angular tolerance relaxed, apex UV bounds extended, concave flag added.

## Unblock cone_15x20_medium (onshape)
Stage 2 detection works perfectly (58 faces, ha=20.56°, err_max=2.4e-6). Full pipeline passes through stage 4.1.

# Stage 2 Extensions: Toroidal Surfaces

## Create ccad test models for tori
Created `filleted_cylinder.lua`, `filleted_hole_block.lua`, `filleted_pipe.lua` with STL+STEP pairs.

## ToroidalHypothesis data structure
Added `ToroidalHypothesis` struct. Plumbed through stages 1, 2, and 3.

## Torus fitting via medial axis / tube center method
Estimate minor radius from normal-line intersections, compute tube centers, fit 3D circle for major circle parameters. Vertex-neighborhood seeding with post-seeding merge. 7 integration tests.

## Torus BFS region growing
Vertex-to-torus distance validation.

## Torus surface selection
Extended stage 2.6 to include toroidal candidates. Added 10 integration tests.

# Stage 3 Extensions for Cone and Torus

## Stage 3.1: Torus OCCT surface creation
Added `Geom_ToroidalSurface` creation from hypothesis parameters.

## Stage 3.2: Torus tangency normals
Added torus analytical normal formulas for tangency detection.

## Stage 3.3: Torus tangent edge curves
Implemented tangent edge curve computation for plane-torus, cylinder-torus, and torus-torus pairs.

## Stage 3.3: Cone tangent edges — not needed
At plane-cone and cylinder-cone junctions, surface normals always differ by the cone half-angle, far exceeding the 2° tangency threshold.

## Stage 3.3: Non-tangent cone and torus edge intersections
Already works for all test models via `GeomAPI_IntSS` with curve selection.

## Stage 3.4: Conical face construction
Already works for all test models. Existing UV-bounds and face creation code handles cones correctly.

## Stage 3.4: Toroidal face construction
Fixed doubly-periodic V-bounds using circular gap algorithm. Fixed filleted_cylinder volume: 13.5% error → <1e-7 relative diff.

# Incorporate Learnings from C++ Implementation

## Port abstracted BFS region growing
Created `src/stage2/bfs.rs` with a `SurfaceType` trait and generic `generate_hypotheses::<S>()` function. Implemented for all 5 surface types. **Regressions**: 8 `block_with_conical_hole` tests fail (missing multi-axis LM rescue for cone seeds).

## Compare C++ vs Rust surface fitting
Evaluated all surface types. Recommendations: keep Rust for all, port sphere Gauss-Newton refinement, analytic Jacobians for cone/torus, torus minor-circle acceptance test, and vertex-position PCA for torus axis.

# Ports from C++ Implementation

## Port sphere Gauss-Newton refinement
After algebraic least-squares sphere fit, runs iterative Gauss-Newton to minimize vertex distance residuals. Up to 20 iterations, convergence at delta_norm < 1e-14, requires ≥6 vertices.

## Port analytic Jacobians for cone LM
Replaced finite-difference Jacobian in ConeLMProblem with closed-form partial derivatives for all 6 parameters (apex position, axis tilt, half-angle).

## Port analytic Jacobians for torus LM
Replaced finite-difference Jacobian in TorusLMProblem with closed-form partial derivatives for all 7 parameters (center, axis tilt, major/minor radii).

## Port torus minor-circle acceptance test
Added angular coverage validation around the torus minor circle (tube cross-section). Computes φ = atan2(axial, radial-R) for each face centroid and requires the angular span ≥ 0.1 rad. Prevents false positives from strips that don't wrap around the tube.

## Port vertex-position PCA for torus axis estimation
Replaced normal covariance with vertex-position PCA in fit_torus_sor(). The smallest eigenvector of the position covariance matrix gives the torus axis — works because axial variance ~ r²/2 while perpendicular variance ~ (R² + r²/2)/2. Also added Kasa circle fit for center location.

## Delete legacy deduce_*_hypotheses() functions
Removed dead code: deduce_planar_hypotheses, deduce_cylindrical_hypotheses, deduce_spherical_hypotheses, deduce_conical_hypotheses, deduce_toroidal_hypotheses. These were superseded by the generic BFS framework.

## Interactive viz features
Added to `--viz-stages` window: initial zoom with 15% margin, help overlay (`?`), face number labels (`N`) with depth-buffer occlusion testing, hover face info (index/area/normal via ray-cast picking), and click-to-skip-seeding (left-click sends `SkipToSeed(face)` to BFS). Integrated egui for all text overlays. Restructured render loop so egui runs after 3D rendering to enable depth-buffer readback.

## Viz bug fixes: occlusion, hover highlight, click-to-skip deadlock
- Replaced unreliable depth-buffer-based face label occlusion with face-normal culling: labels for back-facing faces (normal pointing away from camera) are hidden.
- Added thick bright-yellow wireframe highlight on the hovered face; hovered face's number is drawn larger with a white halo for readability.
- Fixed click-to-skip deadlock: click only sends `SkipToSeed` when `waiting_for_input = true` (BFS is blocked in `recv`), and immediately sets `waiting_for_input = false` so the render loop calls `try_recv()` again and can receive the next overlay — without this the sync channel would deadlock.

## Delete legacy per-type trial BFS functions
Removed 1,305 lines of dead code left over from the pre-generic BFS implementation: `run_cone_trial_bfs()` and `ConeTrialResult` from conical.rs, `run_torus_trial_bfs()` and `TorusTrialResult` from toroidal.rs, `all_vertices_within_*_tolerance()` helpers from conical/toroidal/spherical/planar, `build_vertex_to_faces_map()` and `SphericalSurface::vertex_to_faces()` from spherical.rs, `vertex_to_plane_distance()` from planar.rs, and `viz_custom()` from mod.rs. Cleaned up all `#[allow(dead_code)]` annotations in stage2 and removed unused imports (`VecDeque`, `VizAction`, `VizSender`, `viz_custom`, `viz_face_centroid`, `viz_face_normal`, `REFIT_SKIP_MULTIPLIER`, `build_normal_covariance`) from conical.rs and toroidal.rs.

## Robust cone fitting for pointed cones + prune unselected hypotheses
Rewrote cone fitting algorithm to handle both pointed and truncated cones. For pointed cones (all faces share a single apex vertex), uses apex-vertex seeding with fixed-apex LM refinement. For truncated cones, uses multi-strategy axis estimation from normal differences. Added acceptance test with normal-axis variance threshold and half-angle bounds check to prevent false positive cones on spherical/cylindrical patches. Implemented post-selection pruning of unselected hypotheses with index remapping in `Stage2Output`.

## Move test table to tests/RESULTS.md
Moved the test cases table from DEVELOPMENT_PLAN.md into a new script-generated tests/RESULTS.md. Rewrote scripts/run_test_table.sh to generate the markdown file directly with header, per-directory tables, and summary. The script auto-discovers all STL files under tests/ so new models are picked up automatically. DEVELOPMENT_PLAN.md now references tests/RESULTS.md instead of embedding the table.

## Conical hypothesis visualization
Added cone overlay rendering to the interactive viz. Conical hypotheses (half_angle > 0) now display as translucent cone surfaces instead of cylinders. Added `ConeOverlay` struct, `generate_cone_mesh` function, and `centered_cone_overlay` helper. Updated BFS viz, stage 2.6 competing/winning hypothesis viz, and the render loop to handle cone overlays.

## Toroidal hypothesis visualization
Added torus overlay rendering to the interactive viz. Toroidal hypotheses now display as translucent torus surfaces instead of cylinders at the major radius. Added `TorusOverlay` struct, `generate_torus_mesh` function, and `build_torus_overlay_mesh`. Updated BFS trait, all `viz_overlays` implementations, stage 2.6 competing/winning hypothesis viz, and the render loop.

## Rework cone/cylinder fitting parameters and LM refinement
Overhauled cone and cylinder fitting for robustness: set `REFIT_SKIP_MULTIPLIER` to 10.0 for aggressive refitting, set `seed_vtol_multiplier` to 1.0 (vertices lie precisely on the surface), removed normal-axis consistency check for cones (surface selection resolves overlaps), added vertex recentering around centroid before LM for numerical stability, and increased LM patience to 500 for both cone and cylinder fits. Updated SURFACE_FITTING.md and DEVELOPMENT_PLAN.md to reflect changes.

## Overhaul stage 2 debug output and add --viz-stages-level

Changed default `--vertex-tolerance` from 1e-5 to 5e-5 to better match typical Onshape STL export precision. 

Overhauled BFS debug printing: removed all `[prefix]` tags (e.g. `[spherical-BFS]`), switched to consistent indentation levels (stage summary at 0 indent, per-hypothesis at 2, seed attempts at 4+, per-face BFS at 6+). Split verbosity into 5 levels: 1=stage summaries, 2=per-hypothesis details, 3=per-committed-hypothesis progress/stats, 4=per-seed-attempt details, 5=per-BFS-face exploration.

Added `--viz-stages-level` flag (default 1) to control how often `--viz-stages` pauses: level 1 pauses at committed hypotheses and stage 2.6 selections, level 2 also pauses at accepted seeds, level 3 also pauses at each individual BFS face acceptance. The `viz_level` field is now propagated through `BfsConfig`.

# Session: Test speedup + cone/torus completion (2026-07-06)

## Test suite speedup (2:56 → ~1:15, all 438+ tests green)
- `[profile.dev.package.brepper] opt-level = 3` so tests run optimized crate code.
- Parallelized all --compare DistShapeShape loops across threads (stage 1.2 vertices/centroids, stage 2 hypothesis/selection compares).
- Spherical detection (2.3) runs before conical (2.2); conical seeding skips sphere-committed faces. Kills the O(mesh) bogus cone/cylinder band tiling on fine spheres (pill_fine stage 2.2: 45 s → 0.002 s).
- Cylinder-first fitting with cone-LM skip on exact cylinder fits; LM patience 500 → 150.

## Fixed cone/torus regressions (17 failing tests at session start)
- Restored normal-axis consistency acceptance check for conical hypotheses (plus a per-face variant during BFS expansion).
- Restored angular_tol_multiplier 2.0 for cones (fixes block_with_conical_hole full pipeline); added it for tori (twisted torus quads don't fuse, leaving triangles with doubled normal steps).
- Reworked brittle hypothesis-count tests to assert on stage 2.6 *selected* surfaces.

## Full pipeline for fillet/cone models
- filleted_cylinder, filleted_pipe, filleted_hole_block, block_with_conical_hole, nosecone (ccad + onshape medium), cone_truncated_medium, chamfered_cube ascii now pass --compare through stage 4.1 (new stage-4.1 tests added).
- Torus minor-circle coverage measured from face normals (rejects flat annuli / cylinder bands lying on degenerate torus families).
- Torus face V-bounds: centroids only disambiguate the window, never extend it.
- Post-selection absorption pass merges spurious sub-fit surfaces and tangent-junction bridge facets into larger adjacent surfaces.
- Sphere-cone tangent circles constructed from boundary-vertex averaging perpendicular to the cone axis; reoriented sphere faces compute V bounds from projected boundary circles instead of fixed hemispheres.
- Deterministic fitting: all fit functions take sorted vertex slices instead of HashSets.
- validate_faces wraps OCCT detailed diagnostics in catch_unwind (Standard_NullObject on invalid faces no longer aborts).

## Unblock pipe_elbow (strip-merge torus seeding)
pipe_elbow_10 medium and fine pass --compare through stage 4.1. Strip-merge torus seeds (3-ring neighborhoods spanning ≥2 cylinder hypotheses) recover the elbow torus from the cylinder strips of stage 2.2; cylinder-torus tangent circles try both candidate axes; toroidal UV-bounds classify iso-U vs iso-V boundaries. New tests: elbow torus params match STEP (R=10, r=5), stage-4.1 compare for pipe_elbow_10_medium.

## Multi-shell solids (internal voids)
Stage 1.2 accepts contained negative-volume shells as voids (voids_within_solids stat); stage 3.6 merge_void_solids merges void shells into their container solid with orientation verified by volume. hollow_cube (tests/occt/, generated via OCCT since ccad segfaults) passes stage 4.1 with one 7000 mm³ solid.

## Partial cones + plane-cone tangent generators
cone_tangent_planes: 3 of 3 cones detected (normal-line axis rescue with profile-linearity gate; half-angle cap 80°→85°); plane-cone tangent edges are generator lines through the endpoint vertices; 11/11 edge curves. Stage 3.4 partial-cone face construction issue remains (documented in TODO.md).

## Full-torus seam coverage
full_torus test model (complete donut R=10, r=3) passes the full pipeline: complete tori use natural-bounds face construction (both seams created by OCCT automatically).

## Curvature-guided seeding (stages 1+2)
Per-face curvature field drives ascending seed ordering (flattest regions seed first), and pass-1 junction faces are excluded from gathering. Together these enabled the toroidal-before-conical stage reorder. See DEVELOPMENT_PLAN "Curvature-guided seed ordering". Candidate follow-up (in TODO.md): the multiface_planar / large_sphere masks may now be redundant with the junction mask.

## Detect medium-resolution fillets (stepped_block_toroidal MEDIUM)
All 18 tori + 37 cylinders now selected; the model passes stage 4.1 (dedicated onshape_stepped_block_toroidal_medium_stage41_compare test). Three fixes: (a) the BFS final re-fit keeps the incrementally-maintained fit when a fresh from-scratch LM converges worse — recovered the full-ring R=3.5 pocket-rim fillet torus that refit to 0.87 mm; (b) the cylinder angular-coverage test admits coarse tessellation gaps up to the dihedral limit (medium fillets step ~24° per row, four rows across a 90° tube arc, which the span/3 rule rejected); (c) the default --angular-tolerance was raised from 17.5° to 25° so seed-pair qualification accepts those 24° fillet rows. The coverage allowance is gated twice so it can't admit bogus bands: to cylinders (half_angle≈0), so narrow corner sphere-blend patches can't win selection as 6-face cones (else rounded_cube regresses); and to hypotheses already spanning ≥45°, so thin cylinder bands of a smooth cone stay rejected (else nosecone_truncated over-segments). See LEARNINGS "Medium-resolution tessellation".

## Stage-3 shared-edge topological assembly (DEVELOPMENT_PLAN §3.x rework)
Stage 3.4/3.5 rebuilt around shared edges: every ReconEdge is one TopoDS_Edge used by
both adjacent faces, with per-face pcurves computed analytically (SurfaceChart) and
attached via BRep_Builder; boundary loops are branch-aligned in UV and their winding
classifies wrapped faces; seam edges (iso-curves anchored at shared vertices) and
degenerate pole edges (cone apex, sphere pole) are first-class topology built by one
uniform mechanism; sphere axes point at the region's mean direction so boundaries never
graze UV poles; wires/faces/shells are BRep_Builder-assembled (no MakeWire/MakeFace-
from-wire/ShapeFix in face construction) and shells are closed by construction, so
BRepBuilderAPI_Sewing is gone. Stage 3.3 gained mesh-boundary-chain validation/repair
of edge curves (the chain is ground truth for any junction). Fixes cone_tangent_planes
(apex-cap disks; volume 60% off → 1e-7 relative) and improves stepped_block_toroidal
volume (2.7e-3 → 1.6e-3); ReconEdge.pcurves and BRepVertex.uv_coords are now populated.
Full suite green; 52/57 models pass (only bad/ fail, by design).

## Stage 2.7 NURBS fallback (DEVELOPMENT_PLAN §2.7)
Freeform regions — connected components of leftover single-face planar selections,
grown across smooth edges and absorbing the small accidental analytic fits freeform
surfaces bait (coplanar strip pairs at inflections, cylinder ribbons, bump-peak
cones) — are Floater-mean-value parameterized onto the unit square, resampled at
mesh-vertex-aligned grid stations, fitted with GeomAPI_PointsToBSplineSurface
(interpolation first; OCCT's approximator overshoots its Tol3D near curved patch
boundaries), and validated per-vertex against --surface-tolerance with
own-UV-seeded projections. Stage 3 reconstructs them through a B-spline
SurfaceChart wrapping the exact fitted surface (ShapeAnalysis_Surface with a
local-minimum gap re-check); the seam/winding/pole machinery self-disables.
Planes tangent to a region along their entire shared boundary merge in (blister's
exactly-planar margin); B-splines are never absorption/refinement targets (the
fit doesn't vouch beyond its region). The boundary-chain edge repair gained its
B-spline third arm (line/circle/B-spline, best by chain deviation) — first
demanded by the arch's elliptical junctions. Runs before absorption so accidental
analytic fits can't hide from gathering. All 8 tests/nurbs models pass full
--compare (volumes 2.7e-6…5.2e-3 relative to the exact references) — closed
freeform bands (elliptical_cylinder) are seam-cut, fitted once, and split into
two faces sharing the surface, with the two tagged artificial edges exempted
from the stage 3.1 compare (the periodic-face-with-seam form is the documented
follow-up). 60/65 zoo models pass (5 bad/ excepted); full suite green.
