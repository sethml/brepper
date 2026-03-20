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
