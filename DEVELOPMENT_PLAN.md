# STL to STEP Conversion Utility - Development Plan

## Project Overview

**brepper** (B-Rep from Mesh) - A utility to convert triangulated STL meshes from CAD exports into parametric STEP files with fitted analytic and freeform surfaces.

## Assumptions

- Assumes the vertices of the mesh lie precisely on the CAD surfaces, modulo numerical error/stability.
- Assumes transitions between shapes always have a seam.
- Assumes manifold mesh (each edge shared by exactly 2 faces).
- Assumes consistent face orientation (all normals point outward).

Fails noisily if any of the assumptions are violated. Possible future work: offer a repair mode.

## Dependencies

| Library | Purpose | Version |
|---------|---------|---------|
| OpenCASCADE (OCCT) | STL import, B-Rep modeling, surface fitting, STEP export | latest |

Use the opencascade binding in ../opencascade-rs/crates/opencascade-sys, including porting instructions in PORTING.md and details on the FFI mapping in ../opencascade-rs/crates/opencascade-binding-generator/README.md. It's a reasonably complete binding of the C++ API of the opencascade library which resides in ../opencascade-rs/crates/occt-sys/OCCT, with source code in src/ and documentation in dox/.

### Required OCCT Functionality

This section enumerates the OCCT classes and functions needed, for evaluating Rust binding coverage.

**Stage 1: Mesh Input**
- `RWStl` / `RWStl_Reader` — Read binary/ASCII STL files
- `gp_Pnt` — 3D point representation
- `gp_Dir` / `gp_Vec` — Direction and vector types

**Stage 2: Surface Fitting**
- `gp_Pln` — Plane geometry (normal + distance representation)
- `gp_Cylinder` / `gp_Ax3` — Cylindrical surface geometry (axis + radius)
- `gp_Sphere` — Spherical surface geometry (center + radius)
- `gp_Cone` — Conical surface geometry (axis + half-angle)
- `gp_Torus` — Toroidal surface geometry (for fillets)

**Stage 3: Surface Reconstruction**

*Geometry (infinite surfaces and curves):*
- `Geom_Plane`, `Geom_CylindricalSurface`, `Geom_SphericalSurface`, `Geom_ConicalSurface`, `Geom_ToroidalSurface` — Analytic surface objects
- `Geom_BSplineSurface` — NURBS/B-spline surfaces
- `Geom_Line`, `Geom_Circle`, `Geom_Ellipse` — Analytic curve types
- `Geom_BSplineCurve` — Freeform curves
- `Geom2d_Curve`, `Geom2d_Line`, `Geom2d_Circle`, `Geom2d_BSplineCurve` — 2D curves for pcurves

*Surface fitting:*
- `GeomAPI_PointsToBSplineSurface` — Fit B-spline surface to points

*Surface-surface intersection:*
- `GeomAPI_IntSS` — Compute intersection curves between two surfaces

*Point-curve projection:*
- `GeomAPI_ProjectPointOnCurve` — Project point onto curve (for trimming)
- `ShapeAnalysis_Curve` — Curve analysis utilities

*Topology (B-Rep construction):*
- `TopoDS_Vertex`, `TopoDS_Edge`, `TopoDS_Wire`, `TopoDS_Face`, `TopoDS_Shell`, `TopoDS_Solid`, `TopoDS_Compound` — Topological entities
- `BRep_Builder` — Low-level topology construction
- `BRepBuilderAPI_MakeVertex` — Create vertices
- `BRepBuilderAPI_MakeEdge` — Create edges from curves (with optional pcurves)
- `BRepBuilderAPI_MakeWire` — Assemble edges into wires
- `BRepBuilderAPI_MakeFace` — Create faces from surface + bounding wires
- `BRepBuilderAPI_Sewing` — Stitch faces into shells
- `BRepBuilderAPI_MakeSolid` — Create solids from shells

*Topology exploration:*
- `TopExp_Explorer` — Iterate over sub-shapes
- `TopExp::MapShapes` — Build maps of topology
- `BRep_Tool` — Extract geometry from topology (IsSeam, SameParameter, etc.)

*Shape healing and fixing:*
- `ShapeFix_Shape` — General shape repair
- `ShapeFix_Face` — Face-specific repairs
- `ShapeFix_Wire` — Wire repairs, including `FixEdgeCurves` for pcurve computation
- `ShapeFix_Shell` — Shell orientation fixing
- `ShapeFix_Solid` — Solid creation/repair, `SolidFromShell`

*Shape analysis and validation:*
- `BRepCheck_Analyzer` — Validate B-Rep topology/geometry
- `ShapeAnalysis_Shell` — Shell analysis (e.g., `CheckOrientedShells`)
- `GProp_GProps` / `BRepGProp` — Compute volume, area, center of mass

**Stage 4: Output**
- `STEPControl_Writer` — Write STEP files
- `STEPControl_StepModelType` — Control STEP output mode (`AsIs`, `ManifoldSolidBrep`)
- `Interface_Static` — Set STEP header metadata

*Optional/debugging:*
- `BRepTools::Write` — Write BREP format (OCCT native, for debugging)
- `IGESControl_Writer` — Write IGES format (future enhancement)

## Open Questions

- Links between data structures (face => edge, edge => vertex, etc): Use indices or pointers or references or other?
- Data structure management: Vectors? Reference counting?

*Comment: Indices into vectors are generally the right choice for this domain—they're stable under serialization, easily debuggable, and avoid lifetime complexity. Pointers/references become problematic when hypotheses are deleted or vectors reallocate. Consider using "generation counters" or a slot-map pattern if you need to detect stale indices.*

## Architecture

### Stage 1: Mesh Input & Preprocessing

#### 1.1 Read STL File
Read the input STL file and generate an in-memory representation of the triangle mesh, with fields for future stages to traverse the mesh and fit shapes to sets of mesh faces. Weld vertices by position (with tolerance) to build connectivity.
- **Library**: OCCT `RWStl::ReadFile` (static function, via `rw_stl::read_file_charptr_progressrange_2` in the Rust bindings)
- **Input**: Binary or ASCII STL file
- **Output**: `ConnectedMesh`, storing:
    - `vertices: Vec<MeshVertex>` — Mesh vertices with double-precision 3D coordinates (`x`, `y`, `z: f64`).
    - `faces: Vec<MeshFace>` — Mesh faces with:
        - `vertex_count: u8` — 3 for triangles, 4 for quads. Always 3 after initial STL input; stage 1.3 may merge coplanar triangle pairs into quads.
        - `vertex_indices: [usize; 4]` — Indices of vertices, ordered by right-hand rule.
        - `neighbors: [i32; 4]` — Index of the mesh face across each edge, or -1 if none. Filled in stage 1.2.
        - `normal: Option<[f64; 3]>` — Mesh face normal, computed from vertices in stage 1.2. `None` for degenerate faces. (Changed from `gp_Dir` in original plan — native Rust arrays avoid OCCT binding complexity and keep the mesh data structure self-contained.)
        - `planar_hypothesis: i32` — Index of active planar hypothesis, or `NO_HYPOTHESIS` (-1) if none, or `UNDEDUCED_PLANAR_HYPOTHESIS` (-2) if not yet deduced.
        - `cylindrical_hypothesis: i32` — Index of active cylindrical hypothesis, or `NO_HYPOTHESIS` (-1) if none, or `UNDEDUCED_CYLINDRICAL_HYPOTHESIS` (-2) if not yet deduced.
        - `spherical_hypothesis: i32` — Index of active spherical hypothesis, or `NO_HYPOTHESIS` (-1) if none.
    - `stats: MeshValidationStats` — Statistics described in stage 1.2.

Note: hypotheses vectors (planar, cylindrical, spherical) are stored in `Stage2Output`, not in `ConnectedMesh`. The mesh stores only per-face hypothesis indices. (Changed from original plan which listed hypothesis vectors as part of `ConnectedMesh` — separating them into the stage 2 output keeps each stage's output self-contained and avoids coupling the mesh representation to later stages.)

For now, each face can belong to a single hypothesis of each type. Hopefully that's sufficient since hypotheses should nearly exactly match the vertices, but it's possible that in the future we may need to keep a set of candidate hypothesis indices per type, then select the best one in stage 2.6.

**Vertex welding** uses a spatial hash grid: each vertex is bucketed by integer coordinates `(round(x/tol), round(y/tol), round(z/tol))`. To find matches, adjacent cells (3×3×3 neighborhood) are checked for distance ≤ tolerance. This gives O(1) average lookup. The weld tolerance is derived from `--vertex-tolerance` as `min(1e-9, vertex_tolerance * 0.01)`, so welding is much tighter than the surface-fitting tolerance.

#### 1.2 Mesh Validation
Traverse the faces of the mesh, collecting statistics, validating the geometry, and populating normals and neighbors.
- Collect stats into `MeshValidationStats` and optionally print: `mesh_faces`, `mesh_vertices`, `mesh_edges_open` (edges with only 1 incident face), `mesh_edges_non_manifold` (edges with >2 incident faces), `mesh_edges_inconsistent_orientation`, `mesh_faces_degenerate`, `connected_shells`, `solids`, `voids_within_solids`.
- Compute and populate mesh face normals using Newell's method (handles both triangles and quads with consistent winding).
- Compute mesh face neighbors based on shared mesh edges. For each edge (vertex pair), collect all face uses. For manifold edges (exactly 2 uses), set both faces as neighbors. Also check orientation consistency: properly oriented adjacent faces should traverse the shared edge in opposite directions.
- Count connected shells via BFS over a face adjacency graph. A shell is counted as a "solid" if all its edges are manifold (exactly 2 faces per edge). (Simplified from the original plan which suggested ray casting or signed volume analysis — the current heuristic is sufficient for well-formed meshes. Voids within solids are tracked in the stats struct but not yet computed.)
- With the `--compare` flag: check that all vertices are within `--vertex-tolerance` of the STEP shape, and that all face centroids are within `--surface-tolerance`. Uses `BRepExtrema_DistShapeShape` for bounded distance (respects face trimming). Reports the worst vertex and worst centroid with their distances.

Validation checks, in priority order (first failure is reported):
1. **Degenerate faces**: faces with <3 vertices or zero-area normals.
2. **Non-manifold edges**: edges shared by >2 faces.
3. **Inconsistent orientation**: adjacent faces that traverse a shared edge in the same direction (indicating flipped normals within a shell).
4. **Self-intersections**: intentionally deferred to a later stage where richer geometric predicates are available and O(N²) triangle checks can be avoided.

#### 1.3 Coplanar Triangle Fusion (Triangle → Quad Merging)
Merge adjacent coplanar triangle pairs that share a diagonal edge into quads. This is important because CAD tessellators (including CodeCAD/OCCT and Onshape) typically produce quad-strip patterns on curved surfaces (cylinders, spheres, cones). Each quad is split into two coplanar triangles sharing a diagonal. Without fusion, stage 2.1 groups these into 2-face planar hypotheses, which complicates the "single-face = curved surface candidate" criterion in stages 2.2 and 2.3.

**Algorithm (implemented):**

The algorithm operates in three phases: identify merge pairs, apply merges, then compact the face array.

*Phase 1 — Identify merge pairs:*
- Iterate faces in index order. For each unmerged triangle face `fi`, examine each edge.
- For each manifold edge shared with a neighbor `ni` (where `ni > fi` to avoid processing edges twice), both must be unmerged triangles with valid normals.
- Find the 4th vertex `b` from `ni` (the vertex not on the shared edge). Three vertices (`s0`, `s1`, `a`) define face `fi`'s plane.
- **Coplanarity check**: compute the signed distance from `b` to `fi`'s plane (using `fi`'s unit normal and one of its vertices). If `|distance| <= vertex_tolerance`, the vertex is coplanar within the data's numerical accuracy. (Using `vertex_tolerance` directly is correct since it represents the approximate precision of STL vertex positions — values within this tolerance are indistinguishable from coplanar.)
- **Convexity check**: verify the quad `[s1, a, s0, b]` is strictly convex by checking that the cross product at each of the 4 consecutive vertex triples agrees in direction with the face normal (positive dot product). Non-convex quads would break Newell's method for normal computation.
- If both checks pass, mark both `fi` and `ni` as merged and record the pair. Each triangle participates in at most one merge (greedy: first valid edge wins for each face).

*Phase 2 — Apply merges:*
- For each merge pair `(fi, ni)`: transform `fi` into a quad with `vertex_count = 4`, vertex order `[s1, a, s0, b]` (right-hand rule consistent with `fi`'s winding), and neighbors assembled from both triangles' non-shared edges. Recompute the normal via Newell's method. Update `ni`'s former neighbors to point to `fi`. Mark `ni` as deleted (`vertex_count = 0`).

*Phase 3 — Compact:*
- Remove deleted faces, build an old→new index mapping, remap all neighbor references.
- Update `mesh_faces` stat and record `mesh_quads` count.

**Constraints:**
- Each triangle participates in at most one merge (greedy: first valid merge wins).
- Only merge if the result is a convex quad — this ensures Newell's method produces a correct normal.
- Triangles at surface boundaries (where two curved surfaces meet, or a curved surface meets a plane) will generally NOT be coplanar with their neighbor across that boundary, so they correctly remain as triangles.
- On finely-tessellated curved surfaces (spheres, tori), adjacent quads may be nearly coplanar within `vertex_tolerance`. These additional merges are valid — they reduce face count and create more single-face candidates for curved surface hypothesis seeding.
- Sphere poles: triangles near poles may not have a coplanar neighbor (they border 5 or 6 other triangles in a fan pattern, not a quad strip). These correctly remain as triangles and become single-face planar hypotheses — good candidates for spherical hypothesis seeding.
---

### Stage 2: Surface Fitting

*Comment: The approach described here is "region growing" from seed faces—this works but may produce suboptimal results because the first seed determines the region boundaries. An alternative is RANSAC-style fitting: randomly sample minimal point sets, fit surfaces, count inliers, and keep the best. RANSAC is more robust to the mesh traversal order. The region-growing approach described here also doesn't naturally handle the case where the same surface type appears in disconnected regions (e.g., two separate planar faces with the same orientation). Consider a hybrid: use RANSAC or global clustering first, then refine with region growing.*

**Response: Given vertices that lie within epsilon of precisely on the surfaces involved, hopefully the first seed will grow to the entire region. I think the depth-first search or breadth-first search should effectively find planar surfaces that are partially disconnected. If they're fully disconnected, we'll consider them independent planar surfaces, which is fine.**

#### 2.1 Deduce planar hypotheses
Fit planar hypotheses to all connected sets of coplanar faces. Every face is assigned to exactly one planar hypothesis (including single-face hypotheses for faces on curved surfaces). A planar hypothesis consists of:
- `normal: [f64; 3]` — Unit normal vector pointing outward from the shell/solid.
- `distance: f64` — Signed distance from origin to plane along the normal (plane equation: normal · p = distance).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this plane (unordered set, not right-hand-rule ordered — boundary representation is deferred to stage 3).
- `error_max, error_min: f64` — Maximum (positive) and minimum (most negative) signed distance from any vertex to the plane.
- `error_abs_sum: f64` — Sum of absolute vertex-to-plane distances.

The algorithm uses BFS region growing (not DFS as originally proposed — BFS explores more uniformly and avoids the "narrow corridor drift" problem where DFS could accumulate error along a thin strip before reaching the main planar region):
- Initialize all `planar_hypothesis` to `UNDEDUCED_PLANAR_HYPOTHESIS` (-2).
- For every face index fi:
    - If `face[fi].planar_hypothesis != -2`: continue (already assigned).
    - Create a new planar hypothesis hi. Seed it with face fi's normal and compute the initial plane distance as the average dot product of fi's vertices with the normal.
    - Assign fi to hypothesis hi. Add fi to a BFS queue.
    - BFS loop — pop faces from the queue and examine each neighbor ni:
        - If ni is already assigned to a hypothesis, skip.
        - **Vertex distance check**: compute signed distance from each vertex of ni to the current plane. If all are within `--vertex-tolerance` (default 1e-5 mm), accept directly. If any vertex exceeds `2 * vertex_tolerance`, skip immediately — the face is clearly non-coplanar and re-fitting cannot help (since the plane can shift by at most ~tolerance while keeping existing vertices in tolerance). (No separate angular alignment check is needed — with tight vertex tolerance, any triangle with all vertices in the tolerance band necessarily has a closely aligned normal.)
        - **Re-fit attempt**: if any vertex exceeds tolerance but none exceeds `2 * tolerance`, try re-fitting the plane using area-weighted normal averaging over all current faces plus ni, with the distance as the average projection of all current+new vertices. Then check that **all** vertices (existing and new) are within tolerance of the re-fitted plane. If so, accept the re-fit; otherwise skip ni. (This addresses the comment's concern about previously-accepted faces no longer fitting — the re-fit validation checks all vertices, not just the new face's.)
        - If accepted: assign ni to hypothesis hi, add its vertices to the vertex set, push ni onto the BFS queue.
    - After BFS completes: final re-fit using all collected faces and vertices (area-weighted normal averaging + vertex centroid). Compute error metrics (error_max, error_min, error_abs_sum).
    - Single-face hypotheses are kept (not deleted). They represent faces on curved surfaces that will be assigned to cylindrical/spherical hypotheses in later stages.

Plane fitting uses **area-weighted normal averaging**: each face's unit normal is weighted by the face's triangle area, then the weighted sum is normalized. The plane distance is the mean of `normal · vertex` over all vertices in the hypothesis. This gives larger faces more influence on the plane orientation, which is more robust than unweighted averaging.

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the fitted plane (perpendicular projection).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - (Changed from original plan: the plan described projecting face *edges* onto the hypothesis and checking the enclosed surface region. The centroid-projection approach is simpler and sufficient for validation — it confirms the fitted plane aligns with a real STEP surface at representative interior points. Edge projection would be needed for boundary accuracy validation, which is a stage 3 concern.)

*Comment on 2.1: ~~The DFS vs BFS question is worth considering: DFS can get "trapped" in a narrow corridor and accumulate drift before exploring the main region. BFS explores more uniformly. However, the re-fitting step partially mitigates this.~~ Resolved: BFS was chosen for the implementation. ~~A bigger issue: once you re-fit the plane, previously accepted faces might no longer fit the new plane! Consider a final validation pass that removes faces whose vertices exceed the tolerance after the final fit.~~ Resolved: the re-fit step checks all vertices (existing + new) before accepting, and a final re-fit with error metrics is computed at the end. ~~Also: the "vertices in right-hand-rule order" for the hypothesis is unclear—planar regions aren't simply connected in general (they can have holes), so you'll need a more complex boundary representation.~~ Resolved: vertices are stored as an unordered set; boundary representation is deferred to stage 3.*

#### 2.2 Deduce cylindrical hypotheses
Fit cylindrical hypotheses to connected sets of faces that lie on a common cylinder. Only faces with single-face planar hypotheses (from stage 2.1) are candidates — faces belonging to multi-face planar hypotheses are genuinely flat and get `cylindrical_hypothesis = NO_HYPOTHESIS` (-1). (This works cleanly because stage 1.3 merges coplanar triangle pairs into quads, so curved-surface facets become single quad faces with single-face planar hypotheses.) A cylindrical hypothesis consists of:
- `axis_origin: [f64; 3]` — A point on the cylinder axis.
- `axis_direction: [f64; 3]` — Unit direction vector along the axis.
- `radius: f64` — Radius of the cylinder (always positive).
- `convex: bool` — Whether face normals point away from the axis (convex, like the outside of a pipe) or toward it (concave, like the inside of a hole).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this cylinder.
- `error_max: f64` — Maximum absolute distance from any vertex to the cylinder surface (i.e. max of `| ||v - axis_closest|| - radius |`).
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

The algorithm uses BFS region growing, analogous to stage 2.1 but seeded from pairs of adjacent faces (since a single triangle is always planar and cannot determine curvature):

**Seeding:**
- Initialize all `cylindrical_hypothesis` to `UNDEDUCED` (-2).
- For each face fi where `cylindrical_hypothesis == -2`:
    - If fi belongs to a multi-face planar hypothesis: set `cylindrical_hypothesis = -1`, skip.
    - Search fi's neighbors for a seed partner ni: an adjacent face that is also unassigned (`cylindrical_hypothesis == -2`), also has a single-face planar hypothesis, and has a sufficiently different normal (their cross product magnitude exceeds a minimum threshold, e.g. 0.01, to avoid near-parallel normals that would produce an unreliable axis estimate).
    - If no valid seed partner found: set `cylindrical_hypothesis = -1`, skip. (This face will be considered for spherical or other hypotheses in later stages.)
    - Estimate initial cylinder parameters from the seed pair (fi, ni) — see Cylinder Fitting below.
    - Verify seed validity: check that all vertices of fi and ni are within `--vertex-tolerance` of the estimated cylinder surface. If not, try the next neighbor. If no neighbor produces a valid seed, set `cylindrical_hypothesis = -1` and skip.
    - Determine convexity: compute the vector from the axis to the centroid of fi. If fi's face normal points in the same direction as this vector (positive dot product), the cylinder is convex; otherwise concave.
    - Create a new cylindrical hypothesis, assign both fi and ni, add both to a BFS queue.

**BFS expansion:**
- Pop faces from the queue and examine each neighbor ni:
    - If ni is already assigned a cylindrical hypothesis, skip.
    - If ni belongs to a multi-face planar hypothesis, skip.
    - **Vertex distance check**: for each vertex of ni, compute `| ||v - axis_closest|| - radius |`. If all are within `--vertex-tolerance`, accept directly. If any exceeds `2 * vertex_tolerance`, skip immediately (re-fitting cannot help, same reasoning as for planar hypotheses). If between 1x and 2x, proceed to re-fit attempt.
    - **Convexity check**: compute the vector from the axis to ni's centroid. Check that ni's face normal agrees with the hypothesis convexity (dot product with radial vector has the expected sign). Reject if inconsistent — this prevents merging the inner and outer surfaces of a thin-walled cylinder.
    - **Re-fit attempt**: if any vertex exceeds tolerance but none exceeds 2x, re-fit the cylinder from all current faces plus ni (see Cylinder Fitting below). Check that **all** vertices (existing and new) are within tolerance of the re-fitted cylinder. If so, accept the re-fit; otherwise skip ni.
    - If accepted: assign ni to the hypothesis, add its vertices to the vertex set, push ni onto the BFS queue.
- After BFS completes: final re-fit from all accumulated faces and vertices. Compute error metrics.
- **Centroid validation**: After BFS completes and error metrics are computed, validate that **all** face centroids lie within `--surface-tolerance` of the fitted cylinder surface. This catches spurious fits (e.g. two perpendicular planar faces that happen to fit a cylinder algebraically but whose centroids are far from any actual cylindrical surface). If validation fails, undo assignments and discard the hypothesis.

**Cylinder fitting** (used for seeding and re-fitting):

The fitting proceeds in two steps:

1. **Axis direction estimation from face normals.** On a cylinder, every face normal is perpendicular to the axis: `n · a = 0`. So the axis direction `a` minimizes $\sum w_i (n_i \cdot a)^2$ subject to $\|a\| = 1$, where $w_i$ is the area of face $i$. This is the eigenvector corresponding to the *smallest* eigenvalue of the 3×3 weighted covariance matrix $M = \sum w_i \, n_i \, n_i^T$. For the two-face seed, this simplifies to `a = normalize(n_fi × n_ni)`.

2. **Axis position and radius via 2D circle fitting.** Given the axis direction `a`, choose two orthogonal unit vectors `u`, `w` perpendicular to `a`. Project each vertex onto the 2D plane: `(x_i, y_i) = (v_i · u, v_i · w)`. Fit a circle to the 2D points using the algebraic least-squares method: fit the linear model $x^2 + y^2 + Dx + Ey + F = 0$ via least squares, then `center = (-D/2, -E/2)` and `radius = sqrt(center_x² + center_y² - F)`. The axis origin in 3D is `center_x * u + center_y * w` (the component along `a` is arbitrary).

The normal covariance matrix `M` can be maintained incrementally during BFS (add $w_i \, n_i \, n_i^T$ when accepting a face), but the eigenvector solve and circle fit must be recomputed from scratch when re-fitting, since a change in axis direction invalidates the 2D projection. This is O(n) in vertices but involves only simple arithmetic — no iterative optimization.

**Subtleties and edge cases:**
- **Near-parallel seed normals**: If the cross product of the two seed normals is small, the axis estimate is unreliable. Require `||n_fi × n_ni|| > MIN_CROSS_THRESHOLD` (e.g. 0.01, corresponding to ~0.6° between normals). Skip seed pairs that fail this.
- **Cones vs cylinders**: A cone has varying radius along the axis. If we fit a cylinder to conical faces, vertex distances near the ends will exceed tolerance, so cone faces will be naturally rejected. Cone fitting is deferred.
- **Large-radius cylinders**: A cylinder with very large radius looks locally planar. Such faces will belong to multi-face planar hypotheses from stage 2.1 (since vertices on a nearly-flat cylinder are nearly coplanar within tolerance). So they won't be considered for cylindrical fitting. This is the correct behavior — stage 2.6 will resolve any ambiguity.
- **Cylinder wrapping >180°**: This works naturally with BFS. The convexity check ensures we don't accidentally merge inner and outer surfaces.
- **Multiple cylinders sharing an axis**: E.g. stepped bore holes with different radii. These will form separate hypotheses because vertices at different radii will exceed the vertex tolerance.
- **Seam between cylinder and plane**: Faces at the junction have some vertices on the plane and some on the cylinder. Since the plane vertices aren't on the cylinder surface, the vertex distance check rejects them. This correctly excludes boundary faces.

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the cylinder surface (nearest point: project onto axis, then move radially to radius distance).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - Note: centroid-to-surface distance measures both fitting accuracy and tessellation sagitta (the flat triangle's centroid lies inside/outside the curved surface by approximately $h^2 / 8R$ where $h$ is the chord length). For typical tessellations this is small relative to `--surface-tolerance`, but it could become significant for very coarse meshes on tight-radius cylinders.

*Comment: ~~Cylinders are parameterized by axis (point + direction) and radius—6 DOF total. Minimum 5 points needed for a unique fit, but robust fitting requires more.~~ Addressed: the two-step fitting (normal eigenvector + 2D circle) handles this robustly from 2+ non-parallel faces. ~~Key considerations: (1) A cylindrical patch has principal curvature in one direction only—use this to distinguish from spheres/cones.~~ Addressed: the axis direction estimate will be well-defined for cylinders (normals span a plane perpendicular to axis) but degenerate for spheres (normals span all directions, smallest eigenvalue is not well-separated). ~~(2) "Negative radius" isn't the right framing; instead, track whether the surface normal points toward or away from the axis (convex vs concave).~~ Addressed: `convex` field + convexity check during BFS. ~~(3) For cones: 7 DOF. Cones degenerate to cylinders when half-angle→0, so you might fit cones first and detect near-zero angles.~~ Deferred: cone fitting will be considered as a future enhancement. For now, conical faces will fail the cylinder vertex distance check and remain as single-face planar hypotheses. ~~(4) Watch out for nearly-planar cylindrical patches (large radius)—they may fit planes better.~~ Addressed: large-radius cylinders are captured by multi-face planar hypotheses in stage 2.1 and excluded from cylinder candidates.*

#### 2.3 Deduce spherical hypotheses
TODO. Be sure to handle surfaces with negative curavature correctly - negative radius?

*Comment: Spheres are 4 DOF (center + radius). Minimum 4 non-coplanar points for a unique fit. For "negative curvature" (concave spherical patch), track normal orientation relative to center rather than using negative radius. Key challenge: partial spherical patches are hard to distinguish from cylinders or even planes if the patch is small relative to the radius. Consider requiring a minimum angular extent or using curvature analysis to disambiguate. Also consider toroidal surfaces (donuts, fillets)—they're common in CAD and combine characteristics of cylinders and spheres.*

#### 2.4 Deduce ruled surface hypotheses
TODO Optional. Find mesh which is coplanar on one axis, and model as an extruded curve surface/ruled surface.

*Comment: This is a good idea for capturing linear extrusions and sweeps. A ruled surface is defined by two boundary curves with linear interpolation between them. For extrusions, one "curve" is a point (the surface degenerates to a generalized cylinder). Detection: look for parallel mesh edges that share the same direction. Fitting: project to a plane perpendicular to the ruling direction and fit a 2D curve. Watch out for twisted ruled surfaces (rulings aren't parallel)—these are harder to detect and fit.*

#### 2.5 Deduce NURBS hypotheses
TODO: for groups of adjacent faces which are covered by one- or two-face planar hypotheses and not cylindrical or spherical hypotheses, try to fit a NURBS or b-spline surface to the vertices.

*Comment: NURBS fitting is complex and requires careful consideration: (1) Parameterization: you need to assign (u,v) parameters to each mesh vertex before fitting. Common approaches: conformal mapping, Floater's mean value coordinates, or discrete harmonic mapping. (2) Degree and knot selection: start with bicubic (degree 3×3); knot placement can use chord-length parameterization or be optimized. (3) Regularization: without it, the surface may oscillate. Consider smoothness penalties. (4) An alternative worth considering: use OCCT's `GeomAPI_PointsToBSplineSurface` which handles much of this automatically. (5) For "freeform" regions that are nearly planar, a plane with small tolerance may be preferable to a NURBS that overfits noise.*

#### 2.6 Select surfaces to use for reconstruction
- Iterate until out of valid hypotheses:
    - Select the hypothesis that fits some metric TODO of fitting the most area precisely. Add it to a list of selected surfaces.
    - Mark all faces using that hypothesis used.
    - Delete those faces from all other hypotheses that use them. Delete or mark invalid any hypothesis that ends up with insufficient faces left.
- Every face should be covered by one selected hypothesis.
- With the --compare flag, for each selected surface:
  - Project the bounding edges of the faces that handled by the surface onto the surface.
  - Check that within the projected bounding edge, the surface is within --vertex-tolerance of the surface of a solid from the STEP file.

*Comment: This greedy selection has a potential failure mode: selecting a large-but-poor-fit surface early can fragment remaining regions into pieces too small to fit well. Consider: (1) A quality metric that balances area coverage AND fit quality (e.g., area × (1 - normalized_error)). (2) Penalizing hypotheses that would leave "orphan" faces (faces with no valid remaining hypothesis). (3) Preferring analytic surfaces (plane, cylinder, sphere) over NURBS when fits are comparable, since analytic surfaces are more robust for downstream operations. (4) A backtracking mechanism if selection leads to uncoverable faces. Also: what if a face has no hypothesis at all after this process? This needs explicit handling—possibly flag as error or create a single-face planar patch.*

---

### Stage 3: Surface Reconstruction

*Comment: This stage has the most complexity and is where most CAD reconstruction projects run into trouble. The core challenge is that surface intersections in 3D are numerically delicate—small perturbations in surfaces can cause large changes in intersection curves, or cause intersections to disappear entirely. Consider adding explicit tolerance parameters throughout, and building in diagnostic output for debugging.*

Create a vector of face descriptors, one for each selected surface hypothesis, containing:
    - Reference to the hypothesis.
    - OCCT surface object. (Infinite or bounded?) (Created in 3.1.)
    - OCCT face object. (Created in 3.4.)
    - Vector of indices of adjacent surface face descriptors, ordered topologically (consecutive faces are adjacent to each other). NOTE: if a face connects to this one twice (with some other face(s) in between), it will occur more than once in the vector.
    - Vector of edge wire indices, one for each adjacent surface face. Edge wire i connects this face to adjacent face i.
    - Vector of vertex point indices, once for each pair of adjacent faces. Vertex index i represents the intersection between this face, adjacent face i, and adjacent face i+1%N. If there is only one adjacent face, then this vector is empty; otherwise it contains N points. (This assumes we have a solid - we'll need to adjust this assumption if we ever handle non-solid bodies.)

*Comment: The "one adjacent face" case is actually impossible for a closed solid—every edge has exactly two adjacent faces, and every face has at least 3 edges, so at least 3 neighbors (possibly with repeats). The edge case you probably mean is a face that's topologically a disk vs a face with holes. Also, faces can be adjacent to themselves (e.g., a cylindrical face wrapping around has one edge connecting it to itself). The data structure should handle this explicitly.*

Also create a vector of edge wires, each containing:
    - Indices of two adjacent faces.
    - Indices of two adjacent vertices.
    - Vector<Geom_Curve> for intersections of adjacent faces. TODO: vector or just one curve?
    - For each adjacent face, a Vector<Geom2d_Curve> for intersections in the UV-space of that face's surface. TODO: vector or just one curve?
    - Whether there's a tangency relationship. (More details?)

*Comment: For the "vector or one curve" question: mathematically, two surfaces can intersect in multiple disjoint curves (e.g., a plane cutting through a torus). However, if your mesh connectivity is correct, each edge wire should correspond to exactly one connected component of the intersection, so one curve (which may be a composite/piecewise curve) should suffice. The 2D pcurves are essential for OCCT face construction—make sure they're computed consistently (same parameterization direction, matching endpoints). Consider storing orientation flags to track which direction along the curve corresponds to which adjacent vertex.*
And a vector of vertices, each containing:
    - Indices of N adjacent faces, in topological order. (The faces connected to this vertex.) NOTE: There may be duplicates if a face connects to a vertex in multiple ways.
    - Indices of N adjacent wires, in topological order.
    - 3D point of vertex location.
    - Vector of N 2D points of vertex location in U-V space of the corresponding face.

*Comment: The vertex data structure looks correct. One refinement: the "3D point" should ideally be computed as the intersection of all adjacent surfaces (when possible), not just averaged from mesh vertices, to ensure the B-Rep is geometrically consistent. When three or more surfaces meet at a vertex, over-determination can cause the intersection to fail numerically—you may need least-squares fitting or to accept a small positional tolerance.*

#### 3.1 Create OCCT surface objects
- For each face descriptor:
    - Populate the OCCT surface object with a surface constructed from the hypothesis. Keep track of the mapping from hypothesis to surface index.

*Comment: OCCT surfaces (Geom_Plane, Geom_CylindricalSurface, etc.) are infinite. This is fine, but be aware that some surfaces have natural parameterization bounds (e.g., cylinder's U ∈ [0, 2π]). For consistency with OCCT conventions, ensure: planes use gp_Ax3 with Z as normal; cylinders/cones have Z along axis; spheres have poles at Z extremes. This affects pcurve computation later.*

- Second pass over face descriptors:
    - Populate the adjacent face indices with the indices of adjacent faces as determined from hypotheses sharing a hypothesis mesh vertex index.
    - For each adjacent face:   
        - Look up the adjacent edge wire descriptor, or create one if it doesn't exist yet. Populate it with:
            - Adjacent face indices.
            - Look up or create adjacent vertex descriptor indices. Populate them. For each adjacent vertex descriptor, populate:
                - Adjacent face indices.
                - Adjacent wire indices.

#### 3.2 Detect and create tangency relationships
- Detect edges between faces where there is numerically a very close to tangent relationship. Mark those edges as tangent.
- TODO: should we modify the surfaces to be numerically tangent? This may be challenging - if a surface is numerically close to tangent to two or more adjacent faces, we may need to do some sort of global optimization or iterate to a fixpoint. I suppose we can try to achieve numerical tangency, and if that's not possible, at least ensure that there's intersection along the extend of the shared edge.

*Comment: Tangent detection is critical for fillets and blends. Suggested approach: compute surface normals at several sample points along the shared mesh boundary; if normals agree within tolerance (e.g., < 0.1°), mark as tangent. For enforcing tangency: modifying analytic surfaces is usually wrong (it changes the geometry), but for NURBS you can add tangency constraints to the fit. A more robust approach: accept near-tangency and use a larger intersection tolerance when computing the shared edge. Also consider G2 (curvature) continuity for high-quality fillets—this matters for rendering/machining but may be overkill for your use case.*

#### 3.3 Create OCCT edge wires
- For each edge wire:
    - Compute the intersection between the adjacent faces to create candidate wires - there may be more than one, except for faces with a tangent relationship, special handling may be necessary:
        - For plane tangent to cylinder, create a linear wire where the cylinder's normal matches the plane's normal.
        - For cylinder tangent to cylinder, create a linear wire where the two cylinders' normals match.
        - For cylinder tangent to sphere, create a circular wire where the normals match.
        - ...?

*Comment: The tangent cases are well-identified. Additional cases: sphere tangent to sphere (point contact—degenerate), torus tangent to plane (circle), cone tangent to plane (line through apex or ellipse). For general surface-surface intersection, use OCCT's `GeomAPI_IntSS`. However, IntSS can fail or produce spurious curves for near-tangent surfaces. Consider: (1) Using the mesh boundary as a guide—project mesh boundary vertices onto both surfaces and fit a curve. (2) For analytic surfaces, compute intersections analytically when possible (they have closed-form solutions). (3) Always validate that the computed curve lies on both surfaces within tolerance.*
    - If there are vertices adjacent to this edge, then cut the wire at each vertex. I'm not sure how best to do this, perhaps either cut the wire with one or more adjacent surfaces (what if they're tangent?), or in a separate pass find the intersection of all of the wires at a vertex and cut them there.

*Comment: Cutting at vertices is the right idea but tricky in practice. Recommended approach: (1) Compute all edge curves first (full intersection curves, not yet trimmed). (2) For each B-Rep vertex, compute its 3D position as intersection of three surfaces, or use mesh vertex position as initial guess and project onto each surface. (3) For each edge curve, find the parameter values corresponding to the vertex positions and trim. Use `GeomAPI_ProjectPointOnCurve` for this. The vertex position should be consistent across all edges meeting there—if not, you have a gap that needs tolerance handling.*
    - Take all of the cut wires, and pick a set which is closest to the mesh vertices which lie along the edge, and which chain together to connect the adjacent vertices (or which form a loop, if there are no adjacent vertices).
    - Assemble that set into an OCCT curve to store into the edge wire descriptor.

*Comment: Using mesh vertices as a guide for selecting among multiple intersection curves is a good idea. Be aware that mesh vertices along an edge may not lie exactly on the fitted surfaces (due to fitting error), so use projection rather than direct distance. For seams (edges where a face is adjacent to itself, like a cylinder's wraparound), you won't have two distinct intersection curves—instead, you need to create an iso-parametric curve at a U or V seam of the surface.*
    - Store the vertices at the end of the wire into the adjacent vertex descriptors. Give some kind of error if the vertex is too far from any existing vertex. (TODO: what should we do here? Average them? Store a different vertex coordinate per edge?)
    - TODO: should we represent each wire in 3D space and UV space of each face during the operation? Or compute in one coordinate system and convert to the others after?

*Comment: For the vertex position question: OCCT's B-Rep model requires that all edges meeting at a vertex share the same TopoDS_Vertex (same 3D point). If edges disagree about vertex position, you have a few options: (1) Average and accept the tolerance. (2) Use the vertex with smallest fitting error. (3) Re-fit surfaces with constrained vertex positions. Option 1 is most practical. OCCT TopoDS_Vertex has a tolerance field specifically for this—set it to the maximum deviation. For 3D vs UV: compute in 3D, then derive pcurves using `ShapeAnalysis_Curve::ProjectPointOnCurve` or by evaluating surface inverse. OCCT's BRepBuilderAPI_MakeEdge can create edges with 3D curve + pcurves on both faces simultaneously.*

#### 3.4 Create OCCT faces
- For each face descriptor:
    - Populate OCCT face object via surface bounded by wires extracted from adjacent edge wire descriptors.

*Comment: This step uses `BRepBuilderAPI_MakeFace`. Key considerations: (1) The outer wire must be oriented counter-clockwise when viewed from outside the solid (along the face normal). (2) Inner wires (holes) must be clockwise. (3) Wires must be closed and edges must connect end-to-end within tolerance. (4) If face construction fails, ShapeFix_Face can often repair minor issues. (5) For periodic surfaces (cylinders, spheres), ensure pcurves handle the seam correctly—you may need to add a seam edge explicitly. This is often the most debugging-intensive step.*

#### 3.5 Construct Shells
- Find sets of connected face descriptors via DFS over the face graph (expore from each face to adjacent faces). For each set:
    - Stitch together an OCCT shell.

*Comment: Use `BRepBuilderAPI_Sewing` for this. Key settings: (1) Set sewing tolerance based on your vertex tolerance from earlier. (2) Enable "SameParameterMode" to ensure edge geometry is consistent. (3) After sewing, check `SewedShape()` for the result. Sewing can merge edges that are geometrically close—this is usually desired but verify the topology matches your intent. If sewing produces a compound instead of a shell, faces weren't connected properly. Also verify shell orientability—`ShapeAnalysis_Shell::CheckOrientedShells` can detect Möbius-strip-like errors.*

#### 3.6 Construct Solids
- Convert shells to solid bodies. 
    - TODO: figure out which shells are voids? Does face orientation help here?

*Comment: Yes, face orientation is the key. OCCT convention: face normals point outward from material. For a solid: outer shell normals point out, void shell normals point in (toward the void interior, away from material). To classify: compute signed volume of each shell—positive = outer, negative = inner. Alternatively, pick a point inside the shell and ray-cast to determine if it's inside any other shells. Use `BRepBuilderAPI_MakeSolid` to combine an outer shell with void shells. `ShapeFix_Solid::SolidFromShell` can also create a solid from a single closed shell and orient it correctly. Final validation: `BRepCheck_Analyzer` will verify the solid is valid.*

---

### Stage 4: Output

#### 4.1 Output objects
- Write constructed objects to a STEP file (or potentially other formats).

*Comment: Use `STEPControl_Writer` with `STEPControl_AsIs` mode to preserve your exact geometry. Consider also offering `STEPControl_ManifoldSolidBrep` mode which enforces stricter solid validity. Before export, run `ShapeFix_Shape` as a final cleanup pass (but heed your AGENTS.md warning about investigating root causes of any fixes). Set appropriate STEP header metadata (author, organization, etc.) for traceability. For debugging, also consider outputting intermediate formats: BREP (OCCT native), or individual surfaces/curves to help diagnose reconstruction issues.*

---

## Source Code Organization

Each stage will be represented by a function which takes a ref to a configuration data structure and consumes an input data structure, and returns an output data structure (or error). That data structure is passed to the next stage:

- main:
  - let config = parse flags, load --compare STEP file
  - out1 = stage1(config).unwrap()
  - if config.stage < 2: exit
  - out2 = stage2(congig, out1).unwrap()
  - if config.stage < 3: exit
  - out3 = stage3(config, out2).unwrap()
  - if config.stage < 4: exit
  - stage4(config, out3).unwrap()  # writes output file
- stageN:
  - execute whatever portions of stage1 are allowed by config.stage

Each stage should have a source file stageN.rs, with a definition of that stage's output data structure at the top, the stageN() function, and then whatever other functions are required to implement the stage's functionality.

---

## Implementation Phases

### Foundation
- [x] Project setup (crates, cargo.toml, dependencies)
- [x] Test utility: read an STL and a STEP, compute maximum distance between STL vertices and STEP surfaces, and print it out. Create a script in scripts/ to apply it to all of the stl/step file pairs under tests/ and print out a table of maximum distances.
- [x] Program skeleton: implement main program, flag parsing, and create stub source files for each stage. Create stage output data structures. For portions of the data structures which are unclear at this point, stub them out with comments. Take existing stage 1.1 and stage 1.2 implementations in mesh.rs and reformulate them into stage1.rs and appropriate output data structures.

### Stage 1 Mesh Input
- [x] Stage 1.1: Read STL file into `ConnectedMesh`, including welded vertices and per-face placeholder fields for neighbors, normals, and hypotheses.
- [x] Stage 1.2: Mesh validation pass to compute face normals, edge neighbors, manifold stats, connected shells, and orientation consistency checks.
- [x] Stage 1.3: Coplanar triangle fusion — merge adjacent coplanar triangle pairs into quads. This simplifies surface fitting by ensuring curved-surface facets (cylinder/sphere quad strips) become single faces rather than 2-face planar hypothesis groups.
- [x] Implement tests that all stl/step file pairs in tests/ pass consistency checks and pass when fed to brepper with the --compare flag. Also invent a few file pairs in tests/bad that will fail with the --compare_step flag by editing the location of one or more vertices to fail surface closeness tests. Also invent some bad cases that fail the mesh validation tests in various ways. Ensure that these bad cases fail in the correct way.

### Stage 2: Surface Fitting
- [x] Understand stage 2.1 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.1. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [x] Implement stage 2.1. Make sure --compare passes for all test shapes composed only of planar surfaces.
- [x] Understand stage 2.2 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.2. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [x] Implement stage 2.2: Deduce cylindrical hypotheses.
- [ ] Understand stage 2.3 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.2. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [ ] Implement stage 2.3: Deduce spherical hypotheses.

---

## Testing Strategy

The main testing strategy is to process a set of example stl/step pairs, and use the --compare flag to ensure that the pipeline is working correctly for each. 

1. **Unit Tests**: Individual components (readers, fitters, converters)
2. **Integration Tests**: Full pipeline on known geometries
3. **Regression Tests**: Compare output STEP to reference
4. **Validation**: Round-trip test (STEP → mesh → STEP)

### Test Cases

| Test Case | Input | Expected Output | Status |
|-----------|-------|-----------------|--------|
| Cube | 12 triangles | 6 planes, 1 solid | ✓ Stage 2.1 |
| Wedge | 12 triangles | 6 planes (incl. angled) | ✓ Stage 2.1 |
| T-Shape | 28 triangles | 10 planes | ✓ Stage 2.1 |
| Staircase | 48 triangles | 12 planes | ✓ Stage 2.1 |
| Chamfered Cube | 44 triangles | 26 planes | ✓ Stage 2.1 |
| Stepped Block | Complex planar | Multiple planes | ✓ Stage 2.1 |
| L Bracket | Complex planar | Multiple planes | ✓ Stage 2.1 |
| Cylinder | Tessellated cylinder | 1 cylinder + 2 planes | ✓ Stage 1 |
|| Simple Cylinder (ccad) | 124 triangles | 1 cylinder + 2 planes | ✓ Stage 2.2 |
|| Block with Hole (ccad) | 44 triangles | 6 planes + 1 concave cylinder | ✓ Stage 2.2 |
|| Pipe (ccad) | 244 triangles | 2 cylinders (in/out) + 2 annular planes | ✓ Stage 2.2 |
|| Stepped Cylinder (ccad) | 240 triangles | 2 cylinders + 3 planes | ✓ Stage 2.2 |
|| Two Holes (ccad) | 244 triangles | 6 planes + 2 concave cylinders | ✓ Stage 2.2 |
| Sphere | Tessellated sphere | 1 sphere | ✓ Stage 1 |
| Cone | Tessellated cone | 1 cone + 1 plane | ✓ Stage 1 |
| Fillet | Blended edge | Planes + fillet surface | |
| Complex part | Real CAD export | Matching topology | |

---

## Performance Targets

| Mesh Size | Target Time |
|-----------|-------------|
| < 10K triangles | < 5 seconds |
| 10K - 100K triangles | < 30 seconds |
| 100K - 1M triangles | < 5 minutes |

---

## Future Enhancements

- [ ] Support for OBJ, PLY input formats
- [ ] IGES export option
- [ ] Machine learning for surface type classification
- [ ] Hole detection and filling
- [ ] Feature recognition (holes, pockets, bosses)
