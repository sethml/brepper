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
        - `spherical_hypothesis: i32` — Index of active spherical hypothesis, or `NO_HYPOTHESIS` (-1) if none, or `UNDEDUCED` (-2) if not yet deduced.
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
- With the `--compare_step` flag: check that all vertices are within `--vertex-tolerance` of the STEP shape, and that all face centroids are within `--surface-tolerance`. Uses `BRepExtrema_DistShapeShape` for bounded distance (respects face trimming). Reports the worst vertex and worst centroid with their distances.

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

- With the `--compare_step` flag:
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
    - Search fi's neighbors for a seed partner ni: an adjacent face that is also unassigned (`cylindrical_hypothesis == -2`), also has a single-face planar hypothesis, and has a sufficiently different normal (their cross product magnitude exceeds a minimum threshold, e.g. 0.01, to avoid near-parallel normals that would produce an unreliable axis estimate). Additionally, the dihedral angle between fi and ni must not exceed `--angular-tolerance` (default 17.5°) — this prevents seeding from pairs of faces that meet at too sharp an angle (e.g., adjacent faces of a cube at 90°).
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
    - **Angular tolerance check**: compute the dihedral angle between ni and the current BFS-parent face. If it exceeds `--angular-tolerance`, reject. This limits how fast the surface can curve between adjacent faces.
    - **Re-fit attempt**: if any vertex exceeds tolerance but none exceeds 2x, re-fit the cylinder from all current faces plus ni (see Cylinder Fitting below). Check that **all** vertices (existing and new) are within tolerance of the re-fitted cylinder. If so, accept the re-fit; otherwise skip ni.
    - If accepted: assign ni to the hypothesis, add its vertices to the vertex set, push ni onto the BFS queue.
- After BFS completes: final re-fit from all accumulated faces and vertices. Compute error metrics.
- **Centroid validation**: After BFS completes and error metrics are computed, validate that **all** face centroids lie within `--surface-tolerance` of the fitted cylinder surface. This catches spurious fits (e.g. two perpendicular planar faces that happen to fit a cylinder algebraically but whose centroids are far from any actual cylindrical surface). If validation fails, undo assignments and discard the hypothesis.
- **Minimum face count**: Require at least 3 faces in a cylindrical hypothesis. Any real cylindrical surface from a CAD tessellation will produce at least 3 facets around its circumference. This eliminates spurious 2-face cylinder fits that arise from locally adjacent faces on non-cylindrical curved surfaces (tori, cones, spheres). Combined with centroid validation, this significantly reduces false positives.

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

- With the `--compare_step` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the cylinder surface (nearest point: project onto axis, then move radially to radius distance).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - Note: centroid-to-surface distance measures both fitting accuracy and tessellation sagitta (the flat triangle's centroid lies inside/outside the curved surface by approximately $h^2 / 8R$ where $h$ is the chord length). For typical tessellations this is small relative to `--surface-tolerance`, but it could become significant for very coarse meshes on tight-radius cylinders.

*Comment: ~~Cylinders are parameterized by axis (point + direction) and radius—6 DOF total. Minimum 5 points needed for a unique fit, but robust fitting requires more.~~ Addressed: the two-step fitting (normal eigenvector + 2D circle) handles this robustly from 2+ non-parallel faces. ~~Key considerations: (1) A cylindrical patch has principal curvature in one direction only—use this to distinguish from spheres/cones.~~ Addressed: the axis direction estimate will be well-defined for cylinders (normals span a plane perpendicular to axis) but degenerate for spheres (normals span all directions, smallest eigenvalue is not well-separated). ~~(2) "Negative radius" isn't the right framing; instead, track whether the surface normal points toward or away from the axis (convex vs concave).~~ Addressed: `convex` field + convexity check during BFS. ~~(3) For cones: 7 DOF. Cones degenerate to cylinders when half-angle→0, so you might fit cones first and detect near-zero angles.~~ Deferred: cone fitting will be considered as a future enhancement. For now, conical faces will fail the cylinder vertex distance check and remain as single-face planar hypotheses. ~~(4) Watch out for nearly-planar cylindrical patches (large radius)—they may fit planes better.~~ Addressed: large-radius cylinders are captured by multi-face planar hypotheses in stage 2.1 and excluded from cylinder candidates.*

#### 2.3 Deduce spherical hypotheses

Fit spherical hypotheses to connected sets of faces that lie on a common sphere. All faces are candidates for spherical fitting — faces with cylindrical hypotheses or multi-face planar hypotheses are not excluded, since a face can legitimately belong to multiple surface types (e.g., equator faces of a sphere may also fit a cylinder). Stage 2.6 resolves overlapping hypotheses later. A spherical hypothesis consists of:
- `center: [f64; 3]` — Center of the sphere.
- `radius: f64` — Radius of the sphere (always positive).
- `convex: bool` — Whether face normals point away from the center (convex, like the outside of a ball) or toward it (concave, like a bowl).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this sphere.
- `error_max: f64` — Maximum absolute distance from any vertex to the sphere surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

The algorithm uses BFS region growing, analogous to stages 2.1 and 2.2 but seeded from pairs of adjacent faces with non-parallel normals (analogous to cylindrical seeding). The algebraic sphere fit from 4+ non-coplanar vertices (from 2 adjacent triangles/quads) determines a unique sphere, and validation catches bad fits:

**Seeding:**
- For each face fi where `spherical_hypothesis == UNDEDUCED` (-2):
    - Search fi's neighbors for a seed partner ni: an adjacent face that is also unassigned (`spherical_hypothesis == UNDEDUCED`), and has a sufficiently different normal (cross product magnitude exceeds `MIN_CROSS_THRESHOLD` = 0.01). The dihedral angle between fi and ni must also not exceed `--angular-tolerance` (default 17.5°).
    - If no valid seed partner found: set `spherical_hypothesis = NO_HYPOTHESIS`, skip.
    - Estimate initial sphere parameters from the seed pair (fi, ni) using least-squares sphere fitting (see Sphere Fitting below). The combined vertices of two non-coplanar faces (at least 6 for two triangles, 8 for two quads) provide enough constraints for the 4-DOF sphere fit.
    - Verify seed validity: check that all vertices of fi and ni are within `--vertex-tolerance` of the estimated sphere surface.
    - **Radius sanity check**: reject the seed if the fitted radius exceeds `max_sphere_radius` (see below). This prevents degenerate fits where a very large sphere approximates locally flat chamfer/bevel faces.
    - Determine convexity: the vector from sphere center to fi's centroid should align with fi's normal (positive dot product → convex, negative → concave).
    - Create a new spherical hypothesis, assign both seed faces, add to BFS queue.

*Design note: The original plan called for 3-face seeds with a normals-span-3D eigenvalue check, but this failed on finely tessellated spheres where 3 adjacent faces have nearly identical normals (spanning only ~6°, below the eigenvalue ratio threshold). Pair-based seeding works because the sphere fit itself validates geometry — if the vertices don't lie on a sphere, the fit produces a large radius or high error, which is caught by max_sphere_radius and vertex distance checks.*

**BFS expansion:**
- Pop faces from the queue and examine each neighbor ni:
    - If ni already has a spherical hypothesis, skip.
    - **Vertex distance check**: for each vertex of ni, compute `| ||v - center|| - radius |`. If all within `--vertex-tolerance`, accept directly. If any exceeds `2 * vertex_tolerance`, skip. If between 1x and 2x, re-fit.
    - **Convexity check**: verify ni's normal agrees with the hypothesis convexity (dot product of normal with radial vector has the expected sign).
    - **Angular tolerance check**: compute the dihedral angle between ni and the current BFS-parent face. If it exceeds `--angular-tolerance`, reject.
    - **Re-fit attempt**: re-fit sphere from all current faces plus ni. Check all vertices within tolerance and radius within `max_sphere_radius`. Accept or reject.
    - If accepted: assign ni, add vertices, push onto BFS queue.
- After BFS completes: final re-fit, compute error metrics.
- **Centroid validation**: all face centroids must be within `--surface-tolerance` of the fitted sphere. Discard hypothesis if validation fails.
- **Minimum face count**: require at least 4 faces (a sphere needs at least 4 non-coplanar points for a unique fit).
- **Maximum radius**: `max_sphere_radius = bounding_box_diagonal * 10`, where `bounding_box_diagonal` is the diagonal of the axis-aligned bounding box of all mesh vertices. Compute the bounding box once at stage 2 entry. This limit prevents degenerate sphere fits: on a chamfered cube, corner chamfer faces have normals spanning 3D (e.g., [1,1,1]/√3, [1,1,0]/√2, [1,0,1]/√2) and pass the eigenvalue check, but they're flat — a sphere of radius R ≈ L²/(8·vertex_tolerance) could fit them numerically. For a 14mm model with vertex_tolerance=1e-5, R ≈ 112 km, far exceeding `14mm × 10 = 140mm`. Real spherical features (bearings, domes, lenses) have radius within ~10× the model extent.

**Sphere fitting** (used for seeding and re-fitting):
- Given vertices on a sphere: $|v - c|^2 = r^2$, expand to $v \cdot v - 2 v \cdot c + |c|^2 = r^2$.
- Rearrange: $v \cdot v = 2 v \cdot c - |c|^2 + r^2 = 2 v \cdot c + k$ where $k = r^2 - |c|^2$.
- This is linear in unknowns $(c_x, c_y, c_z, k)$. Solve via least squares: $A x = b$ where $A_{i} = [2v_{ix}, 2v_{iy}, 2v_{iz}, 1]$ and $b_i = v_i \cdot v_i$.
- After solving: center $= (c_x, c_y, c_z)$, radius $= \sqrt{k + |c|^2}$.

**Distinguishing spheres from cylinders:**
- Faces with cylindrical hypotheses are NOT excluded — they may legitimately belong to both a cylinder and a sphere (e.g., equator band). Stage 2.6 disambiguates.
- On cylindrical surfaces, sphere fits from pairs of nearly-parallel-normal faces tend to produce unreliable results (very large or negative r²), which are caught by `max_sphere_radius` and the `r² > 0` check in sphere fitting.
- The pair-based seeding requires non-parallel normals (cross product > 0.01), which naturally excludes seeds from planar regions but allows seeds from curved regions of any type.

**Preventing false positives on chamfers/bevels:**
- On a chamfered or beveled cube, corner and edge chamfer faces can seed sphere fits. However, these faces are flat — the sphere fit produces a degenerate radius vastly larger than the model. The `max_sphere_radius` constraint (10× bounding box diagonal) catches this. For a 14mm model with vertex_tolerance=1e-5, a degenerate sphere fit would have R ≈ 112 km, far exceeding `14mm × 10 = 140mm`.
- The centroid validation check also catches spurious fits where the algebraic fit succeeds but the sphere doesn't actually pass through the face centroids.

**Test models:**
- Simple sphere (ccad): full convex sphere, radius 10 → 1 convex spherical hypothesis.
- Hemisphere (ccad): top half of sphere + flat base → 1 convex spherical hypothesis + 1 multi-face planar.
- Spherical pocket (ccad): block with concave hemispherical cavity → 6 planar + 1 concave spherical.
- Ball on cylinder (ccad): sphere atop cylinder stalk → 1 convex spherical + 1 convex cylindrical + 1 planar.
- Sphere (onshape): full sphere, fine tessellation.
- Dome/hemisphere (onshape): hemisphere with flat base.

*Comment: Spheres are 4 DOF (center + radius). Minimum 4 non-coplanar points for a unique fit. For "negative curvature" (concave spherical patch), track normal orientation relative to center rather than using negative radius. ~~Key challenge: partial spherical patches are hard to distinguish from cylinders or even planes if the patch is small relative to the radius. Consider requiring a minimum angular extent or using curvature analysis to disambiguate.~~ Addressed: the max_sphere_radius constraint prevents degenerate large-radius fits, and pair-based seeding with vertex distance + centroid validation prevents false positives. Also consider toroidal surfaces (donuts, fillets)—they’re common in CAD and combine characteristics of cylinders and spheres. Note: torus faces locally fit spheres, producing sphere hypotheses on torus surfaces (e.g., pipe_elbow). Stage 2.6 will need torus hypothesis support to resolve this.*

#### 2.4 Deduce ruled surface hypotheses
TODO Optional. Find mesh which is coplanar on one axis, and model as an extruded curve surface/ruled surface.

*Comment: This is a good idea for capturing linear extrusions and sweeps. A ruled surface is defined by two boundary curves with linear interpolation between them. For extrusions, one "curve" is a point (the surface degenerates to a generalized cylinder). Detection: look for parallel mesh edges that share the same direction. Fitting: project to a plane perpendicular to the ruling direction and fit a 2D curve. Watch out for twisted ruled surfaces (rulings aren't parallel)—these are harder to detect and fit.*

#### 2.5 Deduce NURBS hypotheses
TODO: for groups of adjacent faces which are covered by one- or two-face planar hypotheses and not cylindrical or spherical hypotheses, try to fit a NURBS or b-spline surface to the vertices.

*Comment: NURBS fitting is complex and requires careful consideration: (1) Parameterization: you need to assign (u,v) parameters to each mesh vertex before fitting. Common approaches: conformal mapping, Floater's mean value coordinates, or discrete harmonic mapping. (2) Degree and knot selection: start with bicubic (degree 3×3); knot placement can use chord-length parameterization or be optimized. (3) Regularization: without it, the surface may oscillate. Consider smoothness penalties. (4) An alternative worth considering: use OCCT's `GeomAPI_PointsToBSplineSurface` which handles much of this automatically. (5) For "freeform" regions that are nearly planar, a plane with small tolerance may be preferable to a NURBS that overfits noise.*

#### 2.6 Select surfaces to use for reconstruction
Assign each mesh face to exactly one "selected surface" for reconstruction. Since stages 2.1–2.3 (and optionally 2.4–2.5) already assign faces to typed hypotheses with minimal overlap, surface selection is a simple per-face priority resolution.

**Per-face priority rule** (highest priority first):
1. **Multi-face planar hypothesis** (face count > 1): the faces are genuinely coplanar. Select the planar hypothesis.
2. **Spherical hypothesis**: if present, select it. (Sphere is more constrained than cylinder — 4 DOF vs 6 — so a valid sphere fit is more specific.)
3. **Cylindrical hypothesis**: if present, select it.
4. **Single-face planar hypothesis**: fallback for faces on surfaces not yet fitted (conical, toroidal, freeform). Stages 2.4–2.5 will replace some of these with ruled or NURBS surfaces.

**Algorithm:**
- For each mesh face, apply the priority rule to determine its selected hypothesis type and index.
- Build a `Vec<SelectedSurface>`, one per unique selected hypothesis, containing:
    - `surface_type: SurfaceType` — enum: Planar, Cylindrical, Spherical (later: Ruled, BSpline).
    - `hypothesis_index: usize` — index into the appropriate hypothesis vector in Stage2Output.
    - `faces: Vec<usize>` — mesh face indices assigned to this surface.
    - `vertices: Vec<usize>` — mesh vertex indices on this surface (union of face vertices).
- Validate: every face is assigned to exactly one selected surface. Error if any face has no valid hypothesis.
- Print a summary: count of selected surfaces by type, total faces covered, any uncovered faces.

This simple priority rule works because a face can have at most one hypothesis of each type. Since stage 2.3 now allows faces to have both cylindrical and spherical hypotheses, a face on the equator of a sphere that also fits a cylinder will be assigned to the sphere (higher priority). Multi-face planar hypotheses still take highest priority since those faces are genuinely coplanar.

*Future work: if the priority rule produces poor results (e.g., a large-radius cylinder incorrectly captured as a multi-face plane), consider fit-error-weighted selection or a greedy area-coverage approach. A more sophisticated metric could balance area × (1 - normalized_error). Backtracking could handle cases where greedy selection fragments remaining regions.*

- With the `--compare_step` flag, for each selected surface:
    - Compute face centroids and project them onto the selected surface.
    - Measure distance from each projected centroid to the nearest surface in the reference STEP file.
    - Report an error if any projected centroid exceeds `--surface-tolerance`.

---

### Stage 3: Surface Reconstruction

Build OCCT B-Rep topology (faces, edges, vertices) from the selected surfaces. This is the most complex stage — surface intersections in 3D are numerically delicate, and small perturbations can cause large changes in intersection curves or cause them to disappear.

**Core data structures:**

A vector of `ReconFace` descriptors, one per selected surface:
- `selected_surface_index: usize` — Index into the selected surfaces from stage 2.6.
- `surface: Handle<Geom_Surface>` — Infinite OCCT surface (created in 3.1).
- `adj_faces: Vec<usize>` — Indices of adjacent ReconFaces, ordered topologically around this face's boundary. A face may appear multiple times if it shares multiple edges with this face (e.g., a cylindrical face adjacent to itself via a seam edge).
- `adj_edges: Vec<usize>` — ReconEdge indices. `adj_edges[i]` is the edge between this face and `adj_faces[i]`.
- `adj_vertices: Vec<usize>` — ReconVertex indices at corners. `adj_vertices[i]` is the vertex between `adj_edges[i]` and `adj_edges[(i+1) % N]`.
- `occt_face: Option<TopoDS_Face>` — Populated in 3.4.

A vector of `ReconEdge` descriptors:
- `face_indices: [usize; 2]` — Indices of the two adjacent ReconFaces.
- `vertex_indices: [usize; 2]` — Indices of the ReconVertices at each end. For closed-loop edges (no vertices), both are `usize::MAX`.
- `curve_3d: Handle<Geom_Curve>` — 3D intersection curve, trimmed to vertex endpoints. One curve per edge suffices — each edge corresponds to one connected component of the surface-surface intersection (the mesh connectivity guarantees this).
- `pcurves: [Handle<Geom2d_Curve>; 2]` — 2D parametric curve on each adjacent face's surface.
- `tangent: bool` — Whether the adjacent surfaces are tangent along this edge (detected in 3.2).
- `mesh_boundary_vertices: Vec<usize>` — Mesh vertex indices along this boundary, used to guide curve selection and validate the intersection curve.

A vector of `ReconVertex` descriptors:
- `point: gp_Pnt` — 3D position, computed by projecting the mesh vertex onto all adjacent surfaces and averaging (rather than raw mesh vertex position, for B-Rep consistency).
- `adj_faces: Vec<usize>` — Indices of ReconFaces meeting at this vertex, in topological order.
- `adj_edges: Vec<usize>` — Indices of ReconEdges meeting at this vertex, in topological order.
- `uv_coords: Vec<(f64, f64)>` — UV coordinates of this vertex on each adjacent face's surface.
- `tolerance: f64` — Maximum deviation from any adjacent surface (used for OCCT's `TopoDS_Vertex` tolerance).

#### 3.1 Create OCCT surface objects and build adjacency graph

**Surface creation:**
- For each selected surface, create the corresponding infinite OCCT surface:
    - Planar → `Geom_Plane` from `gp_Pln` (gp_Ax3 with Z as face normal).
    - Cylindrical → `Geom_CylindricalSurface` from `gp_Cylinder` (axis direction along Z, origin on axis, radius).
    - Spherical → `Geom_SphericalSurface` from `gp_Sphere` (center, radius).
    - (Future: BSpline → `Geom_BSplineSurface` via `GeomAPI_PointsToBSplineSurface`, etc.)
- For consistency with OCCT conventions: planes use gp_Ax3 with Z as normal; cylinders have Z along axis; spheres have poles at Z extremes. This affects pcurve computation in 3.3.

**Adjacency graph construction:**
- Two selected surfaces are adjacent if they share mesh edges (neighboring mesh faces that belong to different selected surfaces).
- For each mesh edge between faces assigned to different selected surfaces:
    - Record the pair of selected surfaces.
    - Collect the mesh vertices along the shared boundary.
- Group shared mesh edges into contiguous boundary segments (connected sequences between the same pair of surfaces). Each segment becomes a ReconEdge.
- Identify B-Rep vertices: mesh vertices where 3+ selected surfaces meet. For each, create a ReconVertex.
    - Compute its 3D position by projecting the mesh vertex position onto all adjacent surfaces and averaging.
    - Record adjacent faces and edges in topological order (walking around the vertex).
- Populate the `adj_faces`, `adj_edges`, and `adj_vertices` vectors of each ReconFace from the constructed graph.

#### 3.2 Detect tangency relationships
For each ReconEdge, determine whether the two adjacent surfaces are tangent along the shared boundary. This matters because `GeomAPI_IntSS` can fail or produce degenerate results for tangent/near-tangent surfaces, requiring special handling in edge curve computation (3.3).

**Detection algorithm:**
- Sample several mesh vertices along the shared boundary (e.g., every 5th boundary vertex, minimum 3 samples).
- At each sample point, compute the surface normal of both adjacent surfaces.
- If all sampled normal pairs agree within a small angle threshold (e.g., 2°), mark the edge as tangent.

**Decision: do not modify surfaces to enforce tangency.** Modifying analytic surfaces would change the geometry. Instead, tangent edges get special handling in edge curve computation (3.3): construct the edge curve directly rather than relying on surface-surface intersection.

For the initial implementation, tangency detection can be deferred (mark all edges as non-tangent) since the early test models (planar + cylindrical + spherical intersections) typically don't have tangent edges. Tangent edges arise from fillets and blends, which will be addressed when those test models are added.

#### 3.3 Create edge curves
For each ReconEdge, compute the 3D intersection curve between the two adjacent surfaces, trim it to the vertex endpoints, and derive pcurves.

**Intersection computation:**
- For non-tangent edges: use `GeomAPI_IntSS` to compute intersection curves. IntSS may return multiple curves (e.g., a plane cutting through a torus); select the one closest to the mesh boundary vertices. For each candidate curve, project the centroids of mesh boundary edges onto the curve and pick the one with smallest total distance.
- For tangent edges: construct curves directly based on surface pair types:
    - Plane tangent to cylinder: line along the cylinder ruling where the cylinder normal matches the plane normal.
    - Cylinder tangent to cylinder: line along the shared ruling direction.
    - Cylinder tangent to sphere: circular arc where the normals match.
    - General tangent case: project mesh boundary vertices onto both surfaces and fit a `Geom_BSplineCurve` through the projected midpoints.

**Trimming to vertex endpoints:**
- For each ReconVertex at an edge endpoint, project the vertex's 3D position onto the intersection curve using `GeomAPI_ProjectPointOnCurve` to get the curve parameter value.
- Trim the curve to the parameter range `[t_start, t_end]`.
- For closed-loop edges (no vertices at endpoints), use the full intersection curve.

**Pcurve computation:**
- Compute all intersections in 3D, then derive pcurves from the 3D curve.
- For analytic curve-on-analytic-surface cases, compute pcurves analytically:
    - Line on plane → `Geom2d_Line`
    - Circle on plane, cylinder, or sphere → `Geom2d_Circle` or `Geom2d_Line` (for iso-parametric curves)
- For other cases: sample points along the 3D curve, project each onto the surface to get UV coordinates, and fit a `Geom2d_BSplineCurve`.
- Alternatively, use `ShapeFix_Wire::FixEdgeCurves` after edge creation to compute missing pcurves automatically.

**Seam edges:**
- A face adjacent to itself (e.g., a full cylinder wrapping around) requires a seam edge at a fixed U or V parameter. Create the seam as an iso-parametric curve on the surface (e.g., U=0 and U=2π for a cylinder). The 3D curve is the same for both parameterizations; the two pcurves differ by the period.

**Vertex position consistency:**
- OCCT requires all edges at a vertex to share the same `TopoDS_Vertex` (same 3D point). When edges computed independently disagree slightly about vertex position, use the averaged position from the ReconVertex and set the vertex tolerance to the maximum deviation.

#### 3.4 Create OCCT faces
For each ReconFace, construct a `TopoDS_Face` from the surface and bounding wires.

**Algorithm:**
- Collect all ReconEdges that border this face. For each, create a `TopoDS_Edge` using `BRepBuilderAPI_MakeEdge`, providing the 3D curve, pcurve on this face's surface, and the two `TopoDS_Vertex` endpoints.
- Order edges into closed wires. Each face has one outer wire and zero or more inner wires (holes). Determine wire connectivity by matching shared vertices.
- Orient wires: outer wire counter-clockwise when viewed from outside the solid (along the face normal), inner wires clockwise. Use the face normal and mesh vertex positions to determine orientation.
- Construct the face: `BRepBuilderAPI_MakeFace` with the surface and outer wire, then add inner wires.
- If face construction fails, attempt repair with `ShapeFix_Face`. Investigate and fix the root cause if repair is needed frequently.
- For periodic surfaces (cylinders, spheres), include the seam edge in the wire and ensure pcurves respect the parameter periodicity.

#### 3.5 Construct shells
Group connected ReconFaces into shells and assemble each group into a `TopoDS_Shell`.

**Algorithm:**
- Find connected components of the ReconFace adjacency graph via BFS/DFS.
- For each connected component, use `BRepBuilderAPI_Sewing` to stitch the `TopoDS_Face` objects into a shell:
    - Set sewing tolerance to `vertex_tolerance`.
    - The sewing operation merges shared edges and ensures consistent edge geometry.
- After sewing, verify the result is a `TopoDS_Shell` (not a `TopoDS_Compound`, which would indicate disconnected faces).
- Validate with `ShapeAnalysis_Shell::CheckOrientedShells` to detect orientation inconsistencies.
- If sewing produces errors, fall back to manual shell construction with `BRep_Builder`.

#### 3.6 Construct solids
Convert closed shells into `TopoDS_Solid` objects.

**Algorithm:**
- For each shell, check if it is closed (all edges shared by exactly 2 faces).
- Determine shell role using signed volume:
    - Compute signed volume via `BRepGProp::VolumeProperties`. Positive volume → outer shell, negative → inner shell (cavity/void).
    - For single-shell solids, use `ShapeFix_Solid::SolidFromShell` which handles orientation automatically.
- For multi-shell solids (outer shell containing voids): combine the outer shell with inner/void shells using `BRepBuilderAPI_MakeSolid`.
- Validate with `BRepCheck_Analyzer`. If validation fails, attempt repair with `ShapeFix_Shape`, but investigate root causes.
- Compute and report final volume using `BRepGProp` for comparison with expected values.

---

### Stage 4: Output

#### 4.1 Output objects
Write constructed solid(s) to a STEP file.

**Algorithm:**
- Create a `STEPControl_Writer`.
- Add each solid (or compound of solids) using `Transfer` with `STEPControl_AsIs` mode to preserve exact geometry.
- Write the output file.
- Optionally write BREP format (OCCT native) for debugging with `BRepTools::Write`.
- With the `--compare_step` flag:
    - Load the reference STEP file.
    - Compare volumes of the output solid vs reference solid using `BRepGProp`.
    - Compute maximum distance between output and reference surfaces using `BRepExtrema_DistShapeShape`.
    - Report pass/fail against `--surface-tolerance`.

---

## Source Code Organization

Each stage will be represented by a function which takes a ref to a configuration data structure and consumes an input data structure, and returns an output data structure (or error). That data structure is passed to the next stage:

- main:
  - let config = parse flags, load --compare_step STEP file
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
- [x] Implement tests that all stl/step file pairs in tests/ pass consistency checks and pass when fed to brepper with the --compare_step flag. Also invent a few file pairs in tests/bad that will fail with the --compare_step flag by editing the location of one or more vertices to fail surface closeness tests. Also invent some bad cases that fail the mesh validation tests in various ways. Ensure that these bad cases fail in the correct way.

### Stage 2: Surface Fitting - planes, cylinders, and spheres, oh my!
- [x] Understand stage 2.1 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.1. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [x] Implement stage 2.1. Make sure --compare_step passes for all test shapes composed only of planar surfaces.
- [x] Understand stage 2.2 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.2. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [x] Implement stage 2.2: Deduce cylindrical hypotheses.
- [x] Understand stage 2.3 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.3. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [x] Implement stage 2.3: Deduce spherical hypotheses.
- [x] Implement stage 2.6: Select surfaces for reconstruction using per-face priority rule.
- [x] Implement the --angular-tolerance flag for cylindrical and spherical surface hypothesis generation.

### Stage 3: Surface Reconstruction
- [ ] Implement stage 3.1: Create OCCT surface objects and build adjacency graph from mesh connectivity.
- [ ] Implement stage 3.2: Detect tangency relationships between adjacent surfaces (initially mark all as non-tangent).
- [ ] Implement stage 3.3: Compute edge curves via surface-surface intersection, trim to vertices, derive pcurves.
- [ ] Implement stage 3.4: Create OCCT faces from surfaces bounded by edge wires.
- [ ] Implement stage 3.5: Stitch faces into shells using BRepBuilderAPI_Sewing.
- [ ] Implement stage 3.6: Construct solids from shells, determine outer/inner shell roles.
- [ ] Revisit stage 3.2: Detect tangency relationships between adjacent surfaces. If there are any problems with models that involve tangent curves, such as tests/onshape/chamfered_cube, then imagine more tests for difficult tangency relationships including cylinder-plane and sphere-cylinder, create those tests, and ensure that tangency detection works correctly.
- [ ] Revisit stage 3.3 if tangency detection was added.

### Stage 4: Output
- [ ] Implement stage 4.1: Write solids to STEP file, validate against reference with --compare_step.

### Stage 2 Extensions: Additional Surface Types
- [ ] Understand stage 2.4 (ruled surfaces) and imagine challenging test shapes. Create test models in tests/ccad/.
- [ ] Implement stage 2.4: Deduce ruled surface hypotheses.
- [ ] Understand stage 2.5 (NURBS/B-spline surfaces) and imagine challenging test shapes. Create test models in tests/ccad/.
- [ ] Implement stage 2.5: Deduce NURBS hypotheses.
- [ ] Imagine test shapes requiring extended surface selection (mixed analytic + freeform). Create test models in tests/ccad/.
- [ ] Revisit stage 2.6: Extend surface selection to handle ruled and NURBS surfaces alongside analytic surfaces.

---

## Testing Strategy

The main testing strategy is to process a set of example stl/step pairs, and use the --compare_step flag to ensure that the pipeline is working correctly for each. 

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
|| Simple Sphere (ccad) | 974 triangles | 1 sphere | ✓ Stage 2.3 |
|| Hemisphere (ccad) | 518 triangles | 1 sphere + 1 plane | ✓ Stage 2.3 |
|| Spherical Pocket (ccad) | 486 triangles | 6 planes + 1 concave sphere | ✓ Stage 2.3 |
|| Ball on Cylinder (ccad) | 764 triangles | 1 sphere + 1 cylinder + 1 plane | ✓ Stage 2.3 |
|| Sphere | Tessellated sphere | 1 sphere | ✓ Stage 2.3 |
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
