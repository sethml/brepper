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
- With the `--compare` flag: check that all vertices are within `--vertex-tolerance` of the STEP shape, and that all face centroids are within `--surface-tolerance`. Uses `BRepExtrema_DistShapeShape` for bounded distance (respects face trimming). Reports the worst vertex and worst centroid with their distances.

Validation checks, in priority order (first failure is reported):
1. **Degenerate faces**: faces with <3 vertices or zero-area normals.
2. **Non-manifold edges**: edges shared by >2 faces.
3. **Inconsistent orientation**: adjacent faces that traverse a shared edge in the same direction (indicating flipped normals within a shell).
4. **Inverted normals**: for each closed shell, compute the signed volume via the divergence theorem (sum of signed tetrahedra formed by each face and the origin: V = (1/6) Σ v0·(v1×v2)). A negative signed volume means all face normals point inward (consistently inverted winding). This catches meshes where orientation is internally consistent (check 3 passes) but globally inverted.
5. **Self-intersections**: intentionally deferred to a later stage where richer geometric predicates are available and O(N²) triangle checks can be avoided.

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

Stages 2.1–2.5 fit analytic surface hypotheses (planes, cylinders, spheres, cones, tori) to groups of mesh faces using BFS region growing with per-surface seeding strategies and fitting math. Full algorithm descriptions are in [SURFACE_FITTING.md](SURFACE_FITTING.md).

#### 2.1 Deduce planar hypotheses
BFS region growing from seed faces, with area-weighted normal averaging and vertex-tolerance acceptance. See [SURFACE_FITTING.md](SURFACE_FITTING.md#21-deduce-planar-hypotheses).

#### 2.2 Deduce cylindrical hypotheses
Multi-seed evaluation with 3-face seeds, BFS expansion with cylinder fit qualification, angular coverage validation. Fitting via normal covariance axis estimation + 2D circle fit + LM refinement. See [SURFACE_FITTING.md](SURFACE_FITTING.md#22-deduce-cylindrical-hypotheses).

#### 2.3 Deduce spherical hypotheses
Vertex-neighborhood seeding, BFS with solid-angle coverage validation. Algebraic least-squares sphere fitting. See [SURFACE_FITTING.md](SURFACE_FITTING.md#23-deduce-spherical-hypotheses).

#### 2.4 Deduce conical hypotheses
Apex-vertex seeding, (h,r) profile fitting and quadric fitting initialization, LM refinement. See [SURFACE_FITTING.md](SURFACE_FITTING.md#24-deduce-conical-hypotheses).

#### 2.5 Deduce toroidal hypotheses
Face-based seeding, medial axis / tube center method with SOR axis fallback, multi-axis LM rescue. See [SURFACE_FITTING.md](SURFACE_FITTING.md#25-deduce-toroidal-hypotheses).

#### 2.6 Select surfaces to use for reconstruction
Assign each mesh face to exactly one "selected surface" for reconstruction. Since stages 2.1–2.5 produce hypotheses that may overlap (a face can belong to one planar, one cylindrical, one spherical, one conical, and one toroidal hypothesis simultaneously), surface selection must resolve these overlaps.


**Greedy area-based selection:**

The algorithm selects hypotheses greedily by total mesh face area, largest first. This naturally resolves conflicts: real surface hypotheses cover many faces with large total area, while bogus hypotheses (e.g., a near-degenerate r=139mm sphere from 3 nearly-coplanar faces) cover minimal area and lose their faces to correctly-fitted surfaces.

**Algorithm:**
1. Compute the geometric area of each mesh face.
2. Build a candidate list of all hypotheses (planar with face count > 1, cylindrical, spherical, conical, toroidal). Each candidate tracks its set of member faces and their total area.

3. **Greedy loop**: while unassigned faces remain:
   a. For each candidate, compute the total area of its still-unassigned faces.
   b. Select the candidate with the largest remaining area.
   c. Assign all its still-unassigned faces to it as a selected surface.
   d. Remove those faces from all other candidates' pools.
4. Assign any remaining unassigned faces to their single-face planar hypothesis (fallback for faces on surfaces not yet fitted: ruled, freeform). Stage 2.8 may replace some of these with NURBS surfaces.

5. Build a `Vec<SelectedSurface>`, one per selected hypothesis, containing:
   - `surface_type: SurfaceType` — enum: Planar, Cylindrical, Spherical, Conical, Toroidal (later: BSpline).

   - `hypothesis_index: usize` — index into the appropriate hypothesis vector in Stage2Output.
   - `faces: Vec<usize>` — mesh face indices assigned to this surface.
   - `vertices: Vec<usize>` — mesh vertex indices on this surface (union of face vertices).
6. Validate: every face is assigned to exactly one selected surface.
7. Print a summary: count of selected surfaces by type, total faces covered.

**Why this works better than a fixed type-priority rule:** A fixed priority (e.g., spherical > cylindrical > multi-face planar) fails when a bogus hypothesis of a high-priority type (3-face r=139mm sphere) competes with a correct hypothesis of a lower-priority type (35-face r=2mm cylinder). Area-based greedy selection naturally picks the correct hypothesis because it covers vastly more area.

**Nerfed tests:**
- (none currently — `onshape_pipe_elbow_stage26` now passes `--compare`, see note below.)
- `onshape_pipe_elbow_stage26_compare`: Passes `--compare` with 290 cylindrical + 0 toroidal selections. The pipe elbow's torus surface is tessellated into many narrow strips whose normals span < 180°, so stage 2.5 torus fitting doesn't fire (normals look locally cylindrical). This is expected — the torus is decomposed into cylindrical patches that individually pass `--compare`.

- With the `--compare` flag, for each selected surface:
    - Compute face centroids and project them onto the selected surface.
    - Measure distance from each projected centroid to the nearest surface in the reference STEP file.
    - For all hypotheses, report an error if distance exceeds `--surface-tolerance`.

#### 2.7 Surface refitting (optional)

After stage 2.6 selects surfaces, some border faces may have been absorbed by a hypothesis that fits within tolerance but isn't the best fit. This step refines assignments and re-fits surfaces to improve accuracy.

**Algorithm:**
1. For each mesh face, check whether it belongs to more than one hypothesis (i.e., it has valid assignments in multiple hypothesis types — planar, cylindrical, spherical, conical, toroidal).

2. For each such face, compute the vertex-to-surface distance for each candidate surface (using the selected surfaces from 2.6, not the raw hypotheses).
3. Reassign the face to the selected surface with the smallest maximum vertex distance.
4. After all reassignments, re-fit each affected surface from its updated face/vertex set to remove any bias introduced by the extra vertices.

This step is expected to be most useful at boundaries between surfaces of different types (e.g., where a cylinder meets a plane), where BFS may have greedily absorbed a border face that technically fits within tolerance but belongs more naturally to the adjacent surface. It may be unnecessary if the greedy area-based selection in 2.6 already produces clean boundaries — implement only if stage 3 reconstruction reveals boundary accuracy problems.

#### 2.8 NURBS fallback

For groups of adjacent faces that remain as single-face planar hypotheses after stages 2.1–2.5 (i.e., they don't fit any analytic surface), try to fit a NURBS or B-spline surface.

**Algorithm:**
1. Identify connected components of single-face planar hypothesis faces that are adjacent to each other in the mesh.
2. For each component with ≥ 4 faces:
   a. Compute UV parameterization of the face vertices (using discrete harmonic mapping or Floater's mean value coordinates).
   b. Use OCCT's `GeomAPI_PointsToBSplineSurface` to fit a B-spline through the parameterized points.
   c. Validate: max vertex-to-surface distance < `--surface-tolerance`.
   d. If valid, create a NURBS hypothesis; if not, leave as single-face planar.
3. Add valid NURBS hypotheses to the surface selection.

*Comment: NURBS fitting is complex: parameterization quality strongly affects fit quality, degree/knot selection matters, and thin features may need special handling. Start with OCCT’s built-in B-spline fitter which handles degree selection and parameterization internally. For “freeform” regions that are nearly planar, a plane with small tolerance is preferable to a NURBS that overfits noise.*

---

### Algorithm Reference

Algorithmic techniques (Unified SOR Framework, Gaussian Map Classification, Medial Axis / Tube Center Method, Quadric Fitting) are documented in [SURFACE_FITTING.md](SURFACE_FITTING.md#algorithm-reference).

---

### Stage 3: Surface Reconstruction

Build OCCT B-Rep topology (faces, edges, vertices) from the selected surfaces. This is the most complex stage — surface intersections in 3D are numerically delicate, and small perturbations can cause large changes in intersection curves or cause them to disappear.

**Core data structures:**

A vector of `FaceDescriptor` structs, one per selected surface:
- `selected_surface_idx: usize` — Index into the selected surfaces from stage 2.6.
- `surface: OwnedPtr<HandleGeomSurface>` — Infinite OCCT surface, type-erased (created in 3.1).
- `adjacent_faces: Vec<usize>` — Indices of adjacent FaceDescriptors, ordered topologically around this face's boundary. A face may appear multiple times if it shares multiple edges with this face (e.g., a cylindrical face adjacent to itself via a seam edge).
- `edge_indices: Vec<usize>` — ReconEdge indices. `edge_indices[i]` is the edge between this face and `adjacent_faces[i]`.
- `vertex_indices: Vec<usize>` — BRepVertex indices at corners. `vertex_indices[i]` is the vertex between `edge_indices[i]` and `edge_indices[(i+1) % N]`.
- (TODO) `occt_face: Option<TopoDS_Face>` — Populated in 3.4.

A vector of `ReconEdge` structs:
- `face_indices: [usize; 2]` — Indices of the two adjacent FaceDescriptors.
- `vertex_indices: [usize; 2]` — Indices of the two BRepVertices at each end. For closed-loop edges (no vertices), both are `usize::MAX`.
- `curve_3d: Option<OwnedPtr<HandleGeomCurve>>` — 3D intersection curve, trimmed to vertex endpoints (populated in 3.3).
- (TODO) `pcurves: [Handle<Geom2d_Curve>; 2]` — 2D parametric curve on each adjacent face's surface.
- `tangent: bool` — Whether the adjacent surfaces are tangent along this edge (detected in 3.2).
- `mesh_boundary_vertices: Vec<usize>` — Mesh vertex indices along this boundary, ordered along the boundary.

A vector of `BRepVertex` structs:
- `point: [f64; 3]` — 3D position (currently raw mesh vertex position; future: project onto adjacent surfaces and average).
- `adjacent_faces: Vec<usize>` — Indices of FaceDescriptors meeting at this vertex, in topological order.
- `adjacent_edges: Vec<usize>` — Indices of ReconEdges meeting at this vertex, in topological order.
- (TODO) `uv_coords: Vec<(f64, f64)>` — UV coordinates of this vertex on each adjacent face's surface.
- (TODO) `tolerance: f64` — Maximum deviation from any adjacent surface.

#### 3.1 Create OCCT surface objects and build adjacency graph

**Surface creation** (`create_occt_surface`):
- For each selected surface, create the corresponding infinite OCCT surface:
    - Planar → `geom::Plane::new_pnt_dir` with `gp::Pnt` at centroid and `gp::Dir` along normal. Wrapped via `to_handle().to_handle_surface()`.
    - Cylindrical → `geom::CylindricalSurface::new_ax3_real` with `gp::Ax3` from axis origin/direction and radius. Wrapped via `to_handle().to_handle_surface()`.
    - Spherical → `geom::SphericalSurface::new_ax3_real` with `gp::Ax3` from center and Z-axis direction. Wrapped via `to_handle().to_handle_surface()`.
    - Conical → `geom::ConicalSurface::new_ax3_real2` with `gp::Ax3` from apex/direction and half-angle + radius at reference height. Wrapped via `to_handle().to_handle_surface()`. The `Ax3` origin is at the cone apex, Z along axis.
    - Toroidal → `geom::ToroidalSurface::new_ax3_real2` with `gp::Ax3` from center/axis direction, major radius $R$, and minor radius $r$. Wrapped via `to_handle().to_handle_surface()`. The `Ax3` origin is at the torus center, Z along rotation axis.
    - (Future: BSpline → `Geom_BSplineSurface` via `GeomAPI_PointsToBSplineSurface`.)

- For consistency with OCCT conventions: planes use gp_Ax3 with Z as normal; cylinders have Z along axis; spheres have poles at Z extremes; cones have Z along axis from apex; tori have Z along rotation axis. This affects pcurve computation in 3.3.


**Adjacency graph construction:**
1. Build face-to-surface map: for each mesh face, record which selected surface it belongs to (`build_face_to_surface_map`).
2. Collect boundary edges: walk all mesh edges; when an edge separates faces belonging to different surfaces, record it as an `UndirectedEdge` keyed by `SurfacePair` (`collect_boundary_edges`).
3. Chain boundary edges: for each surface pair, chain connected boundary edges into ordered vertex sequences (`chain_boundary_edges`). Each chain becomes a `ReconEdge` with the mesh boundary vertices in order.
4. Find corner vertices: mesh vertices where 3+ selected surfaces meet become `BRepVertex` entries (`find_corner_vertices`). Their 3D position is taken from the mesh vertex.
5. Assign vertex indices to edges: for each `ReconEdge`, check if its first/last mesh boundary vertex is a corner vertex and set `vertex_indices` accordingly (or `usize::MAX` for closed-loop edges).
6. Build face descriptors: for each selected surface, create a `FaceDescriptor` with the OCCT surface handle. Walk the boundary edges, ordering them topologically around the face boundary (`walk_edges_around_face`). This uses a rotation-order algorithm: at each vertex, pick the next edge that rotates counter-clockwise around the face. This populates `adjacent_faces`, `edge_indices`, and `vertex_indices` for each face.
7. Build vertex adjacency: for each `BRepVertex`, walk the edges meeting at that vertex in topological order, recording `adjacent_faces` and `adjacent_edges`.
8. Validate topology: check Euler's formula V-E+F=2 for genus-0 models.

- With the `--compare` flag:
    - For each `ReconEdge`, sample the mesh boundary vertices along the edge. For each sample vertex, compute the distance to the nearest edge/vertex in the reference STEP file using `BRepExtrema_DistShapeShape`. Report an error if any sample exceeds `--vertex-tolerance`.
    - For each `BRepVertex`, compute the distance from the vertex's 3D position to the nearest vertex in the reference STEP file. Report an error if any exceeds `--vertex-tolerance`.
    - These checks validate that the adjacency graph's topological features (edges where surfaces meet, corners where 3+ surfaces meet) correspond to real features in the reference model.

#### 3.2 Detect tangency relationships
For each ReconEdge, determine whether the two adjacent surfaces are tangent along the shared boundary. This matters because `GeomAPI_IntSS` can fail or produce degenerate results for tangent/near-tangent surfaces, requiring special handling in edge curve computation (3.3).

**Detection algorithm (implemented):**
- For each ReconEdge, sample mesh boundary vertices along the shared boundary (every 5th vertex, minimum 3 samples, always including first and last).
- At each sample point, compute the outward-facing surface normal of both adjacent surfaces analytically:
    - Planar: the hypothesis normal (constant).
    - Cylindrical: normalized radial direction from axis to point, negated if concave.
    - Spherical: normalized direction from center to point, negated if concave.
    - Conical: radial direction from axis at the point's axial position, rotated by the cone's half-angle away from the axis, negated if concave. Specifically: for point $p$ with axial component $h = (p-A)\cdot\hat{a}$, the normal is the unit vector in the plane of (radial, axis) at angle $(90° - \theta)$ from radial.
    - Toroidal: direction from the nearest tube center point to the surface point, negated if concave. The tube center is $c = C + R \cdot \hat{r}$ where $\hat{r}$ is the radial unit vector from the torus axis to the point.

- Compute the dot product of the two normals. If all sampled pairs agree within 2° (dot > cos(2°) ≈ 0.9994), mark the edge as tangent.
- Points too close to a cylinder axis or sphere center (degenerate) cause the edge to be marked non-tangent.

**Decision: do not modify surfaces to enforce tangency.** Modifying analytic surfaces would change the geometry. Instead, tangent edges get special handling in edge curve computation (3.3): construct the edge curve directly rather than relying on surface-surface intersection.

Tangent edges arise from fillets and blends: part_rounded_cube_10_r2 has 8 plane-cylinder tangent edges, and rounded_cube_10_r2 has 46 tangent edges (plane-cylinder + sphere-cylinder). Toroidal fillets produce tangent edges at both boundaries (torus-plane, torus-cylinder, torus-torus).


- With the `--compare` flag:
    - For each ReconEdge marked as tangent, verify that the corresponding edge in the reference STEP file is also tangent (i.e., the two adjacent STEP surfaces have matching normals along the STEP edge). Report a warning if tangency is detected in the mesh but not in the STEP reference, or vice versa.
    - This is a consistency check, not a geometry check — stage 3.1's --compare already validates that mesh edge boundary vertices lie on STEP edges. Stage 3.2 adds validation that the tangency classification itself is correct.

#### 3.3 Create edge curves
For each ReconEdge, compute the 3D intersection curve between the two adjacent surfaces, trim it to the vertex endpoints, and store in `edge.curve_3d`.

**Intersection computation (implemented):**
- For non-tangent edges: use `GeomAPI_IntSS` to compute intersection curves with tolerance `vertex_tolerance_mm`. IntSS may return multiple curves (e.g., a plane cutting through a torus); `select_closest_curve()` picks the one closest to the mesh boundary vertices by sampling up to 10 evenly-spaced boundary vertices, projecting each onto each candidate curve via `GeomAPI_ProjectPointOnCurve`, and summing distances.
- For tangent edges: construct the curve analytically rather than using `GeomAPI_IntSS` (which fails or produces degenerate results for tangent/near-tangent surfaces). Dispatches by surface pair type:
  - **Plane-cylinder**: `Geom_Line` parallel to the cylinder axis at the tangent point, computed from the cylinder axis, radius, and the plane normal component perpendicular to the axis.
  - **Cylinder-cylinder**: `Geom_Line` along the first cylinder's axis direction.
  - **Sphere-cylinder**: `Geom_Circle` (great circle arc on the sphere). Center and radius come from the sphere hypothesis. The circle plane normal is `normalize(cross(v0 - center, v1 - center))`. Arc selection samples mesh boundary vertices to determine forward vs. reverse arc on the periodic circle.
  - **Plane-cone** (TODO): `Geom_Line` along the cone's generator at the tangent azimuth, from the tangent point toward (or away from) the apex.
  - **Cylinder-cone** (TODO): `Geom_Circle` at the axial position where the cone radius equals the cylinder radius (when coaxial), or `IntSS`-based (non-coaxial).
  - **Plane-torus**: `Geom_Circle` perpendicular to the torus axis. Circle center and radius computed from mesh boundary vertices projected onto the torus axis. Handles both closed loops and arcs.
  - **Cylinder-torus**: `Geom_Circle` perpendicular to the torus axis. Same vertex-averaging algorithm as plane-torus.
  - **Torus-torus**: `Geom_Circle` perpendicular to the torus axis. Same vertex-averaging algorithm.
  - **Sphere-torus** (TODO): `Geom_Circle` — analogous to sphere-cylinder but the junction is a circle on the sphere.
  - **Fallback**: Construct a line from the two vertex endpoints.

  All tangent edges are trimmed to vertex endpoints using the same parameter projection as non-tangent edges. A `validate_tangent_curve` step checks that boundary vertices lie within `surface_tolerance_mm` of the constructed curve.

**Trimming to vertex endpoints (implemented):**
- For each ReconVertex at an edge endpoint, project the vertex's 3D position onto the intersection curve using `GeomAPI_ProjectPointOnCurve` to get the curve parameter value.
- Trim the curve to the parameter range `[t_start, t_end]` using `Geom_TrimmedCurve`.
- For closed-loop edges (vertex_indices == usize::MAX), use the full intersection curve parameter range.
- **Partial arc reconstruction**: When `IntSS` returns a partial arc for a closed-loop edge (e.g., a semicircle when a sphere-plane intersection passes through the sphere's UV poles), reconstruct the full circle from 3 sampled points on the arc using the circumscribed circle formula. This fixes cases where the intersection plane contains the sphere axis, causing `IntSS` to split the circle into two separate arcs in UV space.
- For non-periodic curves, ensure `t_start < t_end` (swap if needed).
- The resulting `TrimmedCurve` is upcast to `HandleGeomCurve` via `.to_handle().to_handle_curve()` and stored in `edge.curve_3d`.

**Pcurve computation (deferred to stage 3.4):**
- Pcurves will be derived when creating `TopoDS_Edge` objects.
- Options: analytic pcurves for simple cases (line on plane, circle on cylinder), or `ShapeFix_Wire::FixEdgeCurves` for general cases.

**Seam edges (future):**
- A face adjacent to itself (e.g., a full cylinder wrapping around) requires a seam edge at a fixed U or V parameter. Create the seam as an iso-parametric curve on the surface.
- Conical faces: full-revolution cones require a seam edge along one generator line (iso-U curve from apex to base).
- Toroidal faces: full-revolution tori (360° around major axis) require a seam at one meridian (iso-U); full-tube tori (360° around minor circle) also need a seam at one parallel (iso-V). A complete torus needs both.


**Vertex position consistency (future):**
- OCCT requires all edges at a vertex to share the same `TopoDS_Vertex` (same 3D point). When edges computed independently disagree slightly about vertex position, use the averaged position from the ReconVertex and set the vertex tolerance to the maximum deviation.

- With the `--compare` flag:
    - For each ReconEdge with a computed `curve_3d`, sample points along the trimmed curve (e.g., 10 evenly-spaced parameter values between start and end). For each sample point, compute the distance to the nearest edge in the reference STEP file using `BRepExtrema_DistShapeShape`. Report an error if any sample exceeds `--vertex-tolerance`.
    - Additionally, verify that the curve endpoints coincide with BRepVertex positions within `--vertex-tolerance`.
    - This validates that the computed intersection curves not only connect the right vertices but also follow the correct geometric path between them.

#### 3.4 Create OCCT faces
For each ReconFace, construct a `TopoDS_Face` from the surface and bounding wires.

**Algorithm (as implemented):**

*Stage 1: Create TopoDS_Edge for each ReconEdge.*
- For each `ReconEdge`, create a `BRepBuilderAPI_MakeEdge` using the trimmed 3D curve and shared `TopoDS_Vertex` endpoints from stage 3.3. For periodic curves (circles from cylinder/sphere intersections), vertex order is matched to the curve parameterization direction: V1 corresponds to `first_parameter()` and V2 to `last_parameter()`. This is critical because `MakeEdge` strips `TrimmedCurve` wrappers and uses `ElCLib::AdjustPeriodic` on the base periodic curve, which always selects the forward (CCW) arc from V1 to V2.
- **Adaptive vertex tolerance**: Before each `MakeEdge`, the distance between each vertex position and its curve endpoint is computed. The OCCT vertex tolerance is updated to at least that distance (with 1% margin), accommodating imprecisely-fitted surface intersection curves that don't pass exactly through mesh vertex positions. This replaces the earlier fixed `vertex_tolerance_mm` approach.

*Stage 2: Create OCCT faces.*
- **Planar faces** use wire-based construction:
  - Group the face's edges into wire loops based on `BRepVertex` connectivity. Edges sharing a vertex index are grouped together; closed-loop edges (vertex_indices both `usize::MAX`) each form a separate wire.
  - Build `BRepBuilderAPI_MakeWire` for each group by sequentially adding `TopoDS_Edge` objects.
  - Identify the outer wire as the group with the most edges (heuristic).
  - Create `BRepBuilderAPI_MakeFace` with the surface handle and outer wire, then add inner wires as holes.

- **Cylindrical, spherical, conical, and toroidal faces** use three approaches:

  - **Full revolution** (all boundary edges are closed loops): UV-bounds construction via `MakeFace::new_handlegeomsurface_real5(surface, umin, umax, vmin, vmax, tol)`. UV bounds are computed from edge projections using a circular gap algorithm that handles periodicity. This automatically creates seam edges for proper B-Rep topology.
  - **Spherical faces with closed-loop edges near poles** (e.g., pill/capsule hemispheres or dome hemispheres cut by a meridional plane): The boundary circles may pass through or near the sphere's UV singularities (poles at V=±π/2 where U is undefined). When any boundary circle passes within 45° of either pole (detected via `ProjectPointOnCurve` + chord-length formula), the sphere surface is recreated with its Z-axis aligned to a reorientation axis. The axis is determined from: (a) the adjacent cylinder's axis direction for tangent edges, or (b) the boundary circle's plane normal for non-tangent edges (computed from 3 sampled points on the curve). This is safe because for all-closed-loop boundaries, the reorientation axis is unique. The hemisphere is selected by projecting a mesh face centroid onto the oriented surface and checking the V-sign. UV-bounds [0,2π] × [0,π/2] or [0,2π] × [-π/2,0] then work correctly since the boundary circle becomes an iso-V curve.
  - **Spherical faces with pole vertices** (e.g., rounded cube corners where cylinders meet at a sphere pole): When a spherical face's UV bounds include V=±π/2, vertex U values at the pole are undefined (singularity). The `compute_uv_bounds_from_edges` function excludes pole U values (within 0.01 radians of ±π/2) from the circular gap algorithm to avoid corrupted U spans. The face is then constructed via UV-bounds `MakeFace::new_handlegeomsurface_real5`, producing a triangular-topology face (3 edges, 2 meridians + 1 parallel; pole collapsed to a vertex). These faces have a harmless BRepCheck edge failure (tolerance ~0) on the first edge, resolved by sewing.
  - **Partial revolution** (open edges with vertices): Wire-based construction with pre-set pcurves. Before building the wire, `BRep_Builder::update_edge` sets `Geom2d_Line` pcurves on each IntSS edge for the periodic surface, mapping the edge's 3D parameter range to (u(t), v(t)) in surface UV space. This ensures MakeFace selects the correct arc and shares edges with adjacent planar faces for sewing.
  - **Conical faces**: Same approaches as cylinder: full revolution via UV-bounds (with circular gap for U and cone apex V-extension), partial via wire+pcurves. The cone apex vertex is handled as a UV singularity (V extended to 0 for non-truncated cones). All test models pass: simple_cone, cone_cylinder, nosecone, block_with_conical_hole.
  - **Toroidal faces**: UV-bounds construction for complete torus patches. The torus UV space is doubly periodic: U ∈ [0,2π) around the major axis, V ∈ [0,2π) around the tube. For V-bounds, the circular gap algorithm handles periodicity (same approach as U-bounds for partial revolution). Multiple sample points per edge curve (3 instead of 1) and mesh face centroid projections provide robust V disambiguation. Fillet patches typically span a small V arc (e.g., π/2 for a 90° fillet).

  - Full spheres (no edges) use `MakeFace::new_handlegeomsurface_real(surface, tolerance)` with natural bounds.

*Stage 3: Validate and compare.*
- Validate each face with `BRepCheck_Analyzer`. For invalid faces, run detailed `BRepCheck_Face` diagnostics (intersect_wires, classify_wires, orientation_of_wires, is_unorientable) and sub-shape analysis (edge tolerance, vertex checks). Faces where the only failure is edges with near-zero tolerance (< `vertex_tolerance_mm`) are treated as harmless (occurs on UV-bounds sphere pole faces, resolved by sewing) and not counted as failures.
- Wire orientation fix: after adding inner wires (holes) to a face, apply `ShapeFix_Face::fix_orientation()` to ensure inner wires are oriented opposite to the outer wire. This eliminates BRepCheck "bad orientation of sub-shape" warnings that previously appeared on all planar faces with holes.
- With `--compare`: for each face, sample a representative mesh face centroid from the surface's mesh faces and compute distance to the reference STEP shape using `BRepExtrema_DistShapeShape`. Report face count comparison. Error if max distance exceeds `--surface-tolerance`.

#### 3.5 Construct shells
Group faces into shells using `BRepBuilderAPI_Sewing`.

**Algorithm (as implemented):**
- Create a `BRepBuilderAPI_Sewing` instance with tolerance set to `vertex_tolerance`.
- Add all `TopoDS_Face` objects to the sewing operation.
- After `Perform()`, extract the sewed shape and iterate over its sub-shapes:
  - `Shell` sub-shapes are used directly.
  - Individual `Face` sub-shapes (e.g., a single-face model like `simple_sphere`) are wrapped in a new `TopoDS_Shell` via `BRep_Builder`.
  - `Solid`, `Compound`, and `Compsolid` sub-shapes are explored recursively for `Shell` children.
- Apply `ShapeFix_Shell::Perform()` (which calls `ShapeFix_Face::Perform()` on each face to re-add pcurves that sewing may have discarded, and `FixFaceOrientation` for edge consistency). Face orientation is NOT post-processed — `FixFaceOrientation` handles edge consistency and `SolidFromShell` (stage 3.6) handles global orientation via `PerformInfinitePoint`. Investigation showed that flipping individual face orientations breaks edge consistency and causes volume computation errors.
- Validate with `ShapeAnalysis_Shell::CheckOrientedShells` and report inconsistencies.
- Report sewing statistics: free edges, multiple edges, contiguous edges.

- With the `--compare` flag:
    - Count shells in the reference STEP file by exploring for `Shell` sub-shapes.
    - Report an error if the sewing produced free edges (edges not shared by any pair of faces), as this indicates incomplete stitching.
    - Report shell count comparison and orientation check results.

#### 3.6 Construct solids
Convert closed shells into `TopoDS_Solid` objects.

**Algorithm:**
- For each shell, use `ShapeFix_Solid::SolidFromShell` which handles orientation automatically.
- Apply `ShapeFix_Shape::Perform()` for comprehensive face/wire/edge/shell/solid fixes.
- Call `BRepLib::SameParameter` with `forced=false` and tight tolerance (`vertex_tolerance_mm`) to ensure pcurve/3D curve agreement without corrupting already-correct pcurves. Using `forced=true` or a large tolerance (e.g. 1.0mm) recomputes all pcurves and can corrupt faces with closely-spaced tangent edges (e.g. sphere faces at sphere-cylinder junctions), causing BRepCheck self-intersection failures.
- Call `BRepLib::UpdateTolerances` and `BRepLib::OrientClosedSolid`.
- Validate with `BRepCheck_Analyzer`. Log warning if validation fails (may indicate upstream face issues from stage 3.4).
- Compute and report volume using `BRepGProp::VolumeProperties`.
- Multi-shell solids (outer shell containing voids) not yet implemented — currently each shell produces one solid.

- With the `--compare` flag:
    - Count solids in the reference STEP file.
    - Compute volumes of both constructed and reference solids. Match by sorted volume. Report warning if relative volume difference exceeds 1%. All current test models achieve volume agreement within ~1e-7 relative difference.
    - Compute `BRepExtrema_DistShapeShape` between each constructed solid and the reference STEP shape. Report error if distance exceeds `--surface-tolerance`. This is the tightest geometric accuracy check.

---

### Stage 4: Output

#### 4.1 Output objects
Write constructed solid(s) to a STEP file.

**Algorithm (implemented):**
- Build a `topo_ds::Compound` containing all solids from stage 3.6.
- Set STEP metadata: AP214 schema, product name "brepper" via `Interface_Static::set_c_val`.
- Create a `STEPControl_Writer`, transfer the compound with `StepModelType::Asis` mode to preserve exact geometry.
- Write the output file. Report error if transfer or write returns non-Done status.
- With the `--compare` flag:
    - Re-read the written STEP file via `STEPControl_Reader` to validate the round-trip.
    - Compare volumes of the output shape vs reference shape using `BRepGProp::VolumeProperties`. Warn on >1% relative difference.
    - Compute maximum distance between output and reference using `BRepExtrema_DistShapeShape`. Hard error if exceeds `--surface-tolerance`.

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

## Testing Strategy

The main testing strategy is to process a set of example stl/step pairs, and use the --compare flag to ensure that the pipeline is working correctly for each. 

1. **Unit Tests**: Individual components (readers, fitters, converters)
2. **Integration Tests**: Full pipeline on known geometries
3. **Regression Tests**: Compare output STEP to reference
4. **Validation**: Round-trip test (STEP → mesh → STEP)

### Test Cases

Each STL test file was run with both debug and release builds of brepper:

```
# Debug build:
./target/debug/brepper <input.stl> -o /dev/null --compare <input.step>
# Release build:
./target/release/brepper <input.stl> -o /dev/null --compare <input.step>
# Special case (rounded_cube_10_r2_coarse): add --angular-tolerance 20
```

Tests last run: 2026-03-19. Files without a corresponding .step have no --compare flag.

#### Bad inputs (expected failures)

| Test File | Triangles | Debug | Release | Notes |
|-----------|-----------|-------|---------|-------|
| tests/bad/cube_on_plane.stl | — | FAIL (0.04s) | FAIL (0.42s) | Stage 1.2 compare: vertex 9mm from STEP |
| tests/bad/cube_shifted.stl | — | FAIL (0.04s) | FAIL (0.04s) | Stage 1.2 compare: vertex 6.9mm from STEP |
| tests/bad/degenerate_face.stl | — | FAIL (0.03s) | FAIL (0.03s) | Validation: 1 degenerate face |
| tests/bad/inconsistent_winding.stl | — | FAIL (0.03s) | FAIL (0.03s) | Validation: 3 flipped edges |
| tests/bad/non_manifold_edge.stl | — | FAIL (0.03s) | FAIL (0.03s) | Validation: 1 non-manifold edge |

#### ccad/generated

| Test File | Triangles | Debug | Release | Notes |
|-----------|-----------|-------|---------|-------|
| ball_on_cylinder.stl | 765 | ✓ (2.04s) | ✓ (1.56s) | |
| block_corner_hole.stl | 44 | ✓ (0.08s) | ✓ (0.08s) | |
| block_with_conical_hole.stl | 312 | FAIL (0.64s) | FAIL (0.53s) | Stage 3.1 compare: edge 5.3mm from STEP |
| block_with_hole.stl | 132 | ✓ (0.45s) | ✓ (0.37s) | |
| cone_cylinder.stl | 314 | ✓ (0.26s) | ✓ (0.23s) | |
| cube.stl | 12 | ✓ (0.05s) | ✓ (0.05s) | |
| filleted_cylinder.stl | 1528 | FAIL (41.02s) | FAIL (26.41s) | Stage 3.5 compare: shell inf from STEP |
| filleted_hole_block.stl | 1484 | FAIL (30.96s) | FAIL (20.58s) | Stage 3.1 compare: edge 0.39mm from STEP |
| filleted_pipe.stl | 2576 | FAIL (39.11s) | FAIL (27.96s) | Stage 3.1 compare: edge 0.59mm from STEP |
| hemisphere.stl | 519 | ✓ (0.94s) | ✓ (0.73s) | |
| nosecone.stl | 713 | ✓ (0.99s) | ✓ (0.77s) | |
| pipe.stl | 244 | ✓ (0.11s) | ✓ (0.11s) | |
| simple_cone.stl | 348 | ✓ (0.28s) | ✓ (0.23s) | |
| simple_cylinder.stl | 124 | ✓ (0.06s) | ✓ (0.05s) | |
| simple_sphere.stl | 976 | ✓ (58.51s) | ✓ (41.04s) | |
| spherical_pocket.stl | 487 | ✓ (4.94s) | ✓ (3.66s) | |
| staircase.stl | 48 | ✓ (0.14s) | ✓ (0.14s) | |
| stepped_cylinder.stl | 240 | ✓ (0.11s) | ✓ (0.11s) | |
| t_shape.stl | 28 | ✓ (0.09s) | ✓ (0.09s) | |
| two_holes.stl | 252 | ✓ (0.50s) | ✓ (0.42s) | |
| wedge.stl | 12 | ✓ (0.05s) | ✓ (0.05s) | |

#### fusion

| Test File | Triangles | Debug | Release | Notes |
|-----------|-----------|-------|---------|-------|
| plate_with_hole_100x50_high.stl | 284 | ✓ (0.62s) | ✓ (0.52s) | |
| plate_with_hole_100x50_low.stl | 112 | ✓ (0.25s) | ✓ (0.22s) | |
| plate_with_hole_100x50_medium.stl | 168 | ✓ (0.35s) | ✓ (0.31s) | |

#### manual

| Test File | Triangles | Debug | Release | Notes |
|-----------|-----------|-------|---------|-------|
| cube.stl | 12 | ✓ (0.05s) | ✓ (0.05s) | |

#### onshape

| Test File | Triangles | Debug | Release | Notes |
|-----------|-----------|-------|---------|-------|
| chamfered_cube_10_c1_medium.ascii.stl | 44 | FAIL (0.06s) | FAIL (0.06s) | Stage 1.2 compare: vertex 4.8e-5mm from STEP (ASCII precision) |
| chamfered_cube_10_c1_medium.stl | 44 | ✓ (0.17s) | ✓ (0.17s) | |
| cone_15x20_medium.stl | 116 | ✓ (0.10s) | ✓ (0.09s) | |
| cone_tangent_planes_medium.stl | 376 | FAIL (3.99s) | FAIL (2.93s) | Stage 3.1 compare: vertex 0.16mm from STEP |
| cylinder_10_r2_fine.stl | 576 | ✓ (0.17s) | ✓ (0.16s) | |
| cylinder_10x30_medium.stl | 232 | ✓ (0.08s) | ✓ (0.08s) | |
| dome_hemisphere_20_fine.stl | 10368 | ✓ (8.33s) | ✓ (7.01s) | |
| l_bracket_simple_medium.stl | 492 | ✓ (0.90s) | ✓ (0.83s) | |
| nosecone_medium.stl | 1572 | FAIL (117.51s) | FAIL (74.87s) | Stage 3.1 compare: edge 0.18mm from STEP |
| part_rounded_cube_10_r2_coarse.stl | 140 | ✓ (0.21s) | ✓ (0.20s) | |
| part_rounded_cube_10_r2_fine.stl | 588 | ✓ (0.71s) | ✓ (0.66s) | |
| pill_coarse.stl | 752 | ✓ (25.68s) | ✓ (16.52s) | |
| pill_fine.stl | 20736 | ✓ (106.29s) | ✓ (74.32s) | |
| pipe_elbow_10_fine.stl | 11232 | FAIL (26.27s) | FAIL (21.78s) | Stage 3.1 compare: edge 7.9mm from STEP |
| pipe_elbow_10_medium.stl | 2320 | FAIL (7.97s) | FAIL (5.85s) | Stage 3.1 adjacency error |
| plate_with_hole_100x50_coarse.stl | 136 | ✓ (0.30s) | ✓ (0.26s) | |
| rounded_cube_10_r2_coarse.stl | 820 | ✓ (3.88s) | ✓ (2.99s) | Uses --angular-tolerance 20 |
| rounded_cube_10_r2_fine.stl | 21466 | ✓ (35.62s) | ✓ (33.19s) | |
| rounded_cube_10_r2_medium.stl | 3548 | ✓ (32.52s) | ✓ (23.44s) | |
| sphere_25_fine.stl | 20448 | ✓ (82.73s) | ✓ (54.04s) | |
| stepped_block_coarse.stl | 52 | ✓ (0.16s) | ✓ (0.16s) | |
| stepped_block_toroidal_fillets_fine.stl | 50962 | FAIL (60.21s) | FAIL (60.04s) | Stage 1.2 compare: vertex 1.2e-5mm from STEP |
| stepped_block_toroidal_fillets_medium.stl | 4576 | FAIL (5.57s) | FAIL (5.51s) | Stage 1.2 compare: vertex 1.1e-5mm from STEP |

#### Summary

- **Total test files**: 52
- **Passing (debug)**: 33 / 52
- **Passing (release)**: 33 / 52
- **Expected failures (bad/)**: 5
- **Failing non-bad tests**: 14 (3 filleted/toroidal, 2 pipe elbow, 2 stepped_block_toroidal_fillets, 1 block_with_conical_hole, 1 cone_tangent_planes, 1 nosecone, 1 chamfered_cube ASCII, 2 compare tolerance marginal)


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
- [ ] Hole detection and filling
- [ ] Feature recognition (holes, pockets, bosses)
- [ ] Curvature-based surface detection: use per-vertex principal curvatures (Rusinkiewicz estimator) as an alternative seeding signal. Principal curvature signatures: plane (k1≈k2≈0), sphere (k1≈k2, constant), cylinder (k1=const, k2≈0), cone (k1 varies, k2≈0), torus (k1=const, k2 varies).
- [ ] Quadric pre-classification: fit 10-coefficient quadric via SVD as a fast O(N) surface classifier before expensive nonlinear fitting. Eigenvalue pattern of the 3×3 quadric matrix directly classifies plane/sphere/cylinder/cone.
- [ ] Adaptive mesh density handling: current BFS tolerance assumes roughly uniform tessellation. For meshes with varying triangle density (e.g., adaptive refinement near features), weight BFS decisions by triangle area.

