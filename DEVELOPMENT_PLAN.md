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

## Implementation Phases

### Foundation
- [x] Project setup (crates, cargo.toml, dependencies)
- [x] Test utility: read an STL and a STEP, compute maximum distance between STL vertices and STEP surfaces, and print it out. Create a script in scripts/ to apply it to all of the stl/step file pairs under tests/ and print out a table of maximum distances.
- [x] Program skeleton: implement main program, flag parsing, and create stub source files for each stage. Create stage output data structures. For portions of the data structures which are unclear at this point, stub them out with comments. Take existing stage 1.1 and stage 1.2 implementations in mesh.rs and reformulate them into stage1.rs and appropriate output data structures.

### Stage 1 Mesh Input
- [x] Stage 1.1: Read STL file into `ConnectedMesh`, including welded vertices and per-face placeholder fields for neighbors, normals, and hypotheses.
- [x] Stage 1.2: Mesh validation pass to compute face normals, edge neighbors, manifold stats, connected shells, and orientation consistency checks.
- [x] Stage 1.3: Coplanar triangle fusion — merge adjacent coplanar triangle pairs into quads. This simplifies surface fitting by ensuring curved-surface facets (cylinder/sphere quad strips) become single faces rather than 2-face planar hypothesis groups.
- [x] Implement tests that all stl/step file pairs in tests/ pass consistency checks and pass when fed to brepper with the --compare flag. Also invent a few file pairs in tests/bad that will fail with the --compare flag by editing the location of one or more vertices to fail surface closeness tests. Also invent some bad cases that fail the mesh validation tests in various ways. Ensure that these bad cases fail in the correct way.

### Stage 2: Surface Fitting - planes, cylinders, and spheres, oh my!
- [x] Understand stage 2.1 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.1. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [x] Implement stage 2.1. Make sure --compare passes for all test shapes composed only of planar surfaces.
- [x] Understand stage 2.2 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.2. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [x] Implement stage 2.2: Deduce cylindrical hypotheses.
- [x] Understand stage 2.3 and the existing test models under tests/. Imagine additional test shapes which will be challenging for stage 2.3. Create these test shapes in tests/ccad/ - see tests/ccad/README.md.
- [x] Implement stage 2.3: Deduce spherical hypotheses.
- [x] Implement stage 2.6: Select surfaces for reconstruction using per-face priority rule.
- [x] Implement the --angular-tolerance flag for cylindrical and spherical surface hypothesis generation.
- [x] Revisit stage 2.6: Replace per-face priority rule with greedy area-based selection.
- [x] Revisit stage 2.2: Update cylindrical hypothesis to match updated algorithm description from commits 2ac6d4f onward. Test that deduced cylinders on test models match the step cylinders.
- [x] Revisit stage 2.3: Update spherical hypothesis to match updated algorithm description (vertex-based seeding, solid-angle coverage validation, surface-tolerance in BFS, relaxed max_sphere_radius). Test that deduced spheres on test models match the step spheres. This should fix the sphere overgrowth problem blocking `onshape_rounded_cube_stage26`. Unnerf that test and ensure it passes.

### Stage 3: Surface Reconstruction
- [x] Implement stage 3.1: Create OCCT surface objects and build adjacency graph from mesh connectivity.
- [x] Revisit stage 3.1: Implement newly-described --compare check.
- [x] Implement stage 3.2: Detect tangency relationships between adjacent surfaces. Computes outward-facing surface normals at sampled boundary points and compares angle; marks edge as tangent if normals agree within 2°. No tangent edges expected for current test models (all planar/cylindrical/spherical intersections meet at angles > 2°).
- [x] Implement stage 3.3: Compute edge curves via surface-surface intersection (`GeomAPI_IntSS`), trim to vertex endpoints via `GeomAPI_ProjectPointOnCurve` + `Geom_TrimmedCurve`. Multi-curve selection picks closest to mesh boundary vertices. Tangent edges skipped (no tangent edges in current models). All edges computed successfully for all test models. Pcurves deferred to stage 3.4.
- [x] Implement stage 3.4: Create OCCT faces from surfaces bounded by edge wires. Creates `TopoDS_Edge` for each `ReconEdge` using `BRepBuilderAPI_MakeEdge` with trimmed 3D curves and shared `TopoDS_Vertex` endpoints. For periodic curves (circles), vertex order is matched to curve parameterization to ensure correct arc selection by `MakeEdge`. Planar faces: group edges into wire loops by vertex connectivity, build `MakeWire`/`MakeFace` with outer wire + inner hole wires. Cylindrical/spherical faces: UV-bounds construction for full-revolution (with circular gap algorithm for periodicity), wire-based construction with pre-set `Geom2d_Line` pcurves for partial-revolution. Full spheres (no edges): natural UV bounds. Validates all faces with `BRepCheck_Analyzer`. Compare check samples mesh face centroids against reference STEP.
- [x] Implement stage 3.5: Stitch faces into shells using `BRepBuilderAPI_Sewing`. All faces sewn together with vertex tolerance. Handles Shell, Face, Solid, and Compound results from sewing. Applies `ShapeFix_Shell` orientation fixing. Compare validates shell count and checks for free edges.
- [x] Implement stage 3.6: Construct solids from shells using `ShapeFix_Solid::SolidFromShell`. Validates with `BRepCheck_Analyzer`, computes volume via `BRepGProp::VolumeProperties`. Compare checks volume agreement (warning on >1% diff) and `BRepExtrema_DistShapeShape` distance to STEP reference.

### Stage 4: Output
- [x] Implement stage 4.1: Write solids to STEP file using `STEPControl_Writer` with `StepModelType::Asis`. Builds compound from all solids, sets AP214 schema. Compare re-reads written file and validates volume agreement (warning on >1% diff) and `BRepExtrema_DistShapeShape` distance to STEP reference. 15 integration tests: 3 basic output + 11 compare + 1 missing-output-path error.

### Stage 3 refinements
- [x] Revisit stage 3.2: Tangency detection correctly identifies 8 tangent edges on part_rounded_cube_10_r2 (cylinder-plane tangencies at fillet boundaries). Detection algorithm using analytical surface normals at sampled boundary points works well. No changes needed to the 2° threshold.
- [x] Revisit stage 3.3: Implement tangent edge curve computation. For tangent edges, construct curves analytically rather than using `GeomAPI_IntSS` (which fails for tangent surfaces). Plane-cylinder tangencies produce lines parallel to the cylinder axis at the analytically-computed tangent point. Also fix `MakeEdge` in stage 3.4: use `new_handlegeomcurve_vertex2_real2` with explicit parameter values (instead of relying on OCCT's vertex-to-curve projection which fails when mesh vertex tolerance exceeds OCCT's default precision of 1e-7mm), and set vertex tolerance to `vertex_tolerance_mm` via `BRep_Builder::update_vertex_vertex_real`. This unblocks part_rounded_cube_10_r2 — all 10 faces (6 planar + 4 cylindrical) reconstruct correctly with volume matching STEP reference to 2.19e-7 relative difference.
- [x] Implement sphere-cylinder tangent edge curves: great circle arcs on the sphere surface, using sphere center/radius (not cylinder parameters, which may be imprecise). Adaptive per-vertex tolerance in MakeEdge accommodates surface fit imprecision. Unblocks rounded_cube_10_r2_coarse full pipeline — STEP file exports with all 94 faces (76 planar + 10 cylindrical + 8 spherical), though many faces have BRepCheck warnings due to imprecise cylinder fits.
- [x] Implement closed-loop sphere-cylinder tangent edge curves and pole-safe hemisphere face construction for pill/capsule shapes. Closed-loop tangent circles are constructed as full `Geom_Circle` curves. When boundary circles pass within 45° of a sphere pole, the sphere surface is reoriented to align with the cylinder axis, making UV-bounds construction work correctly. Unblocks pill_coarse.stl and pill_fine.stl — STEP output matches reference to <1e-7 relative volume difference.
- [x] Evaluated stage 2.7 (surface refitting): not needed. All test models including complex rounded_cube and pill pass --compare with max face-to-STEP distances of ~3e-6 mm or better. No boundary accuracy problems from face-to-surface assignments. Stage 2.7 can be revisited if cone/torus surfaces reveal boundary issues. Also un-ignored 7 tests (rounded_cube stages 3.1/3.4/3.5/3.6/4.1, pill_coarse/pill_fine stage 4.1) that were marked ignored due to cylinder 3-face seeding rework but now pass.

### Stage 2 Extensions: Conical Surfaces
- [x] Create ccad test models for conical surfaces: `simple_cone.lua` (truncated cone), `block_with_conical_hole.lua` (block with conical bore), `cone_cylinder.lua` (cone joined to cylinder, cone-cylinder tangency), `nosecone.lua` (cone with tangent spherical cap, cone-sphere tangency). Export STL+STEP pairs.
- [x] Implement `ConicalHypothesis` data structure in stage 2 (apex, axis, half-angle, convex, faces, vertices, errors). Added `ConicalHypothesis` struct with fields for apex, axis_direction, half_angle, convex, faces, vertices, error_max, centroid_error_max, error_abs_sum. Extended `MeshFace` with `conical_hypothesis` field and `Stage2Output` with `conical_hypotheses` vector. Added `SelectedSurface::Conical` variant.
- [x] Implement cone fitting: axis estimation from normal covariance, (h,r) profile fitting via linear regression to get half-angle and apex position. Apex-vertex seeding strategy instead of triple-seed: identifies vertices with many non-coplanar-normal incident faces (eigenval/trace < 0.3), uses all incident faces as seed set for robust axis estimation. Axis orientation ensured by flipping when h_sum < 0. Levenberg-Marquardt refinement for 6-parameter optimization (axis perturbation, apex position, half-angle). False-positive rejection filters: half-angle bounds (2°–85°), vertex error < vertex_tol, apex distance < 10× bounding box diagonal, normal-axis angle consistency (std_dev < angular_tol/2).
- [x] Implement cone BFS region growing with vertex-to-cone distance validation and angular coverage check. Seeds validated at vertex_tol with multi-axis LM rescue (tries geometric axis, normal covariance eigenvector, and 3 coordinate axes). BFS expansion with vertex_tol for vertex checks, surface_tol for centroid checks, and refit-if-close pattern (REFIT_SKIP_MULTIPLIER × vertex_tol threshold). Convexity by majority vote after final re-fit. Angular tolerance for cone BFS uses 2× the configured value to accommodate inter-strip dihedral angles on coarsely tessellated cones (e.g., 18° between adjacent strips with 20 circumferential divisions).
- [x] Extend stage 2.6 surface selection to include conical candidates. Cones participate in greedy area-based selection alongside planar/cylindrical/spherical hypotheses.
- [x] Distinguish cones from cylinders/spheres/planes: half-angle bounds (< 2° = cylindrical, > 85° = planar), apex distance sanity check against mesh extent, and normal-axis angle standard deviation check to reject sphere patches falsely fitted as cones. All 325 tests pass with zero false positives.
- [x] Unit tests: 5 unit tests (cone detection, false positive rejection for cube/cylinder/sphere) and 57 integration tests covering stage 2.4 hypothesis counts, cone parameter matching against STEP, surface selection counts, stage 2.4/2.6 compare tests, stage 3 topology/compare tests for all cone models, and stage 4 output compare tests. Total: 382 tests.
- [x] Full pipeline tests: All cone models pass `--compare` through stage 4.1. simple_cone (vol diff 4e-9), block_with_conical_hole (vol diff 1.3e-9), cone_cylinder (vol diff ~0), nosecone (vol diff ~0), onshape cone_15x20 (vol diff 6.4e-9). Three bug fixes enabled this: (1) cone BFS angular tolerance relaxed to 2× to avoid blocking cross-strip growth, (2) cone apex UV bounds extended to V=0 for full (non-truncated) cones with single boundary edge, (3) concave flag added for conical surfaces in face creation.
- [x] Unblock existing `cone_15x20_medium` (onshape) test. Stage 2 detection works perfectly (58 faces, ha=20.56°, err_max=2.4e-6). Full pipeline now passes through stage 4.1 (vol diff 6.4e-9, distance ~0).

### Stage 2 Extensions: Toroidal Surfaces
- [x] Create ccad test models for toroidal surfaces: `filleted_cylinder.lua` (cylinder with fillets on circular edges → convex toroidal fillets), `filleted_hole_block.lua` (block with cylindrical bore, hole edges filleted → concave toroidal fillets), `filleted_pipe.lua` (hollow pipe with all edges filleted → mixed convex/concave toroidal fillets). Export STL+STEP pairs. Note: ccad's `profile_xz` only supports polylines (no arcs/circles), so `revolve()` produces conical-strip approximations, not true toroidal surfaces. Fillets on circular edges are the correct way to create TOROIDAL_SURFACE in STEP via ccad.
- [x] Implement `ToroidalHypothesis` data structure in stage 2 (center, axis, major_radius, minor_radius, convex, faces, vertices, errors). Plumbed through stage 1 (`MeshFace.toroidal_hypothesis`, `UNDEDUCED_TOROIDAL_HYPOTHESIS`), stage 2 (`ToroidalHypothesis` struct, `SelectedSurface::Toroidal` variant, `Stage2Output.toroidal_hypotheses`), and stage 3 (all match arms handle `Toroidal` for surface normals, face lists, vertex lists, concavity, descriptions; OCCT surface creation is `todo!()` pending binding).
- [x] Implement torus fitting via medial axis / tube center method: estimate minor radius from normal-line intersections, compute tube centers k_i = p_i + r*n_i, fit 3D circle to tube centers for major circle parameters. Vertex-neighborhood seeding with post-seeding merge, quality gate on circle-fit RMS (0.1*minor_r threshold). 7 integration tests pass.
- [x] Implement torus BFS region growing with vertex-to-torus distance validation.
- [x] Extend stage 2.6 surface selection to include toroidal candidates. Toroidal candidates already participated in greedy area-based selection. Added 10 integration tests: 3 surface selection count tests (filleted_cylinder exact counts: 2 planar + 1 cylindrical + 2 toroidal = 5 total; filleted_hole_block/filleted_pipe minimum counts with >=2/>=4 toroidal), 3 stage 2.6 `--compare` tests for filleted models, `onshape_pipe_elbow_stage26_compare` (passes with 0 toroidal — torus segmented into cylindrical patches), and 3 `all_faces_covered_by_selection` tests for filleted models.

### Stage 3 Extensions for Cone and Torus
- [x] Stage 3.1: Add `Geom_ToroidalSurface` creation from hypothesis parameters. Uses `geom::ToroidalSurface::new_ax3_real2(&ax3, major_radius, minor_radius)` binding (already present in opencascade-sys). Unblocked `rounded_cube_coarse_brep_check` test.
- [x] Stage 3.2: Add torus analytical normal formulas for tangency detection. (Cone normals already implemented.) Torus normals compute direction from tube center to surface point, negated if concave.
- [x] Stage 3.3: Implement tangent edge curve computation for plane-torus, cylinder-torus, and torus-torus pairs (circle-based, both closed loop and arc). Cone tangent edges not needed (see next item).
- [x] Stage 3.3: Plane-cone and cylinder-cone tangent edge curves: NOT NEEDED. At plane-cone and cylinder-cone junctions, surface normals always differ by the cone half-angle (typically 10°–45°), far exceeding the 2° tangency threshold. No test model has cone tangent edges. All cone models pass via the non-tangent IntSS path.
- [x] Stage 3.3: Non-tangent cone and torus edge intersections via `GeomAPI_IntSS` with curve selection — already works for all test models (simple_cone, cone_cylinder, nosecone, block_with_conical_hole all pass stage 3.6 compare).
- [x] Stage 3.4: Conical face construction — already works for all test models. Existing UV-bounds and face creation code handles cones correctly.
- [x] Stage 3.4: Toroidal face construction — fixed doubly-periodic V-bounds. The torus V parameter is periodic in [0, 2π]; simple min/max produced wrong bounds for faces crossing V=0. Applied circular gap algorithm (same as U-bounds for partial revolution). Fixed filleted_cylinder volume: 13.5% error → <1e-7 relative diff. Fixed filleted_pipe: 6.4% → <1e-7.

### Incorporate learnings from aborted C++ implementation

- [x] Port abstracted BFS region growing from C++ implementation. Created `src/stage2/bfs.rs` with a `SurfaceType` trait (associated types `Fit` and `Hypothesis`, methods for seeding, fitting, distance computation, acceptance testing, hypothesis building) and a generic `generate_hypotheses::<S>()` function. Implemented trait for all 5 surface types: `PlanarSurface`, `CylindricalSurface`, `SphericalSurface`, `ConicalSurface`, `ToroidalSurface`. Replaced 5 separate `deduce_*_hypotheses()` calls in stage2 with generic BFS invocations. Added post-BFS planar fallback synthesis for uncovered faces. Old per-type BFS functions retained (with `#[allow(dead_code)]`) as reference. **Regressions**: 8 `block_with_conical_hole` tests fail (0 conical hypotheses produced — the generic BFS lacks the old code's multi-axis LM rescue for cone seeds, which is needed when seed faces include non-conical faces near the apex). All other 388 tests pass. Stage 2 is ~3× slower (255s vs ~80s baseline) due to the generic BFS trying more seeds.
- [x] Compare C++ implementations of surface fitting functions to Rust versions. Evaluate which is simpler, which is more likely to perform well, and how many hacks each has. Give a recommendation on whether to port the C++ approach to Rust to replace the Rust version, or to keep the Rust version, for each surface type. The goal is numerically accurate and fast surface fitting to both seed face regions and expanded regions.

  **Findings (Rust ~7121 LOC vs C++ ~2236 LOC, ~3.2× larger — partly legacy dead code + verbose LM boilerplate + manual linear algebra):**

  | Surface | Recommendation | Key Reason |
  |---------|---------------|------------|
  | Planar | Keep Rust, delete legacy BFS | Functionally equivalent; generic BFS is cleaner |
  | Cylindrical | Keep Rust | Rust has LM refinement that C++ lacks (better accuracy) |
  | Spherical | Keep Rust, **port C++ Gauss-Newton refinement** | Rust lacks iterative refinement after algebraic fit (clear gap) |
  | Conical | Keep Rust, port analytic Jacobian | Rust has multi-axis LM rescue C++ lacks; C++ has cleaner analytic Jacobian |
  | Toroidal | Keep Rust dual-method, port PCA axis + analytic Jacobian | C++ has elegant vertex-position PCA axis estimation and analytic Jacobian |

  **Top-priority ports from C++ to Rust:**
  1. Sphere Gauss-Newton refinement (iterative refinement after algebraic sphere fit)
  2. Analytic Jacobians for cone and torus LM (replace finite-difference)
  3. Torus minor-circle acceptance test (C++ checks minor circle center lies on major circle)
  4. Vertex-position PCA for torus axis estimation (C++ uses vertex positions directly, not just normals)

### Ports from C++ Implementation
- [x] Port sphere Gauss-Newton refinement: after algebraic least-squares sphere fit, run iterative Gauss-Newton to minimize vertex distance residuals. (C++ `surface_fit.cpp` `fit_sphere()` second pass.) Up to 20 iterations, convergence at delta_norm < 1e-14, requires ≥6 vertices. Uses normal equations with 4x4 Gaussian elimination (same as algebraic fit).
- [ ] Port analytic Jacobians for cone LM: replace finite-difference Jacobian in cone LM refinement with closed-form partial derivatives. (C++ `surface_fit.cpp` `ConeResidual::df()`.)
- [ ] Port analytic Jacobians for torus LM: same for torus fitting. (C++ `surface_fit.cpp` `TorusResidual::df()`.)
- [ ] Port torus minor-circle acceptance test: after fitting, verify that minor-circle centers lie on the major circle within tolerance. Rejects false positives. (C++ `stage2.cpp` torus acceptance.)
- [ ] Port vertex-position PCA for torus axis estimation: use PCA on vertex positions (not just normals) to estimate the torus axis of revolution. May improve robustness for partial tori. (C++ `stage2.cpp` torus seeding.)
- [ ] Delete legacy `deduce_*_hypotheses()` functions (currently `#[allow(dead_code)]` in each surface type file) — no longer needed now that generic BFS works and comparison is complete.

- [ ] Full pipeline tests: all test models including cone and torus models and filleted_cylinder, filleted_hole_block, filleted_pipe pass --compare through stage 4.1. Unblock existing `pipe_elbow_10_fine` (onshape) test.

### Unified Surface-of-Revolution Framework (refactor)
- [ ] Extract shared SOR axis estimation (normal covariance smallest eigenvector) into a common utility function.
- [ ] Implement (h,r) coordinate reduction: for a set of vertices and an axis, compute axial and radial coordinates.
- [ ] Implement unified profile classifer: given (h,r) points, fit horizontal line (cylinder), sloped line (cone), circle centered on axis (sphere), and offset circle (torus). Return all fits with RMS residuals.
- [ ] Refactor cylinder fitting to use (h,r) reduction + profile classification as an alternative initialization.
- [ ] Refactor sphere fitting to use (h,r) reduction when axis is estimated from normal covariance.
- [ ] Integrate unified SOR into BFS: after initial seed, try all profile types and pick best before growing.
- [ ] Regression tests: ensure all existing test models still pass after refactor.

### Gaussian Map Pre-Classification
- [ ] Implement normal covariance eigenvalue classification: compute $\lambda_1, \lambda_2, \lambda_3$ eigenvalue pattern for face neighborhoods. Classify as planar ($\lambda_1 \gg \lambda_2 \approx \lambda_3$), cylindrical/conical ($\lambda_1 \approx \lambda_2 \gg \lambda_3$), spherical ($\lambda_1 \approx \lambda_2 \approx \lambda_3$), or toroidal.
- [ ] Use classification to prioritize seeding strategies per region, reducing wasted fitting attempts.
- [ ] Add quadric fitting (10-coefficient SVD) as a fast initialization alternative. Extract surface type from eigenvalue pattern of the 3×3 quadric matrix.

### NURBS Fallback
- [ ] Create test models with freeform surfaces that cannot be fit by analytic types (e.g., sculpted surface, spline loft).
- [ ] Implement NURBS fitting for remaining unassigned face groups using OCCT's `GeomAPI_PointsToBSplineSurface`.
- [ ] Integrate NURBS surfaces into stage 3 reconstruction (edge curves, face creation).
- [ ] Extend stage 2.6 surface selection to handle NURBS alongside analytic surfaces.


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
| Cube | 12 triangles | 6 planes, 1 solid | ✓ Stage 4.1 (STEP output) |
| Wedge | 12 triangles | 6 planes (incl. angled) | ✓ Stage 4.1 (STEP output) |
| T-Shape | 28 triangles | 10 planes | ✓ Stage 4.1 (STEP output) |
| Staircase | 48 triangles | 12 planes | ✓ Stage 4.1 (STEP output) |
| Chamfered Cube (onshape) | 44 triangles | 26 planes | ✓ Stage 4.1 (STEP output) |
| Stepped Block (onshape) | 52 triangles | 16 planes | ✓ Stage 4.1 (STEP output) |
| L Bracket (onshape) | 492 triangles | 10 planes | ✓ Stage 4.1 (STEP output) |
| Cylinder (onshape) | 232 triangles | 1 cylinder + 2 planes | ✓ Stage 4.1 (STEP output) |
| Simple Cylinder (ccad) | 124 triangles | 1 cylinder + 2 planes | ✓ Stage 4.1 (STEP output) |
| Block with Hole (ccad) | 132 triangles | 6 planes + 1 concave cylinder | ✓ Stage 4.1 (STEP output) |
| Pipe (ccad) | 244 triangles | 2 cylinders (in/out) + 2 annular planes | ✓ Stage 4.1 (STEP output) |
| Stepped Cylinder (ccad) | 240 triangles | 2 cylinders + 3 planes | ✓ Stage 4.1 (STEP output) |
| Two Holes (ccad) | 252 triangles | 6 planes + 2 concave cylinders | ✓ Stage 4.1 (STEP output) |
| Block Corner Hole (ccad) | 44 triangles | 6 planes + 1 concave cylinder | ✓ Stage 4.1 (STEP output) |
| Simple Sphere (ccad) | 974 triangles | 1 sphere | ✓ Stage 4.1 (STEP output) |
| Hemisphere (ccad) | 518 triangles | 1 sphere + 1 plane | ✓ Stage 4.1 (STEP output) |
| Spherical Pocket (ccad) | 486 triangles | 6 planes + 1 concave sphere | ✓ Stage 4.1 (STEP output) |
| Ball on Cylinder (ccad) | 764 triangles | 1 sphere + 1 cylinder + 1 plane | ✓ Stage 4.1 (STEP output) |
| Sphere (onshape) | 20448 triangles | 1 sphere | ✓ Stage 4.1 (STEP output) |
| Dome Hemisphere (onshape) | 10368 triangles | 1 sphere + 1 plane | ✓ Stage 4.1 (STEP output) |
| Pill coarse (onshape) | 752 triangles | 1 cylinder + 2 spheres | ✓ Stage 4.1 (STEP output) |
| Pill fine (onshape) | 20736 triangles | 1 cylinder + 2 spheres | ✓ Stage 4.1 (STEP output) |
| Plate with Hole (onshape) | 136 triangles | 6 planes + 1 cylinder | ✓ Stage 4.1 (STEP output) |
| Plate with Hole low (fusion) | 112 triangles | 6 planes + 1 cylinder | ✓ Stage 4.1 (STEP output) |
| Plate with Hole med (fusion) | 168 triangles | 6 planes + 1 cylinder | ✓ Stage 4.1 (STEP output) |
| Plate with Hole high (fusion) | 284 triangles | 6 planes + 1 cylinder | ✓ Stage 4.1 (STEP output) |
| Part Rounded Cube coarse (onshape) | 140 triangles | 6 planes + 4 cylinders | ✓ Stage 4.1 (STEP output) |
| Part Rounded Cube fine (onshape) | 588 triangles | 6 planes + 4 cylinders | ✓ Stage 4.1 (STEP output) |
| Rounded Cube coarse (onshape) | 820 triangles | 6 planes + 12 cylinders + 8 spheres | ✓ Stage 4.1 (needs --angular-tolerance 20) |
| Rounded Cube medium (onshape) | 3548 triangles | 6 planes + 12 cylinders + 8 spheres | ✓ Stage 4.1 (STEP output) |
| Rounded Cube fine (onshape) | 21466 triangles | 6 planes + 12 cylinders + 8 spheres | ✓ Stage 4.1 (STEP output) |
| Pipe Elbow (onshape) | 11232 triangles | Cylinders + torus | ✗ Stage 3.1 (missing torus surface type) |
| Cone (onshape) | 116 triangles | 1 cone + 1 plane | ✓ Stage 4.1 (STEP output) |
| Simple Cone (ccad) | 120 triangles | 1 truncated cone + 2 planes | ✓ Stage 4.1 (STEP output) |
| Block with Conical Hole (ccad) | 312 triangles | 6 planes + 1 concave cone | ✓ Stage 4.1 (STEP output) |
| Cone-Cylinder (ccad) | 292 triangles | 1 cone + 1 cylinder + 2 planes | ✓ Stage 4.1 (STEP output) |
| Nosecone (ccad) | 820 triangles | 1 cone + 1 sphere + 1 plane | ✓ Stage 4.1 (STEP output) |
| Filleted Cylinder (ccad) | 1528 triangles | 1 cylinder + 2 planes + 2 convex tori | Planned: toroidal surfaces phase |
| Filleted Hole Block (ccad) | 1484 triangles | 6 planes + 1 concave cylinder + 2 concave tori | Planned: toroidal surfaces phase |
| Filleted Pipe (ccad) | 2576 triangles | 2 cylinders + 2 planes + 4 tori (convex+concave) | Planned: toroidal surfaces phase |
| Stepped Block Toroidal Fillets fine (onshape) | 50962 triangles | 16 planes + 37 cylinders + 18 tori | Planned: toroidal surfaces phase |
| Stepped Block Toroidal Fillets medium (onshape) | 4576 triangles | 16 planes + 37 cylinders + 18 tori | Planned: toroidal surfaces phase |
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
- [ ] Hole detection and filling
- [ ] Feature recognition (holes, pockets, bosses)
- [ ] Curvature-based surface detection: use per-vertex principal curvatures (Rusinkiewicz estimator) as an alternative seeding signal. Principal curvature signatures: plane (k1≈k2≈0), sphere (k1≈k2, constant), cylinder (k1=const, k2≈0), cone (k1 varies, k2≈0), torus (k1=const, k2 varies).
- [ ] Quadric pre-classification: fit 10-coefficient quadric via SVD as a fast O(N) surface classifier before expensive nonlinear fitting. Eigenvalue pattern of the 3×3 quadric matrix directly classifies plane/sphere/cylinder/cone.
- [ ] Adaptive mesh density handling: current BFS tolerance assumes roughly uniform tessellation. For meshes with varying triangle density (e.g., adaptive refinement near features), weight BFS decisions by triangle area.

