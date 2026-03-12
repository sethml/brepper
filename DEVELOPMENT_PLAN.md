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
Fit cylindrical hypotheses to connected sets of faces that lie on a common cylinder. A cylindrical hypothesis consists of:
- `axis_origin: [f64; 3]` — A point on the cylinder axis.
- `axis_direction: [f64; 3]` — Unit direction vector along the axis.
- `radius: f64` — Radius of the cylinder (always positive).
- `convex: bool` — Whether face normals point away from the axis (convex, like the outside of a pipe) or toward it (concave, like the inside of a hole).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this cylinder.
- `error_max: f64` — Maximum absolute distance from any vertex to the cylinder surface (i.e. max of `| ||v - axis_closest|| - radius |`).
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

The algorithm uses BFS region growing, analogous to stage 2.1 but seeded from triples of nearby faces (since a single triangle is always planar and cannot determine curvature, and two faces provide minimal angular constraint on the cylinder axis). A key challenge is that **any two non-coplanar faces define some cylinder** — the cross product of their normals gives an axis, and a circle fit always produces a radius. The previous algorithm took the first valid seed and committed to it, which meant a bogus 3-face cylinder from an unfortunate seed pair could permanently consume faces that should belong to a correct 20-face cylinder seeded from a different pair. The fix is multi-seed evaluation with 3-face seeds that span more of the cylinder's circumference:

**Pairwise qualification:**
Two faces fi and ni are pairwise qualified to participate in a cylindrical hypothesis if:
- fi and ni are neighbors.
- fi and ni are unassigned (`cylindrical_hypothesis == -2`).
- fi and ni have a sufficiently different normal (cross product magnitude exceeds `MIN_CROSS_THRESHOLD`, e.g. 0.01).
- The dihedral angle between fi and ni does not exceed `--angular-tolerance` (default 17.5°).

**Cylinder fit qualification:**
A face or set of faces is cylinder fit qualified if:
- All vertices are within `--vertex-tolerance` of the cylinder.
- All face centroids are within `--surface-tolerance` of the cylinder.
- **Convexity check**: compute the vector from the axis to the face's centroid. Check that the face's normal agrees with the hypothesis convexity (dot product with radial vector has the expected sign). Reject if inconsistent — this prevents merging the inner and outer surfaces of a thin-walled cylinder.

**Multi-seed evaluation:**
- Initialize all `cylindrical_hypothesis` to `UNDEDUCED` (-2).
- For each face fi where `cylindrical_hypothesis == -2`:
    - For each neighbor n1 of fi where n1 is pairwise qualified with fi:
        - For each neighbor n2 of n1 where n2 is not a neighbor of fi and n2 is pairwise qualified with n1:
            - Produce a cylinder fit to the set of vertices belonging to fi, n1, and n2.
            - If faces fi, n1, and n2 are cylinder fit qualified, accept this set of faces as a valid seed and perform BFS expansion on them using temporary data structures (a local face list and vertex set — do not mutate `mesh.faces[].cylindrical_hypothesis` during trials). The initial BFS queue should contain fi, n1, and n2, so we explore all of their neighbors. Each trial produces a candidate hypothesis with a face count and fitted parameters. After each trial completes with acceptance, compare against the **current best candidate** for this starting face; keep whichever has more area covered. This avoids storing all trial results.
    - If no valid (fi, n1, n2) triple was found, skip fi (leave as `UNDEDUCED` — it may be absorbed by a later face's BFS, or will remain unassigned at the end).
- Accept the **best candidate** as a cylindrical hypothesis: create a cylindrical hypothesis, and for every face in the hypothesis, assign the face's cylindrical hypothesis field.
- Apply **angular coverage validation** (see below) to each trial. If it fails, discard the trial.
- Accept the **best candidate** as a cylindrical hypothesis: create a cylindrical hypothesis, and for every face in the hypothesis, assign the face's cylindrical hypothesis field.
- After all faces have been processed, set any remaining `UNDEDUCED` faces to `NO_HYPOTHESIS`.

Future work (not currently implemented - in the future we'll evaluate whether these optimizations help):
    - **Early termination** of redundant trials: during a trial BFS, if a best candidate already exists, check whether the trial is rediscovering the same cylinder. After the trial adds a face, compare its fitted cylinder parameters against the best candidate's: axis directions nearly parallel (`|a_trial · a_best| > 1 - 1e-6`), radii within `--vertex-tolerance`, and axis-to-axis distance within `--vertex-tolerance`. If parameters match AND every face accepted so far is already in the best candidate's face set, abandon the trial — it can only produce a subset of the existing best. This is the dominant optimization: on a 20-face fillet cylinder with ~3 valid seeds per face, most of the ~60 trials terminate after 2–3 faces instead of exploring all 20.

**BFS expansion** (run once per trial seed):
- Repeatedly pop a face fi from the queue and examine each neighbor ni:
    - If ni is already assigned a cylindrical hypothesis (in the mesh) or claimed by this trial, skip.
    - **Angular tolerance check**: for each of ni's mesh neighbors that is already assigned in this trial, compute the dihedral angle between ni and that neighbor. If any exceeds `--angular-tolerance`, skip ni. Checking all assigned neighbors (not just the BFS parent) provides defense-in-depth against creased NURBS surfaces where BFS might approach a face from a low-angle direction while a high-angle assigned neighbor exists on a different axis. (Note: the full pairwise qualification — including the cross product minimum and unassigned requirement — is NOT applied during BFS expansion. The cross product minimum is a seeding criterion only: on a finely tessellated cylinder, adjacent faces legitimately have near-parallel normals with tiny cross products. The unassigned criterion is N/A since trial-claimed faces are logically assigned within the trial.)
    - If ni fails cylinder fit qualification:
        - If any vertex of ni is further than `2 * vertex_tolerance` from the cylinder or ni's centroid is further than `2 * surface_tolerance` from the cylinder, skip ni. Otherwise, do a re-fit attempt:
        - **Re-fit attempt**: re-fit the cylinder from all current faces plus ni (see Cylinder Fitting below). Check that **all** faces (existing and new) are cylinder fit qualified with the re-fit. If so, accept the re-fit; otherwise reject the re-fit cylinder parameters and skip ni.
    - If accepted: add ni to the trial's face list and vertex set, push ni onto the BFS queue.
- After BFS completes: do a final re-fit from all accumulated faces and vertices. Compute error metrics.

**Angular coverage validation** (applied to each trial):

The problem: any two non-coplanar faces fit *some* cylinder. A bogus seed that doesn't span the cylinder's circumference can grow into a small hypothesis via BFS, consuming faces that should belong to a correct hypothesis from a better seed. The angular coverage check ensures the hypothesis has genuine circumferential support — faces distributed around the cylinder, not clustered on one side.

Algorithm:
1. For each face in the hypothesis, compute the centroid and project it onto the cylinder's angular coordinate θ. Given the cylinder's axis direction `a`, axis origin `o`, and orthogonal basis vectors `u`, `w` perpendicular to `a`: for centroid `c`, compute `radial = c - o - ((c - o) · a) * a`, then `θ = atan2(radial · w, radial · u)`.
2. Sort the θ values (one per face).
3. Compute all N gaps between consecutive sorted θ values, **including the wraparound gap**: `gap_wrap = θ[0] + 2π - θ[N-1]`.
4. The **largest gap** is the "empty arc" — the angular region not covered by any face. The **span** is `2π - largest_gap`.
5. Among the remaining N-1 gaps, the **second-largest gap** must be `≤ span / 3`.

This ensures the hypothesis has at least 3 distinct angular clusters of faces around the cylinder. A bogus cylinder where all faces are side-by-side in a narrow arc will fail: e.g., 5 faces clustered in a 30° arc have span ≈ 30°, but the second-largest internal gap (~7°) may pass. However, such a cluster's centroid validation and vertex tolerance checks should already reject it. The angular coverage check primarily catches the case where a seed grabs faces from 2 sides of a non-cylindrical surface (e.g., two planar faces that happen to fit a cylinder) and BFS adds a few more nearby faces, resulting in a hypothesis with faces clustered in 1–2 angular positions rather than distributed around the circumference.

*Note on applying angular coverage during BFS growth: this was considered but is not useful. Early in BFS, faces are naturally clustered near the seed, so the check would fail prematurely. The vertex distance, convexity, and angular tolerance checks are the correct per-face acceptance criteria during growth. Angular coverage is a global structural validation — it answers "does this completed hypothesis look like a real cylinder?" rather than "should this next face be added?"*

*Note: For partial cylinders (arcs less than 360°), this check still works — it requires 3+ angular clusters within whatever arc is covered, which is the minimum for a reliable cylinder fit. A 90° cylindrical arc with 6 evenly-spaced faces has internal gaps of ~18° and span of ~90°, so second_largest_gap (18°) ≤ span/3 (30°) — passes correctly.*

**Computational cost**: For each starting face with K valid seed partners (neighbors passing pairwise qualification), each partner n1 may have L additional neighbors qualifying as n2, producing up to K×L trial BFS runs per face. Typical K is 1–3 (most faces have 3–4 neighbors, of which 1–2 have sufficiently different normals) and L is 1–2 (n1's neighbors minus fi and minus n1's non-qualifying neighbors). So typically 1–6 trials per face, each of which is a BFS that covers the connected cylinder region. Without early termination, total work is proportional to `sum over faces of (K_i × L_i × region_size_i)`. Since most faces on the same cylinder will discover the same region, the work per cylinder is roughly `N_faces × K × L × N_faces`, dominated by the first few faces that discover the region before their neighbors get committed. The temporary data structures (HashSet for claimed faces, Vec for face list) are lightweight.

**Cylinder fitting** (used for seeding and re-fitting):

The fitting proceeds in three steps:

1. **Axis direction estimation from face normals.** On a cylinder, every face normal is perpendicular to the axis: `n · a = 0`. So the axis direction `a` minimizes $\sum w_i (n_i \cdot a)^2$ subject to $\|a\| = 1$, where $w_i$ is the area of face $i$. This is the eigenvector corresponding to the *smallest* eigenvalue of the 3×3 weighted covariance matrix $M = \sum w_i \, n_i \, n_i^T$. For the three-face seed, this is the full 3×3 eigenvector computation (the cross-product shortcut `a = normalize(n1 × n2)` that works for two-face seeds is not used, since we want the axis estimate that best fits all three faces' normals).

2. **Axis position and radius via 2D circle fitting.** Given the axis direction `a`, choose two orthogonal unit vectors `u`, `w` perpendicular to `a`. Project each vertex onto the 2D plane: `(x_i, y_i) = (v_i · u, v_i · w)`. **Center the 2D coordinates** by subtracting the mean: `x_i' = x_i - mean_x`, `y_i' = y_i - mean_y`. Fit a circle to the centered 2D points using the algebraic least-squares method: fit the linear model $x^2 + y^2 + Dx + Ey + F = 0$ via least squares, then `center' = (-D/2, -E/2)` and `radius = sqrt(center_x'² + center_y'² - F)`. The axis origin in 3D is `(center_x' + mean_x) * u + (center_y' + mean_y) * w` (the component along `a` is arbitrary). Centering is essential for numerical stability: without it, the normal equations matrix becomes ill-conditioned when absolute vertex coordinates are large relative to the arc span (e.g., vertices at (15mm, 12.5mm) spanning 2.5° of a r=2mm cylinder produce sums-of-squares ~1350 that dwarf the r=2 signal, causing wildly wrong radii).
3. **Levenberg-Marquardt refinement.** The initial axis + circle fit provides a good starting point but can have ~0.4° axis error when faces span a narrow arc (<15°). LM refines all 6 parameters simultaneously: `[alpha, beta, qx, qy, qz, radius]` where alpha/beta are tilt angles from the initial axis in two perpendicular directions. The residual for each vertex is `||cross(p - q, a)|| - r` (signed distance to cylinder surface). Numerical Jacobian via central differences (eps=1e-8). Uses `levenberg-marquardt` crate v0.14 with `nalgebra` v0.33 types. Falls back to the step-2 estimate if LM produces a negative radius.

The normal covariance matrix `M` can be maintained incrementally during BFS (add $w_i \, n_i \, n_i^T$ when accepting a face), but the eigenvector solve, circle fit, and LM refinement must be recomputed from scratch when re-fitting, since a change in axis direction invalidates the 2D projection. Steps 1-2 are O(n) in vertices with simple arithmetic; step 3 adds ~5-10 LM iterations (typically <1ms per hypothesis).

**Subtleties and edge cases:**
- **Near-parallel seed normals**: Pairwise qualification requires cross product magnitude > `MIN_CROSS_THRESHOLD` (e.g. 0.01, ~0.6° between normals) for each consecutive pair in the seed chain (fi↔n1 and n1↔n2). This ensures each pair contributes meaningful angular information. Note that fi and n2 need NOT have different normals — on a cylinder, faces at the same circumferential position but different axial positions have identical normals. The angular diversity comes from the intermediate face n1 having a different circumferential position.
- **Cones vs cylinders**: A cone has varying radius along the axis. If we fit a cylinder to conical faces, vertex distances near the ends will exceed tolerance, so cone faces will be naturally rejected. Cone fitting is deferred.
- **Large-radius cylinders**: A cylinder with very large radius looks locally planar. Such faces will belong to multi-face planar hypotheses from stage 2.1 (since vertices on a nearly-flat cylinder are nearly coplanar within tolerance). Multi-face planar faces CAN seed cylinders — the angular coverage validation, minimum face count, and other checks prevent bogus large-radius cylinder fits. Stage 2.6 resolves any ambiguity between cylindrical and planar assignments.
- **Cylinder wrapping >180°**: This works naturally with BFS. The convexity check ensures we don't accidentally merge inner and outer surfaces.
- **Multiple cylinders sharing an axis**: E.g. stepped bore holes with different radii. These will form separate hypotheses because vertices at different radii will exceed the vertex tolerance.
- **Seam between cylinder and plane**: Faces at the junction have some vertices on the plane and some on the cylinder. Since the plane vertices aren't on the cylinder surface, the vertex distance check rejects them. This correctly excludes boundary faces.

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the cylinder surface (nearest point: project onto axis, then move radially to radius distance).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - Note: centroid-to-surface distance measures both fitting accuracy and tessellation sagitta (the flat triangle's centroid lies inside/outside the curved surface by approximately $h^2 / 8R$ where $h$ is the chord length). For typical tessellations this is small relative to `--surface-tolerance`, but it could become significant for very coarse meshes on tight-radius cylinders.

#### 2.3 Deduce spherical hypotheses

Fit spherical hypotheses to connected sets of faces that lie on a common sphere. Faces with cylindrical hypotheses are not excluded — a face can legitimately belong to both a cylinder and a sphere (e.g., equator faces). Stage 2.6 resolves overlapping hypotheses later. A spherical hypothesis consists of:
- `center: [f64; 3]` — Center of the sphere.
- `radius: f64` — Radius of the sphere (always positive).
- `convex: bool` — Whether face normals point away from the center (convex, like the outside of a ball) or toward it (concave, like a bowl).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this sphere.
- `error_max: f64` — Maximum absolute distance from any vertex to the sphere surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

The algorithm uses BFS region growing, analogous to stages 2.1 and 2.2 but seeded from the group of faces surrounding a mesh vertex. A vertex neighborhood on a sphere naturally provides normals spanning 3D (unlike cylinder seeding which only needs normals spanning a 2D plane), making it a more natural seed geometry for the 4-DOF sphere fit. The algebraic sphere fit from 4+ non-coplanar vertices determines a unique sphere, and validation catches bad fits:

**Seeding:**
- For each mesh vertex V, collect the faces incident on V ("surrounding faces"). Exclude any face that has `spherical_hypothesis != UNDEDUCED` or that belongs to a multi-face planar hypothesis (same rationale as cylindrical seeding: their vertices lie on flat surfaces, not spheres). Multi-face planar faces remain UNDEDUCED so BFS can absorb them later.
    - If fewer than 3 non-excluded faces remain around V, skip this vertex (not enough angular diversity for a reliable sphere fit).
    - Compute a spherical fit to all vertices of the remaining surrounding faces using least-squares sphere fitting (see Sphere Fitting below).
    - Skip this vertex if:
        - The fitted radius exceeds `max_sphere_radius` (see below).
        - Any vertex of the seed faces is not within `--vertex-tolerance` of the estimated sphere surface.
        - The dihedral angle between any pair of adjacent seed faces (i.e., seed faces sharing a mesh edge) exceeds `--angular-tolerance` (default 17.5°).
        - The distance from any seed face centroid to the sphere exceeds `--surface-tolerance`.
    - Determine convexity: for each seed face, compute the dot product of the face normal with the vector from sphere center to face centroid. Positive → convex, negative → concave. If convexity is inconsistent across seed faces, skip this vertex.
    - Create a new spherical hypothesis, assign all seed faces, and perform BFS expansion.

*Design note: The original plan called for pair-based seeding (two adjacent faces with non-parallel normals, as in cylindrical), then 3-face seeds with a normals-span-3D eigenvalue check. Pair-based seeding failed on finely tessellated spheres where adjacent faces have nearly identical normals (spanning only ~6°). Vertex-neighborhood seeding is superior for spheres: it naturally collects faces spanning multiple directions around a vertex, providing better constraints for the 4-DOF sphere fit without needing an explicit normals-span check. The sphere fit itself validates geometry — if the vertices don't lie on a sphere, the fit produces a large radius or high error, caught by max_sphere_radius and vertex distance checks.*

**BFS expansion:**
- Pop faces from the BFS queue and examine each neighbor ni:
    - If ni already has a spherical hypothesis, skip.
    - **Vertex distance check**: for each vertex of ni, compute `| ||v - center|| - radius |`. If all are within `--vertex-tolerance`, proceed to remaining checks (no re-fit needed). If any exceeds `2 * vertex_tolerance`, skip (re-fitting cannot help, same reasoning as planar and cylindrical). If between 1x and 2x, flag for re-fit.
    - **Centroid validation**: compute the distance from ni's centroid to the fitted sphere. If beyond `2 * surface_tolerance`, skip. If beyond `surface_tolerance`, flag for re-fit. (Same structure as stage 2.2 cylindrical BFS.)
    - **Convexity check**: verify ni's normal agrees with the hypothesis convexity (dot product of normal with radial vector has the expected sign). If not, skip.
    - **Angular tolerance check**: for each of ni's mesh neighbors that is already assigned to this hypothesis, compute the dihedral angle between ni and that neighbor. If any exceeds `--angular-tolerance`, skip.
    - **Re-fit attempt**: if vertex distance or centroid validation flagged a re-fit, re-fit sphere from all current faces plus ni. If all vertices are within `--vertex-tolerance` and all face centroids within `--surface-tolerance` and radius within `max_sphere_radius`, accept the face with the updated sphere parameters. If not, revert to the previous sphere fit and skip the face.
    - If accepted: assign ni to the spherical hypothesis, and add its neighbors to the BFS queue.
- After BFS completes: final re-fit, compute error metrics.
- **Minimum face count**: require at least 4 faces (a sphere needs at least 4 non-coplanar points for a unique fit).

**Solid-angle coverage validation** (applied to the completed hypothesis, analogous to angular coverage for cylinders):

The problem: sphere BFS can grow along cylinder fillet surfaces, because locally adjacent cylinder faces fit on a sphere within vertex_tolerance. This produces oversized sphere hypotheses spanning a thin strip around a cylinder fillet (observed on `onshape_rounded_cube`). The analog of the cylindrical angular coverage test is a **solid-angle coverage test**: faces on a genuine spherical patch produce centroid-to-center directions that span 3D, while faces on a cylinder fillet produce directions that span only a 1D strip.

Algorithm:
1. For each face in the hypothesis, compute the unit vector from sphere center to face centroid: $d_i = \text{normalize}(\text{centroid}_i - \text{center})$.
2. Compute the area-weighted 3×3 covariance matrix of these unit vectors: $C = \sum w_i \, d_i \, d_i^T$ where $w_i$ is the face's mesh area.
3. Compute eigenvalues of $C$. Let $\lambda_1 \geq \lambda_2 \geq \lambda_3$.
4. Require $\lambda_3 / \lambda_1$ to exceed a minimum threshold (e.g., 0.01). If it fails, the hypothesis directions are nearly coplanar (a strip on the sphere, characteristic of cylinder fillet growth) or nearly collinear (a narrow band), and the hypothesis should be discarded.

*Note: Fillet faces lie along a 1D arc on the sphere's direction space, producing $\lambda_3 \approx 0$. Genuine spherical patches (hemispheres, domes, spherical pockets) have centroid directions spanning 3D, so all eigenvalues are substantial. For partial spherical caps, the directions span at least 2D; very small caps may fail but such small patches are unlikely to be meaningful spherical features.*

*Note: Like cylindrical angular coverage, this should NOT be applied during BFS growth — early in BFS, faces are clustered near the seed and naturally appear low-dimensional. It is a global structural validation for completed hypotheses.*

**Maximum radius**: `max_sphere_radius`. With solid-angle coverage validation and surface-tolerance validation during BFS, the radius limit no longer needs to be tight — those checks prevent the pathological growth patterns that the original 10× bounding box diagonal limit was designed to catch. A large value like `bounding_box_diagonal * 1000` suffices as a numerical guard rail, preventing sphere fits with absurd radii (where floating-point precision of vertex coordinates could make flat faces appear spherical) while not rejecting any plausible real-world spherical feature. Compute the bounding box once at stage 2 entry.

**Sphere fitting** (used for seeding and re-fitting):
- Given vertices on a sphere: $|v - c|^2 = r^2$, expand to $v \cdot v - 2 v \cdot c + |c|^2 = r^2$.
- Rearrange: $v \cdot v = 2 v \cdot c - |c|^2 + r^2 = 2 v \cdot c + k$ where $k = r^2 - |c|^2$.
- This is linear in unknowns $(c_x, c_y, c_z, k)$. Solve via least squares: $A x = b$ where $A_{i} = [2v_{ix}, 2v_{iy}, 2v_{iz}, 1]$ and $b_i = v_i \cdot v_i$.
- After solving: center $= (c_x, c_y, c_z)$, radius $= \sqrt{k + |c|^2}$.

**Distinguishing spheres from cylinders:**
- Faces with cylindrical hypotheses are NOT excluded — they may legitimately belong to both a cylinder and a sphere (e.g., equator band). Stage 2.6 disambiguates.
- On cylindrical surfaces, sphere fits from vertex neighborhoods with nearly-parallel normals tend to produce unreliable results (very large or negative $r^2$), which are caught by `max_sphere_radius` and the $r^2 > 0$ check in sphere fitting. The solid-angle coverage test provides additional protection: faces along a cylinder strip span only 1D in direction space, failing the eigenvalue ratio test.
- Vertex-based seeding with the minimum-3-face requirement and angular diversity naturally avoids seeds from planar regions while allowing seeds from curved regions.

**Preventing false positives on chamfers/bevels:**
- On a chamfered or beveled cube, corner and edge chamfer faces can seed sphere fits. However, these faces are flat — the sphere fit produces a degenerate radius vastly larger than the model. The `max_sphere_radius` constraint catches this even at 1000× bounding box diagonal: for a 14mm model with vertex_tolerance=1e-5, a degenerate sphere fit would have R ≈ 112 km, far exceeding `14mm × 1000 = 14m`.
- The centroid validation check also catches spurious fits where the algebraic fit succeeds but the sphere doesn't actually pass through the face centroids.
- The solid-angle coverage test may also reject chamfer fits where few faces span a limited set of directions.

**Test models:**
- Simple sphere (ccad): full convex sphere, radius 10 → 1 convex spherical hypothesis.
- Hemisphere (ccad): top half of sphere + flat base → 1 convex spherical hypothesis + 1 multi-face planar.
- Spherical pocket (ccad): block with concave hemispherical cavity → 6 planar + 1 concave spherical.
- Ball on cylinder (ccad): sphere atop cylinder stalk → 1 convex spherical + 1 convex cylindrical + 1 planar.
- Sphere (onshape): full sphere, fine tessellation.
- Dome/hemisphere (onshape): hemisphere with flat base.

- With the `--compare` flag:
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the fitted sphere surface (nearest point: normalize the centroid-to-center vector and scale to radius distance).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - Note: same centroid-to-surface sagitta caveat as for cylindrical hypotheses (stage 2.2). For spheres the sagitta is approximately $h^2 / 8R$ where $h$ is the chord length of the facet.

*Comment: Spheres are 4 DOF (center + radius). Minimum 4 non-coplanar points for a unique fit. For "negative curvature" (concave spherical patch), track normal orientation relative to center rather than using negative radius. ~~Key challenge: partial spherical patches are hard to distinguish from cylinders or even planes if the patch is small relative to the radius. Consider requiring a minimum angular extent or using curvature analysis to disambiguate.~~ Addressed: solid-angle coverage validation prevents 1D strip growth (cylinder fillets), surface-tolerance validation in BFS prevents growth along non-spherical surfaces, and max_sphere_radius prevents degenerate large-radius fits on flat faces. Also consider toroidal surfaces (donuts, fillets)—they're common in CAD and combine characteristics of cylinders and spheres. Note: torus faces locally fit spheres, producing sphere hypotheses on torus surfaces (e.g., pipe_elbow). Stage 2.6 will need torus hypothesis support to resolve this.*

#### 2.4 Deduce conical hypotheses

Fit conical hypotheses to connected sets of faces that lie on a common cone. A cone is a surface of revolution (SOR) with a linear profile — it generalizes the cylinder (which has a constant-radius profile). Conical surfaces commonly appear as chamfers, tapered bores, and funnel shapes in CAD models.

A conical hypothesis consists of:
- `apex: [f64; 3]` — Apex (tip) of the cone.
- `axis_direction: [f64; 3]` — Unit direction vector along the axis (from apex outward toward the open end).
- `half_angle: f64` — Half-angle of the cone (in radians). The angle between the axis and any generator line on the surface.
- `convex: bool` — Whether face normals point away from the axis (convex, like the outside of a cone) or toward it (concave, like the inside of a conical hole).
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this cone.
- `error_max: f64` — Maximum absolute distance from any vertex to the cone surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

**Cone parameterization:**
A cone with apex $A$, axis direction $\hat{a}$, and half-angle $\theta$ satisfies: for any point $p$ on the surface, the angle between $(p - A)$ and $\hat{a}$ equals $\theta$. Equivalently, in (h,r) coordinates where $h = (p - A) \cdot \hat{a}$ (signed axial distance from apex) and $r = |(p - A) - h\hat{a}|$ (radial distance from axis), the surface satisfies $r = h \tan\theta$.

The signed distance from a point $p$ to the cone surface (for LM residual) is:
$$d = r_i \cos\theta - h_i \sin\theta$$
where $h_i = (p_i - A) \cdot \hat{a}$ and $r_i = |(p_i - A) - h_i \hat{a}|$.

**Initialization via (h,r) profile fitting:**
This exploits the key insight from the unified SOR framework (see Algorithm Reference): all surfaces of revolution reduce to a 2D profile fitting problem once the axis is estimated.

1. **Axis estimation** from normal covariance (same as cylinder): compute weighted covariance $M = \sum w_i n_i n_i^T$ where $n_i$ are face normals and $w_i$ are face areas. The eigenvector with smallest eigenvalue is the axis direction estimate.
2. **Reduce to (h,r) coordinates**: $h_i = (v_i - \text{centroid}) \cdot \hat{a}$, $r_i = |(v_i - \text{centroid}) - h_i \hat{a}|$.
3. **Fit linear profile** $r = m \cdot h + b$ via weighted linear regression (weights = vertex area).
4. **Extract cone parameters**: half-angle $\theta = \arctan(m)$, apex at distance $-b/m$ from centroid along axis.
5. **Reject if cylindrical**: if $|m| < 0.01$ (nearly constant radius), this is a cylinder, not a cone. Also reject if $\theta > 85°$ (nearly planar).

**Alternative initialization via quadric fitting:**
A cone is a quadric surface: $x^T A x + b^T x + c = 0$ with specific eigenvalue structure. This approach is useful as a cross-check or when the axis is poorly estimated.

1. Build the 10-coefficient quadric feature matrix: for each vertex $p_i = (x,y,z)$, construct row $[x^2, y^2, z^2, xy, xz, yz, x, y, z, 1]$.
2. SVD of the feature matrix. The right singular vector with smallest singular value gives the quadric coefficients.
3. Convert to matrix form: $Q(x) = x^T \mathbf{A} x + b^T x + c$.
4. Extract cone parameters:
   - Apex: $A = -\frac{1}{2} \mathbf{A}^{-1} b$ (gradient of Q is zero at apex).
   - Axis: eigenvector of $\mathbf{A}$ corresponding to the distinct eigenvalue (for a cone, two eigenvalues are equal and one is different).
   - Half-angle: from the eigenvalue ratio.
5. Validate: eigenvalue pattern must match a cone (two eigenvalues approximately equal, one distinct, with specific sign relationships). If the pattern matches a cylinder (one eigenvalue ≈ 0), sphere (all equal), or ellipsoid, skip.

**BFS region growing** (analogous to cylindrical):
- Seed from triples of faces (like cylinders) where normals make a roughly constant angle with the estimated axis, but the radial distance from the axis varies monotonically along the axis.
- Expand via BFS: for each candidate face, check $|r_i - h_i \tan\theta| < \text{vertex\_tolerance}$ for all vertices.
- Centroid validation: face centroid within surface_tolerance of the cone.
- Convexity check: normal direction consistent with hypothesis convexity.
- Angular tolerance: dihedral angle between adjacent faces ≤ angular_tolerance.
- Angular coverage validation (same algorithm as cylinder).
- Minimum 6-face requirement (a cone has 5 independent DOF but 6 faces provides margin).

**Distinguishing cones from cylinders:**
After fitting both a cylinder and a cone to the same face set, compare RMS residuals. If the cone fit is not significantly better than the cylinder fit (< 20% RMS improvement), prefer the simpler cylinder interpretation. This prevents noisy flat-profile data from being misclassified as a very gentle cone.

**LM refinement:**
Optimize $[A_x, A_y, A_z, \alpha, \beta, \theta]$ (6 parameters: apex position 3, axis tilt angles 2, half-angle 1).
Residual: $d_i = r_i \cos\theta - h_i \sin\theta$.

- With the `--compare` flag:
    - Same centroid-projection validation as cylinder/sphere.

*Comment: The (h,r) profile approach is elegant because it handles the cylinder/cone ambiguity naturally: a horizontal line in (h,r) is a cylinder, a sloped line is a cone, and the distinction is just the slope parameter. This is described in the SOR reduction technique. The quadric approach provides an independent validation path.*

#### 2.5 Deduce toroidal hypotheses

Fit toroidal hypotheses to connected sets of faces that lie on a common torus. Toroidal surfaces commonly appear as fillets (rounded edges), blends, and pipe elbows in CAD models. A torus is a surface of revolution with a circular profile offset from the axis.

A toroidal hypothesis consists of:
- `center: [f64; 3]` — Center of the torus (on the axis, in the plane of the major circle).
- `axis_direction: [f64; 3]` — Unit direction vector along the torus axis.
- `major_radius: f64` — Distance from center to the tube centerline ($R$).
- `minor_radius: f64` — Radius of the tube cross-section ($r$). For a fillet of radius $r$, the minor radius is $r$.
- `convex: bool` — Whether the torus patch is on the outer (convex) or inner (concave) surface of the tube.
- `faces: Vec<usize>` — Mesh face indices that fit this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this torus.
- `error_max: f64` — Maximum absolute distance from any vertex to the torus surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

**Torus parameterization:**
A torus with center $C$, axis $\hat{a}$, major radius $R$, minor radius $r$ satisfies: for any point $p$, let $v = p - C$, $\text{axial} = v \cdot \hat{a}$, $\text{radial} = |v - \text{axial} \cdot \hat{a}|$. Then:
$$d = \sqrt{(\text{radial} - R)^2 + \text{axial}^2} - r$$
is the signed distance from $p$ to the torus surface (0 on surface, positive outside tube, negative inside).

**Fitting via medial axis / tube center method:**

Direct nonlinear fitting of 7 torus parameters ($C$, $\hat{a}$, $R$, $r$) is numerically fragile because of the high DOF and the coupling between parameters. The **medial axis approach** linearizes the problem by exploiting the geometric relationship between surface normals and the tube center.

Key insight: for a surface with constant positive principal curvature $\kappa_1 = 1/r$ in one direction, the **tube center** (center of the osculating circle) at each point is $k_i = p_i + r \cdot n_i$ where $n_i$ is the outward surface normal. For a torus, all tube centers lie on the **major circle** — a circle of radius $R$ centered at $C$ in the plane perpendicular to $\hat{a}$. This reduces torus fitting to: (1) estimate $r$, (2) compute tube centers, (3) fit a 3D circle.

**Step 1: Estimate minor radius $r$ from normal-line intersections.**
- For pairs of adjacent (non-parallel) faces in the candidate region, compute the closest approach of their normal lines: $L_1: p_1 + t \cdot n_1$, $L_2: p_2 + s \cdot n_2$.
- Reject pairs with nearly parallel normals ($|n_1 \cdot n_2| > 0.95$).
- The closest-approach distance approximates 0 (for a torus, normal lines converge at the tube center), and the intersection parameter gives $r \approx t \approx s$ (for outward normals on a convex torus).
- Collect many estimates, take the median for robustness against outliers.

**Step 2: Compute tube center points.**
- For each vertex $p_i$ with outward normal $n_i$: $k_i = p_i + r \cdot n_i$ (convex torus) or $k_i = p_i - r \cdot n_i$ (concave torus).
- These points should cluster on the major circle.

**Step 3: Classify spine geometry (detect degeneracies).**
- Compute PCA of the $k_i$ points. Eigenvalue pattern:
  - $\lambda_1 \approx \lambda_2 \gg \lambda_3$: points form a circle → **torus** (proceed).
  - $\lambda_1 \gg \lambda_2 \approx \lambda_3$: points form a line → actually a **cylinder** (delegate to cylinder detection).
  - $\lambda_1 \approx \lambda_2 \approx \lambda_3 \approx 0$: points cluster at a single point → actually a **sphere** (delegate to sphere detection).

**Step 4: Fit 3D circle to tube center points.**
- PCA of $k_i$: smallest eigenvector = circle plane normal = **torus axis direction** $\hat{a}$.
- Project $k_i$ into the circle plane using an orthonormal basis $(u, v)$ perpendicular to $\hat{a}$.
- 2D algebraic circle fit in the projected plane → center $(u_c, v_c)$ and radius $R$.
- Convert back to 3D: torus center $C = \text{mean}(k_i) + u_c \cdot u + v_c \cdot v$, major radius $R$.

**Step 5: LM refinement.**
- Optimize $[C_x, C_y, C_z, \alpha, \beta, R, r]$ (7 parameters: center 3, axis tilt 2, major radius 1, minor radius 1).
- Residual: $d_i = \sqrt{(\text{radial}_i - R)^2 + \text{axial}_i^2} - r$.

**Seeding strategy:**
Toroidal fillets occur at edges between two surfaces where a blend or fillet has been applied. The fillet faces have normals that rotate smoothly as you traverse the fillet. Seeding uses:

- **Edge-based seeding**: For each boundary between two multi-face selected surfaces (from stage 2.6), examine the single-face planar hypothesis faces adjacent to that boundary. These "leftover" faces between two fitted surfaces are strong candidates for fillet/blend toroidal surfaces. Start BFS from clusters of these faces.
- **Normal-rotation detection**: In a seed neighborhood, if normals rotate systematically around the tube cross-section (as measured by normal covariance: the $\lambda_1 \approx \lambda_2 \gg \lambda_3$ eigenvalue pattern, similar to cylinder but with additional normal variation), this suggests a toroidal surface.

**BFS region growing:**
- Vertex-to-torus distance: $|\sqrt{(\text{radial} - R)^2 + \text{axial}^2} - r| < \text{vertex\_tolerance}$.
- Centroid validation: face centroid within surface_tolerance of torus.
- Convexity check: normal direction consistent with hypothesis.
- Angular tolerance: dihedral angle ≤ angular_tolerance.
- Re-fit on marginal faces (same as cylinder/sphere).

**Post-growth validation:**
- Minimum 7-face requirement (torus has 7 DOF).
- Angular coverage around the minor circle: faces should span multiple angular positions around the tube cross-section (analogous to cylinder angular coverage), ensuring the torus fit is well-constrained.

*Comment: The medial axis / tube center method is far more robust for fillet-sized torus patches than direct 7-parameter nonlinear fitting. By reducing the problem to (1) estimate one scalar r, (2) fit a 3D circle to the tube centers, it separates the easy part (circle fitting is linear algebra) from the hard part (estimating the minor radius). The spine classification step (Step 3) also provides a natural way to detect degenerate cases where the "torus" is actually a cylinder or sphere — this happens when the fillet radius approaches the edge length or when the surface is nearly spherical.*

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

This section documents key algorithmic techniques drawn from computational geometry research that inform the surface fitting approach. These are referenced by the stage descriptions above.

#### Unified Surface-of-Revolution (SOR) Framework

All surfaces of revolution — cylinder, cone, sphere, torus — share a common structure: they are generated by rotating a 2D profile curve around an axis. This means they can all be fitted using a unified approach:

1. **Estimate the axis of revolution** from the surface normals: for any SOR, surface normals are perpendicular to the axis (cylinder), make a constant angle with the axis (cone), point toward the axis (sphere pole behavior), or some combination. The axis can be estimated as the smallest eigenvector of the area-weighted normal covariance matrix $M = \sum w_i n_i n_i^T$.

2. **Reduce to (h,r) coordinates**: for each vertex, compute $h = (p - \text{centroid}) \cdot \hat{a}$ (axial position) and $r = |(p - \text{centroid}) - h\hat{a}|$ (radial distance from axis).

3. **Fit candidate profile curves** in the (h,r) plane:
   - **Cylinder**: $r = R$ (constant) — horizontal line. 1 parameter: $R$.
   - **Cone**: $r = m \cdot h + b$ (linear) — sloped line. 2 parameters: $m, b$.
   - **Sphere**: $h^2 + (r - 0)^2 = R^2$ — circle centered on the axis. 2 parameters: $h_0, R$ (center offset and radius). But note: sphere center is ON the axis, so this is $h^2 + r^2 = R^2$ in centered coordinates.
   - **Torus**: $(r - R)^2 + h^2 = r_{\text{minor}}^2$ — circle offset from axis. 2 parameters: $R, r_{\text{minor}}$.

4. **Pick the best profile** by RMS residual. Use BIC or AIC to penalize model complexity (cylinder has 1 DOF, cone has 2, torus has 2).

This framework naturally handles ambiguous cases: a very gentle cone looks like a cylinder (slope $\approx 0$); a sphere with very large radius looks like a plane ($R \to \infty$); a torus with large major radius looks like a cylinder ($R \gg r$).

The current implementation uses separate fitting algorithms for each surface type (stages 2.2–2.5), which is appropriate for incremental development. A future refactoring pass (see Implementation Phases) will unify these into the SOR framework to improve accuracy on ambiguous cases and simplify the codebase.

#### Gaussian Map (Normal Covariance) Classification

The distribution of surface normals (the "Gaussian map") directly reveals the surface type. For a patch of mesh faces, compute the area-weighted covariance matrix of the normals:
$$C = \sum_i w_i (n_i - \bar{n})(n_i - \bar{n})^T$$
where $w_i$ is face area and $\bar{n}$ is the weighted mean normal. The eigenvalues $\lambda_1 \geq \lambda_2 \geq \lambda_3$ classify the surface:

- **Plane**: $\lambda_1 \approx \lambda_2 \approx \lambda_3 \approx 0$ — all normals identical.
- **Cylinder/Cone**: $\lambda_1 \approx \lambda_2 \gg \lambda_3$ — normals lie on a great circle (cylinder) or small circle (cone).
- **Sphere**: $\lambda_1 \approx \lambda_2 \approx \lambda_3$ (all nonzero) — normals distributed isotropically.
- **Torus**: $\lambda_1 > \lambda_2 \gg \lambda_3$ — normals form a band pattern.

The cylinder axis estimation already uses this: the axis is the eigenvector of the uncentered normal covariance $M = \sum w_i n_i n_i^T$ corresponding to the smallest eigenvalue (the direction perpendicular to all normals). Extending this to pre-classify surface type before fitting would improve seeding efficiency.

#### Medial Axis / Tube Center Method

For any surface with a constant principal curvature $\kappa_1 = 1/r$ in one direction, the **medial axis** (locus of centers of curvature) can be computed directly from the surface points and normals:
$$k_i = p_i + r \cdot n_i$$

The geometry of the medial axis reveals the surface type:
- **Cylinder**: medial axis is a line (the cylinder axis) — $k_i$ points are collinear.
- **Sphere**: medial axis is a point (the center) — $k_i$ points cluster.
- **Torus**: medial axis is a circle (the major circle) — $k_i$ points lie on a circle.
- **Cone**: medial axis is a point (the apex) for one direction of curvature.

This is used in stage 2.5 (torus fitting): estimate $r$ from normal-line intersections between adjacent faces, compute tube centers $k_i$, fit a 3D circle to get the major circle, and extract the torus parameters.

#### Quadric Fitting

Any quadric surface (plane, sphere, cylinder, cone, ellipsoid, hyperboloid, paraboloid) can be written as:
$$Q(x, y, z) = Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0$$

Fitting: construct the $N \times 10$ feature matrix and find its null space via SVD (smallest singular vector). The matrix form $x^T \mathbf{A} x + b^T x + c = 0$ encodes the surface type in the eigenvalues of $\mathbf{A}$:
- **Sphere**: $\mathbf{A} = \lambda I$ (three equal eigenvalues).
- **Cylinder**: one eigenvalue is 0 (axis direction), other two equal.
- **Cone**: three nonzero eigenvalues with $\det(\mathbf{A}) \approx 0$ (passes through apex).
- **Plane**: $\mathbf{A} = 0$ (all eigenvalues zero).

Quadric fitting is fast (SVD of an $N \times 10$ matrix, $O(N)$) and provides a good initialization for specialized fitters. It's used as an alternative initialization for cone fitting (stage 2.4).


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
  - **Plane-torus** (TODO): `Geom_Circle` — the torus-plane tangent curve is a circle on the torus at the minor boundary (inner or outer equator). Center on the major circle, radius = R ± r.
  - **Cylinder-torus** (TODO): `Geom_Circle` at the junction where the torus minor circle meets the cylinder. Commonly occurs at fillet-to-cylinder transitions.
  - **Torus-torus** (TODO): `Geom_Circle` where two fillet torus patches meet (e.g., at rounded corners where 3 fillets converge). The shared circle lies on both tori.
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
  - **Conical faces** (TODO): Same three approaches as cylinder (full revolution via UV-bounds, partial via wire+pcurves). The cone apex vertex is a UV singularity (all U values collapse to V=0 or V_apex); handle analogously to sphere poles. For truncated cones (no apex), UV-bounds construction suffices.
  - **Toroidal faces** (TODO): UV-bounds construction for complete torus patches. The torus UV space is doubly periodic: U ∈ [0,2π) around the major axis, V ∈ [0,2π) around the tube. Fillet patches (partial torus) use wire-based construction with pcurves mapping 3D curves to torus UV space. The V range of a fillet is typically a small arc (e.g., 0 to π/2 for a 90° fillet).

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
- [x] Implement cone fitting: axis estimation from normal covariance, (h,r) profile fitting via linear regression to get half-angle and apex position. Apex-vertex seeding strategy instead of triple-seed: identifies vertices with many non-coplanar-normal incident faces (eigenval/trace < 0.3), uses all incident faces as seed set for robust axis estimation. Axis orientation ensured by flipping when h_sum < 0. Levenberg-Marquardt refinement for 6-parameter optimization (axis perturbation, apex position, half-angle). False-positive rejection filters: half-angle bounds (2°–85°), vertex error < surface_tol, apex distance < 10× bounding box diagonal, normal-axis angle consistency (std_dev < angular_tol/2).
- [x] Implement cone BFS region growing with vertex-to-cone distance validation and angular coverage check. Seeds validated with relaxed tolerance (5× surface_tol), BFS expansion with vertex+centroid distance checks, convexity by majority vote after final re-fit. Angular tolerance for cone BFS uses 2× the configured value to accommodate inter-strip dihedral angles on coarsely tessellated cones (e.g., 18° between adjacent strips with 20 circumferential divisions).
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
- [ ] Unit tests: torus fitting on synthetic point sets; BFS growing on fillet meshes.
- [ ] Full pipeline tests: filleted_cylinder, filleted_hole_block, filleted_pipe pass `--compare` through stage 4.1.
- [ ] Unblock existing `pipe_elbow_10_fine` (onshape) test.

### Stage 3 Extensions for Cone and Torus
- [x] Stage 3.1: Add `Geom_ToroidalSurface` creation from hypothesis parameters. Uses `geom::ToroidalSurface::new_ax3_real2(&ax3, major_radius, minor_radius)` binding (already present in opencascade-sys). Unblocked `rounded_cube_coarse_brep_check` test.
- [ ] Stage 3.2: Add torus analytical normal formulas for tangency detection. (Cone normals already implemented.)
- [ ] Stage 3.3: Implement tangent edge curve computation for plane-cone, cylinder-cone, plane-torus, cylinder-torus, and torus-torus pairs.
- [ ] Stage 3.3: Handle non-tangent cone and torus edge intersections via `GeomAPI_IntSS` with curve selection.
- [ ] Stage 3.4: Implement conical face construction (UV-bounds for full revolution, wire+pcurves for partial, apex singularity handling).
- [ ] Stage 3.4: Implement toroidal face construction (UV-bounds for full revolution, wire+pcurves for partial fillet patches, doubly-periodic UV handling).
- [ ] Full pipeline tests: all cone and torus test models pass --compare through stage 4.1.

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

