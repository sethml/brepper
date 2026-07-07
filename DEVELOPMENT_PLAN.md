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
        - `conical_hypothesis: i32` — Index of active conical/cylindrical hypothesis, or `NO_HYPOTHESIS` (-1) if none, or `UNDEDUCED` (-2) if not yet deduced. Cylinders and cones are unified: both use this field with `ConicalHypothesis` (cylinders have `half_angle == 0`).
        - `spherical_hypothesis: i32` — Index of active spherical hypothesis, or `NO_HYPOTHESIS` (-1) if none, or `UNDEDUCED` (-2) if not yet deduced.
    - `stats: MeshValidationStats` — Statistics described in stage 1.2.

Note: hypotheses vectors (planar, conical, spherical, toroidal) are stored in `Stage2Output`, not in `ConnectedMesh`. The mesh stores only per-face hypothesis indices. Cylinders and cones share the same hypothesis vector (`conical_hypotheses`) and per-face field (`conical_hypothesis`), with `half_angle == 0` distinguishing cylinders from true cones.

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
4. **Inverted normals**: for each closed shell, compute the signed volume via the divergence theorem (sum of signed tetrahedra formed by each face and the origin: V = (1/6) Σ v0·(v1×v2)). A negative signed volume is legitimate for an internal void shell (normals point out of the material, i.e. into the cavity): such shells are accepted — and counted in `voids_within_solids` — when their bounding box is contained in some closed positive shell's box. Otherwise a negative volume means all face normals point inward (consistently inverted winding); this catches meshes where orientation is internally consistent (check 3 passes) but globally inverted.
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

**Precision note**: Mesh vertex positions, normals, and all hypothesis parameters use double-precision (`f64`) throughout. Surface fitting arithmetic and BFS distance computations are all in double precision.

#### Generic BFS Hypothesis Generation

Stages 2.1–2.5 share a common algorithm structure for deducing surface hypotheses from the mesh. This section describes the generic algorithm; each stage section below describes how it is parameterized. The Rust implementation uses the `SurfaceType` trait in `bfs.rs` to share code across surface types.

The algorithm operates on a single surface type at a time (planar, conical/cylindrical, spherical, toroidal). For each surface type, the algorithm iterates over all faces that don't yet have a hypothesis of that type, attempting to grow a hypothesis from each face.

**Hypothesis data structure:**
Each surface type defines a hypothesis struct with type-specific parameters (see per-stage sections) plus common fields:
- `faces: Vec<usize>` — Mesh face indices in this hypothesis.
- `vertices: Vec<usize>` — Mesh vertex indices on this surface.
- `error_max: f64` — Maximum absolute distance from any vertex to the fitted surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

Each `MeshFace` has a hypothesis index field per surface type (e.g., `planar_hypothesis`, `conical_hypothesis`, `spherical_hypothesis`, `toroidal_hypothesis`), initialized to `UNDEDUCED` (-2) at the start of each stage. Conical and cylindrical surfaces share the `conical_hypothesis` field (cylinders have `half_angle == 0`).

**Algorithm: `generate_hypotheses<S: SurfaceType>(mesh, config)` — runs the full hypothesis generation for one surface type:**

1. Initialize all hypothesis fields for this surface type to `UNDEDUCED`.
2. For each face `fi` in index order where the face is eligible (hypothesis field == `UNDEDUCED`, plus any surface-type-specific eligibility checks):
    a. **Seed generation**: call `S::generate_seeds()` to find all valid seeds anchored at `fi` (see per-stage sections for seed types). For each candidate seed:
        - Collect the seed faces' vertices.
        - Perform a surface fit from the seed vertices using `S::fit_surface()`.
        - **Seed vertex tolerance**: check all seed face vertices within `--vertex-tolerance` of the fitted surface. Reject seed if any vertex exceeds this.
        - **Seed surface tolerance**: check all seed face centroids within `--surface-tolerance` of the fitted surface. Reject seed if any fails.
        - **Seed angular tolerance**: check dihedral angles between all pairs of adjacent seed faces. Reject seed if any exceeds `--angular-tolerance`.
        - Determine convexity from the seed faces.
    b. **BFS expansion** (per valid seed): using temporary data structures (local face list, vertex set — mesh hypothesis fields are NOT mutated during trials):
        - Initialize BFS queue with all seed faces.
        - Pop faces from the queue. For each neighbor `ni` of the popped face:
            - Skip if `ni` is not eligible for this surface type (already assigned a hypothesis, or fails surface-type-specific eligibility).
            - Skip if `ni` is already claimed by this trial.
            - **Angular tolerance check**: for each of `ni`'s mesh neighbors that is already claimed by this trial, compute the dihedral angle. If any exceeds `--angular-tolerance`, skip `ni`.
            - **Convexity check**: call `S::convexity_check()`. If it fails, skip `ni`.
            - **Surface distance check**: compute distance from each vertex of `ni` to the fitted surface via `S::vertex_distance()`, and from `ni`'s centroid via `S::centroid_distance()`.
                - If all vertices are within `--vertex-tolerance` and centroid within `--surface-tolerance`: accept directly.
                - If any vertex exceeds `REFIT_SKIP_MULTIPLIER * vertex_tolerance` (currently 10.0×) or centroid exceeds `REFIT_SKIP_MULTIPLIER * surface_tolerance`: skip `ni` — too far for refitting to help.
                - Otherwise (between 1× and `REFIT_SKIP_MULTIPLIER`×): attempt a **re-fit** using all current trial faces plus `ni`. If the re-fit satisfies vertex tolerance and surface tolerance for **all** faces (existing and new), accept the re-fit. Otherwise, revert and skip `ni`.
            - If accepted: add `ni` to the trial's face list, add its vertices to the vertex set, push `ni` onto the BFS queue.
        - After BFS completes: do a final re-fit from all accumulated faces and vertices.
    c. **Acceptance tests** (per completed trial): call `S::acceptance_test()` for surface-type-specific post-BFS validation (e.g., angular coverage, solid-angle coverage, minimum face count). If validation fails, discard the trial.
    d. **Best candidate selection**: among all accepted trials for this starting face `fi`, keep the one with the largest total face area. Only two candidates need to be in memory at a time (current best and trial being evaluated).
    e. **Commit**: if a best candidate was found, create a hypothesis and assign it to all of its faces (mutating `mesh.faces[].{type}_hypothesis`). Faces assigned here will be skipped by subsequent starting faces.
3. After all starting faces have been processed: set any remaining `UNDEDUCED` faces to `NO_HYPOTHESIS`.

**Curvature-guided seed ordering**: Start faces are processed in ascending order of a per-face *curvature-difference score* rather than mesh index order. The score is computed once per mesh (in stage 2, shared by all surface stages):

1. **Per-face curvature summary**: for each mesh edge between adjacent faces, the discrete normal curvature along the crossing direction is κ ≈ θ/d (dihedral angle over centroid distance). A face's summary is the maximum κ over its edges — approximately its largest principal curvature.
2. **Neighbor difference**: a face's raw score is the maximum |summary difference| against its edge neighbors. Interior faces of any uniform surface (plane, cylinder, cone, sphere, torus) score ≈ 0 — curvature is locally constant even where it is nonzero — while faces at surface transitions score high. Crucially, this includes **tangent junctions** (fillet↔plane, fillet↔cylinder), which are G1-continuous and therefore invisible to dihedral-angle checks, but G2-discontinuous: curvature jumps from 0 to 1/r.
3. **Diffusion**: 3 Jacobi passes of neighbor averaging (new = ½·self + ½·mean(neighbors)) spread the junction signal a few hops, so ordering prefers deep-interior faces whose seed neighborhoods (1-ring / fan, radius ≈ 1–2) cannot touch a junction.
4. **Junction mask**: separately from the diffused ordering score, faces whose RAW score exceeds 50% of the local curvature magnitude (max of own and neighbors' summaries) are flagged as junction faces. The criterion is local-relative, not mesh-global: sharp creases dominate a global maximum and tangent junctions — the primary target — would fall below any fraction of it.

Effect: interior seeds run first and commit large, correctly-fitted hypotheses; junction faces are usually absorbed by BFS growth before they ever become start faces (committed faces do not seed). **Seed gathering exclusion**: in pass 1, junction faces are also excluded from seed neighborhoods (conical fans, spherical 1-rings, toroidal hop-expansion), so seeds stay on one surface; a second pass retries remaining undeduced faces without the filter so surfaces composed entirely of transition faces (one-row fillets) still seed. Toroidal gathering hop-expands (up to 3 hops, targeting ≥ 16 faces) since small fillet tubes have little interior across the tube but extend along the fillet axis. Ordering and masks are deterministic: ties break by face index.

**Computational cost**: For each starting face, the number of seed trials is typically small (1 for planar/spherical/toroidal, 1–6 for conical/cylindrical). Each trial BFS covers the connected surface region. Since most faces on the same surface will discover the same region, the work is dominated by the first few faces that discover the region before their neighbors get committed. The curvature field costs O(N) plus three O(N) diffusion passes — negligible next to LM fitting.

- With the `--compare` flag (applied to all surface types):
  - For each hypothesis, compute the centroid of each member face.
  - Project each centroid onto the fitted surface (nearest point on the surface).
  - Measure the distance from the projected point to the nearest surface in the reference STEP file using `BRepExtrema_DistShapeShape`.
  - If any projected centroid exceeds `--surface-tolerance`, report an error.
  - Note: centroid-to-surface distance measures both fitting accuracy and tessellation sagitta ($\approx h^2 / 8R$ for curvature radius $R$ and chord length $h$). For typical tessellations this is small relative to `--surface-tolerance`.

#### 2.1 Deduce planar hypotheses

Fit planar hypotheses to all connected sets of coplanar faces. Every face is assigned to exactly one planar hypothesis (including single-face hypotheses for faces on curved surfaces). A planar hypothesis consists of:
- `normal: [f64; 3]` — Unit normal vector pointing outward from the shell/solid.
- `distance: f64` — Signed distance from origin to plane along the normal (plane equation: normal · p = distance).
- `error_max, error_min: f64` — Maximum (positive) and minimum (most negative) signed distance from any vertex to the plane.
- `error_abs_sum: f64` — Sum of absolute vertex-to-plane distances.

Runs the generic BFS hypothesis generation with the following parameterization:

**Seed type**: The starting face `fi` alone. No additional seed faces are needed — a single face defines a plane.

**Surface fit**: Area-weighted normal averaging. Each face's unit normal is weighted by the face's triangle area, then the weighted sum is normalized. The plane distance is the mean of `normal · vertex` over all vertices. This gives larger faces more influence on the plane orientation. See [SURFACE_FITTING.md](SURFACE_FITTING.md) for detailed fitting math.

**Distance metric**: Signed distance from vertex to plane: `d = normal · vertex - distance`. All vertices must satisfy `|d| ≤ vertex_tolerance`.

Because the vertex tolerance is tight for CAD meshes, the angular tolerance check is redundant for planar hypotheses — any face with all vertices in the tolerance band necessarily has a closely aligned normal. It is still applied for consistency with the generic algorithm.

**Acceptance tests**: `face_count >= 2`. Single-face planes are not kept as hypotheses. Faces that end up in a single-face group are left as `NO_HYPOTHESIS`; stage 2.6 synthesizes single-face planar hypotheses for any faces that received no hypothesis from stages 2.1–2.5.

**Execution order note**: Stages execute planar (2.1) → spherical (2.3) →
toroidal (2.5) → conical (2.2) → selection (2.6), keeping the historical
stage numbering for CLI/viz compatibility. The principle is most-constrained
surfaces first: every band of quad rows on a tessellated sphere lies exactly
on some cone or cylinder, and any adjacent pair of quad rows on a torus tube
likewise lies exactly on some cone — so if the conical stage ran first it
would tile spheres and fillets with bogus band hypotheses (slow, and they
starve the later stages' seeding). Sphere and torus acceptance tests are
selective enough to run early (solid-angle coverage; normal-based tube-span
coverage), and cylinders resist bogus torus fits cheaply (tube centers
collapse onto the axis; the (h,r) profile is a constant-radius line).
Conical seeding skips sphere- and torus-committed faces; conical BFS
*expansion* may still absorb them, so tangent junctions work.

#### 2.2 Deduce conical/cylindrical hypotheses

Fit conical and cylindrical hypotheses to connected sets of faces that lie on a common cone or cylinder. Cylinders and cones are unified into a single `ConicalHypothesis` struct (matching the C++ implementation), with `half_angle == 0` for cylinders. A hypothesis consists of:
- `apex: [f64; 3]` — The apex of the cone, or a point on the axis for cylinders.
- `axis_direction: [f64; 3]` — Unit direction vector along the axis.
- `half_angle: f64` — Half-angle of the cone in radians (0 for cylinders).
- `radius_at_apex: f64` — Cylinder radius when `half_angle == 0`. 0 for true cones.
- `convex: bool` — Whether face normals point away from the axis (convex, like the outside of a pipe) or toward it (concave, like the inside of a hole).
- `error_max: f64` — Maximum absolute distance from any vertex to the surface.
- `centroid_error_max: f64` — Maximum absolute distance from any face centroid to the surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

The `is_cylinder()` method returns `true` when `half_angle < 1e-4`, used throughout to branch on cylinder vs. cone behavior.

Runs the generic BFS hypothesis generation with the following parameterization:

**Seed type**: Three faces `(fi, f1, f2)` — the starting face `fi`, a pairwise-qualified neighbor `f1` of `fi`, and another neighbor `f2` of `fi` where `f2` is NOT an edge-wise neighbor of `f1` (so `fi` is between `f1` and `f2` rather than at one end of a chain). Both `f1` and `f2` must be pairwise-qualified with `fi`. This places the seed faces at different circumferential positions around `fi`, providing angular diversity for the fit. Multiple seeds may be generated per starting face (one per valid `(f1, f2)` pair); the best trial by area is kept.

*Pairwise qualification*: Two faces are pairwise qualified if they are neighbors, both unassigned for this surface type (`conical_hypothesis < 0`), have sufficiently different normals (cross product magnitude exceeds `MIN_CROSS_THRESHOLD`, e.g. 0.01), and the dihedral angle does not exceed `--angular-tolerance`. The cross product minimum is a seeding criterion only — during BFS expansion, adjacent faces on a finely tessellated cylinder legitimately have near-parallel normals.

**Surface fit**: Attempts both a cone fit (5-parameter) and a cylinder fit (4-parameter) independently, then selects whichever achieves the lower maximum vertex-to-surface distance. If the cone's fitted half-angle is outside the valid range (2°–85°) or its apex is too far from the mesh, the cone fit is rejected. The fits use: (1) Axis direction from smallest eigenvector of the area-weighted normal covariance matrix $M = \sum w_i n_i n_i^T$. (2) For cylinders: project vertices perpendicular to the axis, fit a 2D circle for axis position and radius. For cones: estimate half-angle and apex from vertex positions along the axis. (3) LM refinement minimizing vertex-to-surface distances, with vertices **recentered** around their centroid for numerical stability (parameters shifted back after fitting), and `patience=500` to ensure convergence well below vertex tolerance. With `--verbose 4` or higher, each seed fit logs both the cone and cylinder fit attempts along with their max vertex/surface errors and rejection reasons. See [SURFACE_FITTING.md](SURFACE_FITTING.md) for detailed fitting math.

**Distance metric**: For cylinders (axis origin `o`, axis direction `a`, radius `r`): distance from point `p` is $| \|p - o - ((p-o)\cdot a)a\| - r |$. For cones (apex, axis, half_angle): distance is the point-to-cone-surface distance along the radial direction.

**BFS-specific checks**:
- **Convexity**: Compute the vector from the axis to the face's centroid. Check that the face's normal agrees with the hypothesis convexity (dot product with radial vector has the expected sign). Reject if inconsistent — this prevents merging the inner and outer surfaces of a thin-walled cylinder.
- **Normal-axis consistency (per face)**: on a cone, every facet normal satisfies |n·axis| = sin(half_angle) exactly (tessellation coarseness rotates facet normals *around* the axis, preserving the dot product). Candidate faces deviating more than 0.05 are rejected — faces on spheres and tori violate this as soon as BFS crosses tube rows, keeping doomed trials small.

**Acceptance tests**:
- **Angular coverage validation**: (1) For each face, compute centroid angular position θ around the axis. (2) Sort θ values. (3) Find the largest gap and compute span = 2π − largest_gap. (4) The second-largest gap must be ≤ max(span / 3, `max_gap`), where `max_gap` is a coarse-tessellation floor. Ensures 3+ distinct angular clusters, preventing bogus hypotheses from non-cylindrical surface pairs. Also requires minimum face count. The `max_gap` floor (= `angular_tol × angular_tol_multiplier`, capped at span/2) is applied **only to cylinders** (half_angle ≈ 0) whose hypothesis **already spans ≥45°**: medium-resolution CAD steps small fillets into a few wide-gap rows (four rows across a 90° tube arc ⇒ 24° gaps > span/3 = 23°), and straight-edge fillets are cylinders. True cones keep the strict span/3 rule — a narrow corner sphere-blend patch fits a 6-face ~50° half-angle cone within tolerance and would otherwise win selection over the correct sphere. The span-≥45° gate keeps the floor from admitting thin cylinder bands of a smooth cone (a slice of a shallow cone is locally almost cylindrical); genuine coarse fillets span a real quarter/half-round (~70°+), spurious bands span ≤~15°.
- **Normal-axis consistency (aggregate)**: reject if the variance of n·axis over member faces exceeds 0.0005 — the faces lie on a sphere or torus, not a cone/cylinder. Without this check, stage 2.2 tiles torus fillets with bogus conical hypotheses, which starves torus seeding (stage 2.4 seeds exclude conical-committed faces).

`seed_vtol_multiplier` is 1.0 for cones/cylinders (vertices lie precisely on the surface, not approximated). `angular_tol_multiplier` is 2.0: tessellated cones with N circumferential divisions have inter-strip dihedral angles of 360°/N, which can exceed the default angular tolerance; vertex-to-surface distance is the primary discriminator. The default `--angular-tolerance` is 25° so that seed-pair qualification (which uses the *raw*, unmultiplied tolerance) accepts the ~24° row steps of medium-resolution fillets.

#### 2.3 Deduce spherical hypotheses

Fit spherical hypotheses to connected sets of faces that lie on a common sphere. Faces with conical/cylindrical hypotheses are not excluded — a face can legitimately belong to both a cylinder and a sphere (e.g., equator faces). Stage 2.6 resolves overlapping hypotheses later. A spherical hypothesis consists of:
- `center: [f64; 3]` — Center of the sphere.
- `radius: f64` — Radius of the sphere (always positive).
- `convex: bool` — Whether face normals point away from the center (convex, like the outside of a ball) or toward it (concave, like a bowl).
- `error_max: f64` — Maximum absolute distance from any vertex to the sphere surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

Runs the generic BFS hypothesis generation with the following parameterization:

**Seed type**: The starting face `fi` plus all faces that share at least one vertex with `fi` (the **1-ring vertex neighborhood**). Only faces that are unassigned for this surface type are included. If fewer than 3 eligible faces remain in the seed, skip this starting face (not enough angular diversity for a reliable sphere fit). Returns a single seed per starting face.

A vertex neighborhood on a sphere naturally provides normals spanning 3D, making it a good seed geometry for the 4-DOF sphere fit.

**Surface fit**: Algebraic least-squares sphere fitting. Linearizes $|v - c|^2 = r^2$ as $2c_x x + 2c_y y + 2c_z z + k = x^2 + y^2 + z^2$ and solves via 4×4 normal equations with Gaussian elimination. Extracts center $(c_x, c_y, c_z)$ and radius $r = \sqrt{k + |c|^2}$. Optionally refined with Gauss-Newton minimizing vertex-to-sphere distances. Maximum radius bounded by `bounding_box_diagonal * 1000` as a numerical guard rail. See [SURFACE_FITTING.md](SURFACE_FITTING.md) for detailed fitting math.

**Distance metric**: $| \|p - \text{center}\| - \text{radius} |$.

**BFS-specific check — convexity**: Verify the face normal agrees with the hypothesis convexity (dot product of normal with vector from center to face centroid has the expected sign).

**Acceptance tests — solid-angle coverage validation**: Computes area-weighted 3×3 covariance matrix of unit direction vectors from sphere center to face centroids. Checks eigenvalue ratio ($\lambda_3 / \lambda_1$) to ensure directions span 3D. Fillet faces along a cylinder produce $\lambda_3 \approx 0$ (1D strip), failing this test. Also requires minimum face count (`MIN_SPHERE_FACES`, at least 4 — a sphere needs at least 4 non-coplanar points for a unique fit).

#### 2.4 Deduce toroidal hypotheses

Fit toroidal hypotheses to connected sets of faces that lie on a common torus. Toroidal surfaces commonly appear as fillets, blends, and pipe elbows. A toroidal hypothesis consists of:
- `center: [f64; 3]` — Center of the torus (on the axis, in the plane of the major circle).
- `axis_direction: [f64; 3]` — Unit direction vector along the torus axis.
- `major_radius: f64` — Distance from center to the tube centerline ($R$).
- `minor_radius: f64` — Radius of the tube cross-section ($r$).
- `convex: bool` — Whether the torus patch is on the outer (convex) or inner (concave) surface of the tube.
- `error_max: f64` — Maximum absolute distance from any vertex to the torus surface.
- `error_abs_sum: f64` — Sum of absolute vertex-to-surface distances.

Runs the generic BFS hypothesis generation with the following parameterization:

**Seed type**: Two seed sources:
- *1-ring seeds* (same as spherical): the starting face `fi` plus the 1-ring vertex neighborhood, for faces not committed to any other hypothesis type.
- *Strip-merge seeds*: for faces committed to a cylindrical hypothesis, a 3-ring edge neighborhood — but only when it spans ≥ 2 distinct conical hypotheses. A coarsely tessellated torus (e.g. a pipe elbow) gets tiled by narrow cylinder strips in stage 2.2, each spanning too little of the major arc for a 1-ring torus fit; a 3-ring neighborhood spans several strips and provides the major-arc curvature. A plain cylinder region is a single hypothesis and is skipped cheaply.

Faces already assigned to spherical hypotheses are excluded from seeds. Seeds with fewer than `MIN_TORUS_FACES` (7) faces or fewer than 9 unique vertices are skipped (torus has 7 DOF, needs sufficient data for a reliable fit). Returns a single seed per starting face.

**Surface fit**: Tries **two methods** and keeps whichever achieves the lowest maximum vertex-to-torus distance: (1) **Medial axis / tube center method**: estimates minor radius from normal-line intersections of adjacent face pairs, computes tube centers $k_i = p_i \pm r \cdot n_i$, fits a 3D circle to the tube centers for axis direction and major radius, then LM refinement of 7 parameters $[C_x, C_y, C_z, \alpha, \beta, R, r]$. (2) **SOR axis method**: estimates axis from PCA of vertex positions (smallest eigenvector), center via Kasa circle fit, radii from (h,r) profile circle fit, then same LM refinement. Rejects if minor_radius ≤ 0, major_radius ≤ 0, or minor_radius > major_radius (self-intersecting). See [SURFACE_FITTING.md](SURFACE_FITTING.md) for detailed fitting math.

**Distance metric**: $|\sqrt{(\text{radial} - R)^2 + \text{axial}^2} - r|$ where $\text{axial} = (p - C) \cdot \hat{a}$ and $\text{radial} = |(p - C) - \text{axial} \cdot \hat{a}|$.

**BFS-specific check — convexity**: Compute the vector from the nearest tube center point to the face centroid. The tube center is $c = C + R \cdot \hat{r}$ where $\hat{r}$ is the radial unit vector from the torus axis to the point. Check that the face normal's dot product with this vector has the expected sign.

**Acceptance tests**: (1) Minimum `MIN_TORUS_FACES` (7) — torus has 7 DOF. (2) **Angular coverage around the major circle** (same azimuthal algorithm as cylindrical/conical). (3) **Minor circle coverage from face normals**: the tube angle of each face is measured from its *normal* (ψ = atan2(n·â, n·r̂)); the span (largest-gap algorithm) must be ≥ 0.3 rad. Normals — not centroid positions — are the decisive signal: a flat annulus or a cylinder band lies exactly on a one-parameter family of degenerate tori and can span a nonzero *positional* tube arc on a large-tube-radius member of that family, but its normals are all parallel (tube span ≈ 0).

`angular_tol_multiplier` is 2.0: torus quads twist (doubly curved), so stage 1.3 coplanar fusion often fails on fillet patches, leaving raw triangles whose normals differ by up to ~2× the tessellation angular step.

#### 2.5 Select surfaces to use for reconstruction
Assign each mesh face to exactly one "selected surface" for reconstruction. Since stages 2.1–2.4 produce hypotheses that may overlap (a face can belong to one planar, one conical/cylindrical, one spherical, and one toroidal hypothesis simultaneously), surface selection must resolve these overlaps.


**Greedy area-based selection:**

The algorithm selects hypotheses greedily by total mesh face area, largest first. This naturally resolves conflicts: real surface hypotheses cover many faces with large total area, while bogus hypotheses (e.g., a near-degenerate r=139mm sphere from 3 nearly-coplanar faces) cover minimal area and lose their faces to correctly-fitted surfaces.

**Algorithm:**
1. Compute the geometric area of each mesh face.
2. Build a candidate list of all hypotheses (planar with face count > 1, conical/cylindrical, spherical, toroidal). Each candidate tracks its set of member faces and their total area.

3. **Greedy loop**: while unassigned faces remain:
   a. For each candidate, compute the total area of its still-unassigned faces.
   b. Select the candidate with the largest remaining area.
   c. Assign all its still-unassigned faces to it as a selected surface.
   d. Remove those faces from all other candidates' pools.
4. Assign any remaining unassigned faces to their single-face planar hypothesis (fallback for faces on surfaces not yet fitted: ruled, freeform). Stage 2.7 may replace some of these with NURBS surfaces.

5. Build a `Vec<SelectedSurface>`, one per selected hypothesis, containing:
   - `surface_type: SurfaceType` — enum: Planar, Cylindrical, Spherical, Conical, Toroidal (later: BSpline).

   - `hypothesis_index: usize` — index into the appropriate hypothesis vector in Stage2Output.
   - `faces: Vec<usize>` — mesh face indices assigned to this surface.
   - `vertices: Vec<usize>` — mesh vertex indices on this surface (union of face vertices).
6. **Absorb redundant surfaces** (`absorb_redundant_surfaces`): merge a selected surface S into an adjacent, larger-area surface T when all of S's vertices lie within tolerance of T's fitted geometry AND all of S's face centroids lie within `--surface-tolerance` of T. Two situations require this:
   - *Spurious sub-fits*: a narrow band of faces on a cone/sphere can fit an "exact" cylinder (e.g. 3 azimuthal strips of a shallow cone have nearly coplanar normals) and get committed before the real surface's BFS reaches those faces. Tolerance: 10× vertex tolerance (exact overlap).
   - *Tangent-junction bridge facets*: a mesh facet spanning a tangent junction (e.g. cone–sphere on the nosecone) has vertices on both sides and can't be absorbed by either surface at vertex tolerance; it becomes an orphan single-face planar surface creating fake edges. Tolerance: `--surface-tolerance` (chordal sag), only for surfaces with ≤ 2 faces.
   The centroid check prevents absorbing flat end caps triangulated as a fan from a boundary vertex (all cap *vertices* lie on the adjacent cone/cylinder's boundary circle, but face interiors do not).
7. **Prune unselected hypotheses**: remove any hypothesis that was not selected from the hypothesis vectors in `Stage2Output`. Remap `SelectedSurface` indices to match the compacted vectors. After this step, `Stage2Output.conical_hypotheses` (etc.) contain only selected hypotheses.
8. Validate: every face is assigned to exactly one selected surface.
9. Print a summary: count of selected surfaces by type, total faces covered.

**Why this works better than a fixed type-priority rule:** A fixed priority (e.g., spherical > cylindrical > multi-face planar) fails when a bogus hypothesis of a high-priority type (3-face r=139mm sphere) competes with a correct hypothesis of a lower-priority type (35-face r=2mm cylinder). Area-based greedy selection naturally picks the correct hypothesis because it covers vastly more area.

**Nerfed tests:**
- (none currently — `onshape_pipe_elbow_stage26` now passes `--compare`, see note below.)
- `onshape_pipe_elbow_stage26_compare`: Passes `--compare` with 290 cylindrical + 0 toroidal selections. The pipe elbow's torus surface is tessellated into many narrow strips whose normals span < 180°, so stage 2.5 torus fitting doesn't fire (normals look locally cylindrical). This is expected — the torus is decomposed into cylindrical patches that individually pass `--compare`.

- With the `--compare` flag, for each selected surface:
    - Compute face centroids and project them onto the selected surface.
    - Measure distance from each projected centroid to the nearest surface in the reference STEP file.
    - For all hypotheses, report an error if distance exceeds `--surface-tolerance`.

#### 2.6 Surface refitting (boundary refinement — implemented)

After selection and absorption, some border faces may have been absorbed by a hypothesis that fits within tolerance but isn't the best fit; the misassignment shifts the mesh region boundary by a row, which shifts reconstructed edge boundary vertices and corner vertices off the true feature curves.

**Algorithm (`refine_boundary_assignments`, as implemented):**
1. Build the face → selected-surface map; walk every face that has a neighbor in a *different* selected surface (region-boundary faces — the original plan keyed on multi-hypothesis membership, but pruning/absorption destroy that; adjacency is the modern equivalent).
2. Compute the face's max vertex distance to its own surface and to each adjacent selected surface.
3. Move the face when an adjacent surface fits at least **2× better**. The margin prevents oscillation for true junction faces, whose vertices lie on both surfaces equally.
4. Iterate to a fixpoint (≤ 10 rounds; boundaries only move locally), then rebuild the affected hypotheses' vertex lists and drop emptied surfaces. Hypothesis *geometry* is not re-fitted: the fits were validated during BFS; only the region map changes.

**Corner vertex refinement** (stage 3, `refine_corner_vertices`, runs after tangency detection): SHARP corner vertices are refined by cyclic projection onto their adjacent fitted surfaces — the raw mesh vertex at a corner is only as accurate as the tessellation puts it (where a fillet terminates, the mesh may have no vertex exactly at the CAD corner). Corners adjacent to any TANGENT edge are excluded: tangent surface pairs intersect along a curve rather than an isolated point, and cyclic projection slides along the tangent direction (observed drifting 0.57 mm). A safety guard keeps the original vertex if the refinement moves more than 2× surface tolerance.

#### 2.7 NURBS fallback

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

*Reconstruction note: when this stage is implemented, do the stage-3 reconstruction rework described in §3.x (proposed) first, as its foundation. NURBS patches have arbitrary trimmed boundaries that would stress the current projection-plus-geometric-sewing assembly exactly where it is most fragile (see the §3.x analysis and the `cone_tangent_planes` case), and the NURBS test models are also what finally give the more-general shared-edge assembler real coverage.*

---

### Algorithm Reference

Surface fitting math details (plane fitting, cylinder fitting, sphere fitting, cone fitting, torus fitting) and algorithmic techniques (Unified SOR Framework, Gaussian Map Classification, Medial Axis / Tube Center Method, Quadric Fitting) are documented in [SURFACE_FITTING.md](SURFACE_FITTING.md).

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
- `point: [f64; 3]` — 3D position: the mesh vertex position, refined for sharp corners by cyclic projection onto the adjacent fitted surfaces (see stage 2.6 refitting notes).
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

- With the `--compare` flag (runs after stage 3.2 tangency detection):
    - For each `ReconEdge`, sample the mesh boundary vertices along the edge. For each sample vertex, compute the distance to the nearest edge/vertex in the reference STEP file using `BRepExtrema_DistShapeShape`. Report an error if any sample exceeds the edge's tolerance: `--vertex-tolerance` for sharp edges, `--surface-tolerance` for tangent edges (an imprecise tessellation can place the region boundary one row off a tangent junction — still within chordal tolerance of both surfaces and harmless, since tangent junction curves are built from the fitted hypotheses, not from these vertices).
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
  - **Sphere-cylinder**: `Geom_Circle` (great circle arc on the sphere). Center and radius come from the sphere hypothesis. The circle plane normal is `normalize(cross(v0 - center, v1 - center))`. Arc selection samples mesh boundary vertices to determine forward vs. reverse arc on the periodic circle. Open arcs use this path for slight cones too (the circle is built from sphere geometry + endpoint vertices, so it is valid for near-cylindrical cones such as fillets fitted with a small half-angle).
  - **Sphere-cone (closed loop)**: the tangent circle is NOT a great circle (it sits at radius R·cos α, offset along the cone axis). Constructed as a `Geom_Circle` perpendicular to the cone axis, with axial position and radius averaged from the edge's mesh boundary vertices (`compute_tangent_edge_curve_axis_circle_closed`) — the same vertex-averaging approach as plane-torus.
  - **Plane-cone** (TODO): `Geom_Line` along the cone's generator at the tangent azimuth, from the tangent point toward (or away from) the apex.
  - **Cylinder-cone** (TODO): `Geom_Circle` at the axial position where the cone radius equals the cylinder radius (when coaxial), or `IntSS`-based (non-coaxial).
  - **Plane-torus**: `Geom_Circle` perpendicular to the torus axis. Circle center and radius computed from mesh boundary vertices projected onto the torus axis. Handles both closed loops and arcs.
  - **Cylinder-torus**: `Geom_Circle` computed by the vertex-averaging algorithm against *both* candidate axes — the torus axis (coaxial fillet junctions) and the cylinder axis (tube-continuation junctions, e.g. a straight pipe leaving a pipe-elbow torus, where the junction is a tube cross-section) — keeping whichever circle best fits the boundary vertices.
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
  - **Spherical faces with closed-loop edges near poles** (e.g., pill/capsule hemispheres or dome hemispheres cut by a meridional plane): The boundary circles may pass through or near the sphere's UV singularities (poles at V=±π/2 where U is undefined). When any boundary circle passes within 45° of either pole (detected via `ProjectPointOnCurve` + chord-length formula), the sphere surface is recreated with its Z-axis aligned to a reorientation axis. The axis is determined from: (a) the adjacent cylinder's axis direction for tangent edges, or (b) the boundary circle's plane normal for non-tangent edges (computed from 3 sampled points on the curve). This is safe because for all-closed-loop boundaries, the reorientation axis is unique. The face's V bounds are computed by projecting samples of each boundary circle onto the oriented surface: for sphere-cylinder tangencies the junction is a great circle (V=0), but for sphere-cone tangencies it sits at a nonzero latitude, so fixed hemisphere bounds would misplace the boundary. With a single boundary circle, the face spans from that circle's latitude to the pole on the side occupied by the mesh (determined by projecting a mesh face centroid and checking the V-sign); with two, it spans between them.
  - **Spherical faces with pole vertices** (e.g., rounded cube corners where cylinders meet at a sphere pole): When a spherical face's UV bounds include V=±π/2, vertex U values at the pole are undefined (singularity). The `compute_uv_bounds_from_edges` function excludes pole U values (within 0.01 radians of ±π/2) from the circular gap algorithm to avoid corrupted U spans. The face is then constructed via UV-bounds `MakeFace::new_handlegeomsurface_real5`, producing a triangular-topology face (3 edges, 2 meridians + 1 parallel; pole collapsed to a vertex). These faces have a harmless BRepCheck edge failure (tolerance ~0) on the first edge, resolved by sewing.
  - **Partial revolution** (open edges with vertices): Wire-based construction with pre-set pcurves. Before building the wire, `BRep_Builder::update_edge` sets `Geom2d_Line` pcurves on each IntSS edge for the periodic surface, mapping the edge's 3D parameter range to (u(t), v(t)) in surface UV space. This ensures MakeFace selects the correct arc and shares edges with adjacent planar faces for sewing.
  - **Conical faces**: Same approaches as cylinder: full revolution via UV-bounds (with circular gap for U and cone apex V-extension), partial via wire+pcurves. The cone apex vertex is handled as a UV singularity (V extended to 0 for non-truncated cones). All test models pass: simple_cone, cone_cylinder, nosecone, block_with_conical_hole.
  - **Toroidal faces**: UV-bounds construction for complete torus patches. The torus UV space is doubly periodic: U ∈ [0,2π) around the major axis, V ∈ [0,2π) around the tube. Boundary edges are classified from per-edge projected samples: *iso-V* circles (fillet junctions: constant tube angle, spanning U) give a face full in U with V from the circular gap algorithm; *iso-U* circles (tube-continuation junctions, e.g. pipe elbows: constant major angle, spanning the full tube) give a face full in V with U from the circular gap algorithm over the edges' U positions, disambiguated by mesh centroid projections. Multiple sample points per edge curve (3 instead of 1) define the V window; mesh face centroid projections only *disambiguate* which of the two candidate arcs the face spans — they never extend the bounds (a cap triangle whose three vertices all lie on the plane-torus junction circle is legitimately absorbed by the torus hypothesis, but its centroid projects slightly past the tangent circle and would otherwise shift the face boundary off the junction, breaking sewing). Fillet patches typically span a small V arc (e.g., π/2 for a 90° fillet).

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
- **Multi-shell solids (internal voids)**: after per-shell solids are built, `merge_void_solids` finds solids whose bounding box is contained in a strictly larger solid's box; each such void's shells are added to the container solid (`BRep_Builder::Add`), trying both shell orientations and keeping the one whose resulting volume equals container − Σ voids. The hollow_cube test model (tests/occt/, generated by `gen-test-models` since ccad cannot express enclosed voids) validates one solid of 7000 mm³ instead of two nested solids.

- With the `--compare` flag:
    - Count solids in the reference STEP file.
    - Compute volumes of both constructed and reference solids. Match by sorted volume. Report warning if relative volume difference exceeds 1%. All current test models achieve volume agreement within ~1e-7 relative difference.
    - Compute `BRepExtrema_DistShapeShape` between each constructed solid and the reference STEP shape. Report error if distance exceeds `--surface-tolerance`. This is the tightest geometric accuracy check.

#### 3.x (proposed) Rework: shared-edge topological assembly instead of geometric sewing

**Status: not started. Deliberately deferred and bound to the NURBS effort (see Future Enhancements). This section records the decision and its trigger so it is not rediscovered from scratch.**

**The current design.** Stage 3.4 builds each face independently and stage 3.5 reconciles them with `BRepBuilderAPI_Sewing`, a *geometric* stitch that merges edges by proximity within a tolerance. Edges *are* created once as shared `TopoDS_Edge` objects (§3.1), but two things defeat that sharing by assembly time:
- `MakeFace`/`MakeWire`/`ShapeFix` copy edges, so the "shared" edge in face A is no longer the identical object in face B.
- The UV-bounds construction paths (full spheres, complete tori, full-revolution and — as of the cone-cap investigation — apex-cap cones) synthesize their *own* natural-bounds boundary edges that were never shared with a neighbour at all.

So sewing is load-bearing, and the `ReconEdge.pcurves: [_; 2]` and `BRepVertex.uv_coords` fields (§3.0) that the original design earmarked for a both-sides-shared representation were never populated — the geometric sewing pass stands in for them. The face-construction + sewing block is ~2100 lines carrying ~150 pole/seam/apex/reorient/ShapeFix/degenerate special-cases.

**Why this is the root of a whole class of fragility.** Every hard case in stage 3 — sphere pole caps, torus seams, cone apex caps, and the tangent-edge/corner-vertex *tolerance relaxations* (the "genuinely ambiguous junction rows") — is ultimately a symptom of trimming faces independently and matching them by geometric tolerance. `cone_tangent_planes` is the sharpest example: its two apex-cap cone faces are pie-slices-with-apex whose two boundary curves lie on the same iso-V line; the wire builder collapses them, and a pure UV-bounds cap gets the *area* right but only matches its neighbours along part of the boundary, so sewing leaves free edges → open shell. The problem spans segmentation → edge topology → trimming → sewing, which is exactly why patching any single layer only half-works.

**The proposed design.** Build the B-rep topologically instead of geometrically:
1. Build each edge once (shared `TopoDS_Edge` + shared vertices, with its 3D curve) — as today.
2. For **both** adjacent faces of every edge, compute and attach the pcurve of that edge on that face's surface (populate `ReconEdge.pcurves`). Poles/seams are handled *only* here, in pcurve computation, uniformly for every periodic surface.
3. Insert degenerate edges at UV poles (cone apex, sphere pole) and seam edges at periodic wraps, shared correctly, as first-class topology rather than per-surface fix-ups.
4. Assemble each face and the shell directly with `BRep_Builder` from the shared edges — the shell is closed *by construction*. Drop `BRepBuilderAPI_Sewing`; it becomes at most a validation check.
5. This makes construction lower-level and more controlled (avoid `ShapeFix`/`MakeFace` edge-copying that breaks sharing). The cone apex cap, sphere pole cap, and torus seam then all reduce to "attach the correct pcurve to the shared edge, plus a degenerate/seam edge" — one uniform mechanism instead of ~150 special-cases.

**Why defer it.** 52/57 models pass; nothing is *blocked*. The only current beneficiary is one unusual model (`cone_tangent_planes`) that already writes valid-ish output with a volume warning. The full payoff arrives with **NURBS**: freeform patches have arbitrary trimmed boundaries and will stress this machinery hard, building NURBS on top of projection+sewing is building on sand, and NURBS is also what finally brings test coverage to validate a more general assembler (today it would be a more-general machine exercised only on the analytic surfaces the current code already handles). So the rework is the *foundation to lay when NURBS begins*, not a standalone task.

**Trigger to pull it forward earlier.** If a model becomes *unfixable* without tighter topology — i.e. the tangent-edge/corner tolerance relaxations stop merely loosening tolerances and start letting a real model fail — do the rework then rather than adding relaxation #N.

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

## Visualization (`--viz-stages`)

The `--viz-stages` flag enables interactive 3D visualization for debugging each pipeline stage. The user specifies a stage pattern (e.g. `2.2` for conical BFS, or `2.*` for all stage 2 sub-stages). The companion `--viz-stages-level` flag (default 1) controls how often the visualization pauses: level 1 pauses at committed hypotheses and stage 2.6 selections, level 2 also pauses at each accepted seed, and level 3 also pauses at each individual BFS face acceptance.

### Architecture

The visualization uses a two-thread architecture with synchronous channel communication:

1. **Main thread**: Owns the window, drives the three-d render loop, and displays overlays.
2. **Pipeline thread**: Runs the brepper pipeline (stages 1–4) in the background.

### Communication pipeline

```
Pipeline thread                          Main/render thread
     │                                         │
     │  VizOverlay (via SyncSender)             │
     ├────────────────────────────────────────►  │
     │                                         │ display overlay, wait for keypress
     │  VizAction (via SyncSender)              │
     │  ◄────────────────────────────────────────┤
     │                                         │
     ▼                                         ▼
```

- `VizSender::show_and_wait(overlay)` sends a `VizOverlay` to the render thread and blocks.
- The render thread displays face highlights, translucent geometry overlays (cylinders, cones, spheres, tori), status text, and optionally moves the camera.
- The user presses **Space** (NextStep), **Shift+Space** (NextSeed/skip to next seed), or **Q** (Quit).
- The render thread sends back a `VizAction` and the pipeline thread resumes.

### Key types (in `src/viz.rs`)

- `VizSender` / `VizReceiver`: Channel endpoints. Created by `viz_channel()`.
- `VizOverlay`: A display frame containing face highlights, cylinder/cone/sphere/torus overlays, line overlays, shape meshes, status text, and optional camera focus.
- `VizAction`: The user's response — `NextStep`, `NextSeed`, `SkipToSeed(face_index)`, or `Quit`.
- `FaceHighlight`, `CylinderOverlay`, `ConeOverlay`, `SphereOverlay`, `TorusOverlay`, `LineOverlay`, `ShapeMeshOverlay`: Overlay primitives.

### Interactive features

The render window (powered by three-d with egui overlay) provides:
- **Camera controls**: Left-drag orbit, right-drag pan, scroll zoom (toward cursor).
- **Face numbers** (`N` key): Thin black labels at face centroids, culled by face-normal direction (back-facing labels are hidden). The hovered face number is displayed larger with a white halo.
- **Hover info**: Ray-cast face picking under the cursor displays face index, vertex count, area, and normal in the bottom-left corner. Thick bright-yellow edges are drawn around the hovered face.
- **Click-to-skip**: Left-clicking a face sends `SkipToSeed(face_index)`, causing the BFS loop to skip ahead to the seed containing that face. Click is only processed while the window is paused waiting for user input (`waiting_for_input = true`) to avoid a deadlock where the BFS thread would block trying to send the next overlay while the render loop also blocks waiting for input.
- **Help overlay** (`?` key): Shows all keyboard and mouse shortcuts.

### Integration with generic BFS

The `SurfaceType` trait provides:
- `viz_overlays(fit, mesh, face_list)`: Returns geometry overlays (e.g., a translucent cylinder for conical fits, a sphere for spherical fits).
- `viz_fit_description(fit)`: Returns a short text description of the fit parameters for the status bar.

The generic `generate_hypotheses()` function in `src/stage2/bfs.rs` sends viz overlays controlled by `--viz-stages-level`:
- **Level 1** (default): Pause only after a hypothesis is committed (best trial selected and assigned to faces); also pause at stage 2.6 selection overlays. Least disruptive — good for checking overall results.
- **Level 2**: Also pause after each accepted seed (before BFS expansion begins). Useful to inspect initial fit quality.
- **Level 3**: Also pause after each individual BFS face acceptance during region growing. Very granular; use for tracing specific expansion decisions.

BFS debug verbosity (controlled by `--verbose`) uses consistent indentation (no `[prefix]` tags):
- **Level 1**: Stage summaries (hypothesis counts, timing).
- **Level 2**: Per-hypothesis details (seed face list, final fit parameters).
- **Level 3**: Per-committed-hypothesis progress (seed/fit counts, profiling stats).
- **Level 4**: Per-seed-attempt details (seed rejected/accepted with error values, fit descriptions).
- **Level 5**: Per-BFS-face exploration details (each face accepted, rejected, or skipped with error values).

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

See [tests/RESULTS.md](tests/RESULTS.md) for the full test results table.

Generated by: `./scripts/run_test_table.sh [path/to/brepper]`


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

