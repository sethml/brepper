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

## Stage 1: Mesh Input & Preprocessing

### 1.1 Read STL File
Read the input STL file and generate an in-memory representation of the triangle mesh, with fields for future stages to traverse the mesh and fit shapes to sets of mesh faces. Weld vertices by position (with tolerance) to build connectivity.
- **Library**: OCCT `TKSTL::RWStl_Reader` from `RWSTL.hxx`
- **Input**: Binary or ASCII STL file
- **Output**: `ConnectedMesh`, storing:
    - Vector of mesh vertices with double-precision 3d coordinates.
    - Vector of mesh faces with:
        - vertex_count;  // 3 or 4 (or 0 if the face is unused).
        - vertex_indices[4];  // Indices of vertices, ordered by right-hand rule. Must be coplanar.
        - neighbors[4];  // Index of the mesh face across each edge, or -1 if none. Filled in stage 1.2.
        - gp_Dir normal;  // Mesh face normal, computed from vertices in stage 1.2.
        - planar_hypothesis;  // Index of active planar hypothesis, or -1 if none, or -2 if not yet deduced.
        - cylindrical_hypothesis;  // Index of active cylindrical hypothesis, or -1 if none.
        - spherical_hypothesis;  // Index of active spherical hypothesis, or -1 if none.
    - Vector of planar hypotheses, cylindrical hypotheses, and spherical hypotheses - defined later.
    - Statistics described in stage 1.2.

For now, each face can belong to a single hypothesis of each type. Hopefully that's sufficient since hypotheses should nearly exactly match the vertices, but it's possible that in the future we may need to keep set of candidate hypothesis indices per type, then select the best one in stage 2.6.

### 1.2 Mesh Validation
Traverse the faces of the mesh, collecting some basic statistics, validating the geometry, and populating normals and neighbors.
- Collect stats and optionally print: number of mesh faces, number of mesh vertices, number of mesh edges with 0 neighbors, number of mesh edges with >1 neighbors, number of connected shells, number of solids, number of voids within solids.
- Compute and populate mesh face normals.
- Compute mesh face neighbors based on shared mesh edges.

At this stage, validate the mesh: Edges with >1 neighbor indicate non-manifold geometry; degenerate triangles (zero area), flipped normals (inconsistent orientation within a shell), and self-intersections.

The "number of solids" and "voids" computation is non-trivial and may require ray casting or signed volume analysis—consider deferring this to after surface reconstruction.

---

## Stage 2: Surface Fitting

*Comment: The approach described here is "region growing" from seed faces—this works but may produce suboptimal results because the first seed determines the region boundaries. An alternative is RANSAC-style fitting: randomly sample minimal point sets, fit surfaces, count inliers, and keep the best. RANSAC is more robust to the mesh traversal order. The region-growing approach described here also doesn't naturally handle the case where the same surface type appears in disconnected regions (e.g., two separate planar faces with the same orientation). Consider a hybrid: use RANSAC or global clustering first, then refine with region growing.*

**Response: Given vertices that lie within epsilon of precisely on the surfaces involved, hopefully the first seed will grow to the entire region. I think the depth-first search or breadth-first search should effectively find planar surfaces that are partially disconnected. If they're fully disconnected, we'll consider them independent planar surfaces, which is fine.**

### 2.1 Deduce planar hypotheses
Fit planar hypotheses to all sets of more than one face which are coplanar. A planar hypothesis consists of:
- gp_Dir normal;  // Vector normal to plane. Points toward the outside of the shell/solid.
- double distance;  // Distance from origin to plane in direction of normal.
- faces;  // Set of mesh face indices which fit this hypothesis.
- vertices;  // Set of mesh face vertex indices which fit this hypothesis, in right-hand-rule order.
- error_max, error_min;  // Maximum (positive) and minimum (negative) distance from a vertex to the plane.
- error_abs_sum;  // Sum of absolute value of distance from vertex to plane.
The algorithm:
- Initialize all planar_hypothesis to -2.
- For every face index fi:
    - If face[fi].planar_hypothesis != -2: continue
    - Push a new planar_hypothesis index hi, initialized with this face.
    - Call explore_neighbors(fi, hi). TODO: this is a depth-first search - should we do a breadth-first search instead?
    - Re-fit the planar hypothesis and update error metrics.
    - TODO: If there is only one face in the planar_hypothesis, should we delete it?
- Function explore_neighbors(fi, hi):
    - Set face[fi].planar_hypothesis to hi.
    - Add fi to planar_hypothesis[hi].faces.
    - Add vertices of this face to planar_hypothesis[hi].vertices.
    - For each neighbor ni of face[fi]:
        - If face[ni].planar_hypothesis != -2, continue.
        - If face[ni] is not sufficiently coplanar with planar_hypothesis[hi], continue. TODO: define parameters for coplanarity, probably an acceptable angular deviation and vertex distance. At this step, accept a vertex distance greater than the parameter, since we can attempt to re-fit if it's out of range.
        - If an vertex of face[ni] has a distance greater from the hypothesis than the acceptable vertex distance:
            - Test re-fitting the planar hypothesis to planar_hypothesis[hi].vertices plus this face's vertices that are not already in planar_hypothesis[hi].vertices.
            - If after re-fitting there is still a vertex with greater error than the acceptable vertex distance, continue.
            - Otherwise, assign the re-fit plane to planar_hypothesis[hi].
        - Call explore_neighbors(ni, hi).

*Comment on 2.1: The DFS vs BFS question is worth considering: DFS can get "trapped" in a narrow corridor and accumulate drift before exploring the main region. BFS explores more uniformly. However, the re-fitting step partially mitigates this. A bigger issue: once you re-fit the plane, previously accepted faces might no longer fit the new plane! Consider a final validation pass that removes faces whose vertices exceed the tolerance after the final fit. Also: the "vertices in right-hand-rule order" for the hypothesis is unclear—planar regions aren't simply connected in general (they can have holes), so you'll need a more complex boundary representation.*

### 2.2 Deduce cylindrical hypotheses
TODO. Optional: Worthwhile to generalize to conic/cylindrical? Be sure to handle surfaces with negative curavature correctly - negative radius?

*Comment: Cylinders are parameterized by axis (point + direction) and radius—6 DOF total. Minimum 5 points needed for a unique fit, but robust fitting requires more. Key considerations: (1) A cylindrical patch has principal curvature in one direction only—use this to distinguish from spheres/cones. (2) "Negative radius" isn't the right framing; instead, track whether the surface normal points toward or away from the axis (convex vs concave). (3) For cones: 7 DOF (axis point, direction, half-angle). Cones degenerate to cylinders when half-angle→0, so you might fit cones first and detect near-zero angles. (4) Watch out for nearly-planar cylindrical patches (large radius)—they may fit planes better.*

### 2.3 Deduce spherical hypotheses
TODO. Be sure to handle surfaces with negative curavature correctly - negative radius?

*Comment: Spheres are 4 DOF (center + radius). Minimum 4 non-coplanar points for a unique fit. For "negative curvature" (concave spherical patch), track normal orientation relative to center rather than using negative radius. Key challenge: partial spherical patches are hard to distinguish from cylinders or even planes if the patch is small relative to the radius. Consider requiring a minimum angular extent or using curvature analysis to disambiguate. Also consider toroidal surfaces (donuts, fillets)—they're common in CAD and combine characteristics of cylinders and spheres.*

### 2.4 Deduce ruled surface hypotheses
TODO Optional. Find mesh which is coplanar on one axis, and model as an extruded curve surface/ruled surface.

*Comment: This is a good idea for capturing linear extrusions and sweeps. A ruled surface is defined by two boundary curves with linear interpolation between them. For extrusions, one "curve" is a point (the surface degenerates to a generalized cylinder). Detection: look for parallel mesh edges that share the same direction. Fitting: project to a plane perpendicular to the ruling direction and fit a 2D curve. Watch out for twisted ruled surfaces (rulings aren't parallel)—these are harder to detect and fit.*

### 2.5 Deduce NURBS hypotheses
TODO: for groups of adjacent faces which are covered by one- or two-face planar hypotheses and not cylindrical or spherical hypotheses, try to fit a NURBS or b-spline surface to the vertices.

*Comment: NURBS fitting is complex and requires careful consideration: (1) Parameterization: you need to assign (u,v) parameters to each mesh vertex before fitting. Common approaches: conformal mapping, Floater's mean value coordinates, or discrete harmonic mapping. (2) Degree and knot selection: start with bicubic (degree 3×3); knot placement can use chord-length parameterization or be optimized. (3) Regularization: without it, the surface may oscillate. Consider smoothness penalties. (4) An alternative worth considering: use OCCT's `GeomAPI_PointsToBSplineSurface` which handles much of this automatically. (5) For "freeform" regions that are nearly planar, a plane with small tolerance may be preferable to a NURBS that overfits noise.*

### 2.6 Select surfaces to use for reconstruction
- Iterate until out of valid hypotheses:
    - Select the hypothesis that fits some metric TODO of fitting the most area precisely. Add it to a list of selected surfaces.
    - Mark all faces using that hypothesis used.
    - Delete those faces from all other hypotheses that use them. Delete or mark invalid any hypothesis that ends up with insufficient faces left.
- Every face should be covered by one selected hypothesis.

*Comment: This greedy selection has a potential failure mode: selecting a large-but-poor-fit surface early can fragment remaining regions into pieces too small to fit well. Consider: (1) A quality metric that balances area coverage AND fit quality (e.g., area × (1 - normalized_error)). (2) Penalizing hypotheses that would leave "orphan" faces (faces with no valid remaining hypothesis). (3) Preferring analytic surfaces (plane, cylinder, sphere) over NURBS when fits are comparable, since analytic surfaces are more robust for downstream operations. (4) A backtracking mechanism if selection leads to uncoverable faces. Also: what if a face has no hypothesis at all after this process? This needs explicit handling—possibly flag as error or create a single-face planar patch.*

---

## Stage 3: Surface Reconstruction

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

### 3.1 Create OCCT surface objects
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

### 3.2 Detect and create tangency relationships
- Detect edges between faces where there is numerically a very close to tangent relationship. Mark those edges as tangent.
- TODO: should we modify the surfaces to be numerically tangent? This may be challenging - if a surface is numerically close to tangent to two or more adjacent faces, we may need to do some sort of global optimization or iterate to a fixpoint. I suppose we can try to achieve numerical tangency, and if that's not possible, at least ensure that there's intersection along the extend of the shared edge.

*Comment: Tangent detection is critical for fillets and blends. Suggested approach: compute surface normals at several sample points along the shared mesh boundary; if normals agree within tolerance (e.g., < 0.1°), mark as tangent. For enforcing tangency: modifying analytic surfaces is usually wrong (it changes the geometry), but for NURBS you can add tangency constraints to the fit. A more robust approach: accept near-tangency and use a larger intersection tolerance when computing the shared edge. Also consider G2 (curvature) continuity for high-quality fillets—this matters for rendering/machining but may be overkill for your use case.*

### 3.3 Create OCCT edge wires
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

### 3.4 Create OCCT faces
- For each face descriptor:
    - Populate OCCT face object via surface bounded by wires extracted from adjacent edge wire descriptors.

*Comment: This step uses `BRepBuilderAPI_MakeFace`. Key considerations: (1) The outer wire must be oriented counter-clockwise when viewed from outside the solid (along the face normal). (2) Inner wires (holes) must be clockwise. (3) Wires must be closed and edges must connect end-to-end within tolerance. (4) If face construction fails, ShapeFix_Face can often repair minor issues. (5) For periodic surfaces (cylinders, spheres), ensure pcurves handle the seam correctly—you may need to add a seam edge explicitly. This is often the most debugging-intensive step.*

### 3.5 Construct Shells
- Find sets of connected face descriptors via DFS over the face graph (expore from each face to adjacent faces). For each set:
    - Stitch together an OCCT shell.

*Comment: Use `BRepBuilderAPI_Sewing` for this. Key settings: (1) Set sewing tolerance based on your vertex tolerance from earlier. (2) Enable "SameParameterMode" to ensure edge geometry is consistent. (3) After sewing, check `SewedShape()` for the result. Sewing can merge edges that are geometrically close—this is usually desired but verify the topology matches your intent. If sewing produces a compound instead of a shell, faces weren't connected properly. Also verify shell orientability—`ShapeAnalysis_Shell::CheckOrientedShells` can detect Möbius-strip-like errors.*

### 3.6 Construct Solids
- Convert shells to solid bodies. 
    - TODO: figure out which shells are voids? Does face orientation help here?

*Comment: Yes, face orientation is the key. OCCT convention: face normals point outward from material. For a solid: outer shell normals point out, void shell normals point in (toward the void interior, away from material). To classify: compute signed volume of each shell—positive = outer, negative = inner. Alternatively, pick a point inside the shell and ray-cast to determine if it's inside any other shells. Use `BRepBuilderAPI_MakeSolid` to combine an outer shell with void shells. `ShapeFix_Solid::SolidFromShell` can also create a solid from a single closed shell and orient it correctly. Final validation: `BRepCheck_Analyzer` will verify the solid is valid.*

---
## Stage 4: Output

### 4.1 Output objects
- Write constructed objects to a STEP file (or potentially other formats).

*Comment: Use `STEPControl_Writer` with `STEPControl_AsIs` mode to preserve your exact geometry. Consider also offering `STEPControl_ManifoldSolidBrep` mode which enforces stricter solid validity. Before export, run `ShapeFix_Shape` as a final cleanup pass (but heed your AGENTS.md warning about investigating root causes of any fixes). Set appropriate STEP header metadata (author, organization, etc.) for traceability. For debugging, also consider outputting intermediate formats: BREP (OCCT native), or individual surfaces/curves to help diagnose reconstruction issues.*

---

## Implementation Phases

### Phase 1: Foundation
- [x] Project setup (crates, cargo.toml, dependencies)
- [x] Test utility: read an STL and a STEP, compute maximum distance between STL vertices and STEP surfaces, and print it out. Create a script in scripts/ to apply it to all of the stl/step file pairs under tests/ and print out a table of maximum distances.

### Phase 2: Stage 1 Mesh Input
- [x] Stage 1.1: Read STL file into `ConnectedMesh`, including welded vertices and per-face placeholder fields for neighbors, normals, and hypotheses.
- [ ] Stage 1.2: Mesh validation pass to compute face normals, edge neighbors, manifold stats, connected shells, and orientation consistency checks.

---

## Testing Strategy

1. **Unit Tests**: Individual components (readers, fitters, converters)
2. **Integration Tests**: Full pipeline on known geometries
3. **Regression Tests**: Compare output STEP to reference
4. **Validation**: Round-trip test (STEP → mesh → STEP)

### Test Cases

| Test Case | Input | Expected Output | Status |
|-----------|-------|-----------------|--------|
| Cube | 12 triangles | 6 planes, 1 solid | ✓ Passing |
| Cylinder | Tessellated cylinder | 1 cylinder + 2 planes | ✓ Passing |
| Sphere | Tessellated sphere | 1 sphere | ✓ Passing |
| Cone | Tessellated cone | 1 cone + 1 plane | ✓ Passing |
| Stepped Block | Complex planar | Multiple planes | ✓ Passing |
| L Bracket | Complex planar | Multiple planes | ✓ Passing |
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
