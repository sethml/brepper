# STL to STEP Conversion Utility - Development Plan

## Project Overview

**brepper** (B-Rep from Mesh) - A utility to convert triangulated STL meshes from CAD exports into parametric STEP files with fitted analytic and freeform surfaces.

## Design Principles

The stage 2 architecture was rewritten twice: from the fit-driven BFS
hypothesis pipeline to partition-then-classify (2026-07-17), and then —
after the fit-free curvature-jump classifier was measured to be
unsalvageable (2026-07-20, see the Stage-2 Retrospective's Resolution) —
to the current **sharp-partition + scoped fit competition** design. These
principles drive it; retired mechanisms and the reasons they died are in
[Abandoned Approaches](#abandoned-approaches).

1. **Exactness is the superpower.** Mesh vertices lie on the true CAD
   surfaces within numerical epsilon. A candidate surface either validates
   every claimed vertex within `--vertex-tolerance` or the claim is wrong —
   validation is decisive, not statistical. We are not solving the
   fit-primitives-in-noisy-scans problem the literature (RANSAC lineage)
   addresses.
2. **The answer's boundaries are already in the mesh.** CAD tessellators emit
   seam edges along B-rep face boundaries — *including tangent junctions*
   (verified 2026-07-17 on the nosecone cone–sphere tangency: 0 of 1572 STL
   triangles straddle the tangency ring; 104 vertices lie exactly on it). So
   the true segmentation is a partition of the mesh face graph whose region
   boundaries are a subset of mesh edges, and every region is connected.
   This is architectural rather than lucky, and it is worth knowing why,
   because the whole pipeline rests on it — see **Mesh provenance and the seam
   assumption** below.
3. **Partition by what is certain; decompose the rest by fit competition.**
   Sharp dihedral edges are exact and cheap — they partition the mesh into
   smooth regions with no calibration (measured: every sharp-edged corpus
   model segments perfectly by them alone). Tangent junctions inside a
   smooth region are NOT reliably detectable fit-free (measured: no
   threshold on local curvature differences separates real junctions from
   real curvature gradients — the classes overlap in that statistic), so
   they are found by fit competition: surfaces fitted from seeds compete on
   validated area, and the exact vertex test places the boundary.
4. **Curvature allocates effort; fitting decides geometry.** The per-face
   jet curvature field orders seeds, prunes and orders fit attempts, seeds
   initial parameters, and gates expensive refits — and never decides
   membership or boundaries. Every curvature error converts to wasted or
   saved time, never to wrong output.
5. **Repair by validation, locally.** Remaining mistakes surface as
   validation failures and are repaired by local fit-validated operations —
   resolve multi-claimed junction rows to the better-fitting side, merge
   fragments a neighbor's surface explains, relabel ±1-row boundary error.
   No global fix-up passes with chordal slop.
6. **Invariants are enforced, not assumed.** Before stage 3: every mesh face
   belongs to exactly one region; every region is mesh-connected; every
   region's surface validates its own faces; results are deterministic.
   (The abandoned pipeline's greedy selection violated connectivity —
   3DBenchy face 5841 reached stage 3.4 claiming three disconnected mesh
   regions.)

## Assumptions

- Assumes the vertices of the mesh lie precisely on the CAD surfaces, modulo numerical error/stability.
- Assumes transitions between shapes always have a seam — verified to hold
  even for tangent transitions (see Design Principles #2 and **Mesh
  provenance** below). This holds for every B-rep CAD tessellator; it does
  *not* carry information for mesh-native sources such as OpenSCAD, which are
  a known, deliberately unhandled gap.
- Assumes manifold mesh (each edge shared by exactly 2 faces).
- Assumes consistent face orientation (all normals point outward).

Fails noisily if any of the assumptions are violated, after stage 1.2's
repair pass has fixed what can be fixed without moving a vertex
(`--strict-mesh` disables repair).

### Mesh provenance and the seam assumption

Design Principle #2 — that a B-rep face boundary is always a chain of mesh
edges, tangent or not — is load-bearing for the entire pipeline, so this
records why it holds and where it stops holding. (Researched 2026-07-28.)

**Why it holds: every B-rep faceter is edge-first.** They all discretize each
edge *once*, shared by both adjacent faces, and then fill each face's
parametric domain out to those fixed boundary polylines. This is not a
stylistic choice a vendor might change — it is how watertightness is
guaranteed. OCCT states it outright: "At the first step all edges from a face
are discretized… At the second step, the faces are tessellated", with "edges
tessellation data … used to form contours of target faces and passed as input
to the triangulation algorithm" and "each face is processed separately". ACIS
is documented as "maintaining the edge consistency between adjacent faces";
Parasolid's vertex-matching machinery exists so that neighbouring face meshes
meet exactly. A tangent edge is an edge in the *topology*, so it is discretized
like any other and no triangle can straddle it.

**What the likely sources are.** For models destined for STL and print-sharing
sites:

| Tool | Kernel | Faceter |
|------|--------|---------|
| Fusion 360, Inventor | ShapeManager (Autodesk's ACIS fork) | quad-tree gridding |
| SolidWorks, Onshape, Solid Edge, Shapr3D, Plasticity | Parasolid | per-face, vendor |
| FreeCAD | Open CASCADE | **Delaunay** (Watson) in UV space |
| Rhino | own | regular quad grid + recursive subdivision |
| Tinkercad, SketchUp, Blender, OpenSCAD | none — mesh/CSG native | n/a |

Note the corpus is Onshape-heavy (Parasolid) plus CodeCAD, so it is *not*
representative of the faceter mix in the wild.

**Delaunay does occur, and the assumption survives it.** FreeCAD — the
fastest-growing free option, and the same mesher brepper itself uses — meshes
with OCCT's `BRepMesh`, which is genuinely Delaunay. But it is a *constrained*
Delaunay whose constraints are the pre-discretized face boundaries, so seams
are preserved exactly as elsewhere. What differs between faceters is the
*interior* pattern: grid-based faceters (ACIS, Rhino, Parasolid) leave a
strongly structured, anisotropic mesh, Delaunay much less so. The practical
consequence is that **boundary preservation is dependable across all of them,
while sampling density and anisotropy are corroborating evidence only** — worth
remembering if the tessellation layout is ever used as a signal (it is not
today; stage 2 uses only dihedral angle and the curvature jet).

**OpenSCAD and other mesh-native sources: a known gap, not handled.** These
have no B-rep at any stage — `cylinder()` is a prism and `sphere()` a
polyhedron from the moment they are created, with CSG run on polyhedra (CGAL,
or the newer Manifold backend). Consequences:

- *The seam assumption holds trivially and uselessly.* Every facet is planar,
  so no triangle can span an orientation change — a stronger guarantee than any
  B-rep faceter gives. But **every** edge is a facet boundary, so nothing
  distinguishes an intended face boundary from a discretization step. The seam
  carries no information; only the curvature jump does (a constant dihedral
  through a fillet, dropping to 0° on the flat it meets).
- *There are no tangencies.* What was meant as a C¹ junction becomes a finite
  dihedral of one discretization step. OpenSCAD's fragment rule is
  `n = ceil(max(min(360/$fa, r*2*PI/$fs), 5))` with defaults `$fa = 12`,
  `$fs = 2`, so `n` caps at 30 and *every* curve carries a dihedral of at least
  12°.
- *Our `--angular-tolerance` of 25° is therefore the arbiter*, with a sharp
  crossover: `$fn ≥ 15` — or, on defaults, radius ≳ 4.5 mm, where
  `ceil(πr) ≥ 15` — reads as smooth and gets an analytic fit, while `$fn ≤ 14`
  reads as sharp and stays a prism. A `$fn = 6` hex nut correctly stays
  hexagonal; the near-ubiquitous `$fn = 100` reads as smooth.
- *Where a fit is attempted it will succeed easily.* OpenSCAD **inscribes** its
  polygons — the reason holes print undersized — so a prism's vertices lie
  exactly on the ideal cylinder of radius r, and vertex validation passes at
  ~zero residual. The binding check is the centroid one against
  `--surface-tolerance`: the sagitta is `r(1 − cos(180/n))`, which is 0.11 mm
  for r = 20, n = 30 (passes 0.4) but 0.55 mm for r = 100, n = 30 (fails, so
  the surface degrades to per-facet planes). Large-radius curves at default
  resolution are thus the expected failure mode.
- *The ambiguity is genuine, not a defect.* A 15-gon prism and a coarsely
  tessellated cylinder are the same mesh; nothing in it says which the designer
  meant. `--angular-tolerance` decides, which makes it a policy dial here
  rather than a measurement.

We are not addressing this yet. When we do, the likely shape is a flag marking
the input as faceted-by-construction, or auto-detection — see TODO.md. Two
signals look promising for detection: an OpenSCAD revolution has a *constant*
dihedral of exactly 360/n around the whole curve (a B-rep tessellation's steps
vary and its density changes at seams), and its vertices sit exactly on the
circumscribed circle rather than straddling the true surface.

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

**Stage 2: Segmentation & Classification**
- `gp_Pln` — Plane geometry (normal + distance representation)
- `gp_Cylinder` / `gp_Ax3` — Cylindrical surface geometry (axis + radius)
- `gp_Sphere` — Spherical surface geometry (center + radius)
- `gp_Cone` — Conical surface geometry (axis + half-angle)
- `gp_Torus` — Toroidal surface geometry (for fillets)
- `GeomAPI_PointsToBSplineSurface` — Fit B-spline surface to points (freeform rung)

**Stage 3: Surface Reconstruction**

*Geometry (infinite surfaces and curves):*
- `Geom_Plane`, `Geom_CylindricalSurface`, `Geom_SphericalSurface`, `Geom_ConicalSurface`, `Geom_ToroidalSurface` — Analytic surface objects
- `Geom_BSplineSurface` — NURBS/B-spline surfaces
- `Geom_Line`, `Geom_Circle`, `Geom_Ellipse` — Analytic curve types
- `Geom_BSplineCurve` — Freeform curves
- `Geom2d_Curve`, `Geom2d_Line`, `Geom2d_Circle`, `Geom2d_BSplineCurve` — 2D curves for pcurves

*Surface-surface intersection:*
- `GeomAPI_IntSS` — Compute intersection curves between two surfaces

*Point-curve projection:*
- `GeomAPI_ProjectPointOnCurve` — Project point onto curve (for trimming)
- `ShapeAnalysis_Curve` — Curve analysis utilities

*Topology (B-Rep construction):*
- `TopoDS_Vertex`, `TopoDS_Edge`, `TopoDS_Wire`, `TopoDS_Face`, `TopoDS_Shell`, `TopoDS_Solid`, `TopoDS_Compound` — Topological entities
- `BRep_Builder` — Low-level topology construction (wires, faces, shells, pcurve/seam/degenerate-edge attachment — the stage-3.4/3.5 assembler)
- `BRepBuilderAPI_MakeVertex` — Create vertices
- `BRepBuilderAPI_MakeEdge` — Create edges from curves with explicit parameters and vertices
- `BRepBuilderAPI_MakeFace` — Natural-bounds faces for complete closed surfaces only (full sphere, complete torus); bounded faces are BRep_Builder-assembled

*Topology exploration:*
- `TopExp_Explorer` — Iterate over sub-shapes
- `TopExp::MapShapes` — Build maps of topology
- `BRep_Tool` — Extract geometry from topology (IsSeam, SameParameter, etc.)

*Shape healing and fixing (stage 3.6 only — face/wire/shell construction uses none):*
- `ShapeFix_Shape` — General shape repair on assembled solids
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
        - `normal: Option<[f64; 3]>` — Mesh face normal, computed from vertices in stage 1.2. `None` for degenerate faces. (Native Rust arrays avoid OCCT binding complexity and keep the mesh data structure self-contained.)
    - `stats: MeshValidationStats` — Statistics described in stage 1.2.

Note: the mesh carries no surface-assignment state. Stage 2's face→region
assignment lives in `Stage2Output` (`face_region: Vec<usize>`), keeping
`ConnectedMesh` a pure stage-1 product. (The abandoned BFS pipeline stored
per-face hypothesis indices of each surface type on `MeshFace`; those fields
are gone.)

**Vertex welding** uses a spatial hash grid: each vertex is bucketed by integer coordinates `(round(x/tol), round(y/tol), round(z/tol))`. To find matches, adjacent cells (3×3×3 neighborhood) are checked for distance ≤ tolerance. This gives O(1) average lookup. The weld tolerance is derived from `--vertex-tolerance` as `min(1e-9, vertex_tolerance * 0.01)`, so welding is much tighter than the surface-fitting tolerance.

#### 1.2 Mesh Validation
Traverse the faces of the mesh, collecting statistics, validating the geometry, and populating normals and neighbors.
- Collect stats into `MeshValidationStats` and optionally print: `mesh_faces`, `mesh_vertices`, `mesh_edges_open` (edges with only 1 incident face), `mesh_edges_non_manifold` (edges with >2 incident faces), `mesh_edges_inconsistent_orientation`, `mesh_faces_degenerate`, `connected_shells`, `solids`, `voids_within_solids`.
- Compute and populate mesh face normals using Newell's method (handles both triangles and quads with consistent winding).
- Compute mesh face neighbors based on shared mesh edges. For each edge (vertex pair), collect all face uses. For manifold edges (exactly 2 uses), set both faces as neighbors. Also check orientation consistency: properly oriented adjacent faces should traverse the shared edge in opposite directions.
- Count connected shells via BFS over a face adjacency graph. A shell is counted as a "solid" if all its edges are manifold (exactly 2 faces per edge).
- **Repair pass** (bypassed by `--strict-mesh`): drop duplicate-vertex and
  coincident faces, stitch T-junction slivers (splitting the long-edge
  neighbor at the interior vertex), unify winding per shell, and reorient
  inverted closed shells — everything fixable without moving a vertex.
- With the `--compare` flag: check that all vertices are within `--vertex-tolerance` of the STEP shape, and that all face centroids are within `--surface-tolerance`. Uses `BRepExtrema_DistShapeShape` for bounded distance (respects face trimming). Reports the worst vertex and worst centroid with their distances.

Validation checks, in priority order (first failure is reported):
1. **Degenerate faces**: faces with <3 vertices or zero-area normals.
2. **Non-manifold edges**: edges shared by >2 faces.
3. **Inconsistent orientation**: adjacent faces that traverse a shared edge in the same direction (indicating flipped normals within a shell).
4. **Inverted normals**: for each closed shell, compute the signed volume via the divergence theorem (sum of signed tetrahedra formed by each face and the origin: V = (1/6) Σ v0·(v1×v2)). A negative signed volume is legitimate for an internal void shell (normals point out of the material, i.e. into the cavity): such shells are accepted — and counted in `voids_within_solids` — when their bounding box is contained in some closed positive shell's box. Otherwise a negative volume means all face normals point inward (consistently inverted winding); this catches meshes where orientation is internally consistent (check 3 passes) but globally inverted.
5. **Self-intersections**: intentionally deferred to a later stage where richer geometric predicates are available and O(N²) triangle checks can be avoided.

#### 1.3 Coplanar Triangle Fusion (Triangle → Quad Merging)
Merge adjacent coplanar triangle pairs that share a diagonal edge into quads. CAD tessellators (including CodeCAD/OCCT and Onshape) typically produce quad-strip patterns on curved surfaces; each quad arrives as two coplanar triangles sharing a diagonal. Fusion halves the face count on such surfaces (a direct stage-2/3 performance win), gives each fused face crossing directions along both grid axes (which sharpens the per-face curvature descriptor in stage 2.1), and lets a real narrow face that tessellates as a single rectangle arrive as one region-forming quad. (The old rationale — making curved-surface faces "single-face planar candidates" for BFS seeding — is gone with the BFS pipeline.)

**Algorithm (implemented):**

The algorithm operates in three phases: identify merge pairs, apply merges, then compact the face array.

*Phase 1 — Identify merge pairs:*
- Iterate faces in index order. For each unmerged triangle face `fi`, examine each edge.
- For each manifold edge shared with a neighbor `ni` (where `ni > fi` to avoid processing edges twice), both must be unmerged triangles with valid normals.
- Find the 4th vertex `b` from `ni` (the vertex not on the shared edge). Three vertices (`s0`, `s1`, `a`) define face `fi`'s plane.
- **Coplanarity check**: compute the signed distance from `b` to `fi`'s plane (using `fi`'s unit normal and one of its vertices). If `|distance| <= vertex_tolerance`, the vertex is coplanar within the data's numerical accuracy.
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
- On finely-tessellated curved surfaces (spheres, tori), adjacent quads may be nearly coplanar within `vertex_tolerance`. These additional merges are valid — they reduce face count further.
- Sphere poles: triangles near poles may not have a coplanar neighbor (they border 5 or 6 other triangles in a fan pattern, not a quad strip). These correctly remain as triangles.

---

### Stage 2: Segmentation & Classification

**Status (2026-07-20): fit-driven pivot — the SECOND revision of the
stage-2 rewrite.** The fit-free curvature-jump edge classifier was retired
after the 2026-07-19 calibration campaign produced decisive evidence
against it (measurements in the Stage-2 Retrospective's *Resolution*;
mechanism autopsy in Abandoned Approaches). The partition now keeps only
the signal that proved exact — sharp dihedral edges — and tangent-junction
discovery moves to fit competition, with the osculating-jet curvature
field demoted from decision-maker to effort-allocation prior. The
transition from the current code state (commit ce85dd0: jet-based 2.1/2.2
with the jump classifier and its hysteresis fixpoint) is **one-shot** — no
compatibility path, no threshold gating: the jump machinery is deleted,
not disabled.

**Progress (Rewrite Roadmap steps 1–7 landed — the rewrite is complete):** the jump classifier is
deleted and the jet demoted to the `JetField` prior (step 1, `edges.rs`);
the whole-region DOF ladder is implemented (step 2, `classify.rs`, stage
2.3(a)); resolution & repair — relabel, merge, connectivity re-split,
invariants — is implemented (step 3, `repair.rs`, stage 2.4); scoped fit
competition decomposes the tangent chains the ladder cannot explain (step 4,
`compete.rs`, stage 2.3(b)); bent-tube/torus promotion competes tori as
first-class candidates (step 5, 2.3(c)); and the freeform rung fits the
leftover components competition cannot explain analytically (step 6,
`bspline.rs`, stage 2.3(e)). Sharp-edged and single-surface-per-region
models flow end to end through stages 3/4 with `--compare` (cube,
block_with_hole, pipe, two_holes, spherical_pocket, simple_cone, hemisphere,
ball_on_cylinder, …), as do the tangent chains: **rounded_cube →
6 planes + 12 fillet cylinders + 8 corner spheres**, **nosecone →
cone+sphere+base**, **pipe_elbow → 2 cylinders + 1 torus**; and the freeform
models: **elliptical_arch → 1 B-spline + 3 walls**, **bspline_wave → 1
B-spline + 5 walls**, **elliptical_cylinder → 1 periodic B-spline band + 2
caps**. Coarse-mesh torus fillets (stepped_block_toroidal, part_rounded_cube)
whose leftover blob neither competition nor freeform can fit still fall back
to a clean per-facet polyhedron. The old BFS pipeline is fully deleted (step
7): the analytic fitters now implement the trimmed `surface::SurfaceType`
trait, driven only by the ladder and competition. The corpus is green; the
models the new pipeline does not yet decompose as the old BFS did are
disabled with blocker comments (see the step-7a list).

**Pipeline:** sharp-edge classification + jet prior (2.1) → smooth-region
extraction (2.2) → per-region surface discovery: whole-region ladder, then
scoped fit competition with model promotion (2.3) → resolution & repair
(2.4) → emit selected surfaces (2.5).

**The division of labor** (the load-bearing rule, Design Principle #4):
the jet prior allocates effort — seed order, rung order, initial
parameters, refit triggers — and **never decides membership or
boundaries**; vertex-tolerance fit validation decides all geometry. Every
jet error therefore converts to wasted or saved time, never to wrong
output.

**Precision note**: Mesh vertex positions, normals, and all fitted-surface
parameters use double-precision (`f64`) throughout.

#### 2.1 Sharp edges + jet curvature prior

Two products, both deterministic, both cheap (O(N·ring), single pass —
the classification fixpoint is gone; nothing iterates on classification
anymore).

**(a) Sharp-edge classification — the only *decision* this stage makes.**
A manifold edge is BOUNDARY when its dihedral angle exceeds
`--angular-tolerance`; everything else is INTERIOR. The default 25° stays
above the tessellator's maximum dihedral step on smooth surfaces (coarse
fillets step ~24°/row). This rule is measured-exact on the corpus: every
sharp-edged analytic model segments perfectly by it alone (cube 6,
block_with_hole 7, pipe 4, two_holes 8, spherical_pocket 7,
chamfered_cube 26, cone_15x20 2, simple_cone 3 — one pass, no knobs
beyond the user-facing tolerance). The fit-free *curvature-jump* rule
that used to accompany it is retired (evidence in the Retrospective's
Resolution; autopsy in Abandoned Approaches): tangent junctions are found
by stage 2.3's fit competition instead.

**(b) The jet curvature prior — data for 2.3, never a decision.** Each
face carries the differential summary of an **anchored osculating jet**
(Cazals & Pouget SGP 2005 — the standard robust curvature estimator for
irregular meshes): a quadratic height field `z = Q(u,w)` fitted to vertex
*positions* over a PCA tangent frame of the face's own + 1-ring vertex
neighborhood. Three properties make it robust to how CAD models are
actually triangulated (each was forced by a measured corpus failure):

- **It uses vertex positions** — exact f32 mesh data — not normals. A long
  skinny strip/sliver triangle's plane is a chord whose normal can sit
  several degrees off the true tangent, so normal-based estimators inherit
  that error; positions are ground truth. The fit is tessellation-layout
  agnostic (aspect-ratio slivers, staggered rows, pole fans, strip meshes).
- **The tangent frame comes from a PCA of own + 1-ring vertices, not the
  facet normal** — the smallest-variance direction of the neighborhood IS
  the local normal and stays accurate on strips. Own vertices dominate the
  frame (weight `OWN_WEIGHT`).
- **The face's own vertices anchor the fit** (weight `OWN_WEIGHT = 50`).
  Tessellators never straddle a B-rep boundary with a triangle, so a face's
  own corners always lie on its own surface. There is **no outlier
  rejection** (rejection forces arbitrary one-side reads that vary face to
  face — measured); sharp-crease neighbors are excluded from the point set,
  and tangent-junction mixing is allowed — junction rows read intermediate
  values, which is fine because no decision ever keys on them.

*Estimator history, for context:* per-edge κ = θ/d failed on irregular
tessellation (the centroid distance doesn't measure arc length in the
rotation direction — the 65:1 sliver fillet of rounded_cube_10_r2_fine
read ~100× low and fragmented ~50×). The neighbor-normal tensor fixed
strips and slivers but deadlocked on balanced-mix junction rings
(stepped_block's blend zones). The jet subsumes both — as a *measurement*;
as a junction *detector* the whole family is dead (see Resolution).

**Per-face prior contents** (`JetField`, on `Segmentation`):
- `(κ1, κ2)` — signed principal curvatures (jet Hessian eigenvalues;
  positive convex). Few-% accurate on region interiors; contaminated on
  junction rows (by design — the anchor mixes).
- Principal directions — the Hessian's eigenvectors lifted to 3D in the
  PCA frame. The κ≈0 direction is a cylinder/cone ruling; the κ1
  direction wraps a tube.
- The PCA normal, and the **medial point** `centroid + n̂/κ1` (defined when
  `|κ1|` clears the noise floor) — a curvature-center / tube-center
  estimate feeding torus initialization and bent-tube promotion.

**Type signatures the prior predicts** (used ONLY to order and prune work
in 2.3): plane (0,0); cylinder (κ,0), κ locally constant; cone (κ,0), κ
drifting along the ruling; sphere (κ,κ); torus (κ1 constant around the
tube, κ2 varying, sign change past the inner equator); anything else →
freeform-ish.

Output: `EdgeClassification` (sharp classes + counts) and `JetField`,
both on `Segmentation`, consumed by 2.3 and by diagnostics/viz.

#### 2.2 Smooth-region extraction

Connected components of the face graph over INTERIOR (non-sharp) edges
(union-find; deterministic numbering by smallest member face index). A
**smooth region** is a maximal set of faces reachable without crossing a
sharp edge. On sharp-edged models it IS the final segmentation; on models
with tangency it is a *chain of tangentially-joined surfaces* that 2.3
decomposes — e.g. the capsule is ONE smooth region (both caps tangent to
the wall), a rounded cube's entire outer shell is ONE smooth region, the
nosecone is TWO (the cone+sphere tangent chain; the base disk), and the
pipe elbow's tube surface is ONE (torus bend tangent to both straight
pipes).

Output (in `Stage2Output.segmentation`):
- `face_region: Vec<usize>` — exactly one region per face (total coverage
  by construction).
- `regions: Vec<Region>` — member faces + member vertices.

Regions are connected by construction. There is no minimum size: a single
facet is a valid region.

**Diagnostics** (report-only at this stage): with `--compare`, sample
midpoints of sharp BOUNDARY edges and check distance to the nearest
reference STEP *edge* (within `--surface-tolerance`); report region counts
vs reference face counts. Smooth-region counts are NOT expected to match
reference face counts on tangent models — that match is 2.3/2.4's job.

#### 2.3 Per-region surface discovery

Decompose each smooth region into its true surfaces. Two mechanisms, in
order; both decide membership exclusively by fit validation.

**Validation (used everywhere, decisive by Design Principle #1):**
- Every claimed vertex within `--vertex-tolerance` of the surface.
- Every claimed face centroid within `--surface-tolerance` (chordal
  sagitta ≈ h²/8R stays well inside for real tessellations).
- Orientation coherence: every face normal agrees in sign with the surface
  normal at its centroid — this determines the hypothesis `convex` flag and
  replaces the old per-face BFS convexity checks.
- Numeric-validity guards retained from the fitters (fit-domain
  constraints, not selection heuristics): cone half-angle in [2°, 85°],
  torus `0 < r ≤ R`, sphere radius ≤ 1000× bounding-box diagonal.
- **Sub-tolerance features**: when the fitted curvature radius is below
  `--surface-tolerance` (3DBenchy has r=0.2mm cylinders against the 0.4mm
  default), centroid checks are vacuous at that scale; validate by vertex
  tolerance alone and mark the region (stage 3's chain validation shares
  the blind spot — see TODO, face-5841 family).

**(a) Whole-region ladder.** Try to explain the ENTIRE region with one
surface, ascending DOF; the first fit that validates wins. Simplest-first
resolves degeneracy (any plane also lies on an infinite-radius cylinder; a
single facet is exactly a plane) without coverage heuristics.

| Rung | Surface | DOF | Fit (see SURFACE_FITTING.md) |
|------|---------|-----|------------------------------|
| 1 | Plane | 3 | Area-weighted normal average |
| 2 | Sphere | 4 | Algebraic LSQ + Gauss-Newton |
| 3 | Cylinder | 5 | Normal-covariance axis + 2D circle fit + LM |
| 4 | Cone | 6 | Two-path (pointed/truncated) + LM; the existing code fits cylinder and cone together and this rung keeps that dual-fit internally, preferring the cylinder when both validate |
| 5 | Torus | 7 | Medial-axis tube-center + SOR axis + LM |

The **jet signature orders and prunes the rungs** (a (κ,κ) region tries
sphere first; (0,0) tries plane; a failed prediction falls through to the
remaining rungs — a misprediction costs one cheap failed fit, never a
wrong region). Fits see the whole region — maximal support, which keeps
the strip problem dead for any region that IS one surface. For large
regions, fit on a deterministic vertex subsample, validate against all
vertices, refit augmented with the worst offenders on a near-miss
(bounded retries). Most regions of most models end here.

**(b) Scoped fit competition** — for regions the ladder cannot explain
(tangent chains). Seed-first multi-type competition, run INSIDE the smooth
region; the sharp boundary cordons the arena, so competition cost is
bounded and cross-region accidents are impossible.

- **Seeds**: faces not yet claimed by any committed hypothesis, in
  jet-interior-first order (faces whose (κ1,κ2) best match their
  neighbors' — the descendant of curvature-guided seed ordering;
  misordering costs only time). Claimed faces never seed.
- **Per seed**: run the surface types' fit-and-expand from that same seed
  (jet-pruned rung set, as in (a)) and commit the hypothesis with the
  **largest validated area** (simplest-DOF tiebreak — a plane ties an
  infinite-radius cylinder). Initial parameters come from the jet:
  r ≈ 1/κ1, sphere center ≈ centroid − n̂/κ1, ruling/axis from the κ≈0
  principal direction, torus tube centers directly from medial points —
  which also repairs a known weakness of the medial-axis torus fitter on
  strip meshes (its pairwise normal-line intersections must reject
  near-parallel normal pairs, |n₁·n₂| > 0.95 — exactly what strips are
  made of; jet κ1 needs no pairing). Expect the multi-axis LM rescue
  ladders to become mostly dead — delete them only once the corpus proves
  it.
- **Expansion**: BFS from the seed; accept a face when the current fit
  validates it (vertex tolerance + orientation coherence). Jet-informed
  candidate ordering grows the clean core first, so refits always run
  from well-conditioned support and junction rows arrive last, when the
  model is already determined.
- **Claims**: expansion may grow through faces claimed by earlier
  hypotheses — junction-row faces whose vertices genuinely lie on two
  surfaces end up claimed by both, which is the true ambiguity; 2.4
  resolves multi-claims by best fit.
- **Soundness** (the sub-region lemma, under exactness): an accidental fit
  is by definition supported on a *sub-region*; the true surface fitted
  from an interior seed expands — with vertex tolerance as the arbiter —
  to its entire true extent. When all types compete from the SAME seed,
  largest-validated-area selects the true surface. Cross-surface leakage
  is bounded to genuine junction rows: two distinct analytic surfaces
  agreeing to vertex tolerance over an extended patch would be
  near-identical. (The abandoned BFS pipeline compared hypotheses born
  from *different* seeds at *different* times — that is why commit order
  mattered there and does not matter here.)

**(c) Growth & model promotion — the degenerate-limit growth model.**
Strip-tessellated curved sweeps (the pipe-elbow bend; extrusions) defeat
seed-local fits: a cylinder fit claims one strip and stalls, because the
neighboring strip sits ~R·Δφ²/2 off the fit — hundreds to thousands of
times vertex tolerance. The cure is growing with a model family whose
degenerate cases are *regular interior points* of the parameter space:

- Parameterize the bent tube by a point on the spine, spine direction,
  tube radius `r`, spine **curvature** `k = 1/R`, and a bending-plane
  azimuth — `k = 0` IS the cylinder, no infinity anywhere (the
  Pratt/Taubin move that makes lines a regular case of circles).
- **Promotion trigger — jet-informed, not blanket**: when growth stalls,
  read the stalled frontier's jet. Tube-continuation signature (κ1 ≈ the
  claimed tube's 1/r, small nonzero κ2) → attempt the promotion refit,
  seeding `k` from a circle/line fit through the **medial points** of the
  claimed support. Incompatible signature ((0,0) plane, (κ,κ) sphere) →
  stop; don't spend the refit. This replaces both the old 10×
  `REFIT_SKIP_MULTIPLIER` window (far too tight for strips) and blanket
  chordal-proximity refitting (too expensive).
- **Regularization**: while support spans one strip, `k` is
  unidentifiable — ridge toward `k = 0`, relaxing as support grows.
  Exactness bounds the over-flexibility risk: flexibility cannot absorb
  faces of a genuinely different surface, because vertex tolerance
  remains the acceptance bar.
- Once enough arc is claimed, hand off to the full torus fitter and
  validate normally. **Model reduction at commit**: collapse `k ≈ 0`
  (spine deviation from straight < vertex tolerance over the claimed
  extent) back to a cylinder.
- Scope guard: a bent-tube-only generalization could lock onto a sphere's
  latitude band and stall (the dual pathology). The sphere rung competing
  from the same seed covers spheres; if band-locking is ever observed
  *winning* a commit, the answer is widening the promotion family (the
  SOR-profile / Dupin-cyclide direction in SURFACE_FITTING.md), never a
  mask.

**(d) Model-conditioned jet bands** (effort gates only). With a hypothesis
in hand, the model PREDICTS each candidate face's (κ1,κ2) — *including
gradients*: a fitted cone predicts (cos α/ρ(h), 0) with ρ varying along
the axis, so the gradient class that broke the fit-free classifier is a
non-issue against a model. Generous bands (≈[⅓×, 3×] of predicted — jet
interiors are few-% accurate, junction rows off by up to ~half) may skip
hopeless refit attempts and order candidates. NEVER to reject a face that
vertex-validates. (This is the old pipeline's one genuinely effective
consistency heuristic — the |n·axis| = sin α normal-axis check — reborn
in curvature space with a model to condition on.)

**Identifiability is not model selection.** Two different questions get
conflated when a fit looks wrong. *Model selection* asks which of several
candidates explains the data best, and is answered by comparing residuals and
coverage. *Identifiability* asks whether the data pins a parameter at all — and
a fit with an unpinned parameter is not a slightly-worse explanation but a free
variable that happened to land somewhere, which no comparison between candidates
can detect. The degenerate tori on `stepped_block_toroidal` are the case in
point: they report R = 20…114 mm on 0.5 mm fillets, and on their own support
they genuinely fit *better* than the cylinder they displace, so every
residual-based or coverage-based rule prefers them. What is wrong with them is
that their support sweeps 3–7° about the axis, which cannot determine a major
radius. `MIN_AXIAL_SWEEP_DEG` states exactly that and nothing else.

Keep the three kinds of guard distinct when adding one, because they have
different justifications and different evidence:

- *Identifiability* — the data does not determine this parameter
  (`MIN_AXIAL_SWEEP_DEG`, `axis_normals_consistent`).
- *Validity* — this surface cannot be represented or exported sanely
  (`MIN/MAX_HALF_ANGLE_DEG` bounding a cone away from a cylinder and a plane,
  `MAX_TORUS_RADIUS_RATIO` keeping a torus off its self-intersecting limit,
  rejecting non-finite fits).
- *Model selection* — which candidate wins. This is the only one that should be
  a comparison, and it is the one that most often gets written as a threshold
  when it should not be.

**Mesh precision is a float32 quantum, not a constant — and it scales with the
model.** The corpus's STLs (Onshape, Fusion, the NURBS set) are *binary*, so
vertices are float32 and carry ~1.19e-7 relative precision. The absolute floor
on vertex accuracy is therefore `max|coordinate| × 1.19e-7`, which is a property
of how far the model sits from the origin as much as of the tessellation:

| model | max coord | quantum | measured stage-1 vertex max dist |
|---|---|---|---|
| rounded_cube fine | 23.1 | 2.76e-6 mm | 1.9e-6 |
| stepped_block_toroidal fine | 164.0 | **1.96e-5 mm** | **1.23e-5** |

Both measurements sit just under their quantum. So the "a mesh vertex lies on
the true CAD surface to ~1e-6 mm" figure this plan leans on repeatedly is not a
constant: it is ~1e-6 for small models and twenty times worse for
stepped_block_toroidal, purely because that model is 164 mm across. (ccad's STLs
are ASCII at 7 significant digits, which works out similar.)

*The consequence for `--vertex-tolerance`.* At its 5e-5 default it means
**18× the quantum on rounded_cube and 2.6× on stepped_block_toroidal** — wildly
different demands dressed as one number. That is a strong candidate explanation
for why stepped_block_toroidal is the model every threshold experiment breaks
on: it runs closest to its own data's precision floor, so small changes flip
which fit wins. It also explains the two surfaces of 55 that sit at 6.9e-5 there
— at 164 mm in float32 they largely cannot do better.

*What separates cleanly in quantum units.* Committed analytic surfaces fit their
own support to well under one quantum, across a 7× range of quanta:
rounded_cube 0.47 median / 0.7 max, pipe_elbow 0.37 / 0.6, hill_with_hole
0.38 / 1.2, stepped_block_toroidal fine 0.31 / 1.0. Against that, the spurious
fits `hill_with_hole` used to commit sat at **80×** its quantum. So "a genuine
fit explains its support to about the precision the data has" is a scale-free
statement that holds across the corpus, where the same claim in millimetres does
not.

*What this does NOT justify.* Making the acceptance tolerance adaptive by
scaling it with the quantum was tried and measured wrong. Loosening
`--vertex-tolerance` on stepped_block_toroidal fine gives 138 / 87 / 148 / 161 /
111 / 76 surfaces at 5e-5 / 1e-4 / 2e-4 / 4e-4 / 8e-4 / 1.5e-3 — chaotic, no
trend, no plateau, exactly the signature of the `GROW_VTOL_FACTOR` sweep. A
scale-adaptive bar would put that model at ~3.5e-4 (preserving rounded_cube's
current 18× ratio), which is *worse* than the default. The quantum is the right
**unit** and the wrong **knob**.

**The flood bar is the vertex tolerance — and the self-calibrating one it
replaced was masking a seeding defect.** From 2026-07-27 the flood scaled its
acceptance bar by each fit's own median residual, clamped between
`--vertex-tolerance` and 30× it. The reasoning was sound as far as it went: a
mesh vertex lies on the true CAD surface to ~1e-6 mm, so a face that belongs
misses its fit by the *fit's* error while a face that does not belong misses it
by a geometric departure, and those are different quantities. `acceptance_tol`,
`residual_scale`, `GROW_VTOL_FACTOR` and `RESIDUAL_SIGMAS` were deleted on
2026-07-28; the bar is now `--vertex-tolerance` itself.

*Why it went.* Three things, none of which is the theory being wrong.

First, **the statistic measured the wrong thing**. `residual_scale` was computed
over the faces the fit was built from and then used to judge faces it was not —
interpolation error spent on extrapolation. A small patch has a *small* residual
(fewer constraints on the same DOF) and *badly determined parameters*, so the
criterion is loosest exactly where it should be tightest.

Second, **the empirical case for the 30× ceiling had gone stale**. It rested on
a 2026-07-21 observation that without it "everything returned NO GROWTH". With
the bar set to exactly the vertex tolerance, 12 of the 14 probe models are
bit-identical — `condition_patch`, the seed-patch fallback, the stale-fit refits
and the area-weighted plane offset had each removed part of the cause.

Third, and decisively, **the response to the constant is non-monotonic**:
stepped_block_toroidal fine gives 138 surfaces at 1×, 90 at 2×, 134 at 5×, 113 at 12× and
the exact 71 at 30×. A bound behaves monotonically with a plateau. A lever does
this. The constant was not buying tolerance, it was choosing which accidental
fit won a contested support.

*What removing it bought.* `bspline_blister` goes from 968 surfaces (959 of them
per-facet) to **7** against a reference of 6 — the spurious analytic patches
competition was carving out of it were what had prevented the freeform fit from
taking the region at all. `hill_with_hole` goes from 47 to **8** against a
reference of 7. Both had been broken and unwatched.

*What it exposed, and what is now blocked on it.* stepped_block_toroidal fine regresses
71 → 138, through **323 seeds that produce no growth** — including *flat* seeds
failing to grow a *plane*, which cannot be a tolerance problem, since a wall's
vertices sit at ~1e-15. The cause is `condition_patch`: it trims a straddling
seed patch until the seed validates but gives up at `MIN_CONDITION_FACES`
*without checking whether it succeeded*, and stepped_block's flat regions are
narrow strips between fillets, so a 32-face BFS patch always straddles and
cannot be trimmed to a clean one. The loose bar hid this by validating the
contaminated fit anyway.

**That floor cannot simply be retuned, and this is measured.** At a conical
conditioning floor of 2, stepped_block_toroidal fine recovers to 75 surfaces with its
cylinders and tori exactly right (37 and 18) — but `filleted_pipe` acquires a
spurious 18-face 5.6° cone that fits at 6e-7. At 6, filleted_pipe is correct and
stepped_block is 136. At 4, both are wrong. Making the floor per-type via a
`SurfaceType::CONDITION_FLOOR` associated constant does not help either, because
the conflict is within a single type (conical). It is the same pathology as the
tolerance it was hiding behind: a lever, not a bound.

**So the next step is not a value.** The seed patch must not straddle in the
first place — `seed_patch` grows a depth-limited BFS neighbourhood and hands
`condition_patch` the job of repairing it afterwards, which is the wrong shape.
Three stepped_block tests carry this as their blocker.

**Snap shared boundaries to the junction, not to a residual comparison.**
Nothing in stage 2 reasoned about where a boundary *should* be — it fell out of
flood plus claim resolution, and every remaining `--compare` failure has been a
boundary sitting a tessellation row off. The seam assumption says the correct
boundary already exists in the mesh as a chain of edges, so this is selection
from a guaranteed-nonempty set rather than a search. `locus.rs` supplies the
criterion and `snap_boundaries` (repair.rs) applies it.

*Use the minimum-distance locus, not the intersection.* At a tangent junction
the surfaces graze rather than cross, and their intersection is catastrophically
ill-conditioned: a plane offset by δ from the torus it should graze meets it in
a circle displaced by √(2rδ), which is how a 4.7e-4 mm plane-fit error became
the recorded 4.4e-2 mm edge error. The locus of minimum distance — the set of
points where the two normal fields agree — has no such amplification, and the
reason is worth remembering: displacing the plane along its own normal leaves
the tangency circle exactly where it was, so the locus does not depend on the
quantity that gets amplified at all. Only a *direction* error moves it, and
linearly. Positions depart quadratically from a tangency while normals turn
linearly; key the boundary off the latter. (This is verified by unit test in
`locus.rs`, not just asserted.)

*The meridian reduction.* Nearly every tangent junction in CAD joins two
surfaces of revolution about a common axis. Reduce the pair to the meridian
half-plane (ρ, h) of that axis, where each surface is a line (plane ⊥ axis,
cylinder, cone) or a circle (sphere, torus) and coaxial tangency is exactly
meridian tangency; the junction is the circle swept by one meridian point q,
and the side test is displacement along the curves' common tangent there. This
subsumes the cases one would otherwise write out separately, and explains why:
at a plane-torus junction ρ is the informative coordinate and h is stationary,
at a cylinder-torus junction it is the other way round, so neither coordinate
works alone and the common tangent is the thing to project onto. Two junctions
are not axisymmetric and are built directly: a plane grazing a cylinder along a
ruling, and the tube junction where a straight fillet (cylinder) continues into
a corner fillet (torus) — there the axes are *perpendicular*, the surfaces share
the tube rather than an axis, and C¹ continuity forces the cylinder's axis to be
tangent to the torus's spine circle, so the junction is the foot of the
perpendicular from the torus centre onto that axis. Note that neither fit says
where its own surface *ends*; the junction has to be read off the pair.

*What makes a move safe.* Three conditions together, none of them a tuned
number: the surface a face would move to must actually fit it (soundness — a
boundary correction must never hand a face to a surface that does not explain
it, and this is what stops a large fused planar quad from being swallowed by a
fillet it happens to touch); the locus must put the whole face on that side,
*including its interior* — a disc triangulated with every vertex on its rim has
all its corners sitting on the junction and on both surfaces at once, so the
corners decide nothing and the centroid must; and the moved set must be
reachable from the current boundary through other moving faces. Orientation —
which side belongs to which surface — cannot be read off the supports globally,
because a locus circle is complete while the junction is only the arc where the
surfaces meet, so a large planar face wanders to the far side of its own
tangency circle while remaining unambiguously the plane's. Read it off the
finely tessellated side, whose extent is comparable to the junction, and only
from faces the other surface cannot explain.

*Loud, not lenient.* A tangent pair whose locus cannot be built or does not
validate is reported and left alone; transversal junctions are out of scope by
construction, being sharp mesh edges that 2.1 already partitions exactly. The
report is what found the rim-triangulated-disc case above. Snap runs *after*
relabel and merge, so it neither fights the approximate criterion nor complains
about slivers merge is about to dissolve; on the corpus nothing is reported.
`boundaries_sit_on_their_junctions` (stage2_competition.rs) pins both counts at
zero, which is a much sharper regression guard than a surface count — a
boundary one row off leaves the counts intact.

*What it retired, and honest accounting.* Relabel's **2× anti-flapping margin**
and **10-round cap** are gone. Both were the price of a criterion that compares
a face against two fits *each contaminated by that face*: it could oscillate, so
the margin held it back and the cap stopped it. In-loop refitting removed the
contamination, and allowing each face to move at most once bounds the
alternation as a statement about convergence rather than a tuned number. Both
were inline literals, so the stage-2 constant count does not fall; `locus.rs`
adds two (`COAXIAL_TOL`, a validity precondition on the meridian reduction, and
`TANGENCY_TOL_RAD`, which is stage 3.2's existing tangency criterion given one
definition shared by both stages), taking the grep count 35 → 37.
**Measured: with the margin and cap gone, snap moves zero faces on every corpus
model — unconstrained relabel now reaches the same boundary on its own, and the
suite passes with snap ablated.** So snap is currently a correctness *guarantee*
rather than a correctness *fix*: the boundaries relabel converges to are the
ones the junction geometry demands, which the `dump_regions locus` diagnostic
confirms corpus-wide (zero wrong-side faces). Note what did *not* happen — snap
is not what made retiring the margin safe. Relabel without the margin passes on
its own; snap supplied the *evidence* that where it converges is where the
junction geometry says it should be, which is a claim nothing could previously
make.

**Kept deliberately, and here is when to revisit it.** An inert pass is exactly
the kind of crust this project keeps having to cut away, so it does not get to
stay on sentiment. Two things should decide it, and neither is true yet. First,
**the corpus is not complete enough to say "inert"**: it is Onshape-heavy
(Parasolid — see *Mesh provenance*), every model in it decomposes well enough
that relabel has an easy job, and the place snap would first earn its keep is a
model where relabel does *not* converge — which today is coarse-mesh territory
that `LEFTOVER_BLOB_MAX` makes unmeasurable anyway. Re-ablate once the corpus
covers other faceters and the coarse path is measurable. Second, **the next
deliberate complexity-reduction pass should put it on the list**: run the
ablation protocol on it like any other candidate, and if it is still moving zero
faces on a corpus that can actually distinguish, delete `snap_boundaries` and
keep `locus.rs` for the diagnostic and the invariant test — the property is
worth enforcing even if nothing needs to act on it. Until one of those, it stays:
it is the exact per-face assignment criterion item C (variational shape
approximation) needs, and deleting it now would only mean writing it again.

**The leftover blob is a boundary problem, not a size problem.** When
competition plus freeform left a connected leftover blob bigger than
`LEFTOVER_BLOB_MAX` (8 faces), stage 2.3 threw away the *entire region's*
decomposition and per-facetted it. The premise was sound — a surface abutting a
jagged per-facet blob reconstructs to an invalid wire, and that fallback took
rounded_cube coarse from 43 BRepCheck failures to 0 — but a local defect
triggered a global catastrophe, which is what made coarse meshes untunable:
several unrelated constants looked load-bearing only because they sat on the
safe side of this cliff.

*Size is the wrong quantity, and this is measured rather than argued.*
`cone_tangent_planes` leaves a blob covering **73%** of its region and
reconstructs with **zero** BRepCheck failures; `bspline_blister` leaves one
covering **43%** and produces sixteen invalid faces. No threshold on the blob —
absolute or relative, before or after dissolution — separates those two, and
both a relative gate and a "does the decomposition explain a minority of the
region" gate were tried and rejected on exactly this evidence.

*What does separate them is where the invalidity comes from.* A face's wire
gains one edge per distinct per-facet neighbour, so the face that goes invalid
is the one whose boundary has been shattered into more pieces than the face has
substance. `dissolve_fragmented_surfaces` (classify.rs) therefore compares each
committed surface's per-facet neighbour count against its own face count and
dissolves only the surfaces that lose — carrying no length and no tuned number,
and leaving the rest of the region alone. It is iterated to a fixpoint, since
dissolving one surface can shatter its neighbour's boundary in turn; in the
limit that reproduces the old whole-region behaviour, but only when every
surface in the region individually earns it.

*Results.* `LEFTOVER_BLOB_MAX` is deleted (stage-2 constant count 37 → 36).
`onshape_cone_tangent_planes_three_cones` is **re-enabled and passing** — its
third cone was never a fitting failure, it was fitted correctly and then thrown
away with the region around it. Corpus decomposition is otherwise identical.

*What this exposed, and did not fix.* `bspline_blister` and `bspline_saddle`
had no coverage past stage 1, and had silently collapsed to an all-per-facet
polyhedron (2137 and 1605 planar surfaces, zero analytic) — a valid solid
carrying no information. They now decompose, and that is pinned by
`freeform_models_decompose_rather_than_collapsing`. They also now build a
handful of invalid faces (4 and 3), which the old discard avoided only by
per-facetting everything. Every one of them is a **single-mesh-face** per-facet
plane or one B-spline — the micro-sliver/degenerate-edge family, not the family
the discard claimed to address — and blister's root cause is upstream regardless:
its 926-face freeform component fails the pole-bound check, so nothing can
explain it. Tracked by an `#[ignore]`d test carrying the whole measurement.

**The jet is a prior, not a gate.** Stage 2.3(b) states this — validation
decides every face — but `seed_patch` was gating its BFS on jet similarity, and
on a small fillet tessellated into a handful of facets the jet is not reliable
enough to gate on: adjacent facets of one exact r = 0.5 arc read κ₁ = 0.023,
0.162 and 2.2. The filter starved the patch below what a 5-DOF fit needs, so
every candidate failed and the arc fell to one plane per facet. It now falls
back to the plain topological neighbourhood whenever the filtered patch is too
small to fit; a straddling patch is what `condition_patch` trims and what flood
validation rejects, so nothing is risked. Wherever the jet is consulted, ask
what happens when it is wrong — the answer must never be "no candidate at all".

**Refit whenever the support changes.** A fit is only a statement about the
faces it was computed from, so any step that changes a surface's face set
invalidates its geometry — and the faces such a step removes are precisely the
junction rows that pulled the fit off true, so a stale fit is not merely
approximate, it carries the error of faces the surface no longer owns. Two
places got this wrong and both are now fixed (2026-07-26): `resolve_claims`
trimmed a hypothesis to the faces it won while keeping the fit from everything
its flood reached, and `relabel_boundaries` iterated without refitting between
rounds. The second is the more interesting one: boundary repair and fitting are
one problem, not two, and comparing a face against two fits that are *both*
contaminated by that same face is what stalls the relabel loop inside its 2×
margin. Alternating assignment and fit is what makes the fixpoint worth
iterating to. (Measured: rounded_cube's r = 2 fillets went from 1.9976–1.9984
at 1.5e-3 mm error to 2.0000 at 1.6e-6; stepped_block_toroidal fine went from
73 surfaces to exactly the reference's 71 and now passes `--compare`.)

Relatedly, a fit's weighting must not inherit the tessellator's sampling
density: `fit_plane` area-weighted its normal but took a plain vertex mean for
the offset, so a large flat face ringed by the many tiny quads of a fillet's
first row was placed at the *midpoint* of the two levels. Area weighting is the
correct discretization of ∫(n·x − d)² dA.

A second degeneracy of the same shape is the *unobserved rotation*
(`MIN_AXIAL_SWEEP_DEG`, compete.rs, 2026-07-26). Each analytic type contains the
next simpler one as a limit — a torus with R → ∞ is a cylinder, a cylinder with
r → ∞ is a plane — and a patch that sweeps only a few degrees about its own axis
never observes the revolution, so that radius is free. A free parameter does not
merely produce a useless fit; it produces a *winning* one, because the richer
model is a super-model of the surface it should have lost to: it covers
everything that surface covers and then floods one row further on the slack,
carrying the largest-validated-area contest. So a candidate whose support does
not sweep far enough about its axis to determine its radius is rejected, and the
simpler surface wins the ground back. (Measured on stepped_block_toroidal:
ten such impostors, sweeping 2.8°–6.7°, against real surfaces at 64.8°–103.4°.)

The normal-axis check itself also survives literally, and as of 2026-07-25 it
gates the conical *fit*, not just candidate faces: where the support is thin
in the axial direction (a one-quad-tall cylinder band between two fillet
tori) the vertex residual barely constrains the axis, and a least-squares
cone tips its axis into the band's own plane and beats the true cylinder on
residual. The facet normals are unambiguous exactly where the positions are
degenerate, so a fit whose own support's median facet violates
|n·axis| = sin α is rejected (`axis_normals_consistent`).

**(e) Freeform leftovers.** Connected components of a smooth region left
unclaimed after competition (nothing analytic explains them) go to the
freeform rung below. Before fitting, **dissolve accidental analytic
islands**: a committed analytic hypothesis that is small (≲16 faces) and
enclosed by the leftover component is provisionally dissolved into it; if
the B-spline fit then validates the whole component *including* the
island's vertices, the island stays freeform, else it is restored. This is
the old NURBS gathering's confetti-absorption guard reborn with a
fit-validated backstop instead of a chordal one. Components below 4 faces,
and anything freeform fitting rejects at every grid density, terminally
degrade to per-facet planes — total coverage can never fail.

**Output**: one hypothesis per committed structure, appended to the typed
vectors in `Stage2Output` (`planar_hypotheses`, `conical_hypotheses` —
cylinders are cones with `half_angle == 0` — `spherical_hypotheses`,
`toroidal_hypotheses`, `bspline_hypotheses`). Hypothesis structs keep
their existing fields: type-specific parameters, `faces`, `vertices`,
`convex`, and error statistics (`error_max`, `centroid_error_max`).

##### Freeform (B-spline) fitting

Entered for the unclaimed leftover components of 2.3(e) — parts of a
smooth region that neither the whole-region ladder nor competition
explains analytically, after the accidental-island dissolution guard. The
fitting machinery carries over intact from the previous architecture
(implemented 2026-07-10 through 2026-07-15, `src/stage2/bspline.rs`); what
is deleted is its old region-*gathering* front end (see Abandoned
Approaches) — the smooth-region partition and competition leftovers now
supply the input. A planar margin of a freeform surface is smoothly
connected to it and lands in the same smooth region automatically, which
is what the old "tangent-plane merge" step existed to fix (competition's
plane rung may claim such a margin; the island-dissolution guard decides,
fit-validated, whether it rejoins the freeform patch).

**Algorithm (`deduce_bspline_surfaces`, per component):**
1. Components planar to fitting exactness (10× vertex tolerance,
   NaN-guarded: the area-weighted normal of a closed region cancels to
   zero) become a multi-face planar hypothesis instead — a plane is
   preferable to a B-spline that overfits noise. Components below 4 faces
   are not fitted (they degrade to per-facet planes).
2. **Parameterize** the region onto the unit square: outer boundary loop pinned by cumulative chord length, 4 corners chosen where the adjacent *outside* region changes (arc-length-quartile fallback when the boundary offers ≠ 4 transitions), interior and hole-loop vertices solved with Floater mean-value weights by Gauss–Seidel (convergence enforced — an unconverged map voids Floater's bijectivity guarantee). Quads are split along the v0–v2 diagonal for weighting and sampling. **Closed bands** (a second boundary loop of extent ≥ half the largest — a small hole floats freely, a comparable second loop means the region wraps, e.g. the wall of a closed prism) are first **seam-cut** into a disk along a shortest vertex path between the two loops (`cut_region_seam`): the path follows FACE perimeter edges only (a path along a quad's split diagonal would tear that quad across the cut), its vertices are duplicated, right-side triangle fans switch to the clones, and the boundary re-chains into one loop whose two cut copies carry a `CUT_OUTSIDE` id — corner detection then lands the 4 corners exactly at the cut endpoints. The clones keep their originals' positions, so the two parameter-extreme grid columns of the around-the-band axis carry identical data (the periodic fit drops the closing one; the clamped fallback fit keeps both and closes exactly through the shared column).
3. **Resample a grid** through the flattened triangulation (uniform-bin point location; numeric misses near corner slivers are clamped to the nearest triangle within one bin — feeding them to the hole membrane extrapolates garbage). Grid stations per axis are the boundary vertices' own pinned positions (stations *between* the vertices of a coarse profile sample chord interiors and pin the fit to the chordal prism — a 1.7% volume deficit on the arch, the mesh's own chordal error), subdivided up to an interior-density target of ~1.5·√(region vertices) (a flat margin tessellates coarsely: hill_with_hole has 8-vertex sides on a 946-vertex region). Station spacing is **conditioned in 3D-arc measure** (`axis_stations`): chord-length-parametrized C2 interpolation amplifies data noise by roughly the square of the adjacent knot-spacing ratio, so a boundary mixing fine facets with one long straight chord (the xyz-cube letter walls: a 21× jump) turned a 0.02 mm-flat ribbon into control points 150 mm off the part while passing through every grid node exactly. Near-coincident stations (< 0.25× the median arc gap — the near-singular end of the same conditioning problem) are merged, and segments are subdivided to a fixpoint where no piece exceeds 2× either neighbor segment's pieces, uniformly in arc rather than UV (the Floater map compresses high-curvature areas). Surviving boundary-vertex stations keep their exact pinned positions; new substations land only inside long segments, where a long boundary chord means the mesh itself is flat and mid-segment samples lie on the true surface. Genuine hole interiors are filled with a Laplacian membrane so the fitter sees smooth data. The grid is flipped if needed so d1u×d1v matches the mesh's outward normals: B-spline faces are always `concave = false`.
4. **Fit**: `GeomAPI_PointsToBSplineSurface::Interpolate` (chord-length, C2 cubic) first — the grid is sampled from an exact-on-surface mesh, and OCCT's least-squares approximator reports `IsDone` while overshooting its `Tol3D` by millimeters near strongly curved patch boundaries. The approximator (tol = surface/8, then surface/2) remains as fallback, and grid density retries at caps 121/181/241. **Closed bands fit periodically**: `Interpolate(…, thePeriodic=true)` with the around-the-band axis as the array's FIRST index (transposed when that axis is v; the outward-orientation flip then reverses the second index and accounts for the transpose's normal negation) and the closing station omitted — OCCT appends the wrap-around copy of the first row itself and calls `SetUPeriodic()`. The periodic path is interpolation-only (the approximator binding carries no periodic flag). The resulting closure is C1, not C2 (position and D1 match across the wrap to ~1e-14; D2 jumps — measured), which nothing downstream requires. The fitted around-axis knots are then rescaled to exactly [0, 2π] via `SetUKnots` — an exact affine reparameterization (~1e-15 drift, measured) — because stage 3.4's wrapped-face machinery hard-codes that period. The U=0 iso is then exactly the cut path (chord-length parameters start at the first supplied row), a free seam gauge.
5. **Validate**: first a **control-point bound check** — every pole must lie within the grid data's bounding box inflated by half its diagonal. An ill-conditioned fit passes through every grid node exactly yet swings wildly between them, and the vertex/centroid validation below samples at data points, making it structurally blind to inter-node excursions; the surface lies in the convex hull of its poles, so poles near the data bound the whole surface near the data (the xyz-cube blow-up put poles 22 diagonals out; legit fits, including periodic bands, stay well inside the margin). Then every region vertex and face centroid within `--surface-tolerance` of the fit, each projection seeded from the vertex's own parameterization value (sequential seeding breaks at mesh row wraps, and neither `NextValueOfUV` nor `ValueOfUV` recovers reliably on many-knot surfaces). On failure at every grid density, the region returns to the 2.4 fallback ladder.
6. **Closed-band result**: the periodic fit produces ONE `BSplineHypothesis` with `u_periodic = true`, reconstructed by stage 3 as one wrapped face with a real seam edge — the form CAD kernels (and the references) write. If the periodic fit fails validation at every grid density, the band retries with the clamped fit and is **seam-split** (`split_cut_hypothesis`) into two `SelectedSurface::BSpline` halves sharing the fitted surface object; their two cut edges are plain edges that a reference model does not have, so the stage 3.1 edge compare will (correctly) flag them — the fallback trades compare-exactness for robustness and is currently exercised by no model.

The `BSplineHypothesis` owns the fitted `Geom_BSplineSurface` (stage 3 reuses the same object via fresh handles) plus a cached, seed-rechecking `ShapeAnalysis_Surface` projector used by the 2.4 repair operations, `--compare`, and stage-3 charts/normals/projections. Stage-3 charts of B-spline faces also carry a parameter-space `home` anchor (the projection of the region's mean face centroid): a projection landing on a coincident-boundary curve of a seam-split band (the fallback) is ambiguous between the two parameter extremes (Newton, the gap re-check, and `ValueOfUV` all cannot tell — both branches are exact), and the branch nearer the face's own interior wins. On a periodic (u_periodic) chart the U-arm of that resolution is replaced by principal-value reduction to [0, 2π) — the same contract as the analytic periodic charts — and the whole wrapped-face machinery (edge-sample branch unwrap, winding classification, seam planning/anchoring, pcurve shifts) applies verbatim. There is no compare exemption: the stage 3.1 edge/vertex compare is fully armed for every edge, including the seam (a real feature both the reconstruction and the reference have).

**Periodic-band junctions (stage 3.3)**: the band's two cap junctions are closed-loop edges (no vertices) and must carry PERIODIC curves — seam planning re-anchors a closed edge with `MakeEdge` over [t, t+period] at a solved seam azimuth, which a closed-but-non-periodic curve refuses. Two mechanisms deliver that (both were needed on elliptical_cylinder): (a) an IntSS curve that comes back closed-but-non-periodic (observed type: `Geom_TrimmedCurve`, exact to ~1e-14) is converted via `GeomConvert::CurveToBSplineCurve` + `SetPeriodic` (exact, measured); (b) an IntSS curve that comes back as a genuine partial arc falls to the boundary-chain repair, whose closed-chain arm interpolates a periodic C2 B-spline through the chain vertices (`GeomAPI_Interpolate` with the periodic flag, closing point not repeated — ~4e-15 chain deviation observed). Partial-arc detection for closed loops keys on the ENDPOINT GAP, not the parameter span: a closed B-spline junction arrives over an arbitrary range (observed [0,1]), and the old span≠2π test replaced the true junction with a bogus circle (14.8 mm off, 12% volume error).

**Periodic closed-band conventions (probed 2026-07-15, `tmp/probe_periodic.rs`; confirmed in `GeomAPI_PointsToBSplineSurface.cxx`)**: periodic along U = the array's FIRST index; supply each row once (OCCT appends the wrap-around copy of row 1 itself; a duplicated closing row would make the system singular); the result already has `SetUPeriodic()` applied. The around-the-band axis is NOT always U in the Floater square (elliptical_cylinder's cut copies land as the u-extreme sides, making v the around axis): `cut_around_axis` detects it from the `CUT_OUTSIDE` loop entries, and `fit_grid` transposes the array so the surface's U is always the around axis — `BSplineHypothesis.u_periodic` is a single flag, no per-axis plumbing. The stage-3.4 wrapped-face machinery (branch unwrap in `sample_edge_uv`, winding classification, `align_loop` offsets, pcurve shifts, seam handling — `grep TWO_PI src/stage3.rs`) hard-codes the 2π period; `SetUKnots` rescales the fitted knots exactly, after which it all applies verbatim once the chart reports `u_periodic()`. Debug probes: `BREPPER_DEBUG_EDGE_UV=1` prints every chart projection during edge sampling; `--stage=3.4 -v` prints per-edge pcurve deviations + UV endpoints per face; the `-vvv` fit-residual-vs-own-grid probe separates fit error from sampler/projection error.

**Implementation notes (from the stage-3 rework, 2026-07-08; all implemented — kept for reference):**
- *The chart contract.* A chart provides `point_uv(p) → (u, v, pole_radius)` and `point_at(u, v)`, plus `u_periodic()`/`v_periodic()`. For a B-spline chart: `pole_radius` is ∞ (like the plane chart — fitted height-field patches have no UV degeneracies), and **both periodic flags must be false for open patches** (`u_periodic()` is true only for the periodic closed-band fit). Every mod-2π unwrap/branch-lift in stage 3.4 is gated on those flags; reducing non-periodic parameter coordinates mod 2π silently corrupts them (a real bug on plane charts once, masked for months — see LEARNINGS). With both flags false, the whole seam/winding/degenerate-pole machinery self-disables and B-spline faces take the plain-loop path end to end.
- *The chart must wrap the exact `Handle(Geom_BSplineSurface)` stored in `FaceDescriptor.surface`* — never a reconstruction. Use one `ShapeAnalysis_Surface` per face: `new_handlegeomsurface` + `value_of_uv(&mut self, pnt, precision)` — note `&mut` (it caches internally). For sampling along an edge curve, prefer `next_value_of_uv(prev_uv, pnt, preci, maxpreci)`. Validate each solution with `gap()` and fall back to a fresh `value_of_uv` when the gap exceeds tolerance.
- *Every `match` on `SelectedSurface` needs a `BSpline` arm* — `create_occt_surface`, `surface_faces`, `surface_normal_at_point`, `project_point_to_selected`, `build_face_chart`, the concave-flag matches, the tangent-edge dispatch in 3.3, and the cosmetic `stype` labels + viz overlays.
- *Fix the fitted patch's normal orientation at fit time.* `GeomAPI_PointsToBSplineSurface` follows the orientation of the input point grid, so build the grid such that d1u×d1v matches the mesh's outward normals — then B-spline faces are always `concave = false`/FORWARD. Choosing the convention at fit time is one line; recovering it later is a per-face heuristic.
- *Edge curves lean on the boundary chain.* `GeomAPI_IntSS` is flakier against B-spline surfaces than between quadrics, and tangent-classified junctions (a fitted patch usually meets its neighbors near-tangentially) may have no analytic construction — both routes end at the mesh-boundary-chain machinery (`repair_edge_curve_from_boundary`, `edge_curve_from_boundary_chain`), which is surface-pair-agnostic and fits line, circle, AND a `GeomAPI_PointsToBSpline` curve through the chain vertices, keeping whichever deviates least.
- *What NURBS will not stress — don't build for it:* seams, winding loops, degenerate pole edges, and edge splitting apply only to periodic/wrapped faces; non-periodic B-spline patches never take those paths.
- *Test models*: `tests/nurbs/` — 8 exact-reference STL/STEP pairs from `gen-test-models` (see tests/nurbs/README.md for the difficulty ladder and the strip-mesh rationale: ruled/extruded models ship full-length strips with NO interior vertices along the ruling, exactly what real CAD exporters emit).

#### 2.4 Resolution & repair

Turns 2.3's claims into a strict, connected, validated partition. All
operations are fit-validated and deterministic (fixed iteration order —
by region/hypothesis index, smallest region first for merges — bounded
rounds).

**Claim resolution** — the boundary placement at tangent junctions: faces
claimed by multiple hypotheses (junction rows) go to the claimant whose
surface fits their vertices best, requiring a **2× margin** to prefer a
different side than the larger claimant (a true junction face fits both
sides equally — the margin prevents arbitrary flapping; residual ties
break to the larger claimant, then lower hypothesis index). Exact vertex
distances, no curvature involved.

**Boundary relabel** (±1-row polish; port of the old
`refine_boundary_assignments`): after resolution, walk region-boundary
faces; move a face when an adjacent region's surface fits its vertices
better; refit every surface a move touched, and iterate to a fixpoint,
allowing each face to move at most once.

**Boundary snap** (`locus.rs`): after relabel and merge, put each *tangent*
junction on the mesh-edge chain the two fitted surfaces' junction geometry
picks out. See "Snap shared boundaries to the junction, not to a residual
comparison" above.

**Merge**: absorb region S into an adjacent larger region T when every
vertex of S lies within `--vertex-tolerance` of T's fitted surface and
every centroid of S within `--surface-tolerance`; refit the union,
revalidate, keep only if valid. Eats junction-circle singletons (faces
whose vertices sit ON a tangent seam — genuinely ambiguous, claimable by
either side), chord-degenerate slivers (a single-facet-wide band of a
coarse cylinder validates as its own facet plane, but its vertices lie
exactly on the neighbor cylinder), and residual same-surface fragments.
This is the abandoned absorb pass's *tight* rule reborn with a principled
criterion — vertex tolerance, at which "on the same surface" is decisive;
never the chordal slop of the abandoned loose rule (0.4mm slop, 794
BRepCheck failures on planar-only 3DBenchy). B-spline surfaces are never
merge TARGETS for faces outside their fitted region: a fit vouches only
for its data (the island-dissolution guard in 2.3(e) is the fit-validated
way into a freeform patch).

**Connectivity re-split + invariants** (hard assertions before stage 3):
resolution and relabeling can disconnect a region — re-split into
mesh-connected components (identical fitted geometry on two non-adjacent
regions is legal). Then assert: every mesh face belongs to exactly one
region (total coverage); every region is mesh-connected; every region's
surface validates the region's own vertices/centroids.

(The old standalone "split repair" is gone: under-segmentation cannot
survive 2.3 — competition decomposes within smooth regions — and a region
nothing explains degrades to per-facet planes via 2.3(e). There is no
mechanism left that re-runs a curvature classifier at lowered thresholds;
that mechanism died with the classifier.)

#### 2.5 Selected surfaces & validation

- Build `Vec<SelectedSurface>` — one per region, the **unchanged stage-3
  input contract**: `surface_type` (Planar, Cylindrical, Spherical, Conical,
  Toroidal, BSpline), `hypothesis_index` into the typed vectors, `faces`,
  `vertices`.
- Prune hypothesis vectors to selected entries and remap indices (as
  before).
- Print a summary: surface counts by type, faces covered, repair-operation
  counts.
- With `--compare`: for each selected surface, project each member face's
  centroid onto the fitted surface and measure the distance to the nearest
  reference STEP surface with `BRepExtrema_DistShapeShape`; error above
  `--surface-tolerance`. (Centroid-to-surface distance includes tessellation
  sagitta ≈ h²/8R — small relative to `--surface-tolerance` for typical
  tessellations.) Also re-run the 2.2 partition diagnostics against the
  final (post-repair) partition, now as hard checks.

#### CLI notes (stage numbering and modes)

- `--stage` sub-stages: **2.1** edge classification, **2.2** region
  extraction, **2.3** classification, **2.4** repair, **2.5** selection.
- `--surface-types=<list>` (new) restricts the ladder and competition rung
  set (default: all). `--surface-types=planar` is **exact planar-only
  mode**: curved smooth regions fail the one-rung ladder, competition's
  plane rung claims only coplanar sets, and the rest degrades to per-facet
  planes via the 2.3(e) terminal fallback — the mode that converts
  3DBenchy with volume matching the official 15.55 cm³ — without any
  chordal-slop absorb by construction.
- `--skip-stage` is retired with the old stage numbering; its
  type-skipping role is `--surface-types`, and its absorb/refinement gates
  have no equivalent (absorb is gone; repair is core, and debugging-level
  disabling of repair operations can be a `-Z`-style debug option if ever
  needed).

#### Implementation hints (read before writing code)

Function-level guidance for whoever implements the pivot, written against
the code state at commit ce85dd0 + this doc change. The transition is
one-shot: delete, don't gate.

**Starting inventory:**
- `src/stage2/edges.rs`: sharp rule + the retired jump classifier
  (hysteresis fixpoint, `FIRE_RATIO`/`UNFIRE_RATIO`) + the jet estimator
  `jet_curvature(origin, z_hint, points) -> Option<(f64, f64)>` with three
  unit tests worth keeping (`jet_recovers_quadratic`,
  `jet_skinny_cylinder_strip`, `jet_anchors_own_surface_at_junction`).
- `src/stage2/regions.rs`: union-find extraction — keep as-is;
  `is_boundary` becomes sharp-only automatically once `EdgeClass::Jump` is
  deleted.
- `src/stage2/mod.rs`: new `stage2()` at the 2.2 frontier + the nerfed
  `stage2_bfs_pipeline` (uncalled) + `Segmentation` on `Stage2Output`.
- `src/bin/dump_regions.rs`: the main diagnostic (modes: summary,
  singletons, adj, probe, geom, leaks, rdesc, jet) — keep it building at
  every step; extend it to dump the `JetField` and per-hypothesis claims.
- `tests/stage2_segmentation.rs`: 14 tests whose region-count expectations
  are pinned to the JUMP classifier's output — **re-pin to smooth-region
  counts by measuring with dump_regions, do not guess** (tangent chains
  merge: capsule 1 smooth region, nosecone 2, pipe-elbow tube 1,
  rounded_cube shell 1; sharp-edged models unchanged: cube 6,
  block_with_hole 7, …).

**edges.rs surgery (roadmap step 1):** delete `EdgeClass::Jump`, the
fixpoint/hysteresis loop, and its constants; keep the sharp rule and the
single-pass jet with sharp-crease exclusion. EXTEND `jet_curvature` to
also return the Hessian's eigenvectors lifted to 3D (they are computed
internally and currently discarded) plus the PCA normal; derive the medial
point `centroid + n̂/κ1` when `|κ1|` clears the noise floor. New output
struct `JetField { k1, k2, dir1, normal, medial: Option<[f64;3]> }` per
face, stored on `Segmentation`.

**Fitters to reuse** (widen visibility as needed — most are `pub(super)`
or private): `fit_plane` (planar.rs), `fit_sphere` (spherical.rs),
`fit_cylinder` (cylindrical.rs), `fit_cone`/`fit_cone_primary`/
`fit_cone_normal_line_rescue` (conical.rs), `fit_torus`/`fit_torus_inner`/
`fit_torus_sor`/`fit_torus_axes_rescue` (toroidal.rs), and the
`vertex_to_{cylinder,sphere,cone,torus}_distance` functions for
validation. With jet initializers the rescue/multi-axis paths should
become dead — delete them only after the corpus proves it. The LM
recentering inside the fitters (vertices recentered about their centroid)
is load-bearing for conditioning; keep it.

**Determinism rules (hard-won — violating them makes fits vary run to
run):** never iterate a `HashSet`/`HashMap` directly — sort first
(floating-point non-associativity changes fits); break all ties by
face/hypothesis index; subsample deterministically (stride, not random).

**Competition mechanics (2.3(b))**: per-face claims as a small vec of
hypothesis ids (multi-claim only occurs on junction rows — keep the
common single-claim case allocation-free); a committed hypothesis records
{type, params, sorted face list, error stats}. Seed loop skips any face
with a claim. Expansion is plain BFS with the vertex-tolerance test;
candidate ordering and refit gating may consult the jet (model-conditioned
bands, §2.3(d)) but the accept/reject decision is the vertex test alone.

**Stage-3 contract**: `SelectedSurface` + the typed hypothesis vectors are
unchanged — stage 3 needs no edits. Before deleting the `MeshFace`
hypothesis index fields in the deletion pass, grep `src/stage3.rs` for
reads of `planar_hypothesis`/`conical_hypothesis`/`spherical_hypothesis`
(`build_face_to_surface_map` should key off `SelectedSurface.faces`;
verify before deleting).

**Deletion pass checklist** (roadmap step 7): `bfs.rs`;
`stage2_bfs_pipeline` + `compute_seed_scores` + `BfsConfig` masks; the
seeding/BFS entry points in the five fitter files; `select_surfaces`,
`absorb_redundant_surfaces`, `prune_unselected_hypotheses` (port
`refine_boundary_assignments` into repair.rs FIRST); `--skip-stage`; the
`MeshFace` hypothesis fields; SURFACE_FITTING.md's seeding/BFS/acceptance
sections get pointers to Abandoned Approaches.

**Cost guardrails**: watch pill_fine and rounded_cube_fine stage-2 wall
time (the old pipeline's pathologies were 429s/479s there when masks were
off; whole-region ladder + one competition per region should land well
under the old masked baseline of 1.2s/4.0s) and count promotion-refit
attempts at `-vvv`.

**Acceptance narratives (what "working" looks like, per model):**
- *nosecone* (smooth regions: cone+sphere chain, base): ladder fails on
  the chain; a sphere-interior seed claims the 644-face sphere, a
  cone-interior seed claims the cone (fitted from large support — its
  curvature gradient is irrelevant to fitting); the junction row is
  multi-claimed and resolved; base is ladder-planar. 4 surfaces.
- *pipe_elbow* (tube surface one smooth region): ladder fails; cylinder
  seeds claim the straight pipes and stall at the bend; jet-triggered
  promotion discovers k=1/R and the bend commits as a torus;
  tube-cross-section junction rows resolve. Torus + 2 cylinders + caps —
  no strip machinery anywhere.
- *rounded_cube* (outer shell one smooth region): competition decomposes
  it into 6 planes + 12 fillet cylinders + 8 corner spheres.
- *tests/nurbs*: leftovers → freeform; exercise the island-dissolution
  guard (accidental cylinder ribbons on the saddle).

**Viz**: 2.1 = sharp `LineOverlay` (red) only — the jump/orange overlay
dies with the classifier; 2.2 = per-region `FaceHighlight` palette; add
2.3 pauses per committed hypothesis reusing the `FaceHighlight` +
type-overlay patterns from `bfs.rs`'s viz block (preserve the per-type
`viz_overlays`/`viz_fit_description` functions as free functions before
deleting `bfs.rs`).

---

### Algorithm Reference

Surface fitting math (plane, cylinder, cone, sphere, torus fitting; LM
refinement; the SOR (h,r) profile framework; medial-axis tube-center method;
quadric fitting) is documented in [SURFACE_FITTING.md](SURFACE_FITTING.md).
Its sections on seed selection, BFS growth parameters, and global acceptance
tests describe the abandoned pipeline (see Abandoned Approaches) — the
fitting math itself carries over unchanged and is what the 2.3 ladder rungs
call.

---

### Stage 3: Surface Reconstruction

Build OCCT B-Rep topology (faces, edges, vertices) from the selected surfaces. This is the most complex stage — surface intersections in 3D are numerically delicate, and small perturbations can cause large changes in intersection curves or cause them to disappear.

The stage-2 rewrite does not change this stage's input contract
(`Vec<SelectedSurface>` + typed hypothesis vectors) or its algorithms. It
strengthens one precondition: selected surfaces are now guaranteed
mesh-connected (post-2.4 invariant), which the old pipeline violated (the
3DBenchy face-5841 3-winding-loop abort).

**Core data structures:**

A vector of `FaceDescriptor` structs, one per selected surface:
- `selected_surface_idx: usize` — Index into the selected surfaces from stage 2.5.
- `surface: OwnedPtr<HandleGeomSurface>` — Infinite OCCT surface, type-erased (created in 3.1).
- `adjacent_faces: Vec<usize>` — Indices of adjacent FaceDescriptors, ordered topologically around this face's boundary. A face may appear multiple times if it shares multiple edges with this face (e.g., a cylindrical face adjacent to itself via a seam edge).
- `edge_indices: Vec<usize>` — ReconEdge indices. `edge_indices[i]` is the edge between this face and `adjacent_faces[i]`.
- `vertex_indices: Vec<usize>` — BRepVertex indices at corners. `vertex_indices[i]` is the vertex between `edge_indices[i]` and `edge_indices[(i+1) % N]`.
- The `TopoDS_Face` objects built from these descriptors live in `Stage3Output.faces` (parallel to `face_descriptors`), populated in 3.4, with `Stage3Output.concave_faces` flagging faces whose natural normal points into the material.

A vector of `ReconEdge` structs:
- `face_indices: [usize; 2]` — Indices of the two adjacent FaceDescriptors.
- `vertex_indices: [usize; 2]` — Indices of the two BRepVertices at each end. For closed-loop edges (no vertices), both are `usize::MAX`.
- `curve_3d: Option<OwnedPtr<HandleGeomCurve>>` — 3D intersection curve, trimmed to vertex endpoints (populated in 3.3).
- `pcurves: [Option<OwnedPtr<HandleGeom2dCurve>>; 2]` — pcurve on each adjacent face's surface, parallel to `face_indices` (populated in 3.4 by the shared-edge assembler).
- `tangent: bool` — Whether the adjacent surfaces are tangent along this edge (detected in 3.2).
- `mesh_boundary_vertices: Vec<usize>` — Mesh vertex indices along this boundary, ordered along the boundary.

A vector of `BRepVertex` structs:
- `point: [f64; 3]` — 3D position: the mesh vertex position, refined for sharp corners by cyclic projection onto the adjacent fitted surfaces (see corner vertex refinement below).
- `adjacent_faces: Vec<usize>` — Indices of FaceDescriptors meeting at this vertex, in topological order.
- `adjacent_edges: Vec<usize>` — Indices of ReconEdges meeting at this vertex, in topological order.
- `uv_coords: Vec<(f64, f64)>` — UV coordinates of this vertex on each adjacent face's surface, parallel to `adjacent_faces` (populated in 3.4).
- Per-vertex tolerances live on the shared `TopoDS_Vertex` objects, raised adaptively in 3.4 to cover each vertex's distance from its curve endpoints.

**Corner vertex refinement** (`refine_corner_vertices`, runs after tangency detection): SHARP corner vertices are refined by cyclic projection onto their adjacent fitted surfaces — the raw mesh vertex at a corner is only as accurate as the tessellation puts it (where a fillet terminates, the mesh may have no vertex exactly at the CAD corner). Corners adjacent to any TANGENT edge are excluded: tangent surface pairs intersect along a curve rather than an isolated point, and cyclic projection slides along the tangent direction (observed drifting 0.57 mm). A safety guard keeps the original vertex if the refinement moves more than 2× surface tolerance, and a refined position landing within max(vertex_tolerance, 1e-6mm) of an edge-adjacent corner's current position is rejected (cyclic projection once pulled both endpoints of an edge onto the same near-parallel-plane intersection fixpoint, collapsing the edge).

#### 3.1 Create OCCT surface objects and build adjacency graph

**Surface creation** (`create_occt_surface`):
- For each selected surface, create the corresponding infinite OCCT surface:
    - Planar → `geom::Plane::new_pnt_dir` with `gp::Pnt` at centroid and `gp::Dir` along normal. Wrapped via `to_handle().to_handle_surface()`.
    - Cylindrical → `geom::CylindricalSurface::new_ax3_real` with `gp::Ax3` from axis origin/direction and radius.
    - Spherical → `geom::SphericalSurface::new_ax3_real` with `gp::Ax3` from center and Z-axis direction.
    - Conical → `geom::ConicalSurface::new_ax3_real2` with `gp::Ax3` from apex/direction and half-angle + radius at reference height. The `Ax3` origin is at the cone apex, Z along axis.
    - Toroidal → `geom::ToroidalSurface::new_ax3_real2` with `gp::Ax3` from center/axis direction, major radius $R$, and minor radius $r$.
    - BSpline → a fresh handle to the exact `Geom_BSplineSurface` object fitted in stage 2.3's freeform rung (never a reconstruction — the charts wrap the same object).
- For consistency with OCCT conventions: planes use gp_Ax3 with Z as normal; cylinders have Z along axis; spheres have poles at Z extremes; cones have Z along axis from apex; tori have Z along rotation axis. This affects pcurve computation.

**Adjacency graph construction:**
1. Build face-to-surface map: for each mesh face, record which selected surface it belongs to (`build_face_to_surface_map` — under the new stage 2 this is a re-export of `face_region`).
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
    - Conical: radial direction from axis at the point's axial position, rotated by the cone's half-angle away from the axis, negated if concave.
    - Toroidal: direction from the nearest tube center point to the surface point, negated if concave.
    - BSpline: D1 cross product on the surface handle, oriented outward.
- Compute the dot product of the two normals. If all sampled pairs agree within 2° (dot > cos(2°) ≈ 0.9994), mark the edge as tangent.
- Points too close to a cylinder axis or sphere center (degenerate) cause the edge to be marked non-tangent.

**Decision: do not modify surfaces to enforce tangency.** Modifying analytic surfaces would change the geometry. Instead, tangent edges get special handling in edge curve computation (3.3): construct the edge curve directly rather than relying on surface-surface intersection.

- With the `--compare` flag: verify that each edge's tangency classification matches the reference STEP (warning on mismatch). This is a consistency check, not a geometry check.

#### 3.3 Create edge curves
For each ReconEdge, compute the 3D intersection curve between the two adjacent surfaces, trim it to the vertex endpoints, and store in `edge.curve_3d`.

**Intersection computation (implemented):**
- For non-tangent edges: use `GeomAPI_IntSS` to compute intersection curves with tolerance `vertex_tolerance_mm`. IntSS may return multiple curves (e.g., a plane cutting through a torus); `select_closest_curve()` picks the one closest to the mesh boundary vertices by sampling up to 10 evenly-spaced boundary vertices, projecting each onto each candidate curve via `GeomAPI_ProjectPointOnCurve`, and summing distances. IntSS construction can throw (wrapped in catch_unwind, routed to the chain fallback).
- For tangent edges: construct the curve analytically rather than using `GeomAPI_IntSS`. Dispatches by surface pair type:
  - **Plane-cylinder**: `Geom_Line` parallel to the cylinder axis at the tangent point.
  - **Cylinder-cylinder**: `Geom_Line` along the first cylinder's axis direction.
  - **Sphere-cylinder**: `Geom_Circle` (great circle arc on the sphere). Arc selection samples mesh boundary vertices to determine forward vs. reverse arc on the periodic circle. Open arcs use this path for slight cones too.
  - **Sphere-cone (closed loop)**: the tangent circle is NOT a great circle (it sits at radius R·cos α, offset along the cone axis). Constructed perpendicular to the cone axis, with axial position and radius averaged from the edge's mesh boundary vertices.
  - **Plane-cone**: `Geom_Line` along the cone's generator at the tangent azimuth.
  - **Cylinder-cone** (TODO): `Geom_Circle` at the axial position where the cone radius equals the cylinder radius (when coaxial), or `IntSS`-based (non-coaxial).
  - **Plane-torus**: `Geom_Circle` perpendicular to the torus axis. Circle center and radius computed from mesh boundary vertices projected onto the torus axis. Handles both closed loops and arcs.
  - **Cylinder-torus**: `Geom_Circle` computed by the vertex-averaging algorithm against *both* candidate axes — the torus axis (coaxial fillet junctions) and the cylinder axis (tube-continuation junctions, e.g. a straight pipe leaving a pipe-elbow torus), keeping whichever circle best fits the boundary vertices.
  - **Torus-torus**: `Geom_Circle` perpendicular to the torus axis. Same vertex-averaging algorithm.
  - **Sphere-torus** (TODO): `Geom_Circle` — analogous to sphere-cylinder but the junction is a circle on the sphere.
  - **Fallback**: Construct a line from the two vertex endpoints.

  All tangent edges are trimmed to vertex endpoints using the same parameter projection as non-tangent edges. A `validate_tangent_curve` step checks that boundary vertices lie within `surface_tolerance_mm` of the constructed curve.

**Trimming to vertex endpoints (implemented):**
- For each BRepVertex at an edge endpoint, project the vertex's 3D position onto the intersection curve using `GeomAPI_ProjectPointOnCurve` to get the curve parameter value.
- Trim the curve to the parameter range `[t_start, t_end]` using `Geom_TrimmedCurve` (`candidate_trim_parameters` handles periodic arc selection, domain clamping for bounded candidates, and rejects degenerate spans).
- For closed-loop edges (vertex_indices == usize::MAX), use the full intersection curve parameter range; closed-loop edges carry PERIODIC curves (converted or chain-interpolated as needed — see the periodic-band junction notes in stage 2.3).
- **Partial arc reconstruction**: When `IntSS` returns a partial arc for a closed-loop edge (e.g., a semicircle when a sphere-plane intersection passes through the sphere's UV poles), reconstruct the full circle from 3 sampled points on the arc.
- For non-periodic curves, ensure `t_start < t_end` (swap if needed).

**Boundary-chain validation and repair (`repair_edge_curve_from_boundary`):**
- The mesh boundary chain lies on the true junction within chordal tolerance, so it is the ground truth for every edge curve regardless of surface pair. After trimming, the curve is validated against a spread of chain vertices; when it strays beyond `--surface-tolerance` — the complementary-arc trap on a periodic junction, a degenerate IntSS result, or a straight-line stand-in for a curved junction — the curve is rebuilt as the circle through three spread chain vertices, re-trimmed with a chain-supported arc (`select_arc_parameters`), and adopted if it fits the chain better.
- The chain is also the **last-resort constructor** (`edge_curve_from_boundary_chain`): when both the IntSS and tangent-analytic constructions fail outright, the edge curve is built from the chain directly (vertex-to-vertex line, upgraded to the chain's circle when curved, periodic B-spline interpolation for closed chains, exact circumscribed circle for 3-vertex loops) and accepted if it fits the chain within `--surface-tolerance`.

**Vertex position consistency:**
- OCCT requires all edges at a vertex to share the same `TopoDS_Vertex`. All edges are built from the shared per-BRepVertex `TopoDS_Vertex` objects; vertex tolerances are raised adaptively to cover each vertex's distance from its curve endpoints (stage 3.4), and pcurve endpoints snap to the shared vertex's UV so wires close in 2D at their vertices.

- With the `--compare` flag:
    - For each ReconEdge with a computed `curve_3d`, sample points along the trimmed curve. For each sample point, compute the distance to the nearest edge in the reference STEP file using `BRepExtrema_DistShapeShape`. Report an error if any sample exceeds `--vertex-tolerance`.
    - Additionally, verify that the curve endpoints coincide with BRepVertex positions within `--vertex-tolerance`.

#### 3.4 Create OCCT faces
For each face, construct a `TopoDS_Face` from the surface and the shared edges — **shared-edge topological assembly**. Every ReconEdge is exactly one `TopoDS_Edge` used by both adjacent faces; each face attaches its own pcurve to the shared edge and assembles its wires and the face directly with `BRep_Builder`. `MakeWire`/`MakeFace`-from-wire/`ShapeFix` (which copy edges and would break the sharing) are not used; one uniform mechanism covers planes, cylinders, cones, spheres, tori, and B-splines alike.

**Algorithm (as implemented):**

*Step 1: Shared vertices and edges.*
- One `TopoDS_Vertex` per `BRepVertex` (tolerance `vertex_tolerance_mm`).
- One `TopoDS_Edge` per `ReconEdge` via `BRepBuilderAPI_MakeEdge` from the trimmed 3D curve with explicit parameters. For periodic curves (circles), vertex order is matched to the curve parameterization direction: V1 corresponds to `first_parameter()` and V2 to `last_parameter()`. This is critical because `MakeEdge` strips `TrimmedCurve` wrappers and uses `ElCLib::AdjustPeriodic` on the base periodic curve, which always selects the forward (CCW) arc from V1 to V2.
- **Adaptive vertex tolerance**: before each `MakeEdge`, each vertex's distance to its curve endpoint is computed and the OCCT vertex tolerance raised to cover it (1% margin).
- Micro-edges whose trimmed curve endpoints coincide within `Precision::Confusion` fall back to the straight vertex-to-vertex edge.

*Step 2: Per-face UV analysis (`scan_face`).*
- Each face gets an analytic **chart** (`SurfaceChart`) mirroring the gp_Ax3 frame of its OCCT surface, with closed-form point→(U,V) and (U,V)→point maps for the five analytic surface types; B-spline faces carry a chart wrapping the exact fitted surface via a cached `ShapeAnalysis_Surface`. All UV work — sampling, windows, orientation tests — is analytic; `GeomAPI_ProjectPointOnSurf` is not used in face construction.
- **Sphere gauge choice**: a sphere's axis is a free gauge, so it is pointed at the mesh region's mean centroid direction from the center. The region's boundary then never passes near a UV pole: poles land strictly inside the region (a polar cap — handled by a degenerate pole edge) or strictly outside.
- Each boundary edge's stored curve (read back from the `TopoDS_Edge`, so `MakeEdge` parameter adjustments are respected) is sampled to branch-continuous UV arrays; U values within `pole_radius_eps` of a chart pole are backfilled from the nearest well-defined sample.
- Boundary edges are grouped into loops by `BRepVertex` connectivity and each loop is **branch-aligned**: consecutive edges' UV endpoints are matched by whole-period shifts. The loop's net **winding** around each periodic direction classifies the face: a loop with winding ±1 wraps the surface. U-jumps across pole junctions are recorded as pole gaps rather than winding.
- More than 2 winding loops on one face is a hard error: ≥3 winding loops cannot bound a *connected* region of a periodic surface. Under the new stage 2 this cannot happen for a valid input (regions are connected by invariant); the check remains as an assertion (it fired on 3DBenchy under the abandoned pipeline, whose greedy selection created disconnected "regions").

*Step 3: Seam planning (`plan_seam_anchors`).*
- A face whose boundary winds a periodic direction needs a **seam edge**; its endpoints must be vertices of the winding loops. Anchors: an existing `BRepVertex` on the loop when there is one; otherwise the (single) closed-loop edge is rebuilt to start at the seam azimuth and gains an anchor vertex at both ends. Because the anchor vertex becomes part of the *shared* edge, every face using that edge sees the same anchor — cross-face consistency is forced by sharing, not by tolerances. Faces are planned in waves so chains of wrapped faces (e.g. plane–cylinder–torus stacks) adopt one consistent seam azimuth; a face whose winding loops carry incompatible pre-determined anchors is a hard error (edge splitting is the recorded future extension).

*Step 4: Per-face assembly (`assemble_face`).*
- **Window**: wrapped faces span `[c, c+2π]` in the wrap direction, bordered by the seam at the anchor coordinate; other periodic directions use the circular-gap window over boundary samples with mesh-centroid majority disambiguation (`circular_window`). Winding loops are rotated to start at their anchor and oriented for a counterclockwise composite; plain loops are lifted into the window.
- **Seam edge**: 3D curve = `Geom_Surface::u_iso`/`v_iso` at the seam coordinate (parameterized exactly by the other coordinate, for every surface type), trimmed between the two anchors (or anchor and pole), built with `MakeEdge` + explicit parameters. It carries **two pcurves** (`BRep_Builder::update_edge` two-pcurve variant): PC1 for its FORWARD occurrence at the window's high border, PC2 for REVERSED at the low border (OCCT resolves FORWARD→PCurve1; verified against a `BRepPrimAPI_MakeCylinder` seam).
- **Degenerate pole edges**: at a cone apex or sphere pole the UV boundary closes across the pole with a degenerate edge — no 3D curve, `degenerated` flag, a constant-V pcurve line spanning the U gap, the pole vertex at both ends. This is the *only* pole mechanism, identical for cones and spheres.
- **Pcurves** (`build_pcurve`): for every (edge, face) incidence, UV samples of the edge's stored curve are interpolated (`Geom2dAPI_Interpolate`) with the curve's own parameters, so SameParameter holds by construction. Endpoint samples snap to the shared vertex's analytic UV (branch-lifted), making wires closed in 2D at their vertices even when a mesh vertex sits slightly off this face's fitted surface. Near-duplicate UV samples are deduplicated, and a fully pointlike projection gets a degree-1 two-pole Geom2d B-spline pcurve with the edge's 3D extent folded into the returned deviation. Seam pcurves pin the wrap coordinate to the window borders exactly. Sampling densifies (33→257 points) until the pcurve reproduces the 3D curve within `vertex_tolerance/10` or the deviation plateaus; the achieved deviation — measured at samples *and* midpoints, *including endpoints* — sets the edge tolerance.
- **Wires and orientation**: the composite wire of a wrapped face is emitted in counterclockwise window order. Plain loops are oriented by the mesh ground truth: the loop that winds around interior mesh-centroid UV points is the outer wire (counterclockwise); loops containing no interior points are holes (clockwise). Wires are `BRep_Builder`-built from the shared edges with per-occurrence orientations; the face is `BRep_Builder::make_face` + wires + `natural_restriction(false)`.
- **Faces with no boundary edges** (full sphere, complete torus) use `MakeFace` natural bounds — OCCT generates their seam/degenerate edges and there are no shared edges to preserve.
- A face is flagged **concave** when its hypothesis is concave: its natural surface normal points into the material, so stage 3.5 adds it to the shell REVERSED.

*Step 5: Validate and compare.*
- Validate each face with `BRepCheck_Analyzer`; invalid faces get detailed `BRepCheck_Face` diagnostics plus a structural dump under `-v`. Note `BRepCheck_Analyzer::IsValid(face)` reflects subshape statuses *in the face's context* (e.g. `InvalidCurveOnSurface`), which is why the edge tolerance accounting above matters.
- With `--compare`: for each face, sample a representative mesh face centroid and compute distance to the reference STEP shape using `BRepExtrema_DistShapeShape`. Report face count comparison. Error if max distance exceeds `--surface-tolerance`.

#### 3.5 Construct shells
Assemble faces into shells directly from the shared-edge topology. Faces already share their `TopoDS_Edge` objects, so shells are **closed by construction** — `BRepBuilderAPI_Sewing` is not used, and `ShapeAnalysis_Shell` remains as validation only.

**Algorithm (as implemented):**
- Group faces into connected components over the ReconEdge adjacency (multi-solid models produce one shell per component, e.g. hollow_cube's outer and void shells).
- For each component, `BRep_Builder::make_shell` and add every face — FORWARD for convex/planar faces, REVERSED for concave ones — then set the shell's closed flag. Orientation is consistent outward by construction: mesh normals point outward (stage 1 validation) and the hypothesis convexity flags encode each surface's relation to the material.
- Validate with `ShapeAnalysis_Shell::CheckOrientedShells` (free edges, orientation inconsistencies) and report.

- With the `--compare` flag: count shells in the reference STEP; error on free edges; report orientation check results.

#### 3.6 Construct solids
Convert closed shells into `TopoDS_Solid` objects.

**Algorithm:**
- For each shell, use `ShapeFix_Solid::SolidFromShell` which handles orientation automatically.
- Apply `ShapeFix_Shape::Perform()` for comprehensive face/wire/edge/shell/solid fixes.
- Call `BRepLib::SameParameter` with `forced=false` and tight tolerance (`vertex_tolerance_mm`) to ensure pcurve/3D curve agreement without corrupting already-correct pcurves. Using `forced=true` or a large tolerance recomputes all pcurves and can corrupt faces with closely-spaced tangent edges, causing BRepCheck self-intersection failures.
- Call `BRepLib::UpdateTolerances` and `BRepLib::OrientClosedSolid`.
- Validate with `BRepCheck_Analyzer`. Log warning if validation fails.
- Compute and report volume using `BRepGProp::VolumeProperties`.
- **Multi-shell solids (internal voids)**: after per-shell solids are built, `merge_void_solids` finds solids whose bounding box is contained in a strictly larger solid's box; each such void's shells are added to the container solid (`BRep_Builder::Add`), trying both shell orientations and keeping the one whose resulting volume equals container − Σ voids.

- With the `--compare` flag:
    - Count solids in the reference STEP file.
    - Compute volumes of both constructed and reference solids. Match by sorted volume. Report a hard **error** if relative volume difference exceeds 1% (a warning until 2026-07-08; a warning-only check let a 46% volume regression slip past a green test suite). Most models achieve volume agreement within ~1e-7 relative difference.
    - Compute `BRepExtrema_DistShapeShape` between each constructed solid and the reference STEP shape. Report error if distance exceeds `--surface-tolerance`.

#### 3.x Design note: shared-edge topological assembly instead of geometric sewing

**Status: IMPLEMENTED (2026-07-08).** Stages 3.4/3.5 above describe the result; this note records the motivation, since it reshaped stage 3's architecture. Stage 3.4 used to build each face independently (`MakeWire`/`MakeFace`/UV-bounds/`ShapeFix` paths) and stage 3.5 reconciled them with `BRepBuilderAPI_Sewing`, a *geometric* stitch that merges edges by proximity within a tolerance. `MakeFace`/`MakeWire`/`ShapeFix` copy edges, so edge sharing was defeated by assembly time and sewing was load-bearing. Every hard case — sphere pole caps, torus seams, cone apex caps, tangent-edge tolerance relaxations — was a symptom of trimming faces independently and matching them by geometric tolerance. The rework made each edge one shared object with per-face pcurves, made seams/degenerate edges first-class topology, and assembled everything with `BRep_Builder`; shells became closed by construction and `cone_tangent_planes` went from 60% volume error to ~1e-7 relative.

**Hard OCCT constraints (established empirically — do not relearn by crashing).**
- `ShapeFix_*` can **segfault on cone faces** (pcurve/edge-curve interaction). The assembler attaches pcurves directly with `BRep_Builder` and never routes faces through ShapeFix. (Stage 3.6's `ShapeFix_Solid`/`ShapeFix_Shape` on assembled solids has been safe.)
- `MakeEdge` on a periodic (circle) curve forces p1 < p2 and takes the forward arc V1→V2, so **vertex order selects which arc**. The same trap applies when *trimming* periodic curves in stage 3.3; the boundary-chain repair is the safety net.
- `BRep_Tool::CurveOnSurface` resolves a seam edge's FORWARD occurrence to PCurve1 and REVERSED to PCurve2; OCCT's own cylinder puts the FORWARD occurrence at the high U border.
- `Geom_Surface::u_iso(u)`/`v_iso(v)` return curves parameterized exactly by the other surface coordinate (verified in ElSLib for cone/sphere/torus/cylinder) — this is what makes seam construction surface-type-agnostic.
- Thread-safety of the OCCT calls used here: LEARNINGS "OCCT Thread Safety" (BRepExtrema, MakeEdge/Face/Vertex, ShapeFix_* are all safe; only STEP Reader/Writer need serializing).

---

### Stage 4: Output

#### 4.1 Output objects
Write constructed solid(s) to a STEP file.

**Algorithm (implemented):**
- Build a `topo_ds::Compound` containing all solids from stage 3.6.
- Set STEP metadata: AP214 schema, product name "brepper" via `Interface_Static::set_c_val`.
- Create a `STEPControl_Writer`, transfer the compound with `StepModelType::Asis` mode to preserve exact geometry.
- Write the output file. Report error if transfer or write returns non-Done status.
- With the `--compare` flag (each check is a hard error as of 2026-07-08):
    - Re-read the written STEP file via `STEPControl_Reader` to validate the round-trip. Error if unreadable (note this makes `-o /dev/null` incompatible with `--compare`; `run_test_table.sh` writes a real temp file).
    - Compare volumes of the output shape vs reference shape using `BRepGProp::VolumeProperties`. Error on >1% relative difference.
    - Compute maximum distance between output and reference using `BRepExtrema_DistShapeShape`. Error if exceeds `--surface-tolerance`.

---

## Rewrite Roadmap

**One-shot from the current code state** (commit ce85dd0: jet-based
2.1/2.2 with the retired jump classifier still in the tree). No
compatibility path, no threshold gating — the jump machinery is deleted in
step 1, and the old-BFS-era integration tests stay red until step 7
re-pins them. Function-level guidance lives in the stage-2 section's
"Implementation hints" (after the CLI notes), including per-model success
criteria under "Acceptance narratives".

1. **Demote the jet, delete the jump** (`edges.rs`): sharp-only
   classification; one-pass `JetField` (κ1, κ2, principal directions,
   medial points) on `Segmentation`; keep the three jet unit tests; re-pin
   `tests/stage2_segmentation.rs` to smooth-region counts (measure with
   `dump_regions`, don't guess); update `dump_regions` and viz (the jump
   orange overlay dies). **[DONE]**
2. **Stage 2.3(a) whole-region ladder** (`classify.rs`): existing fitters,
   terminal per-facet-plane fallback from day one — this un-blocks stages
   3/4 end-to-end and re-arms full-pipeline `--compare` on every
   sharp-edged model immediately. (Jet-pruned rung order deferred as a perf
   optimization; strict DOF order — plane→sphere→conical→torus,
   first-validates-wins — is already efficient via plane-first.) **[DONE]**
3. **Stage 2.4 resolution/merge/relabel + invariants** (`repair.rs`) —
   built BEFORE competition so single-surface models flow end-to-end with
   junction singletons absorbed; ported `refine_boundary_assignments` as
   `relabel_boundaries`. Merge is gated (small S into a larger real T, >2
   faces) so the per-facet fallback storm never chains; no-op and
   invariants-clean on the current corpus, with unit tests exercising
   merge/no-merge/re-split. **[DONE]**
4. **Stage 2.3(b) scoped competition** (`compete.rs`): seeds, per-seed
   rung competition, claims, resolution wiring. Acceptance: nosecone →
   cone+sphere+base (3, matching the reference's 3 faces); rounded_cube
   shell decomposes into 6+12+8; stepped_block. **[DONE for the analytic
   surfaces — plane/sphere/cylinder/cone.]** Notes from implementation:
   growth floods at ~30× vertex-tolerance (a fit from partial support has
   ~1e-4 mm residual, well below the ~0.1 mm a neighbour diverges past a
   tangent junction; every committed hypothesis is refit on full support to
   ~1e-6 mm, and claim resolution corrects any one-row leak); the seed patch
   is depth-limited so it does not straddle a tangent junction between two
   same-class developables (a cone meeting a cylinder) — and, since
   2026-07-25, also bounded in space at `PATCH_RADII` = 4 radii of curvature
   of the seed's own jet, because a hop limit stops bounding the patch once
   the facets are large (three hops walk ~18 mm along a fillet chain swept
   with 6 mm facets, straddling several of its segments); and a region whose
   competition leaves a large per-facet leftover blob (a whole sub-surface
   it could not fit — e.g. a coarse fillet, or a concave torus fillet
   awaiting step 5) is discarded to a clean per-facet polyhedron rather than
   reconstruct invalid wires at the analytic/per-facet boundary. Torus
   fillets (cone_cylinder, filleted_*, stepped_block_toroidal, pipe_elbow)
   correctly fall back to that polyhedron until step 5.
5. **Stage 2.3(c) promotion** (bent-tube `k = 1/R` + jet triggers +
   medial-point spine seeding). Acceptance: pipe_elbow → torus + 2
   cylinders with no strip machinery. **[IMPLEMENTED — awaiting the debug
   session.]** As built, the cheap path sufficed (no `k = 1/R` degenerate
   model, no medial-spine fit yet — both remain contingency if strip/coarse
   meshes defeat this): the torus competes as a first-class candidate in
   `compete.rs` — (a) grown directly from any curved non-sphere-like seed
   (|κ2| < 0.8·|κ1|; recognizes fillet tori and bend tubes generally, not
   just cylinder-adjacent bends), and (b) *promoted* from a cylinder whose
   flood stalled against tube-continuation jet (κ1 ≈ 1/r, |κ2| ≤ 0.8·κ1),
   fitting the existing 7-DOF torus fitter to the pure-bend patch just past
   the stall line — the cylinder fit itself locates the tangency, so the
   patch cannot straddle it (jet-class patches cannot tell bend from
   straight tube: same κ1). Giant-major-radius degenerate fits (a straight
   segment on an absurdly distant torus) are rejected at `R > bbox
   diagonal`, mirroring the ladder acceptance test. Smoke:
   pipe_elbow_10_medium → 2 planes + 2 cylinders + 1 torus, 1 valid solid,
   0 per-facet leftovers; rounded_cube/nosecone competition tests still
   green.
6. **Freeform rewiring** (`bspline.rs`: leftover components +
   island-dissolution guard). Acceptance: tests/nurbs suite green. **[DONE.]**
   As built: `fit_leftover_components` (a NEW front end; the old region
   *gathering* still lives beside it until the step-7 deletion) fits the
   unclaimed leftover components of each competed region. The freeform
   *domain* is the leftover faces PLUS the small (≤ 16-face) analytic islands
   competition committed — a singly-curved patch grows a short accidental
   cylinder band (the arch) and an inflection grows accidental coplanar
   strips (the wave), each *fragmenting* the region into separate leftover
   pieces; growing connected components over that domain reunites them.
   Per component: if the unexplained (leftover) area dominates
   (`absorbed ≤ single`), the islands are absorbed into the fit and, if it
   validates, dissolved (dropped by the caller); otherwise the analytic fits
   are trusted and only the leftover freeform is fitted. Components whose
   leftover part is below `MIN_BSPLINE_FACES` (4), and fits that fail at
   every grid density, degrade to per-facet planes.
   **Revised 2026-07-25** on three measured failures: (i) only *enclosed*
   islands join the domain — one bordering a surface competition explained at
   scale is a real small feature, not an accident adrift in freeform, and
   without this a densely-filleted model (stepped_block_toroidal) pulls its
   whole decomposition into one nonsensical domain; (ii) the islands are also
   absorbed when dropping them *disconnects* the leftovers into pieces that
   are each patch-sized (≥ `MIN_BSPLINE_FACES`) — islands cutting one freeform
   surface into several patches are exactly the accident this rung undoes,
   whereas real analytic surfaces separated by leftover *slivers*
   (part_rounded_cube's four fillet cylinders) must stay committed;
   (iii) the trust-the-analytic branch fits each connected leftover piece
   separately, because `build_region` needs a connected face set — fitting a
   disconnected one produced boundary loops bounding different pieces and the
   spurious error `closed band: no path between boundary loops`. `face_region` (region id
   per mesh face) is the outside-surface key for boundary-corner detection.
   The old `LEFTOVER_BLOB_MAX` whole-region discard is now the *post-freeform*
   safety net: it fires only when a committed analytic surface survives AND
   freeform could not clear a big leftover blob — so a partially-decomposed
   model (fillet tori + an undecomposed blob) keeps its analytic surfaces
   once freeform gives the blob a home, instead of discarding them
   (filleted_hole_block, ccad_cone_cylinder params, filleted_cylinder tori
   all recovered). B-spline surfaces are inert in stage-2.4 repair (never a
   relabel/merge/re-split target — a fit vouches only for its own data).
   Net vs the pre-step-6 baseline, zero regressions: stage2_integration
   132→139, stage3 99→101, stage4 28→31; lib 90, stage1 83, seg 14,
   competition 3 all green. Remaining stage3/4 failures are the documented
   pre-existing set (part_rounded_cube + rounded_cube coarse; ccad
   cone_cylinder/filleted_cylinder/filleted_pipe reference mismatches;
   stepped_block_toroidal coarse/medium LEFTOVER_BLOB).
7. **Deletion pass + corpus re-enable** — split into 7a (corpus re-enable)
   and 7b (deletion refactor).
   - **7a corpus re-enable [DONE].** The stage-2/3/4 integration tests are
     re-pinned to the new pipeline: their old-sub-stage helpers now run the
     full pipeline (`--stage=2`), which greened the sharp-edged/analytic
     count + coverage tests unchanged (cube = 6 planar, block_with_hole =
     6 planar + 1 cylindrical, pipe = 2 + 2 cyl, two_holes = 6 + 2 cyl,
     chamfered_cube = 26 planar, …). Suite is fully green (0 failures);
     the models the new pipeline does not yet decompose as the old BFS did
     are **disabled with a blocker comment** (`#[ignore]`/`ignore =`), per
     the migration policy, NOT nerfed with loosened gates. Disabled
     blockers, all PRE-EXISTING (steps 2–5 classification quality, not
     step 6): (a) **hole → sphere** — a cylindrical hole tessellated as a few
     rings is vertex-ambiguous with a sphere, which simplest-DOF-first
     wrongly committed (plate_with_hole, fusion plate, l_bracket) — **FIXED
     2026-07-23 by jet-informed rung ordering** (`region_is_developable` in
     classify.rs tries the conical rung before the sphere rung for
     developable regions); the 9 fusion `plate_with_hole` tests were
     re-enabled; l_bracket (a second cylinder undetected) and the coarse
     onshape plate (hole faceted, no cylinder region) stay disabled for
     different gaps; (b) **fillet over-segmentation** — a fillet
     chain decomposed into extra tori + a freeform patch with the straight
     body absorbed (filleted_cylinder/hole_block/pipe) — **FIXED 2026-07-24**
     by four competition fixes in compete.rs (magnitude-ordered principal
     curvatures so a concave hole wall reads curved; seed-patch conditioning
     that trims tangent-junction bleed; DOF-aware claim resolution; fat-torus
     `r ≤ 0.7·R` rejection); all 13 filleted tests re-enabled and passing
     `--compare`; (c) **cone/cylinder
     + concave-fillet** decomposition differs from the old BFS
     (cone_cylinder, cone_tangent_planes); (d) **coarse-mesh leftover blob**
     exceeds `LEFTOVER_BLOB_MAX` → per-facet fallback (part_rounded_cube,
     stepped_block_toroidal coarse/medium); (e) **ccad reference
     mismatches** (cone_cylinder/filleted_* stage-4 `--compare`). See TODO.md
     for the current re-enable list.
   - **7b deletion refactor [DONE].** `bfs.rs` is deleted; its `SurfaceType`
     trait — trimmed to the methods the ladder/competition use (`fit_surface`,
     `vertex_distance`, `centroid_distance`, `convexity_check`,
     `determine_convexity`, `make_hypothesis`) — moved to `surface.rs`. Gone:
     `generate_hypotheses`/`BfsConfig`/`Seed`/`HypothesisFaces`;
     `stage2_bfs_pipeline` + `select_surfaces`/`refine_boundary_assignments`/
     `absorb_redundant_surfaces`/`prune_unselected_hypotheses`/
     `compare_selected_surfaces`/`compute_seed_scores` and the old-BFS unit
     tests in `mod.rs`; the old `deduce_bspline_surfaces` gathering front end;
     `--skip-stage`; the `MeshFace.{planar,conical,spherical,toroidal}_hypothesis`
     fields + the `UNDEDUCED_*`/`NO_HYPOTHESIS` constants; and each fitter's
     BFS-only trait methods (`generate_seeds`, `hypothesis_field`,
     `acceptance_test`, `seed_*`, `viz_overlays`, …) and their dead helpers.
     Net −3860/+37 lines across 12 files; behaviour-preserving — the
     integration suites are byte-identical to 7a (stage2 159/26, stage3 101/5,
     stage4 31/6), only the deleted BFS unit tests drop out of the lib count
     (90→74). clippy clean. (Deleting the BFS engine cascaded through the
     shared trait and the MeshFace fields as expected — no clippy-clean
     partial deletion existed.)

Expected wins to verify along the way: pipe_elbow's torus without
strip-merge seeding; nosecone with no absorb pass; stage-2 runtime on
large models collapsing (whole-region fits + one competition per smooth
region; full stage 2 on 3DBenchy was ~15 min under BFS vs ~3 s
planar-only).

---

## Stage-2 Retrospective (2026-07-19, plus the 2026-07-20 Resolution)

Written for whoever picks this up next. Honest account of what the
partition-then-classify rewrite achieved, what it cost, and whether it
looked like the right bet — including vs. the then-Contingency
(seed-first competition, since promoted into §2.3 by the Resolution at
the end of this section). The body below is the 2026-07-19 state, kept
as history; **the Resolution supersedes its recommendations.**

### What exists and works

- **Stages 2.1 + 2.2 are implemented and tested** (`src/stage2/edges.rs`,
  `regions.rs`, `tests/stage2_segmentation.rs` — 14 tests). 2.3–2.5 are
  NOT implemented; `stage2()` stops at the 2.2 frontier and the old BFS
  pipeline sits nerfed in place (`stage2_bfs_pipeline`, uncalled) until
  the deletion pass. **The old stage-2/3/4 integration tests are red by
  design** — do not "fix" them by reviving the old pipeline.
- **The curvature estimator is now principled, and region interiors are
  essentially exact.** Per-face anchored osculating-jet: quadratic
  height-field fitted to vertex *positions* over a PCA frame of the
  own+1-ring point cloud, anchored on the face's own vertices (weight
  50), no outlier rejection, inside a hysteresis fixpoint. On clean
  regions the descriptors are right to a few %: plane (0,0), cylinder
  (1/r, 0), sphere (1/r, 1/r), torus (1/r_tube, k2(ψ) with the correct
  sign change at the inner equator), cone (1/ρ(h), 0). Verified directly
  with `dump_regions` probes on the corpus.
- **The redesign's headline structural win is real.** The pipe-elbow
  torus bend arrives as ONE region with zero strip machinery — the
  strip problem (local seeds can't condition high-DOF fits) is gone,
  because the partition, not the seed, determines fit support.
- Measured state of the corpus: all sharp-edged analytic models exact
  (cube 6, block_with_hole 7, pipe 4, two_holes 8, spherical_pocket 7,
  chamfered_cube 26); pill_fine 3; pill_coarse 3 true + junction
  singletons; pipe_elbow bend one region (+ merged walls); cone_tangent
  11 (was ~160 with earlier estimators); filleted_cylinder 21;
  filleted_pipe/hole_block ~95; nosecone ~130; stepped_block ~260;
  rounded_cube_fine 26 true regions + small confetti; rounded_cube
  medium/coarse ~50-70.

### What was hard, and what it means

- **The estimator was the whole battle** — five variants were tried
  (per-edge κ=θ/d → neighbor-normal tensor → jet+outlier-rejection →
  +weights/thresholds → final jet, no rejection). Two durable lessons:
  (1) use vertex *positions*, not normals — a skinny strip/sliver
  triangle's plane is a chord whose normal is degrees off the true
  tangent, poisoning any normal-based estimate; (2) outlier rejection
  *backfires on tangent junctions* — it forces arbitrary one-side reads
  that vary face to face and fragment the junction. The own-vertex
  anchor alone breaks the junction symmetry. Do not re-add rejection.
- **Junction-circle faces are inherently ambiguous, and no local
  estimator fixes that.** A face whose vertices sit on a tangent seam
  genuinely belongs to the boundary, not to either side; on pill the
  seam-circle neighborhood is even planar by construction, so any
  position fit reads (0,0). These become singleton regions. This is
  EXPECTED and is the 2.4 merge repair's designed input — it is not a
  classifier bug to keep polishing.
- **Similar-radii tangent junctions still miss.** pipe_elbow's outer
  torus↔cylinder junction (σ_min jump 0.067/0.2 = 0.33 < FIRE 0.35)
  doesn't fire, so the cylinder walls merge with the bend. This is the
  documented "similar radii" case for 2.4 split repair, not an oversight.
- **Knob count is a real risk.** The classifier still has FIRE_RATIO,
  UNFIRE_RATIO, FLOOR_COEFF, OWN_WEIGHT, MAX_ROUNDS — fewer and more
  principled than the BFS pipeline's knobs, but the same *kind* of
  machinery. If stage 2.4 accumulates more, this architecture is
  re-inventing the heuristic treadmill it was created to escape.
- The estimator is O(N·ring) per round × ~3 rounds and comfortably fast
  (rounded_cube_fine 10.8k faces: ~0.2 s; pipe_elbow: ~0.02 s).

### Segmentation-first vs. the Contingency — an honest verdict

The Contingency (seed-first multi-type competition) places boundaries by
**fit validation at vertex tolerance** — exact, no curvature estimation
needed, so tangent junctions and junction circles come out clean by
construction. Segmentation-first places boundaries by **curvature
estimation**, which is fundamentally noisy at junction rows and undefined
on junction circles — producing ±1-row boundary error and seam confetti
that the (still unwritten) 2.4 repair loop must absorb.

On the other hand, segmentation-first's **fit support is maximal**:
classification fits see whole regions, which is what dissolves the
strip problem (the Contingency still has it — from any elbow seed the
torus rung is fed a 1-ring and loses to a cylinder strip — unless the
degenerate-limit growth model is built too). And segmentation-first
keeps freeform first-class (a freeform region is a region; the
Contingency leaves freeform as leftover-gathering).

**Verdict: they converge, and the deciding evidence doesn't exist yet.**
The Contingency wins on boundary exactness; segmentation-first wins on
fit conditioning and freeform. What segmentation-first has already
proven: region interiors are exact and the strip problem is gone. What
it has NOT yet proven: that the 2.4 repair loop can absorb the junction
fringe cheaply and without a new pile of knobs. **Build 2.3/2.4 next —
that is the make-or-break test of this architecture.** If repair turns
out fragile or knob-heavy, retreat to the Contingency (its 2.4 repair
loop, invariants, and 2.5 output contract carry over unchanged, and the
degenerate-limit growth model fixes its strip problem); do NOT spend
more effort tuning the edge classifier first.

### What to build next (in order) — SUPERSEDED by the Resolution below

1. **Stage 2.3 classification ladder** (`src/stage2/classify.rs`), plane
   rung first — it enables the coplanar merge and immediately validates
   the repair loop on planar fragmentation (planar caps fragmented at
   junction zones merge trivially). Reuse the existing fitters in
   `planar.rs`/`cylindrical.rs`/`conical.rs`/`spherical.rs`/`toroidal.rs`
   (their fitting/validation functions; their seeding/BFS entry points
   go). Terminal per-facet-plane fallback from day one.
2. **Stage 2.4 repair loop** (`src/stage2/repair.rs`): merge (absorb
   regions whose vertices validate within `--vertex-tolerance` against an
   adjacent larger region's surface — this eats the junction-circle
   singletons and junction-zone confetti), split (re-run 2.1 at lowered
   thresholds inside regions that fail the whole ladder — this is what
   recovers pipe_elbow's outer junction and any residual over-merge),
   boundary relabel (±1 row, the old 2.9 port). Enforce the three
   invariants (total coverage, connectivity, per-region validity).
3. Then wire the freeform rung (bspline.rs, gathering deleted), the
   deletion pass (`bfs.rs`, `stage2_bfs_pipeline`, `--skip-stage`,
   greedy selection, absorb, MeshFace hypothesis fields), and re-enable
   the corpus per Rewrite Roadmap step 6.
4. Useful diagnostics while working: `cargo build --bin dump_regions`
   (modes: summary, singletons, adj, probe, geom, leaks, rdesc, jet);
   `--viz-stages=2.1,2.2` shows sharp (red) / jump (orange) edges and
   region colors.

### Resolution (2026-07-20): pivot decided

The "build 2.3/2.4 next as the make-or-break test" recommendation above
was overtaken by evidence. Direct measurement with `dump_regions` on the
corpus:

- **nosecone_medium**: the sphere is ONE clean 644-face region reading
  σ = (0.080, 0.080) = 1/12.5 — and the cone is **129 fragments whose
  descriptors show a clean gradient**: (0.00, 0.00) at the wide base →
  (0.02, 0) mid-cone → (0.08, 0) near the tip. The estimator measures
  correctly — κ_max = 1/ρ(h) genuinely varies along a cone — and a
  threshold on first differences of curvature fires on real gradients. A
  false-positive *class*, not a bug.
- The same statistic must stay blind at 0.35 to those gradient steps
  while firing at 0.33 for pipe_elbow's genuine torus↔cylinder junction:
  **the two classes overlap in the decision statistic. No threshold
  exists**; calibration cannot fix a statistic that does not separate the
  classes.
- The hysteresis fixpoint hits its 10-round cap (oscillates) on every
  fragmenting model (nosecone 130 regions, cone_truncated 14,
  filleted_pipe 94 regions / 670 jump edges); it converges only where
  jump = 0.
- Five estimator variants, each fixing one corpus class and breaking
  another, is the signature of a signal that is not in the local data at
  junction-row scale. Meanwhile everything decided by SHARP edges alone
  was exact — one round, zero new knobs.

**Decision: scoped pivot** — the architecture in §2 above (second
revision). What carries forward: the sharp-edge partition, region
extraction, the jet estimator (demoted to the §2.1 effort-allocation
prior), viz, `dump_regions`, the segmentation tests, and the invariants.
What changed: tangent-junction discovery moved to fit competition — the
former "Contingency: Seed-First Multi-Type Competition" section was
**promoted into §2.3(b)–(c)** (its standalone section is dissolved; its
sub-region-lemma soundness argument, claim semantics, and the
degenerate-limit growth model all live there now), scoped inside smooth
regions, which fixes its two stated cost concerns (bounded arena, no
cross-sharp-edge accidents). The curvature-jump classifier moved to
Abandoned Approaches. Deciding rationale: in every observed failure mode,
the recovery that works is fit validation at vertex tolerance — curvature
is reliable only where dihedral angles already decide (sharp edges) or in
roles where errors cost time instead of correctness (the prior).

---

## Abandoned Approaches

Retired with the 2026-07-17 stage-2 redesign and the 2026-07-20 fit-driven
pivot. Recorded here because each worked well enough to ship models — but
each existed to compensate for a structural defect elsewhere, and the
compensations kept multiplying (that treadmill is why the architecture
changed, twice). Do not re-add one without identifying which design
principle it would compensate for.

**Fit-free curvature-jump edge classification** (built 2026-07-18/19,
retired 2026-07-20): the attempt to detect tangent (G1-smooth,
G2-discontinuous) junctions from local curvature-descriptor differences,
without fitting. Five estimator variants (per-edge κ = θ/d →
neighbor-normal tensor → jet with outlier rejection → re-weighted jet →
anchored jet without rejection) each fixed one corpus class and broke
another. Terminal evidence (measurements in the Stage-2 Retrospective's
Resolution): the decision statistic — a threshold on local differences of
principal curvatures — cannot separate real curvature *gradients* on one
surface (a cone's κ_max = 1/ρ(h); gradient steps read > 0.35 and
shattered nosecone's cone into 129 fragments while its sphere stayed one
clean 644-face region) from real *junctions* between surfaces with
similar curvatures (pipe_elbow's torus↔cylinder junction reads 0.33,
missed). The classes overlap, so no threshold exists, and the hysteresis
fixpoint oscillated to its round cap on every fragmenting model. What
survives: the sharp-dihedral rule (exact, load-bearing in §2.1(a)) and
the anchored osculating-jet estimator itself — an accurate
interior-curvature field, demoted to the §2.1(b) prior whose errors can
only cost time (Design Principle #4). Do not resurrect a fit-free
junction detector without a statistic that provably separates gradients
from jumps on staggered tessellations.

**Per-type BFS hypothesis generation** (the `SurfaceType` trait framework in
`bfs.rs`: per-face seeding, trial BFS growth with refit-on-near-miss
(`REFIT_SKIP_MULTIPLIER` band), best-trial-by-area commit). The core defect
was *commit-during-generation*: the first hypothesis of a type to claim a
face owned it, so an accidental sub-region fit committed early fragmented
the real surface — a 6-face r=12.4992 cylinder band on the nosecone's
r=12.5 sphere; a 31-face split fragment of nosecone_truncated's real cone
that left a fake corner vertex 80.6 mm off the model until a later pass
rejoined it. Fit-driven growth also *invited* accidental fits (every
quad-row band on a sphere lies exactly on some cone and some torus) and paid
for them at scale. What survives: the fitting math and LM refinement
(SURFACE_FITTING.md), whole-region refit-validate retry, and the convexity
logic (as 2.3's orientation-coherence check).

**Stage execution ordering** (planar → spherical → toroidal → conical →
selection, "most-constrained first"). Existed purely so the conical stage
would not tile spheres and torus tubes with exactly-fitting narrow band
hypotheses, starving the later stages' seeding. With per-region
classification there is no cross-type competition to referee.

**Curvature-guided seed ordering, junction mask, and the large_sphere
mask.** The curvature-difference score (κ = θ/d, neighbor difference,
local-relative 50% threshold) was computed per *face* and used only to order
seeds and exclude junction faces from seed neighborhoods; the large_sphere
mask kept torus seeding off sphere interiors (measured 2026-07-15: without
it, bogus-but-valid torus tilings cost 120–358× in stage-2 runtime —
pill_fine 1.2 s → 429 s). The curvature math is *promoted*, not abandoned:
it is now the 2.1 edge classifier (extended to a signed two-sided
descriptor). The ordering and masking roles are moot — there are no seeds.

**Seed topologies** (conical 3-face fan with optional f3 for 6-vertex
coverage; spherical 1-ring vertex neighborhood; toroidal 1-ring +
3-hop-expansion targeting ≥16 faces; strip-merge torus seeds, already
deleted 2026-07-15 as unreachable). All existed because fits from tiny face
sets are ill-conditioned and had to be grown carefully. Region-level fits
see all faces at once; the conditioning problem is gone (pipe_elbow's torus
is the emblematic case).

**Global acceptance tests**: cylinder/cone angular coverage (second-largest
gap ≤ span/3, plus the coarse-fillet `max_gap` floor gated to cylinders
already spanning ≥45°), sphere solid-angle coverage (eigenvalue ratio
λ3/λ1 of the centroid-direction covariance), torus minor-circle coverage
from face normals (tube span ≥ 0.3 rad), and normal-axis consistency
(per-face |n·axis| deviation ≤ 0.05 during BFS; aggregate variance ≤ 0.0005
at acceptance). Every one is a detector for "this hypothesis is an
accidental sub-region fit on a different surface type" — a class of
hypothesis the new pipeline never constructs, because fits only see maximal
regions and a cone cannot validate an entire sphere. Their rationale
documentation (angular-cluster reasoning, the medium-resolution-fillet gap
analysis) remains valuable context if tiny-region ambiguity ever needs a
diagnostic.

**Greedy area-based surface selection** (old stage 2.6): resolved
overlapping hypotheses (a face could hold one hypothesis per type) largest
total-area first. Obsolete: the partition is disjoint by construction, so
there is nothing to select between. It was also the connectivity-invariant
violator — it kept only a winner's *still-unassigned remnant* with no
connectivity recheck, which is how 3DBenchy face 5841 reached stage 3.4
claiming three disconnected mesh regions (the ≥3-winding-loop abort).

**Redundant-surface absorption** (old stage 2.8, `absorb_redundant_surfaces`):
merged a selected surface into an adjacent larger one under a *tight* rule
(10× vertex tolerance) and a *loose* rule (`--surface-tolerance`, ≤2-face
orphans). The tight rule's targets — spurious sub-fits and split fragments —
cannot form under per-region classification; its sound criterion lives on as
the 2.4 merge. The loose rule was condemned by measurement before the
redesign: its motivating "tangent-junction bridge facet" story was tested
and found false (tessellators emit seam edges even at tangencies), its only
corpus firings absorbed *real* narrow planar faces 0.12–0.27 mm off, and in
planar-only 3DBenchy it introduced up to 0.4 mm of chordal slop (794
BRepCheck face failures, 1% volume loss).

**Per-face multi-hypothesis bookkeeping** (`planar_hypothesis`,
`conical_hypothesis`, `spherical_hypothesis`, `toroidal_hypothesis` fields
on `MeshFace`, with `UNDEDUCED`/`NO_HYPOTHESIS` sentinels). Existed to let
overlapping hypotheses coexist until selection. Replaced by the single
`face_region` map in `Stage2Output`.

**NURBS region gathering** (old stage 2.7 steps 1–3): connected-component
gathering of single-face planar *selections* over smooth edges, absorption
of ≤16-face analytic selections (freeform surfaces bait the analytic
detectors locally — accidental coplanar strips, cylinder ribbons on
saddles), the two protective guards, and the tangent-plane merge of
wholly-smooth-boundary planar neighbors. All of it reconstructed "the
freeform region" from the wreckage of analytic over-detection; the partition
now hands the freeform rung an intact region directly, and exactly-planar
margins of freeform surfaces land in the same region by curvature
continuity.

**Boundary refinement is NOT abandoned** — `refine_boundary_assignments`
(old stage 2.9) is promoted into the 2.4 repair loop unchanged in
mechanism, upgraded in role. Corner vertex refinement likewise survives in
stage 3.

---

## Visualization (`--viz-stages`)

The `--viz-stages` flag enables interactive 3D visualization for debugging pipeline stages. The companion `--viz-stages-level` flag (default 1) controls pause granularity.

Supported stages after the rewrite:
- **2.1**: sharp edges — drawn as red `LineOverlay` segments (offset along
  the adjacent faces' average normal so they sit on top of the mesh), with
  `no_depth_test` available for X-ray inspection. The jump/orange overlay
  dies with the retired classifier. A face-shading mode for the jet prior
  (color by κ1/κ2) is a useful optional addition for calibrating the
  effort gates. (`EdgeHighlight` outlines whole faces rather than
  individual mesh edges, so it stays unpopulated — still TODO.md's delete
  candidate.)
- **2.2**: smooth regions — one `FaceHighlight` per region with a
  deterministic color palette, sharp-boundary `LineOverlay` retained for
  context.
- **2.3**: surface discovery — level 1 pauses at each committed hypothesis
  (claimed faces highlighted + fit overlay/parameters via the per-type
  `viz_overlays` / `viz_fit_description` functions, preserved from
  `bfs.rs` before its deletion); level 2 pauses at each competing rung's
  expansion result; level 3 at promotion attempts.
- **2.4**: resolution & repair — level 1 pauses after convergence; level 2
  at each resolution/merge/relabel operation (multi-claimed junction rows
  are worth a distinct highlight color).
- **3.3** (edge curves), **3.4** (faces), **3.5** (shells), **3.6** (solids)
  — unchanged.

### Architecture

The visualization uses a two-thread architecture with synchronous channel communication:

1. **Main thread**: Owns the window, drives the three-d render loop, and displays overlays.
2. **Pipeline thread**: Runs the brepper pipeline (stages 1–4) in the background.

- `VizSender::show_and_wait(overlay)` sends a `VizOverlay` to the render thread and blocks.
- The render thread displays face highlights, edge highlights, translucent geometry overlays (cylinders, cones, spheres, tori), status text, and optionally moves the camera.
- The user presses **Space** (NextStep), **Shift+Space** (skip to next region/operation), or **Q** (Quit); clicking a face skips ahead to that face's region.

### Key types (in `src/viz.rs`)

- `VizSender` / `VizReceiver`: Channel endpoints. Created by `viz_channel()`.
- `VizOverlay`: A display frame containing face highlights, edge highlights, cylinder/cone/sphere/torus overlays, line overlays, shape meshes, status text, and optional camera focus.
- `VizAction`: The user's response — `NextStep`, `NextSeed`, `SkipToSeed(face_index)`, or `Quit`.
- `FaceHighlight`, `EdgeHighlight`, `CylinderOverlay`, `ConeOverlay`, `SphereOverlay`, `TorusOverlay`, `LineOverlay`, `ShapeMeshOverlay`: Overlay primitives.

The render window (powered by three-d with egui overlay) provides camera controls (orbit/pan/zoom-at-cursor), face-number labels with occlusion (`N`), hover info with ray-cast picking, click-to-skip, and a help overlay (`?`).

Debug verbosity (`--verbose`) uses consistent indentation:
- **Level 1**: Stage summaries (region counts, classification counts, timing).
- **Level 2**: Per-region details (classification result, fit parameters, repair operations).
- **Level 3**: Per-rung fit attempts with error values; repair-loop round summaries.
- **Level 4**: Per-face evidence (curvature descriptors, relabel decisions).

---

## Source Code Organization

Each stage is a function taking a config reference and consuming the
previous stage's output structure:

- main:
  - let config = parse flags, load --compare STEP file
  - out1 = stage1(config)?           # ConnectedMesh
  - out2 = stage2(config, out1)?     # Stage2Output: regions + selected surfaces
  - out3 = stage3(config, out2)?     # Stage3Output: B-Rep topology
  - stage4(config, out3)?            # writes STEP file

Stage 2 module layout (target):

- `src/stage2/mod.rs` — orchestration, `Stage2Output` (`Segmentation`,
  typed hypothesis vectors, `SelectedSurface`), shared helpers.
- `src/stage2/edges.rs` — 2.1 sharp-edge classification + the jet prior
  (`JetField`: κ1/κ2, principal directions, medial points).
- `src/stage2/regions.rs` — 2.2 union-find smooth-region extraction.
- `src/stage2/classify.rs` — 2.3(a) whole-region ladder,
  subsample/validate/refit driver, jet rung ordering.
- `src/stage2/compete.rs` — 2.3(b)–(c) scoped competition: seeds, claims,
  expansion, bent-tube promotion.
- `src/stage2/repair.rs` — 2.4 resolution/merge/relabel, invariants.
- `src/stage2/{planar,cylindrical,conical,spherical,toroidal}.rs` — fitting
  + validation for each surface type (seeding/BFS entry points deleted).
- `src/stage2/bspline.rs` — freeform rung (parameterize/resample/fit/
  validate; gathering deleted; island-dissolution guard).
- `src/stage2/bfs.rs` — deleted with the BFS pipeline.

Data-structure linking uses indices into vectors throughout (stable, easily
debuggable, serialization-friendly; no lifetime complexity).

---

## Testing Strategy

The main testing strategy is to process a set of example stl/step pairs, and use the --compare flag to ensure that the pipeline is working correctly for each.

1. **Unit Tests**: Individual components (readers, fitters, edge classifier, region extraction, repair operations)
2. **Integration Tests**: Full pipeline on known geometries
3. **Regression Tests**: Compare output STEP to reference
4. **Validation**: Round-trip test (STEP → mesh → STEP)

Stage-2 integration tests assert:
- Partition sanity per model at 2.2: **smooth-region** counts (sharp-edged
  models = final counts; tangent models = tangent-chain counts, e.g.
  capsule 1, nosecone 2 — measure with `dump_regions` when re-pinning,
  don't guess) plus the coverage/connectivity/determinism invariants
  (already in `tests/stage2_segmentation.rs`).
- Surface-discovery end states after 2.3/2.4 — the same expectations the
  old suite pinned (cube = 6 planar; block_with_hole = 6 planar +
  1 cylindrical; pipe = 2 planar + 2 cylindrical; ball_on_cylinder =
  1 planar + 1 cylindrical + 1 spherical; hemisphere = 1 planar +
  1 spherical; cone_15x20 = 1 planar + 1 conical; etc.), plus
  fitted-parameter agreement with the reference STEP (axis, radius,
  center) for analytic models.
- Competition/repair behavior on the models that motivated the design
  (nosecone decomposes to 4 surfaces; pipe_elbow's bend promotes to a
  torus; capsule junction rows resolve; junction singletons absorbed).

**Migration reality**: the suite is red during the rewrite (see Rewrite
Roadmap step 7). Re-enable per model; a model stays disabled with a comment
naming its blocker rather than being nerfed with loosened gates.

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

The partition-then-classify design should land large models well inside
these targets: stage 2 does O(E) edge classification, O(N α(N)) extraction,
and one fit-and-validate per region — no per-seed redundancy and no
accidental-fit trials (under BFS, full stage 2 on 3DBenchy's 196k faces ran
~15 minutes; the planar-only path that skips hypothesis generation ran ~3 s,
which is the cost floor the new architecture approaches plus per-region
fitting).

---

## Future Enhancements

- [ ] Support for OBJ, PLY input formats
- [ ] IGES export option
- [ ] Hole detection and filling
- [ ] Feature recognition (holes, pockets, bosses)
- [ ] Richer curvature estimators for 2.1: per-vertex principal curvatures (Rusinkiewicz) would give full curvature tensors instead of the per-edge samples, sharpening the two-sided descriptor on irregular triangulations.
- [ ] Quadric pre-classification: fit a 10-coefficient quadric via SVD as a fast O(N) pre-rung — its eigenvalue pattern classifies plane/sphere/cylinder/cone directly and could order or prune the 2.3 ladder for large regions.
- [ ] Adaptive tessellation handling: `κ_floor` and the curvature descriptors assume locally comparable facet sizes; meshes with strong adaptive refinement may need area-weighted or multi-ring curvature estimates.
- [ ] Repair mode for assumption violations beyond stage 1.2's current repairs (e.g. boundary crack stitching — see TODO.md).
