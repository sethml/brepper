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
  infinite-radius cylinder).
- **Growth is two phases, and only the second is validated.** *This is the
  as-built shape and it differs from the original design; the difference
  matters, so it is spelled out here rather than left to the roadmap notes.*
  1. **Conditioning patch** (`seed_patch` → `condition_patch`): an
     **unvalidated** depth-limited BFS neighbourhood of the seed —
     `MAX_PATCH_DEPTH` 3 hops, `SEED_PATCH_CAP` 32 faces — optionally
     filtered by *jet similarity* (same coarse curvature class; magnitude
     within 1.5×; within `PATCH_RADII` = 4 curvature radii of the seed).
     The initial fit is a least-squares fit to that patch. `condition_patch`
     then trims fit-outlier faces until the seed validates against its own
     patch's fit, bottoming out at a size floor.
  2. **Flood**: BFS from the seed accepting a face only when the current fit
     validates it (vertex tolerance + centroid tolerance + orientation
     coherence), refitting to a fixpoint.
- *Why not the original design.* §2.3(b) first specified that initial
  parameters come **analytically from the jet at the seed face** — r ≈ 1/κ1,
  sphere centre ≈ centroid − n̂/κ1, ruling/axis from the κ≈0 principal
  direction, torus tube centres from medial points. That was abandoned
  because the jet's *magnitude* is not trustworthy at the scale of a small
  fillet: adjacent facets of one exact r = 0.5 arc read κ1 = 0.023, 0.162 and
  2.2, i.e. radii of 43 mm, 6 mm and 0.45 mm for the same feature, so
  "r ≈ 1/κ1" would be wrong by 100×. Fitting a neighbourhood is robust to
  that — but it introduces a failure the original design could not have,
  since there is no patch to straddle when the parameters come from one
  face. See "Why stepped_block_toroidal is hard" above: on a mesh whose
  facet sizes vary 2000:1 the patch cannot stay on its own surface, and the
  jet filter that is supposed to keep it there classifies only 0.17% of that
  model's faces as flat.
- **Ordering**: jet-informed candidate ordering grows the clean core first,
  so refits run from well-conditioned support and junction rows arrive last,
  when the model is already determined.
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

### Rewriting stage 2.3(b): the blind patch is the root cause

*Design proposed 2026-07-28 from a session of measurements; not yet built. Every
claim below has a number behind it, and the dead ends are recorded so they are
not re-derived.*

*(2026-07-29: the missing bootstrap design this section calls for is now
specified — see **The Robustness Rework**, Phase 1 "earned patches". Steps
0/1/3/4 below remain landed and correct; the step-5 deletion happens as part
of that phase.)*

**Root cause.** `grow_from` seeds a hypothesis by fitting to an *unvalidated*
BFS neighbourhood (`seed_patch` → `condition_patch`, depth 3, up to 32 faces).
The patch never survives into the hypothesis — after round 0 it is replaced by
the flood output, so every committed face has passed `face_validates`. What
survives is its influence on the *initial fit parameters*, and that is enough to
produce both observed failures:

1. The seed fails its own contaminated fit → `flood` returns empty → the
   candidate is rejected. 323 such events on `stepped_block_toroidal` fine.
2. Or the contaminated fit is merely *plausible*, and the flood validates a set
   of faces against the wrong surface. Every member genuinely lies within
   tolerance of it.

So validation is working perfectly and is not protecting us: the membership test
is sound, the thing being tested against is wrong. That is also why "prune to
what validates" is nearly a no-op — there is nothing unvalidated to prune.

**What survived this session's measurements, and is what the rework should lean
on.** Validation by vertex distance always behaved. The whole-region DOF ladder
is exact when the region is right (sharp-edged models need no tuning at all).
The per-facet fallback is *geometrically exact* for single-facet surfaces — 14
of this model's 16 planes arrive that way and are correct. Largest-validated-area
selection works (all 18 tori exact, 34/37 cylinders). The junction locus is exact
and verified corpus-wide. Refitting on corrected support is what buys accuracy.

**What did not survive — four local per-face signals, all measured wrong. Do not
re-derive.** Relative *face size* is anti-correlated with junctions (across a
flat|curved junction the median size ratio is 1.26; *within* a plane it is p90
47.5, because Delaunay grading rings a plane with slivers finer than the fillet
facets it abuts). The strip *long edge* is not reliably the ruling (adjacent
strips' long directions differ by up to 76°). A ring-local **shape operator**
from neighbour normals is worse than the jet on 155:1 strips and rank-deficient
on a quad whose four neighbours sample only two directions. And *dimensionless*
curvature κ·extent does not separate either — the plane quad sits at 0.139
against fillet strips at 2.10, inside a continuum.

**Status: steps 0, 1, 3 and 4 are built; step 5 is measured NOT YET POSSIBLE (2026-07-28).**

- **2.3(d), the leftover-component ladder rescue** (`rescue_leftover_components`,
  classify.rs). A leftover connected component *is* a region competition failed
  to explain, so the whole-region DOF ladder runs on each one before freeform
  gets it. No new math was needed: the conical fitter already estimates its axis
  as the smallest eigenvector of the normal covariance, and given an
  uncontaminated face set recovers r = 0.99998 / 0.99999 / 0.50008 against true
  1.0 / 1.0 / 0.5. Competition was simply never handing it a clean face set.
  It can only *add* surfaces — every rung is vertex-tolerance validated, so a
  genuinely freeform component still falls through to 2.3(e).
- **Visited ≠ claimed** (compete.rs). Committing no longer retires the faces a
  hypothesis grew into; every region face seeds exactly once. Faces in the
  *interior* of the explained set are skipped, because a fit grown from the
  middle of a validated surface merely refits that surface — a pragmatic bound,
  not a derived one. Reseeding with no bound at all was measured at over ten
  minutes on `stepped_block_toroidal` fine and abandoned; the bounded form costs
  42 s → 103 s.

- **2.3(d2), junction-constrained growth** (`grow_with_axis` in compete.rs,
  `rescue_with_neighbour_axes` in classify.rs). Borrow the axis direction from
  an adjacent committed surface and grow a cylinder with it fixed, largest
  validated area first, iterating until nothing more grows. With the axis fixed
  the fit is a 3-DOF circle in the perpendicular plane — well conditioned even
  on a short arc — so growth starts from the seed facet alone and there is no
  blind patch to straddle. This is what closes the leftover blob that 2.3(d)
  cannot: a blob spanning several fillets validates no single ladder rung, but
  it is bounded by surfaces that are already right.

Together: that model went **138 → 75 surfaces** (planes 80 → 18 against a
reference 16, cylinders 34 → 38 against 37, **tori 21 → 18, exactly the
reference**, B-splines 3 → 1), with the other 13 probe models bit-identical.

Superseded numbers from the intermediate state: after steps 0–1 alone it was
**138 → 95 surfaces** (cylinders 34 → 36 against a
reference 37, tori 21 → 19 against 18, B-splines 3 → 1, per-facet planes
80 → 39), with the other 13 probe models bit-identical.

- **Step 4, the (d)/(d2) fixpoint.** The two rescues alternate until neither
  commits anything: (d2) splitting a blob can leave a remainder one ladder rung
  explains, and (d) committing a surface gives (d2) another axis to borrow.
  Terminating, since every round strictly shrinks the leftovers. **Measured
  inert on the corpus** — the two already converge in one round on every model
  — so it is a guarantee rather than a fix, kept because the cost is one extra
  pass over an unchanged set.

**Step 5 — deleting `seed_patch`/`condition_patch`/`PATCH_*` — is NOT a cleanup,
and is not yet possible. Measured.** Cutting the blind patch to the seed alone
(`SEED_PATCH_CAP` = 1) collapses `stepped_block_toroidal` fine to **24872
surfaces**, essentially one per facet. The blind patch is still the sole source
of initial parameters for the *free* fit path, which commits ~63 hypotheses
before either rescue runs; (d) and (d2) only clean up after it. Deleting it
requires the constrained path to take over the primary role — every surface
grown on an axis borrowed from a committed neighbour, bootstrapped from
somewhere other than a blind neighbourhood — which is a further design step, not
a removal. The stage-2 constant count therefore stays at 34.

**What the remaining 24-surface gap now is** (measured after steps 0–1, on
`stepped_block_toroidal` fine: 95 against a reference 71). Twenty-three of the
twenty-four surplus surfaces are **single-facet planes** — leftover slivers that
nothing explained — plus one extra torus, one spurious B-spline, one missing
cylinder. Two distinct kinds of sliver, and they want different fixes:

- *Absorbable ones.* Face 6758 sits 7.7e-5 mm from an adjacent committed
  `cyl r=1.0000` and 4.1e-5 from an adjacent single-facet plane — the first is
  1.5× vertex tolerance, the second *inside* it. These are boundary faces that
  stage 2.4's relabel/merge/snap ought to take, and mostly the repair pass does
  reduce 140 committed surfaces to 95.
- *Unfitted curved runs.* Faces 9857/9866 and their neighbours are a chain of
  per-facet planes whose normals rotate progressively — n = [0,−0.280,0.960],
  [0,−0.269,0.963], [0,−0.248,0.969], [0,−0.238,0.971], [0,−0.227,0.974],
  [0,−0.206,0.979] — with adjacent-facet distances of 4.7e-4 to 1.4e-3 mm, far
  above vertex tolerance. That is a fillet arc that was never fitted as a
  cylinder, sitting in per-facet pieces.

**Limitation of 2.3(d) found doing this, and it shapes step 3.** The rescue
fits *one* surface per leftover component. After steps 0–1 the region still
leaves 104 leftovers whose largest connected blob is **49 faces**, and that blob
spans several distinct fillets — so no single ladder rung validates it and the
whole blob falls through. Re-running competition inside it does not help either,
since visited ≠ claimed already seeded every one of those faces. What such a
blob has that competition lacks is *committed neighbours*: it is bounded by
surfaces that are already correct, which is exactly the junction constraint of
step 3. Expect step 3 to be what closes this, not more seeding.

**The design (remaining steps).**

- **Seeding stays universal.** Every unclaimed face is a candidate seed; the
  constraint machinery below is an addition, never a fallback.
- **Visited ≠ claimed.** Growing into a face is a hypothesis, not a conclusion.
  A seed that commits nothing costs one entry from the untried list and nothing
  else; the faces it touched stay available. This is what makes an unreliable
  bootstrap survivable — a 36-face surface gets 36 independent attempts — and it
  makes the plan's long-standing claim that seed misordering "costs only time"
  actually true, where today an early bad hypothesis consumes faces a later,
  better seed needed.
- **The commit gate is IDENTIFIABILITY, not residual.** Measured and important:
  the three spurious tori on `stepped_block_toroidal` fine fit their 7–8 faces
  to **0.52 / 0.72 / 0.87 mesh quanta**, against genuine surfaces at mean 0.39
  and max 0.56. They overlap — a 7-DOF torus on 7 faces is essentially exactly
  determined and fits perfectly by construction. No residual or coverage test
  can reject them, exactly as *Identifiability is not model selection* says. The
  gate must ask whether the support pins the parameters (`MIN_AXIAL_SWEEP_DEG`
  and the conditioning it stands for), and the quantum belongs where it did
  work — as the *unit* for judging a fit whose parameters are already pinned.
- **Expand under junction constraints where a neighbour is already committed.**
  This is what `locus.rs` unlocked. Measured on the three cylinders this model
  loses: each is a fillet-chain segment sitting between committed neighbours of
  *exactly its own tube radius* (cyl r=1.0 / cyl r=1.0; torus r=1.0 / cyl r=1.0;
  torus r=0.5 / cyl r=0.5). Fixing only the axis direction recovers
  **r = 0.99998, 0.99999, 0.50008** from arcs of 65°, 60° and 40°, with maximum
  residuals of 1.0e-5, 1.0e-5 and 5.0e-6 mm — five to ten times *under* vertex
  tolerance. The free 5-DOF fit fails on the **axis**, not the radius, and the
  neighbour supplies exactly that.
- **Iterate seed passes to a fixpoint**, so a segment that cannot be fitted
  freely on pass 1 becomes a nearly-constrained problem on pass 2 once its
  neighbour is committed. Terminating: each pass either commits a face or ends.

**What it deletes.** `seed_patch`, `condition_patch`, `grow_patch`,
`SEED_PATCH_CAP`, `MAX_PATCH_DEPTH`, `PATCH_RADII`, `PATCH_TRIM_ROUNDS`,
`MIN_CONDITION_FACES`, and the jet's role in seeding (leaving it to the ladder's
rung ordering). The flood is sound and stays.

**The open risk, stated narrowly.** Not "a model with no anchor" — seeding stays
universal, so with nothing committed the design degrades to today's behaviour,
which finds 34/37 cylinders and 18/18 tori. The real gap is the **bootstrap for
curved surfaces on coarse meshes**: a plane is exact from one facet, and on fine
meshes the jet is good (four of six sampled fillet strips read κ₁ = 2.0055,
2.0051, 2.0049, 2.0058 against a true 2.0), but the recorded κ₁ = 0.023 / 0.162 /
2.2 spread over one exact r = 0.5 arc was measured on stepped_block **medium**.
That is the one place the blind patch is doing irreplaceable work. Test any
replacement against stepped_block_toroidal medium and part_rounded_cube coarse
before deleting it.

**Why stepped_block_toroidal is hard: the seeding is scale-blind, and its mesh
is not scale-uniform.** This model is one tangent-connected region of 26122
faces (384 sharp edges partition nothing), so all 71 reference surfaces —
16 planes, 37 cylinders, 18 tori — must come out of competition. It is the model
every threshold in stage 2 sits on a cliff for, and the reason is not that its
surfaces are hard to tell apart. They are trivially distinguishable: the fillets
are r = 0.5 and 1 against planes of ~100 mm, a curvature contrast of 4000:1.

*The mesh is wildly non-uniform.* 89% of the faces are fillet facets covering
224 mm²; the planes are 2.5% of the faces covering 13783 mm². Mean facet size is
**0.52 mm on the fillets and 4.6 mm on the flats**, and per-plane it is worse —
a 3232 mm² plane is **12 faces** (16 mm each) and a 550 mm² plane is 9. Facet
*area* varies by ~2000:1 across the model.

*Every neighbourhood mechanism in stage 2 counts faces, not space.* `seed_patch`
grows a depth-3 BFS neighbourhood capped at `SEED_PATCH_CAP = 32` — from a
surface that often has fewer than 32 faces in total (72 of 75 connected flat
components are under 32). Measured: of 200 sampled flat seeds, **none** produced
a patch that stayed on the flat; the median patch is 50% fillet facets and the
worst is 96%. The plane fit is then computed from a patch that is half a
different surface, and `condition_patch` is left to repair it against a floor of
6 and 16 trim rounds.

`PATCH_RADII` exists precisely to bound the patch in space rather than hops —
but it is measured in curvature radii of the seed's own jet, and its doc notes
"flats (jet magnitude ~0) are unbounded, as they should be". So on exactly the
seeds that need it here, the spatial bound is switched off.

*It corrupts the curvature estimate too.* Only **7 of the 16 reference planes**
appear as connected flat components at all. A plane that is one or two quads
ringed by fillets has its jet computed from a neighbourhood dominated by fillet
vertices, so it reads curved before seeding even begins.

*But this is NOT what costs the model its surfaces — correction, 2026-07-28.*
The mechanism above is real and measured, and it is tempting to conclude the
planes are therefore lost. **They are not.** `facet_plane` turns every
unexplained face into its own plane, and for a single-quad plane that fallback
is geometrically *exact* — one quad, one plane, correct normal and offset. The
committed decomposition contains **exactly 16 planar surfaces of meaningful
area totalling 14301 mm²**, matching the reference's 16 planes: fourteen
single-quad "steps" arriving via the fallback, and the two sliver-heavy "sides"
(195 and 193 faces) found analytically by competition. A seed that returns NO
GROWTH on a single-quad plane is harmless, because a per-facet plane and a
committed plane are the same surface there.

The real losses on this model are entirely in the **curved** surfaces: 34
cylinders against a reference 37, 21 tori against 18, 3 spurious B-splines, and
~62 junk single-facet planes of 0.04–0.11 mm² that should belong to neighbouring
fillets. Read the scale-blindness as a mechanism that *will* bite — it is why
the conditioning patch cannot be trusted — not as the explanation for this
model's surface count.

*What follows.* The constants are not badly-chosen values; they are face counts
standing in for lengths, on a mesh where face size varies by 30× across a
junction. That is why their responses are non-monotonic and why loosening or
tightening any of them "changes which accident wins" rather than converging.
**The lead is to bound and filter the seed patch by mesh scale**: a neighbour
whose facet is many times smaller than the seed's is, on a B-rep tessellation,
almost certainly on the other side of a face boundary — the density change the
*Mesh provenance* section lists as corroborating evidence, which here is a 9–30×
signal rather than a subtle one. Provenance-dependent, so it must stay a seeding
prior with validation still deciding every face; but as a prior it is far
stronger on this model than the jet it would supplement.

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

## The Robustness Rework (2026-07-29)

**This section is the rationale behind the phased plan in TODO.md.** TODO.md
carries the step-by-step instructions; this section records why each phase is
shaped the way it is, the design detail for the two new mechanisms (earned
patches, distance-guided exploration), the corpus strategy, what was considered
and deliberately not adopted, and the measurements preserved from the TODO
items the plan retired. It synthesizes two parallel design sessions held
2026-07-29 (one logged at `agent-logs/2026-07-29_14-59_nothinking.md`) plus the
measured record above.

### Goal and thesis

The goal changed: **maximum robustness and minimum fragility over all plausible
CAD-derived meshes**, not score on the current corpus. Fragility is
operationalized as two headline numbers: the stage-2 tuned-constant count
(34 at plan time; `grep -h "^const \|^pub(crate) const \|^pub(super) const "
src/stage2/*.rs | wc -l`) and the **cliff count** — the number of constants
whose ±2× perturbation collapses any model (four unrelated knobs once each
flipped stepped_block medium from 83 surfaces to full per-facet collapse).

The measured fragility concentrates in three structural debts, and the plan
retires them in dependency order:

1. **Evaluation is anecdotal.** 14 hand-watched models, Onshape/Parasolid-heavy;
   three NURBS models rotted silently because nothing watched them. You cannot
   claim generality you cannot measure → Phase 0 builds the instrument first.
2. **Two of three tolerances are guesses.** `--surface-tolerance` defaults to
   0.4 mm against a measured true sag of ~9e-4 mm (~400× too loose), making
   centroid validation nearly vacuous; `--angular-tolerance` is a stand-in for
   a discarded exporter setting. Both are measurable from the mesh itself
   → Phases 0 (diagnostics) and 3 (regime flip).
3. **One component — blind-patch seeding — forces most of the guard pile.**
   `seed_patch`/`condition_patch` fit an unvalidated 32-face BFS neighbourhood;
   the contamination lives on in the fit parameters, producing both the 323
   NO-GROWTH seeds and plausible-but-wrong fits. Roughly two-thirds of stage
   2's constants exist because tangent seams are invisible to it → Phases 1–2
   replace the mechanism, Phase 4 deletes what that makes inert.

Everything else measured principled — stage 1, sharp-edge regions, the DOF
ladder, `fit_validate`, `face_validates`, flood, `resolve_claims`,
relabel+refit, the locus machinery, stage 3 topological assembly, the
`--compare` gates — **is deliberately not touched**.

### Phase 1 design: earned patches (replaces the blind patch)

Replace grab-then-repair with **grow-under-validation**: BFS from the seed,
admitting a frontier face only if it validates against the current fit at
vertex tolerance; when a face fails, attempt **promotion** up the ladder
(plane → cylinder/cone or sphere → torus), refitting accepted+candidate and
requiring *every already-accepted face* to still validate; if no type
validates, the face is across a seam — exclude it and record it as a
**frontier face** (Phase 2 consumes these). The patch can never straddle a
junction by construction, which is the exact property the blind patch lacks.
Commit semantics do not change: the earned (type, params, patch) goes to the
unchanged `grow_from` flood/refit path, and every existing commit gate
(`MIN_AXIAL_SWEEP_DEG`, torus ratio, `MIN_*_FACES`, `resolve_claims`) keeps
its meaning.

What makes promotion workable at tiny support — where a free least-squares fit
is hopeless (sagitta arithmetic: over a 7.5° arc on r = 0.5 the sag is
~0.001 mm against a 2e-5 mm quantum, so a free radius is ~2% ≈ 200× vtol
wrong) — is that **facet normals are exact** and give closed-form,
DOF-reducing initializers:

- **Cylinder — axis = n₁ × n₂.** A planar facet spanning two rulings of a
  cylinder has its normal exactly perpendicular to the axis, at any coarseness,
  so two facets with non-parallel normals give the axis exactly on
  strip/quad-tessellated surfaces. Radius+origin then come from
  `fit_circle_for_axis` (cylindrical.rs), measured well-conditioned on short
  arcs (r = 0.99998 / 0.99999 / 0.50008 from 65°/60°/40° arcs). **Provenance
  caveat:** a general Delaunay triangle's normal is *not* exactly
  perpendicular to the axis; measure axis error per provenance (Phase 0's OCCT
  meshes), and degrade to the normal-covariance axis over the earned patch
  where needed — still clean, still not blind.
- **Sphere — concurrent normal lines.** For any planar facet inscribed in a
  sphere, the line through the facet's **circumcenter** along its normal
  passes exactly through the sphere centre (equidistance). Two facets → centre
  at the near-intersection (common-perpendicular midpoint), radius = distance.
- **Cone** — the existing normal-covariance axis (`smallest_eigenvector_3x3`)
  over the earned patch, guarded by `axis_normals_consistent`, both existing.
- **Torus — meridian coplanarity, the "open case" closed.** Every torus
  normal line lies in the meridian plane through its point, which contains
  the axis — so each facet contributes one coplanarity constraint
  `(pᵢ − a) · (nᵢ × d) = 0` on the axis (point a, direction d; 4 DOF; the
  system is bilinear — solve by alternation or pairwise elimination). Four
  facets pin the axis; then (R, r) is a 2-D circle fit in (ρ, z) profile
  coordinates — `fit_circle_2d` already exists (toroidal.rs). **Caveats:**
  rank-deficient by construction on sphere data (any axis through the centre
  satisfies it) and consistent-but-degenerate on cylinder data, so ladder
  order and a conditioning floor must precede it; grid-aligned quads give
  near-exact meridian normals, Delaunay triangles do not. **Time-boxed**: unit
  test to ~1e-9 on synthetic tori before integrating; if conditioning does not
  pan out, skip it — direct torus seeding stays as today and Phase 1 is not
  blocked.

Two behaviours need more than validation failure:

- **Ambiguity → branch, don't choose.** Two same-height faces of a cylinder
  validate as a sphere centred on the axis (the hole→sphere trap at patch
  scale). Grow both branches; the false one should die within a face or two
  (the next face at a different height is off the sphere by ≈ Δ²/2R ≫ vtol).
  Instrument branch lifetime; persistent branches mean validation is less
  sharp than believed — stop and ask.
- **Bent tube → drift-triggered promotion.** A bend is locally a valid
  cylinder (sag of the bend over a 2 mm patch at R = 14 is ~3.6e-5 mm), so
  admission never fails and promotion never triggers on residual. The signal
  is **systematic fit drift**: along a straight tube successive refits are
  stationary; along a bend the fitted axis precesses with a consistent sign as
  support accumulates. Test drift as sign-consistency of successive axis
  deltas, *not* magnitude (a magnitude bar would be a new tuned constant —
  stop and ask if one seems needed). Success retires `promote_bend` and the
  `BEND_*` constants; failure keeps them, and that is fine.

**Deleted when this lands (after ablation evidence):** `condition_patch`
entirely, `seed_patch`'s jet filtering, `SEED_PATCH_CAP`, `MAX_PATCH_DEPTH`,
`PATCH_RADII`, `PATCH_TRIM_ROUNDS`, `MIN_CONDITION_FACES`, likely
`MAX_TORUS_STARTS` (multiple flood starts exist because the fit patch was
unreliable).

### Reference labeling, and what it settled (Phase 1.0, 2026-07-30)

`src/reference.rs` (`--reference-report`) labels every mesh face with the
reference face it came from and scores the mesh against that ground truth.
Step-by-step numbers are in TODO.md 1.0; what belongs here is why it exists,
why it is nearly free, and the four results that change later phases.

**It is a read-back, not a new computation.** Stage 1's `--compare` already runs
`BRepExtrema_DistShapeShape` from every mesh face centroid to a compound of the
reference's faces. That call already determines which reference face won and at
what (u,v); `min_distance_to_shape` was returning only the distance and dropping
both. `SupportOnShape2` hands the winner back, so ground-truth labeling of the
entire mesh costs one accessor call per face that was already being computed.
The alternative that was on the table — a standalone probe — would have cost far
more than that and measured the wrong object: it reads raw STL triangles, where
the pipeline's mesh has been repaired and quad-fused, which is why
`step_measure_tolerance` needed a `--quad` flag bolted on. **The general lesson
is that a measurement belongs where the computation already happens**, and the
same investigation found the cost of ignoring it: `step_measure_tolerance` still
measures against the reference *solid*, so it silently reports exactly 0 for any
deviation over OCCT's hard-coded 1 µm classifier tolerance (`aTolerance = 0.001`
in `BRepExtrema_DistShapeShape`'s `SolidTreatment`) and is correct below it —
exact where the answer does not matter, zero where it does. That is the same
"never measure a distance against a solid" defect fixed in the pipeline on
2026-07-29, in the one place the fix did not reach. It also refines that entry:
the old stage-1 centroid check was not vacuous on all convex geometry, only
above 1 µm, which is exactly why correct models looked right and nobody noticed.

**1. The closed-form axis error is `perpendicularity error / |n₁ × n₂|`.**
Verified to 0.1% on `pipe_elbow_10` Parasolid medium: a perpendicularity max of
0.0280° over a conditioning p50 of 0.0271 predicts 1.033°, and the measured
error is 1.0342°, identical across all 348 pairs. One relation explains every
cylinder outlier in the sweep, and it separates two causes that a raw error
number confounds — normals that are not perpendicular, and perfect normals
amplified by an ill-conditioned cross product.

**2. The Delaunay caveat did not reproduce, and the reason is architectural.**
The Phase 1 design warned that "a general Delaunay triangle's normal is not
exactly perpendicular to the axis" and told 1.4 to degrade that provenance to a
covariance axis. Measured, OCCT's re-meshes read perpendicularity p50 = p95 =
**0.0000°** and axis error **0.0000°** at every density on five of six models.
The reason is already recorded under "Mesh provenance and the seam assumption":
OCCT's Delaunay is *constrained* by the pre-discretized face boundaries, so a
cylinder's triangles still span two rulings and their normals are still exactly
perpendicular. The provenance fallback is not needed. The single genuine
non-perpendicular case corpus-wide is `l_bracket_simple` Fusion medium
(perpendicularity p50 0.2740°), a small hole carrying 96 facets over 75 mm².

**3. Adjacent facets are the worst possible pairing, and finer meshes are
worse.** `|n₁ × n₂|` for an adjacent pair is `sin(angular step)` — measured p50
0.0432 on stepped_block fine, which is sin(2.48°) and independently recovers
Onshape's 2.5° setting. So the amplification in (1) is ~23× on a fine mesh and
grows as the mesh refines. This inverts the natural assumption that a denser
mesh gives a better closed-form fit, and it is a design constraint on Phase 1.2
rather than a property of any exporter: **pair the seed with the angularly most
distant admitted facet**, where the factor is 1 and the axis error falls to the
perpendicularity error itself (under 0.032° everywhere in the corpus). A 90°
reading, which is the max on three rows, is the fully degenerate pair of two
facets on one ruling strip — no information rather than a wrong answer.

**4. The torus initializer should not be skipped.** Phase 1.1 time-boxes it on
the suspicion that the meridian-coplanarity system is too ill-conditioned to
solve. Its residual — the quantity the solve drives to zero — measures p50
**2.354e-5 mm** over 23 863 torus facets on stepped_block fine, which is the
mesh's own float32 quantum. The constraint is satisfied to the limit of what the
mesh can express.

**The seam assumption now has a direct measurement, and it holds.** Phase 0.5b
measured only half of Design Principle #2 — that mesh *vertices* lie on the
reference surfaces. The other half, that a *facet* never straddles a face
boundary, is what the earned patch depends on and had never been measured.
Labeling gives it for free: 51 of 157 483 facets across 79 meshes disagree about
their own label across their own interior, **0.032%**. The number is an estimate
rather than a bound in either direction — it over-reports where a neighbouring
face curves back within reach of a sample point, and under-reports a straddle
confined to a facet's rim, since the samples cover the middle of a facet — but
both error terms are bounded by something this small, so the conclusion survives
the imprecision even though the count does not. The only recurring cluster is
exactly 12 facets on `stepped_block_toroidal` under OCCT coarse/medium and
Fusion medium/fine, absent under Parasolid and under OCCT fine; consistent with
coarse facets spanning that model's known narrow flat regions, and not chased.

**Segmentation scoring closes the over-segmentation hole**, which this document
recorded as needing its own step. Two area-weighted, threshold-free numbers from
the same labels: `ref coverage` = `Σ_R max_S area(R,S) / Σ area` is 1.0 when no
reference face was split across several of our surfaces, and `surf purity` =
`Σ_S max_R area(S,R) / Σ area` is 1.0 when none of our surfaces spans several
reference faces. They separate the two failure directions that "over-segmented"
was being used for loosely. On the case that motivated it — `rand_filleted_box_0`,
9 surfaces against a reference 26, passing every existing gate at 6.7e-3
relative volume — they read **1.0000 and 0.2201**, correctly diagnosing merging
rather than splitting, where a correct model reads 1.0000/1.0000. This satisfies
the standing rule about metrics: asked what it returns for a grossly wrong
input, it returns something unmistakably different from a right one.

**Freeform reference faces, and the property that actually governs facet-normal
exactness.** The normal report initially covered only analytic reference faces,
which is 32–66% of a NURBS model's area — the rest was announced as unmeasured
by the `analytic area` figure but not scored. Extending it needed only the
footpoint (u,v) that `ParOnFaceS2` already reports from the same distance query,
plus `GeomLProp_SLProps` (which arrived in the bindings on 2026-07-30).
`SLProps` is used rather than a raw `D1U × D1V` because it answers
`IsNormalDefined`, and a B-spline face has poles and seams where the normal
genuinely does not exist; those facets are excluded rather than scored as an
error of some size. Over the 8-model NURBS corpus every freeform facet evaluated
and none landed on a pole.

The result is not the obvious one. Freeform facet normals are 2–3 orders of
magnitude less exact than analytic ones — `bspline_terrain` p50 0.9396°,
`hill_with_hole` 0.9014°, against 0.0000° for planes — but that comparison is
confounded, because the corpus's freeform models are ~100× coarser per facet.
Ruledness separates them better: **`bspline_wave`'s *extrusion* faces read p50
0.0507° at 16.7 mm² per facet, while `bspline_terrain`'s B-spline faces read
0.9396° at 0.53 mm² per facet — the 30× coarser surface is 19× more exact.**

**But ruledness is necessary, not sufficient — the real condition is that the
tessellation be ruling-ALIGNED.** (Raised by the user, then measured; an earlier
version of this section claimed ruledness governs, which is wrong.) The exact
statement: a facet on a cylinder has its normal perpendicular to the axis **iff
two of its vertices share a ruling**. Parametrize the cylinder and project the
vertices along the axis; they land on a circle, the facet's plane is parallel to
the axis iff those projections are collinear, and three points on a circle are
collinear only when two coincide. A tessellator that lays triangles *across* the
ruling family therefore produces skewed normals on a geometrically perfect
cylinder.

Measured over the 12 cylinder-bearing models, 12 275 cylinder facets: facets
with an edge on a ruling are perpendicular to a worst case of **0.0006°**, while
the **25.6%** that straddle three or more rulings reach **1.0528°**. Every
non-zero perpendicularity in the corpus comes from a straddling facet, so
alignment is the whole story — though straddling is necessary rather than
sufficient for a *large* error, since the tilt grows with a facet's angular
spread and shrinks with its axial one. That also explains the extrusion number
above: 0.0507° rather than 0.0000° is itself evidence that those facets are not
aligned, because alignment would force exactness.

Where alignment comes from: quad fusion is a **23× discriminator but not a
guarantee** — 2.7% of fused quads straddle against 62.3% of triangles (a
coplanar quad on a cylinder *usually* spans two rulings, but four points on an
elliptical section are coplanar too). Simple geometry gets alignment for free —
a plain cylinder under OCCT is 120 quads / 0 triangles, 100% aligned, which is
the real reason result (2)'s OCCT re-meshes read 0.0000°, sharper than "the
rulings survive as face boundaries". Straddling concentrates where a tessellator
cannot lay a clean strip: `stepped_block_toroidal` straddles 49–53% at every
OCCT density, and its narrow flat regions between fillets are exactly its known
pathology.

This is *better* news for Phase 1.2 than the aggregate suggested, because the
test needs no ground truth. Given a candidate axis, the minimum angular gap
between a facet's vertices says whether that facet's normal can be trusted, so
the earn loop can **select** the facets it initializes from rather than hoping.
With the amplification law, the rule becomes: pair two **ruling-aligned** facets
that are **angularly distant**. For Phase 5's freeform work the consequence is
unchanged — facet normals cannot be assumed reliable on a doubly-curved patch at
any density, since there are no rulings to align to.

**A numerical fact worth not rediscovering.** An angle assertion near 0° cannot
be tightened past about 1e-6°, because it goes through `acos`, which near an
argument of 1 turns a 1e-16 rounding error into ~1.5e-8 rad. Assert on
`|â × b̂|` instead, which is linear in the error and resolves to machine epsilon.
Measured, not theorized: an exactly-constructed strip cylinder failed a 1e-9°
perpendicularity assertion at one sample out of sixteen, and the same property
passes at 1e-15 in cross-product form. The reported angles in `reference.rs` are
printed to four decimals of a degree, far above that floor, so the report is
unaffected — but do not read exactness below ~1e-6° out of a number from `acos`.

### Phase 2 design: distance-guided exploration

The 26122-face seed sweep costs 103 s on stepped_block fine and its
useful-work fraction is unknown (profile first — do not build this phase's
machinery if the profile says the cost is in fits or floods rather than the
sweep). The principled orderings and bounds, all derived rather than tuned:

- **Distance field**: per-face geodesic distance to the region boundary
  (multi-source Dijkstra over the dual graph, **centroid-gap weights in
  millimetres, not hops** — hop counts lie on graded meshes, which is what
  `PATCH_RADII` existed to patch). Seed **farthest-first** (Gonzalez maximin),
  replacing `neighbor_jet_mismatch` — the jet-based proxy for exactly this
  quantity.
- **Feedback, not detection**: tangent junctions cannot be *detected* (the
  measured-dead pivot; a tangent junction is locally the *smoothest* place on
  the mesh — an anti-edge), but they are *discovered exactly* by Phase 1's
  excluded frontier faces. Add discovered frontiers as new distance sources
  incrementally. The only trustworthy edge map is the one the fits produce.
- **Duplicate abort**: a seed whose reachable ball lies entirely inside one
  committed surface can only refit that surface — skip it. Topological, a
  strict generalization of `interior_of_claims`, no numeric tolerance.
- **Acceptance discipline**: ordering must not change results. Bit-identical
  probe output, or every diff explained; order-dependence is itself the
  fragility being removed — stop and surface it.

### Phase 0 design: the corpus and the instrument

Three provenance sources, cheapest first, all scored against reference STEPs:

1. **OCCT re-mesh of existing references** (automated, free ground truth).
   `gen_test_models` already drives `BRepMesh` with deflection/angle settings
   (`write_tessellated_stl`); add a mode that loads any reference STEP and
   emits STL variants at several settings, plus **off-origin variants**
   (translate the shape ~1 m before both STL and STEP output) so the float32
   quantum question becomes testable. This adds constrained-Delaunay
   provenance to a strip-heavy corpus.
2. **Fusion 360 re-export** (user batch job): import reference STEPs, export
   STL at low/medium/high. Different tessellator, same ground truth.
3. **Wild STEP+STL pairs** (most representative, most care): harvest models
   published with both a STEP/F3D and an STL (Thingiverse API / Printables /
   GrabCAD, scripted or via a browser session). **Pair-verification gate**:
   run stage 1 `--compare`; a pair counts as ground truth only if mesh
   vertices lie on the STEP surfaces to within a few quanta — otherwise
   quarantine as crash-test-only (the STL may be a different revision,
   decimated, or repaired). Record the verification *rate*: it is itself the
   first field measurement of the provenance premise ("vertices lie on the
   true surface") in the wild. Keep licensing clean: store source URL,
   author, license, and hash in a manifest; do not commit files whose license
   does not permit redistribution (gitignore the external corpus, commit the
   manifest).

The instrument: a one-command sweep harness (per-model rows: per-type surface
breakdown, `--compare` distances, volume error, runtime, exit status;
machine-readable; saved baselines + diff mode) and a **fragility meter** that
perturbs/ablates each stage-2 constant and reports which models change.

### Phase 0 as built (2026-07-29)

Steps 0.1, 0.2, 0.3 and 0.6 landed; step-by-step results with numbers are
annotated in TODO.md. What belongs here is the reasoning and the findings that
outlive the steps.

**The corpus was measuring one tessellator and one distance from the origin.**
`tests/remesh/` re-tessellates existing reference STEPs with OCCT's BRepMesh
at three densities plus a ~1 m translation, so ground truth is unchanged and
only provenance varies. That immediately produced the phase's headline result:
`stepped_block_toroidal` fails at **every** OCCT density — 183 surfaces coarse,
88 medium, 97 fine at 35596 triangles — while the Parasolid *fine* mesh of the
same body passes at exactly 71. Fine yields *more* surfaces than medium, so it
is not a resolution problem, and all 18 tori are recovered every time: what
fragments is the planes and cylinders around them. The strip assumption is no
longer a suspicion about the future; it is a measurement, and Phase 1's
initializer caveat ("a general Delaunay triangle's normal is not exactly
perpendicular to the axis") now has a model that demonstrates the cost.

**The seeded random population found what fourteen hand-watched models could
not.** `tests/random/` draws parameters where no constant was ever tuned, and
three of its five families fail on the first sweep:

- Filleted boxes collapse to **0 planes / 0 cylinders / 8 spheres / 1
  B-spline** (reference: 6 + 12 + 8). Their fillet-to-face ratios are
  0.03–0.10; `rounded_cube_10_r2`, the corpus's only filleted box and a model
  that passes exactly, is 0.20. The corpus had been exercising a ratio 2–7×
  larger than a typical part's, and at typical ratios one freeform patch
  swallows the entire flat structure. This is the single most useful new target
  for Phase 1: it is the blind patch's failure mode at scale.
- A torus with r/R = 0.84 becomes **40 surfaces (17 cones + 20 spheres + 2 tori
  + 1 B-spline)** for a one-face reference, while 0.36 and 0.10 are exact.
  `MAX_TORUS_RADIUS_RATIO = 0.7` is rejecting a legitimate solid, not just
  spurious fits. Phase 4.2 lists that guard under "expect these to still fire;
  firing means keep, not tune" — this row is the counter-case to weigh.
- A ball only slightly wider than its post (13.98 on 12.30) splits the sphere
  into three and adds a plane.

**The compare gates cannot see over-segmentation, so per-type counts are not
optional.** `rand_filleted_box_0` passes stage 4.1 at 6.7e-3 relative volume
with 9 surfaces against 26 and *no planes at all*; `rand_cone_chain_1` passes
at 2.0e-3 with a cone replaced by a B-spline. Both report a stage 3.1 edge
distance of exactly 0: that compare is one-sided — each of our edges to its
nearest reference edge — so a decomposition with *fewer* edges satisfies it
trivially, and the 1% volume gate is far too loose to notice. The measurement
protocol's insistence on per-type breakdowns was already recorded (a 3-torus
regression once hid behind an unchanged total of 71); this is a second,
sharper reason. A two-sided edge compare would be a real improvement to the
instrument, but it changes what "pass" means on every model and so belongs in
its own step, not smuggled into Phase 0.

**Pairing must be recorded, not inferred.** A re-mesh's reference lives in
another directory, and searching `tests/` for a matching stem is both ambiguous
(`plate_with_hole_100x50.step` exists under `fusion/` and `onshape/`) and
filesystem-order dependent — the first version of the harness silently scored
the OCCT re-mesh of the Onshape plate against the Fusion reference. The
generator now writes `tests/remesh/sources.tsv` and the harness reads it.
Any future generated corpus should do the same.

**"The strip assumption" was the wrong diagnosis; narrow flat regions are the
right one.** Phase 0.2 read the OCCT re-mesh result — `stepped_block_toroidal`
failing at every Delaunay density while the Parasolid mesh passed — as the price
of the strip assumption. Phase 0.5a's Fusion re-export tested that and it does
not hold. All three tessellators fail the model (best is Parasolid fine at 75
surfaces against a reference 71), so it is not Delaunay-specific; and triangle
count does not explain it either, because OCCT fine has 7.8× the triangles of
Parasolid medium and does *worse* (97 against 84). What the eight measurements
share is where the surplus lives: **planes, from 18 at best to 232 at worst,
while all 18 tori survive every provenance above ~3.7K triangles.** The failure
is narrow flat regions between fillets — precisely what a blind 32-face BFS patch
straddles, and precisely what Phase 1's earned patch is built to fix. Attributing
it to Delaunay would have sent Phase 1 after the wrong property.

Two smaller results from the same experiment. `cone_tangent_planes` decomposes
the same way under Fusion as under Parasolid (4 planes + 2 cones + 1 sphere),
which makes OCCT's 4 planes + 1 cone + 2 B-splines the outlier — so Phase 5c is
looking at a provenance-specific freeform reading, not at the model's nature.
And `nosecone`, which passes on the Onshape mesh at 3 surfaces, fails on the
Fusion mesh at 10 through two separate defects on one flat circular base: about
3000 mm² of it is absorbed into the adjacent cone, and what remains fragments
into six plane surfaces carrying byte-identical parameters. Six coplanar regions
that never merged is a connectivity failure, not a fitting one — the same family
as the "absorbing bridge facets" limit in LEARNINGS.md, but here it is fatal.

**Never measure a distance against a solid — three gates read zero for every
input.** `BRepExtrema_DistShapeShape` returns the **minimum** distance between
two shapes and reports **zero whenever one lies inside the other**. Point against
a solid therefore answers "is it inside?"; shape against shape answers "do they
overlap?". Neither is an accuracy measure, and neither announces itself, because
a wrong answer and a perfect one both print `0.000000e0`. All three of the
pipeline's reference comparisons had it, found together on 2026-07-29 and fixed
in `src/compare.rs`:

- **Stage 1** measured mesh vertices and face centroids against the reference
  solid. It passed a mesh written 1000× too small — a 25 mm sphere spanning
  0.025 mm, wholly inside its reference, distance exactly 0 — and its centroid
  check was vacuous on all convex geometry, since an inscribed tessellation puts
  every facet centroid inside the solid, which is the sag the check exists to
  measure. **One of the two sag truth values this plan quotes was 4.3× low as a
  result**: stepped_block fine measures 3.81e-3 mm, not 8.9e-4. rounded_cube
  fine is unchanged at 9.5e-4, so the old numbers were not uniformly wrong,
  which is exactly what made this hard to notice. Phase 0.7 calibrates its sag
  estimator against these — use the new ones.
- **Stage 3.4** sampled one point per face against the reference solid.
- **Stages 3.6 and 4.1** gated on the distance between our solid and the
  reference solid. Two versions of the same part always overlap, so those gates
  never fired at all. Measured on a decomposition with 9 surfaces against a
  reference 26 and *no planes whatsoever*: both printed the same 0 as a perfect
  model, and only the volume check beside them did any work.

The fix is not simply "use face compounds". For a point that works — a compound
of faces has no inside. But the minimum distance between two *boundaries* that
touch anywhere is still zero, so solid-vs-solid agreement needs a two-sided
quantity: `boundary_deviation` uses `BRepExtrema_ShapeProximity` with an infinite
tolerance, whose proximity value is the thinnest layer containing both shapes.
That now separates the cases it could not before — 1.97e-2 mm for a correct
model against 1.28e-1 mm for the 9-surface one. Its resolution floor is the
triangulation deflection it meshes at, a tenth of the tolerance being gated, so
it sharpens automatically if Phase 3 tightens `--surface-tolerance`.

**The general rule, worth applying to any new check: ask what your metric
returns for a grossly wrong input. If that is indistinguishable from what it
returns for a right one, the metric is wrong.** Both failures here — these gates
and the generator self-check that watched triangle counts while every file it
wrote was unusable — are the same error: verifying a proxy, or verifying nothing,
and reading the silence as success.

**A corpus gate has to be run before the corpus is trusted, not after.** The
first pair-verification pass reported 15 of 15 verified on a Fusion corpus that
was written in metres — every file 1000× too small and useless — because it was
an ad-hoc loop wrapping the stage-1 check described above. `tools/verify_pairs.sh`
is that gate, versioned, and it is now a thin wrapper: stage 1 reports the true
12.4875 mm and fails, so the gate needs no cleverness of its own. It briefly
carried a heuristic for exactly-zero distances, which is deleted — a workaround
for a bug is not worth keeping once the bug is fixed. It does still report the
f32 position bound beside each result, because that reveals coordinate scale for
free: 5.7e-10 mm for the metres sphere against 7.2e-7 mm for a correct mesh of
the same model. This matters well beyond one bad export: 0.5b's wild pairs are
*expected* to include decimated, repaired and rescaled meshes.

The generator side needs the same discipline. Fusion's `STLExportOptions`
inherits ambient state — `unitType` from the Design's default units, which for an
imported STEP was metres, and `printUtility` from the last UI setting — so an
export that does not set *every* property explicitly is not reproducible. That
single fact explains two separate anomalies that cost four runs to untangle:
byte-identical meshes across three refinement presets, and a silent 1000× scale
error. The lesson generalizes past Fusion: when a generator's defaults come from
the environment, pin them all, and verify the artifact rather than a proxy for it.
The batch's own self-check confirmed that triangle counts increased while every
file it wrote was unusable.

**Re-imported STEPs survive as ground truth.** Fusion can heal or reparameterize
geometry on STEP import, which would quietly destroy the pairing. Measured over
15 exports: vertices lie on the *original* reference surfaces at 4.1e-7 to
7.9e-6 mm, so the verification rate for that source is 15/15. That is the first
field application of the 0.5c gate, and it says a commercial CAD round-trip is
safe to treat as ground truth — which matters for 0.5b, where wild pairs are
expected to fail the same gate for unrelated reasons (different revision,
decimated, repaired).

**The provenance premise survives contact with the wild — but a third of
published pairs are not usable ground truth.** Of 18 STEP+STL pairs harvested
from open-source hardware repos (0.5b), 11 correspond and 7 verify outright; 5
describe different objects and 2 have meshes too broken to validate. Where they
correspond the premise holds hard: the *median* vertex sits on the CAD surface at
the float32 quantum, and 4 of the 7 that pass do so at or below one quantum, on
files we did not export.

*The first pass of this reported 6 of 11, and both numbers were wrong.* It
deleted 7 pairs as "different objects" on the strength of `BRepBndLib::Add(shape,
box, useTriangulation = true)` — where the flag means "use a triangulation IF THE
FACE HAS ONE", and a shape straight from STEP has none, so OCCT returned bounds
enclosing the surface poles. On a 30 mm cube that box reads 34.8 mm on one axis.
Two compounding errors followed: a padded reference box compared against an exact
mesh box invented size differences, and the denominator shrank because the
failures had been deleted. One of the deleted pairs verifies at 4.5e-5 mm. **Mesh
a shape before you bound it, keep the URLs and hashes of anything you delete, and
never quote a rate over a corpus you pruned using the same measurement.** Four
durable lessons:

*A max is the wrong statistic for a wild mesh.* Both ToolChanger ducts are exact
to float32 at p95 (6.9e-6 mm against a 9.3e-6 bound) and fail only on a thin tail
reaching 4.3e-4. A single stray vertex — from repair, decimation or a
simplified small feature — is a different fact about a mesh than a systematic
offset, and the max cannot distinguish them. Stage 1's compare therefore reports
**p50 / p95 / max**; the shape of the distribution names the cause on sight.

*Rigid alignment cannot be more precise than the mesh's own sag — so below that
it must not act at all.* A tessellated mesh's bounding box is inset from the true
solid, and asymmetrically, so a transform recovered from bbox centres carries a
sag-scale error. Applying such a transform where none was needed makes a good
pair worse: two ToolChanger ducts went from a vertex p95 of 6.9e-6 mm to 4.4e-4
under a 4.4e-4 mm "correction". Alignment therefore ignores translations below a
noise floor, and the floor is a measured gap rather than a tuned constant — the
corpus's shifts are bimodal with 4000× of empty space between the noise (8e-8 to
1.6e-3 mm) and the real displacements (6.5 mm and up), so any cut inside it gives
the same answer. The skipped magnitude is printed, so a corpus that lands in the
gap announces itself. Above the floor the method recovers exact round-number
placements perfectly (`NucDeck`, 6.52 mm → 4.58e-7 mm) and nothing finer. Going
below would mean fitting the transform to the distance being tested, which is why
alignment stays discrete: bbox centres plus the first of the 24 proper rotations
of a cube whose unsorted extents match. It never sees the distance it helps to
pass, so it cannot make a wrong pair look right.

*Check the file's text precision before theorizing about geometry.* An ASCII STL
at 6 significant digits has ~1e-4 mm resolution on a 50 mm part, which is exactly
the scale of these residuals — a coincidence worth ruling out rather than
assuming. Here it was ruled out: the STEPs carry 15–20 significant digits and the
STLs are binary.

**The cliff count is 9 of 32, and four of the nine are the blind patch.** The
fragility meter (`tools/fragility.sh`) perturbs one constant at a time by ±2×
and reports which models change and which *collapse* — a passing model that
now fails, or a surface count that at least doubles. `MAX_PATCH_DEPTH`,
`MIN_CONDITION_FACES`, `PATCH_TRIM_ROUNDS` and `SEED_PATCH_CAP` are cliffs and
account for most of the collapses; all four are already on Phase 1.5's deletion
list. The thesis above — "one component, blind-patch seeding, forces most of the
guard pile" — was inferred from single-model debugging sessions; it is now a
measurement over 45 models, and Phase 1's success criterion can be stated as
*those four cliffs disappear without new ones appearing*.

Three things the meter found that were not part of the thesis:

- **`NORMAL_AXIS_TOL` (conical.rs, 0.05) is a cliff nobody had flagged**:
  halving it collapses all three `rand_cone_chain` models,
  `block_with_conical_hole` and the OCCT re-mesh of `cone_tangent_planes`. It
  appears on no deletion list and no watch list. Any future work on cone axes
  should treat it as both load-bearing and fragile.
- **`MIN_HALF_ANGLE_DEG` and `PROJ_GAP_RECHECK` are cliffs upward only**
  (2° → 4° kills the random cone chains, whose half-angles are genuinely small;
  0.02 → 0.04 collapses a cone chain and a filleted box). Cliffs are not
  symmetric, so a meter that only tightened or only loosened would have missed
  both.
- **16 of 32 constants are completely insensitive to ±2×.** This is *not* a
  deletion list, and reading it as one would be the same mistake as
  complexity-pricing: a guard that fires far from its own boundary is exactly
  what a correct guard looks like. `MIN_AXIAL_SWEEP_DEG` is measured essential
  (it is what stops a surface of revolution being fitted by support that never
  observes the revolution) *and* measured ±2×-insensitive. Phase 4.2's
  ablations must therefore disable, not perturb — the two experiments answer
  different questions, and the meter answers only "is this value
  load-bearing".

Cost and cadence, since the full matrix costs ~110 minutes: run it at phase
boundaries, not per commit, and use `--only` for the constants a change plausibly
touches. The current implementation is about 3× slower than necessary — it
parallelizes models within one perturbation's sweep, where wall time is pinned by
the single 101 s model at 39% core utilization, when it should parallelize across
perturbations. Phase 2.1's runtime work on that model removes the binding
constraint either way.

**The f32 "quantum" is the wrong quantity for a distance, as recorded — and
here is the right one.** The bound on how far f32 storage can move a vertex off
its true surface is `max over vertices of ‖(ulpₓ, ulp_y, ulp_z)‖ / 2`: each
coordinate's error is at most half an ulp, and the largest projection of that
box onto any unit direction is the norm of its half-widths (Cauchy–Schwarz).
Measured on the 1 m off-origin meshes: the bound is 4.32e-5 mm and the true max
vertex-to-surface distance is 3.93e-5 mm, so the bound is real and tight (0.91),
while `max|coordinate| × 1.19e-7` claims 6.10e-5 mm — 1.4× too loose, because
it ignores that a distance is a projection while also ignoring that all three
coordinates contribute. Both figures are printed in the stage-1 summary (the
naive one for continuity with earlier notes); only the bound is used as a
floor. **The clamp is inert across the whole corpus**, including at 1 m, and
first bites past ~1.02 m where the f32 exponent steps. It is a correctness
floor for a part positioned in a large assembly, not a lever on corpus scores.

### What a tessellator actually bounds (measured 2026-07-30)

`--surface-tolerance` and the exporter's chordal setting are **the same physical
quantity**, so the tolerance is not merely guessable — it is recoverable. Every
tessellator driven so far takes *two* bounds, not one: OCCT's
`IncrementalMesh(shape, deflection, …, angle, …)` and Fusion's `surfaceDeviation`
plus `normalDeviation`. Both subdivide until both bounds hold, so

    achieved sag = min( chordal bound , R · (1 − cos(θ/2)) )

and **which term wins depends on the surface radius**. Measured over
`tests/remesh` (OCCT) and `tests/fusion_reexport` (Fusion) at matched requested
deviations of 0.5 / 0.1 / 0.02 mm and 30° / 15° / 6°, as achieved-sag ÷
requested-deflection:

| model | curvature radius | coarse / medium / fine | binding |
|---|---|---|---|
| `sphere_25` | 12.5 mm | 1.33 / 1.60 / 1.48 | chordal |
| `pipe_elbow_10` | 5 mm tube | 0.27 / 0.29 / 0.27 | angular |
| `cylinder_10x30` | 5 mm | 0.09 / 0.11 / 0.09 | angular |
| `rounded_cube_10_r2` | 2 mm fillets | 0.08 / 0.12 / 0.11 | angular |
| `plate_with_hole` | 5 mm hole | 0.09 / 0.11 / 0.09 | angular |

Three things follow, all measured rather than argued:

- **The angular bound governs sag as θ², exactly.** `cylinder_10x30`'s achieved
  sags are 4.278e-2 / 1.071e-2 / 1.714e-3, ratios of **4.00×** and **6.25×**
  against θ² predictions of (30/15)² = 4.00 and (15/6)² = 6.25.
- **The crossover sits at R ≈ 12–15 mm** for these settings, from
  `R = d / (1 − cos(θ/2))`. `sphere_25` at R = 12.5 mm is the one model on the
  chordal side, and it is the one model that reaches its bound. Fusion shows the
  same pattern at the same crossover, so this is a property of tessellation, not
  of one vendor.
- **A nominal export setting is therefore the wrong number to adopt.** Taking
  `--surface-tolerance` from the exporter's deflection would be up to ~12× too
  loose on exactly the features that matter — a 2 mm fillet achieves 0.08 of the
  requested deflection. Worse, there is no single right value per file: within
  one model the achieved sag varies by an order of magnitude between a 2 mm
  fillet and a 100 mm face.

The chordal-limited rows reading 1.3–1.6× *over* their bound is not the
tessellator overshooting. They are all spheres, and the two quantities differ by
a factor that can be derived exactly: the deflection bounds the **edge** sagitta,
`R θ²/8`, while the gate measures the **facet centroid**. For a triangle whose
vertices lie on a sphere at angular radius `ρ ≈ θ/√3`, the centroid sits
`R(1 − cos ρ) ≈ R θ²/6` off the surface, so

    centroid deviation / edge sagitta = (Rθ²/6)/(Rθ²/8) = 4/3

and `sphere_25_occt_coarse` measures **0.667498 mm against 4/3 × 0.5 = 0.666667**,
a 0.1% match. Medium and fine run 11–20% above 4/3 because OCCT's triangles are
not perfectly equilateral and the reported figure is a max over facets. Divided
out, all of them sit just under the bound as they should.

(This 4/3 is *not* the same factor as the Phase 0.7 estimator's 0.50× on spheres.
That one is 2×, from both principal curvatures displacing a facet centre where
`hθ/8` accounts for one. Two different comparisons; an earlier version of this
note conflated them and quoted ~2× here.)

This is also the explanation for the four coarse Fusion models that stage 1
rejects (0.5a): they are the chordal-limited ones, where the achieved sag really
does reach the requested 0.5 mm and so exceeds the 0.4 mm default.

### Phase 3.2 dependent audit: what `--surface-tolerance` is actually doing

Every one of the 29 uses read on 2026-07-30, before flipping anything. The
headline is that **one constant is standing in for three different length
scales**, and per facet they are:

| quantity | value | on `rounded_cube` fine (R = 2 mm, θ = 2.5°) |
|---|---|---|
| sag | `R(1 − cos θ/2) ≈ Rθ²/8` | 4.8e-4 mm |
| facet width `h` | `Rθ` | 8.7e-2 mm |
| feature scale | radius, wall thickness | ~2 mm |

**`h / sag = 8/θ`** — 183× at θ = 2.5°, 15× at θ = 30°. The ratio between the two
things `stol` is being asked to bound is itself a function of the angular
tolerance, which is why no single constant can track both: the gap between them
swings 12× across the density range.

**1. Disambiguating coincidental vertex agreement** — `face_validates`'s centroid
half, reached from `fit_validate`, both `grow` paths and `merge_regions`. The
load-bearing use, and *not* a fit-quality measure: a flat end cap fan-triangulated
from a boundary vertex has every vertex on the boundary circle, which lies exactly
on the adjacent cone, so vertex-only validation absorbed the whole cap and
`simple_cone` collapsed to one surface. It separates **sag** (a facet that
belongs) from **feature scale** (one that merely coincides). Single threshold:
**no** — measured sag runs 1e-3 mm on fine meshes to 0.667 mm on `sphere_25`
coarse, and 0.667 already exceeds the 0.4 default, which is why four coarse Fusion
models are rejected in stage 1. Derivable exactly from θ and `h`, both per-facet
measurable.

**2. Boundary-placement slack at tangent junctions** — stage 3.1 substitutes
`stol` for `vtol` on edges and vertices flagged tangent, because the region
boundary can land one facet row off the junction. That displacement is **`h`, not
sag**, so this use wants a quantity 183× larger than use 1 on the same mesh.
Not derivable from θ alone (needs a length), but `h` is directly measurable.

**3. Acceptance bars for reference comparison** — stage 1's centroid gate,
stage 3.4 face samples, 3.6 and 4.1 boundary deviation, and the nearest-STEP-face
lookup. "Is our output as close as the mesh permits", where the ceiling on
achievable accuracy *is* the sag. A global value is more defensible here than
anywhere else — a test wants a stable bar — but it still has to scale with the
mesh, or coarse models fail spuriously. They do today.

**4. Reconstruction plausibility bounds** — vertex snapping capped at `2×stol`,
and edge-curve-versus-boundary-chain deviation in four places. Already standing in
for sag by admission; the comment reads "mesh boundary vertices may be slightly
off the analytical tangent curve due to tessellation sagitta. Use
surface_tolerance rather than vertex_tolerance since the sagitta of curved
surfaces can be significant." It cannot say *how much* sagitta, so it uses a
constant.

**5. Freeform fitting tolerances — LIKELY A LIVE BUG, not a tuning question.**
The §2.7 approximator fallback is given `stol/8` and `stol/2`, and the
spline-through-boundary-chain fit `max(stol/10, vtol×10)`. Those points are mesh
*vertices*, which lie on the true surface to ~1e-6 mm, so the fit should be held
to `vtol`. `stol/8` = 0.05 mm is about 5×10⁴ times looser than the data warrants.
Worth checking against Phase 5a's conditioning failures before assuming it is
harmless. Note the boundary-chain site already derives its tolerance with a `vtol`
floor, so the pattern exists in the tree — it is just applied inconsistently.

**6. Effectively inert** — the §2.7 pole-bound margin `0.5·diag + stol`, where the
bounding-box half-diagonal dominates by orders of magnitude. Consistent with the
fragility meter, which found no B-spline constant to be a cliff.

**Conclusion for 3.2.** The pipeline does not need one derived tolerance; it needs
**two per-face quantities plus the one it already has**: `sag(f)` for anything
about a facet's interior (uses 1, 3, 4), `h(f)` for anything about boundary
placement (use 2), and `vtol` for anything about vertices (use 5). All three are
measurable from the mesh with no fit involved — `estimate.rs` already computes the
first two per edge. From θ *alone* neither length is derivable, but `R = h/θ` is
measurable, so one measured length plus θ determines both scales.

### Tolerance sequencing (Phases 0 and 3)

Estimators are built and *reported* in Phase 0 (they double as the
provenance-violation detector: estimator disagreement or a missing dihedral
gap means a decimated/remeshed/faceted-by-construction input — warn now,
degrade gracefully in Phase 6). The **default flip is deferred to Phase 3**,
after the seeding rework stabilizes: a 400× change to `stol` moves behaviour
everywhere (`REFIT_GIVE_UP` is 8×stol; merge gates and `face_validates` use
it), and debugging it concurrently with new seeding would confound both. The
earn loop admits on *vertex* tolerance, so Phase 1 is essentially independent
of the stol regime. Ablation (Phase 4) runs after the flip so deletions are
judged under the final regime.

### Considered and not adopted (with reasons)

- **Recursive subdivision / split-and-merge as the architecture.** Strong
  asymptotics (the tree visits only boundary-containing paths) and its merge
  half already exists (`fit_validate` as union-revalidation), but: fixed
  geometric cuts do not respect seams (boundaries are mesh-edge chains;
  strip-shaped surfaces stay contaminated at every level until cells are
  smaller than the strip), and identifiability *inverts* — a sub-node of a
  real torus legitimately sweeps 5° and fits garbage parameters that validate
  locally, so every sweep guard would need to become "defer, don't reject"
  until merge time, a rework of commit semantics. Phases 1–2 deliver most of
  its benefit at a fraction of the risk. **Held in reserve**; do not start
  without Phase 1–4 measurements in hand and user sign-off.
- **Precomputed junction-evidence field for seeding** (dihedral
  second-difference, diffused; or any jet/dihedral edge map). Any field
  computed *before* fitting is a detector prior — the refuted detector
  returning in a hat. Superseded by the exact, threshold-free fit-derived
  frontier field of Phase 2. (For the record, the raw idea also had a real
  basis: within a tessellated surface the dihedral field is locally
  consistent and *steps* at junctions — but the OCCT-Delaunay and
  coarse-mesh cases erode it, and nothing needs it once frontiers exist.)
- **The "bridge" formulation of type promotion** (grow under an
  exactly-determined uncommittable fit until the new DOF identify). Correct
  shape, but its fragile middle regime — junction faces soaking into
  exactly-determined fits — is mostly obviated by the closed-form
  initializers, which pin the new DOF geometrically at 2–4 faces instead of
  waiting for least-squares identifiability. The earn loop keeps the bridge's
  two sound rules: full-support revalidation on every accretion, and
  identifiability-gated *commits*.
- **VSA (the old item C)** as a planned phase. Earned patches + the existing
  relabel/refit fixpoint already implement the alternation locally, and the
  known-working machinery it would replace (`resolve_claims`, relabel, merge)
  is exactly what a big rewrite puts at risk. Reserve, same gate as
  subdivision. The pending `snap_boundaries` keep/delete decision folds into
  Phase 4 with fresh ablation data (it is the per-face assignment criterion a
  future VSA would need, and the locus diagnostic depends on the module).
- **Retuning anything to close the last corpus gaps first** (stepped fine's 4
  surplus surfaces, medium's ~13). Deliberately deferred: they are boundary-
  and identifiability-residue that Phases 1–4 either fix or isolate, and
  chasing them now re-enters the one-model-at-a-time loop that produced the
  fragility.

### Preserved measurements from retired TODO items

(Everything here is measured; do not re-derive. Sources: 2026-07-26..29
sessions.)

- **Sagitta bound**: no local patch determines a free surface of revolution —
  7.5° arc on r = 0.5 → sag ~0.001 mm vs 2e-5 mm quantum → r ~2% (≈200×
  vtol) wrong. Only DOF reduction helps at small support.
- **Blind-patch facts**: `SEED_PATCH_CAP = 1` collapses stepped_block fine to
  24872 surfaces (the blind patch was the free path's only bootstrap);
  0 of 200 sampled flat seeds produced a patch that stayed on the flat
  (median 50% fillet contamination); 72 of 75 flat components are under the
  32-face cap; a 3232 mm² plane is 12 faces.
- **Borrowed axes**: per-seed borrowed-axis candidates fragment surfaces
  (75 → 100, cylinders 38 → 63, measured twice); (d2) works *only* because it
  selects largest validated area globally before committing. Any borrowed-axis
  use outside (d2) must keep that global selection.
- **Jet reliability**: fine meshes good (κ₁ = 2.005 vs true 2.0); coarse
  garbage (κ₁ = 0.023 / 0.162 / 2.2 across one exact r = 0.5 arc on medium).
  Facet long-edge direction is right ~65% on strips (up to 76° off) —
  corroborating only, never load-bearing.
- **Identifiability, not residual**: impostor cylinders/tori sweep 2.8–6.7°
  vs real surfaces ≥ 64.8°; spurious exactly-determined tori fit at
  0.52/0.72/0.87 quanta vs genuine mean 0.39 / max 0.56 — overlapping, so no
  residual or coverage test rejects them. `candidate_sweep` is vacuous for
  spheres and true cones — a normal-spread identifiability test for those is
  still missing (Phase 5c territory).
- **Exact-normal component exploration** (partial, superseded as a detection
  path): connected components over near-zero-dihedral edges give 342
  components ≥ 10 mm² where the reference has 16 planes, threshold-insensitive
  from 0.001° to 0.2°; top ~16 sum to 14301 mm² vs expected 13783. The
  unexplained tail was suspected large-radius fillets with sub-threshold
  steps. Do not resume as detection; frontier feedback replaces it.
- **hill_with_hole residue**: one junk torus — 17 faces (one over
  `SMALL_SELECTION_ABSORB_FACES = 16`), R = 3.5101 r = 2.3882 (r/R = 0.680,
  just under the 0.7 cap), area 1.533 mm² = **0.095%** of its freeform patch —
  vs `cone_tangent_planes`' real cones at large fractions (1081/163 mm²).
  Relative area is the one untried candidacy quantity. With the decomposition
  forced correct to 7 surfaces it then fails stage 4.1 volume at 6.96e-2
  relative — a separate freeform-accuracy/export issue.
- **Freeform fit failures**: `bspline_blister`'s 926-face component rejected
  by the §2.7 pole-bound check ("control points 2749.5 mm outside the data
  bounding box"), saddle's three components likewise; the check is right, the
  fit conditioning is the bug; nothing recovers, regions fall to per-facet.
  When they do decompose, 4 and 3 invalid micro-sliver faces remain (ignored
  test).
- **Sag truth values** for estimator validation: stage-1 `--compare` centroid
  max ≈ 8.9e-4 mm (stepped fine), 9.5e-4 mm (rounded_cube fine); Onshape
  angular settings 12.5°/6.25°/2.5°; Onshape-coarse worst sag 0.026 mm,
  Fusion-low 0.038 mm.
- **cone_15x20**: its reference STEP is itself faceted (59 planes, not
  analytic) — a corpus-data defect, not a pipeline bug; regenerate or mark
  unusable in Phase 0.

### Map from the retired TODO to the new plan

- Runtime bound (old item 2) → Phase 2 (profile first).
- stepped fine surplus / medium gap (old item 3) → Phase 1 acceptance targets,
  then Phase 4/5 residue.
- NURBS corpus blockers → Phase 5a (blister/saddle fit), Phase 5c
  (hill_with_hole selection).
- B-prime (quantum) → Phase 0 clamp+report, Phase 5b unit-of-quality.
- B (tolerance estimation) → Phase 0 diagnostics + Phase 3 flip.
- C (VSA) and recursive subdivision → Reserve.
- Step 5 (delete the blind patch) → Phase 1 (earned patches are the missing
  bootstrap design).
- snap_boundaries decision → Phase 4 with fresh ablation data.
- Faceted-by-construction / OpenSCAD → Phase 6.
- 3DBenchy, nosecone coverage, stage-1-only audit → Phase 0.
- Radius-refinement (old "robust rather than tuned" item 3) → Phase 4
  diagnostic (radius clouds should collapse with clean supports; if not,
  investigate support disagreement — never clustering, never snapping).

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

`tools/sweep.sh` is the test-results table: one machine-readable row per model
with the per-type breakdown, `--compare` distances, volume error, runtime and
exit status, plus baselines and diffs. It replaced `tests/RESULTS.md` and
`scripts/run_test_table.sh` (deleted 2026-07-29), which produced a
human-readable table with no per-type counts, no baseline to diff against, and
no way to tell an improvement from a regression.

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
