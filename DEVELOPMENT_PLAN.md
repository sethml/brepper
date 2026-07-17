# STL to STEP Conversion Utility - Development Plan

## Project Overview

**brepper** (B-Rep from Mesh) - A utility to convert triangulated STL meshes from CAD exports into parametric STEP files with fitted analytic and freeform surfaces.

## Design Principles

The stage 2 architecture was rewritten (2026-07-17) from a fit-driven BFS
hypothesis pipeline to **partition-then-classify** segmentation. These
principles drive the design; the old pipeline and the reasons it was retired
are recorded in [Abandoned Approaches](#abandoned-approaches).

1. **Exactness is the superpower.** Mesh vertices lie on the true CAD
   surfaces within numerical epsilon. A surface fitted to *all* faces of a
   region either validates every vertex within `--vertex-tolerance` or the
   region is not that surface — validation is decisive, not statistical. We
   are not solving the fit-primitives-in-noisy-scans problem the literature
   (RANSAC lineage) addresses.
2. **The answer's boundaries are already in the mesh.** CAD tessellators emit
   seam edges along B-rep face boundaries — *including tangent junctions*
   (verified 2026-07-17 on the nosecone cone–sphere tangency: 0 of 1572 STL
   triangles straddle the tangency ring; 104 vertices lie exactly on it). So
   the true segmentation is a partition of the mesh face graph whose region
   boundaries are a subset of mesh edges, and every region is connected.
3. **Partition first, fit later.** Segmentation uses fit-free signals
   (dihedral angles, discrete-curvature jumps). Fitting happens once per
   *maximal* region. This makes accidental exact fits — a band of sphere
   faces lying exactly on some cone or torus — structurally irrelevant: a
   sub-region surface can never explain the full region, so it is never even
   hypothesized. The abandoned pipeline spent most of its heuristics (stage
   ordering, seeding masks, coverage acceptance tests) and much of its
   runtime fighting exactly these fits.
4. **Repair by validation, locally.** Segmentation mistakes surface as
   classification failures and are repaired by local operations — split a
   region that fits nothing, merge adjacent regions that fit one surface,
   relabel boundary faces to the better-fitting side. No global fix-up
   passes with chordal slop.
5. **Invariants are enforced, not assumed.** Before stage 3: every mesh face
   belongs to exactly one region; every region is mesh-connected; every
   region's surface validates its own faces; results are deterministic.
   (The abandoned pipeline's greedy selection violated connectivity —
   3DBenchy face 5841 reached stage 3.4 claiming three disconnected mesh
   regions.)

## Assumptions

- Assumes the vertices of the mesh lie precisely on the CAD surfaces, modulo numerical error/stability.
- Assumes transitions between shapes always have a seam — verified to hold
  even for tangent transitions (see Design Principles #2).
- Assumes manifold mesh (each edge shared by exactly 2 faces).
- Assumes consistent face orientation (all normals point outward).

Fails noisily if any of the assumptions are violated, after stage 1.2's
repair pass has fixed what can be fixed without moving a vertex
(`--strict-mesh` disables repair).

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

**Status: this stage is a from-scratch redesign (2026-07-17). The code on
disk still implements the abandoned BFS pipeline until the
[Rewrite Roadmap](#rewrite-roadmap) lands.** There is no migration flag: the
new structure replaces the old outright, and the test corpus is re-enabled
model by model as classification matures.

**Pipeline:** classify mesh edges (2.1) → extract connected regions (2.2) →
classify each region against a ladder of surface types (2.3) → repair by
split/merge/relabel (2.4) → emit selected surfaces (2.5).

**Precision note**: Mesh vertex positions, normals, and all fitted-surface
parameters use double-precision (`f64`) throughout.

#### 2.1 Edge classification

Label every manifold mesh edge INTERIOR or BOUNDARY using fit-free local
signals. Region boundaries are a subset of mesh edges (Design Principles #2),
so this labeling *is* the segmentation; everything downstream refines it.

**Discrete curvature field:**
- For each manifold edge `e` between faces `f1, f2`: the signed dihedral
  angle `θ_e` (positive when the edge is convex — normals diverge; negative
  when concave), the centroid distance `d_e`, and the discrete normal
  curvature along the crossing direction `κ_e = θ_e / d_e`.
- Per-face curvature descriptor: `σ_max(f)` and `σ_min(f)` — the maximum and
  minimum signed `κ_e` over the face's manifold edges. These approximate the
  face's two principal curvatures, sampled along the available crossing
  directions (a fused quad samples both grid axes; a triangle samples three
  directions). This is the junction-mask curvature math from the abandoned
  pipeline, extended from a scalar (max only) to a two-sided descriptor and
  promoted from a seeding hint to the segmentation signal itself.

**Classification rules (per edge):**
1. **Sharp**: `|θ_e| > --angular-tolerance` → BOUNDARY. Catches every
   non-tangent junction. The tolerance must stay above the tessellator's
   maximum dihedral step on smooth surfaces (coarse fillets step ~24°/row),
   so the default remains 25°; recalibrate against the corpus once 2.1 lands.
2. **Curvature jump** (tangent junctions, which are G1-smooth but
   G2-discontinuous): let
   `jump(e) = max(|σ_max(f1) − σ_max(f2)|, |σ_min(f1) − σ_min(f2)|)` and
   `scale(e) = max(|σ_max(f1)|, |σ_max(f2)|)`. BOUNDARY when
   `jump > J · scale` **and** `jump > κ_floor`. `J = 0.5` initially — the
   local-relative criterion carried over from the junction mask, whose
   measured rationale stands: thresholds relative to a mesh-global maximum
   fail because sharp creases dominate it while tangent junctions (the
   primary target) fall below any global fraction.
3. Everything else → INTERIOR.

`κ_floor` is a noise floor so exact-coplanar plane interiors (`κ ≡ 0`) and
fine tessellations of gentle curvature never fire: STL vertices are f32, so
per-edge curvature noise scales like `ε / d_e²` for position noise `ε`; use
`κ_floor = C · vertex_tolerance / d_e²` with `C` calibrated on the corpus.

**Why the two-sided descriptor matters:** `σ_max` alone misses tangent
junctions between surfaces sharing their largest curvature — the capsule
(pill) model's hemisphere↔cylinder junction has `σ_max = 1/r` on both sides,
but `σ_min` jumps from `1/r` (sphere) to `0` (cylinder), which rule 2
catches at full strength. Fillet↔plane (`σ_max`: `1/r` vs 0) and
fillet↔cylinder with sufficiently different radii fire on `σ_max`. A
fillet↔cylinder pair with similar radii (e.g. `1 − r/R < J`) can slip
through — that miss is recovered by stage 2.4's split repair, which re-runs
this classifier at progressively lower `J` inside regions that fail to
classify.

**Expected over-fire:** on coarse cones near the apex, `σ_max` varies
row-to-row (circumferential curvature grows toward the apex), so rule 2 may
split one cone into bands. Harmless: each band classifies as the same cone
and stage 2.4's merge rejoins them on parameter agreement.

Output: `EdgeClassification` keyed by undirected mesh edge, retaining the
class and the curvature evidence (`κ_e`, per-face descriptors) for reuse by
the 2.4 split repair, `--compare` diagnostics, and viz. Deterministic, O(E).

#### 2.2 Region extraction

Connected components of the face graph over INTERIOR edges (union-find).

Output (in `Stage2Output`):
- `face_region: Vec<usize>` — exactly one region per face (total coverage by
  construction).
- `regions: Vec<Region>` — member faces, member vertices, and ordered
  boundary loops (mesh-vertex loops per region, shared with the freeform
  parameterizer and structurally the same boundaries stage 3.1 will chain).

Regions are connected by construction. There is no minimum size: a single
facet is a valid region (and classifies as its own plane — this is the
terminal fallback that guarantees every model converts).

**Diagnostics** (report-only at this stage; repair may still move
boundaries): with `--compare`, sample midpoints of BOUNDARY edges and check
distance to the nearest reference STEP *edge* (within `--surface-tolerance` —
tangent-junction boundaries can legitimately sit one facet row off). Report
counts of suspicious boundary edges and of regions whose faces project onto
reference faces with disagreeing surface geometry (a missed-seam signal).

#### 2.3 Region classification

For each region, fit surface types in **ascending degrees of freedom** and
accept the first fit that validates. Simplest-first resolves degeneracy (any
plane also lies on an infinite-radius cylinder; a single facet is exactly a
plane) without coverage heuristics.

**The ladder:**

| Rung | Surface | DOF | Fit (see SURFACE_FITTING.md) |
|------|---------|-----|------------------------------|
| 1 | Plane | 3 | Area-weighted normal average |
| 2 | Sphere | 4 | Algebraic LSQ + Gauss-Newton |
| 3 | Cylinder | 5 | Normal-covariance axis + 2D circle fit + LM |
| 4 | Cone | 6 | Two-path (pointed/truncated) + LM; the existing code fits cylinder and cone together and this rung keeps that dual-fit internally, preferring the cylinder when both validate |
| 5 | Torus | 7 | Medial-axis tube-center + SOR axis + multi-axis LM rescue |
| 6 | Freeform B-spline | — | See "Freeform (B-spline) fitting" below; entered only via the 2.4 split-probe protocol |

**Validation (per fit, decisive by Design Principle #1):**
- Every region vertex within `--vertex-tolerance` of the surface.
- Every face centroid within `--surface-tolerance` (chordal sagitta ≈ h²/8R
  stays well inside for real tessellations).
- Orientation coherence: every face normal agrees in sign with the surface
  normal at its centroid — this determines the hypothesis `convex` flag and
  replaces the old per-face BFS convexity checks (it is what keeps the inner
  and outer walls of a thin pipe from ever sharing a surface; they are
  separate regions anyway, since the wall rim makes sharp edges).
- Numeric-validity guards retained from the fitters (these are fit-domain
  constraints, not selection heuristics): cone half-angle in [2°, 85°],
  torus `0 < r ≤ R`, sphere radius ≤ 1000× bounding-box diagonal.

**Fits see the whole region.** Initial parameter estimates use all faces
(area-weighted normal covariance, etc.), which conditions far better than
the abandoned pipeline's 3-face seeds — e.g. the pipe-elbow torus, which
previously needed special strip-merge seeding because narrow cylinder strips
had already tiled it, is now fitted once from the full tube and the (h,r)
profile is unambiguous. For large regions, fit on a deterministic vertex
subsample, validate against all vertices, and refit augmented with the worst
offenders on a near-miss (bounded retries).

**Sub-tolerance features**: when the fitted curvature radius is below
`--surface-tolerance` (3DBenchy has fitted cylinders of r=0.2mm against the
0.4mm default), centroid checks are vacuous at that scale; classify by
vertex tolerance alone and mark the region, since stage 3's chain validation
has the same blind spot there (see TODO, face-5841 family).

**Output**: exactly one hypothesis per region, appended to the typed vectors
in `Stage2Output` (`planar_hypotheses`, `conical_hypotheses` — cylinders are
cones with `half_angle == 0` — `spherical_hypotheses`,
`toroidal_hypotheses`, `bspline_hypotheses`). Hypothesis structs keep their
existing fields: type-specific parameters, `faces`, `vertices`, `convex`,
and error statistics (`error_max`, `centroid_error_max`).

##### Freeform (B-spline) fitting

The final rung, reached only after the 2.4 split-probe protocol has ruled
out "this region is really several analytic surfaces with a missed seam".
The fitting machinery carries over intact from the previous architecture
(implemented 2026-07-10 through 2026-07-15, `src/stage2/bspline.rs`); what
is deleted is its old region-*gathering* front end (see Abandoned
Approaches) — the partition now supplies the region. A planar margin of a
freeform surface is curvature-continuous with it and lands in the same
region automatically, which is what the old "tangent-plane merge" step
existed to fix.

**Algorithm (`deduce_bspline_surfaces`, per region):**
1. Regions planar to fitting exactness (10× vertex tolerance, NaN-guarded:
   the area-weighted normal of a closed region cancels to zero) become a
   multi-face planar hypothesis instead — a plane is preferable to a
   B-spline that overfits noise. Regions below 4 faces are not fitted (the
   2.4 fallback ladder handles them).
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

#### 2.4 Repair loop

Runs after initial classification; fixes the two 2.1 failure modes
(over-segmentation and missed seams) plus row-level boundary error, by
validation-driven local operations. Iterate the three operations in a
deterministic order (by region index; smallest region first for merges) to a
fixpoint, bounded rounds.

**Merge** — repairs over-segmentation (coarse-cone banding; any residual
same-surface fragmentation): merge region S into an adjacent larger region T
when every vertex of S lies within `--vertex-tolerance` of T's fitted
surface and every centroid of S within `--surface-tolerance`; refit the
union, revalidate, keep only if valid. This is the abandoned absorb pass's
*tight* rule reborn with a principled criterion (vertex tolerance, at which
"on the same surface" is decisive — never the chordal slop of the abandoned
loose rule, which put absorbed facets up to 0.4mm off their claimed plane
and caused 794 BRepCheck failures on planar-only 3DBenchy). It also folds in
chord-degenerate narrow strips (a single-facet-wide band of a coarse
cylinder validates rung-1 as its own facet plane, but its vertices lie
exactly on the neighbor cylinder, so merge claims it).

**Split** — repairs missed seams (a region that fails every rung is
presumed to be several surfaces): re-run the 2.1 curvature-jump classifier
on the region's internal edges at progressively lower thresholds (J/2, J/4,
…) until the region separates; classify the components. Acceptance
preference, in order:
1. A split whose components ALL classify analytically → accept (this is how
   a fillet↔cylinder junction with similar radii — below the global J — is
   recovered).
2. Otherwise, try the freeform rung on the *whole* region (a genuine
   freeform surface has no internal jump structure worth splitting on; this
   ordering also prevents a smooth-but-merged analytic pair from
   reconstructing as one plausible-but-wrong freeform patch that would erase
   a reference edge).
3. Otherwise accept the best split and recurse on its components. Terminal
   case: a single facet always classifies as a plane, so this ladder cannot
   fail to terminate — the worst case degrades to per-facet planes (exactly
   the planar-only mode that runs 3DBenchy end-to-end), never to a
   conversion failure.

**Boundary relabel** — repairs ±1-row boundary error (the 2.1 classifier and
tessellation can place a tangent boundary one facet row off the true
junction, which would shift stage 3's reconstructed boundary vertices off
the feature curves): walk faces with a neighbor in a different region;
compute the face's max vertex distance to its own surface and to each
adjacent region's surface; move the face when an adjacent surface fits at
least **2× better** (margin prevents oscillation for true junction faces,
whose vertices lie on both surfaces equally). Iterate to a fixpoint (≤ 10
rounds; boundaries only move locally), then rebuild affected regions'
vertex lists. Region *geometry* is not refitted — the fits were validated
against the region interior; only the membership map changes. This is the
previous architecture's `refine_boundary_assignments` promoted from a
fix-up pass to a core mechanism. After relabeling, re-split any region that
became disconnected (each fragment keeps the shared fitted geometry as its
own region — connectivity is the invariant; identical geometry on two
non-adjacent regions is legal).

**Post-2.4 invariants (hard assertions before stage 3):**
- Every mesh face belongs to exactly one region (total coverage).
- Every region is mesh-connected.
- Every region's surface validates the region's own vertices/centroids.

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
- `--surface-types=<list>` (new) restricts the classification ladder
  (default: all). `--surface-types=planar` is **exact planar-only mode**:
  curved regions fail the one-rung ladder and cascade through split repair
  to per-facet planes — the mode that converts 3DBenchy with volume matching
  the official 15.55 cm³ — now without any chordal-slop absorb by
  construction.
- `--skip-stage` is retired with the old stage numbering; its
  type-skipping role is `--surface-types`, and its absorb/refinement gates
  have no equivalent (absorb is gone; repair is core, and debugging-level
  disabling of repair operations can be a `-Z`-style debug option if ever
  needed).

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

The jump from the BFS pipeline to partition-then-classify is immediate — no
compatibility flag, no parallel pipelines, and **no reconciliation of old
and new stage numbers**. The moment 2.1/2.2 land, the new numbering owns
`--stage` and `--viz-stages`, and the old BFS pipeline is disconnected from
`stage2()` — nerfed in place (its code stays in-tree as reference until the
deletion pass, but nothing calls it). `stage2()` runs exactly the new
sub-stages implemented so far and stops at that frontier with no selected
surfaces, so stages 3/4 and most integration tests stay red until 2.3–2.5
land. Models are re-enabled as they pass, easiest first. Suggested order:

1. **Stage 2.1 + 2.2** (edge classifier, region extraction): new
   `src/stage2/edges.rs` + `src/stage2/regions.rs`, plus the report-only
   partition diagnostics. **Viz** (primitives verified against viz.rs):
   boundary edges render as `LineOverlay` segments — arbitrary colored 3D
   lines with a `no_depth_test` option for X-ray inspection — color-coded
   by class (sharp vs curvature-jump); regions render as one
   `FaceHighlight` per region with a deterministic palette (golden-ratio
   hue stepping on region index). `EdgeHighlight` is NOT the right
   primitive — it outlines all edges of listed *faces*, not individual mesh
   edges — so it stays unpopulated and remains TODO.md's delete candidate.
   `--viz-stages=2.1,2.2` needs no parser changes (`parse_stage` accepts
   any major.minor): 2.1 pauses once showing classified edges over the
   mesh, 2.2 pauses once showing colored regions with the boundary overlay
   (`--viz-stages-level` has no finer granularity to offer here; a
   later refinement could pause per classification round at level 2).
   Calibrate `J` and `κ_floor` across the corpus by inspecting partitions
   (viz + diagnostics) against reference STEP face structure *before* any
   fitting exists. This step is also the decision gate for the fallback
   design: if no calibration produces reliable partitions on the corpus,
   switch to
   [seed-first multi-type competition](#contingency-seed-first-multi-type-competition)
   rather than building 2.3/2.4 on bad partitions.
2. **Stage 2.3 analytic ladder**: new `src/stage2/classify.rs` calling the
   existing fitters (`planar.rs`, `cylindrical.rs`, `conical.rs`,
   `spherical.rs`, `toroidal.rs` keep their fitting/validation functions;
   their seeding/BFS entry points go). Terminal per-facet-plane fallback
   from day one, so every model produces output immediately.
3. **Stage 2.4 repair**: `src/stage2/repair.rs` — merge, split (reusing the
   2.1 classifier with a threshold parameter), and boundary relabel (port of
   `refine_boundary_assignments`). Connectivity re-split + the three hard
   invariants.
4. **Freeform rung rewiring**: `bspline.rs` keeps parameterization → grid →
   fit → validate; its gathering front end is deleted; hook it behind the
   2.4 split-probe protocol.
5. **Deletion pass**: `bfs.rs` (generation loop, seed machinery, masks,
   acceptance tests), greedy selection, `absorb_redundant_surfaces`, the
   `MeshFace` hypothesis fields, `--skip-stage`. Update SURFACE_FITTING.md's
   seeding/BFS/acceptance sections with pointers to Abandoned Approaches.
6. **Test migration**: stage-2 integration tests re-target region counts and
   classification results (end-state expectations — cube = 6 planar,
   block_with_hole = 6 planar + 1 cylindrical, etc. — are unchanged);
   re-enable corpus models in roughly this order: pure-planar ccad models →
   analytic primitives → fillet/tangent models (pill, rounded cubes,
   stepped_block_toroidal) → nurbs suite → nosecone/nosecone_truncated →
   pipe_elbow → 3DBenchy.

Expected wins to verify along the way: pipe_elbow's torus recovered without
strip-merge seeding; nosecone with no absorb pass; stepped_block's
single-rectangle walls as real planar regions; stage-2 runtime on large
models collapsing (each region fit O(1) times — no more tiling large spheres
with thousands of band trials; full stage 2 on 3DBenchy was ~15 min under
BFS vs ~3 s planar-only).

---

## Contingency: Seed-First Multi-Type Competition

Not abandoned — considered (2026-07-17) and not chosen. This is the fallback
design if stage 2.1's curvature-jump calibration proves unreliable in
practice (coarse fillets, f32 quantization noise defeating `J`/`κ_floor`
tuning). Recorded so the option and its rationale don't get lost.

**The design**: keep fit-driven BFS growth, but invert the old pipeline's
loops. Iterate over seed faces in curvature-score order (interior faces
first — the existing seed-ordering heuristics); skip faces already claimed
by a committed hypothesis. For each seed, run *every* surface type's
fit-and-expand from that same seed, and commit the hypothesis with the
largest supported area (with a simplest-DOF tiebreak). Expansion may grow
through already-claimed faces — junction faces whose vertices genuinely lie
on two surfaces end up claimed by both — but claimed faces never seed. After
all faces are claimed, resolve overlaps into a strict partition (per-face
best-fit among claimants, the boundary-relabel rule) and re-split any
disconnected remnant.

**Why it is sound, not heuristic** (under the exactness assumption): an
accidental fit is by definition supported on a *sub-region* — the true
surface fitted from an interior seed expands, with vertex tolerance as the
arbiter, to the entire region, so when all types compete on the *same* seed,
largest-area-wins selects the true surface. The old pipeline compared
hypotheses born from different seeds at different times, which is why commit
order mattered. Cross-region leakage is bounded: two distinct analytic
surfaces agreeing to vertex tolerance over an extended patch would be
near-identical, so overlap is confined to genuine junction rows.

**What it fixes without any masks or ordering**: fragmentation (an
interior-seeded surface covers its whole region regardless of prior
commits), stage ordering (per-seed competition replaces it), and the
accidental-fit tiling cost (each region pays ~5 bounded expansions once,
not per-face — the large_sphere pathology disappears without masking). Its
distinctive strength over partition-then-classify: **boundary placement
never depends on discrete curvature estimates** — the fit itself finds the
boundary from exact vertex distances (the capsule's equal-radius
sphere↔cylinder junction needs no two-sided descriptor here), and curvature
ordering errors cost only efficiency, never correctness.

**Why it was not chosen**:
1. Small-seed fit conditioning is unsolved — fits still start from 1-rings/
   fans, so the pipe-elbow torus (strip tessellation; a 1-ring spans no
   major arc) again decomposes into cylinder strips. The degenerate-limit
   growth model below repairs this, at the cost of the machinery it
   describes. Partition-first instead fits the whole region at once and
   needs none of it.
2. Freeform stays second-class: a B-spline patch cannot compete from a seed
   (parameterization needs a bounded region), so freeform surfaces get
   tiled by their largest accidental analytic patches and the old
   leftover-gathering + confetti-absorption heuristics return. (The
   degenerate-limit model does not help here — no analytic family contains
   B-splines.)
3. The residual heuristic surface is larger: overlap resolution + remnant
   re-split, degeneracy guards and the DOF tiebreak (on a plane region a
   near-infinite-radius cylinder *ties* the plane on area), a "good enough
   to suppress seeding" definition, and load-bearing seed ordering.

**Required upgrade if adopted — degenerate-limit growth model** (identified
2026-07-17): grow with a model family whose degenerate cases are *regular
interior points* of the parameter space, so the fit explores as the simpler
surface until the data pulls it off the limit. Concretely for the strip
problem: parameterize the torus by a point on the tube spine, spine
direction, tube radius `r`, spine **curvature** `k = 1/R`, and a
bending-plane azimuth — `k = 0` *is* the cylinder, with no infinity
anywhere (the Pratt/Taubin curvature parameterization that circle fitters
use to make lines a regular case). Growth starts at `k ≈ 0` on one strip;
the refit that admits the neighboring strip makes `k` identifiable and LM
discovers the bend — no seed-size machinery, no discrete promotion event.
Three consequences to design for:
- *The refit trigger, not the fitter, is the actual blocker.* A neighboring
  strip sits ~`R·Δφ²/2` off the current fit — hundreds to thousands of
  times vertex tolerance, far outside the `REFIT_SKIP_MULTIPLIER` (10×)
  window that gates refit attempts today. The general model needs refits
  triggered by *chordal* (surface-tolerance) proximity, or batched per
  stalled frontier; this is where the performance cost concentrates, since
  refits run LM over the whole grown region.
- *The underdetermined phase needs regularization.* While support spans one
  strip, `k` is unidentifiable; ridge the fit toward `k = 0` and relax the
  bias as support grows. The exactness assumption bounds the
  over-flexibility risk: added flexibility cannot absorb faces of a
  genuinely different surface, because within-vertex-tolerance remains the
  acceptance bar and two distinct surfaces agreeing to that tolerance over
  an extended patch would be near-identical.
- *A torus-only generalization re-creates the dual pathology*: a bent-tube
  model that cannot degenerate further locks onto a sphere's latitude band,
  closes the ring, and stalls — spheres get tiled by torus bands again. The
  logical endpoint is growing with a family containing **all** the analytic
  primitives as degenerate cases — the Unified SOR framework of
  SURFACE_FITTING.md (Dupin cyclides are the mathematically complete
  version) — after which classification becomes model *reduction* (collapse
  `k ≈ 0` to cylinder, a degenerate spine to sphere) instead of a ladder of
  fitters. Freeform stays outside any such family.

If adopted, the 2.4 repair loop, the post-stage-2 invariants, and the 2.5
output contract apply to its results unchanged — only 2.1–2.3 would be
replaced by the seed-competition driver (plus the degenerate-limit growth
model above).

---

## Abandoned Approaches

Retired with the 2026-07-17 stage-2 redesign. Recorded here because each
worked well enough to ship models, and some may be worth resurrecting as
*diagnostics* — but all of them existed to compensate for fit-driven greedy
segmentation, and the compensations kept multiplying (that treadmill is why
the architecture changed). Do not re-add one without identifying which
design principle it would compensate for.

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
- **2.1**: edge classification — BOUNDARY edges drawn as `LineOverlay`
  segments color-coded by class (sharp vs curvature-jump), with
  `no_depth_test` available for X-ray inspection. (`EdgeHighlight` outlines
  whole faces rather than individual mesh edges, so it stays unpopulated —
  still TODO.md's delete candidate.)
- **2.2**: regions — one `FaceHighlight` per region with a deterministic
  color palette, boundary `LineOverlay` retained for context.
- **2.3**: per-region classification — level 1 pauses at each classified
  region (fit overlay + parameters via the existing `viz_overlays` /
  `viz_fit_description` per-surface-type functions); level 2 pauses at each
  attempted rung.
- **2.4**: repair — level 1 pauses after the loop converges; level 2 pauses
  at each merge/split/relabel operation.
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

- `src/stage2/mod.rs` — orchestration, `Stage2Output` (`face_region`,
  `regions`, `EdgeClassification`, typed hypothesis vectors,
  `SelectedSurface`), shared helpers.
- `src/stage2/edges.rs` — 2.1 curvature field + edge classifier
  (threshold-parameterized so 2.4's split repair can reuse it).
- `src/stage2/regions.rs` — 2.2 union-find extraction, boundary loops.
- `src/stage2/classify.rs` — 2.3 ladder, subsample/validate/refit driver.
- `src/stage2/repair.rs` — 2.4 merge/split/relabel, invariants.
- `src/stage2/{planar,cylindrical,conical,spherical,toroidal}.rs` — fitting
  + validation for each surface type (seeding/BFS entry points deleted).
- `src/stage2/bspline.rs` — freeform rung (parameterize/resample/fit/
  validate; gathering deleted).
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
- Edge classification / partition sanity per model (with `--compare`: region
  boundary edges near reference edges; region interiors mapping to
  consistent reference surfaces).
- Region classification end states — the same expectations the old suite
  pinned (cube = 6 planar; block_with_hole = 6 planar + 1 cylindrical;
  pipe = 2 planar + 2 cylindrical; ball_on_cylinder = 1 planar +
  1 cylindrical + 1 spherical; hemisphere = 1 planar + 1 spherical;
  cone_15x20 = 1 planar + 1 conical; etc.), plus fitted-parameter agreement
  with the reference STEP (axis, radius, center) for analytic models.
- Repair-loop behavior on the models that motivated it (coarse cones
  re-merged; capsule junction found; stepped_block narrow walls kept as
  real planar regions).

**Migration reality**: the suite is red during the rewrite (see Rewrite
Roadmap step 6). Re-enable per model; a model stays disabled with a comment
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
