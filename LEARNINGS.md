## opencascade-sys API Patterns

The Rust bindings mirror OCCT C++ closely but with automatic name transformations:

### Module and Type Naming
- `BRepBuilderAPI_MakeEdge` → `b_rep_builder_api::MakeEdge`
- `BRep_Tool` → `b_rep::Tool`
- `gp_Pnt` → `gp::Pnt`
- `TopExp_Explorer` → `top_exp::Explorer`
- `GeomAPI_ProjectPointOnSurf` → `geom_api::ProjectPointOnSurf`
- `STEPControl_Reader` → `step_control::Reader`
- `RWStl` → `rw_stl` (free functions)

### Constructors
Named `::new_<compressed_param_types>()`. For example:
```rust
gp::Pnt::new_real3(x, y, z)           // gp_Pnt(double, double, double)
gp::Dir::new_real3(dx, dy, dz)
top_exp::Explorer::new_shape_shapeenum2(shape, from_type, avoid_type)
geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(pt, surf, algo)
```

### Ownership
- All objects are returned as `OwnedPtr<T>` which auto-derefs via `Deref`/`DerefMut`
- `OwnedPtr<T>` is essentially `Box<T>` for OCCT heap objects

### Handles (OCCT `Handle<T>`)
- `Type::to_handle(owned_ptr)` consumes `OwnedPtr<T>` → `OwnedPtr<HandleT>`
- `.get()` dereferences a handle to access the underlying object
- Upcasting: `.to_handle_surface()`, `.to_handle_curve()`, etc.
- Downcasting: `.downcast_to_plane()` returns `Option<OwnedPtr<HandleDerived>>`

### Enums
Typed Rust enums, not raw integers:
```rust
top_abs::ShapeEnum::Face
extrema::ExtAlgo::Grad
```

### Shape Downcasting
`topo_ds::face_shape` is a safe function, but it panics if the shape is not a Face. Always guard with `shape_type()`:
```rust
if shape_ref.shape_type() == top_abs::ShapeEnum::Face {
    let face = topo_ds::face_shape(shape_ref);  // &topo_ds::Shape → &topo_ds::Face
    ...
}
```

### Static Methods → Free Functions
OCCT utility classes (e.g., `BRep_Tool`) become module-level free functions or associated functions on the type.

## Common API Recipes

### Read an STL file → vertex positions
```rust
use opencascade_sys::{rw_stl, message};
let progress = message::ProgressRange::new();
let tri_handle = rw_stl::read_file_charptr_progressrange_2(path, &progress);
let tri = tri_handle.get();  // &Poly_Triangulation
for i in 1..=tri.nb_nodes() {
    let pt = tri.node(i);  // OwnedPtr<gp::Pnt>
    let (x, y, z) = (pt.x(), pt.y(), pt.z());
}
```

### Read a STEP file → TopoDS_Shape
```rust
use opencascade_sys::{step_control, message};
let mut reader = step_control::Reader::new();
reader.read_file_charptr(path);
reader.transfer_roots(&message::ProgressRange::new());
let shape = reader.one_shape();  // OwnedPtr<topo_ds::Shape>
```

### Iterate faces of a shape
```rust
use opencascade_sys::{top_exp, top_abs, topo_ds};
let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
    &shape, top_abs::ShapeEnum::Face, top_abs::ShapeEnum::Shape,
);
while explorer.more() {
    let face = topo_ds::face_shape(explorer.value());
    // ... use face
    explorer.next();
}
```

### Project a point onto a surface (compute distance)
```rust
use opencascade_sys::{b_rep, geom_api, extrema};
let surface = b_rep::Tool::surface_face(face);  // OwnedPtr<HandleGeomSurface>
let projector = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
    &point, &surface, extrema::ExtAlgo::Grad,
);
if projector.is_done() && projector.nb_points() > 0 {
    let dist = projector.lower_distance();
}
```

### Poly_Triangulation index convention
All node/triangle indices are **1-based** (OCCT convention).


### Bounded vs unbounded distance: BRepExtrema vs GeomAPI_ProjectPointOnSurf
`GeomAPI_ProjectPointOnSurf` projects onto the *infinite* underlying `Geom_Surface`, not the bounded face. Use `BRepExtrema_DistShapeShape` (via `b_rep_extrema::DistShapeShape`) instead to compute distance to bounded `TopoDS_Face` shapes. Create a `BRepBuilderAPI_MakeVertex` from the point, then compute `DistShapeShape` between the vertex shape and the target shape. This correctly reports nonzero distance for points that lie on an infinite surface extension but outside the bounded face.

### Stage 2.1 planar hypothesis algorithm
- BFS region growing from seed faces using angular + vertex-distance tolerance.
- Angular threshold of `1e-6` (1 - cos ≈ 0.08°) prevents grouping faces with slightly different normals.
- Re-fitting via area-weighted normal averaging + vertex centroid for plane distance.
- For CAD meshes, vertices lie precisely on surfaces, so fit errors are typically 0 or near-epsilon.
- Coplanar but disconnected faces correctly get separate hypotheses (not a bug).
- The `vertex_tolerance_mm` config parameter controls the maximum vertex-to-plane distance for grouping.

### CodeCAD (ccad) test model generation
- `ccad build` in `tests/ccad/` generates both STL and STEP files from Lua scripts.
- `wedge(dx, dy, dz, ltx)` requires `ltx > 0` (cannot be 0).
- Boolean operations like `union()` work with pre-translated shapes.
- Parts are registered in `project.json`; each entry needs a unique `id`.
- Available primitives: `box()`, `cylinder(diameter, height)`, `sphere(diameter)`, `cone(d1, d2, height)`, `wedge()`, `hex_prism()`.
- Boolean operations: `union()`, `difference()`, `intersection()`.
- Transformations: `translate()`, `rotate_x/y/z()`, `scale()`, `center_xyz()` etc.
- Construction: `extrude(face, height)`, `revolve(profile, angle)`.
- Edge ops: `fillet_all(s, r)`, `chamfer_all(s, d)`.

### CAD tessellation quad-strip patterns
- Tessellated cylinders and spheres from CAD tools (CodeCAD/OCCT, Onshape) use quad-strip patterns where each quad facet splits into 2 coplanar triangles.
- Solution: stage 1.3 merges coplanar triangle pairs back into quads, so curved-surface facets become single faces with single-face planar hypotheses.
- This restores the clean "single-face = curved surface candidate" criterion for stages 2.2 and 2.3.
- Sphere poles are an exception: triangles in the polar fan pattern don't have coplanar neighbors and correctly remain as triangles.

### Stage 2.2 cylindrical hypothesis algorithm
- Two-step fitting: (1) axis direction from smallest eigenvector of area-weighted normal covariance matrix M = Σ wᵢ nᵢ nᵢᵀ, (2) axis position and radius via 2D algebraic circle fit after projecting vertices perpendicular to axis.
- **Critical: center the 2D projected coordinates before circle fitting.** Without centering, the normal equations matrix becomes catastrophically ill-conditioned when absolute vertex coordinates are large relative to the arc span (e.g., vertices at (15, 12.5) spanning 2.5° of a r=2mm cylinder). The sums-of-squares terms (sx2 ≈ 1350) dwarf the signal (radius ≈ 2), causing Cramer's rule to produce wildly wrong radii (e.g., 139mm instead of 2mm). Fix: subtract mean_x/mean_y before fitting, add back when converting 2D center to 3D. This is the standard numerical remedy for polynomial least-squares fitting.
- 3×3 symmetric eigenvalue computation uses depressed cubic (trigonometric solution with cos(θ/3)), plus null-space eigenvector via cross products of (M - λI) rows.
- Multi-seed evaluation with 3-face seeds: for each starting face fi, finds triples (fi, n1, n2) where n1 is a pairwise-qualified neighbor of fi, and n2 is a pairwise-qualified neighbor of n1 that is NOT a neighbor of fi. Runs a trial BFS for each triple using temporary data (no mesh mutation during trials), keeps the best trial (most total area).
- BFS seeds from triples of faces spanning more circumference. Pairwise qualification requires non-parallel normals (cross product > MIN_CROSS_THRESHOLD=0.01) and dihedral angle ≤ angular tolerance. No restrictions on multi-face planar faces.
- Angular coverage validation applied to each trial after BFS: requires 3+ angular clusters of face centroids around the cylinder circumference. Rejects bogus hypotheses where faces are clustered on one side. Early termination is disabled as future work.
- Centroid validation after BFS: each face centroid must be within surface_tolerance of the fitted cylinder. This rejects spurious fits from perpendicular planar faces that algebraically fit a cylinder.
- Minimum 3-face requirement: cylinder hypotheses must contain at least 3 faces. Any real cylindrical surface from a CAD tessellation will have at least 3 facets around its circumference. This eliminates spurious 2-face cylinder fits from locally adjacent faces on tori, cones, and spheres (e.g., cone went from 29 to 0 spurious fits, pipe_elbow from 252 to 51).
- Angular tolerance check: both seed pairs and BFS expansion candidates must have dihedral angle ≤ `--angular-tolerance` (default 17.5°). This prevents fitting cylinders to faces that meet at sharp angles (e.g., cube faces at 90°). CAD tessellators enforce a maximum angular tolerance between adjacent triangles on the same surface.
- Extracting cylinder parameters from STEP: `b_rep_adaptor::Surface::cylinder()` → `gp::Cylinder`, then `.radius()`, `.axis().direction()`, `.axis().location()`. Useful for test verification.

### ccad box() coordinate system
- `box(w, h, d)` places the box in the positive octant: (0,0,0) to (w,h,d). It does NOT center at origin.
- To center a box, use `translate(box(w,h,d), -w/2, -h/2, -d/2)`.
- `cylinder(diameter, height)` is vertical (Z-axis), centered at the XY origin, extends from z=0 to z=height (NOT centered vertically). To make through-holes, translate cylinders vertically so they extend past both sides of the target body.
- `sphere(diameter)` is centered at the origin. To make hemispherical pockets or domes, use `intersection()` or `difference()` with a box.

### Stage 2.3 spherical hypothesis algorithm
- Algebraic least-squares sphere fitting: expand |v-c|²=r² linearly into 4 unknowns (cx, cy, cz, k=r²-|c|²). Normal equations AᵀAx = Aᵀb solved via 4×4 Gaussian elimination with partial pivoting.
- Vertex-neighborhood seeding: for each mesh vertex, collect surrounding non-excluded faces (≥3), fit sphere from all their vertices. Naturally provides angular diversity for the 4-DOF sphere fit. Superior to pair-based seeding which failed on finely tessellated spheres.
- Multi-face planar faces excluded from seeding but absorb-able by BFS, same as cylindrical.
- Solid-angle coverage validation: area-weighted 3×3 covariance matrix of centroid-to-center unit direction vectors. Eigenvalues via depressed cubic (same as cylinder axis). Require λ_min/λ_max ≥ MIN_SPHERE_EIGENVALUE_RATIO (0.01). Rejects strip-like hypotheses from cylinder fillet growth.
- MAX_SPHERE_RADIUS_FACTOR = 1000.0 (was 10.0). With solid-angle coverage + surface-tolerance in BFS, radius limit only needs to prevent numerical precision issues.
- Centroid validation during BFS (not just post-BFS): each face centroid checked against surface_tolerance with 2x hard skip threshold.
- Minimum 4-face requirement: a sphere needs at least 4 non-coplanar points for a unique fit.
- All faces are candidates for BFS expansion (including those with cylindrical hypotheses and multi-face planar hypotheses). Only seeding excludes multi-face planar. Stage 2.6 resolves overlapping hypotheses.
- Torus faces locally fit spheres (e.g., pipe_elbow produces 36 sphere hypotheses). This is expected — stage 2.6 with torus hypothesis support will resolve it.
- Concave spheres (bowls/pockets) have normals pointing toward center — same convexity logic as cylinders.
- Angular tolerance check: same as cylindrical — seed faces and BFS neighbors must have dihedral angle ≤ `--angular-tolerance`. This fixed the spurious sphere detection on the 1mm manual cube.

### Stage 2.6 surface selection
- Greedy area-based selection: iteratively selects the hypothesis with the largest total remaining mesh face area, assigns its unassigned faces, removes those faces from all other candidates, repeats.
- Multi-face planar, cylindrical, and spherical hypotheses are candidates. Single-face planar hypotheses are fallback for unclaimed faces.
- After selection, hypothesis face/vertex lists are updated in-place to reflect only the faces actually assigned.
- Replaced the old per-face priority rule (spherical > cylindrical > multi-face planar > single-face planar) which failed when bogus small-area hypotheses of a high-priority type beat correct large-area hypotheses of a lower-priority type.
- Sphere BFS along cylinder fillets is prevented by solid-angle coverage validation in stage 2.3. Vertex-based seeding also helps by providing better initial seeds.
- Compare tolerance uses `surface_tolerance` for all hypothesis types (stages 2.1-2.3 and 2.6 all use the same tolerance). Previously used `vertex_tolerance` which was too tight for centroid-to-surface distance checks.


### Stage 3.1 adjacency graph construction
- OCCT surface creation: `geom::Plane::new_pnt_dir`, `geom::CylindricalSurface::new_ax3_real`, `geom::SphericalSurface::new_ax3_real`. Type erasure via `.to_handle().to_handle_surface()` to get `OwnedPtr<HandleGeomSurface>`.
- Adjacency is determined from mesh edges: an edge between mesh faces assigned to different selected surfaces is a boundary edge.
- Boundary edges are grouped by `SurfacePair` (canonical ordering), then chained into ordered vertex sequences. Each chain becomes a `ReconEdge`.
- Corner vertices (where 3+ surfaces meet) become `BRepVertex` entries.
- Topological ordering of edges around each face uses rotation-order walk: at each vertex, pick the edge to the face that comes next in counter-clockwise order. For faces without vertices (closed loops like cylinder caps), all edges incident on that face are included.
- Bug encountered: closing vertex in topological walk was missing — after walking N edges around a face, the vertex connecting the last edge back to the first edge must be added separately.
- `wedge(10, 10, 10, 5)` in CodeCAD is a truncated wedge (8 vertices, 6 faces), NOT a triangular prism (6 vertices, 5 faces). The 4th parameter `ltx` is the top face X dimension.
- Euler formula V-E+F=2 is a good sanity check for genus-0 (no holes through the model) manifold solids. Models with through-holes have different Euler characteristics.
- Stage 3.1 --compare check: validates that (a) ReconEdge boundary mesh vertices lie within vertex_tolerance of STEP edges, and (b) BRepVertex positions lie within vertex_tolerance of STEP vertices. Uses `topo_ds::Compound` to build compounds of all STEP edges and vertices separately, then `BRepExtrema_DistShapeShape` against those compounds.
- Building a `topo_ds::Compound`: use `topo_ds::Builder::new()`, `builder.make_compound(&mut compound)`, then `builder.add(compound.as_shape_mut(), shape)`. Also available via `b_rep::Builder` which inherits from `topo_ds::Builder`.
- Shape downcasting: `topo_ds::edge(shape)`, `topo_ds::vertex(shape)`, `topo_ds::face_shape(shape)` etc. These panic on type mismatch — check `shape.shape_type()` first.

### Stage 3.2 tangency detection
- Surface normals at boundary points can be computed analytically for planar/cylindrical/spherical surfaces without needing OCCT UV evaluation. Plane: constant normal. Cylinder: radial direction from axis (negate for concave). Sphere: direction from center (negate for concave).
- Tangency threshold of 2° (cos ≈ 0.9994) works well. Tangent edges arise from fillets (e.g., cylinder tangent to plane at a fillet edge). part_rounded_cube_10_r2 has 8 tangent edges (4 cylinders × 2 touching planes each).

### Stage 3.3 edge curve computation
- `GeomAPI_IntSS::new_handlegeomsurface2_real(S1, S2, Tol)` computes surface-surface intersection curves. Returns 1-indexed curves via `.nb_lines()` and `.line(i)` which returns `&HandleGeomCurve`.
- When IntSS returns multiple curves, select the one closest to mesh boundary vertices by projecting sampled boundary points onto each curve via `GeomAPI_ProjectPointOnCurve` and picking minimum total distance.
- `GeomAPI_ProjectPointOnCurve::new_pnt_handlegeomcurve(&pt, curve)` projects a point onto a curve. `.lower_distance_parameter()` returns the parameter of the closest point.
- `Geom_TrimmedCurve::new_handlegeomcurve_real2(curve, U1, U2)` trims a curve to a parameter range. Convert to handle via `Geom_TrimmedCurve::to_handle(owned).to_handle_curve()`.
- For closed-loop edges (no vertex endpoints, e.g., full circle on a cylinder cap), use the curve's full parameter range `[first_parameter, last_parameter]`.
- `ReconEdge.curve_3d` stores `Option<OwnedPtr<geom::HandleGeomCurve>>`. `OwnedPtr` does not implement `Debug`, so `ReconEdge` needs a manual `Debug` impl.
- All current non-tangent edges compute 100% of edge curves successfully via IntSS (cube 12/12, cylinder 2/2, hemisphere 1/1, ball_on_cylinder 2/2, block_with_hole 15/15, pipe 4/4, spherical_pocket 13/13, chamfered_cube 48/48, part_rounded_cube 16/16).
- **Tangent edge curves**: `GeomAPI_IntSS` fails for tangent surfaces (returns no curves). Instead, construct tangent edge curves analytically. For plane-cylinder tangencies: compute the radial direction from the component of the plane normal perpendicular to the cylinder axis, find the tangent point at `axis_origin + r * radial_unit`, construct a `Geom_Line` through that point along the cylinder axis direction, and trim to the vertex endpoint parameters. Validate that boundary vertices lie within tolerance of the constructed curve.
- **Sphere-cylinder tangent edge curves**: The tangent curve between a sphere and a cylinder (at fillet corners) is a great circle arc on the sphere. Construct a `Geom_Circle` using the sphere's center and radius (NOT the cylinder's, which may be imprecisely fitted). The circle plane normal is `normalize(cross(v0 - center, v1 - center))` where v0, v1 are the endpoint vertex positions. This ensures the circle passes through both vertices to the sphere fit accuracy (~1e-6mm), independent of cylinder parameter quality. Use `gp::Ax2::new_pnt_dir(&center_pnt, &normal_dir)` → `geom::Circle::new_ax2_real(&ax2, radius)` → `.to_handle().to_handle_curve()`. For arc selection on the periodic circle, sample mesh boundary vertices to determine forward vs. reverse arc.
- **Closed-loop sphere-cylinder tangent edge curves** (no vertex endpoints, e.g. pill/capsule shapes): Construct a full `Geom_Circle` with the circle plane normal set to the cylinder axis direction and center at the sphere center. Use the full parameter range [0, 2π] as a `TrimmedCurve`.
- **Adaptive vertex tolerance for MakeEdge**: When surface fits are imprecise (e.g., cylinders fitted from limited mesh data), IntSS intersection curves may not pass exactly through mesh vertex positions. Before creating MakeEdge, compute the distance from each vertex to the curve at its assigned parameter, and update the OCCT vertex tolerance to at least that distance (with 1% margin). This replaces the fixed `vertex_tolerance_mm` (1e-5mm) with position-adaptive tolerances that accommodate surface fit imprecision.

### Stage 3.4 face creation
- **Planar faces**: Create `TopoDS_Edge` from trimmed curve using `BRepBuilderAPI_MakeEdge::new_handlegeomcurve(curve)`. Group edges into wire loops by `BRepVertex` vertex connectivity (union-find style). Build `MakeWire` from edges, `MakeFace` from surface + outer wire + hole wires.
- **Periodic faces (cylinder/sphere) with edges**: Two approaches depending on the edge type:
  - **Full revolution** (all boundary edges are closed loops): Use UV-bounds construction via `MakeFace::new_handlegeomsurface_real5(surface, umin, umax, vmin, vmax, tol)`. This automatically creates seam edges needed for proper B-Rep topology. Circular gap algorithm computes correct u-bounds even when edges cross the 0/2π seam.
  - **Partial revolution** (open edges with vertices): Use wire-based construction WITH pre-set pcurves via `BRep_Builder::update_edge` before calling `MakeFace`. The pcurves (Geom2d_Line) map the edge's 3D parameter range to the surface's UV space. This ensures correct arc selection on the periodic surface and shares edges with adjacent planar faces for proper sewing.
- **Full spheres** (no edges): Use `MakeFace::new_handlegeomsurface_real(surface, tolerance)` which uses natural surface bounds.
- `MakeEdge::edge()` takes `&mut self`, not `&self`. Need `&mut` references when accessing edges.
- `BRepCheck_Analyzer::new_shape_bool(shape, true)` validates faces. `.is_valid()` returns bool.
- The base `geom::Surface` type does NOT have `bounds()` or `first_u_parameter()` etc. in the Rust bindings. These are available on specific subtypes (CylindricalSurface, SphericalSurface, Plane, etc.). For UV bounds, compute them from geometry (mesh vertices, hypothesis parameters) rather than trying to query the surface.
- **Shared vertices**: Create `TopoDS_Vertex` objects from `BRepVertex` positions (`MakeVertex::new_pnt(pt)` then `.vertex().to_owned()`), then set vertex tolerance to `vertex_tolerance_mm` via `BRep_Builder::update_vertex_vertex_real()` (OCCT default tolerance ~1e-7 is too tight for mesh vertices). Use `MakeEdge::new_handlegeomcurve_vertex2_real2(curve, &v1, &v2, fp, lp)` with explicit parameter values to avoid OCCT's vertex-to-curve projection (which fails with `Pointprojectionfailed` when vertex is slightly off the curve). Parameters are always `(first_parameter, last_parameter)` since V1 is always ordered to be the vertex closest to the curve's start point.
- **Pcurves on periodic faces**: For partial-revolution periodic faces, pcurves must be pre-set on IntSS edges using `BRep_Builder::update_edge_edge_handlegeom2dcurve_handlegeomsurface_location_real()` BEFORE creating the MakeFace. Without pre-set pcurves, MakeFace may compute incorrect pcurves that select the wrong arc of the periodic surface.
- **CRITICAL: MakeEdge vertex order for periodic curves**: `BRepBuilderAPI_MakeEdge` strips `TrimmedCurve` wrappers and uses the base curve. For periodic curves (circles), `ElCLib::AdjustPeriodic` forces p1 < p2, always selecting the forward (CCW) arc from V1 to V2. If V1/V2 are in the wrong order relative to the curve parameterization, MakeEdge creates an edge spanning the complementary arc (e.g., 270° instead of 90°). Fix: when creating MakeEdge, ensure V1 is closest to the curve's `first_parameter()` point and V2 to the `last_parameter()` point.
- **BRepCheck_Face diagnostics**: `b_rep_check::Face::new_face(face)` creates a face checker. Call `.minimum()` first to compute all checks, then query `.intersect_wires(false)`, `.classify_wires(false)`, `.orientation_of_wires(false)` (each returns `b_rep_check::Status`), and `.is_unorientable()` (bool).
- **Wire orientation fix for holes**: After adding inner wires (holes) to a `MakeFace`, apply `shape_fix::Face::new_face(make_face.face())` then `.fix_orientation()` to ensure inner wires are oriented opposite to the outer wire. Wrap the result back via `MakeFace::new_face(&fixer.face())`. Without this, OCCT's BRepCheck reports `Badorientationofsubshape` on every planar face with holes.


### Stage 3.3 edge curve computation: tricky cases
- **IntSS degenerate curves for plane-through-cylinder-axis**: When a plane passes through a cylinder's axis, `GeomAPI_IntSS` returns 4 degenerate curves with extreme parameter ranges (~-846M to -846M+47k), none containing the model's actual z-range. `ProjectPointOnCurve` fails on 3 of 4 curves. Fallback: validate IntSS curve by projecting boundary vertices — if distance exceeds 1mm, construct a `Geom_Line` from the two vertex positions instead.
- **`is_periodic()` returns false for IntSS circles**: Circle curves from cylinder×plane⊥axis intersections have `is_periodic()=false` despite being full circles with `first_param=0, last_param=2π`. Detect by checking `(last_param - first_param - 2π).abs() < 1e-6` instead.
- **Closed curve arc selection**: When a curve has param span ≈ 2π (full circle), the default parameter range may select the wrong arc (e.g., 270° instead of 90°). Fix: sample mesh boundary vertices, count which ones lie on the direct vs. complementary arc, and pick the arc with more support.
- **Centroid-based curve selection**: When IntSS returns multiple curves, selecting the one closest to a centroid computed from boundary vertices is more robust than sequential `ProjectPointOnCurve` calls, which can fail on degenerate curves.

### Stage 3.5 shell construction
- `BRepBuilderAPI_Sewing`: Create with `Sewing::new_real(tolerance)` where tolerance = `vertex_tolerance_mm`. Add all faces via `.add(face.as_shape())`, then `.perform(&ProgressRange::new())`.
- The sewed shape can contain mixed types: `Shell` (most common), `Face` (single-face models like `simple_sphere`), `Solid`, `Compound`, `Compsolid`. Iterate sub-shapes by type to extract all shells.
- For `Face` results, wrap in a `TopoDS_Shell` manually: `topo_ds::Builder::new()`, `.make_shell(&mut shell)`, `.add(shell.as_shape_mut(), face.as_shape())`.
- `ShapeAnalysis_Shell::check_oriented_shells(shape, false, false)` returns `true` when BAD edges are found (orientation problems), not when the shell is valid. The return value is inverted from what you might expect.
- `ShapeFix_Shell::new_shell(shell)`, then `.fix_face_orientation(shell, false, false)`, `.perform(&ProgressRange::new())`, `.shell()` returns the fixed shell.
- `top_abs::ShapeEnum::Compsolid` (not `CompSolid`) — watch the capitalization.
- `topo_ds::shell_shape(shape)` downcasts `&Shape` to `&Shell` reference; Shell implements `CppCopyable`/`ToOwned` so you can `.to_owned()` to get `OwnedPtr<Shell>`.
- Sewing statistics: `.nb_free_edges()`, `.nb_multiple_edges()`, `.nb_contigous_edges()` — note the OCCT typo "contigous" (not "contiguous").

### Stage 3.6 solid construction
- `shape_fix::Solid::new()`, then `.solid_from_shell(shell)` creates an `OwnedPtr<topo_ds::Solid>` with automatic orientation handling.
- `topo_ds::solid(shape)` downcasts `&Shape` to `&Solid`; `Solid` is `CppCopyable`/`ToOwned`.
- `b_rep_g_prop::volume_properties_shape_gprops_bool3(shape, &mut gprops, only_closed, skip_shared, use_triangulation)` computes volume. `gprops.mass()` returns the signed volume.
- Volume was historically wrong for models with periodic faces (cylinders/spheres) due to pcurve and vertex ordering issues. These are now resolved — all test models produce volumes matching STEP references to within ~1e-7 relative difference.
- Face orientation is handled correctly by the combination of `ShapeFix_Shell::Perform()` (edge consistency + pcurve repair) and `ShapeFix_Solid::SolidFromShell()` (global orientation via `BRepClass3d_SolidClassifier::PerformInfinitePoint`). Do NOT post-process individual face orientations — flipping individual faces breaks edge consistency, which causes both `BRepGProp::VolumeProperties` and `SolidFromShell` to give wrong results. This was extensively tested: 6/8 models produce perfect volume without any post-processing, vs. various failures when individual faces were flipped.
- `b_rep_check::Analyzer::new_shape_bool(shape, true)` + `.is_valid()` validates a solid.

### Stage 4.1 STEP output
- `step_control::Writer::new()` creates a writer. `writer.transfer_shape_stepmodeltype_bool_progressrange(shape, mode, compgraph, &progress)` transfers shapes. `writer.write(path)` writes the file.
- `step_control::StepModelType::Asis` preserves the exact geometry (no conversion to manifold/faceted form).
- `if_select::ReturnStatus::Retdone` indicates success for both transfer and write.
- `interface::Static::set_c_val(name, val)` sets STEP header metadata. Key parameters: `"write.step.schema"` (e.g., "AP214"), `"write.step.product.name"`.
- OCCT prints transfer statistics to stdout during `transfer_shape`, which cannot be suppressed via the Rust API. This is cosmetic noise in tests.
- For compare validation, re-read the written STEP via `step_control::Reader` to validate the full round-trip (write + read), not just the in-memory shapes.

### Sphere UV pole singularity handling
- Sphere surfaces have UV singularities at poles (V=±π/2) where U is undefined. When a sphere's boundary circle passes through or near a pole, UV-bounds construction fails: the projected pcurve degenerates, producing edges with zero-area faces or incorrect topology.
- **Angular proximity detection**: Use `ProjectPointOnCurve` to find the closest point on the boundary circle to each pole, then compute angular distance via chord-length formula: `cos(angle) = 1 - dist²/(2R²)`. The 45° threshold avoids both false positives and numerical issues near the exact pole.
- **Reorientation is safe for all-closed-loop tangent spheres**: When `has_tangent_closed_loops` is true (all edges are closed loops with no vertices), all tangent cylinder boundaries must share one unique axis. Proof: two tangent great circles from cylinders with different axes would intersect on the sphere in 2 points, creating vertices, which contradicts `all_closed_loops`. So reorienting the sphere to align with the (single) cylinder axis is unambiguous.
- **Failed alternatives**: (1) Wire-based split-at-poles with `ShapeConstruct_ProjectCurveOnSurface` produced correct face areas but wrong volume (cylinder only) due to BRepCheck edge/vertex consistency issues. (2) UV-bounds on the original sphere with mesh-centroid-derived U range produced 6 free edges because UV-bounds boundaries are iso-parametric meridian arcs, geometrically different from the tangent circles.

## Cylinder-Sphere Boundary: BFS Fragmentation

### face_area bug (fixed)
The `face_area` function in stage 2 only computed area from `vertex_indices[0..2]` — correct for triangles but wrong for quads (`vertex_count == 4`) created by stage 1.3 coplanar fusion. For quads, the second triangle (v0, v2, v3) was ignored, computing roughly half the correct area. Fix: check `vertex_count` and add the second triangle's area when it's 4.

### Cylinder fragmentation mechanism (empirically verified)
On rounded_cube_10_r2_fine, one of the 12 quarter-cylinder edges is fragmented into 6 small hypotheses (22+12+9+5+5+34 faces) instead of one 36-face hypothesis. Investigation found:

1. **Fragment faces are NOT sphere faces.** All 108 faces within vertex tolerance of the ideal cylinder have `sph_hyp=-1` (no sphere hypothesis) and `max_cyl_err < 2e-6` against the precise ideal cylinder axis.

2. **Fragmentation is caused by barrier faces.** At 3-4 locations along the cylinder strip, individual faces get `cyl_hyp=-1` (no valid seed partner). These faces cannot grow through committed neighbors or form valid seed pairs (cross product < MIN_CROSS_THRESHOLD with cylinder-side neighbors, sphere-side neighbors already committed by bogus equatorial hypotheses). Each barrier splits the BFS growth region.

3. **Bogus equatorial hypotheses contribute indirectly.** They commit sphere faces that are mesh neighbors of boundary cylinder faces. When BFS reaches a barrier face, its sphere-side neighbors are already committed, and its cylinder-side neighbors have nearly parallel normals. No valid seed → `cyl_hyp=-1` → barrier.

4. **The high vtx_err_max (~1.58e-4) in fragment hypotheses is misleading.** It's computed against the fragment's own poorly-fitted cylinder (fitted from a small portion of the full cylinder). Against the ideal cylinder, all face errors are < 2e-6.

5. **Coordinate precision matters.** Using truncated sphere center coordinates (4 decimal places: 21.1219 instead of 21.1218879546) shifts the ideal cylinder axis enough to make all face errors appear ~1e-5, obscuring the true ~1e-6 errors. Always use full f64 precision for distance measurements.

### Axis drift theory: disproven
Hypothesis: BFS refitting gradually rotates the cylinder axis when sphere faces are added, causing "drift." Empirical testing across rounded_cube_10_r2_{coarse,fine,medium} and pill_coarse showed **zero significant axis drift** on any successful/committed cylinder hypothesis. All axis_dot(seed, refit) values were ≥ 0.999999996. Drift only occurred on tiny (3-4 face) hypotheses already rejected by existing validation (min face count, angular coverage). The axis stability check described in some versions of the development plan would not have solved the actual problem.


### Interactive visualization (viz.rs) architecture
- `three_d::Window::render_loop(self, callback)` consumes `self` — the render loop must own the Window and run on the main thread.
- Channel-based architecture: main thread owns the render loop, pipeline runs in a background thread communicating via `mpsc::sync_channel(0)` for synchronous handoff.
- `VizSender` (pipeline side) calls `show_and_wait()` which blocks until the user presses Space/Shift+Space/Q.
- `VizReceiver` (render side) receives overlays and sends back `VizAction` responses.
- Closures inside `render_loop` must be `'static` — extract helper logic to free functions rather than using closures that capture outer scope references.
- Config has `OwnedPtr<Shape>` which is not Clone — must move Config into the pipeline thread closure, not clone it.