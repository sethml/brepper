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
- 3×3 symmetric eigenvalue computation uses depressed cubic (trigonometric solution with cos(θ/3)), plus null-space eigenvector via cross products of (M - λI) rows.
- BFS seeds from pairs of adjacent faces with non-parallel normals (cross product > MIN_CROSS_THRESHOLD=0.01). Multi-face planar faces are excluded from seeding (both as seed face and seed partner) to prevent bogus large-radius cylinders from nearly-coplanar face pairs on curved surfaces. However, multi-face planar faces remain UNDEDUCED and CAN be absorbed by BFS expansion from nearby legitimate seeds.
- Centroid validation after BFS: each face centroid must be within surface_tolerance of the fitted cylinder. This rejects spurious fits from perpendicular planar faces that algebraically fit a cylinder.
- Minimum 3-face requirement: cylinder hypotheses must contain at least 3 faces. Any real cylindrical surface from a CAD tessellation will have at least 3 facets around its circumference. This eliminates spurious 2-face cylinder fits from locally adjacent faces on tori, cones, and spheres (e.g., cone went from 29 to 0 spurious fits, pipe_elbow from 252 to 51).
- Angular tolerance check: both seed pairs and BFS expansion candidates must have dihedral angle ≤ `--angular-tolerance` (default 17.5°). This prevents fitting cylinders to faces that meet at sharp angles (e.g., cube faces at 90°). CAD tessellators enforce a maximum angular tolerance between adjacent triangles on the same surface.

### ccad box() coordinate system
- `box(w, h, d)` places the box in the positive octant: (0,0,0) to (w,h,d). It does NOT center at origin.
- To center a box, use `translate(box(w,h,d), -w/2, -h/2, -d/2)`.
- `cylinder(diameter, height)` is vertical (Z-axis), centered at the XY origin, extends from z=0 to z=height (NOT centered vertically). To make through-holes, translate cylinders vertically so they extend past both sides of the target body.
- `sphere(diameter)` is centered at the origin. To make hemispherical pockets or domes, use `intersection()` or `difference()` with a box.

### Stage 2.3 spherical hypothesis algorithm
- Algebraic least-squares sphere fitting: expand |v-c|²=r² linearly into 4 unknowns (cx, cy, cz, k=r²-|c|²). Normal equations AᵀAx = Aᵀb solved via 4×4 Gaussian elimination with partial pivoting.
- BFS seeds from pairs of adjacent faces with non-parallel normals (cross product > MIN_CROSS_THRESHOLD=0.01), analogous to cylindrical seeding. Multi-face planar faces excluded from seeding but absorb-able by BFS, same as cylindrical.
- Original plan called for 3-face seeds with normals-span-3D eigenvalue check, but this failed on finely tessellated spheres where 3 adjacent faces span only ~6° and fall below the eigenvalue ratio threshold. Pair-based seeding works because the 4-DOF sphere fit from 6+ vertices validates the geometry.
- `MAX_SPHERE_RADIUS_FACTOR = 2.0`: max_sphere_radius = bounding_box_diagonal × 2. Previously was 10.0, which allowed bogus r≈139mm sphere hypotheses from small arcs on cylinder fillet faces of a 17.3mm-diagonal model.
- `MAX_CYLINDER_RADIUS_FACTOR = 2.0`: analogous cap for cylinder hypotheses. Radius checked at seed acceptance, BFS refit, and final validation.
- Centroid validation after BFS: each face centroid must be within surface_tolerance of the fitted sphere.
- Minimum 4-face requirement: a sphere needs at least 4 non-coplanar points for a unique fit.
- All faces are candidates for BFS expansion (including those with cylindrical hypotheses and multi-face planar hypotheses). Only seeding excludes multi-face planar. Stage 2.6 resolves overlapping hypotheses using per-face priority.
- Torus faces locally fit spheres (e.g., pipe_elbow produces 36 sphere hypotheses). This is expected — stage 2.6 with torus hypothesis support will resolve it.
- Concave spheres (bowls/pockets) have normals pointing toward center — same convexity logic as cylinders.
- Angular tolerance check: same as cylindrical — seed pairs and BFS neighbors must have dihedral angle ≤ `--angular-tolerance`. This fixed the spurious sphere detection on the 1mm manual cube.

### Stage 2.6 surface selection
- Greedy area-based selection: iteratively selects the hypothesis with the largest total remaining mesh face area, assigns its unassigned faces, removes those faces from all other candidates, repeats.
- Multi-face planar, cylindrical, and spherical hypotheses are candidates. Single-face planar hypotheses are fallback for unclaimed faces.
- After selection, hypothesis face/vertex lists are updated in-place to reflect only the faces actually assigned.
- Replaced the old per-face priority rule (spherical > cylindrical > multi-face planar > single-face planar) which failed when bogus small-area hypotheses of a high-priority type beat correct large-area hypotheses of a lower-priority type.
- Sphere-first pipeline (stage >= 2.6): sphere BFS runs before cylinder BFS. Cylinder BFS skips sphere-assigned faces and uses sagitta-based `expansion_tol = max(radius * (1 - cos(angular_tol/2)), vertex_tol)` instead of tight vertex_tol. This allows cylinder BFS to absorb fillet faces that would previously be left as planar fallbacks (e.g., rounded_cube: 36 faces/cylinder vs ~7 before).
- Without sphere-first, cylinder BFS with loose expansion_tol absorbs sphere faces. Without loose expansion_tol, cylinder BFS only captures faces within 10nm of the fit cylinder, which is too strict for sagitta on coarse tessellations.
- Stage gate: sphere-first only at stage >= 2.6, NOT 2.3. At stage 2.3, a cylinder viewed along its axis fits a sphere, so running sphere BFS first absorbs bore-hole faces as false spheres.
- Compare tolerance uses `surface_tolerance` (0.4mm) for all hypothesis types.
- Merge phase was removed from cylinder BFS — it caused cumulative axis drift when merging multiple cylinder hypotheses.


### Stage 3.1 adjacency graph construction
- OCCT surface creation: `geom::Plane::new_pnt_dir`, `geom::CylindricalSurface::new_ax3_real`, `geom::SphericalSurface::new_ax3_real`. Type erasure via `.to_handle().to_handle_surface()` to get `OwnedPtr<HandleGeomSurface>`.
- Adjacency is determined from mesh edges: an edge between mesh faces assigned to different selected surfaces is a boundary edge.
- Boundary edges are grouped by `SurfacePair` (canonical ordering), then chained into ordered vertex sequences. Each chain becomes a `ReconEdge`.
- Corner vertices (where 3+ surfaces meet) become `BRepVertex` entries.
- Topological ordering of edges around each face uses rotation-order walk: at each vertex, pick the edge to the face that comes next in counter-clockwise order. For faces without vertices (closed loops like cylinder caps), all edges incident on that face are included.
- Bug encountered: closing vertex in topological walk was missing — after walking N edges around a face, the vertex connecting the last edge back to the first edge must be added separately.
- `wedge(10, 10, 10, 5)` in CodeCAD is a truncated wedge (8 vertices, 6 faces), NOT a triangular prism (6 vertices, 5 faces). The 4th parameter `ltx` is the top face X dimension.
- Euler formula V-E+F=2 is a good sanity check for genus-0 (no holes through the model) manifold solids. Models with through-holes have different Euler characteristics.