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


### GeomAPI_ProjectPointOnSurf projects onto unbounded surfaces
`GeomAPI_ProjectPointOnSurf` projects onto the *infinite* underlying surface, not the bounded face. When testing distance from a point to a STEP model, a point that is far from a face's bounded area but coincidentally lies on the infinite plane of another face will report distance 0. This means `--compare` vertex-distance checks can silently pass for points that are geometrically far from the actual solid. When creating test files for surface comparison failures, shift vertices away from ALL infinite surface planes, not just in one axis-aligned direction.
