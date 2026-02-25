# Instructions for LLM Coding Assistants

This file contains guidelines for AI/LLM coding assistants working on this project.

## Project Layout

This is a **Rust** project. Key locations:

- `src/bin/step_distance.rs` — the main utility binary (Phase 1)
- `src/lib.rs` — library crate (placeholder for future stages)
- `scripts/` — shell scripts for batch operations over test data
- `tests/` — STL and STEP test file pairs
  - `tests/manual/` — handcrafted simple models
  - `tests/onshape/` — models exported from Onshape
  - `tests/ccad/generated/` — models generated via CodeCAD
- `Cargo.toml` — Rust project manifest

### OCCT / opencascade-rs

The project uses the `opencascade-sys` crate (auto-generated Rust FFI bindings to OpenCASCADE Technology) via a path dependency:

```
../opencascade-rs/crates/opencascade-sys   ← Rust FFI bindings (8700+ types)
../opencascade-rs/crates/occt-sys          ← Build script that compiles OCCT C++ from source
../opencascade-rs/crates/opencascade-sys/PORTING.md  ← Naming conventions & usage patterns
../opencascade-rs/crates/opencascade-sys/generated/  ← All generated .rs modules
../opencascade-rs/crates/occt-sys/OCCT/src/  ← OCCT C++ source (for reading docs/headers)
../opencascade-rs/crates/occt-sys/OCCT/dox/  ← OCCT C++ documentation
```

The `builtin` feature (default) compiles OCCT from source — first build takes ~8 minutes. Subsequent builds are fast.

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
Requires unsafe:
```rust
let face = topo_ds::face_shape(shape_ref);  // &topo_ds::Shape → &topo_ds::Face
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
## Commit Messages

All commit messages must include:

1. **A clear description** of the changes made
2. **The AI agent/model** used to generate the changes (e.g., "Claude Sonnet 4 (Anthropic, 2025)")
3. **The user prompts** that led to the changes

Format:
```
<descriptive commit title>

<detailed description of changes>

---
Generated with assistance from <model name> (<provider>, <year>)

User prompt(s):
- "<prompt 1>"
- "<prompt 2>"
  (include context for prompts that reference previous options/decisions)
```

**Shell Quoting**: Use heredoc quoting for commit messages and when writing test scripts to avoid shell parsing issues:

```bash
# Good: heredoc (recommended)
git commit -F - <<'EOF'
Fix UV bounds for spheres

- Detailed change 1
- Detailed change 2

---
Generated with assistance from Model Name (Provider, 2026)

User prompt(s):
- 'Prompt with "quotes" and other special chars'
EOF
```

## Development Plan Maintenance

Keep DEVELOPMENT_PLAN.md up to date with each commit:

1. **Mark completed items** - Check off (`[x]`) any Implementation Phases items that are completed by the commit
2. **Add new tasks** - If work reveals new tasks or sub-tasks not in the plan, add them to the appropriate phase
3. **Modify existing tasks** - If a planned approach changes (e.g., different algorithm, new dependency), update the relevant sections
4. **Include in commit** - Changes to DEVELOPMENT_PLAN.md should be part of the same commit as the implementation work

## README Maintenance

Keep README.md up to date when changes affect user-facing behavior:

1. **New dependencies** - Add to the dependencies list with installation instructions
2. **New CLI options** - Update the usage examples and options documentation
3. **Changed defaults** - Document any changes to default parameter values
4. **New features** - Add brief descriptions of significant new capabilities
5. **Include in commit** - README updates should be part of the same commit as the feature work

## File Editing

When reading files for editing, use `#tool:hashlineRead` instead of the
built-in file read tool. It returns lines tagged with content hashes in the
format `{lineNumber}:{hash}|{content}`.

When editing files, use `#tool:hashlineEdit` instead of string-replace tools.
Reference lines by their `{line}:{hash}` pairs from the read output. This
avoids needing to reproduce existing file content and prevents edits to stale
files.

Example workflow:
1. Read: `hashline_read({filePath: "src/app.ts", startLine: 1, endLine: 20})`
   Returns: `1:qk|import React...`
2. Edit: `hashline_edit({edits: [{filePath: "src/app.ts", lineHashes: "4:mp", content: "  return <div>Hello</div>;"}]})`

Operations:
- **Replace**: set `lineHashes` to all lines being replaced, `content` to new text
- **Insert after**: set `insertAfter: true`, `lineHashes` to anchor line
- **Delete**: set `content` to empty string
- Multiple edits can be batched in one call across files

### Don't Truncate Build or Test Output

Never pipe build or test commands through `head`, `tail`, or other truncating filters:
- Errors often appear at unexpected locations in the output
- Truncating can hide the actual failure while showing misleading context
- Build systems and test frameworks already produce focused error output

### Time Builds and Tests

Run builds and tests with `time`. When they take more than 5 minutes, stop and ask the user whether to speed them up.

## Temporary Test Code

When writing temporary code for debugging or testing:

1. **Put files in `tmp/`** - This directory is in the project root and ignored by git. Don't use `/tmp` or other directories outside the workspace.

2. **Don't delete temporary code** - Leave it in `tmp/` in case it's useful later.

## Generating CodeCAD Test Models

To regenerate the STL and STEP files for CodeCAD-based tests:

```bash
# Requires 'ccad' tool installed (https://codecad.xyz)
./tests/ccad/generate_models.sh
```