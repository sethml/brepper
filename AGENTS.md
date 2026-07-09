# Instructions for LLM Coding Assistants

This file contains guidelines for AI/LLM coding assistants working on this project. Read ALL of AGENTS.md - no just the first 100 lines.

### Commit Checklist (MANDATORY — verify before every commit)

- `cargo clippy` emits no warnings.
- Build passes with no errors or warnings.

## Project Layout

This is a Rust project. Key locations:

- `README.md` - Program usage. Read it.
- `DEVELOPMENT_PLAN.md` - Development plan. Read it.
- `Cargo.toml` — Rust project manifest
- `src/` — Rust source files
  - `bin/brepper.rs` — the main program
  - `bin/step_measure_tolerance.rs` — utility to measure distances from STL to correpsonding STEP
  - `bin/viewer.rs` — GUI to visualize STL or STEP files
- `scripts/` — shell scripts for batch operations over test data
  - `step_distance_table.sh` - runs `step_measure_tolerance` over all STL/STEP pairs and summarizes statistics
- `tests/` — STL and STEP test file pairs
  - `manual/` — handcrafted simple models
  - `onshape/` — models exported from Onshape
  - `ccad/generated/` — models generated via CodeCAD

### OCCT / opencascade-rs

The project uses the `opencascade-sys` crate (auto-generated Rust FFI bindings to OpenCASCADE Technology) via a path dependency. MANDATORY: All useful C++ OCCT APIs should be available in the binding - if you want to use one that doesn't seem to be bound, complain loudly about it AND STOP.

```
../opencascade-rs/crates/opencascade-sys   ← Rust FFI bindings (8700+ types)
../opencascade-rs/crates/opencascade-sys/PORTING.md  ← Naming conventions & usage patterns
../opencascade-rs/crates/opencascade-sys/generated/  ← All generated .rs modules
../opencascade-rs/crates/occt-sys/OCCT/src/  ← OCCT C++ source (for reading docs/headers)
../opencascade-rs/crates/occt-sys/OCCT/dox/  ← OCCT C++ documentation
```

The `builtin` feature (default) compiles OCCT from source — first build takes ~8 minutes. Subsequent builds are fast.

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

### IMPORTANT: Time Builds and Tests

**Always** prefix build and test commands with `time`. When they take more than 2 minutes ask the user whether to speed them up.

Correct: `time make -j4 2>&1`
Correct: `time ctest --output-on-failure -j4 2>&1`

## Temporary Test Code

When writing temporary code for debugging or testing:

1. **Put files in `tmp/`** - This directory is in the project root and ignored by git. Don't use `/tmp` or other directories outside the workspace.

2. **Don't delete temporary code** - Leave it in `tmp/` in case it's useful later.

## Library Reuse

Before implementing non-trivial functionality (text rendering, file parsing, data structures, etc.), search for existing open-source libraries that provide that functionality. Evaluate whether using a library will be better for code simplicity and maintainability than writing it from scratch. Consider factors like: compatibility with existing dependencies, maintenance status, API complexity, and how well the library fits the use case. Then go with your decision — either use the library or implement it yourself — and briefly explain why.

## Generating CodeCAD Test Models

Test models can be generated with ccad - see `../codecad/docs` for documentation. To regenerate the STL and STEP files for CodeCAD-based tests:

```bash
# Requires 'ccad' tool installed (https://codecad.xyz)
./tests/ccad/generate_models.sh
```
