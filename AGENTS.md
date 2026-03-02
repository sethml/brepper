# Instructions for LLM Coding Assistants

This file contains guidelines for AI/LLM coding assistants working on this project.

## When instructed to iterate:

- Read LEARNINGS.md.
- Read README.md.
- Read DEVELOPMENT_PLAN.md.
- Find the next incomplete item in the "## Implementation Phases" section of DEVELOPMENT_PLAN.md.
- Implement that item. Make sure tests pass. Don't skip important functionality. If there's enough ambiguity in the specification that you're unsure how to proceed, ask the user what to do. Don't be afraid to modify existing data structures or code if it's helpful.
- Update README.md and DEVELOPMENT_PLAN to refect your changes.
- Think about whether items in LEARNINGS.md are obsolete - if so, delete or rewrite them.
- Think about whether you've learned things that may be useful to a future AI agent editing the code - if so, add them to LEARNINGS.md.
- Commit your changes, including any pre-existing staged files.
- Think about whether you encountered code that would be clearer or more correct if refactored or improved. If so, make the improvements, test, and commit.
- Tell the user a limmerick inspired by your work this session.

## Project Layout

This is a **Rust** project. Key locations:

- `README.md` - Program usage. Read it.
- `DEVELOPMENT_PLAN.md` - Development plan. Read it.
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

## Chat Log Export

After committing, give a written summary of what you accomplished in this session, then
export the current chat session to `agent-logs/`. Use this command exactly:
`( ../vscode-export-chat-log/export-chat-log.py --wait && git add agent-logs/ && git commit --amend --no-edit ) &`

## Generating CodeCAD Test Models

To regenerate the STL and STEP files for CodeCAD-based tests:

```bash
# Requires 'ccad' tool installed (https://codecad.xyz)
./tests/ccad/generate_models.sh
```