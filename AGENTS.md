# Instructions for LLM Coding Assistants

This file contains guidelines for AI/LLM coding assistants working on this project.

## IMPORTANT: Surviving Context Compaction

When the user says "iterate", **immediately add the full iteration checklist (steps 1-12) to your todo list before doing any work**. Include steps like "Run clippy", "Commit changes", "Export chat log", "Tell limerick" — not just implementation subtasks. The todo list survives compaction; file contents, verbal instructions, and even the user's original prompt wording do not.

When compacting, include the user's prompt in the conversation summary as guidance. Also include the hashline-edit tool usage instructions.

If you notice a `<conversation-summary>` in your context, that means compaction has occurred — **re-read this file (AGENTS.md) immediately** before continuing work.

## When instructed to iterate:

Follow this checklist IN ORDER. Do not skip steps. Do not reorder.

1. Read LEARNINGS.md, README.md, and DEVELOPMENT_PLAN.md.
2. If given instructions on work to do, follow them. Otherwise find the first item in TODO.md and implement it.
3. Don't worry about running tests for now. We're making changes such that a bunch of tests are broken.
4. Update README.md (if user-facing behavior changed) and DEVELOPMENT_PLAN.md to reflect your changes. When updating DEVELOPMENT_PLAN.md, ensure the algorithm descriptions in the Architecture section accurately describe the code as implemented — if the implementation differs from the original plan, update the plan to match reality.
5. When complete, move the item from TODO.md to the end of DONE.md.
6. Think about whether items in LEARNINGS.md are obsolete - if so, delete or rewrite them.
7. Think about whether you've learned things that may be useful to a future AI agent editing the code - if so, add them to LEARNINGS.md.
8. Commit your changes (see Commit Messages and Commit Checklist below), including any pre-existing staged files.
9. Think about whether you encountered code that would be clearer or more correct if refactored or improved. If so, make the improvements, test, and commit.
10. Run the chat log export command (see Chat Log Export below).
11. Tell the user a limerick inspired by your work this session.
12. Stop. Do not continue to the next iteration unless the user says to.

### Commit Checklist (MANDATORY — verify before every commit)

- `cargo clippy` emits no warnings.
- Build passes with no errors or warnings.
- Commit message uses heredoc format with model name and user prompts (see Commit Messages)
- Chat log exported after commit (see Chat Log Export)

## Project Layout

This is a **Rust** project. Key locations:

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

## Previous Versions

There's an abandoned incomplete C++ implementation of this tool in git branch `brepper-cpp2-abandoned` and directory ../brepper-cpp2. You're in the process of porting some improvements from the C++ implementation to this implementation.

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
2. Edit: `hashline_edit({edits: [{filePath: "src/app.ts", lineHashes: "3:fn,4:mp,5:ze", content: "export function App() {\n  return <div>Goodbye</div>;\n}"}]})`

Operations:
- **Replace**: set `lineHashes` to the lines being replaced (e.g. `"3:fn,4:mp,5:ze"`), `content` to new text
- **Replace large range**: provide at least 3 hashes at each end, e.g. `"167:or,168:dd,169:yj-479:ge,480:dd,481:ea"`
- **Insert after**: set `insertAfter: true`, `lineHashes` to anchor line
- **Insert at file start**: set `insertAfter: true`, `lineHashes` to `"0:"`
- **Delete**: set `content` to empty string
- Multiple edits can be batched in one call across files

### IMPORTANT: Don't Truncate Build or Test Output

**NEVER** pipe build or test commands through `head`, `tail`, `| grep`, or other truncating/filtering commands:
- Errors often appear at unexpected locations in the output
- Truncating can hide the actual failure while showing misleading context
- Build systems and test frameworks already produce focused error output
- If output is too long, read the full output file instead of truncating

### IMPORTANT: Time Builds and Tests

**Always** prefix build and test commands with `time`. When they take more than 2 minutes ask the user whether to speed them up.

Correct: `time make -j4 2>&1`
Correct: `time ctest --output-on-failure -j4 2>&1`
Wrong: `make -j4 2>&1 | tail -10`
Wrong: `ctest -j4 2>&1 | tail -20`

## Temporary Test Code

When writing temporary code for debugging or testing:

1. **Put files in `tmp/`** - This directory is in the project root and ignored by git. Don't use `/tmp` or other directories outside the workspace.

2. **Don't delete temporary code** - Leave it in `tmp/` in case it's useful later.

## Chat Log Export

After **every** commit, give a written summary of what you accomplished in this session, then
export the current chat session to `agent-logs/`. Use this command exactly:
`( ../vscode-export-chat-log/export-chat-log.py --wait && git add agent-logs/ && git commit --amend --no-edit ) &`

## Generating CodeCAD Test Models

Test models can be generated with ccad - see `../codecad/docs` for documentation. To regenerate the STL and STEP files for CodeCAD-based tests:

```bash
# Requires 'ccad' tool installed (https://codecad.xyz)
./tests/ccad/generate_models.sh
```