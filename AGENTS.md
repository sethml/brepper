# Instructions for LLM Coding Assistants

This file contains guidelines for AI/LLM coding assistants working on this project.

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

## API Documentation

When implementing features that use external libraries:

1. **Fetch official documentation first** - Use web fetch tools to retrieve current API documentation from official sources:
   - PCL: https://pointclouds.org/documentation/tutorials/
   - OpenCASCADE: https://dev.opencascade.org/doc/overview/html/
   - Eigen: https://eigen.tuxfamily.org/dox/

2. **Consult header files if needed** - If documentation is insufficient or unclear, read the relevant header files from the installed libraries to understand:
   - Function signatures and parameter types
   - Available overloads
   - Template parameters
   - Expected usage patterns

3. **Avoid assumptions** - Do not guess at API details. When uncertain, gather more context before implementing.

## Project-Specific Context

- This project uses PCL 1.15+, OpenCASCADE 7.9+, and Eigen 3.4+
- Build system is CMake 3.20+
- Target platform is macOS with Homebrew-installed dependencies
- See DEVELOPMENT_PLAN.md for the overall architecture and implementation stages

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

## Code Style Preferences

### Avoid Large Conditional Compilation Blocks

Prefer making dependencies required rather than optional if it would require large `#ifdef` blocks:

- Large conditional blocks lead to hidden build breakages and bugs
- Code paths that aren't regularly compiled tend to bit-rot
- Testing burden doubles when code has multiple compilation configurations
- If a dependency is important enough to use, make it required

**Bad:**
```cpp
#ifdef USE_FEATURE_X
    // 50+ lines of code using feature X
#else
    // 50+ lines of fallback code
#endif
```

**Better:** Make the dependency required, or isolate the feature into a separate optional component that can be tested independently.

Small `#ifdef` blocks (e.g., platform-specific includes, debug logging) are acceptable.

### Command Line Efficiency

When building and running tests, combine commands with `&&` to avoid unnecessary waits:

```bash
# Good: single command, stops on failure, parallel tests
cmake --build build && ctest --test-dir build -j8 --output-on-failure

# Bad: separate commands requiring multiple tool invocations
cmake --build build
ctest --test-dir build --output-on-failure
```

### Run Tests in Parallel

Always run tests with `-j8` (or similar) for parallel execution:

```bash
# Good: parallel test execution (~21 seconds)
ctest --test-dir build -j8 --output-on-failure

# Bad: sequential execution (~98 seconds)
ctest --test-dir build --output-on-failure
```

The test suite is structured with separate `TEST_CASE`s (rather than `SECTION`s) specifically to enable parallel execution via ctest.

### Sanitizer Testing

Periodically run tests with sanitizers enabled to catch memory and concurrency issues:

```bash
# Build with Address Sanitizer + Undefined Behavior Sanitizer
cmake -B build-sanitize -DBUILD_TESTS=ON -DENABLE_ASAN=ON -DENABLE_UBSAN=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-sanitize && ctest --test-dir build-sanitize -j4 --output-on-failure

# Build with Thread Sanitizer (for OpenMP race conditions)
cmake -B build-tsan -DBUILD_TESTS=ON -DENABLE_TSAN=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-tsan && ctest --test-dir build-tsan -j4 --output-on-failure
```

Note: Use `-j4` instead of `-j8` with sanitizers as they increase memory usage significantly.

### Don't Truncate Build or Test Output

Never pipe build or test commands through `head`, `tail`, or other truncating filters:

- Errors often appear at unexpected locations in the output
- Truncating can hide the actual failure while showing misleading context
- Build systems and test frameworks already produce focused error output

```bash
# Bad: might hide the actual error
cmake --build build 2>&1 | head -100

# Good: see all output
cmake --build build
```
