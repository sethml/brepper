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
