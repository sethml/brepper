#!/bin/bash
set -e
# Regenerate OCCT-generated test models: tests/occt/ (models ccad can't make,
# e.g. multi-shell voids) and tests/nurbs/ (exact freeform models targeting
# the DEVELOPMENT_PLAN §2.7 NURBS fallback).
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
cd "$DIR"
cargo build --release --bin gen-test-models
./target/release/gen-test-models tests/occt tests/nurbs
# Normalize STEP timestamps so regeneration doesn't create spurious diffs.
for f in tests/occt/*.step tests/nurbs/*.step; do
    sed -i '' "s/'[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}T[0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}[^']*'/'2000-01-01T00:00:00'/" "$f"
done
echo "Models generated in tests/occt/ and tests/nurbs/"
