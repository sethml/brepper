#!/bin/bash
set -e
# Regenerate OCCT-generated test models (tests/occt/). These models cannot be
# produced via ccad (e.g. multi-shell solids with internal voids).
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
cd "$DIR"
cargo build --release --bin gen-test-models
./target/release/gen-test-models tests/occt
# Normalize STEP timestamps so regeneration doesn't create spurious diffs.
for f in tests/occt/*.step; do
    sed -i '' "s/'[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}T[0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}[^']*'/'2000-01-01T00:00:00'/" "$f"
done
echo "Models generated in tests/occt/"
