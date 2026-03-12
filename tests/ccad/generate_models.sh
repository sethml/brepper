#!/bin/bash
set -e

# Get the directory of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "Generating CodeCAD models..."
if command -v ccad &> /dev/null; then
    ccad build
    # Normalize timestamps in generated STEP files so they don't change on every regeneration.
    # The FILE_NAME line contains a timestamp like '2026-03-11T20:08:38'.
    for f in "$DIR"/generated/*.step; do
        sed -i '' "s/'[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}T[0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}'/'2000-01-01T00:00:00'/" "$f"
    done
    echo "Models generated in $DIR/generated/"
else
    echo "Error: ccad tool not found. Please install CodeCAD to generate models."
    exit 1
fi
