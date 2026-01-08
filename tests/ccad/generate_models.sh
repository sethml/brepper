#!/bin/bash
set -e

# Get the directory of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "Generating CodeCAD models..."
if command -v ccad &> /dev/null; then
    ccad build
    echo "Models generated in $DIR/generated/"
else
    echo "Error: ccad tool not found. Please install CodeCAD to generate models."
    exit 1
fi
