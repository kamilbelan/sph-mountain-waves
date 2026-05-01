#!/bin/bash

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
VARIANTS=("SPH" "PTH")
BASE_DIR="$REPO_ROOT/data/final/stationary"

echo "Starting automated L_inf / L_2 norm plotting..."
echo "==============================================="

for var in "${VARIANTS[@]}"; do
    STD_DIR="${BASE_DIR}/standard/${var}"
    WB_DIR="${BASE_DIR}/well-balanced/${var}"
    OUT_NAME="stationary_norms_${var}"
    
    if [ -d "$STD_DIR" ] && [ -d "$WB_DIR" ]; then
        echo "Processing formulation: $var"
        julia -t 4 --project="$REPO_ROOT" "$REPO_ROOT/scripts/plotting/stationary_norms.jl" "$STD_DIR" "$WB_DIR" "$OUT_NAME"
        echo "-----------------------------------------------"
    else
        echo "Warning: Missing directories for $var. Skipping."
        echo "-----------------------------------------------"
    fi
done

echo "All norm plots generated successfully!"
