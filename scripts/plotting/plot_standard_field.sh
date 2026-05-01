#!/bin/bash

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# Target the standard directory specifically
TARGET_DIR="$REPO_ROOT/data/final/stationary/standard"

echo "Starting automated plotting for Standard schemes..."
echo "==================================================="

# Loop through every subdirectory (PA, PTH, SPH) inside standard/
for run_dir in "$TARGET_DIR"/*; do
    # Check if it is actually a directory
    if [ -d "$run_dir" ]; then
        # Extract the variant name (e.g., "PA") from the path
        variant=$(basename "$run_dir")
        
        # Define the output prefix (e.g., "standard_PA")
        out_name="standard_${variant}"
        
        echo "Processing: $variant"
        echo "Input:  $run_dir"
        echo "Output: $out_name"
        
        # Run the Julia script
        julia -t 4 --project="$REPO_ROOT" "$REPO_ROOT/scripts/plotting/stationary_field.jl" "$run_dir" "$out_name"
        
        echo "---------------------------------------------------"
    fi
done

echo "Batch processing complete!"
