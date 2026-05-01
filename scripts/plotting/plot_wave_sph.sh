#!/bin/bash

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

# List your frame numbers here (e.g., 10 20 30)
FRAMES=(30 35 40 45 50 55)

for XX in "${FRAMES[@]}"; do
    FRAME_FILE=$(printf "frame_%04d.h5" $XX)

    echo "------------------------------------------------"
    echo "Processing SPH: $FRAME_FILE"

    julia -t 4 --project="$REPO_ROOT" "$REPO_ROOT/scripts/plotting/emerging_wave.jl" \
        "$REPO_ROOT/data/final/evolutionary/SPH/" "$FRAME_FILE" SPH
done
