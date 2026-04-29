#!/bin/bash

# List your frame numbers here (e.g., 57 58 90 60)
FRAMES=(64)

for XX in "${FRAMES[@]}"; do
    # Format the number to 4 digits (e.g., 57 -> 0057)
    FRAME_FILE=$(printf "frame_%04d.h5" $XX)
    
    echo "------------------------------------------------"
    echo "Processing SPH: $FRAME_FILE"
    
    julia -t 4 --project=. scripts/sph_tensile_instability.jl \
        data/final/evolutionary/SPH/ "$FRAME_FILE" SPH
done
