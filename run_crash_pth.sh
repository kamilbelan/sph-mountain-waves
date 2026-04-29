#!/bin/bash

# List your frame numbers here
FRAMES=(55 56 57 58 59 60 61 62)

for XX in "${FRAMES[@]}"; do
    FRAME_FILE=$(printf "frame_%04d.h5" $XX)
    
    echo "------------------------------------------------"
    echo "Processing PTH: $FRAME_FILE"
    
    julia -t 4 --project=. scripts/pth_tensile_instability.jl \
        data/final/evolutionary/PTH/ "$FRAME_FILE" PTH
done
