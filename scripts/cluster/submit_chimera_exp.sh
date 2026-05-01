#!/usr/bin/env bash
# Usage: ./submit_chimera_exp.sh experiments/2026_13_04_sponge_med/
#    or: ./submit_chimera_exp.sh 2026_13_04_sponge_med
#
# If the experiment directory contains a .sh file, uses that as the sbatch script.
# Otherwise falls back to the root submit_chimera.sh.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <experiment_dir>"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# accept both "experiments/foo" and just "foo"
input="${1%/}"
if [[ -d "$PROJECT_ROOT/experiments/$input" ]]; then
    exp_dir="$PROJECT_ROOT/experiments/$input"
elif [[ -d "$PROJECT_ROOT/$input" ]]; then
    exp_dir="$PROJECT_ROOT/$input"
else
    echo "Error: Experiment directory not found: $input"
    exit 1
fi

# validate config files
if [[ ! -f "$exp_dir/global_params.toml" ]] || [[ ! -f "$exp_dir/sim_params.toml" ]]; then
    echo "Error: $exp_dir must contain global_params.toml and sim_params.toml"
    exit 1
fi

# find a submit script: use the experiment's .sh if it exists, otherwise fall back to root
SUBMIT_SH=$(find "$exp_dir" -maxdepth 1 -name "*.sh" | head -1)

if [[ -n "$SUBMIT_SH" ]]; then
    echo "Using experiment submit script: $SUBMIT_SH"
else
    SUBMIT_SH="$SCRIPT_DIR/submit_chimera.sh"
    echo "No .sh in experiment dir, falling back to: $SUBMIT_SH"
fi

cd "$PROJECT_ROOT"

# convert to relative paths for SLURM
global_rel="${exp_dir#$PROJECT_ROOT/}/global_params.toml"
sim_rel="${exp_dir#$PROJECT_ROOT/}/sim_params.toml"

echo "Submitting: sbatch $SUBMIT_SH $global_rel $sim_rel"
sbatch "$SUBMIT_SH" "$global_rel" "$sim_rel"

