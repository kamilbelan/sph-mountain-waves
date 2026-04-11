#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <experiment-folder-name> [submit.sh location]"
    echo ""
    echo "Examples:"
    echo "  $0 2026_09_04_diffusion_long"
    echo "  $0 2026_09_04_diffusion_long submit_batch_snehurka.sh"
    exit 1
fi

exp_name="$1"
submit_script="${2:-submit_batch_snehurka.sh}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exp_dir="$SCRIPT_DIR/experiments/$exp_name"

# check that experiment directory exists
if [[ ! -d "$exp_dir" ]]; then
    echo "Error: Experiment directory not found: $exp_dir"
    exit 1
fi

# check that config files exist
global_conf="$exp_dir/global_params.toml"
sim_conf="$exp_dir/sim_params.toml"

if [[ ! -f "$global_conf" ]]; then
    echo "Error: global_params.toml not found at $global_conf"
    exit 1
fi

if [[ ! -f "$sim_conf" ]]; then
    echo "Error: sim_params.toml not found at $sim_conf"
    exit 1
fi

# submit the job
echo "Submitting experiment: $exp_name"
echo "  Global config: $global_conf"
echo "  Sim config:    $sim_conf"
echo ""

cd "$SCRIPT_DIR"
sbatch "$submit_script" "$global_conf" "$sim_conf"
