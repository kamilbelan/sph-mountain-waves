#!/bin/bash -l
#SBATCH --job-name=SPH
#SBATCH --output=logs/SPH-%j.out
#SBATCH --error=logs/SPH-%j.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --partition=express3

#SBATCH --mail-type=ALL
#SBATCH --mail-user=slurm@kamilbelan.anonaddy.com

set -euo pipefail

# ==============================================================================
# USAGE: sbatch submit_batch.sh config/global_params.toml configs/sim_params.toml 
# ==============================================================================
# MUST BE RUN FROM THE PROJECT ROOT! WHERE PROJECT.TOML LIVES


# ==============================================================================
# INPUT 
# ==============================================================================
# check that two config files have been passed as arguments to sbatch

if [ "$#" -ne 2 ]; then
    echo "Error: You must provide two config files."
    echo "Usage: sbatch submit_batch.sh <global.toml> <sweep.toml>"
    exit 1
fi

GLOBAL_CONF=$1
SIM_CONF=$2

# derive a meaningful job name from the experiment folder (e.g. "2026_01_09_highres") and update the job name visible in squeue
JOB_NAME=$(basename "$(dirname "$GLOBAL_CONF")")
scontrol update JobId=$SLURM_JOB_ID Name="$JOB_NAME" 2>/dev/null || true

echo "=== JOB START $(date) on $(hostname) ==="
echo "   Job Name:      $JOB_NAME"
echo "   Global Config: $GLOBAL_CONF"
echo "   Sweep Config:  $SIM_CONF"
echo "   Git commit:    $(git rev-parse --short HEAD)"

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

# set threads to match SLURM request
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "   Running with $JULIA_NUM_THREADS threads"

# navigate to the directory where you submitted the job (ie project root)
# this ensures DrWatson finds Project.toml correctly.
cd $SLURM_SUBMIT_DIR
echo "   Working Directory: $(pwd)"

# define julia binary path
JULIA_BIN=/usr/work/belank/software/julia-1.12.4/bin/julia

# ==============================================================================
# EXECUTION
# ==============================================================================

mkdir -p logs

# run the script, passing the config files as arguments
$JULIA_BIN --sysimage=sph.so --project=. scripts/run_sim.jl "$GLOBAL_CONF" "$SIM_CONF"

echo "=== JOB END $(date) ==="
