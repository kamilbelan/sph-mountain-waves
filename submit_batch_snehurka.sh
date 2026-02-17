#!/bin/bash -l
#SBATCH --job-name=SPH_static
#SBATCH --output=logs/SPH_static-%j.out
#SBATCH --error=logs/SPH_static-%j.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=3:00:00
#SBATCH --partition=express3

#SBATCH --mail-type=ALL
#SBATCH --mail-user=slurm@kamilbelan.anonaddy.com


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

echo "=== JOB START $(date) on $(hostname) ==="
echo "   Global Config: $GLOBAL_CONF"
echo "   Sweep Config:  $SIM_CONF"

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

# run the script, passing the config files as arguments
$JULIA_BIN --project=. scripts/run_sim.jl "$GLOBAL_CONF" "$SIM_CONF"

echo "=== JOB END $(date) ==="
