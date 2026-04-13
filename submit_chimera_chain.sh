#!/usr/bin/env bash
# ==============================================================================
# Chain multiple 12h ffa jobs to achieve long runs (e.g. 72h) without preemption.
#
# USAGE:
#   ./submit_chimera_chain.sh <N_JOBS> <global.toml> <sim.toml> [restart_dir]
#
# EXAMPLES:
#   # Fresh 72h run (6 × 12h jobs):
#   ./submit_chimera_chain.sh 6 experiments/my_exp/global_params.toml experiments/my_exp/sim_params.toml
#
#   # Resume an interrupted chain from an existing run directory:
#   ./submit_chimera_chain.sh 4 experiments/my_exp/global_params.toml experiments/my_exp/sim_params.toml data/sims/EvolWBalancedThetaHopkins/2026-04-13/123456_dr=520/
#
# HOW IT WORKS:
#   Job 1 runs fresh (or resumes from restart_dir if given).
#   Jobs 2..N wait for the previous job via --dependency=afterok.
#   All jobs write to the SAME run directory, continuing from the latest checkpoint.
#
#   The first job creates the run directory and writes its path to a shared file.
#   Subsequent jobs read that file to find the restart_dir.
#
# NOTES:
#   - t_end in sim_params.toml is the TOTAL simulation time, not per-job.
#     Each job runs the time loop from its checkpoint to t_end and checkpoints
#     when the 12h wall time forces SLURM to stop it.
#   - If a job finishes t_end before the wall time, the remaining chained
#     jobs will detect that and exit immediately (no wasted queue time).
# ==============================================================================

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <N_JOBS> <global.toml> <sim.toml> [restart_dir]"
    exit 1
fi

N_JOBS=$1
GLOBAL_CONF=$2
SIM_CONF=$3
RESTART_DIR="${4:-}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SUBMIT_SCRIPT="$SCRIPT_DIR/submit_chimera.sh"

# shared file where job 1 writes the run_dir for subsequent jobs to pick up
# must be on a shared filesystem (not /tmp, which is node-local)
mkdir -p "$SCRIPT_DIR/logs"
CHAIN_FILE="$SCRIPT_DIR/logs/chain_$(date +%Y%m%d_%H%M%S)_$$.txt"
echo "$RESTART_DIR" > "$CHAIN_FILE"

echo "=== Submitting chain of $N_JOBS jobs ==="
echo "    Config: $GLOBAL_CONF / $SIM_CONF"
if [[ -n "$RESTART_DIR" ]]; then
    echo "    Resuming from: $RESTART_DIR"
fi
echo "    Chain file: $CHAIN_FILE"
echo ""

PREV_JOB=""
for i in $(seq 1 "$N_JOBS"); do
    # build sbatch options
    SBATCH_OPTS=""
    if [[ -n "$PREV_JOB" ]]; then
        SBATCH_OPTS="--dependency=afterany:$PREV_JOB"
    fi

    # export CHAIN_FILE so the wrapper can find the restart dir
    JOB_ID=$(sbatch --parsable \
        $SBATCH_OPTS \
        --export=ALL,CHAIN_FILE="$CHAIN_FILE" \
        "$SUBMIT_SCRIPT" "$GLOBAL_CONF" "$SIM_CONF")

    echo "  Job $i/$N_JOBS: $JOB_ID $([ -n "$PREV_JOB" ] && echo "(after $PREV_JOB)" || echo "(first)")"
    PREV_JOB=$JOB_ID
done

echo ""
echo "=== Chain submitted. Monitor with: squeue -u \$USER ==="
echo "    Cancel all: scancel $PREV_JOB  (cancels dependents too)"
