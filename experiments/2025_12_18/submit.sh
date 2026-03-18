#!/bin/bash -l
#SBATCH --job-name=baseline
#SBATCH --output=logs/baseline-%j.out
#SBATCH --error=logs/baseline-%j.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00
#SBATCH --partition=express3

#SBATCH --mail-type=ALL
#SBATCH --mail-user=slurm@kamilbelan.anonaddy.com

# Submit from the project root:
#   sbatch submit_batch_snehurka.sh \
#     experiments/2025_12_18_baseline/global_params.toml \
#     experiments/2025_12_18_baseline/sim_params.toml

sbatch submit_batch_snehurka.sh \
  experiments/2025_12_18_baseline/global_params.toml \
  experiments/2025_12_18_baseline/sim_params.toml
