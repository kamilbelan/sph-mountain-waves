#!/bin/bash -l
#SBATCH --job-name=highres
#SBATCH --output=logs/highres-%j.out
#SBATCH --error=logs/highres-%j.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --partition=express3

#SBATCH --mail-type=ALL
#SBATCH --mail-user=slurm@kamilbelan.anonaddy.com

# Submit from the project root:
#   sbatch submit_batch_snehurka.sh \
#     experiments/2026_01_09_highres/global_params.toml \
#     experiments/2026_01_09_highres/sim_params.toml

sbatch submit_batch_snehurka.sh \
  experiments/2026_01_09_highres/global_params.toml \
  experiments/2026_01_09_highres/sim_params.toml
