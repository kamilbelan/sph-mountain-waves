#!/bin/bash -l
#SBATCH --job-name=example
#SBATCH --output=example.%j.out
#SBATCH --error=example.%j.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=3:00:00
#SBATCH --partition=express3

#SBATCH --mail-type=ALL
#SBATCH --mail-user=belank@karlin.mff.cuni.cz

echo "=== JOB START $(date) on $(hostname) ==="

#module load julia/1.9.1
#module load spack-bin-julia/.1.9.1-gcc9.4.0-uzbijvaynhque2ch

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Running on $(hostname) with $JULIA_NUM_THREADS threads"
#echo "Julia binary: $(which julia)"

cd /usr/users/belank/sph-mountain-waves/jobs/2025_18_12

/usr/work/belank/julia-1.12.4/bin/julia -t $JULIA_NUM_THREADS example_job.jl

echo "=== JOB END $(date) ==="
