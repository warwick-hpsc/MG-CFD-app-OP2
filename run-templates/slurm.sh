#!/bin/bash

#SBATCH --job-name=MG-CFD.<RUN ID>
#SBATCH --output=sbatch.stdout
#SBATCH --error=sbatch.stderr

#SBATCH --nodes=<NODES>
#SBATCH --ntasks=<NTASKS>
#SBATCH --ntasks-per-node=<TPN>
#SBATCH --cpus-per-task=<NTHREADS>
#SBATCH --ntasks-per-core=1
#SBATCH --exclusive
#SBATCH -m block:block

#SBATCH --time=<HOURS>:<MINUTES>:00
#SBATCH --partition=<PARTITION>
#SBATCH --account=<PROJECT CODE>

RUN_CMD="srun --cpu_bind=cores"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
