#!/bin/bash

#SBATCH --job-name=MG-CFD.<RUN ID>
#SBATCH --output=sbatch.stdout
#SBATCH --error=sbatch.stderr

#SBATCH --nodes=<NODES>
#SBATCH --ntasks-per-node=<TPN>
#SBATCH --cpus-per-task=<NTHREADS>
#SBATCH --time=<HOURS>:<MINUTES>:00
#SBATCH --exclusive
#SBATCH --partition=<PARTITION>
#SBATCH --account=<PROJECT CODE>

RUN_CMD="srun -n <NTASKS>"
