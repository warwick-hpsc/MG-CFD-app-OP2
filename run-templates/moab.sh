#!/bin/bash

#MSUB -N <RUN_ID>.MG-CFD.run.sbatch
#MSUB -o sbatch.stdout
#MSUB -e sbatch.stderr

#MSUB -l nodes=<NODES>
#MSUB -l ttc=<TPN>
#MSUB -l walltime=<HOURS>:<MINUTES>:00
#MSUB -n
#MSUB -q <PARTITION>
#MSUB -A <PROJECT CODE>

RUN_CMD="srun --cpus-per-task=<NTHREADS> -n <NTASKS>"