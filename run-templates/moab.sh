#!/bin/bash

#MSUB -N MG-CFD.<RUN ID>
#MSUB -o moab.stdout
#MSUB -e moab.stderr

#MSUB -l nodes=<NODES>
#MSUB -l procs=<NTASKS>
#MSUB -l ttc=<TPN>
#MSUB -n

#MSUB -l walltime=<HOURS>:<MINUTES>:00
#MSUB -q <PARTITION>
#MSUB -A <PROJECT CODE>

RUN_CMD="srun --cpus-per-task=<NTHREADS>"

export OMP_NUM_THREADS=$PBS_NUM_PPN
export OMP_PROC_BIND=true
