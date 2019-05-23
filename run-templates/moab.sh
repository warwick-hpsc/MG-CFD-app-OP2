#!/bin/bash

#MSUB -N MG-CFD.<RUN ID>
#MSUB -o moab.stdout
#MSUB -e moab.stderr

#MSUB -l nodes=<NODES>
#MSUB -l ttc=<TPN>
#MSUB -l walltime=<HOURS>:<MINUTES>:00
#MSUB -n
#MSUB -q <PARTITION>
#MSUB -A <PROJECT CODE>

RUN_CMD="srun --cpus-per-task=<NTHREADS> -n <NTASKS>"
