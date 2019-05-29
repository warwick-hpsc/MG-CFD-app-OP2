#!/bin/bash

#PBS -N MG-CFD.<RUN ID>
#PBS -o pbs.stdout
#PBS -e pbs.stderr

#PBS -l select=<NODES>
#PBS -l walltime=<HOURS>:<MINUTES>:00
#PBS -A <PROJECT CODE>
#PBS -q <PARTITION>

# RUN_CMD="aprun -n <NTASKS> -N <TPN> -d <NTHREADS>"
RUN_CMD="aprun -n <NTASKS> -N <TPN>"

export OMP_NUM_THREADS=<NTHREADS>
export OMP_PROC_BIND=true
