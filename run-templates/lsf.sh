#!/bin/bash

#BSUB -J <RUN_ID>.MG-CFD.run.sbatch
#BSUB -o sbatch.stdout
#BSUB -e sbatch.stderr

#BSUB -nnodes <NODES>
#BSUB -R "span[ptile=<NCPUS_PER_NODE>] affinity[core(<NTHREADS>, same=socket)]"
#BSUB -W <HOURS>:<MINUTES>
#BSUB -x
#BSUB -q <PARTITION>
#BSUB -P <PROJECT CODE>

RUN_CMD="mpirun -n <NTASKS>"