#!/bin/bash

NT=200  # of time steps
PPC=100 # of particles per cell
CPP=64  # of cell per process
DTFACTOR=0.001 #defines a fraction of dt vs. time steps from plasma frequency;
               #must be positive
NPROC=4 # of processors
LHSV=20000 #applied voltage on left-hand side; RHS is grounded; 
PROGNAME=a.out

CMDLINE="mpirun -np $NPROC $PROGNAME -ppc $PPC -ncpp $CPP -nt $NT -dtfactor $DTFACTOR -lhsv $LHSV"
echo $CMDLINE
$CMDLINE
