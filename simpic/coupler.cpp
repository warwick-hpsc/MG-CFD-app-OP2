#include <stdio.h>
#include <mpi.h>
#include "simpic_lib.h"
int main(int argc, char **argv){
    MPI_Init(&argc, &argv);
    MPI_Fint comms_shell = MPI_Comm_c2f(MPI_COMM_WORLD);
    main_simpic(argc, argv, comms_shell);
}