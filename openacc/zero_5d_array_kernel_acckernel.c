//
// auto-generated by op2.py
//

//user function
#include ".././src/Kernels/misc.h"
//user function
//#pragma acc routine
inline void zero_5d_array_kernel_openacc( 
    double* array) {
    for(int j = 0; j < NVAR; j++) {
        array[j] = 0.0;
    }
}

// host stub function
void op_par_loop_zero_5d_array_kernel(char const *name, op_set set,
  op_arg arg0){

  int nargs = 1;
  op_arg args[1];
  const int nk = 1;

  args[0] = arg0;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(nk);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[nk].name      = name;
  OP_kernels[nk].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  zero_5d_array_kernel");
  }

  op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set->size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    #pragma acc parallel loop independent deviceptr(data0)
    for ( int n=0; n<set->size; n++ ){
      zero_5d_array_kernel_openacc(
        &data0[5*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[nk].time     += wall_t2 - wall_t1;
  OP_kernels[nk].transfer += (float)set->size * arg0.size * 2.0f;
}
