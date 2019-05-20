//
// auto-generated by op2.py
//

#include "utils.h"

//user function
__device__ void identify_differences_gpu( 
    const double* test_value,
    const double* master_value,
    double* difference) {













    const double acceptable_relative_difference = 10.0e-9;

    for (int v=0; v<NVAR; v++) {
        double acceptable_difference = master_value[v] * acceptable_relative_difference;
        if (acceptable_difference < 0.0) {
            acceptable_difference *= -1.0;
        }

        if (acceptable_difference < 3.0e-19) {
            acceptable_difference = 3.0e-19;
        }

        double diff = test_value[v] - master_value[v];
        if (diff < 0.0) {
            diff *= -1.0;
        }

        if (diff > acceptable_difference) {
            difference[v] = diff;
        } else {
            difference[v] = 0.0;
        }
    }

}

// CUDA kernel function
__global__ void op_cuda_identify_differences(
  const double *__restrict arg0,
  const double *__restrict arg1,
  double *arg2,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    identify_differences_gpu(arg0+n*5,
                         arg1+n*5,
                         arg2+n*5);
  }
}


//host stub function
void op_par_loop_identify_differences(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];
  const int nk = 22;

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(nk);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[nk].name      = name;
  OP_kernels[nk].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  identify_differences");
  }

  op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set->size > 0) {

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_22
      int nthread = OP_BLOCK_SIZE_22;
    #else
      int nthread = OP_block_size;
    //  int nthread = 128;
    #endif

    int nblocks = 200;

    op_cuda_identify_differences<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (double *) arg1.data_d,
      (double *) arg2.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[nk].time     += wall_t2 - wall_t1;
  OP_kernels[nk].transfer += (float)set->size * arg0.size;
  OP_kernels[nk].transfer += (float)set->size * arg1.size;
  OP_kernels[nk].transfer += (float)set->size * arg2.size * 2.0f;
}
