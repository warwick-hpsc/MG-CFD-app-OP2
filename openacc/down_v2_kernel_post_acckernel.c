//
// auto-generated by op2.py
//

//user function
#include <math.h>
#include "const.h"

int direct_down_v2_kernel_post_stride_OP2CONSTANT;
int direct_down_v2_kernel_post_stride_OP2HOST=-1;
//user function
//#pragma acc routine
inline void down_v2_kernel_post_openacc( 
    const double* residuals1_prolonged,
    const double* residuals1_prolonged_wsum,
    const double* residuals2,
    double* variables2) {


    for (int i=0; i<NVAR; i++) {
        variables2[(i)*direct_down_v2_kernel_post_stride_OP2CONSTANT] += residuals2[(i)*direct_down_v2_kernel_post_stride_OP2CONSTANT] - (residuals1_prolonged[(i)*direct_down_v2_kernel_post_stride_OP2CONSTANT] / (*residuals1_prolonged_wsum));
    }
}

// host stub function
void op_par_loop_down_v2_kernel_post(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(21);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[21].name      = name;
  OP_kernels[21].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  down_v2_kernel_post");
  }

  op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set->size >0) {

    if ((OP_kernels[21].count==1) || (direct_down_v2_kernel_post_stride_OP2HOST != getSetSizeFromOpArg(&arg0))) {
      direct_down_v2_kernel_post_stride_OP2HOST = getSetSizeFromOpArg(&arg0);
      direct_down_v2_kernel_post_stride_OP2CONSTANT = direct_down_v2_kernel_post_stride_OP2HOST;
    }

    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    double* data3 = (double*)arg3.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2,data3)
    for ( int n=0; n<set->size; n++ ){
      down_v2_kernel_post_openacc(
        &data0[n],
        &data1[1*n],
        &data2[n],
        &data3[n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[21].time     += wall_t2 - wall_t1;
  OP_kernels[21].transfer += (float)set->size * arg0.size;
  OP_kernels[21].transfer += (float)set->size * arg1.size;
  OP_kernels[21].transfer += (float)set->size * arg2.size;
  OP_kernels[21].transfer += (float)set->size * arg3.size * 2.0f;
}
