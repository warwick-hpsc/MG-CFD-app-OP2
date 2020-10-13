//
// auto-generated by op2.py
//

//user function
#include "utils.h"

//user function
//#pragma acc routine
inline void identify_differences_openacc( 
    const double* test_value,
    const double* master_value,
    double* difference) {













    const double acceptable_relative_difference = 10.0e-8;

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

// host stub function
void op_par_loop_identify_differences(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(23);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[23].name      = name;
  OP_kernels[23].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  identify_differences");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  if (set_size >0) {


    //Set up typed device pointers for OpenACC

    double* data0 = (double*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double* data2 = (double*)arg2.data_d;
    #pragma acc parallel loop independent deviceptr(data0,data1,data2)
    for ( int n=0; n<set->size; n++ ){
      identify_differences_openacc(
        &data0[5*n],
        &data1[5*n],
        &data2[5*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[23].time     += wall_t2 - wall_t1;
  OP_kernels[23].transfer += (float)set->size * arg0.size;
  OP_kernels[23].transfer += (float)set->size * arg1.size;
  OP_kernels[23].transfer += (float)set->size * arg2.size * 2.0f;
}
