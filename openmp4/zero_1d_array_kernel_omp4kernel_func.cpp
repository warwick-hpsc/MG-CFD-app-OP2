//
// auto-generated by op2.py
//

#include ".././src/Kernels/misc.h"

void zero_1d_array_kernel_omp4_kernel(
  double *data0,
  int dat0size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    double* array = &data0[1*n_op];

    //inline function
    
      *array = 0.0;
    //end inline func
  }

}
