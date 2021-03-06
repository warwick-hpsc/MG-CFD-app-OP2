//
// auto-generated by op2.py
//

#include <math.h>
#include "const.h"

void down_v2_kernel_pre_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    double* weight_sum = &data0[5*n_op];
    double* residual_sum = &data1[1*n_op];

    //inline function
    
      *weight_sum = 0.0;
      residual_sum[VAR_DENSITY] = 0.0;
      residual_sum[VAR_MOMENTUM+0] = 0.0;
      residual_sum[VAR_MOMENTUM+2] = 0.0;
      residual_sum[VAR_MOMENTUM+1] = 0.0;
      residual_sum[VAR_DENSITY_ENERGY] = 0.0;
    //end inline func
  }

}
