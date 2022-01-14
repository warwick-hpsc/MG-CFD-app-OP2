//
// auto-generated by op2.py
//

#include <math.h>
#include "const.h"

void up_post_kernel_omp4_kernel(
  double *data0,
  int dat0size,
  int *data1,
  int dat1size,
  int count,
  int num_teams,
  int nthread,
  int direct_up_post_kernel_stride_OP2CONSTANT){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    double* variable = &data0[n_op];
    const int* up_scratch = &data1[1*n_op];

    //inline function
    
      double avg = (*up_scratch)==0 ? 1.0 : 1.0 / (double)(*up_scratch);
      variable[VAR_DENSITY] *= avg;
      variable[VAR_MOMENTUM+0] *= avg;
      variable[VAR_MOMENTUM+1] *= avg;
      variable[VAR_MOMENTUM+2] *= avg;
      variable[VAR_DENSITY_ENERGY] *= avg;
    //end inline func
  }

}
