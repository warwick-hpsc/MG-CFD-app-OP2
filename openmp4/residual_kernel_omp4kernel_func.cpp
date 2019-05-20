//
// auto-generated by op2.py
//

#include ".././src/Kernels/validation.h"

void residual_kernel_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  int count,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int n_op=0; n_op<count; n_op++ ){
    //variable mapping
    const double* old_variable = &data0[5*n_op];
    const double* variable = &data1[5*n_op];
    double* residual = &data2[5*n_op];

    //inline function
    
      for (int v=0; v<NVAR; v++) {
          residual[v] = variable[v] - old_variable[v];
      }
    //end inline func
  }

}
