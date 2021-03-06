//
// auto-generated by op2.py
//

#include "utils.h"

void identify_differences_omp4_kernel(
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
    const double* test_value = &data0[5*n_op];
    const double* master_value = &data1[5*n_op];
    double* difference = &data2[5*n_op];

    //inline function
    













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
    //end inline func
  }

}
