//
// auto-generated by op2.py
//

#include <math.h>
#include "const.h"

void down_kernel_omp4_kernel(
  double *data0,
  int dat0size,
  double *data1,
  int dat1size,
  double *data2,
  int dat2size,
  int *map3,
  int map3size,
  double *data3,
  int dat3size,
  double *data4,
  int dat4size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data0[0:dat0size],data1[0:dat1size],data2[0:dat2size])\
    map(to:col_reord[0:set_size1],map3[0:map3size],data3[0:dat3size],data4[0:dat4size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map3idx = map3[n_op + set_size1 * 0];

    //variable mapping
    double* variable = &data0[5*n_op];
    const double* residual = &data1[5*n_op];
    const double* coord = &data2[3*n_op];
    const double* residual_above = &data3[5 * map3idx];
    const double* coord_above = &data4[3 * map3idx];

    //inline function
    
      double dx = fabs(coord[0] - coord_above[0]);
      double dy = fabs(coord[1] - coord_above[1]);
      double dz = fabs(coord[2] - coord_above[2]);
      double dm = sqrt(dx*dx + dy*dy + dz*dz);

      variable[VAR_DENSITY]        -= dm* (residual_above[VAR_DENSITY]        - residual[VAR_DENSITY]);
      variable[VAR_MOMENTUM+0]     -= dx* (residual_above[VAR_MOMENTUM+0]     - residual[VAR_MOMENTUM+0]);
      variable[VAR_MOMENTUM+1]     -= dy* (residual_above[VAR_MOMENTUM+1]     - residual[VAR_MOMENTUM+1]);
      variable[VAR_MOMENTUM+2]     -= dz* (residual_above[VAR_MOMENTUM+2]     - residual[VAR_MOMENTUM+2]);
      variable[VAR_DENSITY_ENERGY] -= dm* (residual_above[VAR_DENSITY_ENERGY] - residual[VAR_DENSITY_ENERGY]);
    //end inline func
  }

}
