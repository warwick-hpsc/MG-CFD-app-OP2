//
// auto-generated by op2.py
//

//user function
//************************************************//
// Copyright 2016-2019 University of Warwick

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is furnished
// to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//************************************************//

#ifndef MG_H
#define MG_H

#include <math.h>

#include "const.h"

inline void up_pre_kernel(
    double* variable,
    int* up_scratch)
{
    variable[VAR_DENSITY] = 0.0;
    variable[VAR_MOMENTUM+0] = 0.0;
    variable[VAR_MOMENTUM+1] = 0.0;
    variable[VAR_MOMENTUM+2] = 0.0;
    variable[VAR_DENSITY_ENERGY] = 0.0;
    *up_scratch = 0;
}

inline void up_kernel(
    const double* variable,
    double* variable_above,
    int* up_scratch)
{
    variable_above[VAR_DENSITY]        += variable[VAR_DENSITY];
    variable_above[VAR_MOMENTUM+0]     += variable[VAR_MOMENTUM+0];
    variable_above[VAR_MOMENTUM+1]     += variable[VAR_MOMENTUM+1];
    variable_above[VAR_MOMENTUM+2]     += variable[VAR_MOMENTUM+2];
    variable_above[VAR_DENSITY_ENERGY] += variable[VAR_DENSITY_ENERGY];
    *up_scratch += 1;
}

inline void up_post_kernel(
    double* variable,
    const int* up_scratch)
{
    double avg = (*up_scratch)==0 ? 1.0 : 1.0 / (double)(*up_scratch);
    variable[VAR_DENSITY] *= avg;
    variable[VAR_MOMENTUM+0] *= avg;
    variable[VAR_MOMENTUM+1] *= avg;
    variable[VAR_MOMENTUM+2] *= avg;
    variable[VAR_DENSITY_ENERGY] *= avg;
}

inline void down_kernel(
    double* variable,
    const double* residual,
    const double* coord,
    const double* residual_above,
    const double* coord_above)
{
    double dx = fabs(coord[0] - coord_above[0]);
    double dy = fabs(coord[1] - coord_above[1]);
    double dz = fabs(coord[2] - coord_above[2]);
    double dm = sqrt(dx*dx + dy*dy + dz*dz);

    variable[VAR_DENSITY]        -= dm* (residual_above[VAR_DENSITY]        - residual[VAR_DENSITY]);
    variable[VAR_MOMENTUM+0]     -= dx* (residual_above[VAR_MOMENTUM+0]     - residual[VAR_MOMENTUM+0]);
    variable[VAR_MOMENTUM+1]     -= dy* (residual_above[VAR_MOMENTUM+1]     - residual[VAR_MOMENTUM+1]);
    variable[VAR_MOMENTUM+2]     -= dz* (residual_above[VAR_MOMENTUM+2]     - residual[VAR_MOMENTUM+2]);
    variable[VAR_DENSITY_ENERGY] -= dm* (residual_above[VAR_DENSITY_ENERGY] - residual[VAR_DENSITY_ENERGY]);
}

inline void down_v2_kernel_pre(
    double* weight_sum,
    double* residual_sum)
{
    *weight_sum = 0.0;
    residual_sum[VAR_DENSITY] = 0.0;
    residual_sum[VAR_MOMENTUM+0] = 0.0;
    residual_sum[VAR_MOMENTUM+2] = 0.0;
    residual_sum[VAR_MOMENTUM+1] = 0.0;
    residual_sum[VAR_DENSITY_ENERGY] = 0.0;
}
inline void down_v2_kernel(
    const double* coord2a,
    const double* coord2b,
    const double* coord1a,
    const double* coord1b,
    const double* residuals1a,
    const double* residuals1b,
    double* residuals1a_prolonged,
    double* residuals1b_prolonged,
    double* residuals1a_prolonged_wsum,
    double* residuals1b_prolonged_wsum)
{
    // For each node that has the same coordinates as its MG node parent,
    // the 'prolonged residual' is simply taken directly from the MG node.
    //
    // For each other node N, the 'prolonged residuals' is the weighted average
    // across N's MG node and MG nodes of N's neighbours, requiring an
    // edge-based loop. The weight is 1.0/distance.

    // Process a2:
    double dx_a1a2 = coord2a[0] - coord1a[0];
    double dy_a1a2 = coord2a[1] - coord1a[1];
    double dz_a1a2 = coord2a[2] - coord1a[2];
    if (dx_a1a2 == 0.0 && dy_a1a2 == 0.0 && dz_a1a2 == 0.0) {
        // a2 == a1:
        residuals1a_prolonged[VAR_DENSITY]        = residuals1a[VAR_DENSITY];
        residuals1a_prolonged[VAR_MOMENTUM+0]     = residuals1a[VAR_MOMENTUM+0];
        residuals1a_prolonged[VAR_MOMENTUM+1]     = residuals1a[VAR_MOMENTUM+1];
        residuals1a_prolonged[VAR_MOMENTUM+2]     = residuals1a[VAR_MOMENTUM+2];
        residuals1a_prolonged[VAR_DENSITY_ENERGY] = residuals1a[VAR_DENSITY_ENERGY];
        *residuals1a_prolonged_wsum = 1.0;
    } else {
        // Calculate contribution of a1 -> a2:
        const double idist_a1a2 = 1.0/sqrt(dx_a1a2*dx_a1a2 + dy_a1a2*dy_a1a2 + dz_a1a2*dz_a1a2);
        residuals1a_prolonged[VAR_DENSITY]        += idist_a1a2*residuals1a[VAR_DENSITY];
        residuals1a_prolonged[VAR_MOMENTUM+0]     += idist_a1a2*residuals1a[VAR_MOMENTUM+0];
        residuals1a_prolonged[VAR_MOMENTUM+1]     += idist_a1a2*residuals1a[VAR_MOMENTUM+1];
        residuals1a_prolonged[VAR_MOMENTUM+2]     += idist_a1a2*residuals1a[VAR_MOMENTUM+2];
        residuals1a_prolonged[VAR_DENSITY_ENERGY] += idist_a1a2*residuals1a[VAR_DENSITY_ENERGY];
        *residuals1a_prolonged_wsum += idist_a1a2;

        // Calculate contribution of b1 >- a2:
        double dx_b1a2 = coord1b[0] - coord2a[0];
        double dy_b1a2 = coord1b[1] - coord2a[1];
        double dz_b1a2 = coord1b[2] - coord2a[2];

        const double idist_b1a2 = 1.0/sqrt(dx_b1a2*dx_b1a2 + dy_b1a2*dy_b1a2 + dz_b1a2*dz_b1a2);
        residuals1a_prolonged[VAR_DENSITY]        += idist_b1a2*residuals1b[VAR_DENSITY];
        residuals1a_prolonged[VAR_MOMENTUM+0]     += idist_b1a2*residuals1b[VAR_MOMENTUM+0];
        residuals1a_prolonged[VAR_MOMENTUM+1]     += idist_b1a2*residuals1b[VAR_MOMENTUM+1];
        residuals1a_prolonged[VAR_MOMENTUM+2]     += idist_b1a2*residuals1b[VAR_MOMENTUM+2];
        residuals1a_prolonged[VAR_DENSITY_ENERGY] += idist_b1a2*residuals1b[VAR_DENSITY_ENERGY];
        *residuals1a_prolonged_wsum += idist_b1a2;
    }

    // Process b2:
    double dx_b1b2 = coord2b[0] - coord1b[0];
    double dy_b1b2 = coord2b[1] - coord1b[1];
    double dz_b1b2 = coord2b[2] - coord1b[2];
    if (dx_b1b2 == 0.0 && dy_b1b2 == 0.0 && dz_b1b2 == 0.0) {
        // b2 == b1:
        residuals1b_prolonged[VAR_DENSITY]        = residuals1b[VAR_DENSITY];
        residuals1b_prolonged[VAR_MOMENTUM+0]      = residuals1b[VAR_MOMENTUM+0];
        residuals1b_prolonged[VAR_MOMENTUM+1]      = residuals1b[VAR_MOMENTUM+1];
        residuals1b_prolonged[VAR_MOMENTUM+2]      = residuals1b[VAR_MOMENTUM+2];
        residuals1b_prolonged[VAR_DENSITY_ENERGY] = residuals1b[VAR_DENSITY_ENERGY];
        *residuals1b_prolonged_wsum = 1.0;
    } else {
        // Calculate contribution of b1 -> b2:
        const double idist_b1b2 = 1.0/sqrt(dx_b1b2*dx_b1b2 + dy_b1b2*dy_b1b2 + dz_b1b2*dz_b1b2);
        residuals1b_prolonged[VAR_DENSITY]        += idist_b1b2*residuals1b[VAR_DENSITY];
        residuals1b_prolonged[VAR_MOMENTUM+0]     += idist_b1b2*residuals1b[VAR_MOMENTUM+0];
        residuals1b_prolonged[VAR_MOMENTUM+1]     += idist_b1b2*residuals1b[VAR_MOMENTUM+1];
        residuals1b_prolonged[VAR_MOMENTUM+2]     += idist_b1b2*residuals1b[VAR_MOMENTUM+2];
        residuals1b_prolonged[VAR_DENSITY_ENERGY] += idist_b1b2*residuals1b[VAR_DENSITY_ENERGY];
        *residuals1b_prolonged_wsum += idist_b1b2;

        // Calculate contribution of a1 -> b2:
        double dx_a1b2 = coord1a[0] - coord2b[0];
        double dy_a1b2 = coord1a[1] - coord2b[1];
        double dz_a1b2 = coord1a[2] - coord2b[2];

        const double idist_a1b2 = 1.0/sqrt(dx_a1b2*dx_a1b2 + dy_a1b2*dy_a1b2 + dz_a1b2*dz_a1b2);
        residuals1b_prolonged[VAR_DENSITY]        += idist_a1b2*residuals1b[VAR_DENSITY];
        residuals1b_prolonged[VAR_MOMENTUM+0]     += idist_a1b2*residuals1b[VAR_MOMENTUM+0];
        residuals1b_prolonged[VAR_MOMENTUM+1]     += idist_a1b2*residuals1b[VAR_MOMENTUM+1];
        residuals1b_prolonged[VAR_MOMENTUM+2]     += idist_a1b2*residuals1b[VAR_MOMENTUM+2];
        residuals1b_prolonged[VAR_DENSITY_ENERGY] += idist_a1b2*residuals1b[VAR_DENSITY_ENERGY];
        *residuals1b_prolonged_wsum += idist_a1b2;
    }
}

inline void down_v2_kernel_post(
    const double* residuals1_prolonged,
    const double* residuals1_prolonged_wsum,
    const double* residuals2,
    double* variables2)
{
    // Divide through by weight sum to complete the weighted average started by down_v2_kernel(),
    // then apply the prolonged residual to grid:
    for (int i=0; i<NVAR; i++) {
        variables2[i] += residuals2[i] - (residuals1_prolonged[i] / (*residuals1_prolonged_wsum));
    }
}

#endif
#ifdef VECTORIZE
//user function -- modified for vectorisation
inline void up_pre_kernel_vec( double variable[*][SIMD_VEC], int up_scratch[*][SIMD_VEC], int idx ) {
    variable[VAR_DENSITY][idx] = 0.0;
    variable[VAR_MOMENTUM+0][idx] = 0.0;
    variable[VAR_MOMENTUM+1][idx] = 0.0;
    variable[VAR_MOMENTUM+2][idx] = 0.0;
    variable[VAR_DENSITY_ENERGY][idx] = 0.0;
    up_scratch[0][idx]= 0;
}
#endif

// host stub function
void op_par_loop_up_pre_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;
  //create aligned pointers for dats
  ALIGNED_double       double * __restrict__ ptr0 = (double *) arg0.data;
  __assume_aligned(ptr0,double_ALIGN);
  ALIGNED_int       int * __restrict__ ptr1 = (int *) arg1.data;
  __assume_aligned(ptr1,int_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(15);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: up_pre_kernel\n");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      if (n+SIMD_VEC >= set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      ALIGNED_double double dat0[5][SIMD_VEC];
      ALIGNED_int int dat1[1][SIMD_VEC];
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx0_5 = 5 * arg0.map_data[(n+i) * arg0.map->dim + 0];
        int idx1_1 = 1 * arg0.map_data[(n+i) * arg0.map->dim + 0];

      }
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        up_pre_kernel_vec(
          dat0,
          dat1,
          i);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx0_5 = 5 * arg0.map_data[(n+i) * arg0.map->dim + 0];
        int idx1_1 = 1 * arg0.map_data[(n+i) * arg0.map->dim + 0];

        (ptr0)[idx0_5 + 0] = dat0[0][i];
        (ptr0)[idx0_5 + 1] = dat0[1][i];
        (ptr0)[idx0_5 + 2] = dat0[2][i];
        (ptr0)[idx0_5 + 3] = dat0[3][i];
        (ptr0)[idx0_5 + 4] = dat0[4][i];

        (ptr1)[idx1_1 + 0] = dat1[0][i];

      }
    }

    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map0idx = arg0.map_data[n * arg0.map->dim + 0];

      up_pre_kernel(
        &(ptr0)[5 * map0idx],
        &(ptr1)[1 * map0idx]);
    }
  }

  if (exec_size == 0 || exec_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[15].name      = name;
  OP_kernels[15].count    += 1;
  OP_kernels[15].time     += wall_t2 - wall_t1;
  OP_kernels[15].transfer += (float)set->size * arg0.size;
  OP_kernels[15].transfer += (float)set->size * arg1.size;
  OP_kernels[15].transfer += (float)set->size * arg0.map->dim * 4.0f;
}
