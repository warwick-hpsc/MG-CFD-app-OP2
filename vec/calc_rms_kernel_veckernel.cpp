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

#ifndef VALIDATION_H
#define VALIDATION_H

#include "utils.h"

inline void residual_kernel(
    const double* old_variable, 
    const double* variable, 
    double* residual)
{
    for (int v=0; v<NVAR; v++) {
        residual[v] = variable[v] - old_variable[v];
    }
}

inline void calc_rms_kernel(
    const double* residual, 
    double* rms)
{
    for (int i=0; i<NVAR; i++) {
        *rms += residual[i]*residual[i];
    }
}

inline void identify_differences(
    const double* test_value,
    const double* master_value, 
    double* difference)
{
    // If floating-point operations have been reordered, then a difference
    // is expected due to rounding-errors, but the difference should
    // be smaller than the following:
    //   1 x 10 ^ ( E - 17 + N )
    // Where E = exponent of master value
    //       N = largest difference in exponents of any floating-point
    //           arithmetic operation performed

    // N represents how many of the least-significant base-10 digits
    // of the floating-point mantissa are allowed to differ due to 
    // FP arithmetic reordering. Its value is guessed as 8, as to set it
    // accurately would require a trace of all floating-point operation
    // outputs during the runs.

    const double acceptable_relative_difference = 10.0e-8;

    for (int v=0; v<NVAR; v++) {
        double acceptable_difference = master_value[v] * acceptable_relative_difference;
        if (acceptable_difference < 0.0) {
            acceptable_difference *= -1.0;
        }

        // Ignore any differences smaller than 3e-19:
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

inline void count_non_zeros(
    const double* value, 
    int* count)
{   
    for (int v=0; v<NVAR; v++) {
        if ((*value) > 0.0) {
            (*count)++;
        }
    }
}

inline void count_bad_vals(
    const double* value, 
    int* count)
{   
    #ifdef OPENACC
        // OpenACC compilation is complaining about use of isnan()
    #else
        for (int v=0; v<NVAR; v++) {
            if (isnan(value[v]) || isinf(value[v])) {
                *count += 1;
            }
        }
    #endif
}

#endif

// host stub function
void op_par_loop_calc_rms_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;
  //create aligned pointers for dats
  ALIGNED_double const double * __restrict__ ptr0 = (double *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(14);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  calc_rms_kernel");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      double dat1[SIMD_VEC];
      for ( int i=0; i<SIMD_VEC; i++ ){
        dat1[i] = 0.0;
      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        calc_rms_kernel(
          &(ptr0)[5 * (n+i)],
          &dat1[i]);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        *(double*)arg1.data += dat1[i];
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      calc_rms_kernel(
        &(ptr0)[5*n],
        (double*)arg1.data);
    }
  }

  // combine reduction data
  op_mpi_reduce(&arg1,(double*)arg1.data);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[14].name      = name;
  OP_kernels[14].count    += 1;
  OP_kernels[14].time     += wall_t2 - wall_t1;
  OP_kernels[14].transfer += (float)set->size * arg0.size;
}
