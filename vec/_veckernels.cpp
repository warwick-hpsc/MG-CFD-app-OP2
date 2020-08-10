//
// auto-generated by op2.py
//

#define double_ALIGN 128
#define float_ALIGN 64
#define int_ALIGN 64
#ifdef VECTORIZE
#define SIMD_VEC 4
#define ALIGNED_double __attribute__((aligned(double_ALIGN)))
#define ALIGNED_float __attribute__((aligned(float_ALIGN)))
#define ALIGNED_int __attribute__((aligned(int_ALIGN)))
#ifdef __ICC
  #define DECLARE_PTR_ALIGNED(X, Y) __assume_aligned(X, Y)
#else
  #define DECLARE_PTR_ALIGNED(X, Y)
#endif
#else
#define ALIGNED_double
#define ALIGNED_float
#define ALIGNED_int
#define DECLARE_PTR_ALIGNED(X, Y)
#endif

// global constants

// header
#include "op_lib_cpp.h"

#ifdef PAPI
#include <papi.h>
#endif
void op_par_loop_compute_flux_edge_kernel_instrumented(
  char const *name, op_set set,
  op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3, op_arg arg4
  #ifdef VERIFY_OP2_TIMING
    , double* compute_time_ptr, double* sync_time_ptr
  #endif
  , long* iter_counts_ptr
  #ifdef PAPI
  , long_long* restrict event_counts, int event_set, int num_events
  #endif
);
// user kernel files
#include "initialize_variables_kernel_veckernel.cpp"
#include "zero_5d_array_kernel_veckernel.cpp"
#include "zero_1d_array_kernel_veckernel.cpp"
#include "calculate_cell_volumes_veckernel.cpp"
#include "dampen_ewt_veckernel.cpp"
#include "copy_double_kernel_veckernel.cpp"
#include "calculate_dt_kernel_veckernel.cpp"
#include "get_min_dt_kernel_veckernel.cpp"
#include "compute_step_factor_kernel_veckernel.cpp"
#include "compute_flux_edge_kernel_veckernel.cpp"
#include "compute_bnd_node_flux_kernel_veckernel.cpp"
#include "time_step_kernel_veckernel.cpp"
#include "indirect_rw_kernel_veckernel.cpp"
#include "residual_kernel_veckernel.cpp"
#include "calc_rms_kernel_veckernel.cpp"
#include "count_bad_vals_veckernel.cpp"
#include "up_pre_kernel_veckernel.cpp"
#include "up_kernel_veckernel.cpp"
#include "up_post_kernel_veckernel.cpp"
#include "down_v2_kernel_pre_veckernel.cpp"
#include "down_v2_kernel_veckernel.cpp"
#include "down_v2_kernel_post_veckernel.cpp"
#include "down_kernel_veckernel.cpp"
#include "identify_differences_veckernel.cpp"
#include "count_non_zeros_veckernel.cpp"
