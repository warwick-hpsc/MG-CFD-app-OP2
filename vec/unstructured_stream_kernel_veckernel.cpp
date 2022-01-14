//
// auto-generated by op2.py
//

//user function
#ifndef UNSTRUCTURED_STREAM_H
#define UNSTRUCTURED_STREAM_H

// Unstructured stream kernel
// - performs same data movement as compute_flux_edge() but with minimal arithmetic. 
//   Measures upper bound on performance achievable by compute_flux_edge()
inline void unstructured_stream_kernel(
    const double *variables_a,
    const double *variables_b,
    const double *edge_weight,
    double *fluxes_a, 
    double *fluxes_b)
{
    double ex = edge_weight[0];
    double ey = edge_weight[1];
    double ez = edge_weight[2];

    double p_a, pe_a;
    double3 momentum_a;
    p_a          = variables_a[VAR_DENSITY];
    momentum_a.x = variables_a[VAR_MOMENTUM+0];
    momentum_a.y = variables_a[VAR_MOMENTUM+1];
    momentum_a.z = variables_a[VAR_MOMENTUM+2];
    pe_a         = variables_a[VAR_DENSITY_ENERGY];

    double p_b, pe_b;
    double3 momentum_b;
    p_b          = variables_b[VAR_DENSITY];
    momentum_b.x = variables_b[VAR_MOMENTUM+0];
    momentum_b.y = variables_b[VAR_MOMENTUM+1];
    momentum_b.z = variables_b[VAR_MOMENTUM+2];
    pe_b         = variables_b[VAR_DENSITY_ENERGY];

    double p_a_val  = p_b + ex;
    double pe_a_val = pe_b + ey;
    double mx_a_val = momentum_b.x + ez;
    double my_a_val = momentum_b.y;
    double mz_a_val = momentum_b.z;

    double p_b_val = p_a;
    double pe_b_val = pe_a;
    double mx_b_val = momentum_a.x;
    double my_b_val = momentum_a.y;
    double mz_b_val = momentum_a.z;

    fluxes_a[VAR_DENSITY]  += p_a_val;
    fluxes_a[VAR_MOMENTUM+0] += mx_a_val;
    fluxes_a[VAR_MOMENTUM+1] += my_a_val;
    fluxes_a[VAR_MOMENTUM+2] += mz_a_val;
    fluxes_a[VAR_DENSITY_ENERGY] += pe_a_val;

    fluxes_b[VAR_DENSITY]  += p_b_val;
    fluxes_b[VAR_MOMENTUM+0] += mx_b_val;
    fluxes_b[VAR_MOMENTUM+1] += my_b_val;
    fluxes_b[VAR_MOMENTUM+2] += mz_b_val;
    fluxes_b[VAR_DENSITY_ENERGY] += pe_b_val;
}

#endif
#ifdef VECTORIZE
//user function -- modified for vectorisation
#if defined __clang__ || defined __GNUC__
__attribute__((always_inline))
#endif
inline void unstructured_stream_kernel_vec( const double variables_a[][SIMD_VEC], const double variables_b[][SIMD_VEC], const double edge_weight[][SIMD_VEC], double fluxes_a[][SIMD_VEC], double fluxes_b[][SIMD_VEC], int idx ) {
    double ex = edge_weight[0][idx];
    double ey = edge_weight[1][idx];
    double ez = edge_weight[2][idx];

    double p_a, pe_a;
    double3 momentum_a;
    p_a          = variables_a[VAR_DENSITY][idx];
    momentum_a.x = variables_a[VAR_MOMENTUM+0][idx];
    momentum_a.y = variables_a[VAR_MOMENTUM+1][idx];
    momentum_a.z = variables_a[VAR_MOMENTUM+2][idx];
    pe_a         = variables_a[VAR_DENSITY_ENERGY][idx];

    double p_b, pe_b;
    double3 momentum_b;
    p_b          = variables_b[VAR_DENSITY][idx];
    momentum_b.x = variables_b[VAR_MOMENTUM+0][idx];
    momentum_b.y = variables_b[VAR_MOMENTUM+1][idx];
    momentum_b.z = variables_b[VAR_MOMENTUM+2][idx];
    pe_b         = variables_b[VAR_DENSITY_ENERGY][idx];

    double p_a_val  = p_b + ex;
    double pe_a_val = pe_b + ey;
    double mx_a_val = momentum_b.x + ez;
    double my_a_val = momentum_b.y;
    double mz_a_val = momentum_b.z;

    double p_b_val = p_a;
    double pe_b_val = pe_a;
    double mx_b_val = momentum_a.x;
    double my_b_val = momentum_a.y;
    double mz_b_val = momentum_a.z;

    fluxes_a[VAR_DENSITY][idx]  = p_a_val;
    fluxes_a[VAR_MOMENTUM+0][idx] = mx_a_val;
    fluxes_a[VAR_MOMENTUM+1][idx] = my_a_val;
    fluxes_a[VAR_MOMENTUM+2][idx] = mz_a_val;
    fluxes_a[VAR_DENSITY_ENERGY][idx] = pe_a_val;

    fluxes_b[VAR_DENSITY][idx]  = p_b_val;
    fluxes_b[VAR_MOMENTUM+0][idx] = mx_b_val;
    fluxes_b[VAR_MOMENTUM+1][idx] = my_b_val;
    fluxes_b[VAR_MOMENTUM+2][idx] = mz_b_val;
    fluxes_b[VAR_DENSITY_ENERGY][idx] = pe_b_val;

}
#endif

// host stub function
void op_par_loop_unstructured_stream_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  //create aligned pointers for dats
  ALIGNED_double const double * __restrict__ ptr0 = (double *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr1 = (double *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr2 = (double *) arg2.data;
  DECLARE_PTR_ALIGNED(ptr2,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr3 = (double *) arg3.data;
  DECLARE_PTR_ALIGNED(ptr3,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr4 = (double *) arg4.data;
  DECLARE_PTR_ALIGNED(ptr4,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(12);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: unstructured_stream_kernel\n");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      if (n<set->core_size && n>0 && n % OP_mpi_test_frequency == 0)
        op_mpi_test_all(nargs,args);
      if ((n+SIMD_VEC >= set->core_size) && (n+SIMD_VEC-set->core_size < SIMD_VEC)) {
        op_mpi_wait_all(nargs, args);
      }
      ALIGNED_double double dat0[5][SIMD_VEC];
      ALIGNED_double double dat1[5][SIMD_VEC];
      ALIGNED_double double dat2[3][SIMD_VEC];
      ALIGNED_double double dat3[5][SIMD_VEC];
      ALIGNED_double double dat4[5][SIMD_VEC];
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx0_5 = 5 * arg0.map_data[(n+i) * arg0.map->dim + 0];
        int idx1_5 = 5 * arg0.map_data[(n+i) * arg0.map->dim + 1];
        int idx2_3 = 3 * (n+i);

        dat0[0][i] = (ptr0)[idx0_5 + 0];
        dat0[1][i] = (ptr0)[idx0_5 + 1];
        dat0[2][i] = (ptr0)[idx0_5 + 2];
        dat0[3][i] = (ptr0)[idx0_5 + 3];
        dat0[4][i] = (ptr0)[idx0_5 + 4];

        dat1[0][i] = (ptr1)[idx1_5 + 0];
        dat1[1][i] = (ptr1)[idx1_5 + 1];
        dat1[2][i] = (ptr1)[idx1_5 + 2];
        dat1[3][i] = (ptr1)[idx1_5 + 3];
        dat1[4][i] = (ptr1)[idx1_5 + 4];

        dat2[0][i] = (ptr2)[idx2_3 + 0];
        dat2[1][i] = (ptr2)[idx2_3 + 1];
        dat2[2][i] = (ptr2)[idx2_3 + 2];

        dat3[0][i] = 0.0;
        dat3[1][i] = 0.0;
        dat3[2][i] = 0.0;
        dat3[3][i] = 0.0;
        dat3[4][i] = 0.0;

        dat4[0][i] = 0.0;
        dat4[1][i] = 0.0;
        dat4[2][i] = 0.0;
        dat4[3][i] = 0.0;
        dat4[4][i] = 0.0;

      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        unstructured_stream_kernel_vec(
          dat0,
          dat1,
          dat2,
          dat3,
          dat4,
          i);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx3_5 = 5 * arg0.map_data[(n+i) * arg0.map->dim + 0];
        int idx4_5 = 5 * arg0.map_data[(n+i) * arg0.map->dim + 1];

        (ptr3)[idx3_5 + 0] += dat3[0][i];
        (ptr3)[idx3_5 + 1] += dat3[1][i];
        (ptr3)[idx3_5 + 2] += dat3[2][i];
        (ptr3)[idx3_5 + 3] += dat3[3][i];
        (ptr3)[idx3_5 + 4] += dat3[4][i];

        (ptr4)[idx4_5 + 0] += dat4[0][i];
        (ptr4)[idx4_5 + 1] += dat4[1][i];
        (ptr4)[idx4_5 + 2] += dat4[2][i];
        (ptr4)[idx4_5 + 3] += dat4[3][i];
        (ptr4)[idx4_5 + 4] += dat4[4][i];

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
      int map0idx;
      int map1idx;
      map0idx = arg0.map_data[n * arg0.map->dim + 0];
      map1idx = arg0.map_data[n * arg0.map->dim + 1];

      unstructured_stream_kernel(
        &(ptr0)[5 * map0idx],
        &(ptr1)[5 * map1idx],
        &(ptr2)[3 * n],
        &(ptr3)[5 * map0idx],
        &(ptr4)[5 * map1idx]);
    }
  }

  if (exec_size == 0 || exec_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[12].name      = name;
  OP_kernels[12].count    += 1;
  OP_kernels[12].time     += wall_t2 - wall_t1;
  OP_kernels[12].transfer += (float)set->size * arg0.size;
  OP_kernels[12].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[12].transfer += (float)set->size * arg2.size;
  OP_kernels[12].transfer += (float)set->size * arg0.map->dim * 4.0f;
}
