//
// auto-generated by op2.py
//

#define double_ALIGN 128
#define float_ALIGN 64
#define int_ALIGN 64
#ifdef VECTORIZE
#define SIMD_VEC 8
#define ALIGNED_double __attribute__((aligned(double_ALIGN)))
#define ALIGNED_float __attribute__((aligned(float_ALIGN)))
#define ALIGNED_int __attribute__((aligned(int_ALIGN)))
#ifdef __ICC
  #define DECLARE_PTR_ALIGNED(X, Y) __assume_aligned(X, Y)
#else
  #define DECLARE_PTR_ALIGNED(X, Y)
#endif
#include "../slope/unstructured_stream_kernel_veckernel.h"
#else
#define ALIGNED_double
#define ALIGNED_float
#define ALIGNED_int
#define DECLARE_PTR_ALIGNED(X, Y)
#endif

//user function
#include ".././src/Kernels/unstructured_stream.h"

#ifdef PAPI
#include "papi_funcs.h"
#endif

// host stub function
void op_par_loop_unstructured_stream_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  
  op_par_loop_unstructured_stream_kernel_instrumented(name, set, 
    arg0, arg1, arg2, arg3, arg4
    #ifdef VERIFY_OP2_TIMING
      , NULL, NULL
    #endif
    , NULL
    #ifdef PAPI
    , NULL, 0, 0
    #endif
    );
};

void op_par_loop_unstructured_stream_kernel_instrumented(
  char const *name, op_set set,
  op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3, op_arg arg4
  #ifdef VERIFY_OP2_TIMING
    , double* compute_time_ptr, double* sync_time_ptr
  #endif
  , long* iter_counts_ptr
  #ifdef PAPI
  , long_long* restrict event_counts, int event_set, int num_events
  #endif
  )
{

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc_manytime(12, omp_get_max_threads());
  op_timers_core(&cpu_t1, &wall_t1);
  double non_thread_walltime = 0.0;

  int  ninds   = 2;
  int  inds[5] = {0,0,-1,1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: unstructured_stream_kernel\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_12
    int part_size = OP_PART_SIZE_12;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set->size >0) {
    int nthreads = omp_get_max_threads();

    #ifdef PAPI
      // Init and start PAPI
      long_long* temp_count_stores = NULL;
      if (num_events > 0) {
        temp_count_stores = (long_long*)malloc(sizeof(long_long)*num_events);
        for (int e=0; e<num_events; e++) {
          temp_count_stores[e] = 0;
        }
      }
    #endif

    op_plan *Plan = op_plan_get_stage_upload(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL,0);

    // execute plan
    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all(nargs, args);
      }
      int nblocks = Plan->ncolblk[col];

      // Pause process timing and switch to per-thread timing:
      op_timers_core(&cpu_t2, &wall_t2);
      non_thread_walltime += wall_t2 - wall_t1;
      #pragma omp parallel
      {
        double thr_wall_t1, thr_wall_t2, thr_cpu_t1, thr_cpu_t2;
        op_timers_core(&thr_cpu_t1, &thr_wall_t1);

        int thr = omp_get_thread_num();
        int thr_start = (nblocks * thr) / nthreads;
        int thr_end = (nblocks * (thr+1)) / nthreads;
        if (thr_end > nblocks) thr_end = nblocks;

        #ifdef PAPI
          if (thr == 0) {
            if (num_events > 0) {
              my_papi_start(event_set);
            }
          }
        #endif

        for ( int blockIdx=thr_start; blockIdx<thr_end; blockIdx++ ){
          int blockId  = Plan->blkmap[blockIdx + block_offset];
          int nelem    = Plan->nelems[blockId];
          int offset_b = Plan->offset[blockId];
          #ifndef VECTORIZE
            for ( int n=offset_b; n<offset_b+nelem; n++ ){
              int map0idx = arg0.map_data[n * arg0.map->dim + 0];
              int map1idx = arg0.map_data[n * arg0.map->dim + 1];


              unstructured_stream_kernel(
                &((double*)arg0.data)[5 * map0idx],
                &((double*)arg0.data)[5 * map1idx],
                &((double*)arg2.data)[3 * n],
                &((double*)arg3.data)[5 * map0idx],
                &((double*)arg3.data)[5 * map1idx]);
            }
          #else
            int loop_size = nelem;
            int loop_end = offset_b+nelem;
            int simd_end = offset_b+(loop_size/SIMD_VEC)*SIMD_VEC;
            ALIGNED_double double dat0[5][SIMD_VEC];
            ALIGNED_double double dat1[5][SIMD_VEC];
            ALIGNED_double double dat2[3][SIMD_VEC];
            ALIGNED_double double dat3[5][SIMD_VEC];
            ALIGNED_double double dat4[5][SIMD_VEC];
            for (int o=offset_b ; o < simd_end; o+=SIMD_VEC) {

                #pragma omp simd simdlen(SIMD_VEC)
                // "sl" is SIMD lane:
                for (int sl=0; sl<SIMD_VEC; sl++ ){
                  int n = o+sl;
                  int map0idx = arg0.map_data[n * arg0.map->dim + 0];
                  int map1idx = arg0.map_data[n * arg0.map->dim + 1];

                  #pragma unroll(5)
                  for (int v=0; v<5; v++) {
                      dat0[v][sl] = ((double*)arg0.data)[map0idx*5 + v];
                      dat1[v][sl] = ((double*)arg0.data)[map1idx*5 + v];
                      dat3[v][sl] = 0.0;
                      dat4[v][sl] = 0.0;
                  }

                  for (int v=0; v<3; v++) {
                      dat2[v][sl] = ((double*)arg2.data)[n*3 + v];
                  }
                }

                #pragma omp simd simdlen(SIMD_VEC)
                for (int idx=0; idx<SIMD_VEC; idx++ ){
                  unstructured_stream_kernel_vec(
                    dat0,
                    dat1,
                    dat2,
                    dat3,
                    dat4,
                    idx);
                }

                #pragma omp simd safelen(1)
                for ( int sl=0; sl<SIMD_VEC; sl++ ){
                  int n = o+sl;
                  int map0idx = arg0.map_data[n * arg0.map->dim + 0];
                  int map1idx = arg0.map_data[n * arg0.map->dim + 1];

                  #pragma unroll(5)
                  for (int v=0; v<5; v++) {
                      ((double*)arg3.data)[map0idx*5 + v] += dat3[v][sl];
                      ((double*)arg3.data)[map1idx*5 + v] += dat4[v][sl];
                  }
                }
            }

            // remainder:
            for (int n = simd_end ; n < loop_end; n++) {
              int map0idx = arg0.map_data[n * arg0.map->dim + 0];
              int map1idx = arg0.map_data[n * arg0.map->dim + 1];

              unstructured_stream_kernel(
                &((double*)arg0.data)[5 * map0idx],
                &((double*)arg0.data)[5 * map1idx],
                &((double*)arg2.data)[3 * n],
                &((double*)arg3.data)[5 * map0idx],
                &((double*)arg3.data)[5 * map1idx]);
            }
          #endif
        }

        op_timers_core(&thr_cpu_t2, &thr_wall_t2);
        OP_kernels[12].times[thr]  += thr_wall_t2 - thr_wall_t1;

        #ifdef PAPI
          if (thr == 0) {
            if (num_events > 0) {
              my_papi_stop(event_counts, temp_count_stores, event_set, num_events);
              for (int e=0; e<num_events; e++) temp_count_stores[e] = 0;
            }
          }
        #endif
      }

      // Revert to process-level timing:
      op_timers_core(&cpu_t1, &wall_t1);

      block_offset += nblocks;
    }
    OP_kernels[12].transfer  += Plan->transfer;
    OP_kernels[12].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  non_thread_walltime += wall_t2 - wall_t1;
  OP_kernels[12].name      = name;
  OP_kernels[12].count    += 1;
  OP_kernels[12].times[0] += non_thread_walltime;
}
