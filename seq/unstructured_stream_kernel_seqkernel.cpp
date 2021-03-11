//
// auto-generated by op2.py
//

//user function
#include ".././src/Kernels/unstructured_stream.h"

#ifdef PAPI
#include "papi_funcs.h"
#endif

#include <mpi.h>

// host stub function
void op_par_loop_unstructured_stream_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){
  
  op_par_loop_unstructured_stream_kernel_instrumented(name, set, 
    arg0, arg1, arg2, arg3, arg4
    #ifdef PAPI
    , NULL
    #endif
    );
};

void op_par_loop_unstructured_stream_kernel_instrumented(
  char const *name, op_set set,
  op_arg arg0, op_arg arg1, op_arg arg2, op_arg arg3, op_arg arg4
  #ifdef PAPI
  , long_long** restrict event_counts
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
  op_timing_realloc(12);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: unstructured_stream_kernel\n");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set_size >0) {
    #ifdef MEASURE_MEM_BW
      // Need to ensure that MPI complete before timing. 
      // Not necessary to insert an explicit barrier, at least for 
      // single node benchmarking.
      op_mpi_wait_all(nargs, args);
      op_timers_core(&cpu_t1, &wall_t1);
    #endif

    #ifdef PAPI
      my_papi_start();
    #endif

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        #ifdef PAPI
          my_papi_stop(event_counts);
        #endif
        op_mpi_wait_all(nargs, args);
        #ifdef PAPI
          my_papi_start();
        #endif
      }
      int map0idx = arg0.map_data[n * arg0.map->dim + 0];
      int map1idx = arg0.map_data[n * arg0.map->dim + 1];


      unstructured_stream_kernel(
        &((double*)arg0.data)[5 * map0idx],
        &((double*)arg0.data)[5 * map1idx],
        &((double*)arg2.data)[3 * n],
        &((double*)arg3.data)[5 * map0idx],
        &((double*)arg3.data)[5 * map1idx]);
    }

    #ifdef PAPI
      my_papi_stop(event_counts);
    #endif
  }

  if (set_size == 0 || set_size == set->core_size) {
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
