//
// auto-generated by op2.py
//

//user function
#include ".././src/Kernels/mg.h"

// host stub function
void op_par_loop_down_v2_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  int nargs = 10;
  op_arg args[10];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(19);
  op_timers_core(&cpu_t1, &wall_t1);
  double inner_cpu_t1, inner_cpu_t2, inner_wall_t1, inner_wall_t2;
  double compute_time=0.0, sync_time=0.0;

  if (OP_diags>2) {
    printf(" kernel routine with indirection: down_v2_kernel\n");
  }

  op_timers_core(&inner_cpu_t1, &inner_wall_t1);
  int set_size = op_mpi_halo_exchanges(set, nargs, args);
  op_timers_core(&inner_cpu_t2, &inner_wall_t2);
  sync_time += inner_wall_t2 - inner_wall_t1;

  if (set->size >0) {

    op_timers_core(&inner_cpu_t1, &inner_wall_t1);
    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_timers_core(&inner_cpu_t2, &inner_wall_t2);
        compute_time += inner_wall_t2 - inner_wall_t1;
        op_mpi_wait_all(nargs, args);
        op_timers_core(&inner_cpu_t1, &inner_wall_t1);
        sync_time += inner_wall_t1 - inner_wall_t2;
      }
      int map0idx = arg0.map_data[n * arg0.map->dim + 0];
      int map1idx = arg0.map_data[n * arg0.map->dim + 1];
      int map2idx = arg2.map_data[n * arg2.map->dim + 0];
      int map3idx = arg2.map_data[n * arg2.map->dim + 1];

      down_v2_kernel(
        &((double*)arg0.data)[3 * map0idx],
        &((double*)arg0.data)[3 * map1idx],
        &((double*)arg2.data)[3 * map2idx],
        &((double*)arg2.data)[3 * map3idx],
        &((double*)arg4.data)[5 * map2idx],
        &((double*)arg4.data)[5 * map3idx],
        &((double*)arg6.data)[5 * map0idx],
        &((double*)arg6.data)[5 * map1idx],
        &((double*)arg8.data)[1 * map0idx],
        &((double*)arg8.data)[1 * map1idx]);
    }
    op_timers_core(&inner_cpu_t2, &inner_wall_t2);
    compute_time += inner_wall_t2 - inner_wall_t1;
  }

  op_timers_core(&inner_cpu_t1, &inner_wall_t1);
  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);
  op_timers_core(&inner_cpu_t2, &inner_wall_t2);
  sync_time += inner_wall_t2 - inner_wall_t1;

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[19].name      = name;
  OP_kernels[19].count    += 1;
  OP_kernels[19].time     += wall_t2 - wall_t1;
  OP_kernels[19].transfer += (float)set->size * arg0.size;
  OP_kernels[19].transfer += (float)set->size * arg2.size;
  OP_kernels[19].transfer += (float)set->size * arg4.size;
  OP_kernels[19].transfer += (float)set->size * arg6.size * 2.0f;
  OP_kernels[19].transfer += (float)set->size * arg8.size * 2.0f;
  OP_kernels[19].transfer += (float)set->size * arg0.map->dim * 4.0f;
  OP_kernels[19].transfer += (float)set->size * arg2.map->dim * 4.0f;
}
