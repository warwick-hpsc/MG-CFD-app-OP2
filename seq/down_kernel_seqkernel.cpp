//
// auto-generated by op2.py
//

//user function
#include ".././src/Kernels/mg.h"

// host stub function
void op_par_loop_down_kernel(char const *name, op_set set,
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

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(21);
  op_timers_core(&cpu_t1, &wall_t1);
  double inner_cpu_t1, inner_cpu_t2, inner_wall_t1, inner_wall_t2;
  double compute_time=0.0, sync_time=0.0;

  if (OP_diags>2) {
    printf(" kernel routine with indirection: down_kernel\n");
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
      int map3idx = arg3.map_data[n * arg3.map->dim + 0];


      down_kernel(
        &((double*)arg0.data)[5 * n],
        &((double*)arg1.data)[5 * n],
        &((double*)arg2.data)[3 * n],
        &((double*)arg3.data)[5 * map3idx],
        &((double*)arg4.data)[3 * map3idx]);
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
  OP_kernels[21].name      = name;
  OP_kernels[21].count    += 1;
  OP_kernels[21].time     += wall_t2 - wall_t1;
  OP_kernels[21].transfer += (float)set->size * arg3.size;
  OP_kernels[21].transfer += (float)set->size * arg4.size;
  OP_kernels[21].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[21].transfer += (float)set->size * arg1.size;
  OP_kernels[21].transfer += (float)set->size * arg2.size;
  OP_kernels[21].transfer += (float)set->size * arg3.map->dim * 4.0f;
}
