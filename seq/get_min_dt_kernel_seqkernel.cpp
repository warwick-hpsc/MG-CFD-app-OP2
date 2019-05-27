//
// auto-generated by op2.py
//

//user function
#include ".././src/Kernels/compute_step_factor_kernel.h"

// host stub function
void op_par_loop_get_min_dt_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(7);
  op_timers_core(&cpu_t1, &wall_t1);
  double inner_cpu_t1, inner_cpu_t2, inner_wall_t1, inner_wall_t2;
  double compute_time=0.0, sync_time=0.0;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  get_min_dt_kernel");
  }

  op_timers_core(&inner_cpu_t1, &inner_wall_t1);
  int set_size = op_mpi_halo_exchanges(set, nargs, args);
  op_timers_core(&inner_cpu_t2, &inner_wall_t2);
  sync_time += inner_wall_t2 - inner_wall_t1;

  if (set->size >0) {

    op_timers_core(&inner_cpu_t1, &inner_wall_t1);
    for ( int n=0; n<set_size; n++ ){
      get_min_dt_kernel(
        &((double*)arg0.data)[1*n],
        (double*)arg1.data);
    }
    op_timers_core(&inner_cpu_t2, &inner_wall_t2);
    compute_time += inner_wall_t2 - inner_wall_t1;
  }

  op_timers_core(&inner_cpu_t1, &inner_wall_t1);
  // combine reduction data
  op_mpi_reduce_double(&arg1,(double*)arg1.data);
  op_mpi_set_dirtybit(nargs, args);
  op_timers_core(&inner_cpu_t2, &inner_wall_t2);
  sync_time += inner_wall_t2 - inner_wall_t1;

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[7].name      = name;
  OP_kernels[7].count    += 1;
  OP_kernels[7].time     += wall_t2 - wall_t1;
  OP_kernels[7].transfer += (float)set->size * arg0.size;
}
