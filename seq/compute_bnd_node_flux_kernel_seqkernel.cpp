//
// auto-generated by op2.py
//

//user function
#include ".././src/Kernels/flux.h"

// host stub function
void op_par_loop_compute_bnd_node_flux_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];
  const int nk = 10;

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(nk);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: compute_bnd_node_flux_kernel\n");
  }

  int set_size = op_mpi_halo_exchanges(set, nargs, args);

  if (set->size >0) {

    for ( int n=0; n<set->core_size; n++ ){
      int map2idx = arg2.map_data[n * arg2.map->dim + 0];

      compute_bnd_node_flux_kernel(
        &((int*)arg0.data)[1 * n],
        &((double*)arg1.data)[3 * n],
        &((double*)arg2.data)[5 * map2idx],
        &((double*)arg3.data)[5 * map2idx]);
    }

    if (set_size > set->core_size) {
      op_mpi_wait_all(nargs, args);
      for ( int n=set->core_size; n<set_size; n++ ){
        int map2idx = arg2.map_data[n * arg2.map->dim + 0];

        compute_bnd_node_flux_kernel(
          &((int*)arg0.data)[1 * n],
          &((double*)arg1.data)[3 * n],
          &((double*)arg2.data)[5 * map2idx],
          &((double*)arg3.data)[5 * map2idx]);
      }
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[nk].name      = name;
  OP_kernels[nk].count    += 1;
  OP_kernels[nk].time     += wall_t2 - wall_t1;
  OP_kernels[nk].transfer += (float)set->size * arg2.size;
  OP_kernels[nk].transfer += (float)set->size * arg3.size * 2.0f;
  OP_kernels[nk].transfer += (float)set->size * arg0.size;
  OP_kernels[nk].transfer += (float)set->size * arg1.size;
  OP_kernels[nk].transfer += (float)set->size * arg2.map->dim * 4.0f;
}
