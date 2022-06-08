
#ifndef COMM_AVOID_H
#define COMM_AVOID_H
#ifdef COMM_AVOID


#include "op_lib_mpi.h"
#include "op_mpi_core.h"
#include "op_lib_core.h"

void ca_op_par_loop_test_write_kernel(char const *, op_set,
  op_arg,
  op_arg,
  int, int, int);

void ca_op_par_loop_test_read_kernel(char const *, op_set,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  op_arg,
  int, int, int);

void test_comm_avoid(char const *name, op_dat* p_variables, op_dat p_edge_weights, op_dat p_fluxes, op_map p_edge_to_nodes, op_set set,
                      int nloops, int nchains, int default_variable_index){
   
    int nhalos = nloops; // * nchains;
    int map_index = nhalos;

    int nargs_ex0 = nchains;
    op_arg args_ex0[nchains];

    for(int i = 0; i < nchains; i++){
      args_ex0[i] = op_arg_dat_halo(p_variables[i],0,p_edge_to_nodes,5,"double",OP_READ, nhalos, map_index);
      set_dat_dirty(&args_ex0[i]);
    }

    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    op_timing_realloc(28);
    op_timers_core(&cpu_t1, &wall_t1);
    op_mpi_halo_exchanges_cuda_chained(set, nargs_ex0, args_ex0, nhalos, 1);

    int nargs0 = 2;
    op_arg args0[nchains][2];

    int nargs1 = 5;
    op_arg args1[nchains][5];

    for(int i = 0; i < nchains; i++){
      args0[i][0] = op_arg_dat_halo(p_variables[i],0,p_edge_to_nodes,5,"double",OP_INC, nhalos, map_index);
      args0[i][1] = op_arg_dat_halo(p_variables[i],1,p_edge_to_nodes,5,"double",OP_INC, nhalos, map_index);

      args1[i][0] = op_arg_dat_halo(p_variables[i],0,p_edge_to_nodes,5,"double",OP_READ, nhalos, map_index);
      args1[i][1] = op_arg_dat_halo(p_variables[i],1,p_edge_to_nodes,5,"double",OP_READ, nhalos, map_index);
      args1[i][2] = op_arg_dat_halo(p_edge_weights,-1,OP_ID,3,"double",OP_READ, nhalos, map_index);
      args1[i][3] = op_arg_dat_halo(p_fluxes,0,p_edge_to_nodes,5,"double",OP_INC, nhalos, map_index);
      args1[i][4] = op_arg_dat_halo(p_fluxes,1,p_edge_to_nodes,5,"double",OP_INC, nhalos, map_index);
    }

    // this will do the latency hiding for the first two loops in the chain, since halo extension is done for two levels.
    // other levels will receive 0 as core size.
    int n_lower0 = get_set_core_size(set, nloops - 1);
    int n_lower1 = get_set_core_size(set, nloops);

    for(int i = 0; i < nchains; i++){
      ca_op_par_loop_test_write_kernel("ca_test_write_kernel",set,
                            args0[i][0], args0[i][1], 0, n_lower0, 1);

      ca_op_par_loop_test_read_kernel("ca_test_read_kernel",set,
                            args1[i][0], args1[i][1], args1[i][2], args1[i][3], args1[i][4], 0, n_lower1, 1);
    }

    op_mpi_wait_all_cuda_chained(nargs_ex0, args_ex0);

    int n_upper0 = get_set_size_with_nhalos(set, nloops);
    int n_upper1 = get_set_size_with_nhalos(set, nloops - 1);

    for(int i = 0; i < nchains; i++){
      ca_op_par_loop_test_write_kernel("ca_test_write_kernel",set,
                            args0[i][0], args0[i][1], n_lower0, n_upper0, 0);
      
      ca_op_par_loop_test_read_kernel("ca_test_read_kernel",set,
                            args1[i][0], args1[i][1], args1[i][2], args1[i][3], args1[i][4], n_lower1, n_upper1, 0);
    }

    op_mpi_set_dirtybit_cuda(nargs1, args1[default_variable_index]);

    op_timers_core(&cpu_t2, &wall_t2);
    OP_kernels[28].name      = name;
    OP_kernels[28].count    += 1;
    OP_kernels[28].time     += wall_t2 - wall_t1;
    OP_kernels[28].transfer += (float)set->size * args_ex0[default_variable_index].size;
    OP_kernels[28].transfer += (float)set->size * args_ex0[default_variable_index].size * 2.0f;
    OP_kernels[28].transfer += (float)set->size * args_ex0[default_variable_index].size;
    OP_kernels[28].transfer += (float)set->size * args_ex0[default_variable_index].map->dim * 4.0f;
   
}
#endif
#endif