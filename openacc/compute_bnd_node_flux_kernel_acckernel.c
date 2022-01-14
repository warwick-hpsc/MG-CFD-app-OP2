//
// auto-generated by op2.py
//

//user function
#include <math.h>
#include "inlined_funcs.h"
#include "global.h"
#include "config.h"

int opDat2_compute_bnd_node_flux_kernel_stride_OP2CONSTANT;
int opDat2_compute_bnd_node_flux_kernel_stride_OP2HOST=-1;
int direct_compute_bnd_node_flux_kernel_stride_OP2CONSTANT;
int direct_compute_bnd_node_flux_kernel_stride_OP2HOST=-1;
//user function
//#pragma acc routine
inline void compute_bnd_node_flux_kernel_openacc( 
  const int *g,
  const double *edge_weight,
  const double *variables_b,
  double *fluxes_b) {








    if ((*g) <= 2) {

      #include "flux_boundary.elem_func"

    } else if ((*g) == 3 || ((*g) >= 4 && (*g) <= 7) ) {


      #include "flux_wall.elem_func"
    }

}

// host stub function
void op_par_loop_compute_bnd_node_flux_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(10);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[10].name      = name;
  OP_kernels[10].count    += 1;

  int  ninds   = 2;
  int  inds[4] = {-1,-1,0,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: compute_bnd_node_flux_kernel\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_10
    int part_size = OP_PART_SIZE_10;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {

    if ((OP_kernels[10].count==1) || (opDat2_compute_bnd_node_flux_kernel_stride_OP2HOST != getSetSizeFromOpArg(&arg2))) {
      opDat2_compute_bnd_node_flux_kernel_stride_OP2HOST = getSetSizeFromOpArg(&arg2);
      opDat2_compute_bnd_node_flux_kernel_stride_OP2CONSTANT = opDat2_compute_bnd_node_flux_kernel_stride_OP2HOST;
    }
    if ((OP_kernels[10].count==1) || (direct_compute_bnd_node_flux_kernel_stride_OP2HOST != getSetSizeFromOpArg(&arg1))) {
      direct_compute_bnd_node_flux_kernel_stride_OP2HOST = getSetSizeFromOpArg(&arg1);
      direct_compute_bnd_node_flux_kernel_stride_OP2CONSTANT = direct_compute_bnd_node_flux_kernel_stride_OP2HOST;
    }

    //Set up typed device pointers for OpenACC
    int *map2 = arg2.map_data_d;

    int* data0 = (int*)arg0.data_d;
    double* data1 = (double*)arg1.data_d;
    double *data2 = (double *)arg2.data_d;
    double *data3 = (double *)arg3.data_d;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;
    int set_size1 = set->size + set->exec_size;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      #pragma acc parallel loop independent deviceptr(col_reord,map2,data0,data1,data2,data3)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map2idx;
        map2idx = map2[n + set_size1 * 0];


        compute_bnd_node_flux_kernel_openacc(
          &data0[1 * n],
          &data1[n],
          &data2[map2idx],
          &data3[map2idx]);
      }

    }
    OP_kernels[10].transfer  += Plan->transfer;
    OP_kernels[10].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[10].time     += wall_t2 - wall_t1;
}
