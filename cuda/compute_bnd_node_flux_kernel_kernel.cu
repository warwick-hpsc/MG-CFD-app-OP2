//
// auto-generated by op2.py
//

#include <math.h>
#include "inlined_funcs.h"
#include "global.h"
#include "config.h"

//user function
__device__ void compute_bnd_node_flux_kernel_gpu( 
  const int *g,
  const double *edge_weight,
  const double *variables_b,
  double *fluxes_b) {








    if ((*g) <= 2) {

      
      
      
      
          double p_b = variables_b[VAR_DENSITY];
      
          #ifdef IDIVIDE
          double ip_b = 1.0 / p_b;
          #endif
      
          double pe_b, pressure_b;
          double3 velocity_b, momentum_b;
          double flux_contribution_i_momentum_x_b[NDIM],
                 flux_contribution_i_momentum_y_b[NDIM],
                 flux_contribution_i_momentum_z_b[NDIM],
                 flux_contribution_i_density_energy_b[NDIM];
      
          momentum_b.x = variables_b[VAR_MOMENTUM+0];
          momentum_b.y = variables_b[VAR_MOMENTUM+1];
          momentum_b.z = variables_b[VAR_MOMENTUM+2];
          pe_b = variables_b[VAR_DENSITY_ENERGY];
      
          #ifdef IDIVIDE
          compute_velocity(ip_b, momentum_b, velocity_b);
          #else
          compute_velocity(p_b, momentum_b, velocity_b);
          #endif
      
          double speed_sqd_b = compute_speed_sqd(velocity_b);
          double speed_b = std::sqrt(speed_sqd_b);
          pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);
      
          #ifdef IDIVIDE
          double speed_of_sound_b = compute_speed_of_sound(ip_b, pressure_b);
          #else
          double speed_of_sound_b = compute_speed_of_sound(p_b, pressure_b);
          #endif
      
          compute_flux_contribution(p_b, momentum_b, pe_b,
              pressure_b, velocity_b,
              flux_contribution_i_momentum_x_b,
              flux_contribution_i_momentum_y_b,
              flux_contribution_i_momentum_z_b,
              flux_contribution_i_density_energy_b);
      
          fluxes_b[VAR_DENSITY]        += 0;
          fluxes_b[VAR_MOMENTUM +0]    += edge_weight[0]*pressure_b;
          fluxes_b[VAR_MOMENTUM +1]    += edge_weight[1]*pressure_b;
          fluxes_b[VAR_MOMENTUM +2]    += edge_weight[2]*pressure_b;
          fluxes_b[VAR_DENSITY_ENERGY] += 0;
      

    } else if ((*g) == 3 || ((*g) >= 4 && (*g) <= 7) ) {


      
      
      
      
          double p_b = variables_b[VAR_DENSITY];
      
          #ifdef IDIVIDE
          double ip_b = 1.0 / p_b;
          #endif
      
          double pe_b, pressure_b;
          double3 velocity_b, momentum_b;
          double flux_contribution_i_momentum_x_b[NDIM],
                 flux_contribution_i_momentum_y_b[NDIM],
                 flux_contribution_i_momentum_z_b[NDIM],
                 flux_contribution_i_density_energy_b[NDIM];
      
          momentum_b.x = variables_b[VAR_MOMENTUM+0];
          momentum_b.y = variables_b[VAR_MOMENTUM+1];
          momentum_b.z = variables_b[VAR_MOMENTUM+2];
          pe_b = variables_b[VAR_DENSITY_ENERGY];
      
          #ifdef IDIVIDE
          compute_velocity(ip_b, momentum_b, velocity_b);
          #else
          compute_velocity(p_b, momentum_b, velocity_b);
          #endif
      
          double speed_sqd_b = compute_speed_sqd(velocity_b);
          double speed_b = std::sqrt(speed_sqd_b);
          pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);
      
          #ifdef IDIVIDE
          double speed_of_sound_b = compute_speed_of_sound(ip_b, pressure_b);
          #else
          double speed_of_sound_b = compute_speed_of_sound(p_b, pressure_b);
          #endif
      
          compute_flux_contribution(p_b, momentum_b, pe_b,
                                    pressure_b, velocity_b,
                                    flux_contribution_i_momentum_x_b,
                                    flux_contribution_i_momentum_y_b,
                                    flux_contribution_i_momentum_z_b,
                                    flux_contribution_i_density_energy_b);
      
          double factor_x = 0.5 * edge_weight[0],
                 factor_y = 0.5 * edge_weight[1],
                 factor_z = 0.5 * edge_weight[2];
      
          fluxes_b[VAR_DENSITY] +=
                factor_x*(ff_variable_cuda[VAR_MOMENTUM+0] + momentum_b.x)
              + factor_y*(ff_variable_cuda[VAR_MOMENTUM+1] + momentum_b.y)
              + factor_z*(ff_variable_cuda[VAR_MOMENTUM+2] + momentum_b.z);
      
          fluxes_b[VAR_DENSITY_ENERGY] += 
                factor_x*(ff_flux_contribution_density_energy_cuda[0] + flux_contribution_i_density_energy_b[0])
              + factor_y*(ff_flux_contribution_density_energy_cuda[1] + flux_contribution_i_density_energy_b[1])
              + factor_z*(ff_flux_contribution_density_energy_cuda[2] + flux_contribution_i_density_energy_b[2]);
      
          fluxes_b[VAR_MOMENTUM + 0] += 
                factor_x*(ff_flux_contribution_momentum_x_cuda[0] + flux_contribution_i_momentum_x_b[0])
              + factor_y*(ff_flux_contribution_momentum_x_cuda[1] + flux_contribution_i_momentum_x_b[1])
              + factor_z*(ff_flux_contribution_momentum_x_cuda[2] + flux_contribution_i_momentum_x_b[2]);
      
          fluxes_b[VAR_MOMENTUM + 1] += 
                factor_x*(ff_flux_contribution_momentum_y_cuda[0] + flux_contribution_i_momentum_y_b[0])
              + factor_y*(ff_flux_contribution_momentum_y_cuda[1] + flux_contribution_i_momentum_y_b[1])
              + factor_z*(ff_flux_contribution_momentum_y_cuda[2] + flux_contribution_i_momentum_y_b[2]);
      
          fluxes_b[VAR_MOMENTUM + 2] += 
                factor_x*(ff_flux_contribution_momentum_z_cuda[0] + flux_contribution_i_momentum_z_b[0])
              + factor_y*(ff_flux_contribution_momentum_z_cuda[1] + flux_contribution_i_momentum_z_b[1])
              + factor_z*(ff_flux_contribution_momentum_z_cuda[2] + flux_contribution_i_momentum_z_b[2]);
      
    }


}

// CUDA kernel function
__global__ void op_cuda_compute_bnd_node_flux_kernel(
  const double *__restrict ind_arg0,
  double *__restrict ind_arg1,
  const int *__restrict opDat2Map,
  const int *__restrict arg0,
  const double *__restrict arg1,
  int    block_offset,
  int   *blkmap,
  int   *offset,
  int   *nelems,
  int   *ncolors,
  int   *colors,
  int   nblocks,
  int   set_size) {
  double arg3_l[5];

  __shared__ int    nelems2, ncolor;
  __shared__ int    nelem, offset_b;

  extern __shared__ char shared[];

  if (blockIdx.x+blockIdx.y*gridDim.x >= nblocks) {
    return;
  }
  if (threadIdx.x==0) {

    //get sizes and shift pointers and direct-mapped data

    int blockId = blkmap[blockIdx.x + blockIdx.y*gridDim.x  + block_offset];

    nelem    = nelems[blockId];
    offset_b = offset[blockId];

    nelems2  = blockDim.x*(1+(nelem-1)/blockDim.x);
    ncolor   = ncolors[blockId];

  }
  __syncthreads(); // make sure all of above completed

  for ( int n=threadIdx.x; n<nelems2; n+=blockDim.x ){
    int col2 = -1;
    int map2idx;
    if (n<nelem) {
      //initialise local variables
      for ( int d=0; d<5; d++ ){
        arg3_l[d] = ZERO_double;
      }
      map2idx = opDat2Map[n + offset_b + set_size * 0];


      //user-supplied kernel call
      compute_bnd_node_flux_kernel_gpu(arg0+(n+offset_b)*1,
                                 arg1+(n+offset_b)*3,
                                 ind_arg0+map2idx*5,
                                 arg3_l);
      col2 = colors[n+offset_b];
    }

    //store local variables

    for ( int col=0; col<ncolor; col++ ){
      if (col2==col) {
        arg3_l[0] += ind_arg1[0+map2idx*5];
        arg3_l[1] += ind_arg1[1+map2idx*5];
        arg3_l[2] += ind_arg1[2+map2idx*5];
        arg3_l[3] += ind_arg1[3+map2idx*5];
        arg3_l[4] += ind_arg1[4+map2idx*5];
        ind_arg1[0+map2idx*5] = arg3_l[0];
        ind_arg1[1+map2idx*5] = arg3_l[1];
        ind_arg1[2+map2idx*5] = arg3_l[2];
        ind_arg1[3+map2idx*5] = arg3_l[3];
        ind_arg1[4+map2idx*5] = arg3_l[4];
      }
      __syncthreads();
    }
  }
}


//host stub function
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
  OP_kernels[nk].name      = name;
  OP_kernels[nk].count    += 1;


  int    ninds   = 2;
  int    inds[4] = {-1,-1,0,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: compute_bnd_node_flux_kernel\n");
  }

  //get plan
  #ifdef OP_PART_SIZE_10
    int part_size = OP_PART_SIZE_10;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set->size > 0) {

    op_plan *Plan = op_plan_get(name,set,part_size,nargs,args,ninds,inds);

    //execute plan

    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      #ifdef OP_BLOCK_SIZE_10
      int nthread = OP_BLOCK_SIZE_10;
      #else
      int nthread = OP_block_size;
      #endif

      dim3 nblocks = dim3(Plan->ncolblk[col] >= (1<<16) ? 65535 : Plan->ncolblk[col],
      Plan->ncolblk[col] >= (1<<16) ? (Plan->ncolblk[col]-1)/65535+1: 1, 1);
      if (Plan->ncolblk[col] > 0) {
        op_cuda_compute_bnd_node_flux_kernel<<<nblocks,nthread>>>(
        (double *)arg2.data_d,
        (double *)arg3.data_d,
        arg2.map_data_d,
        (int*)arg0.data_d,
        (double*)arg1.data_d,
        block_offset,
        Plan->blkmap,
        Plan->offset,
        Plan->nelems,
        Plan->nthrcol,
        Plan->thrcol,
        Plan->ncolblk[col],
        set->size+set->exec_size);

      }
      block_offset += Plan->ncolblk[col];
    }
    OP_kernels[nk].transfer  += Plan->transfer;
    OP_kernels[nk].transfer2 += Plan->transfer2;
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[nk].time     += wall_t2 - wall_t1;
}
