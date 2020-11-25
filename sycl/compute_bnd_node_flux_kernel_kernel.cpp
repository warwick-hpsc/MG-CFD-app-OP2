//
// auto-generated by op2.py
//

#include <math.h>

#ifndef INLINED_FUNCS_H
#define INLINED_FUNCS_H
#include "const.h"
#include "structures.h"
#define DEBUGGABLE_ABORT fprintf(stderr, "%s:%d\n", __FILE__, __LINE__); fflush(stderr); fflush(stdout); exit(EXIT_FAILURE);
static inline OP_FUN_PREFIX void zero_array(int nelr, double* variables)
{
    #ifdef OMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int i = 0; i < nelr; i++)
    {
        variables[i] = 0.0;
    }
}
static inline OP_FUN_PREFIX void zero_edges(int nelr, edge* mx, edge* my, edge* mz, edge* p, edge* pe)
{
    #ifdef OMP
    #pragma omp parallel
    {
    #endif
        #ifdef VERITAS
            pdump_start_profile();
        #endif
        #ifdef OMP
        #pragma omp for schedule(static)
        #endif
        #ifdef SIMD
        #pragma simd vectorlength(SIMD)
        #else
        #pragma novector
        #endif
        for(int i = 0; i < nelr; i++)
        {
            mx[i].a = 0.0;
            mx[i].b = 0.0;
            my[i].a = 0.0;
            my[i].b = 0.0;
            mz[i].a = 0.0;
            mz[i].b = 0.0;
            p[i].a = 0.0;
            p[i].b = 0.0;
            pe[i].a = 0.0;
            pe[i].b = 0.0;
        }
        #ifdef VERITAS
            pdump_end_profile();
        #endif
    #ifdef OMP
    }
    #endif
}
static inline OP_FUN_PREFIX void compute_flux_contribution(
    double& density, double3& momentum, 
    double& density_energy, 
    double& pressure, 
    double3& velocity, 
    double* fc_momentum_x, 
    double* fc_momentum_y, 
    double* fc_momentum_z, 
    double* fc_density_energy)
{
    fc_momentum_x[0] = velocity.x*momentum.x + pressure;
    fc_momentum_x[1] = velocity.x*momentum.y;
    fc_momentum_x[2] = velocity.x*momentum.z;
    fc_momentum_y[0] = fc_momentum_x[1];
    fc_momentum_y[1] = velocity.y*momentum.y + pressure;
    fc_momentum_y[2] = velocity.y*momentum.z;
    fc_momentum_z[0] = fc_momentum_x[2];
    fc_momentum_z[1] = fc_momentum_y[2];
    fc_momentum_z[2] = velocity.z*momentum.z + pressure;
    double de_p = density_energy+pressure;
    fc_density_energy[0] = velocity.x*de_p;
    fc_density_energy[1] = velocity.y*de_p;
    fc_density_energy[2] = velocity.z*de_p;
}
static inline OP_FUN_PREFIX double compute_speed_sqd(double3& velocity)
{
    return velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z;
}
static inline OP_FUN_PREFIX double compute_pressure(double& density, double& density_energy, double& speed_sqd)
{
    return (double(GAMMA)-double(1.0))*(density_energy - double(0.5)*density*speed_sqd);
}
#ifdef IDIVIDE
    static inline OP_FUN_PREFIX void compute_velocity(double& inverse_density, double3& momentum, double3& velocity)
    {
        velocity.x = momentum.x*inverse_density;
        velocity.y = momentum.y*inverse_deinlined_funcs.hnsity;
        velocity.z = momentum.z*inverse_density;
    }
    static inline OP_FUN_PREFIX double compute_speed_of_sound(double& inverse_density, double& pressure)
    {
        return cl::sycl::sqrt((double(GAMMA)*pressure)*inverse_density);
    }
#else
    static OP_FUN_PREFIX inline void compute_velocity(double& density, double3& momentum, double3& velocity)
    {
        velocity.x = momentum.x / density;
        velocity.y = momentum.y / density;
        velocity.z = momentum.z / density;
    }
    static inline OP_FUN_PREFIX double compute_speed_of_sound(double& density, double& pressure)
    {
        return cl::sycl::sqrt(double(GAMMA)*pressure / density);
    }
#endif
#endif
#include "global.h"
#include "config.h"

//user function
class compute_bnd_node_flux_kernel_kernel;

//host stub function
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

  op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set->size > 0) {

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL);

    cl::sycl::buffer<double,1> *arg2_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg2.data_d);
    cl::sycl::buffer<double,1> *arg3_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg3.data_d);
    cl::sycl::buffer<int,1> *map2_buffer = static_cast<cl::sycl::buffer<int,1>*>((void*)arg2.map_data_d);
    cl::sycl::buffer<int,1> *arg0_buffer = static_cast<cl::sycl::buffer<int,1>*>((void*)arg0.data_d);
    cl::sycl::buffer<double,1> *arg1_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg1.data_d);
    cl::sycl::buffer<int,1> *blkmap_buffer = static_cast<cl::sycl::buffer<int,1>*>((void*)Plan->blkmap);
    cl::sycl::buffer<int,1> *offset_buffer = static_cast<cl::sycl::buffer<int,1>*>((void*)Plan->offset);
    cl::sycl::buffer<int,1> *nelems_buffer = static_cast<cl::sycl::buffer<int,1>*>((void*)Plan->nelems);
    cl::sycl::buffer<int,1> *ncolors_buffer = static_cast<cl::sycl::buffer<int,1>*>((void*)Plan->nthrcol);
    cl::sycl::buffer<int,1> *colors_buffer = static_cast<cl::sycl::buffer<int,1>*>((void*)Plan->thrcol);
    int set_size = set->size+set->exec_size;
    //execute plan

    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int nthread = SIMD_VEC;

      int nblocks = op2_queue->get_device().get_info<cl::sycl::info::device::max_compute_units>();
      int nblocks2 = Plan->ncolblk[col];
      if (Plan->ncolblk[col] > 0) {
        try {
        op2_queue->submit([&](cl::sycl::handler& cgh) {
          auto ind_arg0 = (*arg2_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
          auto ind_arg1 = (*arg3_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
          auto opDat2Map =  (*map2_buffer).template get_access<cl::sycl::access::mode::read>(cgh);
          auto blkmap    = (*blkmap_buffer).template get_access<cl::sycl::access::mode::read>(cgh);
          auto offset    = (*offset_buffer).template get_access<cl::sycl::access::mode::read>(cgh);
          auto nelems    = (*nelems_buffer).template get_access<cl::sycl::access::mode::read>(cgh);
          auto ncolors   = (*ncolors_buffer).template get_access<cl::sycl::access::mode::read>(cgh);
          auto colors    = (*colors_buffer).template get_access<cl::sycl::access::mode::read>(cgh);

          auto arg0 = (*arg0_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
          auto arg1 = (*arg1_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);

          auto ff_variable_sycl = (*ff_variable_p).template get_access<cl::sycl::access::mode::read>(cgh);
          auto ff_flux_contribution_momentum_x_sycl = (*ff_flux_contribution_momentum_x_p).template get_access<cl::sycl::access::mode::read>(cgh);
          auto ff_flux_contribution_momentum_y_sycl = (*ff_flux_contribution_momentum_y_p).template get_access<cl::sycl::access::mode::read>(cgh);
          auto ff_flux_contribution_momentum_z_sycl = (*ff_flux_contribution_momentum_z_p).template get_access<cl::sycl::access::mode::read>(cgh);
          auto ff_flux_contribution_density_energy_sycl = (*ff_flux_contribution_density_energy_p).template get_access<cl::sycl::access::mode::read>(cgh);

          //user fun as lambda
          auto compute_bnd_node_flux_kernel_gpu = [=]( 
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
                      double speed_b = cl::sycl::sqrt(speed_sqd_b);
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
                      double speed_b = cl::sycl::sqrt(speed_sqd_b);
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
                            factor_x*(ff_variable_sycl[VAR_MOMENTUM+0] + momentum_b.x)
                          + factor_y*(ff_variable_sycl[VAR_MOMENTUM+1] + momentum_b.y)
                          + factor_z*(ff_variable_sycl[VAR_MOMENTUM+2] + momentum_b.z);
                  
                      fluxes_b[VAR_DENSITY_ENERGY] += 
                            factor_x*(ff_flux_contribution_density_energy_sycl[0] + flux_contribution_i_density_energy_b[0])
                          + factor_y*(ff_flux_contribution_density_energy_sycl[1] + flux_contribution_i_density_energy_b[1])
                          + factor_z*(ff_flux_contribution_density_energy_sycl[2] + flux_contribution_i_density_energy_b[2]);
                  
                      fluxes_b[VAR_MOMENTUM + 0] += 
                            factor_x*(ff_flux_contribution_momentum_x_sycl[0] + flux_contribution_i_momentum_x_b[0])
                          + factor_y*(ff_flux_contribution_momentum_x_sycl[1] + flux_contribution_i_momentum_x_b[1])
                          + factor_z*(ff_flux_contribution_momentum_x_sycl[2] + flux_contribution_i_momentum_x_b[2]);
                  
                      fluxes_b[VAR_MOMENTUM + 1] += 
                            factor_x*(ff_flux_contribution_momentum_y_sycl[0] + flux_contribution_i_momentum_y_b[0])
                          + factor_y*(ff_flux_contribution_momentum_y_sycl[1] + flux_contribution_i_momentum_y_b[1])
                          + factor_z*(ff_flux_contribution_momentum_y_sycl[2] + flux_contribution_i_momentum_y_b[2]);
                  
                      fluxes_b[VAR_MOMENTUM + 2] += 
                            factor_x*(ff_flux_contribution_momentum_z_sycl[0] + flux_contribution_i_momentum_z_b[0])
                          + factor_y*(ff_flux_contribution_momentum_z_sycl[1] + flux_contribution_i_momentum_z_b[1])
                          + factor_z*(ff_flux_contribution_momentum_z_sycl[2] + flux_contribution_i_momentum_z_b[2]);
                  
                }
            
            
            };
            
          auto kern = [=](cl::sycl::nd_item<1> item) [[intel::reqd_sub_group_size(SIMD_VEC)]] {
            double arg3_l[5];


            //get sizes and shift pointers and direct-mapped data

            int blocksPerWG = (nblocks2-1)/item.get_group_range(0)+1;
            for ( int idx=item.get_group_linear_id()*blocksPerWG; idx<(item.get_group_linear_id()+1)*blocksPerWG && idx < nblocks2; idx++ ){
              int blockId = blkmap[idx + block_offset];

              int nelem    = nelems[blockId];
              int offset_b = offset[blockId];
              sycl::ONEAPI::sub_group sg = item.get_sub_group();

              int nelems2  = item.get_local_range()[0]*(1+(nelem-1)/item.get_local_range()[0]);
              int ncolor   = ncolors[blockId];


              for ( int n=item.get_local_id(0); n<nelems2; n+=item.get_local_range()[0] ){
                int col2 = -1;
                int map2idx;
                if (n<nelem) {
                  //initialise local variables
                  for ( int d=0; d<5; d++ ){
                    arg3_l[d] = ZERO_double;
                  }
                  map2idx = opDat2Map[n + offset_b + set_size * 0];


                  //user-supplied kernel call
                  compute_bnd_node_flux_kernel_gpu(&arg0[(n+offset_b)*1],
                                                   &arg1[(n+offset_b)*3],
                                                   &ind_arg0[map2idx*5],
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
                  sg.barrier();
                }
              }

            }
          };
          cgh.parallel_for<class compute_bnd_node_flux_kernel_kernel>(cl::sycl::nd_range<1>(nthread*nblocks,nthread), kern);
        });
        }catch(cl::sycl::exception const &e) {
        std::cout << e.what() << std::endl;exit(-1);
        }

      }
      block_offset += Plan->ncolblk[col];
    }
    OP_kernels[10].transfer  += Plan->transfer;
    OP_kernels[10].transfer2 += Plan->transfer2;
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  op2_queue->wait();
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[10].time     += wall_t2 - wall_t1;
}
