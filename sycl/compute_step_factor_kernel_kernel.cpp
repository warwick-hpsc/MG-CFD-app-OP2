//
// auto-generated by op2.py
//

#include <math.h>
#include <cmath>
#include "const.h"
#include "inlined_funcs.h"

//user function
class compute_step_factor_kernel_kernel;

//host stub function
void op_par_loop_compute_step_factor_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  double*arg2h = (double *)arg2.data;
  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(8);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[8].name      = name;
  OP_kernels[8].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  compute_step_factor_kernel\n");
  }

  op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set->size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    allocConstArrays(consts_bytes, "double");
    consts_bytes = 0;
    arg2.data   = OP_consts_h + consts_bytes;
    int arg2_offset = consts_bytes/sizeof(double);
    for ( int d=0; d<1; d++ ){
      ((double *)arg2.data)[d] = arg2h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    mvConstArraysToDevice(consts_bytes, "double");
    cl::sycl::buffer<double,1> *consts = static_cast<cl::sycl::buffer<double,1> *>((void*)OP_consts_d);

    //set SYCL execution parameters
    #ifdef OP_BLOCK_SIZE_8
      int nthread = OP_BLOCK_SIZE_8;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    cl::sycl::buffer<double,1> *arg0_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg0.data_d);
    cl::sycl::buffer<double,1> *arg1_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg1.data_d);
    cl::sycl::buffer<double,1> *arg3_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg3.data_d);
    int set_size = set->size+set->exec_size;
    try {
    op2_queue->submit([&](cl::sycl::handler& cgh) {
      auto arg0 = (*arg0_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
      auto arg1 = (*arg1_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
      auto arg3 = (*arg3_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
      auto consts_d = (*consts).template get_access<cl::sycl::access::mode::read_write>(cgh);

      //user fun as lambda
      auto compute_step_factor_kernel_gpu = [=]( 
            const double* variable,
            const double* volume,
            const double* min_dt,
            double* step_factor) {
            double density = variable[VAR_DENSITY];
        
            double3 momentum;
            momentum.x = variable[VAR_MOMENTUM+0];
            momentum.y = variable[VAR_MOMENTUM+1];
            momentum.z = variable[VAR_MOMENTUM+2];
        
            double density_energy = variable[VAR_DENSITY_ENERGY];
            double3 velocity; compute_velocity(density, momentum, velocity);
            double speed_sqd      = compute_speed_sqd(velocity);
            double pressure       = compute_pressure(density, density_energy, speed_sqd);
            double speed_of_sound = compute_speed_of_sound(density, pressure);
        
            *step_factor = (*min_dt) / (*volume);
        
        };
        
      auto kern = [=](cl::sycl::nd_item<1> item) {

        //process set elements
        for ( int n=item.get_global_linear_id(); n<set_size; n+=item.get_global_range()[0] ){

          //user-supplied kernel call
          compute_step_factor_kernel_gpu(&arg0[n*5],
                               &arg1[n*1],
                               &consts_d[arg2_offset],
                               &arg3[n*1]);
        }

      };
      cgh.parallel_for<class compute_step_factor_kernel_kernel>(cl::sycl::nd_range<1>(nthread*nblocks,nthread), kern);
    });
    }catch(cl::sycl::exception const &e) {
    std::cout << e.what() << std::endl;exit(-1);
    }
    freeConstArrays("double");
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  op2_queue->wait();
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[8].time     += wall_t2 - wall_t1;
  OP_kernels[8].transfer += (float)set->size * arg0.size;
  OP_kernels[8].transfer += (float)set->size * arg1.size;
  OP_kernels[8].transfer += (float)set->size * arg3.size * 2.0f;
}
