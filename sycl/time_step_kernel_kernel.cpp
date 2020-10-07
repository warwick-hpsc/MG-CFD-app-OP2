//
// auto-generated by op2.py
//

#include <math.h>
#include <cmath>
#include "const.h"
#include "inlined_funcs.h"

//user function
class time_step_kernel_kernel;

//host stub function
void op_par_loop_time_step_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  int*arg0h = (int *)arg0.data;
  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(11);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[11].name      = name;
  OP_kernels[11].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  time_step_kernel\n");
  }

  op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set->size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(int));
    allocConstArrays(consts_bytes, "int");
    consts_bytes = 0;
    arg0.data   = OP_consts_h + consts_bytes;
    int arg0_offset = consts_bytes/sizeof(int);
    for ( int d=0; d<1; d++ ){
      ((int *)arg0.data)[d] = arg0h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(int));
    mvConstArraysToDevice(consts_bytes, "int");
    cl::sycl::buffer<int,1> *consts = static_cast<cl::sycl::buffer<int,1> *>((void*)OP_consts_d);

    //set SYCL execution parameters
    #ifdef OP_BLOCK_SIZE_11
      int nthread = OP_BLOCK_SIZE_11;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    cl::sycl::buffer<double,1> *arg1_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg1.data_d);
    cl::sycl::buffer<double,1> *arg2_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg2.data_d);
    cl::sycl::buffer<double,1> *arg3_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg3.data_d);
    cl::sycl::buffer<double,1> *arg4_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg4.data_d);
    int set_size = set->size+set->exec_size;
    try {
    op2_queue->submit([&](cl::sycl::handler& cgh) {
      auto arg1 = (*arg1_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
      auto arg2 = (*arg2_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
      auto arg3 = (*arg3_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
      auto arg4 = (*arg4_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);
      auto consts_d = (*consts).template get_access<cl::sycl::access::mode::read_write>(cgh);

      //user fun as lambda
      auto time_step_kernel_gpu = [=]( 
            const int* rkCycle,
            const double* step_factor,
            double* flux,
            const double* old_variable,
            double* variable) {
            double factor = (*step_factor)/double(RK+1-(*rkCycle));
        
            variable[VAR_DENSITY]        = old_variable[VAR_DENSITY]        + factor*flux[VAR_DENSITY];
            variable[VAR_MOMENTUM+0]     = old_variable[VAR_MOMENTUM+0]     + factor*flux[VAR_MOMENTUM+0];
            variable[VAR_MOMENTUM+1]     = old_variable[VAR_MOMENTUM+1]     + factor*flux[VAR_MOMENTUM+1];
            variable[VAR_MOMENTUM+2]     = old_variable[VAR_MOMENTUM+2]     + factor*flux[VAR_MOMENTUM+2];
            variable[VAR_DENSITY_ENERGY] = old_variable[VAR_DENSITY_ENERGY] + factor*flux[VAR_DENSITY_ENERGY];
        
            flux[VAR_DENSITY]        = 0.0;
            flux[VAR_MOMENTUM+0]     = 0.0;
            flux[VAR_MOMENTUM+1]     = 0.0;
            flux[VAR_MOMENTUM+2]     = 0.0;
            flux[VAR_DENSITY_ENERGY] = 0.0;
        
        };
        
      auto kern = [=](cl::sycl::item<1> item) {

        //process set elements
        int n = item.get_id(0);
        if (n < set_size) {

          //user-supplied kernel call
          time_step_kernel_gpu(&consts_d[arg0_offset],
                               &arg1[n*1],
                               &arg2[n*5],
                               &arg3[n*5],
                               &arg4[n*5]);
        }

      };
      cgh.parallel_for<class time_step_kernel_kernel>(cl::sycl::range<1>(set_size), kern);
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
  OP_kernels[11].time     += wall_t2 - wall_t1;
  OP_kernels[11].transfer += (float)set->size * arg1.size;
  OP_kernels[11].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[11].transfer += (float)set->size * arg3.size;
  OP_kernels[11].transfer += (float)set->size * arg4.size * 2.0f;
}
