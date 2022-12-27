//
// auto-generated by op2.py
//

#include <cmath>
#include "const.h"
#include "structures.h"
#include "global.h"

//user function
class initialize_variables_kernel_kernel;

//host stub function
void op_par_loop_initialize_variables_kernel(char const *name, op_set set,
  op_arg arg0){

  int nargs = 1;
  op_arg args[1];

  args[0] = arg0;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(0);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[0].name      = name;
  OP_kernels[0].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  initialize_variables_kernel\n");
  }

  int exec_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (exec_size > 0) {

    const int direct_initialize_variables_kernel_stride_OP2CONSTANT = getSetSizeFromOpArg(&arg0);
    //set SYCL execution parameters
    #ifdef OP_BLOCK_SIZE_0
      int nthread = OP_BLOCK_SIZE_0;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    double *arg0_d = (double*)arg0.data_d;
    int set_size = set->size+set->exec_size;
    try {
    op2_queue->wait();
    op2_queue->submit([&](cl::sycl::handler& cgh) {
      auto ff_variable_sycl = (*ff_variable_p).template get_access<cl::sycl::access::mode::read>(cgh);

      //user fun as lambda
      auto initialize_variables_kernel_gpu = [=]( 
            double* variables) {
            for(int j = 0; j < NVAR; j++) {
                variables[(j)*direct_initialize_variables_kernel_stride_OP2CONSTANT] = ff_variable_sycl[j];
            }
        
        };
        
      auto kern = [=](cl::sycl::item<1> item) [[intel::reqd_sub_group_size(SIMD_VEC)]] {

        //process set elements
        int n = item.get_id(0);
        if (n < exec_size) {

          //user-supplied kernel call
          initialize_variables_kernel_gpu(&arg0_d[n]);
        }

      };
      cgh.parallel_for<class initialize_variables_kernel_kernel>(cl::sycl::range<1>(set_size), kern);
    });
    }catch(cl::sycl::exception const &e) {
    std::cout << e.what() << std::endl;exit(-1);
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  op2_queue->wait();
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[0].time     += wall_t2 - wall_t1;
  OP_kernels[0].transfer += (float)set->size * arg0.size * 2.0f;
}
