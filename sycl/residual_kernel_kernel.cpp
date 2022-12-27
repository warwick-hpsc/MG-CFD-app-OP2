//
// auto-generated by op2.py
//

#include "utils.h"

//user function
class residual_kernel_kernel;

//host stub function
void op_par_loop_residual_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2){

  int nargs = 3;
  op_arg args[3];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(13);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[13].name      = name;
  OP_kernels[13].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  residual_kernel\n");
  }

  int exec_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (exec_size > 0) {

    const int direct_residual_kernel_stride_OP2CONSTANT = getSetSizeFromOpArg(&arg0);
    //set SYCL execution parameters
    #ifdef OP_BLOCK_SIZE_13
      int nthread = OP_BLOCK_SIZE_13;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    double *arg0_d = (double*)arg0.data_d;
    double *arg1_d = (double*)arg1.data_d;
    double *arg2_d = (double*)arg2.data_d;
    int set_size = set->size+set->exec_size;
    try {
    op2_queue->wait();
    op2_queue->submit([&](cl::sycl::handler& cgh) {

      //user fun as lambda
      auto residual_kernel_gpu = [=]( 
            const double* old_variable,
            const double* variable,
            double* residual) {
            for (int v=0; v<NVAR; v++) {
                residual[(v)*direct_residual_kernel_stride_OP2CONSTANT] = variable[(v)*direct_residual_kernel_stride_OP2CONSTANT] - old_variable[(v)*direct_residual_kernel_stride_OP2CONSTANT];
            }
        
        };
        
      auto kern = [=](cl::sycl::item<1> item) [[intel::reqd_sub_group_size(SIMD_VEC)]] {

        //process set elements
        int n = item.get_id(0);
        if (n < exec_size) {

          //user-supplied kernel call
          residual_kernel_gpu(&arg0_d[n],
                              &arg1_d[n],
                              &arg2_d[n]);
        }

      };
      cgh.parallel_for<class residual_kernel_kernel>(cl::sycl::range<1>(set_size), kern);
    });
    }catch(cl::sycl::exception const &e) {
    std::cout << e.what() << std::endl;exit(-1);
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  op2_queue->wait();
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[13].time     += wall_t2 - wall_t1;
  OP_kernels[13].transfer += (float)set->size * arg0.size;
  OP_kernels[13].transfer += (float)set->size * arg1.size;
  OP_kernels[13].transfer += (float)set->size * arg2.size * 2.0f;
}
