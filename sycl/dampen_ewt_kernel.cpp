//
// auto-generated by op2.py
//

#include <cmath>
#include "const.h"
#include "structures.h"
#include "global.h"

//user function
class dampen_ewt_kernel;

//host stub function
void op_par_loop_dampen_ewt(char const *name, op_set set,
  op_arg arg0){

  int nargs = 1;
  op_arg args[1];

  args[0] = arg0;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(4);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[4].name      = name;
  OP_kernels[4].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  dampen_ewt\n");
  }

  op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (set->size > 0) {

    //set SYCL execution parameters
    #ifdef OP_BLOCK_SIZE_4
      int nthread = OP_BLOCK_SIZE_4;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    cl::sycl::buffer<double,1> *arg0_buffer = static_cast<cl::sycl::buffer<double,1>*>((void*)arg0.data_d);
    int set_size = set->size+set->exec_size;
    try {
    op2_queue->submit([&](cl::sycl::handler& cgh) {
      auto arg0 = (*arg0_buffer).template get_access<cl::sycl::access::mode::read_write>(cgh);

      //user fun as lambda
      auto dampen_ewt_gpu = [=]( 
            double* ewt) {
            ewt[0] *= 1e-7;
            ewt[1] *= 1e-7;
            ewt[2] *= 1e-7;
        
        };
        
      auto kern = [=](cl::sycl::item<1> item) {

        //process set elements
        int n = item.get_id(0);
        if (n < set_size) {

          //user-supplied kernel call
          dampen_ewt_gpu(&arg0[n*3]);
        }

      };
      cgh.parallel_for<class dampen_ewt_kernel>(cl::sycl::range<1>(set_size), kern);
    });
    }catch(cl::sycl::exception const &e) {
    std::cout << e.what() << std::endl;exit(-1);
    }
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  op2_queue->wait();
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[4].time     += wall_t2 - wall_t1;
  OP_kernels[4].transfer += (float)set->size * arg0.size * 2.0f;
}
