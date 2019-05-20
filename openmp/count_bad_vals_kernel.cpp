//
// auto-generated by op2.py
//

//user function
#include ".././src/Kernels/validation.h"

// host stub function
void op_par_loop_count_bad_vals(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int*arg1h = (int *)arg1.data;
  int nargs = 2;
  op_arg args[2];
  const int nk = 14;

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(nk);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  count_bad_vals");
  }

  op_mpi_halo_exchanges(set, nargs, args);
  // set number of threads
  #ifdef _OPENMP
    int nthreads = omp_get_max_threads();
  #else
    int nthreads = 1;
  #endif

  // allocate and initialise arrays for global reduction
  int arg1_l[nthreads*64];
  for ( int thr=0; thr<nthreads; thr++ ){
    for ( int d=0; d<1; d++ ){
      arg1_l[d+thr*64]=ZERO_int;
    }
  }

  if (set->size >0) {

    // execute plan
    #pragma omp parallel for
    for ( int thr=0; thr<nthreads; thr++ ){
      int start  = (set->size* thr)/nthreads;
      int finish = (set->size*(thr+1))/nthreads;
      for ( int n=start; n<finish; n++ ){
        count_bad_vals(
          &((double*)arg0.data)[5*n],
          &arg1_l[64*omp_get_thread_num()]);
      }
    }
  }

  // combine reduction data
  for ( int thr=0; thr<nthreads; thr++ ){
    for ( int d=0; d<1; d++ ){
      arg1h[d] += arg1_l[d+thr*64];
    }
  }
  op_mpi_reduce(&arg1,arg1h);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[nk].name      = name;
  OP_kernels[nk].count    += 1;
  OP_kernels[nk].time     += wall_t2 - wall_t1;
  OP_kernels[nk].transfer += (float)set->size * arg0.size;
}
