//
// auto-generated by op2.py
//

//user function
int direct_calc_rms_kernel_stride_OP2CONSTANT;
int direct_calc_rms_kernel_stride_OP2HOST=-1;
//user function

void calc_rms_kernel_omp4_kernel(
  double *data0,
  int dat0size,
  double *arg1,
  int count,
  int num_teams,
  int nthread,
  int direct_calc_rms_kernel_stride_OP2CONSTANT);

// host stub function
void op_par_loop_calc_rms_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  double*arg1h = (double *)arg1.data;
  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(14);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[14].name      = name;
  OP_kernels[14].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  calc_rms_kernel");
  }

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);

  #ifdef OP_PART_SIZE_14
    int part_size = OP_PART_SIZE_14;
  #else
    int part_size = OP_part_size;
  #endif
  #ifdef OP_BLOCK_SIZE_14
    int nthread = OP_BLOCK_SIZE_14;
  #else
    int nthread = OP_block_size;
  #endif

  double arg1_l = arg1h[0];

  if (set_size >0) {

    if ((OP_kernels[14].count==1) || (direct_calc_rms_kernel_stride_OP2HOST != getSetSizeFromOpArg(&arg0))) {
      direct_calc_rms_kernel_stride_OP2HOST = getSetSizeFromOpArg(&arg0);
      direct_calc_rms_kernel_stride_OP2CONSTANT = direct_calc_rms_kernel_stride_OP2HOST;
    }

    //Set up typed device pointers for OpenMP

    double* data0 = (double*)arg0.data_d;
    int dat0size = getSetSizeFromOpArg(&arg0) * arg0.dat->dim;
    calc_rms_kernel_omp4_kernel(
      data0,
      dat0size,
      &arg1_l,
      set->size,
      part_size!=0?(set->size-1)/part_size+1:(set->size-1)/nthread,
      nthread,
      direct_calc_rms_kernel_stride_OP2CONSTANT);

  }

  // combine reduction data
  arg1h[0] = arg1_l;
  op_mpi_reduce_double(&arg1,arg1h);
  op_mpi_set_dirtybit_cuda(nargs, args);

  if (OP_diags>1) deviceSync();
  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[14].time     += wall_t2 - wall_t1;
  OP_kernels[14].transfer += (float)set->size * arg0.size;
}
