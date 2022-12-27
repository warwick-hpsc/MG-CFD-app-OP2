//
// auto-generated by op2.py
//

#include <math.h>
#include "const.h"

//user function
class down_v2_kernel_kernel;

//host stub function
void op_par_loop_down_v2_kernel(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7,
  op_arg arg8,
  op_arg arg9){

  int nargs = 10;
  op_arg args[10];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;
  args[8] = arg8;
  args[9] = arg9;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(20);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[20].name      = name;
  OP_kernels[20].count    += 1;


  int    ninds   = 5;
  int    inds[10] = {0,0,1,1,2,2,3,3,4,4};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: down_v2_kernel\n");
  }

  //get plan
  #ifdef OP_PART_SIZE_20
    int part_size = OP_PART_SIZE_20;
  #else
    int part_size = OP_part_size;
  #endif

  int exec_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (exec_size > 0) {

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL);

    const int opDat0_down_v2_kernel_stride_OP2CONSTANT = getSetSizeFromOpArg(&arg0);
    const int opDat2_down_v2_kernel_stride_OP2CONSTANT = getSetSizeFromOpArg(&arg2);
    double *ind_arg0 = (double*)arg0.data_d;
    double *ind_arg1 = (double*)arg2.data_d;
    double *ind_arg2 = (double*)arg4.data_d;
    double *ind_arg3 = (double*)arg6.data_d;
    double *ind_arg4 = (double*)arg8.data_d;
    int *opDat0Map = arg0.map_data_d;
    int *opDat2Map = arg2.map_data_d;
    int *blkmap = (int *)Plan->blkmap;
    int *offset = (int *)Plan->offset;
    int *nelems = (int *)Plan->nelems;
    int set_size = set->size+set->exec_size;
    //execute plan

    int block_offset = 0;
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int nthread = 1;

      int nblocks = op2_queue->get_device().get_info<cl::sycl::info::device::max_compute_units>();
      int nblocks2 = Plan->ncolblk[col];
      if (Plan->ncolblk[col] > 0) {
        try {
        op2_queue->wait();
        op2_queue->submit([&](cl::sycl::handler& cgh) {


          //user fun as lambda
          auto down_v2_kernel_gpu = [=]( 
                const double* coord2a,
                const double* coord2b,
                const double* coord1a,
                const double* coord1b,
                const double* residuals1a,
                const double* residuals1b,
                double* residuals1a_prolonged,
                double* residuals1b_prolonged,
                double* residuals1a_prolonged_wsum,
                double* residuals1b_prolonged_wsum) {
            
            
            
            
            
            
            
                double dx_a1a2 = coord2a[(0)*opDat0_down_v2_kernel_stride_OP2CONSTANT] - coord1a[(0)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                double dy_a1a2 = coord2a[(1)*opDat0_down_v2_kernel_stride_OP2CONSTANT] - coord1a[(1)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                double dz_a1a2 = coord2a[(2)*opDat0_down_v2_kernel_stride_OP2CONSTANT] - coord1a[(2)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                if (dx_a1a2 == 0.0 && dy_a1a2 == 0.0 && dz_a1a2 == 0.0) {
            
                    residuals1a_prolonged[(VAR_DENSITY)*opDat0_down_v2_kernel_stride_OP2CONSTANT]        = residuals1a[(VAR_DENSITY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+0)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     = residuals1a[(VAR_MOMENTUM+0)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+1)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     = residuals1a[(VAR_MOMENTUM+1)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+2)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     = residuals1a[(VAR_MOMENTUM+2)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_DENSITY_ENERGY)*opDat0_down_v2_kernel_stride_OP2CONSTANT] = residuals1a[(VAR_DENSITY_ENERGY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    *residuals1a_prolonged_wsum = 1.0;
                } else {
            
                    const double idist_a1a2 = 1.0/cl::sycl::sqrt(dx_a1a2*dx_a1a2 + dy_a1a2*dy_a1a2 + dz_a1a2*dz_a1a2);
                    residuals1a_prolonged[(VAR_DENSITY)*opDat0_down_v2_kernel_stride_OP2CONSTANT]        += idist_a1a2*residuals1a[(VAR_DENSITY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+0)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_a1a2*residuals1a[(VAR_MOMENTUM+0)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+1)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_a1a2*residuals1a[(VAR_MOMENTUM+1)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+2)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_a1a2*residuals1a[(VAR_MOMENTUM+2)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_DENSITY_ENERGY)*opDat0_down_v2_kernel_stride_OP2CONSTANT] += idist_a1a2*residuals1a[(VAR_DENSITY_ENERGY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    *residuals1a_prolonged_wsum += idist_a1a2;
            
                    double dx_b1a2 = coord1b[(0)*opDat2_down_v2_kernel_stride_OP2CONSTANT] - coord2a[(0)*opDat0_down_v2_kernel_stride_OP2CONSTANT];
                    double dy_b1a2 = coord1b[(1)*opDat2_down_v2_kernel_stride_OP2CONSTANT] - coord2a[(1)*opDat0_down_v2_kernel_stride_OP2CONSTANT];
                    double dz_b1a2 = coord1b[(2)*opDat2_down_v2_kernel_stride_OP2CONSTANT] - coord2a[(2)*opDat0_down_v2_kernel_stride_OP2CONSTANT];
            
                    const double idist_b1a2 = 1.0/cl::sycl::sqrt(dx_b1a2*dx_b1a2 + dy_b1a2*dy_b1a2 + dz_b1a2*dz_b1a2);
                    residuals1a_prolonged[(VAR_DENSITY)*opDat0_down_v2_kernel_stride_OP2CONSTANT]        += idist_b1a2*residuals1b[(VAR_DENSITY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+0)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_b1a2*residuals1b[(VAR_MOMENTUM+0)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+1)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_b1a2*residuals1b[(VAR_MOMENTUM+1)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_MOMENTUM+2)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_b1a2*residuals1b[(VAR_MOMENTUM+2)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1a_prolonged[(VAR_DENSITY_ENERGY)*opDat0_down_v2_kernel_stride_OP2CONSTANT] += idist_b1a2*residuals1b[(VAR_DENSITY_ENERGY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    *residuals1a_prolonged_wsum += idist_b1a2;
                }
            
                double dx_b1b2 = coord2b[(0)*opDat0_down_v2_kernel_stride_OP2CONSTANT] - coord1b[(0)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                double dy_b1b2 = coord2b[(1)*opDat0_down_v2_kernel_stride_OP2CONSTANT] - coord1b[(1)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                double dz_b1b2 = coord2b[(2)*opDat0_down_v2_kernel_stride_OP2CONSTANT] - coord1b[(2)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                if (dx_b1b2 == 0.0 && dy_b1b2 == 0.0 && dz_b1b2 == 0.0) {
            
                    residuals1b_prolonged[(VAR_DENSITY)*opDat0_down_v2_kernel_stride_OP2CONSTANT]        = residuals1b[(VAR_DENSITY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+0)*opDat0_down_v2_kernel_stride_OP2CONSTANT]      = residuals1b[(VAR_MOMENTUM+0)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+1)*opDat0_down_v2_kernel_stride_OP2CONSTANT]      = residuals1b[(VAR_MOMENTUM+1)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+2)*opDat0_down_v2_kernel_stride_OP2CONSTANT]      = residuals1b[(VAR_MOMENTUM+2)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_DENSITY_ENERGY)*opDat0_down_v2_kernel_stride_OP2CONSTANT] = residuals1b[(VAR_DENSITY_ENERGY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    *residuals1b_prolonged_wsum = 1.0;
                } else {
            
                    const double idist_b1b2 = 1.0/cl::sycl::sqrt(dx_b1b2*dx_b1b2 + dy_b1b2*dy_b1b2 + dz_b1b2*dz_b1b2);
                    residuals1b_prolonged[(VAR_DENSITY)*opDat0_down_v2_kernel_stride_OP2CONSTANT]        += idist_b1b2*residuals1b[(VAR_DENSITY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+0)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_b1b2*residuals1b[(VAR_MOMENTUM+0)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+1)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_b1b2*residuals1b[(VAR_MOMENTUM+1)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+2)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_b1b2*residuals1b[(VAR_MOMENTUM+2)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_DENSITY_ENERGY)*opDat0_down_v2_kernel_stride_OP2CONSTANT] += idist_b1b2*residuals1b[(VAR_DENSITY_ENERGY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    *residuals1b_prolonged_wsum += idist_b1b2;
            
                    double dx_a1b2 = coord1a[(0)*opDat2_down_v2_kernel_stride_OP2CONSTANT] - coord2b[(0)*opDat0_down_v2_kernel_stride_OP2CONSTANT];
                    double dy_a1b2 = coord1a[(1)*opDat2_down_v2_kernel_stride_OP2CONSTANT] - coord2b[(1)*opDat0_down_v2_kernel_stride_OP2CONSTANT];
                    double dz_a1b2 = coord1a[(2)*opDat2_down_v2_kernel_stride_OP2CONSTANT] - coord2b[(2)*opDat0_down_v2_kernel_stride_OP2CONSTANT];
            
                    const double idist_a1b2 = 1.0/cl::sycl::sqrt(dx_a1b2*dx_a1b2 + dy_a1b2*dy_a1b2 + dz_a1b2*dz_a1b2);
                    residuals1b_prolonged[(VAR_DENSITY)*opDat0_down_v2_kernel_stride_OP2CONSTANT]        += idist_a1b2*residuals1b[(VAR_DENSITY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+0)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_a1b2*residuals1b[(VAR_MOMENTUM+0)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+1)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_a1b2*residuals1b[(VAR_MOMENTUM+1)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_MOMENTUM+2)*opDat0_down_v2_kernel_stride_OP2CONSTANT]     += idist_a1b2*residuals1b[(VAR_MOMENTUM+2)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    residuals1b_prolonged[(VAR_DENSITY_ENERGY)*opDat0_down_v2_kernel_stride_OP2CONSTANT] += idist_a1b2*residuals1b[(VAR_DENSITY_ENERGY)*opDat2_down_v2_kernel_stride_OP2CONSTANT];
                    *residuals1b_prolonged_wsum += idist_a1b2;
                }
            
            };
            
          auto kern = [=](cl::sycl::nd_item<1> item) [[intel::reqd_sub_group_size(SIMD_VEC)]] {


            //get sizes and shift pointers and direct-mapped data

            int blocksPerWG = (nblocks2-1)/item.get_group_range(0)+1;
            for ( int idx=item.get_group_linear_id()*blocksPerWG; idx<(item.get_group_linear_id()+1)*blocksPerWG && idx < nblocks2; idx++ ){
              int blockId = blkmap[idx + block_offset];

              int nelem    = nelems[blockId];
              int offset_b = offset[blockId];


              for ( int n=0; n<nelem; n++ ){
                int map0idx;
                int map1idx;
                int map2idx;
                int map3idx;
                map0idx = opDat0Map[n + offset_b + set_size * 0];
                map1idx = opDat0Map[n + offset_b + set_size * 1];
                map2idx = opDat2Map[n + offset_b + set_size * 0];
                map3idx = opDat2Map[n + offset_b + set_size * 1];


                //user-supplied kernel call
                down_v2_kernel_gpu(&ind_arg0[map0idx],
                                   &ind_arg0[map1idx],
                                   &ind_arg1[map2idx],
                                   &ind_arg1[map3idx],
                                   &ind_arg2[map2idx],
                                   &ind_arg2[map3idx],
                                   &ind_arg3[map0idx],
                                   &ind_arg3[map1idx],
                                   &ind_arg4[map0idx*1],
                                   &ind_arg4[map1idx*1]);
              }

            }
          };
          cgh.parallel_for<class down_v2_kernel_kernel>(cl::sycl::nd_range<1>(nthread*nblocks,nthread), kern);
        });
        }catch(cl::sycl::exception const &e) {
        std::cout << e.what() << std::endl;exit(-1);
        }

      }
      block_offset += Plan->ncolblk[col];
    }
    OP_kernels[20].transfer  += Plan->transfer;
    OP_kernels[20].transfer2 += Plan->transfer2;
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  op2_queue->wait();
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[20].time     += wall_t2 - wall_t1;
}
