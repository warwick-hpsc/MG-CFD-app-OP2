//
// auto-generated by op2.py
//

#include <cmath>
#include "const.h"
#include "structures.h"
#include "global.h"

//user function
class calculate_cell_volumes_kernel;

//host stub function
void op_par_loop_calculate_cell_volumes(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(3);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[3].name      = name;
  OP_kernels[3].count    += 1;


  int    ninds   = 2;
  int    inds[5] = {0,0,-1,1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: calculate_cell_volumes\n");
  }

  //get plan
  #ifdef OP_PART_SIZE_3
    int part_size = OP_PART_SIZE_3;
  #else
    int part_size = OP_part_size;
  #endif

  int exec_size = op_mpi_halo_exchanges_cuda(set, nargs, args);
  if (exec_size > 0) {

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_STAGE_ALL);

    const int opDat0_calculate_cell_volumes_stride_OP2CONSTANT = getSetSizeFromOpArg(&arg0);
    const int direct_calculate_cell_volumes_stride_OP2CONSTANT = getSetSizeFromOpArg(&arg2);
    double *ind_arg0 = (double*)arg0.data_d;
    double *ind_arg1 = (double*)arg3.data_d;
    int *opDat0Map = arg0.map_data_d;
    double *arg2_d = (double*)arg2.data_d;
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
          auto calculate_cell_volumes_gpu = [=]( 
                const double* coords1,
                const double* coords2,
                double* ewt,
                double* vol1,
                double* vol2) {
                double d[NDIM];
                double dist = 0.0;
                for (int i=0; i<NDIM; i++) {
                    d[i] = coords2[(i)*opDat0_calculate_cell_volumes_stride_OP2CONSTANT] - coords1[(i)*opDat0_calculate_cell_volumes_stride_OP2CONSTANT];
                    dist += d[i]*d[i];
                }
                dist = cl::sycl::sqrt(dist);
            
                double area = 0.0;
                for (int i=0; i<NDIM; i++) {
                    area += ewt[(i)*direct_calculate_cell_volumes_stride_OP2CONSTANT]*ewt[(i)*direct_calculate_cell_volumes_stride_OP2CONSTANT];
                }
                area = cl::sycl::sqrt(area);
            
                double tetra_volume = (1.0/3.0)*0.5 *dist *area;
                *vol1 += tetra_volume;
                *vol2 += tetra_volume;
            
                for (int i=0; i<NDIM; i++) {
                    ewt[(i)*direct_calculate_cell_volumes_stride_OP2CONSTANT] = (d[i] / dist) * area;
                }
            
            
            
                for (int i=0; i<NDIM; i++) {
                    ewt[(i)*direct_calculate_cell_volumes_stride_OP2CONSTANT] /= dist;
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
                map0idx = opDat0Map[n + offset_b + set_size * 0];
                map1idx = opDat0Map[n + offset_b + set_size * 1];


                //user-supplied kernel call
                calculate_cell_volumes_gpu(&ind_arg0[map0idx],
                                           &ind_arg0[map1idx],
                                           &arg2_d[n+offset_b],
                                           &ind_arg1[map0idx*1],
                                           &ind_arg1[map1idx*1]);
              }

            }
          };
          cgh.parallel_for<class calculate_cell_volumes_kernel>(cl::sycl::nd_range<1>(nthread*nblocks,nthread), kern);
        });
        }catch(cl::sycl::exception const &e) {
        std::cout << e.what() << std::endl;exit(-1);
        }

      }
      block_offset += Plan->ncolblk[col];
    }
    OP_kernels[3].transfer  += Plan->transfer;
    OP_kernels[3].transfer2 += Plan->transfer2;
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  op2_queue->wait();
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[3].time     += wall_t2 - wall_t1;
}
