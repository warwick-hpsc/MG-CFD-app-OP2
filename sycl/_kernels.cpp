//
// auto-generated by op2.py
//

//global constants
#ifndef MAX_CONST_SIZE
#define MAX_CONST_SIZE 128
#endif
#ifndef SIMD_VEC
#define SIMD_VEC 8
#endif

//header
#include "op_lib_cpp.h"
#include "op_sycl_rt_support.h"
#include "op_sycl_reduction.h"

cl::sycl::buffer<double,1> *smoothing_coefficient_p=NULL;
cl::sycl::buffer<double,1> *ff_variable_p=NULL;
cl::sycl::buffer<double,1> *ff_flux_contribution_momentum_x_p=NULL;
cl::sycl::buffer<double,1> *ff_flux_contribution_momentum_y_p=NULL;
cl::sycl::buffer<double,1> *ff_flux_contribution_momentum_z_p=NULL;
cl::sycl::buffer<double,1> *ff_flux_contribution_density_energy_p=NULL;
cl::sycl::buffer<int,1> *mesh_name_p=NULL;


void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){
  if (!OP_hybrid_gpu) return;
  if (!strcmp(name,"smoothing_coefficient")) {
    smoothing_coefficient_p = static_cast<cl::sycl::buffer<double,1>*>(
        op_sycl_register_const((void*)smoothing_coefficient_p,
            (void*)new cl::sycl::buffer<double,1>((double*)dat,
                cl::sycl::range<1>(dim))));
  }
  else
  if (!strcmp(name,"ff_variable")) {
    ff_variable_p = static_cast<cl::sycl::buffer<double,1>*>(
        op_sycl_register_const((void*)ff_variable_p,
            (void*)new cl::sycl::buffer<double,1>((double*)dat,
                cl::sycl::range<1>(dim))));
  }
  else
  if (!strcmp(name,"ff_flux_contribution_momentum_x")) {
    ff_flux_contribution_momentum_x_p = static_cast<cl::sycl::buffer<double,1>*>(
        op_sycl_register_const((void*)ff_flux_contribution_momentum_x_p,
            (void*)new cl::sycl::buffer<double,1>((double*)dat,
                cl::sycl::range<1>(dim))));
  }
  else
  if (!strcmp(name,"ff_flux_contribution_momentum_y")) {
    ff_flux_contribution_momentum_y_p = static_cast<cl::sycl::buffer<double,1>*>(
        op_sycl_register_const((void*)ff_flux_contribution_momentum_y_p,
            (void*)new cl::sycl::buffer<double,1>((double*)dat,
                cl::sycl::range<1>(dim))));
  }
  else
  if (!strcmp(name,"ff_flux_contribution_momentum_z")) {
    ff_flux_contribution_momentum_z_p = static_cast<cl::sycl::buffer<double,1>*>(
        op_sycl_register_const((void*)ff_flux_contribution_momentum_z_p,
            (void*)new cl::sycl::buffer<double,1>((double*)dat,
                cl::sycl::range<1>(dim))));
  }
  else
  if (!strcmp(name,"ff_flux_contribution_density_energy")) {
    ff_flux_contribution_density_energy_p = static_cast<cl::sycl::buffer<double,1>*>(
        op_sycl_register_const((void*)ff_flux_contribution_density_energy_p,
            (void*)new cl::sycl::buffer<double,1>((double*)dat,
                cl::sycl::range<1>(dim))));
  }
  else
  if (!strcmp(name,"mesh_name")) {
    mesh_name_p = static_cast<cl::sycl::buffer<int,1>*>(
        op_sycl_register_const((void*)mesh_name_p,
            (void*)new cl::sycl::buffer<int,1>((int*)dat,
                cl::sycl::range<1>(dim))));
  }
  else
  {
    printf("error: unknown const name\n"); exit(1);
  }
}

//user kernel files
#include "initialize_variables_kernel_kernel.cpp"
#include "zero_5d_array_kernel_kernel.cpp"
#include "zero_1d_array_kernel_kernel.cpp"
#include "calculate_cell_volumes_kernel.cpp"
#include "dampen_ewt_kernel.cpp"
#include "copy_double_kernel_kernel.cpp"
#include "calculate_dt_kernel_kernel.cpp"
#include "get_min_dt_kernel_kernel.cpp"
#include "compute_step_factor_kernel_kernel.cpp"
#include "compute_flux_edge_kernel_kernel.cpp"
#include "compute_bnd_node_flux_kernel_kernel.cpp"
#include "time_step_kernel_kernel.cpp"
#include "unstructured_stream_kernel_kernel.cpp"
#include "residual_kernel_kernel.cpp"
#include "calc_rms_kernel_kernel.cpp"
#include "count_bad_vals_kernel.cpp"
#include "up_pre_kernel_kernel.cpp"
#include "up_kernel_kernel.cpp"
#include "up_post_kernel_kernel.cpp"
#include "down_v2_kernel_pre_kernel.cpp"
#include "down_v2_kernel_kernel.cpp"
#include "down_v2_kernel_post_kernel.cpp"
#include "down_kernel_kernel.cpp"
#include "identify_differences_kernel.cpp"
#include "count_non_zeros_kernel.cpp"
