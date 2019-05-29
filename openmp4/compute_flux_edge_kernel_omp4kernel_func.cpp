//
// auto-generated by op2.py
//

#include <math.h>
#include "inlined_funcs.h"
#include "global.h"
#include "config.h"

void compute_flux_edge_kernel_omp4_kernel(
  int *map0,
  int map0size,
  double *data2,
  int dat2size,
  double *data0,
  int dat0size,
  double *data3,
  int dat3size,
  int *col_reord,
  int set_size1,
  int start,
  int end,
  int num_teams,
  int nthread){

  #pragma omp target teams num_teams(num_teams) thread_limit(nthread) map(to:data2[0:dat2size]) \
    map(to: smoothing_coefficient_ompkernel)\
    map(to:col_reord[0:set_size1],map0[0:map0size],data0[0:dat0size],data3[0:dat3size])
  #pragma omp distribute parallel for schedule(static,1)
  for ( int e=start; e<end; e++ ){
    int n_op = col_reord[e];
    int map0idx = map0[n_op + set_size1 * 0];
    int map1idx = map0[n_op + set_size1 * 1];

    //variable mapping
    const double *variables_a = &data0[5 * map0idx];
    const double *variables_b = &data0[5 * map1idx];
    const double *edge_weight = &data2[3*n_op];
    double *fluxes_a = &data3[5 * map0idx];
    double *fluxes_b = &data3[5 * map1idx];

    //inline function
    
    double ewt = std::sqrt(edge_weight[0]*edge_weight[0] +
                           edge_weight[1]*edge_weight[1] +
                           edge_weight[2]*edge_weight[2]);

    double p_b = variables_b[VAR_DENSITY];

    #ifdef IDIVIDE
    double ip_b = 1.0 / p_b;
    #endif

    double pe_b, pressure_b;
    double3 velocity_b, momentum_b;
    double flux_contribution_i_momentum_x_b[NDIM],
           flux_contribution_i_momentum_y_b[NDIM],
           flux_contribution_i_momentum_z_b[NDIM],
           flux_contribution_i_density_energy_b[NDIM];

    momentum_b.x = variables_b[VAR_MOMENTUM+0];
    momentum_b.y = variables_b[VAR_MOMENTUM+1];
    momentum_b.z = variables_b[VAR_MOMENTUM+2];
    pe_b = variables_b[VAR_DENSITY_ENERGY];

    #ifdef IDIVIDE
    compute_velocity(ip_b, momentum_b, velocity_b);
    #else
    compute_velocity(p_b, momentum_b, velocity_b);
    #endif

    double speed_sqd_b = compute_speed_sqd(velocity_b);
    double speed_b = std::sqrt(speed_sqd_b);

    pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);

    #ifdef IDIVIDE
    double speed_of_sound_b = compute_speed_of_sound(ip_b, pressure_b);
    #else
    double speed_of_sound_b = compute_speed_of_sound(p_b, pressure_b);
    #endif

    compute_flux_contribution(p_b, momentum_b, pe_b,
        pressure_b, velocity_b,
        flux_contribution_i_momentum_x_b,
        flux_contribution_i_momentum_y_b,
        flux_contribution_i_momentum_z_b,
        flux_contribution_i_density_energy_b);

    double factor_a, factor_b;

    double p_a, pe_a, pressure_a;
    double3 velocity_a, momentum_a;
    double flux_contribution_i_momentum_x_a[NDIM],
           flux_contribution_i_momentum_y_a[NDIM],
           flux_contribution_i_momentum_z_a[NDIM],
           flux_contribution_i_density_energy_a[NDIM];

    p_a = variables_a[VAR_DENSITY];

    #ifdef IDIVIDE
    double ip_a = 1.0 / p_a;
    #endif

    momentum_a.x = variables_a[VAR_MOMENTUM+0];
    momentum_a.y = variables_a[VAR_MOMENTUM+1];
    momentum_a.z = variables_a[VAR_MOMENTUM+2];
    pe_a = variables_a[VAR_DENSITY_ENERGY];

    #ifdef IDIVIDE
    compute_velocity(ip_a, momentum_a, velocity_a);
    #else
    compute_velocity(p_a, momentum_a, velocity_a);
    #endif

    double speed_sqd_a = compute_speed_sqd(velocity_a);
    double speed_a = std::sqrt(speed_sqd_a);
    pressure_a = compute_pressure(p_a, pe_a, speed_sqd_a);

    #ifdef IDIVIDE
    double speed_of_sound_a = compute_speed_of_sound(ip_a, pressure_a);
    #else
    double speed_of_sound_a = compute_speed_of_sound(p_a, pressure_a);
    #endif

    compute_flux_contribution(p_a, momentum_a, pe_a,
                              pressure_a, velocity_a,
                              flux_contribution_i_momentum_x_a,
                              flux_contribution_i_momentum_y_a,
                              flux_contribution_i_momentum_z_a,
                              flux_contribution_i_density_energy_a);

    factor_a = -ewt*smoothing_coefficient_ompkernel*0.5
               *(speed_a + std::sqrt(speed_sqd_b)
               + speed_of_sound_a + speed_of_sound_b);

    factor_b = -ewt*smoothing_coefficient_ompkernel*0.5
               *(speed_b + std::sqrt(speed_sqd_a)
               + speed_of_sound_b + speed_of_sound_a);

    double factor_x = -0.5*edge_weight[0], factor_y = -0.5*edge_weight[1], factor_z = -0.5*edge_weight[2];

    fluxes_a[VAR_DENSITY] +=
        factor_a*(p_a - p_b)
      + factor_x*(momentum_a.x + momentum_b.x)
      + factor_y*(momentum_a.y + momentum_b.y)
      + factor_z*(momentum_a.z + momentum_b.z);

    fluxes_a[VAR_DENSITY_ENERGY] +=
        factor_a*(pe_a - pe_b)
      + factor_x*(flux_contribution_i_density_energy_a[0] + flux_contribution_i_density_energy_b[0])
      + factor_y*(flux_contribution_i_density_energy_a[1] + flux_contribution_i_density_energy_b[1])
      + factor_z*(flux_contribution_i_density_energy_a[2] + flux_contribution_i_density_energy_b[2]);

    fluxes_a[VAR_MOMENTUM + 0] +=
        factor_a*(momentum_a.x - momentum_b.x)
      + factor_x*(flux_contribution_i_momentum_x_a[0] + flux_contribution_i_momentum_x_b[0])
      + factor_y*(flux_contribution_i_momentum_x_a[1] + flux_contribution_i_momentum_x_b[1])
      + factor_z*(flux_contribution_i_momentum_x_a[2] + flux_contribution_i_momentum_x_b[2]);

    fluxes_a[VAR_MOMENTUM + 1] +=
        factor_a*(momentum_a.y - momentum_b.y)
      + factor_x*(flux_contribution_i_momentum_y_a[0] + flux_contribution_i_momentum_y_b[0])
      + factor_y*(flux_contribution_i_momentum_y_a[1] + flux_contribution_i_momentum_y_b[1])
      + factor_z*(flux_contribution_i_momentum_y_a[2] + flux_contribution_i_momentum_y_b[2]);

    fluxes_a[VAR_MOMENTUM + 2] +=
        factor_a*(momentum_a.z - momentum_b.z)
      + factor_x*(flux_contribution_i_momentum_z_a[0] + flux_contribution_i_momentum_z_b[0])
      + factor_y*(flux_contribution_i_momentum_z_a[1] + flux_contribution_i_momentum_z_b[1])
      + factor_z*(flux_contribution_i_momentum_z_a[2] + flux_contribution_i_momentum_z_b[2]);

    fluxes_b[VAR_DENSITY] +=
        factor_b*(p_b - p_a)
      - factor_x*(momentum_a.x + momentum_b.x)
      - factor_y*(momentum_a.y + momentum_b.y)
      - factor_z*(momentum_a.z + momentum_b.z);

    fluxes_b[VAR_DENSITY_ENERGY] +=
        factor_b*(pe_b - pe_a)
      - factor_x*(flux_contribution_i_density_energy_a[0] + flux_contribution_i_density_energy_b[0])
      - factor_y*(flux_contribution_i_density_energy_a[1] + flux_contribution_i_density_energy_b[1])
      - factor_z*(flux_contribution_i_density_energy_a[2] + flux_contribution_i_density_energy_b[2]);

    fluxes_b[VAR_MOMENTUM + 0] +=
        factor_b*(momentum_b.x - momentum_a.x)
      - factor_x*(flux_contribution_i_momentum_x_a[0] + flux_contribution_i_momentum_x_b[0])
      - factor_y*(flux_contribution_i_momentum_x_a[1] + flux_contribution_i_momentum_x_b[1])
      - factor_z*(flux_contribution_i_momentum_x_a[2] + flux_contribution_i_momentum_x_b[2]);

    fluxes_b[VAR_MOMENTUM + 1] +=
        factor_b*(momentum_b.y - momentum_a.y)
      - factor_x*(flux_contribution_i_momentum_y_a[0] + flux_contribution_i_momentum_y_b[0])
      - factor_y*(flux_contribution_i_momentum_y_a[1] + flux_contribution_i_momentum_y_b[1])
      - factor_z*(flux_contribution_i_momentum_y_a[2] + flux_contribution_i_momentum_y_b[2]);

    fluxes_b[VAR_MOMENTUM + 2] +=
        factor_b*(momentum_b.z - momentum_a.z)
      - factor_x*(flux_contribution_i_momentum_z_a[0] + flux_contribution_i_momentum_z_b[0])
      - factor_y*(flux_contribution_i_momentum_z_a[1] + flux_contribution_i_momentum_z_b[1])
      - factor_z*(flux_contribution_i_momentum_z_a[2] + flux_contribution_i_momentum_z_b[2]);
    //end inline func
  }

}
