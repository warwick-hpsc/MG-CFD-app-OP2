//
// auto-generated by op2.py
//

//user function
// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

// #ifndef FLUX_H
// #define FLUX_H

#include <math.h>

#include "inlined_funcs.h"

#include "global.h"
#include "config.h"

// #ifdef VECTORIZE
//user function -- modified for vectorisation
#if defined __clang__ || defined __GNUC__
__attribute__((always_inline))
#endif
inline void compute_bnd_node_flux_kernel_vec( const int g[][SIMD_VEC], const double edge_weight[][SIMD_VEC], const double variables_b[][SIMD_VEC], double fluxes_b[][SIMD_VEC], int idx ) {








    if ((g[0][idx]) <= 2) {

      
      
      
      
      
      
      
      
          double p_b = variables_b[VAR_DENSITY][idx];
      
          #ifdef IDIVIDE
          double ip_b = 1.0 / p_b;
          #endif
      
          double pe_b, pressure_b;
          double3 velocity_b, momentum_b;
          double flux_contribution_i_momentum_x_b[NDIM],
                 flux_contribution_i_momentum_y_b[NDIM],
                 flux_contribution_i_momentum_z_b[NDIM],
                 flux_contribution_i_density_energy_b[NDIM];
      
          momentum_b.x = variables_b[VAR_MOMENTUM+0][idx];
          momentum_b.y = variables_b[VAR_MOMENTUM+1][idx];
          momentum_b.z = variables_b[VAR_MOMENTUM+2][idx];
          pe_b = variables_b[VAR_DENSITY_ENERGY][idx];
      
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
      
          fluxes_b[VAR_DENSITY][idx]        = 0;
          fluxes_b[VAR_MOMENTUM +0][idx]    = edge_weight[0][idx]*pressure_b;
          fluxes_b[VAR_MOMENTUM +1][idx]    = edge_weight[1][idx]*pressure_b;
          fluxes_b[VAR_MOMENTUM +2][idx]    = edge_weight[2][idx]*pressure_b;
          fluxes_b[VAR_DENSITY_ENERGY][idx] = 0;
      

    } else if ((g[0][idx]) == 3 || ((g[0][idx]) >= 4 && (g[0][idx]) <= 7) ) {


      
      
      
      
      
      
      
          double p_b = variables_b[VAR_DENSITY][idx];
      
          #ifdef IDIVIDE
          double ip_b = 1.0 / p_b;
          #endif
      
          double pe_b, pressure_b;
          double3 velocity_b, momentum_b;
          double flux_contribution_i_momentum_x_b[NDIM],
                 flux_contribution_i_momentum_y_b[NDIM],
                 flux_contribution_i_momentum_z_b[NDIM],
                 flux_contribution_i_density_energy_b[NDIM];
      
          momentum_b.x = variables_b[VAR_MOMENTUM+0][idx];
          momentum_b.y = variables_b[VAR_MOMENTUM+1][idx];
          momentum_b.z = variables_b[VAR_MOMENTUM+2][idx];
          pe_b = variables_b[VAR_DENSITY_ENERGY][idx];
      
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
      
          double factor_x = 0.5 * edge_weight[0][idx],
                 factor_y = 0.5 * edge_weight[1][idx],
                 factor_z = 0.5 * edge_weight[2][idx];
      
          fluxes_b[VAR_DENSITY][idx] =
                factor_x*(ff_variable[VAR_MOMENTUM+0] + momentum_b.x)
              + factor_y*(ff_variable[VAR_MOMENTUM+1] + momentum_b.y)
              + factor_z*(ff_variable[VAR_MOMENTUM+2] + momentum_b.z);
      
          fluxes_b[VAR_DENSITY_ENERGY][idx] = 
                factor_x*(ff_flux_contribution_density_energy[0] + flux_contribution_i_density_energy_b[0])
              + factor_y*(ff_flux_contribution_density_energy[1] + flux_contribution_i_density_energy_b[1])
              + factor_z*(ff_flux_contribution_density_energy[2] + flux_contribution_i_density_energy_b[2]);
      
          fluxes_b[VAR_MOMENTUM + 0][idx] = 
                factor_x*(ff_flux_contribution_momentum_x[0] + flux_contribution_i_momentum_x_b[0])
              + factor_y*(ff_flux_contribution_momentum_x[1] + flux_contribution_i_momentum_x_b[1])
              + factor_z*(ff_flux_contribution_momentum_x[2] + flux_contribution_i_momentum_x_b[2]);
      
          fluxes_b[VAR_MOMENTUM + 1][idx] = 
                factor_x*(ff_flux_contribution_momentum_y[0] + flux_contribution_i_momentum_y_b[0])
              + factor_y*(ff_flux_contribution_momentum_y[1] + flux_contribution_i_momentum_y_b[1])
              + factor_z*(ff_flux_contribution_momentum_y[2] + flux_contribution_i_momentum_y_b[2]);
      
          fluxes_b[VAR_MOMENTUM + 2][idx] = 
                factor_x*(ff_flux_contribution_momentum_z[0] + flux_contribution_i_momentum_z_b[0])
              + factor_y*(ff_flux_contribution_momentum_z[1] + flux_contribution_i_momentum_z_b[1])
              + factor_z*(ff_flux_contribution_momentum_z[2] + flux_contribution_i_momentum_z_b[2]);
      
    }


}
// #endif