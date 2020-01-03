// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef FLUX_H
#define FLUX_H

#include <math.h>

#include "inlined_funcs.h"

#include "global.h"
#include "config.h"

inline void compute_bnd_node_flux_kernel(
  const int *g, 
  const float *edge_weight, 
  const float *variables_b, 
  float *fluxes_b)
{
  // if (conf.legacy_mode) {
  //   if (mesh_name == MESH_LA_CASCADE && ((*g)==0 || (*g)==1 || (*g)==2)) {
  //     #include "flux_boundary.elem_func"
  //   } else if (mesh_name == MESH_LA_CASCADE && ((*g)==3 || (*g)==4 || (*g)==5 || (*g)==6)) {
  //     #include "flux_wall.elem_func"
  //   }
  // }

  // else {
    if ((*g) <= 2) {
      // Physical surface.
      #include "flux_boundary.elem_func"
    // } else {
    } else if ((*g) == 3 || ((*g) >= 4 && (*g) <= 7) ) {
      // g==3 => Freestream. Treat as far field.
      // g in range [4,7] => in/out sub/supersonic flow. Also treat as far field.
      #include "flux_wall.elem_func"
    }
  // }
}

inline void compute_flux_edge_kernel(
    const float *variables_a,
    const float *variables_b,
    const float *edge_weight,
    float *fluxes_a, 
    float *fluxes_b)
{
  double ewt = std::sqrt(edge_weight[0]*edge_weight[0] +
                         edge_weight[1]*edge_weight[1] +
                         edge_weight[2]*edge_weight[2]);

  float p_b = variables_b[VAR_DENSITY];

  #ifdef IDIVIDE
  float ip_b = 1.0f / p_b;
  #endif

  float pe_b, pressure_b;
  float3 velocity_b, momentum_b;
  float flux_contribution_i_momentum_x_b[NDIM],
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

  float speed_sqd_b = compute_speed_sqd(velocity_b);
  float speed_b = std::sqrt(speed_sqd_b);

  pressure_b = compute_pressure(p_b, pe_b, speed_sqd_b);

  #ifdef IDIVIDE
  float speed_of_sound_b = compute_speed_of_sound(ip_b, pressure_b);
  #else
  float speed_of_sound_b = compute_speed_of_sound(p_b, pressure_b);
  #endif

  compute_flux_contribution(p_b, momentum_b, pe_b,
      pressure_b, velocity_b,
      flux_contribution_i_momentum_x_b,
      flux_contribution_i_momentum_y_b,
      flux_contribution_i_momentum_z_b,
      flux_contribution_i_density_energy_b);

  float factor_a, factor_b;

  //a
  float p_a, pe_a, pressure_a;
  float3 velocity_a, momentum_a;
  float flux_contribution_i_momentum_x_a[NDIM],
         flux_contribution_i_momentum_y_a[NDIM],
         flux_contribution_i_momentum_z_a[NDIM],
         flux_contribution_i_density_energy_a[NDIM];

  p_a = variables_a[VAR_DENSITY];

  #ifdef IDIVIDE
  float ip_a = 1.0f / p_a;
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

  float speed_sqd_a = compute_speed_sqd(velocity_a);
  float speed_a = std::sqrt(speed_sqd_a);
  pressure_a = compute_pressure(p_a, pe_a, speed_sqd_a);

  #ifdef IDIVIDE
  float speed_of_sound_a = compute_speed_of_sound(ip_a, pressure_a);
  #else
  float speed_of_sound_a = compute_speed_of_sound(p_a, pressure_a);
  #endif

  compute_flux_contribution(p_a, momentum_a, pe_a,
                            pressure_a, velocity_a,
                            flux_contribution_i_momentum_x_a,
                            flux_contribution_i_momentum_y_a,
                            flux_contribution_i_momentum_z_a,
                            flux_contribution_i_density_energy_a);

  //b
  factor_a = -ewt*smoothing_coefficient*0.5f
             *(speed_a + std::sqrt(speed_sqd_b)
             + speed_of_sound_a + speed_of_sound_b);

  factor_b = -ewt*smoothing_coefficient*0.5f
             *(speed_b + std::sqrt(speed_sqd_a)
             + speed_of_sound_b + speed_of_sound_a);

  float factor_x = -0.5f*edge_weight[0], factor_y = -0.5f*edge_weight[1], factor_z = -0.5f*edge_weight[2];

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
}

#endif
