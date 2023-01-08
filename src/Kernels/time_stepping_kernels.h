// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef COMPUTE_STEP_FACTOR_H
#define COMPUTE_STEP_FACTOR_H

#include <math.h>
#include <cmath>

#include "const.h"
#include "inlined_funcs.h"

inline void calculate_dt_kernel(
    const float* variable, 
    const float* volume, 
    float* dt)
{
    float density = variable[VAR_DENSITY];

    float3 momentum;
    momentum.x = variable[VAR_MOMENTUM+0];
    momentum.y = variable[VAR_MOMENTUM+1];
    momentum.z = variable[VAR_MOMENTUM+2];

    float density_energy = variable[VAR_DENSITY_ENERGY];
    float3 velocity; compute_velocity(density, momentum, velocity);
    float speed_sqd      = compute_speed_sqd(velocity);
    float pressure       = compute_pressure(density, density_energy, speed_sqd);
    float speed_of_sound = compute_speed_of_sound(density, pressure);

    *dt = float(0.5f) * (cbrt(*volume) / (sqrt(speed_sqd) + speed_of_sound));
}

inline void get_min_dt_kernel(
    const float* dt, 
    float* min_dt)
{
    if ((*dt) < (*min_dt)) {
        *min_dt = (*dt);
    }
}

inline void compute_step_factor_kernel(
    const float* variable, 
    const float* volume, 
    const float* min_dt, 
    float* step_factor)
{
    float density = variable[VAR_DENSITY];

    float3 momentum;
    momentum.x = variable[VAR_MOMENTUM+0];
    momentum.y = variable[VAR_MOMENTUM+1];
    momentum.z = variable[VAR_MOMENTUM+2];

    float density_energy = variable[VAR_DENSITY_ENERGY];
    float3 velocity; compute_velocity(density, momentum, velocity);
    float speed_sqd      = compute_speed_sqd(velocity);
    float pressure       = compute_pressure(density, density_energy, speed_sqd);
    float speed_of_sound = compute_speed_of_sound(density, pressure);

    // Bring forward a future division-by-volume:
    *step_factor = (*min_dt) / (*volume);
}

inline void time_step_kernel(
    const int* rkCycle,
    const float* step_factor,
    float* flux,
    const float* old_variable,
    float* variable)
{
    float factor = (*step_factor)/float(RK+1-(*rkCycle));

    variable[VAR_DENSITY]        = old_variable[VAR_DENSITY]        + factor*flux[VAR_DENSITY];
    variable[VAR_MOMENTUM+0]     = old_variable[VAR_MOMENTUM+0]     + factor*flux[VAR_MOMENTUM+0];
    variable[VAR_MOMENTUM+1]     = old_variable[VAR_MOMENTUM+1]     + factor*flux[VAR_MOMENTUM+1];
    variable[VAR_MOMENTUM+2]     = old_variable[VAR_MOMENTUM+2]     + factor*flux[VAR_MOMENTUM+2];
    variable[VAR_DENSITY_ENERGY] = old_variable[VAR_DENSITY_ENERGY] + factor*flux[VAR_DENSITY_ENERGY];

    flux[VAR_DENSITY]        = 0.0f;
    flux[VAR_MOMENTUM+0]     = 0.0f;
    flux[VAR_MOMENTUM+1]     = 0.0f;
    flux[VAR_MOMENTUM+2]     = 0.0f;
    flux[VAR_DENSITY_ENERGY] = 0.0f;
}

#endif
