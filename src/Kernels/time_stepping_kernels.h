// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef COMPUTE_STEP_FACTOR_H
#define COMPUTE_STEP_FACTOR_H

#include <math.h>
#include <cmath>

#include "const.h"
#include "inlined_funcs.h"

inline void calculate_dt_kernel(
    const double* variable, 
    const double* volume, 
    double* dt)
{
    double density = variable[VAR_DENSITY];

    double3 momentum;
    momentum.x = variable[VAR_MOMENTUM+0];
    momentum.y = variable[VAR_MOMENTUM+1];
    momentum.z = variable[VAR_MOMENTUM+2];

    double density_energy = variable[VAR_DENSITY_ENERGY];
    double3 velocity; compute_velocity(density, momentum, velocity);
    double speed_sqd      = compute_speed_sqd(velocity);
    double pressure       = compute_pressure(density, density_energy, speed_sqd);
    double speed_of_sound = compute_speed_of_sound(density, pressure);

    *dt = double(0.5) * (cbrt(*volume) / (sqrt(speed_sqd) + speed_of_sound));
}

inline void get_min_dt_kernel(
    const double* dt, 
    double* min_dt)
{
    if ((*dt) < (*min_dt)) {
        *min_dt = (*dt);
    }
}

inline void compute_step_factor_kernel(
    const double* variable, 
    const double* volume, 
    const double* min_dt, 
    double* step_factor)
{
    double density = variable[VAR_DENSITY];

    double3 momentum;
    momentum.x = variable[VAR_MOMENTUM+0];
    momentum.y = variable[VAR_MOMENTUM+1];
    momentum.z = variable[VAR_MOMENTUM+2];

    double density_energy = variable[VAR_DENSITY_ENERGY];
    double3 velocity; compute_velocity(density, momentum, velocity);
    double speed_sqd      = compute_speed_sqd(velocity);
    double pressure       = compute_pressure(density, density_energy, speed_sqd);
    double speed_of_sound = compute_speed_of_sound(density, pressure);

    // Bring forward a future division-by-volume:
    *step_factor = (*min_dt) / (*volume);
}

inline void time_step_kernel(
    const int* rkCycle,
    const double* step_factor,
    double* flux,
    const double* old_variable,
    double* variable)
{
    double factor = (*step_factor)/double(RK+1-(*rkCycle));

    variable[VAR_DENSITY]        = old_variable[VAR_DENSITY]        + factor*flux[VAR_DENSITY];
    variable[VAR_MOMENTUM+0]     = old_variable[VAR_MOMENTUM+0]     + factor*flux[VAR_MOMENTUM+0];
    variable[VAR_MOMENTUM+1]     = old_variable[VAR_MOMENTUM+1]     + factor*flux[VAR_MOMENTUM+1];
    variable[VAR_MOMENTUM+2]     = old_variable[VAR_MOMENTUM+2]     + factor*flux[VAR_MOMENTUM+2];
    variable[VAR_DENSITY_ENERGY] = old_variable[VAR_DENSITY_ENERGY] + factor*flux[VAR_DENSITY_ENERGY];

    flux[VAR_DENSITY]        = 0.0;
    flux[VAR_MOMENTUM+0]     = 0.0;
    flux[VAR_MOMENTUM+1]     = 0.0;
    flux[VAR_MOMENTUM+2]     = 0.0;
    flux[VAR_DENSITY_ENERGY] = 0.0;
}

#endif
