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

#endif
