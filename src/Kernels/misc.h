#ifndef MISC_H
#define MISC_H

#include <cmath>

#include "const.h"
#include "structures.h"
#include "global.h"

inline void initialize_variables_kernel(
    double* variables)
{
    for(int j = 0; j < NVAR; j++) {
        variables[j] = ff_variable[j];
    }
}

inline void zero_5d_array_kernel(
    double* array)
{
    for(int j = 0; j < NVAR; j++) {
        array[j] = 0.0;
    }
}

inline void zero_3d_array_kernel(
    double* array)
{
    for(int j = 0; j < NDIM; j++) {
        array[j] = 0.0;
    }
}

inline void zero_1d_array_kernel(
    double* array)
{
    *array = 0.0;
}

inline void calculate_cell_volumes(
    const double* coords1, 
    const double* coords2, 
    double* ewt,
    double* vol1, 
    double* vol2)
{
    double d[NDIM];
    double dist = 0.0;
    for (int i=0; i<NDIM; i++) {
        d[i] = coords2[i] - coords1[i];
        dist += d[i]*d[i];
    }
    dist = sqrt(dist);

    double area = 0.0;
    for (int i=0; i<NDIM; i++) {
        area += ewt[i]*ewt[i];
    }
    area = sqrt(area);

    double tetra_volume = (1.0/3.0)*0.5 *dist *area;
    *vol1 += tetra_volume;
    *vol2 += tetra_volume;

    // Redirect ewt to be parallel to normal:
    for (int i=0; i<NDIM; i++) {
        ewt[i] = (d[i] / dist) * area;
    }

    // |ewt| currently is face area. Divide through by distance 
    // to produce 'surface vector' with magnitude (area/dm), 
    // for use in flux accumulation:
    for (int i=0; i<NDIM; i++) {
        ewt[i] /= dist;
    }
}

inline void dampen_ewt(
    double* ewt)
{
    // ewt[0] *= 1e-7;
    // ewt[1] *= 1e-7;
    // ewt[2] *= 1e-7;

    ewt[0] *= 5e-8;
    ewt[1] *= 5e-8;
    ewt[2] *= 5e-8;
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
