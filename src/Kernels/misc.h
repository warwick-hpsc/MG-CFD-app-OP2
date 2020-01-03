#ifndef MISC_H
#define MISC_H

#include <cmath>

#include "const.h"
#include "structures.h"
#include "global.h"

inline void initialize_variables_kernel(
    float* variables)
{
    for(int j = 0; j < NVAR; j++) {
        variables[j] = ff_variable[j];
    }
}

inline void zero_5d_array_kernel(
    float* array)
{
    for(int j = 0; j < NVAR; j++) {
        array[j] = 0.0f;
    }
}

inline void zero_3d_array_kernel(
    float* array)
{
    for(int j = 0; j < NDIM; j++) {
        array[j] = 0.0f;
    }
}

inline void zero_1d_array_kernel(
    float* array)
{
    *array = 0.0f;
}

inline void calculate_cell_volumes(
    const float* coords1, 
    const float* coords2, 
    float* ewt,
    float* vol1, 
    float* vol2)
{
    float d[NDIM];
    float dist = 0.0f;
    for (int i=0; i<NDIM; i++) {
        d[i] = coords2[i] - coords1[i];
        dist += d[i]*d[i];
    }
    dist = sqrt(dist);

    float area = 0.0f;
    for (int i=0; i<NDIM; i++) {
        area += ewt[i]*ewt[i];
    }
    area = sqrt(area);

    float tetra_volume = (1.0f/3.0f)*0.5f *dist *area;
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
    float* ewt)
{
    ewt[0] *= 1e-7f;
    ewt[1] *= 1e-7f;
    ewt[2] *= 1e-7f;
}

#endif
