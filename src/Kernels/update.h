// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef UPDATE_H
#define UPDATE_H

inline void update_internal_kernel(
    const float* mx,
    const float* my, 
    const float* mz, 
    const float* p, 
    const float* pe, 
    float* fluxes_a, 
    float* fluxes_b)
{
    fluxes_a[VAR_DENSITY]        += p[0];
    fluxes_a[VAR_MOMENTUM]       += mx[0];
    fluxes_a[VAR_MOMENTUM+1]     += my[0];
    fluxes_a[VAR_MOMENTUM+2]     += mz[0];
    fluxes_a[VAR_DENSITY_ENERGY] += pe[0];

    fluxes_b[VAR_DENSITY]        += p[1];
    fluxes_b[VAR_MOMENTUM]       += mx[1];
    fluxes_b[VAR_MOMENTUM+1]     += my[1];
    fluxes_b[VAR_MOMENTUM+2]     += mz[1];
    fluxes_b[VAR_DENSITY_ENERGY] += pe[1];
}

inline void update_noninternal_kernel(
    const float* mx,
    const float* my, 
    const float* mz, 
    const float* p, 
    const float* pe, 
    float* fluxes_b)
{
    fluxes_b[VAR_DENSITY] += p[0];
    fluxes_b[VAR_MOMENTUM] += mx[0];
    fluxes_b[VAR_MOMENTUM+1] += my[0];
    fluxes_b[VAR_MOMENTUM+2] += mz[0];
    fluxes_b[VAR_DENSITY_ENERGY] += pe[0];
}

#endif