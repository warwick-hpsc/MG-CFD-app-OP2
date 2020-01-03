// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef INLINED_FUNCS_H
#define INLINED_FUNCS_H

// #include <cmath>

#include "const.h"
#include "structures.h"

#define DEBUGGABLE_ABORT fprintf(stderr, "%s:%d\n", __FILE__, __LINE__); fflush(stderr); fflush(stdout); exit(EXIT_FAILURE);

// inline void zero_array(int nelr, float* variables)
static inline OP_FUN_PREFIX void zero_array(int nelr, float* variables)
{
    #ifdef OMP
    #pragma omp parallel for schedule(static)
    #endif
    for(int i = 0; i < nelr; i++)
    {
        variables[i] = 0.0;
    }
}

// inline void zero_edges(int nelr, edge* mx, edge* my, edge* mz, edge* p, edge* pe)
static inline OP_FUN_PREFIX void zero_edges(int nelr, edge* mx, edge* my, edge* mz, edge* p, edge* pe)
{
    #ifdef OMP
    #pragma omp parallel
    {
    #endif
        #ifdef VERITAS
            pdump_start_profile();
        #endif

        #ifdef OMP
        #pragma omp for schedule(static)
        #endif

        #ifdef SIMD
        #pragma simd vectorlength(SIMD)
        #else
        #pragma novector
        #endif
        for(int i = 0; i < nelr; i++)
        {
            mx[i].a = 0.0;
            mx[i].b = 0.0;
            my[i].a = 0.0;
            my[i].b = 0.0;
            mz[i].a = 0.0;
            mz[i].b = 0.0;
            p[i].a = 0.0;
            p[i].b = 0.0;
            pe[i].a = 0.0;
            pe[i].b = 0.0;
        }

        #ifdef VERITAS
            pdump_end_profile();
        #endif

    #ifdef OMP
    }
    #endif
}

// inline void compute_flux_contribution(
static inline OP_FUN_PREFIX void compute_flux_contribution(
    float& density, float3& momentum, 
    float& density_energy, 
    float& pressure, 
    float3& velocity, 
    float* fc_momentum_x, 
    float* fc_momentum_y, 
    float* fc_momentum_z, 
    float* fc_density_energy)
{
    fc_momentum_x[0] = velocity.x*momentum.x + pressure;
    fc_momentum_x[1] = velocity.x*momentum.y;
    fc_momentum_x[2] = velocity.x*momentum.z;

    fc_momentum_y[0] = fc_momentum_x[1];
    fc_momentum_y[1] = velocity.y*momentum.y + pressure;
    fc_momentum_y[2] = velocity.y*momentum.z;

    fc_momentum_z[0] = fc_momentum_x[2];
    fc_momentum_z[1] = fc_momentum_y[2];
    fc_momentum_z[2] = velocity.z*momentum.z + pressure;

    float de_p = density_energy+pressure;
    fc_density_energy[0] = velocity.x*de_p;
    fc_density_energy[1] = velocity.y*de_p;
    fc_density_energy[2] = velocity.z*de_p;
}

static inline OP_FUN_PREFIX float compute_speed_sqd(float3& velocity)
{
    return velocity.x*velocity.x + velocity.y*velocity.y + velocity.z*velocity.z;
}

static inline OP_FUN_PREFIX float compute_pressure(float& density, float& density_energy, float& speed_sqd)
{
    return (float(GAMMA)-float(1.0))*(density_energy - float(0.5)*density*speed_sqd);
}

#ifdef IDIVIDE
    static inline OP_FUN_PREFIX void compute_velocity(float& inverse_density, float3& momentum, float3& velocity)
    {
        velocity.x = momentum.x*inverse_density;
        velocity.y = momentum.y*inverse_deinlined_funcs.hnsity;
        velocity.z = momentum.z*inverse_density;
    }
    static inline OP_FUN_PREFIX float compute_speed_of_sound(float& inverse_density, float& pressure)
    {
        return sqrt((float(GAMMA)*pressure)*inverse_density);
    }
#else
    static OP_FUN_PREFIX inline void compute_velocity(float& density, float3& momentum, float3& velocity)
    {
        velocity.x = momentum.x / density;
        velocity.y = momentum.y / density;
        velocity.z = momentum.z / density;
    }
    static inline OP_FUN_PREFIX float compute_speed_of_sound(float& density, float& pressure)
    {
        return sqrt(float(GAMMA)*pressure / density);
    }
#endif

#endif
