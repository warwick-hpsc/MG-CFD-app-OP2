// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef READ_H
#define READ_H

#include <math.h>

#include "inlined_funcs.h"

#include "global.h"
#include "config.h"
#include "unstructured_stream.h"

inline void test_read_kernel(
    const double *variables_a,
    const double *variables_b,
    const double *edge_weight,
    double *fluxes_a,
    double *fluxes_b)
{
    double ex = edge_weight[0];
    double ey = edge_weight[1];
    double ez = edge_weight[2];

    double p_a, pe_a;
    double3 momentum_a;
    p_a          = variables_a[VAR_DENSITY];
    momentum_a.x = variables_a[VAR_MOMENTUM+0];
    momentum_a.y = variables_a[VAR_MOMENTUM+1];
    momentum_a.z = variables_a[VAR_MOMENTUM+2];
    pe_a         = variables_a[VAR_DENSITY_ENERGY];

    double p_b, pe_b;
    double3 momentum_b;
    p_b          = variables_b[VAR_DENSITY];
    momentum_b.x = variables_b[VAR_MOMENTUM+0];
    momentum_b.y = variables_b[VAR_MOMENTUM+1];
    momentum_b.z = variables_b[VAR_MOMENTUM+2];
    pe_b         = variables_b[VAR_DENSITY_ENERGY];

    double p_a_val  = p_b + ex;
    double pe_a_val = pe_b + ey;
    double mx_a_val = momentum_b.x + ez;
    double my_a_val = momentum_b.y;
    double mz_a_val = momentum_b.z;

    double p_b_val = p_a;
    double pe_b_val = pe_a;
    double mx_b_val = momentum_a.x;
    double my_b_val = momentum_a.y;
    double mz_b_val = momentum_a.z;

    fluxes_a[VAR_DENSITY]  += p_a_val;
    fluxes_a[VAR_MOMENTUM+0] += mx_a_val;
    fluxes_a[VAR_MOMENTUM+1] += my_a_val;
    fluxes_a[VAR_MOMENTUM+2] += mz_a_val;
    fluxes_a[VAR_DENSITY_ENERGY] += pe_a_val;

    fluxes_b[VAR_DENSITY]  += p_b_val;
    fluxes_b[VAR_MOMENTUM+0] += mx_b_val;
    fluxes_b[VAR_MOMENTUM+1] += my_b_val;
    fluxes_b[VAR_MOMENTUM+2] += mz_b_val;
    fluxes_b[VAR_DENSITY_ENERGY] += pe_b_val;

}

#endif
