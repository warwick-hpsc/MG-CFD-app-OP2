#ifndef UNSTRUCTURED_STREAM_H
#define UNSTRUCTURED_STREAM_H

// Unstructured stream kernel
// - performs same data movement as compute_flux_edge() but with minimal arithmetic. 
//   Measures upper bound on performance achievable by compute_flux_edge()
inline void unstructured_stream_kernel(
    const float *variables_a,
    const float *variables_b,
    const float *edge_weight,
    float *fluxes_a, 
    float *fluxes_b)
{
    float ex = edge_weight[0];
    float ey = edge_weight[1];
    float ez = edge_weight[2];

    float p_a, pe_a;
    float3 momentum_a;
    p_a          = variables_a[VAR_DENSITY];
    momentum_a.x = variables_a[VAR_MOMENTUM+0];
    momentum_a.y = variables_a[VAR_MOMENTUM+1];
    momentum_a.z = variables_a[VAR_MOMENTUM+2];
    pe_a         = variables_a[VAR_DENSITY_ENERGY];

    float p_b, pe_b;
    float3 momentum_b;
    p_b          = variables_b[VAR_DENSITY];
    momentum_b.x = variables_b[VAR_MOMENTUM+0];
    momentum_b.y = variables_b[VAR_MOMENTUM+1];
    momentum_b.z = variables_b[VAR_MOMENTUM+2];
    pe_b         = variables_b[VAR_DENSITY_ENERGY];

    float p_a_val  = p_b + ex;
    float pe_a_val = pe_b + ey;
    float mx_a_val = momentum_b.x + ez;
    float my_a_val = momentum_b.y;
    float mz_a_val = momentum_b.z;

    float p_b_val = p_a;
    float pe_b_val = pe_a;
    float mx_b_val = momentum_a.x;
    float my_b_val = momentum_a.y;
    float mz_b_val = momentum_a.z;

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
