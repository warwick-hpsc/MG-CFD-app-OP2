/*
/ File contains copy of just indirect_rw_kernel_vec() from 'vec' folder.
*/

//user function -- modified for vectorisation
#if defined __clang__ || defined __GNUC__
__attribute__((always_inline))
#endif
inline void indirect_rw_kernel_vec( const double variables_a[][SIMD_BLOCK_SIZE], const double variables_b[][SIMD_BLOCK_SIZE], const double *edge_weight, double fluxes_a[][SIMD_BLOCK_SIZE], double fluxes_b[][SIMD_BLOCK_SIZE], int idx ) {
    double ex = edge_weight[0];
    double ey = edge_weight[1];
    double ez = edge_weight[2];

    double p_a, pe_a;
    double3 momentum_a;
    p_a          = variables_a[VAR_DENSITY][idx];
    momentum_a.x = variables_a[VAR_MOMENTUM+0][idx];
    momentum_a.y = variables_a[VAR_MOMENTUM+1][idx];
    momentum_a.z = variables_a[VAR_MOMENTUM+2][idx];
    pe_a         = variables_a[VAR_DENSITY_ENERGY][idx];

    double p_b, pe_b;
    double3 momentum_b;
    p_b          = variables_b[VAR_DENSITY][idx];
    momentum_b.x = variables_b[VAR_MOMENTUM+0][idx];
    momentum_b.y = variables_b[VAR_MOMENTUM+1][idx];
    momentum_b.z = variables_b[VAR_MOMENTUM+2][idx];
    pe_b         = variables_b[VAR_DENSITY_ENERGY][idx];

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

    fluxes_a[VAR_DENSITY][idx]  += p_a_val;
    fluxes_a[VAR_MOMENTUM+0][idx] += mx_a_val;
    fluxes_a[VAR_MOMENTUM+1][idx] += my_a_val;
    fluxes_a[VAR_MOMENTUM+2][idx] += mz_a_val;
    fluxes_a[VAR_DENSITY_ENERGY][idx] += pe_a_val;

    fluxes_b[VAR_DENSITY][idx]  += p_b_val;
    fluxes_b[VAR_MOMENTUM+0][idx] += mx_b_val;
    fluxes_b[VAR_MOMENTUM+1][idx] += my_b_val;
    fluxes_b[VAR_MOMENTUM+2][idx] += mz_b_val;
    fluxes_b[VAR_DENSITY_ENERGY][idx] += pe_b_val;

}
