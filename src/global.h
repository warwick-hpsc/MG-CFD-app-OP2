#ifndef GLOBAL_H
#define GLOBAL_H

extern int mesh_name;

extern float smoothing_coefficient;
extern float ff_variable[5];
extern float ff_flux_contribution_momentum_x[3];
extern float ff_flux_contribution_momentum_y[3];
extern float ff_flux_contribution_momentum_z[3];
extern float ff_flux_contribution_density_energy[3];

extern int levels;
extern int current_level;

#endif
