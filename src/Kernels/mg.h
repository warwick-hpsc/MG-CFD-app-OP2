//************************************************//
// Copyright 2016-2019 University of Warwick

// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
// sell copies of the Software, and to permit persons to whom the Software is furnished 
// to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//************************************************//

#ifndef MG_H
#define MG_H

#include <math.h>

#include "const.h"

inline void up_pre_kernel(
    float* variable, 
    int* up_scratch)
{
    variable[VAR_DENSITY] = 0.0f;
    variable[VAR_MOMENTUM+0] = 0.0f;
    variable[VAR_MOMENTUM+1] = 0.0f;
    variable[VAR_MOMENTUM+2] = 0.0f;
    variable[VAR_DENSITY_ENERGY] = 0.0f;
    *up_scratch = 0;
}

inline void up_kernel(
    const float* variable, 
    float* variable_above, 
    int* up_scratch)
{
    variable_above[VAR_DENSITY]        += variable[VAR_DENSITY];
    variable_above[VAR_MOMENTUM+0]     += variable[VAR_MOMENTUM+0];
    variable_above[VAR_MOMENTUM+1]     += variable[VAR_MOMENTUM+1];
    variable_above[VAR_MOMENTUM+2]     += variable[VAR_MOMENTUM+2];
    variable_above[VAR_DENSITY_ENERGY] += variable[VAR_DENSITY_ENERGY];
    *up_scratch += 1;
}

inline void up_post_kernel(
    float* variable, 
    const int* up_scratch)
{
    float avg = (*up_scratch)==0 ? 1.0f : 1.0f / (float)(*up_scratch);
    variable[VAR_DENSITY] *= avg;
    variable[VAR_MOMENTUM+0] *= avg;
    variable[VAR_MOMENTUM+1] *= avg;
    variable[VAR_MOMENTUM+2] *= avg;
    variable[VAR_DENSITY_ENERGY] *= avg;
}

inline void down_kernel(
    float* variable, 
    const float* residual, 
    const float* coord, 
    const float* residual_above, 
    const float* coord_above)
{
  float dx = fabs(coord[0] - coord_above[0]);
  float dy = fabs(coord[1] - coord_above[1]);
  float dz = fabs(coord[2] - coord_above[2]);
  float dm = sqrt(dx * dx + dy * dy + dz * dz);

  variable[VAR_DENSITY] -=
      dm * (residual_above[VAR_DENSITY] - residual[VAR_DENSITY]);
  variable[VAR_MOMENTUM + 0] -=
      dx * (residual_above[VAR_MOMENTUM + 0] - residual[VAR_MOMENTUM + 0]);
  variable[VAR_MOMENTUM + 1] -=
      dy * (residual_above[VAR_MOMENTUM + 1] - residual[VAR_MOMENTUM + 1]);
  variable[VAR_MOMENTUM + 2] -=
      dz * (residual_above[VAR_MOMENTUM + 2] - residual[VAR_MOMENTUM + 2]);
  variable[VAR_DENSITY_ENERGY] -=
      dm * (residual_above[VAR_DENSITY_ENERGY] - residual[VAR_DENSITY_ENERGY]);
}

inline void down_v2_kernel_pre(
    float* weight_sum, 
    float* residual_sum)
{
    *weight_sum = 0.0f;
    residual_sum[VAR_DENSITY] = 0.0f;
    residual_sum[VAR_MOMENTUM+0] = 0.0f;
    residual_sum[VAR_MOMENTUM+2] = 0.0f;
    residual_sum[VAR_MOMENTUM+1] = 0.0f;
    residual_sum[VAR_DENSITY_ENERGY] = 0.0f;
}
inline void down_v2_kernel(
    const float* coord2a, 
    const float* coord2b, 
    const float* coord1a, 
    const float* coord1b, 
    const float* residuals1a,
    const float* residuals1b,
    float* residuals1a_prolonged, 
    float* residuals1b_prolonged, 
    float* residuals1a_prolonged_wsum,
    float* residuals1b_prolonged_wsum)
{
    // For each node that has the same coordinates as its MG node parent, 
    // the 'prolonged residual' is simply taken directly from the MG node. 
    //
    // For each other node N, the 'prolonged residuals' is the weighted average 
    // across N's MG node and MG nodes of N's neighbours, requiring an 
    // edge-based loop. The weight is 1.0f/distance.

    // Process a2:
    float dx_a1a2 = coord2a[0] - coord1a[0];
    float dy_a1a2 = coord2a[1] - coord1a[1];
    float dz_a1a2 = coord2a[2] - coord1a[2];
    if (dx_a1a2 == 0.0f && dy_a1a2 == 0.0f && dz_a1a2 == 0.0f) {
        // a2 == a1:
        residuals1a_prolonged[VAR_DENSITY]        = residuals1a[VAR_DENSITY];
        residuals1a_prolonged[VAR_MOMENTUM+0]     = residuals1a[VAR_MOMENTUM+0];
        residuals1a_prolonged[VAR_MOMENTUM+1]     = residuals1a[VAR_MOMENTUM+1];
        residuals1a_prolonged[VAR_MOMENTUM+2]     = residuals1a[VAR_MOMENTUM+2];
        residuals1a_prolonged[VAR_DENSITY_ENERGY] = residuals1a[VAR_DENSITY_ENERGY];
        *residuals1a_prolonged_wsum = 1.0f;
    } else {
        // Calculate contribution of a1 -> a2:
        const float idist_a1a2 = 1.0f/sqrt(dx_a1a2*dx_a1a2 + dy_a1a2*dy_a1a2 + dz_a1a2*dz_a1a2);
        residuals1a_prolonged[VAR_DENSITY]        += idist_a1a2*residuals1a[VAR_DENSITY];
        residuals1a_prolonged[VAR_MOMENTUM+0]     += idist_a1a2*residuals1a[VAR_MOMENTUM+0];
        residuals1a_prolonged[VAR_MOMENTUM+1]     += idist_a1a2*residuals1a[VAR_MOMENTUM+1];
        residuals1a_prolonged[VAR_MOMENTUM+2]     += idist_a1a2*residuals1a[VAR_MOMENTUM+2];
        residuals1a_prolonged[VAR_DENSITY_ENERGY] += idist_a1a2*residuals1a[VAR_DENSITY_ENERGY];
        *residuals1a_prolonged_wsum += idist_a1a2;

        // Calculate contribution of b1 >- a2:
        float dx_b1a2 = coord1b[0] - coord2a[0];
        float dy_b1a2 = coord1b[1] - coord2a[1];
        float dz_b1a2 = coord1b[2] - coord2a[2];

        const float idist_b1a2 = 1.0f/sqrt(dx_b1a2*dx_b1a2 + dy_b1a2*dy_b1a2 + dz_b1a2*dz_b1a2);
        residuals1a_prolonged[VAR_DENSITY]        += idist_b1a2*residuals1b[VAR_DENSITY];
        residuals1a_prolonged[VAR_MOMENTUM+0]     += idist_b1a2*residuals1b[VAR_MOMENTUM+0];
        residuals1a_prolonged[VAR_MOMENTUM+1]     += idist_b1a2*residuals1b[VAR_MOMENTUM+1];
        residuals1a_prolonged[VAR_MOMENTUM+2]     += idist_b1a2*residuals1b[VAR_MOMENTUM+2];
        residuals1a_prolonged[VAR_DENSITY_ENERGY] += idist_b1a2*residuals1b[VAR_DENSITY_ENERGY];
        *residuals1a_prolonged_wsum += idist_b1a2;
    }

    // Process b2:
    float dx_b1b2 = coord2b[0] - coord1b[0];
    float dy_b1b2 = coord2b[1] - coord1b[1];
    float dz_b1b2 = coord2b[2] - coord1b[2];
    if (dx_b1b2 == 0.0f && dy_b1b2 == 0.0f && dz_b1b2 == 0.0f) {
        // b2 == b1:
        residuals1b_prolonged[VAR_DENSITY]        = residuals1b[VAR_DENSITY];
        residuals1b_prolonged[VAR_MOMENTUM+0]      = residuals1b[VAR_MOMENTUM+0];
        residuals1b_prolonged[VAR_MOMENTUM+1]      = residuals1b[VAR_MOMENTUM+1];
        residuals1b_prolonged[VAR_MOMENTUM+2]      = residuals1b[VAR_MOMENTUM+2];
        residuals1b_prolonged[VAR_DENSITY_ENERGY] = residuals1b[VAR_DENSITY_ENERGY];
        *residuals1b_prolonged_wsum = 1.0f;
    } else {
        // Calculate contribution of b1 -> b2:
        const float idist_b1b2 = 1.0f/sqrt(dx_b1b2*dx_b1b2 + dy_b1b2*dy_b1b2 + dz_b1b2*dz_b1b2);
        residuals1b_prolonged[VAR_DENSITY]        += idist_b1b2*residuals1b[VAR_DENSITY];
        residuals1b_prolonged[VAR_MOMENTUM+0]     += idist_b1b2*residuals1b[VAR_MOMENTUM+0];
        residuals1b_prolonged[VAR_MOMENTUM+1]     += idist_b1b2*residuals1b[VAR_MOMENTUM+1];
        residuals1b_prolonged[VAR_MOMENTUM+2]     += idist_b1b2*residuals1b[VAR_MOMENTUM+2];
        residuals1b_prolonged[VAR_DENSITY_ENERGY] += idist_b1b2*residuals1b[VAR_DENSITY_ENERGY];
        *residuals1b_prolonged_wsum += idist_b1b2;

        // Calculate contribution of a1 -> b2:
        float dx_a1b2 = coord1a[0] - coord2b[0];
        float dy_a1b2 = coord1a[1] - coord2b[1];
        float dz_a1b2 = coord1a[2] - coord2b[2];

        const float idist_a1b2 = 1.0f/sqrt(dx_a1b2*dx_a1b2 + dy_a1b2*dy_a1b2 + dz_a1b2*dz_a1b2);
        residuals1b_prolonged[VAR_DENSITY]        += idist_a1b2*residuals1b[VAR_DENSITY];
        residuals1b_prolonged[VAR_MOMENTUM+0]     += idist_a1b2*residuals1b[VAR_MOMENTUM+0];
        residuals1b_prolonged[VAR_MOMENTUM+1]     += idist_a1b2*residuals1b[VAR_MOMENTUM+1];
        residuals1b_prolonged[VAR_MOMENTUM+2]     += idist_a1b2*residuals1b[VAR_MOMENTUM+2];
        residuals1b_prolonged[VAR_DENSITY_ENERGY] += idist_a1b2*residuals1b[VAR_DENSITY_ENERGY];
        *residuals1b_prolonged_wsum += idist_a1b2;
    }
}

inline void down_v2_kernel_post(
    const float* residuals1_prolonged, 
    const float* residuals1_prolonged_wsum, 
    const float* residuals2, 
    float* variables2)
{
    // Divide through by weight sum to complete the weighted average started by down_v2_kernel(), 
    // then apply the prolonged residual to grid:
    for (int i=0; i<NVAR; i++) {
        variables2[i] += residuals2[i] - (residuals1_prolonged[i] / (*residuals1_prolonged_wsum));
    }
}

#endif
