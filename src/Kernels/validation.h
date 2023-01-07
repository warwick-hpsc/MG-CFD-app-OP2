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

#ifndef VALIDATION_H
#define VALIDATION_H

#include "utils.h"

inline void residual_kernel(
    const float* old_variable, 
    const float* variable, 
    float* residual)
{
    for (int v=0; v<NVAR; v++) {
        residual[v] = variable[v] - old_variable[v];
    }
}

inline void calc_rms_kernel(
    const float* residual, 
    float* rms)
{
    for (int i=0; i<NVAR; i++) {
        *rms += residual[i]*residual[i];
    }
}

inline void identify_differences(
    const float* test_value,
    const float* master_value, 
    float* difference)
{
    // If floating-point operations have been reordered, then a difference
    // is expected due to rounding-errors, but the difference should
    // be smaller than the following:
    //   1 x 10 ^ ( E - 17 + N )
    // Where E = exponent of master value
    //       N = largest difference in exponents of any floating-point
    //           arithmetic operation performed

    // N represents how many of the least-significant base-10 digits
    // of the floating-point mantissa are allowed to differ due to 
    // FP arithmetic reordering. Its value is guessed as 8, as to set it
    // accurately would require a trace of all floating-point operation
    // outputs during the runs.

    const float acceptable_relative_difference = 10.0e-6f;

    for (int v=0; v<NVAR; v++) {
        float acceptable_difference = master_value[v] * acceptable_relative_difference;
        if (acceptable_difference < 0.0f) {
            acceptable_difference *= -1.0f;
        }

        // Ignore any differences smaller than 3e-19:
        if (acceptable_difference < 3.0e-19f) {
            acceptable_difference = 3.0e-19f;
        }

        float diff = test_value[v] - master_value[v];
        if (diff < 0.0f) {
            diff *= -1.0f;
        }

        if (diff > acceptable_difference) {
            difference[v] = diff;
        } else {
            difference[v] = 0.0f;
        }
    }
}

inline void count_non_zeros(
    const float* value, 
    int* count)
{   
    for (int v=0; v<NVAR; v++) {
        if (value[v] > 0.0f) {
            (*count)++;
        }
    }
}

inline void count_bad_vals(
    const float* value, 
    int* count)
{   
    #if defined(OPENACC) || defined(__HIPSYCL__) || defined(TRISYCL_CL_LANGUAGE_VERSION)
        // OpenACC compilation is complaining about use of isnan()
    #else
        for (int v=0; v<NVAR; v++) {
            if (isnan(value[v]) || isinf(value[v])) {
                *count += 1;
            }
        }
    #endif
}

#endif
