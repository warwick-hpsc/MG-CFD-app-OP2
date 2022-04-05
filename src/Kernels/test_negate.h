// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef NEGATE_H
#define NEGATE_H

#include <math.h>

#include "inlined_funcs.h"

#include "global.h"
#include "config.h"


inline void test_negate_kernel(
    double *variables_a,
    double *variables_b,
    double *flux_a,
    double *flux_b,
    int diff)
{
  // flux_a[VAR_DENSITY] -= (variables_a[VAR_DENSITY] + 1);
  // flux_a[VAR_DENSITY_ENERGY] -= (variables_a[VAR_DENSITY_ENERGY] + 1);
  // flux_a[VAR_MOMENTUM+0] -= (variables_a[VAR_MOMENTUM+0] + 1);
  // flux_a[VAR_MOMENTUM+1] -= (variables_a[VAR_MOMENTUM+1] + 1);
  // flux_a[VAR_MOMENTUM+2] -= (variables_a[VAR_MOMENTUM+2] + 1);

  // flux_b[VAR_DENSITY] -= (variables_b[VAR_DENSITY] + 1);
  // flux_b[VAR_DENSITY_ENERGY] -= (variables_b[VAR_DENSITY_ENERGY] + 1);
  // flux_b[VAR_MOMENTUM+0] -= (variables_b[VAR_MOMENTUM+0] + 1);
  // flux_b[VAR_MOMENTUM+1] -= (variables_b[VAR_MOMENTUM+1] + 1);
  // flux_b[VAR_MOMENTUM+2] -= (variables_b[VAR_MOMENTUM+2] + 1);

  variables_a[VAR_DENSITY] -= diff;
  variables_a[VAR_DENSITY_ENERGY] -= diff;
  variables_a[VAR_MOMENTUM+0] -= diff;
  variables_a[VAR_MOMENTUM+1] -= diff;
  variables_a[VAR_MOMENTUM+2] -= diff;

  variables_b[VAR_DENSITY] -= diff;
  variables_b[VAR_DENSITY_ENERGY] -= diff;
  variables_b[VAR_MOMENTUM+0] -= diff;
  variables_b[VAR_MOMENTUM+1] -= diff;
  variables_b[VAR_MOMENTUM+2] -= diff;

}

#endif
