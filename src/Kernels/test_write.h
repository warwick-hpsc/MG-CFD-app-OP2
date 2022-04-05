// Copyright 2009, Andrew Corrigan, acorriga@gmu.edu
// This code is from the AIAA-2009-4001 paper

#ifndef WRITE_H
#define WRITE_H

#include <math.h>

#include "inlined_funcs.h"

#include "global.h"
#include "config.h"

inline void test_write_kernel(
    double *variables_a,
    double *variables_b)
{
  variables_a[VAR_DENSITY] += 1;
  variables_a[VAR_DENSITY_ENERGY] += 1;
  variables_a[VAR_MOMENTUM+0] += 1;
  variables_a[VAR_MOMENTUM+1] += 1;
  variables_a[VAR_MOMENTUM+2] += 1;

  variables_b[VAR_DENSITY] += 1;
  variables_b[VAR_DENSITY_ENERGY] += 1;
  variables_b[VAR_MOMENTUM+0] += 1;
  variables_b[VAR_MOMENTUM+1] += 1;
  variables_b[VAR_MOMENTUM+2] += 1;

}

#endif
