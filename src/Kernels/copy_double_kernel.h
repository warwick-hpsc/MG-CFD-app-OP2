#ifndef COPY_DOUBLE_KERNEL_H
#define COPY_DOUBLE_KERNEL_H

#include "const.h"

inline void copy_double_kernel(
	const float* variables, 
	float* old_variables)
{
	for (int i=0; i<NVAR; i++) {
		old_variables[i] = variables[i];
	}
}

#endif
