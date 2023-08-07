#ifndef STRUCTURES_H
#define STRUCTURES_H

#if !defined(__HPSYCL__) && !defined(__CUDACC__) && !defined(__HIPCC__)
  struct double3 { double x, y, z; };
#else
  // nvcc will pull in '/usr/include/vector_types.h' which 
  // contains an identical definition of 'double3' (above).
#endif

#define DEBUGGABLE_ABORT fprintf(stderr, "%s:%d\n", __FILE__, __LINE__); fflush(stderr); fflush(stdout); exit(EXIT_FAILURE);
struct edge { double a, b; };
struct edge_neighbour { int a, b; double x, y, z; };

#endif
