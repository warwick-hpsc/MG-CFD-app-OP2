#ifndef STRUCTURES_H
#define STRUCTURES_H

#if !defined(__HIPSYCL__) && !defined(__CUDACC__) && !defined(__HIPCC__)
  struct float3 { double x, y, z; };
#else
  // nvcc will pull in '/usr/include/vector_types.h' which 
  // contains an identical definition of 'float3' (above).
#endif

#define DEBUGGABLE_ABORT fprintf(stderr, "%s:%d\n", __FILE__, __LINE__); fflush(stderr); fflush(stdout); exit(EXIT_FAILURE);
struct edge { float a, b; };
struct edge_neighbour { int a, b; float x, y, z; };

#endif
