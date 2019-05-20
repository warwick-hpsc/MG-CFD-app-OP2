#ifndef STRUCTURES_H
#define STRUCTURES_H

#ifndef __CUDACC__
  struct double3 { double x, y, z; };
#else
  // nvcc will pull in '/usr/include/vector_types.h' which 
  // contains an identical definition of 'double3' (above).
#endif

struct edge { double a, b; };
struct edge_neighbour { int a, b; double x, y, z; };

#endif
