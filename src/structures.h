#include <vector>

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

struct unit{
  char type;//either M for MG-CFD or C for Coupler unit
  int processes;
  char coupling_type; //either S for sliding plane or C for CHT 
  std::vector< std::vector<int> > mgcfd_ranks;//each coupler unit has 2 MG-CFD instances, MG-CFD units will only have themselves
  std::vector< std::vector<int> > coupler_ranks;//MG-CFD instances can have multiple couplers
  std::vector<int> mgcfd_units;//the 2 MG-CFD units coupled by this coupler unit
};

struct locators{//one locator per rank - states type and relative position
  char typelocator;//either M for MG-CFD or C for Coupler unit
  int placelocator;//e.g 1 for 1st coupler/mgcfd unit, 2 for 2nd coupler/mgcfd unit
};
#endif

