#ifndef COMPUTE_NODE_AREA_H
#define COMPUTE_NODE_AREA_H

inline void zero_node_area_kernel(
    float* area)
{
    *area = 0.0;
}

inline void compute_node_area_kernel(
    const float* edge_weights, 
    const float* node1_coords, 
    const float* node2_coords, 
    float* node1_area, 
    float* node2_area)
{
  float ewtx = fabs(edge_weights[0]);
  float ewty = fabs(edge_weights[1]);
  float ewtz = fabs(edge_weights[2]);

  float dx = fabs(node1_coords[0] - node2_coords[0]);
  float dy = fabs(node1_coords[1] - node2_coords[1]);
  float dz = fabs(node1_coords[2] - node2_coords[2]);

  float area = 0.5f * (1.0f / 3.0f) * (dx * ewtx + dy * ewty + dz * ewtz);
  *node1_area += area;
  *node2_area += area;
}

#endif
