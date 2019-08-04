#ifndef COMPUTE_NODE_AREA_H
#define COMPUTE_NODE_AREA_H

inline void zero_node_area_kernel(
    double* area)
{
    *area = 0.0;
}

inline void compute_node_area_kernel(
    const double* edge_weights, 
    const double* node1_coords, 
    const double* node2_coords, 
    double* node1_area, 
    double* node2_area)
{
    double ewtx = cl::sycl::fabs(edge_weights[0]);
    double ewty = cl::sycl::fabs(edge_weights[1]);
    double ewtz = cl::sycl::fabs(edge_weights[2]);

    double dx = cl::sycl::fabs(node1_coords[0] - node2_coords[0]);
    double dy = cl::sycl::fabs(node1_coords[1] - node2_coords[1]);
    double dz = cl::sycl::fabs(node1_coords[2] - node2_coords[2]);

    double area = 0.5 * (1.0/3.0) * (dx*ewtx + dy*ewty + dz*ewtz);
    *node1_area += area;
    *node2_area += area;
}

#endif
