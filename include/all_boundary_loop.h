#ifndef ALL_BOUNDARY_LOOP_H
#define ALL_BOUNDARY_LOOP_H

#include <Eigen/Core>
#include <vector>

// Compute ordered boundary loops for a manifold triangle mesh
// F : #F x 3 matrix of triangle vertex indices
// Returns:
//   vector of loops, where each loop is an ordered list of vertex indices
std::vector<std::vector<int>> all_boundary_loop(const Eigen::MatrixXi& F);

#endif