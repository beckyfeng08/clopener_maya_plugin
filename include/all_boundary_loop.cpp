#include "all_boundary_loop.h"
#include <igl/boundary_loop.h>

std::vector<std::vector<int>> all_boundary_loop(const Eigen::MatrixXi& F)
{
    std::vector<std::vector<int>> loops;

    // libigl computes ordered boundary loops
    igl::boundary_loop(F, loops);

    return loops;
}