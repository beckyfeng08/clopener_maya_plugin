#pragma once

#include <Eigen/Sparse>

namespace closing_flow_detail {

Eigen::MatrixXd dihedral_angles(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::MatrixXi& TT_out);

Eigen::VectorXd discrete_mean_curvature(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F);

Eigen::VectorXi get_all_boundary_vids(
    const Eigen::MatrixXi& F);

Eigen::VectorXd discrete_gaussian_curvature(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F);

Eigen::VectorXi incident_faces(
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& vids,
    int v_incidence = 2);

std::vector<bool> expandRing(
    const Eigen::SparseMatrix<int>& A, 
    const std::vector<bool>& inSet);

}; // namespace closing_flow_detail