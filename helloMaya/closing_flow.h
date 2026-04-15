#pragma once

#include "closing_flow_types.h"
#include "closing_flow_helpers.h"
#include "remesh_botsch.h"
#include <Eigen/Dense>

// Input/output: triangle mesh only (#V x 3, #F x 3)
bool closing_flow(const Eigen::MatrixXd &V_in, const Eigen::MatrixXi &F_in,
                  const ClosingFlowParams &params, Eigen::MatrixXd &V_out,
                  Eigen::MatrixXi &F_out);
