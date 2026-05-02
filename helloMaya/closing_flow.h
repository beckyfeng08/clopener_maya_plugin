#ifndef CLOSING_FLOW_H
#define CLOSING_FLOW_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "closing_flow_types.h"
#include "closing_flow_helpers.h"
#include "remesh_botsch.h"

class ClosingFlow {
public:
    ClosingFlow(const Eigen::MatrixXd& V_in,
                const Eigen::MatrixXi& F_in,
                const ClosingFlowParams& params);

    // Run one iteration. Returns true if the flow should continue,
    // false if it has converged, the active set is empty, or it was cancelled.
    bool step();

    // Accessors for the host (Maya, viewer, etc.)
    const Eigen::MatrixXd& current_V() const { return Vfull_; }
    const Eigen::MatrixXi& current_F() const { return Ffull_; }
    int  iteration()  const { return iter_; }
    bool converged()  const { return converged_; }

    // Let host tweak params between steps if desired
    ClosingFlowParams& params() { return params_; }
    const ClosingFlowParams& params() const { return params_; }

    // Allow host to request early termination (e.g. Esc in Maya)
    void cancel() { converged_ = true; }

    // Timing accessors (parity with original)
    double remesh_seconds_total() const { return remesh_seconds_total_; }

private:
    // Persistent state across iterations
    Eigen::MatrixXd          Vfull_;
    Eigen::MatrixXi          Ffull_;
    Eigen::SparseMatrix<int> A_;
    Eigen::VectorXd          M_;
    Eigen::VectorXi          moving_;
    bool                     recompute_ = true;
    bool                     converged_ = false;
    int                      iter_ = 0;
    ClosingFlowParams        params_;
    double                   remesh_seconds_total_ = 0.0;
};

// Convenience wrapper that preserves the old API.
bool closing_flow(const Eigen::MatrixXd& V_in,
                  const Eigen::MatrixXi& F_in,
                  const ClosingFlowParams& params,
                  Eigen::MatrixXd& V_out,
                  Eigen::MatrixXi& F_out);

#endif