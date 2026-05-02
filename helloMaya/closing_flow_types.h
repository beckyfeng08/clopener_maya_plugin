#pragma once

struct ClosingFlowParams {
    int maxiter = 50;
    double dt = 0.07;
    // double bd = 1.0 / 0.08;
    double bd = 1.0 / 0.08;

    /// Target edge length for Botsch-Kobbelt remeshing. If <= 0, uses
    /// 0.5 * avg_edge_length (adaptive to mesh scale)
    double h = 0.03;

    /// Remesh iterations per flow step
    int remesh_iterations = 1;

    bool opening = false;
    bool always_recompute = false;
    double tol = 1e-5;
    bool quadric_curvature = false;
    bool verbose = false; 
};