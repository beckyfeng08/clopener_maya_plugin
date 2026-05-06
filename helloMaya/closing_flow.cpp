#include "closing_flow.h"
#include "per_face_prin_curvature.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdexcept>

#include <Eigen/Sparse>

#include <igl/adjacency_matrix.h>
#include <igl/avg_edge_length.h>
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/principal_curvature.h>
#include <igl/remove_unreferenced.h>
#include <igl/slice.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/edges.h>
#include <igl/doublearea.h>

namespace cf = closing_flow_detail;

// =============================================================================
// Constructor
// =============================================================================
ClosingFlow::ClosingFlow(const Eigen::MatrixXd& V_in,
    const Eigen::MatrixXi& F_in,
    const ClosingFlowParams& params)
    : Vfull_(V_in), Ffull_(F_in), params_(params)
{
    // Make sure that face and vertex have the correct dimensions (passing in a triangle mesh)
    if (F_in.cols() != 3 || V_in.cols() != 3) {
        throw std::invalid_argument(
            "ClosingFlow: V must be (n,3) and F must be (m,3)");
    }

    avg_edge_ = igl::avg_edge_length(V_in, F_in);
    std::cerr << "avg edge length = " << avg_edge_ << "\n";
    if (params_.verbose) {
        std::cerr << "Shape of Vfull: (" << Vfull_.rows() << ", " << Vfull_.cols() << ")\n";
        std::cerr << "Shape of Ffull: (" << Ffull_.rows() << ", " << Ffull_.cols() << ")\n";
    }

    // recompute_ starts true so the first step() does the full setup
    // (A_, M_, moving_ all get computed in the recompute block below)

    // ---- Convert vertex-index selection to positional selection ----
    // We snapshot the selected vertices' positions now, before any remeshing
    // has changed the vertex set. Step 9 will use spatial proximity to decide
    // which post-remesh vertices belong to the selected region.
    if (params_.selection.size() > 0) {
        const int k = (int)params_.selection.size();
        selection_positions_.resize(k, 3);
        int kept = 0;
        for (int i = 0; i < k; ++i) {
            int idx = params_.selection(i);
            if (idx >= 0 && idx < V_in.rows()) {
                selection_positions_.row(kept++) = V_in.row(idx);
            }
        }
        selection_positions_.conservativeResize(kept, 3);

        // Tolerance: a vertex counts as "in selection" if its squared distance
        // to the nearest original selected vertex is below tol_sq.
        // Using a fraction of the average input edge length keeps this scale-aware.
        double tol = 0.75 * avg_edge_;   // a bit less than one edge
        selection_tol_sq_ = tol * tol;

        if (params_.verbose) {
            std::cerr << "selection: " << kept << " positions stored, "
                << "proximity tolerance = " << tol
                << " (avg edge length = " << avg_edge_ << ")\n";
        }
    }
}

// =============================================================================
// step(): one iteration of the closing/opening flow
// =============================================================================
bool ClosingFlow::step()
{
    if (converged_) return false;

    Eigen::MatrixXd Vprev = Vfull_; // (n, 3)
    Eigen::MatrixXi Fprev = Ffull_; // (m, 3)

    const int nVerts = Vfull_.rows(); // number of total vertices in the mesh
    Eigen::VectorXi frozen(nVerts);   // vertices with 0 flow (dS/dt)

    if (params_.verbose) {
        std::cerr << "=== iter " << iter_ << " ===\n";
        std::cerr << "nVerts: " << nVerts << "\n";

        if (iter_ > 0) {
            std::cerr << "A.rows(): " << A_.rows() << "\n";
            std::cerr << "A.nonZeros(): " << A_.nonZeros() << "\n";
            std::cerr << "moving.size(): " << moving_.size() << "\n";
            if (moving_.size() > 0)
                std::cerr << "moving.maxCoeff(): " << moving_.maxCoeff() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // RECOMPUTE BLOCK: A, M, moving
    // -------------------------------------------------------------------------
    if (recompute_) { // at the very first iteration, or if we are recomputing
        moving_.resize(nVerts);  // resize here instead of declaring

        // Compute A, sparse matrix with adjacency matrix plus identity for self-loops
        igl::adjacency_matrix(Ffull_, A_); // (n, n)
        for (int i = 0; i < Vfull_.rows(); i++) {
            A_.coeffRef(i, i) = 1; // (n, n)
        }

        // mass matrix M
        Eigen::SparseMatrix<double> M_sparse;
        igl::massmatrix(Vfull_, Ffull_, igl::MASSMATRIX_TYPE_VORONOI, M_sparse);

        // diagonal elements and clip
        M_.resize(Vfull_.rows());
        for (int i = 0; i < Vfull_.rows(); i++) {
            M_(i) = std::max(M_sparse.coeff(i, i), 1e-8);
        }

        // Gaussian curvature K
        Eigen::VectorXd K = cf::discrete_gaussian_curvature(Vfull_, Ffull_);

        Eigen::VectorXd H = cf::discrete_mean_curvature(Vfull_, Ffull_);

        Eigen::VectorXd disc = H.array().square() - M_.array() * K.array();
        Eigen::VectorXcd sqrtDisc = disc.cast<std::complex<double>>().array().sqrt();

        Eigen::VectorXcd k_max = H.cast<std::complex<double>>() + sqrtDisc;
        Eigen::VectorXcd k_min = H.cast<std::complex<double>>() - sqrtDisc;

        // Divide by M and take real part
        Eigen::VectorXd curv;
        if (params_.opening) {
            curv = (-k_max.array() / M_.cast<std::complex<double>>().array()).real().matrix();
        }
        else {
            curv = (k_min.array() / M_.cast<std::complex<double>>().array()).real().matrix();
        }

        // Active vertices where curv < -bd AND (if a selection is provided) the
        // vertex is in the user's selection. Empty selection = whole mesh eligible.
        // Closing: k_min < -bd
        // Opening: k_max > bd <-> -k_max < -bd
        const bool use_selection = params_.selection.size() > 0;
        std::vector<bool> in_selection;
        if (use_selection) {
            in_selection.assign(nVerts, false);
            for (int i = 0; i < params_.selection.size(); ++i) {
                int idx = params_.selection(i);
                if (idx >= 0 && idx < nVerts) {
                    in_selection[idx] = true;
                }
            }
            if (params_.verbose) {
                std::cerr << "selection size: " << params_.selection.size()
                    << " (out of " << nVerts << " vertices)\n";
            }
        }

        Eigen::Index nmov = 0, nfroz = 0;
        for (Eigen::Index i = 0; i < nVerts; ++i) {
            bool is_curvy = curv(i) < -params_.bd;
            bool is_allowed = !use_selection || in_selection[i];
            if (is_curvy && is_allowed) {
                moving_(nmov) = static_cast<int>(i);
                nmov++;
            }
            else {
                frozen(nfroz) = static_cast<int>(i);
                nfroz++;
            }
        }
        moving_.conservativeResize(nmov);
        frozen.conservativeResize(nfroz);

        if (params_.verbose) {
            std::cerr << "num active vertices: " << moving_.size() << "\n";
            std::cerr << "num frozen vertices: " << frozen.size() << "\n";
            std::cerr << "active vertices + frozen vertices: "
                << moving_.size() + frozen.size() << "\n";
        }

        if (!params_.always_recompute) {
            recompute_ = false;
        }
    }

    if (moving_.size() == 0) { // we converged, signal stop
        std::cout << "Active set is empty" << std::endl;
        converged_ = true;
        return false;
    }

    // -------------------------------------------------------------------------
    // 2-ring expansion
    // -------------------------------------------------------------------------
    std::vector<bool> isMoving(nVerts, false);
    for (int i = 0; i < moving_.size(); ++i) {
        isMoving[moving_(i)] = true;
    }

    std::vector<bool> ring1 = cf::expandRing(A_, isMoving);
    std::vector<bool> isActive = cf::expandRing(A_, ring1);

    Eigen::Index nfroz2 = 0;
    Eigen::VectorXi frozen_with_2ring(nVerts);
    for (int i = 0; i < nVerts; i++) {
        if (!isActive[i]) {
            frozen_with_2ring(nfroz2) = static_cast<int>(i);
            nfroz2++;
        }
    }
    frozen_with_2ring.conservativeResize(nfroz2);
    if (params_.verbose) {
        std::cerr << "num frozen with 2 ring vertices: " << frozen_with_2ring.size() << "\n";
    }

    if (params_.always_recompute) {
        std::fill(isActive.begin(), isActive.end(), true);
    }

    // ---- Sanity Checks for 2 ring expansion -----
    int nActive = (int)std::count(isActive.begin(), isActive.end(), true);
    if (params_.verbose) {
        std::cerr << "num active (2-ring) vertices: " << nActive << "\n";
        std::cerr << "num moving vertices: " << moving_.size() << "\n";
        std::cerr << "sanity: moving <= active? "
            << (moving_.size() <= (size_t)nActive ? "YES" : "NO") << "\n";
        std::cerr << "nActive after 2-ring: " << nActive << "\n";
        std::cerr << "frozen_with_2ring size: " << frozen_with_2ring.size() << "\n";
        std::cerr << "sanity: nActive + frozen_with_2ring == nVerts? "
            << (nActive + (int)frozen_with_2ring.size() == nVerts ? "YES" : "NO - WARNING") << "\n";
    }

    // -------------------------------------------------------------------------
    // Active and inactive faces
    // -------------------------------------------------------------------------
    std::vector<int> active_vec;
    active_vec.reserve(nVerts);
    for (int i = 0; i < nVerts; ++i) {
        if (isActive[i]) active_vec.push_back(i);
    }
    Eigen::VectorXi active = Eigen::Map<Eigen::VectorXi>(active_vec.data(), (int)active_vec.size());

    // fid_active: face indices where >= 2 vertices are in the active 2-ring
    Eigen::VectorXi fid_active = cf::incident_faces(Ffull_, active);

    if (fid_active.size() == 0) {
        std::cout << "Active set is empty, fid_active.size() == 0" << std::endl;
        converged_ = true;
        return false;
    }

    Eigen::MatrixXi f_active;
    igl::slice(Ffull_, fid_active, 1, f_active);

    // f_inactive: complement of fid_active
    std::vector<bool> isFaceActive(Ffull_.rows(), false);
    for (int i = 0; i < fid_active.size(); i++) {
        isFaceActive[fid_active(i)] = true;
    }

    std::vector<int> fid_inactive_vec;
    fid_inactive_vec.reserve(Ffull_.rows() - fid_active.size());
    for (int f = 0; f < Ffull_.rows(); f++) {
        if (!isFaceActive[f]) fid_inactive_vec.push_back(f);
    }

    Eigen::VectorXi fid_inactive = Eigen::Map<Eigen::VectorXi>(
        fid_inactive_vec.data(), (int)fid_inactive_vec.size());
    Eigen::MatrixXi f_inactive;
    igl::slice(Ffull_, fid_inactive, 1, f_inactive);

    if (params_.verbose) {
        std::cerr << "num active faces (fid_active): " << fid_active.size() << "\n";
        std::cerr << "num inactive faces: " << fid_inactive_vec.size() << "\n";
        std::cerr << "sanity: active + inactive == total? "
            << (fid_active.size() + fid_inactive_vec.size() == (size_t)Ffull_.rows() ? "YES" : "NO") << "\n";
    }

    // -------------------------------------------------------------------------
    // Remove unreferenced vertices from active submesh
    // -------------------------------------------------------------------------
    Eigen::MatrixXd v_active;
    Eigen::MatrixXi f_active_new;
    Eigen::VectorXi I, J;
    igl::remove_unreferenced(Vfull_, f_active, v_active, f_active_new, I, J);
    f_active = f_active_new;

    // fixed_test = setdiff1d(arange(len(v_active)), I[moving])
    std::vector<bool> isMovingLocal(v_active.rows(), false);
    for (int i = 0; i < moving_.size(); ++i) {
        int local_idx = I(moving_(i));
        if (local_idx >= 0)
            isMovingLocal[local_idx] = true;
    }
    std::vector<int> fixed_test_vec;
    fixed_test_vec.reserve(v_active.rows());
    for (int i = 0; i < (int)v_active.rows(); ++i)
        if (!isMovingLocal[i]) fixed_test_vec.push_back(i);
    Eigen::VectorXi fixed_test = Eigen::Map<Eigen::VectorXi>(
        fixed_test_vec.data(), (int)fixed_test_vec.size());

    // Remove unreferenced vertices from inactive submesh
    Eigen::MatrixXd v_inactive;
    Eigen::MatrixXi f_inactive_new;
    Eigen::VectorXi I_inactive, J_inactive;
    igl::remove_unreferenced(Vfull_, f_inactive, v_inactive, f_inactive_new, I_inactive, J_inactive);
    f_inactive = f_inactive_new;

    // V, F are the active submesh
    Eigen::MatrixXd V = v_active;
    Eigen::MatrixXi F = f_active;

    if (params_.verbose) {
        std::cerr << "v_active rows: " << v_active.rows() << "\n";
        std::cerr << "f_active rows: " << f_active.rows() << "\n";
        std::cerr << "I size: " << I.size() << " (should == Vfull.rows() == " << Vfull_.rows() << ")\n";
        std::cerr << "J size: " << J.size() << " (should == v_active.rows() == " << v_active.rows() << ")\n";

        bool I_J_consistent = true;
        for (int i = 0; i < Vfull_.rows(); ++i) {
            if (I(i) >= 0 && J(I(i)) != i) {
                I_J_consistent = false;
                std::cerr << "I/J inconsistency at vertex " << i << "\n";
                break;
            }
        }
        std::cerr << "I/J mapping consistent? " << (I_J_consistent ? "YES" : "NO") << "\n";

        std::cerr << "fixed_test size: " << fixed_test.size() << "\n";
        std::cerr << "isMovingLocal count: "
            << std::count(isMovingLocal.begin(), isMovingLocal.end(), true) << "\n";
        std::cerr << "sanity: fixed_test + moving_local == v_active? "
            << (fixed_test.size() + std::count(isMovingLocal.begin(), isMovingLocal.end(), true)
                == (size_t)v_active.rows() ? "YES" : "NO") << "\n";

        std::cerr << "v_inactive rows: " << v_inactive.rows() << "\n";
        std::cerr << "f_inactive rows: " << f_inactive.rows() << "\n";

        std::cerr << "v_active first vertex: " << v_active.row(0) << "\n";
        std::cerr << "Vfull row J(0):        " << Vfull_.row(J(0)) << "\n";
        std::cerr << "sanity: v_active[0] == Vfull[J[0]]? "
            << (v_active.row(0).isApprox(Vfull_.row(J(0))) ? "YES" : "NO") << "\n";
    }

    // -------------------------------------------------------------------------
    // M restricted to active submesh via J (new->old map)
    // -------------------------------------------------------------------------
    Eigen::VectorXd M_active(v_active.rows());
    for (int i = 0; i < (int)v_active.rows(); ++i)
        M_active(i) = M_(J(i));

    if (params_.verbose) {
        std::cerr << "M_active size: " << M_active.size()
            << " (should == v_active.rows() == " << v_active.rows() << ")\n";
        std::cerr << "M_active min: " << M_active.minCoeff() << " (should be >= 1e-8)\n";
        std::cerr << "M_active max: " << M_active.maxCoeff() << "\n";
        std::cerr << "sanity: M_active size == v_active rows? "
            << (M_active.size() == v_active.rows() ? "YES" : "NO") << "\n";
    }

    // dblA: per-face area weights
    Eigen::VectorXd dblA_vec;
    igl::doublearea(V, F, dblA_vec);
    dblA_vec /= 2.0;

    if (params_.verbose) {
        std::cerr << "dblA_vec size: " << dblA_vec.size()
            << " (should == F.rows() == " << F.rows() << ")\n";
        std::cerr << "dblA_vec min: " << dblA_vec.minCoeff()
            << " (should be > 0, degenerate faces if not)\n";
        std::cerr << "dblA_vec max: " << dblA_vec.maxCoeff() << "\n";
        std::cerr << "sanity: any degenerate faces? "
            << (dblA_vec.minCoeff() <= 0 ? "YES - WARNING" : "NO") << "\n";
    }

    // v_xyz: stack x, y, z coordinates into one long vector [x | y | z]
    int nV_active = (int)V.rows();
    Eigen::VectorXd v_xyz(nV_active * 3);
    v_xyz.segment(0, nV_active) = V.col(0);
    v_xyz.segment(nV_active, nV_active) = V.col(1);
    v_xyz.segment(2 * nV_active, nV_active) = V.col(2);

    if (params_.verbose) {
        std::cerr << "v_xyz size: " << v_xyz.size()
            << " (should == 3 * V.rows() == " << 3 * V.rows() << ")\n";
        std::cerr << "sanity: v_xyz size == 3 * nV_active? "
            << (v_xyz.size() == 3 * nV_active ? "YES" : "NO") << "\n";
        std::cerr << "V(0,:) = " << V.row(0) << "\n";
        std::cerr << "v_xyz[0], v_xyz[nV], v_xyz[2nV] = "
            << v_xyz(0) << ", " << v_xyz(nV_active) << ", " << v_xyz(2 * nV_active) << "\n";
        std::cerr << "sanity: v_xyz stacking correct? "
            << (std::abs(v_xyz(0) - V(0, 0)) < 1e-10 &&
                std::abs(v_xyz(nV_active) - V(0, 1)) < 1e-10 &&
                std::abs(v_xyz(2 * nV_active) - V(0, 2)) < 1e-10 ? "YES" : "NO") << "\n";
    }

    // boundary vertices of active submesh
    Eigen::VectorXi boundary_verts = cf::get_all_boundary_vids(F);
    if (params_.verbose) {
        std::cerr << "boundary_verts size: " << boundary_verts.size() << "\n";
        std::cerr << "sanity: boundary verts < nV_active? "
            << (boundary_verts.size() == 0 || boundary_verts.maxCoeff() < nV_active ? "YES" : "NO") << "\n";
    }

    // fixed = unique(concat(fixed_test, boundary_verts))
    std::vector<int> fixed_vec;
    fixed_vec.reserve(fixed_test.size() + boundary_verts.size());
    for (int i = 0; i < fixed_test.size(); ++i)     fixed_vec.push_back(fixed_test(i));
    for (int i = 0; i < boundary_verts.size(); ++i) fixed_vec.push_back(boundary_verts(i));
    std::sort(fixed_vec.begin(), fixed_vec.end());
    fixed_vec.erase(std::unique(fixed_vec.begin(), fixed_vec.end()), fixed_vec.end());
    Eigen::VectorXi fixed = Eigen::Map<Eigen::VectorXi>(fixed_vec.data(), (int)fixed_vec.size());

    if (params_.verbose) {
        std::cerr << "fixed_test size: " << fixed_test.size() << "\n";
        std::cerr << "boundary_verts size: " << boundary_verts.size() << "\n";
        std::cerr << "fixed size (union): " << fixed.size() << "\n";
        std::cerr << "sanity: fixed <= fixed_test + boundary? "
            << (fixed.size() <= (size_t)(fixed_test.size() + boundary_verts.size()) ? "YES" : "NO") << "\n";
        std::cerr << "sanity: fixed indices in range? "
            << (fixed.size() == 0 ||
                (fixed.minCoeff() >= 0 && fixed.maxCoeff() < nV_active) ? "YES" : "NO") << "\n";
    }

    // fixed_xyz: extend fixed into stacked xyz space
    int nFixed = (int)fixed.size();
    Eigen::VectorXi fixed_xyz(nFixed * 3);
    fixed_xyz.segment(0, nFixed) = fixed;
    fixed_xyz.segment(nFixed, nFixed) = fixed.array() + nV_active;
    fixed_xyz.segment(2 * nFixed, nFixed) = fixed.array() + 2 * nV_active;

    if (params_.verbose) {
        std::cerr << "fixed_xyz size: " << fixed_xyz.size()
            << " (should == 3 * fixed.size() == " << 3 * fixed.size() << ")\n";
        std::cerr << "sanity: fixed_xyz max < 3*nV_active? "
            << (fixed_xyz.maxCoeff() < 3 * nV_active ? "YES" : "NO") << "\n";
    }

    // -------------------------------------------------------------------------
    // Principal curvature directions per face
    // -------------------------------------------------------------------------
    Eigen::MatrixXd PD1, PD2;
    Eigen::VectorXd PC1, PC2;
    pfpc::per_face_prin_curvature(V, F, PD1, PD2, PC1, PC2);
    if (params_.opening)
        PD1 = PD2;

    if (params_.verbose) {
        std::cerr << "PD1 rows: " << PD1.rows() << " cols: " << PD1.cols()
            << " (should be " << F.rows() << " x 3)\n";
        std::cerr << "PD2 rows: " << PD2.rows() << " cols: " << PD2.cols()
            << " (should be " << F.rows() << " x 3)\n";
        std::cerr << "sanity: PD1 rows == F rows? " << (PD1.rows() == F.rows() ? "YES" : "NO") << "\n";
        double pd1_norm_first = PD1.row(0).norm();
        std::cerr << "PD1 row 0 norm: " << pd1_norm_first << " (should be ~1.0)\n";
        std::cerr << "opening mode: " << (params_.opening ? "YES (using PD2)" : "NO (using PD1)") << "\n";
    }

    // proy_matrix = D matrix, per-face dot product with curvature directions
    int nF = F.rows();
    Eigen::SparseMatrix<double> proy_matrix(nF, 3 * nF);
    std::vector<Eigen::Triplet<double>> proy_trips;
    proy_trips.reserve(3 * nF);
    for (int i = 0; i < nF; i++) {
        proy_trips.emplace_back(i, i, PD1(i, 0));  // x block
        proy_trips.emplace_back(i, i + nF, PD1(i, 1));  // y block
        proy_trips.emplace_back(i, i + 2 * nF, PD1(i, 2));  // z block
    }
    proy_matrix.setFromTriplets(proy_trips.begin(), proy_trips.end());

    if (params_.verbose) {
        std::cerr << "proy_matrix shape: (" << proy_matrix.rows() << ", " << proy_matrix.cols() << ")"
            << " (should be nF x 3*nF = " << nF << " x " << 3 * nF << ")\n";
        std::cerr << "sanity proy_matrix rows == nF?    " << (proy_matrix.rows() == nF ? "YES" : "NO") << "\n";
        std::cerr << "sanity proy_matrix cols == 3*nF?  " << (proy_matrix.cols() == 3 * nF ? "YES" : "NO") << "\n";
        double pd1_min_norm = std::numeric_limits<double>::max();
        double pd1_max_norm = 0.0;
        for (int i = 0; i < nF; i++) {
            double n = PD1.row(i).norm();
            pd1_min_norm = std::min(pd1_min_norm, n);
            pd1_max_norm = std::max(pd1_max_norm, n);
        }
        std::cerr << "PD1 row norm range: [" << pd1_min_norm << ", " << pd1_max_norm << "]"
            << " (should be ~1.0 for all rows)\n";
    }

    // G matrix, gradient operator: (3*nF x nV)
    Eigen::SparseMatrix<double> G;
    igl::grad(V, F, G);

    if (params_.verbose) {
        std::cerr << "G shape: (" << G.rows() << ", " << G.cols() << ")"
            << " (should be 3*nF x nV = " << 3 * nF << " x " << nV_active << ")\n";
        std::cerr << "sanity G rows == 3*nF? " << (G.rows() == 3 * nF ? "YES" : "NO") << "\n";
        std::cerr << "sanity G cols == nV?   " << (G.cols() == nV_active ? "YES" : "NO") << "\n";
    }

    // projected_gradient = D * G: (nF x nV)
    Eigen::SparseMatrix<double> projected_gradient = proy_matrix * G;

    if (params_.verbose) {
        std::cerr << "projected_gradient shape: (" << projected_gradient.rows() << ", "
            << projected_gradient.cols() << ")"
            << " (should be nF x nV = " << nF << " x " << nV_active << ")\n";
        std::cerr << "sanity projected_gradient rows == nF?  "
            << (projected_gradient.rows() == nF ? "YES" : "NO") << "\n";
        std::cerr << "sanity projected_gradient cols == nV?  "
            << (projected_gradient.cols() == nV_active ? "YES" : "NO") << "\n";
    }

    // dblA as sparse diagonal matrix
    Eigen::SparseMatrix<double> dblA_sparse(nF, nF);
    std::vector<Eigen::Triplet<double>> area_trips;
    area_trips.reserve(nF);
    for (int i = 0; i < nF; i++)
        area_trips.emplace_back(i, i, dblA_vec(i));
    dblA_sparse.setFromTriplets(area_trips.begin(), area_trips.end());

    if (params_.verbose) {
        std::cerr << "dblA_sparse shape: (" << dblA_sparse.rows() << ", " << dblA_sparse.cols() << ")"
            << " (should be nF x nF = " << nF << " x " << nF << ")\n";
        std::cerr << "dblA min: " << dblA_vec.minCoeff()
            << " (should be > 0, else degenerate faces)\n";
        std::cerr << "sanity no degenerate faces? "
            << (dblA_vec.minCoeff() > 0 ? "YES" : "NO - WARNING") << "\n";
    }

    // int_proy_grad_sq = -Gp^T * dblA * Gp: (nV x nV)
    Eigen::SparseMatrix<double> int_proy_grad_sq =
        -(projected_gradient.transpose() * dblA_sparse * projected_gradient);

    if (params_.verbose) {
        std::cerr << "int_proy_grad_sq shape: (" << int_proy_grad_sq.rows() << ", "
            << int_proy_grad_sq.cols() << ")"
            << " (should be nV x nV = " << nV_active << " x " << nV_active << ")\n";
        std::cerr << "sanity int_proy_grad_sq square? "
            << (int_proy_grad_sq.rows() == int_proy_grad_sq.cols() ? "YES" : "NO") << "\n";
        std::cerr << "sanity int_proy_grad_sq rows == nV? "
            << (int_proy_grad_sq.rows() == nV_active ? "YES" : "NO") << "\n";
    }

    // M as sparse diagonal: block_diag([M, M, M]) shape (3*nV, 3*nV)
    Eigen::SparseMatrix<double> M_diag(nV_active, nV_active);
    std::vector<Eigen::Triplet<double>> m_trips;
    m_trips.reserve(nV_active);
    for (int i = 0; i < nV_active; i++)
        m_trips.emplace_back(i, i, M_active(i));
    M_diag.setFromTriplets(m_trips.begin(), m_trips.end());

    Eigen::SparseMatrix<double> M_block(3 * nV_active, 3 * nV_active);
    std::vector<Eigen::Triplet<double>> mb_trips;
    mb_trips.reserve(3 * nV_active);
    for (int i = 0; i < nV_active; i++) {
        mb_trips.emplace_back(i, i, M_active(i));
        mb_trips.emplace_back(i + nV_active, i + nV_active, M_active(i));
        mb_trips.emplace_back(i + 2 * nV_active, i + 2 * nV_active, M_active(i));
    }
    M_block.setFromTriplets(mb_trips.begin(), mb_trips.end());

    if (params_.verbose) {
        std::cerr << "M_block shape: (" << M_block.rows() << ", " << M_block.cols() << ")"
            << " (should be 3*nV x 3*nV = " << 3 * nV_active << " x " << 3 * nV_active << ")\n";
        std::cerr << "sanity M_block square?       " << (M_block.rows() == M_block.cols() ? "YES" : "NO") << "\n";
        std::cerr << "sanity M_block rows == 3*nV? " << (M_block.rows() == 3 * nV_active ? "YES" : "NO") << "\n";
    }

    // block_diag([int_proy_grad_sq x3]): (3*nV, 3*nV)
    Eigen::SparseMatrix<double> S_block(3 * nV_active, 3 * nV_active);
    std::vector<Eigen::Triplet<double>> s_trips;
    for (int k = 0; k < int_proy_grad_sq.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(int_proy_grad_sq, k); it; ++it) {
            s_trips.emplace_back(it.row(), it.col(), it.value());
            s_trips.emplace_back(it.row() + nV_active, it.col() + nV_active, it.value());
            s_trips.emplace_back(it.row() + 2 * nV_active, it.col() + 2 * nV_active, it.value());
        }
    }
    S_block.setFromTriplets(s_trips.begin(), s_trips.end());

    if (params_.verbose) {
        std::cerr << "S_block shape: (" << S_block.rows() << ", " << S_block.cols() << ")"
            << " (should be 3*nV x 3*nV = " << 3 * nV_active << " x " << 3 * nV_active << ")\n";
        std::cerr << "sanity S_block square?       " << (S_block.rows() == S_block.cols() ? "YES" : "NO") << "\n";
        std::cerr << "sanity S_block rows == 3*nV? " << (S_block.rows() == 3 * nV_active ? "YES" : "NO") << "\n";
    }

    // K = M - tau * G^T * D^T * A * D * G   <-- this is Q
    Eigen::SparseMatrix<double> Q = M_block - 0.01 * params_.dt * S_block;

    if (params_.verbose) {
        std::cerr << "Q shape: (" << Q.rows() << ", " << Q.cols() << ")"
            << " (should be 3*nV x 3*nV = " << 3 * nV_active << " x " << 3 * nV_active << ")\n";
        std::cerr << "sanity Q square?       " << (Q.rows() == Q.cols() ? "YES" : "NO") << "\n";
        std::cerr << "sanity Q rows == 3*nV? " << (Q.rows() == 3 * nV_active ? "YES" : "NO") << "\n";
    }

    // linear = -block_diag(M,M,M) @ v_xyz
    Eigen::VectorXd linear = -(M_block * v_xyz);

    if (params_.verbose) {
        std::cerr << "linear size: " << linear.size()
            << " (should be 3*nV = " << 3 * nV_active << ")\n";
        std::cerr << "sanity linear size == 3*nV? "
            << (linear.size() == 3 * nV_active ? "YES" : "NO") << "\n";
        std::cerr << "v_xyz size: " << v_xyz.size()
            << " (should be 3*nV = " << 3 * nV_active << ")\n";
    }

    // Aeq, Beq: empty equality constraints
    Eigen::SparseMatrix<double> Aeq(0, 3 * nV_active);
    Eigen::VectorXd Beq(0);

    // Fixed values: current positions of fixed vertices in stacked xyz
    Eigen::VectorXd fixed_vals;
    igl::slice(v_xyz, fixed_xyz, fixed_vals);

    if (params_.verbose) {
        std::cerr << "fixed_xyz size: " << fixed_xyz.size()
            << " (should be 3 * fixed.size() = " << 3 * fixed.size() << ")\n";
        std::cerr << "sanity fixed_xyz max < 3*nV? "
            << (fixed_xyz.maxCoeff() < 3 * nV_active ? "YES" : "NO") << "\n";
        std::cerr << "sanity fixed_xyz min >= 0?   "
            << (fixed_xyz.minCoeff() >= 0 ? "YES" : "NO") << "\n";
        bool fixed_xyz_correct = true;
        for (int i = 0; i < fixed.size(); i++) {
            if (fixed_xyz(i) != fixed(i) ||
                fixed_xyz(i + fixed.size()) != fixed(i) + nV_active ||
                fixed_xyz(i + 2 * fixed.size()) != fixed(i) + 2 * nV_active) {
                fixed_xyz_correct = false;
                break;
            }
        }
        std::cerr << "sanity fixed_xyz offsets correct? " << (fixed_xyz_correct ? "YES" : "NO") << "\n";

        std::cerr << "fixed_vals size: " << fixed_vals.size()
            << " (should == fixed_xyz.size() = " << fixed_xyz.size() << ")\n";
        std::cerr << "sanity fixed_vals size == fixed_xyz size? "
            << (fixed_vals.size() == fixed_xyz.size() ? "YES" : "NO") << "\n";
        if (fixed_vals.size() > 0)
            std::cerr << "sanity fixed_vals[0] == v_xyz[fixed_xyz[0]]? "
            << (std::abs(fixed_vals(0) - v_xyz(fixed_xyz(0))) < 1e-10 ? "YES" : "NO") << "\n";
    }

    // -------------------------------------------------------------------------
    // V^(t+1) solution!! (uu reshaped as nx3 matrix U)
    // -------------------------------------------------------------------------
    Eigen::VectorXd uu;
    igl::min_quad_with_fixed_data<double> mqwf_data;
    igl::min_quad_with_fixed_precompute(Q, fixed_xyz, Aeq, true, mqwf_data);
    igl::min_quad_with_fixed_solve(mqwf_data, linear, fixed_vals, Beq, uu);

    Eigen::MatrixXd U(nV_active, 3);
    U.col(0) = uu.segment(0, nV_active);
    U.col(1) = uu.segment(nV_active, nV_active);
    U.col(2) = uu.segment(2 * nV_active, nV_active);

    if (params_.verbose) {
        std::cerr << "uu size: " << uu.size()
            << " (should be 3*nV = " << 3 * nV_active << ")\n";
        std::cerr << "sanity uu size == 3*nV? " << (uu.size() == 3 * nV_active ? "YES" : "NO") << "\n";
        std::cerr << "U shape: (" << U.rows() << ", " << U.cols() << ")"
            << " (should be nV x 3 = " << nV_active << " x 3)\n";
        std::cerr << "sanity U rows == nV? " << (U.rows() == nV_active ? "YES" : "NO") << "\n";
        std::cerr << "sanity U cols == 3?  " << (U.cols() == 3 ? "YES" : "NO") << "\n";

        bool fixed_unchanged = true;
        for (int i = 0; i < fixed.size(); i++) {
            int vi = fixed(i);
            if ((U.row(vi) - v_active.row(vi)).norm() > 1e-6) {
                fixed_unchanged = false;
                std::cerr << "WARNING: fixed vertex " << vi << " moved by "
                    << (U.row(vi) - v_active.row(vi)).norm() << "\n";
                break;
            }
        }
        std::cerr << "sanity fixed vertices unchanged after solve? "
            << (fixed_unchanged ? "YES" : "NO") << "\n";

        Eigen::Vector3d U_min = U.colwise().minCoeff();
        Eigen::Vector3d U_max = U.colwise().maxCoeff();
        Eigen::Vector3d V_min = v_active.colwise().minCoeff();
        Eigen::Vector3d V_max = v_active.colwise().maxCoeff();
        std::cerr << "v_active bbox: [" << V_min.transpose() << "] to [" << V_max.transpose() << "]\n";
        std::cerr << "U bbox:        [" << U_min.transpose() << "] to [" << U_max.transpose() << "]\n";
        double bbox_expansion = (U_max - U_min).norm() / std::max((V_max - V_min).norm(), 1e-8);
        std::cerr << "bbox expansion ratio: " << bbox_expansion
            << " (should be close to 1.0, large values suggest solver diverged)\n";
    }

    // -------------------------------------------------------------------------
    // Remesh
    // -------------------------------------------------------------------------
    int nV_remesh = (int)U.rows();
    // Use percentage of average edge length as the target edge length
    double target_edge_length;
    if (params_.use_relative) {
        target_edge_length = params_.frac_of_avg_edge * avg_edge_;
    }
    else {
        target_edge_length = params_.h;
    }

    Eigen::VectorXd target_vec = Eigen::VectorXd::Constant(nV_remesh, target_edge_length);
    const auto t_remesh_start = std::chrono::steady_clock::now();
    remesh_botsch(U, F, target_vec, params_.remesh_iterations, boundary_verts, true);
    remesh_seconds_total_ +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_remesh_start).count();

    if (params_.verbose) {
        std::cerr << "Finished remeshing active submesh\n";
        std::cerr << "U rows after remesh: " << U.rows() << "\n";
        std::cerr << "F rows after remesh: " << F.rows() << "\n";
    }

    // Compute interior_verts BEFORE merge
    Eigen::VectorXi boundary_verts_remeshed = cf::get_all_boundary_vids(F);
    std::vector<bool> isBoundary(U.rows(), false);
    for (int i = 0; i < boundary_verts_remeshed.size(); i++)
        isBoundary[boundary_verts_remeshed(i)] = true;
    std::vector<int> interior_verts_vec;
    interior_verts_vec.reserve(U.rows());
    for (int i = 0; i < (int)U.rows(); i++)
        if (!isBoundary[i]) interior_verts_vec.push_back(i);
    Eigen::VectorXi interior_verts = Eigen::Map<Eigen::VectorXi>(
        interior_verts_vec.data(), (int)interior_verts_vec.size());

    if (params_.verbose) {
        std::cerr << "boundary verts after remesh: " << boundary_verts_remeshed.size() << "\n";
        std::cerr << "interior verts after remesh: " << interior_verts.size() << "\n";
    }

    // Udup = U
    Eigen::MatrixXd Udup = U;

    // -------------------------------------------------------------------------
    // Merge active and inactive surface
    // -------------------------------------------------------------------------
    int n_inactive_verts = (int)v_inactive.rows();
    int n_inactive_faces = (int)f_inactive.rows();
    int n_active_verts = (int)Udup.rows();
    int n_active_faces = (int)F.rows();

    Eigen::MatrixXd Vfull_new(n_inactive_verts + n_active_verts, 3);
    Vfull_new.topRows(n_inactive_verts) = v_inactive;
    Vfull_new.bottomRows(n_active_verts) = Udup;

    Eigen::MatrixXi Ffull_new(n_inactive_faces + n_active_faces, 3);
    Ffull_new.topRows(n_inactive_faces) = f_inactive;
    Ffull_new.bottomRows(n_active_faces) = F.array() + n_inactive_verts;

    if (params_.verbose) {
        std::cerr << "Vfull_new rows before dedup: " << Vfull_new.rows() << "\n";
        std::cerr << "Ffull_new rows before dedup: " << Ffull_new.rows() << "\n";
    }

    Eigen::VectorXi SVI, SVJ;
    igl::remove_duplicate_vertices(Vfull_new, Ffull_new, 0, Vfull_, SVI, SVJ, Ffull_);

    if (params_.verbose) {
        std::cerr << "Vfull rows after dedup: " << Vfull_.rows() << "\n";
        std::cerr << "Ffull rows after dedup: " << Ffull_.rows() << "\n";
        std::cerr << "SVI size: " << SVI.size() << " (should == Vfull.rows())\n";
        std::cerr << "SVJ size: " << SVJ.size() << " (should == Vfull_new.rows())\n";
        std::cerr << "sanity SVI size == Vfull rows? "
            << (SVI.size() == Vfull_.rows() ? "YES" : "NO") << "\n";
        std::cerr << "sanity SVJ size == Vfull_new rows? "
            << (SVJ.size() == Vfull_new.rows() ? "YES" : "NO") << "\n";
    }

    // ── Step 4: Map interior verts into full mesh ─────────────────────────────
    Eigen::VectorXi interior_active_full(interior_verts.size());
    for (int i = 0; i < interior_verts.size(); i++)
        interior_active_full(i) = SVJ(interior_verts(i) + n_inactive_verts);

    if (params_.verbose) {
        std::cerr << "interior_active_full size: " << interior_active_full.size() << "\n";
        std::cerr << "sanity interior_active_full max < Vfull rows? "
            << (interior_active_full.size() == 0 ||
                interior_active_full.maxCoeff() < Vfull_.rows() ? "YES" : "NO - WARNING") << "\n";
    }

    // ── Step 5: Rebuild A from active region edges ────────────────────────────
    int nVfull_post = (int)Vfull_.rows();
    Eigen::MatrixXi F_global = F.array() + n_inactive_verts;
    Eigen::MatrixXi edges_active;
    igl::edges(F_global, edges_active);

    A_.resize(nVfull_post, nVfull_post);
    std::vector<Eigen::Triplet<int>> a_trips;
    a_trips.reserve(edges_active.rows() * 2 + nVfull_post);
    for (int i = 0; i < edges_active.rows(); i++) {
        int u = SVJ(edges_active(i, 0));
        int v = SVJ(edges_active(i, 1));
        a_trips.emplace_back(u, v, 1);
        a_trips.emplace_back(v, u, 1);
    }
    for (int i = 0; i < nVfull_post; i++)
        a_trips.emplace_back(i, i, 1);
    A_.setFromTriplets(a_trips.begin(), a_trips.end());

    if (params_.verbose) {
        std::cerr << "A rebuilt (active edges via SVJ), size: "
            << A_.rows() << " x " << A_.cols() << "\n";
        std::cerr << "sanity: A.rows() == Vfull.rows()? "
            << (A_.rows() == Vfull_.rows() ? "YES" : "NO - WARNING") << "\n";
    }

    // ── Step 6 ────────────────────────────────────────────────────────────────
    Eigen::VectorXd H_interior_active = cf::discrete_mean_curvature(U, F);
    Eigen::VectorXd K_interior_active = cf::discrete_gaussian_curvature(U, F);

    // ── Step 7 ────────────────────────────────────────────────────────────────
    Eigen::SparseMatrix<double> M_active_sparse;
    igl::massmatrix(U, F, igl::MASSMATRIX_TYPE_VORONOI, M_active_sparse);
    Eigen::VectorXd m_active(U.rows());
    for (int i = 0; i < (int)U.rows(); i++)
        m_active(i) = M_active_sparse.coeff(i, i);

    M_.resize(Vfull_.rows());
    M_.setOnes();
    for (int i = 0; i < interior_verts.size(); i++) {
        int full_idx = interior_active_full(i);
        int local_idx = interior_verts(i);
        M_(full_idx) = M_(full_idx) + m_active(local_idx) - 1.0;
    }
    for (int i = 0; i < M_.size(); i++)
        M_(i) = std::max(M_(i), 1e-8);

    if (params_.verbose) {
        std::cerr << "M updated, size: " << M_.size()
            << " (should == Vfull.rows() == " << Vfull_.rows() << ")\n";
    }

    // ── Step 8 ────────────────────────────────────────────────────────────────
    Eigen::VectorXd m_active_clipped(U.rows());
    for (int i = 0; i < (int)U.rows(); i++)
        m_active_clipped(i) = std::max(m_active(i), 1e-8);

    Eigen::VectorXd disc_active = H_interior_active.array().square()
        - m_active_clipped.array() * K_interior_active.array();
    Eigen::VectorXcd sqrtDisc_active = disc_active.cast<std::complex<double>>().array().sqrt();
    Eigen::VectorXcd k_max_active = H_interior_active.cast<std::complex<double>>() + sqrtDisc_active;
    Eigen::VectorXcd k_min_active = H_interior_active.cast<std::complex<double>>() - sqrtDisc_active;

    Eigen::VectorXd K_interior_active_curv;
    if (params_.opening) {
        K_interior_active_curv = (-k_max_active.array() /
            m_active_clipped.cast<std::complex<double>>().array()).real().matrix();
    }
    else {
        K_interior_active_curv = (k_min_active.array() /
            m_active_clipped.cast<std::complex<double>>().array()).real().matrix();
    }

    // ── Step 9 ────────────────────────────────────────────────────────────────
    // Rebuild moving_ from post-remesh interior vertices that (a) have high
    // curvature AND (b) lie within the user's selected region (if any).
    //
    // (b) is checked geometrically because remeshing has invalidated the
    // original vertex indexing. We compare each candidate's *position* to the
    // snapshot of selected positions taken at construction time.
    const bool use_selection = selection_positions_.rows() > 0;

    // Pre-compute "is this Vfull_ vertex inside the selected region?" once.
    // We only need the answer for vertices in interior_active_full, but
    // computing for all of Vfull_ once and indexing into it is cleaner and
    // the cost is dominated by the kNN lookup itself.
    std::vector<bool> in_selection;
    if (use_selection) {
        in_selection.assign(Vfull_.rows(), false);

        // For each interior_active_full vertex, find squared distance to the
        // nearest stored selection position. We use point_mesh_squared_distance
        // with a degenerate "mesh" of just the selection points (no faces)
        // — but that API needs faces, so we use a direct loop instead.
        // For typical selection sizes (hundreds to low thousands) this is fast.
        //
        // If your selections get huge (10k+) and this becomes a hotspot,
        // build a KD-tree once with nanoflann or igl::octree.
        for (int i = 0; i < interior_active_full.size(); ++i) {
            int full_idx = interior_active_full(i);
            Eigen::RowVector3d p = Vfull_.row(full_idx);

            double min_d2 = std::numeric_limits<double>::max();
            for (int j = 0; j < selection_positions_.rows(); ++j) {
                double d2 = (p - selection_positions_.row(j)).squaredNorm();
                if (d2 < min_d2) {
                    min_d2 = d2;
                    if (min_d2 < selection_tol_sq_) break;  // early out
                }
            }
            if (min_d2 < selection_tol_sq_)
                in_selection[full_idx] = true;
        }
    }

    std::vector<bool> moving_full(Vfull_.rows(), false);
    for (int i = 0; i < interior_verts.size(); i++) {
        int local_idx = interior_verts(i);
        int full_idx = interior_active_full(i);
        bool is_curvy = K_interior_active_curv(local_idx) < -params_.bd;
        bool is_allowed = !use_selection || in_selection[full_idx];
        if (is_curvy && is_allowed)
            moving_full[full_idx] = true;
    }
    std::vector<int> moving_vec;
    moving_vec.reserve(Vfull_.rows());
    for (int i = 0; i < (int)Vfull_.rows(); i++)
        if (moving_full[i]) moving_vec.push_back(i);
    moving_.resize(moving_vec.size());
    for (int i = 0; i < (int)moving_vec.size(); i++)
        moving_(i) = moving_vec[i];

    if (params_.verbose) {
        std::cerr << "moving size after update: " << moving_.size() << "\n";
        if (use_selection) {
            int n_in_sel = (int)std::count(in_selection.begin(), in_selection.end(), true);
            std::cerr << "  vertices inside selection region: " << n_in_sel << "\n";
        }
    }

    // ── Step 10: convergence check ────────────────────────────────────────────
    if ((iter_ + 1) % 10 == 0) {
        Eigen::VectorXd sqr_dist;
        Eigen::VectorXi closest_faces;
        Eigen::MatrixXd closest_points;
        igl::point_mesh_squared_distance(Vfull_, Vprev, Fprev,
            sqr_dist, closest_faces, closest_points);
        double max_dist = sqr_dist.maxCoeff();
        if (params_.verbose) {
            std::cerr << "iter " << iter_ + 1 << " max squared dist from prev: " << max_dist << "\n";
        }
        if (max_dist < params_.tol) {
            if (params_.verbose) {
                std::cerr << "Converged at iter " << iter_ + 1 << "\n";
            }
            ++iter_;
            converged_ = true;
            return false;
        }
    }

    ++iter_;
    return true;
}

// =============================================================================
// Backwards-compatible wrapper
// =============================================================================
bool closing_flow(const Eigen::MatrixXd& V_in,
    const Eigen::MatrixXi& F_in,
    const ClosingFlowParams& params,
    Eigen::MatrixXd& V_out,
    Eigen::MatrixXi& F_out)
{
    if (F_in.cols() != 3 || V_in.cols() != 3)
        return false;

    const auto t_flow_start = std::chrono::steady_clock::now();

    try {
        ClosingFlow flow(V_in, F_in, params);
        for (int i = 0; i < params.maxiter; ++i) {
            if (!flow.step()) break;
        }
        V_out = flow.current_V();
        F_out = flow.current_F();

        const double flow_seconds_total =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_flow_start).count();
        std::cerr << "closing_flow total time: " << flow_seconds_total << " s\n";
        std::cerr << "remeshing total time: " << flow.remesh_seconds_total() << " s\n";

        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "closing_flow failed: " << e.what() << "\n";
        return false;
    }
}