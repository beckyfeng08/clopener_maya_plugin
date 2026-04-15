#include "closing_flow.h"
#include "per_face_prin_curvature.h"

#include <algorithm>
#include <chrono>
#include <iostream>

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

namespace cf = closing_flow_detail;

bool closing_flow(
    const Eigen::MatrixXd &V_in,
    const Eigen::MatrixXi &F_in,
    const ClosingFlowParams &params,
    Eigen::MatrixXd &V_out,
    Eigen::MatrixXi &F_out)
{
    // Make sure that face and vertex have the correct dimensions (passing in a triangle mesh)
    // Make sure that face and vertex have the correct dimensions (passing in a triangle mesh)
    if (F_in.cols() != 3 || V_in.cols() != 3)
        return false;

    // start timer
    // start timer
    const auto t_flow_start = std::chrono::steady_clock::now();
    double remesh_seconds_total = 0;

    // create a copy of input vertices and faces
    // create a copy of input vertices and faces
    Eigen::MatrixXd Vfull = V_in;
    Eigen::MatrixXi Ffull = F_in;

    std::cerr << "Shape of Vfull: (" << Vfull.rows() << ", " << Vfull.cols()
              << ")\n";
    std::cerr << "Shape of Ffull: (" << Ffull.rows() << ", " << Ffull.cols()
              << ")\n";

    bool recompute = true;

    Eigen::SparseMatrix<int> A;
    Eigen::VectorXd M;
    Eigen::VectorXi moving;

    // Start the integration thing over params.maxiter timesteps
    // Start the integration thing over params.maxiter timesteps
    // !!!!!!!!!!!!! NOTE: MAXITER IS HARDCODED TO 1 !!!!!!!!!!!!!
    // When we actually run the algorithm it will be like idk 300 or something big
    for (int iter=0; iter < params.maxiter; iter++) {
        
        Eigen::MatrixXd Vprev = Vfull; // (n, 3)
        Eigen::MatrixXi Fprev = Ffull; // (m, 3)

        const int nVerts = Vfull.rows(); // number of total vertices in the mesh
        // Eigen::VectorXi moving(nVerts); // "active vertices" - those whose min. curvature k < 1/(-r) for closing
        Eigen::VectorXi frozen(nVerts); // vertices with 0 flow (dS/dt)

        // Eigen::SparseMatrix<int> A;
        // Eigen::VectorXd M;

        std::cerr << "=== iter " << iter << " ===\n";
        std::cerr << "nVerts: " << nVerts << "\n";
        
        if (iter > 0) {
            std::cerr << "A.rows(): " << A.rows() << "\n";
            std::cerr << "A.nonZeros(): " << A.nonZeros() << "\n";
            std::cerr << "moving.size(): " << moving.size() << "\n";
            if (moving.size() > 0)
                std::cerr << "moving.maxCoeff(): " << moving.maxCoeff() << "\n";
        }


        
        // if (iter == 1 || recompute) {
        if (recompute) { // at the very first iteration, or if we are recomputing
            moving.resize(nVerts);  // resize here instead of declaring
            
            // Compute A, sparse matrix with adjacency matrix plus identity for self-loops
            
            igl::adjacency_matrix(Ffull, A); // (n, n)
            for (int i = 0; i < Vfull.rows(); i++) {
                A.coeffRef(i, i) = 1; // (n, n)
            }

            // mass matrix M
            Eigen::SparseMatrix<double> M_sparse;
            igl::massmatrix(Vfull, Ffull, igl::MASSMATRIX_TYPE_VORONOI, M_sparse);

            // diagonal elements and clip
            M.resize(Vfull.rows());
            for (int i = 0; i < Vfull.rows(); i++) {
                M(i) = std::max(M_sparse.coeff(i, i), 1e-8);
            }

            // Gaussian curvature K (both methods produce same number vertices - result from one test)
            // Method 1: built-in libigl function
            // Eigen::VectorXd K;
            // igl::gaussian_curvature(Vfull, Ffull, K);
            // Method 2: match python 
            Eigen::VectorXd K = cf::discrete_gaussian_curvature(Vfull, Ffull);
            
            Eigen::VectorXd H = cf::discrete_mean_curvature(Vfull, Ffull);

            Eigen::VectorXd disc = H.array().square() - M.array() * K.array();
            Eigen::VectorXcd sqrtDisc = disc.cast<std::complex<double>>().array().sqrt();

            Eigen::VectorXcd k_max = H.cast<std::complex<double>>() + sqrtDisc;
            Eigen::VectorXcd k_min = H.cast<std::complex<double>>() - sqrtDisc;

            // Divide by M and take real part (matches Python's np.real(...))
            Eigen::VectorXd curv;
            if (params.opening) {
                curv = (-k_max.array() / M.cast<std::complex<double>>().array()).real().matrix();
            } else {
                curv = (k_min.array() / M.cast<std::complex<double>>().array()).real().matrix();
            }


            // if (params.quadric_curvature) {
            //     // Compute principal curvatures using quadric fitting (taken from https://github.com/libigl/libigl/blob/main/tutorial/203_CurvatureDirections/main.cpp)
            //     // Compute curvature directions via quadric fitting. I think this yields better looking results but you end up remeshing more parts so it takes a performance
            //     Eigen::MatrixXd PD1, PD2; // PD1 - V by 3 maximal curvature direction for each vertex, PD2 - " minimal "
            //     igl::principal_curvature(Vfull, Ffull, PD1, PD2, k_max, k_min);
            //     H = 0.5*(k_max+k_min);
            // }




            // Active vertices where curv < -bd
            // Closing: k_min < -bd
            // Opening: k_max > bd <-> -k_max < -bd
            // allocate at most #V entries then shrink at end so we 
            // don't have to loop through all vertices to form frozen 
            // set dtype Eigen::VectorXi for remeshing
            Eigen::Index nmov = 0, nfroz = 0;
            for (Eigen::Index i = 0; i < nVerts; ++i) {
                if (curv(i) < -params.bd) {
                    moving(nmov) = static_cast<int>(i);
                    nmov++;
                } else {
                    frozen(nfroz) = static_cast<int>(i);
                    nfroz++;
                }
            }
            moving.conservativeResize(nmov);
            frozen.conservativeResize(nfroz);
            std::cerr << "num active vertices: " << moving.size() << "\n";
            std::cerr << "num frozen vertices: " << frozen.size() << "\n";
            std::cerr << "active vertices + frozen vertices: " << moving.size() + frozen.size()
                      << "\n";

            if (!params.always_recompute) {
                recompute = false;
            }
        }

        if (moving.size() == 0) { // we converged, break out of loop
            std::cout << "Active set is empty" << std::endl;
            break;
        }

        // ~~~~~~~~~~~~ 2-ring expansion ~~~~~~~~~~~~~~~~~~~
        // Get 2-ring expansion (add a 2-triangle deep buffer out from the moving vertices)
        std::vector<bool> isMoving(nVerts, false);
        for (int i = 0; i < moving.size(); ++i) {
            isMoving[moving(i)] = true;
        }

        std::vector<bool> ring1 = cf::expandRing(A, isMoving);
        std::vector<bool> isActive = cf::expandRing(A, ring1);

        Eigen::Index nfroz2 = 0;
        Eigen::VectorXi frozen_with_2ring(nVerts);
        for (int i=0; i<nVerts; i++) {
            if (!isActive[i]){
                frozen_with_2ring(nfroz2) = static_cast<int>(i);
                nfroz2++;
            }
        }
        frozen_with_2ring.conservativeResize(nfroz2);
        std::cerr << "num frozen with 2 ring vertices: " << frozen_with_2ring.size() << "\n";

        if (params.always_recompute) {
            std::fill(isActive.begin(), isActive.end(), true);
        }

        // ~~~~~~~~~~~~ Sanity Checks for 2 ring expansion ~~~~~~~~~
        int nActive = (int)std::count(isActive.begin(), isActive.end(), true);
        std::cerr << "num active (2-ring) vertices: " << nActive << "\n";
        std::cerr << "num moving vertices: " << moving.size() << "\n";
        std::cerr << "sanity: moving <= active? " << (moving.size() <= (size_t)nActive ? "YES" : "NO") << "\n";
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        std::cerr << "nActive after 2-ring: " << nActive << "\n";
        std::cerr << "frozen_with_2ring size: " << frozen_with_2ring.size() << "\n";
        std::cerr << "sanity: nActive + frozen_with_2ring == nVerts? "
                  << (nActive + (int)frozen_with_2ring.size() == nVerts ? "YES" : "NO - WARNING") << "\n";



            
        // ~~~~~~~~~~~~~~ Active and inactive faces ~~~~~~~~~~~~~
        // Convert isActive bool vector to index vector for incident_faces
        // Speed this up later? (construct active_vec when we are constructing isActive)
        std::vector<int> active_vec;
        active_vec.reserve(nVerts);
        for (int i = 0; i < nVerts; ++i) {
            if (isActive[i]) active_vec.push_back(i);
        }
        Eigen::VectorXi active = Eigen::Map<Eigen::VectorXi>(active_vec.data(), (int)active_vec.size());

        // fid_active: face indices where >= 2 vertices are in the active 2-ring
        Eigen::VectorXi fid_active = cf::incident_faces(Ffull, active);

        if (fid_active.size() == 0) {
            std::cout << "Active set is empty, fid_active.size() == 0" << std::endl;
            break;
        }

        // f_active: index into Ffull with fid_active
        // Eigen::MatrixXi f_active = Ffull(fid_active, Eigen::all);
        Eigen::MatrixXi f_active;
        igl::slice(Ffull, fid_active, 1, f_active);

        // f_inactive: complement of fid_active
        std::vector<bool> isFaceActive(Ffull.rows(), false);
        for (int i = 0; i < fid_active.size(); i++) {
            isFaceActive[fid_active(i)] = true;
        }

        std::vector<int> fid_inactive_vec;
        fid_inactive_vec.reserve(Ffull.rows() - fid_active.size());
        for (int f = 0; f < Ffull.rows(); f++) {
            if (!isFaceActive[f]) fid_inactive_vec.push_back(f);
        }

        Eigen::VectorXi fid_inactive = Eigen::Map<Eigen::VectorXi>(fid_inactive_vec.data(), (int)fid_inactive_vec.size());
        // Eigen::MatrixXi f_inactive = Ffull(fid_inactive, Eigen::all);
        Eigen::MatrixXi f_inactive;
        igl::slice(Ffull, fid_inactive, 1, f_inactive);

        // ~~~~~~~~~~~~~ Sanity checks for faces ~~~~~~~~~~~~~~~~~
        std::cerr << "num active faces (fid_active): " << fid_active.size() << "\n";
        std::cerr << "num inactive faces: " << fid_inactive_vec.size() << "\n";
        std::cerr << "sanity: active + inactive == total? " 
                << (fid_active.size() + fid_inactive_vec.size() == (size_t)Ffull.rows() ? "YES" : "NO") << "\n";
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Remove unreferenced vertices from active submesh
        // I: old->new index map (-1 if unreferenced), J: new->old index map
        Eigen::MatrixXd v_active;
        Eigen::MatrixXi f_active_new;
        Eigen::VectorXi I, J;
        igl::remove_unreferenced(Vfull, f_active, v_active, f_active_new, I, J);
        f_active = f_active_new;

        // fixed_test = setdiff1d(arange(len(v_active)), I[moving])
        // I[moving] maps global moving indices to local active submesh indices
        // fixed_test = all local vertices that are NOT moving
        std::vector<bool> isMovingLocal(v_active.rows(), false);
        for (int i = 0; i < moving.size(); ++i) {
            int local_idx = I(moving(i));
            if (local_idx >= 0)  // -1 means this moving vertex wasn't in the active submesh
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
        igl::remove_unreferenced(Vfull, f_inactive, v_inactive, f_inactive_new, I_inactive, J_inactive);
        f_inactive = f_inactive_new;

        // V, F are the active submesh
        Eigen::MatrixXd V = v_active;
        Eigen::MatrixXi F = f_active;

        // ~~~~~~~~~~~~~~~~~~~~~ Sanity checks for remove unreferenced faces ~~~~~~~~~~

        std::cerr << "v_active rows: " << v_active.rows() << "\n";
        std::cerr << "f_active rows: " << f_active.rows() << "\n";
        std::cerr << "I size: " << I.size() << " (should == Vfull.rows() == " << Vfull.rows() << ")\n";
        std::cerr << "J size: " << J.size() << " (should == v_active.rows() == " << v_active.rows() << ")\n";
        
        // Verify I maps correctly: J[I[v]] == v for all referenced vertices
        bool I_J_consistent = true;
        for (int i = 0; i < Vfull.rows(); ++i) {
            if (I(i) >= 0 && J(I(i)) != i) {
                I_J_consistent = false;
                std::cerr << "I/J inconsistency at vertex " << i << "\n";
                break;
            }
        }
        std::cerr << "I/J mapping consistent? " << (I_J_consistent ? "YES" : "NO") << "\n";

        // After fixed_test computation
        std::cerr << "fixed_test size: " << fixed_test.size() << "\n";
        std::cerr << "isMovingLocal count: " << std::count(isMovingLocal.begin(), isMovingLocal.end(), true) << "\n";
        std::cerr << "sanity: fixed_test + moving_local == v_active? "
                << (fixed_test.size() + std::count(isMovingLocal.begin(), isMovingLocal.end(), true) == (size_t)v_active.rows() ? "YES" : "NO") << "\n";

        // After remove_unreferenced on inactive submesh
        std::cerr << "v_inactive rows: " << v_inactive.rows() << "\n";
        std::cerr << "f_inactive rows: " << f_inactive.rows() << "\n";

        // Sanity check v_active and v_inactive positions are subsets of Vfull
        // Check first and last vertex of v_active exist in Vfull (rough check)
        std::cerr << "v_active first vertex: " << v_active.row(0) << "\n";
        std::cerr << "Vfull row J(0):        " << Vfull.row(J(0)) << "\n";
        std::cerr << "sanity: v_active[0] == Vfull[J[0]]? "
                << (v_active.row(0).isApprox(Vfull.row(J(0))) ? "YES" : "NO") << "\n";

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // M restricted to active submesh via J (new->old map)
        // M was the global clipped mass diagonal; J maps local->global
        // So M_active(i) = M(J(i)) for each local vertex i
        Eigen::VectorXd M_active(v_active.rows());
        for (int i = 0; i < (int)v_active.rows(); ++i)
            M_active(i) = M(J(i));  // M is the global clipped mass vector from earlier
        // sanity check
        std::cerr << "M_active size: " << M_active.size() << " (should == v_active.rows() == " << v_active.rows() << ")\n";
        std::cerr << "M_active min: " << M_active.minCoeff() << " (should be >= 1e-8)\n";
        std::cerr << "M_active max: " << M_active.maxCoeff() << "\n";
        std::cerr << "sanity: M_active size == v_active rows? " << (M_active.size() == v_active.rows() ? "YES" : "NO") << "\n";

        // dblA: per-face area weights (double area / 2 = area)
        Eigen::VectorXd dblA_vec;
        igl::doublearea(V, F, dblA_vec);
        dblA_vec /= 2.0;

        // ~~~~~~~~~~~~~~~~~~~~~ Sanity checks for dblA ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        std::cerr << "dblA_vec size: " << dblA_vec.size() << " (should == F.rows() == " << F.rows() << ")\n";
        std::cerr << "dblA_vec min: " << dblA_vec.minCoeff() << " (should be > 0, degenerate faces if not)\n";
        std::cerr << "dblA_vec max: " << dblA_vec.maxCoeff() << "\n";
        std::cerr << "sanity: any degenerate faces? " << (dblA_vec.minCoeff() <= 0 ? "YES - WARNING" : "NO") << "\n";
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // v_xyz: stack x, y, z coordinates into one long vector [x | y | z]
        int nV_active = (int)V.rows();
        Eigen::VectorXd v_xyz(nV_active * 3);
        v_xyz.segment(0,          nV_active) = V.col(0);
        v_xyz.segment(nV_active,  nV_active) = V.col(1);
        v_xyz.segment(2*nV_active, nV_active) = V.col(2);

        // ~~~~~~~~~~~~~~~~~~ Sanity checks for v_xyz ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        std::cerr << "v_xyz size: " << v_xyz.size() << " (should == 3 * V.rows() == " << 3 * V.rows() << ")\n";
        std::cerr << "sanity: v_xyz size == 3 * nV_active? " << (v_xyz.size() == 3 * nV_active ? "YES" : "NO") << "\n";
        // Check first vertex matches
        std::cerr << "V(0,:) = " << V.row(0) << "\n";
        std::cerr << "v_xyz[0], v_xyz[nV], v_xyz[2nV] = " 
                << v_xyz(0) << ", " << v_xyz(nV_active) << ", " << v_xyz(2*nV_active) << "\n";
        std::cerr << "sanity: v_xyz stacking correct? " 
                << (std::abs(v_xyz(0) - V(0,0)) < 1e-10 && 
                    std::abs(v_xyz(nV_active) - V(0,1)) < 1e-10 && 
                    std::abs(v_xyz(2*nV_active) - V(0,2)) < 1e-10 ? "YES" : "NO") << "\n";
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // boundary vertices of active submesh
        Eigen::VectorXi boundary_verts = cf::get_all_boundary_vids(F);
        // sanity checks:
        std::cerr << "boundary_verts size: " << boundary_verts.size() << "\n";
        std::cerr << "sanity: boundary verts < nV_active? " 
                << (boundary_verts.size() == 0 || boundary_verts.maxCoeff() < nV_active ? "YES" : "NO") << "\n";

        // fixed = unique(concat(fixed_test, boundary_verts))
        std::vector<int> fixed_vec;
        fixed_vec.reserve(fixed_test.size() + boundary_verts.size());
        for (int i = 0; i < fixed_test.size(); ++i)    fixed_vec.push_back(fixed_test(i));
        for (int i = 0; i < boundary_verts.size(); ++i) fixed_vec.push_back(boundary_verts(i));
        std::sort(fixed_vec.begin(), fixed_vec.end());
        fixed_vec.erase(std::unique(fixed_vec.begin(), fixed_vec.end()), fixed_vec.end());
        Eigen::VectorXi fixed = Eigen::Map<Eigen::VectorXi>(fixed_vec.data(), (int)fixed_vec.size());
        // sanity checks:
        std::cerr << "fixed_test size: " << fixed_test.size() << "\n";
        std::cerr << "boundary_verts size: " << boundary_verts.size() << "\n";
        std::cerr << "fixed size (union): " << fixed.size() << "\n";
        std::cerr << "sanity: fixed <= fixed_test + boundary? " 
                << (fixed.size() <= (size_t)(fixed_test.size() + boundary_verts.size()) ? "YES" : "NO") << "\n";
        std::cerr << "sanity: fixed indices in range? " 
                << (fixed.size() == 0 || (fixed.minCoeff() >= 0 && fixed.maxCoeff() < nV_active) ? "YES" : "NO") << "\n";

        // fixed_xyz: extend fixed into stacked xyz space
        // [fixed | fixed + nV | fixed + 2*nV]
        int nFixed = (int)fixed.size();
        Eigen::VectorXi fixed_xyz(nFixed * 3);
        fixed_xyz.segment(0,       nFixed) = fixed;
        fixed_xyz.segment(nFixed,  nFixed) = fixed.array() + nV_active;
        fixed_xyz.segment(2*nFixed, nFixed) = fixed.array() + 2 * nV_active;
        // sanity checks: 
        std::cerr << "fixed_xyz size: " << fixed_xyz.size() << " (should == 3 * fixed.size() == " << 3 * fixed.size() << ")\n";
        std::cerr << "sanity: fixed_xyz max < 3*nV_active? " 
                << (fixed_xyz.maxCoeff() < 3 * nV_active ? "YES" : "NO") << "\n";

        // Principal curvature directions per face
        // Principal curvature directions and values per face
        Eigen::MatrixXd PD1, PD2;
        Eigen::VectorXd PC1, PC2;
        pfpc::per_face_prin_curvature(V, F, PD1, PD2, PC1, PC2);
        if (params.opening)
            PD1 = PD2;
        // Sanity checks
        std::cerr << "PD1 rows: " << PD1.rows() << " cols: " << PD1.cols() << " (should be " << F.rows() << " x 3)\n";
        std::cerr << "PD2 rows: " << PD2.rows() << " cols: " << PD2.cols() << " (should be " << F.rows() << " x 3)\n";
        std::cerr << "sanity: PD1 rows == F rows? " << (PD1.rows() == F.rows() ? "YES" : "NO") << "\n";
        // Check PD1 rows are unit vectors (principal directions should be normalized)
        double pd1_norm_first = PD1.row(0).norm();
        std::cerr << "PD1 row 0 norm: " << pd1_norm_first << " (should be ~1.0)\n";
        std::cerr << "opening mode: " << (params.opening ? "YES (using PD2)" : "NO (using PD1)") << "\n";




        // proy_matrix = D matrix, per-face dot product with curvature directions
        // proy_matrix = hstack([diag(PD1[:,0]), diag(PD1[:,1]), diag(PD1[:,2])])
        // Shape: (nF, 3*nF) — projects gradient onto PD1 direction
        int nF = F.rows();
        Eigen::SparseMatrix<double> proy_matrix(nF, 3 * nF); 
        std::vector<Eigen::Triplet<double>> proy_trips;
        proy_trips.reserve(3 * nF);
        for (int i = 0; i < nF; i++) {
            proy_trips.emplace_back(i, i,        PD1(i, 0));  // x block
            proy_trips.emplace_back(i, i + nF,   PD1(i, 1));  // y block
            proy_trips.emplace_back(i, i + 2*nF, PD1(i, 2));  // z block
        }
        proy_matrix.setFromTriplets(proy_trips.begin(), proy_trips.end());
        // ─── proy_matrix: D ─────────────────────────────────────────────────────────
        std::cerr << "proy_matrix shape: (" << proy_matrix.rows() << ", " << proy_matrix.cols() << ")"
                << " (should be nF x 3*nF = " << nF << " x " << 3*nF << ")\n";
        std::cerr << "sanity proy_matrix rows == nF?    " << (proy_matrix.rows() == nF ? "YES" : "NO") << "\n";
        std::cerr << "sanity proy_matrix cols == 3*nF?  " << (proy_matrix.cols() == 3*nF ? "YES" : "NO") << "\n";
        // PD1 rows should be unit vectors
        double pd1_min_norm = std::numeric_limits<double>::max();
        double pd1_max_norm = 0.0;
        for (int i = 0; i < nF; i++) {
            double n = PD1.row(i).norm();
            pd1_min_norm = std::min(pd1_min_norm, n);
            pd1_max_norm = std::max(pd1_max_norm, n);
        }
        std::cerr << "PD1 row norm range: [" << pd1_min_norm << ", " << pd1_max_norm << "]"
                << " (should be ~1.0 for all rows)\n";

        // G matrix, gradient operator: (3*nF x nV)
        Eigen::SparseMatrix<double> G;
        igl::grad(V, F, G);
        // ─── G: gradient operator ───────────────────────────────────────────────────
        std::cerr << "G shape: (" << G.rows() << ", " << G.cols() << ")"
                << " (should be 3*nF x nV = " << 3*nF << " x " << nV_active << ")\n";
        std::cerr << "sanity G rows == 3*nF? " << (G.rows() == 3*nF ? "YES" : "NO") << "\n";
        std::cerr << "sanity G cols == nV?   " << (G.cols() == nV_active ? "YES" : "NO") << "\n";

        // D*G
        // projected_gradient = proy_matrix @ G: (nF x nV)
        Eigen::SparseMatrix<double> projected_gradient = proy_matrix * G;
        // ─── projected_gradient: D*G ────────────────────────────────────────────────
        std::cerr << "projected_gradient shape: (" << projected_gradient.rows() << ", " << projected_gradient.cols() << ")"
                << " (should be nF x nV = " << nF << " x " << nV_active << ")\n";
        std::cerr << "sanity projected_gradient rows == nF?  " << (projected_gradient.rows() == nF ? "YES" : "NO") << "\n";
        std::cerr << "sanity projected_gradient cols == nV?  " << (projected_gradient.cols() == nV_active ? "YES" : "NO") << "\n";

        // A matrix, diagonal matrix of triangle areas
        // dblA as sparse diagonal matrix
        Eigen::SparseMatrix<double> dblA_sparse(nF, nF);
        std::vector<Eigen::Triplet<double>> area_trips;
        area_trips.reserve(nF);
        for (int i = 0; i < nF; i++)
            area_trips.emplace_back(i, i, dblA_vec(i));
        dblA_sparse.setFromTriplets(area_trips.begin(), area_trips.end());
        // ─── dblA_sparse: A ─────────────────────────────────────────────────────────
        std::cerr << "dblA_sparse shape: (" << dblA_sparse.rows() << ", " << dblA_sparse.cols() << ")"
                << " (should be nF x nF = " << nF << " x " << nF << ")\n";
        std::cerr << "dblA min: " << dblA_vec.minCoeff() << " (should be > 0, else degenerate faces)\n";
        std::cerr << "sanity no degenerate faces? " << (dblA_vec.minCoeff() > 0 ? "YES" : "NO - WARNING") << "\n";


        // G^T * D^T * A * D * G — the Dirichlet energy term (int_proy_grad_sq)
        // int_proy_grad_sq = -Gp^T * dblA * Gp: (nV x nV)
        Eigen::SparseMatrix<double> int_proy_grad_sq =
            -(projected_gradient.transpose() * dblA_sparse * projected_gradient);
        // ─── int_proy_grad_sq: G^T D^T A D G ────────────────────────────────────────
        std::cerr << "int_proy_grad_sq shape: (" << int_proy_grad_sq.rows() << ", " << int_proy_grad_sq.cols() << ")"
                << " (should be nV x nV = " << nV_active << " x " << nV_active << ")\n";
        std::cerr << "sanity int_proy_grad_sq square? "
                << (int_proy_grad_sq.rows() == int_proy_grad_sq.cols() ? "YES" : "NO") << "\n";
        std::cerr << "sanity int_proy_grad_sq rows == nV? "
                << (int_proy_grad_sq.rows() == nV_active ? "YES" : "NO") << "\n";

        // M matrix, mass matrix block diagonal
        // M as sparse diagonal: block_diag([M, M, M]) shape (3*nV, 3*nV)
        Eigen::SparseMatrix<double> M_diag(nV_active, nV_active);
        std::vector<Eigen::Triplet<double>> m_trips;
        m_trips.reserve(nV_active);
        for (int i = 0; i < nV_active; i++)
            m_trips.emplace_back(i, i, M_active(i));
        M_diag.setFromTriplets(m_trips.begin(), m_trips.end());

        // block_diag([M, M, M]): (3*nV, 3*nV)
        Eigen::SparseMatrix<double> M_block(3*nV_active, 3*nV_active);
        std::vector<Eigen::Triplet<double>> mb_trips;
        mb_trips.reserve(3 * nV_active);
        for (int i = 0; i < nV_active; i++) {
            mb_trips.emplace_back(i,              i,              M_active(i));
            mb_trips.emplace_back(i + nV_active,  i + nV_active,  M_active(i));
            mb_trips.emplace_back(i + 2*nV_active, i + 2*nV_active, M_active(i));
        }
        M_block.setFromTriplets(mb_trips.begin(), mb_trips.end());
        // ─── M_block ────────────────────────────────────────────────────────────────
        std::cerr << "M_block shape: (" << M_block.rows() << ", " << M_block.cols() << ")"
                << " (should be 3*nV x 3*nV = " << 3*nV_active << " x " << 3*nV_active << ")\n";
        std::cerr << "sanity M_block square?      " << (M_block.rows() == M_block.cols() ? "YES" : "NO") << "\n";
        std::cerr << "sanity M_block rows == 3*nV? " << (M_block.rows() == 3*nV_active ? "YES" : "NO") << "\n";

        // block_diag([int_proy_grad_sq x3]): (3*nV, 3*nV)
        Eigen::SparseMatrix<double> S_block(3*nV_active, 3*nV_active);
        std::vector<Eigen::Triplet<double>> s_trips;
        for (int k = 0; k < int_proy_grad_sq.outerSize(); k++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(int_proy_grad_sq, k); it; ++it) {
                s_trips.emplace_back(it.row(),              it.col(),              it.value());
                s_trips.emplace_back(it.row() + nV_active,  it.col() + nV_active,  it.value());
                s_trips.emplace_back(it.row() + 2*nV_active, it.col() + 2*nV_active, it.value());
            }
        }
        S_block.setFromTriplets(s_trips.begin(), s_trips.end());
        // ─── S_block ────────────────────────────────────────────────────────────────
        std::cerr << "S_block shape: (" << S_block.rows() << ", " << S_block.cols() << ")"
                << " (should be 3*nV x 3*nV = " << 3*nV_active << " x " << 3*nV_active << ")\n";
        std::cerr << "sanity S_block square?       " << (S_block.rows() == S_block.cols() ? "YES" : "NO") << "\n";
        std::cerr << "sanity S_block rows == 3*nV? " << (S_block.rows() == 3*nV_active ? "YES" : "NO") << "\n";

        // K = M - τ * G^T * D^T * A * D * G                <--- This is Q
        // Q = block_diag(M,M,M) - 0.01 * dt * block_diag(S,S,S)
        Eigen::SparseMatrix<double> Q = M_block - 0.01 * params.dt * S_block;
        // ─── Q = K ──────────────────────────────────────────────────────────────────
        std::cerr << "Q shape: (" << Q.rows() << ", " << Q.cols() << ")"
                << " (should be 3*nV x 3*nV = " << 3*nV_active << " x " << 3*nV_active << ")\n";
        std::cerr << "sanity Q square?       " << (Q.rows() == Q.cols() ? "YES" : "NO") << "\n";
        std::cerr << "sanity Q rows == 3*nV? " << (Q.rows() == 3*nV_active ? "YES" : "NO") << "\n";

        // M * V^t      <---- this is linear
        // linear = -block_diag(M,M,M) @ v_xyz
        Eigen::VectorXd linear = -(M_block * v_xyz);
        // ─── linear = -M*v_xyz ──────────────────────────────────────────────────────
        std::cerr << "linear size: " << linear.size()
                << " (should be 3*nV = " << 3*nV_active << ")\n";
        std::cerr << "sanity linear size == 3*nV? " << (linear.size() == 3*nV_active ? "YES" : "NO") << "\n";
        std::cerr << "v_xyz size: " << v_xyz.size()
                << " (should be 3*nV = " << 3*nV_active << ")\n";

        // Aeq, Beq: empty equality constraints
        Eigen::SparseMatrix<double> Aeq(0, 3 * nV_active);
        Eigen::VectorXd Beq(0);

        // Fixed values: current positions of fixed vertices in stacked xyz
        // Eigen::VectorXd fixed_vals = v_xyz(fixed_xyz);  // requires igl::slice or manual indexing
        Eigen::VectorXd fixed_vals;
        igl::slice(v_xyz, fixed_xyz, fixed_vals);
        // ─── fixed_xyz ──────────────────────────────────────────────────────────────
        std::cerr << "fixed_xyz size: " << fixed_xyz.size()
                << " (should be 3 * fixed.size() = " << 3*fixed.size() << ")\n";
        std::cerr << "sanity fixed_xyz max < 3*nV? "
                << (fixed_xyz.maxCoeff() < 3*nV_active ? "YES" : "NO") << "\n";
        std::cerr << "sanity fixed_xyz min >= 0?   "
                << (fixed_xyz.minCoeff() >= 0 ? "YES" : "NO") << "\n";
        // Check that the three blocks are correctly offset
        bool fixed_xyz_correct = true;
        for (int i = 0; i < fixed.size(); i++) {
            if (fixed_xyz(i) != fixed(i) ||
                fixed_xyz(i + fixed.size()) != fixed(i) + nV_active ||
                fixed_xyz(i + 2*fixed.size()) != fixed(i) + 2*nV_active) {
                fixed_xyz_correct = false;
                break;
            }
        }
        std::cerr << "sanity fixed_xyz offsets correct? " << (fixed_xyz_correct ? "YES" : "NO") << "\n";

        // ─── fixed_vals ─────────────────────────────────────────────────────────────
        std::cerr << "fixed_vals size: " << fixed_vals.size()
                << " (should == fixed_xyz.size() = " << fixed_xyz.size() << ")\n";
        std::cerr << "sanity fixed_vals size == fixed_xyz size? "
                << (fixed_vals.size() == fixed_xyz.size() ? "YES" : "NO") << "\n";
        // Spot check: fixed_vals(0) should equal v_xyz(fixed_xyz(0))
        if (fixed_vals.size() > 0)
            std::cerr << "sanity fixed_vals[0] == v_xyz[fixed_xyz[0]]? "
                    << (std::abs(fixed_vals(0) - v_xyz(fixed_xyz(0))) < 1e-10 ? "YES" : "NO") << "\n";

        // V^(t+1) solution!! (this is vector uu, which we reshape as nx3 matrix U)
        // min_quad_with_fixed solve
        Eigen::VectorXd uu;
        igl::min_quad_with_fixed_data<double> mqwf_data;
        igl::min_quad_with_fixed_precompute(Q, fixed_xyz, Aeq, true, mqwf_data);
        igl::min_quad_with_fixed_solve(mqwf_data, linear, fixed_vals, Beq, uu);

        // Reshape uu from [x|y|z] stacked vector back to (nV x 3) matrix
        // uu = uu.reshape(3, -1).T  in Python
        Eigen::MatrixXd U(nV_active, 3);
        U.col(0) = uu.segment(0,           nV_active);
        U.col(1) = uu.segment(nV_active,   nV_active);
        U.col(2) = uu.segment(2*nV_active, nV_active);


        // ─── after solve: uu and U ──────────────────────────────────────────────────
        std::cerr << "uu size: " << uu.size()
                << " (should be 3*nV = " << 3*nV_active << ")\n";
        std::cerr << "sanity uu size == 3*nV? " << (uu.size() == 3*nV_active ? "YES" : "NO") << "\n";
        std::cerr << "U shape: (" << U.rows() << ", " << U.cols() << ")"
                << " (should be nV x 3 = " << nV_active << " x 3)\n";
        std::cerr << "sanity U rows == nV? " << (U.rows() == nV_active ? "YES" : "NO") << "\n";
        std::cerr << "sanity U cols == 3?  " << (U.cols() == 3 ? "YES" : "NO") << "\n";

        // Fixed vertices should not have moved
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

        // U should not have exploded (rough sanity on bounding box)
        Eigen::Vector3d U_min = U.colwise().minCoeff();
        Eigen::Vector3d U_max = U.colwise().maxCoeff();
        Eigen::Vector3d V_min = v_active.colwise().minCoeff();
        Eigen::Vector3d V_max = v_active.colwise().maxCoeff();
        std::cerr << "v_active bbox: [" << V_min.transpose() << "] to [" << V_max.transpose() << "]\n";
        std::cerr << "U bbox:        [" << U_min.transpose() << "] to [" << U_max.transpose() << "]\n";
        double bbox_expansion = (U_max - U_min).norm() / std::max((V_max - V_min).norm(), 1e-8);
        std::cerr << "bbox expansion ratio: " << bbox_expansion
                << " (should be close to 1.0, large values suggest solver diverged)\n";
        

        int nV_remesh = (int)U.rows();
        Eigen::VectorXd target_vec = Eigen::VectorXd::Constant(nV_remesh, params.h);
        const auto t_remesh_start = std::chrono::steady_clock::now();
        remesh_botsch(U, F, target_vec, params.remesh_iterations, boundary_verts, true);
        remesh_seconds_total +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_remesh_start)
            .count();
        std::cerr << "Finished remeshing active submesh\n";
        std::cerr << "U rows after remesh: " << U.rows() << "\n";
        std::cerr << "F rows after remesh: " << F.rows() << "\n";

        // Compute interior_verts BEFORE merge
        // boundary_verts = get_all_boundary_vids(F)
        // interior_verts = setdiff1d(arange(len(U)), boundary_verts)
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
        std::cerr << "boundary verts after remesh: " << boundary_verts_remeshed.size() << "\n";
        std::cerr << "interior verts after remesh: " << interior_verts.size() << "\n";

        // Udup = U  (Python does a shallow copy here, no perturbation applied)
        Eigen::MatrixXd Udup = U;

        // merge active and inactive surface
        // Vfull = concatenate([v_inactive, Udup], axis=0)
        // Ffull = concatenate([f_inactive, F + len(v_inactive)], axis=0)
        int n_inactive_verts = (int)v_inactive.rows();
        int n_inactive_faces = (int)f_inactive.rows();
        int n_active_verts   = (int)Udup.rows();
        int n_active_faces   = (int)F.rows();

        Eigen::MatrixXd Vfull_new(n_inactive_verts + n_active_verts, 3);
        Vfull_new.topRows(n_inactive_verts)    = v_inactive;
        Vfull_new.bottomRows(n_active_verts)   = Udup;

        Eigen::MatrixXi Ffull_new(n_inactive_faces + n_active_faces, 3);
        Ffull_new.topRows(n_inactive_faces)    = f_inactive;
        Ffull_new.bottomRows(n_active_faces)   = F.array() + n_inactive_verts;

        std::cerr << "Vfull_new rows before dedup: " << Vfull_new.rows() << "\n";
        std::cerr << "Ffull_new rows before dedup: " << Ffull_new.rows() << "\n";

        // Vfull, SVI, SVJ, Ffull = igl.remove_duplicate_vertices(Vfull, Ffull, 0)
        Eigen::VectorXi SVI, SVJ;
        igl::remove_duplicate_vertices(Vfull_new, Ffull_new, 0, Vfull, SVI, SVJ, Ffull);

        std::cerr << "Vfull rows after dedup: " << Vfull.rows() << "\n";
        std::cerr << "Ffull rows after dedup: " << Ffull.rows() << "\n";
        std::cerr << "SVI size: " << SVI.size() << " (should == Vfull.rows())\n";
        std::cerr << "SVJ size: " << SVJ.size() << " (should == Vfull_new.rows())\n";
        std::cerr << "sanity SVI size == Vfull rows? " 
                  << (SVI.size() == Vfull.rows() ? "YES" : "NO") << "\n";
        std::cerr << "sanity SVJ size == Vfull_new rows? " 
                  << (SVJ.size() == Vfull_new.rows() ? "YES" : "NO") << "\n";


        // ── Step 4: Map interior verts into full mesh ────────────────────────────
        // Python: interior_active_full = SVJ[interior_verts + len(v_inactive)]
        // SVJ is now valid since we haven't remeshed yet after dedup
        Eigen::VectorXi interior_active_full(interior_verts.size());
        for (int i = 0; i < interior_verts.size(); i++)
            interior_active_full(i) = SVJ(interior_verts(i) + n_inactive_verts);

        std::cerr << "interior_active_full size: " << interior_active_full.size() << "\n";
        std::cerr << "sanity interior_active_full max < Vfull rows? "
                  << (interior_active_full.size() == 0 ||
                      interior_active_full.maxCoeff() < Vfull.rows() ? "YES" : "NO - WARNING") << "\n";

        // ── Step 5: Rebuild A from active region edges ───────────────────────────
        // Python: E = SVJ[edges(F + len(v_inactive))]
        //         A = sparse matrix from E + identity
        int nVfull_post = (int)Vfull.rows();
        Eigen::MatrixXi F_global = F.array() + n_inactive_verts;
        Eigen::MatrixXi edges_active;
        igl::edges(F_global, edges_active);

        A.resize(nVfull_post, nVfull_post);
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
        A.setFromTriplets(a_trips.begin(), a_trips.end());

        std::cerr << "A rebuilt (active edges via SVJ), size: " << A.rows() << " x " << A.cols() << "\n";
        std::cerr << "sanity: A.rows() == Vfull.rows()? "
                  << (A.rows() == Vfull.rows() ? "YES" : "NO - WARNING") << "\n";


        // Step 6
        // Python: H_interior_active = discrete_mean_curvature(U, F)
        //         K_interior_active = discrete_gaussian_curvature(U, F)
        Eigen::VectorXd H_interior_active = cf::discrete_mean_curvature(U, F);
        Eigen::VectorXd K_interior_active = cf::discrete_gaussian_curvature(U, F);


        // Step 7
        // Python: Mnew = sp.diags(np.ones(len(Vfull)))
        //         M_active = massmatrix(U, F)
        //         m_active = diagonal of M_active
        //         Mnew[interior_active_full, interior_active_full] += m_active[interior_verts] - 1
        //         M = Mnew
        Eigen::SparseMatrix<double> M_active_sparse;
        igl::massmatrix(U, F, igl::MASSMATRIX_TYPE_VORONOI, M_active_sparse);
        Eigen::VectorXd m_active(U.rows());
        for (int i = 0; i < (int)U.rows(); i++)
            m_active(i) = M_active_sparse.coeff(i, i);

        M.resize(Vfull.rows());
        M.setOnes();
        for (int i = 0; i < interior_verts.size(); i++) {
            int full_idx  = interior_active_full(i);
            int local_idx = interior_verts(i);
            M(full_idx) = M(full_idx) + m_active(local_idx) - 1.0;
        }
        for (int i = 0; i < M.size(); i++)
            M(i) = std::max(M(i), 1e-8);

        std::cerr << "M updated, size: " << M.size()
                  << " (should == Vfull.rows() == " << Vfull.rows() << ")\n";


        //  Step 8
        // Python: k_interior_active = H + sqrt(H^2 - M_active*K + 0j)
        //         K_interior_active = real(k[:,1] / safe_diag(M_active))  (closing)
        Eigen::VectorXd m_active_clipped(U.rows());
        for (int i = 0; i < (int)U.rows(); i++)
            m_active_clipped(i) = std::max(m_active(i), 1e-8);

        Eigen::VectorXd disc_active = H_interior_active.array().square()
                                    - m_active_clipped.array() * K_interior_active.array();
        Eigen::VectorXcd sqrtDisc_active = disc_active.cast<std::complex<double>>().array().sqrt();
        Eigen::VectorXcd k_max_active = H_interior_active.cast<std::complex<double>>() + sqrtDisc_active;
        Eigen::VectorXcd k_min_active = H_interior_active.cast<std::complex<double>>() - sqrtDisc_active;

        Eigen::VectorXd K_interior_active_curv;
        if (params.opening) {
            K_interior_active_curv = (-k_max_active.array() /
                m_active_clipped.cast<std::complex<double>>().array()).real().matrix();
        } else {
            K_interior_active_curv = (k_min_active.array() /
                m_active_clipped.cast<std::complex<double>>().array()).real().matrix();
        }

        // Step 9
        // Python: moving = zeros(len(Vfull))
        //         moving_interior_active = is_active(K_interior_active)
        //         moving[interior_active_full] = moving_interior_active[interior_verts]
        //         moving = where(moving)[0]
        std::vector<bool> moving_full(Vfull.rows(), false);
        for (int i = 0; i < interior_verts.size(); i++) {
            int local_idx = interior_verts(i);
            int full_idx  = interior_active_full(i);
            if (K_interior_active_curv(local_idx) < -params.bd)
                moving_full[full_idx] = true;
        }
        std::vector<int> moving_vec;
        moving_vec.reserve(Vfull.rows());
        for (int i = 0; i < (int)Vfull.rows(); i++)
            if (moving_full[i]) moving_vec.push_back(i);
        moving.resize(moving_vec.size());
        for (int i = 0; i < (int)moving_vec.size(); i++)
            moving(i) = moving_vec[i];

        std::cerr << "moving size after update: " << moving.size() << "\n";

        // step 10: convergence check
        // Python: if (iter + 1) % 10 == 0: check point_mesh_squared_distance
        if ((iter + 1) % 10 == 0) {
            Eigen::VectorXd sqr_dist;
            Eigen::VectorXi closest_faces;
            Eigen::MatrixXd closest_points;
            igl::point_mesh_squared_distance(Vfull, Vprev, Fprev, sqr_dist, closest_faces, closest_points);
            double max_dist = sqr_dist.maxCoeff();
            std::cerr << "iter " << iter+1 << " max squared dist from prev: " << max_dist << "\n";
            if (max_dist < params.tol) {
                std::cerr << "Converged at iter " << iter+1 << "\n";
                break;
            }
        }

    }

    const double flow_seconds_total =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t_flow_start)
            .count();
    std::cerr << "closing_flow total time: " << flow_seconds_total << " s\n";
    std::cerr << "remeshing total time: " << remesh_seconds_total << " s\n";

    V_out = std::move(Vfull);
    F_out = std::move(Ffull);
    return true;
}