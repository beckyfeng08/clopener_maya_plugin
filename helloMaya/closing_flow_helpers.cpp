#include "closing_flow_helpers.h"

#include <cmath>

#include <igl/barycenter.h>
#include <igl/edge_flaps.h>
#include <igl/edge_lengths.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/boundary_loop.h>
#include <igl/internal_angles.h>

namespace closing_flow_detail {

    namespace {
        // MSVC does not define M_PI unless _USE_MATH_DEFINES is set before <cmath>;
        // keep a local constant so include order does not matter.
        constexpr double kPi = 3.14159265358979323846264338327950288;
    } // namespace

    // Computes dihedral angles per face per edge (n_faces x 3)
    Eigen::MatrixXd dihedral_angles(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        Eigen::MatrixXi& TT_out)
    {
        // Triangle-triangle adjacency
        Eigen::MatrixXi TT;
        igl::triangle_triangle_adjacency(F, TT);

        // Fix edge convention mismatch - triangle_triangle_adjacency.h
        Eigen::PermutationMatrix<3, 3> perm(3);
        perm.indices() = Eigen::Vector3i(1, 2, 0);
        TT = (TT * perm).eval();
        TT_out = TT;

        // Face normals
        Eigen::MatrixXd FN;
        // Eigen::MatrixXd Z = Eigen::MatrixXd::Constant(F.rows(), 3, 1.0/3.0);
        // igl::per_face_normals(V, F, Z, FN);
        igl::per_face_normals(V, F, FN);

        // Face barycenters
        Eigen::MatrixXd BC;
        igl::barycenter(V, F, BC);

        int nF = F.rows();
        Eigen::MatrixXd angles(nF, 3); // dihedral angle per face per edge

        for (int f = 0; f < nF; f++) { // loop over each face
            for (int e = 0; e < 3; e++) { // loop over each face edge (defined by opposite vertex)
                int adj = TT(f, e);
                if (adj < 0) {
                    angles(f, e) = 0.0;
                    continue;
                }

                // Dot product of the two face normals
                double dot = FN.row(f).dot(FN.row(adj));
                dot = std::max(-1.0, std::min(1.0, dot));

                // angle between normals
                double angle = std::acos(dot);

                // pi - angle, clipped to [0, pi]
                angle = std::max(0.0, kPi - angle);

                // Convexity check
                Eigen::RowVector3d diff = BC.row(adj) - BC.row(f);
                double conv_check = diff.dot(FN.row(f));
                if (conv_check >= 1e-8) {
                    // Non-convex: angle = 2*pi - angle
                    angle = 2.0 * kPi - angle;
                }

                angles(f, e) = angle;
            }
        }

        return angles;
    }

    // Computes discrete mean curvature per vertex
    // 0.5 * SUM over edge lengths times dihedral angle ??
    Eigen::VectorXd discrete_mean_curvature(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F)
    {
        Eigen::MatrixXi TT;
        Eigen::MatrixXd A = dihedral_angles(V, F, TT);  // (F, 3))

        // Edge lengths, reordered to match [2,0,1] column permutation
        Eigen::MatrixXd L_raw;
        igl::edge_lengths(V, F, L_raw);  // (F, 3), columns [01, 12, 20]
        Eigen::MatrixXd L(F.rows(), 3);
        L.col(0) = L_raw.col(2);  // opposite v0 = edge between v1,v2 = col 2 of raw
        L.col(1) = L_raw.col(0);  // opposite v1 = edge between v0,v2 = col 0 of raw  
        L.col(2) = L_raw.col(1);  // opposite v2 = edge between v0,v1 = col 1 of raw

        // Edge flaps for vertex indexing
        Eigen::MatrixXi E, EF, EI;
        Eigen::VectorXi EMAP;
        igl::edge_flaps(F, E, EMAP, EF, EI);

        // EMAP reshape: column major, reorder columns [2,0,1]
        // EMAP is F.rows()*3 long, stored column-major as (F, 3)
        int nF = F.rows();
        Eigen::MatrixXi emap_mat(nF, 3);
        for (int j = 0; j < 3; j++)
            for (int i = 0; i < nF; i++)
                emap_mat(i, j) = EMAP(j * nF + i);

        // Permute columns [2,0,1]
        Eigen::MatrixXi emap(nF, 3);
        emap.col(0) = emap_mat.col(2);
        emap.col(1) = emap_mat.col(0);
        emap.col(2) = emap_mat.col(1);

        // Adjacency matrix is TT

        // Fix boundary: set A to pi where adj < 0 so (pi - A) = 0
        for (int f = 0; f < nF; f++)
            for (int e = 0; e < 3; e++)
                if (TT(f, e) < 0)
                    A(f, e) = kPi;

        // cur = 0.5 * 0.5 * 0.5 * (pi - A) * L, flattened
        // vector of val
        Eigen::VectorXd cur(nF * 3);
        for (int f = 0; f < nF; f++)
            for (int e = 0; e < 3; e++)
                cur(f * 3 + e) = 0.125 * (kPi - A(f, e)) * L(f, e);

        // Accumulate onto vertices via edge endpoints
        Eigen::VectorXd H = Eigen::VectorXd::Zero(V.rows());
        for (int f = 0; f < nF; f++) {
            for (int e = 0; e < 3; e++) {
                int edge_idx = emap(f, e);
                double val = cur(f * 3 + e);
                H(E(edge_idx, 0)) += val;
                H(E(edge_idx, 1)) += val;
            }
        }

        return H;
    }

    Eigen::VectorXi get_all_boundary_vids(const Eigen::MatrixXi& F)
    {
        std::vector<std::vector<int>> loops;
        igl::boundary_loop(F, loops);

        int total = 0;
        for (const auto& lp : loops) total += static_cast<int>(lp.size());

        Eigen::VectorXi bnd(total);
        int idx = 0;
        for (const auto& lp : loops)
            for (int v : lp)
                bnd(idx++) = v;
        return bnd;
    }

    Eigen::VectorXd discrete_gaussian_curvature(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F)
    {
        // Accumulate internal angles per vertex
        Eigen::MatrixXd angles;
        igl::internal_angles(V, F, angles); // (nF, 3)

        Eigen::VectorXd K = Eigen::VectorXd::Constant(V.rows(), 2.0 * kPi);
        for (int f = 0; f < F.rows(); f++)
            for (int e = 0; e < 3; e++)
                K(F(f, e)) -= angles(f, e);

        // Boundary correction
        Eigen::VectorXi bnd = get_all_boundary_vids(F);
        for (int i = 0; i < bnd.size(); i++)
            K(bnd(i)) -= kPi;

        return K;
    }

    Eigen::VectorXi incident_faces(
        const Eigen::MatrixXi& F,
        const Eigen::VectorXi& vids,    // vertex IDs, not face IDs!
        int v_incidence)
    {
        const int nV = F.maxCoeff() + 1;
        const int nF = F.rows();

        // selected[v] = true if v is in vids
        std::vector<bool> selected(nV, false);
        for (int i = 0; i < vids.size(); i++)
            selected[vids(i)] = true;

        // Count how many vertices of each face are selected
        std::vector<int> result;
        for (int f = 0; f < nF; f++)
        {
            int count = 0;
            for (int e = 0; e < F.cols(); e++)
                if (selected[F(f, e)])
                    count++;
            if (count >= v_incidence)
                result.push_back(f);
        }

        // Convert to Eigen vector
        Eigen::VectorXi out(result.size());
        for (int i = 0; i < (int)result.size(); i++)
            out(i) = result[i];
        return out;
    }

    std::vector<bool> expandRing(const Eigen::SparseMatrix<int>& A, const std::vector<bool>& inSet)
    {
        int nVerts = (int)inSet.size();
        std::vector<bool> result(nVerts, false);
        for (int j = 0; j < nVerts; ++j) {
            if (!inSet[j]) continue;
            for (Eigen::SparseMatrix<int>::InnerIterator it(A, j); it; ++it)
                result[it.row()] = true;
        }
        return result;
    }


} // namespace closing_flow_detail