#pragma once

#include "SurfaceMeshDEC.h"
#include "UV.h"
#include <igl/cut_mesh.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/avg_edge_length.h>
#include <igl/ambient_occlusion.h>
#include <igl/facet_components.h>
#include <igl/triangle/scaf.h>
#include <igl/MappingEnergyType.h>
#include <igl/is_border_vertex.h>

struct FaceUVAT {
    std::vector<std::vector<int>> build(double e, double l, Eigen::MatrixXd const &_V, Eigen::MatrixXi const &_F, bool _bijective, bool use_ao, bool use_bb, double bb_x, double bb_y, double bb_coeff)
    {
        V = _V;
        F = _F;
        bijective = _bijective;
        use_bounding_box = use_bb;
        bounding_x = bb_x;
        bounding_y = bb_y;
        bounding_weight = bb_coeff;
        avgl = igl::avg_edge_length(V, F);
        lambda_inv =  1. / ( l * avgl ) * Eigen::VectorXd::Ones( F.rows() );
        cut_memory = Eigen::MatrixXi::Zero(F.rows(), 3);

        if(use_ao)
            useAO();

        v = Eigen::VectorXd::Ones( F.rows() );
        buildMassMatrix();
        buildFaceLaplacian();
        buildFaceIdentity();
        setEpsilon( e );

        UVSolver = ProjectedNewtonSymDir(V, F);
        UVSolver.setBijective(bijective);
        UVSolver.setBoundingBox(use_bounding_box, bounding_x, bounding_y, bounding_weight);

        std::vector<std::vector<int>> cuts;
        X = UVSolver.initiateHarmonicX(I, F, V, cuts);
        return cuts;
    }

    void setEpsilon( double e )
    {
        epsilon         = e * avgl;
        l_1_over_4e     = (1./4./epsilon) * massMatrix * Eigen::VectorXd::Ones( F.rows() );
        l_1_over_4e_Id1 = (1./4./epsilon) * massMatrix * faceIdentity;
        le_ad1_d1       = epsilon * faceLaplacian;
    }

    void useAO() {
        Eigen::VectorXd S;
        Eigen::MatrixXd P = Eigen::MatrixXd(F.rows(), 3);
        Eigen::MatrixXd N = Eigen::MatrixXd(F.rows(), 3);
        double lambda_init = 1. / lambda_inv(0);
        double lambda_inv_init = lambda_inv(0);
        for(int i = 0; i < F.rows(); i++) {
            Eigen::Vector3d v0 = V.row(F(i,0));
            Eigen::Vector3d v1 = V.row(F(i,1));
            Eigen::Vector3d v2 = V.row(F(i,2));
            P.row(i) = (v0+v1+v2) / 3;
            N.row(i) = (v1-v0).cross(v2-v0);
            N.row(i).normalize();
        }
        //igl::ambient_occlusion(V, F, P, N, 128, S);
        S = Eigen::VectorXd(F.rows());
        for(int i = 0; i < F.rows(); i++) {
            S(i) = 1.;
        }

        int truc[59] = {3546, 3527, 3524, 3526, 3567, 3554, 3552,
            3272, 3275, 3306, 3298, 3296, 3299, 3303, 3574, 3265, 3270, 3325, 3273,
            3114, 3106, 3119, 3116, 3078, 3073, 3072, 3074, 3130, 3087,
            2088, 2090, 2059, 2049, 
            2062, 2079, 2064, 2066, 2065,
            2736, 2738, 2746, 2703, 2700, 2702, 2751, 2698, 2696,
            2205, 2713, 2714, 2712, 2719, 2706,
            2228, 2230, 2727, 2723, 2720, 2722};
        for(int i = 0; i < 59; i++) {
            S(truc[i]) = 100.;
        }
        /*
        S(3549) = 100.;
        S(3548) = 100.;
        S(3289) = 100.;
        S(3551) = 100.;
        S(3538) = 100.;*/

        for(int i = 0; i < F.rows(); i++) {
            lambda_inv(i) = lambda_inv(i) * S(i);
        }
        polyscope::getSurfaceMesh("input mesh")->addFaceScalarQuantity( "AO", S );
        /*
        double lambda_inv_mean = lambda_inv.sum() / lambda_inv.rows();
        lambda_inv = lambda_inv / lambda_inv_mean * lambda_inv_init;
        
        double lambda_tot = 0.;
        for(int i = 0; i < lambda_inv.rows(); i++) {
            lambda_tot += 1./ lambda_inv(i);
        }
        double lambda_mean = lambda_tot / lambda_inv.rows();
        lambda_inv = lambda_inv * lambda_mean / lambda_init;
        std::cout << "LAMBDA" << lambda_mean << ", " << lambda_init << std::endl;
        */
    }

    void buildFaceLaplacian() {
        Eigen::MatrixXi TT;
        Eigen::MatrixXi TTi;
        igl::triangle_triangle_adjacency(F, TT, TTi);
        std::vector<Eigen::Triplet<double>> triplets;

        int n = F.rows();
        faceLaplacian = Eigen::SparseMatrix<double>(n, n);
        for(int i = 0; i < n; i++) {
            Eigen::Vector3d ci = V.row(F(i,0)) + V.row(F(i,1)) + V.row(F(i,2));
            ci = ci / 3;
            for(int j = 0; j < 3; j++) {
                if(TT(i,j) >= 0) {
                    Eigen::Vector3d cj = V.row(F(TT(i,j),0)) + V.row(F(TT(i,j),1)) + V.row(F(TT(i,j),2));
                    cj = cj / 3;
                    Eigen::Vector3d e1 = V.row(F(i,j));
                    Eigen::Vector3d e2 = V.row(F(i,(j+1)%3));
                    float coeff = (e1 - e2).norm() / (ci - cj).norm();
                    //float coeff = 1.;
                    triplets.push_back({i, TT(i,j), -coeff});
                    triplets.push_back({i, i, coeff});
                }
            }
        }
        faceLaplacian.setFromTriplets(triplets.begin(), triplets.end());
        faceLaplacian = -faceLaplacian;
    }

    void buildMassMatrix() {
        int n = F.rows();
        massMatrix = Eigen::SparseMatrix<double>(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        for(int i = 0; i < n; i++) {
            Eigen::Vector3d e1 = V.row(F(i,1)) - V.row(F(i,0));
            Eigen::Vector3d e2 = V.row(F(i,2)) - V.row(F(i,0));
            triplets.push_back({i, i, e1.cross(e2).norm()/2});
        }
        massMatrix.setFromTriplets(triplets.begin(), triplets.end());
    }

    void buildFaceIdentity() {
        int n = F.rows();
        faceIdentity = Eigen::SparseMatrix<double>(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        for(int i = 0; i < n; i++) {
            triplets.push_back({i, i, 1});
        }
        faceIdentity.setFromTriplets(triplets.begin(), triplets.end());
    }

    std::vector<std::vector<int>> splitInFacePatches(Eigen::MatrixXi &TT, Eigen::MatrixXi &TTi, double threshold);

    //Use local quantities for patch instead of global vertices index
    //Gives the correct Fpatch and vector for local_to_global translation
    void patchToLocal(std::vector<int> &patch, Eigen::MatrixXi &Fpatch, std::vector<int> &vertices);

    double getFarthestPair(Eigen::MatrixXi &Fpatch, Eigen::MatrixXd &Vpatch, int &i_longest, int &j_longest, std::vector<bool> border_vertices, std::vector<int> vertices);

    Eigen::MatrixXi ATCutAtEdges(Eigen::MatrixXd const&_V, bool retract, bool tree, bool force_border_crossing);

    void cut(Eigen::MatrixXi const& cut_mask, Eigen::MatrixXd const&_V, bool retract);

    void doRetract(Eigen::MatrixXi const&I2);

    void solveOneAlternateStep()
    {
        DGtal::trace.beginBlock("Solving alternate UV AT");
        UVSolver.setX(X);

        DGtal::trace.beginBlock("Solving for u");
        
        X = UVSolver.DGprojectedNewtonSymDir();

        DGtal::trace.endBlock();
        DGtal::trace.beginBlock("Solving for v");
        former_v = v;

        int n = F.rows();
        Eigen::VectorXd guvSquared = Eigen::VectorXd::Zero(n);
        for(int i = 0; i < n; i++) {
            Eigen::Matrix2d Fg = UVSolver.getDeformationGradient(X, i);
            Eigen::Matrix2d Fginv = Fg.inverse();
            guvSquared(i) = 1. / 2. * (Fg.squaredNorm() + Fginv.squaredNorm() - 4.);
        }

        Eigen::SparseMatrix<double> ope_v1 = l_1_over_4e_Id1 - le_ad1_d1;
        ope_v1 += lambda_inv.asDiagonal() * massMatrix * guvSquared.asDiagonal();

        solver_v.compute( ope_v1 );
        DGtal::trace.info() << "Solving V v = l/4e * 1";
        v = solver_v.solve( l_1_over_4e );
        UVSolver.setV(v);


        DGtal::trace.info() << "  => "
            << ((solver_v.info() == Eigen::Success ) ? "OK" : "ERROR")
            << std::endl;
        DGtal::trace.info() << "Min v : " << v.minCoeff() << std::endl;
        printTotalEnergy();
        DGtal::trace.endBlock();
        DGtal::trace.endBlock();
    }

    void printTotalEnergy() {
        int n = F.rows();
        Eigen::VectorXd guvSquared = Eigen::VectorXd::Zero(n);
        double totArea = 0.;
        for(int i = 0; i < n; i++) {
            Eigen::Matrix2d Fg = UVSolver.getDeformationGradient(X, i);
            Eigen::Matrix2d Fginv = Fg.inverse();
            guvSquared(i) = 1. / 2. * (Fg.squaredNorm() + Fginv.squaredNorm() - 4.);
            totArea += getArea(V, F, i);
        }
        Eigen::SparseMatrix<double> ope_v1 = l_1_over_4e_Id1 - le_ad1_d1;
        ope_v1 += lambda_inv.asDiagonal() * massMatrix * guvSquared.asDiagonal();

    }

    double finalStep() {
        Eigen::VectorXd v1 = Eigen::VectorXd::Ones( F.rows() );
        UVSolver.setV(v1);
        UVSolver.setX(X);
        X = UVSolver.DGprojectedNewtonSymDir();

        double res = 0.;
        double totArea = 0.;
        for(int i = 0; i < F.rows(); i++) {
            Eigen::Matrix2d Fg = UVSolver.getDeformationGradient(X, i);
            Eigen::Matrix2d Fginv = Fg.inverse();
            res += getArea(V, F, i) * (Fg.squaredNorm() + Fginv.squaredNorm());
            totArea += getArea(V, F, i);
        }
        res = res / totArea;
        return res;
    }

    /// Computes the norms loo, l2, l1 of (v - former_v), i.e. the
    /// evolution of discontinuity function v.
    ///
    /// @return a tuple (n_infty,n_2,n_1) giving  the loo/l2/l1-norm of (v - former_v)
    std::tuple<double,double,double> diffV() const
    {
        Eigen::VectorXd delta = ( v - former_v ).cwiseAbs();
        const double n_oo = delta.maxCoeff();
        const double n_2  = std::sqrt( delta.squaredNorm() / delta.size() );
        const double n_1  = delta.mean();
        return std::make_tuple( n_oo, n_2, n_1 );
    }

    double avgl;
    Eigen::VectorXd lambda_inv;
    double epsilon;
    bool bijective;
    bool use_bounding_box;
    double bounding_x;
    double bounding_y;
    double bounding_weight;
    Eigen::MatrixXi F;
    Eigen::MatrixXd V;
    Eigen::VectorXd v;
    Eigen::MatrixXd X;
    //Map from original to duplicated vertices
    //Used to find duplicated vertices
    Eigen::VectorXi I;
    Eigen::MatrixXi cut_memory;
    Eigen::VectorXd former_v;
    Eigen::SparseMatrix<double> massMatrix;
    Eigen::SparseMatrix<double> faceLaplacian;
    Eigen::SparseMatrix<double> faceIdentity;
    Eigen::SparseMatrix<double> l_1_over_4e_Id1;
    Eigen::VectorXd l_1_over_4e;
    Eigen::SparseMatrix<double> le_ad1_d1;
    ProjectedNewtonSymDir UVSolver;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_v;
};
