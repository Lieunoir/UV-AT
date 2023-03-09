#pragma once

#include <Eigen/StdVector>
#include "SurfaceMeshDEC.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <igl/boundary_loop.h>
#include <igl/avg_edge_length.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/flip_avoiding_line_search.h>
#include <igl/cut_to_disk.h>
#include <igl/cut_mesh.h>
#include <igl/doublearea.h>
#include <igl/flipped_triangles.h>

extern double getArea(Eigen::MatrixXd const &V, Eigen::MatrixXi const &F, int i);

struct ProjectedNewtonSymDir
{
    std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> deformationHalf;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd X;
    std::vector<std::vector<int>> edges;
    Eigen::VectorXd v;
    double borderE = 1.;
    double borderAlpha = 0.0001;
    double flatV = 0.0000000000001;
    bool bijective = false;
    bool solver_inited = false;
    bool use_bounding_box = false;
    double bounding_x = 0.;
    double bounding_y = 0.;
    double bounding_weight = 1.;

    ProjectedNewtonSymDir() {}

    ProjectedNewtonSymDir(Eigen::MatrixXd const &V, Eigen::MatrixXi const &F) : V(V), F(F) {
        for(int i = 0; i < F.rows(); i++) {
            Eigen::VectorXd V1 = V.row(F(i,1)) - V.row(F(i,0));
            Eigen::VectorXd V2 = V.row(F(i,2)) - V.row(F(i,0));
            double x = V1.dot(V2) / V1.norm();
            double y = sqrt(V2.norm() * V2.norm() - x*x);
            Eigen::Matrix2d A;
            A << 1/V1.norm(), -x/(y*V1.norm()),
                 0.,           1/y;
            deformationHalf.push_back(A);
        }
        v = Eigen::VectorXd::Ones(F.rows());
        borderE = 10. * igl::avg_edge_length(V,F) / 8.;

        igl::boundary_loop(F,edges);
    }

    Eigen::MatrixXd initiateHarmonicX(Eigen::VectorXi &I) {
		Eigen::MatrixXd V2;
        Eigen::MatrixXi F2;
        std::vector<std::vector<int>> cuts;
        return initiateHarmonicX(I, F2, V2, cuts);
    }

    Eigen::MatrixXd initiateHarmonicX(Eigen::VectorXi& I, Eigen::MatrixXi& F2, Eigen::MatrixXd& V2, std::vector<std::vector<int>>& cuts);

    Eigen::Matrix2d getDeformationGradient(Eigen::MatrixXd const &values, int i) {
        Eigen::RowVector2d U1 = values.row(F(i,1)) - values.row(F(i,0));
        Eigen::RowVector2d U2 = values.row(F(i,2)) - values.row(F(i,0));
        Eigen::Matrix2d U;
        U << U1(0), U2(0),
             U1(1), U2(1);
        return U * deformationHalf[i];
    }

    Eigen::Matrix<double, 4, 6> getFGradient(int i) {
        Eigen::Matrix<double, 4, 6> dfdx;
        Eigen::Matrix2d H = deformationHalf[i];
        dfdx << -(H(0,0)+H(1,0)), H(0,0), H(1,0), 0., 0., 0.,
                0., 0., 0., -(H(0,0)+H(1,0)), H(0,0), H(1,0),
                -(H(0,1)+H(1,1)), H(0,1), H(1,1), 0., 0., 0.,
                0., 0., 0., -(H(0,1)+H(1,1)), H(0,1), H(1,1);
        return dfdx;
    }

    Eigen::Vector4d vectorify(Eigen::Matrix2d const &M) {
        int n = M.rows();
        int m = M.cols();
        Eigen::Vector4d res;
        for(int j = 0; j < m; j++) {
            for(int i = 0; i < n; i++) {
                res(i+j*n) = M(i,j);
            }
        }
        return res;
    }

    double symDirEnergy(Eigen::MatrixXd& values) {
        static bool found = false;
        double res = 0.;
        for(int i = 0; i < F.rows(); i++) {
            const Eigen::Matrix2d Fg = getDeformationGradient(values, i);
            const Eigen::Matrix2d Fginv = Fg.inverse();
            res += (v(i) * v(i) + flatV) * getArea(V, F, i) / 2. * (Fg.squaredNorm() + Fginv.squaredNorm() - 4.);
        }

        if(bijective) {
            Eigen::MatrixXd Xv = values;
            double Eb = 0.;
            for(auto eLoop: edges) {
                for(int it1 = 0; it1 < eLoop.size(); it1++) {
                    for(int it2 = 0; it2 < eLoop.size(); it2++) {
                        if(it1 != 0 && it2 != 0) {
                            int i1 = eLoop[it1-1];
                            int i2 = eLoop[it1];
                            int j1 = eLoop[it2-1];
                            int j2 = eLoop[it2];
                            if(i1 != j1 && i1 != j2) {
                                double d = 0.5 * ((Xv.row(j2)-Xv.row(i1)).norm() + (Xv.row(j1)-Xv.row(i1)).norm() - (Xv.row(j1)-Xv.row(j2)).norm());
                                if(d < borderE)
                                    Eb += (borderE / d - 1.) * (borderE / d - 1.);
                            }
                            if(i2 != j1 && i2 != j2) {
                                double d = 0.5 * ((Xv.row(j2)-Xv.row(i2)).norm() + (Xv.row(j1)-Xv.row(i2)).norm() - (Xv.row(j1)-Xv.row(j2)).norm());
                                if(d < borderE)
                                    Eb += (borderE / d - 1.) * (borderE / d - 1.);
                            }
                        }
                    }
                }
            }
            res = res + borderAlpha * Eb;
            std::cout << res << " (" << borderAlpha * Eb << ")" << std::endl;
        }

        if(use_bounding_box) {
            double bounding_energy = 0;
            for(int i = 0; i <  values.rows(); i++) {
                if(values(i,0) > bounding_x)
                    bounding_energy += (values(i,0) - bounding_x)*(values(i,0) - bounding_x);
                if(values(i,0) < 0)   
                    bounding_energy += values(i,0) * values(i,0);
                if(values(i,1) > bounding_y)   
                    bounding_energy += (values(i,1) - bounding_y)*(values(i,1) - bounding_y);
                if(values(i,1) < 0)   
                    bounding_energy += values(i,1) * values(i,1);
            }
            res += bounding_weight * bounding_energy;
        }
        return res;
    }

    Eigen::Vector4d getSymDirEnergyGradient(Eigen::Matrix2d const &svdU, Eigen::Matrix2d const &svdV, Eigen::Vector2d const &Sigma, Eigen::Vector4d const &f, double v) {
        double s1 = Sigma(0);
        double s2 = Sigma(1);
        double I2 = s1*s1 + s2*s2;
        double I3 = s1*s2;
        Eigen::Matrix2d Gtemp;
        Gtemp << Sigma(1), 0.,
                 0., Sigma(0);
        const Eigen::Matrix2d G = svdU * Gtemp * svdV.transpose();
        const Eigen::Vector4d g = vectorify(G);
        return 0.5 * ((1. + 1. / (I3*I3)) * 2. * f - 2. * I2 / (I3*I3*I3) * g);
    }

    Eigen::Matrix2d getTwist(Eigen::Matrix2d const &svdU, Eigen::Matrix2d const &svdV) {
        Eigen::Matrix2d temp;
        temp << 0., -1.,
                1.,  0.;
        return 1. / sqrt(2) *  svdU * temp * svdV.transpose();
    }

    Eigen::Matrix2d getFlip(Eigen::Matrix2d const &svdU, Eigen::Matrix2d const &svdV) {
        Eigen::Matrix2d temp;
        temp << 0., 1.,
                1., 0.;
        return 1. / sqrt(2) *  svdU * temp * svdV.transpose();
    }

    Eigen::Matrix2d getScaling(Eigen::Matrix2d const &svdU, Eigen::Matrix2d const &svdV, int i) {
        Eigen::Matrix2d temp;
        temp << (i == 1 ? 1. : 0.), 0.,
                0., (i == 2 ? 1. : 0.);
        return svdU * temp * svdV.transpose();
    }

    std::array<std::pair<double, Eigen::Vector4d>, 4> evalSymdirEigenSystem(Eigen::Matrix2d const &svdU, Eigen::Vector2d const &Sigma, Eigen::Matrix2d const &svdV, double v) {
        std::array<std::pair<double, Eigen::Vector4d>, 4> res;
        const Eigen::Matrix2d D1 = getScaling(svdU, svdV, 1);
        const Eigen::Vector4d d1 = vectorify(D1);
        const Eigen::Matrix2d D2 = getScaling(svdU, svdV, 2);
        const Eigen::Vector4d d2 = vectorify(D2);
        const Eigen::Matrix2d L = getFlip(svdU, svdV);
        const Eigen::Vector4d l = vectorify(L);
        const Eigen::Matrix2d T = getTwist(svdU, svdV);
        const Eigen::Vector4d t = vectorify(T);
        double s1 = Sigma(0);
        double s2 = Sigma(1);
        double I2 = s1*s1 + s2*s2;
        double I3 = s1*s2;
        res[0] = (std::pair<double, Eigen::Vector4d>(1. + 3. / (s1*s1*s1*s1), d1));
        res[1] = (std::pair<double, Eigen::Vector4d>(1. + 3. / (s2*s2*s2*s2), d2));
        res[2] = (std::pair<double, Eigen::Vector4d>(1. + 1. / (I3*I3) + I2 / (I3*I3*I3), l));
        res[3] = (std::pair<double, Eigen::Vector4d>(1. + 1. / (I3*I3) - I2 / (I3*I3*I3), t));
        return res;
    }

    Eigen::MatrixXd DGprojectedNewtonSymDir(bool singleIt = false);

    void addBorderEGradientPart(int i, int j1, int j2, Eigen::MatrixXd const& Xv, Eigen::VectorXd &b);

    Eigen::VectorXd getBorderEGradient();

    void addBorderEHessianPart(int i, int j1, int j2, Eigen::MatrixXd const& Xv, std::vector<Eigen::Triplet<double>> &triplets);

    Eigen::SparseMatrix<double> getBorderEHessian();

    Eigen::SparseMatrix<double> getBoundingEHessian();

    Eigen::VectorXd getBoundingEGradient();

    void setX(Eigen::MatrixXd const&Xt) {
        X = Xt;
    }

    void setBorderE(double borderEt) {
        borderE = borderEt;
    }

    void setV(Eigen::VectorXd const& _v) {
        v = _v;
    }

    Eigen::MatrixXd getCornerParam() {
        Eigen::MatrixXd Xc = Eigen::MatrixXd(3*F.rows(), 2);
        for(int i = 0; i < F.rows(); i++) {
            for(int j = 0; j < 3; j++) {
                Xc.row(3*i+j) = X.row(F(i,j));
            }
        }
        return Xc;
    }

    void setBijective(bool bij) {
        bijective = bij;
    }

    void setBoundingBox(bool use, double x, double y, double weight) {
        use_bounding_box = use;
        bounding_x = x;
        bounding_y = y;
        bounding_weight = weight;
    }
};

