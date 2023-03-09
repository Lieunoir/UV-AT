#include "UV.h"

double getArea(Eigen::MatrixXd const &V, Eigen::MatrixXi const &F, int i) {
    Eigen::Vector3d V1 = V.row(F(i,1)) - V.row(F(i,0));
    Eigen::Vector3d V2 = V.row(F(i,2)) - V.row(F(i,0));
    return V1.cross(V2).norm() / 2.;
}

Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

Eigen::MatrixXd ProjectedNewtonSymDir::initiateHarmonicX(Eigen::VectorXi& I, Eigen::MatrixXi& F2, Eigen::MatrixXd& V2, std::vector<std::vector<int>>& cuts) {
    igl::cut_to_disk(F, cuts);
    if(cuts.size() == 0) {
        std::vector<int> simpleCut;
        simpleCut.push_back(F(43,0));
        simpleCut.push_back(F(43,1));
        simpleCut.push_back(F(43,2));
        cuts.push_back(simpleCut);
    }

    std::set<std::array<int, 2>> cut_edges;
    for (const auto& cut : cuts) {
        const size_t cut_len = cut.size();
        for (size_t i=0; i<cut_len-1; i++) {
            std::array<int, 2> e{cut[i], cut[i+1]};
            if (e[0] > e[1]) {
                std::swap(e[0], e[1]);
            }
            cut_edges.insert(e);
        }
    }

    const size_t num_faces = F.rows();
    Eigen::MatrixXi cut_mask(num_faces, 3);
    cut_mask.setZero();
    for (size_t i=0; i<num_faces; i++) {
        std::array<int, 2> e0{F(i, 0), F(i, 1)};
        std::array<int, 2> e1{F(i, 1), F(i, 2)};
        std::array<int, 2> e2{F(i, 2), F(i, 0)};
        if (e0[0] > e0[1]) std::swap(e0[0], e0[1]);
        if (e1[0] > e1[1]) std::swap(e1[0], e1[1]);
        if (e2[0] > e2[1]) std::swap(e2[0], e2[1]);

        if (cut_edges.find(e0) != cut_edges.end()) {
            cut_mask(i, 0) = 1;
        }
        if (cut_edges.find(e1) != cut_edges.end()) {
            cut_mask(i, 1) = 1;
        }
        if (cut_edges.find(e2) != cut_edges.end()) {
            cut_mask(i, 2) = 1;
        }
    }

    const Eigen::MatrixXd Vcopy = V;
    igl::cut_mesh(Vcopy, F, cut_mask, V2, F2, I);

    Eigen::VectorXi bnd;
    igl::boundary_loop(F2,bnd);
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V2,bnd,bnd_uv);

    Eigen::VectorXd M;
    igl::doublearea(V, F, M);
    bnd_uv *= sqrt(M.sum() / (2 * igl::PI));
    Eigen::MatrixXd uv_init;
    igl::harmonic(V2,F2,bnd,bnd_uv,1,uv_init);
    if(igl::flipped_triangles(uv_init, F2).size() != 0)
        igl::harmonic(F2,bnd,bnd_uv,1,uv_init);

    F = F2;
    V = V2;
    igl::boundary_loop(F, edges);

    X = uv_init;
    return uv_init;
}

Eigen::MatrixXd ProjectedNewtonSymDir::DGprojectedNewtonSymDir(bool singleIt) {
    int n = V.rows();

    double gradientMax = 10.;
    double precGradient = 0.;
    int iter = 0;
    if(bijective) {
        symDirEnergy(X);
    }

    do {
        iter ++;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(2*n);
        Eigen::SparseMatrix<double> H(2*n, 2*n);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(36 * F.rows());

        for(int i = 0; i < F.rows(); i++) {
            const Eigen::Matrix2d Fg = getDeformationGradient(X, i);
            const Eigen::Vector4d f = vectorify(Fg);

            const Eigen::JacobiSVD<Eigen::Matrix2d> svd(Fg, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const Eigen::Matrix2d svdU = svd.matrixU();
            const Eigen::Matrix2d svdV = svd.matrixV();
            const Eigen::Vector2d Sigma = svd.singularValues();

            const Eigen::Matrix<double, 4, 6> dfdx = getFGradient(i);
            const Eigen::Vector<double, 6> bq = (v(i) * v(i) + flatV) * getArea(V, F, i) * dfdx.transpose() * getSymDirEnergyGradient(svdU, svdV, Sigma, f, v(i));
            b(F(i,0))   += bq(0);
            b(F(i,1))   += bq(1);
            b(F(i,2))   += bq(2);
            b(n+F(i,0)) += bq(3);
            b(n+F(i,1)) += bq(4);
            b(n+F(i,2)) += bq(5);

            const std::array<std::pair<double, Eigen::Vector4d>, 4> eigenSystem = evalSymdirEigenSystem(svdU, Sigma, svdV, v(i));

            Eigen::Matrix4d Hq = Eigen::Matrix4d::Zero();
            for(auto eigenPair : eigenSystem) {
                Hq += std::max(eigenPair.first, 0.) * eigenPair.second * eigenPair.second.transpose();
            }

            const Eigen::Matrix<double, 6, 6> Htemp = (v(i) * v(i) + flatV) * getArea(V, F, i) * dfdx.transpose() * Hq * dfdx;

            for(int k = 0; k < 6; k++) {
                for(int l = 0; l < 6; l++) {
                    triplets.push_back({ F(i,k%3)+(k>=3?n:0), F(i,l%3)+(l>=3?n:0), Htemp(k,l)});
                }
            }
        }
        H.setFromTriplets(triplets.begin(), triplets.end());
        precGradient = gradientMax;
        gradientMax = b.lpNorm<Eigen::Infinity>();
        if(bijective) {
            b += getBorderEGradient() * borderAlpha;
            H += borderAlpha * getBorderEHessian();
        }

        if(use_bounding_box) {
            b += getBoundingEGradient() * bounding_weight;
            H += getBoundingEHessian() * bounding_weight;
        }

        if(!solver_inited || bijective || use_bounding_box) {
            solver_inited = true;
            solver.analyzePattern(H);
        }

        solver.factorize(H);
        Eigen::VectorXd d = solver.solve(b);
        Eigen::MatrixXd D = Eigen::MatrixXd(d.rows() / 2, 2);
        D.col(0) = d.head(d.rows() / 2);
        D.col(1) = d.tail(d.rows() / 2);
        D = X - D;

        std::function<double(Eigen::MatrixXd&)> energy(std::bind(&ProjectedNewtonSymDir::symDirEnergy, this, std::placeholders::_1));
        Eigen::MatrixXd Xprec = X;
        igl::flip_avoiding_line_search(F, X, D, energy);

        if(!use_bounding_box) { 
            Eigen::RowVector2d origin = X.row(0);
            for(int i = 0; i < X.rows(); i++) {
                X.row(i) -= origin;
            }
        }
    } while(gradientMax > 0.01 && abs(gradientMax - precGradient) > 0.00000000001 && !singleIt);

    return X;
}

void ProjectedNewtonSymDir::addBorderEGradientPart(int i, int j1, int j2, Eigen::MatrixXd const& Xv, Eigen::VectorXd &b) {
    int n = V.rows();
    Eigen::Vector2d Xi  = Xv.row(i);
    Eigen::Vector2d Xj1 = Xv.row(j1);
    Eigen::Vector2d Xj2 = Xv.row(j2);
    if(i != j1 && i != j2) {
        double nij1 = (Xj1-Xi).norm();
        double nij2 = (Xj2-Xi).norm();
        double nj1j2 = (Xj1-Xj2).norm();
        double d = 0.5 * ( nij1 + nij2 - nj1j2 );
        double fact = - 2. * borderE * (borderE - d) / (d*d*d);
        if(d < borderE) {
            b(i)    += 0.5 * fact * ((Xi(0)  - Xj1(0)) / nij1 + (Xi(0)  - Xj2(0)) / nij2);
            b(j1)   += 0.5 * fact * ((Xj1(0) - Xi(0))  / nij1 - (Xj1(0) - Xj2(0)) / nj1j2);
            b(j2)   += 0.5 * fact * ((Xj2(0) - Xi(0))  / nij2 - (Xj2(0) - Xj1(0)) / nj1j2);
            b(n+i)  += 0.5 * fact * ((Xi(1)  - Xj1(1)) / nij1 + (Xi(1)  - Xj2(1)) / nij2);
            b(n+j1) += 0.5 * fact * ((Xj1(1) - Xi(1))  / nij1 - (Xj1(1) - Xj2(1)) / nj1j2);
            b(n+j2) += 0.5 * fact * ((Xj2(1) - Xi(1))  / nij2 - (Xj2(1) - Xj1(1)) / nj1j2);
        }
    }
}

Eigen::VectorXd ProjectedNewtonSymDir::getBorderEGradient() {
    int n = V.rows();
    Eigen::VectorXd b = Eigen::VectorXd::Zero(2*n);
    Eigen::MatrixXd Xv = X;

    for(auto eLoop: edges) {
        for(int it1 = 0; it1 < eLoop.size(); it1++) {
            for(int it2 = 0; it2 < eLoop.size(); it2++) {
                if(it1 != 0 && it2 != 0) {
                    int i1 = eLoop[it1-1];
                    int i2 = eLoop[it1];
                    int j1 = eLoop[it2-1];
                    int j2 = eLoop[it2];
                    addBorderEGradientPart(i1, j1, j2, Xv, b);
                    addBorderEGradientPart(i2, j1, j2, Xv, b);
                }
            }
        }
    }

    return b;
}

void ProjectedNewtonSymDir::addBorderEHessianPart(int i, int j1, int j2, Eigen::MatrixXd const& Xv, std::vector<Eigen::Triplet<double>> &triplets) {
    int n = V.rows();
    Eigen::Vector2d Xi  = Xv.row(i);
    Eigen::Vector2d Xj1 = Xv.row(j1);
    Eigen::Vector2d Xj2 = Xv.row(j2);
    if(i != j1 && i != j2) {
        Eigen::VectorXd b = Eigen::VectorXd::Zero(6);
        Eigen::MatrixXd Htemp = Eigen::MatrixXd::Zero(6,6);
        double nij1 = (Xj1-Xi).norm();
        double nij2 = (Xj2-Xi).norm();
        double nj1j2 = (Xj1-Xj2).norm();
        double d = 0.5 * ( nij1 + nij2 - nj1j2 );
        if(d < borderE) {
            double f1 = - 2. * borderE * (borderE - d) / (d*d*d);
            double f2 = 2. * borderE * (d + 3*(borderE - d)) / (d*d*d*d);
            b(0) += 0.5 * ((Xi(0)  - Xj1(0)) / nij1 + (Xi(0)  - Xj2(0)) / nij2);
            b(1) += 0.5 * ((Xi(1)  - Xj1(1)) / nij1 + (Xi(1)  - Xj2(1)) / nij2);
            b(2) += 0.5 * ((Xj1(0) - Xi(0))  / nij1 - (Xj1(0) - Xj2(0)) / nj1j2);
            b(3) += 0.5 * ((Xj1(1) - Xi(1))  / nij1 - (Xj1(1) - Xj2(1)) / nj1j2);
            b(4) += 0.5 * ((Xj2(0) - Xi(0))  / nij2 - (Xj2(0) - Xj1(0)) / nj1j2);
            b(5) += 0.5 * ((Xj2(1) - Xi(1))  / nij2 - (Xj2(1) - Xj1(1)) / nj1j2);
            Eigen::MatrixXd Htemp2 = Eigen::MatrixXd::Zero(4,4);
            Eigen::VectorXd b2 = Eigen::VectorXd::Zero(4);
            b2(0) += Xj1(0) - Xj2(0);
            b2(1) += Xj1(1) - Xj2(1);
            b2(2) += Xj2(0) - Xj1(0);
            b2(3) += Xj2(1) - Xj1(1);
            Htemp2(0,0) = nj1j2;
            Htemp2(1,1) = nj1j2;
            Htemp2(2,2) = nj1j2;
            Htemp2(3,3) = nj1j2;
            Htemp2(0,2) = -nj1j2;
            Htemp2(1,3) = -nj1j2;
            Htemp2(2,0) = -nj1j2;
            Htemp2(3,1) = -nj1j2;
            Htemp2 = - (Htemp2 - b2 * b2.transpose() / nj1j2);
            Htemp += f2 * b * b.transpose();
            Htemp.block<4,4>(2,2) += f1 * Htemp2 / (nj1j2 * nj1j2);
            std::array<int,3> indices = {i, j1, j2};
            for(int k = 0; k < 3; k++) {
                for(int l = 0; l < 3; l++) {
                    triplets.push_back({ indices[k], indices[l], Htemp(2*k,2*l)});
                    triplets.push_back({ indices[k]+n, indices[l], Htemp(2*k+1,2*l)});
                    triplets.push_back({ indices[k], indices[l]+n, Htemp(2*k,2*l+1)});
                    triplets.push_back({ indices[k]+n, indices[l]+n, Htemp(2*k+1,2*l+1)});
                }
            }
        }
    }
}

Eigen::SparseMatrix<double> ProjectedNewtonSymDir::getBorderEHessian() {
    int n = V.rows();
    Eigen::SparseMatrix<double> H(2*n, 2*n);
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::MatrixXd Xv = X;

    for(auto eLoop: edges) {
        for(int it1 = 0; it1 < eLoop.size(); it1++) {
            for(int it2 = 0; it2 < eLoop.size(); it2++) {
                if(it1 != 0 && it2 != 0) {
                    int i1 = eLoop[it1-1];
                    int i2 = eLoop[it1];
                    int j1 = eLoop[it2-1];
                    int j2 = eLoop[it2];
                    addBorderEHessianPart(i1, j1, j2, Xv, triplets);
                    addBorderEHessianPart(i2, j1, j2, Xv, triplets);
                }
            }
        }
    }
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}

Eigen::SparseMatrix<double> ProjectedNewtonSymDir::getBoundingEHessian() {
    int n = V.rows();
    Eigen::SparseMatrix<double> H(2*n, 2*n);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < n; i++) {
        if(X(i,0) > bounding_x)   
            triplets.push_back({ i, i, 2.});
        if(X(i,0) < 0)   
            triplets.push_back({ i, i, 2.});
        if(X(i,1) > bounding_y)   
            triplets.push_back({ i+n, i+n, 2.});
        if(X(i,1) < 0)   
            triplets.push_back({ i+n, i+n, 2.});
    }
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}

Eigen::VectorXd ProjectedNewtonSymDir::getBoundingEGradient() {
    int n = V.rows();
    Eigen::VectorXd b = Eigen::VectorXd::Zero(2*n);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < n; i++) {
        if(X(i,0) > bounding_x)
            b(i) += 2*(X(i,0) - bounding_x);
        if(X(i,0) < 0)   
            b(i) += 2*(X(i,0));
        if(X(i,1) > bounding_y)   
            b(i+n) += 2*(X(i,1) - bounding_y);
        if(X(i,1) < 0)   
            b(i+n) += 2*(X(i,1));
    }
    return b;
}
