#include "AT.h"

// Mostly logic for cutting

struct V_D {
    int index;
    int prec_index;
    double distance;
    bool border_crossed;
};


bool operator>(const V_D& lhs, const V_D& rhs)
{
    return  lhs.distance > rhs.distance;
}

void FaceUVAT::cut(Eigen::MatrixXi const& cut_mask, Eigen::MatrixXd const&_V, bool retract) {
    Eigen::MatrixXi F2;
    Eigen::VectorXi I2;
    Eigen::VectorXi Itemp;
    const Eigen::MatrixXd Xcopy = X;
    igl::cut_mesh(Xcopy, F, cut_mask, X, F2, I2);
    Itemp = I2;
    for(int i = 0; i < I.rows(); i++) {
        Itemp(i) = I(Itemp(i));
    }
    I = Itemp;
    Eigen::MatrixXd V2 = Eigen::MatrixXd(I.rows(), 3);
    for(int i = 0; i < I.rows(); i++) {
        V2.row(i) = _V.row(I(i));
    }
    V = V2;
    F = F2;

    if(retract)
        doRetract(I2);

    UVSolver = ProjectedNewtonSymDir(V, F);
    UVSolver.setBijective(bijective);
    UVSolver.setBoundingBox(use_bounding_box, bounding_x, bounding_y, bounding_weight);

    buildFaceLaplacian();

    UVSolver.setV(v);
    UVSolver.setX(X);
}

std::vector<std::vector<int>> FaceUVAT::splitInFacePatches(Eigen::MatrixXi &TT, Eigen::MatrixXi &TTi, double threshold) {
    int n = F.rows();
    std::vector<bool> explored(n, false);
    // Contains a partition of all the faces
    std::vector<std::vector<int>> patches;

    std::vector<int> patch = std::vector<int>();
    bool skip = true;
    std::queue<int> neighborQueue;
    for(int i = 0; skip && i < n; i++) {
        if(v(i) < threshold) {
            skip = false;
            neighborQueue.push(i);
        }
    }
    while(!skip) {
        while(!neighborQueue.empty()) {
            int i = neighborQueue.front();
            neighborQueue.pop();
            if(!explored[i]) {
                explored[i] = true;
                patch.push_back(i);
                for(int j = 0; j < 3; j++) {
                    if(TT(i,j) != -1 && v(TT(i,j)) < 0.5) {
                        neighborQueue.push(TT(i,j));
                    }
                }
            }
        }
        patches.push_back(patch);
        skip = true;
        for(int i = 0; skip && i < n; i++) {
            if(v(i) < 0.5 && explored[i] == false) {
                patch = std::vector<int>();
                neighborQueue.push(i);
                skip = false;
            }
        }
    }
    return patches;
}

//Use local quantities for patch instead of global vertices index
//Gives the correct Fpatch and vector for local_to_global translation
void FaceUVAT::patchToLocal(std::vector<int> &patch, Eigen::MatrixXi &Fpatch, std::vector<int> &vertices) {
    Fpatch = Eigen::MatrixXi(patch.size(), 3);
    for(int i = 0; i < patch.size(); i++)
        Fpatch.row(i) = F.row(patch[i]);

    vertices = std::vector<int>();
    for(int i = 0; i < patch.size(); i++) {
        for(int j = 0; j < 3; j++) {
            int vertex = Fpatch(i,j);
            bool found = false;
            for(int k = 0; !found && (k < vertices.size()); k++) {
                if(vertices[k] == vertex)
                    found = true;
            }
            if(!found)
                vertices.push_back(vertex);
        }
    }
    std::sort(vertices.begin(), vertices.end());
    for(int j = 0; j < patch.size(); j++) {
        for(int k = 0; k < 3; k++) {
            bool found = false;
            for(int i = 0; !found && (i < vertices.size()); i++) {
                if(vertices[i] == Fpatch(j,k)) {
                    Fpatch(j,k) = i;
                    found = true;
                }
            }
        }
    }
}

double FaceUVAT::getFarthestPair(Eigen::MatrixXi &Fpatch, Eigen::MatrixXd &Vpatch, int &i_longest, int &j_longest, std::vector<bool> border_vertices, std::vector<int> vertices) {
    int nv_patch = Vpatch.rows();
    int nf_patch = Fpatch.rows();
    std::vector<std::vector<int>> adj_list = std::vector<std::vector<int>>(nv_patch, std::vector<int>());
    for(int i = 0; i < nf_patch; i++) {
        for(int j = 0; j < 3; j++) {
            int i_v1 = Fpatch(i,j);
            for(int k = 0; k < 3; k++) {
                int i_v2 = Fpatch(i,k);
                if(k != j) {
                    bool found = false;
                    for(int l = 0; (l < adj_list[i_v1].size()) && !found; l++) {
                        if(adj_list[i_v1][l] == i_v2)
                            found = true;
                    }
                    if (!found) {
                        adj_list[i_v1].push_back(i_v2);
                    }
                }
            }
        }
    }

    //Find shortest path
    std::vector<bool> explored = std::vector<bool>(nv_patch, false);
    std::vector<int> prec_v = std::vector<int>(nv_patch, -1);
    explored[i_longest] = true;
    std::priority_queue<V_D, std::vector<V_D>, std::greater<V_D>> toExplore;
    for(auto v: adj_list[i_longest]) {
        struct V_D v_d;
        v_d.index = v;
        v_d.prec_index = i_longest;
        v_d.distance = 0;
        if (!border_vertices[vertices[i_longest]] || !border_vertices[vertices[v]])
            v_d.distance = (Vpatch.row(v)-Vpatch.row(i_longest)).norm();
        toExplore.push(v_d);
    }

    if(toExplore.empty())
        return 0;

    while(!toExplore.empty()) {
        auto v = toExplore.top();
        toExplore.pop();
        if(!explored[v.index]) {
            prec_v[v.index] = v.prec_index;
            explored[v.index] = true;
            j_longest = v.index;
            for(auto v_2: adj_list[v.index]) {
                struct V_D v_d;
                v_d.index = v_2;
                v_d.prec_index = v.index;
                v_d.distance = v.distance;
                if (!border_vertices[vertices[v.index]] || !border_vertices[vertices[v_2]])
                    v_d.distance = v_d.distance + (Vpatch.row(v.index)-Vpatch.row(v_2)).norm();
                toExplore.push(v_d);
            }
        }
    }

    std::vector<int> path;
    int cur_v = j_longest;
    double length = 0;
    while(cur_v != i_longest) {
        length += (Vpatch.row(cur_v) - Vpatch.row(prec_v[cur_v])).norm();
        cur_v = prec_v[cur_v];
    }
    return length;
}

std::vector<std::vector<bool>> build_border_edges(Eigen::MatrixXi &F, Eigen::MatrixXd &V) {
    std::vector<std::vector<bool>> border_edges = std::vector<std::vector<bool>>(V.rows(), std::vector<bool>(V.rows(), false));
    for(int i = 0; i < F.rows(); i++) {
        for(int j = 0; j < 3; j++) {
            int i1 = F(i, j);
            int i2 = F(i, (j+1)%3);
            border_edges[i1][i2] = !border_edges[i1][i2];
            border_edges[i2][i1] = !border_edges[i2][i1];
        }
    }

    return border_edges;
}

Eigen::MatrixXi FaceUVAT::ATCutAtEdges(Eigen::MatrixXd const&_V, bool retract, bool tree, bool force_border_crossing) {
    //ATUV on distorted patchs with constrained borders
    Eigen::MatrixXi TT;
    Eigen::MatrixXi TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);
    std::vector<std::vector<int>> patches = splitInFacePatches(TT, TTi, 0.5);

    // Build normals
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(V.rows(), 3);
    for(int i = 0; i < F.rows(); i++) {
        double area = getArea(V, F, i);
        Eigen::Vector3d v1 = V.row(F(i,1)) - V.row(F(i,0));
        Eigen::Vector3d v2 = V.row(F(i,2)) - V.row(F(i,0));
        Eigen::Vector3d normal = v1.cross(v2);
        for(int j = 0; j < 3; j++) {
            N.row(F(i,j)) += normal;
        }
    }
    for(int i = 0; i < V.rows(); i++) {
        V.row(i) = V.row(i) / V.row(i).norm();
    }

    std::vector<std::vector<int>> cuts_i;
    std::vector<bool> border_vertices = igl::is_border_vertex(F);
    std::vector<std::vector<bool>> border_edges = build_border_edges(F, V);
    for(auto patch : patches) {
        Eigen::MatrixXi Fpatch;
        std::vector<int> vertices;
        patchToLocal(patch, Fpatch, vertices);
        //Build F, V and X for each patch
        Eigen::MatrixXd Vpatch = Eigen::MatrixXd(vertices.size(), 3);
        Eigen::MatrixXd Xpatch = Eigen::MatrixXd(vertices.size(), 3);
        Eigen::MatrixXd Npatch = Eigen::MatrixXd(vertices.size(), 3);
        for(int i = 0; i < vertices.size(); i++) {
            Vpatch.row(i) = V.row(vertices[i]);
            Xpatch.row(i) = X.row(vertices[i]);
            Npatch.row(i) = N.row(vertices[i]);
        }


        bool border_patch = false;
        if(tree || force_border_crossing) {
            for(int i = 0; i < vertices.size(); i++) {
                if(border_vertices[vertices[i]])
                    border_patch = true;
            }
        }
        if(!border_patch && tree)
            continue;
        bool use_border_crossing = border_patch && force_border_crossing;

        // Compute cuts for real now
        // Compute farthest pair
        int nf_patch = Fpatch.rows();
        int nv_patch = vertices.size();
        double longest_d = 0.;
        int i_longest, j_longest;
        for(int i = 0; i < nv_patch; i++) {
            int j;
            double d = getFarthestPair(Fpatch, Vpatch, i, j, border_vertices, vertices);
            if (d > longest_d) {
                longest_d = d;
                i_longest = i;
                j_longest = j;
            }
        }

        if(longest_d <= 0.)
            continue;

        //Compute adj list
        std::vector<std::vector<int>> adj_list = std::vector<std::vector<int>>(nv_patch, std::vector<int>());
        for(int i = 0; i < nf_patch; i++) {
            for(int j = 0; j < 3; j++) {
                int i_v1 = Fpatch(i,j);
                for(int k = 0; k < 3; k++) {
                    int i_v2 = Fpatch(i,k);
                    if(k != j) {
                        bool found = false;
                        for(int l = 0; (l < adj_list[i_v1].size()) && !found; l++) {
                            if(adj_list[i_v1][l] == i_v2)
                                found = true;
                        }
                        if (!found) {
                            adj_list[i_v1].push_back(i_v2);
                        }
                    }
                }
            }
        }

        std::vector<bool> explored = std::vector<bool>(nv_patch, false);
        std::vector<bool> explored_second = std::vector<bool>(nv_patch, false);
        std::vector<int> prec_v = std::vector<int>(nv_patch, -1);
        std::vector<int> prec_v_second = std::vector<int>(nv_patch, -1);
        explored[i_longest] = true;
        std::priority_queue<V_D, std::vector<V_D>, std::greater<V_D>> toExplore;
        for(auto v: adj_list[i_longest]) {
            struct V_D v_d;
            v_d.index = v;
            v_d.prec_index = i_longest;
            v_d.distance = 0;
            v_d.border_crossed = false;
            if (!border_edges[vertices[i_longest]][vertices[v]])
                v_d.distance = (Vpatch.row(v)-Vpatch.row(i_longest)).norm();
            toExplore.push(v_d);
        }
        if(border_vertices[vertices[j_longest]])
            use_border_crossing = false;

        while((!use_border_crossing && !explored[j_longest]) || (use_border_crossing && !explored_second[j_longest])) {
            auto v = toExplore.top();
            toExplore.pop();
            if((!v.border_crossed && !explored[v.index]) || (v.border_crossed && !explored_second[v.index])) {
                if(!v.border_crossed) {
                    prec_v[v.index] = v.prec_index;
                    explored[v.index] = true;
                }
                else {
                    prec_v_second[v.index] = v.prec_index;
                    explored_second[v.index] = true;
                    //explored_second[v.prec_index] = true;
                }
                for(auto v_2: adj_list[v.index]) {
                    struct V_D v_d;
                    v_d.index = v_2;
                    v_d.prec_index = v.index;
                    v_d.distance = v.distance;
                    if (!border_edges[vertices[v.index]][vertices[v_2]])
                        v_d.distance = v_d.distance + (Vpatch.row(v.index)-Vpatch.row(v_2)).norm();
                    v_d.border_crossed = v.border_crossed || (border_vertices[vertices[v.index]] && !border_vertices[vertices[v_2]]);
                    toExplore.push(v_d);
                }
            }
        }

        std::vector<int> cut;
        int cur_v = j_longest;
        cut.push_back(j_longest);
        if(use_border_crossing) {
            while(!border_vertices[vertices[cur_v]]) {
                cut.push_back(prec_v_second[cur_v]);
                cur_v = prec_v_second[cur_v];
            }
        }
        while(cur_v != i_longest) {
            cut.push_back(prec_v[cur_v]);
            cur_v = prec_v[cur_v];
        }

        //Translate path to global cut
        for(int i = 0; i < cut.size(); i++) {
            cut[i] = vertices[cut[i]];
        }

        cuts_i.push_back(cut);
    }
    Eigen::MatrixXi cut_mask(F.rows(), 3);
    cut_mask.setZero();
    
    int cut_i = 0;
    double max_cut_len = 0.;
    for(int i = 0; i < cuts_i.size(); i++) {
        double cut_len = 0.;
        auto cut = cuts_i[i];
        for(int vi = 0; vi < cut.size()-1; vi++) {
            int v1 = cut[vi];
            int v2 = cut[vi+1];
            cut_len += (V.row(v1) - V.row(v2)).norm();
        }
        if(cut_len > max_cut_len) {
            max_cut_len = cut_len;
            cut_i = i;
        }
    }
    for(auto cut: cuts_i) {
        double cut_len = 0.;
        for(int vi = 0; vi < cut.size()-1; vi++) {
            int v1 = cut[vi];
            int v2 = cut[vi+1];
            cut_len += (V.row(v1) - V.row(v2)).norm();
        }
        if(cut_len > 0.2 * max_cut_len) {
            for(int vi = 0; vi < cut.size()-1; vi++) {
                int v1 = cut[vi];
                int v2 = cut[vi+1];
                if(!border_edges[v1][v2] && !border_edges[v2][v1]) {
                    bool found = false;
                    for(int i = 0; (i < F.rows()) && !found; i++) {
                        for(int j = 0; (j < 3) && !found; j++) {
                            if((F(i,j) == v1 && F(i,(j+1)%3) == v2) || (F(i,j) == v2 && F(i,(j+1)%3) == v1)) {
                                found = true;
                                cut_mask(i, j) = 1;
                                cut_mask(TT(i,j), TTi(i,j)) = 1;
                            }
                        }
                    }
                }
            }
        }
    }
    cut(cut_mask, _V, retract);

    return cut_mask;
}

// Used to seperate edges when using a global bijectivity energy
void FaceUVAT::doRetract(Eigen::MatrixXi const&I2) {
    int originN = 0;
    int newN = I2.rows();
    //We have to check which vertices have been duplicated
    for(int i = 0; i < newN; i++) {
        if(I2(i) > originN)
            originN = I2(i);
    }
    originN += 1;
    std::vector<std::vector<int>> vertexMapping(originN, std::vector<int>());
    for(int i = 0; i < newN; i++) {
        vertexMapping[I2(i)].push_back(i);
    }
    for(int i = 0; i < originN; i++) {
        //This is a duplicated vertex
        if(vertexMapping[i].size() > 1) {
            for(int l = 0; l < vertexMapping[i].size(); l++) {
                int vertexIndex = vertexMapping[i][l];
                //Used to determine which vertices are adjacent on the boundary (=1 vs =2 or =0)
                std::vector<int> vertexCount(newN, 0);
                Eigen::Vector2d direction = {0,0};
                //Used to determine which direction is the interior
                double angle = 0;
                for(int j = 0; j < F.rows(); j++) {
                    for(int k = 0; k < 3; k++) {
                        if(F(j,k) == vertexIndex) {
                            Eigen::Vector2d e1 = X.row(F(j,(k+1)%3)) - X.row(F(j,k));
                            Eigen::Vector2d e2 = X.row(F(j,(k+2)%3)) - X.row(F(j,k));
                            vertexCount[F(j,(k+1)%3)]++;
                            vertexCount[F(j,(k+2)%3)]++;
                            double angleDiff = std::acos(e1.dot(e2) / (e1.norm() * e2.norm()));
                            angle += angleDiff;
                        }
                    }
                }
                bool foundFirst = false;
                Eigen::Vector2d e1;
                Eigen::Vector2d e2;
                for(int j = 0; j < newN; j++) {
                    if(vertexCount[j] == 1) {
                        if(!foundFirst) {
                            e1 = X.row(j) - X.row(vertexIndex);
                            foundFirst = true;
                        } else
                            e2 = X.row(j) - X.row(vertexIndex);
                    }
                }
                direction = 0.5 * (e1 / e1.norm() + e2 / e2.norm());
                if(angle > 3.14159265359)
                    direction = -direction;
                //We limit how far we go in the direction using triangles flips (solve for det = 0)
                for(int j = 0; j < F.rows(); j++) {
                    for(int k = 0; k < 3; k++) {
                        if(F(j,k) == vertexIndex) {
                            Eigen::Vector2d u1 = X.row(F(j,(k+1)%3));
                            Eigen::Vector2d u2 = X.row(F(j,k));
                            Eigen::Vector2d u3 = X.row(F(j,(k+2)%3));
                            double a = direction(0) * (u3(1) - u1(1)) - direction(1) * (u3(0) - u1(0));
                            double b = u1(0) * (u2(1) - u3(1)) + u2(0) * (u3(1) - u1(1)) + u3(0) * (u1(1) - u2(1));
                            if(abs(a) > 0.00000000001) {
                                double t = - b / a;
                                if (t > 0 && t < 1)
                                    direction = 0.9 * t * direction;
                            }
                        }
                    }
                }
                X.row(vertexIndex) += direction;
            }
        }
    }
}
