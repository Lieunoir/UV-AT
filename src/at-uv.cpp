#include "SurfaceMeshBuilder.h"
#include <math.h>
#include <filesystem>
#include <igl/boundary_loop.h>
#include <igl/lscm.h>
#include <igl/upsample.h>

#include "DGtal/base/Common.h"
#include "DGtal/base/Clock.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/SurfaceMesh.h"
#include "DGtal/io/readers/SurfaceMeshReader.h"
#include "DGtal/io/writers/SurfaceMeshWriter.h"

#define WITH_EIGEN
#include "DGtal/math/linalg/EigenSupport.h"

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/utilities/eigen_interop_helpers.h"

#include "polyscope/polyscope.h"
#include "polyscope/pick.h"
#include "polyscope/surface_mesh.h"
#include "UV.h"
#include "AT.h"
#include "CLI11.hpp"

using namespace DGtal;
using namespace Z3i;

Eigen::MatrixXd V;
Eigen::MatrixXi F;
float Epsilon = 0.1;
double get_border_energy(Eigen::MatrixXi const&F_orig, Eigen::MatrixXd const&V, Eigen::MatrixXi const&F, int n);

typedef DGtal::SurfaceMesh< RealPoint, RealVector > SMesh;
typedef DGtal::SurfaceMeshBuilder< RealPoint, RealVector > SMeshBuilder;
typedef DGtal::SurfaceMeshDEC< EigenLinearAlgebraBackend,
        RealPoint, RealVector >     SMeshDEC;
typedef Eigen::Index Index;

SMesh          surfmesh;
SMeshDEC       dec;

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
typedef geometrycentral::surface::ManifoldSurfaceMesh GCSurfaceMesh;
// typedef geometrycentral::surface::SurfaceMesh GCSurfaceMesh;
std::unique_ptr<GCSurfaceMesh>    mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

polyscope::SurfaceMesh *psMesh;

//dec.init( &surfmesh );
SMeshDEC::Form curr_u2;
float Epsilon1 = 0.1;
float Epsilon2 = 0.01;
float Epsilonf = 0.5;
float Lambda   = 1.;
float Alpha    = 100.0;
int   NbIterations = 3;

//Reordering map to associate DEC SurfaceMesh corners id to
// GC/polyscope one.

void registerParam(Eigen::MatrixXd& X, std::string const paramName) {
    polyscope::getSurfaceMesh("input mesh")
        ->addVertexParameterizationQuantity(paramName, X);
    Eigen::MatrixXd Xproxy = Eigen::MatrixXd::Zero(X.rows(), 3);
    Xproxy.col(0) = X.col(0);
    Xproxy.col(1) = X.col(1);
    Xproxy = (Xproxy.rowwise() - Xproxy.colwise().mean());
    polyscope::registerSurfaceMesh("UV mesh", Xproxy, F)->setEdgeWidth(1.0)->setEdgeColor({0.,0.,0.});;
}

void registerCornerParam(Eigen::MatrixXd& X, std::string const paramName) {
    std::vector<std::array<double, 2>> init_uv_GC(3*F.rows());
    for(auto i=0; i < dec.nbCorners(); ++i)
        init_uv_GC[i] = {X(i,0), X(i,1)};
    polyscope::getSurfaceMesh("input mesh")
        ->addParameterizationQuantity(paramName, init_uv_GC);
    Eigen::MatrixXd Xproxy = Eigen::MatrixXd::Zero(X.rows(), 3);
    Xproxy.col(0) = X.col(0);
    Xproxy.col(1) = X.col(1);
    Xproxy = (Xproxy.rowwise() - Xproxy.colwise().mean());
    Eigen::MatrixXi Fp = Eigen::MatrixXi(F.rows(), 3);
    for(int i = 0; i < Fp.rows(); i++) {
        Fp(i, 0) = 3*i;
        Fp(i, 1) = 3*i+1;
        Fp(i, 2) = 3*i+2;
    }
    polyscope::registerSurfaceMesh("UV corner mesh", Xproxy, Fp)->setEdgeWidth(1.0)->setEdgeColor({0.,0.,0.});
}

void registerLoadedCut(Eigen::MatrixXi cut_mask) {
    std::vector<Vector3>               positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    for(int i = 0; i < cut_mask.rows(); i++) {
        for(int j = 0; j < 3; j++) {
            if(cut_mask(i,j) == 1) {
                const auto  cpt = positions.size();
                positions.push_back( geometry->vertexPositions[ F(i,j) ] );
                positions.push_back( geometry->vertexPositions[ F(i,(j+1)%3) ] );
                edgeInds.push_back( { cpt, cpt + 1 } );
            }
        }
    }
    psMesh->addSurfaceGraphQuantity("Loaded cut", positions, edgeInds);
}

void registerNewCut(Eigen::MatrixXi cut_mask) {
    std::vector<Vector3>               positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    for(int i = 0; i < cut_mask.rows(); i++) {
        for(int j = 0; j < 3; j++) {
            if(cut_mask(i,j) == 1) {
                const auto  cpt = positions.size();
                positions.push_back( geometry->vertexPositions[ F(i,j) ] );
                positions.push_back( geometry->vertexPositions[ F(i,(j+1)%3) ] );
                edgeInds.push_back( { cpt, cpt + 1 } );
            }
        }
    }
    psMesh->addSurfaceGraphQuantity("New cut", positions, edgeInds);
}

Eigen::MatrixXi cut_history;
void registerTotalCut(Eigen::MatrixXi cut_mask) {
    for(int i = 0; i < F.rows(); i++) {
        for(int j = 0; j < 3; j++) {
            if(cut_history(i,j) == 1 || cut_mask(i,j) == 1)
                cut_history(i,j) = 1;
        }
    }
    std::vector<Vector3>               positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    for(int i = 0; i < cut_mask.rows(); i++) {
        for(int j = 0; j < 3; j++) {
            if(cut_history(i,j) == 1) {
                const auto  cpt = positions.size();
                positions.push_back( geometry->vertexPositions[ F(i,j) ] );
                positions.push_back( geometry->vertexPositions[ F(i,(j+1)%3) ] );
                edgeInds.push_back( { cpt, cpt + 1 } );
            }
        }
    }
    psMesh->addSurfaceGraphQuantity("Total cut", positions, edgeInds);
}

void registerCurvesFromV1(const SMeshDEC::Form &bv1)
{
    std::vector<Vector3>               positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    for ( auto e = 0; e < bv1.size(); ++e )
        if ( bv1[ e ] != 1.0 )
        {
            const auto vtcs = surfmesh.edgeVertices( e );
            const auto    i = vtcs.first;
            const auto    j = vtcs.second;
            const auto  cpt = positions.size();
            positions.push_back( geometry->vertexPositions[ i ] );
            positions.push_back( geometry->vertexPositions[ j ] );
            edgeInds.push_back( { cpt, cpt + 1 } );
        }
    psMesh->addSurfaceGraphQuantity("Binarized V", positions, edgeInds);
}

// Track stuff to interwind with polyscope
bool faceATUVon = false;
bool faceATUVdone = false;
bool cutDone = false;
bool wholeATUVon = false;
bool wholeCATUVon = false;

FaceUVAT fUVAT;
// Method parameters
double faceATUVe;
int faceATUVi;
bool bijective = false;
bool use_ao = false;
bool prolongate_borders = false;
bool force_border_crossing = true;
bool use_bb = false;
float bb_x = 1.;
float bb_y = 1.;
float bb_coeff = 100.;

void doFacesATUV() {
    double loo, l2, l1;
    fUVAT.solveOneAlternateStep();
    std::tie( loo, l2, l1 ) = fUVAT.diffV();
    DGtal::trace.info() << "|| v^i - v^(i-1) ||_oo = " << loo << std::endl;

    Eigen::VectorXd distorsion = Eigen::VectorXd(F.rows());
    for(int i = 0; i < F.rows(); i++) {
        Eigen::Matrix2d Fg = fUVAT.UVSolver.getDeformationGradient(fUVAT.X, i);
        Eigen::Matrix2d Fginv = Fg.inverse();
        distorsion(i) = log( 1 + 1. / 2. * (Fg.squaredNorm() + Fginv.squaredNorm() - 4.));
    }

    polyscope::getSurfaceMesh("input mesh")->addFaceScalarQuantity( "Face v", fUVAT.v )->setEnabled(true);
    polyscope::getSurfaceMesh("input mesh")->addFaceScalarQuantity( "Distorsion", distorsion );
    Eigen::MatrixXd Xc = fUVAT.UVSolver.getCornerParam();
    registerCornerParam(Xc, "Param with AT");
    polyscope::getSurfaceMesh("UV corner mesh")->addFaceScalarQuantity( "Face v", fUVAT.v )->setEnabled(true);
    polyscope::getSurfaceMesh("UV corner mesh")->addFaceScalarQuantity( "Distorsion", distorsion );

    faceATUVi++;
    if(faceATUVi >= NbIterations) {
        faceATUVi = 0;
        faceATUVe *= Epsilonf;
        if(faceATUVe < Epsilon2) {
            faceATUVon = false;
            faceATUVdone = true;
            Eigen::MatrixXd Xc = fUVAT.UVSolver.getCornerParam();
            registerCornerParam(Xc, "Param with AT");
            polyscope::getSurfaceMesh("UV corner mesh")->addFaceScalarQuantity( "Face v", fUVAT.v )->setEnabled(true);
            polyscope::getSurfaceMesh("UV corner mesh")->addFaceScalarQuantity( "Distorsion", distorsion );
        } else {
            DGtal::trace.info() << "************************** Epsilon " << faceATUVe
                << " ***************************" << std::endl;
            fUVAT.setEpsilon( faceATUVe );
        }
    }
}

void myCallback()
{
    // Select a vertex with the mouse
    if (polyscope::pick::haveSelection()) {
        bool goodSelection = false;
        auto selection = polyscope::pick::getSelection();
        auto selectedSurface = static_cast<polyscope::SurfaceMesh*>(selection.first);
        auto idx = selection.second;

        // Only authorize selection on the input surface and the reconstruction
        auto surf = polyscope::getSurfaceMesh("input mesh");
        goodSelection = goodSelection || (selectedSurface == surf);
        const auto nv = selectedSurface->nVertices();
        const auto nf = selectedSurface->nFaces();
        const auto ne = selectedSurface->nEdges();
        // Validate that it its a face index
        if ( goodSelection )
        {
            if ( idx < nv )
            {
                std::ostringstream otext;
                otext << "Selected vertex = " << idx;
                ImGui::Text( "%s", otext.str().c_str() );
            }
            else if ( idx - nv < nf )
            {
                std::ostringstream otext;
                otext << "Selected face = " << ( idx - nv );
                ImGui::Text( "%s", otext.str().c_str() );
            }
            else if ( idx - nv - nf < ne )
            {
                std::ostringstream otext;
                otext << "Selected edge = " << ( idx - nv - nf );
                ImGui::Text( "%s", otext.str().c_str() );
            }
        }
    }
    ImGuiTabBarFlags tab_bar_flags = ImGuiTabBarFlags_None;
    ImGui::BeginTabBar("UV",ImGuiTabBarFlags_None);
    if (ImGui::BeginTabItem("AT faces"))
    {
        ImGui::SliderFloat("Eps_1 (start fuzzyness)", &Epsilon1, 0.0, 5.0);
        ImGui::SliderFloat("Eps_2 (stop fuzzyness)", &Epsilon2, 0.0, 5.0);
        ImGui::SliderFloat("Lambda (discontinuity length)", &Lambda, 0.0, 1.0);
        ImGui::SliderInt("Nb iterations", &NbIterations, 1, 100);

        ImGui::Separator();
        ImGui::Checkbox("Bijective", &bijective);
        ImGui::Checkbox("Use AO", &use_ao);
        ImGui::Checkbox("Force border crossing", &force_border_crossing);
        ImGui::Separator();
        ImGui::SliderFloat("Bounding box x :", &bb_x, 0.0, 5.0);
        ImGui::SliderFloat("Bounding box y :", &bb_y, 0.0, 5.0);
        ImGui::SliderFloat("Bounding box coeff :", &bb_coeff, 0.0, 5.0);
        ImGui::Checkbox("Use bounding box", &use_bb);
        if ((ImGui::Button("Compute whole") || wholeATUVon) && !faceATUVon) {
            if(wholeATUVon) {
                registerNewCut(fUVAT.ATCutAtEdges(V, bijective, prolongate_borders, force_border_crossing));
                fUVAT.finalStep();

                double border_energy = get_border_energy(F, V, fUVAT.F, fUVAT.V.rows());
                Eigen::MatrixXd Xc = fUVAT.UVSolver.getCornerParam();
                registerCornerParam(Xc, "Param with AT");
                wholeATUVon = false;
                DGtal::trace.endBlock();
            } else {
                DGtal::trace.beginBlock("Computing parameterization");
                faceATUVon = true;
                wholeATUVon = true;
                cutDone = false;
                faceATUVi = 0;
                faceATUVe = Epsilon1;
                DGtal::trace.info() << "************************** Epsilon " << faceATUVe
                    << " ***************************" << std::endl;

                fUVAT.build( Epsilon1, Lambda, V, F, bijective, use_ao, use_bb, bb_x, bb_y, bb_coeff);
            }
        }

        if (ImGui::Button("Compute initial AT-UV") && !faceATUVon) {
            faceATUVon = true;
            faceATUVdone = false;
            cutDone = false;
            faceATUVi = 0;
            faceATUVe = Epsilon1;
            DGtal::trace.info() << "************************** Epsilon " << faceATUVe
                << " ***************************" << std::endl;

            std::vector<std::vector<int>> cuts = fUVAT.build( Epsilon1, Lambda, V, F, bijective, use_ao, use_bb, bb_x, bb_y, bb_coeff );
            std::vector<Vector3>               positions;
            std::vector<std::array<size_t, 2>> edgeInds;
            for(auto cut: cuts) {
                int prec = -1;
                for(auto v: cut) {
                    if(prec != -1) {
                        const auto  cpt = positions.size();
                        positions.push_back( geometry->vertexPositions[ prec ] );
                        positions.push_back( geometry->vertexPositions[ v ] );
                        edgeInds.push_back( { cpt, cpt + 1 } );
                    }
                    prec = v;
                }
            }
            psMesh->addSurfaceGraphQuantity("Initial cut", positions, edgeInds);
        }
        if(faceATUVdone) {
            if (ImGui::Button("Cut")) {
                Eigen::MatrixXi cut_mask = fUVAT.ATCutAtEdges(V, bijective, prolongate_borders, force_border_crossing);
                cutDone = true;
                registerNewCut(cut_mask);
                registerTotalCut(cut_mask);
                Eigen::MatrixXd Xc = fUVAT.UVSolver.getCornerParam();
                registerCornerParam(Xc, "Param with AT");
            }
        }
        if(cutDone) {
            if (ImGui::Button("Refine AT-UV") && !faceATUVon) {
                faceATUVon = true;
                faceATUVi = 0;
                faceATUVe = Epsilon1;
                DGtal::trace.info() << "************************** Epsilon " << faceATUVe
                    << " ***************************" << std::endl;
                fUVAT.setEpsilon( Epsilon1 );
            }
            if(ImGui::Button("Finalize AT-UV")) {
                fUVAT.finalStep();
                double border_energy = get_border_energy(F, V, fUVAT.F, fUVAT.V.rows());

                Eigen::MatrixXd Xc = fUVAT.UVSolver.getCornerParam();
                registerCornerParam(Xc, "Param with AT");
            }
        }
        if(ImGui::Button("Compute UV param")) {
            fUVAT.build( Epsilon1, Lambda, V, F, bijective, use_ao, use_bb, bb_x, bb_y, bb_coeff );
            fUVAT.finalStep();
            double border_energy = get_border_energy(F, V, fUVAT.F, fUVAT.V.rows());

            Eigen::MatrixXd Xc = fUVAT.UVSolver.getCornerParam();
            registerCornerParam(Xc, "Param with AT");
        }
        if(faceATUVon)
            doFacesATUV();
        ImGui::EndTabItem();

    }
    ImGui::EndTabBar();

}

inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

double get_border_length(Eigen::MatrixXi const&F, Eigen::MatrixXd const&V, Eigen::MatrixXi const&F_orig, int n) {
    std::vector<std::vector<bool>> border_edges = std::vector<std::vector<bool>>(n, std::vector<bool>(n, false));
    for(int i = 0; i < F.rows(); i++) {
        for(int j = 0; j < 3; j++) {
            int i1 = F(i, j);
            int i2 = F(i, (j+1)%3);
            border_edges[i1][i2] = !border_edges[i1][i2];
            border_edges[i2][i1] = !border_edges[i2][i1];
        }
    }
    double length = 0;
    for(int i = 0; i < F.rows(); i++) {
        for(int j = 0; j < 3; j++) {
            int i1 = F(i, j);
            int i2 = F(i, (j+1)%3);
            if(border_edges[i1][i2])
                length += (V.row(F_orig(i,j)) - V.row(F_orig(i,(j+1)%3))).norm();
        }
    }
    return length;
}

double get_border_energy(Eigen::MatrixXi const&F_orig, Eigen::MatrixXd const&V, Eigen::MatrixXi const&F, int n) {
    double totArea = 0.;
    double l_orig = get_border_length(F_orig, V, F_orig, V.rows());
    double l_new = get_border_length(F, V, F_orig, n);
    for(int i = 0; i < F.rows(); i++) {
        totArea += getArea(V, F_orig, i);
    }
    
    return 0.5*(l_new - l_orig) / sqrt(totArea / 3.14159265359);
}

void loadUV(char **argv) {
    std::vector<double> us;
    std::vector<double> vs;
    Eigen::MatrixXi F2 = Eigen::MatrixXi(F.rows(), 3);
    std::ifstream infile(argv[2]);
    std::string line;
    int i = 0;
    int i2 = 0;
    while(std::getline(infile, line)) {
        std::istringstream iss(line);
        std::vector<std::string> results((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
        if(results.size() > 0) {
            if(results[0] == "vt") {
                double u = stod(results[1]);
                double v = stod(results[2]);
                us.push_back(u);
                vs.push_back(v);
                //X(i,0) = u;
                //X(i,1) = v;
                i++;
            }
            if(results[0] == "f" && i > 0) {
                std::string delimiter = "/";
                results[1].erase(0, results[1].find(delimiter)+delimiter.length());
                results[2].erase(0, results[2].find(delimiter)+delimiter.length());
                results[3].erase(0, results[3].find(delimiter)+delimiter.length());
                int f1 = stoi(results[1])-1;
                int f2 = stoi(results[2])-1;
                int f3 = stoi(results[3])-1;
                F2(i2,0) = f1;
                F2(i2,1) = f2;
                F2(i2,2) = f3;
                i2++;
            }
        }
    }
    if(i > 0) {
        Eigen::MatrixXd X = Eigen::MatrixXd(i, 2);
        for(int j = 0; j < us.size(); j++) {
            X(j,0) = us[j];
            X(j,1) = vs[j];
        }
        Eigen::MatrixXd X2 = Eigen::MatrixXd::Zero(i, 3);
        X2.col(0) = X.col(0);
        X2.col(1) = X.col(1);

        Eigen::MatrixXd Xcorner = Eigen::MatrixXd(3*F2.rows(), 2);
        for(int j = 0; j < F2.rows(); j++) {
            Xcorner.row(3*j) = X.row(F2(j,0));
            Xcorner.row(3*j+1) = X.row(F2(j,1));
            Xcorner.row(3*j+2) = X.row(F2(j,2));
        }
        polyscope::getSurfaceMesh("input mesh")
            ->addParameterizationQuantity("Loaded param", Xcorner);
        polyscope::registerSurfaceMesh("Loaded UV corner mesh", X2, F2)->setEdgeWidth(1.0)->setEdgeColor({0.,0.,0.});
        Eigen::MatrixXi cut_mask = Eigen::MatrixXi::Zero(F.rows(), 3);
        Eigen::MatrixXi TT;
        Eigen::MatrixXi TTi;
        igl::triangle_triangle_adjacency(F2, TT, TTi);
        Eigen::MatrixXi TT_old;
        Eigen::MatrixXi TTi_old;
        igl::triangle_triangle_adjacency(F, TT_old, TTi_old);
        for(int i = 0; i < F.rows(); i++) {
            for(int j = 0; j < 3; j++) {
                if(TT(i,j) == -1 && TT_old(i,j) != -1)
                    cut_mask(i,j) = 1;
            }
        }
        registerLoadedCut(cut_mask);
        std::cout << "Loaded cut length : " << get_border_energy(F, V, F2, i) << std::endl;
    }
}

void headless_compute(std::ofstream &ofs) {
    bool start = true;
    Clock c;
    long duration;
    while(start || wholeATUVon) {
        if ((start || wholeATUVon) && !faceATUVon) {
            start = false;
            if(wholeATUVon) {
                registerNewCut(fUVAT.ATCutAtEdges(V, bijective, prolongate_borders, force_border_crossing));
                double distorsion_energy = fUVAT.finalStep();
                double border_energy = get_border_energy(F, V, fUVAT.F, fUVAT.V.rows());
                duration = c.stopClock();
                ofs << distorsion_energy << ", " << border_energy << ", " << duration << std::endl;
                Eigen::MatrixXd Xc = fUVAT.UVSolver.getCornerParam();
                registerCornerParam(Xc, "Param with AT");
                wholeATUVon = false;
                DGtal::trace.endBlock();
                return;
            } else {
                DGtal::trace.beginBlock("Computing parameterization");
                c.startClock();
                faceATUVon = true;
                wholeATUVon = true;
                cutDone = false;
                faceATUVi = 0;
                faceATUVe = Epsilon1;
                DGtal::trace.info() << "************************** Epsilon " << faceATUVe
                    << " ***************************" << std::endl;

                fUVAT.build( Epsilon1, Lambda, V, F, bijective, use_ao, use_bb, bb_x, bb_y, bb_coeff);
            }
        }
        if(faceATUVon)
            doFacesATUV();
    }
}

int main(int argc, char **argv)
{
    // Initializing polyscope
    polyscope::init();
    polyscope::state::userCallback = myCallback;

    CLI::App app;
    bool withoutGUI=false;
    app.add_flag("--withoutGUI",withoutGUI,"Disabling the GUI");
    std::string filename;
    app.add_option("-i", filename, "File name");
    std::cout << filename << std::endl;
    CLI11_PARSE(app, argc, argv);
    SMeshBuilder builder;
    bool ok_parse = builder.parseCommandLine( argc, argv );
    if ( ! ok_parse ) return 0;
    builder.buildInput();
    builder.computeNormals();
    builder.fixVertexNormals();

    std::string meshname = std::filesystem::path(builder.filename).stem();
    mkdir("output", 0777);
    mkdir(("output/" + meshname).c_str(), 0777);
    std::ofstream ofs ("output/" + meshname + "/results-uv-at.dat", std::ofstream::out);

    // Initializing corrected DEC
    surfmesh = builder.smesh;
    ofs << surfmesh.nbVertices() << ", " << surfmesh.nbFaces() << std::endl;
    dec.init( &surfmesh );
    dec.requireDiscreteExteriorCalculus();

    // Load Polyscope mesh
    //std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));
    trace.info() << "Build mesh" << std::endl;

    mesh = std::unique_ptr<GCSurfaceMesh>
        ( new GCSurfaceMesh( builder.smesh.allIncidentVertices() ) );
    VertexData< Vector3 > pos( *mesh );
    SMesh::Vertex i = 0;
    for ( auto v : mesh->vertices() )
    {
        auto p = builder.smesh.position( i++ );
        pos[ v ] = Vector3 { p[ 0 ], p[ 1 ], p[ 2 ] };
    }
    trace.info() << "Build geometry" << std::endl;
    geometry = std::unique_ptr<VertexPositionGeometry>
        ( new VertexPositionGeometry( *mesh, pos ) );

    trace.info() << "geometry->requireVertexPositions()" << std::endl;
    geometry->requireVertexPositions();

    V = EigenMap<double, 3, Eigen::RowMajor>(pos);
    F = mesh->getFaceVertexMatrix<int>();
    // Register the mesh with polyscope
    trace.info() << "Register mesh" << std::endl;
    psMesh = polyscope::registerSurfaceMesh
        ( "input mesh",
          V, F);

    trace.info() << "showcase!" << std::endl;

    loadUV(argv);

    cut_history = Eigen::MatrixXi::Zero(F.rows(), 3);

    if(!withoutGUI) {
        polyscope::show();
    } else {
        headless_compute(ofs);
    }
    return 0;
}
