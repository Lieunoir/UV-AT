/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file SurfaceMeshBuilder.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2020/07/11
 *
 * Header file for module SurfaceMeshBuilder.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(SurfaceMeshBuilder_RECURSES)
#error Recursive header files inclusion detected in SurfaceMeshBuilder.h
#else // defined(SurfaceMeshBuilder_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SurfaceMeshBuilder_RECURSES

#if !defined SurfaceMeshBuilder_h
/** Prevents repeated inclusion of headers. */
#define SurfaceMeshBuilder_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <sstream>
#include <string>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/SurfaceMesh.h"
#include "DGtal/shapes/SurfaceMeshHelper.h"
#include "DGtal/io/readers/SurfaceMeshReader.h"
#include "DGtal/io/writers/SurfaceMeshWriter.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/helpers/ShortcutsGeometry.h"

namespace DGtal
{
  namespace po  = boost::program_options;

  /////////////////////////////////////////////////////////////////////////////
  // template class SurfaceMeshBuilder
  /**
     Description of template class 'SurfaceMeshBuilder' <p> \brief Aim:
     Helper class to build a SurfaceMesh from different inputs.

     @tparam TRealPoint an arbitrary model of 3D RealPoint.
     @tparam TRealVector an arbitrary model of 3D RealVector.
  */
  template < typename TRealPoint, typename TRealVector >
  struct SurfaceMeshBuilder
  {
    typedef TRealPoint                              RealPoint;
    typedef TRealVector                             RealVector;
    typedef SurfaceMeshBuilder< RealPoint, RealVector >    Self;
    
    static const Dimension dimension = RealPoint::dimension;
    BOOST_STATIC_ASSERT( ( dimension == 3 ) );

    typedef typename RealVector::Component        Scalar;
    typedef std::vector<Scalar>                   Scalars;
    typedef std::vector< RealPoint >              RealPoints;
    typedef std::vector< RealVector >             RealVectors;
    typedef SurfaceMesh<RealPoint,RealVector>     SMesh;
    typedef SurfaceMeshReader<RealPoint,RealVector> SMeshReader;
    typedef SurfaceMeshWriter<RealPoint,RealVector> SMeshWriter;
    typedef SurfaceMeshHelper<RealPoint,RealVector> SMeshHelper;
    typedef typename SMesh::Size                  Size;
    typedef typename SMesh::Index                 Index;
    typedef typename SMesh::Face                  Face;
    typedef typename SMesh::Edge                  Edge;
    typedef typename SMesh::Vertex                Vertex;
    typedef Z3i::Space                            Space;
    typedef Z3i::KSpace                           KSpace;
    typedef Shortcuts<KSpace>                     SH;
    typedef ShortcutsGeometry<KSpace>             SHG;
    typedef SH::IdxDigitalSurface                 Surface;
    typedef Surface::DigitalSurfaceContainer      Container;
    typedef SH::BinaryImage                       BinaryImage;
    typedef SH::ImplicitShape3D                   ImplicitShape3D;
    typedef SH::DigitalSurface                    DigitalSurface;
    typedef SH::IdxDigitalSurface                 IdxDigitalSurface;
    typedef GradientColorMap<Scalar>              ColorMap;
    typedef typename SMesh::Vertices              Vertices;
    
    //---------------------------------------------------------------------------
  public:
    /// @name Standard services
    /// @{
    
    /// Default destructor.
    ~SurfaceMeshBuilder() = default;
    /// Default copy constructor.
    /// @param other the object to clone
    SurfaceMeshBuilder( const Self& other ) = default;
    /// Default move constructor.
    /// @param other the object to move
    SurfaceMeshBuilder( Self&& other ) = default;
    /// Default assignment constructor.
    /// @param other the object to clone
    /// @return a reference to 'this'.
    Self& operator=( const Self& other ) = default;

    /// @}

    /// Creates the object with the meaningful options.
    SurfaceMeshBuilder();
    /// Parses command-line options.
    bool parseCommandLine( int argc, char** argv );
    /// Build input from data for further processing.
    bool buildInput();
    /// Build inputs from implicit polynomial shapes or 3D .vol image file.
    bool buildPolynomialOrVolInput();
    /// Build inputs from a set of predefined meshes
    bool buildPredefinedMesh();
    /// Build inputs from a polygonal mesh .obj file
    bool buildObjInput();

    /// Compute normals
    void computeNormals();
    /// Fix corrected normals wrt naive normals, so that their inner
    /// product is no smaller than \ref fix_normals_val.
    void fixVertexNormals();

    /// Perturbates positions
    void perturbatePositions();

    //---------------------------------------------------------------------------
  public:
    /// Parses command line given as \a argc, \a argv according to \a
    /// opts and fills in storage map \a vm.
    ///
    /// @param[in] opts the allowed options.
    /// @param[in] argc the size of array \a argv.
    ///
    /// @param[in] argv an array of C-style strings containing the
    /// desired command-line options.
    ///
    /// @param[inout] vm the map options to variables that is updated.
    static bool args2vm( const po::options_description& opts,
			 int argc, char** argv,
			 po::variables_map& vm );
    /**
     * Missing parameter error message.
     *
     * @param param
     */
    static void missingParam( std::string param )
    {
      trace.error() << " Parameter: " << param << " is required.";
      trace.info() << std::endl;
    }

    /// Subdivides the current mesh according to parameter \a subdivide
    /// @param subdivide either "NO", "CENTROID", "CANONIC"
    void subdivideCurrentMesh( std::string subdivide );
    
    /// Subdivides a mesh given as a range of faces and vertex
    /// positions in order to have a triangulated mesh.
    ///
    /// @param[in] input_faces the range of input faces, each face is
    /// a range of indicent vertices.
    ///
    /// @param[in] input_positions the range of input vertex positions.
    ///
    /// @param[in] centroid_subdivision when 'true', subdivides each
    /// face by adding a vertex in the middle of each face, otherwise
    /// makes an umbrella of triangles starting from first vertex of
    /// face.
    ///
    /// @param[out] output_faces the range of output faces, each face is
    /// a range of indicent vertices.
    ///
    /// @param[out] output_positions the range of output vertex positions.
    static void subdivideMesh
    ( const std::vector< Vertices >& input_faces,
      const RealPoints&              input_positions,
      bool                           centroid_subdivision,
      std::vector< Vertices >&       output_faces,
      RealPoints&                    output_positions );

    /// Subdivides a mesh given as a range of faces and vertex
    /// positions in order to have a triangulated mesh. It also
    /// computes the new vertex normals.
    ///
    /// @param[in] input_faces the range of input faces, each face is
    /// a range of indicent vertices.
    ///
    /// @param[in] input_positions the range of input vertex positions.
    ///
    /// @param[in] input_vnormals the range of input vertex normals.
    ///
    /// @param[in] centroid_subdivision when 'true', subdivides each
    /// face by adding a vertex in the middle of each face, otherwise
    /// makes an umbrella of triangles starting from first vertex of
    /// face.
    ///
    /// @param[out] output_faces the range of output faces, each face is
    /// a range of indicent vertices.
    ///
    /// @param[out] output_positions the range of output vertex positions.
    ///
    /// @param[out] output_vnormals the range of output vertex normals.
    static void subdivideMeshWithVertexNormals
    ( const std::vector< Vertices >& input_faces,
      const RealPoints&              input_positions,
      const RealVectors&             input_vnormals,
      bool                           centroid_subdivision,
      std::vector< Vertices >&       output_faces,
      RealPoints&                    output_positions,
      RealVectors&                   output_vnormals );

    /// Subdivides a mesh given as a range of faces and vertex
    /// positions in order to have a triangulated mesh. It also
    /// computes the new face normals.
    ///
    /// @param[in] input_faces the range of input faces, each face is
    /// a range of indicent vertices.
    ///
    /// @param[in] input_positions the range of input vertex positions.
    ///
    /// @param[in] input_fnormals the range of input face normals.
    ///
    /// @param[in] centroid_subdivision when 'true', subdivides each
    /// face by adding a vertex in the middle of each face, otherwise
    /// makes an umbrella of triangles starting from first vertex of
    /// face.
    ///
    /// @param[out] output_faces the range of output faces, each face is
    /// a range of indicent vertices.
    ///
    /// @param[out] output_positions the range of output vertex positions.
    ///
    /// @param[out] output_fnormals the range of output face normals.
    static void subdivideMeshWithFaceNormals
    ( const std::vector< Vertices >& input_faces,
      const RealPoints&              input_positions,
      const RealVectors&             input_fnormals,
      bool                           centroid_subdivision,
      std::vector< Vertices >&       output_faces,
      RealPoints&                    output_positions,
      RealVectors&                   output_fnormals );
    
    //---------------------------------------------------------------------------
  public:
    
    /// Meaningful options for this object.
    po::options_description general_opt;
    /// Parsed command-line
    po::variables_map vm;
    /// command-line as string
    std::string       command_line;
    /// Input filename
    std::string       filename;
    /// Mesh name
    std::string       meshname;
    /// Mesh args
    std::vector< std::string > meshargs;
    /// Gridstep
    double            h;
    /// Subdivision parameters for predefined meshes.
    double            sub_m, sub_n;
    /// Parameters object
    Parameters        params;
    /// The simplified mesh on which all computations are done
    SMesh              smesh;
    
    bool in_obj;
    bool in_vol;
    bool in_mesh;
    bool in_polynomial;

    KSpace                        K;   // Digital space.
    CountedPtr<BinaryImage>       bimage;
    CountedPtr<ImplicitShape3D>   shape;
    CountedPtr<DigitalSurface>    surface;
    CountedPtr<IdxDigitalSurface> idx_surface;
    bool blocky;
    bool dual;

    /// Stores expected normal vectors (if applicable)
    RealVectors   expected_normals;
    /// Stores measured normal vectors (if applicable)
    RealVectors   measured_normals;
    /// Tells how to force corrected normals to have positive dot
    /// product with naive normals
    double        fix_normals_val;
    
  }; // end of class SurfaceMeshBuilder
  
} // namespace DGtal

///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "SurfaceMeshBuilder.ih"
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SurfaceMeshBuilder_h

#undef SurfaceMeshBuilder_RECURSES
#endif // else defined(SurfaceMeshBuilder_RECURSES)

