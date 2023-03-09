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
 * @file SurfaceMeshDEC.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2020/07/16
 *
 * Header file for module SurfaceMeshDEC.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(SurfaceMeshDEC_RECURSES)
#error Recursive header files inclusion detected in SurfaceMeshDEC.h
#else // defined(SurfaceMeshDEC_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SurfaceMeshDEC_RECURSES

#if !defined SurfaceMeshDEC_h
/** Prevents repeated inclusion of headers. */
#define SurfaceMeshDEC_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <sstream>
#include <string>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/SurfaceMesh.h"

namespace DGtal
{
  /////////////////////////////////////////////////////////////////////////////
  // template class SurfaceMeshDEC
  /**
     Description of template class 'SurfaceMeshDEC' <p> \brief Aim:
     Provides methods and operators that builds a discrete exterior
     calculus onto a triangulated surface mesh.

   * @tparam TLinearAlgebraBackend linear algebra backend used (i.e. EigenSparseLinearAlgebraBackend).
     @tparam TRealPoint an arbitrary model of 3D RealPoint.
     @tparam TRealVector an arbitrary model of 3D RealVector.
  */
  template < typename TLinearAlgebraBackend,
             typename TRealPoint, typename TRealVector >
  struct SurfaceMeshDEC
  {
    typedef TLinearAlgebraBackend LinearAlgebraBackend;
    typedef TRealPoint            RealPoint;
    typedef TRealVector           RealVector;
    typedef typename LinearAlgebraBackend::DenseVector::Index  Index;
    typedef typename LinearAlgebraBackend::DenseVector::Scalar Scalar;
    typedef typename LinearAlgebraBackend::DenseVector  DenseVector;
    typedef typename LinearAlgebraBackend::DenseMatrix  DenseMatrix;
    typedef typename LinearAlgebraBackend::SparseMatrix SparseMatrix;
    typedef SurfaceMesh< RealPoint, RealVector >        Mesh;
    typedef typename Mesh::Vertex                       Vertex;
    typedef typename Mesh::Edge                         Edge;
    typedef typename Mesh::Face                         Face;
    typedef typename Mesh::Face                         Corner;
    typedef typename Mesh::Size                         Size;

    typedef DenseVector  Form;
    typedef SparseMatrix LinearOperator;

    /// A triplet (row, col, value), useful for initializaing sparse matrices.
    typedef typename LinearAlgebraBackend::Triplet Triplet;
    /// A range of triplets (row,col,value).
    typedef std::vector<Triplet>                   Triplets;

    // ------------------------- Initialization services ----------------------------
  public:
    /// @name Initialization services
    /// @{
    
    /// Default constructor. The object is invalid.
    SurfaceMeshDEC()
      : myMesh( nullptr ),
        has_corrected_calculus( false ),
        has_discrete_exterior_calculus( false ),
        myNbTriangles( 0 ), myNbQuadrangles( 0 )
    {}
    
    /// Constructor from surface mesh \a smesh.
    /// @param smesh any surface mesh
    SurfaceMeshDEC( ConstAlias< Mesh > smesh )
      : myMesh( nullptr ),
        has_corrected_calculus( false ),
        has_discrete_exterior_calculus( false ),
        myNbTriangles( 0 ), myNbQuadrangles( 0 )
    {
      init( myMesh );
    }

    /// Global initialization from pointer on SurfaceMesh \a pMesh
    /// @param pMesh any pointer on a valid SurfaceMesh.
    void init( const Mesh* pMesh )
    {
      has_corrected_calculus = false;
      has_discrete_exterior_calculus = false;
      if ( pMesh == nullptr )
        {
          clear();
          return;
        }
      myMesh = pMesh;
      initIndexMaps();
      trace.info() << "[SurfaceMeshDEC::init] #V=" << nbCells( 0 )
                   << " #E=" << nbCells( 1 )
                   << " #F=" << nbCells( 2 )
                   << " #C=" << nbCorners() << std::endl;
      initFaceMatrices();
      initDerivativeOperators();
      trace.beginBlock( "Init average operators" );
      computeAverageOperators();
      trace.endBlock();
    }

    /// Clears everything.
    void clear()
    {
      myMesh = nullptr;
      has_corrected_calculus = false;
      has_discrete_exterior_calculus = false;
    }

    /// @}
    
    // ------------------------- Accessors ------------------------------
  public:
    /// @name Accesors
    /// @{
    
    /// @return a pointer to the associated mesh.
    const Mesh* mesh() const
    { return myMesh; }

    /// @return 'true' iff the associated mesh is triangulated
    bool isMeshTriangulated() const
    { return myNbTriangles   == nbCells( 2 ); }
      
    /// @return 'true' iff the associated mesh is quadrangulated
    bool isMeshQuadrangulated() const
    { return myNbQuadrangles == nbCells( 2 ); }

    /// @param d the dimension of cells in 0,1,2
    /// @return the number of \a d-dimensional cells.
    Size nbCells( Dimension d ) const
    {
      ASSERT( myMesh != nullptr );
      switch ( d ) {
      case 0: return myMesh->nbVertices();
      case 1: return myMesh->nbEdges();
      case 2: return myMesh->nbFaces();
      default:
        trace.warning() << "[SurfaceMeshDEC::nbCells] Invalid dimension "
                        << d << std::endl;
        return 0;
      }
    }

    /// @return the number of face corners (i.e. 3 times the number of
    /// faces for a triangulated surface).
    Index nbCorners() const
    {
      return idxFtoC( nbCells( 2 ) );
    }
    
    /// @param f any face index.
    /// @return the index of the first corner of face \a f.
    Corner idxFtoC( const Face f ) const
    { return myFtoCIdx[ f ]; }

    /// @param f any face index.
    /// @return the number of corners of face \a f.
    Size nbCorners( const Face f ) const
    { return myNbCorners[ f ]; }

    /// @param c any corner index.
    /// @return the index of the face containing corner \a c.
    Face idxCtoF( const Corner c ) const
    { return myCtoFIdx[ c ]; }

     /// @param c any face index.
     /// @return the index of the vertex assocaited with the corner \a c.
     Vertex idxCtoV( const Corner c ) const
     { const auto idxF = idxCtoF( c );
       const auto vtcs = mesh()->incidentVertices( idxF );
       const auto idxS = c - idxFtoC( idxF );
       return vtcs[ idxS ];  }

    /// @}

    // ------------------------- Geometric services ---------------------------
  public:
    /// @name Geometric services
    /// @{

    /// @param f any index of triangular or quadrangular face.
    /// @return the area of the given face
    double faceArea( Face f ) const;
    /// @param e any index of edge
    /// @return the length of this edge.
    double edgeLength( Edge e ) const;
    /// @param e any index of edge
    /// @return the length of the dual of this edge.
    double edgeDualLength( Edge e ) const;
    /// @param v any index of vertex
    /// @return (an approximation of) the area of the face dual to vertex v.
    double vertexDualArea( Vertex v ) const;
    
    /// @param f any index of triangular or quadrangular face.
    /// @return the corrected area of the given face
    double faceCorrectedArea( Face f ) const;
    /// @param f any index of triangular face.
    /// @return the corrected area of the given face
    double triangleCorrectedArea( Face f ) const;
    /// @param f any index of quadrangular face.
    /// @return the corrected area of the given face
    double quadrangleCorrectedArea( Face f ) const;

    /// @param f any index of face (triangular, quadrangular or anything).
    /// @return the exact or approximate circumcenter to the given face vertices.
    RealPoint faceCircumcenter( Face f ) const;
    
    /// @param a any point
    /// @param b any point
    /// @param c any point
    /// @return the circumcenter of points \a a, \a b, \a c.
    RealPoint circumcenter( RealPoint a, RealPoint b, RealPoint c ) const;

    /// @param p any point
    /// @param a any point
    /// @param b any point
    /// @param c any point
    
    /// @return 'true' if \a p projected on abc lies inside or on the
    /// boundary of triangle abc.
    bool isInTriangle( RealPoint p,
                       RealPoint a, RealPoint b, RealPoint c ) const;

    /// @param p any point
    /// @param f any face
    /// @return 'true' if \a p projected on face \a f lies inside or on the
    /// boundary of face \a f.
    bool isInFace( RealPoint p, Face f ) const;
        
    /// @}
    
    // ------------------------- Generic Calculus operators --------------------------
  public:
    /// @name Generic Calculus operators
    /// @{
    
    /// @return the topological derivative operator D0: F0 -> F1
    const LinearOperator& D0() const
    { return myD0; }
    /// @return the topological derivative operator D1: F1 -> F2
    const LinearOperator& D1() const
    { return myD1; }

    /// @return the identity operator for 0-forms
    LinearOperator Id0() const
    {
      LinearOperator Id( nbCells( 0 ), nbCells( 0 ) );
      Id.setIdentity();
      return Id;
    }
    /// @return the identity operator for 1-forms
    LinearOperator Id1() const
    {
      LinearOperator Id( nbCells( 1 ), nbCells( 1 ) );
      Id.setIdentity();
      return Id;
    }
    /// @return the identity operator for 2-forms
    LinearOperator Id2() const
    {
      LinearOperator Id( nbCells( 2 ), nbCells( 2 ) );
      Id.setIdentity();
      return Id;
    }

    /// @param n the chosen dimension for values in corners.
    /// @return the identity operator for C-forms
    LinearOperator IdC( const Dimension n = 1 ) const
    {
      LinearOperator Id( nbCells( 2 ) * n, nbCells( 2 ) * n );
      Id.setIdentity();
      return Id;
    }
    
    /// @return a 0-form with all values set to zero
    Form form0() const
    { return Form::Zero( nbCells( 0 ) ); }

    /// @return a 1-form with all values set to zero
    Form form1() const
    { return Form::Zero( nbCells( 1 ) ); }

    /// @return a 2-form with all values set to zero
    Form form2() const
    { return Form::Zero( nbCells( 2 ) ); }

    /// @param n the chosen dimension for values in corners.
    /// @return a Corner-form of given dimension with all values set to zero
    Form formC( const Dimension n = 1 ) const
    { return Form::Zero( nbCorners() * n ); }
    
    /// @}

    // ------------------------- Corrected calculus services ------------------------
  public:
    /// @name Corrected calculus services
    /// @{
    
    /// Defines operators for corrected calculus.
    void requireCorrectedCalculus()
    {
      if ( ! has_corrected_calculus )
        {
          trace.beginBlock( "Init corrected calculus" );
          initCorrectedCalculusOperators();
          has_corrected_calculus = true;
          trace.endBlock();
        }
    }      

    /// @return 'true' if operators for corrected calculus are computed.
    bool hasCorrectedCalculus() const
    { return has_corrected_calculus; }

    /// @return 'true' if operators for corrected calculus are computed.
    void unrequireCorrectedCalculus()
    {
      if ( has_corrected_calculus )
        {
          myCCM0.clear();
          myCCDiagH0.clear();
          myCCDiagDualH2.clear();
          myCCM1.clear();
          myCCDiagDualH1.clear();
          myCCM2.clear();
          myCCDualH0.clear();
          myCCC1IP.clear();
          myCCSharp.clear();
          myCCFlat.clear();
          myCCLL0.clear();
          has_corrected_calculus = false;
        }
    }
        
    /// @return the inner product / mass matrix M0 for 0-form,
    /// i.e. <f,g>_0 := f^T * M0 * g
    const LinearOperator& ccM0() const
    { return myCCM0; }

    /// @return the hodge star for 0-form, which is also the mass matrix M0.
    /// since \f$ \int f \wedge \star_0 g := <f,g>_0 := f^T * M0 * g \f$
    const LinearOperator& ccH0() const
    { return myCCM0; }

    /// @return the diagonal (simplified) hodge star for 0-form
    /// i.e. <f,g>_0 := f^T * ccDiagH0.asDiagonal() * g
    const Form& ccDiagH0() const
    { return myCCDiagH0; }

    /// @return the diagonal hodge star for dual 2-forms.
    const Form& ccDiagDualH2() const
    { return myCCDiagDualH2; }
    
    /// @return the inner product / mass matrix M1 for 1-forms,
    /// i.e. <a,b>_1 := a^T * M1 * b
    const LinearOperator& ccM1() const
    { return myCCM1; }

    /// @return the hodge star for 1-form, which is also the mass matrix M1.
    /// since \f$ \int a \wedge \star_1 b := <a,b>_1 := a^T * M1 * b \f$
    const LinearOperator& ccH1() const
    { return myCCM1; }

    /// @return the diagonal (simplified) Hodge star for 1-forms (use
    /// `asDiagonal()` for operator).
    const Form& ccDiagH1() const
    { return myCCDiagH1; }

    /// @return the diagonal Hodge star for dual 1-forms (use
    /// `asDiagonal()` for operator).
    const Form& ccDiagDualH1() const
    { return myCCDiagDualH1; }

    /// @return the inner product / mass matrix M2 for 2-forms, a
    /// simple diagonal matrix made of inverse of corrected face areas.
    /// i.e. <u,v>_2 := u^T * M2 * v
    const LinearOperator& ccM2() const
    { return myCCM2; }

    /// @return the hodge star for 2-forms, a
    /// simple diagonal matrix made of inverse of corrected face areas, equal to ccM2
    const LinearOperator& ccH2() const
    { return myCCM2; }

    /// @return the hodge star for dual 0-forms
    const LinearOperator& ccDualH0() const
    { return myCCDualH0; }
        
    /// @return the left laplacian LL0 = -1.0 * D0^T * M1 * D0; ie
    /// Lap0 * x = b <=> LL0 * x = M0 * b
    const LinearOperator& ccLL0() const
    { return myCCLL0; }

    /// @return the sharp operator from E to 9T=3C
    const LinearOperator& ccSharp() const
    { return myCCSharp; }

    /// @return the flat operator from 9T=3C to E
    const LinearOperator& ccFlat() const
    { return myCCFlat; }

    /// @return thr C-Vectors inner product matrix.
    const LinearOperator& ccC1IP() const
    { return myCCC1IP; }
    
    
    /// @return The non-conformal edge based Laplacian
    /// (E x E linear operator)
    LinearOperator edgeLaplacian() const;
    /// @return The non-conformal edge inner product
    /// (E x E linear operator)
    LinearOperator edgeInnerProduct() const;
    
    /// @}

    // ------------------ Discrete exterior calculus services ------------------------
  public:
    /// @name Discrete exterior calculus services
    /// @{
    
    /// Defines operators for discrete exterior calculus.
    void requireDiscreteExteriorCalculus()
    {
      if ( ! has_discrete_exterior_calculus )
        {
          trace.beginBlock( "Init discrete exterior calculus" );
          initDiscreteExteriorCalculusOperators();
          has_discrete_exterior_calculus = true;
          trace.endBlock();
        }
    }      

    /// @return 'true' if operators for discrete exterior calculus are computed.
    bool hasDiscreteExteriorCalculus() const
    { return has_discrete_exterior_calculus; }

    /// @return 'true' if operators for discrete exterior calculus are computed.
    void unrequireDiscreteExteriorCalculus()
    {
      if ( has_discrete_exterior_calculus )
        {
          myDECM0.clear();
          myDECM1.clear();
          myDECM2.clear();
          myDECDualH2.clear();
          myDECDualH1.clear();
          myDECDualH0.clear();
          myDECLL0.clear();
          myDECL0.clear();
          myDECFaceGrad0x.clear();
          myDECFaceGrad0y.clear();
          myDECFaceGrad0z.clear();
          has_discrete_exterior_calculus = false;
        }
    }

    /// @return the DEC inner product / mass matrix M0 for 0-form,
    /// i.e. <f,g>_0 := f^T * M0 * g, and equivalently the Hodge star
    /// for 0-forms.
    const LinearOperator& decM0() const
    { return myDECM0; }

    /// @return the DEC hodge star for 0-forms, which is also the mass matrix M0.
    const LinearOperator& decH0() const
    { return myDECM0; }

    /// @return the DEC hodge star for dual 2-forms, which is also the
    /// inverse of mass matrix M0.
    const LinearOperator& decDualH2() const
    { return myDECDualH2; }

    /// @return the DEC inner product / mass matrix M1 for 1-form,
    /// i.e. <a,b>_1 := a^T * M1 * b, and equivalently the Hodge star
    /// for 1-forms.
    const LinearOperator& decM1() const
    { return myDECM1; }

    /// @return the DEC hodge star for 1-forms, which is also the mass matrix M1.
    const LinearOperator& decH1() const
    { return myDECM1; }

    /// @return the DEC hodge star for dual 1-forms, which is also the negative of the
    /// inverse of mass matrix M1.
    const LinearOperator& decDualH1() const
    { return myDECDualH1; }

    /// @return the DEC inner product / mass matrix M2 for 2-forms,
    /// i.e. <u,v>_2 := u^T * M2 * v, and equivalently the Hodge star
    /// for 2-forms.
    const LinearOperator& decM2() const
    { return myDECM2; }

    /// @return the DEC hodge star for 2-forms, which is also the mass matrix M2.
    const LinearOperator& decH2() const
    { return myDECM2; }

    /// @return the DEC hodge star for dual 2-forms, which is also the
    /// inverse of mass matrix M2.
    const LinearOperator& decDualH0() const
    { return myDECDualH0; }

    /// @return the DEC Left Laplacian for 0-forms.
    const LinearOperator& decLL0() const
    { return myDECLL0; }
    
    /// @return the DEC Laplacian for 0-forms.
    const LinearOperator& decL0() const
    { return myDECL0; }

    /// @return the DEC x-gradient per-face operator for 0-forms.
    const LinearOperator& decFaceGrad0x() const
    { return myDECFaceGrad0x; }
    /// @return the DEC y-gradient per-face operator for 0-forms.
    const LinearOperator& decFaceGrad0y() const
    { return myDECFaceGrad0y; }
    /// @return the DEC z-gradient per-face operator for 0-forms.
    const LinearOperator& decFaceGrad0z() const
    { return myDECFaceGrad0z; }
    
    /// @}

    // ------------------ Piecewise constant calculus services ------------------------
  public:
    /// @name Piecewise constant calculus services (pcc)
    /// @{
    
    /// Defines operators for piecewise constant calculus.
    void requirePiecewiseConstantCalculus( bool trivial_metric = false )
    {
      if ( ! has_piecewise_constant_calculus )
        {
          trace.beginBlock( "Init piecewise constant calculus" );
          initPiecewiseConstantCalculusOperators( trivial_metric );
          has_piecewise_constant_calculus = true;
          trace.endBlock();
        }
    }      

    /// @return 'true' if operators for piecewise constant calculus are computed.
    bool hasPiecewiseConstantCalculus() const
    { return has_piecewise_constant_calculus; }

    /// Clears operators for piecewise constant calculus if they were computed.
    void unrequirePiecewiseConstantCalculus()
    {
      if ( has_piecewise_constant_calculus )
        {
          myPCCMU.clear();
          myPCCMV.clear();
          myPCCMW.clear();
          myPCCDivM.clear();
          myPCCGradE.clear();
          myPCCDivE.clear();
          myPCCLapE.clear();
          has_piecewise_constant_calculus = false;
        }
    }

    /// @return the inner product for U-forms, which can be seen as an
    /// inner product for dual 0-forms. It is simply the diagonal of
    /// face areas.
    const LinearOperator& pccMU() const
    { return myPCCMU; }

    /// @return the inner product for V-forms, which can be seen as an
    /// inner product for dual 1-forms.  It is simply the diagonal of
    /// edge lengths.
    const LinearOperator& pccMV() const
    { return myPCCMV; }

    /// @return the inner product for W-forms, which can be seen as an
    /// inner product for corner forms.  It is simply the diagonal of
    /// the length of the line joining corners to the face barycenter.
    const LinearOperator& pccMW() const
    { return myPCCMW; }
    
    /// @return the gradient M operator U -> V for U-forms, which can
    /// be seen as a derivative operator for dual 0-forms.
    LinearOperator pccGradM() const
    { return myD1.transpose(); }

    /// @return the divergence M operator V -> U for V-forms, which
    /// can be seen as anti derivative operator for dual 1-forms. It
    /// is defined as the negative adjoint of \ref pccGradM, using
    /// inner products \ref pccMU and \ref pccMV.
    const LinearOperator& pccDivM() const
    { return myPCCDivM; }

    /// @return the Laplacian M operator U -> U for U-forms, which is
    /// the composition of `pccDivM * pccGradM`.
    LinearOperator pccLapM() const
    { return pccDivM() * pccGradM(); }

    /// @return the gradient E operator V -> W for V-forms, which can
    /// be seen as a kind of derivative operator for dual 1-forms to
    /// corner forms.
    const LinearOperator& pccGradE() const
    { return myPCCGradE; }

    /// @return the divergence E operator W -> V for W-forms, which
    /// can be seen as anti derivative operator for corner forms to dual 1-forms. It
    /// is defined as the negative adjoint of \ref pccGradE, using
    /// inner products \ref pccMV and \ref pccMW.
    const LinearOperator& pccDivE() const
    { return myPCCDivE; }

    /// @return the Laplacian E operator V -> V for V-forms, which is
    /// the composition of `pccDivE * pccGradE`.
    LinearOperator pccLapE() const
    { return myPCCLapE; /* pccDivE() * pccGradE(); */ }
    
    
    // ---------------------------------------------------------------------------
  public:
    /// @name Averaging Calculus operators
    /// @{
    
    /// @return the linear operator that sums the vectorial components
    /// (x,y,z) from 3V and outputs them as scalars (x+y+z) in V
    const LinearOperator& sumV1toV() const
    { return mySumV1toV; }

    /// @return the linear operator that sums the vectorial components
    /// (x,y,z) from 3V and outputs them as three scalars
    /// (x+y+z,x+y+z,x+y+z) in 3V
    const LinearOperator& sumV1toV1() const
    { return mySumV1toV1; }

    /// @return the linear operator that sums the vectorial components
    /// (x,y,z) from 9T=3C and outputs them as scalars (x+y+z) in 3T=C
    const LinearOperator& sumC1toC() const
    { return mySumC1toC; }

    /// @return the linear operator that sums the vectorial components
    /// (x,y,z) from 9T=3C and outputs them as three scalars
    /// (x+y+z,x+y+z,x+y+z) in 9T=3C
    const LinearOperator& sumC1toC1() const
    { return mySumC1toC1; }

    /// @return the sum operator C-Vectors to Face
    const LinearOperator& sumC1toF() const
    { return mySumC1toF; }

    /// @return the sum operator from corners to faces
    const LinearOperator& sumCtoF() const
    { return mySumCtoF; }
    
    /// @return the average operator from corners to faces
    const LinearOperator& avgCtoF() const
    { return myAvgCtoF; }

    /// @return the average operator from vertices to edges
    const LinearOperator& avgVtoE() const
    { return myAvgVtoE; }

    /// @return the average operator from faces to edges
    const LinearOperator& avgFtoE() const
    { return myAvgFtoE; }

    /// @return the average operator from corners to edges
    const LinearOperator& avgCtoE() const
    { return myAvgCtoE; }
    /// @return the average operator C-Vectors -> V-Vectors
    const LinearOperator& avgC1toV1() const
    { return myAvgC1toV1; }
    
    /// @return the operator V -> C
    const LinearOperator& VtoC() const
    { return myVtoC; }
    /// @return the operator V-Vectors -> C-Vectors
    const LinearOperator& V1toC1() const
    { return myV1toC1; }
    
    /// @}

    // ------------------------- Solver services ------------------------------
  public:
    /// @name Solver services
    /// @{

    /// System of the form A x = b, where A is SDP, where you wish to
    /// have Dirichlet boundary conditions u at some places p (`p=1`or `p=0`).
    /// @pre `#row(A) = #col(A) = #col(u)`
    /// @return the linear matrix A' to prefactor.
    LinearOperator dirichletOperator( const LinearOperator& A,
                                      const Form& p ) const;

    /// System of the form A x = b, where A is SDP, where you wish to
    /// have Dirichlet boundary conditions u at some places p.
    /// @pre `#row(A) = #col(A) = #col(u)`
    /// @return the form b' to solve for
    Form dirichletVector( const LinearOperator& A,
                          const Form& b,
                          const Form& p,
                          const Form& u ) const;
    /// System of the form A x = b, where A is SDP, where you wish to
    /// have Dirichlet boundary conditions u at some places p.
    /// @pre `#row(A) = #col(A) = #col(u)`
    /// @return the solution x from the solution xp
    Form dirichletSolution( const Form& xp,
                            const Form& p,
                            const Form& u ) const;
    
    /// @}
    
    // ------------------------- Protected Datas ------------------------------
  protected:
    /// the associated mesh
    const Mesh* myMesh;

    /// Indicates if corrected calculus operators are computed.
    bool has_corrected_calculus;

    /// Indicates if discrete exterior calculus operators are computed
    /// (Caltech sense).
    bool has_discrete_exterior_calculus;

    /// Indicates if piecewise constant calculus operators are computed
    bool has_piecewise_constant_calculus;
        
    /// the number of triangulated faces
    Size myNbTriangles;
    
    /// the number of quadrangulated faces
    Size myNbQuadrangles;
    
    /// the mapping face -> nb corners (f -> 3 if surface is triangulated)
    std::vector< Size >  myNbCorners;

    /// the mapping face -> 1st corner index (f -> 3*f if surface is triangulated)
    std::vector< Corner > myFtoCIdx;
    
    /// the mapping 1st corner -> face index (c -> c/3 if surface is triangulated)
    std::vector< Face >   myCtoFIdx;
    
    // ---------------- calculus protected Datas -----------------------
  protected:
    
    /// the topological derivative operator D0: F0 -> F1
    LinearOperator myD0;
    
    /// the topological derivative operator D1: F1 -> F2
    LinearOperator myD1;
    
    // ---------------- corrected calculus protected Datas -----------------------
  protected:

    /// The inner product / mass matrix M0 for 0-forms, i.e. `<f,g>_0
    /// := f^T * M0 * g`
    LinearOperator myCCM0;
    
    /// The simplified diagonal inner product / mass matrix M0 for
    /// 0-forms, i.e. `<f,g>_0 := f^T * diag( DiagM0 ) * g`
    /// @note Lumped mass matrix, less accurate but easily inversible.
    Form myCCDiagH0;
    
    /// The inverse simplified diagonal inner product / mass matrix M0 for
    /// dual 2-forms
    Form myCCDiagDualH2;
    
    /// The inner product / mass matrix M1 for 1-forms, i.e. `<a,b>_0
    /// := a^T * M1 * b`
    LinearOperator myCCM1;

    /// The simplified Hodge star for 1-forms (use `asDiagonal()`
    /// for operator) (a negative operator).
    Form myCCDiagH1;
    
    /// The simplified Hodge star for dual 1-forms (use `asDiagonal()`
    /// for operator) (a negative operator).
    Form myCCDiagDualH1;
    
    /// The inner product / mass matrix M2 for 2-forms, a diagonal
    /// matrix, or equivalently the Hodge star for 2-forms.
    LinearOperator myCCM2;

    /// The inner product / mass matrix M2 for dual 2-forms, a
    /// diagonal matrix, or equivalently the Hodge star for dual
    /// 2-forms.
    LinearOperator myCCDualH0;

    /// The left laplacian LL0 = -1.0 * D0^T * M1 * D0;
    /// ie Lap0 * x = b <=> LL0 * x = M0 * b
    LinearOperator myCCLL0;

    /// Operator for building a vector field from a 1-form
    LinearOperator myCCSharp;
    /// Operator for building a 1-form from a vector field
    LinearOperator myCCFlat;
    /// C-Vectors inner product matrix.
    LinearOperator myCCC1IP;

    /// Sum operator C to Face
    LinearOperator mySumCtoF;
    /// Sum operator C-Vectors to Corners
    LinearOperator mySumC1toC;
    /// Sum operator C-Vectors to C-Vectors
    LinearOperator mySumC1toC1;
    /// Sum operator V-Vectors to Vertex
    LinearOperator mySumV1toV;
    /// Sum operator V-Vectors to V-Vectors
    LinearOperator mySumV1toV1;
    /// Sum operator C-Vectors to Face
    LinearOperator mySumC1toF;

    /// Average operator Corners -> Faces
    LinearOperator myAvgCtoF;
    /// Average operator Corners -> Edges
    LinearOperator myAvgCtoE;
    /// Average operator Vertices -> Edges
    LinearOperator myAvgVtoE; 
    /// Average operator Faces -> Edges
    LinearOperator myAvgFtoE; 
    /// Average operator C-Vectors -> V-Vectors
    LinearOperator myAvgC1toV1;
    /// Operator V-Vectors -> C-Vectors
    LinearOperator myV1toC1;
    /// Operator V -> C
    LinearOperator myVtoC;

    // ---------------- discrete exterior calculus protected Datas -------------------
  protected:

    /// The DEC inner product / mass matrix M0 for 0-forms, i.e. `<f,g>_0
    /// := f^T * M0 * g`, and equivalently the Hodge star for primal 0-forms
    LinearOperator myDECM0;
    /// The DEC inner product / mass matrix M1 for 1-forms, i.e. `<a,b>_1
    /// := a^T * M1 * b`, and equivalently the Hodge star for primal 1-forms
    LinearOperator myDECM1;
    /// The DEC inner product / mass matrix M2 for 2-forms, i.e. `<u,v>_2
    /// := u^T * M2 * v`, and equivalently the Hodge star for primal 2-forms
    LinearOperator myDECM2;
    /// The DEC Hodge star for dual 2-forms, and the inverse of mass matrix M0.
    LinearOperator myDECDualH2;
    /// The DEC Hodge star for dual 1-forms, and the negative of
    /// inverse of mass matrix M1.
    LinearOperator myDECDualH1;
    /// The DEC Hodge star for dual 0-forms, and the inverse of mass matrix M2.
    LinearOperator myDECDualH0;
    /// The DEC Left Laplacian for 0-forms.
    LinearOperator myDECLL0;
    /// The DEC Laplacian for 0-forms.
    LinearOperator myDECL0;

    /// The x-component constant per-face gradient operator for 0-forms
    /// (triangulated surfaces)
    LinearOperator myDECFaceGrad0x;
    /// The y-component constant per-face gradient operator for 0-forms
    /// (triangulated surfaces)
    LinearOperator myDECFaceGrad0y;
    /// The z-component constant per-face gradient operator for 0-forms
    /// (triangulated surfaces)
    LinearOperator myDECFaceGrad0z;

    // ---------------- piecewise constant calculus protected datas -------------------
  protected:

    /// The inner product / mass matrix MU for U-forms, a diagonal matrix.
    LinearOperator myPCCMU;
    /// The inner product / mass matrix MV for V-forms, a diagonal matrix.
    LinearOperator myPCCMV;
    /// The inner product / mass matrix MW for W-forms, a diagonal matrix.
    LinearOperator myPCCMW;
    /// The divergence M operator V -> U.
    LinearOperator myPCCDivM;
    /// The gradient E operator V -> W.
    LinearOperator myPCCGradE;
    /// The divergence E operator W -> V.
    LinearOperator myPCCDivE;
    /// The Laplacian E operator V -> V.
    LinearOperator myPCCLapE;
    
    // ------------------------- Protected services ------------------------------
  protected:

    void initIndexMaps();
    void initFaceMatrices();
    void initDerivativeOperators();
    void initCorrectedCalculusOperators();
    void initDiscreteExteriorCalculusOperators();
    void initPiecewiseConstantCalculusOperators( bool trivial_metric );
    
    /// Computes the inner product / mass matrix M0 for 0-forms,
    /// i.e. <f,g>_0 := f^T * M0 * g
    void computeCorrectedInnerProduct0();

    /// Computes the inner product / mass matrix M1 for 1-forms,
    /// i.e. <a,b>_0 := a^T * M1 * b
    void computeCorrectedInnerProduct1();

    /// Computes the inner product / mass matrix M2 for 2-forms, a
    /// simple diagonal matrix made of the inverse of the face
    /// corrected areas.
    void computeCorrectedInnerProduct2();

    /// Computes the sharp operator to convert a 1-form into an
    /// embedded vector field.
    void computeCorrectedSharp();

    /// Computes the flat operator to convert an
    /// embedded vector field into a 1-form (normal component is ignored)
    void computeCorrectedFlat();

    /// Computes the C-Vectors inner product matrix
    void computeCornerVectorInnerProduct();
    
    /// Computes the average operators between different cells.
    void computeAverageOperators();

    /// Adds the coefficients used to compute the inner product
    /// between 0-forms within a triangle. Takes into account the
    /// normals at vertices of the mesh.
    void triangleCorrectedInnerProduct0( Triplets& triplets, Face f );

    /// Adds the coefficients used to compute the inner product
    /// between 0-forms within a quadrangle. Takes into account the
    /// normals at vertices of the mesh.
    void quadrangleCorrectedInnerProduct0( Triplets& triplets, Face f );

    /// Adds the coefficients used to compute the inner product
    /// between 0-forms within a triangle. Takes into account the
    /// normals at vertices of the mesh.
    void triangleCorrectedInnerProduct1( Triplets& triplets, Face f );

    /// Adds the coefficients used to compute the inner product
    /// between 0-forms within a quadrangle. Takes into account the
    /// normals at vertices of the mesh.
    void quadrangleCorrectedInnerProduct1( Triplets& triplets, Face f );

    /// @param f any triangular face.
    /// @return the 3x3 matrix representing interaction coefficients between
    /// corners of face \a f when compute inner vector product.
    DenseMatrix triangleCornerVectorInnerProduct  ( const Face f ) const;

    /// @param f any quadrangular face.
    /// @return the 4x4 matrix representing interaction coefficients between
    /// corners of face \a f when compute inner vector product.
    DenseMatrix quadrangleCornerVectorInnerProduct( const Face f ) const;
    
    /// Matrix for transforming a 3D tangent vector rooted at corner c into a
    /// 1-form at its two incident edges. This matrix is invertible.
    DenseMatrix Bc( const Corner c ) const;
    
    /// Matrix for transforming a 3D tangent vector rooted at vertex i in a
    /// 1-form at edge (i,j) and (i,k). This matrix is invertible.
    DenseMatrix Bijk( Vertex i, Vertex j, Vertex k ) const;

    /// @param op any sparse linear operator.
    /// @return its average number of non zero coefficients per row.
    double nnzPerRow( const LinearOperator& op ) const;

    /// @return a dense vector containing the angle at every corners.
    DenseVector cornerAngles() const;
    
    
    // ------------------------- Internal methods ------------------------------
  private:

    /// @param i the index of the corner in the face (between 0 and 2 or 3 ).
    /// @return internal matrices for sharp operator
    DenseMatrix S( const Size i, const Size n ) const
    {
      switch( n ) {
      case 3: switch( i ) {
        case 0: return Sijk;
        case 1: return Sjki;
        default: return Skij;
        }
      case 4: switch( i ) {
        case 0: return Sijkl;
        case 1: return Sjkli;
        case 2: return Sklij;
        default: return Slijk;
        }
      default: trace.error() << "Invalid face" << std::endl;
        return DenseMatrix();
      }
    }

    /// @param i the index of the corner in the face (between 0 and 2 or 3 ).
    /// @return internal matrices for flat operator
    DenseMatrix SI( const Size i, const Size n ) const
    {
      switch( n ) {
      case 3: switch( i ) {
        case 0:  return SIijk;
        case 1:  return SIjki;
        default: return SIkij;
        }
      case 4: switch( i ) {
        case 0:  return SIijkl;
        case 1:  return SIjkli;
        case 2:  return SIklij;
        default: return SIlijk;
        }
      default: trace.error() << "Invalid face" << std::endl;
        return DenseMatrix();
      }
    }

    bool checkNaNInSparseMatrix( const LinearOperator& M, std::string name ) const
    {
      std::string str = "Check NaN in matrix " + name;
      trace.beginBlock( str.c_str() );
      Size nb_nan = 0;
      for ( int k = 0; k < M.outerSize(); ++k )
        for ( typename LinearOperator::InnerIterator it( M, k ); it; ++it )
          {
            if ( std::isnan( it.value() ) )
              {
                // trace.info() << "(" << it.row() << "," << it.col() << ")";
                nb_nan += 1;
              }
          }
      trace.info() << " --> #nan=" << nb_nan << std::endl;
      trace.endBlock();
      return nb_nan == 0;
    }

    bool checkSymmetryInSparseMatrix( const LinearOperator& M, std::string name ) const
    {
      std::string str = "Check Symmetry in matrix " + name;
      trace.beginBlock( str.c_str() );
      double max_diff = 0.0;
      for ( int k = 0; k < M.outerSize(); ++k )
        for ( typename LinearOperator::InnerIterator it( M, k ); it; ++it )
          {
            max_diff = std::max( max_diff,
                                 fabs( it.value() - M.coeff( it.col(), it.row() ) ) );
          }
      trace.info() << " --> #max_diff=" << max_diff << std::endl;
      trace.endBlock();
      return max_diff <= 1e-10;
    }
    
    // ------------------------- Internal data ------------------------------
  private:
    /// Internal matrices for computing inner product on 1-forms (triangular faces).
    DenseMatrix Sijk, Sjki, Skij;
    /// Internal matrices for flat operator on 1-vectors (triangular faces).
    DenseMatrix SIijk, SIjki, SIkij;

    /// Internal matrices for computing inner product on 1-forms (quadrangular faces).
    DenseMatrix Sijkl, Sjkli, Sklij, Slijk;
    /// Internal matrices for flat operator on 1-vectors (quadrangular faces).
    DenseMatrix SIijkl, SIjkli, SIklij, SIlijk;
    
  }; // struct SurfaceMeshDEC

  
} // namespace DGtal

///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "SurfaceMeshDEC.ih"
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SurfaceMeshDEC_h

#undef SurfaceMeshDEC_RECURSES
#endif // else defined(SurfaceMeshDEC_RECURSES)

