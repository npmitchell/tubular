/* =============================================================================================
 *
 *  surface_parameterization.cpp
 *
 *  Find a planar parameterization of a 3D triangular mesh surface representation.
 *  Right now the code only supports parameterization of topological disks.
 *  Depends on CGAL.
 *
 *  by Dillon Cislo
 *  05/11/2019
 *
 *  This is a MEX-File for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <cstdlib>
#include <vector>
#include <iostream>

#include <CGAL/Simple_cartesian.h> 
#include <CGAL/Surface_mesh.h>

// Basic surface parameterization objects
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/MVC_post_processor_3.h>

// Fixed boundary parameterizers
#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>

// Free boundary parameterizers
#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>

// No-boundary parameterizers (topological sphere) **NOT CURRENTLY SUPPORTED**
#include <CGAL/Surface_mesh_parameterization/Orbifold_Tutte_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/orbifold_enums.h>
#include <CGAL/Surface_mesh_parameterization/orbifold_shortest_path.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <boost/function_output_iterator.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                   Kernel;
typedef Kernel::Point_2                                  Point_2;
typedef Kernel::Point_3                                  Point_3;

typedef CGAL::Surface_mesh<Point_3>                      Mesh;
typedef Mesh::Vertex_index                               Vertex_index;
typedef Mesh::Face_index                                 Face_index;

typedef boost::graph_traits<Mesh>::halfedge_descriptor   halfedge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor       face_descriptor;
typedef Mesh::Property_map<vertex_descriptor, Point_2>   UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;

/// Enumeration for type of boundary
enum BORDER_TYPE { FIXED = 1, FREE = 2 };

// Enumeration for the type of fixed boundary parameterization
enum FIXED_BORDER_TYPE { ARC_LENGTH = 1, UNIFORM = 2 };

// Enumeration for the type of fixed boundary target shape
enum FIXED_BORDER_SHAPE { CIRCLE = 1, SQUARE = 2 };

// Enumeration for the fixed boundary parameterization method
enum FIXED_BORDER_METHOD { BARYCENTRIC = 1, AUTHALIC = 2, CONFORMAL = 3, MEAN = 4 };

// Enumeration for the free boundary parameterization method
enum FREE_BORDER_METHOD { LSCM = 1, ARAP = 2 };

///
/// Parameterize a surface using a fixed boundary method
///
void parameterize_fixed_border( SMP::Error_code &err, Mesh &mesh, UV_pmap &uvMap, 
    halfedge_descriptor &bhd, int fixedType, int fixedShape, int paramMethod,
    bool cS, const std::vector<vertex_descriptor> &corners ) {

  // This is so ugly - I hate that I have to do this, but none of the
  // parameterizers are templated....
  typedef SMP::Circular_border_arc_length_parameterizer_3<Mesh> CBALP;
  typedef SMP::Circular_border_uniform_parameterizer_3<Mesh> CBUP;
  typedef SMP::Square_border_arc_length_parameterizer_3<Mesh> SBALP;
  typedef SMP::Square_border_uniform_parameterizer_3<Mesh> SBUP;

  typedef SMP::Barycentric_mapping_parameterizer_3<Mesh, CBALP> B_CBALP;
  typedef SMP::Barycentric_mapping_parameterizer_3<Mesh, CBUP> B_CBUP;
  typedef SMP::Barycentric_mapping_parameterizer_3<Mesh, SBALP> B_SBALP;
  typedef SMP::Barycentric_mapping_parameterizer_3<Mesh, SBUP> B_SBUP;

  typedef SMP::Discrete_authalic_parameterizer_3<Mesh, CBALP> A_CBALP;
  typedef SMP::Discrete_authalic_parameterizer_3<Mesh, CBUP> A_CBUP;
  typedef SMP::Discrete_authalic_parameterizer_3<Mesh, SBALP> A_SBALP;
  typedef SMP::Discrete_authalic_parameterizer_3<Mesh, SBUP> A_SBUP;

  typedef SMP::Discrete_conformal_map_parameterizer_3<Mesh, CBALP> C_CBALP;
  typedef SMP::Discrete_conformal_map_parameterizer_3<Mesh, CBUP> C_CBUP;
  typedef SMP::Discrete_conformal_map_parameterizer_3<Mesh, SBALP> C_SBALP;
  typedef SMP::Discrete_conformal_map_parameterizer_3<Mesh, SBUP> C_SBUP;

  typedef SMP::Mean_value_coordinates_parameterizer_3<Mesh, CBALP> M_CBALP;
  typedef SMP::Mean_value_coordinates_parameterizer_3<Mesh, CBUP> M_CBUP;
  typedef SMP::Mean_value_coordinates_parameterizer_3<Mesh, SBALP> M_SBALP;
  typedef SMP::Mean_value_coordinates_parameterizer_3<Mesh, SBUP> M_SBUP;

  // Run the paramterization method
  switch( paramMethod ) {

    case BARYCENTRIC :

      if ( fixedType == ARC_LENGTH ) {
        if ( fixedShape == CIRCLE ) {
          err = SMP::parameterize( mesh, B_CBALP(), bhd, uvMap );
        } else {
          if ( cS ) {
            SBALP bp( corners[0], corners[1], corners[2], corners[3] );
            err = SMP::parameterize( mesh, B_SBALP( bp ), bhd, uvMap );
          } else {
            err = SMP::parameterize( mesh, B_SBALP(), bhd, uvMap );
          }
        }
      } else {
        if ( fixedShape == CIRCLE ) {
          err = SMP::parameterize( mesh, B_CBUP(), bhd, uvMap );
        } else {
          if ( cS ) {
            SBUP bp( corners[0], corners[1], corners[2], corners[3] );
            err = SMP::parameterize( mesh, B_SBUP( bp ), bhd, uvMap );
          } else {
            err = SMP::parameterize( mesh, B_SBUP(), bhd, uvMap );
          }
        }
      }

      break;

    case AUTHALIC :

      if ( fixedType == ARC_LENGTH ) {
        if ( fixedShape == CIRCLE ) {
          err = SMP::parameterize( mesh, A_CBALP(), bhd, uvMap );
        } else {
          if ( cS ) {
            SBALP bp( corners[0], corners[1], corners[2], corners[3] );
            err = SMP::parameterize( mesh, A_SBALP( bp ), bhd, uvMap );
          } else {
            err = SMP::parameterize( mesh, A_SBALP(), bhd, uvMap );
          }
        }
      } else {
        if ( fixedShape == CIRCLE ) {
          err = SMP::parameterize( mesh, A_CBUP(), bhd, uvMap );
        } else {
          if ( cS ) {
            SBUP bp( corners[0], corners[1], corners[2], corners[3] );
            err = SMP::parameterize( mesh, A_SBUP( bp ), bhd, uvMap );
          } else {
            err = SMP::parameterize( mesh, A_SBUP(), bhd, uvMap );
          }
        }
      }

      break;

    case CONFORMAL :

      if ( fixedType == ARC_LENGTH ) {
        if ( fixedShape == CIRCLE ) {
          err = SMP::parameterize( mesh, C_CBALP(), bhd, uvMap );
        } else {
          if ( cS ) {
            SBALP bp( corners[0], corners[1], corners[2], corners[3] );
            err = SMP::parameterize( mesh, C_SBALP( bp ), bhd, uvMap );
          } else {
            err = SMP::parameterize( mesh, C_SBALP(), bhd, uvMap );
          }
        }
      } else {
        if ( fixedShape == CIRCLE ) {
          err = SMP::parameterize( mesh, C_CBUP(), bhd, uvMap );
        } else {
          if ( cS ) {
            SBUP bp( corners[0], corners[1], corners[2], corners[3] );
            err = SMP::parameterize( mesh, C_SBUP( bp ), bhd, uvMap );
          } else {
            err = SMP::parameterize( mesh, C_SBUP(), bhd, uvMap );
          }
        }
      }

      break;

    case MEAN :

      if ( fixedType == ARC_LENGTH ) {
        if ( fixedShape == CIRCLE ) {
          err = SMP::parameterize( mesh, M_CBALP(), bhd, uvMap );
        } else {
          if ( cS ) {
            SBALP bp( corners[0], corners[1], corners[2], corners[3] );
            err = SMP::parameterize( mesh, M_SBALP( bp ), bhd, uvMap );
          } else {
            err = SMP::parameterize( mesh, M_SBALP(), bhd, uvMap );
          }
        }
      } else {
        if ( fixedShape == CIRCLE ) {
          err = SMP::parameterize( mesh, M_CBUP(), bhd, uvMap );
        } else {
          if ( cS ) {
            SBUP bp( corners[0], corners[1], corners[2], corners[3] );
            err = SMP::parameterize( mesh, M_SBUP( bp ), bhd, uvMap );
          } else {
            err = SMP::parameterize( mesh, M_SBUP(), bhd, uvMap );
          }
        }
      }

      break;

  }

};

///
/// Parameterize a surface using a free boundary method
///
void parameterize_free_border( SMP::Error_code &err, Mesh &mesh, UV_pmap &uvMap, 
    halfedge_descriptor &bhd, int paramMethod,
    bool xS, const std::vector<vertex_descriptor> &xPts ) {

  typedef SMP::Two_vertices_parameterizer_3<Mesh> TVP;

  typedef SMP::LSCM_parameterizer_3<Mesh, TVP> LSCMP;
  typedef SMP::ARAP_parameterizer_3<Mesh, TVP> ARAPP;

  switch( paramMethod ) {

    case LSCM :

      if ( xS ) {
        TVP bp( xPts[0], xPts[1] );
        err = SMP::parameterize( mesh, LSCMP( bp ), bhd, uvMap );
      } else {
        err = SMP::parameterize( mesh, LSCMP(), bhd, uvMap );
      }

      break;

    case ARAP :

      if ( xS ) {
        TVP bp( xPts[0], xPts[1] );
        err = SMP::parameterize( mesh, ARAPP( bp ), bhd, uvMap );
      } else {
        err = SMP::parameterize( mesh, ARAPP(), bhd, uvMap );
      }

      break;

  }

};


///
/// Brief main function to call computational functionalities
///
void mexFunction( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] ) {

  // -------------------------------------------------------------------------------------------
  // INPUT PROCESSING
  // -------------------------------------------------------------------------------------------
  
  // Check for proper number of arguments
  if ( nrhs != 3 ) {
    mexErrMsgIdAndTxt( "MATLAB:surface_parameterization:nargin",
        "SURFACE_PARMETERIZATION requires three input arguments." );
  } else if ( nlhs != 2 ) {
    mexErrMsgIdAndTxt( "MATLAB:surface_parameterization:nargout",
        "SURFACE_PARMETERIZATION requres two output arguments." );
  }

  double *faces = mxGetPr( prhs[0] );         // The face connectivity list
  std::size_t numFaces = mxGetM( prhs[0] );   // The number of faces
  std::size_t sizeFaces = mxGetN( prhs[0] );  // The number of vertices in a single face


  // Check that the input mesh is a triangulation
	if ( sizeFaces != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_parameterization:face_dim",
				"Input mesh must be a triangulation." );
	}


	double *vertex = mxGetPr( prhs[1] ); 	      // The vertex coordinate list
	std::size_t numVertex = mxGetM( prhs[1] );  // The number of vertices
	std::size_t dim = mxGetN( prhs[1] ); 	      // The dimensionality of the vertex list

	// Check the dimensionality of the vertex list
	if ( dim != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_parameterization:vertex_dim",
				"Vertex coordinates must be 3D." );
	}

  // Default parameterization settings
  int idx;

  int borderType = FIXED;
  int fixedType = ARC_LENGTH;
  int fixedShape = CIRCLE;
  int paramMethod = CONFORMAL;
  
  bool cornersSupplied = false;
  std::vector<vertex_descriptor> corners;
  corners.reserve( 4 );
  double *cornersIn;

  bool xPtsSupplied = false;
  std::vector<vertex_descriptor> xPts;
  xPts.reserve( 2 );
  double *xPtsIn;

  if ( ( idx = mxGetFieldNumber( prhs[2], "borderType" ) ) != -1 ) {
    borderType = (int) *mxGetPr(mxGetFieldByNumber( prhs[2], 0, idx ));
  }

  if ( ( idx = mxGetFieldNumber( prhs[2], "fixedType" ) ) != -1 ) {
    fixedType = (int) *mxGetPr(mxGetFieldByNumber( prhs[2], 0, idx ));
  }

  if ( ( idx = mxGetFieldNumber( prhs[2], "fixedShape" ) ) != -1 ) {
    fixedShape = (int) *mxGetPr(mxGetFieldByNumber( prhs[2], 0, idx ));
  }

  if ( ( idx = mxGetFieldNumber( prhs[2], "paramMethod" ) ) != -1 ) {
    paramMethod = (int) *mxGetPr(mxGetFieldByNumber( prhs[2], 0, idx ));
  }

  if ( ( idx = mxGetFieldNumber( prhs[2], "corners" ) ) != -1 ) {
    cornersSupplied = true;
    cornersIn = mxGetPr(mxGetFieldByNumber( prhs[2], 0, idx ));
  }

  if ( ( idx = mxGetFieldNumber( prhs[2], "fixedPoints" ) ) != -1 ) {
    xPtsSupplied = true;
    xPtsIn = mxGetPr(mxGetFieldByNumber( prhs[2], 0, idx ));
  }

  if ( cornersSupplied ) {
    for( int i = 0; i < 4; i++ ) {
      corners.push_back( (vertex_descriptor) (cornersIn[i]-1) );
    }
  }

  if ( xPtsSupplied ) {
    for( int i = 0; i < 2; i++ ) {
      xPts.push_back( (vertex_descriptor) (xPtsIn[i]-1) );
    }
  }

	// Create and populate the polyhedral mesh ---------------------------------------------
	
	// Create vector of 3D point objects
	std::vector<Point_3> points;
	points.reserve( numVertex );
	for( int i = 0; i < numVertex; i++ ) {

		points.push_back( Point_3( vertex[i],
					 vertex[i+numVertex],
					 vertex[i+(2*numVertex)] ) );

	}

	// Create vector of polygon objects
	std::vector< std::vector<std::size_t> > polygons;
	polygons.reserve( numFaces );
	for( int i = 0; i < numFaces; i++ ) {

		std::vector<std::size_t> currentPolygon;
		currentPolygon.reserve( sizeFaces );
		for( int j = 0; j < sizeFaces; j++ ) {

			// NOTE: We subtract 1 from the index to account for
			// MATLAB's 1-indexed array structures
			double index = faces[i+(j*numFaces)]-1.0;
			currentPolygon.push_back( (std::size_t) index );

		}

		polygons.push_back( currentPolygon );

	}

	// Populate the mesh
	Mesh mesh;
  CGAL::Polygon_mesh_processing::orient_polygon_soup( points, polygons );
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, mesh );

  // -------------------------------------------------------------------------------------------
  // MESH PROCESSING
  // -------------------------------------------------------------------------------------------
  
  // A halfedge on the boundary
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border( mesh ).first;

  // The 2D points of the (uv)-parameterization will be written to this property map
  UV_pmap uvMap = mesh.add_property_map<vertex_descriptor, Point_2>
    ("v:uv", Point_2(0,0)).first;

  // Parameterize surface mesh
  SMP::Error_code err;

  if ( borderType == FIXED ) {

    parameterize_fixed_border( err, mesh, uvMap, bhd,
        fixedType, fixedShape, paramMethod,
        cornersSupplied, corners );

  } else if ( borderType == FREE ) {

    parameterize_free_border( err, mesh, uvMap, bhd, paramMethod,
        xPtsSupplied, xPts );

  } else {

    mexErrMsgTxt("Invalid border type provided!");

  }

  if ( mesh.has_garbage() ) { mesh.collect_garbage(); }

  // -------------------------------------------------------------------------------------------
  // OUTPUT PROCESSING
  // -------------------------------------------------------------------------------------------
  
  std::size_t numFacesFinal = mesh.number_of_faces(); // Final number of faces
  std::size_t numVertexFinal = mesh.number_of_vertices(); // Final number of vertices

  plhs[0] = mxCreateDoubleMatrix( numVertexFinal, 2, mxREAL );
  double *uvOut = mxGetPr( plhs[0] );
  
  plhs[1] = mxCreateDoubleMatrix( numFacesFinal, sizeFaces, mxREAL );
  double *facesOut = mxGetPr( plhs[1] );

  // Collect (u,v)-coordinates
  int i = 0;
  BOOST_FOREACH( Vertex_index v, mesh.vertices() ) {

    Point_2 uv = uvMap[v];

    uvOut[i] = (double) uv.x();
    uvOut[i+numVertexFinal] = (double) uv.y();

    i++;

  }

	// Collect face quantities
	i = 0;
	BOOST_FOREACH( Face_index f, mesh.faces() ) {

		// Iterate around the current face to find connectivity
		int j = 0;
		BOOST_FOREACH( Vertex_index v, vertices_around_face(mesh.halfedge(f), mesh) ) {

			// NOTE: We add 1 to the index to account for
			// MATLAB's 1-indexed array structures
			facesOut[i+(j*numFacesFinal)] = (double) v+1.0;
			j++;

		}

		i++;

	}

  return;

};

