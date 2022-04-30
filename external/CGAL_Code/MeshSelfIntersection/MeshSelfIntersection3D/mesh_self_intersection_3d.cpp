/* =============================================================================================
 *
 * 	mesh_self_intersection_3d.cpp
 *
 * 	Determines if a 3D polyhedral mesh contains any self-intersections and returns
 * 	the pairs of intersecting faces.
 *
 * 	by Dillon Cislo
 * 	01/23/2019
 *
 * 	This is a MEX-file for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <vector>
#include <stdexcept>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 	Kernel;
typedef Kernel::Point_3 					Point;
typedef CGAL::Surface_mesh<Point>	 			Mesh;

typedef boost::graph_traits<Mesh>::face_descriptor 		face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor 		vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

///
/// Brief main function to call computational functionalities
///
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------

	// Check for proper number of arguments
	if ( nrhs != 2 ) {
		mexErrMsgIdAndTxt("MATLAB:mesh_self_intersection_3d:nargin",
				"MESH_SELF_INTERSECTION_3D requires two input arguments.");
	} else if ( nlhs != 2 ) {
		mexErrMsgIdAndTxt("MATLAB:mesh_self_intersection_3d:nargout",
				"MESH_SELF_INTERSECTION_3D requires two output arguments.");
	}

	double *faces = mxGetPr( prhs[0] );	 // The face connectivity list
	int numFaces = (int) mxGetM( prhs[0] );	 // The number of faces
	int sizeFaces = (int) mxGetN( prhs[0] ); // The number of vertices in a single face

	double *vertex = mxGetPr( prhs[1] ); 	 // The vertex coordinate list
	int numVertex = (int) mxGetM( prhs[1] ); // The number of vertices
	int dim = (int) mxGetN( prhs[1] ); 	 // The dimensions of the vertex coordinates

	// Check dimensionality of the vertex list
	if ( dim != 3 ) {
		mexErrMsgIdAndTxt("MATLAB:mesh_self_intersection_3d:dimension",
				"Vertex array improperly sized");
	}

	// Create and populate the polyhedral mesh ---------------------------------------------
	
	// Create vector of 3D point objects
	std::vector<Point> points;
	points.reserve( numVertex );
	for( int i = 0; i < numVertex; i++ ) {

		points.push_back( Point( vertex[i],
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
	PMP::orient_polygon_soup( points, polygons );
	PMP::polygon_soup_to_polygon_mesh( points, polygons, mesh );
	
	// -------------------------------------------------------------------------------------
	// SELF-INTERSECTION CHECK
	// -------------------------------------------------------------------------------------
	
	plhs[0] = mxCreateLogicalMatrix(1,1);
	bool *intersecting = mxGetLogicals( plhs[0] ); // Output

	*intersecting = PMP::does_self_intersect( mesh,
		PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)) );

	std::vector< std::pair<face_descriptor, face_descriptor> > intersected_polys;
	PMP::self_intersections( mesh, std::back_inserter(intersected_polys) );

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
	

	int numIntersect = intersected_polys.size();

	if ( numIntersect != 0 ) {

		plhs[1] = mxCreateDoubleMatrix( numIntersect, 2, mxREAL );
		double *int_poly_out = mxGetPr( plhs[1] );
	
		for( int i = 0; i < numIntersect; i++ ) {
	
			// NOTE: We add 1 to the index to account for
			// MATLAB's 1-indexed array structures
			int_poly_out[i] = (double) intersected_polys[i].first+1.0;
			int_poly_out[i+numIntersect] = (double) intersected_polys[i].second+1.0;
	
		}

	} else {

		plhs[1] = mxCreateDoubleMatrix( 0, 0, mxREAL );

	}

	return;

};
