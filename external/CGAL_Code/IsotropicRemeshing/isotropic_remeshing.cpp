/* ==============================================================================================
 *
 * 	isotropic_remeshing.cpp
 *
 * 	Performs an isotropic remeshing of a region of a polygon mesh.  Depends on CGAL
 *
 *  Example Usage 
 *  -------------
 *  [Fnew, Vnew, Fnnew, Vnnew] = ...
 *     isotropic_remeshing(faces, vertices, targetEdgeLen, numIterations)
 * 
 *
 * 	by Dillon Cislo
 * 	02/19/2019
 *
 * 	This is a MEX-File for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <vector>
#include <stdexcept>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <boost/function_output_iterator.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 	Kernel;
typedef Kernel::Point_3 					Point;
typedef Kernel::Vector_3					Vector;

typedef CGAL::Surface_mesh<Point> 				Mesh;
typedef Mesh::Vertex_index 					Vertex_index;
typedef Mesh::Face_index 					Face_index;

typedef boost::graph_traits<Mesh>::face_descriptor 		face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor 		vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor 		halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor 		edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

///
/// A struct holding a vector of edges associated with certain halfedges
///
struct halfedge2edge {

	halfedge2edge( const Mesh &m, std::vector<edge_descriptor> &edges )
		: m_mesh( m ), m_edges( edges ) {}

	void operator()( const halfedge_descriptor &h ) const {
		m_edges.push_back( edge( h, m_mesh ) );
	}

	const Mesh &m_mesh;
	
	std::vector<edge_descriptor> &m_edges;

};

///
/// Brief main function to call computational functionalities
///
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	// Check for proper number of arguments
	if ( nrhs != 4 ) {
		mexErrMsgIdAndTxt( "MATLAB:isotropic_remeshing:nargin",
				"ISOTROPIC_REMESHING requires four input arguments." );
	} else if ( nlhs != 4 ) {
		mexErrMsgIdAndTxt( "MATLAB:isotropic_remeshing:nargout",
				"ISOTROPIC_REMESHING requires four output arguments." );
	}

	double *faces = mxGetPr( prhs[0] );	   // The face connectivity list
	std::size_t numFaces = mxGetM( prhs[0] );  // The number of faces
	std::size_t sizeFaces = mxGetN( prhs[0] ); // The number of vertices in a single face

	double *vertex = mxGetPr( prhs[1] ); 	   // The vertex coordinate list
	std::size_t numVertex = mxGetM( prhs[1] ); // The number of vertices
	std::size_t dim = mxGetN( prhs[1] ); 	   // The dimensionality of the vertex list

	// Check the dimensionality of the vertex list
	if ( dim != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:isotropic_remeshing:vertex_dim",
				"Vertex coordinates must be 3D." );
	}

	// The target edge length of the final output mesh
	double target_edge_length = *mxGetPr( prhs[2] );

	// Check the sign of the target edge length
	if ( target_edge_length <= 0.0 ) {
		mexErrMsgIdAndTxt( "MATLAB:isotropic_remeshing:edge_sign",
				"Target edge length must be positive." );
	}

	// The number of iterations to be run
	int nb_iter = (int) *mxGetPr( prhs[3] );

	// Check the sign of the number of iterations
	if ( nb_iter < 1 ) {
		mexErrMsgIdAndTxt( "MATLAB:isotropic_remeshing:iter_sign",
				"Number of iterations must be positive." );
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
	// MESH PROCESSING
	// -------------------------------------------------------------------------------------
	
	// Split border edges
	std::cout << "Split border... ";

	std::vector<edge_descriptor> border;
	PMP::border_halfedges( mesh.faces(), mesh,
			boost::make_function_output_iterator(halfedge2edge(mesh,border)) );
	PMP::split_long_edges( border, target_edge_length, mesh );

	std::cout << "done." << std::endl;

	// Perform isotropic remeshing
	std::cout << "Start remeshing... ";

	PMP::isotropic_remeshing(
			mesh.faces(),
			target_edge_length,
			mesh,
			PMP::parameters::number_of_iterations(nb_iter)
			.protect_constraints(true) // i.e. protect border
			);

	std::cout << "done." << std::endl;

	// Collect any garbage that may have accumulated in the mesh
	if ( mesh.has_garbage() ) { mesh.collect_garbage(); }

	// Calculate normals
	auto fnormals = mesh.add_property_map<face_descriptor, Vector>
		("f:normals", CGAL::NULL_VECTOR).first;
	auto vnormals = mesh.add_property_map<vertex_descriptor, Vector>
		("v:normals", CGAL::NULL_VECTOR).first;

	PMP::compute_normals( mesh, vnormals, fnormals,
			PMP::parameters::vertex_point_map( mesh.points() ).
			geom_traits( Kernel() ) );
	

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	std::size_t numFacesFinal = mesh.number_of_faces(); // Final number of faces
	std::size_t numVertexFinal = mesh.number_of_vertices(); // Final number of vertices

	// I'm assuming the remeshing process keeps face size the same?
	plhs[0] = mxCreateDoubleMatrix( numFacesFinal, sizeFaces, mxREAL );
	double *facesOut = mxGetPr( plhs[0] );

	plhs[1] = mxCreateDoubleMatrix( numVertexFinal, 3, mxREAL );
	double *vertexOut = mxGetPr( plhs[1] );

	plhs[2] = mxCreateDoubleMatrix( numFacesFinal, 3, mxREAL );
	double *fnOut = mxGetPr( plhs[2] );

	plhs[3] = mxCreateDoubleMatrix( numVertexFinal, 3, mxREAL );
	double *vnOut = mxGetPr( plhs[3] );

	// Collect face quantities
	int i = 0;
	BOOST_FOREACH( Face_index f, mesh.faces() ) {

		// Iterate around the current face to find connectivity
		int j = 0;
		BOOST_FOREACH( Vertex_index v, vertices_around_face(mesh.halfedge(f),mesh) ) {

			// NOTE: We add 1 to the index to account for
			// MATLAB's 1-indexed array structures
			facesOut[i+(j*numFacesFinal)] = (double) v+1.0;
			j++;

		}

		// Face normal vector
		Vector fn = fnormals[f];

		fnOut[i] = fn[0];
		fnOut[i+numFacesFinal] = fn[1];
		fnOut[i+(2*numFacesFinal)] = fn[2];
		
		i++;

	}

	// Collect vertex quantities
	i = 0;
	BOOST_FOREACH( Vertex_index v, mesh.vertices() ) {

		// Vertex coordinates
		Point pp = mesh.point( v );

		vertexOut[i] = pp[0];
		vertexOut[i+numVertexFinal] = pp[1];
		vertexOut[i+(2*numVertexFinal)] = pp[2];

		// Vertex normal vector
		Vector vn = vnormals[v];

		vnOut[i] = vn[0];
		vnOut[i+numVertexFinal] = vn[1];
		vnOut[i+(2*numVertexFinal)] = vn[2];

		i++;

	}

	return;

};


