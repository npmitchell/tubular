/* =============================================================================================
 *
 *  surface_reconstruction.cpp
 *  
 *  This function creates a mesh triangulation from a disordered 3D point cloud using
 *  either the "Advancing front" surface reconstruction method or the "Scale space"
 *  surface reconstruction method. Depends on CGAL.
 *
 *  by Dillon Cislo
 *  12/08/2021
 *
 *  This is a MEX-file for MATLAB
 *  
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/remove_outliers.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/jet_smooth_point_set.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Scale_space_reconstruction_3/Jet_smoother.h>
#include <CGAL/Scale_space_reconstruction_3/Advancing_front_mesher.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
//#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                         Point_3;
typedef Kernel::Vector_3                        Vector_3;
typedef Kernel::Sphere_3                        Sphere_3;
typedef CGAL::Point_set_3<Point_3, Vector_3>    Point_set;
typedef std::array<std::size_t, 3>              Facet;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

enum RECONSTRUCTION_TYPE {

  // Advancing front reconstruction
  ADVANCING_FRONT = 0,

  // Scale space reconstruction
  SCALE_SPACE = 1

};

// Main function
void mexFunction( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------

  // Check for proper number of arguments
  if ( nrhs != 2 ) {
    mexErrMsgIdAndTxt( "MATLAB:surface_reconstruction:nargin",
        "SURFACE_RECONSTRUCTION requires 2 input arguments" );
  } else if ( nlhs != 2 ) {
    mexErrMsgIdAndTxt("MATLAB:surface_reconstruction:nargout",
        "SURFACE_RECONSTRUCTION requries 2 output arguments" );
  }

	// The input point cloud list
	double *pts = mxGetPr( prhs[0] );
	std::size_t numPoints = mxGetM( prhs[0] ); // The number of points
	std::size_t dim = mxGetN( prhs[0] ); // The dimensionality of the vertex list

	// Check the dimensionality of the point cloud
	if ( dim != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_reconstruction:point_dim",
				"Point coordinates must be 3D" );
	}

  // Load input points into point set structure 
  Point_set points;
  for( int i = 0; i < numPoints; i++ ) {

		Point_3 pp = Point_3( pts[i], pts[i+numPoints], pts[i+(2*numPoints)] );
    points.insert( pp );

  }


  // Process input options --------------------------------------------------------------
  
  int idx, tmp;
  
  // The method used for surface reconstruction
  int reconstruction_choice = ADVANCING_FRONT;
  if ( (idx = mxGetFieldNumber( prhs[1], "reconstructionMethod" )) != -1 ) {

    tmp = (int) *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));

    if (tmp == 0) {
      reconstruction_choice = ADVANCING_FRONT;
    } else if (tmp == 1) {
      reconstruction_choice = SCALE_SPACE;
    } else {
      mexErrMsgTxt("Invalid surface reconstruction method supplied!");
    }

  }

  // Whether or not to remove outliers from the point set
  bool remove_outliers = false; 
  if ( (idx = mxGetFieldNumber( prhs[1], "removeOutliers" )) != -1 ) {
    remove_outliers = *mxGetLogicals(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  // The number of neighbors considered in the outlier removal process
  int ro_num_neighbors = 24;
  if ( (idx = mxGetFieldNumber( prhs[1], "numNeighborsRO" )) != -1 ) {
    ro_num_neighbors = (int) *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  if (ro_num_neighbors <= 0) {
    mexErrMsgTxt("Number of neighbors for outlier removal must be positive");
  }

  // The percentage of points to remove in the outlier removal process
  double ro_threshold_percent = 5.0;
  if ( (idx = mxGetFieldNumber( prhs[1], "thresholdPercentRO" )) != -1 ) {
    ro_threshold_percent = *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  if ((ro_threshold_percent <= 0.0) || (ro_threshold_percent > 100.0)) {
    mexErrMsgTxt("Threshold percentage for outlier removal is invalid");
  }

  // The distance scale for outlier removal
  double ro_dist_scale = 2.0;
  if ( (idx = mxGetFieldNumber( prhs[1], "thresholdDistanceRO" )) != -1 ) {
    ro_dist_scale = *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  if (ro_dist_scale <= 0.0) {
    mexErrMsgTxt("Threshold distance for outlier removal is invalid");
  }

  // Whether or not to appyl grid simplification to the point set
  bool simplify_point_set = false;
  if ( (idx = mxGetFieldNumber( prhs[1], "simplifyPointSet" )) != -1 ) {
    simplify_point_set = *mxGetLogicals(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  // The number of neighbors used to determine average point spacing
  int avg_spacing_num_neighbors = 6;
  if ( (idx = mxGetFieldNumber( prhs[1], "numNeighborsSpacing" )) != -1 ) {
    avg_spacing_num_neighbors = (int) *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  if (avg_spacing_num_neighbors <= 0) {
    mexErrMsgTxt("Number of neighbors for determining average spacing must be positive");
  }

  // Compute average spacing using a user specified neighborhood
  double spacing = CGAL::compute_average_spacing<Concurrency_tag>(
      points, avg_spacing_num_neighbors );

  // The scale of the average spacing used to construct the simplification grid
  double sips_spacing_scale = 2.0;
  if ( (idx = mxGetFieldNumber( prhs[1], "spacingScaleSIPS" )) != -1 ) {
    sips_spacing_scale = *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  if (sips_spacing_scale <= 0.0) {
    mexErrMsgTxt("Spacing scale for point set simplification must be positive");
  }

  // Whether or not to smooth the point set
  bool smooth_point_set = false;
  if ( (idx = mxGetFieldNumber( prhs[1], "smoothPointSet" )) != -1 ) {
    smooth_point_set = *mxGetLogicals(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  // The number of neighbors considered in the point smoothing process
  int smps_num_neighbors = 24;
  if ( (idx = mxGetFieldNumber( prhs[1], "numNeighborsSMPS" )) != -1 ) {
    smps_num_neighbors = (int) *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  if (smps_num_neighbors <= 0) {
    mexErrMsgTxt("Number of neighbors for point set smoothing must be positive");
  }

  // Number of jet smoothing iterations applied during scale space reconstruction
  int ss_smooth_iter = 4;
  if ( (idx = mxGetFieldNumber( prhs[1], "smoothIterations" )) != -1 ) {
    ss_smooth_iter = (int) *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  if (ss_smooth_iter <= 0) {
    mexErrMsgTxt("Number of smoothing iterations must be positive");
  }

  // Maximum facet length allowd during scale space reconstruction
  double ss_max_length = std::numeric_limits<double>::infinity();
  if ( (idx = mxGetFieldNumber( prhs[1], "maxFacetLength" )) != -1 ) {
    ss_max_length = *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
  }

  if (ss_max_length <= 0.0) {
    mexErrMsgTxt("Maximum facet length must be positive");
  }


  // ------------------------------------------------------------------------------------
  // POINT SET PROCESSING
  // ------------------------------------------------------------------------------------
  

  // Remove outliers from point set
  if ( remove_outliers ) {

    CGAL::remove_outliers( points,
        ro_num_neighbors,
        points.parameters().threshold_percent( ro_threshold_percent ).
        threshold_distance( ro_dist_scale * spacing ) );

    points.collect_garbage();

  }

  // Simplify point set
  if ( simplify_point_set ) {

    // Simplify using a grid of a use specified scale
    CGAL::grid_simplify_point_set( points, sips_spacing_scale * spacing );

    points.collect_garbage();

  }

  // Smooth the surface
  if ( smooth_point_set ) {
    CGAL::jet_smooth_point_set<Concurrency_tag>( points, smps_num_neighbors);
  }

  // ------------------------------------------------------------------------------------
  // SURFACE RECONSTRUCTION
  // ------------------------------------------------------------------------------------
  
  
  std::vector<Facet> facets; // Output face connectivity list
  std::vector<Point_3> vertices; // Output vertex coordinate list

  if (reconstruction_choice == ADVANCING_FRONT) {

    // The function is called directly using the points raw iterators
    CGAL::advancing_front_surface_reconstruction( points.points().begin(),
        points.points().end(), std::back_inserter(facets) );

    // Copy points for random access
    vertices.reserve( points.size() );
    std::copy( points.points().begin(), points.points().end(),
        std::back_inserter(vertices) );


  } else if (reconstruction_choice == SCALE_SPACE) {

    CGAL::Scale_space_surface_reconstruction_3<Kernel> reconstruct(
        points.points().begin(), points.points().end() );

    // Iteratively apply jet smoothing
    reconstruct.increase_scale( ss_smooth_iter,
        CGAL::Scale_space_reconstruction_3::Jet_smoother<Kernel>() );

    // Mesh with Advancing Front mesher with a user specified
    // maximum face length
    reconstruct.reconstruct_surface(
        CGAL::Scale_space_reconstruction_3::Advancing_front_mesher<Kernel>(
          ss_max_length) );

    // Copy points for output
    vertices.reserve( points.size() );
    std::copy( points.points().begin(), points.points().end(),
        std::back_inserter(vertices) );

    // Copy faces for output
    facets.reserve( reconstruct.number_of_facets() );
    for (const auto &facet : CGAL::make_range( reconstruct.facets_begin(),
          reconstruct.facets_end() ) ) {

      Facet f = {facet[0], facet[1], facet[2]};
      facets.push_back( f );

    }

  } else {

    mexErrMsgTxt("Invalid reconstruction method");

  }

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
  
  std::size_t numVertices = vertices.size();
  std::size_t numFaces = facets.size();

  plhs[0] = mxCreateDoubleMatrix( numFaces, 3, mxREAL );
  double *facesOut = mxGetPr( plhs[0] );

  plhs[1] = mxCreateDoubleMatrix( numVertices, dim, mxREAL );
  double *vertexOut = mxGetPr( plhs[1] );

  // NOTE: We add 1 to the index to accound for
  // MATLAB's 1-indexed array structures
  for ( int i = 0; i < facets.size(); i++ )
    for ( int j = 0; j < 3; j++ )
      facesOut[i+(j*numFaces)] = (double) (facets[i][j]+1.0);

  for ( int i = 0; i < vertices.size(); i++ )
    for ( int j = 0; j < dim; j++ )
      vertexOut[i+(j*numVertices)] = vertices[i][j];

  return;

};
