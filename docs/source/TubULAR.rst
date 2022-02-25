TubULAR
=========


.. toctree::
   :maxdepth: 2


Tube-like sUrface Lagrangian Analysis Resource (TubULAR) class
    

Coordinate Systems
==================
**uv** :  (conformal map onto unit square)
      Conformally mapping the cylinderCutMesh onto the unit square 
      in the plane results in the instantaneous uv coordinate system. 
      Corners of the unit square are taken directly from the cutMesh,
      so are liable to include some overall twist.
**sphi** : (proper length x rectified azimuthal coordinate)
    Cylinder-like coordinate system in which first coordinate is the 
    proper length along the surface and second is a rectified azimuthal
    coordinate. Rectification means that the surface is rotated at 
    each longitudinal coordinate ``u`` so that the surface changes as 
    little as possible over time. Mathematically, ``v`` is taken to 
    ``phi(u)``, where each ``u`` coordinate is offset by a "rotation" 
    about the surface along the ``v`` direction (the circumferential direction). 
    This rectification may be based on surface positions in R^3 (geometric) or 
    based on intensity motion in pullback space (material/Lagrangian) 
    inferred through phase-correlation of the tissue at each ``u`` 
    value (averaged across ``v``).
**sphi_sm** : (proper length x rectified azimuthal coordinate, smoothed in time)
	  Same as sphi, but coordiantes smoothed over time with a tripulse 
	  filter.
**APDV** : (rotated, translated, and scaled coordinate system, like lab frame)
    Physical coordinates used to visualize dynamics. Aligned meshes are in this 
    3D coordinate system.
**xyz** : (data coordinate system)
    Pixel coordinates of the data itself. Raw meshes, cylinder meshes, cutMeshes 
    etc are in this 3D coordinate system.	  

PIV measurement can be automatically computed in pullback space to determine 3D motion. These measurements can be done in a coordinate system of your choice (uv, sphi, or sphi_sm).
	  

APDV coordinate system and Centerline specification
===================================================
QuapSlap allows the user to designate the AP axis, a DV axis, and a
centerline for the segmented object. 
The default APDV coordinate system is automatically determined from the 
elongation axis of the surface and the data axes. 

To customize, you can designate an APDV coordinate system that may be 
biologically relevant by using iLastik training on the
timepoint t0 (which is an attribute of tubi, found as tubi.t0). 
Train in iLastik on anterior (A), posterior (P), background (B), and dorsal (D) location 
in different iLastik channels, making small blobs of high probability near
the anterior end (A), somewhere along the posterior end for (P) so that AP forms the AP axis, and 
anywhere along dorsal for D. Then the centers of mass of
thresholded contiguous high probability are computed for each to
define the coordinate system. By default for this customized option, the ilastik results are
read as ch1=A, ch2=P, ch3=B, ch4=D, but you may specify anteriorChannel, 
posteriorChannel, and dorsalChannel to specify the iLastik
training channel that is used for each.
Name the h5 file output from iLastik as 
..._Probabilities_APDVcoords.h5
For example, dorsal for the gut was chosen at the fused site where 
additional 48YGAL4-expressing muscle-like cells form a seam.
Posterior is at the rear of the yolk, where the endoderm closes, for 
apical surface training. Anterior is at the junction of the midgut 
with the foregut.
Separately, define the AP points for centerline extraction. For most gut
data, the posterior point is different in centerline than it is for AP
centerline specification, so we use a different ILP to train for A and P
for all timepoints of dynamic data. 

Properties
==================
	  
xp : ImSAnE experiment class instance or struct with fields
	expMeta : struct with fields
	fileMeta : struct with fields
dynamic : bool
	true if multiple timepoints, false if fixed
timeInterval : numeric, default=1        
    increment in time between timepoints with 
    indices differing by 1. For example, if
    timePoints are [0,1,2,4] and these are 
    [1,1,2] minutes apart, then timeInterval 
    is 1. 
timeUnits : str (default='min')
	units of the timeInterval (ex 'min')
spaceUnits : str (default = '$\mu$m')
	units of the embedding space (ex '$\mu$m')
imSize : 2x1 int
	size of pullback images to create (default is [a_ratio * 1000, 1000])
dir : struct with fields for different saved data names
	directories where QuapSlap data lives
fileName : struct with fields for different saved data names
	Full filenames (with the absolute path) for various data and results
fileBase : struct with fields for different saved data names
	Relative fileNames to be populated by timestamp ('...%06d...mat')
fullFileBase :  struct with fields for different saved data names
	full path of filenames (like fullfile(tubi.dir.X, tubi.fileBase.X))
ssfactor : int
	subsampling factor for probabilities 
APDV : struct with fields
    'resolution', float 
	    resolution of data in spaceUnits / pixel
    'rot' : 3x3 rotation matrix
		rotation matrix to transform 3d pts from data frame into APDV frame accordint to v' = (rot*v+trans)*resolution
    'trans' : 1x3 float array
		translation vector to transform 3d pts from data frame into APDV frame according to v' = (rot*v+trans)*resolution
flipy  : bool
	whether data is mirror image of lab frame coordinates
nV : int (default=100)
	sampling number along circumferential axis
nU  : int (default=100) 
	sampling number along longitudinal axis
uvexten : string, build during instantiation
	naming extension with nU and nV like '_nU0100_nV0100'
t0 : numeric
	reference time in the experiment, for building rectified sphi coordinates, also default for building pathlines
normalShift : numeric (default=0)                 
	shift to apply to meshes in pixel space along normal direction
a_fixed : numeric (default=1)
	aspect ratio for fixed-width pullbacks in uv and sphi coordinates
phiMethod : str specifier ('3dcurves' or 'texture')
    method for stabilizing the ``v`` axis: either use the geometry of the surface to minimize
    surface deformation from timepoint to timepoint, so that v(t) for a given value of ``u`` is 
    very near v(t-1) at the same value of ``u``. 
    This is a method for determining Phi map in pullback mesh creation, with 
    the full map from embedding to pullback being [M'=(Phi)o(Z)o(M)], where M is a nearly-conformal mapping to the
    plane, Z changes the spacing of the longitudinal coordinates to reflect the average proper length traversed from 
    the 'ring' of u=const to the 'ring' of u=const+du.
    This string specifier must be '3dcurves' (geometric phi stabilization) 
    or 'texture' (optical flow phi stabilization)
endcapOptions : struct with fields 
	adist : numeric
		distance around anterior point A which is removed from tubular mesh (sliced off)
	pdist : numeric
		distance around anterior point P which is removed from tubular mesh (sliced off) 
	tref : numeric
		timestamp for reference time used to define th point on the endcap
	    at which we cut the cylinder mesh into a cylinderCutMesh (a topological disk/square). 
	    This "dorsal" point for other timepoints are identified by pointmatching.
	Additional fields allowed 
plotting : struct with fields 
    'preview' : bool (default=false)
		display intermediate results
    'save_ims' : bool (default=true)
		save images along the way
    'xyzlim_um_buff' : 3x2 float
		xyzlimits in um in RS coord sys with buffer
    'xyzlim_raw' : 3x2 float
		xyzlimits in pixels
    'xyzlim_pix' : 3x2 float
		xyzlimits in pixels, in rotated and translated coordinates
    'xyzlim_um' : 3x2 float
		xyzlimits in um in rotated, translated, and scaled coordinate system
    'colors' : Nx3 RGB values
		color cycle for tubi
apdvPts = struct('anteriorPts', [], ...
    'posteriorPts', [], ... 
    'antPts_sm', [], ...
    'postPts_sm', [], ... 
    'dorsalPts', [], ... 
    'antPts_rs', [], ... 
    'postPts_rs', [], ... 
    'dorsPts_rs', [])
apdvOptions : optional struct with fields 
	options used for finding APDV frame
currentTime : numeric
    timestamp, assigned whenever time is set
currentMesh = struct with fields
    'rawMesh' : struct with fields 
		original mesh found by surface detection
    'alignedMesh': struct with fields    
		APDV rotated and scaled mesh (raw mesh in APDV coordinates)     
    'cylinderMesh', [],
		original mesh with endcaps cut off
    'cylinderMeshClean', [], 
		cylinder mesh with "ears" removed (ears give difficulty in mapping to the plane)
    'cutMesh', [],
		cylinder mesh with a seam given by cutPath
    'cutPath', [], 
		vertex indices of the cutMesh along which the periodic seam is cut
    'uvcutMesh', [],
		rectilinear cutMesh in (u,v) from Dirichlet map result to rectangle 
    'spcutMesh', [],
		rectilinear cutMesh in (s,phi) 'surface Lagrangian' parameterization
    'spcutMeshSm', [],
		rectilinear cutMesh in (s,phi) smoothed in time
    'spcutMeshSmRS', [],
		rectilinear cutMesh in (s,phi) smoothed in time with rotated scaled embedding
    'spcutMeshSmRSC', [],
		rectilinear cutMesh as closed cylinder (topological annulus), in (s,phi) smoothed, with rotated scaled embedding
    'ricciMesh' : struct with fields 
		ricci flow result pullback mesh, topological annulus          
currentCline = struct('mss', [], ...
    'mcline', [], ...
    'avgpts', []) ;
data = struct('adjustlow', 0, ...
    'adjusthigh', 0, ...
    'axisOrder', [1 2 3], ...
    'ilastikOutputAxisOrder', 'cxyz') options for scaling and transposing image intensity data
currentData = struct('IV', [], ...
    'adjustlow', 0, ...
    'adjusthigh', 0 )           image intensity data in 3d and scaling
currentVelocity : struct with field
	'piv3d' : struct with current timepoint's PIV information in both 2d and 3d
piv : struct with fields for all timepoints' PIV information
    'imCoords' : str specifier (default='sp_sme')
		image coord system for measuring PIV / optical flow) ;
    'Lx' :  int
		width of image, in pixels (x coordinate)
    'Ly' : int 
		height of image, in pixels (y coordinate)
    'raw' : struct
		raw PIV results from disk/PIVLab
    'smoothed' : 
		smoothed PIV results after gaussian blur
    'smoothing_sigma' numeric (default=1 ) 
		sigma of gaussian smoothing on PIV, in units of PIV sampling grid pixels
velocityAverage : struct with fields
  vsmM : (#timePoints-1) x (nX*nY) x 3 float array
      3d velocities at PIV evaluation coordinates in um/dt rs
  vfsmM : (#timePoints-1) x (2*nU*(nV-1)) x 3 float array
      3d velocities at face barycenters in um/dt rs
  vnsmM : (#timePoints-1) x (nX*nY) float array
      normal velocity at PIV evaluation coordinates in um/dt rs
  vvsmM : (#timePoints-1) x (nU*nV) x 3 float array
      3d velocities at (1x resolution) mesh vertices in um/min rs
  v2dsmM : (#timePoints-1) x (nX*nY) x 2 float array
      2d velocities at PIV evaluation coordinates in pixels/ min
  v2dsmMum : (#timePoints-1) x (nX*nY) x 2 float array
      2d velocities at PIV evaluation coordinates in scaled pix/min, but 
      proportional to um/min (scaled by dilation of map)
    'v3d'
		3D velocities in embedding space [pix/dt]
    'v2d'
		2D tangential velocities in pullback
    'v2dum'
		2D tangential velocity scaled by speed in true embedding space
    'vn'
		normal velocity in spaceUnits per timeInterval timeUnits
    'vf'
		velocity vielf on face barycenters after Lagrangian avg
    'vv', []) ;  
		velocity field on vertices after Lagrangian avg
cleanCntrlines : 
	centerlines in embedding space after temporal averaging
smoothing : struct with fields
    'lambda' : float (default=0.00)
	      diffusion const for field smoothing on mesh, if nonzero
    'lambda_mesh' : float (default=0.00)      
		diffusion const for vertex smoothing of mesh itself, if nonzero
    'nmodes' : int (default=7)
		number of low freq modes to keep per DV hoop
    'zwidth' : int (default=1) 
		half-width of tripulse filter applied along zeta/z/s/u direction in pullback space, in units of du/dz/ds/dzeta
pathlines : struct with fields
    't0' : numeric
		time at which pathlines are rectilinear (from which time the points are advected).
		This is a timestamp (not an index) at which pathlines form regular grid in space.
    'refMesh'
		reference mesh for pathline advection
    'piv', []
		Lagrangian pathlines from piv coords
    'vertices', [], 
		Lagrangian pathlines from mesh vertices
    'vertices3d',
		Lagrangian pathlines from mesh vertices
    'faces', [],
		Lagrangian pathlines from mesh face barycenters
    'beltrami', struct with fields for beltrami coefficient evaluated along pathlines
       	'mu_material', [], ...
        'mu_material_filtered', [], ...
        'mu_material_vertices', [], ...
        'fitlerOptions', struct with fields describing filtering     
currentStrain : struct with fields
    'pathline' :
		strain from pathlines
    struct('t0Pathlines', [], ...   
		t=0 timepoint for pathlines in question
    'strain' :
		strain evaluated along pathlines
    'beltrami' : 
		beltrami coefficient for pathlines
	

Full Contents:
==================

.. automodule:: @TubULAR
    :show-inheritance:
    :members:

	
	
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
