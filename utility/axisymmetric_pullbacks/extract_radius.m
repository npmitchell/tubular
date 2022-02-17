%% 
% Extract radius by taking ratio of geodesics to A and P from each point.
% Denote the distance to A as da and the distance to P as dp. 
% Equate the ratio da / (da + dp) to the pathlength, s. 
% Match to centerline point with pathlength s. Measure Euclidean distance.
%
% By Noah Mitchell
%==========================================================================

clear; close all; clc;

%% Parameters
overwrite = false ;
resave_ims = false ;
save_ims = true ;
normal_shift = 10 ;
a_fixed = 2 ;
preview = false ;
washout2d = 0.5 ;
washout3d = 0.5 ;

%% Add paths
addpath(genpath('/mnt/data/code/gptoolbox'));
addpath('/mnt/data/code/gut_matlab/') ;
addpath_recurse('/mnt/data/code/imsaneV1.2.3/external/') ;
addpath_recurse('/mnt/data/code/gut_matlab/plotting/') ;

%% Define some colors
blue = [0 0.4470 0.7410] ;
orange = [0.8500 0.3250 0.0980] ;
yellow = [0.9290, 0.6940, 0.1250] ;
purple = [0.4940, 0.1840, 0.5560] ;
green = [0.4660, 0.6740, 0.1880] ;
sky = [0.3010, 0.7450, 0.9330] ;
red = [0.6350, 0.0780, 0.1840] ;
brick = [0.800000 0.250000 0.330000] ;
light_green =[0.560000 0.930000 0.560000] ;
light_gray = [0.830000 0.830000 0.830000] ;
bwr = diverging_cmap([0:0.01:1], 1, 2) ;


%% Initialize Some Directory Definitions ==================================

% The top level data directory
meshDir = [ dataDir ...
    'msls_output_prnun5_prs1_nu0p00_s0p10_pn2_ps4_l1_l1' ];

% The dataset name base for the AD points
ADBase = '/mesh_apical_stab_%06d/adorsal';

% The dataset name base for the PD points
PDBase = '/mesh_apical_stab_%06d/pdorsal';

% The folder where the pullback images will be saved
nshift = strrep(sprintf('%03d', normal_shift), '-', 'n') ;
imFolder = fullfile(meshDir, ['PullbackImages_' nshift 'step'] ) ;
imFolder_e = [imFolder '_extended'] ;
imFolder_r = [imFolder '_relaxed'] ;
imFolder_re = [imFolder '_relaxed_extended'] ;
cutFolder = fullfile(meshDir, 'cutMesh') ;
cutMeshImagesDir = fullfile(cutFolder, 'images') ;
cylCutDir = fullfile(meshDir, 'cylindercut') ;
cylCutMeshOutDir = fullfile(cylCutDir, 'cleaned') ;
% The file name base for the cylinder meshes
cylinderMeshCleanBase = fullfile( cylCutMeshOutDir, ...
    'mesh_apical_stab_%06d_cylindercut_clean.ply' );


% The file constaing the AD/PD points
dpFile = fullfile( cylCutDir, 'ap_boundary_dorsalpts.h5' );

tomake = {imFolder, imFolder_e, imFolder_r, imFolder_re,...
    pivDir, cutFolder, cutMeshImagesDir, cylCutMeshOutDir} ;
for i = 1:length(tomake)
    dir2make = tomake{i} ;
    if ~exist( dir2make, 'dir' )
        mkdir(dir2make);
    end
end

%% Load rotation, translation, resolution
rot = dlmread(fullfile(meshDir, 'rotation_APDV.txt')) ;
xyzlim = dlmread(fullfile(meshDir, 'xyzlim_APDV_um.txt'), ',', 1, 0) ;
trans = dlmread(fullfile(meshDir, 'translation_APDV.txt'));
resolution = dlmread(fullfile(meshDir, 'resolution.txt'), ',', 1, 0) ;
if resolution(:) == resolution(1)
    resolution = resolution(1) ;
else
    error('Have not handled case for anisotropic resolution')
end

%% Iterate Through Time Points to Compute radii ========================

for t = xp.fileMeta.timePoints
    % Load the cutMesh
    cutMeshfn = fullfile(cutFolder, [fileNameBase, '_cutMesh.mat']) ;
    cutMeshfn = sprintf(cutMeshfn, t) ;
    load(cutMeshfn, 'cutMesh')

    % % Create global graph
    % GG = makeGraph(faceIn, vertexIn) ;
    % % Get min geodesic distance from the anterior cut point at the same phi
    % % value
    % apath = shortestpath( GG, cp1, cp2 )' ;
    % adist = 
    % % Get min geodesic distance from the anterior cut point at the same phi
    % % value
    % ppath = shortestpath( GG, cp1, cp2 )' ;
    % pdist = 
    
    
    % Consider all vertices
    
    % Get starting points (anterior)
    mesh3dfn =  sprintf( cylinderMeshCleanBase, t ) ;
    outadIDxfn = fullfile(cylCutMeshOutDir, 'apIDx.h5') ;
    outpdIDxfn = fullfile(cylCutMeshOutDir, 'pdIDx.h5') ;
    
    
    % Load the mesh, the anterior and posterior boundary IDs
    mesh = read_ply_mod(mesh3dfn) ;
    adIDx = h5read(outadIDxfn, ['/' sprintf('%06d', t)]) ;
    pdIDx = h5read(outpdIDxfn, ['/' sprintf('%06d', t)]) ;
    start_points = adIDx ;

    options = {} ;
    DD = perform_fast_marching_mesh(mesh.v, mesh.f, start_points, options) ;

    % perform_fast_marching_mesh - launch the Fast Marching algorithm on a 3D mesh.
    %
    %   [D,S,Q] = perform_fast_marching_mesh(vertex, faces, start_points, options)
    %
    %   vertex, faces: a 3D mesh
    %   start_points(i) is the index of the ith starting point .
    %
    %   D is the distance function to the set of starting points.
    %   S is the final state of the points : -1 for dead (ie the distance
    %       has been computed), 0 for open (ie the distance is only a temporary
    %       value), 1 for far (ie point not already computed). Distance function
    %       for far points is Inf.
    %   Q is the index of the closest point. Q is set to 0 for far points.
    %       Q provide a Voronoi decomposition of the domain. 
    
    
    
end
