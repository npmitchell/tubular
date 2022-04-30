%% Convert .ply mesh to .off mesh and save

clear; close all; clc;
addpath('../ply_codes/')

% Prepare path for meshes ===========================================
% path = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/data/48Ygal4-UAShisRFP/2016021015120_objFiles/';
path = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/simulation/simple_sequence/';
path = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/data/48Ygal4-UAShistRFP/201901021550_folded_2part/plys_depth10/'
% meshes = dir([path, 'chiral_tube*.ply']) ;
meshes = dir([path, 'mesh_apical*.ply']) ;

outdir = fullfile(path, 'off_meshes/') ;

% Create output dirs
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Prepare for iteration
% chisums is the integral chirality
dmyk = 1
dt = 1;
todo = 1:dt:min(length(meshes), 71) ;


% Load each PLY mesh and convert to OFF
for ii = todo
    % Load the mesh
    meshfn = fullfile(meshes(ii).folder, meshes(ii).name) ;
    disp(['Loading ', meshfn])
    [tri, pts] = ply_read(meshfn, 'tri') ;
    tri = tri' ;
    pts = pts' ;
    basename = split(meshes(ii).name, '.ply') ;
    outfn = fullfile(outdir, [basename{1}, '.off']) ;
    
    triang = triangulation( tri, pts);
    write_mesh_off(outfn, triang)
end
disp('done')
    