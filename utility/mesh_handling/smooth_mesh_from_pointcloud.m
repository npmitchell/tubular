function [V, F] = smooth_mesh_from_pointcloud(PU, PNU)
% Smoothing procedure for zebrafish meshes
%
% Parameters
% ----------
% PU : Nx3 float
%   pointcloud points from read point cloud from file as:
%  [PU, ~, ~, ~, PNU, ~] = readOBJ(pointCloudFileName);
% PNY : Nx3 float
%   normal directions for the pointcloud
%
% Returns
% -------
% V : Qx3 float
%   smoothed mesh vertex coordinates
% F : Qx3 int
%   face connectivity list
%
%
% 

poisson_surface_reconstruction(PU, PNU, 'poisson_mesh.off');

% Read in result from file and delete original file
[V, F, ~, ~, ~] = readOFF('poisson_mesh.off');
delete('poisson_mesh.off');

[V, F] = smooth_remesh(V, F) ;
