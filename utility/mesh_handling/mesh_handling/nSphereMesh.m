function [ pts, tri ] = nSphereMesh(npts, type)
% nSphereMesh(generation, type)
% Generate a nearly spherical mesh of at least npts vertices by bisecting 
%   a platonic solid
% 
% Parameters
% ----------
% n : int
%   minimum allowed number of vertices on the sphere
% type : string specifier ('oct', 'tet', 'ico', default='ico')
%   type of Platonic solid to begin with for bisecting
% 
% Returns
% -------
% pts : Nx3 float array
%   vertices of approx spherical mesh
% tri : #faces x 3 int array
%   connectivity list indexing into P for each face
% 
% NPM 2020

if nargin < 2
    type = 'ico' ;
end

nout = 0 ;
generation = 0 ;
while nout < npts
    [ pts, tri ] = generateSphereMesh( generation, type );
    nout = size(pts, 2) ;
    generation = generation + 1 ;
    disp(['sphere has ' num2str(nout) ' pts'])
end
pts = pts' ;
disp(['done with ' num2str(nout) ' pts on the sphere'])
