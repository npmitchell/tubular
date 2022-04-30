function [cutMesh, mesh] = buildCylinderMesh(nU, nV) 
% buildCylinderMesh(nU, nV)
%   create a model cylinder mesh with a pullback space in the unit square
%   and an embedding space in R^3 of a cylinder of unit radius and unit
%   length
% 
% Parameters
% ----------
% nU : int
% nV : int
% 
% Returns
% -------
% cutMesh : cylinder mesh with seam, struct with fields
%   u : 2D vertex positions in unit square
%   v : 3D vertex positions (as cut cylinder)
%   nU : int
%   nV : int
% mesh : closed cylinder mesh, struct with fields
%   u : 2D vertex positions in unit square
%   v : 3D vertex positions (as closed cylinder)
%   nU : int
%   nV : int
%
% NPMitchell 2020

zz = linspace(0, 1, nU) ;
phi = linspace(0, 2*pi, nV) ;
[zz, phi] = meshgrid(zz, phi) ;
zz = zz' ;
phi = phi' ;
uv = [ zz(:), phi(:) ./ (2*pi) ] ;
vv = [ zz(:), cos(phi(:)), sin(phi(:)) ] ; 
faces = defineFacesRectilinearGrid(uv, nU, nV) ;
cutMesh.f = faces ;
cutMesh.v = vv ;
cutMesh.u = uv ;
cutMesh.nU = nU ;
cutMesh.nV = nV ;
mesh = glueCylinderCutMeshSeam(cutMesh) ;

