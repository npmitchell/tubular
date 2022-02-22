function [vfield3d, jacobian_2d_to_3d] = ...
    pushVectorField2Dto3DMesh(vfield2d, v2d, v3d, ff, fieldfaces)
%PUSHVECTORFIELD2DTO3DMESH(vfield2d, v2d, v3d, ff, fieldfaces)
% Push a vector field defined on faces of a mesh from a 2d mesh to its 3d
% embedding.
% Acts on contravariant vectors, NOT on covariant vectors.
%
% Parameters
% ----------
% vfield2d : N x 2 float array
%   vector field defined in 2d, evaluated at faces indexed by fieldfaces
% v2d : #vertices x 2 float array
%   mesh vertices in 2d domain space
% v3d : #vertices x 3 float array
%   mesh vertices in 3d target/embedding space
% ff : #faces x 3 int array
%   connectivity list
% fieldfaces : length(vfield2d, 1) 
%   face indices at which vector field is defined
%
% Returns
% -------
% vfield2d : #faces x 2 float array
%   The vector field mapped to the 2d mesh
% jacobian_2d_to_3d : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element
% 
%
% Dillon Cislo, Noah Mitchell 2020

jacobian_2d_to_3d = jacobian2Dto3DMesh(v2d, v3d, ff) ;

if nargin < 5
    if size(vfield2d, 1) == size(ff, 1)
        fieldfaces = 1:size(ff, 1) ;
    else
        error('First pass fieldfaces, obtainable via pointLocation(tri)')
    end
end

% Push Vector Fields Forward to 3D 
vfield3d = zeros(size(vfield2d, 1), 3);
for f = 1:length(fieldfaces)
    vfield3d(f,:) = jacobian_2d_to_3d{fieldfaces(f)} * vfield2d(f,:)';
end
