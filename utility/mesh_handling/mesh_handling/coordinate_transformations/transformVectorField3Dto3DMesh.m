function [vfield3d, jacobian_3d_to_3d] = transformVectorField3Dto3DMesh( ...
    old_vfield3d, oldv3d, newv3d, ff, fieldfaces)
%TRANSFORMVECTORFIELD3DTO3DMESH( old_vfield3d, oldv3d, newv3d, ff,
%fieldfaces)
% Transform a vector field defined on faces of a 3d mesh to the faces of an
% updated 3d mesh with the same face connectivity
%
% Parameters
% ----------
% old_vfield3d : N x 3 float array
%   vector field defined in the initial 3d configuration,
%   evaluated at faces indexed by fieldfaces
% oldv3d : #vertices x 3 float array
%   mesh vertices in 3d domain space
% newv3d : #vertices x 3 float array
%   mesh vertices in 3d image space
% ff : #faces x 3 int array
%   connectivity list
% fieldfaces : length(vfield3d, 1) 
%   face indices at which vector field is defined
%
% Returns
% -------
% vfield3d : N x 3 float array
%   The vector field mapped to the updated 3d mesh configuration
% jacobian_3d_to_3d : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element
% 
%
% Dillon Cislo, Noah Mitchell 2020

jacobian_3d_to_3d = jacobian3Dto3DMesh( newv3d, oldv3d, ff );

% Transform vector fields
vfield3d = zeros(size(old_vfield3d));
for f = 1:size(fieldfaces, 1)
    vfield3d(f,:) = ...
        (jacobian_3d_to_3d{fieldfaces(f)} * old_vfield3d(f,:).').';
end

end

