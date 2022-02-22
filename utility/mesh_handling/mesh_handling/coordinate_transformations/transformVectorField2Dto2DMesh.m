function [vfield2d, jacobian_2d_to_2d] = ...
    transformVectorField2Dto2DMesh( ...
    old_vfield2d, oldv2d, newv2d, ff, fieldfaces )
%TRANSFORMVECTORFIELD2DTO2DMESH( old_vfield2d, oldv2d, newv2d, ff,
%fieldfaces)
% Transform a vector field defined on faces of a 2d mesh to the faces of an
% updated 2d mesh with the same face connectivity
%
% Parameters
% ----------
% old_vfield2d : N x 2 float array
%   vector field defined in the initial 2d configuration,
%   evaluated at faces indexed by fieldfaces
% oldv2d : #vertices x 2 float array
%   mesh vertices in 2d domain space
% newv2d : #vertices x 2 float array
%   mesh vertices in 2d image space
% ff : #faces x 2 int array
%   connectivity list
% fieldfaces : length(vfield2d, 1) 
%   face indices at which vector field is defined
%
% Returns
% -------
% vfield2d : N x 2 float array
%   The vector field mapped to the updated 3d mesh configuration
% jacobian_2d_to_2d : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element
% 
%
% Dillon Cislo, Noah Mitchell 2020

jacobian_2d_to_2d = jacobian2Dto2DMesh( newv2d, oldv2d, ff );

% Transform vector fields
vfield2d = zeros(size(old_vfield2d));
for f = 1:size(fieldfaces, 1)
    vfield2d(f,:) = ...
        (jacobian_2d_to_2d{fieldfaces(f)} * old_vfield2d(f,:).').';
end

end


