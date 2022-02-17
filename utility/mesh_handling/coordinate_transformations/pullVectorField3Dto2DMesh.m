function [vfield2d, jac3d_to_2d] = ...
    pullVectorField3Dto2DMesh(vfield3d, v2d, v3d, ff, fieldfaces)
%PULLVECTORFIELD3DTO2DMESH(vfield2d, v2d, v3d, ff, fieldfaces)
% Pullback a vector field defined on faces of a mesh from a 3d mesh 
% embeddingto to its image.
%
% Parameters
% ----------
% vfield3d : N x 3 float array
%   vector field defined in 3d, evaluated at faces in the embedding mesh 
%   indexed by fieldfaces
% v2d : #vertices x 2 float array
%   mesh vertices in 2d target/image space
% v3d : #vertices x 3 float array
%   mesh vertices in 3d embedding/domain space 
% ff : #faces x 3 int array
%   connectivity list
% fieldfaces : length(vfield2d, 1) 
%   face indices at which vector field is defined
%
% Returns
% -------
% vfield2d : #fieldfaces x 2 float array
%   The vector field mapped to the 2d mesh
% jac3d_to_2d : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element
%
% Dillon Cislo, Noah Mitchell 2020

jac3d_to_2d = jacobian3Dto2DMesh(v2d, v3d, ff) ;

if nargin < 5
    fieldfaces = 1:size(ff, 1) ;
end

% Pullback Vector Fields to Domain of Parameterization ====================
vfield2d = zeros(size(fieldfaces, 1), 2);
for f = 1:length(fieldfaces)
    vfield2d(f,:) = jac3d_to_2d{fieldfaces(f)} * vfield3d(f,:)';
    
    if any(any(isnan(jac3d_to_2d{fieldfaces(f)})))
        debugMsg(1, ['jacobian is singular for face ' ...
            num2str(fieldfaces(f))])
    end
    
end


