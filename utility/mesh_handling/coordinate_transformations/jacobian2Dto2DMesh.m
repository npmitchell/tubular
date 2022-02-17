function J_2D_to_2D = jacobian2Dto2DMesh(newv2d, oldv2d, ff)
%JACOBIAN2DTO2DMESH(oldV2D, newV3D, F) Construct jacobian for mapping
%2d->2d
%
% Parameters
% ----------
% newv2d : #vertices x 2 float array
%   The updated vertex locations in 2d (image space)
% oldv2d : #vertices x 2 float array
%   The initial vertex locations in 2d (domain space)
% ff : #faces x 3 int array
%   connectivity list for triangulated meshes
% 
% Returns
% -------
% J_2D_To_2D : length(F) x 1 cell array
%   A cell array containing the jacobian for each face as each element
% 
% Dillon Cislo, Noah Mitchell 2021

J_2D_to_2D = cell(size(ff,1), 1);

for f = 1:size(ff,1)
    
    % (3x1) Edge vectors of the current face in the initial configuration
    ej0 = [(oldv2d(ff(f,1),:)-oldv2d(ff(f,3),:)).'; 0];
    ek0 = [(oldv2d(ff(f,2),:)-oldv2d(ff(f,1),:)).'; 0];
    
    % Calculate the unit normal of the initial face and its double area
    n0 = cross(ej0,ek0,1);
    dblA0 = sqrt(sum(n0.^2));
    n0 = n0 ./ dblA0;
    
    % (3x1) Edge vectors of the current face in the updated configuration
    ej = [(newv2d(ff(f,1),:)-newv2d(ff(f,3),:)).'; 0];
    ek = [(newv2d(ff(f,2),:)-newv2d(ff(f,1),:)).'; 0];
    
    % The (3x3) Jacobian matrix for the current face
    J = (ek * cross(n0,ej0,1).' - ej * cross(n0,ek0,1).') ./ dblA0;
    
    % Extract only the (2x2) components
    J_2D_to_2D{f} = J(1:2, 1:2);
    
end


end

