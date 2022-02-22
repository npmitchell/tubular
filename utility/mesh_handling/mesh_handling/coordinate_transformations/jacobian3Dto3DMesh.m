function J_3D_to_3D = jacobian3Dto3DMesh(newv3d, oldv3d, ff)
%JACOBIAN3DTO3DMESH(oldv3d, newv3d, F) Construct jacobian for mapping
%3d->3d
%
% Parameters
% ----------
% newv3d : #vertices x 3 float array
%   The updated vertex locations in 3d (image space)
% oldv3d : #vertices x 3 float array
%   The initial vertex locations in 3d (domain space)
% ff : #faces x 3 int array
%   connectivity list for triangulated meshes
% 
% Returns
% -------
% J_3D_To_3D : length(F) x 1 cell array
%   A cell array containing the jacobian for each face as each element
% 
% Dillon Cislo, Noah Mitchell 2021

J_3D_to_3D = cell(size(ff,1), 1);

for f = 1:size(ff,1)
    
    % Edge vectors of the current face in the initial configuration
    ej0 = (oldv3d(ff(f,1),:)-oldv3d(ff(f,3),:)).';
    ek0 = (oldv3d(ff(f,2),:)-oldv3d(ff(f,1),:)).';
    
    % Calculate the unit normal of the initial face and its double area
    n0 = cross(ej0,ek0,1);
    dblA0 = sqrt(sum(n0.^2));
    n0 = n0 ./ dblA0;
    
    % Edge vectors of the current face in the updated configuration
    ej = (newv3d(ff(f,1),:)-newv3d(ff(f,3),:)).';
    ek = (newv3d(ff(f,2),:)-newv3d(ff(f,1),:)).';
    
    % The Jacobian matrix for the current face
    J = (ek * cross(n0,ej0,1).' - ej * cross(n0,ek0,1).') ./ dblA0;
    
    J_3D_to_3D{f} = J;
    
end


end

