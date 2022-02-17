function J_2D_To_3D = jacobian2Dto3DMesh(v2d, v3d, ff)
%JACOBIAN2DTO3DMESH(v2d, v3d, ff) Construct jacobian for mapping 2d->3d
% Acts on contravariant vectors, NOT on covariant vectors
%
% Parameters
% ----------
% v2d : #vertices x 2 float array
%   The vertex locations in 2d (domain space)
% v3d : #vertices x 3 float array
%   The vertex locations in 3d (image space)
% ff : #faces x 3 int array
%   connectivity list for triangulated meshes
% 
% Returns
% -------
% J_2D_To_3D : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element
% 
% Dillon Cislo, NPMitchell 2020

% Some convenience variables
uv12 = v2d(ff(:,2), :) - v2d(ff(:,1), :);
uv13 = v2d(ff(:,3), :) - v2d(ff(:,1), :);

xyz12 = v3d(ff(:,2), :) - v3d(ff(:,1), :);
xyz13 = v3d(ff(:,3), :) - v3d(ff(:,1), :);

% This just returns the second element of an array
get2 = @(A) A(2);

% I don't really have a physical explanation for this operator, but I'm
% sure there is one. It returns the second element of the wedge product of
% two column vectors. Here {u,v} are each (2x1) column vectors
JOp = @(u,v) get2(u * v' - v * u');

J_2D_To_3D = cell(size(ff,1), 1);

for f = 1:size(ff,1)
    
    % Calculate the each component of the operator
    JF = zeros([3 2]);
    
    JF(1,1) = -JOp([uv12(f,2); xyz12(f,1)], [uv13(f,2); xyz13(f,1)]);
    JF(2,1) = -JOp([uv12(f,2); xyz12(f,2)], [uv13(f,2); xyz13(f,2)]);
    JF(3,1) = -JOp([uv12(f,2); xyz12(f,3)], [uv13(f,2); xyz13(f,3)]);
    JF(1,2) = JOp([uv12(f,1); xyz12(f,1)], [uv13(f,1); xyz13(f,1)]);
    JF(2,2) = JOp([uv12(f,1); xyz12(f,2)], [uv13(f,1); xyz13(f,2)]);
    JF(3,2) = JOp([uv12(f,1); xyz12(f,3)], [uv13(f,1); xyz13(f,3)]);
    
    JD = JOp(uv12(f,:)', uv13(f,:)');
    
    JF = JF ./ JD;
    
    J_2D_To_3D{f} = JF;
     
end