function norms = faceNormals(ff, xyz)
%FACE_NORMALS(ff, xyz)
%
% Parameters
% ----------
% ff : #faces x 3 int array
%   indices into xyz of each face, consistently oriented (assumed)
% xyz : #vertices x 3 float array
% 
%
% Returns
% -------
% normals : #faces x 3 float array
%   "outward" facing normal
%   
% 
% NPMitchell 2020

v21 = xyz(ff(:, 2), :) - xyz(ff(:, 1), :) ;
v31 = xyz(ff(:, 3), :) - xyz(ff(:, 1), :) ;
nx = -v21(:, 2) .* v31(:, 3) + v21(:, 3) .* v31(:, 2) ;
ny = -v21(:, 3) .* v31(:, 1) + v21(:, 1) .* v31(:, 3) ;
nz = -v21(:, 1) .* v31(:, 2) + v21(:, 2) .* v31(:, 1) ;
norms = [nx, ny, nz] ;
% normalize the normals
norms = norms ./ sqrt(sum(norms.^2, 2)) ;
