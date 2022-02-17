function [newV3d_mapSpace, fieldfaces, bc] = ...
    barycentricMap3DtoND(oldFaces, oldV3d, oldV_mapSpace, newV3d, options)
%[newV3d_mapSpace, fieldfaces, bc] = barycentricMap3DtoND(oldFaces, ...
%   oldV3d, oldV_mapSpace, newV3d, options)
% 
% Parameters
% ----------
% oldFaces : M x 3 int array
%   The old (sampling) mesh connectivity list, indexing into the old vertex 
%   array (should index into both oldV3d and oldV_mapSpace the same way)
% oldV3d : P x 2 float array
%   The old (sampling) mesh vertex locations in 3d sampling space
% oldV_mapSpace : P x 2 or 3 float array
%   The mesh vertex locations in mapping space (NDimensional)
% newV3d : N x 3 float array
%   The points to map using barycentric coordinates, in 3d sampling space
% 
%
% Returns
% -------
% newV3d_mapSpace : #newVtx x D float array
%   mapped points in a new space (D-dimensional space, same space as
%   oldV_mapSpace). 
% fieldfaces : N x 1 int array
%   indices into faces of the faces on which uv live. These are where 
%   the supplied field uv is defined on the mesh, and where the output
%   pts live on the 3d mesh.
% tr0 : triangulation object
%   triangulation of the input faces and 2d vertices (v2d)
%
% NPMitchell 2021
%

epsVtx = 1e-5 ;
preview = false ;
if nargin > 3
    if isfield(options, 'epsVtx')
        epsVtx = options.epsVtx;
    end
    if isfield(options, 'preview')
        preview = options.preview;
    end
end

[sqrD, qF, qV] = point_mesh_squared_distance(newV3d, oldV3d, oldFaces ) ;
% check that after remeshing, pushing the remeshed vertices onto the
% original mesh requires almost no motion in 3d
try
    assert(all(abs(sqrt(sqrD) - vecnorm(qV - newV3d, 2, 2))) < epsVtx) 
    assert(all(vecnorm(qV - newV3d, 2, 2) < epsVtx))
catch
    % Check visually
    clf
    trisurf(triangulation(oldFaces, oldV3d), ...
        'edgecolor', 'none', 'facecolor', 'k', 'faceVertexCData', 0) ; 
    hold on ; 
    scatter3(newV3d(:, 1), newV3d(:, 2), ...
        newV3d(:, 3), 10, sqrD, 'filled', 'markeredgecolor', 'none')
    axis equal
    colorbar

    error(['Supplied pts are far from the original ', ...
        'surface by more than ' num2str(epsVtx) ': max deviation is ' num2str(max(vecnorm(qV - newV3d, 2, 2) ))])
    
end

% Get barycentric coordinates via oldMesh, and find newV3d_mapSpace 
% NOTE: USE THE PROJECTED VERTEX LOCATIONS qV
bc = barycentric_coordinates(qV, oldV3d(oldFaces(qF, 1), :), ...
    oldV3d(oldFaces(qF, 2), :), oldV3d(oldFaces(qF, 3), :)) ;

if any(any(bc < 0))
    error('Some negative barycentric coordinates! Why?')
end
newV3d_mapSpace = bc(:, 1) .* oldV_mapSpace(oldFaces(qF, 1), :) + ...
     bc(:, 2) .* oldV_mapSpace(oldFaces(qF, 2), :) + ...
     bc(:, 3) .* oldV_mapSpace(oldFaces(qF, 3), :) ;

 if nargout > 1
     % qF are "query faces"
     fieldfaces = qF ;
 end
 
 