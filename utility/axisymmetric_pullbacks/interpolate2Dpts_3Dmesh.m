function [pts, fieldfaces, tr0, baryc0] = ...
    interpolate2Dpts_3Dmesh(faces, v2d, v3d, uv)
%INTERPOLATE2DPTS_3DMESH Map points in 2D to 3D using mesh faces  
%   Map points in 2D space living in 2D representation of a mesh to 3D
%   space living on the 3D representation of the same mesh, using
%   barycentric coordinates. 
%
% Parameters
% ----------
% faces : M x 3 int array
%   The mesh connectivity list, indexing into the vertex array(s)
% v2d : P x 2 float array
%   The mesh vertex locations in 2d
% v3d : P x 3 float array
%   The mesh vertex locations in 3d
% uv : N x 2 float array
%   The points to map to 3D using barycentric coordinates
%
% Returns
% -------
% pts : N x 3 float array
%   the 3d point coordinates, living on 3D-representation of the mesh
% fieldfaces : N x 1 int array
%   indices into faces of the faces on which uv live. These are where 
%   the supplied field uv is defined on the mesh, and where the output
%   pts live on the 3d mesh.
% tr0 : triangulation object
%   triangulation of the input faces and 2d vertices (v2d)
% baryc0 : N x 3 float array
%   barycentric coordinates of the query points in mesh faces
%
% See also
% --------
% barycentericMap2d(faces, v2d, vmap, uv) -- essentially the same
% function but allows 2d to 2d
%
% NPMitchell 2019

% If uv is a single point, make sure it is the right shape
if isequal(size(uv), [2 1])
    uv = reshape(uv, [1, 2]) ;
end

tr0 = triangulation(faces, v2d) ;
[fieldfaces, baryc0] = pointLocation(tr0, uv) ; 

% Note: if we define mesh then we can use baryc to recover the position as:
% uu_check = [sum(baryc0 .* mesh.u(mesh.f(fieldfaces, :), 1)', 2), ...
%            sum(baryc0 .* mesh.u(mesh.f(fieldfaces, :), 2)', 2)] ;


% Handle case where there are NaNs
bad = find(isnan(fieldfaces)) ;
baryc0(isnan(fieldfaces), :) = 0 ; 
fieldfaces(isnan(fieldfaces)) = 1 ;  

% Interpolate the position in 3D given relative position within 2D
% triangle.
% x123(i) is the x coords of the elements of triangle t_contain(i)
vxa = v3d(:, 1) ;
vya = v3d(:, 2) ;
vza = v3d(:, 3) ;

if size(vxa, 1) ~= size(v2d, 1)
    disp(['size of v3d supplied = ' num2str(size(v3d))])
    disp(['size of v2d supplied = ' num2str(size(v2d))])
    error('Length of v3d and v2d must match (same #pts)')
end

% Map to the faces
tria = tr0.ConnectivityList(fieldfaces, :) ;

% Modulo the vertex IDs: trisa are the triangle vertex IDs
% trisa = mod(tria, size(vxa, 1)) ;
% trisa(trisa == 0) = size(vxa, 1) ;
% x123a = vxa(tria) ;
% y123a = vya(tria) ;
% z123a = vza(tria) ;

% Multiply the vertex positions by relative weights.
% Note that baryc gives weights of the three vertices of triangle
% t_contain(i) for pointLocation x0(i), y0(i)
if sum(size(tria)) == 4
    pts = [sum(baryc0 .* vxa(tria)', 2), ...
        sum(baryc0 .* vya(tria)', 2), ...
        sum(baryc0 .* vza(tria)', 2) ] ;
else
    pts = [sum(baryc0 .* vxa(tria), 2), ...
        sum(baryc0 .* vya(tria), 2), ...
        sum(baryc0 .* vza(tria), 2) ] ;
end

% Handle case where there are NaNs
pts(bad, :) = NaN  ;
fieldfaces(bad) = NaN ;

%%%
% Could attempt to fix the bad indices by interpolating
% Decided not to do this -- If behavior like this is desired, then focus 
% on mapping bad pts into the domain of validity in other code before 
% calling this function.
% if ~isempty(bad)
%     Xai = scatteredInterpolant(v2d(:, 1), v2d(:, 2), vxa, 'linear', 'nearest') ;
%     Yai = scatteredInterpolant(v2d(:, 1), v2d(:, 2), vya, 'linear', 'nearest') ;
%     Zai = scatteredInterpolant(v2d(:, 1), v2d(:, 2), vza, 'linear', 'nearest') ;
%     pts(bad, :) = [Xai(uv(bad, 1), uv(bad, 2)), ...
%         Yai(uv(bad, 1), uv(bad, 2)), Zai(uv(bad, 1), uv(bad, 2))] ;
%     fieldfaces(bad) = 
% end
%%

end

