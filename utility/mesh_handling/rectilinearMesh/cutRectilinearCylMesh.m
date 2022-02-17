function cutMesh = cutRectilinearCylMesh(mesh, options)
%cutRectilinearMesh(mesh, options)
%
% Parameters
% ----------
% mesh : struct, closed cylinder mesh with fields
%   nU : int
%   v : (nU*(nV-1)) x 3 float array, optional
%       3d vertices of the mesh embedding
%   u : (nU*(nV-1)) x 2 float array, optional
%       2d vertices of the rectilinear mesh in pullback space
%   f : #faces x 3 int array
%       indices into v (or equivalently into u) of mesh connectivity
%       (faces)
%   vn : (nU*nV) x 3 float array, optional
%       vertex normals 
% options : optional struct with fields
%   ignoreRectangularConstraint : bool
%       allow the assumption that 2d (pullback) coordinates are ordered as
%       rectilinear mesh, without checking this   
% 
% Returns 
% -------
% cutMesh : struct with same fieldnames as input mesh, with properties
%   nU : int
%   v : (nU*nV) x 3 float array
%       3d vertices of the mesh embedding
%   u : (nU*nV) x 2 float array
%       2d vertices of the rectilinear mesh in pullback space
%   f : #faces x 3 int array
%       indices into v (or equivalently into u) of mesh connectivity
%       (faces)
% If input mesh contains field vn, then output will also:
%   vn : (nU*nV) x 3 float array
%       vertex normals 
%
% See also
% --------
% 
% NPMitchell 2020

if nargin > 1
    if isfield(options, 'vmax')
        vmax = options.vmax ;
    else
        vmax = 1 ;
    end
    if isfield(options, 'ignoreRectangularConstraint')
        ignoreRectangularConstraint = options.ignoreRectangularConstraint ;
    elseif isfield(options, 'ignoreRectilinearConstraint')
        ignoreRectangularConstraint = options.ignoreRectilinearConstraint ;
    else
        ignoreRectangularConstraint = false ;
    end
else
    vmax = 1 ;
    ignoreRectangularConstraint = false ;
end

nU = mesh.nU ;
nV = length(mesh.v(:, 1)) / nU + 1;
cutMesh = mesh ;

% Duplicate the first row as the last row
if isfield(mesh, 'v')
    cutMesh.v(nU*(nV-1) + 1:nU*nV, :) = mesh.v(1:nU, :) ;
end

% Duplicate the first row as last in pullback space
if isfield(mesh, 'u')
    cutMesh.u(nU*(nV-1) + 1:nU*nV, 1) = mesh.u(1:nU, 1) ;
    cutMesh.u(nU*(nV-1) + 1:nU*nV, 2) = mesh.u(1:nU, 2) + vmax ;
end

% Redefine faces
if ignoreRectangularConstraint
    cutMesh.f = defineFacesRectilinearGrid([], nU, nV) ;    
else
    cutMesh.f = defineFacesRectilinearGrid(mesh.u, nU, nV) ;
end
assert(all(size(cutMesh.f) == size(mesh.f)))

% Duplicate normals for added points
if isfield(mesh, 'vn')
    cutMesh.vn((nV-1)*nU + 1:nU*nV, :) = cutMesh.vn(1:nU, :) ;
end
