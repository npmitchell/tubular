function [cutMeshClosed, glue2cut] = glueRectCylinderCutMeshSeam(cutMesh)
%glueRectCylinderCutMeshSeam(cutMesh)
% close a rectilinear cutMesh by gluing the seam back together (convert
% from topological disk to topological cylinder).
% 
% Parameters
% ----------
% cutMesh : struct
%   the cylinderCutMesh object with fields
%       nU : int
%       nV : int
%       v : (nU*nV) x 3 float array
%       f : (nU-1)*(nV-1) x 3 int array
%       vn : (nU*nV) x 3 float array [optional]
%       u : nU*nV 
%
% Returns
% -------
% cutMeshClosed : struct with fields
%       nU : int
%       nV : int
%       v : (nU*nV) x 3 float array
%       f : (nU-1)*(nV-1) x 3 int array
%       vn : (nU*nV) x 3 float array [optional]
%       u : nU*(nV-1) x 3 float array 
%
% glue2cut : #vertices in cutMesh x 1 int
%   glue2cut(i) gives the index in output glued mesh of the ith vertex in 
%   the cut mesh, so cutVtx = glueVtx(glue2cut)
% 
% See also
% ---------
% closeRectilinearCylMesh.m --> uses encaps as well
% glueCylinderCutMeshSeam.m --> allows for non-rectlinear meshes (todo)
%
% NPMitchell 2020
%

% Ensure vertices are list of vectors, not a rectilinear structure -- 
% that is, shape is [nU*nV, 3], not [nU,nV,3]
nU = cutMesh.nU ;
nV = cutMesh.nV ;
% Make the u and v fields are list of vectors not gridded structures 
if ~size(cutMesh.u, 2) == 2 || any(size(cutMesh.u) == nU)
    cutMesh.u = reshape(cutMesh.u, [nU * nV, 2]) ;
end
if ~size(cutMesh.v, 2) == 3 || any(size(cutMesh.v) == nU)
    cutMesh.v = reshape(cutMesh.v, [nU * nV, 3]) ;
end

cutMeshClosed = struct();
cutMeshClosed.v = cutMesh.v(1:end-nU, :) ;
% Redefine the last row of the cutMesh to be the same as the first
% Take nU*(nV-1):nU*nV --> 0:nU-1 --> 1:nU
cutMeshClosed.f = mod(cutMesh.f, (nV-1)*nU + 1) ;
% Since modulo sets repeated vertex 1,2,3 -> 0,1,2, we adjust indices by 1.
cutMeshClosed.f(cutMesh.f > (nV-1)*nU) = cutMeshClosed.f(cutMesh.f > (nV-1)*nU) + 1 ;
if isfield(cutMesh, 'vn')
    cutMeshClosed.vn = cutMesh.vn(1:end-nU, :) ;
else  
    cutMeshClosed.vn = per_vertex_normals(cutMeshClosed.v, ...
        cutMeshClosed.f, 'Weighting', 'angle') ;
end
cutMeshClosed.u = cutMesh.u(1:end-nU, :) ;
cutMeshClosed.nU = nU ;
cutMeshClosed.nV = nV ;

if nargout > 1
    % produce a map from the glued mesh back to the cut mesh indices, 
    % so that cutVtx = glueVtx(glue2cut)
    glue2cut = [ 1:(length(cutMesh.u(:, 1))-nU), 1:nU ] ;
end
