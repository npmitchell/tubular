function varargout = removeMeshEars(varargin)
%REMOVEMESHEARS Remove vertices that are connected to only one face.
%
%   [V, F] = removeMeshEars(V, F)
%   [V, F] = removeMeshEars(MESH)
%   [V, F, verticesKept, facesKept] = removeMeshEars(V, F)
%   [V, F, verticesKept, facesKept] = removeMeshEars(MESH)
%   Remove vertices that are connected to only one face. This removes also
%   "pending" faces.
%   Note that if the mesh has boundary, this may remove some regular faces
%   located on the boundary.
%
%   Example
%   removeMeshEars
%
%   See also
%     meshes3d, ensureManifoldMesh
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2019-01-08,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2019 INRA - Cepia Software Platform.
% Edits NPMitchell 2021 for more output information

[vertices, faces] = parseMeshData(varargin{:});

nVertices = size(vertices, 1);

% for each vertex, determine the number of faces it belongs to
vertexDegree = zeros(nVertices, 1);
for iv = 1:nVertices
    vertexDegree(iv) = sum(sum(faces == iv, 2) > 0);
end

% remove vertices with degree 1
inds = find(vertexDegree == 1);
[vertices, faces, facesKept] = removeMeshVertices(vertices, faces, inds);


%% Format output
if nargout < 3
    varargout = formatMeshOutput(nargout, vertices, faces);
elseif nargout == 3
    verticesKept = find(vertexDegree ~= 1) ;
    varargout = {vertices, faces, verticesKept} ;
elseif nargout == 4
    verticesKept = find(vertexDegree ~= 1) ;
    varargout = {vertices, faces, verticesKept, facesKept} ;
else
    error('improper number of output arguments')
end