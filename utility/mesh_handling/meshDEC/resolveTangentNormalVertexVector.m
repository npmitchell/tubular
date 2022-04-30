function [v0n, v0t, vertexnormals] = ...
    resolveTangentNormalVertexVector(faces, vertices, vectors, varargin)
%resolveTangentNormalVertexVector(faces, vertices, vectors, varargin)
% Resolve a 3d vector field defined on VERTICES into the tangential and 
% normal components. 
% Similar to resolveTangentNormalVector() but on vertices.
% Also like resolveTangentNormalVertexVelocities, but with fewer outputs.
%
% Parameters
% ----------
% faces : #faces x 3 int array
%   indices into vertices2d of the mesh faces
% vertices : #vertices x 3 float array
%   mesh embedding in 3d
% vectors : N x 3 float array
%   vectors in 3d evaluated at points living in faces specified by 
%   fieldvertices, or by faces in order supplied if varargin 'fieldvertices' not
%   supplied.
% varargin : optional keyword arguments
%   'fieldvertices' : N x 1 int array, optional (default=1:nFaces)
%       indices into faces in which velocities v0 are evaluated
%   'vertexnormals' : #faces x 3 float array
%       normal vectors for each face 
%
% Returns
% -------
% v0n : N x 1 float array
%   normal velocities at each evaluation point
% v0t : N x 3 float array
%   original scale velocities in tangent plane (as 3d vectors)
% vertexnormals #fieldvertices x 3 float array
%   the normal vectors on each face where the vector field is defined
%
% NPMitchell 2020

% Obtain face normals, either from varargin or compute them
if ~isempty(varargin)
    for i = 1:length(varargin)
        if isa(varargin{i},'double') 
            continue;
        end
        if isa(varargin{i},'logical')
            continue;
        end

        if ~isempty(regexp(varargin{i},'^[Vv]ertex[Nn]ormals','match'))
            vertexnormals = varargin{i+1} ;
        end
        if ~isempty(regexp(varargin{i},'^[Ff]ield[Vv]ertices','match'))
            fieldvertices = varargin{i+1} ;
        end
    end
end
if ~exist('vertexnormals', 'var')
    % Option 1 : compute using matlab built-in
    vertexnormals = per_vertex_normals(vertices, faces, 'Weighting', 'angle') ;
    % Option 2 : project faces onto the already-computed normals
    % aux_alternate_velocity_projection
end

% If fieldvertices is not supplied, then assume that vector(i, :) is
% associated with face(i, :).
if ~exist('fieldvertices', 'var')
    fieldvertices = 1:size(vertices, 1) ;
else
    try
        assert(size(vectors, 1) == size(vertices, 1))
    catch
        error('#vectors and #vertices must match unless fieldvertices supplied')
    end
end

%% Take dot product of flow fields with normals in 3D
v0n = dot(vertexnormals(fieldvertices, :), vectors, 2) ;
% Subtract the normal velocity to obtain tangential velocity
v0t = vectors - v0n .* vertexnormals(fieldvertices, :) ;
