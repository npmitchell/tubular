function [v0n, v0t, v0t2d, jac, facenormals, g_ab, dilation] = ...
    resolveTangentNormalVelocities(faces, vertices, v0, fieldfaces, ...
    vertices2d, varargin)
%RESOLVETANGENTNORMALVELOCITIES(faces, vertices, v0, vertices2d, fieldfaces, varargin)
% Resolve a 3d vector field into the tangential and normal components of a
% 2d mesh in its embedding and in its 2d image (flattened mesh)
%
% Parameters
% ----------
% faces : #faces x 3 int array
%   indices into vertices2d of the mesh faces
% vertices : #vertices x 3 float array
%   mesh embedding in 3d
% v0 : N x 3 float array
%   velocities in 3d evaluated at points living in faces specified by 
%   fieldfaces
% vertices2d : #vertices x 2 float array
%   mesh embedding in 2d
% fieldfaces : N x 1 int array 
%   indices into faces in which velocities v0 are evaluated
% vertices2d : #vertices x 2 float array
%   coordinates in 2d pullback, required if more than two output variables
% varargin : keyword arguments
%   'facenormals' : #faces x 3 float array, normal vectors for each face
%
% Returns
% -------
% v0n : N x 1 float array
%   normal velocities at each evaluation point
% v0t : N x 3 float array
%   original scale velocities in tangent plane (as 3d vectors)
% v0t2d : #fieldfaces x 2 float array
%   The vector field mapped to the 2d mesh, scaled + distorted by jacobian
% jac : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element -- 
%   transformation from 3d to 2d
% facenormals #fieldfaces x 3 float array
%   the normal vectors on each face where the vector field is defined
% g_ab : #fieldfaces x 2 x 2 
%   the metric tensor on each field face
% dilation : #fieldfaces x 1 float array
%   dilation factor that maps the magnitude of vf in 3d to the
%   magnitude of the vf in 2d pullback
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

        if ~isempty(regexp(varargin{i},'^[Ff]ace[Nn]ormals','match'))
            facenormals = varargin{i+1} ;
        end
    end
else
    % Option 1 : compute using matlab built-in
    facenormals = faceNormal( triangulation(faces, vertices) );
    % Option 2 : project faces onto the already-computed normals
    % aux_alternate_velocity_projection
end

%% Take dot product of flow fields with normals in 3D
v0n = dot(facenormals(fieldfaces, :), v0, 2) ;
% Subtract the normal velocity to obtain tangential velocity
v0t = v0 - v0n .* facenormals(fieldfaces, :) ;

%% Compute the tangential velocities in plane
if nargout > 2
    % u is 3d, w is 2d. jac takes u->w, jjac takes w->u
    [v0t2d, jac] = pullVectorField3Dto2DMesh(v0t, vertices2d, ...
        vertices, faces, fieldfaces) ;
    
    % % Inspect result
    % bc = barycenter(vertices2d, faces) ;
    % tri(triangulation(faces(fieldfaces, :), vertices2d(:, 1), vertices2d(:, 2)))
    % hold on;
    % quiver(bc(fieldfaces, 1), bc(fieldfaces, 2), v0t2d(:, 1), v0t2d(:, 2), 0)
    % bad = find(isnan(v0t2d(:, 1))) ;
    % plot(bc(bad, 1), bc(bad, 2), 'r.')
end

if nargout > 4
    % Compute dilation factor that maps the magnitude of vf in 3d to the
    % magnitude of the vf in 2d pullback
    g_ab = zeros(size(fieldfaces, 1), 2, 2);
    dilation = zeros(size(fieldfaces, 1), 1) ;
    for f = 1:length(fieldfaces)
        qg = jac{fieldfaces(f)} * jac{fieldfaces(f)}' ;
        g_ab(f, :, :) =  qg ;
        dilation(f) = sqrt(det(qg)) ;
    end
end


