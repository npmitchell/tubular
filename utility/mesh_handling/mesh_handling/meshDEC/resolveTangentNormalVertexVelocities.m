function [v0n, v0t, v0t2d, jac, vertexnormals, g_ab, dilation] = ...
    resolveTangentNormalVertexVelocities(faces, vertices, v0, fieldvertices, ...
    vertices2d, varargin)
%RESOLVETANGENTNORMALVELOCITIES(faces, vertices, v0, vertices2d, fieldvertices, varargin)
% Resolve a 3d vector field into the tangential and normal components of a
% 2d mesh in its embedding and in its 2d image (flattened mesh)
% NOTE: unlike similar functions for faces, this one returns v0t2d in units
% of the embedding length scale! For velocities in the units of the
% pullback, divide output v0t2d by dilation. 
% 
% TODO: fix v0t2d.field2Faces, which currently seems to be nonsense.
% TODO: handle more general periodic meshes (periodic in 1D)
%
% Parameters
% ----------
% faces : #faces x 3 int array
%   indices into vertices2d of the mesh faces
%   If cutPath is supplied, this must be the closed 3d mesh faces
% vertices : #vertices x 3 float array
%   mesh embedding in 3d
%   If cutPath is supplied, this must be the closed 3d mesh vertices
% v0 : N x 3 float array
%   velocities in 3d evaluated at points living in faces specified by 
%   fieldvertices
% vertices2d : #vertices x 2 float array
%   mesh embedding in 2d
% fieldvertices : N x 1 int array 
%   indices into faces in which velocities v0 are evaluated
% vertices2d : #vertices x 2 float array
%   coordinates in 2d pullback, required if more than two output variables
% varargin : keyword arguments
%   'vertexnormals' : #faces x 3 float array, normal vectors for each face
%   'cutPath' : cutPath of mesh if mesh is periodic, for mesh averaging
%   operators
%
% Returns
% -------
% v0n : N x 1 float array
%   normal velocities at each evaluation point
% v0t : N x 3 float array
%   original scale velocities in tangent plane (as 3d vectors)
% v0t2d : Struct with fields:
%   uv2vertices : #fieldvertices x 2 float array
%       The vector field mapped to the 2d mesh in the local pullback frame IN
%       THE EMBEDDING --> so in units of the embedding space!
%       This one is trustworthy.
%   field2Faces : this one seems not to work at present
% jac : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element -- 
%   transformation from 3d to 2d
% vertexnormals #fieldvertices x 3 float array
%   the normal vectors on each face where the vector field is defined. Note
%   that if this is passed as input, the same is returned.
% g_ab : #fieldvertices x 2 x 2 
%   the metric tensor on each field vertex
% dilation : struct with fields
%   #faces x 1 float array 
%       dilation factor that maps the magnitude of vf in 3d to the
%       magnitude of the vf in 2d pullback
%   #vertices x 1 float array 
%       dilation factor that maps the magnitude of vf in 3d to the
%       magnitude of the vf in 2d pullback, pushed to vertices
%
% NPMitchell 2021

preview = false ;

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
            vertexnormals = varargin{i+1} ;
        end
        if ~isempty(regexp(varargin{i},'^[Cc]ut[Pp]ath','match'))
            cutPath = varargin{i+1} ;
        end
        if ~isempty(regexp(varargin{i},'^[Pp]review','match'))
            preview = varargin{i+1} ;
        end
    end
end
if ~exist('cutPath', 'var')
    cutPath = [] ;
end

if ~exist('vertexnormals', 'var')
    % Option 1 : compute using matlab built-in
    vertexnormals = per_vertex_normals(vertices, faces, 'Weighting', 'angle');
    % Option 2 : project faces onto the already-computed normals
    % aux_alternate_velocity_projection
end

%% Take dot product of flow fields with normals in 3D
v0n = dot(vertexnormals(fieldvertices, :), v0, 2) ;
% Subtract the normal velocity to obtain tangential velocity
v0t = v0 - v0n .* vertexnormals(fieldvertices, :) ;

%% Compute the tangential velocities in plane
if nargout > 2
    [V2F, F2V] = meshAveragingOperators(faces, vertices) ;
    
    % Now make cutMesh so that when we unwrap, the seam is not mapped to a
    % set of long flipped triangles connecting the periodic dimension!
    
    % Todo: update to allow for other mesh local topologies! 
    % cylinderCutMesh(faces, vertices,[], cutPath(:, 1), cutPath(:, 2))
    nU = size(vertices2d, 1) - size(vertices, 1) ;
    nV = size(vertices2d, 1) / nU ;
    isInteger = (nV == round(nV)) ;
    if isInteger && all(all(cutPath == [1:nU; nU*(nV-1)+1:nU*nV]')) 
        mesh = struct() ;
        mesh.f = faces ; 
        mesh.v = vertices ;
        mesh.nU = nU ;
        mesh.nV = nV ;
        mesh.u = vertices2d(1:nU*(nV-1), :) ;
        cutMesh = cutRectilinearCylMesh(mesh) ;
    else
        error('Update code to allow for other connectivities!')
    end
    
    % Take xhat to 3d, yhat to 3d, then resolve v0t on each vector. 
    % This avoids averaging the field, instead averaging the hat vectors!
    xhats = ones(size(faces, 1), 2) ;
    xhats(:, 2) = 0 ;
    x3d = pushVectorField2Dto3DMesh(xhats, vertices2d, cutMesh.v, ...
        cutMesh.f, 1:size(faces, 1)) ;
    yhats = ones(size(faces, 1), 2) ;
    yhats(:, 1) = 0 ;
    y3d = pushVectorField2Dto3DMesh(yhats, vertices2d, cutMesh.v, ...
        cutMesh.f, 1:size(faces, 1)) ;
    
    % Check it
    % bc3 = barycenter(cutMesh.v, cutMesh.f) ;
    % close all
    % quiver3(bc3(:, 1), bc3(:, 2), bc3(:, 3), ...
    %     x3d(:, 1), x3d(:, 2), x3d(:, 3), 1)
    % hold on;    
    % quiver3(bc3(:, 1), bc3(:, 2), bc3(:, 3), ...
    %     y3d(:, 1), y3d(:, 2), y3d(:, 3), 1)
    % axis equal
    
    assert(all(vecnorm(x3d, 2, 2) > 0))
    assert(all(vecnorm(y3d, 2, 2) > 0))
    
    x3dV = F2V * x3d ;
    y3dV = F2V * y3d ;
    
    x3dV = normalizerow(x3dV) ;
    y3dV = normalizerow(y3dV) ;
    
    v0t2d_norms2Vertices = [dot(x3dV, v0t, 2), dot(y3dV, v0t, 2) ] ;
    
    %% Other approach is to move the field itself onto the faces
    % This doesn't seem to work
    % todo: fix it?
    [v0t2d_field2Faces, jac] = pullVectorField3Dto2DMesh(V2F * v0t, vertices2d, ...
        cutMesh.v, cutMesh.f) ;
    v0t2d_field2Faces = F2V * v0t2d_field2Faces ;
    
    v0t2d = struct() ;
    v0t2d.uv2vertices = v0t2d_norms2Vertices ;
    v0t2d.field2faces = v0t2d_field2Faces ;
    
    % Dilation
    g_ab = zeros(size(faces, 1), 2, 2);
    dilation = zeros(size(faces, 1), 1) ;
    for f = 1:length(faces)
        qg = jac{f} * jac{f}' ;
        g_ab(f, :, :) =  qg ;
        dilation(f) = sqrt(det(qg)) ;
    end
    % Return dilation on both faces AND vertices
    dd = struct() ;
    dd.faces = dilation ;
    dd.vertices = F2V * dilation ;
    dilation = dd ;
    v0t2d_field2Faces = v0t2d_field2Faces ./ dilation.vertices ;
    
    if preview
        % Check dilation
        d2d = dilation.vertices ;
        d2d(nU*(nV-1)+1:nU*nV) = d2d(1:nU) ;
        trisurf(cutMesh.f, vertices2d(:, 1), vertices2d(:, 2), 0*vertices2d(:, 1), d2d, 'edgecolor', 'none')
        view(2) ;
        axis equal;

        % preview agreement/discrepancy
        diffmag = vecnorm(v0t2d_norms2Vertices - v0t2d_field2Faces, 2, 2) ;
        alignmag = dot(v0t2d_norms2Vertices, v0t2d_field2Faces, 2) ;
        tmpV = reshape(v0t2d_norms2Vertices, [nU, nV-1, 2]) ;
        Vnorm = reshape(vecnorm(v0t2d_norms2Vertices, 2, 2), [nU, nV-1]) ;
        Fnorm = reshape(vecnorm(v0t2d_field2Faces, 2, 2), [nU, nV-1]) ;
        Vang = atan2(v0t2d_norms2Vertices(:, 2), v0t2d_norms2Vertices(:, 1)) ;
        Fang = atan2(v0t2d_field2Faces(:, 2), v0t2d_field2Faces(:, 1)) ;

        %% Check it
        clf
        quiver(mesh.u(:, 1), mesh.u(:, 2), ...
            v0t2d_norms2Vertices(:, 1), v0t2d_norms2Vertices(:, 2), 1)
        hold on;
        quiver(mesh.u(:, 1), mesh.u(:, 2), ...
            v0t2d_field2Faces(:, 1), v0t2d_field2Faces(:, 2), 1)
        figure
        plot(Vang(:), Fang(:), '.')
        plot(Vnorm(:), Fnorm(:), '.')
        close all

        %% Check each component in 3d
        maxval1 = rms1d(abs(v0t2d_field2Faces(:))) ;
        maxval2 = rms1d(abs(v0t2d_norms2Vertices(:))) ;
        maxval = max(maxval1, maxval2) ;

        subplot(2, 2, 1)
        trisurf(triangulation(faces, vertices), ...
            v0t2d_norms2Vertices(:, 1), 'edgecolor', 'none')
        title('Ratio of resolved tangent norms to tangent norms')
        colorbar
        caxis([-maxval, maxval])
        colormap(bwr)
        axis equal ; view([0,0])
        title('$(v_\perp)_u$', 'interpreter', 'latex')
        subplot(2, 2, 2)
        trisurf(triangulation(faces, vertices), ...
            v0t2d_norms2Vertices(:, 2), 'edgecolor', 'none')
        title('Ratio of resolved tangent norms to tangent norms')
        colorbar
        caxis([-maxval, maxval])
        colormap(bwr)
        axis equal ; view([0,0])
        title('$(v_\perp)_v$', 'interpreter', 'latex')

        %% Check 
        subplot(2, 2, 3)
        trisurf(triangulation(faces, vertices), ...
            v0t2d_field2Faces(:, 1), 'edgecolor', 'none')
        title('Ratio of resolved tangent norms to tangent norms')
        colorbar
        caxis([-maxval, maxval])
        colormap(bwr)
        axis equal ; view([0,0])
        title('$(v_\perp)_u$', 'interpreter', 'latex')
        subplot(2, 2, 4)
        trisurf(triangulation(faces, vertices), ...
            v0t2d_field2Faces(:, 2), 'edgecolor', 'none')
        title('Ratio of resolved tangent norms to tangent norms')
        colorbar
        caxis([-maxval, maxval])
        colormap(bwr)
        axis equal ; view([0,0])
        title('$(v_\perp)_v$', 'interpreter', 'latex')


        %% Check in 3d
        trisurf(triangulation(faces, vertices), ...
            Vnorm(:) ./ vecnorm(v0t, 2, 2), 'edgecolor', 'none')
        title('Ratio of resolved tangent norms to tangent norms')
        colorbar
        maxval = max(abs(Vnorm(:) ./ vecnorm(v0t, 2, 2) - 1)) ;
        caxis([1-maxval, 1+maxval])
        colormap(bwr)
        axis equal ;

        plotPolarField(Vnorm, Vang) ;

    end
    
    % % Inspect result
    % bc = barycenter(vertices2d, faces) ;
    % tri(triangulation(faces(fieldvertices, :), vertices2d(:, 1), vertices2d(:, 2)))
    % hold on;
    % quiver(bc(fieldvertices, 1), bc(fieldvertices, 2), v0t2d(:, 1), v0t2d(:, 2), 0)
    % bad = find(isnan(v0t2d(:, 1))) ;
    % plot(bc(bad, 1), bc(bad, 2), 'r.')

    % Compute dilation factor that maps the magnitude of vf in 3d to the
    % magnitude of the vf in 2d pullback
    
    % u is 3d, w is 2d. jac takes u->w, jjac takes w->u
    % jac = jacobian3Dto2DMesh(vertices2d, vertices3d, faces) ;
        
end


