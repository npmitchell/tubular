function [cutMeshClosed, glue2cut] = glueCylinderCutMeshSeam(cutMesh)
%GLUECUTMESHSEAM(cutMesh)
% close a cutMesh by gluing the seam back together 
% 2020: handled rectilinear case
% 2021: generalized to non-rectilinear cutMesh
% 
% Parameters
% ----------
% cutMesh : struct
%   the cylinderCutMesh object with fields
%       nU : int (ignored if pathPairs is supplied)
%       nV : int (ignored if pathPairs is supplied)
%       v : (nU*nV) x 3 float array
%       f : (nU-1)*(nV-1) x 3 int array
%       vn : (nU*nV) x 3 float array [optional]
%       u : nU*nV
%       pathPairs
% 
% Returns
% -------
% cutMeshClosed : struct
% glue2cut : #vertices in cutMesh x 1 int
%   glue2cut(i) gives the index in output glued mesh of the ith vertex in 
%   the cut mesh, so cutVtx = glueVtx(glue2cut)
% 
% See also
% ---------
% closeRectilinearCylMesh.m -- just the rectilinear case
%
% NPMitchell 2020-2021
%


% Rectangular case
if isfield(cutMesh, 'nU') && isfield(cutMesh, 'nV')
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

    if nargout > 1
        % produce a map from the glued mesh back to the cut mesh indices, 
        % so that cutVtx = glueVtx(glue2cut)
        glue2cut = [ 1:(length(cutMesh.u(:, 1))-nU), 1:nU ] ;
    end
else
    try
        assert(isfield(cutMesh, 'pathPairs'))
    catch
        error('pathPairs must be supplied if nU and nV are not / if mesh is not triangulated rectangle in lattice structure')
    end

    %--------------------------------------------------------------------------
    % Modify the face list on the right hand side
    %--------------------------------------------------------------------------
    
    % Iterate over cut path vertices
    for i = 1:length(cutMesh.pathPairs(:, 1))

        % The vertex ID of the current cut path vertex
        % pick one side of the cut (call it "right" side) 
        vv = cutMesh.pathPairs(i, 1);

        % The duplicate ID of the current vertex ("left" side)
        vP = cutMesh.pathPairs(i, 2);

        % Replace the vertex ID in ALL faces
        cutMesh.f( cutMesh.f == vP ) = vv;

    end
    
    % Build output
    cutMeshClosed = struct() ;
    
    % Remove the duplicate "left" copy of vertices
    [ cutMeshClosed.f, cutMeshClosed.v, unreferenced, oldVertexIDx, glue2cut_part1 ] = ...
        remove_unreferenced_vertices_from_mesh(cutMesh.f, cutMesh.v) ;
    assert(isempty(setdiff(unreferenced, cutMesh.pathPairs(:, 2))))
    
    % vertex normals
    if isfield(cutMesh, 'vn')
        cutMeshClosed.vn = cutMesh.vn(oldVertexIDx, :) ;
    else  
        cutMeshClosed.vn = per_vertex_normals(cutMeshClosed.v, ...
            cutMeshClosed.f, 'Weighting', 'angle') ;
    end
    cutMeshClosed.u = cutMesh.u(oldVertexIDx, :) ;
    
    % finish glue2cut
    if nargout > 1
        error('finish here: add the vertices that are gone using pathPairs indices and oldVertexIDx and such')
        glue2cut = glue2cut_part1 ;
    end
    

end
