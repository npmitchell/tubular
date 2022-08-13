function [newFace, newFaceIdx, newVertex, newVertexIdx] = clip_mesh_mod(face, vertex)
% remove dangling triangles and ears from a mesh
%
% [newface, newfaceIdx] = clip_mesh(face)
% 
% Input Parameters:
%   face:         Connectivity list of triangulation
%   vertex:       Coordinates of triangulation points
%
% Output Parameters:
%   newFace:      Connectivity list of clipped triangulation 
%   newFaceIdx:   Indices of new faces in old triangulation
%   newVertex:    Coordinates of clipped triangulation points
%   newVertexIdx: Indices of new vertices in old triangulation
% 
% By Dillon Cislo

% Make certain that face is a (nf x 3) array
if size(face,1)<size(face,2)
    face=face';
end

% Lone points are defined as vertices belonging to ears (boundary triangles
% sharing only one edge with the mesh) or vertices belonging to triangles
% sharing only one vertex with the edge of the mesh
num_lone_points = Inf;

% Make temporary copies of the input parameters
f = face;
v = vertex;

while( num_lone_points ~= 0 )
    
    % Iterate through the vertices to find the number of lone points
    num_lone_points = 0;
    for i = 1:size(v,1)
        % The number of faces of which vertex i is a member
        vf = sum( any( f == i, 2 ) );
        if vf == 1
            num_lone_points = num_lone_points + 1;
        end
    end
    
    % Iterate through the vertices again to extract the IDs of the lone
    % points
    lone_points = zeros(num_lone_points,1);
    num_lone_points = 0;
    for i = 1:size(v,1)
        % The number of faces of which vertex i is a member
        vf = sum( any( f == i, 2 ) );
        if vf == 1
            num_lone_points = num_lone_points + 1;
            lone_points(num_lone_points) = i;
        end
    end
    
    % Remove the lone vertices at the specified index values
    newVertex = v;
    newVertex(lone_points, :) = [];
    
    % Find the new index for each of the new vertices
    [~, newVertexIdx] = ismember(v, newVertex, 'rows');
    
    % Find any faces that contained the lone vertices and remove them
    newFace = f;
    [~, loc_v] = ismember(lone_points, newFace);
    [row_v, ~] = ind2sub(size(newFace), loc_v);
    row_v = unique(row_v);
    newFace(row_v, :) = [];
    
    % Now update the vertex indices in the connectivity list
    newFace = newVertexIdx(newFace);
    
    % Update the temporary copies
    f = newFace;
    v = newVertex;
    
end

% Find the indices of the new vertices/faces in the old triangulation
[~, newVertexIdx] = ismember( vertex, newVertex, 'rows' );
[~, newFaceIdx] = ismember( face, newFace, 'rows' );

end


    
