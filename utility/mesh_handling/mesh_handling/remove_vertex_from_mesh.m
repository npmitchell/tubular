function [ newFace, newVertex, oldVertexIDx, newVertexIDx ] = ...
    remove_vertex_from_mesh( face, vertex, rmVIDx )
%REMOVE_VERTEX_FROM_MESH Removes a subset of the vertices and the
%associated faces from a mesh triangulation.
%   INPUT PARAMETERS:
%       - face:         #Fx3 face connectivity list
%       - vertex:       #VxD vertex coordinate list
%       - rmVIDx:       #Nx1 list of vertex IDs to remove
%
%   OUTPUT PARAMETERS:
%       - newFace:      #F'x3 updated face connectivity list
%       - newVertex:    #V'xD updated vertex coordinate list
%       - oldVertexIDx: #V'x1 list of the vertex IDs of the updated
%                       vertices in the old vertex list, so that, for ex,
%                       vn_new = old_mesh.vn(oldVertexIdx, :)
%       - newVertexIDx: #V x1 list of vertex IDs of the old vertices in
%                       the new vertex list, so that
%                       kept_old_vertices = new_mesh.v(newVertexIDx, :)
%
% By Dillon Cislo, modified NPMitchell 2019

% Alternate version by NPMitchell:
% 
%     fv.vertices(inside) = [] ;
%     % Now remove rows from faces list and lower all indices above the
%     % removed points
%     pts_to_remove = find(inside) ;
%     lower = 0 * fv.faces(:) ;
%     rows_to_remove = 0 * fv.faces(:, 1) ;
%     for jj = 1:length(pts_to_remove)
%         pt = pts_to_remove(jj) ;
%         lower(fv.faces(:) > pt) = lower(fv.faces(:) > pt) + 1;
%         % find rows where pt is a member
%         rowrmjj = sum(ismember(fv.faces, pt), 2) > 0 ;
%         rows_to_remove = rows_to_remove || rowrmjj ;
%         

% Remove the vertices at the specified index values
newVertex = vertex;
newVertex(rmVIDx, :) = [];

% Find the new index for each of the new vertices
[~, newVertexIDx] = ismember( vertex, newVertex, 'rows');

% Find any faces that contained the indicated vertices and remove them
newFace = face;

loc_v = false(size(face));
for i = 1:numel(rmVIDx)
    loc_v = loc_v | ( newFace == rmVIDx(i) );
end
loc_v = find( loc_v );

[row_v, ~] = ind2sub(size(newFace), loc_v);
row_v = unique(row_v);
newFace(row_v, :) = [];

% Now update the vertex indices in the connectivity list
newFace = newVertexIDx(newFace);

% Find the vIDs of the new vertices in the old vertex list
[ ~, oldVertexIDx] = ismember( newVertex, vertex, 'rows' );

end

