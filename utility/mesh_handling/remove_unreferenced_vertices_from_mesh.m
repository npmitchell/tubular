function [ face, vertex, unreferenced, oldVertexIDx, newVertexIDx ] = ...
    remove_unreferenced_vertices_from_mesh( face, vertex )
%REMOVE_UNREFERENCED_VERTICES_FROM_MESH 
%   Removes a subset of the vertices which have no associated faces from 
%   a mesh triangulation.
%
%   INPUT PARAMETERS:
%       - face:         #Fx3 face connectivity list
%       - vertex:       #VxD vertex coordinate list
%
%   OUTPUT PARAMETERS:
%       - face:         #F'x3 updated face connectivity list
%       - vertex:       #V'xD updated vertex coordinate list
%       - unreferenced: #N x1 list of vertex IDs to remove
%       - oldVertexIDx: #V'x1 list of the vertex IDs of the updated
%                       vertices in the old vertex list, so that, for ex,
%                       vn_new = old_mesh.vn(oldVertexIdx, :)
%       - newVertexIDx: #V x1 list of vertex IDs of the old vertices in
%                       the new vertex list, so that
%                       kept_old_vertices = new_mesh.v(newVertexIDx, :)
%
% NPMitchell 2020

unreferenced = [] ;
oldVertexIDx = 1:size(vertex, 1) ;
if ~isempty(setdiff(1:length(vertex), face(:)))
    unreferenced = setdiff(1:size(vertex, 1), face(:)) ;
    [ face, vertex, oldVertexIDx, newVertexIDx] = ...
        remove_vertex_from_mesh( face, vertex, unreferenced ) ;
end
        