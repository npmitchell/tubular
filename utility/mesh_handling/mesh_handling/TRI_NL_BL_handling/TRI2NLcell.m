function NLcell = TRI2NLcell(TRI, pts)
% TRI2NL(TRI) Convert a triangulation object to Neighbor Cell array (NLcell)
% Note: This is slower than TRI2NLcell(TRI, pts) 
% 
% Parameters
% ----------
% TRI : #faces x 3 int array
%   triangle array, with TRI(i, :) being indices into pts of ith face
% pts : N x 3 float array
%   3D positions of points. Row i is the position of the ith particle.
%
%
% Returns
% -------
% NLcell : #pts x 1 cell array
%   NLcell{i} contains the indices (into triangulationObj.Points) of the
%   neighbors of the vertex i
%
% NPMitchell 2019

% 1. Get all the triangles attached to a particular vertex in the
% triangulation.  
attachedTriangles = vertexAttachments(triangulation(TRI, pts));

NLcell = cell(size(pts, 1), 1) ;
for i = 1: size(NLcell,1)
    % 2. Use the connectivity list to get the vertex indices of all these
    % triangles
    verticesOfTI = TRI(attachedTriangles{i},:);
    % 3. Find all the unique vertices and remove the current vertex
    NLcell{i} = setdiff(unique(verticesOfTI), i);
end

