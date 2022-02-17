function NL = TRI2NL(TRI, pts)
% triang2NL(triangulationObj) Convert a triangulation list (TRI) to Neighbor List (NL)
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
% NL : array of dimension #pts x max(#neighbors)
%   The NL(i,:) contains indices for the neighbors for the ith point, 
%   buffered by zeros if a particle does not have the maximum # nearest 
%   neighbors.
%
% NPMitchell 2019

% 1. Get all the triangles attached to a particular vertex in the
% triangulation.  
attachedTriangles = vertexAttachments(triangulation(TRI, pts));

NLcell = cell(size(pts, 1), 1) ;
for i = 1: size(NLcell,1)
    % 2. Use the connectivity list to get the vertex indices of all these
    % triangles
    verticesOfTI = TRI(attachedTriangles{i}, :);
    % 3. Find all the unique vertices and remove the current vertex
    NLcell{i} = setdiff(unique(verticesOfTI), i);
end

% Convert cell to matrix
% First get max cell element length (max #NNeighbors)
maxnn = 0 ;
for kk = 1:length(NLcell)
    maxnn = max(maxnn, length(NLcell{kk})) ;
end
NL = zeros(length(NLcell), maxnn) ;
for kk=1:length(NLcell)
    NL(kk, 1:length(NLcell{kk})) = NLcell{kk} ;
end
