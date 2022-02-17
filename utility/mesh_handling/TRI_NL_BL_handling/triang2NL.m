function NL = triang2NL(triangulationObj)
% triang2NL(triangulationObj) Convert a triangulation object to Neighbor List (NL)
% NOTE: This is actually slower than TRI2NL(TRI, pts). 
% 
% Parameters
% ----------
% triangulationObj : triangulation object 
%   a MATLAB triangulation object with properties Points, ConnectivityList
%   created, for ex, via triangulation(faces, vertices)
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
attachedTriangles = vertexAttachments(triangulationObj);

% For some reason copying is faster than indexing into the trngln object
TRI = triangulationObj.ConnectivityList ;

NLcell = cell(size(triangulationObj.Points, 1), 1) ;
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
