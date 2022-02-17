function NLcell = triang2NLcell(triangulationObj)
% TRI2NL(TRI) Convert a triangulation object to Neighbor Cell array (NLcell)
% Note: This is slower than TRI2NLcell(TRI, pts) 
% 
% Parameters
% ----------
% triangulationObj : triangulation object 
%   a MATLAB triangulation object with properties Points, ConnectivityList
%   created, for ex, via triangulation(faces, vertices)
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
attachedTriangles = vertexAttachments(triangulationObj);
TRI = triangulationObj.ConnectivityList ;

NLcell = cell(size(triangulationObj.Points, 1), 1) ;
for i = 1: size(NLcell,1)
    % 2. Use the connectivity list to get the vertex indices of all these
    % triangles
    verticesOfTI = TRI(attachedTriangles{i},:);
    % 3. Find all the unique vertices and remove the current vertex
    NLcell{i} = setdiff(unique(verticesOfTI), i);
end

