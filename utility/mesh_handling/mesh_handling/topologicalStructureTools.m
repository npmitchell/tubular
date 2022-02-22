function [eIDx, feIDx, bulkEdgeIDx] = topologicalStructureTools(tri)
% TOPOLOGICALSTRUCTURETOOLS(tri)
%   Compute topological structure tools of triangulation: bond list, face
%   bond IDs, and label of whether bonds are in bulk or edge
% 
% Parameters
% ----------
% tri : triangulation object 
%    triangulation object with properties Points, ConnectivityList
% 
% Returns 
% -------
% eIDx : #bonds x 2 int
%   bond vertex IDs
% feIDx : #faces x 3 int
%   face bond IDs. feIDx(i, :) gives the bond indices (indices into eIDx) 
%   of face i. 
% bulkEdgeIDx : #bonds x 1 int
%   label of whether in bulk (0) or on edge (1)
%
% Dillon Cislo & NPMitchell 2020

% Construct Edge-Face Correspondence -------------------------------------- 
% A array of the fIDs of the faces attached to a particular edge
% If an edge is a border edges ( i.e., only attached to a single face ),
% then that fID is listed twice for dimensional consistency
resizeCell = @(x) repmat( x, 1, 1+mod(numel(x),2) );
edgeFace = tri.edgeAttachments( tri.edges );
edgeFace = cell2mat(cellfun( resizeCell, edgeFace, 'Uni', false ));

% A list of bulk edge IDs
bulkEdgeIDx = diff(edgeFace, 1, 2) ~= 0;

% Construct Face-Edge Correspondence --------------------------------------
F = tri.ConnectivityList ;
e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
e3IDx = sort( [ F(:,2), F(:,1) ], 2 );

eIDx = tri.edges;
eIDx = sort( eIDx, 2 );

[~, e1IDx] = ismember( e1IDx, eIDx, 'rows' );
[~, e2IDx] = ismember( e2IDx, eIDx, 'rows' );
[~, e3IDx] = ismember( e3IDx, eIDx, 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];
