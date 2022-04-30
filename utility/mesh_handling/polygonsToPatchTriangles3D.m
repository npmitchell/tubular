function [ff, vertices, faceMemberIDs] = ...
    polygonsToPatchTriangles3D(vertices, cellIDs, options)
%polygonsToPatchTriangles3D(vertices, faceIDs, options)
% Prepare vertices and faceIDs for plotting via
%         patch('Faces',ff,'Vertices',vv,...
%             'FaceVertexCData',colorV,'FaceColor','flat', ...
%             'Edgecolor', 'none');
%
% This is guaranteed to look ok only if polygons are convex wrt their
% centroids.
%
%
% Parameters
% ----------
% vertices : (#pts in all polygons x D float) or (#vertices x D float)
%   If ~vertexBasedTiling, vertices(i, :) is a vertex of polygon
%   faceIDs(i).
%   If vertexBasedTiling, vertices(i, :) gives a vertex in Ddims that is
%   shared by multiple polygons
% cellIDs : (#pts in all polygons x 1 int array) or (#polygons x 1 cell of
%       int array with vertex indices)
%   If ~vertexBasedTiling, cellIDs(i) states which polygon vertices(i) is 
%   a member of. 
%   If vertexBasedTiling, cellIDs{i} gives the indices of vertices which
%   compose the polygon i when traversed.
% options : optional struct with fields
%   vertexBasedTiling (default=false)
%   centroids : (optional, otherwise computed)
%       centroids of each polygon
%
%
%
% Default options
vertexBasedTiling = false ;

% Unpack options
if isfield(options, 'vertexBasedTiling')
    vertexBasedTiling = options.vertexBasedTiling ;
end
if isfield(options, 'centroids')
    centroids = options.centroids ;
end
if isfield(options, 'keep')
    keep = options.keep ;
end


% Prepare face indexing
if vertexBasedTiling
    % vertices are shared by multiple polygons
    nCells = length(cellIDs) ;
    ff = zeros(nCells * 7, 3) ;
    faceMemberIDs = nan(nCells * 7, 1) ;
    if ~exist('keep', 'var')
        keep = 1:nCells ;
    end
else    
    assert(size(cellIDs, 1) == size(vertices, 1)) ;
    
    % get total # polygon vertices
    ntris = size(cellIDs, 1) ;
    ff = zeros(ntris, 3) ;
    faceMemberIDs = nan(ntris, 1) ;
    nCells = max(cellIDs) ;
    
    % prepare centroids -- estimate as mean vertex position for each cell
    if ~exist('centroids', 'var')
        centroids = nan(max(cellIDs), size(vertices, 2) ) ;
        for qq = unique(cellIDs)        
            centroids(qq, :) = mean(vertices(cellIDs==qq, :)) ;
        end
    end
    if ~exist('keep', 'var')
        keep = 1:nCells ;
    end
end

nVertices = size(vertices, 1) ;
dmyk = 1 ;
if ~vertexBasedTiling
    % Corrected segmenation is boundary-based, not vertex-based tiling
    % (so each cell is a polygon whose boundary we know)
    for cid = 1:nCells
        if ismember(cid, keep)
            face = find(cellIDs == cid) ;
            if ~isempty(face)
                for vid = 1:length(face)
                    if vid < length(face)
                        addface = [face(vid), face(vid+1), nVertices + cid] ;
                    else
                        addface = [face(vid), face(1), nVertices + cid] ;
                    end
                    ff(dmyk, :) = addface ;
                    faceMemberIDs(dmyk) = cid ;
                    dmyk = dmyk + 1 ;
                end
            end
        end
    end
else
    % Uncorrected segmentation is vertex-based tiling
    % with vertices shared by cells known by their faces.
    for cid = 1:nCells
        if ismember(cid, keep)
            face = cellIDs{cid} ;
            if ~isempty(face)
                for vid = 1:length(face)
                    if vid < length(face)
                        addface = [face(vid), face(vid+1), nVertices + cid] ;
                    else
                        addface = [face(vid), face(1), nVertices + cid] ;
                    end
                    ff(dmyk, :) = addface ;
                    faceMemberIDs(dmyk) = cid ;
                    dmyk = dmyk + 1 ;
                end
            end
        end
    end
end
ff = ff(1:dmyk-1, :) ;
faceMemberIDs = faceMemberIDs(1:dmyk-1) ;

% Augment vertices with centroids for patch triangles
vertices = [ vertices; centroids] ;

