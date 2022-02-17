function [ff, vv, qualityV, orig] = polygons2triangulation(polygons, ...
    vertices, centroids, options)
% [ff, vv, qualityV] = polygons2triangulation(polygons, vertices, centroids, options)
% Convert sufficiently convex polygons into triangular faces for plotting
% as a single patch object
%
% Parameters
% ----------
% polygons : #polygons x 1 or 1 x #polygons cell array of polygon indices
%   polygons{i} is the vertex id list of vertices that make polygon i
% vertices : 
% options : optional struct with fields
%   qualities : #polygons x 1
%   keep : #polygons x 1 bool array
%       keep(i) indicates whether to include polygon i in output 
%       triangulation
%
% Returns
% -------

nCells = length(polygons) ;
keep = 1:nCells ;
has_qualities = false ;
if isfield(options, 'keep')
    keep = options.keep ;
end
if isfield(options, 'qualities')
    qualities = options.qualities ;
    has_qualities = true ;
    quality_is_vector = size(qualities, 2) > 1 ;
elseif isfield(options, 'quality')
    qualities = options.quality ;
    has_qualities = true ;
    quality_is_vector = size(qualities, 2) > 1 ;
end

% Preallocate arrays
nVertices = size(vertices, 1) ;
ff = zeros(nCells * 7, 3) ;
if has_qualities
    qualityV = nan(nCells * 7, 1) ;
end

% Index for stuffing vectors
dmyk = 1 ;
for cid = 1:nCells
    if ismember(cid, keep)
        face = polygons{cid} ;
        if ~isempty(face)
            for vid = 1:length(face)
                if vid < length(face)
                    addface = [face(vid), face(vid+1), nVertices + cid] ;
                else
                    addface = [face(vid), face(1), nVertices + cid] ;
                end
                ff(dmyk, :) = addface ;
                if has_qualities
                    if quality_is_vector
                        qualityV(dmyk, :) = qualities(cid, :) ;
                    else
                        qualityV(dmyk) = qualities(cid) ;
                    end
                end
                dmyk = dmyk + 1 ;
            end
        end
    end
end

% Truncate the face list and add centroids to vertex list
ff = ff(1:dmyk-1, :) ;
vv = [vertices; centroids] ;
if has_qualities
    qualityV = qualityV(1:dmyk-1) ;
end