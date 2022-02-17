function [path, d, edgepath] = nearestpath(G, vertexIn, s, t, p, MethodFlag)
% NEARESTPATH Compute shortest path between two nodes
%
%   PATH = NEARESTPATH(G, S, T, P, options) computes the shortest path 
%   starting at node S and ending at node T that is piecwise nearest to 
%   another curve in 3D space, P. 
%
%   If the graph is weighted (that is, G.Edges
%   contains a Weight variable), then those weights are used as the
%   distances along the edges in the graph. Otherwise, all edge distances
%   are taken to be 1. PATH contains all nodes on the shortest path. PATH
%   is a vector of numeric node IDs if S and T are node indices, a string 
%   vector if S and T are string node names, or a cell array of character 
%   vectors if S and are character vector node names. If node T is
%   unreachable from node S, then PATH is empty.
%
%   [PATH,D] = NEARESTPATH(G,S,T,PATH) also returns the length of the shortest
%   path, D. If node T is unreachable from node S, then D is Inf.
%
%   [PATH,D,EDGEPATH] = NEARESTPATH(G,S,T) also returns the edges on the
%   path from node S to node T.
%
%   [PATH,D] = NEARESTPATH(...,'Method',METHODFLAG) optionally specifies
%   the method to use in computing the shortest path.
%   METHODFLAG can be:
%
%         'auto'  -  Uses 'unweighted' if no weights are set, and
%                    'positive' otherwise. This method is the default.
%
%   'unweighted'  -  Treats all edge weights as 1.
%
%     'positive'  -  Requires all edge weights to be positive.
%
%   Example:
%       % Create and plot a graph. Compute the shortest path from node 7
%       % to node 8.
%       s = [1 1 2 3 3 4 4 6 6 7 8 7 5];
%       t = [2 3 4 4 5 5 6 1 8 1 3 2 8];
%       G = graph(s,t);
%       plot(G)
%       path = nearestpath(G,7,8)
%
%   Example:
%       % Create and plot a weighted graph. Compute and highlight
%       % the shortest path from node 3 to node 8.
%       s = [1 3 1 2 2 6 6 7 7 3  3 9  9  4 12 11 11  8];
%       t = [2 1 4 5 6 7 8 5 8 9 10 5 10 11  4 10 12 12];
%       weights = [10 10 10 10 10 1 1 1 1 1 1 1 1 1 1 1 1 1];
%       G = graph(s,t,weights);
%       p = plot(G,'EdgeLabel',G.Edges.Weight);
%       path = nearestpath(G,3,8)
%       highlight(p, path,'EdgeColor','red')
%
%   See also SHORTESTPATH, SHORTESTPATHTREE, DISTANCES, DIGRAPH/SHORTESTPATH
% 
%   NPMitchell 2019

if nargin < 6
    MethodFlag = 'auto' ;
end

paths = cell(length(p) - 1, 1) ;
if nargout > 2
    edgepaths = cell(length(p) - 1, 1) ;
end
if nargout > 1
    dists = cell(length(p) - 1, 1) ;
end

for ii = 1:(length(p)-1)
    % Obtain starting point
    if ii == 1
        % This is the starting point of the entire path
        si = s;
    else
        % Find a node in G nearest p(ii, :)
        si = pointMatch(p(ii, :), vertexIn) ;
    end
    
    % Obtain ending point
    if ii == length(p) - 1
        % This is the termination point of the final path
        ti = t ;
    else
        % Find a node in G nearest p(ii, :)
        ti = pointMatch(p(ii + 1, :), vertexIn) ;
    end
    
    % Compute this mini path
    varargin
    if nargout > 2
        [paths{ii}, dists{ii}, edgepaths{ii}] = shortestpath(G, si, ti, 'Method', MethodFlag) ;
    elseif nargout > 1
        [paths{ii}, dists{ii}] = shortestpath(G, si, ti, 'Method', MethodFlag) ;
    else
        paths{ii} = shortestpath(G, si, ti, 'Method', MethodFlag) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collate the paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preallocate
path = [] ;
if nargout > 1
    d = [] ;
    if nargout > 2
        edgepath = [] ;
    end
end
% collect paths
for ii = 1:length(paths)
    if ii == 1
        path = paths{ii} ;
    else
        path = [path(1:end-1), paths{ii}] ;
    end
    
    % Also collate other output if desired
    if nargout > 1
        d = d + dists{ii} ;
        if nargout > 2
            edgepath = [edgepath, edgepaths{ii}] ;
        end
    end
end

path = path' ;

