function polygons = Cdat2polygons(Cdat, xy, BL, NL, options)
% Convert a struct Cdat of cell information to a list of polygons with
% ordered vertices, one for each cell.
%
% Note that this is, in the general case, an NP hard problem. 
% If there are both cyclic components (ie triangles) to
% the graph of bonds, then we can't used a mixed-sign weight of shortest
% path in a graph. Could use bellmanFordShortestPaths() to find
% negative cycles of a digraph, but here I find all possible paths
% and take the longest. 
%
% Todo: in walking method, initially choose bond direction that has largest
% angle of all possibilities. Handle case where the initial bond direction
% is wrong (is that possible?) so that we always pick the right angle,
% whether max or min.
% Todo: Consider removing bond self-intersections before proceeding.
%
% Parameters
% ----------
% Cdat : cell data struct with fields
%   
% xy : #vertices x 2 numeric array of vertex positions
% BL : #bonds x 2 int array
% options : optional struct with optional fields
%   roi : [xmin, xmax; ymin, ymax] or polygonal region
%       region of interest in which to accept cell centroids
%   method : 'graph' or 'walking' (string specifier, default='graph')
%       method to use to extract polygons
%   debug : bool (default=false)
%       debug the current function acting on inputs (more visualization &
%       output)
%   maxNumEdges : int (default=20)
%       maximum number of edges a valid cell can have
%   pausetime : numeric value (default=0.3) ;   
%       pausing time in seconds for viewing badly behaved cell polygons
%
%
% Returns
% -------
% polygons : #cells x 1 struct of cells
%   polygons{i} gives the ordered n vertex indices of cell i, so that the
%   handedness curls along the z dimension of the 2d tesselation of cells
% 
% NPMitchell 2021

% Default options & unpacking
method = 'graph' ;  % 'walking' or 'graph'
debug = false ;     % whether to debug the algorithm
maxNumEdges = 20 ;  % maximum number of edges a valid cell can have
pausetime = 0.3 ;   % pausing time in seconds for viewing badly behaved cell polygons
if nargin < 5
    options = struct() ;
end

if isfield(options, 'method')
    method = options.method ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'maxNumEdges')
    maxNumEdges = options.maxNumEdges ;
end
if isfield(options, 'pausetime')
    pausetime = options.pausetime ;
end
if isfield(options, 'roi')
    roi = options.roi ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot it to check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if debug
%     disp('plotting for debug')
%     clf
%     plotBonds(xy, BL)
% end

% Compute bond lengths for graph approach        
bL = vecnorm(xy(BL(:, 2), :) - xy(BL(:, 1), :), 2, 2) ;
        
for cid = 1:length(Cdat)
    if mod(cid, 100) == 1
        disp(['Constructing polygon for cell ' num2str(cid) '/' num2str(length(Cdat))])
    end
    
    % what are the vertices involved in polygon cid? (ie this cell)
    ci = Cdat(cid).nverts ;

    % % Preview the polygon to close
    % if debug 
    %     plot(xy(ci, 1), xy(ci, 2), '.')
    % end
    
    if numel(roi) == 4
        inROI = Cdat(cid).centroid.coord(1) > roi(1, 1) & ...
            Cdat(cid).centroid.coord(1) < roi(1, 2) & ...
            Cdat(cid).centroid.coord(2) > roi(2,1) & ...
            Cdat(cid).centroid.coord(2) < roi(2,2) ;  
    else
        error('handle more general roi here')
    end
    
    % Get path along connections    
    if length(ci) < maxNumEdges && length(ci) > 2 && inROI
        if strcmpi(method, 'walking')
            % start with the first one
            [rows, cols] = find(BL == ci(1)) ;
            % find which bond to move along
            candidates = intersect(setdiff(unique(BL(rows, :)), ci(1)), ci) ;
            % walk along one of these
            nextpt = candidates(1) ;
            pg = [ci(1), nextpt] ;
            closed = false ;
            while ~closed

                % next point candidates
                [rows, cols] = find(BL == nextpt) ;
                cand0 = intersect(setdiff(unique(BL(rows, :)), nextpt), ci) ;
                candidates = setdiff(cand0, pg) ;

                if debug
                    disp('stepping...')
                    clf
                    plotBonds(xy, BL)
                    plot(xy(ci, 1), xy(ci, 2), '.')
                    plot(xy(pg(1:end-1), 1), xy(pg(1:end-1), 2), '^')
                    hold on;
                    plot(xy(candidates, 1), xy(candidates, 2), 'o')
                    plot(xy(nextpt, 1), xy(nextpt, 2), 's')
                    % xlim([0, 140])
                    % ylim([700, 900])
                    pause(0.1)
                end

                if length(candidates) == 1
                    nextpt = candidates(1) ;
                    pg = [pg, nextpt] ;
                elseif length(candidates) > 1
                    % the neighboring polygons include at least two of the 
                    % current polygon vertices

                    % Compute angle to each candidate and accept maximum
                    u = xy(nextpt, :) - xy(pg(end-1), :) ;
                    vcands = xy(candidates, :) - xy(nextpt, :) ;
                    
                    % Could compute angle using 3d vector
                    % u(:, 3) = 0 ;
                    % vcands(:, 3) = 0 ;
                    th = zeros(length(candidates), 1) ;
                    for cc = 1:length(candidates)
                        vcand = vcands(cc, :) ;
                        x1 = u(1); y1 = u(2) ;
                        x2 = vcand(1); y2 = vcand(2) ;
                        th(cc) = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
                        % Could compute angle using 3d vector
                        % th(cc) = atan2d(norm(cross(u,vcand)), dot(u,vcand));
                    end
                    % Choosing maximum theta as next id
                    [~, maxid] = max(th) ;
                    nextpt = candidates(maxid) ;
                    pg = [pg, nextpt] ;
                elseif isempty(candidates) 
                    % This is acceptible if the current polygon is about to close
                    candidates = setdiff(cand0, pg(2:end)) ;
                    if candidates(1) == ci(1) && length(candidates) == 1
                        nextpt = ci(1) ;
                    else
                        error('This should not happen. Can we close the polygon?')
                    end
                end

                if nextpt == ci(1) 
                    closed = true ;
                end
            end
        else
            % Graph approach: 
            % Obtain all existing bonds between vertices in ci
            combs = nchoosek(ci, 2) ;
            pairsV = BL(ismember(BL, combs,'rows'), :) ;
            weights0 = bL(ismember(BL, combs,'rows'), :) ;
            [uniqueV, iV, i0] = unique(pairsV(:)) ;
            % pairs0 = pairsV(iV) ;
            % pairsV = pairs0(i0) ;
            pairs0 = reshape(i0, size(pairsV)) ;

            % Sever the first pair that has a degree of 2
            G0 = graph(pairs0(:, 1), pairs0(:, 2), 1) ;
            bID = 1 ;
            severed  = false ;
            clearvars pairs
            while ~severed 
                if bID > size(pairs0, 1)
                    disp('WARNING: did not find a pair of degree 2. Cannot be a polygon')
                    severed = true ;
                else
                    bond2sever = pairs0(bID, :) ;
                    if degree(G0, bond2sever(1)) == 2
                        keep = setdiff(1:size(pairs0, 1), bID) ;
                        pairs = pairs0(keep, :) ;
                        weights = weights0(keep) ;
                        severed = true ;
                    else
                        disp('skipping bond since not deg=2')
                        bID = bID + 1 ;
                    end
                end
            end
            
            % Severed graph (one bond cut) --> if pairs exist, then we have
            % a good graph with one cut
            if exist('pairs', 'var')
                G1 = graph(pairs(:, 1), pairs(:, 2), weights) ;
                % G1 = digraph([pairs(:, 1); pairs(:, 2)], [pairs(:, 2); pairs(:, 1)], [weights; weights]) ;
                % Subplot each graph: one should be ring-like, one should be
                % arc-like
                if debug
                    plot(G0); hold on;
                    plot(G1) ; %, 'EdgeLabel', G1.Edges.Weight)
                    waitfor(gcf)
                end

                startpt = bond2sever(1) ;
                endpt = bond2sever(2) ;        
                % pg0 = shortestpath(G1, startpt, endpt) ;

                % Find all paths between nodes with pathbetweennodes
                % refs: 
                % https://www.mathworks.com/matlabcentral/answers/12283-all-simple-paths-problem
                % https://github.com/kakearney/pathbetweennodes-pkg
                nnode = max(pairs(:));
                nedge = 2*size(pairs, 1);
                adj = sparse([pairs(:, 1); pairs(:, 2)], ...
                    [pairs(:, 2); pairs(:, 1)], ones(nedge,1), nnode, nnode);
                pth = pathbetweennodes(adj, startpt, endpt) ;

                % convert back to vertex ids
                if isempty(pth)
                    skipCell = true ;
                elseif length(pth) == 1
                    pg0 = pth{1} ;
                    skipCell = false ;
                else
                    % Sum up pathlengths
                    % Consider each path
                    plens = zeros(length(pth), 1) ;
                    for pathID = 1:length(pth)
                        ppath = pth{pathID} ;
                        plen = 0 ;
                        Vpath = uniqueV(ppath)' ;
                        % Consider each step in the path, add up Euclidean length
                        for ss = 1:length(ppath)
                            if ss < length(ppath)
                                step = Vpath(ss:ss+1) ;
                            elseif ss == length(ppath)
                                step = [Vpath(ss) Vpath(1)] ;   
                            else
                                error('Path cannot be longer than itself!')
                            end
                            plen = plen + vecnorm(xy(step(2), :) - xy(step(1), :), 2, 2) ;
                        end
                        plens(pathID) = plen ;
                    end
                    % Choose the longest path
                    [~, maxID] = max(plens) ;
                    try
                        % longest length path
                        pg0 = pth{maxID} ;
                        % success! this is a good cell
                        skipCell = false ;
                    catch
                        error('what went wrong with choosing max length path?')
                    end
                end
            else
                try 
                    assert(all(degree(G0, pairs0(:, 1)) > 2) || ...
                        all(degree(G0, pairs0(:, 1)) == 1))
                    disp('Could not cut a bond -- all vertices are deg>2')
                    skipCell = true ;
                catch
                    error('Unclear what is failing. Debug here.')
                end
            end
            
            % Convert back to vertex indices
            % If we have a path, then re-index original vertices for this
            % path
            if ~skipCell
                % we have a path
                try
                    % re-index with original vertex indices rather than
                    % simplified N indices for an N-gon
                    pg = uniqueV(pg0)' ;
                    % This is a good cell
                    skipCell = false ;
                catch
                    disp('found path but could not index...')

                    clf
                    plotBonds(xy, BL)
                    hold on;
                    plot(xy(ci, 1), xy(ci, 2), '.')
                    plot(xy(pairsV(:), 1), xy(pairsV(:), 2), '^')

                    for tid = 1:length(pairsV)
                        plot([xy(pairsV(tid, 1), 1), xy(pairsV(tid, 2), 1)], ...
                            [xy(pairsV(tid, 1), 2), xy(pairsV(tid, 2), 2)], 'k--') ;
                        pause(pausetime)
                        xlim([min(xy(pairsV(:), 1)) - 1, max(xy(pairsV(:), 1)) + 1])
                        ylim([min(xy(pairsV(:), 2)) - 1, max(xy(pairsV(:), 2)) + 1])
                    end
                    disp('Ignoring this cell')
                    skipCell = true ;
                end
            end
        end

        % Check whether angle differences sum to 2pi or -2pi
        
        if ~skipCell
            th = zeros(length(pg), 1) ;
            for cc = 1:length(pg)
                if cc < length(pg) - 1
                    step1 = xy(pg(cc+1), :) - xy(pg(cc), :) ;
                    step2 = xy(pg(cc+2), :) - xy(pg(cc+1), :) ;
                elseif cc < length(pg)
                    step1 = xy(pg(cc+1), :) - xy(pg(cc), :) ;
                    step2 = xy(pg(1), :) - xy(pg(cc+1), :) ;
                elseif cc == length(pg)
                    step1 = xy(pg(1), :) - xy(pg(cc), :) ;
                    step2 = xy(pg(2), :) - xy(pg(1), :) ;
                else
                    error('Bad indexing of polygon pg')
                end
                % step1 = [step1 0] ;
                % step2 = [step2 0] ;
                % th(cc) = atan2(norm(cross(step1,step2)), dot(step1,step2));

                x1 = step1(1) ; y1 = step1(2) ;
                x2 = step2(1) ; y2 = step2(2) ;
                th(cc) = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
            end

            % Path should be closed 
            try
                assert((abs(sum(th)) - 360) < 1e-2)
            catch
                disp(['polygon does not seem closed: sum(theta) = ' num2str(sum(th))])
                clf
                plotBonds(xy, BL)
                hold on;
                plot(xy(ci, 1), xy(ci, 2), '.')
                for tempid = 1:length(pg)
                    plot(xy(pg(tempid), 1), xy(pg(tempid), 2), '^')
                    xlim([min(xy(pg, 1)) - 1, max(xy(pg, 1)) + 1])
                    ylim([min(xy(pg, 2)) - 1, max(xy(pg, 2)) + 1])
                    pause(pausetime)
                end
                skipCell = true ;
                disp('polygon does not seem to be closed!')
            end

            if ~skipCell 
                % Flip if negative
                if sum(th) < 0
                    pg = flipud(pg) ;
                end 

                polygons{cid} = pg ;

                if debug
                    clf
                    plotBonds(xy, BL)
                    hold on;
                    plot(xy(ci, 1), xy(ci, 2), '.')
                    for tempid = 1:length(pg)
                        plot(xy(pg(tempid), 1), xy(pg(tempid), 2), '^')
                        pause(0.1)
                        xlim([min(xy(pg, 1)) - 1, max(xy(pg, 1)) + 1])
                        ylim([min(xy(pg, 2)) - 1, max(xy(pg, 2)) + 1])
                    end
                end
            else
                polygons{cid} = [] ;
            end
        else
            polygons{cid} = [] ;
        end
        
        
    else
        if ~inROI
            disp('Cell not in ROI')    
        else
            disp(['Cell has either < 3 vertices or more than ' num2str(maxNumEdges)])
        end
        polygons{cid} = [] ;
    end
    
    
end

end