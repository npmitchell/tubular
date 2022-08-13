function GG = unpackManualIlastikGroundTruthH5(h5fileBase, timePoints, ...
    imOutDir, rawImFileBase, Options)
% Convert ilastik manual tracking output (series of H5 files, one for each
% timepoint) into a digraph with nodes corresponding to centroids of
% tracked segmented objects and edges corresponding to associations between
% nodes. Segment labels are string enumerating the cells by an ID.
%
% Parameters
% ----------
% h5fileBase : str
%   full path with %06d for tidx-stamp (time index with zero 
%   indexing as of 2021) stamp for output h5 files from ilastik manual
%   tracking step.
% timePoints : N x 1 int array
%   timePoints associate with the series of tracked frames such that
%   timePoints(tidx) is the timestamp of the tidx'th frame, where (as of
%   2021) tidx is zero-indexed and sprintf(h5fileBase, tidx) is the path to
%   the tidx'th frame's tracking data.
% imOutDir : str or empty (optional)
%   if nonempty, output directory for color image of tracking result
%   illustrating tracks by color on top of raw data. 
%   As of 2021, discontinuous tracks are not followed and merged tracks are 
%   not separated in a differentiated fashion.
% rawImFileBase : str or empty (optional if plotting images is desired)
%   full path base for tidx-stamp output image for tracking result overlaid
%   on raw data
% Options : struct with fields:
%   allowSplitting : bool (default = false)
%   allowMerging : bool (default = false) 
%   maxN_for_plot : int (default = 5000)
%       maximum number of distinct colors in the output images
%
% Returns
% --------
% GG : MATLAB digraph 
%   digraph with each node corresponding to a XY location (in pixels)
%   initially associated here with a segmented object that is tracked.
%   Edges connect tracked objects between subsequent timeframes, including
%   mergers and splits, but not (as of 2021) across multiple timeframes if
%   an object disappears or is occluded.
%   
% NPMitchell 2021

if ~strcmp(h5fileBase(end-2:end), '.h5')
    h5fileBase = [h5fileBase '.h5'] ;
end

if nargin < 2
    timePoints = 1:length(dir(h5fileBase)) ;
    try
        assert(length(timePoints) > 0)
    catch
        error('No matching files found')
    end
end
if nargin < 3
    imOutDir = [] ;
end

% Process Options
maxN_for_plot = 5000 ;
allowSplitting = false ;
allowMerging = false ;
if nargin < 5
    if isfield(Options, 'maxN_for_plot')
        maxN_for_plot = Options.maxN_for_plot ;
    end
    if isfield(Options, 'allowSplitting')
        allowSplitting = Options.allowSplitting ;
    end
    if isfield(Options, 'allowMerging')
        allowMerging = Options.allowMerging ;
    end
end
if ~isempty(imOutDir)
    jets = jet(maxN_for_plot) ;
    jetshuffle = jets(randperm(maxN_for_plot), :) ;
end

mfmoveList = [] ;

% countN1 is the "base" index of all previously-accounted-for nodes
% up until the PREVIOUS timepoint
countN1 = 0 ;
for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    disp(['condidering t = ' num2str(tp)])
    tracksfn_tt = sprintf(h5fileBase, tidx-1) ;
    
    im = h5read(tracksfn_tt, '/segmentation/labels') ;
    im = squeeze(im(1, :, :)) ;
    
    % Query all tracked objects among those segmented
    if tidx > 1
        moves = h5read(tracksfn_tt, '/tracking/Moves')' ;
        try
            appr = h5read(tracksfn_tt, '/tracking/Appearances')' ;
            disp('Detected appearances')
        catch
            appr = [] ;
        end
        % cell label disappeared in current file
        try
            disa = h5read(tracksfn_tt, '/tracking/Disappearances')' ;
        catch
            disa = [] ;
        end
        % ancestor (previous file) | descendant 1 (current file) | descendant 2 (current file)
        try
            splt = h5read(tracksfn_tt, '/tracking/Splits')' ;
        catch
            splt = [] ;
        end
        % descendant/child (in current file) | number of objects in child
        try
            merg = h5read(tracksfn_tt, '/tracking/Mergers')' ;
        catch
            merg = [] ;
        end
        
        try 
            mfmoves = h5read(tracksfn_tt, '/tracking/MultiFrameMoves') ;
            disp('Found multi-frame moves')
            if ~isempty(mfmoveList)
                mfmoveList = [mfmoveList; cat(2, mfmoves', tp * ones(size(mfmoves', 1), 1))];
            else
                mfmoveList = cat(2, mfmoves', tp * ones(size(mfmoves', 1), 1)) ;
            end
        catch
            mfmoves = [] ;
        end
        
        regp = regionprops(im, 'Centroid') ;
        com1 = vertcat(regp.Centroid) ;
        
        %%%%%%%%%%%%%%%%%%%%
        % Thread tracks
        %%%%%%%%%%%%%%%%%%%%
        
        % SCHEMATIC OF DIGRAPH
        % ......... nodes for t<t-2
        % ||||||||| <- edges added before last time
        % ......... nodes [1,countN1] ( for time t-2 )
        % | / | | | <- edges added last time
        % ..  . . . nodes [countN1+1, count0] (for time t-1)
        % | \ |  \  <-- these are the edges we are adding in newEdges
        % .  ..   . nodes [count0+1, maxNodeNumber] (for time t)
        
        
        % First timepoint is special since no endpoints into it. Handle
        % here --> define allCellTimes
        if tidx == 2
            % Assign FIRST TIMEPOINT track endpts
            spt = com0(moves(:, 1), :) ;
            allCellCenters = spt ;
            numCellsT = size(spt, 1) ;
            allCellTimes = repmat(timePoints(1), numCellsT, 1) ;
            allCellTimes0 = allCellTimes ;
            % Local segment indices of first timepoint are prevEndptSegIdx
            prevEndptSegIdx = moves(:, 1) ;
            prevSegmentLabels = (1:length(prevEndptSegIdx))' ;
            allSegmentLabelsArr = prevSegmentLabels ;
        end
        % sptIdx -> indices into segment IDs of previous timepoint
        sptIdx = [moves(:, 1),  -1 * ones(size(appr)) ] ; 
        
        % All current IDs from moves, appearances, etc
        % eptIdx -> indices into segment IDs of current timepoint
        eptIdx = [moves(:, 2), appr] ;
                
        % Prevent mergers if desired
        % if length(unique(eptIdx)) < length(eptIdx) && ~allowMerging
        %     disp('Cutting out mergers...')
        %     [ eptIdx, index4unique, ~] = unique(eptIdx, 'stable') ;
        %    sptIdx = sptIdx(index4unique) ;
        % end
        
        % Remove any appearances from the startpt list
        sptIdx(sptIdx < 0) = [] ;
        
        % convert IDs into COM --> Upix
        ept = com1(eptIdx, :) ;
        numCellsT = size(ept, 1) ;
                
        % count0 is the "base" index of all previously-accounted-for nodes
        % count0 counts up until the current timepoint
        count0 = length(allCellTimes) ;
                
        % Create edges of graph connecting the global Node # of start
        % points in previous timepoint to global Node # of the end points
        % in current timepoint
        globalIdxStart = zeros(size(sptIdx)) ;
        globalIdxEnd = zeros(size(sptIdx)) ;
        donePrevEndptIdx = [] ;
        newTrackId = max(allSegmentLabelsArr(:)) + 1 ;
        newSegmentLabels = zeros(size(sptIdx)) ;
                
        for qq = 1:length(sptIdx) 
            % initialize the global IDx of the endpoint
            currentGlobalEndPtIdx = count0 + qq ;
            sqq = sptIdx(qq) ;
            % Find spt Index in previous endpoints Index in order to just
            % keep the segments that were used as Nodes
            % CONVERT SEGMENTATION IDx TO NODE IDx
            match = find(prevEndptSegIdx == sqq) ;
            if length(match) == 1
                % There is a single match, so connect the current idx to
                % this single previous endpoint (now start point)
                
                % This is the node index: number of nodes before previous
                % timepoint (countN1) + index of node in previous 
                % timepoint's endpoint node list (match).
                globalIdxStart(qq) = match + countN1 ;
                globalIdxEnd(qq) = currentGlobalEndPtIdx ;
                try
                    newSegmentLabels(qq) = prevSegmentLabels(match) ;
                catch
                    error('Could not index into prevSegmentLabels')
                end
                
                
            elseif ~isempty(match) && ~all(ismember(match, donePrevEndptIdx))
                % There are multiple matches --> the current idx split from
                % a number of matches. 
                match2keep = ~ismember(match, donePrevEndptIdx) ;
                match = match(match2keep) ;
                if ~isempty(match)
                    % Mark that we have used this previous endpoint already
                    donePrevEndptIdx = [ donePrevEndptIdx, match(1) ] ;
                    newSegmentLabels(qq) = prevSegmentLabels(match(1)) ;
                    
                    % update current Global endpoint index
                    globalIdxStart(qq) = match(1) + countN1 ;
                    globalIdxEnd(qq) = currentGlobalEndPtIdx ;
                else
                    error('This should not happen?')
                    if ~allowSplitting
                        disp('No unused parent. Disallowing splitting...')
                        % assign new track label to segment label array
                        newSegmentLabels(qq) = newTrackId ;
                        newTrackId = newTrackId + 1 ;   
                    else
                        disp('No unused parent in sight -- presumably a multi-split event? Creating new track...')
                        % assign new track label to segment label array
                        newSegmentLabels(qq) = newTrackId ;
                        newTrackId = newTrackId + 1 ;
                    end
                end
            else
                if ~isempty(match)
                    disp('No splitting allowed, creating new track...')
                else
                    disp('Could not find start point in previous endpoint list -- likely multi-Frame Move')
                end
                newSegmentLabels(qq) = newTrackId ;
                newTrackId = newTrackId + 1 ;
            end
        end
        if length(eptIdx) > length(sptIdx) 
            % assign newSegmentLabels for appearances
            error('here')
            
        end
        
        % HEre assert that there are no zeros in globalIdxStart
        assert(isequal(find(globalIdxStart == 0), find(globalIdxEnd == 0)) )
        
        % Now we ascribe the global node IDs (time-cell ID)
        assert(length(sptIdx) == length(eptIdx)) ;
        % globalIdxEnd = ((1:size(sptIdx, 1)) + count0)' ;
        newEdges = [globalIdxStart, globalIdxEnd] ;
        
        % Remove any rows that have null matches (not sure why these exist)
        [row, ~] = find(newEdges == 0) ;
        if ~isempty(row)
            disp('Removing rows with null matches in newEdges')
            newEdges(row, :) = [] ;
        end
        
        if tidx == 2
            edgesT = newEdges ;
        else
            edgesT = [edgesT; newEdges] ;
        end
        try
            assert(all(edgesT(:) > 0))
        catch
            error('Some edges are assigned to be zero nodes!')
        end
        
        % Update the global cell center list
        % NOTE: globalIdxEnd is ordered as [sinks; appearances/new tracks]
        assert(size(allCellCenters, 1) <= min(newEdges(:, 2)) - 1)
        allCellCenters = [allCellCenters; ept];
        
        % Update the times in the global list
        allCellTimes0 = allCellTimes ;
        allCellTimes =  [allCellTimes; repmat(tp, numCellsT, 1)];
        
        % Make sure #nodes == #labels == #timepoints
        assert(length(newSegmentLabels) == numCellsT)
        assert(size(ept, 1) == numCellsT)
    
        % check min/max of edge start/endpts
        % Note: this assertion is wrong because we could have had a split
        % that emanated from node globalIdxEnd - 1 and removed.
        % assert(max(globalIdxStart) == min(globalIdxEnd) - 1 - length(appr))
        
        % The minimum index of any source node in the edges that we are
        % adding in this timepoint is larger than countN1, which was the
        % highest index of all nodes before the previous timepoint. New
        % edges connect the previous timepoint (countN1, count0] 
        % to the current timepoint (with nodes  > count0)
        assert(min(newEdges(:, 1)) >= countN1 + 1)
        
        % SCHEMATIC OF DIGRAPH
        % ......... nodes for t<t-2
        % ||||||||| <- edges added before last time
        % ......... nodes [1,countN1] ( for time t-2 )
        % | / | | | <- edges added last time
        % ..  . . . nodes [countN1+1, count0] (for time t-1)
        % | \ |  \  <-- these are the edges we are adding in newEdges
        % .  ..   . nodes [count0+1, maxNodeNumber] (for time t)
        
        % Update the segment labels (track IDs)
        allSegmentLabelsArr = [allSegmentLabelsArr; newSegmentLabels ] ;
        
        % assign current centers of mass to previous array (for next tp)
        com0 = com1 ;
        % current timepoint's endpoint segmentation indices are now old
        prevEndptSegIdx = eptIdx ;
        % current timepoint's segment labels (track IDs) are now old
        prevSegmentLabels = newSegmentLabels ;
        % countN1 is the "base" index of all previously-accounted-for nodes
        % up until the PREVIOUS timepoint
        countN1 = count0 ;
        
        assert(numel(eptIdx) == numel(newSegmentLabels)) ;
        
    else
        regp = regionprops(im, 'Centroid') ;
        com0 = vertcat(regp.Centroid) ;
        disp('Waiting for second tp to assign object tracks in first tp')
        im0 = im ;
        
        % cell label disappeared in current file
        try
            disa = h5read(tracksfn_tt, '/tracking/Disappearances')' ;
        catch
            disa = [] ;
        end
        
        % descendant/child (in current file) | number of objects in child
        try
            merg = h5read(tracksfn_tt, '/tracking/Mergers')' ;
        catch
            merg = [] ;
        end
    end
    
    % Output image of tracks if output image directory specified
    if ~isempty(imOutDir)
        % Plot the tracks as colors on raw image or black canvas
        if tidx > 1
            if tidx == 2
                if ~isempty(rawImFileBase)
                    raw0 = imread(sprintf(rawImFileBase, timePoints(1))) ;
                else
                    raw0 = 0* im0 ;
                end
                im2 = 0*im0 ;
                for qq = 1:length(sptIdx)
                    im2(im0 == sptIdx(qq)) = newSegmentLabels(qq) ;
                end
                if any(newSegmentLabels > maxN_for_plot)
                    im2 = mod(im2, maxN_for_plot) ;
                end
                rgb0 = label2rgb(im2, jetshuffle, [0,0,0]);
                zeroIdx = all(rgb0 == 0, 3) ;
                red = squeeze(rgb0(:, :, 1)) ;
                grn = squeeze(rgb0(:, :, 2)) ;
                blu = squeeze(rgb0(:, :, 3)) ;
                red(zeroIdx) = raw0(zeroIdx) ;
                grn(zeroIdx) = raw0(zeroIdx) ;
                blu(zeroIdx) = raw0(zeroIdx) ;
                rgb0 = cat(3, red, grn, blu) ;

                if ~exist(imOutDir, 'dir')
                    mkdir(imOutDir)
                end
                imwrite(rgb0, fullfile(imOutDir, ...
                    sprintf('tracks_color_%06d.png', timePoints(1))))
            end

            % permit no raw image being passed to function -- black bg
            if ~isempty(rawImFileBase)
                raw = imread(sprintf(rawImFileBase, tp)) ;
            else
                raw = 0 * im ;
            end
            im2 = 0*im ;
            for qq = 1:length(eptIdx)
                im2(im == eptIdx(qq)) = newSegmentLabels(qq) ;
            end
            if any(newSegmentLabels > maxN_for_plot)
                im2 = mod(im2, maxN_for_plot) ;
            end
            
            rgb = label2rgb(im2, jetshuffle, [0,0,0]);
            zeroIdx = all(rgb == 0, 3) ;
            red = squeeze(rgb(:, :, 1)) ;
            grn = squeeze(rgb(:, :, 2)) ;
            blu = squeeze(rgb(:, :, 3)) ;
            red(zeroIdx) = raw(zeroIdx) ;
            grn(zeroIdx) = raw(zeroIdx) ;
            blu(zeroIdx) = raw(zeroIdx) ;
            rgb = cat(3, red, grn, blu) ;

            imwrite(rgb, fullfile(imOutDir, ...
                sprintf('tracks_color_%06d.png', tp)))
        end
    end
end

% Convert segment Label array into cell
nCells = size(allCellTimes, 1) ;
allSegmentLabels = cell(nCells, 1) ;
ndigits = num2str(length(num2str(max(allSegmentLabelsArr)))) ;
for qq = 1:length(allSegmentLabelsArr)
    allSegmentLabels{qq} = sprintf(...
        ['%0' ndigits 'd'], allSegmentLabelsArr(qq)) ;
end
allGenerations = ones([nCells, 1]);

% A table containing node properties with which to construct the digraph
nodeTable = table( allCellCenters, allCellTimes, ...
    allSegmentLabels, allGenerations, ...
    'VariableNames', {'UPix', 'T', 'Segment', 'Generation'} );


% Construct the graph and update node properties
GG = digraph;
GG = GG.addnode(nodeTable);

for eID = 1:size(edgesT, 1)
    
   % Update the edge list in the digraph
   GG = addedge(GG, edgesT(eID, 1), edgesT(eID, 2));
    
end

% 
% % Create edgeTable in an independent fashion
% % Consider each timepoint and find nodes that link same segmentLabels
% for tidx = 1:length(timePOints)
%     tp = timePoints(tidx) ;
%     unique(GG.Nodes.T)
%     
% end



%% Assign Basic Lineage IDs ===============================================

%--------------------------------------------------------------------------
% Determine the 'progenitor' cells which have no parents in the tracking
% graph structure
%--------------------------------------------------------------------------

% Find the in-degree of each node
inDeg = indegree(GG);

assert(isequal(unique(inDeg), [0; 1]), 'Invalid tracking structure');

progCells = find(inDeg == 0);

clear inDeg

%--------------------------------------------------------------------------
% Add a lineage ID field to each node of the tracking graph
%--------------------------------------------------------------------------

lineageID = -ones(size(GG.Nodes,1),1);
lineageID(progCells) = 1:numel(progCells);

GG.Nodes.LineageID = lineageID;

clear lineageID

%--------------------------------------------------------------------------
% Update the lineage ID field of each node using a depth-first search
%--------------------------------------------------------------------------

for i = 1:numel(progCells)
    
    % Get the IDs of all nodes descended from the current progenitor cell
    descNodes = dfsearch(GG, progCells(i), ...
        { 'discovernode', 'edgetonew' } );
    descNodes = descNodes.Node;
    descNodes = descNodes( ~isnan(descNodes) );
    
    descNodes( ismember(descNodes, progCells) ) = [];
    
    % Udpate the lineage ID
    GG.Nodes(descNodes,:).LineageID = repmat(i, numel(descNodes), 1);
    
end

clear descNodes

assert( ~any(GG.Nodes.LineageID < 0), 'Lineage improperly assigned!' );
disp('DONE');