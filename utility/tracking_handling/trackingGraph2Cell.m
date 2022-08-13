function tracks = trackingGraph2Cell(GG, timePoints) 
%
% Order track by tidx=1:length(timePoints) (NOTE: not timePoints) as index
%
% Returns
% -------
% tracks : #lineages x 1 cell array of #timepoints x 2 float array
%   non-dividing tracks for each object through time with Node # in GG for
%   reference: tracks{i}(j, :) is [Ux, Uy, lineageID] for the ith object 
%   in the jth timepoint, timePoints(j).
%


nLineages = max(GG.Nodes.LineageID) ;
tracks = cell(nLineages, 1) ;
nodeUsed = false(size(GG.Nodes, 1)) ;

% Unpack each lineage into track
for ii = 1:nLineages
    disp(['Creating lineage ' num2str(ii)])
    trackii = zeros(length(timePoints), 3) ;
    for tidx = 1:length(timePoints)
        disp(['t = ', num2str(timePoints(tidx))])
        % Find all objects at this timepoint from a given lineage ii
        TID = find(GG.Nodes.T == timePoints(tidx));
        LID = find(GG.Nodes.LineageID == ii) ;
        match = intersect(TID, LID) ;
        match = setdiff(match, find(nodeUsed)) ;
        
        % If match is unique add Upix and NodeID to lineage
        if length(match) == 1
            % Find coordinates for this match
            Uxy = GG.Nodes(match, :).UPix ;
            trackii(tidx, :) = [Uxy(1), Uxy(2), match] ;
            nodeUsed(match) = true ;
        elseif length(match) > 1
            disp(['More than one match for lineage ' num2str(ii)])
            % Take the closer lineage in Upix space
            Umatches = GG.Nodes(match, :).UPix ;
            Uprev = GG.Nodes(trackii(tidx-1, 3), :).UPix ;
            if isnan(trackii(tidx-1, 3))
                disp([ 'Previous match was missing, ' , ...
                    'choosing first match for lineage ', ...
                    'continuation arbitrarily'])
                match = match(1) ;
            else
                disp(' -> Find nearest match for lineage continuation')
                [~, minID] = min(vecnorm(Umatches - Uprev, 2, 2)) ;
                match = match(minID) ;
            end
        end
        
        % Populate track with match    
        if isempty(match)
            disp(['No match for lineage ' num2str(ii) ' at t=' num2str(timePoints(tidx))])
            trackii(tidx, :) = [NaN, NaN, NaN] ;
        else
            % Find coordinates for this match
            Uxy = GG.Nodes(match, :).UPix ;
            trackii(tidx, :) = [Uxy(1), Uxy(2), match] ;
            nodeUsed(match) = true ;
        end
    
    end
    
    % Append track to track collection cell array
    tracks{ii} = trackii  ;

end
    
    
    
    
    
    