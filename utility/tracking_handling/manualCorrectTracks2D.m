function [tracks, trackGraph] = manualCorrectTracks2D(currentTracks, ...
    fileBase, timePoints, trackOutfn, tracks2Correct)
% Manually correct object tracking in 2D grayscale image sequence.
%
% Parameters
% ----------
% currentTracks : cell array or empty
%   existing tracks to build off of
% timePoints : length #timepoints numeric array
%   timepoints in which to track
% trackOutfn : str path
%   where to periodically save the results as .mat
% 
% Returns
% -------
% tracks : nTracks x 1 cell array of #timepoints x 2 float arrays
%   positions of each track over time. NaNs are permitted for frames with
%   no object (disappearances)
% trackGraph : digraph object
%   digraph representation of the corrected tracking results with Nodes 
%   and Edges
% 
% NPMitchel 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script version options
% subdir = 'muscle_normalShiftn10_p05_n50_s1p00_lambda0p005_maxProj';
% imDir = fullfile('./', subdir, 'muscle_imagestack_LUT') ;
% trackOutfn = fullfile('./muscle_tracks.mat') ;
% 
% timePoints = 1:60 ;
% nTracks = 300 ;
% fileBase = fullfile(imDir, 'Time_%06d_c1_stab_pbspsm_LUT.tif') ;
% 
% %--------
% % inserted here
% load(fullfile(QS.dir.tracking, './endoderm_tracks.mat', 'tracks')
% %--------
% 
% currentTracks = tracks ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Unpack inputs
pausetime = 0.5 ;
tracks = currentTracks ;
bluecolor = [0 , 0.4470, 0.7410 ];
orange = [ 0.8500,  0.3250 , 0.0980 ];
green = [ 0.2, 0.9, 0.2] ;
lwidth = 3 ;
markerSize = 20 ;
nTracks = length(tracks) ;
recap = false ;

%
if nargin < 5
    tracks2Correct = 1:nTracks ;
end

%% Consider each object
% tracks2Correct = [101:200]  
nTracks = length(tracks) ;
fig = [] ;
figCount = 0 ;      
for dmyii = 1:length(tracks2Correct) 
    ii = tracks2Correct(dmyii) ;
    close all
    % Consider each timepoint
    tidx = 1 ;
    Xlim = [];  
    Ylim = [] ;
    trackii = currentTracks{ii} ;
    keepTracking = true ;
    while keepTracking
        tp = timePoints(tidx) ;

        figCount = figCount + 1 ;
        
        % Initialize the Figure
        if figCount > 20 
            close all
            figCount = 0 ;
        end
        if isempty(fig) || ~ishandle(fig)
            fig = figure('units', 'normalized', 'outerposition', [0.5 0 0.5 1]);
        end

         if isempty(Xlim)
            if ~isnan(trackii(1, 1))
                Xlim = trackii(1, 1) + [-200, 200] ;
                Ylim = trackii(1, 2) + [-100, 100] ;
            end
        end

        [ax0, ax1, ax2] = ...
            plotCurrentTrackForCorrection(ii, tidx, timePoints, fileBase, trackii, ...
            currentTracks, Xlim, Ylim) ;

        % Now acquire or get more info about nearby timepoints
        msg = 'Press <space>: accept, e: edit, p: play, o: play from start, a: t-1, s: t+1' ;
        fprintf(msg)
        sgtitle(msg)
        pause
        currkey=get(gcf, 'CurrentKey'); 


        while ismember(currkey, {'d', 'g', 'p', 'o', 'a', 's', 'j'}) 

            switch currkey
                case {'d'}
                    keepTracking = false ;
                    currkey = 'space' ;
                case {'g'}
                    tidx = input('Go to timepoint index: ') ;
                    if tidx > length(timePoints) 
                        tidx = length(timePoints) ;
                    elseif tidx < 1 
                        tidx = 1 ;
                    end
                    currkey = 'x' ;
                case {'r'}
                    % Rapid click for 5 frames
                    Xlim = get(ax2, 'XLim');
                    Ylim = get(ax2, 'YLim');

                    nRapid = 5 ;
                    axes(ax0) ; cla ;
                    axes(ax1) ; cla ;
                    axes(ax2) ; cla ;
                    qq = 0 ;

                    while tidx < length(timePoints) + 1 && qq < nRapid
                        [ax0, minX, maxX, minY, maxY] = rapidClickCurrentTrack(ii, ...
                            tidx, timePoints, fileBase, trackii, ...
                            tracks, Xlim, Ylim) ;
                        msg = 'Click with crosshairs: acquire / <escape>: no detection' ;
                        sgtitle(msg)
                        disp(msg)
                        [xx,yy] = ginput(1) ;
                        currkey=get(gcf, 'CurrentKey'); 
                        if strcmpi(currkey, 'escape') || strcmpi(currkey, 'backspace')
                            try
                                trackii(tidx, :) = [NaN, NaN, 0] ;
                            catch
                                trackii(:, 3) = NaN ;
                                trackii(tidx, :) = [NaN, NaN, 0] ;
                            end
                        else
                            try
                                trackii(tidx, :) = [xx+minX,yy+minY, 0] ;
                            catch
                                trackii(:, 3) = NaN ;
                                trackii(tidx, :) = [xx+minX,yy+minY, 0] ;
                            end
                        end
                        tidx = tidx + 1 ;
                        qq = qq + 1 ;
                    end
                    % Xlim = Xlim + [-20, 20] ;
                    % Ylim = Ylim + [-20, 20] ;
                    [ax0, ax1, ax2] = plotCurrentTrack(ii, tidx, timePoints, fileBase, trackii, ...
                        currentTracks, Xlim, Ylim) ;

                    anyEdits = true ;
                    currkey = 'x' ;
                case {'a', 's'}
                    % Decrememt/Increment timepoint
                    if strcmpi(currkey, 'a')
                        tidx = tidx - 1;
                        if tidx < 1
                            tidx = length(timePoints) ;
                        end
                    elseif strcmpi(currkey, 's')
                        tidx = tidx + 1;
                        
                        if tidx > length(timePoints)
                            tidx = 1 ;
                        end
                    end
                    [ax0, ax1, ax2] = ...
                        plotCurrentTrackForCorrection(ii, tidx, timePoints, fileBase, trackii, ...
                        currentTracks, Xlim, Ylim) ;

                    % replay or move on to acquire
                    fprintf('Press any key to acquire or p to replay')
                    pause
                    Xlim = get(ax2, 'XLim');
                    Ylim = get(ax2, 'YLim');
                    currkey=get(gcf,'CurrentKey'); 
                    
                case {'j'}
                    msg = 'Joining track to another' ;
                    disp(msg)
                    sgtitle(msg)

                    % Look for nearby nuclei
                    nearby = zeros(nTracks, 2) ;
                    for trackID = 1:nTracks
                        nearby(trackID, :) = tracks{trackID}(tidx, 1:2) ;
                    end
                    selfXY = nearby(ii, :) ;
                    if any(isnan(selfXY))
                        selfXY = tracks{ii}(tidx-1, 1:2) ;
                    end
                    farAway = 1000 ;
                    dists = vecnorm(nearby - selfXY, 2, 2) ;
                    dists(ii) = farAway ;
                    [~, minID] = nanmin(dists) ;
                    
                    % Show gluing procedure ---------------------------
                    
                    % Plot current lineage in orange
                    im = imread(sprintf(fileBase, tp));
                    axes(ax2)
                    imshow(im) ;
                    hold on;
                    plot(trackii(tidx-1, 1), ...
                        trackii(tidx-1, 2), 'o', ...
                        'markerSize', markerSize, ...
                        'lineWidth', lwidth, 'color', orange)
                    plot(trackii(1:tidx-1, 1), ...
                        trackii(1:tidx-1, 2), '-', ...
                        'markerSize', markerSize, ...
                        'lineWidth', lwidth, 'color', orange)
                    
                    plot(tracks{minID}(tidx, 1), ...
                        tracks{minID}(tidx, 2), 'o', ...
                        'markerSize', markerSize, ...
                        'lineWidth', lwidth, 'color', green)
                    hold on;
                    plot(tracks{minID}(tidx:end, 1), ...
                        tracks{minID}(tidx:end, 2), '-', ...
                        'markerSize', markerSize, ...
                        'lineWidth', lwidth, 'color', green)
                    set(gca, 'XLim', Xlim , 'YLim', Ylim);
                    
                    msg = ['Join to track ' num2str(minID) '? <space>: accept, e: edit track to join'] ;
                    disp(msg)
                    sgtitle(msg)
                    pause
                    currkey=get(gcf,'CurrentKey'); 
                    if strcmpi(currkey, 'space')
                        % Accept this gluing
                        trackii(tidx:end, :) = tracks{minID}(tidx:end, :) ;
                    elseif strcmpi(currkey, 'e')
                        % Click on another nuclei to glue to that track
                        axes(ax2)
                        imshow(im) ;
                        hold on;
                        plot(nearby(trackID, 1), nearby(trackID, 2), 's')
                        
                        msg = 'Click on another nuclei to glue to its track' ;
                        sgtitle(msg)
                        disp(msg)
                        pause
                        [nx,ny] = ginput(1) ;
                        dists = vecnorm(nearby - [nx,ny], 2, 2) ;
                        dists(ii) = farAway ;
                        [~, minID] = nanmin(dists) ;                       
                        
                        % Glue track
                        trackii(tidx:end, :) = tracks{minID}(tidx:end, :) ;
                    else
                        msg = 'bad key press, ignoring...' ;
                        disp(msg)
                        sgtitle(msg)
                    end
                    
                    
                case {'p', 'o'}
                    % Play nearby timepoints to help identify cell center                
                    Xlim = get(gca, 'XLim');
                    Ylim = get(gca, 'YLim');

                    if strcmp(currkey, 'p')
                        tidx2play = max(1, tidx-4):min(tidx+1, length(timePoints)) ;
                    elseif strcmp(currkey, 'o')
                        tidx2play = 1:tidx ;
                    end
                    times2play = timePoints(tidx2play) ;
                    axes(ax2) ;
                    for kk = 1:length(times2play)
                        ttemp = times2play(kk) ;
                        im1 = imread(sprintf(fileBase, ttemp));

                        cla
                        imshow(im1) ;
                        hold on;
                        % Plot other lineages in blue
                        otherIDs = setdiff(1:nTracks, [ii]) ;
                        for trackID = 1:otherIDs
                            plot(tracks{trackID}(tidx2play(kk), 1), ...
                                tracks{trackID}(tidx2play(kk), 2), 'o', ...
                                'markerSize', markerSize, ...
                                'lineWidth', lwidth, 'color', bluecolor)
                        end
                        % Plot current lineage in orange
                        plot(trackii(tidx2play(kk), 1), ...
                            trackii(tidx2play(kk), 2), 'o', ...
                            'markerSize', markerSize, ...
                            'lineWidth', lwidth, 'color', orange)
                        if tidx > 1
                            set(gca, 'XLim', Xlim , 'YLim', Ylim);
                        end
                        title(['playing track ' num2str(ii) ': t=' num2str(ttemp)])

                        pause(pausetime)
                    end

                    % Re-gain original image onto axis
                    im = imread(sprintf(fileBase, tidx));
                    cla()
                    imshow(im) ;
                    title(['Track ' num2str(ii) ': t=' num2str(tp)])
                    set(gca, 'XLim', Xlim , 'YLim', Ylim);

                    % replay or move on to acquire
                    fprintf('Press any key to acquire or p to replay')
                    pause
                    currkey=get(gcf,'CurrentKey'); 
                otherwise
                    % disp('done with adjustments')

            end  % end of the switch

        end % end of the while loop

        if strcmpi(currkey, 'e')
            % Acquire the XY coordinate for this timepoint
            msg = 'EDIT POSITION: <space> to click with crosshairs: acquire / <escape>: no detection' ;
            sgtitle(msg)
            disp(msg)
            pause
            [xx,yy] = ginput(1) ;
            currkey=get(gcf,'CurrentKey'); 
            if strcmpi(currkey, 'escape') || strcmpi(currkey, 'backspace')
                trackii(tidx, :) = [NaN, NaN, 0] ;
            else
                trackii(tidx, :) = [xx,yy, trackii(tidx, 3)] ;
            end
        else
            msg = 'position accepted' ;
            sgtitle(msg)
            disp(msg)                
            % pause(0.1)
        end

        Xlim = get(gca, 'XLim');
        Ylim = get(gca, 'YLim');

        % prep for next timepoint
        % Save tracks so far
        tracks{ii} = trackii ;
        tidx = tidx + 1 ;


        % Check if we should exit
        if tidx > length(timePoints) 
            msg = 'All done with detections? [y/n]' ;
            sgtitle(msg)
            disp(msg)
            okstring = input(msg, 's') ;
            notOk = contains(lower(okstring), 'n') ;
                       
            keepTracking = (tidx < length(timePoints) + 1) || notOk ;
            tidx = length(timePoints) ;
        end
    end
    tracks{ii} = trackii ;
    currentTracks = tracks ;
    
    %% Save tracks
    inputString = input('Save tracks [y/n]?', 's') ;
    if contains(lower(inputString), 'y') 
        disp(['saving completed tracks to : ' trackOutfn])
        save(trackOutfn, 'tracks')
    end
    close all
    
    %% Instant replay
    if recap
        if ~ishandle(fig)
            fig = figure('units', 'normalized', 'outerposition', [0.5 0 0.5 1]);
        end
        times2play = timePoints ;
        for kk = 1:length(times2play)
            ttemp = times2play(kk) ;
            im1 = imread(sprintf(fileBase, ttemp));
            imshow(im1) ;
            hold on;

            for iprev = 1:ii-1
                trackprev = currentTracks{iprev} ;
                plot(trackprev(kk, 1), trackprev(kk, 2), 'o', 'color', ...
                    bluecolor, 'markerSize', markerSize, 'lineWidth', lwidth)
            end

            tracknow = currentTracks{ii} ;
            plot(tracknow(kk, 1), tracknow(kk, 2), 'o', 'color', orange, ...
                'markerSize', markerSize, 'lineWidth', lwidth)

            hold off ;
            
            title(['playing tracks 1-' num2str(ii) ': t=' num2str(ttemp)])

            pause(0.01 * pausetime)
        end
    end
    
    %% Review track
    
    dmyii = dmyii + 1 ;
end

% Optionally convert to digraph
if nargout > 1
    trackGraph = trackingCell2Graph(tracks) ;
end

