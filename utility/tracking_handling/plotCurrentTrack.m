function [ax0, ax1, ax2, minX, maxX, minY, maxY] = plotCurrentTrack(ii, ...
    tidx, timePoints, fileBase, trackii, ...
    currentTracks, Xlim, Ylim)
%plotCurrentTrack(ii, tidx, timePoints, fileBase, trackii, ...
%         currentTracks, bluecolor, orange, lwidth, markerSize, Xlim, Ylim)
%
% Parameters
% ----------
% ii : track index
% tidx : timestamp index
% timePoints : times over which to index
% fileBase : 
%   sprintf(fileBase, timePoints(tidx)) is the current image in
%   which to track
% trackii : #timePoints x 2 float array
%   current track 
% currentTracks : cell array of #timePoints x 2 float arrays
%   all cell tracks acquired until now
% bluecolor
%
%
%
%
% Returns
% -------
% [ax0, ax1, ax2] : axis handles 
%   ax0: axis handle for axis 1 (prev timepoint) 
%   ax1: axis handle for merge RGB overlay 
%   ax2: axis handle for current timepoint
% minX, maxX, minY, maxY
%
%
%
%
% NPMitchell 2021

%-----------------------------
% Plotting options
%-----------------------------
bluecolor = [0 , 0.4470, 0.7410 ];
orange = [ 0.8500,  0.3250 , 0.0980 ];
green = [ 0.2000    0.9000    0.2000 ] ;
lwidth = 3 ;
markerSize = 20 ;
mode = 0 ;  % full (1) or crop (0)

% Process some input
nTracks = length(currentTracks) ;
otherIDs = setdiff(1:nTracks, ii) ;

%-----------------------------
% Image loading
%-----------------------------
tp0 = timePoints(max(1, tidx-1)) ;
tp = timePoints(tidx) ;

% Overlay previous and current image and crop each
im0 = imread(sprintf(fileBase, tp0));
im = imread(sprintf(fileBase, tp));

if ~mode
    if ~isempty(Xlim) && ~isempty(Ylim)
        minX = max(1, round(Xlim(1))-20) ;
        maxX = min(size(im, 2), round(Xlim(2))+20) ;
        minY = max(1, round(Ylim(1))-20) ;
        maxY = min(size(im, 1), round(Ylim(2))+20) ;
        imCrop = im(minY:maxY, minX:maxX) ;
        im0 = im0(minY:maxY, minX:maxX) ;
        exten1 = ' [do NOT click here]';
        exten2 = ' [acquire on THIS image]';
    else
        imCrop = im ;
        minX = 0 ;
        minY = 0 ;
        maxX = Inf ;
        maxY = Inf ;
        exten1 = ' [acquire on any image]';
        exten2 = ' [acquire on any image]';
    end
    rgb = cat(3, im0, imCrop, imCrop) ;
else
    rgb = cat(3, im0, im, im) ;
    exten1 = ' [acquire on any image]';
    exten2 = ' [acquire on any image]';
end


%-----------------------------
% First axis
%-----------------------------
ax0 = subtightplot(2, 2, 1) ;
imshow(im0) ;
hold on;
% Show lineage through previous timepoint -- SELF -- OLD
if tidx > 1
    plot(trackii(1:tidx-1, 1)-minX, trackii(1:tidx-1, 2)-minY, '-', 'color', orange, 'lineWidth', lwidth)
    plot(trackii(tidx-1, 1)-minX, trackii(tidx-1, 2)-minY, 'o', 'color', orange, 'markerSize', markerSize, 'lineWidth', lwidth)
else
    % This is the first timepoint so show the current lineage starting
    % position in orange -- SELF -- CURRENT
    plot(trackii(tidx, 1)-minX, trackii(tidx, 2)-minY, 'o', 'color', orange, 'markerSize', markerSize, 'lineWidth', lwidth)
end

% Build nearby id's to show
nearbyDone  = nan(length(otherIDs), 2) ;
kk = 1 ;
for otherID = otherIDs(otherIDs<ii)
    otherTrack = currentTracks{otherID} ;
    xx = otherTrack(max(1, tidx-1), 1) ;
    yy = otherTrack(max(1, tidx-1), 2) ;
    if xx > minX && xx < maxX && yy > minY && yy < maxY
        nearbyDone(kk, :) = [xx-minX, yy-minY] ;
        kk = kk + 1 ;
    end
end
nearbyDone = nearbyDone(1:kk-1, :) ;
% Show other tracks DONE in blue -- OTHERS -- CURRENT
plot(nearbyDone(:, 1), nearbyDone(:,2),...
    's', 'color', bluecolor, 'markerSize', markerSize, 'lineWidth', 1) 


% Show other tracks NOT DONE in green -- OTHERS -- CURRENT
if ii < nTracks
    kk = 1 ;
    % try
        nearbyNext  = nan(length(otherIDs), 2) ;
        for otherID = otherIDs(otherIDs>ii)
            otherTrack = currentTracks{otherID} ;
            if isempty(otherTrack)
                disp('plotCurrentTrack: Track is empty!')
            else
                xx = otherTrack(max(1, tidx-1), 1) ;
                yy = otherTrack(max(1, tidx-1), 2) ;
                if xx > minX && xx < maxX && yy > minY && yy < maxY
                    nearbyNext(kk, :) = [xx-minX, yy-minY] ;
                    kk = kk + 1 ;
                end
            end
        end
        nearbyNext = nearbyNext(1:kk-1, :) ;

        % Show other tracks NOT DONE in green -- OTHERS -- CURRENT
        plot(nearbyNext(:, 1), nearbyNext(:, 2),...
            's', 'color', green, 'markerSize', markerSize, 'lineWidth', 1) 
    % catch
    %     disp('Could not display track next')
    % end
end

hold off;
clearvars nearbyNext nearbyDone


hold off;

title(['t=' num2str(tp0) exten1 ])

%-----------------------------
% Second axis
%-----------------------------
ax1 = subtightplot(2, 2, 2) ;
imshow(rgb);
hold on;

% Build nearby id's to show in CURRENT time
% Past tracks
nearbyDone  = nan(length(otherIDs), 2) ;
kk = 1 ;
for otherID = otherIDs
    otherTrack = currentTracks{otherID} ;
    if isempty(otherTrack)
        disp('plotCurrentTrack: no other tracks to plot!')
    else
        xx = otherTrack(tidx, 1) ;  
        yy = otherTrack(tidx, 2) ;
        if xx > minX && xx < maxX && yy > minY && yy < maxY
            nearbyDone(kk, :) = [xx-minX, yy-minY] ;
            kk = kk + 1 ;
        end
    end
end
nearbyDone = nearbyDone(1:kk-1, :) ;

% OTHERS -- CURRENT: Plot other tracks' CURRENT positions on overlay
plot(nearbyDone(:, 1), nearbyDone(:, 2), 's', ...
    'color', bluecolor, 'markerSize', markerSize, 'lineWidth', 1)

% Future tracks
if ii < nTracks
    nearbyNext  = nan(length(otherIDs), 2) ;
    kk = 1 ;
    for otherID = otherIDs(otherIDs>ii)
        otherTrack = currentTracks{otherID} ;
        if isempty(otherTrack)
            disp('plotCurrentTrack: no other tracks to plot')
        else
            xx = otherTrack(tidx, 1) ;
            yy = otherTrack(tidx, 2) ;
            if xx > minX && xx < maxX && yy > minY && yy < maxY
                nearbyNext(kk, :) = [xx-minX, yy-minY] ;
                kk = kk + 1 ;
            end
        end
    end
    nearbyNext = nearbyNext(1:kk-1, :) ;
    % OTHERS -- CURRENT: Plot other tracks' CURRENT positions on overlay
    plot(nearbyNext(:, 1), nearbyNext(:, 2), 's', ...
        'color', green, 'markerSize', markerSize, 'lineWidth', 1)
end

% SELF -- OLD
if tidx > 1
    plot(trackii(tidx-1, 1)-minX, trackii(tidx-1, 2)-minY, 'o', ...
            'color', orange,  'markerSize', markerSize, 'lineWidth', lwidth)
end
hold off;
title(['old Self, current Others, ' exten1])

%-----------------------------
% Third axis
%-----------------------------
ax2 = subtightplot(2, 1, 2) ;
imshow(im) ;
hold on;
% Plot other tracks (not current track)
plot(nearbyDone(:, 1)+minX, nearbyDone(:, 2)+minY, 's', ...
    'color', bluecolor, 'markerSize', markerSize, 'lineWidth', 1)
% Future tracks
if ii < nTracks
    plot(nearbyNext(:, 1)+minX, nearbyNext(:, 2)+minY, 's', ...
        'color', green, 'markerSize', markerSize, 'lineWidth', 1)
end

% Plot current position if it exists for this track
if ~isnan(trackii(tidx, 1))
    plot(trackii(tidx, 1), trackii(tidx, 2), 'o', 'color', ...
        orange, 'markerSize', markerSize, 'lineWidth', 1)
end
hold off;
if isempty(Xlim)
    Xlim = get(gca, 'XLim');
    Ylim = get(gca, 'YLim');
end
if tidx > 1
    set(gca, 'XLim', Xlim , 'YLim', Ylim);
end
title(['Track ' num2str(ii) ': t=' num2str(tp) exten2])
