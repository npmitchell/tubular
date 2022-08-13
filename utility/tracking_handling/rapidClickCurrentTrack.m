function [ax0, minX, maxX, minY, maxY] = rapidClickCurrentTrack(ii, ...
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
tp = timePoints(tidx) ;

% Overlay previous and current image and crop each
im = imread(sprintf(fileBase, tp));

if ~mode
    if ~isempty(Xlim) && ~isempty(Ylim)
        minX = max(1, round(Xlim(1))) ;
        maxX = min(size(im, 2), round(Xlim(2))) ;
        minY = max(1, round(Ylim(1))) ;
        maxY = min(size(im, 1), round(Ylim(2))) ;
        imCrop = im(minY:maxY, minX:maxX) ;
    else
        imCrop = im ;
        minX = 0 ;
        minY = 0 ;
        maxX = size(im, 2) ;
        maxY = size(im, 1) ;
    end
end
exten = ' [acquire on THIS image]';


%-----------------------------
% Only axis
%-----------------------------
cla
ax0 = gca ;
imshow(imCrop);
hold on;

% Build nearby id's to show in CURRENT time
% Past tracks
nearbyDone  = nan(length(otherIDs), 2) ;
kk = 1 ;
for otherID = otherIDs
    otherTrack = currentTracks{otherID} ;
    if ~isempty(otherTrack)
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
        xx = otherTrack(tidx, 1) ;
        yy = otherTrack(tidx, 2) ;
        if xx > minX && xx < maxX && yy > minY && yy < maxY
            nearbyNext(kk, :) = [xx-minX, yy-minY] ;
            kk = kk + 1 ;
        end
    end
    nearbyNext = nearbyNext(1:kk-1, :) ;
    % OTHERS -- CURRENT: Plot other tracks' CURRENT positions on overlay
    plot(nearbyNext(:, 1), nearbyNext(:, 2), 's', ...
        'color', green, 'markerSize', markerSize, 'lineWidth', 1)
end

% Plot current position if it exists for this track
if ~isnan(trackii(tidx, 1))
    plot(trackii(tidx, 1), trackii(tidx, 2), 'o', 'color', ...
        orange, 'markerSize', markerSize, 'lineWidth', 1)
end
hold off;
title(['Track ' num2str(ii) ': t=' num2str(tp) exten])
