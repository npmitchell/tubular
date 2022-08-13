function [k] = flipThroughStackFindLayerOverlay(allstack, overlay, ...
    title_preamble, axis, bigstep, fig, k0, xlims, ylims)
%FLIPTHROUGHSTACKFINDLAYER Find a layer flipping through stack with a
%   binary overlay atop the image, as done for segmentation inspection.
%   Go through 3d data and change the last index of the stack to find a
%   desired layer. When found, press Enter. Returns the layer index 
%   visible when Enter/Return was pressed. 
%
%   Small steps: <--->
%   Big steps: ^
%              |
%              v
%
% Parameters
% ----------
% allstack : N x M x P images
%   The stack of images to flip through
% imoverlay : N x M binary image
%   The overlay (for ex, binary skeleton or mask)
% title_preamble : str
%   Title preceding the layer index 
% axis : int
%   Which axis to flip through
% 
% Returns
% -------
% k : int
%   The layer index visible when Enter/Return is pressed.
% 
% NPMitchell 2019

if nargin < 4
    axis = 3 ;
end
if nargin < 5
    bigstep = 10 ;
    fig = gcf ;
end
if nargin < 6
    fig = gcf ;
end

if nargin < 7 
    % Start with first frame as initial setting
    k = round(0.5 * size(allstack, axis)) ;
else
    if isempty(k0)
        % Start with first frame as initial setting
        k = round(0.5 * size(allstack, axis)) ;
    else
        k = k0 ;
    end
end
if nargin < 8
    xlims = [] ;
end
if nargin < 9
    ylims = [] ;
end

% Flip through the images to find the right stack layer
ax = axes('Parent', fig);
max_k = size(allstack, axis) ;


% get intensity limits
if axis == 3
    slice = (squeeze(allstack(:, :, k))) ;
elseif axis == 2
    slice = (squeeze(allstack(:, k, :))) ;
elseif axis == 1
    slice = (squeeze(allstack(k, :, :))) ;
end
[f,x] = ecdf(slice(:));
f1 = find(f>0, 1, 'first');
f2 = find(f<0.999, 1, 'last');
imin = double(x(f1)) ;
imax = double(x(f2)) ;

pressed_enter = false ;
while ~pressed_enter
    if axis == 3
        slice = (squeeze(allstack(:, :, k))) ;
    elseif axis == 2
        slice = (squeeze(allstack(:, k, :))) ;
    elseif axis == 1
        slice = (squeeze(allstack(k, :, :))) ;
    end
    slice = mat2gray(slice, [imin imax]) ;
    bim = imoverlay(slice, overlay, 'blue');
    imshow(bim)
    if ~isempty(xlims)
        xlim(xlims)
    end
    if ~isempty(ylims)
        ylim(ylims)
    end
    % title('')
    % hold(ax, 'on');
    title([title_preamble ': # ',num2str(k),'']);
    % hold(ax, 'off');
    was_a_key = waitforbuttonpress ;
    disp('key pressed')
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
        pressed_enter = true ;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
        k = k + 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
        k = k - 1; 
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'uparrow')
        k = k + bigstep;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'downarrow')
        k = k - bigstep;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'z')
        ax = gca ;
        funits = get(fig,'Units');
        aunits = get(ax,'Units');  % Store current units for axes and figure in order to restore them at the end of the routine.
        tmp = waitforbuttonpress ;
        BOUND = rbbox;       % Get the boundary of the box
        set(fig,'Units','pixels'); % Set the figure units to 'pixels' 
        set(ax,'Units','pixels');  % Set the axis units to 'pixels' 
        P = get(ax, 'Position');  % Retrieve the axes position.
        XLimit = get(ax, 'Xlim');  YLimit = get(ax, 'Ylim');
        % Get the coordinates of the lower left corner of the box
        % Its height and width, and the X-Y coordinates of the 4 corners
        DeltaX = XLimit(2)-XLimit(1);  DeltaY = YLimit(2)-YLimit(1);
        LeftDist = BOUND(1)-P(1);  UpDist = BOUND(2)-P(2);
        % Defining some useful quantities which will be used often
        XLow = XLimit(1)+DeltaX*LeftDist/P(3);
        XHigh = XLimit(1)+DeltaX*(LeftDist+BOUND(3))/P(3);
        YLow = YLimit(1)+DeltaY*UpDist/P(4);
        YHigh = YLimit(1)+DeltaY*(UpDist+BOUND(4))/P(4);
        % This is the X-Y information about the corners of the RBBOX
        xlims = [XLow, XHigh]; ylims = [YLow, YHigh];
        % ZOOM IN: Set the axes limits according to the zoom box
        set(ax,'XLim',xlims,'YLim',ylims);
        set(fig,'Units',funits);
        set(ax,'Units',aunits);  % Reset the figure and axis units to their original settings
        
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'i')
        imax = max(0, imax * 0.9) ;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'k')
        imax = max(0, imax * 1.1) ;
    end
    
    % Regulate the bounds of k (no negative layers)
    if k < 1
        k = 1 ;
    elseif k > max_k
        k = max_k ;
    end
    
    xlims = xlim() ;
    ylims = ylim() ;
end
return 

