function [k] = flipThroughStackFindLayer(allstack, title_preamble, axis, bigstep, fig, normalize)
%FLIPTHROUGHSTACKFINDLAYER Find a layer flipping through stack
%   Go through 3d data and change the last index of the stack to find a
%   desired layer. When found, press Enter. Returns the layer index 
%   visible when Enter/Return was pressed.
%
% Parameters
% ----------
% allstack : N x M x P images
%   The stack of images to flip through
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

if nargin < 3
    axis = 3 ;
    bigstep = 10 ;
    fig = gcf ;
    normalize = false ;
elseif nargin < 4
    bigstep = 10 ;
    fig = gcf ;
    normalize = false ;
elseif nargin < 5
    fig = gcf ;
    normalize = false ;
end

% Flip through the images to find the right stack layer
ax = axes('Parent', fig);
max_k = size(allstack, axis) ;

% Start with first frame as initial setting
k = 1 ;

% get intensity limits
if normalize
    if axis == 3
        trace = (squeeze(allstack(:, :, k))) ;
    elseif axis == 2
        trace = (squeeze(allstack(:, k, :))) ;
    elseif axis == 1
        trace = (squeeze(allstack(k, :, :))) ;
    end
    [f,x] = ecdf(trace(:));
    f1 = find(f>0, 1, 'first');
    f2 = find(f<0.999, 1, 'last');
    imin = double(x(f1)) ;
    imax = double(x(f2)) ;
else
    imin = double(min(allstack(:))) ;
    imax = double(max(allstack(:))) ;
end

pressed_enter = false ;
while k <= max_k && ~pressed_enter
    if axis == 3
        trace = (squeeze(allstack(:, :, k))) ;
    elseif axis == 2
        trace = (squeeze(allstack(:, k, :))) ;
    elseif axis == 1
        trace = (squeeze(allstack(k, :, :))) ;
    end
    trace = mat2gray(trace, [imin imax]) ;
    imshow(trace)
    hold(ax, 'on');
    title(ax, [title_preamble ': # ',num2str(k),'']);
    hold(ax, 'off');
    if exist('Xlim', 'var')
        set(ax, 'xlim', Xlim) 
        set(ax, 'ylim', Ylim) 
    end
    was_a_key = waitforbuttonpress;
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
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'i')
        imax = max(0, imax * 0.9) ;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'k')
        imax = max(0, imax * 1.1) ;
    end
    
    Xlim = xlim ;
    Ylim = ylim ;
    
    % Regulate the bounds of k (no negative layers)
    if k < 1
        k = 1 ;
    elseif k > max_k
        k = max_k ;
    end
end
return 

