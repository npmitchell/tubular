function [fig, ax1, ax2, imhandle, shandle] = heatmap_on_image(im, ...
    xfield, yfield, cfield, options)
% HEATMAP_ON_IMAGE(im, xfield, yfield, cfield, options)
%   Plot a scalar field as heatmap on image
%   This function has some unresolved bugs 2019-10-06: please fix
%
% Parameters
% ----------
% im : 2d image 
% xfield : x values for heatmap
% yfield : y values for heatmap
% cfield : scalar field values evaluated at xfield, yfield
% options : struct optionally containing the following
%   colorstyle -- 'pcolor', 'imagesc', or 'scatter'
%   flipy      -- whether to flip the heatmap field in Y
%   alpha      -- opacity value (0, 1) 
%   clims      -- color limits [min max]
%   cmap       -- colormap
%
% Returns
% -------
% fig : figure instance
% ax1 : the image axis
% ax2 : the heatmap axis
% 
% NPMitchell 2019

if isfield(options, 'colorstyle')
    colorstyle = options.colorstyle ;
else
    colorstyle = 'imagesc' ;
end
if isfield(options, 'flipy')
    flipy = options.flipy ;
else
    flipy = true ;
end
if isfield(options, 'alpha')
    alphaVal = options.alpha ;
else
    alphaVal = 0.5 ;
end
if isfield(options, 'clims')
    clims = options.clims ;
else
    clims = [-1 1];
end
if isfield(options, 'cmap')
    cmap = options.cmap ;
else
    cmap = [] ;
end

% Background image
fig = figure ;
ax1 = axes;
imhandle = imshow(im); 
colormap(ax1,'gray');
set(ax1,'ydir','normal');

% Foreground image
ax2=axes;
tmp = gcf ;
set(tmp.Number, 'currentaxes', ax2)
if strcmp(colorstyle, 'pcolor')
    % Use pcolor. Check that this is correct
    shandle = pcolor(xx, yy, aniso) ;   
    shandle.FaceColor = 'interp' ;
    shandle.EdgeColor = 'none' ;
    set(shandle, 'facealpha', alphaVal)
elseif strcmp(colorstyle, 'imagesc')
    % Instead use imagesc
    shandle = imagesc(xfield, yfield, cfield) ;
    alpha(ax2, alphaVal)
else
    shandle = scatter(xfield, yfield, cfield, 'filled', 'edgecolor', 'none') ;
    alpha(ax2, alphaVal)
end
    
caxis(clims) ;
set(ax2,'color','none','visible','off');

if ~isempty(cmap)
    colormap(ax2, cmap);
end
if flipy
    set(ax2,'ydir','normal');
end

% Link the axes
linkaxes([ax1 ax2])
axis equal

% Make both axes invisible
axis off
set(ax1, 'Visible', 'off')
% Get the current axis size in case things get disrupted
originalPos_ax1 = get(ax1, 'Position') ;
originalPos_ax2 = get(ax2, 'Position') ;
title(['t=' num2str(t) ' min'],...
    'FontWeight','Normal')
set(ax2, 'Position', originalPos_ax2);

end

