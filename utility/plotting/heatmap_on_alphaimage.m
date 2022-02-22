function h2 = heatmap_on_alphaimage(im, xfield, yfield, cfield, options)
% HEATMAP_ON_IMAGE Plot a scalar field as heatmap on image with opacity
% This function has some unresolved bugs 2019-10-06: please fix
% Note: xfield increases from top to bottom of image, yfield from L to R
% 
% 
% Parameters
% ----------
% im : 2d image 
% xfield : x values for heatmap, can be ndgrid, linspace, or scattered (not meshgrid)
% yfield : y values for heatmap, can be ndgrid, linspace, or scattered (not meshgrid)
% cfield : scalar field values evaluated at xfield, yfield
% options : struct optionally containing the following
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


if isfield(options, 'flipy')
    flipy = options.flipy ;
else
    flipy = true ;
end
if isfield(options, 'clims')
    clims = options.clims ;
else
    if isfield(options, 'caxis')
        clims = options.caxis;
    else
        clims = [-1 1];
    end
end
if isfield(options, 'cmap')
    cmap = options.cmap ;
else
    cmap = [] ;
end
if isfield(options, 'gridded')
    if options.gridded
        scattered = false ;
    else
        scattered = true ;
    end
    check_for_grid = false ;
else
    check_for_grid = true ;
end


% Get x and y linspace for the image
pp = 1:size(im, 1) ;  % #rows
qq = 1:size(im, 2) ;  % #columns
[xg, yg] = ndgrid(pp, qq) ;
% error('break')

% Check if the input data is gridded
if check_for_grid
    dx = gradient(xfield) ;
    dy = gradient(yfield) ;
    if all(dx == dx(1)) && all(dy == dy(1))
        scattered = false ;
    else 
        scattered = true ;
    end
end

% Background image is in imshow then used as alpha
if scattered
    gF = scatteredInterpolant(xfield, yfield, cfield);
    vF = gF(xg, yg) ;
else
    if ~all(size(xfield) == size(cfield))
        % Here we assume since the size is not right, we make a grid
        % This ensures that the image is not transposed. 
        % Note that meshgrid would be transposed.
        [xfield, yfield] = ndgrid(xfield, yfield) ;
    end
    gF = griddedInterpolant(xfield, yfield, cfield);
    vF = gF(xg, yg) ;
    
    if any(size(options.alpha) > 1)
        if all(size(options.alpha) == size(cfield))
            alphF = griddedInterpolant(xfield, yfield, options.alpha) ;
            opacityF = alphF(xg, yg) ;
        else
            size(options.alpha)
            size(cfield)
            error('Have not coded for this alpha shape yet')
        end
    else
        opacityF = 1  ;
    end
end

% PLOT IT
% h1 = imshow(im) ;
h2 = imagesc(qq, pp, vF) ;
if any(size(options.alpha) > 1)
    % Use the product of the image and the opacity field F as alpha values
    opacity = double(im) .* opacityF ;
    set(h2, 'AlphaData', uint8(opacity))
else
    set(h2, 'AlphaData', im)
end

caxis(clims)

if ~isempty(cmap)
    colormap(cmap);
end

% Flip the image in Y if flipy is true
if flipy
    set(gca,'ydir','normal');
end

end

