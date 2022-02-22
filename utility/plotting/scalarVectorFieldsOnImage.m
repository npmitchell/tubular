function [h1, h2, h3] = scalarVectorFieldsOnImage(im, xxs, yys, sf, ...
    xxv, yyv, vx, vy, options)
%SCALARVECTORFIELDSONIMAGE(im, xx, yy, vx, vy, options)
%   Plot both a scalar field (sf) and vector field (vx,vy) on an image im. 
%   Coloring of the overlaid field is given by the vector field phase 
%   (or nematic order) and intensity is set by scalar field. Optionally,
%   overlay a quiver plot on the colored image with a desired subsampling
%   of the velocity field. Fields may be defined on vertices or on faces of
%   a mesh implied by xxs, yys, together with options.faces. 
%
% im : QxRx3 RGB float or int array
%   the image onto which we overlay the scalar and vector fields.
%   If not RGB, the image will adopt the current colormap.
% xxs : N x 1 float array
%   x values of scalar field evaluation points
% yys : M x 1 float array
%   y values of scalar field evaluation points
% sf : N x M float array
%   scalar field defined at (xxs(i), yys(j)) for all i=1:N, j=1:M
% xxv : P x 1 float array
%   x values of vector field evaluation points
% yyv : Q x 1 float array
%   y values of vector field evaluation points
% vx : P*Q x 1 float array
%   velocity in x direction at (xxv, yyv)
% vy : P*Q x 1 float array
%   velocity in y direction at (xxv, yyv)
% options : struct with fields 
%   faces : #faces x 3 (optional, if sf is defined on faces)
%       connectivity list of vertices xxs,yys if sf is defined on faces
%   style : str ('diverging' (default) 'phase' 'nematic')
%       style of overlaid scalar field, default is diverging
%       If 'phase' or 'nematic', use 'options.angle' to pass vector angles.
%   angle : #faces x 3 (optional, if style == 'phase')
%       scalar field for phasemap. Note that sf signals opacity if style ==
%       'phase' and angle will replace sf for color if supplied.
%       Otherwise, sf is used for both.
%   sscale : float (optional, used only if style == 'phase', default=max(abs(sf(:))))
%       scalar field maximum for opacity/alpha. If zero, color limits are
%       set to min(sf) and max(sf)
%   alpha : float (optional, used if style ~= 'phase', default=0.8)
%       opacity of overlaid scalar field
%   outfn : str
%       path to save image if given
%   label : str
%       colorbar label. Default is '$|v|$ [$\mu$m / min]' 
%   qsubsample : int (default=10)
%       subsampling factor of the quiver field
%   overlay_quiver : bool (default=true)
%       whether to show the quiverplot overlay
%   qscale : float
%       overall multiplication of the quiver magnitude for overlay
%       (then plotted at true value, "scale=0" in quiver())
%   outfn : str
%       output filename for figure as png 
%   figWidth : int (optional, default = 16) 
%       figure width in cm
%   figHeight : int (optional, default = 10) 
%       figure height in cm
%   cticks : numeric 1d array (optional)
%       colorbar tick values, if specified. Otherwise default.
%
% Returns
% -------
% h1 : handle for imshow
% h2 : handle for imagesc
% h3 : handle for quiverplot
%
% Example usage
% ------------- 
% % Create example image
% im = peaks ;
% [xxs, yys] = meshgrid(1:size(im, 1), 1:size(im, 2)) ;
% % Create example velocity field on the same field (here xxv=xxs, yyv=yys)
% vx = im ;
% vy = reshape(im', size(im)) ;
% sf = sqrt(vx.^2 + vy.^2) ;        % scalar field is velocity magnitude
% options.angle = atan2(vy, vx) ;  % polar field
% % Supply 1d x and y lists
% xyf.x = xxs(1, :) ;
% xyf.y = yys(:, 1)' ;
% options.visibility = 'on' ;
% options.overlay_quiver = false ;
%  scalarVectorFieldsOnImage(ones(size(im)), xs, ys, sf, ...
%     xs, ys, vx(:), vy(:), options)
%
% See also
% --------
% VECTORFIELDHEATPHASEONIMAGE, 
% 
% NPMitchell 2020

%% Default options
labelstr = '' ;
overlay_quiver = true ;
qsubsample = 10 ;
qscale = 10 ;
sscale = max(abs(sf(:))) ;
alphaVal = 0.8 ;
style = 'phase' ;  % 'diverging' (scalar), 'nematic', or 'phase' (polar)
% figure parameters
figWidth = 16 ; % cm
figHeight = 10 ; % cm
visibility = 'off' ;

%% Unpack options
if isfield(options, 'style') 
    style = options.style ;
end
if isfield(options, 'sscale') 
    sscale = options.sscale ;
end
if isfield(options, 'label')
    labelstr = options.label ;
end
if isfield(options, 'overlay_quiver')
    overlay_quiver = options.overlay_quiver ;
end
if isfield(options, 'qsubsample')
    qsubsample = options.qsubsample ;
end
if isfield(options, 'qscale') 
    qscale = options.qscale ;
end
if isfield(options, 'alpha') 
    alphaVal = options.alpha ;
end
if isfield(options, 'figWidth') 
    figWidth = options.figWidth ;
end
if isfield(options, 'figHeight') 
    figHeight = options.figHeight ;
end
if isfield(options, 'visibility')
    if strcmpi(options.visibility, 'off')
        visibility = 'off' ;
    elseif strcmpi(options.visibility, 'on')
        visibility = 'on' ;
    else
        error("options.visibility must be set to 'on' or 'off'")
    end
end

% 
% vangle = reshape(mod(atan2(vy, -vx), 2* pi), gridsz) ;
% speed = reshape(vecnorm([v2dsm_ii(:, 1), v2dsm_ii(:, 2)], 2, 2), gridsz);

%% Set up the figure
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', visibility) ;

% If grayscale image is passed, convert to RGB
if length(size(im)) < 3
    im = cat(3, im, im, im) ;
end
h1 = imshow(im) ;
hold on;

%% Determine if the given scalar field is defined on vertices or faces
sf_on_faces = false ;
if isfield(options, 'faces')
    if length(sf) == size(options.faces, 1) 
        sf_on_faces = true ;
    else
        error('faces passed do not match scalar field size')
    end
end

if sf_on_faces
    % Plot sf on faces
    disp('scalarVectorFieldsOnImage: plotting sf on faces')
    h2 = patch( 'Faces', options.faces, 'Vertices', [xxs(:), yys(:)],...
            'FaceVertexCData', sf, ...
            'FaceColor', 'flat', 'EdgeColor', 'none') ;    
    
    % If we use a phase style, then we obtain alpha info from sf
    % Otherwise, we apply a uniform alpha given by alphaVal
    if strcmpi(style, 'phase') || strcmpi(style, 'nematic')
        set(h2, 'AlphaData', sf / sscale)
    elseif strcmpi(style, 'diverging')
        alpha(alphaVal) ;
        if sscale > 0
            caxis([-sscale, sscale])
        else
            caxis([-max(abs(sf(:))), max(abs(sf(:)))])
        end
    end
else
    % The scalar field is defined on VERTICES (same as #pts given)
    % check that size of sf is compatible with xxs and yys
    same_size = length(sf(:)) == length(xxs(:)) && length(sf(:)) == length(yys(:)) ;
    commensurate = length(sf(:)) == (length(xxs(:)) * length(yys(:))) ;
    boundxy = length(xxs(:)) == 2 && length(yys(:)) == 2 ;
    if ~same_size && ~commensurate && ~boundxy
        error('Scalar field sf must be commensurate size with xxs, yys or with faces')
    end
    
    % Plot on vertices
    disp('scalarVectorFieldsOnImage: plotting sf on vertices')
    if isfield(options, 'angle')
        sfangle = options.angle ;
    elseif isfield(options, 'angles')
        sfangle = options.angles ;
    else
        sfangle = sf ;
    end
    h2 = imagesc(xxs, yys, sfangle) ;
    if strcmpi(style, 'phase') || strcmpi(style, 'nematic')
        set(h2, 'AlphaData', sf / sscale)
    elseif strcmpi(style, 'diverging')
        set(h2, 'alphaData', 0.5)
        if sscale > 0
            caxis([-sscale, sscale])
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADD QUIVER PLOT TO AXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if overlay_quiver   
    if qsubsample > 1
        disp('scalarVectorFieldsOnImage: subsampling quiver')
        wwv = length(xxv) ;
        hhv = length(yyv) ;
        vx = reshape(vx, [wwv, hhv]) ;
        vy = reshape(vy, [wwv, hhv]) ;
        QX = imresize(vx, [wwv / qsubsample, hhv / qsubsample], 'bicubic') ;
        QY = imresize(vy, [wwv / qsubsample, hhv / qsubsample], 'bicubic') ;
        xq = 1:qsubsample:wwv ;
        yq = 1:qsubsample:hhv ;
        [xg, yg] = meshgrid(xxv(xq), yyv(yq)) ;

        h3 = quiver(xg(:), yg(:), qscale * QX(:), qscale * QY(:), 0, 'k', 'LineWidth', 1.2) ;
    else
        size(xxv)
        size(vx)
        size(vy)
        h3 = quiver(xxv(:), yyv(:), qscale * vx(:), qscale * vy(:), 0, 'k', 'LineWidth', 1.2) ;
    end
else
    h3 = [] ;
end

%% Add title (optional)
if isfield(options, 'title')
    title(options.title, 'Interpreter', 'latex')
end

%% Add the colorbar in the style set in options struct
if strcmp(style, 'phase')
    %%%%%%%%%%%%%%%%%%%
    % Phasemap (0, 2*pi)
    %%%%%%%%%%%%%%%%%%%
    colormap phasemap
    caxis([0, 2*pi])
    if isfield(options, 'ylim')
        ylim(options.ylim)
        % ylim([size(im, 2) * 0.25, size(im, 2) * 0.75])
    end
    set(gca, 'Position', [0 0.11 0.85 0.8]) ;
    % Add phasebar
    phasebar('location', [0.87, 0.7, 0.1, 0.1]) ;
    % Add colorbar
    cax = axes('Position',[.9 .3 .02 .3]) ;
    [~, yyq] = meshgrid(0:4, 0:100) ;
    imshow(fliplr(yyq/max(yyq(:))))
    axis on
    yyaxis right
    ylabel(labelstr, 'color', 'k', ...
        'Interpreter', 'Latex')
    yticks([0 1])
    yticklabels({'0', num2str(sscale)})
    xticks([])
    yyaxis left
    yticks([])
    cax.YAxis(1).Color = 'k';
    cax.YAxis(2).Color = 'k';
elseif strcmp(style, 'nematic')
    %%%%%%%%%%%%%%%%%%%
    % Phasemap (0, pi)
    %%%%%%%%%%%%%%%%%%%
    colormap phasemap
    caxis([0, pi])
    if isfield(options, 'ylim')
        ylim(options.ylim)
        % ylim([size(im, 2) * 0.25, size(im, 2) * 0.75])
    end
    set(gca, 'Position', [0 0.11 0.85 0.8]) ;
    % Add phasebar
    phasebar('colormap', phasemap, ...
        'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic') ;
    % Add colorbar
    cax = axes('Position',[.9 .3 .02 .3]) ;
    [~, yyq] = meshgrid(0:4, 0:100) ;
    imshow(fliplr(yyq/max(yyq(:))))
    axis on
    yyaxis right
    ylabel(labelstr, 'color', 'k', ...
        'Interpreter', 'Latex')
    yticks([0 1])
    yticklabels({'0', num2str(sscale)})
    xticks([])
    yyaxis left
    yticks([])
    cax.YAxis(1).Color = 'k';
    cax.YAxis(2).Color = 'k';
elseif strcmp(style, 'diverging')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diverging colormap (-sscale, sscale)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colormap bwr
    if isfield(options, 'ylim')
        ylim(options.ylim)
    end
    set(gca, 'Position', [0 0.11 0.85 0.8]) ;
    % Add colorbar
    c = colorbar('Position',[.9 .333 .02 .333]) ;
    % ylabel(cax, labelstr, 'color', 'k', ...
    %     'Interpreter', 'Latex')
    
    % Make colorbar share the alpha of the image
    % Manually flush the event queue and force MATLAB to render the colorbar
    % necessary on some versions
    drawnow
    % Get the color data of the object that correponds to the colorbar
    cdata = c.Face.Texture.CData;
    % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
    cdata(end,:) = uint8(alphaVal * cdata(end,:));
    % Ensure that the display respects the alpha channel
    c.Face.Texture.ColorType = 'truecoloralpha';
    % Update the color data with the new transparency information
    c.Face.Texture.CData = cdata;
    c.Label.Interpreter = 'latex' ;
    c.Label.String = labelstr ;
    if isfield(options, 'cticks')
        c.Ticks = options.cticks ;
    end
else
    error('have not coded for this style yet')
end

%% Save the image to disk if outfn is supplied
if isfield(options, 'outfn')
    disp(['scalarVectorFieldsOnImage: saving ' options.outfn])
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
    saveas(fig, options.outfn) ;   
    close all
end
