function [fig, h1, h2, h3] = ...
    scalarVectorFieldsOnSurface(faces, vertices, sf, ...
                                xxv, yyv, zzv, vx, vy, vz, options)
%SCALARVECTORFIELDSONIMAGE(im, xx, yy, vx, vy, options)
%   Plot both a scalar field (sf) as heatmap and vector field (vx,vy) as 
%   quiver arrows.
%
% Parameters
% ----------
% faces : F x 3 int array
%   indices into vertices of F mesh faces
% vertices : V x 3 float array
%   coordinates of V mesh vertices
% sf : Fx1 or Vx1 float array 
%   scalar field to plot as heatmap, defined on faces or vertices
% xxv : Q x 1 float array
%   x coordinates of velocity evaluation points
% yyv : Q x 1 float array
%   y coordinates of velocity evaluation points
% zzv : Q x 1 float array
%   z coordinates of velocity evaluation points
% vx : Q x 1 float array
%   velocity in x direction
% vy : Q x 1 float array
%   velocity in y direction
% vz : Q x 1 float array
%   velocity in z direction
% qopts : struct with fields
%   style : str ('diverging' 'phase')
%       style of overlaid scalar field, default is diverging
%   sscale : float (optional, used only if style == 'phase')
%       scalar field maximum for opacity/alpha
%   alpha : float (optional, used if style ~= 'phase')
%       opacity of overlaid scalar field
%   outfn : str
%       path to save image if given
%   label : str
%       colorbar label. Default is '$|v|$ [$\mu$m / min]' 
%   qsubsample : int
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
% 
%
% NPMitchell 2020

%% Default options
labelstr = '' ;
interpreter = 'latex' ;
qscale = 10 ;
sscale = 0 ;
alphaVal = 0.8 ;
style = 'diverging' ;
qsubsample = 1 ;
xlabelstr = '' ;
ylabelstr = '' ;
zlabelstr = '' ;
titlestr = '' ;
% figure parameters
figWidth = 16 ;  % cm
figHeight = 10 ; % cm
lw = 1.2 ;       % linewidth
axposition = [0 0.11 0.85 0.8] ;
axisOff = true ;

% Unpack options
if isfield(options, 'style') 
    style = options.style ;
end
if isfield(options, 'sscale') 
    sscale = options.sscale ;
end
if isfield(options, 'label')
    labelstr = options.label ;
end
if isfield(options, 'title')
    titlestr = options.title ;
end
if isfield(options, 'qscale') 
    qscale = options.qscale ;
end
if isfield(options, 'qsubsample')
    qsubsample = options.qsubsample ;
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
if isfield(options, 'linewidth') 
    lw = options.linewidth ;
end
if isfield(options, 'axPosition') 
    axposition = options.axPosition ;
end
if isfield(options, 'axisOff') 
    axisOff = options.axisOff ;
end

% Set up the figure
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;

% Add the scalar field defined on faces OR on vertices
if size(sf, 1) == size(vertices, 1)
    h1 = patch( 'Faces', faces, 'Vertices', vertices, ...
        'FaceVertexCData', sf, ...
        'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'SpecularStrength', 0.1, 'DiffuseStrength', 0.1, ...
        'AmbientStrength', 0.8 );
elseif size(sf, 1) == size(faces, 1)
    % Add the scalar field defined on faces
    if strcmp(style, 'phasemap')
        cmap = phasemap ;
        colors = mapValueToColor(sf, [0, 2*pi], cmap) ;
    elseif strcmp(style, 'diverging')
        cmap = bwr ;
        if sscale > 0
            colors = mapValueToColor(sf, [-sscale, sscale], cmap) ;
        else
            colors = mapValueToColor(sf, [min(sf(:)), max(sf(:))], cmap) ;
        end
    end
    
    % plot in 3d as patches
    h1 = patch( 'Faces', faces, 'Vertices', vertices, ...
        'FaceVertexCData', colors, 'FaceColor', 'flat', ...
        'EdgeColor', 'none', ...
        'SpecularStrength', 0.1, 'DiffuseStrength', 0.1, ...
        'AmbientStrength', 0.8 );
else
    error('Colors for patch object must be #faces x 3 or #verts x 3')
end
hold on;
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUIVER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if qscale > 0
    if qsubsample > 1
        h3 = quiver3(xxv(1:qsubsample:end), ...
            yyv(1:qsubsample:end), zzv(1:qsubsample:end), ...
            qscale * vx(1:qsubsample:end), ...
            qscale * vy(1:qsubsample:end), ...
            qscale * vz(1:qsubsample:end), 0, 'k',...
            'LineWidth', lw) ;
    else
        h3 = quiver3(xxv(:), yyv(:), zzv(:), ...
            qscale * vx(:), qscale * vy(:), qscale * vz(:), 0, 'k',...
            'LineWidth', lw) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View and axis limits
if isfield(options, 'view') 
    view(options.view(1), options.view(2)) ;
end
if isfield(options, 'xlim') 
    xlim(options.xlim) ;  
end
if isfield(options, 'ylim') 
    ylim(options.ylim) ;
end
if isfield(options, 'zlim')
    zlim(options.zlim) ;
end

%% Labels and plot title
if ~isempty(xlabelstr)
    xlabel(xlabel, 'Interpreter', interpreter)
end
if ~isempty(ylabelstr)
    ylabel(ylabel, 'Interpreter', interpreter)
end
if ~isempty(zlabelstr)
    zlabel(zlabel, 'Interpreter', interpreter)
end
% Add title (optional)
if ~isempty(titlestr)
    title(titlestr, 'Interpreter', interpreter) 
end
hold on;


%% Add the colorbar in the style set in options struct
if strcmp(style, 'phase')
    disp('setting scalar field to phase style')
    %%%%%%%%%%%%%%%%%%%
    % Phasemap
    %%%%%%%%%%%%%%%%%%%
    colormap phasemap
    if isfield(options, 'ylim')
        ylim(options.ylim)
    end
    caxis([0, 2*pi])
    set(gca, 'Position', axposition) ;
    % Add phasebar
    phasebar('location', [0.87, 0.7, 0.1, 0.1]) ;
    % Add colorbar
    cax = axes('Position',[.9 .3 .02 .3]) ;
    [~, yyq] = meshgrid(0:4, 0:100) ;
    imshow(fliplr(yyq/max(yyq(:))))
    axis on
    yyaxis right
    ylabel(labelstr, 'color', 'k', ...
        'Interpreter', interpreter)
    yticks([0 1])
    yticklabels({'0', num2str(sscale)})
    xticks([])
    yyaxis left
    yticks([])
    cax.YAxis(1).Color = 'k';
    cax.YAxis(2).Color = 'k';
elseif strcmp(style, 'diverging')
    disp('setting scalar field to diverging style')
    colormap bwr
    if isfield(options, 'ylim')
        ylim(options.ylim)
    end
    set(gca, 'Position', axposition) ;
    if axisOff
        axis off
    end
    
    % Set color axis limits
    if sscale > 0
        caxis([-sscale, sscale])
    end
    
    % Add colorbar
    c = colorbar('Position',[.9 .333 .02 .333]) ;
    % ylabel(cax, labelstr, 'color', 'k', ...
    %     'Interpreter', interpreter)
    
    % Make colorbar share the alpha of the image
    % Manually flush the event queue and force MATLAB to render the colorbar
    % necessary on some versions
    drawnow
    % Get the color data of the object that correponds to the colorbar
    cdata = c.Face.Texture.CData;
    % Change the 4th channel (alpha channel) to vary
    cdata(end,:) = uint8(alphaVal * cdata(end,:));
    % Ensure that the display respects the alpha channel
    c.Face.Texture.ColorType = 'truecoloralpha';
    % Update the color data with the new transparency information
    c.Face.Texture.CData = cdata;
    c.Label.Interpreter = interpreter ;
    c.Label.String = labelstr ;
    if isfield(options, 'cticks')
        c.Ticks = options.cticks ;
    end
else
    error(['have not coded for this style: ' style])
end


%% Save the image if outfn is supplied
if isfield(options, 'outfn')
    disp(['scalarVectorFieldsOnImage: saving ' options.outfn])
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
    saveas(fig, options.outfn) ;   
    close all
end
