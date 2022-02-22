function [h1, cb] = scalarFieldOnSurface(faces, vertices, sf, options)
%SCALARVECTORFIELDSONIMAGE(faces, vertices, sf, options)
%   Plot a scalar field (sf) as heatmap on a 3d mesh surface
%
% faces : M x 3 float array
%   mesh vertex connectivity list (triangulated faces)
% vertices : N x 3 float array
%   mesh vertices in 3D
% sf : Nx1 or Mx1 float array
%   scalar field to plot on surface
% qopts : struct with fields
%   style : str ('diverging' 'phase', default='diverging')
%       style of overlaid scalar field
%   sscale : float
%       scalar field maximum for clims/opacity
%   alpha : float (optional, used if style ~= 'phase')
%       opacity of overlaid scalar field
%   outfn : str
%       path to save image if given
%   label : str
%       colorbar label. Default is '$|v|$ [$\mu$m / min]' 
%   outfn : str
%       output filename for figure as png 
%   figWidth : int (optional, default = 16) 
%       figure width in cm
%   figHeight : int (optional, default = 10) 
%       figure height in cm
%   cbarPosition : 4 x 1 float array
%       position of the colorbar
%
% Returns
% -------
% h1 : handle for trisurf object
% cb : handle for colorbar object
%
% NPMitchell 2020

%% Default options
labelstr = '' ;
interpreter = 'latex' ;
sscale = max(abs(sf(:))) ;
alphaVal = 0.8 ;
style = 'diverging' ;
xlabelstr = '' ;
ylabelstr = '' ;
zlabelstr = '' ;
titlestr = '' ;
% figure parameters
figWidth = 16 ;  % cm
figHeight = 10 ; % cm
axposition = [0 0.11 0.85 0.8] ;
axisOff = true ;
if strcmp(style, 'diverging')
    cbar_position = [0.8327, 0.1095, 0.0374, 0.8167] ; 
    % formerly:  [.9 .333 .02 .333] ;
elseif strcmp(style, 'phasemap')
    cbar_position = [.9 .3 .02 .3] ;
end

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
if isfield(options, 'title')
    titlestr = options.title ;
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
if isfield(options, 'axPosition') 
    axposition = options.axPosition ;
end
if isfield(options, 'axisOff') 
    axisOff = options.axisOff ;
end
if isfield(options, 'cbarPosition')
    cbar_position = options.cbarPosition ;
end

% Set up the figure
if isfield(options, 'fig')
    fig = options.fig ;
end
if isfield(options, 'ax')
    set(fig, 'CurrentAxes', options.ax) ;
elseif isfield(options, 'axis')
    set(fig, 'CurrentAxes', options.axis) ;
elseif isfield(options, 'fig')
    % pass here.
    % we will next do: ax = gca ;
else
    close all
    fig = figure('units', 'normalized', ...
        'outerposition', [0 0 1 1], 'visible', 'off') ;
end
ax = gca ;

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
elseif strcmp(style, 'positive')
    cmap = parula ;
    colors = mapValueToColor(sf, [0, sscale], cmap) ;
end

% % check it
% close all
% set(gcf, 'visible', 'on')
% histogram(sf)
% waitfor(gcf)
% plot(colors')
% waitfor(gcf)
% error('here')

if size(colors, 1) == size(vertices, 1)
    h1 = patch( 'Faces', faces, 'Vertices', vertices, ...
        'FaceVertexCData', colors, ...
        'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'SpecularStrength', 0.1, 'DiffuseStrength', 0.1, ...
        'AmbientStrength', 0.8 );
elseif size(colors, 1) == size(faces, 1)
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

%% View and axis limits
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

% Labels and plot title
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
    cax = axes('Position', cbar_position) ;
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
elseif strcmp(style, 'diverging') || strcmp(style, 'positive')
    if strcmp(style, 'diverging')
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
    elseif strcmp(style, 'positive')
        disp('setting scalar field to positive style')
        colormap parula
        if isfield(options, 'ylim')
            ylim(options.ylim)
        end
        set(gca, 'Position', axposition) ;
        if axisOff
            axis off
        end

        % Set color axis limits
        if sscale > 0
            caxis([0, sscale])
        end
    end
    
    % Add colorbar
    cb = colorbar('Position', cbar_position) ;
    
    % ylabel(cax, labelstr, 'color', 'k', ...
    %     'Interpreter', interpreter)
    
    % Make colorbar share the alpha of the image
    % Manually flush the event queue and force MATLAB to render the colorbar
    % necessary on some versions
    drawnow
    % Get the color data of the object that correponds to the colorbar
    cdata = cb.Face.Texture.CData;
    % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
    cdata(end,:) = uint8(alphaVal * cdata(end,:));
    % Ensure that the display respects the alpha channel
    cb.Face.Texture.ColorType = 'truecoloralpha';
    % Update the color data with the new transparency information
    cb.Face.Texture.CData = cdata;
    cb.Label.Interpreter = interpreter ;
    cb.Label.String = labelstr ;
else
    error('have not coded for this style yet')
end

%% Save the image if outfn is supplied
if isfield(options, 'outfn')
    if ischar(options.outfn)
        disp(['scalarVectorFieldsOnImage: saving ' options.outfn])
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
        saveas(fig, options.outfn) ;   
        close all
    else
        % save one view for each filename
        for kk = 1:size(options.outfn, 2)
            outfn = options.outfn{kk} ;
            if kk == 1
                view(0, 0)
            elseif kk == 2
                view(0, 90)
            elseif kk == 3
                view(0, 180)
            elseif kk == 4
                view(0, 270)
            elseif kk == 5
                view(90, 0)
            elseif kk == 6
                view(180, 0)
            else
                error('have not defined this view yet')
            end
            disp(['scalarVectorFieldsOnImage: saving ' outfn])
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperPosition', [0 0 figWidth figHeight]);  
            saveas(fig, outfn) ;  
        end
        close all
        
    end
end
