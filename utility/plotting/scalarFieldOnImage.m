function [h1, h2] = scalarFieldOnImage(im, xy_or_fxy, sf, alphaVal, ...
    scale, labelOptions, varargin)
%SCALARFIELDONIMAGE(im, xx, yy, field, alphaVal, scale, label)
% Plot a scalar field over an image, colored by magnitude, with const alpha
% The heatmap style may be diverging, phasemap, positive, or negative.
%
% Example Usage
% -------------
% labelOpts.label = '$|v|$ [$\mu$m/min]' ;
% scalarFieldOnImage(im, xx, yy, reshape(vecnorm(vsm_ii, 2, 2), gridsz),...
%     alphaVal, vtscale, labelOpts, 'Style', 'Positive')
%
% Parameters
% ----------
% im : PxQ numeric array
%   RGB or grayscale image
% xy_or_fxy : N x 2 float array, or struct with fields (faces/f, xy/v) 
%   xy coordinates of the field evaluation locations, or struct with faces
%   and vertex locations for drawing patches if field is defined on faces
%   The name of the fields are a bit flexible: faces can be faces or f,
%   vertices can be xy or v or pts or vertices
% sf : NxM float array
%   the scalar field to plot as heatmap
% alphaVal : float
%   the opacity of the heatmap
% scale : float
%   maximum absolute value for the field to be saturated in colormap
% labelOptions: struct with fields
%   label : str
%       colorbar label, interpreted through Latex by default
%   title : str
%       optional title, interpreted through Latex by default
%   xlabel : str
%       optional title, interpreted through Latex by default
%   ylabel : str
%       optional title, interpreted through Latex by default
% varargin : keyword arguments (optional, default='diverging') 
%   options for the plot, with names
%   'style' : 'phasemap' or 'diverging' or 'positive' or 'negative'
%   'interpreter' : 'Latex', 'default'/'none'
% 
%
% Returns
% -------
% h1 : handle for imshow
% h2 : handle for imagesc
%
% See also
% --------
% VECTORFIELDHEATPHASEONIMAGE, SCALARVECTORFIELDSONIMAGE
%
% NPMitchell 2020

%% Default label Options
xlabelstr = '' ;
ylabelstr = '' ;
titlestr = '' ;
label = '' ;

%% Unpack labelOptions
if isfield(labelOptions, 'xlabel')
    xlabelstr = labelOptions.xlabel ;
end
if isfield(labelOptions, 'ylabel')
    ylabelstr = labelOptions.ylabel ;
end
if isfield(labelOptions, 'title')
    titlestr = labelOptions.title ;
end
if isfield(labelOptions, 'label')
    label = labelOptions.label ;
end
if isfield(labelOptions, 'colormap')
    cmap = labelOptions.colormap ;
elseif isfield(labelOptions, 'cmap')
    cmap = labelOptions.cmap ;
end

%% Unpack options for style (diverging, positive, negative) and cmap
style = 'diverging' ;     % default is diverging
label_interpreter = 'Latex' ; % for colorbar label
if ~isempty(varargin)
    for i = 1:length(varargin)
        if isa(varargin{i},'double') 
            continue;
        end
        if isa(varargin{i},'logical')
            continue;
        end
        if ~isempty(regexp(varargin{i},'^[Ss]tyle','match'))
            stylestr = varargin{i+1} ;
        elseif ~isempty(regexp(varargin{i},'^[Ii]nterpreter','match'))
            label_interpreter = varargin{i+1} ;
        end
    end
    if ~isempty(regexp(stylestr,'^[Pp]hasemap','match'))
        style = 'phasemap' ;
    elseif ~isempty(regexp(stylestr,'^[Dd]iverging','match'))
        style = 'diverging' ;
    elseif ~isempty(regexp(stylestr,'^[Pp]ositive','match'))
        style = 'positive' ;
    elseif ~isempty(regexp(stylestr,'^[Nn]egative','match'))
        style = 'negative' ;
    end
    
end

% Show the image
h1 = imshow(im) ; hold on;
if isnumeric(xy_or_fxy)
    % Overlay the scalar field defined on vertices/xy points
    % Are xy supplied as two 1d arrays spanning a 2d outer product space?
    if size(xy_or_fxy, 1) == size(sf, 1) * size(sf, 2)
        h2 = imagesc(xy_or_fxy(:, 1), xy_or_fxy(:, 2), sf) ;
    else
        h2 = imagesc(xy_or_fxy(:, 1), xy_or_fxy(:, 2), sf) ;
    end
    if strcmp(style, 'phasemap')
        caxis(gca, [0, 2*pi]) ;
        if exist('cmap', 'var')
            colormap(cmap)
        else
            colormap(bwr) ;
        end
    elseif strcmp(style, 'diverging')
        if scale > 0
            caxis(gca, [-scale, scale]) ;
        else
            caxis(gca, [min(sf(:)), max(sf(:))])
        end
        if exist('cmap', 'var')
            colormap(cmap)
        else
            colormap(bwr) ;
        end
    elseif strcmp(style, 'positive')
        if scale > 0
            caxis(gca, [0, scale]) ;
        else
            caxis(gca, [0, max(sf(:))])
        end
        if exist('cmap', 'var')
            colormap(cmap) ;
        end
    elseif strcmp(style, 'negative')
        if scale > 0
            caxis(gca, [-scale, 0]) ;
        else
            caxis(gca, [min(sf(:)), 0]) ;
        end
        if exist('cmap', 'var')
            colormap(cmap) ;
        end
    end
elseif isa(xy_or_fxy, 'struct')
    % Add the scalar field defined on faces
    % Unpack xy_or_fxy
    if isfield(xy_or_fxy, 'faces')
        FF = xy_or_fxy.faces ;
    elseif isfield(xy_or_fxy, 'f')
        FF = xy_or_fxy.f ;
    else
        error('Face list for patches must be supplied as f or faces')
    end
    if isfield(xy_or_fxy, 'xy')
        V2D = xy_or_fxy.xy ;
    elseif isfield(xy_or_fxy, 'vertices')
        V2D = xy_or_fxy.vertices ;
    elseif isfield(xy_or_fxy, 'vertex')
        V2D = xy_or_fxy.vertex ;
    elseif isfield(xy_or_fxy, 'v')
        V2D = xy_or_fxy.v ;
    elseif isfield(xy_or_fxy, 'pts')
        V2D = xy_or_fxy.pts ;
    else
        error(['2d vertices of patches must be supplied as ', 
                'xy, v, vertex, vertices, or pts'])
    end
    
    % Create colors to paint patches
    if strcmp(style, 'phasemap')
        % Phasemap style
        cmap = phasemap ;
        colormap phasemap
        colors = mapValueToColor(sf, [0, 2*pi], cmap) ;
    elseif strcmp(style, 'nematic')
        % Phasemap style
        cmap = phasemap ;
        colormap phasemap
        colors = mapValueToColor(sf, [0, pi], cmap) ;
    elseif strcmp(style, 'diverging')
        % Diverging style
        if exist('cmap', 'var')
            colormap(cmap) ;
        else
            cmap = bwr ;
            colormap bwr
        end
        if scale > 0
            colors = mapValueToColor(sf, [-scale, scale], cmap) ;
            caxis(gca, [-scale, scale])
        else
            colors = mapValueToColor(sf, [min(sf(:)), max(sf(:))], cmap) ;
            caxis(gca, [min(sf(:)), max(sf(:))])
        end
    elseif strcmp(style, 'positive')
        % positive only
        if exist('cmap', 'var')
            colormap(cmap) ;
        else
            cmap = parula ;
            colormap parula
        end
        if scale < 0
            scale = max(sf(:)) ; 
        end
        colors = mapValueToColor(sf, [0, scale], cmap) ;
        caxis(gca, [0, scale])
    elseif strcmp(style, 'negative')
        % negative only
        if exist('cmap', 'var')
            colormap(cmap) ;
        else
            cmap = parula ;
            colormap parula
        end
        if scale < 0
            scale = min(sf(:)) ; 
        end
        colors = mapValueToColor(sf, [-scale, 0], cmap) ;
        caxis(gca, [-scale, 0])
    end
    h1 = patch( 'Faces', FF, 'Vertices', V2D, ...
        'FaceVertexCData', colors, 'FaceColor', 'flat', ...
        'EdgeColor', 'none') ;
end
alpha(alphaVal) ;

%% Labels
if ~isempty(titlestr)
    if ~strcmp(label_interpreter, 'default') && ~strcmp(label_interpreter, 'none')
        title(titlestr, 'Interpreter', label_interpreter) ;
    else
        title(titlestr) ;
    end
end
if ~isempty(xlabelstr)
    if ~strcmp(label_interpreter, 'default') && ~strcmp(label_interpreter, 'none')
        xlabel(xlabelstr, 'Interpreter', label_interpreter) ;
    else
        xlabel(xlabelstr) ;
    end
end
if ~isempty(ylabelstr)
    if ~strcmp(label_interpreter, 'default') && ~strcmp(label_interpreter, 'none')
        ylabel(ylabelstr, 'Interpreter', label_interpreter) ;
    else
        ylabel(ylabelstr) ;
    end
end

% Axis position
if isfield(labelOptions, 'axPosition')
    set(gca, 'Position', labelOptions.axPosition) ;
end


%% Colorbar settings
if isfield(labelOptions, 'cbPosition')
    c = colorbar('Position', labelOptions.cbPosition) ;
else
    c = colorbar();
end
% Manually flush the event queue and force MATLAB to render the colorbar
% necessary on some versions
drawnow
% Get the color data of the object that correponds to the colorbar
cdata = c.Face.Texture.CData;
% Change the 4th channel (alpha channel) to a fraction of initial value 
cdata(end,:) = uint8(alphaVal * cdata(end,:));
% Ensure that the display respects the alpha channel
c.Face.Texture.ColorType = 'truecoloralpha';
% Update the color data with the new transparency information
c.Face.Texture.CData = cdata;
if ~strcmp(label_interpreter, 'default') && ~strcmp(label_interpreter, 'none')
    c.Label.Interpreter = label_interpreter ;
end
c.Label.String = label ;
