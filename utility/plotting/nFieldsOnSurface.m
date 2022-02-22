function [axs, cbs, meshHandles] = ...
    nFieldsOnSurface(meshes, fields, options)
%[axs, cbs, meshHandles] = nFieldsOnSurface(meshes, fields, options)
%
% Inputs
% ------
% meshes : struct with fields f and v 
%          OR 2x1 cell array of faces and vertices
%          OR triangulation object with fields ConnectivityList and Points
%   The meshes on which to plot the scalar fields. Can be 2D or 3D.
% fields : cell array of (#vertices x 1) or (#faces x 1) float array, 
%   or length2 cell of magnitude and angle for nematic fields
%   The fields to plot on the surfaces.
% options : optional struct with optional fields
%   clim : numeric or 2x1 numeric
%       colorlimit, by default set to [-clim, clim] if single value
%   clims : cell of length-1 or length-2 numeric arrays, overwrites clim 
%       colorlimit for each field, by default set to [-clim1, clim1] if 
%       single value
%   axs : axis instances for each field
%   edgecolor : color specifier (default='none')
%       edge color for trisurf
%   cmap : colormap for both axes (default='blueblackred')
%       if a single colormap governs all axes, supply here
%   cmaps : cell array of colormaps
%       colormap for each axis, overwrites cmap
%   labels : cell of strings
%       titles for each subplot
%   makeCbar : nfields x 1 bool array or cell array
%       whether to make colorbar for each axis
%   masterCbar : bool 
%       make a single colorbar ruling all axes
%   xyzlims : 3x2 numeric or (nfields)x1 cell array of 3x2 numeric arrays
%       the xyz limits for all panels as [xmin,xmax;ymin,ymax;zmin,zmax]
%   xlim : length 2 numeric array or length nfields cell array of length 2
%       numeric arrays
%       xlim values if xyzlims not provided
%   ylim : length 2 numeric array or length nfields cell array of length 2
%       numeric arrays
%       ylim values if xyzlims not provided
%   zlim : length 2 numeric array or length nfields cell array of length 2
%       numeric arrays
%       zlim values if xyzlims not provided
%   view : length nfields cell array of length 2 arrays
%       viewing angle for each panel
%   style : 'diverging', 'positive', 'negative'
%       color limit style for scalar fields
%   
%   
% Returns
% -------
% [axs, cbs, meshHandles] : cell arrays of handles for figure objects
%
% Example usage
% -------------
% options.labels = {'$\Re \mu$', '$\rho/2$', '', ''} ;
% options.makeCbar = [false, false, true, true];
% options.axisOff = 'true' ;
% options.view = {[0,0], [0,0],[0,90],[0,90]} ;
% xyzlim = [ -22.4450, 260.0300 ;  -76.5010, 65.4510 ; -59.8050, 80.7680];
% options.xyzlims = {xyzlim,xyzlim, [0,1;0,1;-1,1], [0,1;0,1;-1,1]} ;
% [axs, cbs] = ...
%     nFieldsOnSurface({m3, m3, m2, m2}, ...
%     {real(mu), imag(mu), real(mu), imag(mu)}, options) ;
% expandSecondAxesRow(axs, -0.05)
% sgtitle('deformation at t=0')  
%
% NPMitchell 2020

%% Default options
nfields = length(fields) ;
axisOff = false ;
phasebarPosition = [0.82, 0.12, 0.1, 0.135] ;
visible = 'off' ;
style = 'diverging' ;
interpreter = 'latex' ;
nematic_or_polar = 'nematic' ;
if nargin < 3
    options = struct() ;
end
if isfield(options, 'visible')
    visible = options.visible ;
end
if isfield(options, 'subplotGrouping')
    subplotGrouping = options.subplotGrouping ;
else
    subplotGrouping = [] ;
end
if isfield(options, 'polarStyle')
    nematic_or_polar = lower(options.polarStyle) ;
    if ~strcmpi(nematic_or_polar, 'nematic') && ...
            ~strcmpi(nematic_or_polar, 'polar')
        error('polarStyle must be <nematic> or <polar>')
    end
end

%% Define axes
set(gcf, 'visible', visible)
if isfield(options, 'axs')
    axs = options.axs ;
else
    axs = cell(nfields, 1) ;
    if ~isempty(subplotGrouping)
        for qq = 1:nfields
            axs{qq} = subplot(subplotGrouping(1), subplotGrouping(2), qq) ;
        end
    else    
        if nfields < 4
            for qq = 1:nfields
                axs{qq} = subplot(1, nfields, qq) ;
            end
        elseif nfields == 4
            for qq = 1:nfields
                axs{qq} = subplot(2, 2, qq) ;
            end
        elseif nfields < 7
            for qq = 1:nfields
                axs{qq} = subplot(2, 3, qq) ;
            end
        else
            error('code for default axis arrangement for this number of axes')
        end
    end
end

% Define climits
if isfield(options, 'clim')
    clim = options.clim ;
elseif isfield(options, 'climit')
    clim = options.climit ;
end
if isfield(options, 'clims')
    clims = options.clims ;
end
if isfield(options, 'style')
    style = options.style ;
end

% Define xyzlims
if isfield(options, 'xlim')
    if isa(options.xlim, 'cell')
        try
            assert(numel(options.xlim) == nfields)
        catch
            error('number of elements in options.xlim ~= nfields')
        end
        xlimits = options.xlim ;
    else
        xlimit = options.xlim ;
    end
end
if isfield(options, 'ylim')
    if isa(options.ylim, 'cell')
        try
            assert(numel(options.ylim) == nfields)
        catch
            error('number of elements in options.ylim ~= nfields')
        end
        ylimits = options.ylim ;
    else
        ylimit = options.ylim ;
    end
end
if isfield(options, 'zlim')
    if isa(options.zlim, 'cell')
        try
            assert(numel(options.zlim) == nfields)
        catch
            error('number of elements in options.zlim ~= nfields')
        end
        zlimits = options.zlim ;
    else
        zlimit = options.zlim ;
    end
end
if isfield(options, 'xyzlims')
    if isa(options.xyzlims, 'cell')
        try
            assert(numel(options.xyzlims) == nfields)
        catch
            error('number of elements in options.xyzlim ~= nfields')
        end
        xlimits = {options.xyzlims{1}(1, :)} ;
        ylimits = {options.xyzlims{1}(2, :)} ;
        zlimits = {options.xyzlims{1}(3, :)} ;
        for qq = 2:nfields
            xlimits = cat(1, xlimits, options.xyzlims{qq}(1, :)) ;
            ylimits = cat(1, ylimits, options.xyzlims{qq}(2, :)) ;
            zlimits = cat(1, zlimits, options.xyzlims{qq}(3, :)) ;
        end
    else
        xlimit = options.xyzlims(1, :) ;
        ylimit = options.xyzlims(2, :) ;
        zlimit = options.xyzlims(3, :) ;
    end
end

% Unpack view/views
if isfield(options, 'views')
    views = options.views ;
elseif isfield(options, 'view')
    views = options.view ;
end

% axis off?
if isfield(options, 'axisOff')
    axisOff = options.axisOff ;
end
if isfield(options, 'axisOn')
    axisOff = ~options.axisOn ;
end

% Define edgecolor
edgecolor = 'none' ;
if isfield(options, 'edgecolor')
    edgecolor = options.edgecolor ;
end

% Define colormaps
if isfield(options, 'colormap')
    cmap = options.colormap ;
elseif isfield(options, 'cmap')
    cmap = options.cmap ;
else
    caxis([-1, 1])
    cmap = blueblackred(256) ;
end
if isfield(options, 'colormaps')
    cmaps = options.colormaps ;
elseif isfield(options, 'cmaps')
    cmaps = options.cmaps ;
end

if isfield(options, 'labels')
    labels = options.labels ;
else
    labels = [] ;
end
if isfield(options, 'interpreter')
    interpreter = options.interpreter ;
end

if isfield(options, 'makeCbar')
    % grab boolean (or numeric with 0 and 1) array of whether to make a
    % colorbar for each axis
    makeCbar = options.makeCbar ;
    if isa(makeCbar, 'cell')
        makeCbar = cell2mat(options.makeCbar) ;
    end
    assert(length(makeCbar) == 1 || length(makeCbar) == nfields)
else
    if isfield(options, 'masterCbar')
        masterCbar = options.masterCbar ;
        if ~masterCbar 
            makeCbar = true(nfields, 1) ;
        else
            makeCbar = false(nfields, 1) ;
        end
    else
        makeCbar = true(nfields, 1) ;
    end
end

if isfield(options, 'phasebarPosition')
   phasebarPosition = options.phasebarPosition ;
end

if isfield(options, 'masterCbar')
    masterCbar = options.masterCbar ;
    if masterCbar
        if isfield(options, 'masterCbarPosition')
            masterCbarPosition = options.masterCbarPosition ;
        else
            masterCbarPosition = [0.39, 0.05, 0.22, 0.0341];
        end
    end
else
    masterCbar = false ;
end


%% Panels
for qq = 1:nfields
    % Unpack meshes if there are multiple
    if isa(meshes, 'cell')
        if isa(meshes{qq}, 'struct') || ...
                isa(meshes{qq}, 'cell') || ...
                isa(meshes{qq}, 'triangulation')
            mesh = meshes{qq} ;
        else
            mesh = meshes ;
        end
    else
        mesh = meshes ;
    end
    
    set(gcf,'CurrentAxes', axs{qq})
    
    % Check if scalar field
    if ~isa(fields{qq}, 'cell')
        % Scalar field
        if isa(mesh, 'struct')
            if size(mesh.v, 2) == 3
                meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                    mesh.v(:, 2), mesh.v(:, 3), ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            elseif size(mesh.v, 2) == 2
                meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                    mesh.v(:, 2), 0*mesh.v(:, 1), ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            end
        elseif isa(mesh, 'cell')    
            if size(mesh{2}, 2) == 3
                meshHandles{qq} = trisurf(mesh{1}, mesh{2}(:, 1), ...
                    mesh{2}(:, 2), mesh{2}(:, 3), ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            elseif size(mesh{2}, 2) == 2
                meshHandles{qq} = trisurf(mesh{1}, mesh{2}(:, 1), ...
                    mesh{2}(:, 2), 0*mesh{2}(:, 1), ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            end
        elseif isa(mesh, 'triangulation') 
            if size(mesh.Points, 2) == 3
                meshHandles{qq} = trisurf(mesh, ...
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            elseif size(mesh.Points, 2) == 2
                meshHandles{qq} = trisurf(mesh.ConnectivityList, ...
                    mesh.Points(:, 1), mesh.Points(:, 2), 0*mesh.Points(:, 1), ... 
                    'FaceVertexCData', fields{qq}(:), 'edgecolor', edgecolor) ;
            end    
        end
        axis equal
        
        if makeCbar(qq)
            cbs{qq} = colorbar('location', 'southOutside') ;
        end
        
        if exist('clims', 'var')
            if numel(clims{qq}) == 1
                caxis([-clims{qq}, clims{qq}])
            elseif numel(clims{qq}) == 2
                caxis([clims{qq}(1), clims{qq}(2)])
            end
        elseif exist('clim', 'var')
            if numel(clim) == 1
                caxis([-clim, clim])
            elseif numel(clim) == 2
                caxis([clim(1), clim(2)])
            else
                error('clim provided has > 2 elements')
            end
        elseif strcmpi(style, 'diverging')
            caxis([-max(abs(fields{qq})), max(abs(fields{qq}))])
        elseif strcmpi(style, 'positive')
            caxis([0, max(abs(fields{qq}))])
        elseif strcmpi(style, 'negative')
            caxis([-max(abs(fields{qq})), 0])
        else
            error(['unrecognized style=' style])
        end
        
        if makeCbar(qq) && isfield(options, 'cbarlabels')
            % xlabel(cbs{qq}, options.cbarlabels{qq}, 'interpreter', 'latex') 
            %set(cbs{qq}{2}, 'xlabel', options.cbarlabels{qq})
            ylabel(cbs{qq}, options.cbarlabels{qq}, 'interpreter', interpreter)
            
        end
        
        if ~isempty(labels)
            if length(labels) > qq - 1
                if ~isempty(labels{qq})
                    title(labels{qq}, 'Interpreter', 'Latex')   
                end
            end
        end
        if isfield('cmaps', 'var')
            colormap(axs{qq}, cmaps{qq})
        else
            colormap(axs{qq}, cmap)
        end
    else
        % Nematic field or polar field
        % Unpack field
        dev = fields{qq}{1} ;
        theta = fields{qq}{2} ;
        
        if exist('clims', 'var')
            if numel(clims{qq}) == 1
                clim_dev = clims{qq} ;
            elseif numel(clims{qq}) == 2
                clim_dev = clims{qq}(2) ;
            elseif numel(clims{qq}) == 0
                clim_dev = max(abs(dev(:))) ;
            end
        elseif exist('clim', 'var')
            if numel(clim) == 1
                clim_dev = clim ;
            elseif numel(clim) == 2
                clim_dev = clim(2) ;
            elseif numel(clims{qq}) == 0
                clim_dev = max(abs(dev(:))) ;
            else
                error('clim provided has > 2 elements')
            end
        else
            clim_dev = max(abs(dev(:))) ;
        end
        caxis([0, clim_dev])
        
        % Intensity from dev and color from the theta
        pm256 = phasemap(256) ;
        if strcmpi(nematic_or_polar, 'nematic')
            indx = max(1, round(mod(2*theta(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
        else
            % Range of colormap is -pi to pi, so add pi to index at zero!
            indx = max(1, round(mod(theta(:)+pi, 2*pi)*size(pm256, 1)/(2 * pi))) ;
        end
        colors = pm256(indx, :) ;
        colors = min(dev(:) / clim_dev, 1) .* colors ;
        
        if isa(mesh, 'struct')
            if size(mesh.v, 2) == 3
                meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                    mesh.v(:, 2), mesh.v(:, 3), ...
                    'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
            elseif size(mesh.v, 2) == 2
                meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                    mesh.v(:, 2), 0*mesh.v(:, 1), ...
                    'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
            else
                error('mesh.v must be 2d or 3d list of vector coordinates')
            end
        elseif isa(mesh, 'cell')    
            if size(mesh{2}, 2) == 3
                meshHandles{qq} = trisurf(mesh{1}, mesh{2}(:, 1), ...
                    mesh{2}(:, 2), mesh{2}(:, 3), ...
                    'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
            elseif size(mesh.v, 2) == 2
                meshHandles{qq} = trisurf(mesh.f, mesh.v(:, 1), ...
                    mesh.v(:, 2), 0*mesh.v(:, 1), ...
                    'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
            else
                error('mesh.v must be 2d or 3d list of vector coordinates')
            end
        end
        axis equal
        if ~isempty(labels)
            title(labels{qq}, 'Interpreter', 'Latex')   
        end

        % Colorbar and phasewheel
        if makeCbar(qq)
            colormap(gca, phasemap)
            cbs{qq} = cell(2, 1) ;
            cbs{qq}{1} = phasebar('colormap', phasemap, ...
                'location', phasebarPosition, 'style', nematic_or_polar) ;
            shrink = max(0.6 - 0.1 * (mod(nfields, 3)-2), 0.1) ;
            % axis off
            % view(2)
            cbs{qq}{2} = colorbar('location', 'southOutside') ;
            drawnow
            axpos = get(axs{qq}, 'position') ;
            cbpos = get(cbs{qq}{2}, 'position') ;
            set(cbs{qq}{2}, 'position', [cbpos(1), cbpos(2), cbpos(3)*shrink, cbpos(4)])
            set(axs{qq}, 'position', axpos) ;
            hold on;
            
            % label the colorbar
            if isfield(options, 'cbarlabels')
                %set(cbs{qq}{2}, 'xlabel', options.cbarlabels{qq})
                ylabel(cbs{qq}{2}, options.cbarlabels{qq}, 'interpreter', interpreter)
            end
        end
        
        % Set colorlimits
        if exist('clims', 'var')
            if numel(clims{qq}) == 1
                caxis([0, clims{qq}])
            elseif numel(clims{qq}) == 2
                caxis([clims{qq}(1), clims{qq}(2)])
            end
        elseif exist('clim', 'var')
            if numel(clim) == 1
                caxis([0, clim])
            elseif numel(clim) == 2
                caxis([clim(1), clim(2)])
            else
                error('clim provided has > 2 elements')
            end
        end
        colormap(gca, gray) 
        
        if ~isempty(labels)
            title(labels{qq}, 'Interpreter', 'Latex')   
        end
    end
    
    % Set view from "views"
    if exist('views', 'var')
        if isa(views, 'cell')
            if ~isempty(views{qq})
                view(views{qq}(1), views{qq}(2))
            end
        elseif numel(views) == 2
            view(views(1), views(2))
        else
            error('Number of elements in options.view is not 2 or nfields')
        end
    end
    
    % Are axes off?
    if axisOff
        axis off
    end
    
    % Set xyzlims
    if exist('xlimits', 'var')
        xlim(xlimits{qq}) ;
    elseif exist('xlimit', 'var')
        xlim(xlimit) ;
    end
    if exist('ylimits', 'var')
        ylim(ylimits{qq}) ;
    elseif exist('ylimit', 'var')
        ylim(ylimit) ;
    end
    if exist('zlimits', 'var')
        zlim(zlimits{qq}) ;
    elseif exist('zlimit', 'var')
        zlim(zlimit) ;
    end
end

% Master colorbar
if masterCbar
    axpos = get(gca, 'position') ;
    if exist('cbs', 'var')
        cbs{length(cbs)+1} = colorbar('location', 'southOutside') ;
    else
        cbs{1} = colorbar('location', 'southOutside') ;
    end
    set(gca, 'position', axpos) ;
    set(cbs{length(cbs)}, 'position', masterCbarPosition) ;
end

end