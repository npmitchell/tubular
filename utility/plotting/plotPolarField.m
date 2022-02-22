function [meshHandle, cbs] = plotPolarField(mag, theta, options)
% plotPolarField()
% 
%
% Parameters
% ----------
%
% Returns
% -------
%
% NPMitchell 2021

clim_mag = max(mag(:)) ;
mesh = [] ;
edgecolor = 'none' ;
label = '' ;
makeCbar = true ;
axisOff = false ;
if nargin > 2
    if isfield(options, 'clim_mag')
        clim_mag = options.clim_mag ;
    elseif isfield(options, 'clim')
        clim_mag = options.clim ;
    elseif isfield(options, 'climit')
        clim_mag = options.climit ;    
    end
    if isfield(options, 'mesh')
        mesh = options.mesh ;
    end
    if isfield(options, 'edgecolor')
        edgecolor = options.edgecolor ;
    end
    if isfield(options, 'label')
        label = options.label ;
    end
    if isfield(options, 'makeCbar')
        makeCbar = options.makeCbar ;
    end
    if isfield(options, 'axisOff')
        axisOff = options.axisOff ;
    elseif isfield(options, 'axisOn')
        axisOff = ~ options.axisOn ;
    end
end

if isempty(mesh)
    % OLD version:
    % [uu, vv] = meshgrid(linspace(0, 1, size(mag, 1)), ...
    %     linspace(0, 1, size(mag, 2))) ;
    % mesh.v = [uu(:), vv(:), 0*uu(:)] ;
    
    % NEW version: 
    [uu, vv] = meshgrid(1:size(mag, 2), 1:size(mag, 1)) ;
    mesh.v = [uu(:), vv(:), 0*uu(:)] ;
    mesh.f = defineFacesRectilinearGrid([], size(mag, 1), size(mag, 2)) ;
end
if ~isfield(mesh, 'f')
    mesh.f = defineFacesRectilinearGrid([], size(mag, 1), size(mag, 2)) ;
end
is2d = all(mesh.v(:, 3) == 0) ;


caxis([0, clim_mag])

% Intensity from mag and color from the theta
pm256 = phasemap(256) ;
% NOTE: The phasebar ranges from -pi to pi, so add pi to put zero in the
% middle of the colormap!
indx = max(1, round(mod(theta(:)+pi, 2*pi)*size(pm256, 1)/(2 * pi))) ;
colors = pm256(indx, :) ;
colors = min(mag(:) / clim_mag, 1) .* colors ;


if isa(mesh, 'struct')
    meshHandle = trisurf(mesh.f, mesh.v(:, 1), ...
        mesh.v(:, 2), mesh.v(:, 3), ...
        'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
elseif isa(mesh, 'cell')    
    meshHandle = trisurf(mesh{1}, mesh{2}(:, 1), ...
        mesh{2}(:, 2), mesh{2}(:, 3), ...
        'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
end

axis equal
if ~isempty(label)
    title(label, 'Interpreter', 'Latex')   
end

% Colorbar and phasewheel
if makeCbar
    colormap(gca, phasemap)
    cbs = cell(2, 1) ;
    cbs{1} = phasebar('colormap', phasemap, ...
        'location', [0.76, 0.05, 0.12, 0.135], 'style', 'polar') ;
    shrink = 0.6 ;
    
    if axisOff
        axis off
    end
    
    cbs{2} = colorbar('location', 'southOutside') ;
    drawnow
    axpos = get(gca, 'position') ;
    cbpos = get(cbs{2}, 'position') ;
    set(cbs{2}, 'position', [cbpos(1), cbpos(2), cbpos(3)*shrink, cbpos(4)])
    set(gca, 'position', axpos) ;
    hold on;

    % label the colorbar
    if isfield(options, 'cbarlabels')
        %set(cbs{qq}{2}, 'xlabel', options.cbarlabels{qq})
        ylabel(cbs{qq}{2}, options.cbarlabels{qq}, 'interpreter', interpreter)
    end
end

% Set colorlimits
if exist('clim', 'var')
    error('why clim? should use clim_mag.')
    if numel(clim) == 1
        caxis([0, clim])
    elseif numel(clim) == 2
        caxis([clim(1), clim(2)])
    else
        error('clim provided has > 2 elements')
    end
else
    caxis([0, clim_mag])
end
colormap(gca, gray) 

if is2d
    view(2)
end

