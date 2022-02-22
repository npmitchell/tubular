function [ax1, ax2, cb1, cb2, pbar, meshTr, meshDev] = ...
    plotTraceDeviator(mesh, tre, dev, theta, options)
%plotTraceDeviator(mesh, tre, dev, theta, options)
%
% Inputs
% ------
% mesh : 
% tre : #vertices x 1 or #faces x 1 float array 
% dev : #vertices x 1 or #faces x 1 float array
% theta : #vertices x 1 or #faces x 1 float array
% options : struct with fields
%   clim_trace
%   clim_deviator
%   
% 
% Returns
% -------
%
% Example usage
% -------------
% labels = {'$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon] $', ...
%     '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$'} ;
% options.labels = labels ;
% plotTraceDeviator(mesh, tre, dev, theta, options)
% sgtitle('strain t=0')  
%
% NPMitchell 2020

% Default options
clim_trace = max(abs(0.5 * tre(:))) ; 
clim_deviator = max(abs(dev(:))) ;
if all(dev(:) == clim_deviator)
    clim_deviator = 1 ;
end
    
% Define axes
if isfield(options, 'axis1')
    ax1 = options.axis1 ;
else
    ax1 = subplot(1, 2, 1) ;
end
if isfield(options, 'axis2')
    ax2 = options.axis2 ;
else
    ax2 = subplot(1, 2, 2) ;
end

% Define climits
if isfield(options, 'clim_trace')
    clim_trace = options.clim_trace ;
end
if isfield(options, 'clim_deviator')
    clim_deviator = options.clim_deviator ;
end

% Define
edgecolor = 'none' ;
if isfield(options, 'edgecolor')
    edgecolor = options.edgecolor ;
end

% Define colormaps
if isfield(options, 'colormap_trace')
    cmapTr = options.colormap_trace ;
end
if isfield(options, 'colormap_deviator')
    cmapDev = options.colormap_deviator ;
else
    cmapDev = phasemap(256) ;
end

if isfield(options, 'labels')
    labels = options.labels ;
else
    labels = [] ;
end

% Panel 1
set(gcf,'CurrentAxes', ax1)
% indx = min(length(cmapTr), max(1, round(0.5 * tre(:) / cliom)*size(cmapDev, 1))) ;
% colors = cmapTr(indx, :) ;
meshTr = trisurf(mesh.f, mesh.v(:, 1), ...
    mesh.v(:, 2), mesh.v(:, 3), ...
    'FaceVertexCData', 0.5 * tre(:), 'edgecolor', edgecolor) ;
% daspect([1,1,1])
axis equal
cb1 = colorbar('location', 'southOutside') ;
caxis([-clim_trace, clim_trace])
if ~isempty(labels)
    title(labels{1}, 'Interpreter', 'Latex')   
end
if exist('cmapTr', 'var')
    colormap(cmapTr)
else
    colormap(blueblackred(256)) ;
end
% axis off
% view(2)

% Panel 2 
set(gcf,'CurrentAxes',ax2)
% Intensity from dev and color from the theta
indx = max(1, round(mod(2*theta(:), 2*pi)*size(cmapDev, 1)/(2 * pi))) ;
colors = cmapDev(indx, :) ;
colors = min(dev(:) / clim_deviator, 1) .* colors ;
meshDev = trisurf(mesh.f, mesh.v(:, 1), ...
    mesh.v(:, 2), mesh.v(:, 3), ...
    'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
% daspect([1,1,1])
axis equal
if ~isempty(labels)
    title(labels{2}, 'Interpreter', 'Latex')   
end

% Colorbar and phasewheel
colormap(gca, phasemap)
pbar = phasebar('colormap', phasemap, ...
    'location', [0.82, 0.12, 0.1, 0.135], 'style', 'nematic') ;
% axis off
% view(2)
cb2 = colorbar('location', 'southOutside') ;
drawnow
axpos = get(ax2, 'position') ;
cbpos = get(cb2, 'position') ;
set(cb2, 'position', [cbpos(1), cbpos(2), cbpos(3)*0.6, cbpos(4)])
set(ax2, 'position', axpos) ;
hold on;
caxis([0, clim_deviator])
colormap(gca, gray)

end