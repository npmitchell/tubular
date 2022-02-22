function [axs, cbs, meshes] = ...
    threeScalarFieldsOnSurface(mesh, sf1, sf2, sf3, options)
%threeScalarFieldsOnSurface(mesh, sf1, sf2, sf3, options)
%
% Inputs
% ------
% mesh : struct with fields f and v OR 2x1 cell array of faces and vertices
%   the mesh on which to plot the scalar fields
% sf1 : #vertices x 1 or #faces x 1 float array 
% sf2 : #vertices x 1 or #faces x 1 float array
% sf3 : #vertices x 1 or #faces x 1 float array
% options : optional struct with optional fields
%   clim : numeric or 2x1 numeric
%       colorlimit, by default set to [-clim, clim] if single value
%   clim1 : numeric, overwrites clim for axis1
%       colorlimit, by default set to [-clim1, clim1] if single value
%   clim2 : numeric, overwrites clim for axis2
%       colorlimit, by default set to [-clim2, clim2] if single value
%   clim3 : numeric, overwrites clim for axis3
%       colorlimit, by default set to [-clim3, clim3] if single value
%   axis1 : axis instance for sf1
%   axis2 : axis instance for sf2
%   axis3 : axis instance for sf3
%   edgecolor : color specifier
%   cmap : colormap for both axes
%   cmap1 : colormap for axis1, overwrites cmap
%   cmap2 : colormap for axis2, overwrites cmap
%   cmap3 : colormap for axis3, overwrites cmap
%   labels : cell of strings
%       titles for each subplot
% 
% Returns
% -------
% [axs, cbs, meshes] : handles for figure objects
%
% Example usage
% -------------
% labels = {'$\Re \mu$', '$\rho/2$', '$\Re\mu - \rho/2$'} ;
% options.labels = labels ;
% twoScalarFieldsOnSurface(mesh, sf1, sf2, sf3, options)
% sgtitle('deformation at t=0')  
%
% NPMitchell 2020

% Default options
clim1 = max(sf1(:)) ; 
clim2 = max(sf2(:)) ;
clim3 = max(sf3(:)) ;
    
% Define axes
set(gcf, 'visible', 'off')
if isfield(options, 'axis1')
    ax1 = options.axis1 ;
else
    ax1 = subplot(1, 3, 1) ;
end
if isfield(options, 'axis2')
    ax2 = options.axis2 ;
else
    ax2 = subplot(1, 3, 2) ;
end
if isfield(options, 'axis3')
    ax3 = options.axis3 ;
else
    ax3 = subplot(1, 3, 3) ;
end

% Define climits
if isfield(options, 'clim')
    clim1 = options.clim ;
    clim2 = options.clim ;
    clim3 = options.clim ;
end
if isfield(options, 'clim1')
    clim1 = options.clim1 ;
elseif ~isfield(options, 'clim')
    clim1 = max(abs(sf1(:))) ;
end
if isfield(options, 'clim2')
    clim2 = options.clim2 ;
elseif ~isfield(options, 'clim')
    clim2 = max(abs(sf2(:))) ;
end
if isfield(options, 'clim3')
    clim3 = options.clim3 ;
elseif ~isfield(options, 'clim')
    clim3 = max(abs(sf3(:))) ;
end

% Define
edgecolor = 'none' ;
if isfield(options, 'edgecolor')
    edgecolor = options.edgecolor ;
end

% Define colormaps
if isfield(options, 'colormap')
    cmap = options.colormap ;
else
    caxis([-1, 1])
    cmap = blueblackred(256) ;
end
if isfield(options, 'colormap1')
    cmap1 = options.colormap1 ;
else
    cmap1 = cmap ;
end
if isfield(options, 'colormap2')
    cmap2 = options.colormap2 ;
else
    cmap2 = cmap ;
end
if isfield(options, 'colormap3')
    cmap3 = options.colormap3 ;
else
    cmap3 = cmap ;
end

if isfield(options, 'labels')
    labels = options.labels ;
else
    labels = [] ;
end

%% Panel 1
set(gcf,'CurrentAxes', ax1)
if isa(mesh, 'struct')
    mesh1 = trisurf(mesh.f, mesh.v(:, 1), ...
        mesh.v(:, 2), mesh.v(:, 3), ...
        'FaceVertexCData', sf1(:), 'edgecolor', edgecolor) ;
elseif isa(mesh, 'cell')    
    mesh1 = trisurf(mesh{1}, mesh{2}(:, 1), ...
        mesh{2}(:, 2), mesh{2}(:, 3), ...
        'FaceVertexCData', sf1(:), 'edgecolor', edgecolor) ;
end
axis equal
cb1 = colorbar('location', 'southOutside') ;
caxis([-clim1, clim1])
if ~isempty(labels)
    title(labels{1}, 'Interpreter', 'Latex')   
end
colormap(cmap1)
% axis off
% view(2)

%% Panel 2 
set(gcf,'CurrentAxes', ax2)
if isa(mesh, 'struct')
    mesh2 = trisurf(mesh.f, mesh.v(:, 1), ...
        mesh.v(:, 2), mesh.v(:, 3), ...
        'FaceVertexCData', sf2(:), 'edgecolor', edgecolor) ;
elseif isa(mesh, 'cell')    
    mesh2 = trisurf(mesh{1}, mesh{2}(:, 1), ...
        mesh{2}(:, 2), mesh{2}(:, 3), ...
        'FaceVertexCData', sf2(:), 'edgecolor', edgecolor) ;
end
axis equal
cb2 = colorbar('location', 'southOutside') ;
caxis([-clim2, clim2])
if ~isempty(labels)
    title(labels{2}, 'Interpreter', 'Latex')   
end
colormap(cmap2)

%% Panel 3
set(gcf,'CurrentAxes', ax3)
if isa(mesh, 'struct')
    mesh3 = trisurf(mesh.f, mesh.v(:, 1), ...
        mesh.v(:, 2), mesh.v(:, 3), ...
        'FaceVertexCData', sf3(:), 'edgecolor', edgecolor) ;
elseif isa(mesh, 'cell')    
    mesh3 = trisurf(mesh{1}, mesh{2}(:, 1), ...
        mesh{2}(:, 2), mesh{2}(:, 3), ...
        'FaceVertexCData', sf3(:), 'edgecolor', edgecolor) ;
end
axis equal
cb3 = colorbar('location', 'southOutside') ;
caxis([-clim3, clim3])
if ~isempty(labels)
    title(labels{3}, 'Interpreter', 'Latex')   
end
colormap(cmap3)

%% Package output
axs = [ax1, ax2, ax3] ;
cbs = [cb1, cb2, cb3] ;
meshes = [mesh1, mesh2, mesh3] ;

end