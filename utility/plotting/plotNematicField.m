function [meshHandle, cbs] = plotNematicField(mag, theta, options)
% [meshHandle, cbs] = plotNematicField(mag, theta, options)
% 
% Plot nematic field on a mesh using a phasemap
%
% Parameters
% ----------
% options : optional struct with fields
%   clim_mag : float, colorlimit for magnitude
%   mesh : mesh on which to plot nematic field 
%   edgecolor : color spec
%   label/title : string
%   xyzlim/xyzlims
%   xlim : length 2 numeric, overwrites xyzlim(1, :)
%   ylim : length 2 numeric, overwrites xyzlim(2, :)
%   zlim : length 2 numeric, overwrites xyzlim(3, :)
%   cbarlabel : string 
%   makeCbar : bool 
%   axisOff : bool
%   axisOn : bool, ignored if axisOff specified
%   qsub / qsubsampling / quiverSubsampling : int, subsampling factor for
%       plotting nematic rods as quiverplot
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
interpreter = 'latex';
addQuiver = false ;
qsub = max(1, round(length(mag) / 200)) ;
qScale = 1 ;
cbarlabel = [] ;
if nargin > 2
    if isfield(options, 'clim_mag')
        clim_mag = options.clim_mag ;
    end
    if isfield(options, 'mesh')
        mesh = options.mesh ;
        if isempty(mesh) || ~isfield(mesh, 'f')
            reverseYAxis = true ;
        else
            reverseYAxis = false ;
        end
    else
        reverseYAxis = true ;
    end
    if isfield(options, 'reverseYAxis')
        reverseYAxis = options.reverseYAxis ;
    end
    if isfield(options, 'edgecolor')
        edgecolor = options.edgecolor ;
    end
    if isfield(options, 'title')
        label = options.title ;
    elseif isfield(options, 'label')
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
    if isfield(options, 'interpreter')
        interpreter = options.interpreter ;
    end
    if isfield(options, 'visible')
        set(gcf, 'visible', options.visible)
    end
    if isfield(options, 'addQuiver')
        addQuiver = options.addQuiver ;
    elseif isfield(options, 'quiver')
        addQuiver = options.quiver ;
    end
    if isfield(options, 'qsub')
        qsub = options.qsub ;
    elseif isfield(options, 'qsubsampling')
        qsub = options.qsubsampling ;
    elseif isfield(options, 'quiverSubsampling')
        qsub = options.quiverSubsampling ;
    end
    if isfield(options, 'qScale')
        qScale =  options.qScale ;
    end
    if isfield(options, 'cbarlabel')
        cbarlabel = options.cbarlabel ;
    end
else
    options = struct() ;
    reverseYaxis = true ;
end

if isempty(mesh)
    % [uu, vv] = meshgrid(linspace(0, 1, size(mag, 1)), ...
    %     linspace(0, 1, size(mag, 2))) ;
    [uu, vv] = meshgrid(1:size(mag, 2), 1:size(mag, 1)) ;
    mesh.v = [uu(:), vv(:), 0*uu(:)] ;
    mesh.f = defineFacesRectilinearGrid([], size(mag, 1), size(mag, 2)) ;
end
if ~isfield(mesh, 'f')
    mesh.f = defineFacesRectilinearGrid([], size(mag, 1), size(mag, 2)) ;
end
if size(mesh.v, 2) > 2
    is2d = all(mesh.v(:, 3) == 0) ;
elseif size(mesh.v, 2) == 2
    is2d = true ;
else
    error('mesh vertices must be 2d or 3d')
end

caxis([0, clim_mag])

% Intensity from mag and color from the theta
pm256 = phasemap(256)  ;
pm256 = pm256 / max(pm256(:)) ;
indx = max(1, round(mod(2*theta(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
colors = pm256(indx, :) ;
colors = min(mag(:) / clim_mag, 1) .* colors ;


if isa(mesh, 'struct')
    if is2d
        meshHandle = trisurf(mesh.f, mesh.v(:, 1), ...
            mesh.v(:, 2), 0 * mesh.v(:, 2), ...
            'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
    else
        meshHandle = trisurf(mesh.f, mesh.v(:, 1), ...
            mesh.v(:, 2), mesh.v(:, 3), ...
            'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
    end
elseif isa(mesh, 'cell')
    if is2d 
        meshHandle = trisurf(mesh{1}, mesh{2}(:, 1), ...
            mesh{2}(:, 2), 0 * mesh{2}(:, 1), ...
            'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
    else
        meshHandle = trisurf(mesh{1}, mesh{2}(:, 1), ...
            mesh{2}(:, 2), mesh{2}(:, 3), ...
            'FaceVertexCData', colors, 'edgecolor', edgecolor) ;
    end
end

%% QUIVER
% Add quiver to the plot
if addQuiver
    hold on;
    idx = 1:qsub:size(mesh.f, 1) ;
    bc = barycenter(mesh.v, mesh.f) ;
    
    if size(theta, 1) == size(mesh.f, 1)
        % Add quiver to subsampled faces
        if is2d
            % Add quiver to 2D faces
            Qf = mag .* [cos(theta), sin(theta), 0*theta] ;
            plot_on_faces = true ;
        else
            % Add quiver to 3D faces
            Qf = pushVectorField2Dto3DMesh(mag .* [cos(theta), sin(theta)], ...
                mesh.u, mesh.v, mesh.f, 1:length(mesh.f)) ;
            plot_on_faces = true ;
        end
    elseif size(theta, 1) == size(mesh.v, 1)
        if is2d
            % Just plot the bars on the vertices
            Qf = mag .* [cos(theta), sin(theta), 0*theta] ;
            plot_on_faces = false ;
        else
            % Push nematic to faces
            [V2F, ~] = meshAveragingOperators(mesh.f, mesh.v) ;
            thetaF = V2F * theta ;
            QfF = [cos(thetaF), sin(thetaF)] ;
            Qf = pushVectorField2Dto3DMesh(QfF, ...
                mesh.u, mesh.v, mesh.f, 1:length(mesh.f)) ;
            % normalize before multiplying by mag
            Qf = Qf ./ vecnorm(Qf, 2, 2) ;
            Qf = (V2F * mag) .* Qf ;
            % normalize
            plot_on_faces = true ;
        end
    else
        error('Could not add quiver since theta is not #faces or #vertices x 1 array')
    end
    
    if plot_on_faces
        % Place rods at 3D face barycenters
        quiver3(bc(idx, 1), bc(idx, 2), bc(idx, 3), ...
            Qf(idx, 1), Qf(idx, 2), Qf(idx, 3), qScale, 'w', 'ShowArrowHead', 'off') 
        quiver3(bc(idx, 1), bc(idx, 2), bc(idx, 3), ...
            -Qf(idx, 1), -Qf(idx, 2), -Qf(idx, 3), qScale, 'w', 'ShowArrowHead', 'off')
    else
        quiver(mesh.v(idx, 1), mesh.v(idx, 2), mesh.v(idx, 3),...
            Qf(idx, 1), Qf(idx, 2), Qf(idx, 3), qScale, 'w', 'ShowArrowHead', 'off') 
        quiver3(mesh.v(idx, 1), mesh.v(idx, 2), mesh.v(idx, 3), ...
            -Qf(idx, 1), -Qf(idx, 2), -Qf(idx, 3), qScale, 'w', 'ShowArrowHead', 'off')
            
    end
end

%% FORMATTING
axis equal
if ~isempty(label)
    title(label, 'Interpreter', 'Latex') ;
end

% Colorbar and phasewheel
if makeCbar
    colormap(gca, phasemap)
    cbs = cell(2, 1) ;
    cbs{1} = phasebar('colormap', phasemap, ...
        'location', [0.76, 0.05, 0.12, 0.135], ...
        'style', 'nematic') ;
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
    if isfield(options, 'cbarlabel')
        ylabel(cbs{2}, cbarlabel, 'interpreter', interpreter)
    end
    caxis([0, clim_mag])
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

% Axis labels
if exist('options', 'var')
    if isfield(options, 'xlabel')
        xlabel(options.xlabel, 'interpreter', interpreter)
    end
    if isfield(options, 'ylabel')
        ylabel(options.ylabel, 'interpreter', interpreter)
    end
    if isfield(options, 'zlabel')
        zlabel(options.zlabel, 'interpreter', interpreter)
    end
end

% Axis labels
if exist('options', 'var')
    if isfield(options, 'xyzlim')
        xlim(options.xyzlim(1, :))
        xlim(options.xyzlim(2, :))
        xlim(options.xyzlim(3, :))
    elseif isfield(options, 'xyzlims')
        xlim(options.xyzlims(1, :))
        xlim(options.xyzlims(2, :))
        xlim(options.xyzlims(3, :))
    end
    if isfield(options, 'xlim')
        xlim(options.xlim)
    end
    if isfield(options, 'ylim')
        ylim(options.ylim)
    end
    if isfield(options, 'zlim')
        zlim(options.zlim)
    end
end

% View from top if 2d
if is2d
    view(2)
end

% If no mesh supplied, treat the array like an image (like imagesc)
if reverseYAxis
    set(gca,'YDir','reverse')
end

