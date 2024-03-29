function plotStrainRateTimePoint(tubi, tp, options)
%plotStrainRateTimePoint(tubi, tp, options)
%   Plot the traceful and traceless components of the strain rate tensor
%   defined on each face, using vertex-based results loaded from disk. 
%   Note that the vertex-based results are more heavily smoothed (via
%   spectral filtering) than the face-based results on disk.
%
% Parameters
% ----------
% tubi : TubULAR class instance
% tp : int 
%   timepoint in units of (1/tubi.timeInterval) * tubi.timeUnits
% options: struct with fields
%   overwrite : bool
%       Overwrite existing results on disk
%   mesh : struct with fields f,v
%       current Mesh, allowed to be supplied to allow speedup (prevents
%       delay from loading cutMesh from disk)
%   cutMesh : struct with fields f,v, pathPairs, nU, nV
%       current cutMesh, allowed to be supplied to allow speedup (prevents
%       delay from loading cutMesh from disk)
%   clim_trace : 2x1 numeric
%       color limits for the colormap/bar of the isotropic strain rate
%   clim_deviatoric : 2x1 numeric
%       color limits for the colormap/bar of the deviatoric strain rate
%   samplingResolution : '1x' or 'resampled'
%       currently only 1x resolution is allowed
%   averagingStyle : 'Lagrangian' or 'none'
%       currently we have restricted the velocity averaging style in
%       tubular to either be Lagrangian averaging or none
%   debug : bool
%       if true, show intermediate results
%   plot_comparison : bool
%       plot a comparison between trace of measured strain rate and
%       measured DEC divergence minus H*normal velocity, where H is the
%       mean curvature. 
% Additional optional fields of options that are already specified by tubi
% metadata but can be overwritten by including them in options are:
%     lambda = tubi.smoothing.lambda ;
%     lambda_mesh = tubi.smoothing.lambda_mesh ;
%     nmodes = tubi.smoothing.nmodes ;
%     zwidth = tubi.smoothing.zwidth ;
%   
% 
% NPMitchell 2020

tidx = tubi.xp.tIdx(tp) ;

%% Unpack required params
lambda = tubi.smoothing.lambda ;
lambda_mesh = tubi.smoothing.lambda_mesh ;
nmodes = tubi.smoothing.nmodes ;
zwidth = tubi.smoothing.zwidth ;

if nargin < 3
    options = struct() ;
else
    if isfield(options, 'lambda')
        lambda = options.lambda ;
    end
    if isfield(options, 'lambda')
        lambda_mesh = options.lambda_mesh ;
    end
    if isfield(options, 'nmodes')
        nmodes = options.nmodes ;
    end
    if isfield(options, 'zwidth')
        zwidth = options.zwidth ;
    end
end
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
debug = false ;
plot_comparison = true ;

%% Parameters
overwrite = false ;
clim_trace = 0.05 ;
clim_deviatoric = 0.05 ;
averagingStyle = 'Lagrangian' ;
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'mesh')
    mesh = options.mesh ;
end
if isfield(options, 'cutMesh')
    cutMesh = options.cutMesh ;
end
if isfield(options, 'clim_trace')
    clim_trace = options.clim_trace ;
end
if isfield(options, 'clim_deviatoric')
    clim_deviatoric = options.clim_deviatoric ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'plot_comparison')
    plot_comparison = options.plot_comparison ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
else
    error("Could not parse samplingResolution: set to '1x'")
end


%% Unpack tubi
t0 = tubi.t0set() ;
tubi.getXYZLims ;
xyzlim = tubi.plotting.xyzlim_um ;


%% Output directory

% build tubi.dir.strainRate.smoothing anew to allow for periods in the name
l_lmesh = strrep(sprintf('lambda%0.03f_lmesh%0.3f_modes%02dw%02d', ...
    lambda, lambda_mesh, nmodes, zwidth), '.', 'p') ;

egImDir = ...
    fullfile(tubi.dir.uvCoord, 'strainRate', l_lmesh) ;
measurementDir = ...
    fullfile(tubi.dir.uvCoord, 'strainRate', l_lmesh, 'measurements') ;
if ~exist(egImDir, 'dir')
    mkdir(egImDir) ;
end
if ~exist(measurementDir, 'dir')
    mkdir(measurementDir) ;
end

buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;

%% load from tubi
if doubleResolution
    nU = tubi.nU * 2 - 1 ;
    nV = tubi.nV * 2 - 1 ;
else
    nU = tubi.nU ;
    nV = tubi.nV ;    
end

%% load the metric strain
% Define metric strain filename        
if ~isfield(options, 'tre') || ~isfield(options, 'dev') || ...
        ~isfield(options, 'theta')
    estrainFn = fullfile(measurementDir, sprintf(tubi.fileBase.strainRate, tp)) ;
    disp(['Loading strainrate results from disk: ' estrainFn])
    load(estrainFn, 'strainrate', 'tre_vtx', 'dev_vtx', 'theta_vtx')
    tre = reshape(tre_vtx(:, 1:end-1), [size(tre_vtx, 1) * (size(tre_vtx,2)-1), 1]) ; 
    dev = reshape(dev_vtx(:, 1:end-1), [size(tre_vtx, 1) * (size(tre_vtx,2)-1), 1]) ; 
    theta = reshape(theta_vtx(:, 1:end-1), [size(tre_vtx, 1) * (size(tre_vtx,2)-1), 1]) ; 
else
    tre = options.tre ;
    dev = options.dev ;
    theta = options.theta ;
end

%% Prepare directories for images
dirs2make = { egImDir, fullfile(egImDir, 'strainRate3d'), ...
    fullfile(egImDir, 'strainRate2d') } ;
for ii = 1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bbr256 = blueblackred(256) ;
% clf
% set(gcf, 'visible', 'off')
% imagesc([-1, 0, 1; -1, 0, 1])
% caxis([0, 1])
% pos256 = bluewhitered(256) ;
close all
pm256 = phasemap(256) ;

%% Plot the metric components on trisurf
% denom = sqrt(tg(:, 1, 1) .* tg(:, 2, 2)) ;
% NOTE: \varepsilon --> ${\boldmath${\varepsilon}$}$
labels = {'$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon] $', ...
    '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$'} ;
time_in_units = (tp - t0) * tubi.timeInterval ;
tstr = [': $t=$', sprintf('%03d', time_in_units), ' ', tubi.timeUnits ];

%% consider each metric element & plot in 3d
fn = fullfile(egImDir, 'strainRate3d', sprintf([tubi.fileBase.spcutMeshSmRSC '.png'], tp));
if ~exist(fn, 'file') || overwrite
    clf
    set(gcf, 'visible', 'off') ;
    for qq = 1:2
        % For each view (dorsal, ventral, left, right)
        % for pp = 1:4
        subplot(1, 2, qq) ;
        if qq == 1
            % trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            %     'FaceVertexCData', 0.5* tre, 'edgecolor', 'none')
            trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                0.5* tre(:), 'edgecolor', 'none')
            caxis([-clim_trace, clim_trace])
            colormap(gca, bbr256)
            colorbar('location', 'southOutside') ;      
            % ylabel(cb, labels{qq}, 'Interpreter', 'Latex')

        else
            % Intensity from dev and color from the theta
            indx = max(1, round(mod(2*theta, 2*pi)*size(pm256, 1)/(2 * pi))) ;
            colors = pm256(indx, :) ;
            colors = min(dev(:) / clim_deviatoric, 1) .* colors ;
            trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                'FaceVertexCData', colors, 'edgecolor', 'none')

            % Colorbar and phasewheel
            colormap(gca, phasemap) ;
            phasebar('colormap', phasemap, ...
                'location', [0.82, 0.1, 0.1, 0.135], 'style', 'nematic') ;
            ax = gca ;
            get(gca, 'position') ;
            cb = colorbar('location', 'southOutside') ;
            drawnow
            axpos = get(ax, 'position') ;
            cbpos = get(cb, 'position') ;
            set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3)*0.6, cbpos(4)])
            set(ax, 'position', axpos) ;
            hold on;
            caxis([0, clim_deviatoric])
            colormap(gca, gray)
        end

        axis equal
        xlim(xyzlim(1, :))
        ylim(xyzlim(2, :))
        zlim(xyzlim(3, :))
        axis off
        title(labels{qq}, 'Interpreter', 'Latex')      

        % Save images;
        % dorsal
        % if pp == 1
        %     view(0, 90)
        % ventral
        % elseif pp == 2
        %     view(0, -90)
        % left
        % elseif pp == 3
        %     view(0, 0)
        % right
        % elseif pp == 4
        %     view(0, 180)
        % end
        % end

        % left view
        view(0, 0) 
    end
    sgtitle(['strain rate, ', tstr], 'Interpreter', 'latex') 

    % Save the image
    disp(['Saving 3d strain image: ' fn])
    saveas(gcf, fn) ;
    clf

end

%% Now plot in 2d
close all
set(gcf, 'visible', 'off') ;
fn = fullfile(egImDir, 'strainRate2d', ...
        sprintf([tubi.fileBase.spcutMeshSm '.png'], tp));
if ~exist(fn, 'file') || overwrite
    % Panel 1
    subplot(1, 2, 1) ;
    % trisurf(cutMesh.f, ...
    %     cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
    %     cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
    %     'FaceVertexCData', 0.5* tre, 'edgecolor', 'none')
    trisurf(cutMesh.f, ...
        cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
        cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
        0.5* tre_vtx(:), 'edgecolor', 'none')
    daspect([1,1,1])
    cb = colorbar('location', 'southOutside') ;

    caxis([-clim_trace, clim_trace])
    title(labels{1}, 'Interpreter', 'Latex')   
    colormap(bbr256)
    axis off
    view(2)

    % Panel 2 
    subplot(1, 2, 2) ;
    % Intensity from dev and color from the theta
    indx = max(1, round(mod(2*theta_vtx(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
    colors = pm256(indx, :) ;
    colors = min(dev_vtx(:) / clim_deviatoric, 1) .* colors ;
    trisurf(cutMesh.f, cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
        cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
        'FaceVertexCData', colors, 'edgecolor', 'none')
    daspect([1,1,1]) ;
    title(labels{2}, 'Interpreter', 'Latex')   

    % Colorbar and phasewheel
    colormap(gca, phasemap)
    phasebar('colormap', phasemap, ...
        'location', [0.82, 0.12, 0.1, 0.135], 'style', 'nematic') ;
    axis off
    view(2)
    ax = gca ;
    get(gca, 'position') ;
    cb = colorbar('location', 'southOutside') ;
    drawnow
    axpos = get(ax, 'position') ;
    cbpos = get(cb, 'position') ;
    set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3)*0.6, cbpos(4)])
    set(ax, 'position', axpos) 
    hold on;
    caxis([0, clim_deviatoric])
    colormap(gca, gray)

    % Save the image
    sgtitle(['strain rate, ', tstr], 'Interpreter', 'latex') 
    saveas(gcf, fn) ;
    clf
end    
close all


%% Compare trace to trace of gdot determined via kinematics
close all
set(gcf, 'visible', 'off') ;
fn = fullfile(egImDir, 'strainRate2d', ...
        sprintf(['compare_' tubi.fileBase.spcutMeshSm '.png'], tp));
if (~exist(fn, 'file') || overwrite) && plot_comparison
    % Load gdot trace from kinematics
    % l_lmesh_lerr = strrep(sprintf('lambda%0.03f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d', ...
    %     lambda, lambda_mesh, nmodes, zwidth), '.', 'p') ;
    % fn_gdot = fullfile(tubi.dir.metricKinematics.root, l_lmesh_lerr, ...
    %    'measurements', sprintf('gdot_vertices_%06d.mat', tp)) ;
    % load(fn_gdot, 'gdot')
    
    % This doesn't work if there are dots on the path
    fn_gdot = sprintf(tubi.fullFileBase.metricKinematics.gdot, tp) ;
    load([fn_gdot, '.mat'], 'gdot')

    % Panel 1
    subplot(1, 2, 1)
    trisurf(cutMesh.f, ...
        cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
        cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
        0.5*tre_vtx(:), 'edgecolor', 'none')
    daspect([1,1,1])
    colorbar('location', 'southOutside') ;
    caxis([-clim_trace, clim_trace])
    title(labels{1}, 'Interpreter', 'Latex')   
    colormap(bbr256)
    axis off
    view(2)

    % Panel 2 
    subplot(1, 2, 2)
    % Comparison 1/2 * Tr[g^{-1}gdot]
    trisurf(cutMesh.f, cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
        cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
        0.5 * gdot, 'edgecolor', 'none')
    daspect([1,1,1])
    colorbar('location', 'southOutside') ;
    caxis([-clim_trace, clim_trace])
    title('$\frac{1}{2}$Tr$[g^{-1}\dot{g}]$', 'Interpreter', 'Latex')   
    colormap(bbr256)
    axis off
    view(2)

    % Save the image
    sgtitle(['comparison $\frac{1}{2}$Tr$[g^{-1}\dot{g}]$, ', tstr], 'Interpreter', 'latex') 
    
    disp(['Saving 2d strain image: ' fn])
    saveas(gcf, fn) ;
    clf
end    
close all


%% Compare to freshly computed
close all
set(gcf, 'visible', 'off') ;
fn = fullfile(egImDir, 'strainRate2d', ...
        sprintf(['compare_fresh_' tubi.fileBase.spcutMeshSm '.png'], tp));
if (~exist(fn, 'file') || overwrite) && debug
    DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
    H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;

    %% OPTION1: Could compute divergence anew on smoothed mesh
    % fvel = squeeze(face_vels(tidx, :, :)) ;
    % divv3d_test = DEC.divergence(fvel) ;

    %% OPTION2: Load divv from disk and convert using smoothed mesh
    % [~, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
    if strcmp(averagingStyle, 'Lagrangian')
        if doubleResolution
            dec_tp = load(sprintf(tubi.fullFileBase.decAvg2x, tp)) ;
        else
            dec_tp = load(sprintf(tubi.fullFileBase.decAvg, tp)) ;
        end
    else
        if doubleResolution
            dec_tp = load(sprintf(tubi.fullFileBase.decSimAvg2x, tp)) ;
        else
            dec_tp = load(sprintf(tubi.fullFileBase.decSimAvg, tp)) ;
        end
    end

    % Smooth divergence(v) [divv]
    divv3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
                            lambda, 'implicit', dec_tp.divs.raw) ;

    % veln = divv / (2H) ; 
    % veln_pred = DEC.divergence(fvel) ./ (2.0 * H) ;
    % veln is normal velocity field (velocities along vertex normals)
    % mesh.vn are vertex normals

    % Smooth the velocities in space using gptoolbox
    vx = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 1)) ;
    vy = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 2)) ;
    vz = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 3)) ;
    vxs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
        lambda, 'implicit', vx') ;
    vys = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
        lambda, 'implicit', vy') ;
    vzs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
        lambda, 'implicit', vz') ;

    % Actual normal velocity -- currently using ORIGINAL mesh normals
    veln3d = sum(cat(2, vxs, vys, vzs) .* mesh.vn, 2) ;

    % Predict the divergence
    % H2vn3d = 2 * H3d .* veln3d ;

    % Extend to have another row for 2d map
    divv2d = divv3d ;
    divv2d(nU*(nV-1) + 1:nU*nV) = divv3d(1:nU) ;
    H2d = H3d ;
    H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;
    veln2d = veln3d ;
    veln2d(nU*(nV-1)+1:(nU*nV)) = veln3d(1:nU) ;
    H2vn2d = 2 * H2d .* veln2d ;

    % The difference
    gdot2d = divv2d - H2vn2d ;
    % gdot3d = divv3d - H2vn3d ;

    % Panel 1
    subplot(1, 2, 1) ;
    trisurf(cutMesh.f, ...
        cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
        cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
        'FaceVertexCData', tre, 'edgecolor', 'none')
    daspect([1,1,1])
    cb = colorbar('location', 'southOutside') ;
    caxis([-clim_trace, clim_trace])
    title('Tr$[g^{-1}\dot{g}]$', 'Interpreter', 'Latex')   
    colormap(bbr256)
    axis off
    view(2)

    % Panel 2 
    subplot(1, 2, 2) ;
    % Comparison 1/2 * Tr[g^{-1}gdot]
    trisurf(cutMesh.f, cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
        cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
        gdot2d(:), 'edgecolor', 'none')
    daspect([1,1,1])
    colorbar('location', 'southOutside') ;
    caxis([-clim_trace, clim_trace])
    title('$\frac{1}{2}$Tr$[g^{-1}\dot{g}]$', 'Interpreter', 'Latex')   
    colormap(bbr256)
    axis off
    view(2)

    % Save the image
    sgtitle(['comparison fresh $\frac{1}{2}$Tr$[g^{-1}\dot{g}]$, ', tstr], 'Interpreter', 'latex') 
    saveas(gcf, fn) ;
    clf
end