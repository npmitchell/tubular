function plotPathlineStrainRateTimePoint(QS, tp, options)
%plotPathlineStrainRateTimePoint(QS, tp, options)
%   Plot the traceful and traceless components of the strain rate tensor
%   defined at pathlines moving in the domain of parameterization over
%   time.
%   Here, plot the trajectories as immobile in a grid. Uses vertex 
%   pathlines. 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% tp : int 
%   timepoint in units of (1/QS.timeInterval) * QS.timeUnits
% options: struct with fields
%   
% 
% NPMitchell 2020

tidx = QS.xp.tIdx(tp) ;

%% Unpack required params
lambda = options.lambda ;
lambda_mesh = options.lambda_mesh ;
t0Pathline = options.t0Pathline ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
debug = false ;

%% Parameters
overwrite = false ;
clim_trace = 0.05 ;
clim_deviatoric = 0.05 ;
clim_strain = 0.5 ;
averagingStyle = 'Lagrangian' ;
plot_comparison = false ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'plot_comparison')
    plot_comparison = options.plot_comparison ;
end
if isfield(options, 'mesh')
    mesh = options.mesh ;
else
    mesh = [] ;
end
if isfield(options, 'cutMesh')
    cutMesh = options.cutMesh ;
else
    cutMesh = [] ;
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

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
    sresStr = 'doubleRes_' ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
t0 = QS.t0set() ;
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
% Output directory
egImDir = sprintf( QS.dir.strainRate.pathline.root, t0Pathline) ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

% Check if ALL plot files already exist
% FLOWS
fn_strainRate2d = fullfile(egImDir, 'strainRate2d', ...
    sprintf([QS.fileBase.spcutMeshSm '.png'], tp));
fn_strain2d = fullfile(egImDir, 'strain2d', ...
    sprintf([QS.fileBase.spcutMeshSm '.png'], tp));
fn_strainRateCompare = fullfile(egImDir, 'strainRate2d_compareDEC', ...
    sprintf(['compare_' QS.fileBase.spcutMeshSm '.png'], tp));
redo_images = ~exist(fn_strainRate2d, 'file') || ...
              ~exist(fn_strain2d, 'file') || ...
              (~exist(fn_strainRateCompare, 'file') && ...
                plot_comparison) || ...
              overwrite ;
          
% Create the figures if they are not on disk
if redo_images
    %% load the metric strain
    % Define metric strain filename        
    if ~isfield(options, 'tre') || ~isfield(options, 'dev') || ...
            ~isfield(options, 'theta') || ...
            ~isfield(options, 'strain_tr') || ~isfield(options, 'strain_dv') || ...
            ~isfield(options, 'strain_th')
        estrainFn = fullfile(sprintf(QS.dir.strainRate.pathline.measurements, ...
            t0Pathline), sprintf(QS.fileBase.strainRate, tp)) ;
        disp(['Loading strainrate results from disk: ' estrainFn])
        load(estrainFn, 'strainrate', 'tre', 'dev', 'theta')
    else
        tre = options.tre ;
        dev = options.dev ;
        theta = options.theta ;
        % strain_tr = options.strain_tr ;
        % strain_dv = options.strain_dv ;
        % strain_th = options.strain_th ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Timepoint specific operations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    titlestr = ['$t=$' sprintf('%03d', tp - t0) ' ' QS.timeUnits ] ;
    % Load meshes if not supplied    
    if isempty(mesh) && redo_images
        % If they don't all already exist in RAM, load the meshes
        % Load current mesh
        if doubleResolution
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC2x, tp)) ;
            mesh = tmp.spcutMeshSmRSC2x ;
        else
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
            mesh = tmp.spcutMeshSmRSC ;
        end

        % Smooth the mesh with lambda_mesh
        if lambda_mesh > 0 
            disp('smoothing mesh vertices before computations')
            tri = triangulation(mesh.f, mesh.v) ;
            fbndy = tri.freeBoundary ;
            fbndy = fbndy(:, 1) ;
            mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda_mesh, 'implicit', mesh.v) ;
        end

        % Load cutMesh also
        if isempty(cutMesh)
            % Load cutMesh
            if doubleResolution
                tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS2x, tp)) ;
                cutMesh = tmp.spcutMeshSmRS2x ;
                clearvars tmp
            else
                tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
                cutMesh = tmp.spcutMeshSmRS ;
                clearvars tmp
            end
        end
    end


    %% Prepare directories for images
    dirs2make = { fullfile(egImDir, 'strain2d'), ...
        fullfile(egImDir, 'strainRate2d'), ...
        fullfile(egImDir, 'strainRate2d_compareDEC'), ...
        fullfile(egImDir, 'strainRate2d_compareDECnew')} ;
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
    % bwr256 = bluewhitered(256) ;
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
    time_in_units = (tp - t0) * QS.timeInterval ;
    tstr = [': $t=$', sprintf('%03d', time_in_units), QS.timeUnits ];

    %% Now plot STRAIN RATE in 2d
    close all
    set(gcf, 'visible', 'off') ;
    if ~exist(fn_strainRate2d, 'file') || overwrite
        % Panel 1
        subplot(1, 2, 1)
        trisurf(cutMesh.f, ...
            cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
            'FaceVertexCData', 0.5 * tre(:), 'edgecolor', 'none')
        daspect([1,1,1])
        cb = colorbar('location', 'southOutside') ;
        caxis([-clim_trace, clim_trace])
        title(labels{1}, 'Interpreter', 'Latex')   
        colormap(bbr256)
        axis off; view(2)

        % Panel 2 
        subplot(1, 2, 2)
        % Intensity from dev and color from the theta
        indx = max(1, round(mod(2*theta(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
        colors = pm256(indx, :) ;
        colors = min(dev(:) / clim_deviatoric, 1) .* colors ;
        trisurf(cutMesh.f, cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
            'FaceVertexCData', colors, 'edgecolor', 'none')
        daspect([1,1,1])
        title(labels{2}, 'Interpreter', 'Latex')   

        % Colorbar and phasewheel
        colormap(gca, phasemap)
        phasebar('colormap', phasemap, ...
            'location', [0.82, 0.12, 0.1, 0.135], 'style', 'nematic') ;
        axis off
        view(2)
        ax = gca ;
        cb = colorbar('location', 'southOutside') ;
        drawnow
        axpos = get(ax, 'position') ;
        cbpos = get(cb, 'position') ;
        set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3)*0.6, cbpos(4)])
        set(ax, 'position', axpos) ;
        hold on;
        caxis([0, clim_deviatoric])
        colormap(gca, gray)

        % Save the image
        sgtitle(['strain rate, ', tstr], 'Interpreter', 'latex') 
        disp(['saving fig: ' fn_strainRate2d])
        saveas(gcf, fn_strainRate2d) ;
        clf
    end    
    close all

% %% Now plot STRAIN in 2d
% close all
% set(gcf, 'visible', 'off') ;
% if ~exist(fn_strain2d, 'file') || overwrite
%     % Panel 1 --> trace
%     subplot(1, 2, 1)
%     trisurf(cutMesh.f, ...
%         cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
%         cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
%         'FaceVertexCData', 0.5 * strain_tr(:), 'edgecolor', 'none')
%     daspect([1,1,1])
%     cb = colorbar('location', 'southOutside') ;
% 
%     caxis([-clim_strain, clim_strain])
%     title(labels{1}, 'Interpreter', 'Latex')   
%     colormap(bbr256)
%     axis off
%     view(2)
% 
%     % Panel 2 --> deviator 
%     subplot(1, 2, 2)
%     % Intensity from dev and color from the theta
%     indx = max(1, round(mod(2*strain_th(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
%     colors = pm256(indx, :) ;
%     colors = min(strain_dv(:) / clim_strain, 1) .* colors ;
%     trisurf(cutMesh.f, cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
%         cutMesh.u(:, 2), 0*cutMesh.u(:, 1), ...
%         'FaceVertexCData', colors, 'edgecolor', 'none')
%     daspect([1,1,1])
%     title(labels{2}, 'Interpreter', 'Latex')   
% 
%     % Colorbar and phasewheel
%     colormap(gca, phasemap)
%     phasebar('colormap', phasemap, ...
%         'location', [0.82, 0.12, 0.1, 0.135], 'style', 'nematic') ;
%     axis off
%     view(2)
%     ax = gca ;
%     cb = colorbar('location', 'southOutside') ;
%     drawnow
%     axpos = get(ax, 'position') ;
%     cbpos = get(cb, 'position') ;
%     set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3)*0.6, cbpos(4)])
%     set(ax, 'position', axpos) ;
%     hold on;
%     caxis([0, clim_strain])
%     colormap(gca, gray)
% 
%     % Save the image
%     sgtitle(['strain rate, ', tstr], 'Interpreter', 'latex') 
%     saveas(gcf, fn_strain2d) ;
%     clf
% end    
% close all

    close all
end
