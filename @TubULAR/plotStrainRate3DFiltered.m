function plotStrainRate3DFiltered(QS, options)
%plotStrainRate3DFiltered(QS, options)
%   load spatially-smoothed 
%   epsilon = 1/2 (\nabla_i v_j + \nabla_j v_i) - vN b_{ij} from disk,
%   smooth with time filter and plot in 3d
%   
% Parameters
% ----------
% QS : QuapSlap class object instance
% options : struct with fields
%   lambda : float
%   lambda_mesh : float
%   overwrite : bool 
%   preview : bool 
%   averagingStyle : str specifier ('Lagrangian' or 'Simple')
%   samplingResolution : '1x' or '2x'
%   debug : bool 
%   nTimePoints
% 
% Returns 
% -------
%
% Saves to disk
% -------------
%
%
% NPMitchell 2020

%% Default options
lambda = 0.01 ;
lambda_mesh = 0.0 ;
overwrite = false ;
preview = true ;
averagingStyle = 'Lagrangian' ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
debug = false ;
nTimePoints = 7 ;
clim_trace = 0.05 ;
clim_deviatoric = 0.05 ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'nmodes')
    nmodes = options.nmodes ;
end
if isfield(options, 'zwidth')
    zwidth = options.zwidth ;
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
if isfield(options, 'nTimePoints')
    nTimePoints = options.nTimePoints ;
end
if isfield(options, 'clim_trace')
    clim_trace = options.clim_trace ;
end
if isfield(options, 'clim_deviatoric')
    clim_deviatoric = options.clim_deviatoric ;
end

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
t0 = QS.t0set() ;
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
% Output directory
egImDir = strrep(sprintf( ...
    QS.dir.strainRate.smoothing, lambda, lambda_mesh, nmodes, zwidth), '.', 'p') ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;

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

% Pre-assign timepoints with velocities
tpts = QS.xp.fileMeta.timePoints(1:end-1) ;
tp2do = [tpts(1:20:end), setdiff(tpts, tpts(1:50:end))] ;

% Build metric from mesh
for tp = tp2do
    disp(['t = ' num2str(tp)])

    tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
    mesh = tmp.spcutMeshSmRSC ;
    tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
    cutMesh = tmp.spcutMeshSmRS ;

    % DEBUG
    % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
    % mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) * mesh.ar ;
    mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
    cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
    
    glueMesh = glueCylinderCutMeshSeam(cutMesh) ;

    % Load current mesh +/- nTimePoints
    first = true ;
    tp2doFilter = tp-nTimePoints:tp+nTimePoints ;
    for tpsFilter = tp2doFilter
        tpj = min(max(tpsFilter, QS.xp.fileMeta.timePoints(1)), ...
            QS.xp.fileMeta.timePoints(end-1)) ;
        tidx = QS.xp.tIdx(tpj) ;
        
        %% load the metric strain
        % Define metric strain filename        
        estrainFn = fullfile(strrep(sprintf(QS.dir.strainRate.measurements, ...
            lambda, lambda_mesh), '.', 'p'), ...
            sprintf(QS.fileBase.strainRate, tpj)) ;
        disp(['Loading strainrate results from disk: ' estrainFn])
        % save(estrainFn, 'strainrate', 'tre', 'dev', 'theta', 'theta_pb', ...
        %     'dvij', 'gg', 'bb', 'lambda', 'lambda_mesh', 'readme', ...
        %     'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
        %     'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v', ...
        %     'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
        %     'strainrate_vtx', 'tre_vtx', 'dev_vtx', 'theta_vtx', ...
        %     'gg_vtx', 'bb_vtx', 'dx_faces', 'dy_faces', 'dx_vtx', 'dy_vtx')
        % load(estrainFn, 'tre_vtx', 'dev_vtx', 'theta_vtx', 'dx_vtx', 'dy_vtx')
        load(estrainFn, 'tre', 'dev', 'theta')
        
        if first
            tres = tre ;
            devs = dev ;
            thetas = theta ;
            first = false ;
        else
            tres = cat(2, tres, tre) ;
            devs = cat(2, devs, dev) ;
            thetas = cat(2, thetas, theta) ;
        end
        
    end
    
    % Build TIME-AVERAGED tre, dev, theta
    % average over time
    tre = mean(tres, 2) ;
    % Note that this is not actually dv averaging: averaging over time
    [dev, theta] = QS.dvAverageNematic(devs, thetas) ;
    
    % Transfer to vertices
    [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.u) ;
    tre = reshape(F2V * tre, [nU, nV-1]) ;
    dev = reshape(F2V * dev, [nU, nV-1]) ;
    theta = reshape(F2V * theta, [nU, nV-1]) ;
    
    % Space-filter the data using spectral filter
    filtOpts = struct('nmodes', nmodes, 'xwidth', zwidth) ;
    tre = modeFilterQuasi1D(tre, filtOpts) ;
    ctheta = modeFilterQuasi1D(dev .* cos(2 * theta), filtOpts) ;
    stheta = modeFilterQuasi1D(dev .* sin(2 * theta), filtOpts) ;
    dev = sqrt(ctheta.^2 + stheta.^2) ;
    theta = 0.5 * mod(atan2(stheta, ctheta), 2*pi) ;
    
    % Restrict to single cover and unravel for vertices
    % tre = tre_avg(1:nU*(nV-1)) ;
    % theta = theta_avg(1:nU*(nV-1)) ;
    % dev = dev_avg(1:nU*(nV-1)) ;
    
    % denom = length(tp2doFilter) ;
    % vmean = vmean / denom ;
    
    %% Prepare directories for images
    dirs2make = { egImDir, fullfile(egImDir, 'strainRate3d_filtered'), ...
        fullfile(egImDir, 'strainRate2d_filtered') } ;
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
    time_in_units = (tp - t0) * QS.timeInterval ;
    tstr = [': $t=$', sprintf('%03d', time_in_units), ' ', QS.timeUnits ];

    %% consider each metric element & plot in 3d
    fn = fullfile(egImDir, 'strainRate3d_filtered', sprintf([QS.fileBase.spcutMeshSmRSC '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        clf
        set(gcf, 'visible', 'off') ;
        for qq = 1:2
            % For each view (dorsal, ventral, left, right)
            % for pp = 1:4
            subplot(1, 2, qq) ;
            if qq == 1
                trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                    'FaceVertexCData', 0.5* tre(:), 'edgecolor', 'none')
                % trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                %     0.5* tre, 'edgecolor', 'none')
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
        disp(['Saving figure: ' fn])
        saveas(gcf, fn) ;
        clf

    end

    %% Now plot in 2d
    close all
    set(gcf, 'visible', 'off') ;
    fn = fullfile(egImDir, 'strainRate2d_filtered', ...
            sprintf([QS.fileBase.spcutMeshSm '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        % Panel 1
        subplot(1, 2, 1) ;
        tre2d = tre ;
        tre2d(:, nV) = tre2d(:, 1) ;
        trisurf(cutMesh.f, ...
            cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
            cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
            'FaceVertexCData', 0.5* tre2d(:), 'edgecolor', 'none')
        
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
        theta2d = theta ;
        theta2d(:, nV) = theta(:, 1) ;
        dev2d = dev ;
        dev2d(:, nV) = dev2d(:, 1) ;
        indx = max(1, round(mod(2*theta2d, 2*pi)*size(pm256, 1)/(2 * pi))) ;
        colors = pm256(indx, :) ;
        colors = min(dev2d(:) / clim_deviatoric, 1) .* colors ;
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
        disp(['Saving figure: ' fn])
        saveas(gcf, fn) ;
        clf
    end    
    close all

end

