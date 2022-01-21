function plotStrainRate(QS, options)
%plotStrainRateTimePoint(QS, tp, options)
%   Plot the traceful and traceless components of the strain rate tensor
%   defined on each face, both for individual timepoints, and 
%   also over time as kymographs
%
% Parameters
% ----------
% QS : QuapSlap class instance
% tp : int 
%   timepoint in units of (1/QS.timeInterval) * QS.timeUnits
% options: struct with fields
%   
% Returns
% -------
% 
% 
% NPMitchell 2020

% Default parameters
overwrite = false ;
clim_trace = 0.05 ;
clim_deviatoric = 0.05 ;
averagingStyle = 'Lagrangian' ;
skipTimePoint = false ;
lambda = QS.smoothing.lambda ;
lambda_mesh = QS.smoothing.lambda_mesh ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;

%% Unpack required params
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
debug = false ;

%% Parameters
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
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
if isfield(options, 'skipTimePoint')
    skipTimePoint = options.skipTimePoint ;
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
srlambdaDir = fullfile(QS.dir.strainRate.root, ...
    [strrep(sprintf('lambda%0.3f_lmesh%0.3f', ...
    lambda, lambda_mesh), '.', 'p'), ...
    sprintf('_modes%02dw%02d', nmodes, zwidth)]) ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
nU = QS.nU ;
nV = QS.nV ;
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Prepare directories for images
dirs2make = { srlambdaDir, ...
    fullfile(srlambdaDir, 'strainRate3d'), ...
    fullfile(srlambdaDir, 'strainRate2d') } ;
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
bwr256 = bluewhitered(256) ;
bbr256 = blueblackred(256) ;
% clf
% set(gcf, 'visible', 'off')
% imagesc([-1, 0, 1; -1, 0, 1])
% caxis([0, 1])
% pos256 = bluewhitered(256) ;
close all
pm256 = phasemap(256) ;

%% Collate data from each quarter of the gut
if strcmp(averagingStyle, 'simple')
    sKDir = fullfile(QS.dir.strainRateSimple, ...
        [strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f'], ...
        lambda, lambda_mesh), '.', 'p'), ...
        sprintf('_modes%02dw%02d', nmodes, zwidth)]);
else
    sKDir = fullfile(QS.dir.strainRate.root, ...
        [strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f'], ...
        lambda, lambda_mesh), '.', 'p'), ...
        sprintf('_modes%02dw%02d', nmodes, zwidth)]);
end

%% Load or collate kymograph data
datdir = fullfile(sKDir, 'measurements') ;
apKymoFn = fullfile(datdir, 'apKymographStrainRate.mat') ;
lKymoFn = fullfile(datdir, 'leftKymographStrainRate.mat') ;
rKymoFn = fullfile(datdir, 'rightKymographStrainRate.mat') ;
dKymoFn = fullfile(datdir, 'dorsalKymographStrainRate.mat') ;
vKymoFn = fullfile(datdir, 'ventralKymographStrainRate.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if files_exist
    load(apKymoFn, 'tr_apM', 'dv_apM', 'th_apM')
    load(lKymoFn, 'tr_lM', 'dv_lM', 'th_lM')
    load(rKymoFn, 'tr_rM', 'dv_rM', 'th_rM')
    load(dKymoFn, 'tr_dM', 'dv_dM', 'th_dM')
    load(vKymoFn, 'tr_vM', 'dv_vM', 'th_vM')
else
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        tidx = QS.xp.tIdx(tp) ;

        %% Define metric strain filename        
        estrainFn = fullfile(srlambdaDir, 'measurements', ...
            sprintf(QS.fileBase.strainRate, tp)) ;
        disp(['t=' num2str(tp) ': Loading strainrate results from disk: ' estrainFn])
        load(estrainFn, 'strainrate', 'tre', 'dev', 'theta', ...
            'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
            'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v', ...
            'lambda', 'lambda_mesh')

        %% Plot strain rate for this timepoint
        if ~skipTimePoint
            tpOpts.lambda = lambda ;
            tpOpts.lambda_mesh = lambda_mesh ;
            tpOpts.overwrite = overwrite ;
            QS.plotStrainRateTimePoint(tp, tpOpts)
        end

        %% Store in matrices    
        % dv averaged
        tr_apM(tidx, :) = tre_ap ;
        dv_apM(tidx, :) = dev_ap ;
        th_apM(tidx, :) = theta_ap ;

        % left quarter
        tr_lM(tidx, :) = tre_l ;
        dv_lM(tidx, :) = dev_l ;
        th_lM(tidx, :) = theta_l ;

        % right quarter
        tr_rM(tidx, :) = tre_r ;
        dv_rM(tidx, :) = dev_r ;
        th_rM(tidx, :) = theta_r ;

        % dorsal quarter
        tr_dM(tidx, :) = tre_d ;
        dv_dM(tidx, :) = dev_d ;
        th_dM(tidx, :) = theta_d ;

        % ventral quarter
        tr_vM(tidx, :) = tre_v ;
        dv_vM(tidx, :) = dev_v ;
        th_vM(tidx, :) = theta_v ;

    end
    
    %% Save kymographs
    save(apKymoFn, 'tr_apM', 'dv_apM', 'th_apM')
    save(lKymoFn, 'tr_lM', 'dv_lM', 'th_lM')
    save(rKymoFn, 'tr_rM', 'dv_rM', 'th_rM')
    save(dKymoFn, 'tr_dM', 'dv_dM', 'th_dM')
    save(vKymoFn, 'tr_vM', 'dv_vM', 'th_vM')
end

%% Store kymograph data in cell arrays
tresK = {tr_apM, tr_lM, tr_rM, tr_dM, tr_vM} ;
devsK = {dv_apM, dv_lM, dv_rM, dv_dM, dv_vM} ;
thetasK = {th_apM, th_lM, th_rM, th_dM, th_vM} ;

%% Now plot different measured quantities as kymographs
% Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
dvDir = fullfile(sKDir, 'avgDV') ;
lDir = fullfile(sKDir, 'avgLeft') ;
rDir = fullfile(sKDir, 'avgRight') ;
dDir = fullfile(sKDir, 'avgDorsal') ;
vDir = fullfile(sKDir, 'avgVentral') ;
outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
titleadd = {': circumferentially averaged', ...
    ': left side', ': right side', ': dorsal side', ': ventral side'} ;
tps = QS.xp.fileMeta.timePoints(1:end-1) - t0 ;
titles = {'dilation, $\frac{1}{2}\textrm{Tr}[g^{-1}\dot{\varepsilon}]$ ',...
          'shear, $||\varepsilon-\frac{1}{2}\mathrm{Tr}\left[\mathbf{g}^{-1}\dot{\varepsilon}\right)\bf{g}||$'} ;  
for qq = 1:length(outdirs)
    %% Prep the output directory for this averaging
    odir = outdirs{qq} ;
    if ~exist(odir, 'dir')
        mkdir(odir)
    end

    %% Plot Trace kymograph  
    % Check if images already exist on disk    
    label = '$\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]$' ;
    name = 'tre' ;
    fn = fullfile(odir, [ name '.png']) ;
    fn_zoom = fullfile(odir, [name '_zoom_early.png']) ;

    close all; set(gcf, 'visible', 'off')
    clim_zoom = clim_trace * 0.5 ;
    if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
        trK = tresK{qq} ;
        imagesc((1:nU)/nU, tps, trK)
        caxis([-clim_trace, clim_trace])
        colormap(bbr256)
        
        % Titles 
        title([titles{1}, titleadd{qq}], 'Interpreter', 'Latex')
        ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
        xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')

        % Plot fold identifications
        hold on;
        fons1 = max(1, fons(1)) ;
        fons2 = max(1, fons(2)) ;
        fons3 = max(1, fons(3)) ;
        plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
        plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
        plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

        tidx0 = QS.xp.tIdx(t0) ;
        cb = colorbar() ;
        ylabel(cb, label, 'Interpreter', 'Latex')  
        
        % Save 
        fn = fullfile(odir, [name '.png']) ;
        export_fig(fn, '-png', '-nocrop', '-r200')   

        % Zoom in on small values
        caxis([-clim_zoom, clim_zoom])
        colormap(bbr256)
        fn = fullfile(odir, [name '_zoom.png']) ;
        disp(['saving ', fn])
        export_fig(fn, '-png', '-nocrop', '-r200')   
        % Zoom in on early times
        ylim([min(tps), max(fons) + 10])
        caxis([-clim_zoom, clim_zoom])
        colormap(bbr256)
        fn = fullfile(odir, [name '_zoom_early.png']) ;
        disp(['saving ', fn])
        export_fig(fn, '-png', '-nocrop', '-r200')   
    end
    

    %% Plot DV-averaged/quarter-avgeraged kymograph -- Deviator   
    % denom = sqrt(tg(:, 1, 1) .* tg(:, 2, 2)) ;
    % NOTE: \varepsilon --> ${\boldmath${\varepsilon}$}$
    label = '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$' ;
    name = 'dev' ;
    
    % Check if images already exist on disk
    fn = fullfile(odir, [ name '.png']) ;
    fn_zoom = fullfile(odir, [name '_zoom_early.png']) ;
    zoomstrs = {'', '_zoom'} ;
    
    if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
        for pp = 1:2
            zoomstr = zoomstrs{pp} ;        
            if pp == 1
                clim = clim_deviatoric ;
            else
                clim = clim_deviatoric * 0.5 ;
            end
            
            % Unpack what to plot (averaged kymographs, vary averaging region)
            devK = devsK{qq} ;
            thetaK = thetasK{qq} ;

            % Map intensity from dev and color from the theta
            indx = max(1, round(mod(2*thetaK(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
            colors = pm256(indx, :) ;
            devKclipped = min(devK / clim, 1) ;
            colorsM = devKclipped(:) .* colors ;
            colorsM = reshape(colorsM, [size(devK, 1), size(devK, 2), 3]) ;

            % Plot the kymograph
            close all
            set(gcf, 'visible', 'off')
            imagesc((1:nU)/nU, tps, colorsM)
            caxis([-clim, clim])
            % Add folds to plot
            hold on;
            fons1 = max(1, fons(1)) ;
            fons2 = max(1, fons(2)) ;
            fons3 = max(1, fons(3)) ;
            plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
            plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
            plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

            % Colorbar and phasewheel
            colormap(gca, phasemap)
            phasebar('colormap', phasemap, ...
                'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic') ;
            ax = gca ;
            cb = colorbar('location', 'eastOutside') ;
            drawnow
            axpos = get(ax, 'position') ;
            cbpos = get(cb, 'position') ;
            set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3), cbpos(4)*0.6])
            set(ax, 'position', axpos) 
            hold on;
            caxis([0, clim])
            colormap(gca, gray)

            % title and save
            title([titles{2}, titleadd{qq}], 'Interpreter', 'Latex')
            ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
            ylabel(cb, label, 'Interpreter', 'Latex')  
            fn = fullfile(odir, [ name zoomstr '.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   

            if pp == 2
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                fn = fullfile(odir, [name zoomstr '_early.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end
disp('done with plotting strain rate')


