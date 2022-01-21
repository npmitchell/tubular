function plotPathlineStrain(QS, options)
% plotPathlineStrain(QS, options)
%   Plot the strain (from piv pathlines) along pathlines as kymographs and 
%   correlation plots. These are computed via finite differencing of 
%   pathline mesh metrics. That is, vertices are advected along piv flow 
%   and projected into 3D from pullback space to embedding, then the
%   metrics g_{ij} of those meshes are measured and differenced. 
% 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs         : bool
%   plot_kymographs_cumsum  : bool
%   plot_correlations       : bool 
%   plot_fold_strainRate    : bool
%   plot_lobe_strainRate    : bool
%   plot_fold_strain        : bool
%   plot_lobe_strain        : bool
% 
% Returns
% -------
% <none>
%
% NPMitchell 2020

error('not finished, since this is adapted from StrainRate. Mostly there, but do check things')

%% Default options 
overwrite = false ;
plot_kymographs = true ;
plot_kymographs_strain = true ;
plot_fold_strainRate = true ;
plot_lobe_strainRate = true ;
plot_fold_strain = true ;
plot_lobe_strain = true ;
maxWFrac = 0.03 ;  % Maximum halfwidth of feature/fold regions as fraction of ap length
t0 = QS.t0set() ;
t0Pathline = t0 ;
featureIDOptions = struct() ;

sRate_trace_label = '$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\dot{\varepsilon}]$';
sRate_deviator_label = ...
    '$||\dot{\varepsilon}-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\dot{\varepsilon}\right]\bf{g}||$' ;
% sRate_theta_label = '$\theta_{\mathrm{Dev}[\varepsilon]}$'; 
strain_trace_label = '$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon]$';
strain_trace_label_norm = '$||\mathrm{Tr} [\varepsilon]||$';
strain_deviator_label = ...
    '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$' ;
strain_deviator_label_short = ...
    '$||$Dev$\left[\varepsilon\right]||$' ;
% strain_theta_label = '$\theta_{\mathrm{Dev}[\dot{\varepsilon}]}$'; 

%% Parameter options
lambda = QS.smoothing.lambda ;
lambda_mesh = QS.smoothing.lambda_mesh ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;
climit = 0.2 ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
%% parameter options
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
end

% smoothing parameter options
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

% plotting and sampling
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'climitWide')
    climitWide = options.climitWide ;
else
    climitWide = climit * 3; 
end
if isfield(options, 'maxWFrac')
    maxWFrac = options.maxWFrac ;
end
if isfield(options, 'featureIDOptions')
    featureIDOptions = options.featureIDOptions ;
end

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_kymographs_strain')
    plot_kymographs_strain = options.plot_kymographs_strain ;
end
if isfield(options, 'plot_fold_strainRate')
    plot_fold_strainRate = options.plot_fold_strainRate ;
end
if isfield(options, 'plot_lobe_strainRate')
    plot_lobe_strainRate = options.plot_lobe_strainRate ;
end
if isfield(options, 'plot_fold_strain')
    plot_fold_strain = options.plot_fold_strain ;
end
if isfield(options, 'plot_lobe_strain')
    plot_lobe_strain = options.plot_lobe_strain ;
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
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

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
bluecolor = QS.plotting.colors(1, :) ;
orangecolor = QS.plotting.colors(2, :) ;
yellowcolor = QS.plotting.colors(3, :) ;
purplecolor = QS.plotting.colors(4, :) ;
greencolor = QS.plotting.colors(5, :) ;
graycolor = QS.plotting.colors(8, :) ;
browncolor = QS.plotting.colors(9, :) ;

%% Choose colors
trecolor = yellowcolor ;
Hposcolor = greencolor ;
Hnegcolor = purplecolor ;
Hsz = 3 ;  % size of scatter markers for mean curvature

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

%% Test incompressibility of the flow on the evolving surface
% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - t0 ;

% Output directory is inside metricKinematics dir
mKPDir = sprintf(QS.dir.strainRate.pathline.root, t0Pathline) ;
datdir = fullfile(mKPDir, 'measurements') ;
% Data for kinematics on meshes (defined on vertices) [not needed here]
% mdatdir = fullfile(mKDir, 'measurements') ;

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];

%% Compute or load all timepoints
apKymoFn = fullfile(datdir, 'apKymographPathlineStrainRate.mat') ;
lKymoFn = fullfile(datdir, 'leftKymographPathlineStrainRate.mat') ;
rKymoFn = fullfile(datdir, 'rightKymographPathlineStrainRate.mat') ;
dKymoFn = fullfile(datdir, 'dorsalKymographPathlineStrainRate.mat') ;
vKymoFn = fullfile(datdir, 'ventralKymographPathlineStrainRate.mat') ;
apSKymoFn = fullfile(datdir, 'apKymographPathlineStrain.mat') ;
lSKymoFn = fullfile(datdir, 'leftKymographPathlineStrain.mat') ;
rSKymoFn = fullfile(datdir, 'rightKymographPathlineStrain.mat') ;
dSKymoFn = fullfile(datdir, 'dorsalKymographPathlineStrain.mat') ;
vSKymoFn = fullfile(datdir, 'ventralKymographPathlineStrain.mat') ;
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
    % preallocate for kymos
    ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
    dv_apM = zeros(ntps, nU) ;   % dv averaged
    tr_apM = zeros(ntps, nU) ;   
    th_apM = zeros(ntps, nU) ;   
    dv_lM = zeros(ntps, nU) ;    % left averaged
    tr_lM = zeros(ntps, nU) ;
    th_lM = zeros(ntps, nU) ;
    dv_rM = zeros(ntps, nU) ;    % right averaged
    tr_rM = zeros(ntps, nU) ;
    th_rM = zeros(ntps, nU) ;
    dv_dM = zeros(ntps, nU) ;    % dorsal averaged
    tr_dM = zeros(ntps, nU) ;
    th_dM = zeros(ntps, nU) ;
    dv_vM = zeros(ntps, nU) ;    % ventral averaged
    tr_vM = zeros(ntps, nU) ;
    th_vM = zeros(ntps, nU) ;

    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        close all
        disp(['t = ' num2str(tp)])
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        srfn = fullfile(datdir, sprintf('strainRate_%06d.mat', tp))   ;

        % Load timeseries measurements
        load(srfn, 'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
            'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v')
   
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
        
        %% Save kymograph data
        save(apKymoFn, 'tr_apM', 'dv_apM', 'th_apM')
        save(lKymoFn, 'tr_lM', 'dv_lM', 'th_lM')
        save(rKymoFn, 'tr_rM', 'dv_rM', 'th_rM')
        save(dKymoFn, 'tr_dM', 'dv_dM', 'th_dM')
        save(vKymoFn, 'tr_vM', 'dv_vM', 'th_vM')
    end
end


%% PathlineIntegratedStrain
% % Now load or collate kymograph data for accumulated strains
% files_exist = exist(apSKymoFn, 'file') && ...
%     exist(lSKymoFn, 'file') && exist(rSKymoFn, 'file') && ...
%     exist(dSKymoFn, 'file') && exist(vSKymoFn, 'file') ;
% if files_exist
%     load(apSKymoFn, 'str_apM', 'sdv_apM', 'sth_apM')
%     load(lSKymoFn, 'str_lM', 'sdv_lM', 'sth_lM')
%     load(rSKymoFn, 'str_rM', 'sdv_rM', 'sth_rM')
%     load(dSKymoFn, 'str_dM', 'sdv_dM', 'sth_dM')
%     load(vSKymoFn, 'str_vM', 'sdv_vM', 'sth_vM')
% else
%     % preallocate for kymos
%     ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
%     sdv_apM   = zeros(ntps, nU) ;   % dv averaged
%     str_apM   = zeros(ntps, nU) ;   
%     sth_apM   = zeros(ntps, nU) ;   
%     sdv_lM   = zeros(ntps, nU) ;    % left averaged
%     str_lM = zeros(ntps, nU) ;
%     sth_lM = zeros(ntps, nU) ;
%     sdv_rM   = zeros(ntps, nU) ;    % right averaged
%     str_rM = zeros(ntps, nU) ;
%     sth_rM = zeros(ntps, nU) ;
%     sdv_dM   = zeros(ntps, nU) ;    % dorsal averaged
%     str_dM = zeros(ntps, nU) ;
%     sth_dM = zeros(ntps, nU) ;
%     sdv_vM   = zeros(ntps, nU) ;    % ventral averaged
%     str_vM = zeros(ntps, nU) ;
%     sth_vM = zeros(ntps, nU) ;
% 
%     for tp = QS.xp.fileMeta.timePoints(1:end-1)
%         close all
%         disp(['t = ' num2str(tp)])
%         tidx = QS.xp.tIdx(tp) ;
% 
%         % Check for timepoint measurement on disk
%         sfn = fullfile(datdir, sprintf('strain_%06d.mat', tp))   ;
%         load(sfn, 'strain_tr_ap', 'strain_tr_l', 'strain_tr_r', ...
%             'strain_tr_d', 'strain_tr_v', ...
%             'strain_th_ap', 'strain_th_l', 'strain_th_r', ....
%             'strain_th_d', 'strain_th_v', ...
%             'strain_dv_ap', 'strain_dv_l', 'strain_dv_r', ...
%             'strain_dv_d', 'strain_dv_v') ;
% 
%         %% Store accumulated strain in matrices
%         % dv averaged
%         str_apM(tidx, :) = strain_tr_ap ;
%         sdv_apM(tidx, :) = strain_dv_ap ;
%         sth_apM(tidx, :) = strain_th_ap ;
% 
%         % left quarter
%         str_lM(tidx, :) = strain_tr_l ;
%         sdv_lM(tidx, :) = strain_dv_l ;
%         sth_lM(tidx, :) = strain_th_l ;
% 
%         % right quarter
%         str_rM(tidx, :) = strain_tr_r ;
%         sdv_rM(tidx, :) = strain_dv_r ;
%         sth_rM(tidx, :) = strain_th_r ;
% 
%         % dorsal quarter
%         str_dM(tidx, :) = strain_tr_d ;
%         sdv_dM(tidx, :) = strain_dv_d ;
%         sth_dM(tidx, :) = strain_th_d ;
% 
%         % ventral quarter
%         str_vM(tidx, :) = strain_tr_v ;
%         sdv_vM(tidx, :) = strain_dv_v ;
%         sth_vM(tidx, :) = strain_th_v ;
%     end
%     
%     %% Save kymograph data
%     save(apSKymoFn, 'str_apM', 'sdv_apM', 'sth_apM')
%     save(lSKymoFn, 'str_lM', 'sdv_lM', 'sth_lM')
%     save(rSKymoFn, 'str_rM', 'sdv_rM', 'sth_rM')
%     save(dSKymoFn, 'str_dM', 'sdv_dM', 'sth_dM')
%     save(vSKymoFn, 'str_vM', 'sdv_vM', 'sth_vM')
% end

%% Load mean curvature kymograph data
try
    metricKDir = QS.dir.metricKinematics.root ;
    dirs = dir(fullfile(metricKDir, ...
        strrep(sprintf('lambda%0.3f_lmesh%0.3f_lerr*_modes%02dw%02d', ...
        lambda, lambda_mesh, nmodes, zwidth), '.', 'p'))) ; 
    try
        assert(length(dirs) == 1)
    catch
        disp('More than one matching metricKDir found:')
        for qq = 1:length(dirs)
            fullfile(metricKDir, dirs(qq).name)
        end
        error('Please ensure that the metricKDir is unique')
    end
    metricKDir = fullfile(metricKDir, dirs(1).name) ;
    loadDir = fullfile(metricKDir, sprintf('pathline_%04dt0', t0Pathline), ...
        'measurements') ;
    apKymoMetricKinFn = fullfile(loadDir, 'apKymographMetricKinematics.mat') ;
    load(apKymoMetricKinFn, 'HH_apM')
    lKymoMetricKinFn = fullfile(loadDir, 'leftKymographMetricKinematics.mat') ;
    load(lKymoMetricKinFn, 'HH_lM')
    rKymoMetricKinFn = fullfile(loadDir, 'rightKymographMetricKinematics.mat') ;
    load(rKymoMetricKinFn, 'HH_rM')
    dKymoMetricKinFn = fullfile(loadDir, 'dorsalKymographMetricKinematics.mat') ;
    load(dKymoMetricKinFn, 'HH_dM')
    vKymoMetricKinFn = fullfile(loadDir, 'ventralKymographMetricKinematics.mat') ;
    load(vKymoMetricKinFn, 'HH_vM')
catch
    error('Run QS.plotPathlineMetricKinematics() before QS.plotPathlineStrainRate()')
end
    
%% Load/compute featureIDs
featureIDs = QS.getPathlineFeatureIDs('vertices', featureIDOptions) ;

%% Store kymograph data in cell arrays
trsK = {0.5*tr_apM, 0.5*tr_lM, 0.5*tr_rM, 0.5*tr_dM, 0.5*tr_vM} ;
dvsK = {dv_apM, dv_lM, dv_rM, dv_dM, dv_vM} ;
thsK = {th_apM, th_lM, th_rM, th_dM, th_vM} ;
HHsK = {HH_apM, HH_lM, HH_rM, HH_dM, HH_vM} ;

%% Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
dvDir = fullfile(mKPDir, 'avgDV') ;
lDir = fullfile(mKPDir, 'avgLeft') ;
rDir = fullfile(mKPDir, 'avgRight') ;
dDir = fullfile(mKPDir, 'avgDorsal') ;
vDir = fullfile(mKPDir, 'avgVentral') ;
outdirs = {dvDir, lDir, rDir, dDir, vDir} ;

%% Now plot different measured quantities as kymographs
clim_zoom = climit / 3 ;
if plot_kymographs
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        trK = trsK{qq} ;
        dvK = dvsK{qq} ;
        thK = thsK{qq} ;
        
        titles = {['dilation rate, ' sRate_trace_label],...
            ['shear rate, ', sRate_deviator_label]} ;
        labels = {[sRate_trace_label ' ' unitstr], ...
            [sRate_deviator_label ' ' unitstr]} ;
        names = {'dilation', 'deviator'} ;
        
        %% Plot strainRate DV-averaged Lagrangian pathline kymographs 
        % Check if images already exist on disk

        % Consider both wide color limit and narrow
        zoomstr = {'_wide', '', '_zoom'} ;
        climits = {climit, 0.5*climit, 0.25*climit} ;
        
        for pp = 1:length(climits)
            %% Plot STRAIN RATE traceful DV-averaged pathline kymograph
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{1} zoomstr{pp} '.png']) ;
            fn_early = fullfile(odir, [names{1} zoomstr{pp} '_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_early, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                imagesc((1:nU)/nU, tps, trK)
                caxis([-climits{pp}, climits{pp}])
                colormap(bbr256)
                
                % Plot fold identifications
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                t1ones = ones(size(tps(fons1:end))) ;
                t2ones = ones(size(tps(fons2:end))) ;
                t3ones = ones(size(tps(fons3:end))) ;
                tidx0 = QS.xp.tIdx(t0) ;

                % OPTION 1: use identified div(v) < 0
                plot(featureIDs(1) * t1ones / nU, tps(fons1:end))
                plot(featureIDs(2) * t2ones / nU, tps(fons2:end))
                plot(featureIDs(3) * t3ones / nU, tps(fons3:end))
                
                % % Add folds to plot
                % hold on;
                % tidx0 = QS.xp.tIdx(t0) ;
                % % Which is the first fold (at t0)?
                % [~, minID] = min(fons) ;
                % thisFoldTimes = tps(fons(minID):end) ;
                % t_ones = ones(size(thisFoldTimes)) ;
                % plot(folds.folds(tidx0, minID) * t_ones / nU, thisFoldTimes)
                % % Now plot the later folds at their Lagrangian locations
                % % determined by div(v) being very negative.
                % laterID = setdiff(1:length(fons), minID) ;
                % for ii = laterID
                %     thisFoldTimes = tps(fons(ii):end) ;
                %     t_ones = ones(size(thisFoldTimes)) ;
                %     plot(featureIDs(ii) * t_ones / nU, thisFoldTimes)
                % end
                
                % Titles 
                title([titles{1}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                
                % title and save
                ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                
                % Zoom in on early times
                ylim([min(tps), min(max(fons) + 10, max(tps))])
                disp(['saving ', fn_early])
                export_fig(fn_early, '-png', '-nocrop', '-r200')  
                clf
            end

            %% DEVIATOR -- strain rate
            fn = fullfile(odir, [ names{2} zoomstr{pp} '.png']) ;
            fn_early = fullfile(odir, [names{2} zoomstr{pp} '_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_early, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                % Map intensity from dev and color from the theta
                indx = max(1, round(mod(2*thK(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
                colors = pm256(indx, :) ;
                devKclipped = min(dvK / climits{pp}, 1) ;
                colorsM = devKclipped(:) .* colors ;
                colorsM = reshape(colorsM, [size(dvK, 1), size(dvK, 2), 3]) ;
                imagesc((1:nU)/nU, tps, colorsM)
                caxis([0, climits{pp}])

                % Add folds to plot
                hold on;

                % Plot fold identifications
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                t1ones = ones(size(tps(fons1:end))) ;
                t2ones = ones(size(tps(fons2:end))) ;
                t3ones = ones(size(tps(fons3:end))) ;

                % OPTION 1: use identified div(v) < 0
                plot(featureIDs(1) * t1ones / nU, tps(fons1:end))
                plot(featureIDs(2) * t2ones / nU, tps(fons2:end))
                plot(featureIDs(3) * t3ones / nU, tps(fons3:end))

                % Colorbar and phasewheel
                colormap(gca, phasemap)
                phasebar('colormap', phasemap, ...
                    'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic') ;
                ax = gca ;
                get(gca, 'position')
                cb = colorbar('location', 'eastOutside') ;
                drawnow
                axpos = get(ax, 'position') ;
                cbpos = get(cb, 'position') ;
                set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3), cbpos(4)*0.6])
                set(ax, 'position', axpos) 
                hold on;
                caxis([0, climits{pp}])
                colormap(gca, gray)

                % title and save
                title([titles{2}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
                    
                % title and save
                ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                
                % Zoom in on early times
                ylim([min(tps), min(max(fons) + 10, max(tps))])
                disp(['saving ', fn_early])
                export_fig(fn_early, '-png', '-nocrop', '-r200')   
            end
        end
    end
end

%% Compare kymo of cumulative STRAIN to summed strainRate raw
if plot_kymographs_strain
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        rtrK = trsK{qq} ;  % strainRate
        rdvK = dvsK{qq} ;
        rthK = thsK{qq} ;
        
        labels = {strain_trace_label, strain_deviator_label} ;
        titles = {['dilation, ' strain_trace_label],...
            ['shear, ', strain_deviator_label]} ;
        names = {'strainRate_trace_t0', 'strain_deviator_t0'} ;
        
        %% Plot STRAIN DV-averaged Lagrangian pathline kymographs 
        % Check if images already exist on disk

        % Consider both wide color limit and narrow
        zoomstr = {'_wide', '', '_zoom'} ;
        climits = {climitWide, 0.5*(climit + climitWide), climit} ;
        
        for pp = 1:length(climits)
            %% Plot STRAIN traceful DV-averaged pathline kymograph
            % Check if filenames exist on disk
            fn = fullfile(odir, [ 'check_' names{1} zoomstr{pp} '.png']) ;
            fn_early = fullfile(odir, [ 'check_' names{1} zoomstr{pp} '_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_early, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                subplot(2, 1, 1)
                imagesc((1:nU)/nU, tps, rtrK)
                caxis([-climits{pp}, climits{pp}])
                colormap(bbr256)
                
                subplot(2, 1, 2)
                tidx0 = QS.xp.tIdx(t0) ;
                imagesc((1:nU)/nU, tps, cumsum(rtrK, 1) - sum(rtrK(1:tidx0, :), 1))
                caxis([-climits{pp}, climits{pp}])
                colormap(bbr256)

                % Plot fold identifications
                for subplotID = 1:2
                    subplot(2, 1, subplotID)
                    hold on;
                    fons1 = max(1, fons(1)) ;
                    fons2 = max(1, fons(2)) ;
                    fons3 = max(1, fons(3)) ;
                    t1ones = ones(size(tps(fons1:end))) ;
                    t2ones = ones(size(tps(fons2:end))) ;
                    t3ones = ones(size(tps(fons3:end))) ;

                    % OPTION 1: use identified div(v) < 0
                    plot(featureIDs(1) * t1ones / nU, tps(fons1:end))
                    plot(featureIDs(2) * t2ones / nU, tps(fons2:end))
                    plot(featureIDs(3) * t3ones / nU, tps(fons3:end))
                end
                
                % Title and save
                subplot(2, 1, 1)
                title([titles{1}, titleadd{1}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
                % Second axis labels/title
                subplot(2, 1, 2)
                title(['summed strain rate', titleadd{1}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on early times
                subplot(2, 1, 1)
                tp_onset = QS.xp.fileMeta.timePoints(fons) - t0 ;
                ylim([max(min(tps), -30), min(max(tp_onset) + 10, max(tps))])
                subplot(2, 1, 2)
                ylim([min(tps), min(max(fons) + 10, max(tps))])
                caxis([-climits{pp}, climits{pp}])
                colormap(bbr256)
                tmp = strsplit(fn_early, filesep) ;
                disp(['saving ', tmp{end}])
                export_fig(fn_early, '-png', '-nocrop', '-r200')   
            end
        
%             %% Plot STRAIN deviator DV-averaged pathline kymograph 
%             % Check if images already exist on disk
%             fn = fullfile(odir, [ 'check_' names{2} zoomstr{pp} '.png']) ;
%             fn_early = fullfile(odir, ['check_' names{2} zoomstr{pp} '_early.png']) ;
%             if ~exist(fn, 'file') || ~exist(fn_early, 'file') || overwrite
% 
%                 close all
%                 set(gcf, 'visible', 'off')
%                 % Map intensity from dev and color from the theta
%                 subplot(2, 1, 1) ;
%                 popts.x = (1:nU)/nU ;
%                 popts.y = tps ; 
%                 popts.colormap = pm256 ;
%                 popts.climit = climits{pp} ;
%                 vectorHeatMap(s_dvK, 2 * s_thK, popts) ;
% 
%                 % Do raw sum
%                 [raw_sdK, raw_stK] = QS.averageNematic(dvK, thK) ;
%                 vectorHeatMap(dvK, thK, popts) ;
%                 
%                 for subplotID = 1:2
%                     % Plot fold identifications
%                     hold on;
%                     fons1 = max(1, fons(1)) ;
%                     fons2 = max(1, fons(2)) ;
%                     fons3 = max(1, fons(3)) ;
%                     t1ones = ones(size(tps(fons1:end))) ;
%                     t2ones = ones(size(tps(fons2:end))) ;
%                     t3ones = ones(size(tps(fons3:end))) ;
%                     tidx0 = QS.xp.tIdx(t0) ;
% 
%                     % OPTION 1: use identified div(v) < 0
%                     plot(featureIDs(1) * t1ones / nU, tps(fons1:end))
%                     plot(featureIDs(2) * t2ones / nU, tps(fons2:end))
%                     plot(featureIDs(3) * t3ones / nU, tps(fons3:end))
%                 end
% 
%                 % Titles 
%                 subplot(2, 1, 1)
%                 title([titles{2}, titleadd{2}], 'Interpreter', 'Latex')
%                 ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
%                 xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
%                 % Second axis labels/title
%                 subplot(2, 1, 2)
%                 title(['summed strain rate', titleadd{1}], 'Interpreter', 'Latex')
%                 ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
%                 xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
% 
%                 % Colorbar and phasewheel
%                 colormap(gca, phasemap)
%                 phasebar('colormap', phasemap, ...
%                     'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic') ;
%                 ax = gca ;
%                 get(gca, 'position')
%                 cb = colorbar('location', 'eastOutside') ;
%                 drawnow
%                 axpos = get(ax, 'position') ;
%                 cbpos = get(cb, 'position') ;
%                 set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3), cbpos(4)*0.6])
%                 set(ax, 'position', axpos) 
%                 hold on;
%                 caxis([0, climits{pp}])
%                 colormap(gca, gray)
% 
%                 % title and save
%                 ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
%                 disp(['saving ', fn])
%                 export_fig(fn, '-png', '-nocrop', '-r200')   
% 
%                 % Zoom in on early times
%                 ylim([min(tps), min(max(fons) + 10, max(tps))])
%                 disp(['saving ', fn_early])
%                 export_fig(fn, '-png', '-nocrop', '-r200')   
%             end
        end
    end
end

%% Kymographs of cumulative STRAIN along pathlines 
if plot_kymographs_strain
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;
    
    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        trK = strK{qq} ;
        
        labels = {strain_trace_label, strain_deviator_label} ;
        titles = {['dilation, ', strain_trace_label],...
            ['shear, ' strain_deviator_label]};
        names = {'strain_trace_t0', 'strain_deviator_t0'} ;
        
        %% Plot STRAIN DV-averaged Lagrangian pathline kymographs 
        % Check if images already exist on disk

        % Consider both wide color limit and narrow
        zoomstr = {'_wide', '', '_zoom'} ;
        climits = {climitWide, 0.5*(climit + climitWide), climit} ;
        
        for pp = 1:length(climits)
            %% Plot STRAIN traceful DV-averaged pathline kymograph
            % Check if filenames exist on disk
            fn = fullfile(odir, [ names{1} zoomstr{pp} '.png']) ;
            fn_early = fullfile(odir, [ names{1} zoomstr{pp} '_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_early, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                imagesc((1:nU)/nU, tps, trK)
                caxis([-climits{pp}, climits{pp}])
                colormap(bbr256)

                % Plot fold identifications
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                t1ones = ones(size(tps(fons1:end))) ;
                t2ones = ones(size(tps(fons2:end))) ;
                t3ones = ones(size(tps(fons3:end))) ;
                tidx0 = QS.xp.tIdx(t0) ;
                
                % OPTION 1: use identified div(v) < 0
                plot(featureIDs(1) * t1ones / nU, tps(fons1:end))
                plot(featureIDs(2) * t2ones / nU, tps(fons2:end))
                plot(featureIDs(3) * t3ones / nU, tps(fons3:end))

                % title and save
                title([titles{1}, titleadd{1}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on early times
                ylim([min(tps), min(max(fons) + 10, max(tps))])
                caxis([-climits{pp}, climits{pp}])
                colormap(bbr256)
                tmp = strsplit(fn_early, filesep) ;
                disp(['saving ', tmp{end}])
                export_fig(fn_early, '-png', '-nocrop', '-r200')   
            end
        
            %% Plot STRAIN deviator DV-averaged pathline kymograph 
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{2} zoomstr{pp} '.png']) ;
            fn_early = fullfile(odir, [names{2} zoomstr{pp} '_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_early, 'file') || overwrite

                close all
                set(gcf, 'visible', 'off')
                % Map intensity from dev and color from the theta
                indx = max(1, round(mod(2*sthK{qq}(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
                colors = pm256(indx, :) ;
                devKclipped = min(sdvK{qq} / climits{pp}, 1) ;
                colorsM = devKclipped(:) .* colors ;
                colorsM = reshape(colorsM, [size(sdvK{qq}, 1), size(sdvK{qq}, 2), 3]) ;
                imagesc((1:nU)/nU, tps, colorsM)

                % Plot fold identifications
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                t1ones = ones(size(tps(fons1:end))) ;
                t2ones = ones(size(tps(fons2:end))) ;
                t3ones = ones(size(tps(fons3:end))) ;
                tidx0 = QS.xp.tIdx(t0) ;

                % OPTION 1: use identified div(v) < 0
                plot(featureIDs(1) * t1ones / nU, tps(fons1:end))
                plot(featureIDs(2) * t2ones / nU, tps(fons2:end))
                plot(featureIDs(3) * t3ones / nU, tps(fons3:end))

                % Titles 
                title([titles{2}, titleadd{2}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')

                % Colorbar and phasewheel
                colormap(gca, phasemap)
                phasebar('colormap', phasemap, ...
                    'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic') ;
                ax = gca ;
                get(gca, 'position')
                cb = colorbar('location', 'eastOutside') ;
                drawnow
                axpos = get(ax, 'position') ;
                cbpos = get(cb, 'position') ;
                set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3), cbpos(4)*0.6])
                set(ax, 'position', axpos) 
                hold on;
                caxis([0, climits{pp}])
                colormap(gca, gray)

                % title and save
                ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on early times
                ylim([min(tps), min(max(fons) + 10, max(tps))])
                disp(['saving ', fn_early])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end

%% Plot 1D curves for region around each fold
% Sample divv near each fold
foldw = 0.05 ;
endw = 0.10 ;
cut = [round(endw*nU), featureIDs(1)-round(foldw*nU), ...
        featureIDs(1)+round(foldw*nU), featureIDs(2)-round(foldw*nU), ...
        featureIDs(2)+round(foldw*nU), featureIDs(3)-round(foldw*nU), ... 
        featureIDs(3)+round(foldw*nU), round((1-endw) * nU)] ;
lobes = { cut(1):cut(2), cut(3):cut(4), cut(5):cut(6), cut(7):cut(8) } ;
avgStrings = {'dv-averaged', 'left side', 'right side', ...
    'dorsal side', 'ventral side'} ;
avgLabel = {'dv', 'left', 'right', 'dorsal', 'ventral'} ;
titleRateFoldBase = 'Lagrangian pathline strain rate in folds, ' ;
titleRateLobeBase = 'Lagrangian pathline strain rate in lobes, ' ;
titleFoldBase = 'Lagrangian pathline strain in folds, ' ;
titleLobeBase = 'Lagrangian pathline strain in lobes, ' ;

if length(featureIDs) == 3
    foldYlabels = {'anterior fold', 'middle fold', 'posterior fold'} ;
    lobeYlabels = {'lobe 1', 'lobe 2', 'lobe 3', 'lobe 4'} ;
    foldYlabelsStrainRate = {...
        {'anterior fold'; 'strain rate, $\dot{\varepsilon}$'}, ...
        {'middle fold'; 'strain rate, $\dot{\varepsilon}$'}, ...
        {'posterior fold'; 'strain rate, $\dot{\varepsilon}$'}} ;
    lobeYlabelsStrainRate = {...
        {'lobe 1'; 'strain rate, $\dot{\varepsilon}$'}, ...
        {'lobe 2'; 'strain rate, $\dot{\varepsilon}$'}, ...
        {'lobe 3'; 'strain rate, $\dot{\varepsilon}$'}, ...
        {'lobe 4'; 'strain rate, $\dot{\varepsilon}$'}};
    foldYlabelsStrain = {...
        {'anterior fold'; 'strain, $\varepsilon$'}, ...
        {'middle fold'; 'strain, $\varepsilon$'}, ...
        {'posterior fold'; 'strain, $\varepsilon$'}} ;
    lobeYlabelsStrain = {...
        {'lobe 1'; 'strain, $\varepsilon$'}, ...
        {'lobe 2'; 'strain, $\varepsilon$'}, ...
        {'lobe 3'; 'strain, $\varepsilon$'}, ...
        {'lobe 4'; 'strain, $\varepsilon$'}};
    foldYlabelsStrainRateComposite = {...
        ['anterior fold, ' sRate_trace_label], ...
        ['anterior fold, ' sRate_deviator_label_short], ...
        ['middle fold, ' sRate_trace_label], ...
        ['middle fold, ' sRate_deviator_label_short], ...
        ['posterior fold, ' sRate_trace_label], ...
        ['posterior fold, ' sRate_deviator_label_short]} ;
    foldYlabelsStrainComposite = {...
        ['anterior fold, ' strain_trace_label], ...
        ['anterior fold, ' strain_deviator_label_short], ...
        ['middle fold, ' strain_trace_label], ...
        ['middle fold, ' strain_deviator_label_short], ...
        ['posterior fold, ' strain_trace_label], ...
        ['posterior fold, ' strain_deviator_label_short]} ;
    % foldYlabelsStrainRatio = {...
    %     'anterior fold, $\frac{||\mathrm{Dev}\left[\varepsilon\right]||}{ ||\mathrm{Tr} [\varepsilon]||}$', ...
    %     'middle fold, $\frac{||\mathrm{Dev}\left[\varepsilon\right]||}{ ||\mathrm{Tr} [\varepsilon]||}$', ...
    %     'posterior fold, $\frac{||\mathrm{Dev}\left[\varepsilon\right]||}{ ||\mathrm{Tr} [\varepsilon]||}$'} ;
    lobeYlabelsStrainRateComposite = {...
        ['lobe 1, ' sRate_trace_label], ...
        ['lobe 1, ' sRate_deviator_label_short], ...
        ['lobe 2, ' sRate_trace_label], ...
        ['lobe 2, ' sRate_deviator_label_short], ...
        ['lobe 3, ' sRate_trace_label], ...
        ['lobe 3, ' sRate_deviator_label_short], ...
        ['lobe 4, ' sRate_trace_label], ...
        ['lobe 4, ' sRate_deviator_label_short]} ;
    lobeYlabelsStrainComposite = {...
        ['lobe 1, ' strain_trace_label], ...
        ['lobe 1, ' strain_deviator_label_short], ...
        ['lobe 2, ' strain_trace_label], ...
        ['lobe 2, ' strain_deviator_label_short], ...
        ['lobe 3, ' strain_trace_label], ...
        ['lobe 3, ' strain_deviator_label_short], ...
        ['lobe 4, ' strain_trace_label], ...
        ['lobe 4, ' strain_deviator_label_short]} ;
    % lobeYlabelsStrainRatio = {...
    %   'lobe 1 $\frac{||\mathrm{Dev}\left[\varepsilon\right]||}{ ||\mathrm{Tr} [\varepsilon]||}$', ...
    %   'lobe 2, $\frac{||\mathrm{Dev}\left[\varepsilon\right]||}{ ||\mathrm{Tr} [\varepsilon]||}$', ...
    %   'lobe 3, $\frac{||\mathrm{Dev}\left[\varepsilon\right]||}{ ||\mathrm{Tr} [\varepsilon]||}$', ...
    %   'lobe 4, $\frac{||\mathrm{Dev}\left[\varepsilon\right]||}{ ||\mathrm{Tr} [\varepsilon]||}$'} ;
else
    error('How many featureIDs? Allow variability here')
end
for qq = 1:length(trsK)
    %% STRAIN RATE
    trK = trsK{qq} ;
    dvK = dvsK{qq} ;
    thK = thsK{qq} ;
    HHK = HHsK{qq} ;
    
    % FOLDS: Explore a range of widths for plotting
    for width = 1:round(0.03 * nU)
        % Define the regions that are considered folds in Lagrangian
        % longitudinal coordinates 
        foldRegions = cell(length(featureIDs), 1) ;
        for jj = 1:length(featureIDs)
            foldRegions{jj} = (featureIDs(jj)-width):(featureIDs(jj)+width) ;
        end
        
        %% fold kinematics -- Strain rate, one axis for each fold
        fn = fullfile(outdirs{qq}, ...
            [sprintf('fold_strainRate_w%03d_', 2*width+1), ...
            avgLabel{qq}, '.png']) ;
        fn_withH = fullfile(outdirs{qq}, ...
            [sprintf('fold_strainRate_w%03d_', 2*width+1), ...
            avgLabel{qq}, '_withH.png']) ;
        if plot_fold_strainRate
            plotOpts.overwrite = overwrite ;
            plotOpts.H_on_yyaxis = true ; 
            
            aux_plotPathlineStrainFeatures_subpanels(QS, fn, fn_withH, ...
            featureIDs, width, nU, tps, trK, dvK, thK, HHK, ...
            avgStrings{qq}, titleRateFoldBase, ...
            foldYlabelsStrainRate, sRate_trace_label, sRate_deviator_label, ...
            trecolor, Hposcolor, Hnegcolor, Hsz, pm256, plotOpts)
        end
        
        %% fold kinematics -- Strain rate, all folds on one axis
        fnbase = fullfile(outdirs{qq}, ...
            [sprintf('fold_strainRate_compare_w%03d_', 2*width+1), ...
            avgLabel{qq}]) ;
        fns.fn = [fnbase '.png'] ;
        fns.early = [fnbase '_early.png'] ;
        fns.norms = [fnbase '_norms.png'] ;
        fns.norms_early = [fnbase '_norms_early.png'] ;
        fns.withH = [fnbase '_withH.png'] ;
        plotOpts.overwrite = overwrite ;
        plotOpts.phasecmap = pm256 ;
        measurements.fons = fons - t0 ;
        measurements.regions = foldRegions ;
        data.timepoints = tps ;
        data.tr = trK ;
        data.dv = dvK ;
        data.th = thK ;
        data.HH = HHK ;
        labels.avgString = [avgStrings{qq}, ', ', ...
            '$w_{\textrm{fold}}=', ...
            num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'];
        labels.titleBase = titleFoldBase ;
        labels.ylabel = 'fold strain rate, $\langle\dot{\varepsilon}\rangle$' ;
        labels.ylabel_ratios = 'fold strain rate, $\langle||\mathrm{Tr}[\dot{\varepsilon}]||\rangle/\langle||\mathrm{Dev}[\dot{\varepsilon}]||\rangle$' ;
        labels.legend = foldYlabelsStrainRateComposite ;
        labels.legend_ratios = foldYlabels ;
        labels.trace = strain_trace_label ;
        labels.deviator = strain_deviator_label ;
        aux_plotPathlineStrainRegions(QS, ...
            fns, measurements, data, labels, plotOpts)
    end
       
    %% Plot lobe Kinematics -- STRAIN RATE, one axis per lobe
    if plot_lobe_strainRate
        fn = fullfile(outdirs{qq}, ...
            ['lobe_strainRate_' avgLabel{qq} '.png']) ;
        fn_withH = fullfile(outdirs{qq}, ...
            ['lobe_strainRate_' avgLabel{qq} '_withH.png']) ;
        aux_plotPathlineStrainLobes_subpanels(QS, ...
            fn, fn_withH, ...
            lobes, tps, trK, dvK, thK, HHK, avgStrings{qq}, ...
            titleRateLobeBase, lobeYlabelsStrainRate, ...
            sRate_trace_label, sRate_deviator_label, trecolor,...
            Hposcolor, Hnegcolor, Hsz, pm256, overwrite) 
    end
    
    %% Lobe STRAIN RATE -- all lobes on one axis
    if plot_lobe_strainRate
        fnbase = fullfile(outdirs{qq}, 'lobe_strainRate_compare') ;
        fns.fn = [fnbase '.png'] ;
        fns.early = [fnbase '_early.png'] ;
        fns.norms = [fnbase '_norms.png'] ;
        fns.norms_early = [fnbase '_norms_early.png'] ;
        fns.withH = [fnbase '_withH.png'] ;
        plotOpts.overwrite = overwrite ;
        plotOpts.phasecmap = pm256 ;
        measurements.fons = fons - t0 ;
        measurements.regions = lobes ;
        data.timepoints = tps ;
        data.tr = trK ;
        data.dv = dvK ;
        data.th = thK ;
        data.HH = HHK ;
        labels.avgString = avgStrings{qq} ;
        labels.titleBase = titleLobeBase ;
        labels.ylabel = 'lobe strain rate, $\langle\dot{\varepsilon}\rangle$' ;
        labels.ylabel_ratios = 'lobe strain rate, $\langle||\mathrm{Tr}[\dot{\varepsilon}]||\rangle/\langle||\mathrm{Dev}[\dot{\varepsilon}]||\rangle$' ;
        labels.legend = lobeYlabelsStrainRateComposite ;
        labels.legend_ratios = lobeYlabels ;
        labels.trace = strain_trace_label ;
        labels.deviator = strain_deviator_label ;
        aux_plotPathlineStrainRegions(QS, ...
            fns, measurements, data, labels, plotOpts)     
    end
    
    %% Fold kinematics -- STRAIN, one fold on each axis
    trK = strK{qq} ;
    dvK = sdvK{qq} ;
    thK = sthK{qq} ;
    HHK = HHsK{qq} ;
    for width = 1:round(0.05 * nU)
        % Define the regions that are considered folds in Lagrangian
        % longitudinal coordinates 
        foldRegions = cell(length(featureIDs), 1) ;
        for jj = 1:length(featureIDs)
            foldRegions{jj} = (featureIDs(jj)-width):(featureIDs(jj)+width) ;
        end

        % Define filenames
        fn = fullfile(outdirs{qq}, ...
            [sprintf('fold_strain_w%03d_', 2*width+1), ...
            avgLabel{qq}, '.png']) ;
        fn_withH = fullfile(outdirs{qq}, ...
            [sprintf('fold_strain_w%03d_', 2*width+1), ...
            avgLabel{qq}, '_withH.png']) ;
        
        if plot_fold_strain
            plotOpts.overwrite = overwrite ;
            plotOpts.H_on_yyaxis = true ; 
            aux_plotPathlineStrainFeatures_subpanels(QS, fn, fn_withH, ...
            featureIDs, width, nU, tps, trK, dvK, thK, HHK, ...
            avgStrings{qq}, titleFoldBase, foldYlabelsStrain, ...
            strain_trace_label, strain_deviator_label, ...
            trecolor, Hposcolor, Hnegcolor, Hsz, pm256, plotOpts)
        end
        
        %% Plot STRAIN -- all folds on one axis
        if plot_fold_strain
            fnbase = fullfile(outdirs{qq}, ...
                [sprintf('fold_strain_compare_w%03d_', 2*width+1), ...
                avgLabel{qq}]) ;
            fns.fn = [fnbase '.png'] ;
            fns.early = [fnbase '_early.png'] ;
            fns.norms = [fnbase '_norms.png'] ;
            fns.norms_early = [fnbase '_norms_early.png'] ;
            fns.withH = [fnbase '_withH.png'] ;
            plotOpts.overwrite = overwrite ;
            plotOpts.phasecmap = pm256 ;
            measurements.fons = fons - t0 ;
            measurements.regions = foldRegions ;
            data.timepoints = tps ;
            data.tr = trK ;
            data.dv = dvK ;
            data.th = thK ;
            data.HH = HHK ;
            labels.avgString = [avgStrings{qq}, ', ', ...
                '$w_{\textrm{fold}}=', ...
                num2str(100*(2*width + 1)/ nU), '$\%$\, L_\zeta$'];
            labels.titleBase = titleFoldBase ;
            labels.ylabel = 'fold strain, $\langle\varepsilon\rangle$' ;
            labels.ylabel_ratios = 'fold strain, $\langle||\mathrm{Dev}[\varepsilon]||\rangle/\langle||\mathrm{Tr}[\varepsilon]||\rangle$' ;
            labels.legend = foldYlabelsStrainComposite ;
            labels.legend_ratios = foldYlabels ;
            labels.trace = strain_trace_label ;
            labels.deviator = strain_deviator_label ;
            aux_plotPathlineStrainRegions(QS, ...
                fns, measurements, data, labels, plotOpts)
        end
    end
    
    %% Lobe kinematics -- STRAIN, one axis per lobe
    fn = fullfile(outdirs{qq}, sprintf('lobe_strain.png')) ;
    fn_withH = fullfile(outdirs{qq}, sprintf('lobe_strain_withH.png')) ;
    if plot_lobe_strain
        plotOpts.overwrite = overwrite ;
        aux_plotPathlineStrainLobes_subpanels(QS, fn, fn_withH, ...
            lobes, tps, trK, dvK, thK, HHK, ...
            avgStrings{qq}, titleLobeBase, lobeYlabelsStrain, ...
            strain_trace_label, strain_deviator_label, ...
            trecolor, Hposcolor, Hnegcolor, Hsz, pm256, plotOpts)
    end

    %% lobe kinematics -- STRAIN, all lobes on one axis
    if plot_lobe_strain
        fnbase = fullfile(outdirs{qq}, 'lobe_strain_compare') ;
        fns.fn = [fnbase '.png'] ;
        fns.early = [fnbase '_early.png'] ;
        fns.norms = [fnbase '_norms.png'] ;
        fns.norms_early = [fnbase '_norms_early.png'] ;
        fns.withH = [fnbase '_withH.png'] ;
        plotOpts.overwrite = overwrite ;
        plotOpts.phasecmap = pm256 ;
        measurements.fons = fons - t0 ;
        measurements.regions = lobes ;
        data.timepoints = tps ;
        data.tr = trK ;
        data.dv = dvK ;
        data.th = thK ;
        data.HH = HHK ;
        labels.avgString = avgStrings{qq} ;
        labels.titleBase = titleLobeBase ;
        labels.ylabel = 'lobe strain, $\langle\varepsilon\rangle$' ;
        labels.ylabel_ratios = 'lobe strain, $\langle||\mathrm{Dev}[\varepsilon]||\rangle/\langle||\mathrm{Tr}[\varepsilon]||\rangle$' ;
        labels.legend = lobeYlabelsStrainComposite ;
        labels.legend_ratios = lobeYlabels ;
        labels.trace = strain_trace_label ;
        labels.deviator = strain_deviator_label ;
        aux_plotPathlineStrainRegions(QS, ...
            fns, measurements, data, labels, plotOpts)     
    end

end

disp('done')



