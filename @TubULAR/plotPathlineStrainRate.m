function plotPathlineStrainRate(tubi, options)
% plotPathlineStrainRate(tubi, options)
%   Plot the strain rate along pathlines as kymographs and 
%   correlation plots. These are computed via finite differencing of 
%   pathline mesh metrics. That is, vertices are advected along piv flow 
%   and projected into 3D from pullback space to embedding, then the
%   metrics g_{ij} of those meshes are measured and differenced. 
% 
%
% Parameters
% ----------
% tubi : QuapSlap class instance
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

%% Default options 
overwrite = false ;
plot_kymographs = true ;
t0 = tubi.t0set() ;
t0Pathline = t0 ;

sRate_trace_label = '$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\dot{\varepsilon}]$';
sRate_deviator_label = ...
    '$||\dot{\varepsilon}-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\dot{\varepsilon}\right]\bf{g}||$' ;
sRate_deviator_label_short = '$||$Dev$\left[\dot{\varepsilon}\right]||$' ;

%% Parameter options
lambda = tubi.smoothing.lambda ;
lambda_mesh = tubi.smoothing.lambda_mesh ;
nmodes = tubi.smoothing.nmodes ;
zwidth = tubi.smoothing.zwidth ;
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
else 
    error("Could not parse samplingResolution: set to '1x'")
end

%% Unpack tubi
tubi.getXYZLims ;
xyzlim = tubi.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
fons = tubi.t0set() - tubi.xp.fileMeta.timePoints(1) ;

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
bluecolor = tubi.plotting.colors(1, :) ;
orangecolor = tubi.plotting.colors(2, :) ;
yellowcolor = tubi.plotting.colors(3, :) ;
purplecolor = tubi.plotting.colors(4, :) ;
greencolor = tubi.plotting.colors(5, :) ;
graycolor = tubi.plotting.colors(8, :) ;
browncolor = tubi.plotting.colors(9, :) ;

%% Choose colors
trecolor = 'k' ;  % yellowcolor ;
Hposcolor = greencolor ;
Hnegcolor = purplecolor ;
Hsz = 3 ;  % size of scatter markers for mean curvature

%% load from tubi
if doubleResolution
    nU = tubi.nU * 2 - 1 ;
    nV = tubi.nV * 2 - 1 ;
else
    nU = tubi.nU ;
    nV = tubi.nV ;    
end

%% Test incompressibility of the flow on the evolving surface
% We relate the normal velocities to the divergence / 2 * H.
tps = tubi.xp.fileMeta.timePoints(1:end-1) - t0 ;

% Output directory is inside metricKinematics dir
mKPDir = sprintf(tubi.dir.strainRate.pathline.root, t0Pathline) ;
datdir = fullfile(mKPDir, 'measurements') ;
% Data for kinematics on meshes (defined on vertices) [not needed here]
% mdatdir = fullfile(mKDir, 'measurements') ;

% Unit definitions for axis labels
unitstr = [ '[1/' tubi.timeUnits ']' ];

%% Compute or load all timepoints
apKymoFn = fullfile(datdir, 'apKymographPathlineStrainRate.mat') ;
lKymoFn = fullfile(datdir, 'leftKymographPathlineStrainRate.mat') ;
rKymoFn = fullfile(datdir, 'rightKymographPathlineStrainRate.mat') ;
dKymoFn = fullfile(datdir, 'dorsalKymographPathlineStrainRate.mat') ;
vKymoFn = fullfile(datdir, 'ventralKymographPathlineStrainRate.mat') ;
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
    ntps = length(tubi.xp.fileMeta.timePoints(1:end-1)) ;
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

    for tp = tubi.xp.fileMeta.timePoints(1:end-1)
        close all
        disp(['t = ' num2str(tp)])
        tidx = tubi.xp.tIdx(tp) ;

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

%% Load mean curvature kymograph data
% This information should have been saved during
% plotPathlineMetricKinematics() earlier.
try
    metricKDir = tubi.dir.metricKinematics.root ;
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
    error('Run tubi.plotPathlineMetricKinematics() before tubi.plotPathlineStrainRate()')
end
    
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

for odirId = 1:length(outdirs)
    if ~exist(outdirs{odirId}, 'dir')
        mkdir(outdirs{odirId})
    end
end

%% Now plot different measured quantities as kymographs
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
        names = {'dilation_rate', 'deviator_rate'} ;
        
        %% Plot strainRate DV-averaged Lagrangian pathline kymographs 
        % Check if images already exist on disk

        % Consider both wide color limit and narrow
        zoomstr = {'_wide', '', '_zoom'} ;
        climits = {climit, 0.5*climit, 0.25*climit} ;
        
        for pp = 1:length(climits)
            %% Plot STRAIN RATE traceful DV-averaged pathline kymograph
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{1} zoomstr{pp} '.png']) ;
            if ~exist(fn, 'file')  || overwrite
                close all
                set(gcf, 'visible', 'off')
                imagesc((1:nU)/nU, tps, trK)
                caxis([-climits{pp}, climits{pp}])
                colormap(bbr256)
                
                % Titles 
                title([titles{1}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                
                % title and save
                ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                
                clf
            end

            %% DEVIATOR -- strain rate
            fn = fullfile(odir, [ names{2} zoomstr{pp} '.png']) ;
            if ~exist(fn, 'file')  || overwrite
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
                ylabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
                    
                % title and save
                ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                clf
            end
        end
    end
end

disp('done')



