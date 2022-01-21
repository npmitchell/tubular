function plotPathlineMetricKinematics(QS, options)
% plotPathlineMetricKinematics(QS, options)
%   Plot the metric Kinematics along pathlines as kymographs and 
%   correlation plots.
%   Out-of-plane motion is v_n * 2H, where v_n is normal velocity and H is
%   mean curvature.
%   In-plane motion considered here is div(v_t) where v_t is tangential
%   velocity on the curved surface.
%   The difference div(v_t) - vn*2H = 1/2 Tr[g^{-1} dot{g}], which is a 
%   measure of isotropic metric change over time (dot is the partial 
%   derivative wrt time). 
% 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs         : bool
%   plot_kymographs_cumsum  : bool
%   plot_correlations       : bool 
%   plot_fold_kinematics    : bool
%   plot_lobe_kinematics    : bool
% 
% Returns
% -------
% <none>
%
% NPMitchell 2020

%% Default options 
overwrite = false ;
plot_kymographs = true ;
plot_kymographs_cumsum = false ;
plot_kymographs_cumprod = false ;
plot_correlations = true ;
plot_fold_kinematics = true ;
plot_lobe_kinematics = true ;
plot_cumsum_cumprod = false ;
plot_ap = true ;
plot_left = true ;
plot_right = true ;
plot_dorsal = true ;
plot_ventral = true ;
maxWFrac = 0.03 ;
% Load time offset for first fold, t0
t0 = QS.t0set() ;  
t0Pathline = QS.t0 ;

%% Parameter options
lambda = QS.smoothing.lambda ;
lambda_mesh = QS.smoothing.lambda_mesh ;
lambda_err = QS.smoothing.lambda_err ;
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
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_err')
    lambda_err = options.lambda_err ;
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
if isfield(options, 't0Pathline')
    t0 = options.t0Pathline ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'climit_err')
    climit_err = options.climit_err ;
else
    % default is equal to climit
    climit_err = climit ;
end
if isfield(options, 'climit_veln')
    climit_veln = options.climit_veln ;
else
    % default is equal to 3 um/min
    climit_veln = 3 ;
end
if isfield(options, 'climit_H')
    climit_H = options.climit_H ;
else
    % default is equal to 2*climit
    climit_H = 0.2 ;
end
if isfield(options, 'climit_radius')
    climit_radius = options.climit_radius ;
else
    climit_radius = 0 ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
% if isfield(options, 'climitWide')
%     climitWide = options.climitWide ;
% else
%     climitWide = climit * 3; 
% end
if isfield(options, 'maxWFrac')
    maxWFrac = options.maxWFrac ;
end

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_kymographs_cumsum')
    plot_kymographs_cumsum = options.plot_kymographs_cumsum ;
end
if isfield(options, 'plot_kymographs_cumprod')
    plot_kymographs_cumprod = options.plot_kymographs_cumprod ;
end
if isfield(options, 'plot_cumsum_cumprod')
    plot_cumsum_cumprod = options.plot_cumsum_cumprod ;
end
if isfield(options, 'plot_correlations')
    plot_correlations = options.plot_correlations ;
end
if isfield(options, 'plot_ap')
    plot_ap = options.plot_ap ;
end
if isfield(options, 'plot_left')
    plot_left = options.plot_left ;
end
if isfield(options, 'plot_right')
    plot_right = options.plot_right ;
end
if isfield(options, 'plot_dorsal')
    plot_dorsal = options.plot_dorsal ;
end
if isfield(options, 'plot_ventral')
    plot_ventral = options.plot_ventral ;
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
mKDir = fullfile(QS.dir.metricKinematics.root, ...
    strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
    lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));

%% Get fold locations for pathlines
try
    folds = load(sprintf(QS.fileName.pathlines.featureIDs, t0Pathline)) ;
catch
    disp('Feature IDs for pathlines not stored on disk. Defining...')
    folds = QS.measurePathlineFeatureIDs ;
end
foldNonPathline = load(QS.fileName.fold) ;
fons = foldNonPathline.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
bwr256 = bluewhitered(256) ;
bluecolor = QS.plotting.colors(1, :) ;
orangecolor = QS.plotting.colors(2, :) ;
yellowcolor = QS.plotting.colors(3, :) ;
purplecolor = QS.plotting.colors(4, :) ;
greencolor = QS.plotting.colors(5, :) ;
graycolor = QS.plotting.colors(8, :) ;
browncolor = QS.plotting.colors(9, :) ;

%% Choose colors
divvcolor = bluecolor ;
H2vncolor = orangecolor ;
gdotcolor = yellowcolor ;
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
tps = QS.xp.fileMeta.timePoints(1:end-1) - t0;

% preallocate for cumulative error
ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
HH_apM   = zeros(ntps, nU) ;   % dv averaged
divv_apM = zeros(ntps, nU) ;
veln_apM = zeros(ntps, nU) ;
gdot_apM = zeros(ntps, nU) ;
radi_apM = zeros(ntps, nU) ;
HH_lM   = zeros(ntps, nU) ;    % left averaged
divv_lM = zeros(ntps, nU) ;
veln_lM = zeros(ntps, nU) ;
gdot_lM = zeros(ntps, nU) ;
radi_lM = zeros(ntps, nU) ;
HH_rM   = zeros(ntps, nU) ;    % right averaged
divv_rM = zeros(ntps, nU) ;
veln_rM = zeros(ntps, nU) ;
gdot_rM = zeros(ntps, nU) ;
radi_rM = zeros(ntps, nU) ;
HH_dM   = zeros(ntps, nU) ;    % dorsal averaged
divv_dM = zeros(ntps, nU) ;
veln_dM = zeros(ntps, nU) ;
gdot_dM = zeros(ntps, nU) ;
radi_dM = zeros(ntps, nU) ;
HH_vM   = zeros(ntps, nU) ;    % ventral averaged
divv_vM = zeros(ntps, nU) ;
veln_vM = zeros(ntps, nU) ;
gdot_vM = zeros(ntps, nU) ;
radi_vM = zeros(ntps, nU) ;

% Output directory is inside metricKinematics dir
mKPDir = fullfile(mKDir, sprintf('pathline_%04dt0', t0Pathline)) ;
datdir = fullfile(mKPDir, 'measurements') ;
% Data for kinematics on meshes (defined on vertices) [not needed here]
% mdatdir = fullfile(mKDir, 'measurements') ;

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];
Hunitstr = [ '[1/' QS.spaceUnits ']' ];
vunitstr = [ '[' QS.spaceUnits '/' QS.timeUnits ']' ];

% colormap prep    
close all
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
close all

% Compute or load all timepoints
apKymoFn = fullfile(datdir, 'apKymographMetricKinematics.mat') ;
lKymoFn = fullfile(datdir, 'leftKymographMetricKinematics.mat') ;
rKymoFn = fullfile(datdir, 'rightKymographMetricKinematics.mat') ;
dKymoFn = fullfile(datdir, 'dorsalKymographMetricKinematics.mat') ;
vKymoFn = fullfile(datdir, 'ventralKymographMetricKinematics.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if files_exist
    load(apKymoFn, 'HH_apM', 'gdot_apM', 'divv_apM', 'veln_apM', ...
        'H2vn_apM', 'radius_apM')
    load(lKymoFn, 'HH_lM', 'gdot_lM', 'divv_lM', 'veln_lM', ...
        'H2vn_lM', 'radius_lM')
    load(rKymoFn, 'HH_rM', 'gdot_rM', 'divv_rM', 'veln_rM', ...
        'H2vn_rM', 'radius_rM')
    load(dKymoFn, 'HH_dM', 'gdot_dM', 'divv_dM', 'veln_dM', ...
        'H2vn_dM', 'radius_dM')
    load(vKymoFn, 'HH_vM', 'gdot_vM', 'divv_vM', 'veln_vM', ...
        'H2vn_vM', 'radius_vM')
else
    disp('Collating data into kymographs')
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        close all
        disp(['t = ' num2str(tp)])
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        Hfn = fullfile(datdir, sprintf('HH_series_%06d.mat', tp))   ;
        efn = fullfile(datdir, sprintf('gdot_series_%06d.mat', tp)) ;
        dfn = fullfile(datdir, sprintf('divv_series_%06d.mat', tp)) ;
        nfn = fullfile(datdir, sprintf('veln_series_%06d.mat', tp)) ;
        H2vnfn = fullfile(datdir, sprintf('H2vn_series_%06d.mat', tp)) ;

        % Load timeseries measurements
        load(Hfn, 'HH', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
        load(efn, 'gdot', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
        load(dfn, 'divv', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
        load(nfn, 'veln', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
        load(H2vnfn, 'H2vn', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 
        load(radifn, 'radius', 'radius_ap', 'radius_l', 'radius_r', 'radius_d', 'radius_v') 

        %% Store in matrices
        % dv averaged
        HH_apM(tidx, :) = HH_ap ;
        gdot_apM(tidx, :) = gdot_ap ;
        divv_apM(tidx, :) = divv_ap ;
        veln_apM(tidx, :) = veln_ap ;
        H2vn_apM(tidx, :) = H2vn_ap ;
        radi_apM(tidx, :) = radius_ap ;

        % left quarter
        HH_lM(tidx, :) = HH_l ;
        gdot_lM(tidx, :) = gdot_l ;
        divv_lM(tidx, :) = divv_l ;
        veln_lM(tidx, :) = veln_l ;
        H2vn_lM(tidx, :) = H2vn_l ;
        radi_lM(tidx, :) = radius_l ;

        % right quarter
        HH_rM(tidx, :) = HH_r ;
        gdot_rM(tidx, :) = gdot_r ;
        divv_rM(tidx, :) = divv_r ;
        veln_rM(tidx, :) = veln_r ;
        H2vn_rM(tidx, :) = H2vn_r ;
        radi_rM(tidx, :) = radius_r ;

        % dorsal quarter
        HH_dM(tidx, :) = HH_d ;
        gdot_dM(tidx, :) = gdot_d ;
        divv_dM(tidx, :) = divv_d ;
        veln_dM(tidx, :) = veln_d ;
        H2vn_dM(tidx, :) = H2vn_d ;
        radi_dM(tidx, :) = radius_d ;

        % ventral quarter
        HH_vM(tidx, :) = HH_v ;
        gdot_vM(tidx, :) = gdot_v ;
        divv_vM(tidx, :) = divv_v ;
        veln_vM(tidx, :) = veln_v ;
        H2vn_vM(tidx, :) = H2vn_v ;
        radi_vM(tidx, :) = radius_v ;
    end
    
    % Save the DV-averaged kymographs
    disp('Saving compiled kymograph data to disk')
    save(apKymoFn, 'HH_apM', 'gdot_apM', 'divv_apM', ...
        'veln_apM', 'H2vn_apM', 'radius_apM')
    save(lKymoFn, 'HH_lM', 'gdot_lM', 'divv_lM', 'veln_lM', 'H2vn_lM', 'radius_lM')
    save(rKymoFn, 'HH_rM', 'gdot_rM', 'divv_rM', 'veln_rM', 'H2vn_rM', 'radius_rM')
    save(dKymoFn, 'HH_dM', 'gdot_dM', 'divv_dM', 'veln_dM', 'H2vn_dM', 'radius_dM')
    save(vKymoFn, 'HH_vM', 'gdot_vM', 'divv_vM', 'veln_vM', 'H2vn_vM', 'radius_vM')
end

%% Store kymograph data in cell arrays
HHsK = {HH_apM, HH_lM, HH_rM, HH_dM, HH_vM} ;
gdotsK = {gdot_apM, gdot_lM, gdot_rM, gdot_dM, gdot_vM} ;
divvsK = {divv_apM, divv_lM, divv_rM, divv_dM, divv_vM} ;
velnsK = {veln_apM, veln_lM, veln_rM, veln_dM, veln_vM} ;
H2vnsK = {H2vn_apM, H2vn_lM, H2vn_rM, H2vn_dM, H2vn_vM} ;
radisK = {radius_apM, radius_lM, radius_rM, radius_dM, radius_vM} ;

%% Obtain location of folds / features
featureOpts = struct() ;
featureIDs = QS.getPathlineFeatureIDs('vertices', featureOpts) ;

%% Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
dvDir = fullfile(mKPDir, 'avgDV') ;
lDir = fullfile(mKPDir, 'avgLeft') ;
rDir = fullfile(mKPDir, 'avgRight') ;
dDir = fullfile(mKPDir, 'avgDorsal') ;
vDir = fullfile(mKPDir, 'avgVentral') ;
outdirs = {dvDir, lDir, rDir, dDir, vDir} ;

%% Now plot different measured quantities as kymographs
do_plots = [plot_ap, plot_left, plot_right, plot_dorsal, plot_ventral] ;
if plot_kymographs
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = find(do_plots)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        gdotK = gdotsK{qq} ;
        divvK = divvsK{qq} ;
        velnK = velnsK{qq} ;
        H2vnK = H2vnsK{qq} ;
        radiK = radisK{qq} ;
        m2plot = {gdotK, HHK, divvK, velnK, H2vnK, radiK} ;
        titles = {'$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'mean curvature, $H$', ...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $v_n$', ...
            'normal motion, $v_n 2 H$', ...
            'radius'} ;
        labels = {['$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['mean curvature, $H$ ' Hunitstr], ...
            ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $v_n$ ' vunitstr] , ...
            ['normal motion, $v_n 2 H $ ' unitstr], ...
            ['radius [' QS.spaceUnits ']']} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn', 'radi'} ;
        climits = [climit, climit_H, climit, climit_veln, climit_err, climit_radius] ;

        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{pp} '.png']) ;
            fn_zoom = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            
            if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                imagesc((1:nU)/nU, tps, m2plot{pp})
                if climits(pp) > 0
                    caxis([-climits(pp), climits(pp)])
                end
                colormap(bwr256)

                % Plot feature/fold identifications
                hold on;
                for ii = 1:length(fons)
                    fonsi = max(1, fons(ii)) ;
                    tones = ones(size(tps(fonsi:end))) ;
                    plot(featureIDs(ii) * tones / nU, tps(fonsi:end))
                end
                
                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                fn = fullfile(odir, [ names{pp} '.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on small values
                caxis([-climits(pp)/3, climits(pp)/3])
                colormap(bwr256)
                fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                caxis([-climits(pp)/3, climits(pp)/3])
                colormap(bwr256)
                fn = fullfile(odir, [names{pp} '_zoom_early.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end

%% Kymographs of cumulative sums along pathlines
do_plots = [plot_ap, plot_left, plot_right, plot_dorsal, plot_ventral] ;
if plot_kymographs_cumsum
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = find(do_plots)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        HHK = HHsK{qq} ;
        if t0Pathline > 1
            gdotK = cumsum(gdotsK{qq}, 1) - sum(gdotsK{qq}(1:(t0Pathline-1),:), 1) ;
            divvK = cumsum(divvsK{qq}, 1) - sum(divvsK{qq}(1:(t0Pathline-1),:), 1) ;
            velnK = cumsum(velnsK{qq}, 1) - sum(velnsK{qq}(1:(t0Pathline-1),:), 1) ;
            H2vnK = cumsum(H2vnsK{qq}, 1) - sum(H2vnsK{qq}(1:(t0Pathline-1),:), 1) ;
        else
            gdotK = cumsum(gdotsK{qq}, 1) ;
            divvK = cumsum(divvsK{qq}, 1) ;
            velnK = cumsum(velnsK{qq}, 1) ;
            H2vnK = cumsum(H2vnsK{qq}, 1) ;
        end
        m2plot = {gdotK, divvK, velnK, H2vnK} ;
        titles = {'$\int_{-\infty}^t \textrm{d}t \,  \frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'divergence of flow, $\int_{-\infty}^t \textrm{d}t \,  \nabla \cdot \mathbf{v}$', ...
            'normal velocity, $\int_{-\infty}^t \textrm{d}t \,  v_n$', ...
            'normal motion, $\int_{-\infty}^t \textrm{d}t \,  v_n 2 H$'} ;
        labels = {'$\int_{-\infty}^t \textrm{d}t \, \frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]$ ', ...
            '$\int_{-\infty}^t \textrm{d}t \,  \nabla \cdot \mathbf{v}$ ', ...
            'normal velocity, $\int_{-\infty}^t \textrm{d}t \,  v_n$ ', ...
            'normal motion, $\int_{-\infty}^t \textrm{d}t \,  v_n 2 H $ '} ;
        names = {'Igdot', 'Idivv', 'Iveln', 'IH2vn'} ;
        climits = [climit, climit, climit_veln, climit] ;
        climits = climits * 4; 
        
        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{pp} '.png']) ;
            fn_zoom = fullfile(odir, [names{pp} '_zoom.png']) ;
            fn_zoom_early = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([-climits(pp), climits(pp)])
                colormap(bwr256)
                % Add features/folds to plot
                hold on;
                for ii = 1:length(fons)
                    fonsi = max(1, fons(ii)) ;
                    tones = ones(size(tps(fonsi:end))) ;
                    plot(featureIDs(ii) * tones / nU, tps(fonsi:end))
                end
                
                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                tmp = strsplit(fn, filesep) ;
                disp(['saving ', tmp{end}, ': ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on small values
                caxis([-climits(pp)/2, climits(pp)/2])
                colormap(bwr256)
                tmp = strsplit(fn_zoom, filesep) ;
                disp(['saving ', tmp{end}])
                export_fig(fn_zoom, '-png', '-nocrop', '-r200')   
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                caxis([-climits(pp)/2, climits(pp)/2])
                colormap(bwr256)
                tmp = strsplit(fn_zoom_early, filesep) ;
                disp(['saving ', tmp{end}])
                export_fig(fn_zoom_early, '-png', '-nocrop', '-r200')   
            end
        end
    end
end

%% Kymographs of cumulative products along pathlines 
if plot_kymographs_cumprod
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = find(do_plots)
        % Prep the output directory for this averaging
        odir = fullfile(outdirs{qq}, 'cumprods') ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        gdotK_pos = cumprod(1 + gdotsK{qq}(tps>eps, :), 1) ;
        gdotK_neg = flipud(cumprod(flipud(1 ./ (1 + gdotsK{qq}(tps<eps, :))), 1)) ;
        gdotK = cat(1, gdotK_neg, gdotK_pos) ;
        divvK_pos = cumprod(1 + divvsK{qq}(tps>eps, :), 1) ;
        divvK_neg = flipud(cumprod(flipud(1 ./ (1 + divvsK{qq}(tps<eps, :))), 1)) ;
        divvK = cat(1, divvK_neg, divvK_pos) ;
        H2vnK_pos = cumprod(1 + H2vnsK{qq}(tps>eps, :), 1) ;
        H2vnK_neg = flipud(cumprod(flipud(1 ./ (1 + H2vnsK{qq}(tps<eps, :))), 1)) ;
        H2vnK = cat(1, H2vnK_neg, H2vnK_pos) ;
        m2plot = {gdotK, divvK, H2vnK} ;
        titles = {'$\bar{\epsilon} \equiv \Pi_{0}^{t}\, \left(1+\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]\right)$',...
            'divergence of flow, $\Pi_{0}^{t}\,\left[ 1 + \nabla \cdot \mathbf{v}(\tau)\right]$', ...
            'normal motion, $\Pi_{0}^{t}\, \left[ 1 + v_n(\tau) 2 H(\tau)\right]$'} ;
        labels = {'$\bar{\epsilon}$', ...
            '$ \Pi_{0}^{t}\, \left[ 1 + \nabla \cdot \mathbf{v}(\tau) \right]$' , ...
            'normal motion, $\Pi_{0}^{t}\, \left[ 1 + v_n(\tau) 2 H(\tau)\right] $' } ;
        names = {'Pgdot_t0', 'Pdivv_t0', 'PH2vn_t0'} ;
        climits = [climit, climit, climit] ;
        climits = climits * 3; 
        
        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{pp} '.png']) ;
            fn_zoom = fullfile(odir, [names{pp} '_zoom_early.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([1-climits(pp), 1+climits(pp)])
                colormap(bwr256)
                % Add folds to plot
                hold on;
                fons1 = max(1, fons(1)) ;
                fons2 = max(1, fons(2)) ;
                fons3 = max(1, fons(3)) ;
                t1ones = ones(size(tps(fons1:end))) ;
                t2ones = ones(size(tps(fons2:end))) ;
                t3ones = ones(size(tps(fons3:end))) ;
                tidx0 = QS.xp.tIdx(t0) ;
                plot(folds.folds(tidx0, 1) * t1ones / nU, tps(fons1:end))
                plot(folds.folds(tidx0, 2) * t2ones / nU, tps(fons2:end))
                plot(folds.folds(tidx0, 3) * t3ones / nU, tps(fons3:end))

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on small values
                caxis([1-climits(pp)/3, 1+climits(pp)/3])
                colormap(bwr256)
                fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                tmp = strsplit(fn, filesep) ;
                disp(['saving ', tmp{end}])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                % Zoom in on early times
                ylim([min(tps), max(fons) + 10])
                caxis([1-climits(pp)/3, 1+climits(pp)/3])
                colormap(bwr256)
                tmp = strsplit(fn_zoom, filesep) ;
                disp(['saving ', tmp{end}])
                export_fig(fn_zoom, '-png', '-nocrop', '-r200')   
            end
        end
    end
end


%% Metric Kinematic Correlations
% Plot both all time and select times
timeSpans = {tps, tps(tps < max(fons) + 11)} ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation between terms div(v) and 2Hvn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
padN = round(0.1 * nU) ;
cols = padN:nU - padN ;
corrDir = fullfile(mKPDir, 'correlations_DVquadrants_smoothed') ;

if plot_correlations
    if ~exist(corrDir, 'dir')
        mkdir(corrDir)
    end
    for sigma = 0
        outputFileNames = {fullfile(corrDir, ...
            sprintf('correlation_sigma%02d_alltime_div_2Hvn', sigma)), ...
            fullfile(corrDir, ...
            sprintf('correlation_sigma%02d_earlytimes_div_2Hvn', sigma))} ;
        if ~exist([outputFileNames{1} '.png'], 'file') || ...
                ~exist([outputFileNames{2} '.png'], 'file') || overwrite
            alphaVal = 0.6 ;
            sz = 10 ;
            cmap = parula ;
            close all
            set(gcf, 'visible', 'off')
            % Consider each timespan (early or entire series)
            for tspanIdx = 1:2
                fnout = outputFileNames{tspanIdx} ;
                timeSpan_i = timeSpans{tspanIdx} ;
                ntspan = length(timeSpan_i) ;
                titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
                markers = QS.plotting.markers ;
                colors = mapValueToColor(1:ntspan, [1, ntspan], cmap) ;
                close all
                sphCollection = cell(4, 1) ;
                sposCollection = cell(4, 1) ;
                sphCollection2 = cell(4, 1) ;
                sposCollection2 = cell(4, 1) ;
                sphCollection3 = cell(4, 1) ;
                sposCollection3 = cell(4, 1) ;
                for qq = 1:4  % consider left, right, dorsal, ventral
                    disp(['qq = ', num2str(qq), ': ', titles{qq}])
                    divv = divvsK{qq + 1}(:, cols) ;    
                    H2vn = H2vnsK{qq + 1}(:, cols) ;

                    % Optional: smooth here
                    if sigma > 0
                        divv = imgaussfilt(divv, sigma);            
                        H2vn = imgaussfilt(H2vn, sigma);  
                    end

                    % Check the smoothing on kymographs
                    figure(2) ;
                    sphCollection2{qq} = subplot(2, 2, qq) ;
                    imagesc(cols/nU, tps, divv); 
                    caxis([-climit, climit])
                    colormap(bwr256)
                    figure(3) ;
                    sphCollection3{qq} = subplot(2, 2, qq) ;
                    imagesc(cols/nU, tps, H2vn); 
                    caxis([-climit, climit])
                    colormap(bwr256)

                    figure(1) ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(divv(row, :), H2vn(row, :), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        figure(1)
                        xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                                'Interpreter', 'Latex') ;
                        figure(2)
                        xlabel('ap position, $\zeta/L$', ...
                            'Interpreter', 'Latex') ;
                        figure(3)
                        xlabel('ap position, $\zeta/L$', ...
                            'Interpreter', 'Latex') ;
                    end
                    figure(1)
                    axis equal
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    xmin = max(-2 * climit, xlims(1)) ;
                    xmax = min( 2 * climit, xlims(2)) ;
                    ymin = max(-2 * climit, ylims(1)) ;
                    ymax = min( 2 * climit, ylims(2)) ;

                    % Label the y axis if on the left column
                    ylabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
                    title(titles{qq}, 'Interpreter', 'Latex')

                    figure(2)
                    ylabel(['time [' QS.timeUnits, ']'], 'Interpreter', 'Latex') ;
                    title(titles{qq}, 'Interpreter', 'Latex')
                    figure(3)
                    ylabel(['time [' QS.timeUnits, ']'], 'Interpreter', 'Latex') ;
                    title(titles{qq}, 'Interpreter', 'Latex')

                    % Grab axis position
                    sposCollection2{qq} = get(sphCollection2{qq}, 'Position');
                    sposCollection3{qq} = get(sphCollection3{qq}, 'Position');
                end

                % Adjust xlim and ylim for all four panels in Figure 1
                for qq = 1:4
                    axes(sphCollection{qq})
                    % Add dashed y=x line
                    xlim([xmin, xmax])
                    ylim([ymin, ymax])
                    plot(2 * climit * [-1, 1], 2 * climit * [-1, 1], 'k--')
                    sposCollection{qq} = get(sphCollection{qq}, 'Position');
                end

                % Move subplots left a bit for colorbar space
                for qq = 1:length(sphCollection)
                    spos = sposCollection{qq} ;
                    wh = min(spos(3)-0.05, spos(4)) ;
                    if mod(qq, 2) == 1
                        axes(sphCollection{qq})   
                        set(sphCollection{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                        axis equal
                        set(sphCollection2{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                        set(sphCollection3{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                    else
                        axes(sphCollection{qq})   
                        set(sphCollection{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                        axis equal
                        set(sphCollection2{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                        set(sphCollection3{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                    end
                end

                % master titles (suptitles)
                figure(1) ;
                sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ',...
                    '$\sigma=$', num2str(sigma), ' ', QS.timeUnits], ...
                    'Interpreter', 'Latex') ;

                figure(2) ;
                sgtitle(['$\nabla \cdot \bf{v}_\parallel,$ $\sigma=$', ...
                    num2str(sigma), ' ', QS.timeUnits], ...
                        'Interpreter', 'Latex') ;
                figure(3) ;
                sgtitle(['$2Hv_n$ $\sigma=$', ...
                    num2str(sigma), ' ', QS.timeUnits], 'Interpreter', 'Latex') ;

                % Add colorbar
                figure(1) ;
                c = colorbar('Position',[.9 .333 .02 .333]) ;
                % Make colorbar share the alpha of the image
                % Manually flush the event queue and force MATLAB to render the colorbar
                % necessary on some versions
                drawnow
                % Get the color data of the object that correponds to the colorbar
                cdata = c.Face.Texture.CData;
                % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
                cdata(end,:) = uint8(alphaVal * cdata(end,:));
                % Ensure that the display respects the alpha channel
                c.Face.Texture.ColorType = 'truecoloralpha';
                % Update the color data with the new transparency information
                c.Face.Texture.CData = cdata;
                c.Label.Interpreter = 'Latex' ;
                c.Label.String = ['time [' QS.timeUnits ']'] ;
                c.Ticks = [0, 1] ;
                c.TickLabels = [tps(1), max(timeSpan_i)] ;

                figure(2) ;
                c = colorbar('Position',[.9 .333 .02 .333]) ;
                figure(3) ;
                c = colorbar('Position',[.9 .333 .02 .333]) ;

                % Save figure
                figure(1)
                saveas(gcf, [fnout '.png']) ;
                figure(2)
                saveas(gcf, [fnout '_kymo_divv.png']) ;
                figure(3)
                saveas(gcf, [fnout '_kymo_H2vn.png']) ;
                close all
                set(gcf, 'visible', 'off')
            end
        end
        disp('done with correlation plots betweeen divv and H2vn')
    end
end

%% Plot 1D curves for region around each fold
% Sample divv/H2vn/gdot in each lobe
if plot_cumsum_cumprod
    foldw = 0.05 ;
    endw = 0.10 ;

    if length(featureIDs) == 1
        cut = [round(endw*nU), featureIDs(1)-round(foldw*nU), ...
                featureIDs(1)+round(foldw*nU), round((1-endw) * nU)] ;
        lobes = { cut(1):cut(2), cut(3):cut(4) } ;
    elseif length(featureIDs) == 2
        cut = [round(endw*nU), featureIDs(1)-round(foldw*nU), ...
                featureIDs(1)+round(foldw*nU), featureIDs(2)-round(foldw*nU), ...
                featureIDs(2)+round(foldw*nU),  round((1-endw) * nU)] ;
        lobes = { cut(1):cut(2), cut(3):cut(4), cut(5):cut(6) } ;
    elseif length(featureIDs) == 3
        cut = [round(endw*nU), featureIDs(1)-round(foldw*nU), ...
                featureIDs(1)+round(foldw*nU), featureIDs(2)-round(foldw*nU), ...
                featureIDs(2)+round(foldw*nU), featureIDs(3)-round(foldw*nU), ... 
                featureIDs(3)+round(foldw*nU), round((1-endw) * nU)] ;
        lobes = { cut(1):cut(2), cut(3):cut(4), cut(5):cut(6), cut(7):cut(8) } ;
    else
        error('Code for this number of features here')
    end


    avgStrings = {'dv-averaged', 'left side', 'right side', ...
        'dorsal side', 'ventral side'} ;
    do_plots = [plot_ap, plot_left, plot_right, plot_dorsal, plot_ventral] ;
    avgLabel = {'dv', 'left', 'right', 'dorsal', 'ventral'} ;
    titleFoldBase = 'Lagrangian metric kinematics along folds, ' ;
    titleLobeBase = 'Lagrangian metric kinematics along lobes, ' ;

    % todo: modify if different # of features
    foldYlabels = {'anterior fold', 'middle fold', 'posterior fold'} ;
    lobeYlabels = {'lobe 1', 'lobe 2', 'lobe 3', 'lobe 4'} ;
    for qq = find(do_plots)
        divv = divvsK{qq} ;
        H2vn = H2vnsK{qq} ;
        HH = HHsK{qq} ;
        if ~exist(outdirs{qq}, 'dir')
            mkdir(outdirs{qq})
        end

        %% Plot fold Kinematics -- instantaneous data
        % Explore a range of widths for plotting
        for width = 1:round(maxWFrac * nU)
            fn = fullfile(outdirs{qq}, ...
                [sprintf('fold_kinematics_w%03d_', 2*width+1), ...
                avgLabel{qq}, '.png']) ;
            fn_withH = fullfile(outdirs{qq}, ...
                [sprintf('fold_kinematics_w%03d_', 2*width+1), ...
                avgLabel{qq}, '_withH.png']) ;

            aux_plotPathlineMetricKinematicsFolds_subpanels(QS, ...
                fn, fn_withH, ...
                featureIDs, width, nU, tps, divv, H2vn, HH, titleFoldBase, ...
                foldYlabels, avgStrings{qq}, divvcolor, H2vncolor, ...
                Hposcolor, Hnegcolor, Hsz, overwrite)
        end

        %% Plot lobe Kinematics -- instantaneous data
        fn = fullfile(outdirs{qq}, ...
            ['lobe_kinematics_' avgLabel{qq} '.png']) ;
        fn_withH = fullfile(outdirs{qq}, ...
            ['lobe_kinematics_' avgLabel{qq} '_withH.png']) ;
        aux_plotPathlineMetricKinematicsLobes_subpanels(QS, fn, fn_withH, ...
            lobes, tps, divv, H2vn, HH, lobeYlabels, avgStrings{qq}, ...
            titleLobeBase, divvcolor, H2vncolor, ...
            Hposcolor, Hnegcolor, Hsz, overwrite)

        %% Fold kinematics -- cumprod gdot
        cumsum_cumprod = {'cumsum', 'cumprod'} ;

        for spij = 1:2
            sumprod = cumsum_cumprod{spij} ;
            odir = fullfile(outdirs{qq}, [sumprod 's']) ;
            if ~exist(odir, 'dir')
                mkdir(odir)
            end
            for width = 1:round(maxWFrac * nU)
                %% Plot kinematics on separate axis for each fold
                % ms2plot = {'gdot', 'divv', 'H2vn'} ;    
                fn = fullfile(odir, ...
                    ['fold_kinematics_', sumprod, ...
                    sprintf('_w%03d_', 2*width+1), ...
                    avgLabel{qq}, '.png']) ;
                fn_withH = fullfile(odir, ...
                    ['fold_kinematics_', sumprod, ...
                    sprintf('w%03d_', 2*width+1), ...
                    avgLabel{qq}, '_withH.png']) ;
                if plot_fold_kinematics
                    aux_plotPathlineMetricKinematicsFolds_integrated_subpanels(QS, ...
                        fn, featureIDs, width, nU, tps, divv, H2vn, titleFoldBase, ...
                        foldYlabels, avgStrings{qq}, ...
                        divvcolor, H2vncolor, gdotcolor, overwrite, sumprod)
                end

                %% Plot all fold gdots on one axis
                ms2plot = {'gdot', 'divv', 'H2vn'} ;
                for jj = 1:length(ms2plot)
                    m2plot = ms2plot{jj} ; 

                    fn = fullfile(odir, ...
                        ['fold_kinematics_', m2plot, '_', sumprod, ...
                        sprintf('_compare_w%03d_', 2*width+1), ...
                        avgLabel{qq}, '.png']) ;

                    fn_withH = fullfile(odir, ...
                        ['fold_kinematics_', m2plot, '_', sumprod, ...
                        sprintf('_compare_w%03d_', 2*width+1), ...
                        avgLabel{qq}, '_withH.png']) ;
                    aux_plotPathlineMetricKinematicsFolds_integrated(QS, m2plot, fn, fn_withH, ...
                        featureIDs, width, nU, tps, divv, H2vn, HH, foldYlabels, ...
                        avgStrings{qq}, Hsz, overwrite, sumprod) 
                end
            end

            %% Lobe kinematics -- cumprod gdot
            if plot_lobe_kinematics
                fn = fullfile(outdirs{qq}, ...
                    ['lobe_kinematics_' sumprod '.png']) ;
                fn_withH = fullfile(outdirs{qq}, ...
                    ['lobe_kinematics_' sumprod '_withH.png']) ;
                aux_plotPathlineMetricKinematicsLobes_integrated_subplots(QS, ...
                    fn, fn_withH, ...
                    lobes, tps, divv, H2vn, lobeYlabels, avgStrings{qq}, ...
                    titleLobeBase, overwrite, sumprod)
            end

            %% Plot all lobes' kinematic quantitiy (gdot/H2vn/divv) on one axis
            if  plot_lobe_kinematics
                ms2plot = {'gdot', 'divv', 'H2vn'} ;
                for jj = 1:length(ms2plot)
                    m2plot = ms2plot{jj} ; 

                    fn = fullfile(odir, ...
                        ['lobe_kinematics_', m2plot, '_', ...
                        sumprod, '_compare.png']) ;
                    fn_withH = fullfile(odir, ...
                        ['lobe_kinematics_' m2plot '_', ...
                        sumprod, '_compare_withH.png']) ;
                    aux_plotPathlineMetricKinematicsLobes_integrated(QS, m2plot, fn, fn_withH, ...
                        lobes, tps, divv, H2vn, HH, lobeYlabels, ...
                        avgStrings{qq}, Hsz, overwrite, sumprod)
                end
            end
        end
    end
end

disp('done')



