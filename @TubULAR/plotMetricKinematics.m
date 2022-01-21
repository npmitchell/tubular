function plotMetricKinematics(QS, options)
% plotMetricKinematics(QS, options)
%   Plot the metric Kinematics as kymographs and correlation plots
%   Out-of-plane motion is v_n * 2H, where v_n is normal velocity and H is
%   mean curvature.
%   In-plane motion considered here is div(v_t) where v_t is tangential
%   velocity on the curved surface.
%   The difference div(v_t) - vn*2H = 1/2 Tr[g^{-1} dot{g}], which is a 
%   measure of isotropic metric change over time (dot is the partial
%   derivative wrt time). 
% 
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs : bool
%   plot_kymographs_cumsum : bool
%   plot_correlations : bool
%   plot_gdot_correlations : bool
%   plot_gdot_decomp : bool
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
overwrite_timePoints = false ;
plot_kymographs = true ;
plot_correlations = true ;
plot_gdot_correlations = false ;
plot_raw_scatter_correlations = false ;
plot_spaceMaps = false ;
plot_gdot_decomp = true ;
plot_flows = true ;
plot_factors = true ;
plot_Hgdot = false ;

%% Parameter options
lambda = QS.smoothing.lambda ; 
lambda_err = QS.smoothing.lambda_err ;
lambda_mesh = QS.smoothing.lambda_mesh ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;
climit = 0.2 ;
climit_err = 0.2 ;
climit_veln = climit * 10 ;
climit_H = climit * 2 ;
climit_radius = 0 ;
displayTimepointInRange = 0.25 ;
figResolutionStr = '-r600' ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwrite_timePoints')
    overwrite_timePoints = options.overwrite_timePoints ;
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
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'climit_err')
    climit_err = options.climit_err ;
end
if isfield(options, 'climit_veln')
    climit_veln = options.climit_veln ;
end
if isfield(options, 'climit_H')
    climit_H = options.climit_H ;
end
if isfield(options, 'climit_radius')
    climit_radius = options.climit_radius ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'plot_flows')
    plot_flows = options.plot_flows ;
end
if isfield(options, 'plot_factors')
    plot_factors = options.plot_factors ;
end
if isfield(options, 'plot_Hgdot')
    plot_Hgdot = options.plot_Hgdot ;
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
if isfield(options, 'plot_correlations')
    plot_correlations = options.plot_correlations ;
end
if isfield(options, 'plot_raw_correlations')
    plot_raw_correlations = options.plot_raw_correlations ;
else
    plot_raw_correlations = plot_correlations ;
end
if isfield(options, 'plot_raw_scatter_correlations')
    plot_raw_scatter_correlations = options.plot_raw_scatter_correlations ;
end
if isfield(options, 'plot_spaceMaps')
    plot_spaceMaps = options.plot_spaceMaps ;
end
if isfield(options, 'plot_gdot_correlations')
    plot_gdot_correlations = options.plot_gdot_correlations ;
end
if isfield(options, 'plot_gdot_decomp')
    plot_gdot_decomp = options.plot_gdot_decomp ;
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
folds = load(QS.fileName.fold) ;

%% Colormap
close all
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;

%% Load time offset for first fold, t0
QS.t0set() ;
tfold = QS.t0 ;
% timepoints at which features/folds begin
fons = folds.fold_onset - tfold ;
% timepoint indices at which features/folds begin
fonsIdx = folds.fold_onset - QS.xp.fileMeta.timePoints(1) + 1 ;

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
tps = (QS.xp.fileMeta.timePoints(1:end-1) - tfold) * QS.timeInterval ;
vtimePoints = QS.xp.fileMeta.timePoints(1:end-1) ;

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
outdir = fullfile(mKDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];
Hunitstr = [ '[1/' QS.spaceUnits ']' ];
vunitstr = [ '[' QS.spaceUnits '/' QS.timeUnits ']' ];
runitstr = [ '[' QS.spaceUnits ']' ] ;
    
% Compute or load all timepoints
for tp = QS.xp.fileMeta.timePoints(1:end-1)
    close all
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Check for timepoint measurement on disk
    Hfn = fullfile(outdir, sprintf('HH_vertices_%06d.mat', tp))   ;
    efn = fullfile(outdir, sprintf('gdot_vertices_%06d.mat', tp)) ;
    dfn = fullfile(outdir, sprintf('divv_vertices_%06d.mat', tp)) ;
    nfn = fullfile(outdir, sprintf('veln_vertices_%06d.mat', tp)) ;
    H2vnfn = fullfile(outdir, sprintf('H2vn_vertices_%06d.mat', tp)) ;
    rfn = fullfile(outdir, sprintf('radius_vertices_%06d.mat', tp)) ;

    % Load timeseries measurements
    load(Hfn, 'HH_filt', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
    load(efn, 'gdot_filt', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
    load(dfn, 'divv_filt', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
    load(nfn, 'veln_filt', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
    load(H2vnfn, 'H2vn_filt', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 
    load(rfn, 'radius_filt', 'radius_ap', 'radius_l', 'radius_r', 'radius_d', 'radius_v') 

    % separate 2d/3d data
    H2d = HH_filt ;
    H3d = HH_filt(:, 1:nV-1) ;
    gdot2d = gdot_filt ;
    gdot3d = gdot_filt(:, 1:nV-1);
    divv2d = divv_filt ;
    divv3d = divv_filt(:, 1:nV-1) ;
    veln2d = veln_filt ;
    veln3d = veln_filt(:, 1:nV-1) ;
    H2vn2d = H2vn_filt ;
    H2vn3d = H2vn_filt(:, 1:nV-1) ;
    radi2d = radius_filt ;
    radi3d = radius_filt(:, 1:nV-1) ;
    
    %% Plot results
    % % operational plotting options
    pOptions.overwrite = overwrite_timePoints ;
    pOptions.plot_flows = plot_flows ;
    pOptions.plot_Hgdot = plot_Hgdot ;
    pOptions.plot_factors = plot_factors ;
    % parameter plotting options
    pOptions.doubleResolution = doubleResolution; 
    pOptions.lambda = lambda ;
    pOptions.lambda_err = lambda_err ;
    pOptions.lambda_mesh = lambda_mesh ;
    pOptions.nmodes = nmodes ;
    pOptions.zwidth = zwidth ;
    pOptions.H2vn2d = H2vn2d ;
    pOptions.divv2d = divv2d ;
    pOptions.gdot2d = gdot2d ;
    pOptions.veln2d = veln2d ;
    pOptions.radi2d = radi2d ;
    pOptions.H2d = H2d ;
    pOptions.H2vn3d = H2vn3d ;
    pOptions.divv3d = divv3d ;
    pOptions.gdot3d = gdot3d ;
    pOptions.veln3d = veln3d ;
    pOptions.radi3d = radi3d ;
    pOptions.H3d = H3d ;
    pOptions.cutMesh = [] ;
    pOptions.mesh = [] ;
    pOptions.climit = climit ;
    pOptions.climit_err = climit ;
    pOptions.climit_veln = climit_veln ;
    pOptions.climit_H = climit_H ;
    QS.plotMetricKinematicsTimePoint(tp, pOptions)
    
    %% Store in matrices
    HH_M(tidx, :, :) = HH_filt ;
    gdot_M(tidx, :, :) = gdot_filt ;
    divv_M(tidx, :, :) = divv_filt ;
    veln_M(tidx, :, :) = veln_filt ;
    H2vn_M(tidx, :, :) = H2vn_filt ;
    
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

%% Store kymograph data in cell arrays
HHsK = {HH_apM, HH_lM, HH_rM, HH_dM, HH_vM} ;
gdotsK = {gdot_apM, gdot_lM, gdot_rM, gdot_dM, gdot_vM} ;
divvsK = {divv_apM, divv_lM, divv_rM, divv_dM, divv_vM} ;
velnsK = {veln_apM, veln_lM, veln_rM, veln_dM, veln_vM} ;
H2vnsK = {H2vn_apM, H2vn_lM, H2vn_rM, H2vn_dM, H2vn_vM} ;
radisK = {radi_apM, radi_lM, radi_rM, radi_dM, radi_vM} ;

%% Now plot different measured quantities as kymographs
if plot_kymographs
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
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
            ['radius ' runitstr]} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn', 'radius'} ;
        climits = [climit, climit_H, climit, climit_veln, climit, climit_radius] ;
        cmaps = {bwr256, bwr256, bwr256, bwr256, bwr256, 'parula'};
            
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
                colormap(cmaps{pp})
                % Add folds to plot
                hold on;
                for ii = 1:length(fonsIdx)
                    fonsi = max(1, fonsIdx(ii) ) ; 
                    plot(double(folds.folds(fonsi:end-1, ii)) / double(nU), ...
                        tps(fonsi:end))            
                end

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                fn = fullfile(odir, [ names{pp} '.png']) ;
                disp(['saving kymograph ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   

                % Zoom in on small values
                if climits(pp) > 0
                    caxis([-climits(pp)/3, climits(pp)/3])
                    fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                    disp(['saving kymograph detail ', fn])
                    export_fig(fn, '-png', '-nocrop', '-r200')   
                    % Zoom in on early times
                    ylim([min(tps), min(max(tps), max(fons) + 10)])
                    caxis([-climits(pp)/3, climits(pp)/3])
                    fn = fullfile(odir, [names{pp} '_zoom_early.png']) ;
                    disp(['saving kymograph detail: ', fn])
                    export_fig(fn, '-png', '-nocrop', '-r200')  
                end
            end
        end
    end
end

%% Metric Kinematic Correlations -- vertex-by-vertex
% Plot both all time and select times
timeSpans = {tps} ;
% define pre-folding time chunk(s)
DTime = 30 ;    % 30 minute chunks
if tps(1) < 0
    nchunks = ceil(abs(tps(1)/DTime)) ;
    for qq = 1:nchunks
        tspanStart = max(min(tps), -(nchunks-qq+1)*DTime) ;
        tspanEnd = -(nchunks-qq)*DTime ;
        timeSpans{1+qq} = tspanStart:QS.timeInterval:tspanEnd ;
    end
else
    qq = 0 ;
end
% Now add time chunks after initial folding t0
if max(tps) > 0
    nchunks = ceil(max(tps)/DTime) ;
    for pp = 1:nchunks
        tspanStart = (pp-1)*DTime ;
        tspanEnd = min(max(tps), pp*DTime) ;
        timeSpans{1+qq+pp} = tspanStart:QS.timeInterval:tspanEnd ;
    end
end
% Include one tpsan of -1hr, 1hr, 2hr after t0
hrTSpans = {max(min(tps),-60):QS.timeInterval:0, 0:QS.timeInterval:60, 60:QS.timeInterval:min(max(tps),120)} ;
timeSpans{end+1} = hrTSpans{1} ;
timeSpans{end+1} = hrTSpans{2} ;
timeSpans{end+1} = hrTSpans{3} ;
% Include one tpsan of -1.5hr, 1.5hr, 3hr after t0
lhrTSpans = {max(min(tps),-90):QS.timeInterval:0, 0:QS.timeInterval:90, ...
    90:QS.timeInterval:min(max(tps), 180), ...
    max(min(tps), 0):QS.timeInterval:min(max(tps),75)} ;
% Include one tpsan of -30m:30m, 30m:90min, 90min:150min wrt t0
lhrTSpans = {max(min(tps),-30):QS.timeInterval:min(max(tps), 30), ...
    min(max(tps), 30):QS.timeInterval:90, ...
    max(min(tps), 90):QS.timeInterval:min(max(tps), 150), ...
    0:QS.timeInterval:90} ;
timeSpans{end+1} = lhrTSpans{1} ;
timeSpans{end+1} = lhrTSpans{2} ;
timeSpans{end+1} = lhrTSpans{3} ;
timeSpans{end+1} = lhrTSpans{4} ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation between terms div(v) and 2Hvn for each point (no DV
% averaging).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
padN = round(0.1 * nU) ;
cols = padN:nU - padN ;
if strcmp(samplingResolution, '1x')
    q0 = round(nV * 0.125) ;
    q1 = round(nV * 0.375) ;
    q2 = round(nV * 0.625) ;
    q3 = round(nV * 0.875) ;
end
left = q0:q1 ;
ventral = q1:q2 ;
right = q2:q3 ;
dorsal = [q3:nV, 1:q1] ;
rows = {left, right, dorsal, ventral} ;
    
%% SPATIAL VERSION -- ON MESH, SIMPLE-AVERAGED IN TIME
% Predefine file names
mapDir = fullfile(mKDir, 'spaceMap') ;
outputFileNames = {fullfile(mapDir, 'correlation_alltime_div_2Hvn')} ;
files_exist = true ;
for qq = 2:length(timeSpans)
    outputFileNames{qq} = ...
        fullfile(mapDir, ...
        sprintf('spaceMapKinematics_times%02d_div_2Hvn', qq)) ;
    files_exist = files_exist && exist([outputFileNames{qq}, '_zoom.pdf'], 'file') ;
end
if plot_spaceMaps && (~files_exist || overwrite)
    if ~exist(mapDir, 'dir')
        mkdir(mapDir)
    end
    close all
    set(gcf, 'visible', 'off')
    % Consider each timespan (early or entire series)
    for tspanIdx = 2:length(timeSpans)
        fnout = outputFileNames{tspanIdx} ;
        timeSpan_i = timeSpans{tspanIdx} ;
        
        minTP = max(min(timeSpan_i(1)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        maxTP = max(min(timeSpan_i(end)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        minTidx = QS.xp.tIdx(minTP) ;
        maxTidx = QS.xp.tIdx(maxTP) ;
        tidx_i = QS.xp.tIdx(minTP):QS.xp.tIdx(maxTP) ;
        displayTidx = minTidx + round((maxTidx - minTidx) * displayTimepointInRange) ;
        tidx2plot = min(displayTidx, length(QS.xp.fileMeta.timePoints)) ;
        close all
        
        % Get middle mesh for this time range if not ALLTime
        if tspanIdx == 1
            QS.setTime(QS.t0set()) ;
            mesh = QS.getCurrentSPCutMeshSmRS ;
        else
            % Match displayTimepointRange in range
            QS.setTime(QS.xp.fileMeta.timePoints(tidx2plot)) ;
            mesh = QS.getCurrentSPCutMeshSmRS ;
        end
        m2d = mesh ;
        m2d.v = m2d.u ;
        m2d.v(:, 1) = m2d.v(:, 1) / max(m2d.v(:, 1)) ;
        
        % Get mean divv and mean H2vn for these times
        mdivv = squeeze(mean(divv_M(tidx_i, :, :), 1)) ;
        mH2vn = squeeze(mean(H2vn_M(tidx_i, :, :), 1)) ;
        halfgdot = mdivv - mH2vn ;
        pOptions.labels = {'$2Hv_n$', '$\nabla \cdot v_t$', ...
            '$\frac{1}{2}\mathrm{Tr}\left[\mathbf{g}^{-1} \dot{\mathbf{g}}\right]=\nabla \cdot v_t-2Hv_n$', ...
            '', '', ''} ;
        pOptions.cbarlabels = {'', '', '', ...
            '$2Hv_n$', '$\nabla \cdot v_t$', ...
            '$\frac{1}{2}\mathrm{Tr}\left[\mathbf{g}^{-1} \dot{\mathbf{g}}\right]$'} ;
        pOptions.makeCbar = [false, false, false, true, true, true] ;
        % pOptions.cmap = bwr ;
        % pOptions.cmap = twilight_shifted_trueWhite(256) ;
        pOptions.cmap = twilight_shifted_mod(256) ;
        pOptions.axisOn = false ;
        pOptions.clim = [-climit, climit] ;
        pOptions.visible = 'On' ;
        pOptions.view = {[0,0], [0,0], [0,0], [0,90], [0,90], [0,90]} ;
        pOptions.xlim = {xyzlim(1, :), xyzlim(1, :), xyzlim(1, :), ...
            [0,1], [0,1], [0,1]};
        pOptions.zlim = {xyzlim(3, :), xyzlim(3, :), xyzlim(3, :), ...
            [0,1], [0,1], [0,1]};
        [axs, cbs, meshHandles]  = ...
            nFieldsOnSurface({mesh, mesh, mesh, m2d, m2d, m2d}, ...
            {mH2vn(:), mdivv(:), halfgdot(:), ...
            mH2vn(:), mdivv(:), halfgdot(:)}, pOptions) ;
        
        % Enlarge lower 2d axes to fit 3d embedding axis size
        pos = {} ;
        for ii = 1:3
            set(gcf, 'CurrentAxes', axs{ii})
            pos{ii} = get(gca, 'position') ;
        end
        for ii = 4:length(axs)
            set(gcf, 'CurrentAxes', axs{ii})
            p2 = get(gca, 'position') ;
            set(gca, 'position', [pos{ii-3}(1), p2(2)-pos{ii-3}(4)*0.1, pos{ii-3}(3), pos{ii-3}(4)])
        end
        
        % Correlation for 10%-90% AP axis
        lowU = 0.1 * max(m2d.u(:, 1)) ;
        hiU = 0.9 * max(m2d.u(:, 1)) ;
        selectSpace = (m2d.u(:, 1) > lowU) & (m2d.u(:, 1) < hiU) ;
        RR = corrcoef(mH2vn(selectSpace), mdivv(selectSpace)) ;
        rho = RR(1, 2) ;
        corrString = [' $\rho=$' sprintf('%0.3f', rho)] ;
        
        % master titles (suptitles)
        if tspanIdx > 1
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ', ...
                num2str(min(timeSpan_i)) '$<t<$' num2str(max(timeSpan_i)) ...
                ' ' QS.timeUnits corrString], ...
                'Interpreter', 'Latex') ;
        else
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$' corrString], ...
                'Interpreter', 'Latex') ;
        end
        
        % Save figure
        disp(['Saving spaceMap: ' fnout])
        print(gcf, [fnout '.png'], '-dpng', figResolutionStr) ;
        saveas(gcf, [fnout '.pdf']) ;
        
        % Zoom in on color
        for ii = 1:length(axs)
            set(gcf, 'CurrentAxes', axs{ii})
            caxis([-climit / 2, climit / 2])
        end
        % Save figure
        print(gcf, [fnout '_zoom.png'], '-dpng', figResolutionStr) ;
        saveas(gcf, [fnout '_zoom.pdf']) ;
        
        % Zoom out on color
        for ii = 1:length(axs)
            set(gcf, 'CurrentAxes', axs{ii})
            caxis([-climit * 2, climit * 2])
        end
        % Save figure
        print(gcf, [fnout '_out.png'], '-dpng', figResolutionStr) ;
        saveas(gcf, [fnout '_out.pdf']) ;
        
        close all
        clearvars sphCollection
        set(gcf, 'visible', 'off')
    end
    disp('done with spaceMaps of divv, H2vn, and difference')
end

%% HISTCOUNTS VERSION -- vertex-by-vertex
corrDir = fullfile(mKDir, 'correlations') ;
outputFileNames = {fullfile(corrDir, 'histc_correlation_alltime_div_2Hvn.png')} ;
files_exist = exist(outputFileNames{1}, 'file') ;
for qq = 2:length(timeSpans)
    outputFileNames{qq} = ...
        fullfile(corrDir, sprintf('histc_correlation_times%02d_div_2Hvn.png', qq)) ;
    files_exist = files_exist && exist(outputFileNames{qq}, 'file') ;
end 
if plot_raw_correlations && (~files_exist || overwrite)
    xedges = linspace(-2 * climit, 2 * climit, 50) ;
    yedges = xedges ;
    if ~exist(corrDir, 'dir')
        mkdir(corrDir)
    end
    cmap = parula ;
    close all
    set(gcf, 'visible', 'off')
    % Consider each timespan (early or entire series)
    
    for tspanIdx = 1:length(timeSpans)
        fnout = outputFileNames{tspanIdx} ;
        timeSpan_i = timeSpans{tspanIdx} ;
        minTP = max(min(timeSpan_i(1)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        maxTP = max(min(timeSpan_i(end)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        tidx_i = QS.xp.tIdx(minTP):QS.xp.tIdx(maxTP) ;
        
        % ntspan = length(timeSpan_i) ;
        titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
        close all
        sphCollection = cell(4, 1) ;
        sposCollection = cell(4, 1) ;
        nnmax = 0 ;
        for qq = 1:4  % consider left, right, dorsal, ventral
            disp(['qq = ', num2str(qq), ': ', titles{qq}])
            divv = divv_M(:, rows{qq}, cols) ;    
            H2vn = H2vn_M(:, rows{qq}, cols) ;
            
            sphCollection{qq} = subplot(2, 2, qq) ;
            divv_tp = divv(tidx_i, :, :) ;
            H2vn_tp = H2vn(tidx_i, :, :) ;
            nn = histcounts2(divv_tp(:), H2vn_tp(:), ...
                xedges, yedges) ;
            nnmax = max(nnmax, max(nn(:))) ;
        end
        nnmax = round(nnmax /1000) * 1000 ; 
        
        for qq = 1:4  % consider left, right, dorsal, ventral
            disp(['qq = ', num2str(qq), ': ', titles{qq}])
            divv = divv_M(:, rows{qq}, cols) ;    
            H2vn = H2vn_M(:, rows{qq}, cols) ;
            
            sphCollection{qq} = subplot(2, 2, qq) ;
            divv_tp = divv(tidx_i, :, :) ;
            H2vn_tp = H2vn(tidx_i, :, :) ;
            nn = histcounts2(divv_tp(:), H2vn_tp(:), ...
                xedges, yedges) ;
            imagesc(xedges, yedges, log10(nn)') ;
            caxis([0, log10(nnmax)])
            set(gca,'YDir','normal')
            hold on ;

            % Label the x axis if on the bottom row
            if qq > 2
                xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                        'Interpreter', 'Latex') ;
            end
            axis equal

            % xlim([max(-2 * climit, xlims(1)), min(2 * climit, xlims(2))])
            % ylim([max(-2 * climit, ylims(1)), min(2 * climit, ylims(2))])
            % xlims = get(gca, 'xlim') ;
            % ylims = get(gca, 'ylim') ;
            % leftdot = max(xlims(1), ylims(1)) ;
            % rightdot = min(xlims(2), ylims(2)) ;
            % plot([leftdot, rightdot], [leftdot, rightdot], 'k--')

            % Label the y axis if on the left column
            ylabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
            title(titles{qq}, 'Interpreter', 'Latex')

            % Grab axis position
            sposCollection{qq} = get(sphCollection{qq}, 'Position');
        end

        % Adjust xlim and ylim for all four panels in Figure 1
        for qq = 1:4
            axes(sphCollection{qq})
            % Add dashed y=x line
            xlim([-2*climit, 2*climit])
            ylim([-2*climit, 2*climit])
            plot(2 * climit * [-1, 1], 2 * climit * [-1, 1], 'w--')
            sposCollection{qq} = get(sphCollection{qq}, 'Position');
        end
        
        % Move subplots left a bit for colorbar space
        for qq = 1:length(sphCollection)
            spos = sposCollection{qq} ;
            if ~exist('wh', 'var')
                wh = min(spos(3)-0.05, spos(4)) ;
            end
            if qq == 1
                yrow1 = spos(2) ;
            elseif qq == 2
                yrow1 = min(yrow1, spos(2)) ;
            elseif qq == 3
                yrow2 = spos(2) ;
            elseif qq == 4
                yrow2 = min(yrow2, spos(2)) ;
            end
        end
        for qq = 1:length(sphCollection)
            spos = sposCollection{qq} ;
            if qq < 3
                yrow = yrow1 ;
            else
                yrow = yrow2 ;
            end
            if mod(qq, 2) == 1
                set(sphCollection{qq}, 'Position', [spos(1)-0.01, yrow, wh, wh])    
            else
                set(sphCollection{qq}, 'Position', [spos(1)-0.06, yrow, wh, wh])    
            end
        end

        % master titles (suptitles)
        if tspanIdx > 1
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ', ...
                num2str(min(timeSpan_i)) '$<t<$' num2str(max(timeSpan_i)) ...
                ' ' QS.timeUnits], ...
                'Interpreter', 'Latex') ;
        else
            sgtitle('$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$', ...
                'Interpreter', 'Latex') ;
        end
        
        % Add colorbar
        c = colorbar('Position', [.85 .333 .02 .333]) ;
        % % Make colorbar share the alpha of the image
        % % Manually flush the event queue and force MATLAB to render the colorbar
        % % necessary on some versions
        % drawnow
        % % Get the color data of the object that correponds to the colorbar
        % cdata = c.Face.Texture.CData;
        % % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
        % cdata(end,:) = uint8(alphaVal * cdata(end,:));
        % % Ensure that the display respects the alpha channel
        % c.Face.Texture.ColorType = 'truecoloralpha';
        % % Update the color data with the new transparency information
        % c.Face.Texture.CData = cdata;
        c.Label.Interpreter = 'Latex' ;
        c.Label.String = 'counts' ;
        c.Ticks = linspace(1, round(log10(nnmax)), round(log10(nnmax)));
        c.TickLabels = logspace(1, round(log10(nnmax)), round(log10(nnmax))) ;
        
        % Save figure
        disp(['saving histcount correlations: ' fnout])
        saveas(gcf, fnout) ;
        close all
        clearvars sphCollection
        set(gcf, 'visible', 'off')
    end
end

%% SCATTERPLOT VERSION
% Predefine file names
corrDir = fullfile(mKDir, 'correlations') ;
outputFileNames = {fullfile(corrDir, 'correlation_alltime_div_2Hvn.png')} ;
files_exist = exist(outputFileNames{1}, 'file')  ;
for qq = 2:length(timeSpans)
    outputFileNames{qq} = ...
        fullfile(corrDir, sprintf('correlation_times%02d_div_2Hvn.png', qq)) ;
    files_exist = files_exist && exist(outputFileNames{qq}, 'file') ;
end

if plot_raw_scatter_correlations && (~files_exist || overwrite)
    if ~exist(corrDir, 'dir')
        mkdir(corrDir)
    end
    alphaVal = 0.05 ;
    sz = 2 ;
    cmap = phasemap(256) ;
    close all
    set(gcf, 'visible', 'off')
    % Consider each timespan (early or entire series)
    for tspanIdx = 1:length(timeSpans)
        fnout = outputFileNames{tspanIdx} ;
        timeSpan_i = timeSpans{tspanIdx} ;
        
        minTP = max(min(timeSpan_i(1)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        maxTP = max(min(timeSpan_i(end)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        tidx_i = QS.xp.tIdx(minTP):QS.xp.tIdx(maxTP) ;
        
        ntspan = length(timeSpan_i) ;
        titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
        markers = QS.plotting.markers ;
        colors = mapValueToColor(1:ntspan, [1, ntspan], cmap) ;
        close all
        sphCollection = cell(4, 1) ;
        sposCollection = cell(4, 1) ;
        for qq = 1:4  % consider left, right, dorsal, ventral
            disp(['qq = ', num2str(qq), ': ', titles{qq}])
            divv = divv_M(:, rows{qq}, cols) ;    
            H2vn = H2vn_M(:, rows{qq}, cols) ;
            
            sphCollection{qq} = subplot(2, 2, qq) ;
            for row = tidx_i
                disp(['row = ', num2str(row)])
                divv_tp = divv(row, :, :) ;
                H2vn_tp = H2vn(row, :, :) ;
                set(gcf, 'visible', 'off')
                scatter(divv_tp(:), H2vn_tp(:), sz, ...
                    markers{qq}, 'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', colors(row-tidx_i(1)+1, :), ...
                    'MarkerEdgeAlpha', alphaVal) ;
                hold on ;
            end

            % Label the x axis if on the bottom row
            if qq > 2
                xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                        'Interpreter', 'Latex') ;
            end
            axis equal
            % Add dashed y=x line
            xlims = get(gca, 'xlim') ;
            ylims = get(gca, 'ylim') ;
            
            % Get xmin, xmax, ymin, ymax for the ALL TIMES plot
            if tspanIdx == 1
                xmin = max(-2 * climit, xlims(1)) ;
                xmax = min( 2 * climit, xlims(2)) ;
                ymin = max(-2 * climit, ylims(1)) ;
                ymax = min( 2 * climit, ylims(2)) ;
            end

            % xlim([max(-2 * climit, xlims(1)), min(2 * climit, xlims(2))])
            % ylim([max(-2 * climit, ylims(1)), min(2 * climit, ylims(2))])
            % xlims = get(gca, 'xlim') ;
            % ylims = get(gca, 'ylim') ;
            % leftdot = max(xlims(1), ylims(1)) ;
            % rightdot = min(xlims(2), ylims(2)) ;
            % plot([leftdot, rightdot], [leftdot, rightdot], 'k--')

            % Label the y axis if on the left column
            ylabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;

            % Grab axis position
            sposCollection{qq} = get(sphCollection{qq}, 'Position');
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
            if ~exist('wh', 'var')
                wh = min(spos(3)-0.05, spos(4)) ;
            end
            if qq == 1
                yrow1 = spos(2) ;
            elseif qq == 2
                yrow1 = min(yrow1, spos(2)) ;
            elseif qq == 3
                yrow2 = spos(2) ;
            elseif qq == 4
                yrow2 = min(yrow2, spos(2)) ;
            end
        end
        for qq = 1:length(sphCollection)
            spos = sposCollection{qq} ;
            if qq < 3
                yrow = yrow1 ;
            else
                yrow = yrow2 ;
            end
            if mod(qq, 2) == 1
                set(sphCollection{qq}, 'Position', [spos(1)-0.01, yrow, wh, wh])    
            else
                set(sphCollection{qq}, 'Position', [spos(1)-0.06, yrow, wh, wh])    
            end
        end

        % master titles (suptitles)
        if tspanIdx > 1
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ', ...
                num2str(min(timeSpan_i)) '$<t<$' num2str(max(timeSpan_i)) ...
                ' ' QS.timeUnits], ...
                'Interpreter', 'Latex') ;
        else
            sgtitle('$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$', ...
                'Interpreter', 'Latex') ;
        end

        % Add colorbar
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
        c.TickLabels = [min(timeSpan_i), max(timeSpan_i)] ;

        % Save figure
        disp(['saving correlation plot' fnout])
        saveas(gcf, fnout) ;
        close all
        clearvars sphCollection
        set(gcf, 'visible', 'off')
        
    end
    disp('done with correlation plots betweeen divv and H2vn')
end

%% Metric Kinematic Correlations for each hour, one single axis only (all quadrants)
collateTSpans = {0:min(max(tps), 120)} ;
for qq = 1:length(hrTSpans)
    collateTSpans{qq+1} = hrTSpans{qq} ;
end
for qq = 1:length(lhrTSpans)
    collateTSpans{qq+1} = lhrTSpans{qq} ;
end
padN = round(0.1 * nU) ;
cols = padN:nU - padN ;
corrDVqDir = fullfile(mKDir, 'correlations_DVquadrants_collapsed') ;
files_exist = true ;
for qq = 1:length(collateTSpans)
    outputFileNames{qq} = ...
        fullfile(corrDVqDir, sprintf('correlation_allQuad_times%02d_div_2Hvn', qq)) ;
    files_exist = files_exist && ...
        exist([outputFileNames{qq} '.png'], 'file') && ...
        exist([outputFileNames{qq} '_data.mat'], 'file');
end

if plot_correlations && (~files_exist || overwrite)
    if ~exist(corrDVqDir, 'dir')
        mkdir(corrDVqDir)
    end    
    alphaVal = 0.2 ;
    sz = 4 ;
    cmap = phasemap(256) ;
    close all
    set(gcf, 'visible', 'off')
    % Consider each hr long timespan (early or entire series)
    for tspanIdx = 1:length(collateTSpans)
        fnout = outputFileNames{tspanIdx} ;
        timeSpan_i = collateTSpans{tspanIdx} ;
        
        minTP = max(min(timeSpan_i(1)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        maxTP = max(min(timeSpan_i(end)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        tidx_i = QS.xp.tIdx(minTP):QS.xp.tIdx(maxTP) ;
        
        markers = QS.plotting.markers ;
        colors = mapValueToColor(1:5, [1, 5], cmap) ;
        close all
        allX = [] ;
        allY = [] ;
        for qq = 1:4  % consider left, right, dorsal, ventral
            disp(['qq = ', num2str(qq), ': ', titles{qq}])
            divv = divvsK{qq + 1}(:, cols) ;    
            H2vn = H2vnsK{qq + 1}(:, cols) ;
            xx = divv(tidx_i, :) ;
            yy = H2vn(tidx_i, :) ;
            scatter(xx(:), yy(:), sz, ...
                markers{qq}, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', colors(qq, :), ...
                'MarkerEdgeAlpha', alphaVal) ;
            hold on ;
            allX = [allX(:); xx(:)] ;
            allY = [allY(:); yy(:)] ;
            
            % Label the x axis if on the bottom row
            xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                    'Interpreter', 'Latex') ;
            ylabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
        end
        
        % Add dashed y=x line
        xlims = get(gca, 'xlim') ;
        ylims = get(gca, 'ylim') ;
        axis equal
        xlim([max(-2 * climit, xlims(1)), min(2 * climit, xlims(2))])
        ylim([max(-2 * climit, ylims(1)), min(2 * climit, ylims(2))])
        xlims = get(gca, 'xlim') ;
        ylims = get(gca, 'ylim') ;
        leftdot = max(xlims(1), ylims(1)) ;
        rightdot = min(xlims(2), ylims(2)) ;
        plot([leftdot, rightdot], [leftdot, rightdot], 'k--')
         
        % % master titles (suptitles)
        % figure(1) ;
        % % Add colorbar
        % c = colorbar() ;
        % % Make colorbar share the alpha of the image
        % % Manually flush the event queue and force MATLAB to render the colorbar
        % % necessary on some versions
        % drawnow
        % % Get the color data of the object that correponds to the colorbar
        % colormap(cmap)
        % cdata = c.Face.Texture.CData;
        % % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
        % cdata(end,:) = uint8(alphaVal * cdata(end,:));
        % % Ensure that the display respects the alpha channel
        % c.Face.Texture.ColorType = 'truecoloralpha';
        % % Update the color data with the new transparency information
        % c.Face.Texture.CData = cdata;
        % c.Label.Interpreter = 'Latex' ;
        % c.Label.String = ['time [' QS.timeUnits ']'] ;
        % c.Ticks = [0, 1] ;
        % c.TickLabels = [tps(1), max(timeSpan_i)] ;

        % Label correlation
        % Correlation for 10%-90% AP axis
        RR = corrcoef(allX(:), allY(:)) ;
        rho = RR(1, 2) ;
        corrString = [' $\rho=$' sprintf('%0.3f', rho)] ;
        
        sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ', ...
            num2str(min(timeSpan_i)) '$<t<$' num2str(max(timeSpan_i)) ...
            ' ' QS.timeUnits corrString], ...
            'Interpreter', 'Latex') ;

        % Save figure
        figure(1)
        saveas(gcf, [fnout '.png']) ;
        % disp('saving pdf version...')
        % saveas(gcf, [fnout '.pdf']) ;
        disp('done')
        close all
        set(gcf, 'visible', 'off')
        
        % Save plot data as mat
        save([fnout '_data.mat'], 'allX', 'allY', 'timeSpan_i', 'tidx_i', 'rho')
    end
    disp('done with collapsed correlation plots betweeen divv and H2vn')
end




%% Metric Kinematic Correlations for each hour, one single axis only 
% (all quadrants) -- ONE MARKER ONLY
padN = round(0.1 * nU) ;
cols = padN:nU - padN ;
corrDVqDir = fullfile(mKDir, 'correlations_DVquadrants_collapsed_singleMarker') ;
files_exist = true ;
for qq = 1:length(collateTSpans)
    outputFileNames{qq} = ...
        fullfile(corrDVqDir, sprintf('correlation_allQuad_times%02d_div_2Hvn_singleMarker', qq)) ;
    files_exist = files_exist && exist([outputFileNames{qq} '.png'], 'file') ;
end

if plot_correlations && (~files_exist || overwrite)
    if ~exist(corrDVqDir, 'dir')
        mkdir(corrDVqDir)
    end    
    alphaVal = 0.2 ;
    sz = 4 ;
    cmap = phasemap(256) ;
    close all
    set(gcf, 'visible', 'off')
    % Consider each hr long timespan (early or entire series)
    for tspanIdx = 1:length(collateTSpans)
        fnout = outputFileNames{tspanIdx} ;
        timeSpan_i = collateTSpans{tspanIdx} ;
        
        minTP = max(min(timeSpan_i(1)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        maxTP = max(min(timeSpan_i(end)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        tidx_i = QS.xp.tIdx(minTP):QS.xp.tIdx(maxTP) ;
        
        markers = QS.plotting.markers ;
        colors = mapValueToColor(1:5, [1, 5], cmap) ;
        close all
        allX = [] ;
        allY = [] ;
        for qq = 1:4  % consider left, right, dorsal, ventral
            disp(['qq = ', num2str(qq), ': ', titles{qq}])
            divv = divvsK{qq + 1}(:, cols) ;    
            H2vn = H2vnsK{qq + 1}(:, cols) ;
            xx = divv(tidx_i, :) ;
            yy = H2vn(tidx_i, :) ;
            scatter(xx(:), yy(:), sz, ...
                'o', 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', 'k', ...
                'MarkerEdgeAlpha', alphaVal) ;
            hold on ;
            allX = [allX(:); xx(:)] ;
            allY = [allY(:); yy(:)] ;
            
            % Label the x axis if on the bottom row
            xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                    'Interpreter', 'Latex') ;
            ylabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
        end
        
        % Add dashed y=x line
        xlims = get(gca, 'xlim') ;
        ylims = get(gca, 'ylim') ;
        axis equal
        xlim([max(-2 * climit, xlims(1)), min(2 * climit, xlims(2))])
        ylim([max(-2 * climit, ylims(1)), min(2 * climit, ylims(2))])
        xlims = get(gca, 'xlim') ;
        ylims = get(gca, 'ylim') ;
        leftdot = max(xlims(1), ylims(1)) ;
        rightdot = min(xlims(2), ylims(2)) ;
        plot([leftdot, rightdot], [leftdot, rightdot], 'k--')
         
        % % master titles (suptitles)
        % figure(1) ;
        % % Add colorbar
        % c = colorbar() ;
        % % Make colorbar share the alpha of the image
        % % Manually flush the event queue and force MATLAB to render the colorbar
        % % necessary on some versions
        % drawnow
        % % Get the color data of the object that correponds to the colorbar
        % colormap(cmap)
        % cdata = c.Face.Texture.CData;
        % % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
        % cdata(end,:) = uint8(alphaVal * cdata(end,:));
        % % Ensure that the display respects the alpha channel
        % c.Face.Texture.ColorType = 'truecoloralpha';
        % % Update the color data with the new transparency information
        % c.Face.Texture.CData = cdata;
        % c.Label.Interpreter = 'Latex' ;
        % c.Label.String = ['time [' QS.timeUnits ']'] ;
        % c.Ticks = [0, 1] ;
        % c.TickLabels = [tps(1), max(timeSpan_i)] ;

        % Label correlation
        % Correlation for 10%-90% AP axis
        RR = corrcoef(allX(:), allY(:)) ;
        rho = RR(1, 2) ;
        corrString = [' $\rho=$' sprintf('%0.3f', rho)] ;
        
        sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ', ...
            num2str(min(timeSpan_i)) '$<t<$' num2str(max(timeSpan_i)) ...
            ' ' QS.timeUnits corrString], ...
            'Interpreter', 'Latex') ;

        % Save figure
        figure(1)
        saveas(gcf, [fnout '.png']) ;
        disp('saving pdf version...')
        saveas(gcf, [fnout '.pdf']) ;
        disp('done')
        close all
        set(gcf, 'visible', 'off')
        
        % Save plot data as mat
        save([fnout '_data.mat'], 'allX', 'allY', 'timeSpan_i', 'tidx_i', 'rho')
    end
    disp('done with collapsed correlation plots betweeen divv and H2vn')
end

%% Metric Kinematic Correlations for each timeSpan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation between terms div(v) and 2Hvn for each quadrant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
padN = round(0.1 * nU) ;
cols = padN:nU - padN ;

corrDVqDir = fullfile(mKDir, 'correlations_DVquadrants') ;
% Warning: not truly Lagrangian averaged if sigma > 0. 
files_exist = exist(outputFileNames{1}, 'file')  ;
for qq = 2:length(timeSpans)
    outputFileNames{qq} = ...
        fullfile(corrDVqDir, sprintf('correlation_times%02d_div_2Hvn.png', qq)) ;
    files_exist = files_exist && exist(outputFileNames{qq}, 'file') ;
end

if plot_correlations && (~files_exist || overwrite)
    if ~exist(corrDVqDir, 'dir')
        mkdir(corrDVqDir)
    end    
    alphaVal = 0.2 ;
    sz = 4 ;
    cmap = phasemap(256) ;
    close all
    set(gcf, 'visible', 'off')
   
    % Consider each timespan (early or entire series)
    for tspanIdx = 1:length(timeSpans)
        fnout = outputFileNames{tspanIdx} ;
        timeSpan_i = timeSpans{tspanIdx} ;
        ntspan = length(timeSpan_i) ;
        
        minTP = max(min(timeSpan_i(1)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        maxTP = max(min(timeSpan_i(end)/QS.timeInterval+tfold, max(vtimePoints)), min(vtimePoints)) ;
        tidx_i = QS.xp.tIdx(minTP):QS.xp.tIdx(maxTP) ;
        
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
        
        allX = [] ;
        allY = [] ;
        for qq = 1:4  % consider left, right, dorsal, ventral
            disp(['qq = ', num2str(qq), ': ', titles{qq}])
            divv = divvsK{qq + 1}(:, cols) ;    
            H2vn = H2vnsK{qq + 1}(:, cols) ;

            % Optional: smooth here
            % if sigma > 0
            %     divv = imgaussfilt(divv, sigma);            
            %     H2vn = imgaussfilt(H2vn, sigma);  
            % end

            % Check the quadrant kymographs
            figure(2) ;
            sphCollection2{qq} = subplot(2, 2, qq) ;
            imagesc(cols/nU, timeSpan_i, divv(tidx_i, :)); 
            caxis([-climit, climit])
            colormap(bwr256)
            figure(3) ;
            sphCollection3{qq} = subplot(2, 2, qq) ;
            imagesc(cols/nU, timeSpan_i, H2vn(tidx_i, :)); 
            caxis([-climit, climit])
            colormap(bwr256)

            figure(1) ;
            sphCollection{qq} = subplot(2, 2, qq) ;
            for row = tidx_i
                disp(['row = ', num2str(row)])
                scatter(divv(row, :), H2vn(row, :), sz, ...
                    markers{qq}, 'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', colors(row-tidx_i(1)+1, :), ...
                    'MarkerEdgeAlpha', alphaVal) ;
                hold on ;
                x2add = divv(row,:) ;
                y2add = H2vn(row, :) ;
                allX = [allX(:); x2add(:) ] ;
                allY = [allY(:); y2add(:) ] ;
            end

            % Label the x axis if on the bottom row
            if qq > 2
                figure(1)
                xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                        'Interpreter', 'Latex') ;
                figure(2)
                xlabel(['ap position, $\zeta/L$'], ...
                    'Interpreter', 'Latex') ;
                figure(3)
                xlabel(['ap position, $\zeta/L$'], ...
                    'Interpreter', 'Latex') ;
            end
            figure(1)
            axis equal
            % Add dashed y=x line
            xlims = get(gca, 'xlim') ;
            ylims = get(gca, 'ylim') ;
            xlim([max(-2 * climit, xlims(1)), min(2 * climit, xlims(2))])
            ylim([max(-2 * climit, ylims(1)), min(2 * climit, ylims(2))])
            xlims = get(gca, 'xlim') ;
            ylims = get(gca, 'ylim') ;
            leftdot = max(xlims(1), ylims(1)) ;
            rightdot = min(xlims(2), ylims(2)) ;
            plot([leftdot, rightdot], [leftdot, rightdot], 'k--')

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
            sposCollection{qq} = get(sphCollection{qq}, 'Position');
            sposCollection2{qq} = get(sphCollection2{qq}, 'Position');
            sposCollection3{qq} = get(sphCollection3{qq}, 'Position');
        end

        % Move subplots left a bit for colorbar space
        for qq = 1:length(sphCollection)
            spos = sposCollection{qq} ;
            wh = min(spos(3)-0.05, spos(4)) ;
            if mod(qq, 2) == 1
                set(sphCollection{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                set(sphCollection2{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                set(sphCollection3{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
            else
                set(sphCollection{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                set(sphCollection2{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                set(sphCollection3{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
            end
        end
        
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

        
        % Label correlation
        % Correlation for 10%-90% AP axis
        RR = corrcoef(allX(:), allY(:)) ;
        rho = RR(1, 2) ;
        corrString = [' $\rho=$' sprintf('%0.3f', rho)] ;
        
        if tspanIdx > 1
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ', ...
                num2str(min(timeSpan_i)) '$<t<$' num2str(max(timeSpan_i)) ...
                ' ' QS.timeUnits corrString], ...
                'Interpreter', 'Latex') ;
        else
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$ ' corrString], ...
                'Interpreter', 'Latex') ;
        end

        figure(2) ;
        c = colorbar('Position',[.9 .333 .02 .333]) ;
        if tspanIdx > 1
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ', ...
                num2str(min(timeSpan_i)) '$<t<$' num2str(max(timeSpan_i)) ...
                ' ' QS.timeUnits], ...
                'Interpreter', 'Latex') ;
        else
            sgtitle('$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$', ...
                'Interpreter', 'Latex') ;
        end
        figure(3) ;
        c = colorbar('Position',[.9 .333 .02 .333]) ;
        if tspanIdx > 1
            sgtitle(['$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$, ', ...
                num2str(min(timeSpan_i)) '$<t<$' num2str(max(timeSpan_i)) ...
                ' ' QS.timeUnits], ...
                'Interpreter', 'Latex') ;
        else
            sgtitle('$2Hv_n$ vs $\nabla \cdot \bf{v}_\parallel$', ...
                'Interpreter', 'Latex') ;
        end

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
    disp('done with correlation plots betweeen divv and H2vn')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlations between gdot and each term in the sum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputFileNames = cell(2, 1) ;
outputFileNames{1} = {fullfile(mKDir, 'correlation_alltime_div_gdot'), ...
                   fullfile(mKDir, 'correlation_earlytimes_div_gdot')} ;
outputFileNames{2} = {fullfile(mKDir, 'correlation_alltime_2Hvn_gdot'), ...
                   fullfile(mKDir, 'correlation_earlytimes_2Hvn_gdot')} ;
files_exist = exist(outputFileNames{1}{1}, 'file') && ...
    exist(outputFileNames{1}{2}, 'file') && ...
    exist(outputFileNames{2}{1}, 'file') && ...
    exist(outputFileNames{2}{2}, 'file') ;

if plot_gdot_correlations && (~files_exist || overwrite)
    alphaVal = 0.1 ;
    sz = 4 ;
    cmap = phasemap(256) ;
    close all
    set(gcf, 'visible', 'off')
    for pairIdx = 1:2
        for tspanIdx = 1:2
            fnout = outputFileNames{pairIdx}{tspanIdx} ;
            timeSpan_i = timeSpans{tspanIdx} ;
            ntspan = length(timeSpan_i) ;
            titles = {'left lateral', 'right lateral', 'dorsal', 'ventral'} ;
            markers = QS.plotting.markers ;
            colors = mapValueToColor(1:ntspan, [1, ntspan], cmap) ;
            close all
            cols = round(nV * [0.2, 0.85]) ;
            sphCollection = cell(4, 1) ;
            sposCollection = cell(4, 1) ;
            for qq = 1:4  % consider left, right, dorsal, ventral
                disp(['qq = ', num2str(qq), ': ', titles{qq}])
                if pairIdx == 1
                    divv = divvsK{qq + 1} ;
                    gdot = gdotsK{qq + 1} ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(divv(row, cols), gdot(row, cols), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        xlabel(['$\nabla \cdot \bf{v}_\parallel$ ' unitstr], ...
                                'Interpreter', 'Latex') ;
                    end
                    axis equal
                    % Add dashed y=x line
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    leftdot = max(xlims(1), ylims(1)) ;
                    rightdot = min(xlims(2), ylims(2)) ;
                    plot([leftdot, rightdot], [leftdot, rightdot], 'k--')
                else
                    H2vn = H2vnsK{qq + 1} ;
                    gdot = gdotsK{qq + 1} ;
                    sphCollection{qq} = subplot(2, 2, qq) ;
                    for row = 1:ntspan
                        disp(['row = ', num2str(row)])
                        scatter(H2vn(row, cols), gdot(row, cols), sz, ...
                            markers{qq}, 'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', colors(row, :), ...
                            'MarkerEdgeAlpha', alphaVal) ;
                        hold on ;
                    end

                    % Label the x axis if on the bottom row
                    if qq > 2
                        xlabel(['$2Hv_n$ ' unitstr], 'Interpreter', 'Latex') ;
                    end
                    axis equal
                    % Add dashed y=x line
                    xlims = get(gca, 'xlim') ;
                    ylims = get(gca, 'ylim') ;
                    leftdot = max(xlims(1), -ylims(2)) ;
                    rightdot = min(xlims(2), -ylims(1)) ;
                    plot([leftdot, rightdot], [-leftdot, -rightdot], 'k--')
                end

                % Label the y axis if on the left column
                if qq == 1 || qq == 3
                    ylabel(['$\frac{1}{2}\textrm{Tr}[g^{-1} \dot{g}]$ ' unitstr], ...
                            'Interpreter', 'Latex')
                end
                title(titles{qq}, 'Interpreter', 'Latex')

                % Grab axis position
                sposCollection{qq} = get(sphCollection{qq}, 'Position');
            end

            % Move subplots left a bit for colorbar space
            for qq = 1:length(sphCollection)
                spos = sposCollection{qq} ;
                wh = min(spos(3)-0.05, spos(4)) ;
                if mod(qq, 2) == 1
                    set(sphCollection{qq}, 'Position', [spos(1)-0.01, spos(2), wh, wh])
                else
                    set(sphCollection{qq}, 'Position', [spos(1)-0.06, spos(2), wh, wh])
                end
            end

            % Add colorbar
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

            % Save figure
            saveas(gcf, [fnout '.png']) ;
            saveas(gcf, [fnout '.pdf']) ;
            close all
            set(gcf, 'visible', 'off')
        end
    end
    disp('done')
end


%% Metric Kinematics -- plot isotropic component only (part of decomposition)
if plot_gdot_decomp
    % Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
    dvDir = fullfile(mKDir, 'avgDV') ;
    lDir = fullfile(mKDir, 'avgLeft') ;
    rDir = fullfile(mKDir, 'avgRight') ;
    dDir = fullfile(mKDir, 'avgDorsal') ;
    vDir = fullfile(mKDir, 'avgVentral') ;
    outdirs = {dvDir, lDir, rDir, dDir, vDir} ;
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    gdotK0 = gdotsK{1} ;
    isogrowth = sum(gdotK0, 2) ;
    
    %% Plot as 1d Curve
    plot(tps, isogrowth); 
    xlabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
    ylabel('Isotropic component of growth') 
    saveas(gcf, fullfile(odir, 'average_growth.png'))  
    
    %% Plot quadrant contributions
    for qq = 1:length(outdirs)
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
        m2plot = {gdotK, HHK, divvK, velnK, H2vnK} ;
        titles = {'$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
            'mean curvature, $H$', ...
            'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
            'normal velocity, $v_n$', ...
            'normal motion, $v_n 2 H$'} ;
        labels = {['$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
            ['mean curvature, $H$ ' Hunitstr], ...
            ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
            ['normal velocity, $v_n$ ' vunitstr] , ...
            ['normal motion, $v_n 2 H $ ' unitstr]} ;
        names = {'gdot', 'HH', 'divv', 'veln', 'H2vn'} ;
        climits = [climit, climit_H, climit, climit_veln, climit_err] ;

        %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
        for pp = 1:length(m2plot)            
            fn = fullfile(odir, [ names{pp} '.png']) ;
            if ~exist(fn, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                imagesc((1:nU)/nU, tps, m2plot{pp})
                caxis([-climits(pp), climits(pp)])
                colormap(bwr256)
                % Add folds to plot
                hold on;
                fons1 = max(1, fonsIdx(1)) ;
                fons2 = max(1, fonsIdx(2)) ;
                fons3 = max(1, fonsIdx(3)) ;
                plot(folds.folds(fons1:end-1, 1) / nU, tps(fons1:end))
                plot(folds.folds(fons2:end-1, 2) / nU, tps(fons2:end))
                plot(folds.folds(fons3:end-1, 3) / nU, tps(fons3:end))

                % title and save
                title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
                disp(['saving kinetic terms as plots: ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
            end
        end
    end
end
