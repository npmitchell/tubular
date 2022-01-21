function plotMetricKinematicsInEulerianFrame(QS, options)
% plotMetricKinematicsInEulerianFrame(QS, options)
%   Plot the metric Kinematics as kymographs with absolute AP position in
%   Eulerian frame on the x axis.
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
%   
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;

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

%% Operational options

if isfield(options, 'plot_kymographs_cumsum')
    plot_kymographs_cumsum = options.plot_kymographs_cumsum ;
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
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
close all
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;

%% Load time offset for first fold, t0
QS.t0set() ;
tfold = QS.t0 ;

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
tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;

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
outdir = fullfile(mKDir, 'measurementsEulerianParameterization') ;
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


%% First grab x positions in Eulerian space 
nU = QS.nU ;
nV = QS.nV ;
q0 = round(nV * 0.125) ;
q1 = round(nV * 0.375) ;
q2 = round(nV * 0.625) ;
q3 = round(nV * 0.875) ;
left = q0:q1 ;
ventral = q1:q2 ;
right = q2:q3 ;
dorsal = [q3:nV, 1:q1] ;

% preallocate
ntps = length(QS.xp.fileMeta.timePoints) ;
nAP = QS.nU ;
xDV_allTime = zeros(ntps, nAP) ;
xD_allTime = zeros(ntps, nAP) ;
xV_allTime = zeros(ntps, nAP) ;
xL_allTime = zeros(ntps, nAP) ;
xR_allTime = zeros(ntps, nAP) ;
tpsGrid = zeros(ntps, nAP) ;
for tidx = 1:ntps
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    QS.setTime(tp)
    QS.loadCurrentSPCutMeshSmRS() ;
    mm = QS.currentMesh.spcutMeshSmRS ;
    
    % Find dosal vertex ap position
    vx = reshape(mm.v(:,1), [nU, nV]) ;
    xDV = mean(vx(:, 1:end-1), 2) ;
    xL = mean(vx(:, left), 2) ;
    xR = mean(vx(:, right), 2) ;
    xD = mean(vx(:, dorsal), 2) ;
    xV = mean(vx(:, ventral), 2) ;
    xDV_allTime(tidx, :) = xDV ;
    xD_allTime(tidx, :) = xD ;
    xV_allTime(tidx, :) = xV ;
    xL_allTime(tidx, :) = xL ;
    xR_allTime(tidx, :) = xR ;
    tpsGrid(tidx, :) = (tp - QS.t0) * QS.timeInterval ;
end

%% Now plot different measured quantities as kymographs

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
    H2vnK = H2vnsK{qq} ;
    m2plot = {gdotK, HHK, divvK, H2vnK} ;
    titles = {'$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]=\nabla\cdot\mathbf{v}_\parallel-v_n 2H$',...
        'mean curvature, $H$', ...
        'divergence of flow, $\nabla \cdot \mathbf{v}$', ...
        'bending, $v_n 2 H$'} ;
    labels = {['$\frac{1}{2}\textrm{Tr}[g^{-1}\dot{g}]$ ' unitstr], ...
        ['mean curvature, $H$ ' Hunitstr], ...
        ['$\nabla \cdot \mathbf{v}$ ' unitstr], ...
        ['normal motion, $v_n 2 H $ ' unitstr]} ;
    
    names = {'gdot', 'HH', 'divv',  'H2vn'} ;
    climits = [climit, climit_H, climit, climit] ;
    cmaps = {bwr256, bwr256, bwr256, bwr256};

    %% Plot gdot/HH/divv/veln/H2vn DV-averaged kymograph
    for pp = 1:length(m2plot)

        % Check if images already exist on disk
        fn = fullfile(odir, [ names{pp} '.png']) ;
        fn_zoom = fullfile(odir, [names{pp} '_zoom_early.png']) ;

        if ~exist(fn, 'file') || ~exist(fn_zoom, 'file') || overwrite
            close all
            set(gcf, 'visible', 'off')
            
            % Here interpolate the matrix onto Eulerian frame near folds
            if qq == 1
                si = scatteredInterpolant(xDV_allTime, tpsGrid, m2plot{pp}) ;
            elseif qq == 2
                si = scatteredInterpolant(xL_allTime, tpsGrid, m2plot{pp}) ;
            elseif qq == 3
                si = scatteredInterpolant(xR_allTime, tpsGrid, m2plot{pp}) ;
            elseif qq == 4
                si = scatteredInterpolant(xD_allTime, tpsGrid, m2plot{pp}) ;
            elseif qq == 5
                si = scatteredInterpolant(xV_allTime, tpsGrid, m2plot{pp}) ;
            else
                error('unrecognized averaging DV/D/V/L/R')
            end
            values = si(xfixed, tpsGrid) ;
            error('stop')
            
            imagesc(xfixed, tps, values)
            if climits(pp) > 0
                caxis([-climits(pp), climits(pp)])
            end
            colormap(cmaps{pp})
            % Add folds to plot
            hold on;
            for ii = 1:length(fons)
                fonsi = max(1, fons(ii) ) ; 
                plot(double(folds.folds(fonsi:end-1, ii)) / double(nU), ...
                    tps(fonsi:end))            
            end

            % title and save
            title([titles{pp}, titleadd{qq}], 'Interpreter', 'Latex')
            ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
            xlabel('ap position [$\mu$m]', 'Interpreter', 'Latex')
            cb = colorbar() ;
            ylabel(cb, labels{pp}, 'Interpreter', 'Latex')  
            fn = fullfile(odir, [ names{pp} '.png']) ;
            disp(['saving ', fn])
            export_fig(fn, '-png', '-nocrop', '-r200')   

            % Zoom in on small values
            if climits(pp) > 0
                caxis([-climits(pp)/3, climits(pp)/3])
                fn = fullfile(odir, [names{pp} '_zoom.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                % Zoom in on early times
                ylim([min(tps), min(max(tps), max(fons-tfold) + 10)])
                caxis([-climits(pp)/3, climits(pp)/3])
                fn = fullfile(odir, [names{pp} '_zoom_early.png']) ;
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')  
            end
        end
    end
end


