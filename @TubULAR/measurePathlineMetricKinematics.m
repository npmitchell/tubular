function measurePathlineMetricKinematics(QS, options)
% measurePathlineMetricKinematics(QS, options)
%   Query the metric Kinematics along lagrangian pathlines.
%   Plot results as kymographs and correlation plots.
%   Out-of-plane motion is v_n * 2H, where v_n is normal velocity and H is
%   mean curvature.
%   In-plane motion considered here is div(v_t) where v_t is tangential
%   velocity on the curved surface.
%   The difference div(v_t) - vn*2H = Tr[g^{-1} dot{g}], which is a measure
%   of isotropic metric change over time (dot is the partial derivative wrt
%   time). 
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
plot_kymographs = true ;
plot_kymographs_cumsum = true ;
plot_correlations = true ;
plot_gdot_correlations = false ;
plot_gdot_decomp = true ;

%% Parameter options
lambda = QS.smoothing.lambda ; 
lambda_mesh = QS.smoothing.lambda_mesh ;
lambda_err = QS.smoothing.lambda_err ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;
climit = 0.2 ;
climit_err = 0.2 ;
climit_veln = climit * 10 ;
climit_H = climit * 2 ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = "Lagrangian" ;
QS.t0set() ;
t0Pathline = QS.t0 ;

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
elseif isfield(options, 'lambda_error')
    lambda_err = options.lambda_error ;
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
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_kymographs_cumsum')
    plot_kymographs_cumsum = options.plot_kymographs_cumsum ;
end
if isfield(options, 'plot_correlations')
    plot_correlations = options.plot_correlations ;
end
if isfield(options, 'plot_gdot_correlations')
    plot_gdot_correlations = options.plot_gdot_correlations ;
end
if isfield(options, 'plot_gdot_decomp')
    plot_gdot_decomp = options.plot_gdot_decomp ;
end
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
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
if strcmp(averagingStyle, 'Lagrangian')
    mKDir = fullfile(QS.dir.metricKinematics.root, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
        lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));
else
    mKDir = fullfile(QS.dir.metricKinematicsSimple, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02'], ...
        lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));
end

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

%% Build timepoint list so that we first do every 10, then fill in details
lastIdx = length(QS.xp.fileMeta.timePoints) - 1 ;
coarseIdx = 1:10:lastIdx ;
fineIdx = setdiff(1:lastIdx, coarseIdx) ;
allIdx = [coarseIdx, fineIdx ] ;
tp2do = QS.xp.fileMeta.timePoints(allIdx) ;

% DONE WITH PREPARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pathlines to build Kymographs along pathlines
QS.loadPullbackPathlines(t0Pathline, 'vertexPathlines')
vP = QS.pathlines.vertices ;

% Output directory is inside metricKinematics dir
mKPDir = fullfile(mKDir, sprintf('pathline_%04dt0', t0Pathline)) ;
outdir = fullfile(mKPDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
% Data for kinematics on meshes (defined on vertices)
mdatdir = fullfile(mKDir, 'measurements') ;

% Load Lx, Ly by loadingPIV
QS.loadPIV()

% Discern if piv measurements are done on a double covering or the meshes
if strcmp(QS.piv.imCoords(end), 'e')
    doubleCovered = true ;
end

% Compute or load all timepoints
for tp = tp2do
    close all
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;
    QS.setTime(tp) ;
    
    % Check for timepoint measurement on disk
    Hfn = fullfile(outdir, sprintf('HH_pathline%04d_%06d.mat', t0Pathline, tp))   ;
    efn = fullfile(outdir, sprintf('gdot_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    dfn = fullfile(outdir, sprintf('divv_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    nfn = fullfile(outdir, sprintf('veln_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    rfn = fullfile(outdir, sprintf('radius_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    H2vnfn = fullfile(outdir, sprintf('H2vn_pathline%04d_%06d.mat', t0Pathline, tp)) ;
    files_missing = ~exist(Hfn, 'file') || ~exist(efn, 'file') || ...
         ~exist(dfn, 'file') || ~exist(nfn, 'file') || ...
          ~exist(H2vnfn, 'file') || ~exist(rfn, 'file') ;
    
    if overwrite || files_missing
        disp('Computing pathline metric kinematics...')
        % Load timeseries measurements defined on mesh vertices
        HfnMesh = fullfile(mdatdir, sprintf('HH_vertices_%06d.mat', tp))   ;
        efnMesh = fullfile(mdatdir, sprintf('gdot_vertices_%06d.mat', tp)) ;
        dfnMesh = fullfile(mdatdir, sprintf('divv_vertices_%06d.mat', tp)) ;
        nfnMesh = fullfile(mdatdir, sprintf('veln_vertices_%06d.mat', tp)) ;
        rfnMesh = fullfile(mdatdir, sprintf('radius_vertices_%06d.mat', tp)) ;
        H2vnfnMesh = fullfile(mdatdir, sprintf('H2vn_vertices_%06d.mat', tp)) ;

        try
            load(HfnMesh, 'HH_filt')
            load(efnMesh, 'gdot_filt')
            load(dfnMesh, 'divv_filt')
            load(nfnMesh, 'veln_filt') 
            load(H2vnfnMesh, 'H2vn_filt') 
            load(rfnMesh, 'radius_filt') 
            HH = HH_filt ;
            gdot = gdot_filt ;
            divv = divv_filt ;
            veln = veln_filt ;
            H2vn = H2vn_filt ;
            radius = radius_filt ;
        catch
            msg = 'Run QS.measureMetricKinematics() ' ;
            msg = [msg 'with lambdas=(mesh,lambda,err)=('] ;
            msg = [msg num2str(lambda_mesh) ','] ;
            msg = [msg num2str(lambda) ','] ;
            msg = [msg num2str(lambda_err) ')'] ;
            msg = [msg ' before running ', ...
                    'QS.measurePathlineMetricKinematics()'] ;
            error(msg)
        end
        % Interpolate from vertices onto pathlines
        xx = vP.vX(tidx, :, :) ;
        yy = vP.vY(tidx, :, :) ;
        XY = [xx(:), yy(:)] ;
        Lx = vP.Lx(tidx) ;
        Ly = vP.Ly(tidx) ;
        options.Lx = Lx ;
        options.Ly = Ly ;
        XY = QS.doubleToSingleCover(XY, Ly) ;
        HH = QS.interpolateOntoPullbackXY(XY, HH, options) ;
        gdot = QS.interpolateOntoPullbackXY(XY, gdot, options) ;
        divv = QS.interpolateOntoPullbackXY(XY, divv, options) ;
        veln = QS.interpolateOntoPullbackXY(XY, veln, options) ;
        H2vn = QS.interpolateOntoPullbackXY(XY, H2vn, options) ;
        radius = QS.interpolateOntoPullbackXY(XY, radius, options) ;
                
        % OPTION 1: simply reshape, tracing each XY dot to its t0Pathline
        % grid coordinate
        HH = reshape(HH, [nU, nV]) ;
        gdot = reshape(gdot, [nU, nV]) ;
        divv = reshape(divv, [nU, nV]) ;
        veln = reshape(veln, [nU, nV]) ;
        H2vn = reshape(H2vn, [nU, nV]) ;
        radius = reshape(radius, [nU, nV]) ;
        
        %% OPTION 2: the following regrids onto original XY coordinates,
        % rendering the process of following pathlines moot. 
        % Average into AP bins and take mean along 1/4 DV hoop arcs
        % if doubleCovered
        %     vminmax = [0.25 * Ly, 0.75 * Ly] ;
        % else
        %     vminmax = [1, Ly] ;
        % end
        %
        % Note the transposition: to plot as APDV, imshow(m')
        % HH = binData2dGrid([XY, HH], [1,Lx], vminmax, nU, nV) ;
        % gdot = binData2dGrid([XY, gdot], [1,Lx], vminmax, nU, nV) ;
        % divv = binData2dGrid([XY, divv], [1,Lx], vminmax, nU, nV) ;
        % veln = binData2dGrid([XY, veln], [1,Lx], vminmax, nU, nV) ;
        % H2vn = binData2dGrid([XY, H2vn], [1,Lx], vminmax, nU, nV) ;
           
        % Average along DV -- do not ignore last row at nV since not quite
        % redundant in this version of the algorithm -- is that true?
        HH_ap = nanmean(HH, 2) ;
        gdot_ap = nanmean(gdot, 2) ;
        divv_ap = nanmean(divv, 2) ;
        veln_ap = nanmean(veln, 2) ;
        H2vn_ap = nanmean(H2vn, 2) ;
        radius_ap = nanmean(radius, 2) ;
        
        % quarter bounds
        q0 = round(nV * 0.125) ;
        q1 = round(nV * 0.375) ;
        q2 = round(nV * 0.625) ;
        q3 = round(nV * 0.875) ;
        left = q0:q1 ;
        ventral = q1:q2 ;
        right = q2:q3 ;
        dorsal = [q3:nV, 1:q1] ;
        
        % left quarter
        HH_l = nanmean(HH(:, left), 2) ;
        gdot_l = nanmean(gdot(:, left), 2) ;
        divv_l = nanmean(divv(:, left), 2) ;
        veln_l = nanmean(veln(:, left), 2) ;
        H2vn_l = nanmean(H2vn(:, left), 2) ;
        radius_l = nanmean(radius(:, left), 2) ;
        
        % right quarter
        HH_r = nanmean(HH(:, right), 2) ;
        gdot_r = nanmean(gdot(:, right), 2) ;
        divv_r = nanmean(divv(:, right), 2) ;
        veln_r = nanmean(veln(:, right), 2) ;
        H2vn_r = nanmean(H2vn(:, right), 2) ;
        radius_r = nanmean(radius(:, right), 2) ;
        
        % dorsal quarter
        HH_d = nanmean(HH(:, dorsal), 2) ;
        gdot_d = nanmean(gdot(:, dorsal), 2) ;
        divv_d = nanmean(divv(:, dorsal), 2) ;
        veln_d = nanmean(veln(:, dorsal), 2) ;
        H2vn_d = nanmean(H2vn(:, dorsal), 2) ;
        radius_d = nanmean(radius(:, dorsal), 2) ;
        
        % ventral quarter
        HH_v = nanmean(HH(:, ventral), 2) ;
        gdot_v = nanmean(gdot(:, ventral), 2) ;
        divv_v = nanmean(divv(:, ventral), 2) ;
        veln_v = nanmean(veln(:, ventral), 2) ;
        H2vn_v = nanmean(H2vn(:, ventral), 2) ;
        radius_v = nanmean(radius(:, ventral), 2) ;
        
        % Save results
        save(Hfn, 'HH', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
        save(efn, 'gdot', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
        save(dfn, 'divv', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
        save(nfn, 'veln', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
        save(H2vnfn, 'H2vn', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v')
        save(rfn, 'radius', 'radius_ap', 'radius_l', 'radius_r', ...
            'radius_d', 'radius_v')
    end
end
disp('done with measuring pathline metric kinematics')

%% Combine DV-averaged profiles into kymographs
apKymoFn = fullfile(outdir, 'apKymographMetricKinematics.mat') ;
lKymoFn = fullfile(outdir, 'leftKymographMetricKinematics.mat') ;
rKymoFn = fullfile(outdir, 'rightKymographMetricKinematics.mat') ;
dKymoFn = fullfile(outdir, 'dorsalKymographMetricKinematics.mat') ;
vKymoFn = fullfile(outdir, 'ventralKymographMetricKinematics.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if ~files_exist || overwrite
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        close all
        disp(['t = ' num2str(tp)])
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        Hfn = fullfile(outdir, sprintf('HH_pathline%04d_%06d.mat', t0Pathline, tp))   ;
        efn = fullfile(outdir, sprintf('gdot_pathline%04d_%06d.mat', t0Pathline, tp)) ;
        dfn = fullfile(outdir, sprintf('divv_pathline%04d_%06d.mat', t0Pathline, tp)) ;
        nfn = fullfile(outdir, sprintf('veln_pathline%04d_%06d.mat', t0Pathline, tp)) ;
        rfn = fullfile(outdir, sprintf('radius_pathline%04d_%06d.mat', t0Pathline, tp)) ;
        H2vnfn = fullfile(outdir, sprintf('H2vn_pathline%04d_%06d.mat', t0Pathline, tp)) ;

        % Load timeseries measurements
        load(Hfn, 'HH', 'HH_ap', 'HH_l', 'HH_r', 'HH_d', 'HH_v')
        load(efn, 'gdot', 'gdot_ap', 'gdot_l', 'gdot_r', 'gdot_d', 'gdot_v')
        load(dfn, 'divv', 'divv_ap', 'divv_l', 'divv_r', 'divv_d', 'divv_v')
        load(nfn, 'veln', 'veln_ap', 'veln_l', 'veln_r', 'veln_d', 'veln_v') 
        load(H2vnfn, 'H2vn', 'H2vn_ap', 'H2vn_l', 'H2vn_r', 'H2vn_d', 'H2vn_v') 
        load(rfn, 'radius', 'radius_ap', 'radius_l', 'radius_r', 'radius_d', 'radius_v') 

        %% Store in matrices
        % dv averaged
        HH_apM(tidx, :) = HH_ap ;
        gdot_apM(tidx, :) = gdot_ap ;
        divv_apM(tidx, :) = divv_ap ;
        veln_apM(tidx, :) = veln_ap ;
        H2vn_apM(tidx, :) = H2vn_ap ;
        radius_apM(tidx, :) = radius_ap ; 

        % left quarter
        HH_lM(tidx, :) = HH_l ;
        gdot_lM(tidx, :) = gdot_l ;
        divv_lM(tidx, :) = divv_l ;
        veln_lM(tidx, :) = veln_l ;
        H2vn_lM(tidx, :) = H2vn_l ;
        radius_lM(tidx, :) = radius_l ; 

        % right quarter
        HH_rM(tidx, :) = HH_r ;
        gdot_rM(tidx, :) = gdot_r ;
        divv_rM(tidx, :) = divv_r ;
        veln_rM(tidx, :) = veln_r ;
        H2vn_rM(tidx, :) = H2vn_r ;
        radius_rM(tidx, :) = radius_r ; 

        % dorsal quarter
        HH_dM(tidx, :) = HH_d ;
        gdot_dM(tidx, :) = gdot_d ;
        divv_dM(tidx, :) = divv_d ;
        veln_dM(tidx, :) = veln_d ;
        H2vn_dM(tidx, :) = H2vn_d ;
        radius_dM(tidx, :) = radius_d ; 

        % ventral quarter
        HH_vM(tidx, :) = HH_v ;
        gdot_vM(tidx, :) = gdot_v ;
        divv_vM(tidx, :) = divv_v ;
        veln_vM(tidx, :) = veln_v ;
        H2vn_vM(tidx, :) = H2vn_v ;
        radius_vM(tidx, :) = radius_v ; 
    end
    
    disp('Saving DV-averaged kymograph data')
    % Save the DV-averaged kymographs
    save(apKymoFn, 'HH_apM', 'gdot_apM', 'divv_apM', 'veln_apM', ...
        'H2vn_apM', 'radius_apM')
    save(lKymoFn, 'HH_lM', 'gdot_lM', 'divv_lM', 'veln_lM', ...
        'H2vn_lM', 'radius_lM')
    save(rKymoFn, 'HH_rM', 'gdot_rM', 'divv_rM', 'veln_rM', ...
        'H2vn_rM', 'radius_rM')
    save(dKymoFn, 'HH_dM', 'gdot_dM', 'divv_dM', 'veln_dM', ...
        'H2vn_dM', 'radius_dM')
    save(vKymoFn, 'HH_vM', 'gdot_vM', 'divv_vM', 'veln_vM', ...
        'H2vn_vM', 'radius_vM')
    
end
disp('done measuring pathline metric kinematics (div(vt) and 2Hvn, etc)')