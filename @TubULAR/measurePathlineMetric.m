function measurePathlineStrain(QS, options)
% measurePathlineStrain(QS, options)
%   Integrate the metric strain rate along Lagrangian pathlines.
%   Allow for median filtering along Lagrangian pathlines to avoid 
%   spurious spikes in accumulated strain.
%   Plot results in 2d and/or 3d for each timepoint.
%   
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   overwrite : bool, overwrite data on disk
%   overwriteImages : bool, overwrite images of results on disk
%   plot_comparison : bool, plot a comparison with DEC traceful dilation
%   median_filter_strainRates : bool
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
overwriteImages = false ;
plot_dzdp = false ;
plot_comparison = false ;
median_filter_strainRates = false ;
climitInitial = 0.05 ;
climitRamp = 0.005 ;
climitRatio = 1 ;

%% Parameter options
lambda_mesh = 0.002 ;
lambda = 0.01 ; 
debug = false ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = "Lagrangian" ;
LagrangianFrame = 'zetaphi' ;
% Load time offset for first fold, t0 -- default pathline t0
QS.t0set() ;
t0 = QS.t0 ;
% By default, t0Pathline = t0 (see below)

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'plot_comparison')
    plot_comparison = options.plot_comparison ;
end
if isfield(options, 'plot_dzdp')
    plot_dzdp = options.plot_dzdp ;
end

%% parameter options
if isfield(options, 'climitInitial')
    climitInitial = options.climitInitial ;
end
if isfield(options, 'climitRamp')
    climitRamp = options.climitRamp ;
end
if isfield(options, 'climitRatio')
    climitRatio = options.climitRatio ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
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
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
else
    t0Pathline = t0 ;
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
close all

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - t0;

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];
vunitstr = [ '[' QS.spaceUnits '/' QS.timeUnits ']' ];
    
% DONE WITH PREPARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pathlines to build metric from advected mesh faces
QS.loadPullbackPathlines(t0Pathline, 'vertexPathlines')
vP = QS.pathlines.vertices ;

% Output directory is inside pathline dir
outdir = QS.dir.pathlines.fundForms ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Load Lx, Ly by loadingPIV. 
QS.loadPIV()
Xpiv = QS.piv.raw.x ;
Ypiv = QS.piv.raw.y ;

% Also need velocities to advect mesh
% QS.loadVelocityAverage('vv')
% vPIV = QS.velocityAverage.vv ;

% Discern if piv measurements are done on a double covering or the meshes
if strcmp(QS.piv.imCoords(end), 'e')
    doubleCovered = true ;
end

%% INTEGRATE STRAINRATE INTO STRAIN ON PATHLINES
% Compute or load all timepoints
load(sprintf(QS.fileName.pathlines.v3d, t0), 'v3dPathlines')
vX3rs = v3dPathlines.vXrs ;
vY3rs = v3dPathlines.vYrs ;
vZ3rs = v3dPathlines.vZrs ;
t0 = v3dPathlines.t0 ;
tIdx0 = v3dPathlines.tIdx0 ;

ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
for tidx = 1:ntps
    % Identify current timepoint
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    
    % Do the fund forms for the lagrangian frame already exist?
    ffn = sprintf(QS.fullFileBase.pathline.fundForms, t0, tp) ;
    
    if ~exist(ffn, 'file')
        disp(['t = ' num2str(tp)])
        QS.setTime(tp) ;

        % Load Lagrangian advected vertices as mesh for this timepoint
        % load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp), 'spcutMeshSmRS') ;
        xx = vX3rs(tidx, :) ;
        yy = vY3rs(tidx, :) ;
        zz = vZ3rs(tidx, :) ;
        v3d = [ xx(:), yy(:), zz(:) ] ;
        mesh = struct() ;
        mesh.v = v3d ;
        mesh.u = lagrangianFrame ;
        pathPairs = [ (1:nU)', (nV-1)*nU + (1:nU)' ] ;
        mesh.pathPairs = pathPairs ;

        %% Tile the mesh and the associated velocity vectors to triple cover
        tileCount = [1, 1] ;
        [ TF, TV2D, TV3D, ~ ] = tileAnnularCutMesh( mesh, tileCount ) ;
        % The face list should now be 3x the original size
        assert(size(TF, 1) == 3 * size(cutMesh.f, 1)) ;

        %% Compute the fundamental forms
        [gg, ~] = constructFundamentalForms(TF, TV3D, TV2D) ;
        % Use tiled, open mesh (glued seam) to compute the second fundamental form
        [~, bb] = constructFundamentalForms(TF, TV3D, TV2D) ;

        fundForms.gg = gg ;
        fundForms.bb = bb ;
    else
        disp(['FundForms for t = ' num2str(tp) ' already on disk'])
    end
end
disp('done with integrated pathline strain calculations')


%% Combine DV-averaged gg and bb into kymographs
apKymoFn = fullfile(outdir, 'apKymographLagrangianMetric.mat') ;
lKymoFn = fullfile(outdir, 'leftKymographLagrangianMetric.mat') ;
rKymoFn = fullfile(outdir, 'rightKymographLagrangianMetric.mat') ;
dKymoFn = fullfile(outdir, 'dorsalKymographLagrangianMetric.mat') ;
vKymoFn = fullfile(outdir, 'ventralKymographLagrangianMetric.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if ~files_exist || overwrite || true
    disp('Compiling kymograph data to save to disk...')
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        srfn = fullfile(outdir, sprintf('strain_%06d.mat', tp))   ;

        % Load timeseries measurements
        load(srfn, 'strain_tr_ap', 'strain_tr_l', 'strain_tr_r', ...
            'strain_tr_d', 'strain_tr_v', ...
            'strain_th_ap', 'strain_th_l', 'strain_th_r', ....
            'strain_th_d', 'strain_th_v', ...
            'strain_dv_ap', 'strain_dv_l', 'strain_dv_r', ...
            'strain_dv_d', 'strain_dv_v') ;

        %% Store accumulated strain in matrices
        % dv averaged
        str_apM(tidx, :) = strain_tr_ap ;
        sdv_apM(tidx, :) = strain_dv_ap ;
        sth_apM(tidx, :) = strain_th_ap ;

        % left quarter
        str_lM(tidx, :) = strain_tr_l ;
        sdv_lM(tidx, :) = strain_dv_l ;
        sth_lM(tidx, :) = strain_th_l ;

        % right quarter
        str_rM(tidx, :) = strain_tr_r ;
        sdv_rM(tidx, :) = strain_dv_r ;
        sth_rM(tidx, :) = strain_th_r ;

        % dorsal quarter
        str_dM(tidx, :) = strain_tr_d ;
        sdv_dM(tidx, :) = strain_dv_d ;
        sth_dM(tidx, :) = strain_th_d ;

        % ventral quarter
        str_vM(tidx, :) = strain_tr_v ;
        sdv_vM(tidx, :) = strain_dv_v ;
        sth_vM(tidx, :) = strain_th_v ;
    end
    
    %% Save kymographs
    disp('Saving kymograph data files for Lagrangian pathlines')
    save(apKymoFn, 'str_apM', 'sdv_apM', 'sth_apM')
    disp(['Saved kymograph data to: ' apKymoFn])
    save(lKymoFn, 'str_lM', 'sdv_lM', 'sth_lM')
    disp(['Saved kymograph data to: ' lKymoFn])
    save(rKymoFn, 'str_rM', 'sdv_rM', 'sth_rM')
    disp(['Saved kymograph data to: ' rKymoFn])
    save(dKymoFn, 'str_dM', 'sdv_dM', 'sth_dM')
    disp(['Saved kymograph data to: ' dKymoFn])
    save(vKymoFn, 'str_vM', 'sdv_vM', 'sth_vM')
    disp(['Saved kymograph data to: ' vKymoFn])
    disp('done with strain kymograph data saving')
else
    disp('strain kymograph data already on disk')    
end

