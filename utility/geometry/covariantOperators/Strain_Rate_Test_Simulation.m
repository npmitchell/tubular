%% Make a necking tube, compute the strain rate and lagrangian strain
% 
% addpath('/mnt/data/code/gut_matlab/addpath_recurse/')
% addpath('/mnt/data/code/gut_matlab/plotting/')
% addpath('/mnt/data/code/gut_matlab/plotting/export_fig/')
% addpath('/mnt/data/code/gut_matlab/mesh_handling/')
% addpath_recurse('/mnt/data/code/gptoolbox/')
clear all
nU = 100 ;
nV = 100 ;
velType = 'linearCompression' ;  % 'zero' 'linearCompression' 'linearExtension'
vnscale = 5 ;
vtscale = 5 ;

%% PATHS ==================================================================
if nU == nV 
    meshDir = ['/mnt/data/simulations/strain_tube/necking_10tp_n', ...
        num2str(nU), '_', velType, '/'] ;
else
    meshDir = ['/mnt/data/simulations/strain_tube/necking_10tp_nU', ...
        num2str(nU), '_nV', num2str(nV), '_', velType, '/'] ;
end
dataDir = meshDir ;
if ~exist(meshDir, 'dir')
    mkdir(meshDir)
end
cd('/mnt/data/code/gut_matlab/master_pipeline/')
aux_paths_and_colors
cd(dataDir)

%% QS DEFINITION
opts.meshDir = meshDir ;
% set up shell Experiment class instance
xp = project.Experiment(meshDir, meshDir);
% A filename base template - to be used throughout this script
fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = 'data_%03d.tif';
fileMeta.nChannels          = 1 ;
fileMeta.timePoints         = 0:9 ;
fileMeta.stackResolution    = [1,1,1];
fileMeta.swapZT             = 1;
first_tp = 1 ;
expMeta                     = struct();
expMeta.channelsUsed        = 1;
expMeta.channelColor        = 1;
expMeta.description         = 'necking tube';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();
detectOptions = struct( 'channel', 1, ...
        'ssfactor', 1, ...
        'niter', 1,...
        'niter0', 1, ...
        'lambda1', 1, ...
        'lambda2', 1, ...
        'nu', 1, ...
        'smoothing', 1, ...
        'post_nu', 1, ...
        'post_smoothing', 1, ...
        'exit_thres', 1, ...
        'foreGroundChannel', 1, ...
        'fileName', 'test_data', ...
        'mslsDir', meshDir, ...
        'ofn_ls', 'ls', ...
        'ofn_ply', 'mesh_tmp',...
        'ofn_smoothply', 'mesh', ...
        'ms_scriptDir', 'none', ...
        'timepoint', xp.currentTime, ...
        'zdim', 1, ...
        'pre_nu', 1, ...
        'pre_smoothing', 1, ...
        'mlxprogram', 'none', ...
        'init_ls_fn', 'none', ... % set to none to load prev tp
        'run_full_dataset', 'none',... % projectDir, ... % set to 'none' for single tp
        'radius_guess', 'none', ...
        'dset_name', 'exported_data',...
        'center_guess', 'none',... % xyz of the initial guess sphere ;
        'save', false, ... % whether to save images of debugging output
        'plot_mesh3d', false, ...
        'dtype', 'h5',...
        'mask', 'none',...
        'mesh_from_pointcloud', false, ...
        'prob_searchstr', 'none', ...
        'preilastikaxisorder', 'xyzc', ... 
        'ilastikaxisorder', 'cxyz', ... 
        'physicalaxisorder', 'xyzc', ... 
        'include_boundary_faces', true, ...
        'smooth_with_matlab', 0) ;
xp.setDetectOptions( detectOptions );
opts.flipy = false ;
opts.timeInterval = 1 ;
opts.timeUnits = '$\tau$' ;
opts.spaceUnits = '$\mu$m' ;
opts.nV = nV ;
opts.nU = nU ;
opts.normalShift = 0 ;
opts.a_fixed = 1.0 ;
opts.adjustlow = 1.00 ;                  %  floor for intensity adjustment
opts.adjusthigh = 99.9 ;                 % ceil for intensity adjustment (clip)
% opts.adjustlow = 0 ;                  %  floor for intensity adjustment
% opts.adjusthigh = 0 ;                 % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;
options.lambda_mesh = 0.0 ;
options.lambda = 0.0 ;
options.lambda_err = 0.0 ;
disp('defining QS')
QS = QuapSlap(xp, opts) ;
disp('done')

%% Build cutMesh timeseries
overwrite = false ;
x = linspace(-1, 1, nU) ;
theta = linspace(0, 2*pi, nV)' ;
xx = (ones(size(theta)) * x)' ;
% Create a dummy y array
y = linspace(0, 1, length(theta)) ;
yy = ones(size(x))' * y ;
uv = [xx(:), yy(:)] ;
uv = [0.5 * uv(:, 1) + 0.5, uv(:, 2)] ;
faces = defineFacesRectilinearGrid(uv, length(x), length(theta)) ;

% The tube has longitudinal profile 1-tmax*gaussian()
tmax = 0.8 ;
sigma = 0.2 ;
mu = 0. ;
for tp = QS.xp.fileMeta.timePoints 
    tidx = QS.xp.tIdx(tp);
    disp(['tp=' num2str(tp)])
    
    % Mesh filenames for output
    cmfn = sprintf(QS.fullFileBase.cutMesh, tp) ;
    meshfn = sprintf(QS.fullFileBase.mesh, tp) ;
    spcmfn = sprintf(QS.fullFileBase.spcutMeshSm, tp) ;
    spcmRSfn = sprintf(QS.fullFileBase.spcutMeshSmRS, tp) ;
    spcMSmRSCfn = sprintf(QS.fullFileBase.spcutMeshSmRSC, tp) ;
        
    if ~exist(cmfn, 'file') || ~exist(meshfn, 'file') || ...
            ~exist(spcmfn, 'file') || overwrite

        tt = tp / length(QS.xp.fileMeta.timePoints) ;
        y = 1 - tt * exp(-(x-mu).^2 / sigma.^2) ;
        radius = ones(size(theta)) * y ;

        yy = (cos(theta) * ones(size(x))) .* radius ;
        zz = (sin(theta) * ones(size(x))) .* radius ;
        yy = yy' ;
        zz = zz' ;
        cutMesh = struct() ;
        cutMesh.v = [xx(:), yy(:), zz(:)] * 100 ;
        cutMesh.u = uv ;
        cutMesh.f = faces ;
        cutMesh.pathPairs = [ (1:nU)', ((nU*(nV-1)+1):nU*nV)' ] ;
        cutMesh.nU = nU ;
        cutMesh.nV = nV ;

        % Plot the cutMesh -- 3D
        subplot(2, 2, 1)
        set(gcf, 'visible', 'on')
        h = trisurf(faces, xx, yy, zz, uv(:, 1), 'EdgeColor', 'none') ; 
        colormap parula
        axis equal
        axis off
        [AO,C,l] = apply_ambient_occlusion(h, 'Factor', 0.4, 'Samples', 1000, 'AddLights', true) ;
        camlight(90,  20) 
        view(2)
        % perspective 3D
        subplot(2, 2, 3:4)
        set(gcf, 'visible', 'on')
        h = trisurf(faces, xx, yy, zz, uv(:, 1), 'EdgeColor', 'none') ; 
        colormap parula
        axis equal
        axis off
        [AO,C,l] = apply_ambient_occlusion(h, 'Factor', 0.4, 'Samples', 1000, 'AddLights', true) ;
        camlight(90,  20) 

        % Plot the cutMesh -- 2D
        subplot(2, 2, 2)
        trisurf(cutMesh.f, uv(:, 1), uv(:, 2), 0*uv(:, 2), uv(:, 1), ...
            'EdgeColor', 'none')
        axis equal
        xlabel('u')
        ylabel('v')
        view(2) ;
        saveas(gcf, [cmfn(1:end-4) '.png'])   

        % Glue cylinder cutMesh into mesh
        mesh = glueCylinderCutMeshSeam(cutMesh) ;
        disp(['Saving glued mesh: ' meshfn])
        plywrite_with_normals(meshfn, mesh.f, mesh.v, mesh.vn)

        % Save cutMesh
        cutMesh.vn = mesh.vn(1:nU*(nV-1), :) ;
        cutMesh.vn(nU*(nV-1)+1:nU*nV, 3) = cutMesh.vn(1:nU, 3) ;
        save(cmfn, 'cutMesh')
        
        % Save cutMesh to smoothed/spcutMeshSm
        spcutMeshSm = cutMesh ;
        save(spcmfn, 'spcutMeshSm') ;
        
        % Save cutMesh to smoothed/spcutMeshSmRS
        spcutMeshSmRS = spcutMeshSm ;
        save(spcmRSfn, 'spcutMeshSmRS') ;
        
        % Save glued cylinder cutMesh into spcMSmRSC
        spcutMeshSmRSC = mesh ;
        save(spcMSmRSCfn, 'spcutMeshSmRSC') ;
    end
end

%% Build 3D velocities -- No Velocity =====================================

overwrite = false ;

pivfn = fullfile(QS.dir.piv, 'piv_results.mat') ;
if ~exist(pivfn, 'file') || overwrite
    pivfntmp = fullfile(QS.dir.piv, 'piv_results_template.mat') ;
    load(pivfntmp, 'piv')
    ntps = length(QS.xp.fileMeta.timePoints) ;
    x = cell(ntps, 1) ;
    y = cell(ntps, 1) ;
    u_original = cell(ntps, 1) ;
    v_original = cell(ntps, 1) ;
    u_filtered = cell(ntps, 1) ;
    v_filtered = cell(ntps, 1) ;
    for tidx = 1:ntps
        tp = QS.xp.fileMeta.timePoints(tidx);
        disp(['tp=' num2str(tp)])

        x{tidx} = piv.x{1} ;
        y{tidx} = piv.y{1} ;
        small = 1e-7 ;
        if strcmpi(velType, 'zero')
            u_filtered{tidx} = small * piv.u_filtered{1} ;
            v_filtered{tidx} = small * piv.v_filtered{1} ;
            u_original{tidx} = small * piv.u_original{1} ;
            v_original{tidx} = small * piv.v_original{1} ;
        elseif strcmpi(velType, 'linearcompression')
            meanx = mean(piv.x{1}(:)) ;
            noiseu = small * (rand(size(piv.u_filtered{1}))-0.5) ;
            noisev = small * (rand(size(piv.v_filtered{1}))-0.5) ;
            maxVel = meanx * 0.5 ;
            minVel = -maxVel ; 
            ramp = 0.1 * (meanx - x{tidx}) ;
            u_filtered{tidx} = ramp + noiseu ;
            v_filtered{tidx} = noisev ;
            u_original{tidx} = ramp + noiseu ;
            v_original{tidx} = noisev ;
        elseif strcmpi(velType, 'linearextension')
            meanx = mean(piv.x{1}(:)) ;
            noiseu = small * (rand(size(piv.u_filtered{1}))-0.5) ;
            noisev = small * (rand(size(piv.v_filtered{1}))-0.5) ;
            maxVel = meanx * 0.5 ;
            minVel = -maxVel ; 
            ramp = max(minVel, min(maxVel, meanx - x{tidx})) ;
            u_filtered{tidx} = ramp + noiseu ;
            v_filtered{tidx} = noisev ;
            u_original{tidx} = ramp + noiseu ;
            v_original{tidx} = noisev ;
        end
    end
    disp(['Saving ' pivfn]) 
    save(pivfn, 'x', 'y', 'u_original', 'v_original', ...
        'u_filtered', 'v_filtered') ;
end

%% Mark fold at t=0
QS.setTime(1) ;
QS.loadCurrentCutMesh() ;
fold_onset = [ 0 ] ;
onest = ones(length(QS.xp.fileMeta.timePoints), 1)  ;
folds = round(nU*0.5) * onest ;
ntps = length(QS.xp.fileMeta.timePoints) ;
ssmax = zeros(ntps, 1) ;
ssfold = zeros(ntps, 1) ;
rssmax = zeros(ntps, 1) ;
rssfold = zeros(ntps, 1) ;
for tidx = 1:ntps
    meshv = QS.currentMesh.cutMesh.v ;
    meshx = meshv(:, 1) ;
    ssmax(tidx) = (max(meshx) - min(meshx)) ;
    ssfold(tidx) =  0.5 * ssmax(tidx) ;
    % Max proper length along surface
    rssmax(tidx) = sum(vecnorm(meshv(2:nU, :) - meshv(1:nU-1, :), 2, 2)) ;
    rssfold(tidx) = 0.5 * sum(vecnorm(meshv(2:folds(tidx), :) ...
        - meshv(1:folds(tidx)-1, :), 2, 2)) ;
end
save(QS.fileName.fold, 'fold_onset', 'folds', 'ssmax', 'ssfold', ...
    'rssmax', 'rssfold')
QS.loadFeatures()

%% Measure velocities =============================================
disp('Making map from pixel to xyz to compute velocities in 3d for smoothed meshes...')
options = struct() ;
options.overwrite = false ;
options.preview = false ;
options.show_v3d_on_data = false ;
options.save_ims = true ;
options.vnscale = vnscale ;
options.vtscale = vtscale ;
QS.measurePIV3d(options) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lagrangian dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pullback pathline time averaging of velocities
options = struct() ;
options.overwrite = false ;
options.preview = false ;
options.twidth = 0 ;
options.plotOptions = struct() ;
options.plotOptions.vnscale = vnscale ;
options.plotOptions.vtscale = vtscale ;
QS.timeAverageVelocities(options)
%% Velocity plots for pathline time averaging 
options.overwrite = true ;
options.plot_vxyz = false ;
options.invertImage = true ;
options.averagingStyle = 'Lagrangian'; 
options.samplingResolution = '1x'; 
options.vnscale = vnscale ;
options.vtscale = vtscale ;
QS.plotTimeAvgVelocities(options)
%% Divergence and Curl (Helmholtz-Hodge) for Lagrangian
options = struct() ;
options.overwrite = true ;
options.samplingResolution = '1x' ;
options.averagingStyle = 'Lagrangian' ;
options.clipDiv = 1e9 ;
options.clipRot = 1e9 ;
options.sscaleDiv = 3.0 ;
options.sscaleRot = 0.5 ;
QS.helmholtzHodge(options) ;

%% Compressibility & kinematics for Lagrangian
options = struct() ;
options.overwrite = true ;
options.samplingResolution = '1x'; 
QS.measureMetricKinematics(options)
%% Metric Kinematics Kymographs & Correlations
options = struct() ;
options.overwrite = true ;
options.overwrite_timePoints = true ;
options.plot_Hgdot = false ;
options.plot_flows = true ;
options.plot_factors = true ;
options.plot_kymographs = true ;
options.plot_kymographs_cumsum = true ;
options.plot_kymographs_cumprod = true ;
options.plot_correlations = true ;
options.plot_raw_correlations = true ;
options.plot_gdot_correlations = false ;
options.plot_gdot_decomp = false ;
options.climit = 0.2 * vnscale;
options.climit_err = 0.2 * vnscale ;
options.climit_veln = vnscale;
options.climit_H = 0.5;
options.climit_radius = 0 ;
QS.plotMetricKinematics(options)

%% Pullback pathlines connecting Lagrangian grids
options = struct() ;
options.overwrite = true ;
options.preview = false ;
options.debug = false ; 
QS.measurePullbackPathlines(options)
% Query velocities along pathlines
options = struct() ;
options.overwrite = true ;
options.preview = false ;
QS.measurePathlineVelocities(options)
% plot the pathline velocities 
options = struct() ;
options.overwrite = true ;
options.gridTopology = 'triangulated' ;
QS.plotPathlineVelocities(options)

%% Measure Pathline Kinematics
options = struct() ;
options.overwrite = true ;
QS.measurePathlineMetricKinematics(options)
%% Plot Pathline Kinematics
options = struct() ;
options.overwrite = true ;
options.plot_kymographs = true ;
options.plot_kymographs_cumsum = true ;
options.plot_kymographs_cumprod = false ;
options.plot_correlations = false ;
options.plot_fold_kinematics = true ;
options.plot_lobe_kinematics = true ;
options.climit = 0.10 ;
options.plot_ap = true ;
options.plot_left = false ;
options.plot_right = false ;
options.plot_dorsal = false ;
options.plot_ventral = false ;
options.climit = 0.2 * vnscale;
options.climit_err = 0.2 * vnscale ;
options.climit_veln = vnscale;
options.climit_H = 0.5;
options.climit_radius = 0 ;
QS.plotPathlineMetricKinematics(options)
%% Strain rate (epsilon = 1/2 (djvi+divj) -vn bij)
options = struct() ;
options.overwrite = true ;
options.overwriteImages = true ;
options.preview = false ;
QS.measureStrainRate(options) 
%% Kymograph strain rates
options = struct() ;
options.overwrite = true ;
options.skipTimePoint = true ;
options.clim_trace = 0.05 ;
options.clim_deviatoric = 0.05 ;
QS.plotStrainRate(options)

%% Measure strain rate along pathlines
options = struct() ;
options.overwrite = false ;
options.overwriteImages = true ;
options.plot_dzdp = false ;
QS.measurePathlineStrainRate(options)
% Measure strain along pathlines
options = struct() ;
options.overwrite = true ;
options.overwriteImages = true ;
options.plot_dzdp = false ;
options.median_filter_strainRates = false ;
options.climitInitial = 0.05 ;
options.climitRamp = 0.01 ;
options.climitRatio = 1 ;
QS.measurePathlineStrain(options)
% Pathline strain rate plots
options = struct() ;
options.overwrite = true ;
options.plot_kymographs = false ;
options.plot_kymographs_strain = true ;
options.plot_fold_strainRate = false ;
options.plot_lobe_strainRate = false ;
options.plot_fold_strain = true ;
options.plot_lobe_strain = true ;
options.climit = 0.05 ;
options.climitWide = 1.0 ;
QS.plotPathlineStrainRate(options)
