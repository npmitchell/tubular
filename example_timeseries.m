%% EXAMPLE MASTER PIPELINE FOR TIME SERIES DATA (3D + time)
% by NPMitchell & Dillon Cislo

% TO DO:
% ------
%  - redefine APDV coords as viewingCoordinates (global ref frame)
%  - install instructions for morphsnakes as an option
%  - iLastik instructions (optional)
%  - imsane instructions (optional)
%  - do we need gptoolbox? so far just for laplacian_smooth()
%  - DEC tutorials --> STANDALONE AND INTEGRATED TUTORIALS
%  - Ricci flow tutorials or remove this

% This is a pipeline to analyze dynamic tube-like surfaces in 3D data.
% A tube-like surface is one that is either cylindrical in topology or
% elongated and spherical in topology. The pipeline will work for spheres,
% but two 'endcaps' will be truncated in order to transform it to a
% cylindrical topology.


%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
% cd /mnt/data/organoids/tube/
% cd /mnt/data/code/tubular/example/ ;
cd /mnt/data/tubular_test/fly_midgut/ ;

dataDir = cd ;

%% ADD PATHS TO THIS ENVIRONMENT ==========================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
addpath(fileparts(origpath))
addpath(fullfile('utility', 'addpath_recurse'))
addpath_recurse('utility')
addpath_recurse('/mnt/data/code/gptoolbox')
addpath('TexturePatch')
addpath('DEC')
addpath(fullfile('utility','plotting'))
addpath(fullfile('utility','plotting'))
% go back to the data
cd(dataDir)

%% DEFINE NEW MASTER SETTINGS
if ~exist('./masterSettings.mat', 'file')
    % Metadata about the experiment
    stackResolution = [.2619 .2619 .2619] ;
    nChannels = 1 ;
    channelsUsed = 1 ;
    timePoints = 1:81; %86:211 ;
    ssfactor = 4 ;
    % whether the data is stored inverted relative to real position
    flipy = false ; 
    timeInterval = 2 ;  % physical interval between timepoints
    timeUnits = 'min' ; % physical unit of time between timepoints
    spaceUnits = '$\mu$m' ; % physical unit of time between timepoints
    scale = 1.0 ;      % scale for conversion to 16 bit
    file32Base = 'TP%d_Ch0_Ill0_Ang0,45,90,135,180,225,270,315.tif'; 
    fn = 'Time_%06d_c1_stab';
    fn_prestab = 'Time_%06d_c1.tif';
    set_preilastikaxisorder = 'xyzc' ;
    swapZT = 1 ;
    masterSettings = struct('stackResolution', stackResolution, ...
        'nChannels', nChannels, ...
        'channelsUsed', channelsUsed, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeInterval', timeInterval, ...
        'timeUnits', timeUnits, ...
        'spaceUnits', spaceUnits, ...
        'scale', scale, ...
        'file32Base', file32Base, ...
        'fn', fn,...
        'fn_prestab', fn_prestab, ...
        'swapZT', swapZT, ...
        'set_preilastikaxisorder', set_preilastikaxisorder, ...
        'nU', 150, ...  % 150 for mef2 data with posterior midgut loop
        'nV', 100); 
    disp('Saving masterSettings to ./masterSettings.mat')
    if exist('./masterSettings.mat', 'file')
        ui = input('This will overwrite the masterSettings. Proceed (Y/n)?', 's') ;
        if ~isempty(ui) && (strcmp(ui(1), 'Y') || strcmp(ui(1), 'y'))
            save('./masterSettings.mat', 'masterSettings')
            loadMaster = false ;
        else
            disp('Loading masterSettings from disk instead of overwriting')
            loadMaster = true ;
        end
    else
        save('./masterSettings.mat', 'masterSettings')
        loadMaster = false ;
    end
else
    loadMaster = true ;
end

if loadMaster
    % LOAD EXISTING MASTER SETTINGS
    disp('Loading masterSettings from ./masterSettings.mat')
    load('./masterSettings.mat', 'masterSettings')
    % Unpack existing master settings
    stackResolution = masterSettings.stackResolution ;
    nChannels = masterSettings.nChannels ;
    channelsUsed = masterSettings.channelsUsed ;
    timePoints = masterSettings.timePoints ;
    ssfactor = masterSettings.ssfactor ;
    % whether the data is stored inverted relative to real position
    flipy = masterSettings.flipy ; 
    timeInterval = masterSettings.timeInterval ;  % physical interval between timepoints
    timeUnits = masterSettings.timeUnits ; % physical unit of time between timepoints
    spaceUnits = masterSettings.spaceUnits ; % unit of distance of full resolution data pixels ('$\mu$m')
    fn = masterSettings.fn ;
    set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
    swapZT = masterSettings.swapZT ;
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;
end
dir16bit = fullfile(dataDir) ;

%% I. INITIALIZE ImSAnE PROJECT ===========================================
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
cd(dir16bit)
dataDir = cd ;
projectDir = dataDir ;
% [ projectDir, ~, ~ ] = fileparts(matlab.desktop.editor.getActiveFilename); 
cd(projectDir);
if projectDir(end) ~= '/'
    projectDir = [projectDir '/'];
end

% Start by creating an experiment object, optionally pass on the project
% directory (otherwise it will ask), and change into the directory of the
% data.  This serves as a front-end for data loading, detection, fitting
% etc.
xp = project.Experiment(projectDir, dataDir);

% Set file and experiment meta data
%
% We assume on individual image stack for each time point, labeled by time.
%  To be able to load the stack, we need to tell the project wehre the data
%  is, what convention is assumed for the file names, available time
%  points, and the stack resolution.  Options for modules in ImSAnE are
%  organized in MATLAB structures, i.e a pair of field names and values are
%  provided for each option.
%
% The following file metadata information is required:
% * 'directory'         , the project directory (full path)
% * 'dataDir'           , the data directory (full path)
% * 'filenameFormat'    , fprintf type format spec of file name
% * 'timePoints'        , list of itmes available stored as a vector
% * 'stackResolution'   , stack resolution in microns, e.g. [0.25 0.25 1]
%
% The following file metadata information is optional:
% * 'imageSpace'        , bit depth of image, such as uint16 etc., defined
%                         in Stack class
% * 'stackSize'         , size of stack in pixels per dimension 
%                         [xSize ySize zSize]
% * 'swapZT'            , set=1 if time is 3rd dimension and z is 4th

% A filename base template - to be used throughout this script
fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [fn, '.tif'];
fileMeta.nChannels          = nChannels;
fileMeta.timePoints         = timePoints ;
fileMeta.stackResolution    = stackResolution;
fileMeta.swapZT             = masterSettings.swapZT;

% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating  the sample, which time point
% to use for fitting, etc.
%
% The following project metadata information is required:
% * 'channelsUsed'      , the channels used, e.g. [1 3] for RGB
% * 'channelColor'      , mapping from element in channels used to RGB = 123
% * 'dynamicSurface'    , Not implemented yet, future plan: boolean, false: static surface
% * 'detectorType'      , name of detector class, e.g. radielEdgeDetector
%                         ,(user threshholded), fastCylinderDetector
% * 'fitterType'        , name of fitter class
%
% The following project meta data information is optional:
% * 'description'     , string describing the data set set experiments metadata, 
%                                such as a description, and if the surface is dynamic,
%                                or requires drift correction of the sample.
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 

% first_tp is also required, which sets the tp to do individually.
first_tp = 1 ;
expMeta                     = struct();
expMeta.channelsUsed        = channelsUsed ;
expMeta.channelColor        = 1;
expMeta.description         = 'Drosophila gut';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

% Now set the meta data in the experiment.
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

clear fileMeta expMeta

%% LOAD THE FIRST TIME POINT ==============================================
xp.setTime(xp.fileMeta.timePoints(1)) ;
% xp.loadTime(xp.fileMeta.timePoints(first_tp));
% xp.rescaleStackToUnitAspect();

%% SET DETECTION OPTIONS ==================================================
% Load/define the surface detection parameters
msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
if exist(msls_detOpts_fn, 'file') 
    load(msls_detOpts_fn, 'detectOptions')
else
    channel = 1;
    foreGroundChannel = 1;
    zdim = 2 ;
    ssfactor=4 ;
    niter=35 ;
    niter0=35 ;
    ofn_ply='mesh_ms_' ;
    ofn_ls='msls_' ;
    ofn_smoothply = 'mesh_' ;
    ms_scriptDir='/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/' ;
    pre_nu=-5 ;
    pre_smoothing=0 ;
    lambda1=1 ;
    lambda2=1 ;
    exit_thres=0.000005 ;
    smoothing=0.5 ;
    nu=0 ;
    post_nu=2 ;
    post_smoothing=3 ;    
    init_ls_fn = 'msls_initguess.h5' ;
    mlxprogram = fullfile(meshlabCodeDir, ...
        'laplace_surface_rm_resample30k_reconstruct_LS3_1p2pc_ssfactor4.mlx') ;
    radius_guess = 40 ;
    center_guess = '200,75,75' ;
    dtype = 'h5' ;
    mask = 'none' ;
    prob_searchstr = '_stab_Probabilities.h5' ;
    preilastikaxisorder= set_preilastikaxisorder; ... % axis order in input to ilastik as h5s. To keep as saved coords use xyzc
    ilastikaxisorder= 'cxyz'; ... % axis order as output by ilastik probabilities h5
    imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
    include_boundary_faces = true ;
    smooth_with_matlab = -1;
    
    % Name the output mesh directory --------------------------------------
    mslsDir = [fullfile(projectDir, 'msls_output') filesep];

    % Surface detection parameters ----------------------------------------
    detectOptions = struct( 'channel', channel, ...
        'ssfactor', ssfactor, ...
        'niter', niter,...
        'niter0', niter0, ...
        'lambda1', lambda1, ...
        'lambda2', lambda2, ...
        'nu', nu, ...
        'smoothing', smoothing, ...
        'post_nu', post_nu, ...
        'post_smoothing', post_smoothing, ...
        'exit_thres', exit_thres, ...
        'foreGroundChannel', foreGroundChannel, ...
        'fileName', sprintf( fn, xp.currentTime ), ...
        'mslsDir', mslsDir, ...
        'ofn_ls', ofn_ls, ...
        'ofn_ply', ofn_ply,...
        'ms_scriptDir', ms_scriptDir, ...
        'timepoint', xp.currentTime, ...
        'zdim', zdim, ...
        'ofn_smoothply', ofn_smoothply, ...
        'pre_nu', pre_nu, ...
        'pre_smoothing', pre_smoothing, ...
        'mlxprogram', mlxprogram, ...
        'init_ls_fn', init_ls_fn, ... % set to none to load prev tp
        'run_full_dataset', run_full_dataset,... % projectDir, ... % set to 'none' for single tp
        'radius_guess', radius_guess, ...
        'dset_name', 'exported_data',...
        'center_guess', center_guess,... % xyz of the initial guess sphere ;
        'save', true, ... % whether to save images of debugging output
        'plot_mesh3d', false, ...
        'dtype', dtype,...
        'mask', mask,...
        'mesh_from_pointcloud', false, ...
        'prob_searchstr', prob_searchstr, ...
        'preilastikaxisorder', preilastikaxisorder, ... 
        'ilastikaxisorder', ilastikaxisorder, ... 
        'physicalaxisorder', imsaneaxisorder, ... 
        'include_boundary_faces', include_boundary_faces, ...
        'smooth_with_matlab', smooth_with_matlab, ...
        'pythonVersion', '3') ;

    % save options
    if exist(msls_detOpts_fn, 'file')
        disp('Overwriting detectOptions --> renaming existing as backup')
        backupfn1 = [msls_detOpts_fn '_backup1'] ;
        if exist(backupfn1, 'file')
            backupfn2 = [msls_detOpts_fn '_backup2'] ; 
            system(['mv ' backupfn1 ' ' backupfn2])
        end
        system(['mv ' msls_detOpts_fn ' ' backupfn1])
    end
    disp('Saving detect Options to disk')
    save(msls_detOpts_fn, 'detectOptions') ;
end

% Overwrite certain parameters for script structure
detectOptions.fileName = sprintf( fn, xp.currentTime ) ;
mslsDir = detectOptions.mslsDir ;

% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );
disp('done')

%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
for tt = xp.fileMeta.timePoints
    if ~exist(fullfile(projectDir, [sprintf(fn, tt) '.h5']), 'file')
        disp(['Did not find file: ', fullfile(projectDir, [sprintf(fn, tt) '.h5'])])
        xp.loadTime(tt);
        xp.rescaleStackToUnitAspect();
        % make a copy of the detectOptions and change the fileName
        detectOpts2 = detectOptions ;
        detectOpts2.fileName = sprintf( fn, xp.currentTime ) ;
        xp.setDetectOptions( detectOpts2 );
        xp.detector.prepareIlastik(xp.stack);
        disp(['done outputting downsampled data h5: tp=' num2str(tt) ' for surface detection'])
    else
        disp(['h5 ' num2str(tt) ' was already output, skipping...'])
    end
end    
disp('Open with ilastik if not already done')

%% TRAIN NON-STABILIZED DATA IN ILASTIK TO IDENTIFY SURFACE ==============
% Open ilastik, train pre-stab h5s until probabilities and uncertainty are 
% satisfactory, then run on stab images.

%% Create MorphSnakes LevelSets from the Probabilities output of ilastik ==
% Now detect all surfaces
if strcmp(detectOptions.run_full_dataset, projectDir)
    disp('Running dataset mode, in which all surfaces are extracted serially using morphsnakes')
    xp.setTime(xp.fileMeta.timePoints(1));
    detectOpts2 = detectOptions ;
    detectOpts2.fileName = sprintf( fn, xp.currentTime ) ;
    xp.setDetectOptions( detectOpts2 );
    xp.detectSurface();
else
    % Use morphsnakes to extract surfaces for all timepoints INDIVIDUALLY 
    assert(strcmp(detectOptions.run_full_dataset, 'none'))
    for tp = xp.fileMeta.timePoints
        xp.setTime(tp);
        
        % make a copy of the detectOptions and change the fileName
        detectOpts2 = detectOptions ;
        detectOpts2.timepoint = xp.currentTime ;
        detectOpts2.fileName = sprintf( fn, xp.currentTime );
        xp.setDetectOptions( detectOpts2 );
        xp.detectSurface();
        
        % For next time, use the output mesh as an initial mesh, which will
        % be searched for if init_ls_fn is set to 'none'.
        detectOpts2.init_ls_fn = 'none' ;
    end
end

%% Now we have 3d data volumes and surfaces. Define a TubULAR object. 
% To visualize data on these surfaces and compute how these surfaces deform
% we now define TubULAR object.
nU = masterSettings.nU ;
nV = masterSettings.nV ;
opts.meshDir = mslsDir ;
opts.flipy = flipy ;
opts.timeInterval = timeInterval ;
opts.timeUnits = timeUnits ;
opts.spaceUnits = spaceUnits ;
opts.nU = nU ;
opts.nV = nV ;
opts.t0 = xp.fileMeta.timePoints(1) ;
opts.normalShift = 10 ;
opts.a_fixed = 2.0 ;
opts.adjustlow = 1.00 ;         % floor for intensity adjustment
opts.adjusthigh = 99.9 ;        % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;
opts.lambda_mesh = 0.002 ;
opts.lambda = 0.01 ;
opts.lambda_err = 0.01 ;
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done')

%% Inspect all meshes in 3D
for tp = xp.fileMeta.timePoints
    % Load the mesh
    meshfn = sprintf( tubi.fullFileBase.mesh, tp ) ;    
    mesh = read_ply_mod(meshfn) ;
    % Plot the mesh in 3d. Color here by Y coordinate
    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        mesh.v(:, 3), 'edgecolor', 'none', 'Facealpha', 0.5)
    title(['t=' num2str(tp)])
    axis equal
    view(2)
    pause(0.1)
end

%% Obtain APDV coordinates of the surface. 
% There are two options for obtaining these coordinates. 
%   1. Automatically determine A and P by the extremal points of the
%   surface mesh along the elongated axis of the mesh, and define DV as
%   pointing perpendicular to this.
%   2. Train in iLastik for an anterior spot in 3d (A), a posterior spot in
%   3d (P), and a spot which is dorsal to the line connecting A and P. Any
%   dorsal point is fine, as long as it points dorsal to the AP axis
%   defined by A and P. See picture below.
%
%     example:
%                 D
%                 |
%                 |
%       A -----------------P
%
% Here we use option 2. We must prepare APDV ilastik training first outside
% MATLAB.
% Train on anterior (A), posterior (P), background (B), and 
% dorsal anterior (D) location in different iLastik channels by having a
% blob centered on a point that you wish to identify as each label (in 3D).
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
% Training for dorsal (D) is only needed at the reference time point, t0,
% because that's the only one that's used. 
%
% A dorsal blob for the gut is marked at the site where the gut closes,
% with 48YGAL4-expressing cells form a seam.
% Posterior is at the rear of the yolk, where the endoderm closes, for 
% apical surface training.
% Anterior is at the junction of the midgut with the foregut.

%% 3. align_meshes_APDV or load transformations if already done
% Compute APDV coordinate system
alignAPDVOpts = struct() ;
tubi.computeAPDVCoords(alignAPDVOpts) ;

% Compute the APD COMs
apdvOpts = struct() ;
[acom_sm, pcom_sm] = tubi.computeAPDCOMs(apdvOpts) ;

% Align the meshes APDV & plot them
tubi.alignMeshesAPDV(alignAPDVOpts) ;

disp('done')

%% PLOT ALL TEXTURED MESHES IN 3D =========================================

% Establish texture patch options
% Get limits and create output dir
metadat = struct() ;
metadat.reorient_faces = false ;            % set to true if some mesh normals may be inverted (requires gptoolbox if true)
metadat.normal_shift = tubi.normalShift ;   % normal push, in pixels, along normals defined in data XYZ space
metadat.texture_axis_order = [1 2 3] ;      % texture space sampling. If the surface and dataspace have axis permutation, enter that here

% Psize is the linear dimension of the grid drawn on each triangular face
Options.PSize = 5 ;
Options.EdgeColor = 'none';
Options.numLayers = [0, 0];  % how many layers to MIP over/bundle into stack, as [outward, inward]
Options.layerSpacing = 2 ;    % layers are 2 pixels apart

% Plot on surface for all TP 
options = metadat ;
tubi.plotSeriesOnSurfaceTexturePatch(options, Options)

%% EXTRACT CENTERLINES
% Skip if already done
% Note: these just need to be 'reasonable' centerlines for topological
% checks on the orbifold cuts. Therefore, use as large a resolution ('res')
% as possible that still forms a centerline passing through the mesh
% surface, since the centerline computed here is just for constraining the 
% mapping to the plane.
cntrlineOpts.overwrite = false ;         % overwrite previous results
cntrlineOpts.overwrite_ims = false ;     % overwrite previous results
cntrlineOpts.weight = 0.1;               % for speedup of centerline extraction. Larger is less precise
cntrlineOpts.exponent = 1.0 ;            % how heavily to scale distance transform for speed through voxel
cntrlineOpts.res = 4.0 ;                 % resolution of distance tranform grid in which to compute centerlines
cntrlineOpts.preview = false ;           % preview intermediate results
cntrlineOpts.reorient_faces = false ;    % not needed for our well-constructed meshes
cntrlineOpts.dilation = 0 ;              % how many voxels to dilate the segmentation inside/outside before path computation
% Note: this can take about 400s per timepoint for res=2.0, so use as big a 
%   res value as possible.
%
tubi.extractCenterlineSeries(cntrlineOpts)
disp('done with centerlines')

%% Identify anomalies in centerline data
idOptions.ssr_thres = 15 ;  % distance of sum squared residuals in um as threshold
tubi.generateCleanCntrlines(idOptions) ;
disp('done with cleaning up centerlines')

%% Cylinder cut mesh --> transforms a topological sphere into a topological cylinder
% Look for options on disk. If not saved, define options.
if ~exist(tubi.fileName.endcapOptions, 'file')
    endcapOpts = struct( 'adist_thres', 20, ...  % 20, distance threshold for cutting off anterior in pix
                'pdist_thres', 14, ...  % 15-20, distance threshold for cutting off posterior in pix
                'tref', tubi.t0) ;  % reference timepoint at which time dorsal-most endcap vertices are defined
    tubi.setEndcapOptions(endcapOpts) ;
    % Save the options to disk
    tubi.saveEndcapOptions() ;
else
    % load endcapOpts
    tubi.loadEndcapOptions() ;
    endcapOpts = tubi.endcapOptions ;
end

methodOpts.save_figs = true ;   % save images of cutMeshes along the way
methodOpts.preview = false  ;     % display intermediate results
tubi.sliceMeshEndcaps(endcapOpts, methodOpts) ;

%% Clean Cylinder Meshes
% This removes "ears" from the endcaps of the tubular meshes (cylindrical
% meshes)
cleanCylOptions = struct() ;
tubi.cleanCylMeshes(cleanCylOptions)
disp('done cleaning cylinder meshes')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORBIFOLD -> begin populating tubi.dir.mesh/gridCoords_nUXXXX_nVXXXX/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate Through Time Points to Create Pullbacks ========================
for tt = tubi.xp.fileMeta.timePoints
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
    
    %----------------------------------------------------------------------
    % Create the Cut Mesh
    %----------------------------------------------------------------------
    cutMeshfn = sprintf(tubi.fullFileBase.cutMesh, tt) ;
    cutPathfn = sprintf(tubi.fullFileBase.cutPath, tt) ;
    if ~exist(cutMeshfn, 'file') || ~exist(cutPathfn, 'file')
        if exist(cutMeshfn, 'file')
            disp('Overwriting cutMesh...') ;
        else
            disp('cutMesh not found on disk. Generating cutMesh... ');
        end
        options = struct() ;
        tubi.generateCurrentCutMesh(options)
        disp('Saving cutP image')
        % Plot the cutPath (cutP) in 3D
        tubi.plotCutPath(tubi.currentMesh.cutMesh, tubi.currentMesh.cutPath)
        compute_pullback = true ;
    else
        fprintf('Loading Cut Mesh from disk... ');
        tubi.loadCurrentCutMesh()
        compute_pullback = ~isempty(tubi.currentMesh.cutPath) ;
    end
    
    spcutMeshOptions = struct() ;
    spcutMeshOptions.t0_for_phi0 = tubi.t0set() ;  % which timepoint do we define corners of pullback map
    spcutMeshOptions.save_phi0patch = false ;
    spcutMeshOptions.iterative_phi0 = false ;
    spcutMeshOptions.smoothingMethod = 'none' ;
    tubi.plotting.preview = false ;
    tubi.generateCurrentSPCutMesh([], spcutMeshOptions) ;
    
    % Compute the pullback if the cutMesh is ok
    if compute_pullback || ~exist(sprintf(tubi.fullFileBase.im_sp, tt), 'file')
        tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
    else
        disp('Skipping computation of pullback')
    end
    clear Options IV
        
end
disp('Done with generating spcutMeshes and cutMeshes')

%% Inspect coordinate system charts using (s,phi) coordinate system ('sp')
options = struct() ;
options.coordSys = 'sp' ;
tubi.coordSystemDemo(options)

%% OPTIONAL: COMPUTE MESH SURFACE AREA AND VOLUME =========================
options = struct() ;
tubi.measureSurfaceAreaVolume(options)
disp('done')

%% OPTIONAL: COMPUTE WRITHE OF MEANCURVE CENTERLINES ======================
options = struct() ;
tubi.measureWrithe(options)
disp('done')

%% OPTIONAL: Plot fancy "cross-section" view of centerlines
options = struct() ;
tubi.plotClineXSections(options)

%% Smooth the sphi grid meshes in time ====================================
options = struct() ;
options.width = 4 ;  % width of kernel, in #timepoints, to use in smoothing meshes
tubi.smoothDynamicSPhiMeshes(options) ;

%% Plot the time-smoothed meshes
tubi.plotSPCutMeshSmRS(options) ;

% Inspect coordinate system charts using smoothed meshes
options = struct()
options.coordSys = 'spsm' ;
tubi.coordSystemDemo(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Redo Pullbacks with time-smoothed meshes ===============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Skip if already done
disp('Create pullback using S,Phi coords with time-averaged Meshes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = tubi.xp.fileMeta.timePoints
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = tubi.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    tubi.setTime(tt) ;
    
    % Establish custom Options for MIP --> choose which pullbacks to use
    pbOptions = struct() ;
    pbOptions.numLayers = [0 0] ; % how many onion layers over which to take MIP
    pbOptions.generate_spsm = true ;
    pbOptions.generate_sp = false ;
    tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
end

%% TILE/EXTEND SMOOTHED IMAGES IN Y AND RESAVE =======================================
% Skip if already done
options = struct() ;
options.coordsys = 'spsm' ;
tubi.doubleCoverPullbackImages(options)
disp('done')

%% PERFORM PIV ON PULLBACK MIPS ===========================================
% % Compute PIV in PIVLab
% % ---------------------
% % Open PIVLab
% % Select all frames in meshDir/PullbackImages_XXXstep_sphi/smoothed_extended/
% % Select Sequencing style 1-2, 2-3, ... 
% %
% % Below are settings in PIVlab that work well for this dataset:
% % Image Preprocessing (used to select all, but now:)
% %  --> Enable CLAHE with 20 pix
% %  --> DO NOT Enable highpass with 15 pix
% %  --> DO NOT Enable Intensity capping
% %  --> Wiener2 denoise filter with 3 pix
% %  --> DO NOT Auto constrast stretch
% % PIV settings: 
% %  --> 128 (32 step), 64 (32 step), 32 (16 step), 16 (8 step)
% %  --> Linear window deformation interpolator
% %  --> 5x repeated correlation 
% %  --> Disable auto-correlation
% % Post-processing
% %  --> Standard deviation filter: 7 stdev
% %  --> Local median filter: thres=5, eps=0.1
% %  --> Interpolate missing data
% % Export 
% %  --> File > Save > MAT file

%% Measure velocities =============================================
disp('Making map from pixel to xyz to compute velocities in 3d for smoothed meshes...')
options = struct() ;
options.show_v3d_on_data = false ;
tubi.measurePIV3d(options) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lagrangian dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pullback pathline time averaging of velocities
options = struct() ;
tubi.timeAverageVelocities(options)
% Velocity plots for pathline time averaging 
options.plot_vxyz = false ;
options.invertImage = true ;
options.averagingStyle = 'Lagrangian'; 
tubi.plotTimeAvgVelocities(options)
% Divergence and Curl (Helmholtz-Hodge) for Lagrangian
options = struct() ;
options.averagingStyle = 'Lagrangian' ;
options.lambda = 0 ;
options.lambda_mesh = 0 ; 
tubi.helmholtzHodge(options) ;

% Compressibility & kinematics for Lagrangian
options = struct() ;
tubi.measureMetricKinematics(options)

%% Metric Kinematics Kymographs & Correlations -- Bandwidth Filtered
options = struct() ;
tubi.plotMetricKinematics(options)

%% Pullback pathlines connecting Lagrangian grids
options = struct() ;
tubi.measurePullbackPathlines(options)

%% Query velocities along pathlines
options = struct() ;
tubi.measurePathlineVelocities(options)
% plot the pathline velocities 
options = struct() ;
options.gridTopology = 'triangulated' ;
tubi.plotPathlineVelocities(options)

% Measure Pathline Kinematics
options = struct() ;
tubi.measurePathlineMetricKinematics(options)

% Plot Pathline Kinematics
options = struct() ;
tubi.plotPathlineMetricKinematics(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create ricci mesh at t0 to measure Beltrami coefficient in pathlines
options = struct() ;
options.climit = 1 ;
options.coordSys = 'ricci' ;
tubi.measureBeltramiCoefficient(options) ;

%% Generate all Beltramis from all Riccis & plot aspect ratio over time
options = struct() ;
tubi.computeRicciMeshes(options) 

%% Strain rate (epsilon = 1/2 (djvi+divj) -vn bij)
options = struct() ;
tubi.measureStrainRate(options) 

%% Plot time-averaged strain rates in 3d on mesh
options = struct() ;
tubi.plotStrainRate3DFiltered(options) 

%% Kymograph strain rates
options = struct() ;
options.clim_trace = 0.05 ;
options.clim_deviatoric = 0.05 ;
tubi.plotStrainRate(options)

% Measure strain rate along pathlines
options = struct() ;
options.overwriteImages = false ;
options.plot_dzdp = false ;
tubi.measurePathlineStrainRate(options)

% Pathline strain rate plots
options = struct() ;
options.plot_kymographs = false ;
options.climit = 0.05 ;
options.climitWide = 1.0 ;
tubi.plotPathlineStrainRate(options)

%% Measure strain along pathlines -- note this is from pathlines, not integrating rates
options = struct() ;
options.plot_dzdp = false ;
options.climitInitial = 0.05 ;
options.climitRamp = 0.01 ;
options.climitRatio = 1 ;
tubi.measurePathlineStrain(options)
tubi.plotPathlineStrain(options)

%% Measure coarse-grained bond contraction and dilation in zeta=s/L, phi
% todo: consider putting this in Ricci map frame instead of (s,phi) frame
options = struct() ;
tubi.measureDxDyStrainFiltered(options) ;


%% Extra functionality to demonstrate: PULLBACK STACKS ====================
% Skip if already done
disp('Create pullback stack using S,Phi coords with time-averaged Meshes');
% Load options
overwrite = true ;
optionfn = fullfile(tubi.dir.im_r_sme_stack, 'spcutMeshSmStackOptions.mat') ;
if ~exist(optionfn, 'file') || overwrite
    spcutMeshSmStackOptions.layer_spacing = 0.5 / tubi.APDV.resolution ; % pixel resolution roughly matches xy
    spcutMeshSmStackOptions.n_outward = 20 ;
    spcutMeshSmStackOptions.n_inward = 40 ;
    spcutMeshSmStackOptions.smoothIter = 0 ;
    spcutMeshSmStackOptions.preSmoothIter = 35 ;
    spcutMeshSmStackOptions.imSize = 500 ;
    % Save options
    save(optionfn, 'spcutMeshSmStackOptions')
else
    load(optionfn, 'smSPCutMeshStackOptions')
end
spcutMeshSmStackOptions.overwrite = overwrite ;
tubi.generateSPCutMeshSmStack(spcutMeshSmStackOptions)
