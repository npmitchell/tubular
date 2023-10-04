function makeH5SeriesPrestabForTraining(masterSettings)
% NPMitchell 
%
% This is a pipeline to create h5s before stabilization, to be used on all
% training.

%% Unpack arguments
if strcmp(masterSettings.fn_prestab(end-3:end), '.tif')
    fn = masterSettings.fn_prestab(1:end-4) ;
else
    fn = masterSettings.fn_prestab ;
end
nChannels = masterSettings.nChannels ;
channelsUsed = masterSettings.channelsUsed ;    
timePoints = masterSettings.timePoints ;
ssfactor = masterSettings.ssfactor ;
stackResolution = masterSettings.stackResolution ;
dataDir = masterSettings.dir16bit_prestab ;
preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
swapZT = masterSettings.swapZT ;

%% I. INITIALIZE

%% INITIALIZE ImSAnE PROJECT ==============================================
% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
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
% Set required additional information on the files.
%
% We assume on individual image stack for each time point, labeled by time.
%  To be able to load the stack, we need to tell the project wehre the data
%  is, what convention is assumed for the file names, available time
%  points, and the stack resolution.  Options for modules in ImSAnE are
%  organized in MATLAB structures, i.e a pair of field names and values are
%  provided for each option.
%
% The following file metadata information is required:
%
% * 'directory'         , the project directory (full path)
% * 'dataDir'           , the data directory (full path)
% * 'filenameFormat'    , fprintf type format spec of file name
% * 'timePoints'        , list of itmes available stored as a vector
% * 'stackResolution'   , stack resolution in microns, e.g. [0.25 0.25 1]
%
% The following file metadata information is optional:
%
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
fileMeta.swapZT             = swapZT;

% Set required additional information on the experiment. A verbal data set
% description, Jitter correct by translating  the sample, which time point
% to use for fitting, etc.
%
% The following project metadata information is required:
%
% * 'channelsUsed'      , the channels used, e.g. [1 3] for RGB
% * 'channelColor'      , mapping from element in channels used to RGB = 123
% * 'dynamicSurface'    , Not implemented yet, future plan: boolean, false: static surface
% * 'detectorType'      , name of detector class, e.g. radielEdgeDetector
%                         ,(user threshholded), fastCylinderDetector
% * 'fitterType'        , name of fitter class
%
% The following project meta data information is optional:
%
% * 'description'     , string describing the data set set experiments metadata, 
%                                such as a description, and if the surface is dynamic,
%                                or requires drift correction of the sample.
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 

% first_tp is also required, which sets the tp to do individually.
first_tp = timePoints(1) ;
expMeta                     = struct();
expMeta.channelsUsed        = channelsUsed ;
expMeta.channelColor        = 1:length(channelsUsed) ;
expMeta.description         = 'Apical membrane in Drosophila gut';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = first_tp;
expMeta.detectorType        = 'surfaceDetection.morphsnakesDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

% Now set the meta data in the experiment.
xp.setFileMeta(fileMeta);
xp.setExpMeta(expMeta);
xp.initNew();

clear fileMeta expMeta

%% LOAD THE FIRST TIME POINT ==============================================
xp.setTime(xp.fileMeta.timePoints(1));

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other options
%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: difference from CAAX
mlxprogram = 'surface_rm_resample20k_reconstruct_LS3_ssfactor4.mlx';
mslsDir = projectDir;

% Mesh marching options
run_full_dataset = 'none' ; 

% Load/define the surface detection parameters
msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
if exist(msls_detOpts_fn, 'file') 
    load(msls_detOpts_fn)
    disp(['WARNING: overwriting preilastikaxisorder to saved value in ' msls_detOpts_fn])
    disp(['WARNING: overwriting ssfactor to saved value in ' msls_detOpts_fn])
else
    % these are dummy options -- do not save them
    channel = channelsUsed ;
    foreGroundChannel = 1;
    ssfactor = ssfactor;
    niter = 25 ;
    niter0 = 25 ;
    ofn_ply = 'mesh_apical_ms_stab_' ; 
    ofn_ls = 'msls_apical_stab_' ;
    ofn_smoothply = 'mesh_apical_stab_' ;
    
    ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/' ;
    if ~exist(ms_scriptDir, 'dir')    
        ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/' ;
    end
    lambda1 = 1 ;
    lambda2 = 1 ;
    exit_thres = 0.000001 ;
    smoothing = 0.10 ;
    nu = 0.00 ;
    pre_nu = -5 ;
    pre_smoothing = 1 ;
    post_nu = 2;
    post_smoothing = 4 ;
    zdim = 2 ;
    init_ls_fn = 'msls_apical_stab_init.h5';
    mlxprogram = './surface_rm_resample20k_reconstruct_LS3_1p2pc_ssfactor4.mlx' ;
    radius_guess = 40 ;
    center_guess = '100,100,100' ;
    dtype = 'h5' ;
    mask = 'none' ;
    prob_searchstr = '_stab_Probabilities.h5' ;
    ilastikaxisorder= 'cxyz'; ... % axis order as output by ilastik probabilities h5
    imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
    include_boundary_faces = true ;
    pythonVersion = '3' ;
end

% Surface detection parameters --------------------------------------------
detectOptions = struct( 'channel', channel, ...
    'ssfactor', ssfactor, ...
    'niter', niter,...
    'niter0', niter0, ...
    'lambda1', lambda1, ...
    'lambda2', lambda2, ...
    'pressure', nu, ...
    'tension', smoothing, ...
    'post_pressure', post_nu, ...
    'post_tension', post_smoothing, ...
    'exit_thres', exit_thres, ...
    'foreGroundChannel', foreGroundChannel, ...
    'fileName', sprintf( fn, xp.currentTime ), ...
    'mslsDir', mslsDir, ...
    'ofn_ls', ofn_ls, ...
    'ofn_ply', ofn_ply,...
    'ms_scriptDir', ms_scriptDir, ...
    'timepoint', xp.currentTime, ...
    'zdim', zdim, ...
    'pre_pressure', pre_nu, ...
    'pre_tension', pre_smoothing, ...
    'ofn_smoothply', ofn_smoothply, ...
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
    'physicalaxisorder', imsaneaxisorder, ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
    'preilastikaxisorder', preilastikaxisorder, ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
    'ilastikaxisorder', ilastikaxisorder, ... % axis order as output by ilastik probabilities h5. To keep as saved coords use xyzc
    'include_boundary_faces', include_boundary_faces, ...
     'smooth_with_matlab', -1, ...
     'pythonVersion', pythonVersion) ;

% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );

%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
for t = xp.fileMeta.timePoints
    if ~exist(fullfile(projectDir, [sprintf(fn, t) '.h5']), 'file') && ...
            ~exist(fullfile(projectDir, ['../' sprintf(fn, t) '.h5']), 'file')
        disp(['Did not find file: ', fullfile(projectDir, [sprintf(fn, t) '.h5'])])
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();
        % make a copy of the detectOptions and change the fileName
        detectOpts2 = detectOptions ;
        detectOpts2.fileName = sprintf( fn, xp.currentTime ) ;
        xp.setDetectOptions( detectOpts2 );
        xp.detector.prepareIlastik(xp.stack);
        disp(['done outputting downsampled data h5: tp=' num2str(t) ' for surface detection'])
    else
        disp(['h5 ' num2str(t) ' was already output, skipping...'])
    end
end    
disp('Open pre-stabilized h5s with ilastik if not already done')

