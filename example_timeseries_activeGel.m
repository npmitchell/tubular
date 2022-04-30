%% EXAMPLE MASTER PIPELINE FOR TIME SERIES DATA (3D + time)
% by NPMitchell & Dillon Cislo

% This is a pipeline to analyze a dynamic tube-like interface in 3D data.
% A tube-like surface is one that is either cylindrical in topology or
% elongated and spherical in topology. If the initial mesh is a spherical
% topology, then two 'endcaps' will be truncated in order to transform it 
% to a cylindrical topology.
%
% Here we analyze a liquid-phase-separated microtubule-kinesin gel with
% DNA in one phase and microtubules in the other phase.

% TO DO:
% ------
%  - install instructions for morphsnakes as an option
%  - iLastik instructions (optional)
%  - imsane instructions -- tutorial with incorporation 
%  - do we need gptoolbox? so far just for laplacian_smooth() --> run
%  without, for vertex_normals use meshAveragingOperators, barycenter
%  - Make sure Ricci is not default
%  - imsane integralDetector vs morphsnakesDetector
%  - piv 
%  - add clicker for A and P and D
%  - ask Seb what heart tube resolution is

%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
% cd /mnt/data/tubular_test/activeGel/202203221330_region2_ATayarGel/timesequence/ ;
% cd /Volumes/WOON/tubular_test/activeGel/202203221330_region2_ATayarGel/timesequence/ ;
% cd ~/Desktop/tubular_test/activeGel/202203221330_region2_ATayarGel/timesequence/ ;
% cd /mnt/data/tubular_test/activeGel/202203221330_region3_ATayarGel/timesequence/ ;
 cd /mnt/data/tubular_test/activeGel/202203221330_region4_ATayarGel/imagesequence/ ;
% cd /mnt/data/tubular_test/activeGel/202203221330_region5_ATayarGel/imagesequence/ ;
% cd /mnt/data/tubular_test/activeGel/202203221330_region5_ATayarGel/imagesequence/ ;
% cd /mnt/data/tubular_test/activeGel/20220323_passive_region1_ATayarGel/imagesequence/ ;
% cd /mnt/data/tubular_test/activeGel/20220323_passive_region3_ATayarGel/imagesequence/ ;
% cd /mnt/data/tubular_test/activeGel/20220323_passive_region4_ATayarGel/imagesequence/ ;
 cd /mnt/data/tubular_test/activeGel/20220323_passive_region5_ATayarGel/imagesequence/ ;


dataDir = cd ;

%% ADD PATHS TO THIS ENVIRONMENT ==========================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
addpath(fileparts(origpath))
addpath(genpath('utility'))
addpath(genpath('external'))
addpath(genpath('/mnt/data/code/gut_matlab/gptoolbox/'))
addpath(genpath('RicciFlow_MATLAB'))
addpath(genpath('TexturePatch'))
% addpath_recurse('/mnt/data/code/gptoolbox')
% addpath('TexturePatch')
addpath('DEC')
% addpath(fullfile('utility','plotting'))
% go back to the data
cd(dataDir)

%% DEFINE NEW MASTER SETTINGS
if ~exist('./masterSettings.mat', 'file')  
    % Metadata about the experiment
    stackResolution = [1.1 1.1 2.2] ;  % resolution in spaceUnits per pixel
    nChannels = 2 ;             % how many channels is the data (ex 2 for GFP + RFP)
    channelsUsed = [1,2] ;          % which channels are used for analysis
    timePoints = 0:24;       % timepoints to include in the analysis
    ssfactor = 4 ;              % subsampling factor
    flipy = false ;             % whether the data is stored inverted relative to real position in lab frame
    timeInterval = 30 ;          % physical interval between timepoints
    timeUnits = 's' ;         % physical unit of time between timepoints
    spaceUnits = '$\mu$m' ;     % physical unit of time between timepoints
    fn = 'passivereg5_T%02d';        % filename string pattern
    set_preilastikaxisorder = 'xyzc' ; % data axis order for subsampled h5 data (ilastik input)
    swapZT = 0 ;                % whether to swap the z and t dimensions
    masterSettings = struct('stackResolution', stackResolution, ...
        'nChannels', nChannels, ...
        'channelsUsed', channelsUsed, ...
        'timePoints', timePoints, ...
        'ssfactor', ssfactor, ...
        'flipy', flipy, ...
        'timeInterval', timeInterval, ...
        'timeUnits', timeUnits, ...
        'spaceUnits', spaceUnits, ...
        'fn', fn,...
        'swapZT', swapZT, ...
        'set_preilastikaxisorder', set_preilastikaxisorder, ...
        'nU', 100, ...  
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1: Surface detection using TubULAR's getSurface method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup a working directory for the project, where extracted surfaces,
% metadata and debugging output will be stored.  Also specifiy the
% directory containing the data.
cd(dir16bit)
dataDir = dir16bit ;
projectDir = dataDir ;

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

% Set required additional information on the experiment.
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
expMeta.description         = 'microtubule-kinesin gel';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
expMeta.fitTime             = fileMeta.timePoints(first_tp);
expMeta.detectorType        = 'surfaceDetection.integralDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

%% SET DETECTION OPTIONS ==================================================
% Load/define the surface detection parameters
overwriteDetectOptions = true ;
msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
if exist(msls_detOpts_fn, 'file') && ~overwriteDetectOptions
    load(msls_detOpts_fn, 'detectOptions')
else
    outputfilename_ply='mesh_ls_' ;
    outputfilename_ls='ls_' ;
    outputfilename_smoothply = 'mesh_' ;
    ms_scriptDir='/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/' ;   
    init_ls_fn = 'msls_initguess' ;
    meshlabCodeDir = '/mnt/data/code/meshlab_codes/';
    mlxprogram = fullfile(meshlabCodeDir, ...
        'laplace_surface_rm_resample30k_reconstruct_LS3_1p2pc_ssfactor4.mlx') ;
    prob_searchstr = '_stab_Probabilities.h5' ;
    preilastikaxisorder = set_preilastikaxisorder; ... % axis order in input to ilastik as h5s. To keep as saved coords use xyzc
    ilastikaxisorder= 'cxyz'; ... % axis order as output by ilastik probabilities h5
    imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc
    
    % Name the output mesh directory --------------------------------------
    mslsDir = [fullfile(projectDir, 'msls_output') filesep];

    % Surface detection parameters ----------------------------------------
    detectOptions = struct( 'channel', 1, ...
        'ssfactor', 4, ...
        'niter', 40,...
        'niter0', 100, ...
        'enforceSingleComponent', true, ...
        'pre_pressure', 0, ...
        'pressure', -0.01, ...
        'tension', 0.5, ...
        'post_pressure', 0, ...
        'exit_thres', 1e-7, ...
        'foreGroundChannel', 1, ...
        'fileName', fn, ...
        'mslsDir', mslsDir, ...
        'ofn_ls', outputfilename_ls, ...
        'ofn_ply', outputfilename_ply,...
        'ms_scriptDir', ms_scriptDir, ...
        'timepoint', timePoints(1), ...
        'zdim', 2, ...
        'ofn_smoothply', outputfilename_smoothply, ...
        'mlxprogram', mlxprogram, ...
        'init_ls_fn', init_ls_fn, ... % set to none to load prev tp
        'run_full_dataset', projectDir,... % projectDir, ... % set to 'none' for single tp
        'radius_guess', 40, ...
        'dset_name', 'exported_data',...
        'save', true, ... % whether to save images of debugging output
        'center_guess', '30,30,10',... % xyz of the initial guess sphere ;
        'plot_mesh3d', false, ...
        'dtype', 'mat',...
        'mask', 'none',...
        'mesh_from_pointcloud', false, ...
        'prob_searchstr', prob_searchstr, ...
        'physicalaxisorder', imsaneaxisorder, ... 
        'preilastikaxisorder', preilastikaxisorder, ... 
        'ilastikaxisorder', ilastikaxisorder, ... 
        'include_boundary_faces', true, ...
        'smooth_with_matlab', 0.075);  % set this to >0 to use matlab laplacian filter instead of meshlab

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
mslsDir = detectOptions.mslsDir ;


%% Collect options into xp experiment struct
xp = struct('expMeta', expMeta, 'fileMeta', fileMeta, ...
    'detectOptions', detectOptions) ;


%% Now we have 3d data volumes and surfaces. Define a TubULAR object. 
% To visualize data on these surfaces and compute how these surfaces deform
% we now define TubULAR object.
nU = masterSettings.nU ;
nV = masterSettings.nV ;
opts.meshDir = mslsDir ;        % Directory where meshes reside
opts.flipy = flipy ;            % Set to true if data volume axes are inverted in chirality wrt physical lab coordinates
opts.timeInterval = timeInterval ; % Spacing between adjacent timepoints in units of timeUnits 
opts.timeUnits = timeUnits ;    % units of time, so that adjacent timepoints are timeUnits * timeInterval apart
opts.spaceUnits = spaceUnits ;  % Units of space in LaTeX, for ex '$mu$m' for micron
opts.nU = nU ;                  % How many points along the longitudinal axis to sample surface
opts.nV = nV ;                  % How many points along the circumferential axis to sample surface
opts.t0 = xp.fileMeta.timePoints(1) ;   % reference timepoint used to define surface-Lagrangian and Lagrangian measurements
opts.normalShift = 0 ;         % Additional dilation acting on surface for texture mapping
opts.a_fixed = 2.0 ;            % Fixed aspect ratio of pullback images. Setting to 1.0 is most conformal mapping option.
opts.adjustlow = [304, 221] ;         % floor for intensity adjustment
opts.adjusthigh = [22000, 2500] ; %3300       % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;   % Method for following surface in surface-Lagrangian mapping [(s,phi) coordinates]
opts.lambda_mesh = 0.00 ;       % Smoothing applied to the mesh before DEC measurements
opts.lambda = 0.0 ;             % Smoothing applied to computed values on the surface
opts.lambda_err = 0.0 ;         % Additional smoothing parameter, optional
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')



%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
tubi.prepareIlastik()

%% TRAIN NON-STABILIZED DATA IN ILASTIK TO IDENTIFY SURFACE ==============
% Open ilastik, train on h5s until probabilities and uncertainty are 
% satisfactory for extracting a mesh. For example, here we train on the
% membrane (channel 1) and the yolk (channel 2), so that a level set will
% enclose the yolk but not escape through the membrane.

%% Create MorphSnakes LevelSets from the Probabilities output of ilastik ==
% Now detect all surfaces
disp('Surfaces are extracted serially using active contours')
overwrite = false  ;
tubi.getMeshes(overwrite);

%% Inspect all meshes in 3D
for tp = xp.fileMeta.timePoints
    % Load the mesh
    meshfn = sprintf( tubi.fullFileBase.mesh, tp ) ;    
    mesh = read_ply_mod(meshfn) ;
    % Plot the mesh in 3d. Color here by Y coordinate
    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        mesh.v(:, 3), 'edgecolor', 'none', 'Facealpha', 0.5)
    title(['t=' num2str(tp)])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    view(2)
    pause(0.5)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2: TubULAR -- surface parameterization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Define global orientation frame (for viewing in canonical frame)
% Compute APDV coordinate system
alignAPDVOpts = struct() ;
alignAPDVOpts.overwrite = false ;
tubi.computeAPDVCoords(alignAPDVOpts) ;

%% PLOT ALL TEXTURED MESHES IN 3D (OPTIONAL: this is SLOW) ================
% Establish texture patch options
metadat = struct() ;
metadat.reorient_faces = false ;            % set to true if some mesh normals may be inverted (requires gptoolbox if true)
metadat.normal_shift = tubi.normalShift ;   % normal push, in pixels, along normals defined in data XYZ space
metadat.texture_axis_order = [1 2 3] ;      % texture space sampling. If the surface and dataspace have axis permutation, enter that here
metadat.blackFigure = true ;
metadat.channel = 2;
metadat.perspective_angle = [-20,-20] ;
metadat.perspective_angle = [70,-50] ;
% metadat.perspective_angle = [130,-50] ;
Options.PSize = 5 ;          % Psize is the linear dimension of the grid drawn on each triangular face. Set PSize > 1 for refinement of texture on each triangle of the surface triangulation. Higher numbers are slower but give more detailed images.
Options.numLayers = [10, 2];  % how many layers to MIP over/bundle into stack, as [outward, inward]
Options.layerSpacing = 1 ;   % Distance between layers over which we take MIP, in pixels, 
% metadat.plot_time_points = [8] ;

metadat.plot_dorsal = false ;
metadat.plot_left = false ;
metadat.plot_right = false ;
metadat.plot_ventral = false ;
metadat.plot_perspective = true ;
metadat.overwrite = true ;
 % Options.Imax = 2000.0 ;
% tubi.data.adjusthigh = [2500, 400] ;
% tubi.data.adjusthigh = [2500, 400] ;

% Plot on surface for all timepoints 
tubi.plotSeriesOnSurfaceTexturePatch(metadat, Options)

%% 
for tp = tubi.xp.fileMeta.timePoints
    tubi.setTime(tp)
    tubi.maskCurrentDataWithMesh()
end

%% 
error('done with pipeline')

%% Select the endcaps for the centerline computation (A and P) and a point
% along which we will form a branch cut for mapping to the plane (D).
apdvOpts = struct() ;
apdvOpts.overwrite = false ;
[apts_sm, ppts_sm] = tubi.computeAPDpoints(apdvOpts) ;

%% Align the meshes in the APDV global frame & plot them
tubi.alignMeshesAPDV(alignAPDVOpts) ;

disp('done')

%% EXTRACT CENTERLINES
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
idOptions.ssr_thres = 15 ;  % distance of sum squared residuals in um as threshold for removing spurious centerlines
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

methodOpts.overwrite = true ;
methodOpts.save_figs = true ;   % save images of cutMeshes along the way
methodOpts.preview = false  ;     % display intermediate results
tubi.sliceMeshEndcaps(endcapOpts, methodOpts) ;

%% Clean Cylinder Meshes
% This removes "ears" from the endcaps of the tubular meshes (cylindrical
% meshes)
cleanCylOptions = struct() ;
cleanCylOptions.overwrite = false ;
tubi.cleanCylMeshes(cleanCylOptions)
disp('done cleaning cylinder meshes')
    

