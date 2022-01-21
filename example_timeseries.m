%% EXAMPLE MASTER PIPELINE FOR TIME SERIES DATA (3D + time)
% by NPMitchell & Dillon Cislo

% TO DO:
% ------
%  - redefine APDV coords as viewingCoordinates (global ref frame)
%  - install instructions for morphsnakes as an option
%  - iLastik instructions (optional)
%  - imsane instructions (optional)
%  - do we need gptoolbox?
%  - DEC tutorials --> STANDALONE AND INTEGRATED TUTORIALS


%
% This is a pipeline to take the surface of the growing Drosophila gut and
% conformally map patches of it to the unit disk

% ============= %
%  CONVENTIONS
% ============= % 
% >> APDV ilastik training <<
% ---------------------------
% Train on anterior (A), posterior (P), background (B), and 
% dorsal anterior (D) location in different iLastik channels.
% Default order is 1-A, 2-P, 3-bg, 4-D. 
% Train on pre-stabilized H5s, NOT on stabilized H5s. 
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% By default, name the h5 file output from iLastik as 
%    QS.fullFileBase.apdProb --> [QS.fileBase.name '_Probabilities_APDVcoords.h5'] 
% and/or 
%    QS.fullFileBase.apCenterlineProb --> [QS.fileBase.name '_Probabilities_apcenterline.h5'] 
% 
% The APDV coordinate system -- this is where global rotation is defined -- 
% is computed via 
%    apdvOpts.aProbFileName = QS.fullFileBase.apdProb ;
%    apdvOpts.pProbFileName = QS.fullFileBase.apCenterlineProb ;
%    QS.computeAPDVCoords(alignAPDVOpts) ;
% However, the APD (anterior, posterior, dorsal) centers of mass (COMs) for
% centerline extraction can be independent. 
% Forexample, if the original segmentation used a third channel which 
% happen to be the posterior end, we can specify the posterior COM (P) via:
%    apdvOpts.aProbFileName = QS.fullFileBase.apdProb ; % filename pattern for the apdv training probabilities
%    apdvOpts.pProbFileName = QS.fullFileBase.prob ;
%    apdvOpts.posteriorChannel = 3 ;
%    [acom_sm, pcom_sm] = QS.computeAPDCOMs(apdvOpts) ;
%
% >> Centerline extraction <<
% ---------------------------
% STEP 1: THE INITIAL CENTERLINE IS COMPUTED
% Centerlines are first computed using a fast marching procedure on a 
% distance transform (DT) of the mesh data in dataspace. The startpoint and
% endpoint of this initial centerline can be specified using APDVcoords
% training (and the resulting COM extraction) from the previous section or
% else can be unique to the centerline extraction. This independence is
% useful, for example, in convoluted organs like the gut, where the 
% centerline of the muscle data curves away from the AP axis in the
% posterior where the midgut turns back on itself to connect to the
% hindgut.
%
% For redundant use of APDV coords for centerline, use the option:
%
% For independent training use the following convention:
% Train on anterior (A), posterior (P), background (B), and 
% dorsal anterior (D) location in different iLastik channels.
% Default order is 1-A, 2-P, 3-bg, 4-D. 
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
%
% STEP 2: PRECISE CENTERLINE
% The initial centerline is used to constrain the winding number of the
% orbifold mapping. The second centerline is a series of average 3d points,
% one for each DV ring of the quasi-axisymmetric pullback. 
%

%% Clear workspace ========================================================
% We start by clearing the memory and closing all figures
clear; close all; clc;
cd /mnt/data/code/tubular/example/ ;

dataDir = cd ;

%% PATHS ==================================================================
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
aux_paths_and_colors
% go back to the data
cd(dataDir)

%% DEFINE NEW MASTER SETTINGS
if overwrite_masterSettings || ~exist('./masterSettings.mat', 'file')
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
    t0_for_phi0 = 1 ;
    tidx0_for_stab = 10 ;
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
        't0_for_phi0', t0_for_phi0, ... % 40 for mef2 single channel, 110 for CAAX excellent
        'tidx0_for_stab', tidx0_for_stab, ...
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
    scale = masterSettings.scale ;      % scale for conversion to 16 bit
    file32Base = masterSettings.file32Base ; 
    fn = masterSettings.fn ;
    fn_prestab = masterSettings.fn_prestab ;
    set_preilastikaxisorder = masterSettings.set_preilastikaxisorder ;
    swapZT = masterSettings.swapZT ;
    t0_for_phi0 = masterSettings.t0_for_phi0 ;
    tidx0_for_stab = masterSettings.tidx0_for_stab ;
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;
end
dir32bit = fullfile(dataDir, 'deconvolved_32bit') ;
dir16bit = fullfile(dataDir, 'deconvolved_16bit') ;
dir16bit_prestab = fullfile(dir16bit, 'data_pre_stabilization') ;

%% END OF EXPERIMENT METADATA =============================================
% =========================================================================
% =========================================================================
%% -IIV. make MIPs for 32bit images
% Skip if already done
mipDir = fullfile(dir32bit, 'mips32bit') ;
Options.overwrite_mips = overwrite_mips ;
Options.scale = scale ;
Options.allowMissingFiles = false ;
makeMips(timePoints, dir32bit, file32Base, mipDir, Options)

%  -IV. convert 32 to 16bit images
% Skip if already done
convert32to16bit(timePoints, scale, dir32bit, dir16bit_prestab,...
    file32Base, fn_prestab)

% % Rename stab to prestab -- hack
% % fns = fullfile('./deconvolved_16bit/Time*stab')
% % for qq=1:length(fns)
% %     command = ['mv ' fullfile(fns.folder, fns.name) fullfile(dir16bit, fns.name
% % end

% -III. make MIPs for 16bit images
% Skip if already done
mipDir = fullfile(dir16bit_prestab, 'mips') ;
Options.overwrite_mips = false ;
Options.scale = -1 ; % do NOT rescale intensities during intensity projection
makeMips(timePoints, dir16bit_prestab, fn_prestab, mipDir, Options)

%  -II. stabilize images, based on script stabilizeImagesCorrect.m
% Skip if already done
% name of directory to check the stabilization of mips
mips_stab_check = fullfile(mipDir, 'stab_check') ;
mipoutdir = fullfile(mipDir, 'mips_stab') ;
% Choose bit depth as typename
typename = 'uint16' ;
% Give file names for I/O
fileNameIn = fullfile(dir16bit_prestab, fn_prestab) ;
fileNameOut = fullfile(dir16bit, [fn '.tif']) ;
rgbName = [fn '.png'] ;
typename = 'uint16' ; 
Options = struct(); 
Options.overwrite_mips = false ;
Options.overwrite_tiffs = false ;
Options.im_intensity = 1 ; % 0.01 ;
Options.imref_intensity = 1 ; % 0.005 ;
Options.forceNoDx = false ;
stabilizeImages(fileNameIn, fileNameOut, rgbName, typename, ...
    timePoints, timePoints, tidx0_for_stab, ...
    mipDir, mipoutdir, mips_stab_check, Options)


%   -I. master_gut_timeseries_prestab_for_training.m
% Skip if already done
cd(dir16bit)
dataDir = cd ;
masterSettings.dir16bit_prestab = dir16bit_prestab ;
makeH5SeriesPrestabForTraining(masterSettings)
cd(dir16bit)

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

%% SET DETECT OPTIONS =====================================================
% Must run this section for later functionality.
% Mesh extraction options
if run_full_dataset_ms
    run_full_dataset = projectDir ; 
else
    run_full_dataset = 'none' ;
end

% Load/define the surface detection parameters
msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat') ;
if exist(msls_detOpts_fn, 'file') && ~overwrite_detOpts
    load(msls_detOpts_fn, 'detectOptions')
else
    channel = 1;
    foreGroundChannel = 1;
    % ssfactor = 4;
    % niter = 15 ;
    % niter0 = 15 ;
    % ofn_ply = 'mesh_' ; % mesh_apical_ms_stab_' ; 
    % ofn_ls = 'msls_' ; % mesh_apical_stab_' ;
    % ofn_smoothply = 'mesh_apical_stab_' ; % mesh_apical_stab_' ;
    % lambda1 = 1 ;
    % lambda2 = 1 ;
    % exit_thres = 0.0001 ;
    % if contains(projectDir, 'LifeActRuby') || contains(projectDir, 'CAAXmCherry') 
    %     nu = 4 ;  % For volumetric
    %     smoothing = 0.20 ;
    % else
    %     nu = 0.00 ;
    %     smoothing = 0.10 ;
    % end
    % pre_nu = -4 ;
    % pre_smoothing = 0 ;
    % post_nu = 3 ;
    % post_smoothing = 1 ;
    zdim = 2 ;
    
    % for caax_great
    ssfactor=4
    niter=35
    niter0=35
    ofn_ply='mesh_apical_ms_stab_'
    ofn_ls='msls_'
    ofn_smoothply = 'mesh_stab_' ; % mesh_apical_stab_' ;
    ms_scriptDir='/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper/'
    pre_nu=-5
    pre_smoothing=0
    lambda1=1
    lambda2=1
    exit_thres=0.000005
    smoothing=0.5
    nu=0
    post_nu=2
    post_smoothing=3
    
    init_ls_fn = 'msls_initguess.h5';
    % mlxprogram = fullfile(meshlabCodeDir, ...
    %     'laplace_refine_HCLaplace_LaplaceSPreserve_QuadEdgeCollapse60kfaces.mlx') ;
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
    
    % Name the output mesh directory ------------------------------------------
    % msls_exten = ['_prnu' strrep(strrep(num2str(pre_nu, '%d'), '.', 'p'), '-', 'n')];
    % msls_exten = [msls_exten '_prs' strrep(num2str(pre_smoothing, '%d'), '.', 'p') ];
    % msls_exten = [msls_exten '_nu' strrep(num2str(nu, '%0.2f'), '.', 'p') ];
    % msls_exten = [msls_exten '_s' strrep(num2str(smoothing, '%0.2f'), '.', 'p') ];
    % msls_exten = [msls_exten '_pn' num2str(post_nu, '%d') '_ps',...
    %     num2str(post_smoothing)];
    % msls_exten = [msls_exten '_l' num2str(lambda1) '_l' num2str(lambda2) ];
    meshDir = [projectDir 'msls_output' filesep];
    % meshDir = [meshDir msls_exten '/'] ;

    % Surface detection parameters --------------------------------------------
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
        'mslsDir', meshDir, ...
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
detectOptions.run_full_dataset = run_full_dataset ;
detectOptions.ms_scriptDir = ms_scriptDir ;
meshDir = detectOptions.mslsDir ;
% These are now in QS.
% meshFileBase = [ofn_smoothply '%06d'] ;
% alignedMeshBase = [ofn_smoothply '%06d_APDV_um'] ;

% Set detect options ------------------------------------------------------
xp.setDetectOptions( detectOptions );
disp('done')

%% CREATE THE SUBSAMPLED H5 FILE FOR INPUT TO ILASTIK =====================
% skip if already done

for tt = xp.fileMeta.timePoints
    if ~exist(fullfile(projectDir, 'stabilized_h5s', [sprintf(fn, tt) '.h5']), 'file')
        disp(['Did not find file: ', fullfile(projectDir, 'stabilized_h5s', [sprintf(fn, tt) '.h5'])])
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


%if ~exist([projectDir sprintf(fn, xp.currentTime) '.h5'], 'file')
 %   xp.detector.prepareIlastik(xp.stack);
  %  disp('done outputting downsampled data h5 for surface detection')
%else
 %   disp('h5 was already output, skipping...')
%end
%disp('Open with ilastik if not already done')

%% TRAIN NON-STABILIZED DATA IN ILASTIK TO IDENTIFY APICAL/YOLK ===========
% Skip if already done.
% Open ilastik, train pre-stab h5s until probabilities and uncertainty are 
% satisfactory, then run on stab images.

%% Note that old h5's used to be different order. To convert, do
% Skip 
% if false
%     tmp = h5read(fullfile(meshDir, init_ls_fn), '/implicit_levelset');
%     % OLD ORDER: yxz for implicit levelset
%     tmp2 = permute(tmp, [3,2,1]);
%     h5create(fullfile(meshDir, init_ls_fn), '/implicit_levelset', size(tmp2), 'Datatype', 'int8')
%     h5write(fullfile(meshDir, init_ls_fn), '/implicit_levelset', tmp2)
% end

%% Create MorphSnakesLevelSet from the Probabilities from ilastik ========
% Skip if already done
% Now detect all surfaces
detectOptions.run_full_dataset = 'none' ;  % projectDir ; % 'none' ;  % override here
if strcmp(detectOptions.run_full_dataset, projectDir)
    % assert(run_full_dataset_ms)
    disp('Running dataset mode')
    xp.setTime(xp.fileMeta.timePoints(1));
    detectOpts2 = detectOptions ;
    detectOpts2.fileName = sprintf( fn, xp.currentTime ) ;
    detectOpts2.nu = 4 ;
    detectOpts2.niter0 = 5 ;
    xp.setDetectOptions( detectOpts2 );
    xp.detectSurface();
else
    assert(~run_full_dataset_ms)
    assert(strcmp(detectOptions.run_full_dataset, 'none'))
    % Morphosnakes for all remaining timepoints INDIVIDUALLY ==============
    for tp = xp.fileMeta.timePoints(59:59)
        %try
            xp.setTime(tp);
            % xp.loadTime(tp) ;
            % xp.rescaleStackToUnitAspect();

            % make a copy of the detectOptions and change the fileName
            detectOpts2 = detectOptions ;
            % detectOpts2.post_smoothing = 1 ;
            detectOpts2.timepoint = xp.currentTime ;
            detectOpts2.fileName = sprintf( fn, xp.currentTime );
            detectOpts2.mlxprogram = fullfile(meshlabCodeDir, ...
                  'surface_rm_resample30k_reconstruct_LS3_ssfactor4.mlx') ;
            %  _octree12.mlx') ;
            % detectOpts2.mlxprogram = fullfile(meshlabCodeDir, ...
            %     'laplace_surface_rm_resample25k_reconstruct_LS3_1p2pc_ssfactor4_vip10.mlx') ;
            % detectOpts2.mlxprogram = fullfile(meshlabCodeDir, ...
            %      'laplace_surface_rm_resample25k_reconstruct_LS3_wu13_ssfactor4.mlx') ;
            xp.setDetectOptions( detectOpts2 );
            xp.detectSurface();
            % For next time, use the output mesh as an initial mesh
            detectOpts2.init_ls_fn = 'none' ;
        % catch
        %        error('here')
        %     disp('Could not create mesh -- skipping for now')
        %     % On next timepoint, use the tp previous to current time
        %     detectOptions.init_ls_fn = [detectOptions.ofn_ls, ...
        %             num2str(tp - 1, '%06d' ) '.' detectOptions.dtype] ;
        % end
    end
end

%% Define QuapSlap object
nU = masterSettings.nU ;
nV = masterSettings.nV ;
opts.meshDir = meshDir ;
opts.flipy = flipy ;
opts.timeInterval = timeInterval ;
opts.timeUnits = timeUnits ;
opts.spaceUnits = spaceUnits ;
opts.nU = nU ;
opts.nV = nV ;
opts.normalShift = 10 ;
opts.a_fixed = 2.0 ;
opts.adjustlow = 1.00 ;         % floor for intensity adjustment
opts.adjusthigh = 99.9 ;        % ceil for intensity adjustment (clip)
opts.phiMethod = 'curves3d' ;
opts.lambda_mesh = 0.002 ;
opts.lambda = 0.01 ;
opts.lambda_err = 0.01 ;
disp('defining QS')
QS = QuapSlap(xp, opts) ;
disp('done')

%% Make some mips of shallow stacks
adjustIV = false ; 
% QS.makeMIPs(1, {400:490, 510:550, 570:630, 800:880}, [], adjustIV)
stacks = {400:490, 510:550, 570:630, 800:880} ;
QS.makeMIPs(1, stacks, [], adjustIV)

%% HACK -- transpose meshes
% xyzc -> zyxc
% for tp = 50:60
%     disp(['transposing axes for mesh t=' num2str(tp)])
%     meshfn = sprintf(QS.fullFileBase.mesh, tp) ;  
%     msh = read_ply_mod(meshfn) ;
%     msh.v = msh.v(:, [2, 3,1]) ;
%     msh.vn = msh.vn(:, [2,3,1]) ;
%     plywrite_with_normals(meshfn, msh.f, msh.v, msh.vn)
% end
% disp('done with hack to transpose meshes')



%% Inspect all meshes in 3D
% Skip if already done

% Make an output directory for the quick-and-dirty inspection
for tp = xp.fileMeta.timePoints(48:end)
    % Load the mesh
    meshfn = sprintf( QS.fullFileBase.mesh, tp ) ;    
    mesh = read_ply_mod(meshfn) ;
    % Plot the mesh in 3d. Color here by Y coordinate
    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        mesh.v(:, 3), 'edgecolor', 'none', 'Facealpha', 0.5)
    % saveas(gcf, fullfile(outputdir, sprintf('inspect_%04d.png', tp)))
    title(['t=' num2str(tp)])
    axis equal
    view(2)
    pause(0.1)
end



%% APDV ilastik training
% Train on anterior (A), posterior (P), background (B), and 
% dorsal anterior (D) location in different iLastik channels.
% Perform on pre-stabilized H5s, NOT on stabilized H5s. 
% anteriorChannel, posteriorChannel, and dorsalChannel specify the iLastik
% training channel that is used for each specification.
% Name the h5 file output from iLastik as ..._Probabilities_apcenterline.h5
% Train for anterior dorsal (D) only at the first time point, because
% that's the only one that's used. 
% Dorsal anterior for the gut is at the fused site where additional 
% 48YGAL4-expressing muscle-like cells form a seam.
% Posterior is at the rear of the yolk, where the endoderm closes, for 
% apical surface training.
% Anterior is at the junction of the midgut with the foregut.

%% 3. align_meshes_APDV or load transformations if already done
% Skip if already done.
overwrite_alignAPDVOpts = false ;

% Try to load the results to test if we must do calculation 
try
    [rot, trans] = QS.getRotTrans() ;
    [xyzlim_raw, xyzlim, xyzlim_um, xyzlim_buff] = getXYZLims(QS) ;
    assert(exist(QS.fileName.startendPt, 'file') ~= 0)
    redo_alignmesh = false ;
    
    % check all timepoints
    tidx = 1; 
    while ~redo_alignmesh && tidx <= length(xp.fileMeta.timePoints)
        tt = xp.fileMeta.timePoints(tidx) ;
        searchfn = sprintf(QS.fullFileBase.alignedMesh, tt);
        if ~strcmp(searchfn(end-3:end), '.ply')
            searchfn = [searchfn '.ply'] ;
        end
        if ~exist(searchfn, 'file')
            redo_alignmesh = true ;
            disp(['did not find aligned mesh for tp = ' num2str(tt)])
        end
        tidx = tidx + 1;
    end
    
catch
    redo_alignmesh = true ;
end

if redo_alignmesh || overwrite_APDVMeshAlignment || overwrite_APDVCOMs
    disp('Calculating/Overwriting APDVMeshAlignment')
    % OUTPUTS OF THIS SECTION
    % -------
    % xyzlim.txt 
    %   xyzlimits of raw meshes in units of full resolution pixels (ie not
    %   downsampled)
    % xyzlim_APDV.txt 
    %   xyzlimits of rotated and translated meshes in units of full resolution 
    %   pixels (ie not downsampled)
    % xyzlim_APDV_um.txt 
    %   xyz limits of rotated and translated meshes (APDV coord sys) in microns
    % rotation_APDV.txt
    %   rotation matrix to align mesh to APDV frame
    % translation_APDV.txt
    %   translation vector to align mesh to APDV frame
    % xyzlim.txt 
    %   raw bounding box in original frame (not rotated), in full res pixels
    % xyzlim_APDV.txt
    %   bounding box in rotated frame, in full resolution pixels
    % xyzlim_APDV_um.txt
    %   bounding box in rotated frame, in microns
    % apdv_coms_rs.h5
    %   Centers of mass for A,aux_paths_and_colors P, and D in microns in rotated APDV coord system
    % 
    % Notes
    % -----
    % vertices are in units of pixels (at full resolution)
    % To take mesh to rotated + translated mesh in physical units, apply:
    %         xs = mesh.vertex.z ;
    %         ys = mesh.vertex.y ;
    %         zs = mesh.vertex.x ;
    %         vtx_rs = (rot * vtx' + trans)' * resolution
    %         
    save_figs = true ;  % save images of cntrline, etc, along the way
    preview = false ;  % display intermediate results, for debugging
    dorsal_thres = 0.5 ;  % threshold for extracting Dorsal probability cloud 
    buffer = 5 ;  % extra space in meshgrid of centerline extraction, to ensure mesh contained in volume
    plot_buffer = 40;  % extra space in plots, in um
    weight = 0.1;  % for speedup of centerline extraction. Larger is less precise
    normal_step = 0.5 ;  % how far to move normally from ptmatched vtx if a/pcom is not inside mesh
    eps = 0.01 ;  % value for DT outside of mesh in centerline extraction
    meshorder = 'xyz' ;  % ordering of axes in loaded mesh wrt iLastik output
    anteriorChannel = 1;  % which channel of APD training is anterior
    posteriorChannel = 2;  % which channel of APD training is posterior 
    dorsalChannel = 4 ;  % which channel of APD training is dorsal

    clearvars opts
    optsfn = fullfile(projectDir, 'alignAPDV_Opts.mat') ;
    if exist(optsfn, 'file') && ~overwrite_alignAPDVOpts
        disp('Loading options from disk')
        load(optsfn, 'alignAPDVOpts', 'apdvOpts')
    else
        disp('No alignAPDV_Opts on disk or overwriting, defining')
        apdvOpts.smwindow = 30 ;
        apdvOpts.dorsal_thres = dorsal_thres ;
        apdvOpts.buffer = buffer ;  
        apdvOpts.plot_buffer = plot_buffer ;
        apdvOpts.anteriorChannel = anteriorChannel ;
        apdvOpts.posteriorChannel = posteriorChannel ;
        apdvOpts.dorsalChannel = dorsalChannel ;% filename pattern for the apdv training probabilities
        apdvOpts.tpref = t0_for_phi0 ;
        apdvOpts.ilastikOutputAxisOrder = 'cyxz' ;
        
        alignAPDVOpts.weight = weight ;
        alignAPDVOpts.normal_step = normal_step ;
        alignAPDVOpts.eps = eps ;
        alignAPDVOpts.meshorder = meshorder ;
        alignAPDVOpts.anteriorChannel = anteriorChannel ;
        alignAPDVOpts.posteriorChannel = posteriorChannel ;
        alignAPDVOpts.dorsalChannel = dorsalChannel ;
        alignAPDVOpts.dorsal_thres = 0.5 ;
        alignAPDVOpts.tref = t0_for_phi0 ;
        alignAPDVOpts.ilastikOutputAxisOrder = 'cyxz' ;
        
        % alignAPDVOpts.fn = fn ;  % filename base

        % save the optionsc
        save(optsfn, 'alignAPDVOpts', 'apdvOpts')
    end

    % Add script-instance-specific options
    apdvOpts.overwrite = overwrite_APDVCOMs ;  % recompute APDV coms from training
    apdvOpts.check_slices = false ;
    apdvOpts.preview = preview ;
    apdvOpts.preview_com = false ;
    
    %% Compute APDV coordinate system
    alignAPDVOpts.overwrite = false ;
    QS.computeAPDVCoords(alignAPDVOpts) ;
    
    % Compute the APD COMs
    apdvOpts.overwrite = overwrite_APDVCOMs ;
    [acom_sm, pcom_sm] = QS.computeAPDCOMs(apdvOpts) ;
    
    % Align the meshes APDV & plot them
    alignAPDVOpts.overwrite_ims = overwrite_alignedMeshIms ;  % overwrite images even if centerlines are not overwritten
    alignAPDVOpts.overwrite = overwrite_APDVCOMs || overwrite_APDVMeshAlignment ; % recompute APDV rotation, translation
    QS.alignMeshesAPDV(alignAPDVOpts) ;
else
    % Display APDV COMS over time
    acom_sm = h5read(QS.fileName.apdv, '/acom_sm') ;
    pcom_sm = h5read(QS.fileName.apdv, '/pcom_sm') ;
    acoms = h5read(QS.fileName.apdv, '/acom') ;
    pcoms = h5read(QS.fileName.apdv, '/pcom') ;
    dcom = dlmread(QS.fileName.dcom) ;
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        % Plot the APDV points
        clf
        plot3(acom_sm(tidx, 1), acom_sm(tidx, 2), acom_sm(tidx, 3), 'ro')
        hold on;
        plot3(acoms(tidx, 1), acoms(tidx, 2), acoms(tidx, 3), 'r.')
        plot3(pcom_sm(tidx, 1), pcom_sm(tidx, 2), pcom_sm(tidx, 3), 'b^')
        plot3(pcoms(tidx, 1), pcoms(tidx, 2), pcoms(tidx, 3), 'b.')
        plot3(dcom(1, 1), dcom(1, 2), dcom(1, 3), 'cs')
        axis equal
        title(['t = ', num2str(tp)]) 
        pause(0.01)
    end
    disp('Already done')
end
disp('done')
clearvars normal_step 


%% MAKE MASKED DATA FOR PRETTY VIDEO ======================================
% Skip if already done
% Generate masks for isolating intensity data of the gut alone, for making
% pretty videos without fiducial markers and other spurious fluorescent
% bits. 
% Save output h5s trained on stabilized h5s from iLastik as 
%   -->   <QS.fileBase.name>_Probabilities_mask3d.h5
% and, optionally, a probability cloud to exclude 
% (applied post extraction of previous mask)
%   -->   <QS.fileBase.name>_Probabilities_maskDorsal.h5
% options = struct() ;
% options.maskMask = false ;
% QS.generateMaskedData(options)

%% MAKE ORIENTED MASKED DATA FOR PRETTY VIDEO =============================
% % Skip if already done
% QS.alignMaskedDataAPDV()


%% MUSCLE: PLOT ALL TEXTURED MESHES IN 3D =========================================
% Skip if already done
overwrite_TextureMeshOpts = false ;

% Get limits and create output dir
% Establish texture patch options
metafn = fullfile(QS.dir.texturePatchIm, 'metadat_muscle.mat') ;
if ~exist(metafn, 'file') || overwrite_TextureMeshOpts
    [~,~,~,xyzbuff] = QS.getXYZLims() ;
    xyzbuff(:, 1) = xyzbuff(:, 1) - 20 ; 
    xyzbuff(:, 2) = xyzbuff(:, 2) + 20 ; 
    % Define & Save metadata
    metadat.xyzlim = xyzbuff ;                  % xyzlimits
    metadat.reorient_faces = false ;            % if some normals are inverted
    metadat.normal_shift = QS.normalShift ;             % normal push, in pixels, along normals defined in data XYZ space
    metadat.texture_axis_order = QS.data.axisOrder ;    % texture space sampling
    metadat.smoothing_lambda = 0.006 ;   
    metadat.texture_shift = 0 ;
    
    % Psize is the linear dimension of the grid drawn on each triangular face
    Options = struct() ;
    Options.PSize = 5 ;
    Options.EdgeColor = 'none';
    QS.getRotTrans() ;
    Options.Rotation = QS.APDV.rot ;
    Options.Translation = QS.APDV.trans ;
    Options.Dilation = QS.APDV.resolution ;
    % outward, inward --> reasonable values are ~[3,3]
    
    % ENDODERM
    % Options.numLayers = [1, -5];  % at layerSpacing 2, 2 marches ~0.5 um 
    % Options.layerSpacing = 2 ;
    
    % MUSCLE
    metadat.normal_shift = 35 ;  % normal push, in pixels, along normals defined in data XYZ space
    Options.numLayers = [13, 0];  % at layerSpacing 2, 2 marches ~0.5 um 
    Options.layerSpacing = 2 ;
    
    
    % Save it
    save(metafn, 'metadat', 'Options')
else
    disp('Loading options for muscle texturepatch')
    load(metafn, 'metadat', 'Options')
end


% MUSCLE ADJUST METADAT OPTIONS
metadat.normal_shift = 35 ;  % normal push, in pixels, along normals defined in data XYZ space
Options.numLayers = [13, 0];  % at layerSpacing 2, 2 marches ~0.5 um 
Options.layerSpacing = 2 ;
metadat.smoothing_lambda = 0.006 ;   

% Save it
save(metafn, 'metadat', 'Options')
    
% Use first timepoint's intensity limits throughout
QS.setDataLimits(QS.xp.fileMeta.timePoints(1), 1.0, 99.99)

% Plot on surface for all TP 
options = metadat ;
options.overwrite = false ;
options.plot_dorsal = true ;
options.plot_ventral = true ;
options.plot_right = true ;
options.plot_left = true ;
options.timePoints = QS.xp.fileMeta.timePoints ;
options.figOutDir = fullfile(QS.dir.texturePatchIm, 'muscle_layer') ;

options.plot_perspective = true ;
options.texture_axis_order = [2 1 3] ; % for Hand Hand HistGFP
QS.plotSeriesOnSurfaceTexturePatch(options, Options)
clearvars Options



%% ENDODERM: PLOT ALL TEXTURED MESHES IN 3D =========================================
% Skip if already done
overwrite_TextureMeshOpts = false ;

% Get limits and create output dir
% Establish texture patch options
metafn = fullfile(QS.dir.texturePatchIm, 'metadat.mat') ;
if ~exist(metafn, 'file') || overwrite_TextureMeshOpts
    [~,~,~,xyzbuff] = QS.getXYZLims() ;
    xyzbuff(:, 1) = xyzbuff(:, 1) - 20 ; 
    xyzbuff(:, 2) = xyzbuff(:, 2) + 20 ; 
    % Define & Save metadata
    metadat.xyzlim = xyzbuff ;                  % xyzlimits
    metadat.reorient_faces = false ;            % if some normals are inverted
    metadat.normal_shift = QS.normalShift ;             % normal push, in pixels, along normals defined in data XYZ space
    metadat.texture_axis_order = [2 1 3] ; % for Hand Hand HistGFP,  % texture space sampling
        
    % Psize is the linear dimension of the grid drawn on each triangular face
    Options.PSize = 5 ;
    Options.EdgeColor = 'none';
    QS.getRotTrans() ;
    Options.Rotation = QS.APDV.rot ;
    Options.Translation = QS.APDV.trans ;
    Options.Dilation = QS.APDV.resolution ;
    % outward, inward --> reasonable values are ~[3,3]
    Options.numLayers = [1, -5];  % at layerSpacing 2, 2 marches ~0.5 um 
    Options.layerSpacing = 2 ;
    
    % Save it
    save(metafn, 'metadat', 'Options')
else
    disp('Loading options for endoderm texturepatch')
    load(metafn, 'metadat', 'Options')
end

% Use first timepoint's intensity limits throughout
QS.clearTime()
QS.setDataLimits(QS.xp.fileMeta.timePoints(1), 1.0, 99.9)
% Plot on surface for all TP 
options = metadat ;
options.overwrite = false ;
options.plot_dorsal = true ;
options.plot_ventral = true ;
options.plot_right = true ;
options.plot_left = true ;
options.plot_perspective = true ;
QS.plotSeriesOnSurfaceTexturePatch(options, Options)





%% Morphsnakes demo evolution visualizeMeshEvolution orthoviews data planes
options = struct() ;
options.plot_growth = false ;
options.plot_evolution = true ;
options.growth_t0 = 45 ;
options.overwrite = false ;
options.preview = false ;
QS.data.adjustlow = 1 ;
QS.data.adjusthigh = 99 ;
options.timePoints = xp.fileMeta.timePoints;
options.plot_texturePatches = false ;
options.lambda = 0.001 %35 ;
options.normalShift = 10 ; 
options.viewAngles = [-20, 20] ;
options.xMin = -51 ;
options.xMax = 268 ;
options.yMin = -90 ;
options.yMax = 90 ;
options.zMin = -75 ;
options.zMax = 80 ;
options.y0 = -15 ;
QS.visualizeMeshEvolution(options)

%% EXTRACT CENTERLINES
% Skip if already done
% Note: these just need to be 'reasonable' centerlines for topological
% checks on the orbifold cuts.
exponent = 1.0 ;
res = 4.0 ; 
cntrlineOpts.overwrite = overwrite_centerlines ;     % overwrite previous results
cntrlineOpts.overwrite_ims = overwrite_centerlineIms ;     % overwrite previous results
cntrlineOpts.weight = 0.1;               % for speedup of centerline extraction. Larger is less precise
cntrlineOpts.exponent = exponent ;       % how heavily to scale distance transform for speed through voxel
cntrlineOpts.res = res ;                 % resolution of distance tranform grid in which to compute centerlines
cntrlineOpts.preview = false ;           % preview intermediate results
cntrlineOpts.reorient_faces = false ;    % not needed for our well-constructed meshes
cntrlineOpts.dilation = 0 ;        % how many voxels to dilate the segmentation inside/outside before path computation
cntrlineOpts.skipErrors = false ;
cntrlineOpts.permuteDataAxisOrder = 'xyz' ;
cntrlineOpts.timePoints = fliplr(QS.xp.fileMeta.timePoints)
% Note: this takes about 400s per timepoint for res=2.0
%
QS.extractCenterlineSeries(cntrlineOpts)
disp('done with centerlines')


%% Fix flip in Y for centerlines
% aux_fix_flip

%% Identify anomalies in centerline data
% Skip if already done
overwrite_idAnomClines = true ;
idOptions.ssr_thres = 15 ;  % distance of sum squared residuals in um as threshold
idOptions.overwrite = overwrite_idAnomClines ;
QS.generateCleanCntrlines(idOptions) ;


%% Cylinder cut mesh
% Skip if already done
overwrite_endcapOpts = true ;
if overwrite_endcapOpts || ~exist(QS.fileName.endcapOptions, 'file')
    endcapOpts = struct( 'adist_thres', 20, ...  % 20, distance threshold for cutting off anterior in pix
                'pdist_thres', 14, ...  % 15-20, distance threshold for cutting off posterior in pix
                'tref', 1) ;  % reference timepoint at which time dorsal-most endcap vertices are defined
    QS.setEndcapOptions(endcapOpts) ;
    % Save the options to disk
    QS.saveEndcapOptions() ;
else
    % load endcapOpts
    QS.loadEndcapOptions() ;
    endcapOpts = QS.endcapOptions ;
end

clearvars methodOpts
methodOpts.overwrite = overwrite_endcapOpts ;  % recompute sliced endcaps
methodOpts.save_figs = true ;   % save images of cutMeshes along the way
methodOpts.preview = false  ;     % display intermediate results
QS.sliceMeshEndcaps(endcapOpts, methodOpts) ;

%% Clean Cylinder Meshes
% May skip if already done
cleanCylOptions.overwrite = true ;
cleanCylOptions.save_ims = true ;
QS.cleanCylMeshes(cleanCylOptions)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORBIFOLD -> begin populating Qs.dir.mesh/gridCoords_nUXXXX_nVXXXX/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
% Overwriting options
overwrite_pullbacks = false ;
overwrite_cutMesh = false ;
overwrite_spcutMesh = false ;
% Plotting params
washout2d = 0.5 ;

%% Iterate Through Time Points to Create Pullbacks ========================
% Skip if already done
% outcutfn = fullfile(cutFolder, 'cutPaths_%06d.txt') ;
overwrite_spcutMesh = false ;
overwrite_cutMesh = false ;
tp2do = [t0_for_phi0:max(xp.fileMeta.timePoints), ...
    fliplr(min(xp.fileMeta.timePoints):(t0_for_phi0-1))] ;

for tt = tp2do
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    QS.setTime(tt) ;
    
    %----------------------------------------------------------------------
    % Create the Cut Mesh
    %----------------------------------------------------------------------
    cutMeshfn = sprintf(QS.fullFileBase.cutMesh, tt) ;
    cutPathfn = sprintf(QS.fullFileBase.cutPath, tt) ;
    if ~exist(cutMeshfn, 'file') || ~exist(cutPathfn, 'file') || overwrite_cutMesh
        if exist(cutMeshfn, 'file')
            disp('Overwriting cutMesh...') ;
        else
            disp('cutMesh not found on disk. Generating cutMesh... ');
        end
        options = struct() ;
        options.preview = false ;
        options.t0 = t0_for_phi0 ;
        QS.generateCurrentCutMesh(options)
        disp('Saving cutP image')
        % Plot the cutPath (cutP) in 3D
        QS.plotCutPath(QS.currentMesh.cutMesh, QS.currentMesh.cutPath)
        compute_pullback = true ;
    else
        fprintf('Loading Cut Mesh from disk... ');
        QS.loadCurrentCutMesh()
        compute_pullback = ~isempty(QS.currentMesh.cutPath) ;
    end
    
    spcutMeshOptions = struct() ;
    spcutMeshOptions.t0_for_phi0 = t0_for_phi0 ;
    spcutMeshOptions.overwrite = overwrite_spcutMesh ;
    spcutMeshOptions.save_phi0patch = false ;
    spcutMeshOptions.iterative_phi0 = true ;
    spcutMeshOptions.smoothingMethod = 'none' ;
    QS.plotting.preview = false ;
    QS.generateCurrentSPCutMesh([], spcutMeshOptions) ;
    
    % Compute the pullback if the cutMesh is ok
    if compute_pullback || ~exist(sprintf(QS.fullFileBase.im_sp, tt), 'file')
        pbOptions.overwrite = false ;
        pbOptions.generate_uv = false ;
        pbOptions.generate_sphi = true ;
        pbOptions.generate_uphi = false ;
        pbOptions.generate_relaxed = false ;
        QS.generateCurrentPullbacks([], [], [], pbOptions) ;
    else
        disp('Skipping computation of pullback')
    end
    clear Options IV
        
    % Save SMArr2D (vertex positions in the 2D pullback) -----------------
    % disp(['Saving meshStack to disk: ' mstckfn])
    % save(mstckfn, 'meshStack') ;
    % 
    % %% Save SMArr2D (vertex positions in the 2D pullback) -----------------
    % disp(['Saving spmeshStack to disk: ' spmstckfn])
    % if generate_sphi_coord
    %     save(spmstckfn, 'spmeshStack') ;
    % end
end
clearvars t cutP spcutMesh spvn3d ss pp uv tileCount slin plin plotfn
clearvars IVloaded IV uphi sphi
disp('Done with generating spcutMeshes and cutMeshes')

%% FLIP Y for spcutMesh mcline
% DO ONLY ONCE IF NEEDED (no longer needed for datasets processed after
% 12-03-2020
% for tidx = 1:length(QS.xp.fileMeta.timePoints)
%     tp = QS.xp.fileMeta.timePoints(tidx) ;
%     load(sprintf(QS.fullFileBase.spcutMesh, tp), 'spcutMesh') ;
%     spcutMesh.mcline(:, 2) = - spcutMesh.mcline(:, 2) ;
%     spcutMesh.avgpts(:, 2) = - spcutMesh.avgpts(:, 2) ;
%     save(sprintf(QS.fullFileBase.spcutMesh, tp), 'spcutMesh') ;
% end


%% Preview results ========================================================
check = false ;
if check
    aux_preview_results
end

%% FIND THE FOLDS SEPARATING COMPARTMENTS =================================
% Skip if already done
options = struct() ;
options.overwrite = true ;
options.preview = true ;
options.first_tp_allowed = [-1, 120, -1] ;  % enforce that no folds before this tp
options.guess123 = [.29 .55 .80] ;
options.max_wander = 5 ;
QS.identifyFolds(options)
disp('done')

%% Circularlity of constrictions -- todo?
% options = struct() ;
% QS.measureFoldRadiiVariance(options)

%% COMPUTE MESH SURFACE AREA AND VOLUME ===================================
% Skip if already done
% Note: doing this after fold identification so that t0 is defined for
% plotting purposes
options = struct() ;
options.overwrite = true ;
QS.measureSurfaceAreaVolume(options)
disp('done')

% RECOMPUTE WRITHE OF MEANCURVE CENTERLINES ==============================
% Skip if already done
options = struct() ;
options.overwrite = true ;
QS.measureWrithe(options)
disp('done')

%% Plot fancy "cross-section" view of centerlines
options = struct() ;
QS.plotClineXSections(options)

%% Compute surface area and volume for each compartment ===================
% Skip if already done
options = struct() ;
options.overwrite = true ;
QS.measureLobeDynamics(options) ;

%% plot length, area, and volume for each lobe ============================
% Skip if already done
options = struct() ;
options.overwrite = true ;
QS.plotLobes(options)

% Plot motion of avgpts & DVhoops at folds in yz plane over time =========
% Skip if already done
overwrite_lobeims = true ;
QS.plotConstrictionDynamics(overwrite_lobeims) ;
disp('done')

%% SMOOTH MEAN CENTERLINE RADIUS ==========================================
% Skip if already done
% todo: rework this section so that it takes place after smoothing meshes
% overwrite = false ;
% aux_smooth_avgptcline_radius_before_mesh_smoothing(overwrite, ...
%     QS.xp.fileMeta.timePoints, QS.fullFileBase.spcutMesh, ...
%     QS.fullFileBase.clineDVhoop, radiusImDir, ...
%     QS.APDV.rot, QS.APDV.trans, QS.APDV.resolution, ...
%     QS.plotting.xyzlim_um, QS.nU, QS.nV)

%% Smooth the sphi grid meshes in time ====================================
% Skip if already done
options = struct() ;
options.overwrite = false ;
options.width = 4 ;  % before 2020-12-09, this was 6.
QS.smoothDynamicSPhiMeshes(options) ;

%% Fix one field of spcutMeshSm / spcutMeshSmRS / spcutMeshSmRSC -- 2020-12-21
% for tidx = 1:length(QS.xp.fileMeta.timePoints)
%     tt = QS.xp.fileMeta.timePoints(tidx) ;
%     
%     spcutMeshBase = QS.fullFileBase.spcutMesh ;
%     spcutMeshSmBase = QS.fullFileBase.spcutMeshSm ;
%     spcutMeshSmRSBase = QS.fullFileBase.spcutMeshSmRS ;
%     spcutMeshSmRSCBase = QS.fullFileBase.spcutMeshSmRSC ;
% 
%     % Load the files
%     load(sprintf(spcutMeshBase, tt), 'spcutMesh') ;
%     load(sprintf(spcutMeshSmBase, tt), 'spcutMeshSm') ;
%     load(sprintf(spcutMeshSmRSBase, tt), 'spcutMeshSmRS') ;
%     load(sprintf(spcutMeshSmRSCBase, tt), 'spcutMeshSmRSC') ;
%     
%     spcutMeshSm.u = spcutMesh.sphi ;
%     spcutMeshSmRS.u = spcutMesh.sphi ;
%     spcutMeshSmRSC.u = spcutMesh.sphi(1:nU*(nV-1), :) ;
%     assert(size(spcutMeshSmRSC.u, 1) == size(spcutMeshSmRSC.v, 1))
%     
%     % Resave the files
%     disp(['Saving ' sprintf(spcutMeshSmBase, tt)])
%     save(sprintf(spcutMeshSmBase, tt), 'spcutMeshSm') ;
%     save(sprintf(spcutMeshSmRSBase, tt), 'spcutMeshSmRS') ;
%     save(sprintf(spcutMeshSmRSCBase, tt), 'spcutMeshSmRSC') ;
%             
% end

%% Plot the time-smoothed meshes
% Skip if already done
options.overwrite = false ;
QS.plotSPCutMeshSmRS(options) ;

% Images for publication/presentation on method & coordinate system
% Skip if already done
% Create coordinate system charts visualization using smoothed meshes

% QS.coordSystemDemo()

% COMPUTE MEAN AND GAUSSIAN CURVATURES OF SMOOTHED MESHES
% Skip if already done
options = struct() ;
options.overwrite = false ;
QS.measureCurvatures(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Redo Pullbacks with time-smoothed meshes ===============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Skip if already done
disp('Create pullback using S,Phi coords with time-averaged Meshes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = fliplr(QS.xp.fileMeta.timePoints)
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = QS.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
    QS.setTime(tt) ;
    
    % OPTION 1: Keep constant luminosity throughout, modify default 
    % intensity limits.
    % if tidx == 1        
    %     % Use first timepoint's intensity limits throughout
    %     QS.setDataLimits(QS.xp.fileMeta.timePoints(1), 1.0, 99.995)
    % end
    % QS.data.adjustlow 
    % QS.data.adjusthigh
    
    % OPTION 2 : Adjust intensity to scale from timepoint to timepoint
    adjustlow = 1.00 ;         % floor for intensity adjustment
    adjusthigh = 99.9 ;        % ceil for intensity adjustment (clip)
    QS.data.adjustlow = adjustlow ;
    QS.data.adjusthigh = adjusthigh ;
    
    % Establish custom Options for MIP
    pbOptions = struct() ;
    pbOptions.overwrite = false ;
    pbOptions.numLayers = [0 0] ; % previously [7, 7] ;  % previously [5,5]
    pbOptions.layerSpacing = 0.75 ;
    pbOptions.generate_rsm = true ;
    pbOptions.generate_spsm = false ;
    pbOptions.generate_sphi = false ;
    pbOptions.generate_uvprime = false ;
    pbOptions.generate_ruvprime = false ;
    QS.data.adjustlow = adjustlow ;
    QS.data.adjusthigh = adjusthigh ;
    QS.generateCurrentPullbacks([], [], [], pbOptions) ;
end

%% TILE/EXTEND SMOOTHED IMAGES IN Y AND RESAVE =======================================
% Skip if already done
options = struct() ;
options.overwrite = false ;
% options.coordsys = 'spsm' ;
% QS.doubleCoverPullbackImages(options)
options.coordsys = 'rsm' ;
QS.doubleCoverPullbackImages(options)
disp('done')

%% Compute whether pullback is isothermal -> metric images
options = struct() ;
options.overwrite = false ;
options.makeRawMetricComponentFigures = false ;
options.lambda_mesh = 0.002 ;
options.coordSys = 'ricci' ; % 'spsm_rs'
QS.plotMetric(options) ;
% !!! todo: continue here for diagnostic
% Note: demo_FundFormNematic

%% Compare 2nd fund form director to radon of image
options = struct() ;
options.method = 'pullback' ; 
options.coordSys = 'spsmre' ;
QS.measurePolarity(options) ;
% inside measurePolarity: comparePolarityHopfDifferential()

%% Add radii to spcutMeshSm's in post
% 
% for tidx = 1:length(QS.xp.fileMeta.timePoints)
%     tp = QS.xp.fileMeta.timePoints(tidx) ;
%     % Load this timepoint spcutMeshSm
%     QS.setTime(tp)
%     spcutMeshSm = QS.getCurrentSPCutMeshSm() ;
%     spcutMeshSmRS = QS.getCurrentSPCutMeshSmRS() ;
%     spcutMeshSmRSC = QS.getCurrentSPCutMeshSmRSC() ;
%     
%     % Make avgpts in pixel space (not RS)
%     fprintf('Resampling uvgrid3d curves in pix...\n')
%     nU = spcutMeshSm.nU ;
%     nV = spcutMeshSm.nV ;
%     curves3d_pix = reshape(spcutMeshSm.v, [nU, nV, 3]) ;
%     c3d_dsv_pix = zeros(size(curves3d_pix)) ;  % in units of pix
%     avgpts_pix = zeros(nU, 3) ;
%     radius_pix = zeros(nU, nV) ;
%     for i=1:nU
%         % Note: no need to add the first point to the curve
%         % since the endpoints already match exactly in 3d and
%         % curvspace gives a curve with points on either
%         % endpoint (corresponding to the same 3d location).
%         c3d_dsv_pix(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d_pix(i, :, :)), nV, true) ;
%         if vecnorm(squeeze(c3d_dsv_pix(i, 1, :)) - squeeze(c3d_dsv_pix(i, end, :))) > 1e-7
%             error('endpoints do not join! Exiting')
%         end
%         % Drop the final endpoint in the mean pt determination
%         avgpts_pix(i, :) = mean(squeeze(c3d_dsv_pix(i, 1:end-1, :)), 1) ; 
%         radius_pix(i, :) = vecnorm(squeeze(curves3d_pix(i, :, :)) - avgpts_pix(i, :), 2, 2) ;
%     end
% 
%     % uvpcutMesh.raw.avgpts_pix = avgpts_pix ;
%     % uvpcutMesh.raw.radius_pix = radius_pix ;
%     spcutMeshSm.avgpts_um = QS.xyz2APDV(avgpts_pix) ;
%     spcutMeshSm.radius_um = radius_pix * QS.APDV.resolution ;
%     
%     % Add to RS and RSC versions
%     spcutMeshSmRS.avgpts_um = spcutMeshSm.avgpts_um ;
%     spcutMeshSmRS.radius_um = spcutMeshSm.radius_um ;
%     % RSC
%     spcutMeshSmRSC.avgpts_um = spcutMeshSm.avgpts_um ;
%     spcutMeshSmRSC.radius_um = spcutMeshSm.radius_um(:, 1:end-1) ;
%     
%     % Save spcutMeshSm / RS / RSC
%     save(sprintf(QS.fullFileBase.spcutMeshSm, QS.currentTime), ...
%         'spcutMeshSm')
%     save(sprintf(QS.fullFileBase.spcutMeshSmRS, QS.currentTime), ...
%         'spcutMeshSmRS')
%     save(sprintf(QS.fullFileBase.spcutMeshSmRSC, QS.currentTime), ...
%         'spcutMeshSmRSC')
%     
%     if mod(tidx, 10) == 0 
%         clf; set(gcf, 'visible', 'on')
%         trisurf(triangulation(spcutMeshSmRSC.f, spcutMeshSmRSC.v), 'edgecolor', 'none')
%         axis equal 
%         title(num2str(tidx))
%         pause(1)
%         cla
%         trisurf(triangulation(spcutMeshSmRS.f, spcutMeshSmRS.v), 'edgecolor', 'k')
%         axis equal
%         title([num2str(tidx) ' open'])
%         pause(0.0001)
%     end
% end

%% Cell Segmentation
options = struct() ;
options.overwrite = false ;
options.overwriteImages = false;
options.timePoints = [93:15:263] ;
QS.generateCellSegmentation2D(options) 

options = struct() ;
options.overwrite = false ;
options.overwriteImages = true;
options.timePoints = [93:15:263] ;
QS.processCorrectedCellSegmentation2D(options) 

options = struct() ;
options.timePoints = [93:15:212] ;
options.overwrite = false ;
options.overwriteImages = false ;
options.correctedSegmentation = true ;
QS.generateCellSegmentation3D(options) 
options = struct() ;
options.correctedSegmentation = true ;
options.timePoints = [93:15:212] ;
options.overwrite = false  ;
QS.plotSegmentationStatisticsLobes(options)

%%
options = struct() ;
options.timePoints = [93:15:263] ;
options.overwrite = false ;
QS.generateCellSegmentationPathlines3D(options)
% QS.estimateIntercalationRate(options)

%% 
options = struct() ;
QS.visualizeSegmentationPatch(options)

%% Measure Cell density
% Skip if already done
options = struct() ;
options.overwrite = false ;
options.preview = false ;
QS.measureCellDensity('nuclei', options)
QS.plotCellDensity(options) ;

%% Cell density kymograph
% Skip if already done
options = struct() ;
options.overwrite = false ;
options.timePoints = 0:85 ;
QS.plotCellDensityKymograph(options)

%% CREATE PULLBACK STACKS =================================================
% Skip if already done
disp('Create pullback stack using S,Phi coords with time-averaged Meshes');
% Load options
overwrite = true ;
optionfn = fullfile(QS.dir.im_r_sme_stack, 'spcutMeshSmStackOptions.mat') ;
if ~exist(optionfn, 'file') || overwrite
    spcutMeshSmStackOptions.layer_spacing = 0.5 / QS.APDV.resolution ; % pixel resolution roughly matches xy
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
QS.generateSPCutMeshSmStack(spcutMeshSmStackOptions)

%% TRAIN IN ILASTIK ON MIDGUT TISSUE TO GET THICKNESS
% Read thickness training output, train with first channel fg, second bg,
% third optional folded-over tissue. Name output <filename>_thickness.h5

%% MEASURE THICKNESS ======================================================
% Measure thickness from images of imFolder_sp_r_sme_stack (extendedLUT_smoothed)
thicknessOptions = struct() ;
QS.measureThickness(thicknessOptions)

%% DUMP OR LOAD HERE [break]
clearvars fold_ofn fig1exist fig2exist id idx mcline nsmM prevcline tmp
clearvars fig1exist fig1fn fig2exist fp1fn fp2fn fp3fn fp4fn IV 
clearvars n_inward n_outward TV2D TV3D 
clearvars layer_spacing 
dumpfn = fullfile(meshDir, 'orbifold_dump_before_piv.mat') ;
save(dumpfn)
load(dumpfn)
clearvars dumpfn


%% PERFORM PIV ON PULLBACK MIPS ===========================================
% % Compute PIV in PIVLab
% % ---------------------
% % Open PIVLab
% % Select all frames in meshDir/PullbackImages_010step_sphi/smoothed_extended/
% % Select Sequencing style 1-2, 2-3, ... 
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
disp('Loading PIV results...')
tmp = load(fullfile(QS.dir.piv, 'piv_results.mat')) ;


%% Measure velocities =============================================
disp('Making map from pixel to xyz to compute velocities in 3d for smoothed meshes...')
options = struct() ;
options.overwrite = true ;
options.preview = false ;
options.show_v3d_on_data = false ;
options.save_ims = true ;
QS.measurePIV3d(options) ;


%% AUTOCORRELATIONS
% overwrite_autocorrelations = false ;
% do_acorr = false ;
% redo_acorr = ~exist(fullfile(pivSimAvgDir, 'autocorr_velocities.png'), 'file') ;
% if (redo_acorr || overwrite_autocorrelations) && do_acorr
%     aux_autocorrelations
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First do very simpleminded averaging of velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options.overwrite = false ;
% options.plot_vxyz = true ;
% QS.timeAverageVelocitiesSimple('1x', options)
% %% VELOCITY PLOTS
% options.overwrite = true ;
% options.plot_vxyz = false ;
% options.invertImage = true ;
% QS.plotTimeAvgVelSimple(options)
% %% Divergence and Curl (Helmholtz-Hodge)
% options = struct() ;
% options.overwrite = false ;
% options.samplingResolution = '1x' ; 
% options.averagingStyle = 'simple' ;
% QS.helmholtzHodge(options) ;
% %% Measure Compressibility (div(v), 2*vn*H, and gdot)
% options = struct() ;
% options.overwrite = false ;
% options.plot_Hgdot = false ;
% options.plot_flows = true ;
% options.plot_factors = true ;
% options.plot_kymographs = true ;
% options.plot_kymographs_cumsum = true ;
% options.plot_correlations = true ;
% options.plot_gdot_correlations = false ;
% options.plot_gdot_decomp = true ;
% options.lambda_mesh = 0.002 ;
% options.lambda = 0.02 ;
% options.lambda_err = 0.03 ;
% options.samplingResolution = '1x'; 
% QS.measureMetricKinematics(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOUBLE RESOLUTION Simple average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measure velocities at 2*nU x 2*nV resolution
% options = struct() ;
% options.overwrite = false ;
% options.preview = false ;
% options.show_v3d_on_data = false ;
% options.save_ims = true ;
% QS.generateSPCutMeshSm2x(options.overwrite) ;
% QS.measurePIV3d2x(options) ;
% %% time average the velocities
% options.overwrite = false ;
% options.plot_vxyz = true ;
% QS.timeAverageVelocitiesSimple('2x', options) ;
% %% Velocity plots for 2x
% options.overwrite = false ;
% options.plot_vxyz = false ;
% options.invertImage = true ;
% options.samplingResolution = '2x'; 
% QS.plotTimeAvgVelSimple(options)
% %% Divergence and Curl (Helmholtz-Hodge) for 2x
% options = struct() ;
% options.overwrite = true ;
% options.samplingResolution = '2x' ;
% QS.helmholtzHodgeSimple(options) ;
% %% doubleResolution compressibility
% options = struct() ;
% options.overwrite = false ;
% options.plot_Hgdot = false ;
% options.plot_flows = true ;
% options.plot_factors = true ;
% options.plot_kymographs = true ;
% options.plot_kymographs_cumsum = true ;
% options.plot_correlations = true ;
% options.plot_gdot_correlations = false ;
% options.plot_gdot_decomp = true ;
% options.samplingResolution = '2x'; 
% QS.measureMetricKinematics(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lagrangian dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pullback pathline time averaging of velocities
options = struct() ;
options.overwrite = true ;
options.preview = false ;
QS.timeAverageVelocities(options)
% Velocity plots for pathline time averaging 
options.overwrite = false ;
options.plot_vxyz = false ;
options.invertImage = true ;
options.averagingStyle = 'Lagrangian'; 
options.samplingResolution = '1x'; 
QS.plotTimeAvgVelocities(options)
% Divergence and Curl (Helmholtz-Hodge) for Lagrangian
options = struct() ;
options.overwrite = false ;
options.samplingResolution = '1x' ;
options.averagingStyle = 'Lagrangian' ;
options.lambda = 0 ;
options.lambda_mesh = 0 ; 
QS.helmholtzHodge(options) ;

% Compressibility & kinematics for Lagrangian
options = struct() ;
options.overwrite = false ;
options.samplingResolution = '1x'; 
options.lambda_mesh = 0 ;
options.lambda_err = 0.00 ;
options.lambda = 0.00 ;
options.nmodes = 7 ;  %% bandwidth filtering
options.zwidth = 2 ; 
QS.measureMetricKinematics(options)

%% Metric Kinematics Kymographs & Correlations -- Bandwidth Filtered
options = struct() ;
options.overwrite = false ;
options.overwrite_timePoints = false ;
options.plot_Hgdot = false ;
options.plot_flows = false ;
options.plot_factors = false ;
options.plot_kymographs = true ;
options.plot_kymographs_cumsum = true ;
options.plot_kymographs_cumprod = true ;
options.plot_correlations = true ;
options.plot_spaceMaps = true ;
options.plot_raw_correlations = true ;
options.plot_gdot_correlations = false ;
options.plot_gdot_decomp = false ;
options.lambda = 0.00 ;
options.lambda_err = 0.00 ;
options.lambda_mesh = 0 ;
options.nmodes = 7 ;  %% bandwidth filtering
options.zwidth = 2 ; 
options.climit = 0.3 ;
QS.plotMetricKinematics(options)

%% Measure shear stresses against gradients in pressure
% eta \nabla_i d^i_j = \nalba_j p, 
% where d^i_j = 1/2(\nabla_i v_j + \nalba_j v_i) - b_ij vn
options = struct() ;
options.overwrite = true ;
options.zwidth = 0 ;
QS.measureStokesForces(options) ;
options.plot_kymographs = true ;
QS.plotStokesForces(options) ;  % !!! continue here



%% Pullback pathlines connecting Lagrangian grids
options = struct() ;
options.overwrite = false ;
options.preview = false ;
options.debug = false ; 
QS.measurePullbackPathlines(options)

%% Query velocities along pathlines
options = struct() ;
options.overwrite = false ;
options.preview = false ;
QS.measurePathlineVelocities(options)
% plot the pathline velocities 
options = struct() ;
options.overwrite = false ;
options.gridTopology = 'triangulated' ;
QS.plotPathlineVelocities(options)

% Measure Pathline Kinematics
options = struct() ;
options.overwrite = false ;
options.lambda = 0 ;
options.lambda_error = 0 ;
options.lambda_mesh = 0 ;
options.nmodes = 7 ;
options.zwidth = 1 ; 
QS.measurePathlineMetricKinematics(options)
% todo: resume this BLOCK !!!

% Plot Pathline Kinematics
options = struct() ;
options.overwrite = true ;
options.plot_kymographs = true ;
options.plot_kymographs_cumsum = false ;
options.plot_kymographs_cumprod = false ;
options.plot_correlations = true ;
options.plot_fold_kinematics = true ;
options.plot_lobe_kinematics = true ;
options.climit = 0.30 ;
QS.plotPathlineMetricKinematics(options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create ricci mesh at t0 to measure Beltrami coefficient in pathlines
options = struct() ;
options.overwrite = false ;
options.climit = 1 ;
options.coordSys = 'ricci' ;
QS.measureBeltramiCoefficient(options) ;

%% Generate all Beltramis from all Riccis & plot aspect ratio over time
options = struct() ;
options.overwrite = false ;
QS.computeRicciMeshes(options) %!!!

%% Generate all Beltramis from all Riccis & Plot aspect ratio for isothermal PB over time
% Note: this really shouldn't be necessary, as we show in detail now
options = struct() ;
QS.measureBeltramiCoefficientPullbackToPullback(options) ;

%% Check mu by computing for all times and get the same thing as 
% embedding back to same Reference tp

for tp = QS.xp.fileMeta.timePoints
    disp(['t = ', num2str(tp)])
    QS.setTime(tp)
    try
        rmesh = QS.loadCurrentRicciMesh() ;
    catch
        disp(['could not load Ricci mesh for t=' num2str(tp)])
    end
    
    % Compute mu from mapping from PB(0) to PB(1)
    
end


%% Compare to linearized description
options = struct() ;
options.overwrite = false ;
options.overwriteImages = false ;
options.coordSys = 'ricci' ;
QS.compareBeltramiToLinearizedConstriction(options) ;

%% Compare to nonlinear isothermal description
options = struct() ;
options.overwrite = true ;
options.overwriteImages = false ;
options.coordSys = 'ricci' ;
QS.compareBeltramiToConstriction(options) ;

%% UVPrime option (ignore)
master_uvprime_module
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Measure twist (d v_phi / d zeta)
options = struct() ;
options.overwrite = false ;
QS.measureTwist(options)

%% Measure required stress pattern
options = struct() ;
options.overwrite = false ;
QS.measureStressPattern(options) ;

%% Measure metric strain rate (epsilon=gdot/2)
% option = struct() ;
% options.overwrite = true ;
% options.preview = false ;
% QS.measureMetricStrainRate(options) 

%% Strain rate (epsilon = 1/2 (djvi+divj) -vn bij)
options = struct() ;
options.overwrite = false ;
options.overwriteImages = false ;
options.preview = false ;
options.lambda = 0.000 ;
options.lambda_error = 0 ;
options.lambda_mesh = 0 ;
options.nmodes = 7 ;
options.zwidth = 1 ; 
QS.measureStrainRate(options) 

%% Plot time-averaged strain rates in 3d on mesh
options = struct() ;
% error('This does not look great -- tweak how the averaging is done to be Lagrangian?')
options.overwrite = false ;
options.preview = false ;
QS.plotStrainRate3DFiltered(options) 

%% Kymograph strain rates
options = struct() ;
options.overwrite = false ;
options.skipTimePoint = true ;
options.clim_trace = 0.05 ;
options.clim_deviatoric = 0.05 ;
QS.plotStrainRate(options)

% Measure strain rate along pathlines
options = struct() ;
options.overwrite = false ;
options.overwriteImages = false ;
options.plot_dzdp = false ;
QS.measurePathlineStrainRate(options)
% Integrate strain rates along pathlines (noisy)
% options = struct() ;
% options.overwrite = true ;
% options.overwriteImages = true ;
% options.plot_dzdp = false ;
% options.median_filter_strainRates = false ;
% options.climitInitial = 0.05 ;
% options.climitRamp = 0.01 ;
% options.climitRatio = 1 ;
% QS.measurePathlineIntegratedStrain(options)

% Pathline strain rate plots
options = struct() ;
options.overwrite = false ;
options.plot_kymographs = false ;
options.plot_fold_strainRate = true ;
options.plot_lobe_strainRate = true ;
options.climit = 0.05 ;
options.climitWide = 1.0 ;
QS.plotPathlineStrainRate(options)

%% (REPEAT FROM ABOVE) Pullback pathlines connecting Lagrangian grids
% options = struct() ;
% options.overwrite = false ;
% options.preview = false ;
% options.debug = false ; 
% QS.measurePullbackPathlines(options)
%% Measure strain along pathlines -- note this is from pathlines, not integrating rates
options = struct() ;
options.overwrite = true ;
options.overwriteImages = false ;
options.plot_dzdp = false ;
options.climitInitial = 0.05 ;
options.climitRamp = 0.01 ;
options.climitRatio = 1 ;
QS.measurePathlineStrain(options)
QS.plotPathlineStrain(options)

%% Measure coarse-grained bond contraction and dilation in zeta=s/L, phi
% todo: consider putting this in Ricci map frame instead of (s,phi) frame
options = struct() ;
options.overwrite = true ;
options.overwriteImages = true ;
QS.measureDxDyStrainFiltered(options) ;

%% Output gg and bb in fixed frame (t0 Lagrangian zeta phi frame) -- this is done in measurePathlineStrain()
% options = struct() ;
% options.overwrite = false ;
% measurePathlineMetric(QS, options)

%% Simulate this experiment
options = struct() ;
QS.simulateNES(options)

%% Measure surface area growh of lobes and folds in Eulerian frame
options = struct() ;
options.overwrite = false ;
options.preview = false ;
QS.measureEulerianMetricDynamics(options)

%% Plot texture map in 3d with velocities
for i = 1:size(vM, 1)
    t = time(i) ;
    vtqfn0 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_perspective.png']) ;
    vtqfn1 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_dorsal.png']) ;
    vtqfn2 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_ventral.png']) ;
    vtqfn3 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_lateralL.png']) ;
    vtqfn4 = fullfile(pivSimAvgImQDir, [sprintf('%04d', time(i)) '_lateralR.png']) ;

    % Plot the tangential velocity as quiver on top of the embedding,
    % colored by normal velocity (all in one in embedding!)
    if ~exist(vtqfn0, 'file') || overwrite_vsm_plots
        disp(['Creating ' vtqfn])
        tic

        v2dsmum_ii = squeeze(v2dsmMum(i, :, :)) ;
        vnsm_ii = squeeze(vnsmM(i, :, :)) ;

        % Load the data for the current time point ------------------------
        xp.setT=ime(t) ;
        % Load 3D data for coloring mesh pullback
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();

        fig = figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off') ;

        % Psize is the linear dimension of the grid drawn on each triangular face
        Options.PSize = 5;
        Options.EdgeColor = 'none';
        Options.Rotation = rot ;
        Options.FaceScalarField = vn_interp ;

        % Raw stack data
        IV = xp.stack.image.apply();
        IV = imadjustn(IV{1});

        % First args are physical vertices, then texture faces (same number as 
        % physical faces, but connectivity can be different), then texture
        % vertices, which can also be different. The vertices are in row, column, 
        % page format, so x and y get flipped. IV is the texture volume.
        % Options.PSize 
        [ patchIm, imref, zeroID, MIP, SIP ] = ...
            texture_patch_3d( mesh.f, mesh.v, ...
            mesh.f, mesh.v(:, [2 1 3]), IV, Options );

        % Add quiver to image
        % Downsample by a factor of qsubsample
        xx = piv.x{i}(1, :)' ;
        yy = piv.y{i}(:, 1) ;
        ww = length(xx) ;
        hh = length(yy) ;
        vx = reshape(v2dsmum_ii(:, 1), [ww, hh]) ;
        vy = reshape(v2dsmum_ii(:, 2), [ww, hh]) ;
        QX = imresize(vx, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
        QY = imresize(vy, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
        xq = 1:qsubsample:ww ;
        yq = 1:qsubsample:hh ;
        [xg, yg] = meshgrid(xq, yq) ;
        h2 = quiver(xg(:), yg(:), QX(:), QY(:), 0) ;

        % Format the figure and axis
        xlim(xyzmin(1, :))
        ylim(xyzmin(2, :))
        zlim(xyzmin(3, :))
        title(['$t = $' num2str(t) ' min'], 'Interpreter', 'Latex', 'Color', 'white') 
        xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
        ylabel('DV position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
        zlabel('lateral position [$\mu$m]', 'Interpreter', 'Latex', 'Color', 'white')
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

        % Make background black
        set(gca,'Color','k')
        set(gcf, 'InvertHardCopy', 'off');
        set(gcf, 'Color', 'k')

        % Make tick labels white
        labeltype = {'XTickLabel', 'YTickLabel', 'ZTickLabel'} ;
        for q = 1:2
            % get the current tick labeks
            ticklabels = get(gca, labeltype{q});
            % prepend a color for each tick label
            ticklabels_new = cell(size(ticklabels));
            for i = 1:length(ticklabels)
                ticklabels_new{i} = ['\color{white} ' ticklabels{i}];
            end
            % set the tick labels
            set(gca, labeltype{q}, ticklabels_new);
        end

        % Capture all five views (perspective + DVLR
        disp(['saving figure...' num2str(t, '%06d')])
        saveas(fig, vtqfn0)
        view(0, 90)
        % DORSAL
        saveas(fig, vtqfn1)
        view(0, 270)
        % VENTRAL
        saveas(fig, vtqfn2)

        % Lateral views
        view(0, 0)
        % make white z tick labels
        q = 3;
        % get the current tick labeks
        ticklabels = get(gca, labeltype{q});
        % prepend a color for each tick label
        ticklabels_new = cell(size(ticklabels));
        for i = 1:length(ticklabels)
            ticklabels_new{i} = ['\color{white} ' ticklabels{i}];
        end
        % set the tick labels
        set(gca, labeltype{q}, ticklabels_new);

        % LEFT
        saveas(fig, vtqfn3)
        view(0, 180)
        % RIGHT
        saveas(fig, vtqfn4)
        close all
        toc
        clear Options IV
    end

end

% Check the orientation of the phasebar
% imshow(im)
% hold on;
clf
xsize = 2000 ;
ysize = 2000 ;
step = 100 ;
imsz = [ xsize ysize ] ;
[xx, yy] = meshgrid(1:step:xsize, 1:step:ysize) ;
ucheck = xx - xsize * 0.5 ;
vcheck = yy - ysize * 0.5 ;
vangle = mod(atan2(vcheck, -ucheck), 2* pi) ;
imshow(im * washout2d + max(im) * (1-washout2d)) ;
hold on ;
h2 = imagesc(xx(:), yy(:), vangle) ;
hold on ;
quiver(xx, yy, ucheck, vcheck, 5) ;
colormap phasemap
phasebar
caxis([0, 2*pi])
axis on
set(gcf, 'visible', 'on')
waitfor(gcf)


%% Simpleminded streamlines from velocity scaled by dilation
% % Build 3d grid of positions and velocities
% % Assume that x0, y0 are fixed for all time
% ntimes = length(piv.x) ;
% x0 = piv.x{i} ;
% y0 = piv.y{i} ;
% xlen = size(v2dsmM, 1) ;
% ylen = size(v2dsmM, 2) ;
% vdat = v2dsmM(:) ;
% vxdat = reshape(vdat(:, 1), [ntimes, xlen, ylen]) ;
% vydat = reshape(vdat(:, 2), [ntimes, xlen, ylen]) ;
% vzdat = ones(size(vydat)) ; % we march through time at 1 index / timestep
% % Define positions we track through the streamlines
% startx = x0(1:200:end) ;
% starty = y0(1:200:end) ;
% startz = zeros(size(starty)) ;
% streamline(x0, y0, z0, vxdat, vydat, vzdat, startx, starty, startz)
% view(2)
% save(gcf, fullfile(pivDir, 'streamlines.png')) 
% error('break')

