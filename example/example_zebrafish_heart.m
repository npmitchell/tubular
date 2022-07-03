%% TubULAR Analysis Pipeline: Zebrafish Heart (3D + time)
% by NPMitchell and Dillon Cislo

% This is a pipeline to analyize the dynamic tube-like surface of a growing
% zebrafish heart in (3+1)D. During the developmental stage at which this
% data was procured, the zebrafish heart is a tube with cylindrical
% topology. Our surface extraction technology, however, renderes the fitted
% surfaces as elongated and spherical in topology. This pipeline will
% demonstarate:
%
% (1) How to restore the cylindrical topology by the removal of 'endcaps'
% (2) How to slice open and unwrap the tube-like surface into the plane
% (3) How to construct time-dependent 'surface Lagrangian' coordinates
% (4) How to apply the DEC to analyze flow on the shape-shifting surface
%
% TODO:
%   - WE REALLY NEED THE EXPERIMENT METADATA (time interval, space units
%   etc.)
%   - Your axis order choices are crap
%   - Should probably just make a new MeshLab script that would deprecate
%   the mesh simplification step post surface detection

% Add necessary directories to the path
tubularDir = '/mnt/data/code/tubular';
addpath(genpath(tubularDir));

% Add optional external code to the path
% rmpath(genpath('mnt/data/code/gptooolbox'));
addpath(genpath('/mnt/data/code/gptoolbox'));

%% TubULAR Pipeline Initialization ========================================

% We start by clearing the memory and closing all figures
clear; close all; clc;

% The directory containing the zebrafish heart data
dataDir = '/mnt/data/tubular_test/zebrafish_heart/';

% The directory where project files will be generated and saved
[projectDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
cd(projectDir);

% Define TubULAR master settings
overwriteSettings = true;
if (~exist(fullfile(projectDir, 'masterSettings.mat'), 'file') || overwriteSettings)

    stackResolution = [.3524 .3524 2];  % resolution in spaceUnits per pixel
    nChannels = 1;                      % how many channels is the data (ex 2 for GFP + RFP)
    channelsUsed = 1;                   % which channels are used for analysis
    timePoints = 1:30; % 1:50;          % timepoints to include in the analysis
    ssfactor = 4;                       % subsampling factor
    flipy = false ;                     % whether the data is stored inverted relative to real position in lab frame
    timeInterval = 1;                   % physical interval between timepoints
    timeUnits = 'min';                  % physical unit of time between timepoints
    spaceUnits = '$\mu$m';              % physical unit of time between timepoints
    fn = 'Stack_Repeat_014_Time_%03d';	% filename string pattern
    set_preilastikaxisorder = 'xyzc';	% data axis order for subsampled h5 data (ilastik input)
    swapZT = true;                      % whether to swap the z and t dimensions
    
    masterSettings = struct( ...
        'stackResolution', stackResolution, ...
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
        'nV', 100 );
    
    disp(['Saving master settings to ' ...
        fullfile(projectDir, 'masterSettings.mat')]);
    
    save(fullfile(projectDir, 'masterSettings.mat'), 'masterSettings');
    
    clear stackResolution nChannels channelsUsed timePoints
    clear ssfactor flipy timeInterval timeUnits spaceUnits fn
    clear set_preilastikaxisorder swapZT
    
    disp('Saving masterSettings to ./masterSettings.mat')
    
else
    
    disp(['Loading master settings from ' ...
        fullfile(projectDir, 'masterSettings.mat')]);
    
   load(fullfile(projectDir, 'masterSettings.mat'), 'masterSettings');
    
    
end

%% ************************************************************************
% *************************************************************************
%      PART 1: SURFACE DETECTION USING ImSAnE's 'integralDetector'
% *************************************************************************
% *************************************************************************

% Add ImSAnE to the path
addpath(genpath('/mnt/data/code/imsane_for_git/imsane'));
rmpath('/mnt/data/code/imsane_for_git/imsane/external/CGAL_Code/IsotropicRemeshing');

%% Initialize ImSAnE Project ==============================================

% We start by creating the experiment object, which holds the project
% metadata and serves as a front-end for a number of tasks, including data
% loading, detection, fitting, etc.  When instantiating the experiment
% object we must indicate the path of the project directory and the data
% directory, either by passing them as arguments to the constructor or by
% picking them manually from a dialog box

xp = project.Experiment(projectDir, dataDir);

% Set File Metadata -------------------------------------------------------
% First we set the metadata pertaining to the raw data files in the
% structure 'fileMeta'.  ImSAnE assumes that individual time points are
% saved as separate image stacks and that filenames in a time series are
% identical up to a integer specifying the timepoint.
%
% The following file metadata information is required:
%
% * 'directory'         , the project directory (full path)
% * 'dataDir'           , the data directory (full path)
% * 'filenameFormat'    , fprintf type format spec of file name
% * 'timePoints'        , list of times available stored as a vector
% * 'stackResolution'   , stack resolution in microns, e.g. [0.25 0.25 1]
%
% The following file metadata information is optional:
%
% * 'imageSpace'        , bit depth of image, such as uint16 etc., defined
%                         in Stack class
% * 'stackSize'         , size of stack in pixels per dimension 
%                         [xSize ySize zSize]
% * 'swapZT'            , set=1 if time is 3rd dimension and z is 4th
% * 'nChannels'         , The number of channels in the raw data

fileMeta                    = struct();
fileMeta.dataDir            = dataDir;
fileMeta.filenameFormat     = [ masterSettings.fn, '.tif' ];
fileMeta.nChannels          = masterSettings.nChannels;
fileMeta.timePoints         = masterSettings.timePoints;
fileMeta.stackResolution    = masterSettings.stackResolution;
fileMeta.swapZT             = masterSettings.swapZT;

xp.setFileMeta(fileMeta);

% Set Experiment Metadata -------------------------------------------------
% Next we set additional information regarding our experiment as fields in
% the 'expMeta' structure.
%
% The following project metadata information is required:
%
% * 'channelsUsed'      , the channels used, e.g. [1 3] for RGB
% * 'channelColor'      , mapping from element in channels used to RGB=123
% * 'dynamicSurface'    , Boolean, false: static surface
% * 'detectorType'      , name of detector class
% * 'fitterType'        , name of fitter class
%
% The following project meta data information is optional:
%
% * 'description'     , string describing the data set
% * 'jitterCorrection', Boolean, false: No fft based jitter correction 

expMeta                     = struct();
expMeta.channelsUsed        = masterSettings.channelsUsed;
expMeta.channelColor        = 1;
expMeta.description         = 'A beating Zebrafish heart';
expMeta.dynamicSurface      = 1;
expMeta.jitterCorrection    = 0;
expMeta.fitTime             = fileMeta.timePoints(1);
expMeta.detectorType        = 'surfaceDetection.morphsnakesDetector';
expMeta.fitterType          = 'surfaceFitting.meshWrapper';

xp.setExpMeta(expMeta);

% Initialize New Experiment -----------------------------------------------
% Finally we call initNew(), which reads the stack size from the first 
% available time point, then initializes fitter and detector and creates 
% fitOptions and detectOptions based on their defaults.

xp.initNew();

clear fileMeta expMeta

%% Load First Time Point ==================================================
% First we load the data for the fittingtime point. This is not necessary
% for surface detection. However, it is required to properly store
% information generated in subsequent steps (i.e. generating the surface of
% interest object) and we might as well get it out of the way now and also
% demonstrate some visualization functionality offered by ImSAnE

xp.loadTime( xp.expMeta.fitTime );
xp.rescaleStackToUnitAspect();

% The data in xp.stack is stored as a 'Stack' object. The easiest way
% to look at a cross section of the data is using getSlice

% imshow( xp.stack.getSlice( 'z', 90 ), [] );

%% Set Surface Detection Objects ==========================================
% We now attempt to detect the surface of interest.  In this pipeline, we
% use the 'integralDetector'.  This detector is essentially a wrapper for
% the morphological snakes external code module.  The input to this module
% is a sub-sampled Ilastik pixel probability map.  Given this map, the
% method attempts to segmemt the image volume into distinct regions by
% minimizing the following (schematic) energy functional
%
% Etot = Est + Ep + Ein + Eout
%
% Est amounts to a surface tension, Ep amounts to a pressure, Ein tries to
% make all voxels inside the segmented region have similar intensities, and
% Eout tries to make all of the voxels outside the segmented region have
% similar intensities. It is assumed that the region of interest is a
% closed volume.
%
% See ImSAnE's 'surfaceDetection.integralDetector' for option documentation

msls_detOpts_fn = fullfile(projectDir, 'msls_detectOpts.mat');

if exist(msls_detOpts_fn, 'file')
    
    load(msls_detOpts_fn, 'detectOptions');
    
else
    
    % The name of the meshlab script that generates a mesh from the point
    % cloud
    mlxprogram = 'laplace_surface_rm_resample30k_reconstruct_LS3_1p2pc_ssfactor4.mlx';
    mlxprogram = fullfile( '/mnt/data/code/meshlab_codes/', mlxprogram );
    
    detectOptions = struct();
    detectOptions.channel = masterSettings.channelsUsed(1);
    detectOptions.ssfactor = masterSettings.ssfactor;
    detectOptions.niter = 200;
    detectOptions.niter0 = 400;
    detectOptions.lambda1 = 1;
    detectOptions.lambda2 = 1;
    detectOptions.pressure = 0.1;
    detectOptions.tension = 2;
    detectOptions.pre_pressure = 0;
    detectOptions.pre_tension = 0;
    detectOptions.post_pressure = 0;
    detectOptions.post_tension = 0;
    detectOptions.exit_thres = 1e-6;
    detectOptions.foreGroundChannel = masterSettings.channelsUsed(1);
    detectOptions.fileName = sprintf( masterSettings.fn, xp.expMeta.fitTime );
    detectOptions.mslsDir = fullfile( projectDir, 'MorphSnakesOutput' );
    detectOptions.ofn_ls = 'msls_DP_';
    detectOptions.ofn_ply = 'mesh_DP_ms_';
    detectOptions.ms_scriptDir = '/mnt/data/code/morphsnakes_wrapper/morphsnakes_wrapper';
    detectOptions.timepoint = xp.expMeta.fitTime;
    detectOptions.zdim = 3;
    detectOptions.ofn_smoothply = 'mesh_DP_';
    detectOptions.mlxprogram = fullfile( projectDir, mlxprogram );
    detectOptions.init_ls_fn = 'empty_string';
    detectOptions.run_full_dataset = false;
    detectOptions.radius_guess = 200;
    detectOptions.dset_name = 'exported_data';
    detectOptions.center_guess = '30,125,125';
    detectOptions.save = true;
    detectOptions.plot_mesh3d = false;
    detectOptions.dtype = 'h5';
    detectOptions.mask = 'none';
    detectOptions.mesh_from_pointcloud = true;
    detectOptions.prob_searchstr = '_Probabilities.h5';
    detectOptions.preilastikaxisorder = masterSettings.set_preilastikaxisorder;
    detectOptions.ilastikaxisorder = 'cxyz'; %'cxyz';
    detectOptions.physicalaxisorder = 'yxzc';
    detectOptions.include_boundary_faces = true;
    detectOptions.smooth_with_matlab = 0.01;
    detectOptions.pythonVersion = '';
    
end

xp.setDetectOptions( detectOptions );

%% Batch Ilastik Pre-Processing ===========================================
% In this section we load the formatted time series data and render the
% .tif files into sub-sampled .h5 files appropriate for segmentation in
% Ilastik
%
% -------------------------------------------------------------------------
% NOTE: This section is slow - it requires us to load and format the raw
% data stack of each individual time point.  If you are just following
% along with the tutorial, skip this section - it has already been done for
% you!
%--------------------------------------------------------------------------

% Create the directory that will hold the Ilastik prep files
if ~exist(fullfile(projectDir, 'prepFiles'), 'dir')
    mkdir(fullfile(projectDir, 'prepFiles'));
end 

for t = xp.fileMeta.timePoints
    
    prepFileName = fullfile( projectDir, sprintf( [ 'prepFiles/' ...
        masterSettings.fn ], t ) );
    
    if ~exist(prepFileName, 'file')
        
        xp.loadTime(t);
        xp.rescaleStackToUnitAspect();
        
        detectOptions.fileName = prepFileName;
        detectOptions.timepoint = xp.currentTime;
        
        xp.setDetectOptions(detectOptions);
        xp.detector.prepareIlastik(xp.stack);
        
    end
    
end

clear prepFileName

%% Batch Creation of Ilastik Prediction ===================================
% At this point the user should switch to Ilastik and train to create a
% pixel probability map using the sub-sampled prep files.  Don't work too
% hard here! You dont have to train on every single time point.  Training
% can be performed on a subset of the time series (or even on a single time
% point if you're bold) and used to batch process the rest!
%
% The rest of this pipeline assumes that the pixel probability maps are
% stored in a folder 'projectDir/probFiles'.  Make sure to set the export
% preferences in Ilastik accordingly!

if ~exist(fullfile(projectDir, 'probFiles'), 'dir')
   mkdir(fullfile(projectDir, 'probFiles'));
end 

%--------------------------------------------------------------------------
% NOTE: This section is slow - it requires us to complete an Ilastik
% training workflow and batch process each time point in the data set.  If
% you are just following along with the tutorial, skip this section - it
% has already been done for you!
%--------------------------------------------------------------------------

%% Create MorphSnakes Level Sets from Ilastik Probabilities Output ========

% Detect the Surface at the First Time Point ------------------------------

% The name of the probability file used for detection
probFileName = fullfile( projectDir, ...
    sprintf( [ 'probFiles/' masterSettings.fn ], ...
    xp.fileMeta.timePoints(1) ) );

% Update detector options
detectOptions.fileName = probFileName;
detectOptions.timepoint = xp.fileMeta.timePoints(1);

xp.setDetectOptions( detectOptions );

% Calling detectSurface() runs the morphsnakes detection protocol
xp.detectSurface()

% Detect the Surfaces at all Subsequent Time Points -----------------------

for t = xp.fileMeta.timePoints(2:end)
    
    tidx = xp.tIdx(t);
    
    % The name of the probability file used for detection
    probFileName = fullfile( projectDir, ...
        sprintf( [ 'probFiles/' masterSettings.fn ], t ) );
    
    % Update detector options
    detectOptions.fileName = probFileName;
    detectOptions.timepoint = t;
    
    % Use the previous time points surface as an initial guess
    detectOptions.init_ls_fn = sprintf([detectOptions.ofn_ls '%06d.h5'], ...
        xp.fileMeta.timePoints(tidx-1));
    
    xp.setDetectOptions( detectOptions );
    
    % Calling detectSurface() runs the morphsnakes detection protocol
    xp.detectSurface();
     
end

clear probFileName

%% Generate Point Cloud/Mesh from Level Set ===============================
% This section extracts a simplied point set from the volumetric level set
% corresponding to the cylindrical surface of the heart and then constructs
% a watertight mesh using Poisson surface reconstruction. The
close all; clc;

% Create a directory to hold the output point cloud files
if ~exist(fullfile(projectDir, 'objFiles'),'dir')
    mkdir(fullfile(projectDir, 'objFiles'));
end

% Create a directory to hold the output mesh files
if ~exist(fullfile(projectDir, 'meshFiles'),'dir')
    mkdir(fullfile(projectDir, 'meshFiles'));
end

for t = xp.fileMeta.timePoints()
    
    fprintf('Now processing timepoint T=%d\n', t);
    
    pointCloudFileName = fullfile( projectDir, ...
        sprintf('objFiles/pointCloud_T%03d.obj', t) );
    
    meshFileName = fullfile( projectDir, ...
        sprintf('meshFiles/mesh_T%03d.off', t) );
    
    if ~exist(pointCloudFileName, 'file')
        
        %------------------------------------------------------------------
        % Extract the implicit level set as a 3D binary array
        %------------------------------------------------------------------
        
        % The file name of the current time point
        implicitLSFile = fullfile( detectOptions.mslsDir, ...
            sprintf( 'msls_DP_%06d.h5', t ) );
        
        % The 3D binary array
        bwLS = h5read( implicitLSFile, '/implicit_levelset' );
        
        % Keep only the largest connected component from the implicit
        % levelset
        bwCC = bwconncomp(bwLS);
        [~, maxID] = max(cellfun(@numel, bwCC.PixelIdxList));
        bwLS = false(size(bwLS));
        bwLS(bwCC.PixelIdxList{maxID}) = true;
        
        % Fill in holes in the volume
        bwLS = imfill( bwLS, 'holes' );
        
        % Re-assess level set quality
        bwCC = bwconncomp(bwLS);
        if( bwCC.NumObjects ~= 1 )
            warning('Level set is multiply connected!');
        end
        
        bwCC = bwconncomp(~bwLS);
        if ( bwCC.NumObjects ~= 1 )
            warning( 'Level set complement is multiply connected!');
        end
        
        % Skeletonize each cross section of the binary volume to produce a
        % minimal point set
        bwLS0 = permute(bwLS, [2 3 1]);
        bwLS = false(size(bwLS0));
        for i = 1:size(bwLS,3)
            progressbar(i,size(bwLS,3));
            bwLS(:,:,i) = bwskel(bwLS0(:,:,i));
        end
        bwCC = bwconncomp(bwLS);
        
        % Extract the (x,y,z)-locations of the level set points
        % (in pixel space)
        clear P
        [ P(:,1), P(:,2), P(:,3) ] = ind2sub( size(bwLS), ...
            vertcat(bwCC.PixelIdxList{:}) );
        
        % Rescale boundary pixel locations to the size of the image stack
        P = (P-1) * detectOptions.ssfactor + 1;
        
        % Downsample/denoise the point set
        PC = pointCloud(P);
        PC = pcdownsample(PC, 'nonuniformGridSample', 18); %150
        % PC = pcdenoise(PC, 'NumNeighbors', 24, 'Threshold', 2.5); % 24
        P = PC.Location;
        
        fprintf('Number of points = %d\n', size(P,1));
        
        %------------------------------------------------------------------
        % Estimate Point Set Normals
        %------------------------------------------------------------------
        
        param = struct();
        param.estimation_procedure = 2; % Uses PCA normal estimation
        param.number_of_neighbors = 18;
        orient_neighbors = 12;
        
        % P0 = P;
        [ PN, P, oriented ] = point_set_normals( P, param, ...
            orient_neighbors );
        
        % Check to see if any points were removed during normal orientation
        if oriented
            warning(['Point set normal orientation procedure', ...
                'removed points from the input set!']);
        end
        
        % View results ----------------------------------------------------
        % ssf = false( size(P,1), 1 );
        % sst = false( size(P,1), 1 );
        % sst(1:5:end) = true;
        % % sst = true( size(P, 1), 1 );
        %
        % scatter3(  P([sst; ssf; ssf]), P([ssf; sst; ssf]), ...
        %     P([ssf; ssf; sst]), 'filled' );
        % hold on
        % quiver3( P([sst; ssf; ssf]), P([ssf; sst; ssf]), ...
        %     P([ssf; ssf; sst]), ...
        %     PN([sst; ssf; ssf]), PN([ssf; sst; ssf]), ...
        %     PN([ssf; ssf; sst]), ...
        %     1, 'LineWidth', 2 );
        % axis equal
        % hold off
        
        %------------------------------------------------------------------
        % Upsample Point Set
        %------------------------------------------------------------------
        
        numPoints = 3000; % The number of output points
        sharpnessAngle = 90; % Controls the sharpness of the results
        edgeSensitvity = 0; % Controls sampling density near inferred edges
        neighborRadius = eps; % Initial neighbor radius
        
        [PU, PNU] = upsample_point_set( P, PN, numPoints, sharpnessAngle, ...
            edgeSensitvity, neighborRadius );
        
        % View results ----------------------------------------------------
        % ssf = false( size(PU,1), 1 );
        % sst = true( size(PU, 1), 1 );
        % % sst = false( size(P,1), 1 );
        % % sst(1:50:end) = true;
        %
        % scatter3(  PU([sst; ssf; ssf]), PU([ssf; sst; ssf]), ...
        %     PU([ssf; ssf; sst]), 'filled' );
        % hold on
        % quiver3( PU([sst; ssf; ssf]), PU([ssf; sst; ssf]), ...
        %     PU([ssf; ssf; sst]), ...
        %     PNU([sst; ssf; ssf]), PNU([ssf; sst; ssf]), ...
        %     PNU([ssf; ssf; sst]), ...
        %     1, 'LineWidth', 2 );
        % axis equal
        % hold off
        
        %------------------------------------------------------------------
        % Output Point Cloud/Mesh Files
        %------------------------------------------------------------------
        
        % Output point cloud file
        OBJ = struct();
        OBJ.vertices = PU;
        OBJ.vertices_normal = PNU;
        OBJ.objects(1).type='f';
        OBJ.objects(1).data.vertices = [];
        OBJ.objects(1).data.normal = [];
        write_wobj(OBJ, pointCloudFileName );
        
    else
        
        % Read in point cloud from file
        [PU, ~, ~, ~, PNU, ~] = readOBJ(pointCloudFileName);

    end
    
    if ~exist(meshFileName, 'file')
        
        %------------------------------------------------------------------
        % Generate Mesh From Point Cloud
        %------------------------------------------------------------------
        
        poisson_surface_reconstruction(PU, PNU, 'poisson_mesh.off');
        
        % Read in result from file and delete original file
        [V, F, ~, ~, ~] = readOFF('poisson_mesh.off');
        delete('poisson_mesh.off');
        
        % Clean mesh surface ----------------------------------------------
        
        % Isotropically remesh the surface
        tar_length = 10;
        num_iter = 5;
        protect_constraints = false;
        [F, V, ~, ~] = isotropic_remeshing( F, V, ...
            tar_length, num_iter, protect_constraints );
        
        % Attempt to remove localized mesh spikes by Laplacian relaxation
        V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
            'uniform', [], 2, 'implicit', 1000);
        
        % Try to remove self-intersections
        [intersects, ~] = mesh_self_intersection_3d(F, V);
        intCount = 0;
        while intersects
            
            [V, F] = clean_mesh(V, F, 'MinDist', 0, 'MinArea', 0, ...
                'MinAngle', 0, 'SelfIntersections', 'remove', ...
                'SmallTriangles', 'remove');
            
            % Peform a another isotropic remeshing
            [F, V, ~, ~] = isotropic_remeshing(F, V, ...
                tar_length, num_iter, protect_constraints);
            
            % Another round of spike relaxation
            V = relax_mesh_spikes(F, V, deg2rad(60), pi/2, ...
                'uniform', [], 2, 'implicit', 1000);
            
            [intersects, ~] = mesh_self_intersection_3d(F, V);
            
            intCount = intCount + 1;
            if intCount > 20
                error('Unable to remove self-intersections');
            end
            
        end
        
        % Another round of isotropic remeshing
        [F, V, ~, ~] = isotropic_remeshing(F, V, ...
            tar_length, num_iter, protect_constraints);
        
        % Smooth the entire mesh fixing the boundary
        V = laplacian_smooth(V, F, 'cotan', [], 0.025, 'implicit', V, 10);
        
        % Peform a final isotropic remeshing
        [F, V, ~, ~] = isotropic_remeshing(F, V, ...
            tar_length, num_iter, protect_constraints);
        
        % Mesh quality checks ---------------------------------------------
        
        E = edges(triangulation(F, V));
        
        numBdy = numel(DiscreteRicciFlow.compute_boundaries(F));
        if (numBdy ~= 0)
            error( ['Mesh at time point T = %d has %d ' ...
                'boundary components'], t, numBdy );
        end
        
        eulerChi = size(F,1) + size(V,1) - size(E,1);
        if (eulerChi ~= 2)
            error( ['Mesh at time point T = %d is not ' ...
                'a topological sphere'], t );
        end
        
        [intersects, ~] = mesh_self_intersection_3d(F, V);
        if intersects
            error( ['Mesh at time point T = %d contains ' ...
                'self-intersections'], t );
        end
        
        % hold on
        % trisurf(triangulation(F, V), 'FaceColor', 0.8 * ones([1 3]));
        % hold off
        
        %------------------------------------------------------------------
        % Output Mesh Files
        %------------------------------------------------------------------
        
        writeOFF(meshFileName, V, F);
        
    end
    
end

clear pointCloudFileName meshFileName implicitLSFile
clear bwLS bwCC maxID bwLS0 P PC param P0 oriented sst ssf
clear numPoints sharpnessAngle edgeSensitivity neighborRadius OBJ
clear tar_length num_iter protect_constraintsintCount

clear E numBdy eulerChi intersects
clear PU PNU F V

%% View Point Cloud Intersection with Data ================================
close all; clc;

% Choose time point to be visualized
t = xp.currentTime;

% Load raw image data
if ((t ~= xp.currentTime) || isempty(xp.stack))
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
end

% Load point cloud from file
pointCloudFileName = fullfile( projectDir, ...
    sprintf('objFiles/pointCloud_T%03d.obj', t) );

v3D = readOBJ(pointCloudFileName);

% Choose the slice to be visualized ---------------------------------------
sliceNum = 100;
sliceDim = 'z';

if strcmpi(sliceDim, 'x')
    sliceDimNum = 1;
    plotDim = [2 3];
elseif strcmpi(sliceDim, 'y')
    sliceDimNum = 2;
    plotDim = [3 1];
elseif strcmpi(sliceDim, 'z')
    sliceDimNum = 3;
    plotDim = [1 2];
else
    error('Invalid slice dimension');
end

% Find point cloud intersection with the image plane ----------------------
incThresh = 1;
inPlane = abs(v3D(:, sliceDimNum) - sliceNum) <= incThresh;

% Generate figure ---------------------------------------------------------
imshow( xp.stack.getSlice( sliceDim, sliceNum ), [] );

hold on

% scatter3(v3D(inPlane, 1), v3D(inPlane, 2), v3D(inPlane, 3), 'filled', 'r');
scatter(v3D(inPlane, plotDim(1)), v3D(inPlane, plotDim(2)), 'filled', 'r');

hold off

clear t tidx pointCloudFileName v3D sliceNum sliceDim sliceDimNum
clear incThresh inPlane plotDim

%% View Mesh Intersection With Data (Uses 'gptoolbox') ====================
close all; clc;

% Choose time point to be visualized
t = xp.currentTime;

% Load raw image data
if ((t ~= xp.currentTime) || isempty(xp.stack))
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
end

% Load mesh from file
meshFile = fullfile( projectDir, sprintf('meshFiles/mesh_T%03d.off', t) );
[V, F, ~, ~, ~] = readOFF(meshFile);
% V = V(:, [3 2 1]); % Re-order axes for visualization purposes

% mesh = read_ply_mod('MorphSnakesOutput/mesh_DP_000001.ply');
% trisurf(triangulation(mesh.f, mesh.v)); axis equal
% [F, V, ~, ~] = isotropic_remeshing(mesh.f, mesh.v, 10, 5);
% oldV = V;
% V = oldV(:, [2 1 3]);

% Choose the slice to be visualized ---------------------------------------
sliceNum = 500; % round(1345/2);
sliceDim = 'x';

if strcmpi(sliceDim, 'x')
    sliceDimNum = 1;
    plotDim = [2 3];
elseif strcmpi(sliceDim, 'y')
    sliceDimNum = 2;
    plotDim = [3 1];
elseif strcmpi(sliceDim, 'z')
    sliceDimNum = 3;
    plotDim = [1 2];
else
    error('Invalid slice dimension');
end

% Calculate the contour intersection of the mesh with the plane -----------
[U, UF, ~, ~, ~, ~] = slice_isolines(V, F, V(:, sliceDimNum), ...
    sliceNum, 'Manifold', true);
UE = edges(triangulation(UF, U));

inSlice = find(abs(U(:, sliceDimNum) - sliceNum) < 1e-10);
C = U(inSlice, :);
CE = UE(all(ismember(UE, inSlice),2), :);
CE = changem(CE, 1:size(C,1), inSlice);

% Generate figure ---------------------------------------------------------
imshow( xp.stack.getSlice( sliceDim, sliceNum ), [] );
hold on

% scatter3(C(:,1), C(:,2), C(:,3), 'filled', 'r');

for i = 1:size(CE,1)
    
    plot([C(CE(i,1), plotDim(1)); C(CE(i,2), plotDim(1))], ...
        [C(CE(i,1), plotDim(2)); C(CE(i,2), plotDim(2))], ...
        '-or', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    
end
    
hold off

% trisurf(triangulation(UF, U), 'FaceColor', 0.8*[1 1 1]);
% hold on
% 
% scatter3(C(:,1), C(:,2), C(:,3), 'filled', 'r');
% 
% for i = 1:size(CE,1)
%     
%     plot3([C(CE(i,1),1); C(CE(i,2),1)], ...
%         [C(CE(i,1),2); C(CE(i,2),2)], ...
%         [C(CE(i,1),3); C(CE(i,2),3)], '-ob', 'LineWidth', 2);
%     
% end
%     
% hold off
% 
% xlabel('x'); ylabel('y'); zlabel('z');
% axis equal

clear t tidx  meshFIleC CE UE C plotDim
clear AC GC U UF inSlice sliceNum sliceDim sliceDimNum

%% View the Data on the Mesh Surface ======================================
close all; clc

% Choose time point to be visualized
t = xp.currentTime;

% Load raw image data
if ((t ~= xp.currentTime) || isempty(xp.stack))
    xp.loadTime(t);
    xp.rescaleStackToUnitAspect();
end

% Load mesh from file
meshFile = fullfile( projectDir, sprintf('meshFiles/mesh_T%03d.off', t) );
[V, F, ~, ~, ~] = readOFF(meshFile);

% The raw data image stack
IV = xp.stack.image.apply();
IV = double(IV{1})/65535;
% IV{1} = imadjustn(IV{1});
% IV{2} = imadjustn(IV{2});
% IV{3} = zeros(size(IV{1})); % RGB requires 3 channels

% Texture patch options
textureOptions = struct();
textureOptions.EdgeColor = 'none';
textureOptions.numLayers = [2, 2];
textureOptions.layerSpacing = 2;

texture_patch_3d( F, V, F, V(:, [2 1 3]), IV, textureOptions);
axis equal
colormap bone

clear IV textureOptions t V F meshFile

%% Inspect All Meshes in 3D ===============================================
close all; clc;

for t = xp.fileMeta.timePoints()
    
    % Load the mesh
    meshFile = fullfile( projectDir, sprintf('meshFiles/mesh_T%03d.off', t) );
    [V, F, ~, ~, ~] = readOFF(meshFile);
    
    % Plot the mesh in 3d. Color here by Z coordinate
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), ...
        V(:, 3), 'EdgeColor', 'none', 'FaceAlpha', 0.5)
    
    hold on
    scatter3(V(any(V<=0.5,2),1), V(any(V<=0.5,2),2), V(any(V<=0.5,2),3), ...
        'filled', 'r');
    hold off
    
    title(['t=' num2str(t)])
    axis equal
    view([-80 0]);
    % view(1)
    pause(0.1)
    
end

clear meshFile V F

%% ************************************************************************
% *************************************************************************
%               PART 2: TubULAR -- SURFACE PARAMETERIZATION
% *************************************************************************
% *************************************************************************

%% Set Up TubULAR Directories =============================================

if ~exist('TubULAR_Results', 'dir')
    
    mkdir(fullfile(projectDir, 'TubULAR_Results'));
    
    % Re-write meshes into the format used by TubULAR
    for t = xp.fileMeta.timePoints()
        
        % Load the mesh
        meshFile = fullfile( projectDir, ...
            sprintf('meshFiles/mesh_T%03d.off', t) );
        [V, F, ~, ~, ~] = readOFF(meshFile);
        
        VN = per_vertex_normals(V, F, 'Weighting', 'angle');
        VN = normalizerow(VN);
        
        outputMeshFile = fullfile( projectDir, ...
            sprintf(['TubULAR_Results/'...
            detectOptions.ofn_smoothply '%06d.ply'], t) );
        
        plywrite_with_normals(outputMeshFile, F, V, VN);
        
    end
    
    clear t meshFile V F VN outputMeshFile
    
end


%% Initialize the TubULAR Object ==========================================
close all; clc;

tubOpts = struct();
tubOpts.meshDir = fullfile(projectDir, 'TubULAR_Results');	% Director where meshes produced in previous part reside. All TubULAR results will be stored relative to this directory
tubOpts.flipy = masterSettings.flipy;               % Set to true if data volume axes are inverted in chirality wrt physical lab coordinates               
tubOpts.timeInterval = masterSettings.timeInterval; % Spacing between adjacent timepoints in units of timeUnits
tubOpts.timeUnits = masterSettings.timeUnits;       % Units of time, so that adjacent timepoints are timeUnits * timeInterval apart
tubOpts.spaceUnits = masterSettings.spaceUnits;     % Units of space in LaTeX, for ex '$mu$m' for micron
tubOpts.nU = masterSettings.nU;         % How many points along the longitudinal axis to sample surface
tubOpts.nV = masterSettings.nV;         % How many points along the circumferential axis to sample surface
tubOpts.t0 = xp.fileMeta.timePoints(1);	% Reference timepoint used to define surface-Lagrangian and Lagrangian measurements
tubOpts.normalShift = 2;        % Additional dilation acting on surface for texture mapping
tubOpts.a_fixed = 2.0;          % Fixed aspect ratio of pullback images. Setting to 1.0 is most conformal mapping option.
tubOpts.adjustlow = 1.00;       % Floor for intensity adjustment
tubOpts.adjusthigh = 99.9;      % ceil for intensity adjustment (clip)
tubOpts.phiMethod = 'curved3d';	% Method for following surface in surface-Lagrangian mapping [(s,phi) coordinates]
tubOpts.lambda_mesh = 0;        % Smoothing applied to the mesh before DEC measurements
tubOpts.lambda = 0;             % Smoothing applied to computed values on the surface
tubOpts.lambda_err = 0;         % Additional smoothing parameter, optional

disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, tubOpts) ;
disp('done defining TubULAR instance')

clear tubOpts

%% Inspect All Meshes in 3D (TubULAR Style) ===============================

for tp = xp.fileMeta.timePoints
    
    % Load the mesh
    meshfn = sprintf( tubi.fullFileBase.mesh, tp );    
    mesh = read_ply_mod(meshfn);
    
    % Plot the mesh in 3D
    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        mesh.v(:, 3), 'Edgecolor', 'none', 'Facealpha', 0.5);
    title(['t=' num2str(tp)]);
    axis equal
    view([-80 0]);
    
    pause(0.1)
    
end

clear tp meshfn mesh

%% Obtain APDV Coordinates of the Surface =================================
% There are ftwo options for obtaining these coordinates. 
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
% Here we use option 1.
% Define global orientation frame (for viewing in canonical frame)
% Compute APDV coordinate system
alignAPDVOpts = struct() ;
alignAPDVOpts.overwrite = false ;
tubi.computeAPDVCoords(alignAPDVOpts) ;

%% Determine A-P Points/Endcaps From Point Cloud ==========================
close all; clc;

figure('units', 'normalized', ...
    'outerposition', [0.5 0 0.5 1],  'color', [1 1 1]);

% Store point LOCATIONS
apts = zeros(numel(xp.fileMeta.timePoints), 3);
ppts = zeros(numel(xp.fileMeta.timePoints), 3);

% Store encap vertex IDs
allAntEndCaps = cell(numel(xp.fileMeta.timePoints), 3);
allPosEndCaps = cell(numel(xp.fileMeta.timePoints), 3);

for t = xp.fileMeta.timePoints()
    
    tidx = xp.tIdx(t);
    fprintf('Processing time point T = %d\n', t);
    
    % Read in point cloud from file
    pointCloudFileName = fullfile( projectDir, ...
        sprintf('objFiles/pointCloud_T%03d.obj', t) );
    [P, ~, ~, ~, ~, ~] = readOBJ(pointCloudFileName);
    
    % Read in mesh from file
    meshFileName = sprintf( tubi.fullFileBase.mesh, t );
    mesh = read_ply_mod(meshFileName);
    V = mesh.v; F = mesh.f;
    E = edges(triangulation(F,V));
    
    % Determine vertex neighbor 1-rings
    allNIDx = cell(size(V,1), 1);
    for i = 1:size(V,1)
        nIDx = unique(E(any(E == i, 2), :));
        nIDx(nIDx == i) = [];
        allNIDx{i} = nIDx(:);
    end
    
    % Format Point Cloud for Endcap Extraction ----------------------------
    % Point cloud is matched to mesh vertices. The point cloud then
    % undergoes several rounds of topological erosion/dilation to smooth
    % out the boundaries of contiguous regions not included in the point
    % set
    
    % Point match the point cloud to mesh vertices
    PVIDx = unique(knnsearch(V, P));
    
    inPC = ismember((1:size(V,1)).', PVIDx);
    
    % Dilate the point set - points that are adjacent to points contained
    % within the set, but are not currently contained within the set are
    % added to the set
    dilateIter = 5;
    for i = 1:dilateIter
        
        % Find all points adjacent to points in the point set
        adjIDx = vertcat(allNIDx{inPC});
        
        newIDx = ~inPC & ismember((1:size(V,1)).', adjIDx);
        inPC(newIDx) = true;

    end
    
    % Erode the point set - boundary points, that are included in the point
    % set, but adjacent to points that are not included in the point set,
    % are removed from the point set
    erodeIter = 5;
    for i = 1:erodeIter
        
        % Find all points adjacent to points that are NOT in the point set
        adjIDx = vertcat(allNIDx{~inPC});
        
        bdyIDx = inPC & ismember((1:size(V,1)).', adjIDx);
        inPC(bdyIDx) = false;

    end
    
    PVIDx = find(inPC);
    
    % Determine A-P Points and Endcaps ------------------------------------
    % This method works by removing all points with distances below a user
    % specified threshold and then denoting the two largest remaining mesh
    % components as the respective endcaps. This method will FAIL for
    % pathological cases with strange connectivity (e.g. the two
    % prospective endcap regions are connected by a thin band of vertices).
    % Empirical testing for this data set show that the anterior cap ALWAYS
    % has more vertices than the posterior cap - if this were not true the
    % user could use spatial orientation information to delineate which cap
    
    % Determine the distance of each mesh vertex from the nearest vertex
    % matched to a point cloud point
    D = heat_geodesic(V, F, PVIDx, [], 'Legacy', true);

    Dthresh = 10;
    [FP, VP, ~, ~] = remove_vertex_from_mesh(F, V, D < Dthresh);
    FP(any(FP == 0, 2), :) = [];
    [FP, VP, ~, ~, ~] = remove_unreferenced_vertices_from_mesh(FP, VP);
    
    [~, VA, AIDx, ~] = remove_isolated_mesh_components(FP, VP);
    [FP, VP, ~, ~] = remove_vertex_from_mesh(FP, VP, AIDx);
    FP(any(FP == 0, 2), :) = [];
    [FP, VP, ~, ~, ~] = remove_unreferenced_vertices_from_mesh(FP, VP);
    [~, VP, ~, ~] = remove_isolated_mesh_components(FP, VP);
    
    AIDx = knnsearch(V, VA);
    PIDx = knnsearch(V, VP);
    
    % The A-P points are simply taken to be the locations of maximum
    % distance within each patch
    [~, AID] = max(D(AIDx)); AID = AIDx(AID);
    [~, PID] = max(D(PIDx)); PID = PIDx(PID);
    
    apts(tidx, :) = V(AID, :);
    ppts(tidx, :) = V(PID, :);
    
    allAntEndCaps{tidx} = AIDx(:);
    allPosEndCaps{tidx} = PIDx(:);
    
    % View Results --------------------------------------------------------
 
    trisurf(triangulation(F, V), D, ...
        'FaceColor', 'interp', 'EdgeColor', 'k');
    hold on
    scatter3(V(PVIDx,1), V(PVIDx,2), V(PVIDx,3), 'filled', 'k');
    scatter3(V(AIDx,1), V(AIDx,2), V(AIDx,3), 'filled', 'g')
    scatter3(V(PIDx,1), V(PIDx,2), V(PIDx,3), 'filled', 'c')
    scatter3(V(AID,1), V(AID,2), V(AID,3), 'filled', 'r')
    scatter3(V(PID,1), V(PID,2), V(PID,3), 'filled', 'm')
    hold off
    axis equal
    colorbar
    view([-85, -45]);
    title(sprintf('T = %d\n', t));
    
    pause(0.1);

end

clear t tidx pointCloudFileName P PN meshFileName V F E
clear allNIDx i nIDx PVIDx inPC dilateIter adjIDx newIDx
clear erodeIter bdyIDx D Dthresh FP VP FA VA AIDx PIDx AID PID ts

%% Save A-P Points for Centerline Computation =============================

apdvOpts = struct();
apdvOpts.overwrite = false;
apdvOpts.custom_apts = apts;
apdvOpts.custom_ppts = ppts;
apdvOpts.smwindow = -1;

[apts_sm, ppts_sm] = tubi.computeAPDpoints(apdvOpts);

%% Align Meshes in the APDV Global Frame and Plot Them ====================
fprintf('Aligning/plotting meshes... ');
tubi.alignMeshesAPDV(alignAPDVOpts);
fprintf('Done\n');

%% Plot All Textured Meshes in 3D (OPTIONAL: this is SLOW) ================
% Establish texture patch options
metadat = struct() ;
metadat.reorient_faces = false ;            % set to true if some mesh normals may be inverted (requires gptoolbox if true)
metadat.normal_shift = tubi.normalShift ;   % normal push, in pixels, along normals defined in data XYZ space
metadat.texture_axis_order = [1 2 3] ;      % texture space sampling. If the surface and dataspace have axis permutation, enter that here
textureOptions.PSize = 5 ;          % Psize is the linear dimension of the grid drawn on each triangular face. Set PSize > 1 for refinement of texture on each triangle of the surface triangulation. Higher numbers are slower but give more detailed images.
textureOptions.numLayers = [0, 0];  % how many layers to MIP over/bundle into stack, as [outward, inward]
textureOptions.layerSpacing = 2 ;   % Distance between layers over which we take MIP, in pixels, 

% Plot on surface for all timepoints 
tubi.plotSeriesOnSurfaceTexturePatch(metadat, textureOptions)

%% Extract Centerlines ====================================================
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
cntrlineOpts.preview = true ;           % preview intermediate results
cntrlineOpts.reorient_faces = false ;    % not needed for our well-constructed meshes
cntrlineOpts.dilation = 0 ;              % how many voxels to dilate the segmentation inside/outside before path computation
% Note: this can take about 400s per timepoint for res=2.0, so use as big a 
%   res value as possible.

fprintf('Extracting centerlines... ');
tubi.extractCenterlineSeries(cntrlineOpts)
fprintf('Done\n');

%% Identify Anomalies in Centerline Data ==================================

idOptions.ssr_thes = 15; % Distance of sum squared residuals in 'spaceUnits'
                         % as a threshold for removing spurious centerlines
                         
fprintf('Cleaning centerlines... ');
tubi.generateCleanCntrlines(idOptions);
fprintf('Done\n');

%% Generate Cylindrical Mesh ==============================================
% Transforms a topological sphere into a topological cylinder

% Look for options on disk. If not saved, define options
overwrite_endcap_options = false;
if ~exist(tubi.fileName.endcapOptions, 'file') || overwrite_endcap_options
    
    endcapOpts = struct( ...
        'custom_aidx', allAntEndCaps, ...  % Points that will be removed to create the anterior endcap
        'custom_pidx', allPosEndCaps, ...  % Points that will be removed to create the posterior endcap
        'tref', tubi.t0);   % Reference timepoint at which time dorsal-most endcap vertices are defined
    
    tubi.setEndcapOptions(endcapOpts);
    
    % Save the options to disk
    tubi.saveEndcapOptions();
    
else
    
    % load endcapOpts
    tubi.loadEndcapOptions();
    endcapOpts = tubi.endcapOptions;
    
end

methodOpts.overwrite = true;
methodOpts.save_figs = true;    % save images of cutMeshes along the way
methodOpts.preview = false;     % display intermediate results
tubi.sliceMeshEndcaps(endcapOpts, methodOpts) ;
   
clear overwrite_endcap_options

%% Clean Cylinder Meshes ==================================================
% This removes "ears" from the endcaps of the tubular meshes (cylindrical
% meshes)

cleanCylOptions = struct() ;
cleanCylOptions.overwrite = false ;

fprintf('Cleaning cylinder meshes... ');
tubi.cleanCylMeshes(cleanCylOptions)
fprintf('Done\n');

%% Generate Orbifold Cut Meshes ===========================================
% Begin populating tubi.dir.mesh/gridCoords_nUXXXX_nVXXXX/

overwrite = true;

% Iterate Through Time Points to Create Pullbacks
for tt = tubi.xp.fileMeta.timePoints
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = xp.tIdx(tt);
    
    % Load the data for the current time point
    tubi.setTime(tt) ;
    
    %----------------------------------------------------------------------
    % Create the Cut Mesh
    %----------------------------------------------------------------------
    cutMeshfn = sprintf(tubi.fullFileBase.cutMesh, tt) ;
    cutPathfn = sprintf(tubi.fullFileBase.cutPath, tt) ;
    if ~exist(cutMeshfn, 'file') || ~exist(cutPathfn, 'file') || overwrite
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
        pbOptions = struct() ;
        tubi.generateCurrentPullbacks([], [], [], pbOptions) ;
    else
        disp('Skipping computation of pullback')
    end
        
end
disp('Done with generating spcutMeshes and cutMeshes')

%% Inspect Coordinate System Charts Using (s,phi) Coordinate system ('sp')
options = struct() ;
options.coordSys = 'sp' ;
tubi.coordSystemDemo(options)