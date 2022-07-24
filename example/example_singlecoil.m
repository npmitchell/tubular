%% Example pipeline for surface parameterization (including centerline extraction)
% on a single timepoint.

% First navigate to the location where this script resides:
origpath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(origpath))
cd ../


%% Add paths
dataDir = '/path/to/single_coil/' ;
cd(origpath)
addpath(origpath) ;
addpath(genpath('utility'))
addpath(genpath('external'))
addpath('TexturePatch')
addpath('DECLab')
addpath('RicciFlow_MATLAB')
% go back to the data


%%
outdir = dataDir ;
if ~exist(outdir, 'dir')
    mkdir(outdir) ;
end
cd(dataDir)

%% Define metadata
if ~exist(fullfile(dataDir, 'xp.mat'), 'file')
    %% DEFINE NEW MASTER SETTINGS
    if ~exist('./masterSettings.mat', 'file')
        % Metadata about the experiment
        stackResolution = [1, 1, 1] ;  % resolution in spaceUnits per pixel
        nChannels = 1 ;             % how many channels is the data (ex 2 for GFP + RFP)
        channelsUsed = 1 ;          % which channels are used for analysis
        timePoints = 0;       % timepoints to include in the analysis
        ssfactor = 4 ;              % subsampling factor
        flipy = false ;             % whether the data is stored inverted relative to real position in lab frame
        timeInterval = 1 ;          % physical interval between timepoints
        timeUnits = 'min' ;         % physical unit of time between timepoints
        spaceUnits = [char(956) 'm'] ;     % physical unit of time between timepoints
        fn = 'Time_%06d';        % filename string pattern
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
    %% PART 1: Define the metadata for the project
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd(dir16bit)
    dataDir = cd ;
    projectDir = dataDir ;

    % A filename base template - to be used throughout this script
    fileMeta                    = struct();
    fileMeta.dataDir            = dataDir;
    fileMeta.filenameFormat     = [fn, '.tif'];
    fileMeta.nChannels          = nChannels;
    fileMeta.timePoints         = timePoints ;
    fileMeta.stackResolution    = stackResolution;
    fileMeta.swapZT             = masterSettings.swapZT;

    % first_tp is also required, which sets the tp to do individually.
    first_tp = 1 ;
    expMeta                     = struct();
    expMeta.channelsUsed        = channelsUsed ;
    expMeta.channelColor        = 1;
    expMeta.description         = 'example tube';
    expMeta.dynamicSurface      = 1;
    expMeta.jitterCorrection    = 0;  % 1: Correct for sample translation
    expMeta.fitTime             = fileMeta.timePoints(first_tp);

    %% SET DETECTION OPTIONS ==================================================
    % Load/define the surface detection parameters
    detOpts_fn = fullfile(projectDir, 'detectOpts.mat') ;
    if exist(detOpts_fn, 'file')
        disp('loading detectOptions')
        load(detOpts_fn, 'detectOptions')
    else
        outputfilename_ply='mesh_ls_' ;
        outputfilename_ls='ls_' ;
        outputfilename_smoothply = 'mesh_' ;
        init_ls_fn = 'ls_initguess' ;
        prob_searchstr = '_stab_Probabilities.h5' ;
        preilastikaxisorder = set_preilastikaxisorder; ... % axis order in input to ilastik as h5s. To keep as saved coords use xyzc
        ilastikaxisorder= 'cxyz'; ... % axis order as output by ilastik probabilities h5
        imsaneaxisorder = 'xyzc'; ... % axis order relative to mesh axis order by which to process the point cloud prediction. To keep as mesh coords, use xyzc

        % Name the output mesh directory --------------------------------------
        meshDir = [fullfile(projectDir, 'mesh_output') filesep];

        % Surface detection parameters ----------------------------------------
        detectOptions = struct( 'channel', 1, ...
            'ssfactor', 1, ...
            'niter', 35,...
            'niter0', 160, ...
            'pre_pressure', -5, ...
            'pre_tension', 0, ...
            'pressure', 0, ...
            'tension', 0.5, ...
            'post_pressure', 2, ...
            'post_tension', 3, ...
            'exit_thres', 1e-7, ...
            'foreGroundChannel', 1, ...
            'fileName', sprintf( fn, 0 ), ...
            'meshDir', meshDir, ...
            'ofn_ls', outputfilename_ls, ...
            'ofn_ply', outputfilename_ply,...
            'timepoint', 0, ...
            'zdim', 2, ...
            'ofn_smoothply', outputfilename_smoothply, ...
            'init_ls_fn', init_ls_fn, ... % set to none to load prev tp
            'run_full_dataset', projectDir,... % projectDir, ... % set to 'none' for single tp
            'radius_guess', 40, ...
            'dset_name', 'exported_data',...
            'center_guess', '200,75,75',... % xyz of the initial guess sphere ;
            'save', true, ... % whether to save images of debugging output
            'plot_mesh3d', false, ...
            'dtype', 'h5',...
            'mask', 'none',...
            'mesh_from_pointcloud', false, ...
            'prob_searchstr', prob_searchstr, ...
            'preilastikaxisorder', preilastikaxisorder, ... 
            'ilastikaxisorder', ilastikaxisorder, ... 
            'physicalaxisorder', imsaneaxisorder, ... 
            'include_boundary_faces', true, ...
            'smooth_with_matlab', 0.01) ;

        % save options
        if exist(detOpts_fn, 'file')
            disp('Overwriting detectOptions --> renaming existing as backup')
            backupfn1 = [detOpts_fn '_backup1'] ;
            if exist(backupfn1, 'file')
                backupfn2 = [detOpts_fn '_backup2'] ; 
                system(['mv ' backupfn1 ' ' backupfn2])
            end
            system(['mv ' detOpts_fn ' ' backupfn1])
        end
        disp('Saving detect Options to disk')
        save(detOpts_fn, 'detectOptions') ;
    end

    % Overwrite certain parameters for script structure
    meshDir = detectOptions.meshDir ;

    %% Define Experiment as struct
    xp = struct('fileMeta', fileMeta, ...
        'expMeta', expMeta, 'detectOptions', detectOptions) ;
    disp('done')
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 2: TubULAR -- surface parameterization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now we have 3d data volumes and surfaces. Define a TubULAR object. 
    % To visualize data on these surfaces and compute how these surfaces deform
    % we now define TubULAR object.
    nU = masterSettings.nU ;
    nV = masterSettings.nV ;
    opts = struct() ;
    opts.meshDir = meshDir ;        % Directory where meshes reside
    opts.flipy = flipy ;            % Set to true if data volume axes are inverted in chirality wrt physical lab coordinates
    opts.timeInterval = timeInterval ; % Spacing between adjacent timepoints in units of timeUnits 
    opts.timeUnits = timeUnits ;    % units of time, so that adjacent timepoints are timeUnits * timeInterval apart
    opts.spaceUnits = spaceUnits ;  % Units of space in LaTeX, for ex '$mu$m' for micron
    opts.nU = nU ;                  % How many points along the longitudinal axis to sample surface
    opts.nV = nV ;                  % How many points along the circumferential axis to sample surface
    opts.normalShift = 10 ;         % Additional dilation acting on surface for texture mapping
    opts.a_fixed = 1.0 ;            % Fixed aspect ratio of pullback images. Setting to 1.0 is most conformal mapping option.
    opts.adjustlow = 1.00 ;         % floor for intensity adjustment
    opts.adjusthigh = 99.9 ;        % ceil for intensity adjustment (clip)
    opts.phiMethod = 'curves3d' ;   % Method for following surface in surface-Lagrangian mapping [(s,phi) coordinates]
    opts.lambda_mesh = 0.00 ;       % Smoothing applied to the mesh before DEC measurements
    opts.lambda = 0.0 ;             % Smoothing applied to computed values on the surface
    opts.lambda_err = 0.0 ;         % Additional smoothing parameter, optional
    opts.zwidth = 1 ;
    opts.nmodes = 7 ;
    % opts.t0 = xp.fileMeta.timePoints(1) ;   % reference timepoint used to define surface-Lagrangian and Lagrangian measurements
    % opts.t0 = 123 ;
    % opts.t0 = 37 ;
    % opts.t0 = 1 ;

    disp('saving xp struct and opts to disk')
    save(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')
else
    disp('loading xp struct from disk')
    load(fullfile(dataDir, 'xp.mat'), 'xp', 'opts')
end

%% TubULAR class instance
disp('defining TubULAR class instance (tubi= tubular instance)')
tubi = TubULAR(xp, opts) ;
disp('done defining TubULAR instance')


%% Define global orientation frame (for viewing in canonical frame)
% Compute APDV coordinate system
alignAPDVOpts = struct() ;
alignAPDVOpts.overwrite = false ;
tubi.computeAPDVCoords(alignAPDVOpts) ;

%% Select the endcaps for the centerline computation (A and P) and a point
% along which we will form a branch cut for mapping to the plane (D).
apdvOpts = struct() ;
apdvOpts.overwrite = false ;
apdvOpts.autoAP = true ;
[apts_sm, ppts_sm] = tubi.computeAPDpoints(apdvOpts) ;

%% Align the meshes in the APDV global frame & plot them
alignAPDVOpts.overwrite = true ;
alignAPDVOpts.forceEndpointsInside = true ;
alignAPDVOpts.normal_step = 2 ;
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
cntrlineOpts.res = 1 ;                 % resolution of distance tranform grid in which to compute centerlines
cntrlineOpts.preview = false ;           % preview intermediate results
cntrlineOpts.reorient_faces = false ;    % not needed for our well-constructed meshes
cntrlineOpts.dilation = 0 ;              % how many voxels to dilate the segmentation inside/outside before path computation

% Note: this can take about 400s per timepoint for res=2.0, so use as big a 
%   res value as possible.
%
tubi.generateFastMarchingCenterlines(cntrlineOpts)
disp('done with centerlines')

%% Identify anomalies in centerline data
idOptions.ssr_thres = 15 ;  % distance of sum squared residuals in um as threshold for removing spurious centerlines
tubi.cleanFastMarchingCenterlines(idOptions) ;
disp('done with cleaning up centerlines')

%% Cylinder cut mesh --> transforms a topological sphere into a topological cylinder
% Look for options on disk. If not saved, define options.
if ~exist(tubi.fileName.endcapOptions, 'file') 
    endcapOpts = struct( 'adist_thres', 1.5, ...  % 20, distance threshold for cutting off anterior in pix
                'pdist_thres', 1.5, ...  % 15-20, distance threshold for cutting off posterior in pix
                'tref', tubi.t0) ;  % reference timepoint at which time dorsal-most endcap vertices are defined
    tubi.setEndcapOptions(endcapOpts) ;
    % Save the options to disk
    disp('saving endcap options to disk')
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORBIFOLD -> begin populating tubi.dir.mesh/gridCoords_nUXXXX_nVXXXX/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overwrite = false ;
% Iterate Through Time Points to Create Pullbacks ========================
for tt = tubi.xp.fileMeta.timePoints
    disp(['NOW PROCESSING TIME POINT ', num2str(tt)]);
    tidx = tubi.xp.tIdx(tt);
    
    % Load the data for the current time point ------------------------
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
        compute_pullback = false ;
    end
    
    uvcutMesh = tubi.getCurrentUVCutMesh() ;
    
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

%% Plot new centerline
cylmesh = tubi.getCurrentCylinderMeshClean() ;
mesh = tubi.getCurrentCutMesh() ;
[ss, ~, cline] = tubi.getCurrentClineDVhoop() ;
cline = tubi.APDV2xyz(cline) ;
fmcline = tubi.getCleanFastMarchingCenterlines() ;
fmcline = fmcline{1}(:, 2:4) ;
mesh_sma = mesh_sm ;
mesh_sma.v = mesh_sm.v + 10 * mesh_sm.vn ;

%% Plot result
falph = 0.2 ;
[colors, ~] = define_colors ;
sky = colors(6, :) ;
disp('plotting resulting mesh')
clf
trisurf(triangulation(mesh_sm.f, mesh_sma.v(:, [3,1,2])), 'edgecolor', 'none',...
    'facecolor', sky, 'facealpha', falph)
axis equal
hold on;
plot3(sz, sx, sy, 'k-', 'linewidth', 2)
lighting gouraud
material dull
camlight 
axis off
set(gcf, 'color', 'white')
export_fig(gcf, fullfile(outdir, 'mesh_and_skeleton_bwskel_10step.png'), '-nocrop', '-r600')

% Plot result -- fmarching
[colors, names] = define_colors ;
sky = colors(6, :) ;
disp('plotting resulting mesh')

clf
trisurf(triangulation(mesh.f, mesh.v(:, [3,1,2])), 'edgecolor', 'none',...
    'facecolor', sky, 'facealpha', falph)
axis equal
hold on;
% plot3(sz, sx, sy, '.-')
plot3(fmcline(:, 3), fmcline(:, 1), fmcline(:, 2), 'k-', 'linewidth', 2)
lighting gouraud
material dull
camlight 
axis off
set(gcf, 'color', 'white')
export_fig(gcf, fullfile(outdir, 'mesh_and_skeleton_fastMarching.png'), '-nocrop', '-r600')

% plot result -- DV hoops centerline
clf
trisurf(triangulation(mesh.f, mesh.v(:, [3,1,2])), 'edgecolor', 'none',...
    'facecolor', sky, 'facealpha', falph)
axis equal
hold on;
% plot3(sz, sx, sy, '.-')
plot3(cline(:, 3), cline(:, 1), cline(:, 2), 'k-', 'linewidth', 2)
lighting gouraud
material dull
camlight 
axis off
set(gcf, 'color', 'white')
export_fig(gcf, fullfile(outdir, 'mesh_and_skeleton_avgptsDVhoop.png'), '-nocrop', '-r600')

% Show associations with mesh vertices
clf
spcutMesh = tubi.getCurrentSPCutMesh() ;
trisurf(triangulation(spcutMesh.f, spcutMesh.v(:, [3,1,2])), 'edgecolor', 'none',...
    'facevertexcdata', spcutMesh.sphi(:, 1), 'facealpha', falph)
axis equal
hold on;
% plot3(sz, sx, sy, '.-')
cmap = viridis(length(cline)) ;
for ii = 1:length(cline)
    color = cmap(ii, :) ;
    scatter3(cline(ii, 3), cline(ii, 1), cline(ii, 2), 10, color, 'filled', ...
        'markeredgecolor', 'none')
end
lighting gouraud
material dull
camlight 
axis off
set(gcf, 'color', 'white')
export_fig(gcf, fullfile(outdir, 'mesh_and_skeleton_avgptsDVhoop_associations.png'), '-nocrop', '-r600')

% associations v2

clf
spcutMesh = tubi.getCurrentSPCutMesh() ;
trisurf(triangulation(spcutMesh.f, spcutMesh.v(:, [3,1,2])), 'edgecolor', 'none',...
     'facevertexcdata', spcutMesh.sphi(:, 1), 'facealpha', 0.0)
 
axis equal
hold on;
% plot3(sz, sx, sy, '.-')
cmap = viridis(length(cline)) ;
for ii = fliplr(1:3:length(cline))
    color = cmap(ii, :) ;
    for jj = 1:5:nV
        xyz = spcutMesh.v(ii+(jj-1)*nV, :) ;
        xx = xyz(1); yy = xyz(2); zz = xyz(3) ;
        plot3([cline(ii, 3), zz], [cline(ii, 1), xx], [cline(ii, 2), yy], ...
            '-', 'color', color)
    end
end
lighting gouraud
material dull
camlight 
axis off
set(gcf, 'color', 'white')
export_fig(gcf, fullfile(outdir, 'mesh_and_skeleton_avgptsDVhoop_associations2.png'), '-nocrop', '-r600')


%% PLOT ALL TEXTURED MESHES IN 3D (OPTIONAL: this is SLOW) ================
% Establish texture patch options
metadat = struct() ;
metadat.reorient_faces = false ;            % set to true if some mesh normals may be inverted (requires gptoolbox if true)
metadat.normal_shift = tubi.normalShift ;   % normal push, in pixels, along normals defined in data XYZ space
metadat.texture_axis_order = [1 2 3] ;      % texture space sampling. If the surface and dataspace have axis permutation, enter that here
Options.PSize = 5 ;          % Psize is the linear dimension of the grid drawn on each triangular face. Set PSize > 1 for refinement of texture on each triangle of the surface triangulation. Higher numbers are slower but give more detailed images.
Options.numLayers = [0, 0];  % how many layers to MIP over/bundle into stack, as [outward, inward]
Options.layerSpacing = 2 ;   % Distance between layers over which we take MIP, in pixels, 
metadat.blackFigure = false ;

% Plot on surface for all timepoints 
tubi.plotSeriesOnSurfaceTexturePatch(metadat, Options)