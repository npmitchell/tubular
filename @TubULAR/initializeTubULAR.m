function initializeTubULAR(tubi, xp, opts)
%initializeTubULAR(tubi, xp, opts)
%   Hidden method for instantiating TubULAR class
%
% Parameters
% ----------
% tubi : TubULAR object whose properties to fill in with this function
% xp : Imsane Experiment class instance belonging to tubi or struct
%    if imsane class instance, tubular uses imsane fields to populate
%    metadata. Otherwise this needs to be a struct with fields
%      fileMeta
%      expMeta
%    See docs for entries in these fields and their descriptions.
% opts : struct with fields
%   flipy : bool
%       Set to true if lab coordinates are mirrored along some axis wrt to 
%       data coords.
%   useBioformats : bool
%       Set to true if data should be loaded using the Bioformats importer,
%       if possible
%   meshDir : str
%       path to where meshes are stored and output will be placed
%   timeUnits : str
%       units of time, for ex 'min'
%   spaceUnits : str
%       units of space, for ex '$mu$m' for microns
%   nU : int
%       resolution in the longitudinal direction, in number of sampling
%       points per proper length of the mesh
%   nV : int
%       resolution in the circumferential direction, in number of sampling
%       points per circumference
%   lambda : optional float 
%       smoothing applied to fields on a mesh before computation of DEC
%       derivatives
%   lambda_mesh : optional float 
%       smoothing (diffusion constant) applied to mesh vertices before
%       computation of DEC fields
%   lambda_err : optional float
%       additional smoothing for fields derived from DEC fields
%   
%
% NPMitchell 2020-2022

%% PROPERTIES

% Build timepoint index lookup table if not already baked into xp as an
% ImSAnE experiment class instance:
if isa(xp, 'struct')
    if ~isfield(xp, 'tIdx')
        xp.tIdx = containers.Map(xp.fileMeta.timePoints, ...
            1:length(xp.fileMeta.timePoints)) ;
    end
end

makeDirs = true;
tubi.xp = xp ;
tubi.flipy = opts.flipy ;
meshDir = opts.meshDir ;
tubi.timeUnits = opts.timeUnits ;
tubi.spaceUnits = opts.spaceUnits ;
tubi.fileBase.fn = xp.fileMeta.filenameFormat ;
tubi.ssfactor = xp.detectOptions(1).ssfactor ;
tubi.nU = opts.nU ;
tubi.nV = opts.nV ;
if isfield(opts, 'timeInterval')
    tubi.timeInterval = opts.timeInterval ;
elseif isfield(opts, 'timeinterval')
    tubi.timeInterval = opts.timeinterval ;
else
    tubi.timeInterval = 1 ;
end

% Define reference time t0
tubi.fileName.t0 = fullfile(meshDir, 't0_referenceTime.txt') ;
if isfield(opts, 't0')
    tubi.t0 = opts.t0 ;
    t0supplied_in_opts = true ;
else
    t0supplied_in_opts = false ;
    if exist(tubi.fileName.t0, 'file')
        tubi.t0 = dlmread(tubi.fileName.t0, ',', 1, 0) ;
    else
        tubi.t0 = tubi.xp.fileMeta.timePoints(1) ;
    end
end
dynamic = length(tubi.xp.fileMeta.timePoints) > 1 ;
tubi.dynamic = dynamic ; 

if isfield(opts, 'normalShift')
    tubi.normalShift = opts.normalShift ;
end
if isfield(opts, 'a_fixed')
    tubi.a_fixed = opts.a_fixed ;
end
if isfield(opts, 'imSize')
    tubi.imSize = opts.imSize ;
else
    tubi.imSize = [tubi.a_fixed * 1000, 1000] ;
end
if isfield(opts, 'adjustlow')
    tubi.data.adjustlow = opts.adjustlow ;
end
if isfield(opts, 'adjusthigh')
    tubi.data.adjusthigh = opts.adjusthigh ;
end
if isfield(opts, 'axisOrder')
    tubi.data.axisOrder = opts.axisOrder ;
end
if isfield(opts, 'ilastikOutputAxisOrder')
    tubi.data.ilastikOutputAxisOrder = opts.ilastikOutputAxisOrder ;
end
if isfield(opts, 'makeDirs')
    makeDirs = opts.makeDirs;
end
if isfield(opts, 'useBioformats')
    tubi.useBioformats = opts.useBioformats;
else
    tubi.useBioformats = true;
end
    
% Assign which pullback coordsys is used for velocimetry
if isfield(opts, 'pivPullback')
    tubi.pivPullback = pivPullback ;
end

%% NAMING
uvexten = sprintf('_nU%04d_nV%04d', tubi.nU, tubi.nV) ;

% Naming extension
tubi.uvexten = uvexten ;

% APDV coordinate system
tubi.APDV.resolution = min(xp.fileMeta.stackResolution) ;
tubi.APDV.rot = [] ;
tubi.APDV.trans = [] ;

% Plotting settings
tubi.plotting.colors = define_colors() ;
tubi.plotting.markers = {'o', 's', '^', 'v', '*', '>', '<'} ;

% Directories of the tubi object
% Meshes and measurements before gridding into pullback coords
tubi.dir.data = xp.fileMeta.dataDir ;
tubi.fileName.apdvOptions = fullfile(tubi.dir.data, 'alignAPDV_Options.mat') ;
tubi.dir.mesh = meshDir ;
tubi.dir.alignedMesh = fullfile(meshDir, 'alignedMesh') ;
tubi.dir.cntrline = fullfile(meshDir, 'centerline') ;
tubi.dir.cylinderMesh = fullfile(meshDir, 'cylinderMesh') ;
tubi.dir.cutMesh = fullfile(meshDir, 'cutMesh') ;
tubi.dir.rawRicciMesh = fullfile(meshDir, 'rawRicciMesh') ;
tubi.dir.cylinderMeshClean = fullfile(tubi.dir.cylinderMesh, 'cleaned') ;
tubi.dir.texturePatchIm = fullfile(meshDir, 'images_texturepatch') ;
tubi.dir.mip = fullfile(meshDir, 'mips', 'dim%d_pages%04dto%04d') ;
% Mips
tubi.fileBase.mip = ['mip_' tubi.timeStampStringSpec '.tif'] ;
tubi.fullFileBase.mip = fullfile(tubi.dir.mip, ['mip_' tubi.timeStampStringSpec '.tif']) ;

% After gridding into (u,v) / (zeta,phi) pullback coords
uvDir = fullfile(tubi.dir.mesh, sprintf('gridCoords_nU%04d_nV%04d', tubi.nU, tubi.nV)) ;
tubi.dir.uvCoord = uvDir ;

%% Since Windows machines have issues with sprintf due to the presence of 
% backslashes in filenames. We therefore use num2str() with the timestamp
% format specifier. 
timeStampFormatSpec = '%\d+[a-z]';
tmp = regexp(tubi.fileBase.fn, timeStampFormatSpec, 'match');
tubi.timeStampStringSpec = tmp{1};

%% fileBases
tubi.fileBase.name = xp.fileMeta.filenameFormat(1:end-4) ;
tubi.fullFileBase.name = fullfile(tubi.dir.data, tubi.fileBase.name) ;
if ~strcmp(tubi.fullFileBase.name(end-4:end), '.tif') || ...
    ~strcmp(tubi.fullFileBase.name(end-5:end), '.tiff') 
    tubi.fullFileBase.name = [tubi.fullFileBase.name '.tif'] ;
end
% Note: default is of the form: tubi.fileBase.mesh = 'mesh_%06d.ply' ;
try
    tubi.fileBase.mesh = xp.detector.options.ofn_smoothply ;
    if strcmpi(tubi.fileBase.mesh(end-3:end), '.ply')
        tubi.fileBase.mesh = tubi.fileBase.mesh(1:end-4) ;
    end
    if length(tubi.xp.fileMeta.timePoints) > 1
        if ~contains(tubi.fileBase.mesh, '%')
            tubi.fileBase.mesh = [tubi.fileBase.mesh tubi.timeStampStringSpec] ;
        end
    end
catch
    if contains(xp.detectOptions.ofn_smoothply, '%') && ...
             contains(xp.detectOptions.ofn_smoothply, 'd')
        if contains(xp.detectOptions.ofn_smoothply, '.ply')
            tubi.fileBase.mesh = xp.detectOptions.ofn_smoothply(1:end-4) ;
        else
            tubi.fileBase.mesh = xp.detectOptions.ofn_smoothply ;
        end
    else
        if contains(xp.detectOptions.ofn_smoothply, '.ply')
            tubi.fileBase.mesh = [xp.detectOptions.ofn_smoothply(1:end-4) tubi.timeStampStringSpec] ;
        else
            tubi.fileBase.mesh = [xp.detectOptions.ofn_smoothply tubi.timeStampStringSpec] ;
        end
            
    end
     
end
tubi.fileBase.alignedMesh = ...
    [tubi.fileBase.mesh '_APDV_um'] ;
tubi.fileBase.apdProb = [tubi.fileBase.name '_Probabilities_APDVcoords.h5']  ;
tubi.fileBase.apCenterlineProb = [tubi.fileBase.name '_Probabilities_apcenterline.h5']  ;
tubi.fileBase.prob = [tubi.fileBase.name '_Probabilities.h5'] ; 
tubi.fileBase.centerlineXYZ = ...
    [tubi.fileBase.mesh '_centerline_exp1p0_res*.txt' ] ;
tubi.fileBase.centerlineAPDV = ...
    [tubi.fileBase.mesh '_centerline_scaled_exp1p0_res*.txt' ] ;
tubi.fileBase.cylinderMesh = ...
    [tubi.fileBase.mesh '_cylindercut.ply'] ;
tubi.fileBase.apBoundary = ['ap_boundary_indices_' tubi.timeStampStringSpec '.mat'];
tubi.fileBase.cylinderKeep = ['cylinderMesh_keep_indx_' tubi.timeStampStringSpec '.mat'] ;
tubi.fileName.apBoundaryDorsalPts = 'ap_boundary_dorsalpts.h5' ;

%% Metric
tubi.dir.metric = struct() ;
tubi.dir.metric.data = fullfile(tubi.dir.uvCoord, 'metric', '%s_lmesh%0.3f') ;
tubi.fullFileBase.metric = fullfile(tubi.dir.metric.data, [tubi.fileBase.name '.mat']) ;
tubi.dir.metric.g_images2d = fullfile(tubi.dir.metric.data, 'g_images2d') ;
tubi.dir.metric.b_images2d = fullfile(tubi.dir.metric.data, 'b_images2d') ;
tubi.dir.metric.g_images3d = fullfile(tubi.dir.metric.data, 'g_images3d') ;
tubi.dir.metric.b_images3d = fullfile(tubi.dir.metric.data, 'b_images3d') ;

%% define string for smoothing params
if isfield(opts, 'lambda')
    tubi.smoothing.lambda = opts.lambda ;
end
if isfield(opts, 'lambda_mesh')
    tubi.smoothing.lambda_mesh = opts.lambda_mesh ;
end
if isfield(opts, 'lambda_err')
    tubi.smoothing.lambda_err = opts.lambda_err ;
end
if isfield(opts, 'nmodes')
    tubi.smoothing.nmodes = opts.nmodes ;
end
if isfield(opts, 'zwidth')
    tubi.smoothing.zwidth = opts.zwidth ;
end
l_lmesh_lerr =  strrep(sprintf('lambda%0.03f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d', ...
    tubi.smoothing.lambda, tubi.smoothing.lambda_mesh, ...
    tubi.smoothing.lambda_err, tubi.smoothing.nmodes, tubi.smoothing.zwidth), '.', 'p') ;
l_lmesh = strrep(sprintf('lambda%0.03f_lmesh%0.3f_modes%02dw%02d', ...
    tubi.smoothing.lambda, tubi.smoothing.lambda_mesh, ...
    tubi.smoothing.nmodes, tubi.smoothing.zwidth), '.', 'p') ;

%% Dynamics -- strain rates, kinematics
if dynamic
    % Metric strain dirs
    tubi.dir.metricKinematics = struct() ;
    tubi.dir.metricKinematics.root = fullfile(uvDir, 'metricKinematics') ;
    tubi.dir.metricKinematics.smoothing = fullfile(uvDir, 'metricKinematics', ...
        l_lmesh_lerr) ;
    tubi.dir.metricKinematics.measurements = ...
        fullfile(uvDir, 'metricKinematics', l_lmesh_lerr, 'measurements') ;
    tubi.fullFileBase.metricKinematics = struct() ;
    tubi.fullFileBase.metricKinematics.divv = ...
        fullfile(tubi.dir.metricKinematics.measurements, ['divv_vertices_' ...
        tubi.timeStampStringSpec]) ;
    tubi.fullFileBase.metricKinematics.H2vn = ...
        fullfile(tubi.dir.metricKinematics.measurements, ['H2vn_vertices_' ...
        tubi.timeStampStringSpec]) ;
    tubi.fullFileBase.metricKinematics.gdot = ...
        fullfile(tubi.dir.metricKinematics.measurements, ['gdot_vertices_' ...
        tubi.timeStampStringSpec]) ;

    % Metric Kinematics along pathlines
    tubi.dir.metricKinematics.pathline = struct() ;
    tubi.dir.metricKinematics.pathline.root = ...
        fullfile(uvDir, 'metricKinematics', l_lmesh_lerr, 'pathline_%04dt0') ;
    tubi.dir.metricKinematics.pathline.measurements = ...
        fullfile(tubi.dir.metricKinematics.pathline.root, 'measurements') ;
    tubi.fileBase.metricKinematics = struct() ;
    tubi.fileBase.metricKinematics.pathline.Hfn = ...
        ['HH_pathline%04d_' tubi.timeStampStringSpec '.mat'] ;
    tubi.fileBase.metricKinematics.pathline.gdot = ...
        ['gdot_pathline%04d_' tubi.timeStampStringSpec '.mat'] ;
    tubi.fileBase.metricKinematics.pathline.divv = ...
        ['divv_pathline%04d_' tubi.timeStampStringSpec '.mat'] ;
    tubi.fileBase.metricKinematics.pathline.veln = ...
        ['veln_pathline%04d_' tubi.timeStampStringSpec '.mat'] ;
    tubi.fileBase.metricKinematics.pathline.radius = ...
        ['radius_pathline%04d_' tubi.timeStampStringSpec '.mat'] ;
    tubi.fileBase.metricKinematics.pathline.kymographs = struct() ;
    tubi.fileBase.metricKinematics.pathline.kymographs.ap = 'apKymographMetricKinematics.mat' ;
    tubi.fileBase.metricKinematics.pathline.kymographs.d = 'dKymographMetricKinematics.mat' ;
    tubi.fileBase.metricKinematics.pathline.kymographs.v = 'vKymographMetricKinematics.mat' ;
    tubi.fileBase.metricKinematics.pathline.kymographs.l = 'lKymographMetricKinematics.mat' ;
    tubi.fileBase.metricKinematics.pathline.kymographs.r = 'rKymographMetricKinematics.mat' ;

    % Metric deformation ('metric strain' or 'meshMetric)
    tubi.dir.gstrain = fullfile(uvDir, 'metricStrain') ;
    tubi.dir.gstrainRate = fullfile(tubi.dir.gstrain, 'rateMetric') ;
    tubi.dir.gstrainRateIm = fullfile(tubi.dir.gstrainRate, 'images') ;
    tubi.fileBase.gstrainRate = ['gstrainRate_' tubi.timeStampStringSpec '.mat'] ;
    tubi.fullFileBase.gstrainRate = fullfile(tubi.dir.gstrainRate, ...
        tubi.fileBase.gstrainRate) ; 

    %% Strain Rate
    tubi.dir.strainRate = struct() ;
    tubi.dir.strainRate.root = fullfile(uvDir, 'strainRate') ;
    tubi.dir.strainRate.smoothing = fullfile(uvDir, 'strainRate', l_lmesh) ;
    tubi.dir.strainRate.measurements = ...
        fullfile(uvDir, 'strainRate', l_lmesh, 'measurements') ;
    tubi.fileBase.strainRate = ['strainRate_' tubi.timeStampStringSpec '.mat'] ;
    tubi.fullFileBase.strainRate = fullfile(tubi.dir.strainRate.measurements, ...
        tubi.fileBase.strainRate) ; 
    %% Pathline-based strain measurement --> from pathline path
    % see pathline section

    %% Strain rates along pathlines -- to be filled in with smoothing and t0
    tubi.dir.strainRate.pathline = struct() ;
    tubi.dir.strainRate.pathline.root = ...
        fullfile(uvDir, 'strainRate', l_lmesh, 'pathline_%04dt0') ;
    tubi.dir.strainRate.pathline.measurements = ...
        fullfile(tubi.dir.strainRate.pathline.root, 'measurements') ;
        
    %% Forces from stokes force balance    
    % Metric strain dirs
    tubi.dir.stokesForces = struct() ;
    tubi.dir.stokesForces.root = fullfile(uvDir, 'stokesForces') ;
    tubi.dir.stokesForces.smoothing = fullfile(uvDir, 'stokesForces', ...
        l_lmesh_lerr) ;
    tubi.dir.stokesForces.measurements = ...
        fullfile(uvDir, 'stokesForces', l_lmesh_lerr, 'measurements') ;
    tubi.fullFileBase.stokesForces = struct() ;
    tubi.fullFileBase.stokesForces.Lapv = ...
        fullfile(tubi.dir.metricKinematics.measurements, ['Lapv_vertices_' tubi.timeStampStringSpec]) ;
    tubi.fullFileBase.stokesForces.Kv = ...
        fullfile(tubi.dir.metricKinematics.measurements, ['Kv_vertices_' tubi.timeStampStringSpec]) ;
    tubi.fullFileBase.stokesForces.gradP = ...
        fullfile(tubi.dir.metricKinematics.measurements, ['gradP_vertices_' tubi.timeStampStringSpec]) ;

end

% shorten variable names for brevity
clineDir = tubi.dir.cntrline ;


%% Clean Cylinder Mesh
tubi.fileName.aBoundaryDorsalPtsClean = ...
    fullfile(tubi.dir.cylinderMeshClean, 'adIDx.h5') ;
tubi.fileName.pBoundaryDorsalPtsClean = ...
    fullfile(tubi.dir.cylinderMeshClean, 'pdIDx.h5') ;

%% cutMesh
tubi.fullFileBase.cutPath = fullfile(tubi.dir.cutMesh, ...
    ['cutPaths_' tubi.timeStampStringSpec '.txt']) ;

%% fileNames
nshift = strrep(sprintf('%03d', tubi.normalShift), '-', 'n') ;
shiftstr = ['_' nshift 'step'] ;
tubi.fileName.rot = fullfile(meshDir, 'rotation_APDV.txt') ;
tubi.fileName.trans = fullfile(meshDir, 'translation_APDV.txt') ;
tubi.fileName.xyzlim_raw = fullfile(meshDir, 'xyzlim_raw.txt') ;
tubi.fileName.xyzlim_pix = fullfile(meshDir, 'xyzlim_APDV.txt') ;
tubi.fileName.xyzlim_um = ...
    fullfile(meshDir, 'xyzlim_APDV_um.txt') ;
tubi.fileName.xyzlim_um_buff = ...
    fullfile(meshDir, ['xyzlim_APDV_um' shiftstr '.txt']) ;
% fileNames for APDV and cylinderMesh
tubi.fullFileBase.apdProb = fullfile(tubi.dir.data, tubi.fileBase.apdProb) ;
tubi.fullFileBase.apCenterlineProb = fullfile(tubi.dir.data, tubi.fileBase.apCenterlineProb) ;
tubi.fullFileBase.prob = fullfile(tubi.dir.data, tubi.fileBase.prob) ;
tubi.fileName.apdv = ...
    fullfile(clineDir, 'apdv_pts_for_centerline.h5') ;
tubi.fileName.dpt = fullfile(meshDir, 'dpt_for_rot.txt') ;
tubi.fileName.startendPt = fullfile(clineDir, 'startendpt.h5') ;
tubi.fileName.cleanFMCenterlines = ...
    fullfile(clineDir, 'centerlines_anomalies_fixed.mat') ;
tubi.fileName.apBoundaryDorsalPts = ...
    fullfile(tubi.dir.cylinderMesh, 'ap_boundary_dorsalpts.h5') ;
tubi.fileName.endcapOptions = ...
    fullfile(tubi.dir.cylinderMesh, 'endcapOptions.mat') ;
tubi.fileName.apdBoundary = ...
    fullfile(tubi.dir.cylinderMesh, 'ap_boundary_dorsalpts.h5') ;

%% FileNamePatterns
tubi.fullFileBase.mesh = ...
    fullfile(tubi.dir.mesh, [tubi.fileBase.mesh '.ply']) ;
tubi.fullFileBase.alignedMesh = ...
    fullfile(tubi.dir.alignedMesh, [tubi.fileBase.alignedMesh '.ply']) ;
% fileNames for centerlines
tubi.fullFileBase.centerlineXYZ = ...
    fullfile(clineDir, tubi.fileBase.centerlineXYZ) ;
tubi.fullFileBase.centerlineAPDV = ...
    fullfile(clineDir, tubi.fileBase.centerlineAPDV) ;
tubi.fullFileBase.cylinderMesh = ...
    fullfile(tubi.dir.cylinderMesh, tubi.fileBase.cylinderMesh) ;
tubi.fullFileBase.apBoundary = ...
    fullfile(tubi.dir.cylinderMesh, tubi.fileBase.apBoundary) ;
% tubi.fullFileBase.apBoundaryDorsalPts = tubi.fileName.apBoundaryDorsalPts ;
tubi.fullFileBase.cylinderKeep = ...
    fullfile(tubi.dir.cylinderMesh, tubi.fileBase.cylinderKeep) ;
tubi.fileBase.cylinderMeshClean = [tubi.fileBase.mesh '_cylindercut_clean.ply'] ;
tubi.fullFileBase.cylinderMeshClean = ...
    fullfile(tubi.dir.cylinderMesh, 'cleaned',...
    tubi.fileBase.cylinderMeshClean) ;   

%% Define cutMesh directories
% cutMesh = fullfile(meshDir, 'cutMesh') ;
% cutMeshBase = fullfile(cutMesh, [tubi.fileBase.name, '_cutMesh.mat']) ;
imFolderBase = fullfile(uvDir, ['PullbackImages' shiftstr] ) ;
ricciMeshDir = fullfile(uvDir, ['ricci_cutMesh' shiftstr], 'noResampling') ;
ricciMeshDirWithResampling = fullfile(uvDir, ['ricci_cutMesh' shiftstr], ...
    'withResampling') ;
uvmeshDir = fullfile(uvDir, ['uv_cutMesh' shiftstr]) ;
sphiDir = fullfile(uvDir, ['sphi_cutMesh' shiftstr]) ;
if dynamic
    sphiSmDir = fullfile(sphiDir, 'smoothed') ;
    sphiSmRSDir = fullfile(sphiDir, 'smoothed_rs') ;
    % sphiSmRSImDir = fullfile(sphiSmRSDir, 'images') ;
    % sphiSmRSPhiImDir = fullfile(sphiSmRSImDir, 'phicolor') ;
    sphiSmRSCDir = fullfile(sphiDir, 'smoothed_rs_closed') ;
else
    sphiRSDir = fullfile(sphiDir, 'rs') ;
    sphiRSCDir = fullfile(sphiDir, 'rs_closed') ;
end

%% Images of pullbacks in mesh coordinate systems
imFolder_sp = [imFolderBase '_sphi'] ;
imFolder_spe = fullfile(imFolder_sp, 'extended') ;
imFolder_up = [imFolderBase '_uphi'] ;
imFolder_upe = fullfile(imFolder_up, 'extended') ;
imFolder_ricci = [ imFolderBase '_ricci' ] ;
imFolder_ricci_e = fullfile(imFolder_ricci, 'extended') ;
imFolder_re_stack = fullfile([imFolderBase, '_sphi_relaxed'], ...
    'extended_stack') ;  
imFolder_spe_stack = fullfile([imFolderBase, '_sphi'], ...
    'extended_stack') ;  
imFolder_pivPathline = [imFolderBase, '_pivPathlines_', tubi.timeStampStringSpec, 't0'] ;

% time-averaged meshes
if dynamic
    imFolder_spsm = fullfile(imFolder_sp, 'smoothed') ;
    imFolder_spsme = fullfile(imFolder_sp, 'smoothed_extended') ;  % raw LUT, no histeq
    imFolder_spsmeLUT = fullfile(imFolder_sp, 'smoothed_extended_LUT') ;  % with histeq
    imFolder_rsm = fullfile([imFolderBase, '_sphi_relaxed'], 'smoothed');
    imFolder_rsme = fullfile([imFolderBase, '_sphi_relaxed'], 'smoothed_extended') ;
    imFolder_rsme_stack = fullfile([imFolderBase, '_sphi_relaxed'], ...
        'smoothed_extended_stack') ;  % with histeq?
end

% Folder for curvature measurements
if dynamic
    KHSmDir = fullfile(sphiSmRSCDir, 'curvature') ;
    KSmDir = fullfile(KHSmDir, 'gauss') ;
    HSmDir = fullfile(KHSmDir, 'mean') ;
else
    KHSmDir = fullfile(sphiRSCDir, 'curvature') ;
    KSmDir = fullfile(KHSmDir, 'gauss') ;
    HSmDir = fullfile(KHSmDir, 'mean') ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Port into tubi
% (already above) tubi.dir.cutMesh = fullfile(meshDir, 'cutMesh') ;
tubi.dir.uvcurMesh = uvmeshDir ;
tubi.dir.spcutMesh = sphiDir ;
if dynamic
    tubi.dir.spcutMeshSm = sphiSmDir ;
    tubi.dir.spcutMeshSmRS = sphiSmRSDir ;
    tubi.dir.spcutMeshSmRSC = sphiSmRSCDir ;
else
    tubi.dir.spcutMeshRS = sphiRSDir ;
    tubi.dir.spcutMeshRSC = sphiRSCDir ;
end
tubi.dir.ricci = struct() ;
tubi.dir.ricci.data = ricciMeshDir ;
tubi.dir.ricci.mesh = fullfile(ricciMeshDir, 'meshes') ;
tubi.dir.ricci.solution = fullfile(ricciMeshDir, 'ricciSolutions') ;
tubi.dir.ricci.mu = fullfile(ricciMeshDir, 'beltramiCoefficients') ;

tubi.dir.ricci.dataWithResampling = ricciMeshDirWithResampling ;
tubi.dir.ricci.meshWithResampling = fullfile(ricciMeshDirWithResampling, 'meshes') ;
tubi.dir.ricci.solutionWithResampling = fullfile(ricciMeshDirWithResampling, 'ricciSolutions') ;
tubi.dir.ricci.muWithResampling = fullfile(ricciMeshDirWithResampling, 'beltramiCoefficients') ;

tubi.dir.rawRicci = struct() ;
rawRicciMeshDir = tubi.dir.rawRicciMesh ;
tubi.dir.rawRicci.data = tubi.dir.rawRicciMesh ;
tubi.dir.rawRicci.mesh = fullfile(rawRicciMeshDir, 'meshes') ;
tubi.dir.rawRicci.solution = fullfile(rawRicciMeshDir, 'ricciSolutions') ;
tubi.dir.rawRicci.mu = fullfile(rawRicciMeshDir, 'beltramiCoefficients') ;

%% Raw Ricci mesh coordinates (truly conformal via Ricci flow from cleaned cylinder mesh)
tubi.fileBase.rawRicciMesh = ['rawRicciMesh_%04diter_' tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.rawRicciMesh = fullfile(tubi.dir.rawRicci.mesh, tubi.fileBase.rawRicciMesh) ;
tubi.fileBase.rawRicciSolution = ['rawRicciSolution_%04diter_' tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.rawRicciSolution = fullfile(tubi.dir.rawRicci.solution, tubi.fileBase.rawRicciSolution) ;
tubi.fileBase.rawRicciMu = ['rawRicciMesh_mus_%04diter_' tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.rawRicciMu = fullfile(tubi.dir.rawRicci.mu, tubi.fileBase.rawRicciMu) ;

%% Ricci mesh coordinates (truly conformal via Ricci flow from reparameterized gridCoords mesh (spsm uv or sp))
tubi.fileBase.ricciMesh = ['ricciMesh_%04diter_' tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.ricciMesh = fullfile(tubi.dir.ricci.mesh, tubi.fileBase.ricciMesh) ;
tubi.fileBase.ricciSolution = ['ricciSolution_%04diter_' tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.ricciSolution = fullfile(tubi.dir.ricci.solution, tubi.fileBase.ricciSolution) ;
tubi.fileBase.ricciMu = ['ricciMesh_mus_%04diter_' tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.ricciMu = fullfile(tubi.dir.ricci.mu, tubi.fileBase.ricciMu) ;

tubi.fileBase.ricciMeshWithResampling = ['ricciMeshWithResampling_%04diter_' ...
    tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.ricciMeshWithResampling = fullfile(tubi.dir.ricci.meshWithResampling, tubi.fileBase.ricciMeshWithResampling) ;
tubi.fileBase.ricciSolutionWithResampling = ['ricciSolutionWithResampling_%04diter_' ...
    tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.ricciSolutionWithResampling = fullfile(tubi.dir.ricci.solutionWithResampling, tubi.fileBase.ricciSolutionWithResampling) ;
tubi.fileBase.ricciMuWithResampling = ['ricciMesh_musWithResampling_%04diter_' ...
    tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.ricciMuWithResampling = ...
    fullfile(tubi.dir.ricci.muWithResampling, tubi.fileBase.ricciMuWithResampling) ;

%% centerlines
tubi.dir.clineDVhoop = ...
    fullfile(tubi.dir.uvCoord, ...
    ['centerline_from_DVhoops' shiftstr]) ;
tubi.dir.writhe =  fullfile(tubi.dir.clineDVhoop, 'writhe') ;

%% Images
tubi.dir.im_uv = [imFolderBase '_uv'] ;
tubi.dir.im_uve = [imFolderBase '_uv_extended'] ;
tubi.dir.im_r = [imFolderBase '_sphi_relaxed'] ;
tubi.dir.im_re = fullfile(tubi.dir.im_r, 'extended') ;
tubi.dir.im_sp = imFolder_sp ;
tubi.dir.im_spe = imFolder_spe ;
tubi.dir.im_spe_stack = imFolder_spe_stack ;
tubi.dir.im_up = imFolder_up ;
tubi.dir.im_upe = imFolder_upe ;
if dynamic
    tubi.dir.im_sp_sm = imFolder_spsm ;
    tubi.dir.im_sp_sme = imFolder_spsme ;
    tubi.dir.im_sp_smeLUT = imFolder_spsmeLUT ;
    tubi.dir.im_r_sm = imFolder_rsm ;
    tubi.dir.im_r_sme = imFolder_rsme ;
    tubi.dir.im_r_sme_stack = imFolder_rsme_stack ;
    tubi.dir.im_pivPathlines = imFolder_pivPathline ;
end
tubi.dir.im_ricci = imFolder_ricci ;
tubi.dir.im_ricci_e = imFolder_ricci_e ;
tubi.dir.im_re_stack = imFolder_re_stack ;

tubi.dir.curvatures = KHSmDir ;
tubi.dir.meanCurvature = HSmDir ;
tubi.dir.gaussCurvature = KSmDir ;
tubi.fileBase.curvatures = ['guass_mean_curvature_' tubi.timeStampStringSpec '.mat'] ;
tubi.fullFileBase.curvatures = fullfile(tubi.dir.curvatures, tubi.fileBase.curvatures);
tubi.fullFileBase.cutMesh = ...
    fullfile(tubi.dir.cutMesh, [tubi.fileBase.name, '_cutMesh.mat']) ;
tubi.fullFileBase.phi0fit = ...
    fullfile(tubi.dir.spcutMesh, ['phi0s_' tubi.timeStampStringSpec '_%02d.png']) ; 
tubi.fullFileBase.clineDVhoop = ...
    fullfile(tubi.dir.clineDVhoop,...
    ['centerline_from_DVhoops_' tubi.timeStampStringSpec '.mat']);

%% filenames for writhe of the parameterized mesh
tubi.fileName.writhe = fullfile(tubi.dir.writhe, ...
    ['writhe_sphi' uvexten '_avgpts.mat']) ;

%% FEATURES
% Define any features within specialized scripts
% features =
featuresDir = fullfile(uvDir, 'features') ;
tubi.dir.features = featuresDir ;
tubi.fileName.features = struct() ;
tubi.fileName.features.folds = fullfile(featuresDir, ...
  ['fold_locations_sphi' uvexten '_avgpts.mat']) ;
tubi.fileName.features.features = ...
  fullfile(featuresDir, ['lobe_dynamics' uvexten '.mat']) ;

%% uvcutMesh
tubi.fullFileBase.uvcutMesh = ...
    fullfile(uvmeshDir, ['uvcutMesh_' tubi.timeStampStringSpec '.mat']) ;
tubi.fileBase.uvcutMesh = ['uvcutMesh_' tubi.timeStampStringSpec] ;

%%  spcutMesh and pullbacks
tubi.fullFileBase.spcutMesh = ...
    fullfile(sphiDir, ['spcutMesh_' tubi.timeStampStringSpec '.mat']) ;
tubi.fileBase.spcutMesh = ['spcutMesh_' tubi.timeStampStringSpec] ;
if dynamic
    tubi.fullFileBase.spcutMeshSm = ...
        fullfile(sphiSmDir, ['spcutMeshSm_' tubi.timeStampStringSpec '.mat']) ;
    tubi.fileBase.spcutMeshSm = ['spcutMeshSm_' tubi.timeStampStringSpec] ;
    tubi.fullFileBase.spcutMeshSmRS = ...
        fullfile(sphiSmRSDir, ['spcutMeshSmRS_' tubi.timeStampStringSpec '.mat']) ;
    tubi.fileBase.spcutMeshSmRS = ['spcutMeshSmRS_' tubi.timeStampStringSpec] ;
    
    tubi.fullFileBase.spcutMeshSmRSC = ...
        fullfile(sphiSmRSCDir, ['spcMSmRSC_' tubi.timeStampStringSpec '.mat']) ;
    tubi.fullFileBase.spcutMeshSmRSCPLY = ...
        fullfile(sphiSmRSCDir, ['spcMSmRSC_' tubi.timeStampStringSpec '.ply']) ;
    tubi.fileBase.spcutMeshSmRSC = ['spcMSmRSC_' tubi.timeStampStringSpec] ;
else
    tubi.fullFileBase.spcutMeshRS = ...
        fullfile(sphiRSDir, ['spcutMeshRS_' tubi.timeStampStringSpec '.mat']) ;
    tubi.fileBase.spcutMeshRS = ['spcutMeshRS_' tubi.timeStampStringSpec] ;
    
    tubi.fullFileBase.spcutMeshRSC = ...
        fullfile(sphiRSCDir, ['spcMRSC_' tubi.timeStampStringSpec '.mat']) ;
    tubi.fullFileBase.spcutMeshRSCPLY = ...
        fullfile(sphiRSCDir, ['spcMRSC_' tubi.timeStampStringSpec '.ply']) ;
    tubi.fileBase.spcutMeshRSC = ['spcMRSC_' tubi.timeStampStringSpec] ;
end

% uv
tubi.fileBase.im_uv = [tubi.fileBase.name, '_pbuv.tif'] ;
tubi.fullFileBase.im_uv = fullfile(tubi.dir.im_uv, tubi.fileBase.im_uv) ;
% relaxed sp
tubi.fileBase.im_r = [tubi.fileBase.name, '_pbr.tif'] ;
tubi.fullFileBase.im_r = fullfile(tubi.dir.im_r, tubi.fileBase.im_r) ;
tubi.fileBase.im_re = [tubi.fileBase.name, '_pbre.tif'] ;
tubi.fullFileBase.im_re = fullfile(tubi.dir.im_re, tubi.fileBase.im_re) ;
tubi.fileBase.im_sp = [tubi.fileBase.name, '_pbsp.tif'] ;
tubi.fullFileBase.im_sp = fullfile(tubi.dir.im_sp, tubi.fileBase.im_sp);
tubi.fileBase.im_spe = [tubi.fileBase.name, '_pbspe.tif'] ;
tubi.fullFileBase.im_spe = fullfile(tubi.dir.im_spe, tubi.fileBase.im_spe);
tubi.fileBase.im_speLUT = [tubi.fileBase.name, '_pbspe_LUT.tif'] ;
tubi.fullFileBase.im_speLUT = fullfile(tubi.dir.im_spe, tubi.fileBase.im_speLUT);
tubi.fileBase.im_up = [tubi.fileBase.name, '_pbup.tif'] ;
tubi.fullFileBase.im_up = fullfile(tubi.dir.im_up, tubi.fileBase.im_up) ;
tubi.fileBase.im_ricci = [tubi.fileBase.name, '_ricci.tif'] ;
tubi.fullFileBase.im_ricci = fullfile(tubi.dir.im_ricci, tubi.fileBase.im_ricci) ;

if dynamic
    tubi.fileBase.im_pivPathlines = [tubi.fileBase.name, '_pivPathlines.tif'] ;
    tubi.fullFileBase.im_pivPathlines = fullfile(tubi.dir.im_pivPathlines, tubi.fileBase.im_pivPathlines) ;
end

 %% Smoothed pullbacks (pb)
if dynamic
    
    tubi.fileBase.im_sp_sm = [tubi.fileBase.name, '_pbspsm.tif'] ;
    tubi.fileBase.im_sp_sme = [tubi.fileBase.name, '_pbspsme.tif'] ;
    tubi.fileBase.im_r_sm = [tubi.fileBase.name, '_pbrsm.tif'] ;
    tubi.fileBase.im_r_sme = [tubi.fileBase.name, '_pbrsme.tif'] ;

    tubi.fullFileBase.im_sp_sm = ...
        fullfile(tubi.dir.im_sp_sm, tubi.fileBase.im_sp_sm);
    tubi.fullFileBase.im_sp_sme = ...
        fullfile(tubi.dir.im_sp_sme, tubi.fileBase.im_sp_sme);
    tubi.fullFileBase.im_r_sm = ...
        fullfile(tubi.dir.im_r_sm, tubi.fileBase.im_r_sm);
    tubi.fullFileBase.im_r_sme = ...
        fullfile(tubi.dir.im_r_sme, tubi.fileBase.im_r_sme);
   
end


%% DYNAMICS
if dynamic
    %% PIV
    % By default, we use sp_sme as the piv coordinate system, but could be
    % other.
    tubi.dir.piv = struct() ;
    tubi.dir.piv.root = fullfile(uvDir, 'piv') ;  
    tubi.dir.piv.v3d = fullfile(tubi.dir.piv.root, 'piv3d') ;
    tubi.dir.piv.vt2d = fullfile(tubi.dir.piv.root, 'vt2d') ;
    tubi.dir.piv.vn2d = fullfile(tubi.dir.piv.root, 'vn2d') ;
    tubi.dir.piv.dilation = fullfile(tubi.dir.piv.root, 'dilation') ;
    tubi.fileName.pivRaw = struct() ;
    tubi.fileName.pivRaw.raw = fullfile(tubi.dir.piv.root, 'piv_results.mat') ;
    tubi.fullFileBase.piv3d = fullfile(tubi.dir.piv.v3d, 'piv3d_%04d.mat') ;
    
    % pathlines data
    tubi.dir.pathlines = struct() ;
    tubi.dir.pathlines.data = fullfile(tubi.dir.piv.root, 'pathlines', 't0_%04d') ; % HARDCODED TIMESTAMP FORMAT DELIMITER HERE FOR T0
    tubi.dir.pathlines.XY = fullfile(tubi.dir.pathlines.data, 'images_XY') ;
    tubi.dir.pathlines.XYZ = fullfile(tubi.dir.pathlines.data, 'images_XYZ') ;
    tubi.dir.pathlines.vXY = fullfile(tubi.dir.pathlines.data, 'images_vXY') ;
    tubi.dir.pathlines.fXY = fullfile(tubi.dir.pathlines.data, 'images_fXY') ;
    tubi.dir.pathlines.v3d = fullfile(tubi.dir.pathlines.data, 'images_v3d') ;
    tubi.dir.pathlines.f3d = fullfile(tubi.dir.pathlines.data, 'images_f3d') ;
    tubi.dir.pathlines.strain_images = fullfile(tubi.dir.pathlines.data, 'images_strain') ;
    tubi.dir.pathlines.strain = fullfile(tubi.dir.pathlines.data, 'strain') ;
    tubi.dir.pathlines.dxdyFiltered = fullfile(tubi.dir.pathlines.data, 'DxDyStrainFiltered') ;

    % pathline fileNames
    tubi.fileName.pathlines = struct() ;
    pdir = tubi.dir.pathlines.data ;
    tubi.fileName.pathlines.featureIDs = fullfile(pdir, 'featureIDs.mat') ;
    tubi.fileName.pathlines.XY = fullfile(pdir, 'piv_pathlines_XY.mat') ;
    tubi.fileName.pathlines.XYZ = fullfile(pdir, 'piv_pathlines_XYZ.mat') ;
    tubi.fileName.pathlines.vXY = fullfile(pdir, 'piv_pathlines_vXY.mat') ;
    tubi.fileName.pathlines.v3d = fullfile(pdir, 'piv_pathlines_v3d.mat') ;
    tubi.fileName.pathlines.fXY = fullfile(pdir, 'piv_pathlines_fXY.mat') ;
    tubi.fileName.pathlines.f3 = fullfile(pdir, 'piv_pathlines_f3d.mat') ;
    tubi.fileName.pathlines.refMesh = fullfile(pdir, 'refMesh.mat') ;
    % Ricci flow
    tubi.dir.pathlines.quasiconformal = fullfile(pdir, 'quasiconformal') ;    % for ricci comparison to refMesh ricci pullback
    tubi.dir.pathlines.ricci = struct() ;
    tubi.dir.pathlines.ricci.data = fullfile(pdir, 'ricci') ;             % for ricci flow on advected vertices themselves
    tubi.dir.pathlines.ricci.mu = fullfile(pdir, 'ricci', 'beltramiCoefficientsInstantaneous') ; 
    tubi.dir.pathlines.ricci.quasiconformal = fullfile(pdir, 'ricci', 'beltramiCoefficientsMaterial') ; 
    tubi.dir.pathlines.ricci.solution = fullfile(pdir, 'ricci', 'ricciSolutions') ; 
    tubi.dir.pathlines.ricci.mesh = fullfile(pdir, 'ricci', 'ricciMeshes') ; 
    % Strain and indentation
    tubi.dir.pathlines.radius = fullfile(pdir, 'images_radii_vertices') ;
    tubi.dir.pathlines.indentation = fullfile(pdir, 'images_indentation') ;
    tubi.dir.pathlines.kymographs = fullfile(pdir, 'kymographs') ;

    % pathline fileNames for conformal and kymo
    tubi.fileName.pathlines.quasiconformal = ...
        fullfile(tubi.dir.pathlines.quasiconformal, 'mu_v3dv2d.mat') ;
    tubi.fileName.pathlines.radius = fullfile(pdir, 'pathline_radii.mat') ;
    tubi.fileName.pathlines.indentation = fullfile(pdir, 'pathline_indentation.mat') ;
    tubi.fileName.pathlines.kymographs = struct() ;
    tubi.fileName.pathlines.kymographs.radius = ...
        fullfile(tubi.dir.pathlines.kymographs, 'radiusKymographs.mat') ;
    tubi.fileName.pathlines.kymographs.indentation = ...
        fullfile(tubi.dir.pathlines.kymographs, 'indentationKymographs.mat') ;
    tubi.fileName.pathlines.kymographs.mu = ...
        fullfile(tubi.dir.pathlines.kymographs, 'muKymographs.mat') ;

    %% Pathline-based strain measurement --> from pathline path
    tubi.fileBase.strain = ['strain_' tubi.timeStampStringSpec '.mat'] ;
    tubi.fullFileBase.pathlines.strain = fullfile(tubi.dir.pathlines.strain, ...
        tubi.fileBase.strain) ; 

    tubi.fullFileBase.pathlines.ricciQuasiconformal = fullfile(tubi.dir.pathlines.ricci.quasiconformal, ...
        ['mu_%04diter_' tubi.timeStampStringSpec '_v2dv2d.mat']) ;
    tubi.fullFileBase.pathlines.ricciMesh = fullfile(tubi.dir.pathlines.ricci.mesh, ...
        tubi.fileBase.ricciMesh) ;
    tubi.fullFileBase.pathlines.ricciSolution = fullfile(tubi.dir.pathlines.ricci.solution, ...
        tubi.fileBase.ricciSolution) ;
    tubi.fullFileBase.pathlines.ricciMu = fullfile(tubi.dir.pathlines.ricci.mu, ...
        tubi.fileBase.ricciMu) ;

    % pathline fundamental forms
    % fdir = tubi.dir.pathlines.fundForms ;
    % tubi.fileName.pathlines.fundForms = fullfile(fdir, 'fundForms.mat') ;
    % pathline velocities
    tubi.fileName.pathlines.velocities = struct() ; 
    pdir = tubi.dir.pathlines.data ;
    tubi.dir.pathlines.velocities = fullfile(pdir, 'velocities') ;
    pvdir = tubi.dir.pathlines.velocities ;
    tubi.fileName.pathlines.velocities.v3d = fullfile(pvdir, 'vM.mat')  ;
    tubi.fileName.pathlines.velocities.vn  = fullfile(pvdir, 'vnM.mat') ;
    tubi.fileName.pathlines.velocities.vv  = fullfile(pvdir, 'vvM.mat') ;
    tubi.fileName.pathlines.velocities.vf  = fullfile(pvdir, 'vfM.mat') ;
    tubi.fileName.pathlines.velocities.v2d = fullfile(pvdir, 'v2dM.mat') ;
    tubi.fileName.pathlines.velocities.v2dum = fullfile(pvdir, 'v2dMum.mat') ;
    tubi.fileName.pathlines.velocities.v3dsm   = fullfile(pvdir, 'vsmM.mat') ;
    tubi.fileName.pathlines.velocities.vnsm    = fullfile(pvdir, 'vnsmM.mat') ;
    tubi.fileName.pathlines.velocities.vvsm    = fullfile(pvdir, 'vvsmM.mat') ;
    tubi.fileName.pathlines.velocities.vfsm    = fullfile(pvdir, 'vfsmM.mat') ;
    tubi.fileName.pathlines.velocities.v2dsm   = fullfile(pvdir, 'v2dsmM.mat') ;
    tubi.fileName.pathlines.velocities.v2dsmum = fullfile(pvdir, 'v2dsmMum.mat') ;

    % pathline fileBases
    tubi.fileBase.pathlines = struct() ;
    tubi.fileBase.pathlines.strain = fullfile(tubi.dir.pathlines.strain, tubi.fileBase.strain) ;
    tubi.fileBase.pathlines.dxdyFiltered = fullfile(tubi.dir.pathlines.dxdyFiltered, ['dxdy_' tubi.fileBase.strain]) ;

    %% Velocities -- Lagrangian averaging
    tubi.dir.piv.avg = fullfile(tubi.dir.piv.root, 'lagrangianAvg') ;
    % tubi.dir.piv.avgCurl = fullfile(tubi.dir.piv.avg, 'curl') ;
    % tubi.dir.piv.avgDvg = fullfile(tubi.dir.piv.avg, 'dvg') ;
    tubi.fileName.pivAvg = struct() ;
    tubi.fileName.pivAvg.v2dum = fullfile(tubi.dir.piv.avg, 'v2dMum_avg.mat') ;
    tubi.fileName.pivAvg.v2d = fullfile(tubi.dir.piv.avg, 'v2dM_avg.mat') ;
    tubi.fileName.pivAvg.vn  = fullfile(tubi.dir.piv.avg, 'vnM_avg.mat') ;
    tubi.fileName.pivAvg.v3d = fullfile(tubi.dir.piv.avg, 'vM_avg.mat')  ;
    tubi.fileName.pivAvg.vv  = fullfile(tubi.dir.piv.avg, 'vvM_avg.mat') ;
    tubi.fileName.pivAvg.vf  = fullfile(tubi.dir.piv.avg, 'vfM_avg.mat') ;
    tubi.fileName.pivAvg.vv2d  = fullfile(tubi.dir.piv.avg, 'vv2dM_avg.mat') ;
    tubi.fileName.pivAvg.vf2d  = fullfile(tubi.dir.piv.avg, 'vf2dM_avg.mat') ;
    tubi.fileName.pivAvg.twist  = fullfile(tubi.dir.piv.avg, 'twist_dvphi_ds_from_pivAvg.mat') ;
    % Helmholtz-Hodge and DEC -- Lagrangian averaging
    tubi.dir.piv.avgDEC = struct() ;
    tubi.dir.piv.avgDEC.data   = fullfile(tubi.dir.piv.avg, 'dec') ;
    tubi.dir.piv.avgDEC.lap2d   = fullfile(tubi.dir.piv.avg, 'dec_lap2d') ;
    tubi.dir.piv.avgDEC.lap3d   = fullfile(tubi.dir.piv.avg, 'dec_lap3d') ;
    tubi.dir.piv.avgDEC.lap3dTexture = fullfile(tubi.dir.piv.avg, 'dec_lap3dTexture') ;
    tubi.dir.piv.avgDEC.div2d  = fullfile(tubi.dir.piv.avg, 'dec_div2d') ;
    tubi.dir.piv.avgDEC.div3d  = fullfile(tubi.dir.piv.avg, 'dec_div3d') ;
    tubi.dir.piv.avgDEC.div3dTexture = fullfile(tubi.dir.piv.avg, 'dec_div3dTexture') ;
    tubi.dir.piv.avgDEC.rot2d  = fullfile(tubi.dir.piv.avg, 'dec_rot2d') ;
    tubi.dir.piv.avgDEC.rot3d  = fullfile(tubi.dir.piv.avg, 'dec_rot3d') ;
    tubi.dir.piv.avgDEC.rot3dTexture = fullfile(tubi.dir.piv.avg, 'dec_rot3dTexture') ;
    tubi.dir.piv.avgDEC.harm2d = fullfile(tubi.dir.piv.avg, 'dec_harm2d') ;
    tubi.dir.piv.avgDEC.harm3d = fullfile(tubi.dir.piv.avg, 'dec_harm3d') ;
    tubi.fullFileBase.decAvg = fullfile(tubi.dir.piv.avgDEC.data, ...
                                  [tubi.fileBase.name '_dec.mat'] ) ;

end

%% PCA and mode analysis
tubi.dir.PCAoverTime = fullfile(tubi.dir.mesh, 'PCAoverTime') ;

%% LBS decomposition analysis
tubi.dir.LBSoverTime = fullfile(tubi.dir.mesh, 'LBSoverTime') ;

%% Ensure directories
if makeDirs
    dirs2make = struct2cell(tubi.dir) ;
    disp('creating dirs: ')
    disp(tubi.dir)
    for ii=1:length(dirs2make)
        dir2make = dirs2make{ii} ;
        if isa(dir2make, 'struct')
            dirfields = struct2cell(dir2make) ;
            for qq = 1:length(dirfields)
                dir2make = dirfields{qq} ;
                if isa(dir2make, 'struct')
                    dirfieldsSub = struct2cell(dir2make) ;
                    for pp = 1:length(dirfieldsSub)
                        dir2makeSub = dirfieldsSub{pp} ;
                        if ~exist(dir2makeSub, 'dir') && ~contains(dir2makeSub, '%') 
                            mkdir(dir2makeSub)
                        end
                    end
                else
                    if ~exist(dir2make, 'dir') && ~contains(dir2make, '%')
                        mkdir(dir2make)
                    end
                end
            end
        else
            if ~exist(dir2make, 'dir')
                try
                    mkdir(dir2make)
                catch
                    warning(['Could not make directory: ' dir2make])
                    mkdir(dir2make)
                end
            end
        end
    end
end


%% Segmentation / tracking -- could be extra functionality
tubi.dir.segmentation = fullfile(tubi.dir.mesh, 'cellSegmentation') ;
tubi.dir.tracking = fullfile(tubi.dir.mesh, 'cellTracking') ;

tubi.fullFileBase.segmentation = ...
     fullfile(tubi.dir.segmentation, ...
     [tubi.fileBase.name, '_Probabilities.h5']) ;
tubi.fileBase.segmentation2d = [tubi.fileBase.name, '_segmentation2d'] ;
tubi.fileBase.segmentation3d = [tubi.fileBase.name, '_segmentation3d'] ;
tubi.fullFileBase.segmentation2d = fullfile(tubi.dir.segmentation, 'seg2d', ...
     [tubi.fileBase.segmentation2d '.mat']) ;
tubi.fullFileBase.segmentation3d = fullfile(tubi.dir.segmentation, 'seg3d', ...
     [tubi.fileBase.segmentation3d '.mat']) ;
 
% corrected segmentations
tubi.fullFileBase.segmentation2dCorrected = fullfile(tubi.dir.segmentation, 'seg2d_corrected_%s', ...
     [tubi.fileBase.segmentation2d '.mat']) ;
tubi.fullFileBase.segmentation3dCorrected = fullfile(tubi.dir.segmentation, 'seg3d_corrected', ...
    [tubi.fileBase.segmentation3d '.mat']) ;
tubi.fullFileBase.segmentation2dCorrectedBinary = fullfile(tubi.dir.segmentation, 'seg2d_corrected_%s', ...
     'binary_maps', [tubi.fileBase.segmentation2d '_binary_map.png']) ;

% % nuclei only for voronoi measurement (may be identified through membrane training as non-membrane regions)
% tubi.fullFileBase.cellProbabilities = ...
%      fullfile(tubi.dir.cellProbabilities, ...
%      [tubi.fileBase.name, '_pbrsme_Probabilities.h5']) ;
% tubi.fileBase.cellID = [tubi.fileBase.name, '_pbrsme_cells'] ;
% tubi.fullFileBase.cellID = fullfile(tubi.dir.cellID, ...
%      [tubi.fileBase.cellID '.mat']) ;



%% Save t0 if supplied
if t0supplied_in_opts
    disp(['Writing t0 to disk: ' tubi.fileName.t0])
    write_txt_with_header(tubi.fileName.t0, tubi.t0, 't0, the reference timestamp')
end








