function initializeQuapSlap(QS, xp, opts)
%initializeQuapSlap(QS, xp, opts)
%   Hidden method for instantiating QuapSlap class
%
% Parameters
% ----------
% QS : QuapSlap object whose properties to fill in
% xp : Imsane Experiment class instance belonging to QS
% opts : struct with fields
%   xp : ImSAnE class object instance
%   flipy : bool
%   meshDir : str
%   timeUnits : str
%   spaceUnits : str
%   nU : int
%   nV : int
%   lambda : optional float 
%   lambda_mesh : optional float 
%   lambda_err : optional float
%   
%
% NPMitchell 2020

%% PROPERTIES
QS.xp = xp ;
QS.flipy = opts.flipy ;
meshDir = opts.meshDir ;
QS.timeUnits = opts.timeUnits ;
QS.spaceUnits = opts.spaceUnits ;
QS.fileBase.fn = xp.fileMeta.filenameFormat ;
QS.ssfactor = xp.detectOptions(1).ssfactor ;
QS.nU = opts.nU ;
QS.nV = opts.nV ;
if isfield(opts, 'timeInterval')
    QS.timeInterval = opts.timeInterval ;
elseif isfield(opts, 'timeinterval')
    QS.timeInterval = opts.timeinterval ;
else
    QS.timeInterval = 1 ;
end
dynamic = length(QS.xp.fileMeta.timePoints) > 1 ;
QS.dynamic = dynamic ; 

QS.normalShift = opts.normalShift ;
QS.a_fixed = opts.a_fixed ;
if isfield(opts, 'adjustlow')
    QS.data.adjustlow = opts.adjustlow ;
end
if isfield(opts, 'adjusthigh')
    QS.data.adjusthigh = opts.adjusthigh ;
end
if isfield(opts, 'axisOrder')
    QS.data.axisOrder = opts.axisOrder ;
end
if isfield(opts, 'ilastikOutputAxisOrder')
    QS.data.ilastikOutputAxisOrder = opts.ilastikOutputAxisOrder ;
end

% Assign which pullback coordsys is used for velocimetry
if isfield(opts, 'pivPullback')
    QS.pivPullback = pivPullback ;
end

%% NAMING
uvexten = sprintf('_nU%04d_nV%04d', QS.nU, QS.nV) ;

% Naming extension
QS.uvexten = uvexten ;

% APDV coordinate system
QS.APDV.resolution = min(xp.fileMeta.stackResolution) ;
QS.APDV.rot = [] ;
QS.APDV.trans = [] ;

% Plotting settings
QS.plotting.colors = define_colors() ;
QS.plotting.markers = {'o', 's', '^', 'v', '*', '>', '<'} ;

% Directories of the QS object
% Meshes and measurements before gridding into pullback coords
QS.dir.data = xp.fileMeta.dataDir ;
QS.fileName.apdv_options = fullfile(QS.dir.data, 'alignAPDV_Opts.mat') ;
QS.dir.mesh = meshDir ;
QS.dir.maskedData = fullfile(meshDir, 'maskedData') ;
QS.dir.alignedMesh = fullfile(meshDir, 'alignedMesh') ;
QS.dir.cntrline = fullfile(meshDir, 'centerline') ;
QS.dir.cylinderMesh = fullfile(meshDir, 'cylinderMesh') ;
QS.dir.cutMesh = fullfile(meshDir, 'cutMesh') ;
QS.dir.rawRicciMesh = fullfile(meshDir, 'rawRicciMesh') ;
QS.dir.cylinderMeshClean = fullfile(QS.dir.cylinderMesh, 'cleaned') ;
QS.dir.texturePatchIm = fullfile(meshDir, 'images_texturepatch') ;
QS.dir.mip = fullfile(meshDir, 'mips', 'dim%d_pages%04dto%04d') ;

% Mips
QS.fullFileBase.mip = fullfile(QS.dir.mip, 'mip_%06d.tif') ;

% After gridding into (u,v) / (zeta,phi) pullback coords
uvDir = fullfile(QS.dir.mesh, sprintf('gridCoords_nU%04d_nV%04d', QS.nU, QS.nV)) ;
QS.dir.uvCoord = uvDir ;


%% fileBases
QS.fileBase.name = xp.fileMeta.filenameFormat(1:end-4) ;
QS.fileBase.mesh = ...
    [xp.detector.options.ofn_smoothply '%06d'] ;
QS.fileBase.alignedMesh = ...
    [QS.fileBase.mesh '_APDV_um'] ;
QS.fileBase.apdProb = [QS.fileBase.name '_Probabilities_APDVcoords.h5']  ;
QS.fileBase.apCenterlineProb = [QS.fileBase.name '_Probabilities_apcenterline.h5']  ;
QS.fileBase.prob = [QS.fileBase.name '_Probabilities.h5'] ; 
QS.fileBase.centerlineXYZ = ...
    [QS.fileBase.mesh '_centerline_exp1p0_res*.txt' ] ;
QS.fileBase.centerlineAPDV = ...
    [QS.fileBase.mesh '_centerline_scaled_exp1p0_res*.txt' ] ;
QS.fileBase.cylinderMesh = ...
    [QS.fileBase.mesh '_cylindercut.ply'] ;
QS.fileBase.apBoundary = 'ap_boundary_indices_%06d.mat';
QS.fileBase.cylinderKeep = 'cylinderMesh_keep_indx_%06d.mat' ;
QS.fileName.apBoundaryDorsalPts = 'ap_boundary_dorsalpts.h5' ;

%% Metric
QS.dir.metric = struct() ;
QS.dir.metric.data = fullfile(QS.dir.uvCoord, 'metric', '%s_lmesh%0.3f') ;
QS.fullFileBase.metric = fullfile(QS.dir.metric.data, [QS.fileBase.name '.mat']) ;
QS.dir.metric.g_images2d = fullfile(QS.dir.metric.data, 'g_images2d') ;
QS.dir.metric.b_images2d = fullfile(QS.dir.metric.data, 'b_images2d') ;
QS.dir.metric.g_images3d = fullfile(QS.dir.metric.data, 'g_images3d') ;
QS.dir.metric.b_images3d = fullfile(QS.dir.metric.data, 'b_images3d') ;

%% define string for smoothing params
if isfield(opts, 'lambda')
    QS.smoothing.lambda = opts.lambda ;
end
if isfield(opts, 'lambda_mesh')
    QS.smoothing.lambda_mesh = opts.lambda_mesh ;
end
if isfield(opts, 'lambda_err')
    QS.smoothing.lambda_err = opts.lambda_err ;
end
if isfield(opts, 'nmodes')
    QS.smoothing.nmodes = opts.nmodes ;
end
if isfield(opts, 'zwidth')
    QS.smoothing.zwidth = opts.zwidth ;
end
l_lmesh_lerr =  strrep(sprintf('lambda%0.03f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d', ...
    QS.smoothing.lambda, QS.smoothing.lambda_mesh, ...
    QS.smoothing.lambda_err, QS.smoothing.nmodes, QS.smoothing.zwidth), '.', 'p') ;
l_lmesh = strrep(sprintf('lambda%0.03f_lmesh%0.3f_modes%02dw%02d', ...
    QS.smoothing.lambda, QS.smoothing.lambda_mesh, ...
    QS.smoothing.nmodes, QS.smoothing.zwidth), '.', 'p') ;

%% Dynamics -- strain rates, kinematics
if dynamic
    % Metric strain dirs
    QS.dir.metricKinematics = struct() ;
    QS.dir.metricKinematics.root = fullfile(uvDir, 'metricKinematics') ;
    QS.dir.metricKinematics.smoothing = fullfile(uvDir, 'metricKinematics', ...
        l_lmesh_lerr) ;
    QS.dir.metricKinematics.measurements = ...
        fullfile(uvDir, 'metricKinematics', l_lmesh_lerr, 'measurements') ;
    QS.fullFileBase.metricKinematics = struct() ;
    QS.fullFileBase.metricKinematics.divv = ...
        fullfile(QS.dir.metricKinematics.measurements, 'divv_vertices_%06d') ;
    QS.fullFileBase.metricKinematics.H2vn = ...
        fullfile(QS.dir.metricKinematics.measurements, 'H2vn_vertices_%06d') ;
    QS.fullFileBase.metricKinematics.gdot = ...
        fullfile(QS.dir.metricKinematics.measurements, 'gdot_vertices_%06d') ;

    % Metric Kinematics along pathlines
    QS.dir.metricKinematics.pathline = struct() ;
    QS.dir.metricKinematics.pathline.root = ...
        fullfile(uvDir, 'metricKinematics', l_lmesh_lerr, 'pathline_%04dt0') ;
    QS.dir.metricKinematics.pathline.measurements = ...
        fullfile(QS.dir.metricKinematics.pathline.root, 'measurements') ;
    QS.fileBase.metricKinematics = struct() ;
    QS.fileBase.metricKinematics.pathline.Hfn = 'HH_pathline%04d_%06d.mat' ;
    QS.fileBase.metricKinematics.pathline.gdot = 'gdot_pathline%04d_%06d.mat' ;
    QS.fileBase.metricKinematics.pathline.divv = 'divv_pathline%04d_%06d.mat' ;
    QS.fileBase.metricKinematics.pathline.veln = 'veln_pathline%04d_%06d.mat' ;
    QS.fileBase.metricKinematics.pathline.radius = 'radius_pathline%04d_%06d.mat' ;
    QS.fileBase.metricKinematics.pathline.kymographs = struct() ;
    QS.fileBase.metricKinematics.pathline.kymographs.ap = 'apKymographMetricKinematics.mat' ;
    QS.fileBase.metricKinematics.pathline.kymographs.d = 'dKymographMetricKinematics.mat' ;
    QS.fileBase.metricKinematics.pathline.kymographs.v = 'vKymographMetricKinematics.mat' ;
    QS.fileBase.metricKinematics.pathline.kymographs.l = 'lKymographMetricKinematics.mat' ;
    QS.fileBase.metricKinematics.pathline.kymographs.r = 'rKymographMetricKinematics.mat' ;

    % Metric deformation ('metric strain' or 'meshMetric)
    QS.dir.gstrain = fullfile(uvDir, 'metricStrain') ;
    QS.dir.gstrainRate = fullfile(QS.dir.gstrain, 'rateMetric') ;
    QS.dir.gstrainRateIm = fullfile(QS.dir.gstrainRate, 'images') ;
    QS.dir.gstrainMesh = fullfile(QS.dir.gstrain, 'meshMetric') ;
    QS.dir.gstrainMeshIm = fullfile(QS.dir.gstrainMesh, 'images') ;
    QS.fileBase.gstrainMesh = 'gstrainMesh_%06d.mat' ;
    QS.fileBase.gstrainRate = 'gstrainRate_%06d.mat' ;
    QS.fullFileBase.gstrainMesh = fullfile(QS.dir.gstrainMesh, ...
        QS.fileBase.gstrainMesh) ; 
    QS.fullFileBase.gstrainRate = fullfile(QS.dir.gstrainRate, ...
        QS.fileBase.gstrainRate) ; 

    %% Strain Rate
    QS.dir.strainRate = struct() ;
    QS.dir.strainRate.root = fullfile(uvDir, 'strainRate') ;
    QS.dir.strainRate.smoothing = fullfile(uvDir, 'strainRate', l_lmesh) ;
    QS.dir.strainRate.measurements = ...
        fullfile(uvDir, 'strainRate', l_lmesh, 'measurements') ;
    QS.fileBase.strainRate = 'strainRate_%06d.mat' ;
    QS.fullFileBase.strainRate = fullfile(QS.dir.strainRate.measurements, ...
        QS.fileBase.strainRate) ; 
    %% Pathline-based strain measurement --> from pathline path
    % see pathline section

    %% Strain rates along pathlines -- to be filled in with smoothing and t0
    QS.dir.strainRate.pathline = struct() ;
    QS.dir.strainRate.pathline.root = ...
        fullfile(uvDir, 'strainRate', l_lmesh, 'pathline_%04dt0') ;
    QS.dir.strainRate.pathline.measurements = ...
        fullfile(QS.dir.strainRate.pathline.root, 'measurements') ;
    
    
    %% Forces from stokes force balance    
    % Metric strain dirs
    QS.dir.stokesForces = struct() ;
    QS.dir.stokesForces.root = fullfile(uvDir, 'stokesForces') ;
    QS.dir.stokesForces.smoothing = fullfile(uvDir, 'stokesForces', ...
        l_lmesh_lerr) ;
    QS.dir.stokesForces.measurements = ...
        fullfile(uvDir, 'stokesForces', l_lmesh_lerr, 'measurements') ;
    QS.fullFileBase.stokesForces = struct() ;
    QS.fullFileBase.stokesForces.Lapv = ...
        fullfile(QS.dir.metricKinematics.measurements, 'Lapv_vertices_%06d') ;
    QS.fullFileBase.stokesForces.Kv = ...
        fullfile(QS.dir.metricKinematics.measurements, 'Kv_vertices_%06d') ;
    QS.fullFileBase.stokesForces.gradP = ...
        fullfile(QS.dir.metricKinematics.measurements, 'gradP_vertices_%06d') ;

end

% shorten variable names for brevity
clineDir = QS.dir.cntrline ;


%% Clean Cylinder Mesh
QS.fileName.aBoundaryDorsalPtsClean = ...
    fullfile(QS.dir.cylinderMeshClean, 'adIDx.h5') ;
QS.fileName.pBoundaryDorsalPtsClean = ...
    fullfile(QS.dir.cylinderMeshClean, 'pdIDx.h5') ;

%% cutMesh
QS.fullFileBase.cutPath = fullfile(QS.dir.cutMesh, 'cutPaths_%06d.txt') ;

%% fileNames
nshift = strrep(sprintf('%03d', QS.normalShift), '-', 'n') ;
shiftstr = ['_' nshift 'step'] ;
QS.fileName.rot = fullfile(meshDir, 'rotation_APDV.txt') ;
QS.fileName.trans = fullfile(meshDir, 'translation_APDV.txt') ;
QS.fileName.xyzlim_raw = fullfile(meshDir, 'xyzlim_raw.txt') ;
QS.fileName.xyzlim_pix = fullfile(meshDir, 'xyzlim_APDV.txt') ;
QS.fileName.xyzlim_um = ...
    fullfile(meshDir, 'xyzlim_APDV_um.txt') ;
QS.fileName.xyzlim_um_buff = ...
    fullfile(meshDir, ['xyzlim_APDV_um' shiftstr '.txt']) ;
% fileNames for APDV and cylinderMesh
QS.fullFileBase.apdProb = fullfile(QS.dir.data, QS.fileBase.apdProb) ;
QS.fullFileBase.apCenterlineProb = fullfile(QS.dir.data, QS.fileBase.apCenterlineProb) ;
QS.fullFileBase.prob = fullfile(QS.dir.data, QS.fileBase.prob) ;
QS.fileName.apdv = ...
    fullfile(clineDir, 'apdv_coms_from_training.h5') ;
QS.fileName.dcom = fullfile(meshDir, 'dcom_for_rot.txt') ;
QS.fileName.startendPt = fullfile(clineDir, 'startendpt.h5') ;
QS.fileName.cleanCntrlines = ...
    fullfile(clineDir, 'centerlines_anomalies_fixed.mat') ;
QS.fileName.apBoundaryDorsalPts = ...
    fullfile(QS.dir.cylinderMesh, 'ap_boundary_dorsalpts.h5') ;
QS.fileName.endcapOptions = ...
    fullfile(QS.dir.cylinderMesh, 'endcapOptions.mat') ;
QS.fileName.apdBoundary = ...
    fullfile(QS.dir.cylinderMesh, 'ap_boundary_dorsalpts.h5') ;

%% FileNamePatterns
QS.fullFileBase.mesh = ...
    fullfile(QS.dir.mesh, [QS.fileBase.mesh '.ply']) ;
QS.fullFileBase.alignedMesh = ...
    fullfile(QS.dir.alignedMesh, [QS.fileBase.alignedMesh '.ply']) ;
% fileNames for centerlines
QS.fullFileBase.centerlineXYZ = ...
    fullfile(clineDir, QS.fileBase.centerlineXYZ) ;
QS.fullFileBase.centerlineAPDV = ...
    fullfile(clineDir, QS.fileBase.centerlineAPDV) ;
QS.fullFileBase.cylinderMesh = ...
    fullfile(QS.dir.cylinderMesh, QS.fileBase.cylinderMesh) ;
QS.fullFileBase.apBoundary = ...
    fullfile(QS.dir.cylinderMesh, QS.fileBase.apBoundary) ;
% QS.fullFileBase.apBoundaryDorsalPts = QS.fileName.apBoundaryDorsalPts ;
QS.fullFileBase.cylinderKeep = ...
    fullfile(QS.dir.cylinderMesh, QS.fileBase.cylinderKeep) ;
QS.fullFileBase.cylinderMeshClean = ...
    fullfile(QS.dir.cylinderMesh, 'cleaned',...
    [QS.fileBase.mesh '_cylindercut_clean.ply']) ;            

%% Define cutMesh directories
% cutMesh = fullfile(meshDir, 'cutMesh') ;
% cutMeshBase = fullfile(cutMesh, [QS.fileBase.name, '_cutMesh.mat']) ;
imFolderBase = fullfile(uvDir, ['PullbackImages' shiftstr] ) ;
ricciMeshDir = fullfile(uvDir, ['ricci_cutMesh' shiftstr], 'noResampling') ;
ricciMeshDirWithResampling = fullfile(uvDir, ['ricci_cutMesh' shiftstr], ...
    'withResampling') ;
sphiDir = fullfile(uvDir, ['sphi_cutMesh' shiftstr]) ;
if dynamic
    sphiSmDir = fullfile(sphiDir, 'smoothed') ;
    sphiSmRSDir = fullfile(sphiDir, 'smoothed_rs') ;
    % sphiSmRSImDir = fullfile(sphiSmRSDir, 'images') ;
    % sphiSmRSPhiImDir = fullfile(sphiSmRSImDir, 'phicolor') ;
    sphiSmRSCDir = fullfile(sphiDir, 'smoothed_rs_closed') ;
    sphiSmDir2x = fullfile(sphiDir, 'smoothed_doubleResolution') ;
    sphiSmRSDir2x = fullfile(sphiDir, 'smoothed_rs_doubleResolution') ;
    sphiSmRSCDir2x = fullfile(sphiDir, 'smoothed_rs_closed_doubleResolution') ;
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
imFolder_uvprime = [imFolderBase '_uvprime'] ;
imFolder_uvprime_e = fullfile(imFolder_uvprime, 'extended') ;
imFolder_re_stack = fullfile([imFolderBase, '_sphi_relaxed'], ...
    'extended_stack') ;  
imFolder_spe_stack = fullfile([imFolderBase, '_sphi'], ...
    'extended_stack') ;  

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

% Lobe/fold identification paths
lobeDir = fullfile(uvDir, 'lobes') ;
foldHoopImDir = fullfile(lobeDir, 'constriction_hoops') ;
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
% Port into QS
% (already above) QS.dir.cutMesh = fullfile(meshDir, 'cutMesh') ;
QS.dir.spcutMesh = sphiDir ;
if dynamic
    QS.dir.spcutMeshSm = sphiSmDir ;
    QS.dir.spcutMeshSmRS = sphiSmRSDir ;
    QS.dir.spcutMeshSmRSC = sphiSmRSCDir ;
else
    QS.dir.spcutMeshRS = sphiRSDir ;
    QS.dir.spcutMeshRSC = sphiRSCDir ;
end
QS.dir.ricci = struct() ;
QS.dir.ricci.data = ricciMeshDir ;
QS.dir.ricci.mesh = fullfile(ricciMeshDir, 'meshes') ;
QS.dir.ricci.solution = fullfile(ricciMeshDir, 'ricciSolutions') ;
QS.dir.ricci.mu = fullfile(ricciMeshDir, 'beltramiCoefficients') ;

QS.dir.ricci.dataWithResampling = ricciMeshDirWithResampling ;
QS.dir.ricci.meshWithResampling = fullfile(ricciMeshDirWithResampling, 'meshes') ;
QS.dir.ricci.solutionWithResampling = fullfile(ricciMeshDirWithResampling, 'ricciSolutions') ;
QS.dir.ricci.muWithResampling = fullfile(ricciMeshDirWithResampling, 'beltramiCoefficients') ;

QS.dir.rawRicci = struct() ;
rawRicciMeshDir = QS.dir.rawRicciMesh ;
QS.dir.rawRicci.data = QS.dir.rawRicciMesh ;
QS.dir.rawRicci.mesh = fullfile(rawRicciMeshDir, 'meshes') ;
QS.dir.rawRicci.solution = fullfile(rawRicciMeshDir, 'ricciSolutions') ;
QS.dir.rawRicci.mu = fullfile(rawRicciMeshDir, 'beltramiCoefficients') ;

%% UVprime coordinates (twist-adjusted conformal maps of sp-smoothed meshes)
QS.dir.uvpcutMesh = fullfile(uvDir, ['uvprime_cutMesh' shiftstr]) ;
QS.fileBase.uvpcutMesh = 'uvpcutMesh_%06d' ;
QS.fullFileBase.uvpcutMesh = fullfile(QS.dir.uvpcutMesh, ...
    [QS.fileBase.uvpcutMesh '.mat']) ;

%% Raw Ricci mesh coordinates (truly conformal via Ricci flow from cleaned cylinder mesh)
QS.fileBase.rawRicciMesh = 'rawRicciMesh_%04diter_%06d.mat' ;
QS.fullFileBase.rawRicciMesh = fullfile(QS.dir.rawRicci.mesh, QS.fileBase.rawRicciMesh) ;
QS.fileBase.rawRicciSolution = 'rawRicciSolution_%04diter_%06d.mat' ;
QS.fullFileBase.rawRicciSolution = fullfile(QS.dir.rawRicci.solution, QS.fileBase.rawRicciSolution) ;
QS.fileBase.rawRicciMu = 'rawRicciMesh_mus_%04diter_%06d.mat' ;
QS.fullFileBase.rawRicciMu = fullfile(QS.dir.rawRicci.mu, QS.fileBase.rawRicciMu) ;

%% Ricci mesh coordinates (truly conformal via Ricci flow from reparameterized gridCoords mesh (spsm uv or sp))
QS.fileBase.ricciMesh = 'ricciMesh_%04diter_%06d.mat' ;
QS.fullFileBase.ricciMesh = fullfile(QS.dir.ricci.mesh, QS.fileBase.ricciMesh) ;
QS.fileBase.ricciSolution = 'ricciSolution_%04diter_%06d.mat' ;
QS.fullFileBase.ricciSolution = fullfile(QS.dir.ricci.solution, QS.fileBase.ricciSolution) ;
QS.fileBase.ricciMu = 'ricciMesh_mus_%04diter_%06d.mat' ;
QS.fullFileBase.ricciMu = fullfile(QS.dir.ricci.mu, QS.fileBase.ricciMu) ;

%% double resolution cutMeshes
if dynamic
    QS.dir.spcutMeshSm2x = sphiSmDir2x ;
    QS.dir.spcutMeshSmRS2x = sphiSmRSDir2x ;
    QS.dir.spcutMeshSmRSC2x = sphiSmRSCDir2x ;
end

%% centerlines
QS.dir.clineDVhoop = ...
    fullfile(QS.dir.uvCoord, ...
    ['centerline_from_DVhoops' shiftstr]) ;
QS.dir.writhe =  fullfile(QS.dir.clineDVhoop, 'writhe') ;

%% Images
QS.dir.im_uv = [imFolderBase '_uv'] ;
QS.dir.im_uve = [imFolderBase '_uv_extended'] ;
QS.dir.im_r = [imFolderBase '_sphi_relaxed'] ;
QS.dir.im_re = fullfile(QS.dir.im_r, 'extended') ;
QS.dir.im_r_uvprime = [imFolderBase '_uvprime_relaxed'] ;
QS.dir.im_sp = imFolder_sp ;
QS.dir.im_spe = imFolder_spe ;
QS.dir.im_spe_stack = imFolder_spe_stack ;
QS.dir.im_up = imFolder_up ;
QS.dir.im_upe = imFolder_upe ;
if dynamic
    QS.dir.im_sp_sm = imFolder_spsm ;
    QS.dir.im_sp_sme = imFolder_spsme ;
    QS.dir.im_sp_smeLUT = imFolder_spsmeLUT ;
    QS.dir.im_r_sm = imFolder_rsm ;
    QS.dir.im_r_sme = imFolder_rsme ;
    QS.dir.im_r_sme_stack = imFolder_rsme_stack ;
end
QS.dir.im_uvprime = imFolder_uvprime ;
QS.dir.im_uvprime_e = imFolder_uvprime_e ;
QS.dir.im_ricci = imFolder_ricci ;
QS.dir.im_ricci_e = imFolder_ricci_e ;
QS.dir.im_re_stack = imFolder_re_stack ;
QS.dir.segmentation = fullfile(QS.dir.mesh, 'cellSegmentation') ;
QS.dir.tracking = fullfile(QS.dir.mesh, 'cellTracking') ;

if dynamic
    QS.dir.cellProbabilities = fullfile(QS.dir.im_r_sme, 'cellProbabilities') ;
    QS.dir.cellID = fullfile(QS.dir.im_r_sme, 'cellID') ;
else
    QS.dir.cellProbabilities = fullfile(QS.dir.im_re, 'cellProbabilities') ;
    QS.dir.cellID = fullfile(QS.dir.im_re, 'cellID') ;    
end

QS.dir.lobe = lobeDir ;
QS.dir.foldHoopIm = foldHoopImDir ;
QS.dir.curvatures = KHSmDir ;
QS.dir.meanCurvature = HSmDir ;
QS.dir.gaussCurvature = KSmDir ;
QS.fullFileBase.curvatures = fullfile(QS.dir.curvatures, 'gauss_mean_curvature_%06d.mat') ;
QS.fullFileBase.cutMesh = ...
    fullfile(QS.dir.cutMesh, [QS.fileBase.name, '_cutMesh.mat']) ;
QS.fullFileBase.phi0fit = ...
    fullfile(QS.dir.spcutMesh, 'phi0s_%06d_%02d.png') ; 
QS.fullFileBase.clineDVhoop = ...
    fullfile(QS.dir.clineDVhoop,...
    'centerline_from_DVhoops_%06d.mat');

%% filenames for lobe dynamics/statics
QS.fileName.fold = fullfile(lobeDir, ...
    ['fold_locations_sphi' uvexten '_avgpts.mat']) ;
QS.fileName.lobeDynamics = ...
    fullfile(lobeDir, ['lobe_dynamics' uvexten '.mat']) ;
QS.fileName.writhe = fullfile(QS.dir.writhe, ...
    ['writhe_sphi' uvexten '_avgpts.mat']) ;

%%  spcutMesh and pullbacks
QS.fullFileBase.spcutMesh = ...
    fullfile(sphiDir, 'spcutMesh_%06d.mat') ;
QS.fileBase.spcutMesh = 'spcutMesh_%06d' ;
if dynamic
    QS.fullFileBase.spcutMeshSm = ...
        fullfile(sphiSmDir, 'spcutMeshSm_%06d.mat') ;
    QS.fileBase.spcutMeshSm = 'spcutMeshSm_%06d' ;
    QS.fullFileBase.spcutMeshSmRS = ...
        fullfile(sphiSmRSDir, 'spcutMeshSmRS_%06d.mat') ;
    QS.fileBase.spcutMeshSmRS = 'spcutMeshSmRS_%06d' ;
    
    QS.fullFileBase.spcutMeshSmRSC = ...
        fullfile(sphiSmRSCDir, 'spcMSmRSC_%06d.mat') ;
    QS.fullFileBase.spcutMeshSmRSCPLY = ...
        fullfile(sphiSmRSCDir, 'spcMSmRSC_%06d.ply') ;
    QS.fileBase.spcutMeshSmRSC = 'spcMSmRSC_%06d' ;
    
    %% double resolution spcutMeshSm's
    QS.fullFileBase.spcutMeshSm2x = ...
        fullfile(sphiSmDir2x, 'spcutMeshSm2x_%06d.mat') ;
    QS.fileBase.spcutMeshSm2x = 'spcutMeshSm2x_%06d' ;
    QS.fullFileBase.spcutMeshSmRS2x = ...
        fullfile(sphiSmRSDir2x, 'spcutMeshSmRS2x_%06d.mat') ;
    QS.fileBase.spcutMeshSmRS2x = 'spcutMeshSmRS2x_%06d' ;
    QS.fullFileBase.spcutMeshSmRSC2x = ...
        fullfile(sphiSmRSCDir2x, 'spcMSmRSC2x_%06d.mat') ;
    QS.fullFileBase.spcutMeshSmRSCPLY2x = ...
        fullfile(sphiSmRSCDir2x, 'spcMSmRSC2x_%06d.ply') ;
    QS.fileBase.spcutMeshSmRSC2x = 'spcMSmRSC2x_%06d' ;
else
    QS.fullFileBase.spcutMeshRS = ...
        fullfile(sphiRSDir, 'spcutMeshRS_%06d.mat') ;
    QS.fileBase.spcutMeshRS = 'spcutMeshRS_%06d' ;
    
    QS.fullFileBase.spcutMeshRSC = ...
        fullfile(sphiRSCDir, 'spcMRSC_%06d.mat') ;
    QS.fullFileBase.spcutMeshRSCPLY = ...
        fullfile(sphiRSCDir, 'spcMRSC_%06d.ply') ;
    QS.fileBase.spcutMeshRSC = 'spcMRSC_%06d' ;
end

% uv
QS.fileBase.im_uv = [QS.fileBase.name, '_pbuv.tif'] ;
QS.fullFileBase.im_uv = fullfile(QS.dir.im_uv, QS.fileBase.im_uv) ;
% relaxed sp
QS.fileBase.im_r = [QS.fileBase.name, '_pbr.tif'] ;
QS.fullFileBase.im_r = fullfile(QS.dir.im_r, QS.fileBase.im_r) ;
QS.fileBase.im_re = [QS.fileBase.name, '_pbre.tif'] ;
QS.fullFileBase.im_re = fullfile(QS.dir.im_re, QS.fileBase.im_re) ;
QS.fileBase.im_sp = [QS.fileBase.name, '_pbsp.tif'] ;
QS.fullFileBase.im_sp = fullfile(QS.dir.im_sp, QS.fileBase.im_sp);
QS.fileBase.im_spe = [QS.fileBase.name, '_pbspe.tif'] ;
QS.fullFileBase.im_spe = fullfile(QS.dir.im_spe, QS.fileBase.im_spe);
QS.fileBase.im_speLUT = [QS.fileBase.name, '_pbspe_LUT.tif'] ;
QS.fullFileBase.im_speLUT = fullfile(QS.dir.im_spe, QS.fileBase.im_speLUT);
QS.fileBase.im_up = [QS.fileBase.name, '_pbup.tif'] ;
QS.fullFileBase.im_up = fullfile(QS.dir.im_up, QS.fileBase.im_up) ;
QS.fileBase.im_ricci = [QS.fileBase.name, '_ricci.tif'] ;
QS.fullFileBase.im_ricci = fullfile(QS.dir.im_ricci, QS.fileBase.im_ricci) ;


 %% Smoothed pullbacks (pb)
if dynamic
    if length(QS.xp.expMeta.channelsUsed) > 1
        
        QS.fileBase.im_sp_sm = [ 'Time_%06d_stab_pbspsm.tif'] ;
        QS.fileBase.im_sp_sme = [ 'Time_%06d_stab_pbspsme.tif'] ;
        QS.fileBase.im_r_sm = [ 'Time_%06d_stab_pbrsm.tif'] ;
        QS.fileBase.im_r_sme = [ 'Time_%06d_stab_pbrsme.tif'] ;
    else
        QS.fileBase.im_sp_sm = [QS.fileBase.name, '_pbspsm.tif'] ;
        QS.fileBase.im_sp_sme = [QS.fileBase.name, '_pbspsme.tif'] ;
        QS.fileBase.im_r_sm = [QS.fileBase.name, '_pbrsm.tif'] ;
        QS.fileBase.im_r_sme = [QS.fileBase.name, '_pbrsme.tif'] ;
    end
    QS.fullFileBase.im_sp_sm = ...
        fullfile(QS.dir.im_sp_sm, QS.fileBase.im_sp_sm);
    QS.fullFileBase.im_sp_sme = ...
        fullfile(QS.dir.im_sp_sme, QS.fileBase.im_sp_sme);
    QS.fullFileBase.im_r_sm = ...
        fullfile(QS.dir.im_r_sm, QS.fileBase.im_r_sm);
    QS.fullFileBase.im_r_sme = ...
        fullfile(QS.dir.im_r_sme, QS.fileBase.im_r_sme);
   
end

% uvprime
QS.fileBase.im_uvprime = [QS.fileBase.name, '_pbuvprime.tif'] ;
QS.fullFileBase.im_uvprime = fullfile(QS.dir.im_uvprime, QS.fileBase.im_uvprime) ;
QS.fileBase.im_uvprime_e = [QS.fileBase.name, '_pbuvpe.tif'] ;
QS.fullFileBase.im_uvprime_e = fullfile(QS.dir.im_uvprime_e, QS.fileBase.im_uvprime_e) ;
QS.fileBase.im_r_uvprime = [QS.fileBase.name, '_pbruvprime.tif'] ;
QS.fullFileBase.im_r_uvprime = fullfile(QS.dir.im_r_uvprime, QS.fileBase.im_r_uvprime) ;
% ricci
QS.fileBase.im_uvprime_e = [QS.fileBase.name, '_pbuvpe.tif'] ;
QS.fullFileBase.im_uvprime_e = fullfile(QS.dir.im_uvprime_e, QS.fileBase.im_uvprime_e) ;
QS.fileBase.im_r_uvprime = [QS.fileBase.name, '_pbruvprime.tif'] ;
QS.fullFileBase.im_r_uvprime = fullfile(QS.dir.im_r_uvprime, QS.fileBase.im_r_uvprime) ;

%% Cells segmentation / nuclei
QS.fullFileBase.segmentation = ...
     fullfile(QS.dir.segmentation, ...
     [QS.fileBase.name, '_Probabilities.h5']) ;
QS.fileBase.segmentation2d = [QS.fileBase.name, '_segmentation2d'] ;
QS.fileBase.segmentation3d = [QS.fileBase.name, '_segmentation3d'] ;
QS.fullFileBase.segmentation2d = fullfile(QS.dir.segmentation, 'seg2d', ...
     [QS.fileBase.segmentation2d '.mat']) ;
QS.fullFileBase.segmentation3d = fullfile(QS.dir.segmentation, 'seg3d', ...
     [QS.fileBase.segmentation3d '.mat']) ;
 
% corrected segmentations
QS.fullFileBase.segmentation2dCorrected = fullfile(QS.dir.segmentation, 'seg2d_corrected_%s', ...
     [QS.fileBase.segmentation2d '.mat']) ;
QS.fullFileBase.segmentation3dCorrected = fullfile(QS.dir.segmentation, 'seg3d_corrected', ...
    [QS.fileBase.segmentation3d '.mat']) ;
QS.fullFileBase.segmentation2dCorrectedBinary = fullfile(QS.dir.segmentation, 'seg2d_corrected_%s', ...
     'binary_maps', [QS.fileBase.segmentation2d '_binary_map.png']) ;
 
% nuclei only for voronoi measurement (could be ID'd through membrane training)
QS.fullFileBase.cellProbabilities = ...
     fullfile(QS.dir.cellProbabilities, ...
     [QS.fileBase.name, '_pbrsme_Probabilities.h5']) ;
QS.fileBase.cellID = [QS.fileBase.name, '_pbrsme_cells'] ;
QS.fullFileBase.cellID = fullfile(QS.dir.cellID, ...
     [QS.fileBase.cellID '.mat']) ;
 
%% Thickness
QS.dir.thickness = fullfile(QS.dir.uvCoord, 'thickness', ...
    '%s', 'stack_%02d_%02d_%0.2fum') ; 
% --> coordSys, n_outward, n_inward, stepSize (layerSpacing * QS.APDV.resolution)
QS.fullFileBase.thickness = struct() ;
QS.fullFileBase.thickness.ims2d = fullfile(QS.dir.thickness, 'images2d', ...
    [QS.fileBase.name '_thickness_2d.png']) ;
QS.fullFileBase.thickness.ims3d = fullfile(QS.dir.thickness, 'images3d', ...
    [QS.fileBase.name '_thickness_3d.png']) ;
QS.fullFileBase.thickness.data = fullfile(QS.dir.thickness, ...
    [QS.fileBase.name '_thickness.mat']) ;

%% DYNAMICS
if dynamic
    %% PIV
    % By default, we use sp_sme as the piv coordinate system, but could be
    % other.
    QS.dir.piv = struct() ;
    QS.dir.piv.root = fullfile(uvDir, 'piv') ;  
    QS.dir.piv.v3d = fullfile(QS.dir.piv.root, 'piv3d') ;
    QS.dir.piv.vt2d = fullfile(QS.dir.piv.root, 'vt2d') ;
    QS.dir.piv.vn2d = fullfile(QS.dir.piv.root, 'vn2d') ;
    QS.dir.piv.dilation = fullfile(QS.dir.piv.root, 'dilation') ;
    QS.fileName.pivRaw = struct() ;
    QS.fileName.pivRaw.raw = fullfile(QS.dir.piv.root, 'piv_results.mat') ;
    QS.fullFileBase.piv3d = fullfile(QS.dir.piv.v3d, 'piv3d_%04d.mat') ;
    
    QS.dir.pivMultiChannel = struct() ;
    QS.dir.pivMultiChannel.root = fullfile(uvDir, 'pivMultiChannel') ;  
    QS.dir.pivMultiChannel.v3d = ...
        fullfile(QS.dir.pivMultiChannel.root, 'piv3d') ;  
    QS.dir.pivMultiChannel.vt2d = ...
        fullfile(QS.dir.pivMultiChannel.root, 'vt2d') ;
    QS.dir.pivMultiChannel.vn2d = ...
        fullfile(QS.dir.pivMultiChannel.root, 'vn2d') ;
    QS.dir.pivMultiChannel.dilation = ...
        fullfile(QS.dir.pivMultiChannel.root, 'dilation') ;
    QS.dir.pivMultiChannel.relativeMotion = ...
        fullfile(QS.dir.pivMultiChannel.root, 'relativeMotion') ;
    QS.fileName.pivRawMultiChannel = struct() ;
    QS.fileName.pivRawMultiChannel.raw = ...
        fullfile(QS.dir.pivMultiChannel.root, 'piv_results_ch%d.mat') ;
    QS.fullFileBase.pivMultiChannel = struct() ;
    QS.fullFileBase.pivMultiChannel.v3d = ...
        fullfile(QS.dir.pivMultiChannel.v3d, 'piv3d_ch_%06d.mat') ;
    
    % pathlines data
    QS.dir.pathlines = struct() ;
    QS.dir.pathlines.data = fullfile(QS.dir.piv.root, 'pathlines', 't0_%04d') ; 
    QS.dir.pathlines.XY = fullfile(QS.dir.pathlines.data, 'images_XY') ;
    QS.dir.pathlines.XYZ = fullfile(QS.dir.pathlines.data, 'images_XYZ') ;
    QS.dir.pathlines.vXY = fullfile(QS.dir.pathlines.data, 'images_vXY') ;
    QS.dir.pathlines.fXY = fullfile(QS.dir.pathlines.data, 'images_fXY') ;
    QS.dir.pathlines.v3d = fullfile(QS.dir.pathlines.data, 'images_v3d') ;
    QS.dir.pathlines.f3d = fullfile(QS.dir.pathlines.data, 'images_f3d') ;
    QS.dir.pathlines.strain_images = fullfile(QS.dir.pathlines.data, 'images_strain') ;
    QS.dir.pathlines.strain = fullfile(QS.dir.pathlines.data, 'strain') ;
    QS.dir.pathlines.dxdyFiltered = fullfile(QS.dir.pathlines.data, 'DxDyStrainFiltered') ;

    % pathline fileNames
    QS.fileName.pathlines = struct() ;
    pdir = QS.dir.pathlines.data ;
    QS.fileName.pathlines.featureIDs = fullfile(pdir, 'featureIDs.mat') ;
    QS.fileName.pathlines.XY = fullfile(pdir, 'piv_pathlines_XY.mat') ;
    QS.fileName.pathlines.XYZ = fullfile(pdir, 'piv_pathlines_XYZ.mat') ;
    QS.fileName.pathlines.vXY = fullfile(pdir, 'piv_pathlines_vXY.mat') ;
    QS.fileName.pathlines.v3d = fullfile(pdir, 'piv_pathlines_v3d.mat') ;
    QS.fileName.pathlines.fXY = fullfile(pdir, 'piv_pathlines_fXY.mat') ;
    QS.fileName.pathlines.f3 = fullfile(pdir, 'piv_pathlines_f3d.mat') ;
    QS.fileName.pathlines.refMesh = fullfile(pdir, 'refMesh.mat') ;
    % Ricci flow
    QS.dir.pathlines.quasiconformal = fullfile(pdir, 'quasiconformal') ;    % for ricci comparison to refMesh ricci pullback
    QS.dir.pathlines.ricci = struct() ;
    QS.dir.pathlines.ricci.data = fullfile(pdir, 'ricci') ;             % for ricci flow on advected vertices themselves
    QS.dir.pathlines.ricci.mu = fullfile(pdir, 'ricci', 'beltramiCoefficientsInstantaneous') ; 
    QS.dir.pathlines.ricci.quasiconformal = fullfile(pdir, 'ricci', 'beltramiCoefficientsMaterial') ; 
    QS.dir.pathlines.ricci.solution = fullfile(pdir, 'ricci', 'ricciSolutions') ; 
    QS.dir.pathlines.ricci.mesh = fullfile(pdir, 'ricci', 'ricciMeshes') ; 
    % Strain and indentation
    QS.dir.pathlines.radius = fullfile(pdir, 'images_radii_vertices') ;
    QS.dir.pathlines.indentation = fullfile(pdir, 'images_indentation') ;
    QS.dir.pathlines.kymographs = fullfile(pdir, 'kymographs') ;

    % pathline fileNames for conformal and kymo
    QS.fileName.pathlines.quasiconformal = ...
        fullfile(QS.dir.pathlines.quasiconformal, 'mu_v3dv2d.mat') ;
    QS.fileName.pathlines.radius = fullfile(pdir, 'pathline_radii.mat') ;
    QS.fileName.pathlines.indentation = fullfile(pdir, 'pathline_indentation.mat') ;
    QS.fileName.pathlines.kymographs = struct() ;
    QS.fileName.pathlines.kymographs.radius = ...
        fullfile(QS.dir.pathlines.kymographs, 'radiusKymographs.mat') ;
    QS.fileName.pathlines.kymographs.indentation = ...
        fullfile(QS.dir.pathlines.kymographs, 'indentationKymographs.mat') ;
    QS.fileName.pathlines.kymographs.mu = ...
        fullfile(QS.dir.pathlines.kymographs, 'muKymographs.mat') ;

    %% Pathline-based strain measurement --> from pathline path
    QS.fileBase.strain = 'strain_%06d.mat' ;
    QS.fullFileBase.pathlines.strain = fullfile(QS.dir.pathlines.strain, ...
        QS.fileBase.strain) ; 

    QS.fullFileBase.pathlines.ricciQuasiconformal = fullfile(QS.dir.pathlines.ricci.quasiconformal, ...
        'mu_%04diter_%06d_v2dv2d.mat') ;
    QS.fullFileBase.pathlines.ricciMesh = fullfile(QS.dir.pathlines.ricci.mesh, ...
        QS.fileBase.ricciMesh) ;
    QS.fullFileBase.pathlines.ricciSolution = fullfile(QS.dir.pathlines.ricci.solution, ...
        QS.fileBase.ricciSolution) ;
    QS.fullFileBase.pathlines.ricciMu = fullfile(QS.dir.pathlines.ricci.mu, ...
        QS.fileBase.ricciMu) ;

    % pathline fundamental forms
    % fdir = QS.dir.pathlines.fundForms ;
    % QS.fileName.pathlines.fundForms = fullfile(fdir, 'fundForms.mat') ;
    % pathline velocities
    QS.fileName.pathlines.velocities = struct() ; 
    pdir = QS.dir.pathlines.data ;
    QS.dir.pathlines.velocities = fullfile(pdir, 'velocities') ;
    pvdir = QS.dir.pathlines.velocities ;
    QS.fileName.pathlines.velocities.v3d = fullfile(pvdir, 'vM.mat')  ;
    QS.fileName.pathlines.velocities.vn  = fullfile(pvdir, 'vnM.mat') ;
    QS.fileName.pathlines.velocities.vv  = fullfile(pvdir, 'vvM.mat') ;
    QS.fileName.pathlines.velocities.vf  = fullfile(pvdir, 'vfM.mat') ;
    QS.fileName.pathlines.velocities.v2d = fullfile(pvdir, 'v2dM.mat') ;
    QS.fileName.pathlines.velocities.v2dum = fullfile(pvdir, 'v2dMum.mat') ;
    QS.fileName.pathlines.velocities.v3dsm   = fullfile(pvdir, 'vsmM.mat') ;
    QS.fileName.pathlines.velocities.vnsm    = fullfile(pvdir, 'vnsmM.mat') ;
    QS.fileName.pathlines.velocities.vvsm    = fullfile(pvdir, 'vvsmM.mat') ;
    QS.fileName.pathlines.velocities.vfsm    = fullfile(pvdir, 'vfsmM.mat') ;
    QS.fileName.pathlines.velocities.v2dsm   = fullfile(pvdir, 'v2dsmM.mat') ;
    QS.fileName.pathlines.velocities.v2dsmum = fullfile(pvdir, 'v2dsmMum.mat') ;

    % pathline fileBases
    QS.fileBase.pathlines = struct() ;
    QS.fileBase.pathlines.strain = fullfile(QS.dir.pathlines.strain, QS.fileBase.strain) ;
    QS.fileBase.pathlines.dxdyFiltered = fullfile(QS.dir.pathlines.dxdyFiltered, ['dxdy_' QS.fileBase.strain]) ;


    %% UVPrime pathlines -- alternative to spsm piv, not very conformal, so not very useful
    % By default, we also include uvp_sme as a coordinate system for
    % quasiconformal measurements --> note that these are less Lagrangian than
    % sp_sme or up_sme, so they are treated as independent from the main
    % pipeline in which we measure velocities
    QS.dir.piv_uvprime = fullfile(uvDir, 'piv_uvprime_e') ; 

    QS.dir.pathlines_uvprime = struct() ;
    QS.dir.pathlines_uvprime.data = fullfile(QS.dir.piv_uvprime, 'pathlines', 't0_%04d') ; 
    pdir = QS.dir.pathlines_uvprime.data ;
    QS.dir.pathlines_uvprime.XY = fullfile(pdir, 'images_uvprime_XY') ;
    QS.dir.pathlines_uvprime.XYZ = fullfile(pdir, 'images_uvprime_XYZ') ;
    QS.dir.pathlines_uvprime.vXY = fullfile(pdir, 'images_uvprime_vXY') ;
    QS.dir.pathlines_uvprime.fXY = fullfile(pdir, 'images_uvprime_fXY') ;
    QS.dir.pathlines_uvprime.v3d = fullfile(pdir, 'images_uvprime_v3d') ;
    QS.dir.pathlines_uvprime.f3d = fullfile(pdir, 'images_uvprime_f3d') ;
    QS.dir.pathlines_uvprime.strain_images = fullfile(pdir, 'images_strain') ;
    QS.dir.pathlines_uvprime.strain = fullfile(pdir, 'strain') ;
    QS.dir.pathlines_uvprime.quasiconformal = fullfile(pdir, 'quasiconformal') ;
    QS.dir.pathlines_uvprime.radius = fullfile(pdir, 'images_radii_vertices') ;
    QS.dir.pathlines_uvprime.indentation = fullfile(pdir, 'images_indentation') ;
    QS.dir.pathlines_uvprime.kymographs = fullfile(pdir, 'kymographs') ;
    QS.fileName.pathlines_uvprime = struct() ;
    QS.fileName.pathlines_uvprime.featureIDs = fullfile(pdir, 'featureIDs.mat') ;
    QS.fileName.pathlines_uvprime.XY = fullfile(pdir, 'piv_pathlines_uvprime_XY.mat') ;
    QS.fileName.pathlines_uvprime.XYZ = fullfile(pdir, 'piv_pathlines_uvprime_XYZ.mat') ;
    QS.fileName.pathlines_uvprime.vXY = fullfile(pdir, 'piv_pathlines_uvprime_vXY.mat') ;
    QS.fileName.pathlines_uvprime.v3d = fullfile(pdir, 'piv_pathlines_uvprime_v3d.mat') ;
    QS.fileName.pathlines_uvprime.fXY = fullfile(pdir, 'piv_pathlines_uvprime_fXY.mat') ;
    QS.fileName.pathlines_uvprime.f3 = fullfile(pdir, 'piv_pathlines_uvprime_f3d.mat') ;
    QS.fileName.pathlines_uvprime.refMesh = fullfile(pdir, 'refMesh.mat') ;
    QS.fileName.pathlines_uvprime.quasiconformal = ...
        fullfile(QS.dir.pathlines_uvprime.quasiconformal, 'mu_v3dv2d.mat') ;
    QS.fileName.pathlines_uvprime.radius = fullfile(pdir, 'pathline_radii.mat') ;
    QS.fileName.pathlines_uvprime.indentation = fullfile(pdir, 'pathline_indentation.mat') ;
    QS.fileName.pathlines_uvprime.kymographs = struct() ;
    QS.fileName.pathlines_uvprime.kymographs.radius = ...
        fullfile(QS.dir.pathlines.kymographs, 'radiusKymographs.mat') ;
    QS.fileName.pathlines_uvprime.kymographs.indentation = ...
        fullfile(QS.dir.pathlines.kymographs, 'indentationKymographs.mat') ;
    QS.fileName.pathlines_uvprime.kymographs.mu = ...
        fullfile(QS.dir.pathlines.kymographs, 'muKymographs.mat') ;


    %% Velocities -- Lagrangian averaging
    QS.dir.piv.avg = fullfile(QS.dir.piv.root, 'lagrangianAvg') ;
    % QS.dir.piv.avgCurl = fullfile(QS.dir.piv.avg, 'curl') ;
    % QS.dir.piv.avgDvg = fullfile(QS.dir.piv.avg, 'dvg') ;
    QS.fileName.pivAvg = struct() ;
    QS.fileName.pivAvg.v2dum = fullfile(QS.dir.piv.avg, 'v2dMum_avg.mat') ;
    QS.fileName.pivAvg.v2d = fullfile(QS.dir.piv.avg, 'v2dM_avg.mat') ;
    QS.fileName.pivAvg.vn  = fullfile(QS.dir.piv.avg, 'vnM_avg.mat') ;
    QS.fileName.pivAvg.v3d = fullfile(QS.dir.piv.avg, 'vM_avg.mat')  ;
    QS.fileName.pivAvg.vv  = fullfile(QS.dir.piv.avg, 'vvM_avg.mat') ;
    QS.fileName.pivAvg.vf  = fullfile(QS.dir.piv.avg, 'vfM_avg.mat') ;
    % Helmholtz-Hodge and DEC -- Lagrangian averaging
    QS.dir.piv.avgDEC = struct() ;
    QS.dir.piv.avgDEC.data   = fullfile(QS.dir.piv.avg, 'dec') ;
    QS.dir.piv.avgDEC.lap2d   = fullfile(QS.dir.piv.avg, 'dec_lap2d') ;
    QS.dir.piv.avgDEC.lap3d   = fullfile(QS.dir.piv.avg, 'dec_lap3d') ;
    QS.dir.piv.avgDEC.lap3dTexture = fullfile(QS.dir.piv.avg, 'dec_lap3dTexture') ;
    QS.dir.piv.avgDEC.div2d  = fullfile(QS.dir.piv.avg, 'dec_div2d') ;
    QS.dir.piv.avgDEC.div3d  = fullfile(QS.dir.piv.avg, 'dec_div3d') ;
    QS.dir.piv.avgDEC.div3dTexture = fullfile(QS.dir.piv.avg, 'dec_div3dTexture') ;
    QS.dir.piv.avgDEC.rot2d  = fullfile(QS.dir.piv.avg, 'dec_rot2d') ;
    QS.dir.piv.avgDEC.rot3d  = fullfile(QS.dir.piv.avg, 'dec_rot3d') ;
    QS.dir.piv.avgDEC.rot3dTexture = fullfile(QS.dir.piv.avg, 'dec_rot3dTexture') ;
    QS.dir.piv.avgDEC.harm2d = fullfile(QS.dir.piv.avg, 'dec_harm2d') ;
    QS.dir.piv.avgDEC.harm3d = fullfile(QS.dir.piv.avg, 'dec_harm3d') ;
    QS.fullFileBase.decAvg = fullfile(QS.dir.piv.avgDEC.data, ...
                                  [QS.fileBase.name '_dec.mat'] ) ;

    %% Velocities -- no averaging
    QS.dir.pivRaw = fullfile(QS.dir.piv.root, 'noAvg') ;
    QS.fileName.pivRaw.v2dum = fullfile(QS.dir.pivRaw, 'v2dMum_simpletimeavg.mat') ;
    QS.fileName.pivRaw.v2d = fullfile(QS.dir.pivRaw, 'v2dM_simpletimeavg.mat') ;
    QS.fileName.pivRaw.vn  = fullfile(QS.dir.pivRaw, 'vnM_simpletimeavg.mat') ;
    QS.fileName.pivRaw.v3d = fullfile(QS.dir.pivRaw, 'vM_simpletimeavg.mat') ;
    QS.fileName.pivRaw.vv  = fullfile(QS.dir.pivRaw, 'vvM_simpletimeavg.mat') ;
    QS.fileName.pivRaw.vf  = fullfile(QS.dir.pivRaw, 'vfM_simpletimeavg.mat');
    % Helmholtz-Hodge and DEC -- simple/surface-Lagrangian averaging
    QS.dir.pivRawDEC = struct() ;
    QS.dir.pivRawDEC.data   = fullfile(QS.dir.pivRaw, 'dec') ;
    QS.dir.pivRawDEC.lap2d   = fullfile(QS.dir.pivRaw, 'dec_lap2d') ;
    QS.dir.pivRawDEC.lap3d   = fullfile(QS.dir.pivRaw, 'dec_lap3d') ;
    QS.dir.pivRawDEC.lap3dTexture = fullfile(QS.dir.pivRaw, 'dec_lap3dTexture') ;
    QS.dir.pivRawDEC.div2d  = fullfile(QS.dir.pivRaw, 'dec_div2d') ;
    QS.dir.pivRawDEC.div3d  = fullfile(QS.dir.pivRaw, 'dec_div3d') ;
    QS.dir.pivRawDEC.div3dTexture = fullfile(QS.dir.pivRaw, 'dec_div3dTexture') ;
    QS.dir.pivRawDEC.rot2d  = fullfile(QS.dir.pivRaw, 'dec_rot2d') ;
    QS.dir.pivRawDEC.rot3d  = fullfile(QS.dir.pivRaw, 'dec_rot3d') ;
    QS.dir.pivRawDEC.rot3dTexture = fullfile(QS.dir.pivRaw, 'dec_rot3dTexture') ;
    QS.dir.pivRawDEC.harm2d = fullfile(QS.dir.pivRaw, 'dec_harm2d') ;
    QS.dir.pivRawDEC.harm3d = fullfile(QS.dir.pivRaw, 'dec_harm3d') ;
    QS.fullFileBase.decRaw = fullfile(QS.dir.pivRawDEC.data, ...
                                  [QS.fileBase.name '_dec.mat'] ) ;

    %% Velocities -- simple/surface-Lagrangian averaging
    QS.dir.pivSimAvg = fullfile(QS.dir.piv.root, 'simpleAvg') ;
    % QS.dir.pivSimAvgCurl = fullfile(QS.dir.pivSimAvg, 'rot') ;
    % QS.dir.pivSimAvgDvg = fullfile(QS.dir.pivSimAvg, 'dvg') ;
    QS.fileName.pivSimAvg = struct() ;
    QS.fileName.pivSimAvg.v2dum = fullfile(QS.dir.pivSimAvg, 'v2dMum_simpletimeavg.mat') ;
    QS.fileName.pivSimAvg.v2d = fullfile(QS.dir.pivSimAvg, 'v2dM_simpletimeavg.mat') ;
    QS.fileName.pivSimAvg.vn  = fullfile(QS.dir.pivSimAvg, 'vnM_simpletimeavg.mat') ;
    QS.fileName.pivSimAvg.v3d = fullfile(QS.dir.pivSimAvg, 'vM_simpletimeavg.mat') ;
    QS.fileName.pivSimAvg.vv  = fullfile(QS.dir.pivSimAvg, 'vvM_simpletimeavg.mat') ;
    QS.fileName.pivSimAvg.vf  = fullfile(QS.dir.pivSimAvg, 'vfM_simpletimeavg.mat');
    % Helmholtz-Hodge and DEC -- simple/surface-Lagrangian averaging
    QS.dir.pivSimAvgDEC = struct() ;
    QS.dir.pivSimAvgDEC.data   = fullfile(QS.dir.pivSimAvg, 'dec') ;
    QS.dir.pivSimAvgDEC.lap2d   = fullfile(QS.dir.pivRaw, 'dec_lap2d') ;
    QS.dir.pivSimAvgDEC.lap3d   = fullfile(QS.dir.pivRaw, 'dec_lap3d') ;
    QS.dir.pivSimAvgDEC.lap3dTexture = fullfile(QS.dir.pivRaw, 'dec_lap3dTexture') ;
    QS.dir.pivSimAvgDEC.div2d  = fullfile(QS.dir.pivSimAvg, 'dec_div2d') ;
    QS.dir.pivSimAvgDEC.div3d  = fullfile(QS.dir.pivSimAvg, 'dec_div3d') ;
    QS.dir.pivSimAvgDEC.div3dTexture = fullfile(QS.dir.pivSimAvg, 'dec_div3dTexture') ;
    QS.dir.pivSimAvgDEC.rot2d  = fullfile(QS.dir.pivSimAvg, 'dec_rot2d') ;
    QS.dir.pivSimAvgDEC.rot3d  = fullfile(QS.dir.pivSimAvg, 'dec_rot3d') ;
    QS.dir.pivSimAvgDEC.rot3dTexture = fullfile(QS.dir.pivSimAvg, 'dec_rot3dTexture') ;
    QS.dir.pivSimAvgDEC.harm2d = fullfile(QS.dir.pivSimAvg, 'dec_harm2d') ;
    QS.dir.pivSimAvgDEC.harm3d = fullfile(QS.dir.pivSimAvg, 'dec_harm3d') ;
    QS.fullFileBase.decSimAvg = fullfile(QS.dir.pivSimAvgDEC.data, ...
                                  [QS.fileBase.name '_dec.mat'] ) ;

    %% Double resolution
    QS.dir.piv3d2x = fullfile(QS.dir.piv.root, 'piv3dDoubleRes') ;
    QS.dir.pivt2d2x = fullfile(QS.dir.piv.root, 'vt2dDoubleRes') ;
    QS.dir.pivn2d2x = fullfile(QS.dir.piv.root, 'vn2dDoubleRes') ;
    QS.dir.pivdilation2x = fullfile(QS.dir.piv.root, 'dilationDoubleRes') ;
    QS.fullFileBase.piv3d2x = fullfile(QS.dir.piv3d2x, 'piv3dDoubleRes_%04d.mat') ;
    %% 2x Velocities -- Lagrangian averaging
    QS.dir.piv.avg2x = fullfile(QS.dir.piv.root, 'lagrangianAvgDoubleRes') ;
    QS.fileName.pivAvg2x = struct() ;
    QS.fileName.pivAvg2x.v2dum = fullfile(QS.dir.piv.avg2x, 'v2dMum_timeavg2x.mat') ;
    QS.fileName.pivAvg2x.v2d = fullfile(QS.dir.piv.avg2x, 'v2dM_timeavg2x.mat') ;
    QS.fileName.pivAvg2x.vn  = fullfile(QS.dir.piv.avg2x, 'vnM_timeavg2x.mat') ;
    QS.fileName.pivAvg2x.v3d = fullfile(QS.dir.piv.avg2x, 'vM_timeavg2x.mat') ;
    QS.fileName.pivAvg2x.vv  = fullfile(QS.dir.piv.avg2x, 'vvM_timeavg2x.mat') ;
    QS.fileName.pivAvg2x.vf  = fullfile(QS.dir.piv.avg2x, 'vfM_timeavg2x.mat') ;
    QS.dir.piv.avgDEC2x = struct() ;
    QS.dir.piv.avgDEC2x.data   = fullfile(QS.dir.piv.avg2x, 'dec') ;
    QS.dir.piv.avgDEC2x.div2d  = fullfile(QS.dir.piv.avg2x, 'dec_div2d') ;
    QS.dir.piv.avgDEC2x.div3d  = fullfile(QS.dir.piv.avg2x, 'dec_div3d') ;
    QS.dir.piv.avgDEC2x.div3dTexture = fullfile(QS.dir.piv.avg2x, 'dec_div3dTexture') ;
    QS.dir.piv.avgDEC2x.rot2d  = fullfile(QS.dir.piv.avg2x, 'dec_rot2d') ;
    QS.dir.piv.avgDEC2x.rot3d  = fullfile(QS.dir.piv.avg2x, 'dec_rot3d') ;
    QS.dir.piv.avgDEC2x.rot3dTexture = fullfile(QS.dir.piv.avg2x, 'dec_rot3dTexture') ;
    QS.dir.piv.avgDEC2x.harm2d = fullfile(QS.dir.piv.avg2x, 'dec_harm2d') ;
    QS.dir.piv.avgDEC2x.harm3d = fullfile(QS.dir.piv.avg2x, 'dec_harm3d') ;
    QS.fullFileBase.decAvg2x = fullfile(QS.dir.piv.avgDEC2x.data, ...
                                  [QS.fileBase.name '_dec.mat'] ) ;
    %% 2x Velocities -- simple/surface-Lagrangian averaging
    QS.dir.pivSimAvg2x = fullfile(QS.dir.piv.root, 'simpleAvgDoubleRes') ;
    QS.fileName.pivSimAvg2x = struct() ;
    QS.fileName.pivSimAvg2x.v2dum = ...
        fullfile(QS.dir.pivSimAvg2x, 'v2dMum_simpletimeavg2x.mat') ;
    QS.fileName.pivSimAvg2x.v2d = ...
        fullfile(QS.dir.pivSimAvg2x, 'v2dM_simpletimeavg2x.mat') ;
    QS.fileName.pivSimAvg2x.vn = ...
        fullfile(QS.dir.pivSimAvg2x, 'vnM_simpletimeavg2x.mat') ;
    QS.fileName.pivSimAvg2x.v3d = ...
        fullfile(QS.dir.pivSimAvg2x, 'vM_simpletimeavg2x.mat') ;
    QS.fileName.pivSimAvg2x.vv = ...
        fullfile(QS.dir.pivSimAvg2x, 'vvM_simpletimeavg2x.mat') ;
    QS.fileName.pivSimAvg2x.vf = ...
        fullfile(QS.dir.pivSimAvg2x, 'vfM_simpletimeavg2x.mat') ;
    QS.dir.pivSimAvgDEC2x = struct() ;
    QS.dir.pivSimAvgDEC2x.data   = fullfile(QS.dir.pivSimAvg2x, 'dec') ;
    QS.dir.pivSimAvgDEC2x.div2d  = fullfile(QS.dir.pivSimAvg2x, 'dec_div2d') ;
    QS.dir.pivSimAvgDEC2x.div3d  = fullfile(QS.dir.pivSimAvg2x, 'dec_div3d') ;
    QS.dir.pivSimAvgDEC2x.div3dTexture = ...
        fullfile(QS.dir.pivSimAvg2x, 'dec_div3dTexture') ;
    QS.dir.pivSimAvgDEC2x.rot2d  = fullfile(QS.dir.pivSimAvg2x, 'dec_rot2d') ;
    QS.dir.pivSimAvgDEC2x.rot3d  = fullfile(QS.dir.pivSimAvg2x, 'dec_rot3d') ;
    QS.dir.pivSimAvgDEC2x.rot3dTexture = ...
        fullfile(QS.dir.pivSimAvg2x, 'dec_rot3dTexture') ;
    QS.dir.pivSimAvgDEC2x.harm2d = fullfile(QS.dir.pivSimAvg2x, 'dec_harm2d') ;
    QS.dir.pivSimAvgDEC2x.harm3d = fullfile(QS.dir.pivSimAvg2x, 'dec_harm3d') ;
    QS.fullFileBase.decSimAvg2x  = fullfile(QS.dir.pivSimAvgDEC2x.data, ...
                                  [QS.fileBase.name '_dec.mat'] ) ;

    %% Eulerian kinematics
    QS.dir.eulerianKinematics = fullfile(QS.dir.uvCoord, 'eulerianKinematics') ;
end

%% Ensure directories
dirs2make = struct2cell(QS.dir) ;
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
                    if ~exist(dir2makeSub, 'dir') && ~contains(dir2makeSub, '%04d') ...
                            && ~contains(dir2makeSub, '%0.3f') ...
                            && ~contains(dir2makeSub, '%06d') ...
                            && ~contains(dir2makeSub, '%s')
                        mkdir(dir2makeSub)
                    end
                end
            else
                if ~exist(dir2make, 'dir') && ~contains(dir2make, '%04d') ...
                        && ~contains(dir2make, '%0.3f') ...
                        && ~contains(dir2make, '%06d') ...
                        && ~contains(dir2make, '%s')
                    mkdir(dir2make)
                end
            end
        end
    else
        if ~exist(dir2make, 'dir')
            mkdir(dir2make)
        end
    end
end