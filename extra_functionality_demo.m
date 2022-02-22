
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 5: Extra functionality to demonstrate, just for fun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% OPTIONAL: Make 3D pullback stacks, rather than MIPs ====================
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

%% Measure coarse-grained bond contraction and dilation in zeta=s/L, phi
% todo: consider putting this in Ricci map frame instead of (s,phi) frame
options = struct() ;
tubi.measureDxDyStrainFiltered(options) ;