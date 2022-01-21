function measureCurvatures(QS, options) 
%measureCurvatures(QS, options) 
%   Measure and plot mean and gaussian curvature for smoothed meshes
%   (QS.currentMesh.spcutMeshSm) for each timepoint
% 
% Parameters
% ----------
% QS : QuapSlap class instance
%   Current QuapSlap object
% options : struct with fields
%   overwrite : bool, default=false
%       overwrite results on disk
% 
% NPMitchell 2020


%% Unpack QS
KSmDir = QS.dir.gaussCurvature ;
HSmDir = QS.dir.meanCurvature ;
KHSmDir = QS.dir.curvatures ;
[~,~,~,xyzlim] = QS.getXYZLims() ;

%% Unpack options
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end

%% Ensure directories exist for images
dirs2do = {fullfile(KSmDir, 'latleft'), ...
    fullfile(KSmDir, 'dorsal'), ...
    fullfile(KSmDir, 'latright'), ...
    fullfile(KSmDir, 'ventral'), ...
    fullfile(HSmDir, 'latleft'), ...
    fullfile(HSmDir, 'dorsal'), ...
    fullfile(HSmDir, 'latright'), ...
    fullfile(HSmDir, 'ventral')} ;

for qq=1:length(dirs2do)
    mkdir(dirs2do{qq})
end

%% Compute the curvatures and save images
getderivatives=0;
disp('Computing/Loading Curvatures...')
for tt = QS.xp.fileMeta.timePoints
    Kfns = {fullfile(KSmDir, 'latleft', sprintf('gausscurv_latleft_%06d.png', tt)), ...
        fullfile(KSmDir, 'dorsal', sprintf('gausscurv_dorsal_%06d.png', tt)), ...
        fullfile(KSmDir, 'latright', sprintf('gausscurv_latright_%06d.png', tt)), ...
        fullfile(KSmDir, 'ventral', sprintf('gausscurv_ventral_%06d.png', tt)) };
    Hfns = {fullfile(HSmDir, 'latleft', sprintf('meancurv_latleft_%06d.png', tt)), ...
        fullfile(HSmDir, 'dorsal', sprintf('meancurv_dorsal_%06d.png', tt)), ...
        fullfile(HSmDir, 'latright', sprintf('meancurv_latright_%06d.png', tt)),...
        fullfile(HSmDir, 'ventral', sprintf('meancurv_ventral_%06d.png', tt))} ;
    
    KHfn = sprintf(QS.fullFileBase.curvatures, tt) ;
    % Note this used to say:
    % KHfn = fullfile(KHSmDir, sprintf('gauss_mean_curvature_%06d.mat', tt)) ;
    
    if ~exist(KHfn, 'file') || overwrite
        % load glued / closed mesh for curvature computation
        load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tt), 'spcutMeshSmRSC') ;
        fv.faces = spcutMeshSmRSC.f ;
        fv.vertices = spcutMeshSmRSC.v ;
        PrincipalCurvatures = GetCurvatures( fv, getderivatives);
        gaussCurv = (PrincipalCurvatures(1,:).*PrincipalCurvatures(2,:))';
        meanCurv = -0.5 * (PrincipalCurvatures(1,:) + PrincipalCurvatures(2,:))';

        opts.outfn = Kfns ;
        opts.label = 'Gaussian curvature, $K$' ;
        opts.sscale = 0.01 ;
        opts.cbarPosition = [.85 .333 .02 .333] ;
        opts.xlim = xyzlim(1, :) ;
        opts.ylim = xyzlim(2, :) ;
        opts.zlim = xyzlim(3, :) ;
        scalarFieldOnSurface(fv.faces, fv.vertices, gaussCurv, opts)
        opts.outfn = Hfns ;
        opts.sscale = 0.3 ;
        opts.label = 'mean curvature, $H$' ;
        scalarFieldOnSurface(fv.faces, fv.vertices, meanCurv, opts)
        % Save the curvatures as a mat file
        save(KHfn, 'meanCurv', 'gaussCurv')
    else
        disp('Curvatures already computed & saved on disk')
    end
end