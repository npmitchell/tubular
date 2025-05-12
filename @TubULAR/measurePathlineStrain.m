function measurePathlineStrain(tubi, options)
% measurePathlineStrain(tubi, options)
%   Compute strain from integrated pathlines deforming mesh vertices.
%   Measurements are taken with respect to fixed Lagrangian frame. 
%   Plot results in 2d and/or 3d for each timepoint.
%   
%
% Parameters
% ----------
% tubi : QuapSlap class instance
% options : struct with fields
%   overwrite : bool, overwrite data on disk
%   overwriteImages : bool, overwrite images of results on disk
%   plot_comparison : bool, plot a comparison with DEC traceful dilation
%   median_filter_strainRates : bool
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
overwriteImages = false ;
plot_dzdp = false ;
plot_comparison = false ;
median_filter_strainRates = false ;
climitInitial = 0.05 ;
climitRamp = 0.005 ;
climitRatio = 1 ;
preview = false ;
notClosedTube = 0; 

%% Parameter options
lambda_mesh = 0.002 ;
lambda = 0.01 ; 
debug = false ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = "Lagrangian" ;
LagrangianFrame = 'zetaphi' ;
% Load time offset for first fold, t0 -- default pathline t0
tubi.t0set() ;
t0 = tubi.t0 ;
% By default, t0Pathline = t0 (see below)

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'plot_comparison')
    plot_comparison = options.plot_comparison ;
end
if isfield(options, 'plot_dzdp')
    plot_dzdp = options.plot_dzdp ;
end

%% parameter options
if isfield(options, 'climitInitial')
    climitInitial = options.climitInitial ;
end
if isfield(options, 'climitRamp')
    climitRamp = options.climitRamp ;
end
if isfield(options, 'climitRatio')
    climitRatio = options.climitRatio ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
else
    t0Pathline = t0 ;
end
if isfield(options, 'notClosedTube')
    notClosedTube = options.notClosedTube; 
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
else
    error("Could not parse samplingResolution: set to '1x'")
end

%% Unpack tubi
tubi.getXYZLims ;
xyzlim = tubi.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
fons = tubi.t0set() - tubi.xp.fileMeta.timePoints(1) ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
bbr256 = blueblackred(256) ;
close all

%% load from tubi
if doubleResolution
    nU = tubi.nU * 2 - 1 ;
    nV = tubi.nV * 2 - 1 ;
else
    nU = tubi.nU ;
    nV = tubi.nV ;    
end

% We relate the normal velocities to the divergence / 2 * H.
tps = tubi.xp.fileMeta.timePoints(1:end-1) - t0;

% Unit definitions for axis labels
unitstr = [ '[1/' tubi.timeUnits ']' ];
vunitstr = [ '[' tubi.spaceUnits '/' tubi.timeUnits ']' ];
    
% DONE WITH PREPARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pathlines to build metric from advected mesh faces
tubi.loadPullbackPathlines(t0Pathline, 'vertexPathlines')
vP = tubi.pathlines.vertices ;

% Output directory is inside pathline dir
outdir = sprintf(tubi.dir.pathlines.strain, t0) ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Load Lx, Ly by loadingPIV. 
tubi.loadPIV()
Xpiv = tubi.piv.raw.x ;
Ypiv = tubi.piv.raw.y ;

% Also need velocities to advect mesh
% tubi.loadVelocityAverage('vv')
% vPIV = tubi.velocityAverage.vv ;

% Discern if piv measurements are done on a double covering or the meshes
if strcmp(tubi.piv.imCoords(end), 'e')
    doubleCovered = true ;
end

%% INTEGRATE STRAINRATE INTO STRAIN ON PATHLINES
% Compute or load all timepoints
pb = tubi.getPullbackPathlines([], 'vertexPathlines3d') ;
% load(sprintfm(tubi.fileName.pathlines.v3d, t0), 'v3dPathlines')
% load(sprintfm(tubi.fileName.pathlines.vXY, t0), 'vertexPathlines')
vX3rs = pb.vertices3d.vXrs ;
vY3rs = pb.vertices3d.vYrs ;
vZ3rs = pb.vertices3d.vZrs ;
t0 = pb.vertices3d.t0 ;
tIdx0 = pb.vertices3d.tIdx0 ;

if ~all(isfinite(vX3rs(:)))
    error('Some vertex 3D positions are infinite. Check interpolation')
end

% Define reference mesh
refMeshFn = fullfile(sprintfm(tubi.dir.pathlines.data, t0Pathline), ...
        'refMesh.mat') ;
if exist(refMeshFn, 'file') || overwrite
    load(refMeshFn, 'refMesh')
else
    refMesh = struct() ; 
    vX = vertexPathlines.vX(tIdx0, :) ; 
    vY = vertexPathlines.vY(tIdx0, :) ;
    vXY = [vX(:), vY(:)] ;
    Lx = vertexPathlines.Lx ;
    Ly = vertexPathlines.Ly ;
    refMesh.f = defineFacesRectilinearGrid(vXY, tubi.nU, tubi.nV) ;
    refMesh.u = tubi.XY2uv([Lx(tIdx0), Ly(tIdx0)], vXY, 1, 1) ;
    x0 = vX3rs(tIdx0, :) ;
    y0 = vY3rs(tIdx0, :) ;
    z0 = vZ3rs(tIdx0, :) ;
    refMesh.v = [ x0(:), y0(:), z0(:) ] ; 
    pathPairs = [ (1:nU)', (nV-1)*nU + (1:nU)' ] ;
    refMesh.pathPairs = pathPairs ;
    refMesh.nU = tubi.nU ;
    refMesh.nV = tubi.nV ;
    assert(numel(refMesh.u(:, 1)) == tubi.nU * tubi.nV)

    % Save reference Mesh as lagrangian frame
    save(refMeshFn, 'refMesh') ;
end

% Set the rotated scaled vertices as field "v" for
% inducedStrainPeriodicMesh() call
refMesh.v = refMesh.vrs ;

ntps = length(tubi.xp.fileMeta.timePoints(1:end-1)) ;
tidx2do = 1:40:ntps ;
tidx0 = tubi.xp.tIdx(tubi.t0set()) ;
tidx2do = [tidx2do setdiff(tidx0:10:ntps, tidx2do)] ;
tidx2do = [tidx2do setdiff(1:ntps, tidx2do)] ;
for tidx = tidx2do
    % Identify current timepoint
    tp = tubi.xp.fileMeta.timePoints(tidx) ;
    
    % Do the fund forms for the lagrangian frame already exist?
    % ffn = sprintfm(tubi.fullFileBase.pathline.fundForms, t0, tp) ;
    
    ffn = sprintfm(tubi.fullFileBase.pathlines.strain, t0, tp) ;
    if ~exist(ffn, 'file') || overwrite 
        disp(['t = ' num2str(tp)])
        tubi.setTime(tp) ;

        % Load Lagrangian advected vertices as mesh for this timepoint
        % load(sprintfm(tubi.fullFileBase.spcutMeshSmRS, tp), 'spcutMeshSmRS') ;
        xx = vX3rs(tidx, :) ;
        yy = vY3rs(tidx, :) ;
        zz = vZ3rs(tidx, :) ;
        v3d = [ xx(:), yy(:), zz(:) ] ;
        mesh = struct() ;
        mesh.f = refMesh.f ;
        mesh.v = v3d ;
        mesh.u = refMesh.u ;
        mesh.pathPairs = refMesh.pathPairs ;
        mesh.nU = tubi.nU ;
        mesh.nV = tubi.nV ;
        
        [strain, tre, dev, theta, outputStruct] = ...
            inducedStrainPeriodicMesh(mesh, refMesh, options) ;
        
        %   fundForms : struct with fields
        %       gg : first fundamental form on each face for tiled mesh
        %       bb : second fundamental form on each face for tiled mesh
        %   bondDxDy : struct with fields
        %       dx : length of dx in embedding space / length of dx in pullback space
        %       dy : length of dy in embedding space / length of dy in pullback space
        %       dxTiled : length of dx in embedding space / length of dx in
        %           pullback space for tiled mesh
        %       dyTiled : length of dy in embedding space / length of dy in
        %           pullback space for tiled mesh
        %   theta_pb : #tiledFaces x 1 float array
        %       elongation angle, theta, in pullback space (differs from
        %       theta by scaling of the eigenvectors by (dx, dy) returned in
        %       bondDxDy
        %   faceIDs : #faces x 1 int array
        %       indices of tiled faces that return the faces of the
        %       original input mesh, so that strain_orig = strain(faceIDs)
        fundForms = outputStruct.fundForms ;
        bondDxDy = outputStruct.bondDxDy ;
        % The following are already stored in bondDxDy:
        %  dx_strain = outputStruct.dx_strain ;
        %  dy_strain = outputStruct.dy_strain ;
        %  dbonds_ref = outputStruct.dbonds_ref ;
        %  dbonds_def = outputStruct.dbonds_def ;
        theta_pb = outputStruct.theta_pb ;
        faceIDs = outputStruct.faceIDs ;
        
        % Ratio of initial to 
        fractionalAreaChange = (doublearea(mesh.v, mesh.f) - ...
            doublearea(refMesh.v, refMesh.f)) ./ doublearea(refMesh.v, refMesh.f) ;
        
        
        % SMOOTHED version
        e11 = squeeze(fundForms.g1(faceIDs, 1, 1)-fundForms.g0(faceIDs, 1, 1)) ;
        e12 = squeeze(fundForms.g1(faceIDs, 1, 2)-fundForms.g0(faceIDs, 1, 2)) ;
        e21 = squeeze(fundForms.g1(faceIDs, 2, 1)-fundForms.g0(faceIDs, 2, 1)) ;
        e22 = squeeze(fundForms.g1(faceIDs, 2, 2)-fundForms.g0(faceIDs, 2, 2)) ;

        % Smooth all four components with laplacian
        [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        
        % Force no imaginary values in V2F and F2V
        if notClosedTube
            smoothing_options = struct('widthX', 3, 'nmodes', 0) ; 
            e11 = F2V * e11; e12 = F2V * e12; e21 = F2V * e21; e22 = F2V * e22; 
            dA0 = F2V * fractionalAreaChange; 
        else
            F2V = real(F2V) ;
            smoothing_options = struct('widthX', 3, 'nmodes', 5) ; 
            e11 = modeFilterQuasi1D(reshape(F2V * e11, [nU, nV]), smoothing_options) ;
            e12 = modeFilterQuasi1D(reshape(F2V * e12, [nU, nV]), smoothing_options) ;
            e21 = modeFilterQuasi1D(reshape(F2V * e21, [nU, nV]), smoothing_options) ;
            e22 = modeFilterQuasi1D(reshape(F2V * e22, [nU, nV]), smoothing_options) ;
            dA0 = modeFilterQuasi1D(reshape(F2V * fractionalAreaChange, [nU, nV]), smoothing_options) ;
        end
        
        % e11 = laplacian_smooth(mesh.v,mesh.f,'cotan',[], lambda,'explicit',F2V*e11,1000) ;
        % e12 = laplacian_smooth(mesh.v,mesh.f,'cotan',[], lambda,'explicit',F2V*e12,1000) ;
        % e21 = laplacian_smooth(mesh.v,mesh.f,'cotan',[], lambda,'explicit',F2V*e21,1000) ;
        % e22 = laplacian_smooth(mesh.v,mesh.f,'cotan',[], lambda,'explicit',F2V*e22,1000) ;
        strain_sm_vertices = zeros(numel(e11(:)), 2, 2) ;
        strain_sm_vertices(:, 1, 1) = e11(:) ;
        strain_sm_vertices(:, 1, 2) = e12(:) ;
        strain_sm_vertices(:, 2, 1) = e21(:) ;
        strain_sm_vertices(:, 2, 2) = e22(:) ;
        fractionalAreaChange_sm_vertices = dA0 ;
        
        % % Check result
        % subplot(2, 2, 1)
        % trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
        %     'facevertexcdata', strain_sm_vertices(:, 1, 1)); axis equal; view([0,0]) ;
        % subplot(2, 2, 2)
        % trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
        %     'facevertexcdata', strain_sm_vertices(:, 1, 2)); axis equal; view([0,0]) ;
        % subplot(2, 2, 3)
        % trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
        %     'facevertexcdata', strain_sm_vertices(:, 2, 1)); axis equal; view([0,0]) ;
        % subplot(2, 2, 4)
        % trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
        %     'facevertexcdata', strain_sm_vertices(:, 2, 2)); axis equal; view([0,0]) ;
        % sgtitle(['t=' num2str(tp)])
        % pause(0.1)
        
        % Smooth the strain tensors
        strain_sm = zeros(numel(faceIDs), 2, 2) ;
        strain_sm(:, 1, 1) = V2F*e11(:) ;
        strain_sm(:, 1, 2) = V2F*e12(:) ;
        strain_sm(:, 2, 1) = V2F*e21(:) ;
        strain_sm(:, 2, 2) = V2F*e22(:) ;
        fractionalAreaChange_sm = V2F*fractionalAreaChange_sm_vertices(:) ; 

        tre_sm = zeros(numel(faceIDs), 1) ;
        dev_sm = zeros(numel(faceIDs), 1) ;
        theta_sm = zeros(numel(faceIDs), 1) ;
        for qq = 1:numel(faceIDs)
            fid = faceIDs(qq) ;
            eq = squeeze(strain_sm(qq, :, :)) ;
            gq = squeeze(fundForms.g0(fid, :, :)) ;
            dxTiled_ii = bondDxDy.dxTiled(fid) ;
            dyTiled_ii = bondDxDy.dyTiled(fid) ;
            try
                [tre_sm(qq), dev_sm(qq), theta_sm(qq), theta_pb_sm(qq)] = ...
                    traceDeviatorPullback(eq, gq, dxTiled_ii, dyTiled_ii) ;
            catch
                tre_sm(qq) = NaN ;
                dev_sm(qq) = NaN ;
                theta_sm(qq) = NaN ;
                theta_pb_sm(qq) = NaN ;
            end
        end
        
        % Check result
        if preview
            subplot(2, 2, 1)
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
                'facevertexcdata', tre_sm) ; axis equal; view([0,0]) ; colorbar
            title('tre')
            subplot(2, 2, 2)
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
                'facevertexcdata', dev_sm) ; axis equal; view([0,0]) ; colorbar
            title('dev')
            subplot(2, 2, 3)
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
                'facevertexcdata', fractionalAreaChange_sm) ; axis equal; view([0,0]) ; colorbar
            title('(A-A_0)/A_0')
            subplot(2, 2, 4)
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none', ...
                'facevertexcdata', theta_sm) ; axis equal; view([0,0]) ; colorbar
            title('theta')
            pause(0.1)
        end
        
        save(ffn, 'strain', 'tre', 'dev', 'theta', 'fundForms', 'bondDxDy', ...
            'theta_pb', 'faceIDs', ...
            'fractionalAreaChange', ...
            'strain_sm', 'strain_sm_vertices', ...
            'tre_sm', 'dev_sm', 'theta_sm', 'theta_pb_sm', ...
            'fractionalAreaChange_sm', 'fractionalAreaChange_sm_vertices', ...
            'smoothing_options') ;
    else
        disp(['Strain + fundForms for t = ' num2str(tp) ' already on disk'])
    end

    %% Plot strain for this timepoint
    options.overwrite = overwriteImages ;
    tubi.plotPathlineStrainTimePoint(tp, options)
end

disp('done with integrated pathline strain calculations')


%% Combine DV-averaged gg and bb into kymographs
apKymoFn = fullfile(outdir, 'apKymographLagrangianMetric.mat') ;
lKymoFn = fullfile(outdir, 'leftKymographLagrangianMetric.mat') ;
rKymoFn = fullfile(outdir, 'rightKymographLagrangianMetric.mat') ;
dKymoFn = fullfile(outdir, 'dorsalKymographLagrangianMetric.mat') ;
vKymoFn = fullfile(outdir, 'ventralKymographLagrangianMetric.mat') ;

% Create mesh averaging operators
glueRefMesh = glueCylinderCutMeshSeam(refMesh) ;
[~, F2V] = meshAveragingOperators(glueRefMesh.f, glueRefMesh.v) ;

files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if ~files_exist || overwrite || true
    disp('Compiling kymograph data to save to disk...')
    for tp = tubi.xp.fileMeta.timePoints(1:end-1)
        tidx = tubi.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        srfn = fullfile(outdir, sprintf('strain_%06d.mat', tp))   ;

        % Load timeseries measurements
        load(srfn, 'tre', 'dev', 'theta', 'faceIDs') ;
        tre = F2V * tre(faceIDs) ;
        dev = F2V * dev(faceIDs) ;
        theta = F2V * theta(faceIDs) ;
        % Double first row as last
        tre(length(tre)+1:length(tre)+nU) = tre(1:nU) ;
        dev(length(dev)+1:length(dev)+nU) = dev(1:nU) ;
        theta(length(theta)+1:length(theta)+nU) = theta(1:nU) ;
        
        %% OPTION 1: simply reshape, tracing each XY pathline pt to its t0
        % % grid coordinate
        tre = reshape(tre, [nU, nV]) ;
        dev = reshape(dev, [nU, nV]) ;
        theta = reshape(theta, [nU, nV]) ;
        
        %% Average strainRATE along pathline DV hoops
        % Average along DV -- ignore last redudant row at nV
        [dv_ap, th_ap] = ...
            tubi.dvAverageNematic(dev(:, 1:nV-1), theta(:, 1:nV-1)) ;
        tr_ap = mean(tre(:, 1:nV-1), 2) ;
        
        % quarter bounds
        q0 = round(nV * 0.125) ;
        q1 = round(nV * 0.375) ;
        q2 = round(nV * 0.625) ;
        q3 = round(nV * 0.875) ;
        left = q0:q1 ;
        ventral = q1:q2 ;
        right = q2:q3 ;
        dorsal = [q3:nV, 1:q1] ;
        
        % left quarter
        [dv_l, th_l] = ...
            tubi.dvAverageNematic(dev(:, left), theta(:, left)) ;
        tr_l = mean(tre(:, left), 2) ;
        
        % right quarter
        [dv_r, th_r] = ...
            tubi.dvAverageNematic(dev(:, right), theta(:, right)) ;
        tr_r = mean(tre(:, right), 2) ;
        
        % dorsal quarter
        [dv_d, th_d] = ...
            tubi.dvAverageNematic(dev(:, dorsal), theta(:, dorsal)) ;
        tr_d = mean(tre(:, dorsal), 2) ;
        
        % ventral quarter
        [dv_v, th_v] = ...
            tubi.dvAverageNematic(dev(:, ventral), theta(:, ventral)) ;
        tr_v = mean(tre(:, ventral), 2) ;
        
        %% Store accumulated strain in matrices
        % dv averaged
        str_apM(tidx, :) = tr_ap ;
        sdv_apM(tidx, :) = dv_ap ;
        sth_apM(tidx, :) = th_ap ;

        % left quarter
        str_lM(tidx, :) = tr_l ;
        sdv_lM(tidx, :) = dv_l ;
        sth_lM(tidx, :) = th_l ;

        % right quarter
        str_rM(tidx, :) = tr_r ;
        sdv_rM(tidx, :) = dv_r ;
        sth_rM(tidx, :) = th_r ;

        % dorsal quarter
        str_dM(tidx, :) = tr_d ;
        sdv_dM(tidx, :) = dv_d ;
        sth_dM(tidx, :) = th_d ;

        % ventral quarter
        str_vM(tidx, :) = tr_v ;
        sdv_vM(tidx, :) = dv_v ;
        sth_vM(tidx, :) = th_v ;
    end
    
    %% Save kymographs
    disp('Saving kymograph data files for Lagrangian pathlines')
    save(apKymoFn, 'str_apM', 'sdv_apM', 'sth_apM')
    disp(['Saved kymograph data to: ' apKymoFn])
    save(lKymoFn, 'str_lM', 'sdv_lM', 'sth_lM')
    disp(['Saved kymograph data to: ' lKymoFn])
    save(rKymoFn, 'str_rM', 'sdv_rM', 'sth_rM')
    disp(['Saved kymograph data to: ' rKymoFn])
    save(dKymoFn, 'str_dM', 'sdv_dM', 'sth_dM')
    disp(['Saved kymograph data to: ' dKymoFn])
    save(vKymoFn, 'str_vM', 'sdv_vM', 'sth_vM')
    disp(['Saved kymograph data to: ' vKymoFn])
    disp('done with strain kymograph data saving')
else
    disp('strain kymograph data already on disk')    
end

