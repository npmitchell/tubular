function measureStokesForces(QS, options)
%[gP_apM, Kv_apM, Lv_apM, veln_apM] = measureStokesForces(QS, options)
%   Measure \nabla_j p / eta = (Laplace + K) v_i, meaning the laplacian of
%   the flow plus gaussian curvature times the flow give the force in the
%   tissue.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%   alphaVal : float
%       the opacity of the heatmap to overlay
%   invertImage : bool
%       invert the data pullback under the velocity map
%
% Returns 
% -------
%
% Saves to disk
% -------------
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
preview = false ;
plot_flows = true ;
plot_factors = true ;
lambda = QS.smoothing.lambda ; 
lambda_mesh = QS.smoothing.lambda_mesh ;
lambda_err = QS.smoothing.lambda_err ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;
climit = 0.1 ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = 'Lagrangian' ;

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
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'lambda_err')
    lambda_err = options.lambda_err ;
end
if isfield(options, 'nmodes')
    nmodes = options.nmodes ;
end
if isfield(options, 'zwidth')
    zwidth = options.zwidth ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
    sresStr = 'doubleRes_' ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;

if strcmp(averagingStyle, 'simple')
    sFDir = fullfile(QS.dir.stokesForcesSimple.root, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
        lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));
else
    sFDir = fullfile(QS.dir.stokesForces.root, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f_lerr%0.3f_modes%02dw%02d'], ...
        lambda, lambda_mesh, lambda_err, nmodes, zwidth), '.', 'p'));
end

folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
bwr256 = bluewhitered(256) ;

%% Load vertex-based velocity measurements
if strcmp(averagingStyle, 'Lagrangian')
    if doubleResolution
        vvsmMfn = fullfile(QS.dir.pivAvg2x, 'vvM_avg2x.mat')  ;
        tmp = load(vvsmMfn) ;
        vertex_vels = tmp.vvsmM ;
    else
        vvsmMfn = QS.fileName.pivAvg.vv ;
        % Note that this is fullfile(QS.dir.pivAvg, 'vvM_avg.mat')  
        tmp = load(vvsmMfn) ;
        vertex_vels = tmp.vvsmM ;   
    end
elseif strcmp(averagingStyle, 'simple')
    if doubleResolution
        vvsmMfn = fullfile(QS.dir.pivSimAvg2x, 'vvM_simpletimeavg2x.mat')  ;
        tmp = load(vvsmMfn) ;
        vertex_vels = tmp.vvsmM ;
    else
        vvsmMfn = fullfile(QS.dir.pivSimAvg, 'vvM_simpletimeavg.mat')  ; 
        tmp = load(vvsmMfn) ;
        vertex_vels = tmp.vvsmM ;   
    end
end

%% Load time offset for first fold, t0
QS.t0set() ;
tfold = QS.t0 ;

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

%% Test incompressibility of the flow on the evolving surface
% preallocate for cumulative error
ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
% Kv_apM   = zeros(ntps, nU) ;   % dv averaged
% Lv_apM = zeros(ntps, nU) ;
% veln_apM = zeros(ntps, nU) ;
% gP_apM = zeros(ntps, nU) ;
% Kv_lM   = zeros(ntps, nU) ;    % left averaged
% Lv_lM = zeros(ntps, nU) ;
% veln_lM = zeros(ntps, nU) ;
% gP_lM = zeros(ntps, nU) ;
% Kv_rM   = zeros(ntps, nU) ;    % right averaged
% Lv_rM = zeros(ntps, nU) ;
% veln_rM = zeros(ntps, nU) ;
% gP_rM = zeros(ntps, nU) ;
% Kv_dM   = zeros(ntps, nU) ;    % dorsal averaged
% Lv_dM = zeros(ntps, nU) ;
% veln_dM = zeros(ntps, nU) ;
% gP_dM = zeros(ntps, nU) ;
% Kv_vM   = zeros(ntps, nU) ;    % ventral averaged
% Lv_vM = zeros(ntps, nU) ;
% veln_vM = zeros(ntps, nU) ;
% gP_vM = zeros(ntps, nU) ;

% Build timepoint list so that we first do every 10, then fill in details
lastIdx = length(QS.xp.fileMeta.timePoints) - 1 ;
veryCoarseIdx = 1:50:lastIdx ;
coarseIdx = setdiff(1:10:lastIdx, veryCoarseIdx) ;
fineIdx = setdiff(1:lastIdx, [veryCoarseIdx, coarseIdx]) ;
allIdx = [170, veryCoarseIdx, coarseIdx, fineIdx ] ;
tp2do = QS.xp.fileMeta.timePoints(allIdx) ;

% Output directory is inside metricKinematics dir
outdir = fullfile(sFDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
fig2dDir = fullfile(sFDir, 'images_stokesForces2d') ;
if ~exist(fig2dDir, 'dir')
    mkdir(fig2dDir)
end
fig3dDir = fullfile(sFDir, 'images_stokesForces3d') ;
if ~exist(fig3dDir, 'dir')
    mkdir(fig3dDir)
end
    
% Compute or load all timepoints
for tp = tp2do
    close all
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Check for timepoint measurement on disk
    Kfn = fullfile(outdir, sprintf('KK_vertices_%06d.mat', tp)) ;
    Pfn = fullfile(outdir, sprintf('gradP_vertices_%06d.mat', tp)) ;
    Lfn = fullfile(outdir, sprintf('Lapv_vertices_%06d.mat', tp)) ;

    redo_comp = overwrite || ~exist(Kfn, 'file') || ~exist(Pfn, 'file') ;
    redo_comp = redo_comp || ~exist(Lfn, 'file')  ;

    if redo_comp 
        disp('computing inferred stokes forces')
        tic 
        if doubleResolution
            % Load current mesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC2x, tp)) ;
            mesh = tmp.spcutMeshSmRSC2x ;

            % Load cutMesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS2x, tp)) ;
            cutMesh = tmp.spcutMeshSmRS2x ;
            clearvars tmp
        else
            % Load current mesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
            mesh = tmp.spcutMeshSmRSC ;

            % Load cutMesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
            cutMesh = tmp.spcutMeshSmRS ;
            clearvars tmp
        end
        
        % Compute mean curvature
        % Smooth the mesh with lambda_mesh
        mesh_sm = mesh ; 
        if lambda_mesh > 0 
            disp('smoothing mesh vertices before computations')
            tri = triangulation(mesh.f, mesh.v) ;
            fbndy = tri.freeBoundary ;
            fbndy = fbndy(:, 1) ;
            mesh_sm.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda_mesh, 'implicit', mesh.v) ;
        end
        
        % Make DEC on ORIGINAL MESH, use smoothed mesh for lacplacian 
        % smoothing of fields on mesh
        DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
        K3d = load(sprintf(QS.fullFileBase.curvatures, tp), 'gaussCurv') ;
        K3d = K3d.gaussCurv ;
        
        %% Check units of Gauss curv
        if preview
            figure;
            trisurf(triangulation(mesh.f, mesh.v), K3d, 'edgecolor', 'none')
            caxis([-max(abs(K3d(:))), max(abs(K3d(:)))]) ;
            colorbar ; axis equal
            colormap(bwr)
            title('Gauss curvature')
            pause(1)
            close all

            K3dgrid = reshape(K3d, [nU, nV-1]) ;
            Kap = mean(K3dgrid, 2) ;
            plot(Kap)
            title('Gauss curvature (ap)')
            pause(1)
        end
        
        %% OPTION1: Could compute Laplacian on original mesh 
        % (ie mesh is smoothed in time, not space, loaded from disk)
        % fvel = squeeze(face_vels(tidx, :, :)) ;
        % Lv3d = DEC.laplacian(fvel) ;
        
        %% OPTION2: Load Lapv from disk and convert using smoothed mesh
        % [~, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        if strcmp(averagingStyle, 'Lagrangian')
            if doubleResolution
                dec_tp = load(sprintf(QS.fullFileBase.decAvg2x, tp)) ;
            else
                dec_tp = load(sprintf(QS.fullFileBase.decAvg, tp)) ;
            end
        else
            error('Are you sure you want to be using simple averaging?')
            if doubleResolution
                dec_tp = load(sprintf(QS.fullFileBase.decSimAvg2x, tp)) ;
            else
                dec_tp = load(sprintf(QS.fullFileBase.decSimAvg, tp)) ;
            end
        end
        
        % Smooth Laplacian(v) [lapv]
        Lv3d = dec_tp.lapvs.raw ;
        if lambda > 0
            Lv3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
                                    lambda, 'implicit', Lv3d) ;
        end
        
        % Smooth the velocities in space using gptoolbox
        vx = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 1)) ;
        vy = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 2)) ;
        vz = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, 3)) ;
        if lambda > 0
            vxs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
                lambda, 'implicit', vx') ;
            vys = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
                lambda, 'implicit', vy') ;
            vzs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
                lambda, 'implicit', vz') ;
        else
            vxs = vx' ; 
            vys = vy' ; 
            vzs = vz' ;   
        end
        
        % Compute difference / sum?
        Kv3d = K3d .* [vxs, vys, vzs] ;
        % Extend to have another row for 2d map
        Lv2d = Lv3d ;
        Lv2d((nU*(nV-1) + 1):(nU*nV), :) = Lv3d(1:nU, :) ;
        Kv2d = Kv3d ;
        Kv2d((nU*(nV-1) + 1):(nU*nV), :) = Kv3d(1:nU, :) ;
        
        % Smooth the difference 
        if lambda_err > 0 
            gP3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
                lambda_err, 'implicit', Lv3d + Kv3d) ;
            % expand error to 2d map
            gP2d = gP3d ;
            gP2d(nU*(nV-1)+1:nU*nV, :) = gP3d(1:nU, :) ;
        else
            gP2d = Lv2d + Kv2d ;
            gP3d = Lv3d + Kv3d ;
        end
    
        %% Store data on disk
        Kv = reshape(Kv2d, [nU,nV,3]) ;
        Lv = reshape(Lv2d, [nU,nV,3]) ;
        gP = reshape(gP2d, [nU,nV,3]) ;
        
        % Bandpass filter modes
        if nmodes > 0
            filterOptions.nmodes = nmodes ;
            filterOptions.zwidth = zwidth ;
            Kv_filt = Kv(:, 1:nV-1, :) ;
            gP_filt = gP(:, 1:nV-1, :) ;
            Lv_filt = gP(:, 1:nV-1, :) ;
            Kv_filt(:, :, 1) = modeFilterQuasi1D(Kv(:, 1:nV-1, 1), filterOptions) ; 
            Kv_filt(:, :, 2) = modeFilterQuasi1D(Kv(:, 1:nV-1, 2), filterOptions) ; 
            Kv_filt(:, :, 3) = modeFilterQuasi1D(Kv(:, 1:nV-1, 3), filterOptions) ;   
            gP_filt(:, :, 1) = modeFilterQuasi1D(gP(:, 1:nV-1, 1), filterOptions) ; 
            gP_filt(:, :, 2) = modeFilterQuasi1D(gP(:, 1:nV-1, 2), filterOptions) ; 
            gP_filt(:, :, 3) = modeFilterQuasi1D(gP(:, 1:nV-1, 3), filterOptions) ; 
            Lv_filt(:, :, 1) = modeFilterQuasi1D(Lv(:, 1:nV-1, 1), filterOptions) ;   
            Lv_filt(:, :, 2) = modeFilterQuasi1D(Lv(:, 1:nV-1, 2), filterOptions) ;   
            Lv_filt(:, :, 3) = modeFilterQuasi1D(Lv(:, 1:nV-1, 3), filterOptions) ;   
            
            % 3d version does not have periodic row
            Kv3d_filt = Kv_filt ;
            gP3d_filt = gP_filt ;
            Lv3d_filt = Lv_filt ;
            
            % Append periodic row
            Kv_filt(:, nV, :) = Kv_filt(:, 1, :) ;
            gP_filt(:, nV, :) = gP_filt(:, 1, :) ;
            Lv_filt(:, nV, :) = Lv_filt(:, 1, :) ;
            
            % 2d version has periodic row
            Kv2d_filt = Kv_filt ;
            gP2d_filt = gP_filt ;
            Lv2d_filt = Lv_filt ;
        end
        
        % Average along DV -- ignore last redudant row at nV
        Kv_ap = mean(Kv_filt(:, 1:nV-1, :), 2) ;
        gP_ap = mean(gP_filt(:, 1:nV-1, :), 2) ;
        Lv_ap = mean(Lv_filt(:, 1:nV-1, :), 2) ;
        
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
        Kv_l = mean(Kv_filt(:, left, :), 2) ;
        gP_l = mean(gP_filt(:, left, :), 2) ;
        Lv_l = mean(Lv_filt(:, left, :), 2) ;
        
        % right quarter
        Kv_r = mean(Kv_filt(:, right, :), 2) ;
        gP_r = mean(gP_filt(:, right, :), 2) ;
        Lv_r = mean(Lv_filt(:, right, :), 2) ;
        
        % dorsal quarter
        Kv_d = mean(Kv_filt(:, dorsal, :), 2) ;
        gP_d = mean(gP_filt(:, dorsal, :), 2) ;
        Lv_d = mean(Lv_filt(:, dorsal, :), 2) ;
        
        % ventral quarter
        Kv_v = mean(Kv_filt(:, ventral, :), 2) ;
        gP_v = mean(gP_filt(:, ventral, :), 2) ;
        Lv_v = mean(Lv_filt(:, ventral, :), 2) ;
        
        %% Save timeseries measurements
        disp('Saving timeseries measurements')
        save(Kfn, 'Kv', 'Kv_filt', 'Kv_ap', 'Kv_l', 'Kv_r', 'Kv_d', 'Kv_v')
        save(Pfn, 'gP', 'gP_filt', 'gP_ap', 'gP_l', 'gP_r', 'gP_d', 'gP_v')
        save(Lfn, 'Lv', 'Lv_filt', 'Lv_ap', 'Lv_l', 'Lv_r', 'Lv_d', 'Lv_v')
        save_lambdas = true ;
        toc
        
    else
        % Load timeseries measurements
        load(Kfn, 'Kv_filt', 'Kv_ap', 'Kv_l', 'Kv_r', 'Kv_d', 'Kv_v')
        load(Pfn, 'gP_filt', 'gP_ap', 'gP_l', 'gP_r', 'gP_d', 'gP_v')
        load(Lfn, 'Lv_filt', 'Lv_ap', 'Lv_l', 'Lv_r', 'Lv_d', 'Lv_v')
        
        Kv2d_filt = Kv_filt ;
        Kv3d_filt = Kv_filt(:, 1:nV-1, :) ;
        gP2d_filt = gP_filt ;
        gP3d_filt = gP_filt(:, 1:nV-1, :);
        Lv2d_filt = Lv_filt ;
        Lv3d_filt = Lv_filt(:, 1:nV-1, :) ;
        
        % Ensure that mesh and cutMesh are not assumed to be some other
        % timepoint's meshes
        
        if doubleResolution
            % Load current mesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC2x, tp)) ;
            mesh = tmp.spcutMeshSmRSC2x ;

            % Load cutMesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS2x, tp)) ;
            cutMesh = tmp.spcutMeshSmRS2x ;
            clearvars tmp
        else
            % Load current mesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
            mesh = tmp.spcutMeshSmRSC ;

            % Load cutMesh
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
            cutMesh = tmp.spcutMeshSmRS ;
            clearvars tmp
        end
        
        save_lambdas = false ;
    end    

    %% Save lambdas used to disk
    if save_lambdas    
        lambs = [lambda, lambda_mesh, nmodes, zwidth] ;
        header = ['laplacian smoothing and spectral filtering parameters: lambda for ', ...
            'velocity & div(v), lambda for mesh), '...
            '#modes, zwidth of spectral filter'] ;
        filename = fullfile(sFDir, 'lambdas.txt') ;
        write_txt_with_header(filename, lambs, header)
    end
    
    %% Plot results
    % operational plotting options
    pOptions.overwrite = overwrite ;
    % parameter plotting options
    pOptions.doubleResolution = doubleResolution; 
    pOptions.lambda = lambda ;
    pOptions.lambda_err = lambda_err ;
    pOptions.lambda_mesh = lambda_mesh ;
    pOptions.nmodes = nmodes ;
    pOptions.zwidth = zwidth ;
    pOptions.cutMesh = cutMesh ;
    pOptions.mesh = mesh ;
    pOptions.climit = climit ;
    
    fig3dfn = fullfile(fig3dDir, sprintf('stokesForces3d_%06d.png', tp)) ;
    fig2dfn = fullfile(fig2dDir, sprintf('stokesForces2d_%06d.png', tp)) ;
    fig3dfnZ = fullfile(fig3dDir, sprintf('stokesForces3d_zoom_%06d.png', tp)) ;
    fig2dfnZ = fullfile(fig2dDir, sprintf('stokesForces2d_zoom_%06d.png', tp)) ;

    
    % QS.plotStokesForcesTimePoint(tp, pOptions)
    
    % function here:
    if ~exist(fig3dfn, 'file') || ~exist(fig3dfnZ, 'file') || ...
       ~exist(fig2dfn, 'file') || ~exist(fig2dfnZ, 'file') || ...
       overwrite || overwriteImages
       
        % prep 2D mesh
        m2 = cutMesh ;
        m2.v = m2.u ;
        m2.v(:, 1) = m2.v(:, 1) / max(m2.v(:, 1)) ;
        cutPath = [1:nU; nU*(nV-1)+1:nU*nV]' ;

        % Push onto tangential and normal components (vertex-based)
        [Lvn, Lvt, Lvt2dS, jac, vertexnormals, g_ab, dilation] = ...
            resolveTangentNormalVertexVelocities(mesh.f, mesh.v, ...
            reshape(Lv3d_filt, [nU*(nV-1), 3]), 1:size(mesh.v, 1), ...
            m2.v, 'cutPath', cutPath) ;
        Lvt2d = Lvt2dS.uv2vertices ;
        [Kvn, Kvt, Kvt2dS] = ...
            resolveTangentNormalVertexVelocities(mesh.f, mesh.v, ...
            reshape(Kv3d_filt, [nU*(nV-1), 3]), 1:size(mesh.v, 1), ...
            m2.v, 'vertexnormals', vertexnormals, 'cutPath', cutPath) ;
        Kvt2d = Kvt2dS.uv2vertices ;
        [gPn, gPt, gPt2dS] = ...
            resolveTangentNormalVertexVelocities(mesh.f, mesh.v, ...
            reshape(gP3d_filt, [nU*(nV-1), 3]), 1:size(mesh.v, 1), ...
            m2.v, 'vertexnormals', vertexnormals, 'cutPath', cutPath) ;
        gPt2d = gPt2dS.uv2vertices ;

        % Magnitudes and angles
        Lvtmag = vecnorm(Lvt, 2, 2) ;
        Lvtang = atan2(Lvt2d(:, 2), Lvt2d(:, 1)) ;
        Kvtmag = vecnorm(Kvt, 2, 2) ;
        Kvtang = atan2(Kvt2d(:, 2), Kvt2d(:, 1)) ;
        gPtmag = vecnorm(gPt, 2, 2) ;
        gPtang = atan2(gPt2d(:, 2), gPt2d(:, 1)) ;

        % prep fields
        Lvn2d = Lvn ;
        Lvn2d(nU*(nV-1)+1:nU*nV) = Lvn(1:nU) ;
        Lvtmag2d = Lvtmag ;
        Lvtmag2d(nU*(nV-1)+1:nU*nV) = Lvtmag(1:nU) ;
        Lvtang2d = Lvtang ;
        Lvtang2d(nU*(nV-1)+1:nU*nV) = Lvtang(1:nU) ;

        Kvn2d = Kvn ;
        Kvn2d(nU*(nV-1)+1:nU*nV) = Kvn(1:nU) ;
        Kvtmag2d = Kvtmag ;
        Kvtmag2d(nU*(nV-1)+1:nU*nV) = Kvtmag(1:nU) ;
        Kvtang2d = Kvtang ;
        Kvtang2d(nU*(nV-1)+1:nU*nV) = Kvtang(1:nU) ;

        % gP as 2d unraveled
        gPn2d = gPn ;
        gPn2d(nU*(nV-1)+1:nU*nV) = gPn(1:nU) ;
        gPtmag2d = gPtmag ;
        gPtmag2d(nU*(nV-1)+1:nU*nV) = gPtmag(1:nU) ;
        gPtang2d = gPtang ;
        gPtang2d(nU*(nV-1)+1:nU*nV) = gPtang(1:nU) ;

        % 3D plot
        if ~exist(fig3dfn, 'file') || overwriteImages
            opts = struct() ;
            opts.labels = {'$\big($d$\star$d$v_{\parallel}^\flat\big)_\perp$', ...
                '$Kv_\perp$', ...
                '$(\nabla p)_\perp$', ...
                '$\big($d$\star$d$v_{\parallel}^\flat\big)_\perp$', ...
                 '$Kv_\parallel$', ...
                '$(\nabla p)_\parallel$'} ;
            opts.xyzlim = xyzlim ;
            opts.view = [0, 0] ;
            opts.clim = climit ;
            opts.polarStyle = 'polar'; 
            opts.axisOff = true ;
            opts.makeCbar = false(6, 1) ;
            opts.masterCbar = true ;
            [axs, cbs, meshHandles] = ...
                nFieldsOnSurface({mesh,mesh,mesh,mesh,mesh,mesh}, ...
                {Lvn, Kvn, gPn, ...
                {Lvtmag, Lvtang}, {Kvtmag, Kvtang}, {gPtmag, gPtang}}, ...
                opts) ;
            % Extra colorbars
            extraAxPos = [0.39, 0.55, 0.22, 0.0341] ;
            phasePos = [0.905, 0.05, 0.07, 0.135] ;
            phasebar('colormap', phasemap, ...
                 'location', phasePos, 'style', 'polar') ;
            axpos = get(axs{1}, 'position') ;
            set(gcf, 'CurrentAxes', axs{1})
            tmp = colorbar('location', 'southOutside') ;
            set(axs{1}, 'position', axpos) ;
            set(tmp, 'position', extraAxPos) ;

            disp(['Saving 3d mesh stokes figure: ' fig3dfn])
            export_fig(fig3dfn, '-transparent', '-nocrop', '-r200')
        end
        
        %% 2d plots
        opts.view = [0, 90] ;
        opts.makeCbar = false(6, 1) ;
        nFieldsOnSurface({m2,m2,m2, m2,m2,m2}, ...
            {Lvn2d, Kvn2d, gPn2d, ...
            {Lvtmag2d, Lvtang2d}, {Kvtmag2d, Kvtang2d}, {gPtmag2d, gPtang2d}}, ...
            opts)
        phasebar('colormap', phasemap, ...
             'location', phasePos, 'style', 'polar') ;
        axpos = get(axs{1}, 'position') ;
        set(gcf, 'CurrentAxes', axs{1})
        tmp = colorbar('location', 'southOutside') ;
        set(axs{1}, 'position', axpos) ;
        set(tmp, 'position', extraAxPos) ;

        disp(['Saving 2d mesh stokes figure: ' fig2dfn])
        export_fig(fig2dfn, '-transparent', '-nocrop', '-r200')

        % 3D plot ZOOM
        opts.clim = climit * 0.3 ;
        opts.makeCbar = false(6, 1) ;
        nFieldsOnSurface({mesh,mesh,mesh,mesh,mesh,mesh}, ...
            {Lvn, Kvn, gPn, ...
            {Lvtmag, Lvtang}, {Kvtmag, Kvtang}, {gPtmag, gPtang}}, ...
            opts)
        phasebar('colormap', phasemap, ...
             'location', phasePos, 'style', 'polar') ;
        axpos = get(axs{1}, 'position') ;
        set(gcf, 'CurrentAxes', axs{1})
        tmp = colorbar('location', 'southOutside') ;
        set(axs{1}, 'position', axpos) ;
        set(tmp, 'position', extraAxPos) ;

        disp(['Saving 3d zoom stokes fig: ' fig3dfnZ])
        export_fig(fig3dfnZ, '-transparent', '-nocrop', '-r200')

        % 2D plot ZOOM
        opts.view = [0, 90] ;
        nFieldsOnSurface({m2,m2,m2, m2,m2,m2}, ...
            {Lvn2d, Kvn2d, gPn2d, ...
            {Lvtmag2d, Lvtang2d}, {Kvtmag2d, Kvtang2d}, {gPtmag2d, gPtang2d}}, ...
            opts)
        phasebar('colormap', phasemap, ...
             'location', phasePos, 'style', 'polar') ;
        axpos = get(axs{1}, 'position') ;
        set(gcf, 'CurrentAxes', axs{1})
        tmp = colorbar('location', 'southOutside') ;
        set(axs{1}, 'position', axpos) ;
        set(tmp, 'position', extraAxPos) ;

        disp(['Saving 2d zoom stokes fig: ' fig2dfnZ])
        export_fig(fig2dfnZ, '-transparent', '-nocrop', '-r200')
    end
    
end
