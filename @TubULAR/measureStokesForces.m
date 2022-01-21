function measureStokesForces(QS, options)
%[gP_apM, Kv_apM, Lv_apM, veln_apM] = measureStokesForces(QS, options)
%   Measure \nabla_j p / eta = (Laplace + K) v_i, meaning the laplacian of
%   the flow plus gaussian curvature times the flow give the force in the
%   tissue.
%   
%   Note: this version has tangential basis decomposition for Laplacian
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
coordSys = 'spsme' ;
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
    nU = QS.nU ;
    nV = QS.nV ;
    piv3dFileBase = QS.fullFileBase.piv3d ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
    sresStr = 'doubleRes_' ;    
    nU = 2 * QS.nU - 1 ;
    nV = 2 * QS.nV - 1 ;
    piv3dFileBase = QS.fullFileBase.piv3d2x ;

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
lastIdx = ntps - 1 ;
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
    Kfn = fullfile(outdir, sprintf('K2v_vertices_%06d.mat', tp)) ;
    Pfn = fullfile(outdir, sprintf('gradP_vertices_%06d.mat', tp)) ;
    Lfn = fullfile(outdir, sprintf('Lapv_vertices_%06d.mat', tp)) ;
    grad2Trfn = fullfile(outdir, sprintf('grad2Tr_vertices_%06d.mat', tp)) ;
    gvnbM2Hfn = fullfile(outdir, sprintf('gvnbM2H_vertices_%06d.mat', tp)) ;
    
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
        H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
        H2d = extendFieldAlongCutMeshSeam(H3d, cutMesh) ;
        
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
        
        % Actual normal velocity -- currently using ORIGINAL mesh normals
        veln3d = sum(cat(2, vxs, vys, vzs) .* mesh.vn, 2) ;
        veln2d = veln3d ;
        veln2d(nU*(nV-1)+1:nU*nV) = veln2d(1:nU) ;

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
        % if strcmp(averagingStyle, 'Lagrangian')
        %     if doubleResolution
        %         dec_tp = load(sprintf(QS.fullFileBase.decAvg2x, tp)) ;
        %     else
        %         dec_tp = load(sprintf(QS.fullFileBase.decAvg, tp)) ;
        %     end
        % else
        %     error('Are you sure you want to be using simple averaging?')
        %     if doubleResolution
        %         dec_tp = load(sprintf(QS.fullFileBase.decSimAvg2x, tp)) ;
        %     else
        %         dec_tp = load(sprintf(QS.fullFileBase.decSimAvg, tp)) ;
        %     end
        % end
        % 
        % % Smooth Laplacian(v) [lapv]
        % Lv3d = dec_tp.lapvs.raw ;

        
        %% NOTE: here we do tangential basis decomposition of Laplacian
        QS.getVelocityAverage('vf', 'v2dum') ;
        vf = QS.velocityAverage.vf ;
        % v2dum = QS.velocityAverage.v2dum ;
        
        % Obtain smoothed velocities on all faces
        vfsm = squeeze(vf(tidx, :, :)) ;
        % v2dsmum_ii = squeeze(v2dum(tidx, :, :)) ;

        % Use current time's tiled smoothed mesh
        % Note: vfsmM is in um/min rs
        piv3d = load(sprintf(piv3dFileBase, tp)) ;
        piv3d = piv3d.piv3dstruct ;
        FF = piv3d.m0f ;   % #facesx3 float: mesh connectivity list
        % V2D = piv3d.m0XY ; % Px2 float: 2d mesh vertices in pullback image pixel space
        % v3drs = QS.xyz2APDV(piv3d.m0v3d) ;
        
        % Grab cutMesh from piv3d. Could grab from disk instead...
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
        
        cutM.f = FF ;
        cutM.u = cutMesh.u ;
        V2D = cutM.u ;
        cutM.v = cutMesh.v ; % v3drs ;
        cutM.nU = nU ;
        cutM.nV = nV ;
        
        % Now resolve the vector field for decomposition
        [v0n, v0t, v0t2d, jac3d_to_2d, ~, ~, dilation] = ...
            resolveTangentNormalVelocities(cutM.f, cutM.v, vfsm, ...
            1:length(FF), cutM.u ) ;

        m2d = cutM ;
        m2d.v = cat(2, cutM.u, 0*cutM.u(:, 1)) ;
        % DEC2d = DiscreteExteriorCalculus(m2d.f, m2d.v) ;
        % DEC2d.laplacianVectorField(v0t, cutM.u, cutM.v)
        
        % -------------------------------------
        % DEC 1-form method for Laplacian term
        % -------------------------------------
        % This is -(\delta \diff v^{\flat})^\sharp 
        % Same as -g^{ab} (\delta d v^{\flat})_b as 3d vector (ie sharped)
        %  > flat vector on faces to one form (flatDP) 
        %  > 1form to 2form (d1) > dual zero form (hd2) 
        %  > dual one form (dd0) 
        %  > star (inv(hd1)) brings to primal one form.
        %  > Lastly, sharp DP maps primal dual form to tangent vectors on
        %       faces -- note this is -delta in the notation of Arroyo &
        %       DeSimone 2009 PRE equation 18.
        nFaces = size(mesh.f, 1) ;
        [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        Lvf = DEC.sharpPD * inv(DEC.hd1) * DEC.dd0 * DEC.hd2 * DEC.d1 * DEC.dualVectorToPrimal1Form(v0t) ;
        Lvf = reshape(Lvf, [nFaces, 3]) ;
        Lv3d = F2V * Lvf ;
        % Smooth results
        if lambda > 0
            Lv3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
                                    lambda, 'implicit', Lv3d) ;
        end
        Lv2d = extendFieldAlongCutMeshSeam(Lv3d, cutMesh) ;
        
        [~, ~, Lv3d_uv, ~, ~, ~, ~] = ...
            resolveTangentNormalVelocities(cutM.f, cutM.v, Lvf, ...
            1:length(FF), cutM.u ) ;
        Lv3d_uv = F2V * Lv3d_uv ;
        Lv_magf = vecnorm(Lvf, 2, 2) ;
        Lv_angf = atan2(Lv3d_uv(:, 2), Lv3d_uv(:, 1)) ;
        Lv_mag3d = vecnorm(Lv3d, 2, 2) ;
        Lv_ang3d = atan2(Lv3d_uv(:, 2), Lv3d_uv(:, 1)) ;
        Lv_mag2d = extendFieldAlongCutMeshSeam(Lv_mag3d, cutMesh) ;
        Lv_ang2d = extendFieldAlongCutMeshSeam(Lv3d_ang, cutMesh) ;
        
        % Check it
        opts = struct() ;
        opts.mesh = mesh ;
        opts.climit = 0.02 ;
        plotPolarField(Lv_mag3d, Lv_ang3d, opts)
        opts = struct() ;
        opts.mesh = m2d ;
        opts.climit = 1 ;
        plotPolarField(Lv_mag2d, Lv_ang2d, opts)
        
        % ------------------------------
        % Component-by-component method
        % ------------------------------
        % vecs = v0t ;
        % gr = grad(V2D, cutM.f) ;
        % euev = gr * cutM.v ;
        % eu = euev(1:length(cutM.f), :) ;
        % ev = euev(length(cutM.f)+1:2*length(cutM.f), :) ;
        % 
        % % without normalization 
        % % --> note that we have uvecs and vvecs from v0t2d already.
        % % wu = dot(eu, vecs, 2) ;
        % % wv = dot(ev, vecs, 2) ;
        % % uvecs_num = wu - (wv .* dot(eu,ev,2)) ./ vecnorm(ev, 2, 2).^2 ;
        % % uvecs_den = vecnorm(eu, 2, 2).^2 - dot(eu, ev, 2).^2 ./ vecnorm(ev, 2, 2).^2 ;
        % % uvecs = uvecs_num ./ uvecs_den ;
        % % vvecs = (wu - uvecs .* vecnorm(eu, 2, 2).^2) ./ dot(eu, ev, 2) ;
        % % assert(max(abs(v0t2d(:, 1) - uvecs(:))) < 1e-15)
        % 
        % Lv1 = DEC.laplacian(F2V * v0t2d(:, 1)) ;
        % Lv2 = DEC.laplacian(F2V * v0t2d(:, 2)) ;
        % 
        % Lv_uv = [Lv1, Lv2] ;
        % 
        % % Push to 3d
        % Lv3d = Lv1 .* (F2V * eu) + Lv2 .* (F2V * ev) ;
        % 
        % % Plot in 2d and 3d
        % Lvmag = vecnorm(Lv3d, 2, 2) ; % sqrt(Lv1.^2 + Lv2.^2) ./  (F2V*dilation') ;
        % Lvang = atan2(Lv2, Lv1) ;
        % Lvmag2D = Lvmag ;
        % Lvang2D = Lvang ;
        % Lvmag2D(nU*(nV-1)+1:nU*nV) = Lvmag2D(1:nU) ;
        % Lvang2D(nU*(nV-1)+1:nU*nV) = Lvang2D(1:nU) ;
        % opts = struct() ;
        % opts.axisOff = true ;
        % opts.labels = {'$\nabla^2 v_\parallel$', '$\nabla^2 v_\parallel$'} ;
        % opts.polarStyle = 'polar' ;
        % opts.view = {[0,0], [0, 90]} ;
        % opts.climit = 2*rms1d(Lvmag) ;
        % nFieldsOnSurface({mesh, m2d}, {{Lvmag, Lvang}, {Lvmag2D, Lvang2D}}, opts)
        % set(gcf, 'visible', 'on')
        % saveas(gcf, fullfile(fig2dDir, 'factors_Lv.png'))
        % 
        % %% Checks: uvecs = v0t2d, and vvecs = v0t2d 
        % % figure ;
        % % subplot(2, 3, 1)
        % % trisurf(triangulation(m2d.f, m2d.v), uvecs, 'edgecolor', 'none')
        % % view(2); axis equal; colorbar
        % % title('uvecs')
        % % subplot(2, 3, 2)
        % % trisurf(triangulation(m2d.f, m2d.v), v0t2d(:, 1), 'edgecolor', 'none')
        % % view(2); axis equal; colorbar
        % % title('v0t2d')
        % % subplot(2, 3, 3)
        % % trisurf(triangulation(m2d.f, m2d.v), uvecs ./ v0t2d(:, 1), 'edgecolor', 'none')
        % % title('difference')
        % % view(2); axis equal; colorbar
        % % 
        % % % Check 2
        % % subplot(2, 3, 4)
        % % trisurf(triangulation(m2d.f, m2d.v), vvecs, 'edgecolor', 'none')
        % % view(2); axis equal; axis off; colorbar
        % % title('uvecs')
        % % subplot(2, 3, 5)
        % % trisurf(triangulation(m2d.f, m2d.v), v0t2d(:, 2), 'edgecolor', 'none')
        % % view(2); axis equal ; axis off; colorbar
        % % title('v0t2d')
        % % subplot(2, 3, 6)
        % % trisurf(triangulation(m2d.f, m2d.v), vvecs ./ v0t2d(:, 2), 'edgecolor', 'none')
        % % title('difference')
        % % view(2); axis equal; axis off; colorbar
        
        % -----------------------------------------------
        % 3D "method" for laplacian term -- likely wrong
        % ------------------------------------------------
        % This is too simple-minded: 
        % Lv3d = DEC.laplacianVectorField(v0t, cutM.u) ;
       
        % done -------------------------------------------
        
        %% TERM 4
        % Gaussian curvature term
        % Tangential component of velocities only, but on vertices
        v0t_vtx = F2V * v0t ;
        Kv3d = 2 * K3d .* v0t_vtx ;
        
        %% TERM 2
        % grad trace term
        efn = sprintf(QS.fullFileBase.strainRate, tp) ;
        load(efn, 'tre_vtx')
        tre_vtx3d = tre_vtx(:) ;
        tre_vtx3d = tre_vtx3d(1:nU*(nV-1)) ;
        
        % QS.fullFileBase.metricKinematics.gdot, tp);
        gradtr2_faces = 2 * DEC.gradient(tre_vtx3d) ;
        gradtr2_3d = F2V * gradtr2_faces ;
        
        
        %% TERM 3 -- normal velocity (vn) term %%
        % on faces
        % ---------
        % vn = load(QS.fileName.pivAvg.vn) ;
        % Note we could use veln3d from earlier, but then we'd have to push
        % to faces.
        % vn_faces = v0n ;
        % Hfaces = V2F * H3d ;
        % stuff = vn_faces .* (bb - 2 * Hfaces .* gg) ;
        nFaces = length(cutM.f) ;
        [gg, bb] = constructFundamentalForms(cutM.f, cutM.v, cutM.u) ;
        ggInv = zeros(length(gg), 2, 2) ;
        bbInv = zeros(length(bb), 2, 2) ;
        for faceId = 1:length(gg)
            ggInv(faceId, :, :) = inv(gg{faceId}) ;
            bbInv(faceId, :, :) = inv(bb{faceId}) ;
        end
        % assert(all(ggInv(:, 1, 2) == ggInv(:, 2, 1)))
        H3d_faces = V2F * H3d ;
        bM2H_faces = (bbInv - 2 * H3d_faces .* ggInv) ;
        
        % 2-dimensional basis gradient to get basis vectors        
        gr = grad(V2D, cutM.f) ;
        % euev = gr * cutM.v ;
        % eu = euev(1:length(cutM.f), :) ;
        % ev = euev(length(cutM.f)+1:2*length(cutM.f), :) ;
        
        gradvn = (gr * veln2d) ;
        gradvn = [gradvn(1:nFaces), gradvn(nFaces+1:end)] ;
        gvnbM2H_faces = zeros(nFaces, 2) ;
        for faceId = 1:nFaces
            % matrix product of grad(vn) with (b -2Hg) on this face
            gvnbM2H_faces(faceId, :) = ...
                -2 * gradvn(faceId, :) * squeeze(bM2H_faces(faceId, :, :)) ;
        end      
        
        % euhatF = eu ./ vecnorm(eu, 2, 2) ;
        % evhatF = ev ./ vecnorm(ev, 2, 2) ;
        % euhatV = (F2V * eu) ./ vecnorm((F2V * eu), 2, 2) ;
        % evhatV = (F2V * ev) ./ vecnorm((F2V * ev), 2, 2) ;
        
        % gvnbM2H3d_faces = gvnbM2H_faces(:, 1) .* euhatF + gvnbM2H_faces(:, 2) .* evhatF;
        % gvnbM2H_Fvtx_uv = F2V * gvnbM2H_faces ;
        % gvnbM2H_Fvtx = F2V * gvnbM2H3d_faces ;
        
        % gvnbM2H_faces is CONTRAVARIANT!
        gvnbM2H_3d = pushVectorField2Dto3DMesh(gvnbM2H_faces, V2D, v3d, faces, 1:size(faces,1)) ;
        
        if lambda > 0
            gvnbM2H_3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
                                    lambda, 'implicit', gvnbM2H_3d) ;
        end
        gvnbM2H_2d = extendFieldAlongCutMeshSeam(gvnbM2H_mag3d, cutMesh) ;
        
        % polar coordinate version
        gvnbM2H_mag3d = vecnorm(gvnbM2H_3d, 2, 2) ;
        gvnbM2H_mag2d = extendFieldAlongCutMeshSeam(gvnbM2H_mag3d, cutMesh) ;
        
        [~, ~, gvnbM2H_uv, ~, ~, ~, ~] = ...
            resolveTangentNormalVelocities(cutM.f, cutM.v, gvnbM2H_2d, ...
            1:length(FF), cutM.u ) ;
        gvnbM2H_ang2d = atan2(gvnbM2H_uv(:, 2), gvnbM2H_uv(:, 1)) ;
        gvnbM2H_ang3d = gvnbM2H_ang3d(1:nU*(nV-1), :) ;
        
        % % on vertices
        % % -----------
        % % Create INVERSE fundamental forms
        % % push bb and gg onto vertices
        % bbIvtx = zeros(nU*(nV-1), 2, 2) ;
        % bbIvtx(:, 1, 1) = F2V * bbInv(:, 1, 1) ;
        % bbIvtx(:, 1, 2) = F2V * bbInv(:, 1, 2) ;
        % bbIvtx(:, 2, 1) = F2V * bbInv(:, 2, 1) ;
        % bbIvtx(:, 2, 2) = F2V * bbInv(:, 2, 2) ;
        % ggIvtx = zeros(nU*(nV-1), 2, 2) ;
        % ggIvtx(:, 1, 1) = F2V * ggInv(:, 1, 1) ;
        % ggIvtx(:, 1, 2) = F2V * ggInv(:, 1, 2) ;
        % ggIvtx(:, 2, 1) = F2V * ggInv(:, 2, 1) ;
        % ggIvtx(:, 2, 2) = F2V * ggInv(:, 2, 2) ;
        % bM2H = (bbIvtx - 2 * H3d .* ggIvtx) ;
        % gradvnV = F2V * gradvn ;
        % if lambda > 0
        %     gradvnV = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
        %                             lambda, 'implicit', gradvnV) ;
        %     bM2H = laplacian_smooth(mesh.v, mesh.f, 'cotan', [],...
        %                             lambda, 'implicit', bM2H) ;      
        % end
        % gvnbM2H_vtx = zeros(nU*(nV-1), 2) ;
        % for vId = 1:size(ggIvtx, 1)
        %     gvnbM2H_vtx(vId, :) = ...
        %         -2 * gradvnV(vId, :) * squeeze(bM2H(vId, :, :)) ;
        % end
        % % for comparison, look at vtx calculation compared to face calc 
        % % pushed to vtx
        % gvnbM2H = gvnbM2H_vtx(:, 1) .* euhatV + gvnbM2H_vtx(:, 2) .* evhatV;
        % 

        % polar coordinate version
        gvnbM2H_vtx_mag3d = vecnorm(gvnbM2H, 2, 2) ;
        gvnbM2H_vtx_mag2d = gvnbM2H_vtx_mag3d ;
        gvnbM2H_vtx_mag2d(nU*(nV-1)+1:nU*nV) = gvnbM2H_vtx_mag2d(1:nU) ;
        gvnbM2H_vtx_ang3d = atan2(gvnbM2H_vtx(:, 2), gvnbM2H_vtx(:, 1)) ;
        gvnbM2H_vtx_ang2d = gvnbM2H_vtx_ang3d ;
        gvnbM2H_vtx_ang2d(nU*(nV-1)+1:nU*nV) = gvnbM2H_vtx_ang2d(1:nU) ;
        
        % clf
        % quiver3(mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        %     gvnbM2H(:, 1), gvnbM2H(:, 2), gvnbM2H(:, 3), 1)
        opts = struct() ;
        opts.polarStyle = 'polar' ;
        opts.view = {[0,0], [0, 90], [0,0], [0, 90]} ;
        opts.axisOff = true ;
        opts.visible = 'on' ;
        Fclim = rms1d(gvnbM2H_Fvtx_mag3d(:)) ;
        % Vclim = 0.1 * rms1d(gvnbM2H_vtx_mag3d(:)) ;
        opts.labels = {'$-2 \nabla v_n \cdot (b - 2H g)^\sharp$', ...
            '$-2\nabla v_n \cdot (b - 2H g)^\sharp$', ...
            '$-2\nabla v_n \cdot (b - 2H g)^\sharp$', ...
            '$-2\nabla v_n \cdot (b - 2H g)^\sharp$'} ;
        opts.clims = {Fclim, Fclim, Fclim, Fclim} ;
        nFieldsOnSurface({mesh, m2d, mesh, m2d}, ...
            {{gvnbM2H_Fvtx_mag3d, gvnbM2H_Fvtx_ang3d}, ...
            {gvnbM2H_Fvtx_mag2d, gvnbM2H_Fvtx_ang2d}, ...
            {gvnbM2H_vtx_mag3d, gvnbM2H_vtx_ang3d}, ...
            {gvnbM2H_vtx_mag2d, gvnbM2H_vtx_ang2d}}, opts) ;
        set(gcf, 'visible', 'on')
        saveas(gcf, fullfile(fig2dDir, 'term3_factors_term3full.png'))
        
        % Plot each factor in the vn term
        [gradvnV_mag3d, gradvnV_ang3d] = cartesian2polar(gradvnV) ;
        gradvnV_mag2d = extendFieldAlongCutMeshSeam(gradvnV_mag3d, cutMesh) ;
        gradvnV_ang2d = extendFieldAlongCutMeshSeam(gradvnV_ang3d, cutMesh) ;
        ggIvtx2d = extendFieldAlongCutMeshSeam(ggIvtx, cutMesh) ;
        bbIvtx2d = extendFieldAlongCutMeshSeam(bbIvtx, cutMesh) ;
        Hgvtx2d = extendFieldAlongCutMeshSeam(-2 * H3d .* ggIvtx, cutMesh) ;
        opts = struct() ;
        opts.polarStyle = 'polar' ;
        
        % GRAD vn
        opts.view = {[0,0], [0, 90]} ;
        opts.axisOff = true ;
        opts.visible = 'on' ;
        opts.labels = {'$\nabla v_n$', '$\nabla v_n$'} ;
        opts.subplotGrouping = [1,2] ;
        opts.climit = 2*rms1d(gradvnV_mag3d(:)) ;
        nFieldsOnSurface({mesh, m2d, ...
            }, ...
            {{gradvnV_mag3d, gradvnV_ang3d}, ...  % grad v_n
            {gradvnV_mag2d, gradvnV_ang2d}, ...
            }, opts) ;
        saveas(gcf, fullfile(fig2dDir, 'term3_factors_gradvn.png'))
        
        
        % b^{ab} plot
        opts.view = {[0,0], [0, 90], [0,0], [0, 90], [0,0], [0, 90]} ;
        opts.axisOff = true ;
        opts.visible = 'on' ;
        opts.clim = 2 * rms1d(bbIvtx2d(:)) ;
        opts.labels = {'$b^{11}$', '$b^{11}$', '$b^{12}$', '$b^{12}$', ...
            '$b^{22}$', '$b^{22}$'} ;
        opts.subplotGrouping = [2,4] ;
        % opts.climit = max(gvnbM2H_vtx_mag3d(:)) ;
        nFieldsOnSurface({mesh, m2d,mesh, m2d,mesh, m2d, ...
            }, ...
            {bbIvtx(:, 1, 1), bbIvtx2d(:, 1, 1), ... % bb
            bbIvtx(:, 1, 2), bbIvtx2d(:, 1, 2), ... % bb2
            bbIvtx(:, 2, 2), bbIvtx2d(:, 2, 2), ... % bb4
            }, opts) ;
        saveas(gcf, fullfile(fig2dDir, 'term3_factors_bb.png'))
        
        
        % g^{ab} plot
        clf
        opts = struct() ;
        opts.view = {[0,0], [0, 90], [0,0], [0, 90], [0,0], [0, 90]} ;
        opts.axisOff = true ;
        opts.visible = 'on' ;
        opts.clim = 2 * rms1d(ggIvtx(:)) ;
        opts.labels = {'$g^{11}$', '$g^{11}$', '$g^{12}$',  '$g^{12}$', ...
                 '$g^{22}$', '$g^{22}$'} ;
        opts.subplotGrouping = [2,4] ;
        % opts.climit = max(gvnbM2H_vtx_mag3d(:)) ;
        nFieldsOnSurface({mesh, m2d,mesh, m2d,mesh, m2d, ...
            }, ...
            {ggIvtx(:, 1, 1), ggIvtx2d(:, 1, 1), ... % gg
            ggIvtx(:, 1, 2), ggIvtx2d(:, 1, 2), ... % gg2
            ggIvtx(:, 2, 2), ggIvtx2d(:, 2, 2), ... % gg4
            }, opts)
        saveas(gcf, fullfile(fig2dDir, 'term3_factors_gg.png'))
        
        % g^{ab} plot
        opts = struct() ;
        opts.view = {[0,0], [0, 90]} ;
        opts.axisOff = true ;
        opts.visible = 'on' ;
        opts.clim = 2 * rms1d(H3d(:)) ;
        opts.labels = {'$H$', '$H$'} ;
        opts.subplotGrouping = [1, 2] ;
        % opts.climit = max(gvnbM2H_vtx_mag3d(:)) ;
        nFieldsOnSurface({mesh, m2d}, {H3d, H2d}, opts)
        saveas(gcf, fullfile(fig2dDir, 'term3_factors_HH.png'))
        
        % 2*H*g^{ab} plot
        opts.view = {[0,0], [0, 90], [0,0], [0, 90], [0,0], [0, 90]} ;
        opts.axisOff = true ;
        opts.visible = 'on' ;
        H2g3d = -2 * H3d .* ggIvtx ;
        opts.clim = 2 * rms1d(H2g3d(:)) ;
        opts.labels = {'$-2Hg^{11}$', '$-2Hg^{11}$', '$-2Hg^{12}$',  '$-2Hg^{12}$', ...
                 '$-2Hg^{22}$', '$-2Hg^{22}$'} ;
        opts.subplotGrouping = [2,4] ;
        % opts.climit = max(gvnbM2H_vtx_mag3d(:)) ;
        nFieldsOnSurface({mesh, m2d,mesh, m2d,mesh, m2d, ...
            }, ...
            {H2g3d(:, 1, 1), Hgvtx2d(:, 1, 1), ... % gg
            H2g3d(:, 1, 2), Hgvtx2d(:, 1, 2), ... % gg2
            H2g3d(:, 2, 2), Hgvtx2d(:, 2, 2), ... % gg4
            }, opts)
        saveas(gcf, fullfile(fig2dDir, 'term3_factors_H2g.png'))        
        
        % Extend to have another row for 2d map
        Lv2d = extendFieldAlongCutMeshSeam(Lv3d, cutMesh) ;
        Kv2d = extendFieldAlongCutMeshSeam(Kv3d, cutMesh) ;
        gvnbM2H_3d = gvnbM2H ;
        gradtr2_2d = extendFieldAlongCutMeshSeam(gradtr2_3d, cutMesh) ;
        gvnbM2H_2d = extendFieldAlongCutMeshSeam(gvnbM2H_3d, cutMesh) ;
        
        % Smooth the forces (difference/sum) 
        gP3d = Lv3d + Kv3d + gradtr2_3d + gvnbM2H_3d ;
        
        if lambda_err > 0 
            gP3d = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
                lambda_err, 'implicit', gP3d) ;
            % expand error to 2d map
            gP2d = extendFieldAlongCutMeshSeam(gP3d, cutMesh) ;
        else
            gP2d = Lv2d + Kv2d + gradtr2_2d + gvnbM2H_2d ;
        end
        
        %% Take u,v components of each term and the result
        % --> already done for gvnbM2H_vtx
        gvnbM2H_uv3d = gvnbM2H_vtx ;
        gvnbM2H_uv2d = extendFieldAlongCutMeshSeam(gvnbM2H_uv3d, cutMesh) ;
        Kv_uv3d = [dot(Kv3d, (F2V*eu), 2), dot(Kv3d, (F2V*ev), 2)] ;
        Kv_uv2d = extendFieldAlongCutMeshSeam(Kv_uv3d, cutMesh) ;
        gradtr2_uv3d = [dot(Kv3d, (F2V*eu), 2), dot(Kv3d, (F2V*ev), 2)] ;
        gradtr2_uv2d = extendFieldAlongCutMeshSeam(gradtr2_uv3d, cutMesh) ;
        
        % Check that we get the same thing by subtracting in uv as
        % projecting gP3d along uv 
        gP_uv3d = [dot(gP3d, (F2V*eu), 2), dot(gP3d, (F2V*ev), 2)] ;
        gP_uv2d = extendFieldAlongCutMeshSeam(gP_uv3d, cutMesh) ;
        
        %% Plot final forces
        clf ;
        opts = struct() ;
        opts.view = {[0,0], [0, 90], [0,0], [0, 90], [0,0], [0, 90]} ;
        opts.axisOff = true ;
        opts.visible = 'on' ;
        opts.polarStyle = 'polar' ;
        opts.clim = 2*rms1d(Lvmag) ;
        opts.labels = {'$2\nabla_b \varepsilon^{ab}$', ...
            '$2\nabla_b \varepsilon^{ab}$'} ;
        opts.subplotGrouping = [1, 2] ;
        % opts.climit = max(gvnbM2H_vtx_mag3d(:)) ;
        nFieldsOnSurface({mesh, m2d}, ...
            {{vecnorm(gP3d, 2, 2), atan2(gP_uv3d(:, 2), gP_uv3d(:, 1))}, ...
            {vecnorm(gP2d, 2, 2), atan2(gP_uv2d(:, 2), gP_uv2d(:, 1))}}, opts)
        saveas(gcf, fullfile(fig2dDir, 'gPuv.png'))
        
        
        %% Store data on disk
        % Lv_uv2d = reshape(Lv_uv2d, [nU, nV, 2]) ;     
        % Kv_uv2d = reshape(Kv_uv2d, [nU, nV, 2]) ;     
        % gradtr2_uv2d = reshape(gradtr2_uv2d, [nU, nV, 2]) ;
        % gvnbM2H_uv2d = reshape(gvnbM2H_uv2d, [nU, nV, 2]) ;
        % gP_uv2d = reshape(gP_uv2d, [nU, nV, 2]) ;
        % 3d vectors for filtering
        Lv = reshape(Lv2d, [nU,nV,3]) ;
        grad2Tr = reshape(gradtr2_2d, [nU, nV, 3]) ;
        gvnbM2H = reshape(gvnbM2H_2d, [nU, nV, 3]) ;
        Kv = reshape(Kv2d, [nU,nV,3]) ;
        gP = reshape(gP2d, [nU,nV,3]) ;
        
        % Bandpass filter modes
        if nmodes > 0
            filterOptions.nmodes = nmodes ;
            filterOptions.zwidth = zwidth ;
            Lv_filt = Lv(:, 1:nV-1, :) ;
            grad2Tr_filt = grad2Tr(:, 1:nV-1, :) ;
            gvnbM2H_filt = gvnbM2H(:, 1:nV-1, :) ;
            Kv_filt = Kv(:, 1:nV-1, :) ;
            gP_filt = gP(:, 1:nV-1, :) ;
            Lv_filt(:, :, 1) = modeFilterQuasi1D(Lv_filt(:,:,1), filterOptions) ;   
            Lv_filt(:, :, 2) = modeFilterQuasi1D(Lv_filt(:,:,2), filterOptions) ;   
            Lv_filt(:, :, 3) = modeFilterQuasi1D(Lv_filt(:,:,3), filterOptions) ;   
            grad2Tr_filt(:, :, 1) = modeFilterQuasi1D(grad2Tr_filt(:,:,1), filterOptions) ;   
            grad2Tr_filt(:, :, 2) = modeFilterQuasi1D(grad2Tr_filt(:,:,2), filterOptions) ;   
            grad2Tr_filt(:, :, 3) = modeFilterQuasi1D(grad2Tr_filt(:,:,3), filterOptions) ;   
            gvnbM2H_filt(:, :, 1) = modeFilterQuasi1D(gvnbM2H_filt(:,:,1), filterOptions) ;   
            gvnbM2H_filt(:, :, 2) = modeFilterQuasi1D(gvnbM2H_filt(:,:,2), filterOptions) ;   
            gvnbM2H_filt(:, :, 3) = modeFilterQuasi1D(gvnbM2H_filt(:,:,3), filterOptions) ;   
            Kv_filt(:, :, 1) = modeFilterQuasi1D(Kv(:, 1:nV-1, 1), filterOptions) ; 
            Kv_filt(:, :, 2) = modeFilterQuasi1D(Kv(:, 1:nV-1, 2), filterOptions) ; 
            Kv_filt(:, :, 3) = modeFilterQuasi1D(Kv(:, 1:nV-1, 3), filterOptions) ;   
            gP_filt(:, :, 1) = modeFilterQuasi1D(gP(:, 1:nV-1, 1), filterOptions) ; 
            gP_filt(:, :, 2) = modeFilterQuasi1D(gP(:, 1:nV-1, 2), filterOptions) ; 
            gP_filt(:, :, 3) = modeFilterQuasi1D(gP(:, 1:nV-1, 3), filterOptions) ; 
            
            % 3d version does not have periodic row
            Lv3d_filt = Lv_filt ;
            grad2Tr3d_filt = grad2Tr_filt ;
            gvnbM2H3d_filt = gvnbM2H_filt ;
            Kv3d_filt = Kv_filt ;
            gP3d_filt = gP_filt ;
            
            % Append periodic row
            Lv_filt(:, nV, :) = Lv_filt(:, 1, :) ;
            grad2Tr_filt(:, nV, :) = grad2Tr_filt(:, 1, :) ;
            gvnbM2H_filt(:, nV, :) = gvnbM2H_filt(:, 1, :) ;
            Kv_filt(:, nV, :) = Kv_filt(:, 1, :) ;
            gP_filt(:, nV, :) = gP_filt(:, 1, :) ;
            
            % 2d version has periodic row
            Lv2d_filt = Lv_filt ;
            grad2Tr2d_filt = grad2Tr_filt ;
            gvnbM2H2d_filt = gvnbM2H_filt ;
            Kv2d_filt = Kv_filt ;
            gP2d_filt = gP_filt ;
        end
        
        % Average along DV -- ignore last redudant row at nV
        Lv_ap = mean(Lv_filt(:, 1:nV-1, :), 2) ;
        grad2Tr_ap = mean(grad2Tr_filt(:, 1:nV-1, :), 2) ;
        gvnbM2H_ap = mean(gvnbM2H_filt(:, 1:nV-1, :), 2) ;
        Kv_ap = mean(Kv_filt(:, 1:nV-1, :), 2) ;
        gP_ap = mean(gP_filt(:, 1:nV-1, :), 2) ;
        
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
        Lv_l = mean(Lv_filt(:, left, :), 2) ;
        grad2Tr_l = mean(grad2Tr_filt(:, left, :), 2) ;
        gvnbM2H_l = mean(gvnbM2H_filt(:, left, :), 2) ;
        Kv_l = mean(Kv_filt(:, left, :), 2) ;
        gP_l = mean(gP_filt(:, left, :), 2) ;
        
        % right quarter
        Lv_r = mean(Lv_filt(:, right, :), 2) ;
        grad2Tr_r = mean(grad2Tr_filt(:, right, :), 2) ;
        gvnbM2H_r = mean(gvnbM2H_filt(:, right, :), 2) ;
        Kv_r = mean(Kv_filt(:, right, :), 2) ;
        gP_r = mean(gP_filt(:, right, :), 2) ;
        
        % dorsal quarter
        Lv_d = mean(Lv_filt(:, dorsal, :), 2) ;
        grad2Tr_d = mean(grad2Tr_filt(:, dorsal, :), 2) ;
        gvnbM2H_d = mean(gvnbM2H_filt(:, dorsal, :), 2) ;
        Kv_d = mean(Kv_filt(:, dorsal, :), 2) ;
        gP_d = mean(gP_filt(:, dorsal, :), 2) ;
        
        % ventral quarter
        Lv_v = mean(Lv_filt(:, ventral, :), 2) ;
        grad2Tr_v = mean(grad2Tr_filt(:, ventral, :), 2) ;
        gvnbM2H_v = mean(gvnbM2H_filt(:, ventral, :), 2) ;
        Kv_v = mean(Kv_filt(:, ventral, :), 2) ;
        gP_v = mean(gP_filt(:, ventral, :), 2) ;
        
        %% Save timeseries measurements
        disp('Saving timeseries measurements')
        save(Lfn, 'Lv', 'Lv_filt', 'Lv_ap', 'Lv_l', 'Lv_r', 'Lv_d', 'Lv_v')
        save(grad2Trfn, 'grad2Tr', 'grad2Tr_filt', 'grad2Tr_ap', 'grad2Tr_l', 'grad2Tr_r', 'grad2Tr_d', 'grad2Tr_v')
        save(gvnbM2Hfn, 'gvnbM2H', 'gvnbM2H_filt', 'gvnbM2H_ap', 'gvnbM2H_l', 'gvnbM2H_r', 'gvnbM2H_d', 'gvnbM2H_v')
        save(Kfn, 'Kv', 'Kv_filt', 'Kv_ap', 'Kv_l', 'Kv_r', 'Kv_d', 'Kv_v')
        save(Pfn, 'gP', 'gP_filt', 'gP_ap', 'gP_l', 'gP_r', 'gP_d', 'gP_v')
        save_lambdas = true ;
        toc
        
    else
        % Load timeseries measurements
        load(Lfn, 'Lv_filt', 'Lv_ap', 'Lv_l', 'Lv_r', 'Lv_d', 'Lv_v')
        load(grad2Trfn, 'grad2Tr_filt', 'grad2Tr_ap', 'grad2Tr_l', 'grad2Tr_r', 'grad2Tr_d', 'grad2Tr_v')
        load(gvnbM2Hfn, 'gvnbM2H_filt', 'gvnbM2H_ap', 'gvnbM2H_l', 'gvnbM2H_r', 'gvnbM2H_d', 'gvnbM2H_v') 
        load(Kfn, 'Kv_filt', 'Kv_ap', 'Kv_l', 'Kv_r', 'Kv_d', 'Kv_v')
        load(Pfn, 'gP_filt', 'gP_ap', 'gP_l', 'gP_r', 'gP_d', 'gP_v')
        
        Lv2d_filt = Lv_filt ;
        Lv3d_filt = Lv_filt(:, 1:nV-1, :) ;
        grad2Tr2d_filt = grad2Tr_filt ;
        grad2Tr3d_filt = grad2Tr_filt(:, 1:nV-1, :) ;
        gvnbM2H2d_filt = gvnbM2H_filt ;
        gvnbM2H3d_filt = gvnbM2H_filt(:, 1:nV-1, :) ;
        Kv2d_filt = Kv_filt ;
        Kv3d_filt = Kv_filt(:, 1:nV-1, :) ;
        gP2d_filt = gP_filt ;
        gP3d_filt = gP_filt(:, 1:nV-1, :);
        
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
        [grad2Trn, grad2Trt, grad2Trt2dS] = ...
            resolveTangentNormalVertexVelocities(mesh.f, mesh.v, ...
            reshape(grad2Tr3d_filt, [nU*(nV-1), 3]), 1:size(mesh.v, 1), ...
            m2.v, 'vertexnormals', vertexnormals, 'cutPath', cutPath) ;
        grad2Trt2d = grad2Trt2dS.uv2vertices ;
        [gvnbM2Hn, gvnbM2Ht, gvnbM2Ht2dS] = ...
            resolveTangentNormalVertexVelocities(mesh.f, mesh.v, ...
            reshape(gvnbM2H3d_filt, [nU*(nV-1), 3]), 1:size(mesh.v, 1), ...
            m2.v, 'vertexnormals', vertexnormals, 'cutPath', cutPath) ;
        gvnbM2Ht2d = gvnbM2Ht2dS.uv2vertices ;
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
        grad2Trtmag = vecnorm(grad2Trt, 2, 2) ;
        grad2Trtang = atan2(grad2Trt2d(:, 2), grad2Trt2d(:, 1)) ;
        gvnbM2Htmag = vecnorm(gvnbM2Ht, 2, 2) ;
        gvnbM2Htang = atan2(gvnbM2Ht2d(:, 2), gvnbM2Ht2d(:, 1)) ;
        Kvtmag = vecnorm(Kvt, 2, 2) ;
        Kvtang = atan2(Kvt2d(:, 2), Kvt2d(:, 1)) ;
        gPtmag = vecnorm(gPt, 2, 2) ;
        gPtang = atan2(gPt2d(:, 2), gPt2d(:, 1)) ;

        % prep fields
        Lvn2d = extendFieldAlongCutMeshSeam(Lvn, cutMesh) ;
        Lvtmag2d = extendFieldAlongCutMeshSeam(Lvtmag, cutMesh) ;
        Lvtang2d = extendFieldAlongCutMeshSeam(Lvtang, cutMesh) ;
        
        grad2Trn2d = extendFieldAlongCutMeshSeam(grad2Trn, cutMesh) ;
        grad2Trtmag2d = extendFieldAlongCutMeshSeam(grad2Trtmag, cutMesh) ;
        grad2Trtang2d = extendFieldAlongCutMeshSeam(grad2Trtang, cutMesh) ;
        
        gvnbM2Hn2d = extendFieldAlongCutMeshSeam(gvnbM2Hn, cutMesh) ;
        gvnbM2Htmag2d = extendFieldAlongCutMeshSeam(gvnbM2Htmag, cutMesh) ;
        gvnbM2Htang2d = extendFieldAlongCutMeshSeam(gvnbM2Htang, cutMesh) ;
        
        Kvn2d = extendFieldAlongCutMeshSeam(Kvn, cutMesh) ;
        Kvtmag2d = extendFieldAlongCutMeshSeam(Kvtmag, cutMesh) ;
        Kvtang2d = extendFieldAlongCutMeshSeam(Kvtang, cutMesh) ;

        % gP as 2d unraveled
        gPn2d = extendFieldAlongCutMeshSeam(gPn, cutMesh) ;
        gPtmag2d = extendFieldAlongCutMeshSeam(gPtmag, cutMesh) ;
        gPtang2d = extendFieldAlongCutMeshSeam(gPtang, cutMesh) ;
        
        % 3D plot
        if ~exist(fig3dfn, 'file') || overwriteImages
            opts = struct() ;
            opts.labels = {'$\big($d$\star$d$v_{\parallel}^\flat\big)_\perp$', ...
                '$2\nabla(\textrm{Tr}\varepsilon)$', ...
                '$-2\nabla_b v_n(b^{ab}-2Hg^{ab})$', ...
                '$2Kv_\perp$', ...
                '$\big($d$\star$d$v_{\parallel}^\flat\big)_\perp$', ...
                '$2\nabla(\textrm{Tr}\varepsilon)$', ...
                '$-2\nabla_b v_n(b^{ab}-2Hg^{ab})$', ...
                 '$2Kv_\parallel$'} ;
            opts.xyzlim = xyzlim ;
            opts.view = [0, 0] ;
            opts.clim = 40 ;
            opts.polarStyle = 'polar'; 
            opts.axisOff = true ;
            opts.makeCbar = false(8, 1) ;
            opts.masterCbar = true ;
            opts.subplotGrouping = [2, 4] ;
            [axs, cbs, meshHandles] = ...
                nFieldsOnSurface({mesh,mesh,mesh,mesh, mesh,mesh,mesh,mesh}, ...
                {Lvn, grad2Trn, gvnbM2Hn, Kvn, ...
                {Lvtmag, Lvtang}, ...
                {grad2Trtmag, grad2Trtang}, ...
                {gvnbM2Htmag, gvnbM2Htang}, ...
                {Kvtmag, Kvtang}}, opts) ;
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
        nFieldsOnSurface({m2,m2,m2,m2, m2,m2,m2,m2,}, ...
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
