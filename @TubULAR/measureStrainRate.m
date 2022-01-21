function measureStrainRate(QS, options)
%measureStrainRate(QS, options)
%   Compute epsilon = 1/2 (\nabla_i v_j + \nabla_j v_i) - vN b_{ij} 
%   The result on faces is smoothed insofar as the velocities are smoothed
%   with options.lambda and the mesh (governing b_{ij}) is smoothed with
%   options.labmda_mesh. 
%   The result on vertices is additionally smoothed via a quasi-1D spectral 
%   filter governed by options.nmodes and options.zwidth. 
%   
% Parameters
% ----------
% QS : QuapSlap class object instance
% options : struct with fields
%   lambda : float
%   lambda_mesh : float
%   overwrite : bool 
%   overwriteImages : bool
%   preview : bool 
%   averagingStyle : str specifier ('Lagrangian' or 'Simple')
%   samplingResolution : '1x' or '2x'
%   debug : bool 
% 
% Returns 
% -------
% (none)
%
% Saves to disk
% -------------
% strainRate results in QS.dir.strainRate.measurements/strainRate_%06d.mat
%   -> one result per timepoint, with spectral smoothing applied to
%   vertex-based results
% 
% NPMitchell 2020

%% Default options
lambda = QS.smoothing.lambda ;
lambda_mesh = QS.smoothing.lambda_mesh ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;
overwrite = false ;
overwriteImages = false ;
preview = true ;
averagingStyle = 'Lagrangian' ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
debug = false ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
elseif isfield(options, 'overwrite_ims')
    overwriteImages = options.overwrite_ims ;
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
if isfield(options, 'nmodes')
    nmodes = options.nmodes ;
end
if isfield(options, 'zwidth')
    zwidth = options.zwidth ;
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


%% ------------------------------------------------------------------------
% Construct 2D mesh corresponding to the planar domain of
% parameterization
%--------------------------------------------------------------------------
% % Load from file
% path = fullfile(NESpath, 'NES_Examples') ;
% mesh = read_ply_mod(fullfile(path, 'tube_simple_h1p00_R0p00_w1p00.ply')) ;
% rmID = [length(mesh.v)-1, length(mesh.v)] ;
% [F, V] = remove_vertex_from_mesh(mesh.f, mesh.v, rmID) ;
% % Center the mesh around the x axis
% midx = 0.5 * (max(mesh.v(:, 1)) + min(mesh.v(:, 1))) ;
% V(:, 1) = V(:, 1) - midx ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nU = QS.nU ;
nV = QS.nV ;

% Load vertex-based velocity measurements
QS.getVelocityAverage('vv', 'vf')
vertex_vels = QS.velocityAverage.vv ;
face_vels = QS.velocityAverage.vf ;

% Pre-assign timepoints with velocities
tpts = QS.xp.fileMeta.timePoints(1:end-1) ;
tp2do = [tpts(1:10:end), setdiff(tpts, tpts(1:10:end))] ;

% Build metric from mesh
for tp = tp2do
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Load current mesh
    tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
    mesh = tmp.spcutMeshSmRSC ;
    
    % DEBUG
    % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
    % mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) * mesh.ar ;
    mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
    clearvars tmp

    % Define strain rate filename        
    estrainFn = fullfile(strrep(sprintf( ...
        QS.dir.strainRate.measurements, lambda, lambda_mesh, nmodes, zwidth), '.', 'p'), ...
        sprintf('strainRate_%06d.mat', tp));
        
    % Compute the strain rate if not on disk
    if ~exist(estrainFn, 'file') || overwrite
        if exist(estrainFn, 'file')
            disp('Overwriting strain rate on disk')
        else
            disp('Computing strain rate anew')
        end
        
        % Smooth the mesh vertices
        if lambda_mesh > 0 
            tri = triangulation(mesh.f, mesh.v) ;
            fbndy = tri.freeBoundary ;
            fbndy = fbndy(:, 1) ;
            mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda_mesh, 'implicit', mesh.v) ;
            % check smoothed mesh
            % trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none')
        end
        % Define top struct tools
        [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        
        % Smooth the velocities in space using gptoolbox
        if lambda > 0 
            vraw = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, :)) ;
            vs = laplacian_smooth(mesh.v, mesh.f, 'cotan', [], ...
                lambda, 'implicit', vraw) ;   
            % Push vectors onto faces
            vf = V2F * vs ;
        else
            vf = squeeze(face_vels(tidx, :, :)) ;
            
            % check it
            % tmp = squeeze(face_vels(...
            %     tidx, size(mesh.f, 1)+1:2*size(mesh.f, 1), :)) ;
            % 
            % tmpLoad = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
            % tmpMesh = tmpLoad.spcutMeshSmRS ;
            % [TF, TU] = tileAnnularCutMesh(tmpMesh, [1,1]) ;
            % bcs = barycenter(TU, TF) ;
            % scatter(bcs(:, 1), bcs(:, 2), tmp(:, 1))
        end
        
        % %% Checking -- debug
        % if debug
        %     % Obtain mean curvature H for checking against trace(b_ij)
        %     DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
        %     H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
        %     H2d = H3d ;
        %     H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;
        % 
        %     Hf = V2F * H3d ;
        %     divv = DEC.divergence(vf) ;
        %     divf = V2F * divv ;
        %     vss = F2V * vf ;
        % end
        
        %% Convert to 2D mesh
        mesh.nU = QS.nU ;
        cutMesh = cutRectilinearCylMesh(mesh) ;
        
        % Compute the strain rate tensor
        disp('Computing strainRates via covariantDerivative') 
        % tre : traceful dilation
        % dev : deviatoric magnitude
        % theta : angle of elongation
        srmopts = struct() ;
        srmopts.mesh = mesh ;
        srmopts.debug = debug ;
        
        % Compute strain rate by making cutMesh a triple cover so that no
        % boundary effects are present. Periodic in Y --> obtain pathPairs
        % Assumes rectilinear mesh structure.
        pathPairs = [ (1:nU)', (nV-1)*nU + (1:nU)' ] ;
        cutMesh.pathPairs = pathPairs ;
        [strainrate, tre, dev, theta, outStruct] = ...
            strainRatePeriodicMesh(cutMesh, vf, srmopts) ;
        gg = outStruct.fundForms.gg ;
        bb = outStruct.fundForms.bb ;
        dx_faces = outStruct.bondDxDy.dx ;
        dy_faces = outStruct.bondDxDy.dy ;
        theta_pb = outStruct.theta_pb ;
        dvij = outStruct.dvij ;
        
        % Metric strain -- separate trace and deviatoric strain comp, angle
        epsilon_zz = zeros(size(strainrate, 1), 1) ;  % strain rate zeta zeta
        epsilon_zp = zeros(size(strainrate, 1), 1) ;  % strain rate zeta phi
        epsilon_pz = zeros(size(strainrate, 1), 1) ;  % strain rate phi zeta
        epsilon_pp = zeros(size(strainrate, 1), 1) ;  % strain rate phi phi
        g_zz = zeros(size(strainrate, 1), 1) ;  % metric tensor zeta zeta
        g_zp = zeros(size(strainrate, 1), 1) ;  % metric tensor zeta phi
        g_pz = zeros(size(strainrate, 1), 1) ;  % metric tensor phi zeta
        g_pp = zeros(size(strainrate, 1), 1) ;  % metric tensor phi phi
        b_zz = zeros(size(strainrate, 1), 1) ;  % 2nd fund form zeta zeta
        b_zp = zeros(size(strainrate, 1), 1) ;  % 2nd fund form zeta phi
        b_pz = zeros(size(strainrate, 1), 1) ;  % 2nd fund form phi zeta
        b_pp = zeros(size(strainrate, 1), 1) ;  % 2nd fund form phi phi
        % eigv1 = zeros(size(strainrate, 1), 1) ;  % eigenvalue of smaller direction
        % eigv2 = zeros(size(strainrate, 1), 1) ;  % eigenvalue of larger deformation direction (along theta)
        for qq = 1:size(strainrate, 1)
            eq = strainrate{qq} ;
            gq = gg{qq} ;
            bq = bb{qq} ;
            
            %% Full strain rate for mesh averaging onto vertices
            epsilon_zz(qq) = eq(1, 1) ;
            epsilon_zp(qq) = eq(1, 2) ;
            epsilon_pz(qq) = eq(2, 1) ;
            epsilon_pp(qq) = eq(2, 2) ;
            %% Full metric tensor for mesh averaging onto vertices
            g_zz(qq) = gq(1, 1) ;
            g_zp(qq) = gq(1, 2) ;
            g_pz(qq) = gq(2, 1) ;
            g_pp(qq) = gq(2, 2) ;
            %% 2nd fundamental form for mesh averaging onto vertices
            b_zz(qq) = bq(1, 1) ;
            b_zp(qq) = bq(1, 2) ;
            b_pz(qq) = bq(2, 1) ;
            b_pp(qq) = bq(2, 2) ;
        end
                
        %% Collate results as DVavg, L, R, D, V
        % Find trace and deviator on vertices instead of faces
        epsilon_zz_vtx = F2V * epsilon_zz ;
        epsilon_zp_vtx = F2V * epsilon_zp ;
        epsilon_pz_vtx = F2V * epsilon_pz ;
        epsilon_pp_vtx = F2V * epsilon_pp ;
        
        assert(all(epsilon_zp_vtx == epsilon_pz_vtx))
        % Smooth with spectral mode filter if nmodes > 0 
        if nmodes > 0
            filtOpts = struct('nmodes', nmodes, 'zwidth', zwidth) ;
            epsilon_zz_vtx = modeFilterQuasi1D(reshape(epsilon_zz_vtx, [nU, nV-1]), filtOpts) ;
            epsilon_zp_vtx = modeFilterQuasi1D(reshape(epsilon_zp_vtx, [nU, nV-1]), filtOpts) ;
            epsilon_pz_vtx = modeFilterQuasi1D(reshape(epsilon_pz_vtx, [nU, nV-1]), filtOpts) ;
            epsilon_pp_vtx = modeFilterQuasi1D(reshape(epsilon_pp_vtx, [nU, nV-1]), filtOpts) ;
            epsilon_zz_vtx = epsilon_zz_vtx(:) ;
            epsilon_zp_vtx = epsilon_zp_vtx(:) ;
            epsilon_pz_vtx = epsilon_pz_vtx(:) ;
            epsilon_pp_vtx = epsilon_pp_vtx(:) ;
        end
        
        g_zz_vtx = F2V * g_zz ;
        g_zp_vtx = F2V * g_zp ;
        g_pz_vtx = F2V * g_pz ;
        g_pp_vtx = F2V * g_pp ;
        b_zz_vtx = F2V * b_zz ;
        b_zp_vtx = F2V * b_zp ;
        b_pz_vtx = F2V * b_pz ;
        b_pp_vtx = F2V * b_pp ;
        dx_vtx = F2V * dx_faces ;
        dy_vtx = F2V * dy_faces ;
        
        % pre-allocate theta and theta_pullback
        theta_vtx = 0 * epsilon_zz_vtx ;
        theta_pb_vtx = 0 * epsilon_zz_vtx ;
        for qq = 1:size(epsilon_zz_vtx, 1)
            %% Traceful dilation
            eq = [epsilon_zz_vtx(qq), epsilon_zp_vtx(qq); ...
                  epsilon_pz_vtx(qq), epsilon_pp_vtx(qq)] ;
            gq = [g_zz_vtx(qq), g_zp_vtx(qq); ...
                  g_pz_vtx(qq), g_pp_vtx(qq)] ;
            
            % traceful component -- 1/2 Tr[g^{-1} gdot] = Tr[g^{-1} eps] 
            try
                [tre_vtx(qq), dev_vtx(qq), ...
                    theta_vtx(qq), theta_pb_vtx(qq)] = ...
                    traceDeviatorPullback(eq, gq, dx_vtx(qq), dy_vtx(qq)) ;
            catch
                error('here')
                
            end
        end
        % modulo is not necessary
        % theta_vtx = mod(theta_vtx, pi) ;
        
        %% Store measurements on vertices in grouped arrays
        strainrate_vtx = [epsilon_zz_vtx, epsilon_zp_vtx, ...
            epsilon_pz_vtx, epsilon_pp_vtx] ;
        gg_vtx = [g_zz_vtx, g_zp_vtx, ...
            g_pz_vtx, g_pp_vtx] ;
        bb_vtx = [b_zz_vtx, b_zp_vtx, ...
            b_pz_vtx, b_pp_vtx] ;
        tre_vtx((nU * (nV-1) + 1):nU*nV) = tre_vtx(1:nU) ;
        dev_vtx((nU * (nV-1) + 1):nU*nV) = dev_vtx(1:nU) ;
        theta_vtx((nU * (nV-1) + 1):nU*nV) = theta_vtx(1:nU) ;
        tre_vtx = reshape(tre_vtx, [nU,nV]) ;
        dev_vtx = reshape(dev_vtx, [nU,nV]) ;
        theta_vtx = reshape(theta_vtx, [nU,nV]) ;
        
        % Average along DV -- ignore last redudant row at nV
        [dev_ap, theta_ap] = ...
            QS.dvAverageNematic(dev_vtx(:, 1:nV-1), theta_vtx(:, 1:nV-1)) ;
        tre_ap = mean(tre_vtx(:, 1:nV-1), 2) ;
        
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
        [dev_l, theta_l] = ...
            QS.dvAverageNematic(dev_vtx(:, left), theta_vtx(:, left)) ;
        tre_l = mean(tre_vtx(:, left), 2) ;
        
        % right quarter
        [dev_r, theta_r] = ...
            QS.dvAverageNematic(dev_vtx(:, right), theta_vtx(:, right)) ;
        tre_r = mean(tre_vtx(:, right), 2) ;
        
        % dorsal quarter
        [dev_d, theta_d] = ...
            QS.dvAverageNematic(dev_vtx(:, dorsal), theta_vtx(:, dorsal)) ;
        tre_d = mean(tre_vtx(:, dorsal), 2) ;
        
        % ventral quarter
        [dev_v, theta_v] = ...
            QS.dvAverageNematic(dev_vtx(:, ventral), theta_vtx(:, ventral)) ;
        tre_v = mean(tre_vtx(:, ventral), 2) ;
        
        % save the metric strain
        readme.strainrate = 'strain rate on faces, epsilon=1/2(nabla_i v_j + nabla_j v_i) - vn b_ij' ;
        readme.tre = 'Tr[g^{-1} epsilon] scalar field on faces -- multiply by 1/2 to compare to dev';
        readme.dev = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) magnitude, scalar field on faces';
        readme.theta = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector in embedding space';
        readme.dvij = 'symmetrized covariant derivative 0.5 * (d_i v_j + d_j v_i)';
        readme.gg = 'metric tensor on faces';
        readme.bb = 'second fundamental form on faces';
        readme.lambda = 'Laplacian smoothing on velocities' ;
        readme.lambda_mesh = 'Laplacian smoothing on mesh vertices' ;
        readme.strainrate_vtx = 'strain rate on vertices, epsilon=1/2(nabla_i v_j + nabla_j v_i) - vn b_ij' ;
        readme.gg_vtx = 'metric tensor on vertices';
        readme.bb_vtx = 'second fundamental form on vertices';
        readme.tre_vtx = 'Tr[g^{-1} epsilon], on mesh vertices' ;
        readme.dev_vtx = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) on mesh vertices';
        readme.theta_vtx = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector on mesh vertices';
        readme.dev_ap = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged circumferentially';
        readme.dev_l = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on left quarter, on vertices';
        readme.dev_r = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on right quarter, on vertices';
        readme.dev_d = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on dorsal quarter, on vertices';
        readme.dev_v = 'sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on ventral quarter, on vertices';
        readme.theta_ap = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged circumferentially, on vertices';
        readme.theta_l = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on left quarter, on vertices';
        readme.theta_r = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on right quarter, on vertices';
        readme.theta_d = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on dorsal quarter, on vertices';
        readme.theta_v = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on ventral quarter, on vertices';
        readme.tre_ap = 'Tr[g^{-1} epsilon], averaged circumferentially, on vertices';
        readme.tre_l = 'Tr[g^{-1} epsilon], averaged on left quarter, on vertices';
        readme.tre_r = 'Tr[g^{-1} epsilon], averaged on right quarter, on vertices';
        readme.tre_d = 'Tr[g^{-1} epsilon], averaged on dorsal quarter, on vertices';
        readme.tre_v = 'Tr[g^{-1} epsilon], averaged on ventral quarter, on vertices';
        readme.note = 'The pullback space is taken to range from zeta=[0, 1] and phi=[0, 1]' ; 
        disp(['saving ', estrainFn])
        save(estrainFn, 'strainrate', 'tre', 'dev', 'theta', 'theta_pb', ...
            'dvij', 'gg', 'bb', 'lambda', 'lambda_mesh', 'readme', ...
            'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v', ...
            'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
            'strainrate_vtx', 'tre_vtx', 'dev_vtx', 'theta_vtx', ...
            'gg_vtx', 'bb_vtx', 'dx_faces', 'dy_faces', 'dx_vtx', 'dy_vtx')
        
        % Save info histogram as debug check
        % clf
        % plot(dev(:) .* cos(theta(:)), dev(:) .* sin(theta(:)), '.')
        % xlabel('deviator$[\varepsilon]_\zeta$', 'interpreter', 'latex')
        % ylabel('deviator$[\varepsilon]_\phi$', 'interpreter', 'latex')  
        % axis equal
        % xlim([-0.1, 0.1])
        % ylim([-0.1, 0.1])
        % saveas(gcf, fullfile(strrep(sprintf( ...
        %     QS.dir.strainRate.measurements, lambda, lambda_mesh), '.', 'p'), ...
        %     sprintf('strainRate_check_%06d.png', tp)));
    else
        % Convert to 2D mesh
        mesh.nU = QS.nU ;
        cutMesh = cutRectilinearCylMesh(mesh) ;
    end        
    
    % Plot the result
    options.overwrite = overwriteImages ;
    options.mesh = mesh ;
    options.cutMesh = cutMesh ;
    options.lambda = lambda ;
    options.lambda_mesh = lambda_mesh ;
    options.nmodes = nmodes ;
    options.zwidth = zwidth ;
    options.debug = debug ;
    QS.plotStrainRateTimePoint(tp, options)
end
disp('done with measuring strain rates')

