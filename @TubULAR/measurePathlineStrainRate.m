function measurePathlineStrainRate(QS, options)
% measurePathlineStrainRate(QS, options)
%   Query the metric strain rate along lagrangian pathlines.
%   Plot results as kymographs.
%   
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs : bool
%   plot_kymographs_cumsum : bool
%   plot_gdot_correlations : bool
%   plot_gdot_decomp : bool
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
overwriteImages = false ;
plot_comparison = false ;

%% Parameter options
lambda_mesh = QS.smoothing.lambda_mesh ;
lambda = QS.smoothing.lambda ; 
debug = false ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = "Lagrangian" ;
% Load time offset for first fold, t0 -- default pathline t0
QS.t0set() ;
t0 = QS.t0 ;
% By default, t0Pathline = t0 (see below)

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'plot_comparison')
    plot_comparison = options.plot_comparison ;
end
%% parameter options
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
else
    % default lambda_mesh is equal to lambda 
    lambda_mesh = lambda ;
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
if strcmp(averagingStyle, 'Lagrangian')
    sKDir = QS.dir.strainRate.smoothing;
else
    error('Have not implemented strain rate measurements based on simple averaging')
end
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
bwr256 = bluewhitered(256) ;

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - t0;

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];
vunitstr = [ '[' QS.spaceUnits '/' QS.timeUnits ']' ];
    
% DONE WITH PREPARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pathlines to build Kymographs along pathlines
QS.loadPullbackPathlines(t0Pathline, 'vertexPathlines')
vP = QS.pathlines.vertices ;

% Output directory is inside StrainRate dir
sKPDir = fullfile(sKDir, sprintf('pathline_%04dt0', t0Pathline)) ;
outdir = fullfile(sKPDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
% Data for kinematics on meshes (defined on vertices)
mdatdir = fullfile(sKDir, 'measurements') ;

% Load Lx, Ly by loadingPIV. 
QS.loadPIV()
Xpiv = QS.piv.raw.x ;
Ypiv = QS.piv.raw.y ;

% Also need velocities to advect mesh
% QS.loadVelocityAverage('vv')
% vPIV = QS.velocityAverage.vv ;

% Discern if piv measurements are done on a double covering or the meshes
if strcmp(QS.piv.imCoords(end), 'e')
    doubleCovered = true ;
end

% Compute or load all timepoints
strain = zeros(size(vP.vX, 2) * size(vP.vX, 3), 4) ; 
strain(:, 2) = 1e4 ;
strain(:, 3) = 1e4 ;
for tp = QS.xp.fileMeta.timePoints(1:end-1)
    close all
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;
    QS.setTime(tp) ;
    
    % Check for timepoint measurement on disk, on mesh vertices 
    estrainFn = fullfile(outdir, sprintf('strainRate_%06d.mat', tp)) ;
    
    if overwrite || ~exist(estrainFn, 'file') 
        % Load timeseries measurements defined on mesh vertices 
        srfnMesh = fullfile(mdatdir, sprintf('strainRate_%06d.mat', tp)) ;
        try
            load(srfnMesh, 'strainrate_vtx', 'gg_vtx', 'dx_vtx', 'dy_vtx') 
        catch
            msg = 'Run QS.measurePathlineStrainRate() ' ;
            msg = [msg 'with lambdas=(mesh,lambda,err)=('] ;
            msg = [msg num2str(lambda_mesh) ','] ;
            msg = [msg num2str(lambda) ','] ;
            msg = [msg ' before running ', ...
                    'QS.measurePathlineStrainRate()'] ;
                
            load(srfnMesh, 'strainrate_vtx', 'gg_vtx', 'dx_vtx', 'dy_vtx') 
            error(msg)
        end
        %% Interpolate StrainRate from vertices onto pathlines
        xx = vP.vX(tidx, :, :) ;
        yy = vP.vY(tidx, :, :) ;
        XY = [xx(:), yy(:)] ;
        Lx = vP.Lx(tidx) ;
        Ly = vP.Ly(tidx) ;
        options.Lx = Lx ;
        options.Ly = Ly ;
        XY = QS.doubleToSingleCover(XY, Ly) ;
        
        %% Recall strain rate at gridded vertices
        % strainrate from vertices to pathlines
        exx = strainrate_vtx(:, 1) ;
        exy = strainrate_vtx(:, 2) ;
        eyx = strainrate_vtx(:, 3) ;
        eyy = strainrate_vtx(:, 4) ;
        exx(nU*(nV-1)+1:nU*nV) = exx(1:nU) ;
        exy(nU*(nV-1)+1:nU*nV) = exy(1:nU) ;
        eyx(nU*(nV-1)+1:nU*nV) = eyx(1:nU) ;
        eyy(nU*(nV-1)+1:nU*nV) = eyy(1:nU) ;
        ezz = QS.interpolateOntoPullbackXY(XY, exx, options) ;
        ezp = QS.interpolateOntoPullbackXY(XY, exy, options) ;
        epz = QS.interpolateOntoPullbackXY(XY, eyx, options) ;
        assert(all(ezp == epz))
        epp = QS.interpolateOntoPullbackXY(XY, eyy, options) ;
        strainrate = [ezz, ezp, epz, epp] ;
        
        % Metric from vertices to pathlines
        gxx = gg_vtx(:, 1) ;
        gxy = gg_vtx(:, 2) ;
        gyx = gg_vtx(:, 3) ;
        gyy = gg_vtx(:, 4) ;
        gxx(nU*(nV-1)+1:nU*nV) = gxx(1:nU) ;
        gxy(nU*(nV-1)+1:nU*nV) = gxy(1:nU) ;
        gyx(nU*(nV-1)+1:nU*nV) = gyx(1:nU) ;
        gyy(nU*(nV-1)+1:nU*nV) = gyy(1:nU) ;
        gzz = QS.interpolateOntoPullbackXY(XY, gxx, options) ;
        gzp = QS.interpolateOntoPullbackXY(XY, gxy, options) ;
        gpz = QS.interpolateOntoPullbackXY(XY, gyx, options) ;
        gpp = QS.interpolateOntoPullbackXY(XY, gyy, options) ;
        assert(all(gzp == gpz))
        gg = [gzz, gzp, gpz, gpp] ;
        
        % dx and dy scalings
        dx = dx_vtx ;
        dy = dy_vtx ;
        dx(nU*(nV-1)+1:nU*nV) = dx(1:nU) ;
        dy(nU*(nV-1)+1:nU*nV) = dy(1:nU) ;
        dz = QS.interpolateOntoPullbackXY(XY, dx, options) ;
        dp = QS.interpolateOntoPullbackXY(XY, dy, options) ;
        
        %% Save image of dz and dp
        % dzdpDir = fullfile(sprintf(QS.dir.strainRate.pathline.root, ...
        %         lambda, lambda_mesh), 'dzdp') ;
        % if plot_dzdp && (~exist(
        %     % Check the ratio of lengths compared to pullback lengths
        %     subplot(2, 2, 1)
        %     trisurf(triangulation(mesh1.f, mesh1.v), dz, 'edgecolor', 'none')
        %     title('d$\zeta^{3D}/$d$\zeta^{2D}$', 'interpreter', 'latex')
        %     axis equal
        %     subplot(2, 2, 2)
        %     trisurf(triangulation(mesh1.f, mesh1.v), dp, 'edgecolor', 'none')
        %     title('d$\phi^{3D}/$d$\phi^{2D}$', 'interpreter', 'latex')
        %     axis equal
        %     subplot(2, 2, 3)
        %     meandz = mean(reshape(dz, [nU,nV]), 2) ;
        %     plot((1:nU)/nU, meandz, '.')
        %     if all(meandz < 600)
        %         ylim([0, 800])
        %     end
        %     ylabel('$\langle$d$\zeta\rangle_{dv}$', 'interpreter', 'latex')
        %     xlabel('ap position, $\zeta$', 'interpreter', 'latex')
        %     % save this test image
        %     saveas(gcf, fullfile(dzdpDir, sprintf('dzdp_%06d.png', tidx))) ;
        % end
        
        %% Compute the strain rate trace and deviator at the pathlines
        for qq = 1:size(exx, 1)
            eq = [ezz(qq), ezp(qq); ...
                  epz(qq), epp(qq)] ;
            gq = [gzz(qq), gzp(qq); ...
                  gpz(qq), gpp(qq)] ;
            
            % traceful component is 1/2 Tr[g^{-1} gdot] = Tr[g^{-1} eps] 
            % deviatoric component is eps - 1/2 (traceful component) g
            [tre(qq), dev(qq), theta(qq)] = ...
                traceDeviatorPullback(eq, gq, dz(qq), dp(qq)) ;
        end
                
        %% OPTION 1: simply reshape, tracing each XY pathline pt to its t0
        % % grid coordinate
        tre = reshape(tre, [nU, nV]) ;
        dev = reshape(dev, [nU, nV]) ;
        theta = reshape(theta, [nU, nV]) ;
        
        %% Average strainRATE along pathline DV hoops
        % Average along DV -- ignore last redudant row at nV
        [dev_ap, theta_ap] = ...
            QS.dvAverageNematic(dev(:, 1:nV-1), theta(:, 1:nV-1)) ;
        tre_ap = mean(tre(:, 1:nV-1), 2) ;
        
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
            QS.dvAverageNematic(dev(:, left), theta(:, left)) ;
        tre_l = mean(tre(:, left), 2) ;
        
        % right quarter
        [dev_r, theta_r] = ...
            QS.dvAverageNematic(dev(:, right), theta(:, right)) ;
        tre_r = mean(tre(:, right), 2) ;
        
        % dorsal quarter
        [dev_d, theta_d] = ...
            QS.dvAverageNematic(dev(:, dorsal), theta(:, dorsal)) ;
        tre_d = mean(tre(:, dorsal), 2) ;
        
        % ventral quarter
        [dev_v, theta_v] = ...
            QS.dvAverageNematic(dev(:, ventral), theta(:, ventral)) ;
        tre_v = mean(tre(:, ventral), 2) ;
        
        %% Check result
        if debug
            close all
            scatter(1:nU, strain_th_ap, 20 * strain_dv_ap / max(strain_dv_ap),...
                strain_dv_ap, 'filled')
            pause(2)
        end
        
        % save the metric strain
        readme.strainrate = 'strain rate tensor interpolated onto pathlines' ;
        readme.gg = 'metric tensor interpolated onto pathlines' ;
        readme.tre = 'Tr[g^{-1} epsilon]';
        readme.dev = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) on mesh vertices';
        readme.theta = 'arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector on mesh vertices';
        readme.dev_ap = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged circumferentially';
        readme.dev_l = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged on left quarter, on vertices';
        readme.dev_r = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged on right quarter, on vertices';
        readme.dev_d = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged on dorsal quarter, on vertices';
        readme.dev_v = 'sqrt( Tr[g^{-1} deviator[epsilon] g^{-1} deviator[epsilon]] ) averaged on ventral quarter, on vertices';
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
        
        readme.note = 'Evaluated for Lagrangian paths. The pullback space is taken to range from zeta=[0, 1] and phi=[0, 1]' ; 
        disp(['saving ', estrainFn])
        save(estrainFn, 'strainrate', 'gg', 'dz', 'dp', 'readme', ...
            'tre', 'dev', 'theta', ...
            'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
            'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v')
    else
        disp('strainRate already on disk')
        % Load mesh
        % disp('Loading mesh')
        % tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
        % mesh1 = tmp.spcutMeshSmRS ;
        % clearvars tmp
        % % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
        % mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) * mesh1.ar ;
    end
    
    % Plot the result
    plotOpts.overwrite = overwriteImages ;
    plotOpts.cutMesh = [] ;
    plotOpts.lambda = lambda ;
    plotOpts.lambda_mesh = lambda_mesh ;
    plotOpts.debug = debug ;
    plotOpts.t0Pathline = t0Pathline ;
    plotOpts.plot_comparison = plot_comparison ;
    disp('plotting pathline strain rate for this timepoint')
    QS.plotPathlineStrainRateTimePoint(tp, plotOpts)
end
disp('done with measuring pathline strain rate')


%% Combine DV-averaged profiles into kymographs
apKymoFn = fullfile(outdir, 'apKymographPathlineStrainRate.mat') ;
lKymoFn = fullfile(outdir, 'leftKymographPathlineStrainRate.mat') ;
rKymoFn = fullfile(outdir, 'rightKymographPathlineStrainRate.mat') ;
dKymoFn = fullfile(outdir, 'dorsalKymographPathlineStrainRate.mat') ;
vKymoFn = fullfile(outdir, 'ventralKymographPathlineStrainRate.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if ~files_exist 
    disp('Compiling kymograph data to save to disk...')
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        close all
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        srfn = fullfile(outdir, sprintf('strainRate_%06d.mat', tp))   ;

        % Load timeseries measurements
        load(srfn, 'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
            'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v') ;

        %% Store in matrices
        % dv averaged
        tr_apM(tidx, :) = tre_ap ;
        dv_apM(tidx, :) = dev_ap ;
        th_apM(tidx, :) = theta_ap ;

        % left quarter
        tr_lM(tidx, :) = tre_l ;
        dv_lM(tidx, :) = dev_l ;
        th_lM(tidx, :) = theta_l ;

        % right quarter
        tr_rM(tidx, :) = tre_r ;
        dv_rM(tidx, :) = dev_r ;
        th_rM(tidx, :) = theta_r ;

        % dorsal quarter
        tr_dM(tidx, :) = tre_d ;
        dv_dM(tidx, :) = dev_d ;
        th_dM(tidx, :) = theta_d ;

        % ventral quarter
        tr_vM(tidx, :) = tre_v ;
        dv_vM(tidx, :) = dev_v ;
        th_vM(tidx, :) = theta_v ;

    end
    
    %% Save kymographs
    disp('Saving kymograph data files for Lagrangian pathlines')
    save(apKymoFn, 'tr_apM', 'dv_apM', 'th_apM')
    disp(['Saved kymograph data to: ' apKymoFn])
    save(lKymoFn, 'tr_lM', 'dv_lM', 'th_lM')
    disp(['Saved kymograph data to: ' lKymoFn])
    save(rKymoFn, 'tr_rM', 'dv_rM', 'th_rM')
    disp(['Saved kymograph data to: ' rKymoFn])
    save(dKymoFn, 'tr_dM', 'dv_dM', 'th_dM')
    disp(['Saved kymograph data to: ' dKymoFn])
    save(vKymoFn, 'tr_vM', 'dv_vM', 'th_vM')
    disp(['Saved kymograph data to: ' vKymoFn])
    disp('done with strainRate kymograph data saving')
else
    disp('strainRate kymograph data already on disk')    
end
    