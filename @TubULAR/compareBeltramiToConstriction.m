function compareBeltramiToConstriction(QS, options)
% Use results of QS.measureBeltramiCoefficient() to compare quasiconformal
% Beltrami coefficient to indentation
%
% Parameters
% ----------
%
% Returns
% -------
%
% NPMitchell 2020

%% Default options
coordSys = 'ricci' % could be 'uvprime' or 'ricci'
if isempty(QS.t0)
    t0Pathlines = QS.t0set() ;
else
    t0Pathlines = QS.t0 ;
end
t0 = QS.t0 ;
overwrite = false ;
if nargin < 2 
    options = struct() ;
end
climit = 0.6 ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end
if isfield(options, 't0Pathlines')
    t0Pathlines = options.t0Pathlines ;
else
    disp('Using default t0 for pathlines')
end
if isfield(options, 'climit')
    climit = options.climit ;
end

% Load pathline velocities in UVPrime coordinates
if strcmpi(coordSys, 'ricci')
    v3d = load(sprintf(QS.fileName.pathlines.v3d, t0Pathlines));
    load(sprintf(QS.fileName.pathlines.refMesh, t0Pathlines), 'refMesh') ;
elseif strcmpi(coordSys, 'uvprime')
    v3d = load(sprintf(QS.fileName.pathlines_uvprime.v3d, t0Pathlines));
    load(sprintf(QS.fileName.pathlines_uvprime.refMesh, t0Pathlines), 'refMesh') ;
end
v3x = v3d.v3dPathlines.vXrs ;
v3y = v3d.v3dPathlines.vYrs ;
v3z = v3d.v3dPathlines.vZrs ;

% Load refMesh
u2d = [refMesh.u(:, 1)/max(refMesh.u(:, 1)), refMesh.u(:, 2), 0*refMesh.u(:, 2)] ;
mesh2 = struct();
mesh2.f = refMesh.f ;
mesh2.v = u2d ;
refMesh2plot = mesh2 ;

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Compute dzeta -- material frame distance between DV hoops
%     % (dzeta along the surface from each hoop to the next)
%     % The distance from one hoop to another is the
%     % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
%     refMesh.dzeta = reshape(vecnorm(...
%         diff(reshape(refMesh.vrs, [refMesh.nU, refMesh.nV, 3])), 2, 3), ...
%         [refMesh.nU-1, refMesh.nV]) ;
%     refMesh.dzeta_mean = nanmean(refMesh.dzeta, 2) ;
%         
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % save refMesh
%     save(sprintf(QS.fileName.pathlines.refMesh, t0Pathlines), 'refMesh')
    
umax_ss = sum(refMesh.dzeta_mean) ;

% Get xyzlim
[~,~,~,xyzlim] = QS.getXYZLims() ;

indent_computed = false ;
loaded_mu = false ;

recompute = true ;
if ~overwrite
    try
        error('')
        recompute = false ;
    catch 
        disp('Could not load from disk, computing...')
    end
end

if recompute        
    % Compute full expression for mu and plot in 3d and 2d 
    % Load radii for vertices of pathline meshes at all timepoints

    if strcmpi(coordSys, 'ricci')
        radfn = sprintf(QS.fileName.pathlines.radius, t0Pathlines) ;
        v3dfn = sprintf(QS.fileName.pathlines.v3d, t0Pathlines) ;
        reffn = sprintf(QS.fileName.pathlines.refMesh, t0Pathlines) ;
        mu_all = load(sprintf(QS.fileName.pathlines.quasiconformal, t0Pathlines)) ;
    elseif strcmpi(coordSys, 'uvprime')
        radfn = sprintf(QS.fileName.pathlines_uvprime.radius, t0Pathlines) ;
        v3dfn = sprintf(QS.fileName.pathlines_uvprime.v3d, t0Pathlines) ;
        reffn = sprintf(QS.fileName.pathlines_uvprime.refMesh, t0Pathlines) ;
        fntmp = fullfile(sprintf(...
            QS.dir.pathlines_uvprime.quasiconformal, t0Pathlines), ...
            'mu_mesh.png') ;
        mu_all = load(sprintf(QS.fileName.pathlines_uvprime.quasiconformal, t0Pathlines)) ;
    end
    
    load(radfn, 'vRadiusPathlines') ;
    % Load vertices in 3d for pathline meshes at all timepoints
    load(v3dfn, 'v3dPathlines') ;
    % Load refMesh
    load(reffn, 'refMesh') ;
    mu = bc_metric(refMesh.f, refMesh.u_ricci, refMesh.v, 3) ;
    
    % Check original mu
    close all
    m2plot = triangulation(refMesh.f, ...
        [refMesh.u(:, 1)/max(refMesh.u(:, 1)), ...
        refMesh.u(:, 2), 0*refMesh.u(:,1)]) ;
    opts.view = [0, 90] ;
    opts.labels = {'$\Re\mu$', '$\Im \mu$'} ;
    [axs, cbs, meshHandles] =  ...
        nFieldsOnSurface(m2plot, {real(mu), imag(mu)}, opts) ;
    for qq = 1:2
        set(gcf,'CurrentAxes',axs{qq})
        xlabel('$u$', 'interpreter', 'latex') ;
        ylabel('$v$', 'interpreter', 'latex') ;
        axis equal; caxis(max(abs(real(mu))) * [-1, 1])
        colormap blueblackred
        grid off; axis tight ;
    end
    sgtitle(['$\mu$ for reference mesh'], 'Interpreter', 'latex')
    fntmp = fullfile(sprintf(...
        QS.dir.pathlines.quasiconformal, t0Pathlines), 'mu_t0.png') ;
    saveas(gcf, fntmp)
        
    % Load true mu Beltrami coefficient for comparison
    mu_all = mu_all.mu_material_filtered ;
    nU = refMesh.nU ;
    nV = refMesh.nV ;
    
    % Create output directories
    if strcmpi(coordSys, 'ricci')
        predDir3d = fullfile(sprintf(...
            QS.dir.pathlines.quasiconformal, t0Pathlines), ...
            'images_pred3d') ;
        predDir2d = fullfile(sprintf(...
            QS.dir.pathlines.quasiconformal, t0Pathlines), ...
            'images_pred2d') ;
    else
        predDir3d = fullfile(sprintf(...
            QS.dir.pathlines_uvprime.quasiconformal, t0Pathlines), ...
            'images_pred3d') ;
        predDir2d = fullfile(sprintf(...
            QS.dir.pathlines_uvprime.quasiconformal, t0Pathlines), ...
            'images_pred2d') ;
    end
    if ~exist(predDir3d, 'dir')
        mkdir(predDir3d)
    end
    if ~exist(predDir2d, 'dir')
        mkdir(predDir2d)
    end
    
    % Get radius at t=0
    tidx0 = QS.xp.tIdx(t0Pathlines) ;
    radv0 = reshape(vRadiusPathlines.radii(tidx0, :, :), [nU*nV, 1]) ;
    
    % Consider all timepoints
    timePoints = QS.xp.fileMeta.timePoints ;
    firstSet = 1:10:length(timePoints) ;
    tidx2do = [160, firstSet] ; %, setdiff(1:length(timePoints), firstSet)] ;
    for tidx = tidx2do
        tp = QS.xp.fileMeta.timePoints(tidx) ;

        %% Mesh for this timepoint is made from vertex pathlines
        mesh = refMesh ;
        xx = v3dPathlines.vXrs(tidx, :, :) ;
        yy = v3dPathlines.vYrs(tidx, :, :) ;
        zz = v3dPathlines.vZrs(tidx, :, :) ;
        mesh.v = [xx(:), yy(:), zz(:)] ;
        
        %% Compute rzeta and rphi, as is done in measurePathlineStrain
        % Find dx and dy -- lengths of projected / lengths of pullback
        % [~, dbonds, ~] = labelRectilinearMeshBonds(mesh) ;
        % dx = vecnorm(dbonds.realSpace.u, 2, 2) ./ ...
        %     vecnorm(dbonds.baseSpace.u, 2, 2) ;
        % dy = vecnorm(dbonds.realSpace.v, 2, 2) ./ ...
        %     vecnorm(dbonds.baseSpace.v, 2, 2) ;
        
        % Finite gradient operator
        radv = reshape(vRadiusPathlines.radii(tidx, :, :), [nU*nV, 1]) ;
        radv = radv ./ radv0 ;
        
        dec = DiscreteExteriorCalculus(refMesh.f, ...
            [refMesh.u_ricci(:, 1), refMesh.u_ricci(:, 2), 0*refMesh.u_ricci(:, 1)]) ;
        dxdy = dec.gradient(radv(:)) ;
        dx = dxdy(:, 1) ;
        dy = dxdy(:, 2) ;
        
        % Average radii onto faces
        glueMesh = glueCylinderCutMeshSeam(mesh) ;
        [V2F, F2V] = meshAveragingOperators(glueMesh.f, glueMesh.v) ;
        % radii_glued = vRadiusPathlines.radii(tidx, :, 1:end-1) ;
        dxv = real(F2V * dx) ;
        dyv = real(F2V * dy) ;
        dxv(nU*(nV-1)+1:nU*nV) = dxv(1:nU) ;
        dyv(nU*(nV-1)+1:nU*nV) = dyv(1:nU) ;
        
        % RHS mu on faces
        % radf = V2F * radii_glued(:) ;
        % num = 1 - radf.^2 + dxf.^2 - dyf.^2 + 2 * 1j * dxf .* dyf ;
        % denom = 1 + radf.^2 + dxf.^2 + dyf.^2 + 2*sqrt(radf.^2 .*(1+dxf.^2) + dyf.^2) ;
        % mu_pred_faces = num ./ denom ;
        
        % RHS mu on vertices
        num = 1 - radv.^2 + dxv.^2 - dyv.^2 + 2.0 * 1j * dxv .* dyv ;
        denom = 1 + radv.^2 + dxv.^2 + dyv.^2 + 2*sqrt(radv.^2 .*(1+dxv.^2) + dyv.^2) ;
        mu_pred_vtx = num ./ denom ;
        % test = (1 - radv.^2) ./ (1 + radv.^2 + 2*radv) ;
        
        %% Check dxv and dyv
        % trisurf(triangulation(meshScaled.f, meshScaled.v), dxv, 'edgecolor', 'none')
        % caxis([-max(abs(dxv)), max(abs(dxv))])
        % colormap bwr; colorbar;
        % axis equal
        % 
        % trisurf(triangulation(meshScaled.f, meshScaled.v), dyv, 'edgecolor', 'none')
        % caxis([-max(abs(dyv)), max(abs(dyv))])
        % colormap bwr; colorbar;
        % colorbar()
        % axis equal
        
        %% Check the result
        subplot(2,3,1)
        imagesc(reshape(dxv, [nU,nV]))
        caxis([-max(abs(dxv)), max(abs(dxv))])
        title('$r_\zeta$', 'interpreter', 'latex') ; colorbar() ; axis off
        
        subplot(2,3,2)
        imagesc(reshape(dyv, [nU,nV]))
        caxis([-max(abs(dyv)), max(abs(dyv))])
        title('$r_\phi$', 'interpreter', 'latex') ; colorbar() ; axis off
        
        subplot(2,3,3)
        imagesc(reshape(radv, [nU,nV]))
        caxis([-max(abs(radv)), max(abs(radv))])
        title('$r$', 'interpreter', 'latex') ; colorbar() ; axis off
        
        subplot(2,3,4)
        imagesc(reshape(real(num), [nU,nV]))
        caxis([-max(abs(num)), max(abs(num))])
        title('numerator', 'interpreter', 'latex') ; colorbar() ; axis off
        
        subplot(2,3,5)
        imagesc(reshape(real(denom), [nU,nV]))
        caxis([-max(abs(denom)), max(abs(denom))])
        title('denom', 'interpreter', 'latex') ; colorbar() ; axis off
        
        subplot(2,3,6)
        imagesc(reshape(real(mu_pred_vtx), [nU,nV]))
        set(gcf, 'visible', 'on')
        colormap bwr
        title('$\mu$', 'interpreter', 'latex') ; colorbar() ; axis off
        caxis(max(abs(real(mu_pred_vtx))) * [-1, 1])
        pause(0.1)
        close all
        
        % Mode filter prediction
        filterOptions.widthX = 3 ;
        filterOptions.nmodesY = 5 ;
        filterOptions.preview = false ;
        % truncate final u row for mode filtering
        mu_pred_vtx = mu_pred_vtx(1:nU*(nV-1)) ;
        mu_pred_vtx = reshape(mu_pred_vtx, [nU, nV-1]) ;
        muMVfilt_re = modeFilterQuasi1D(real(mu_pred_vtx), filterOptions) ;
        muMVfilt_im = modeFilterQuasi1D(imag(mu_pred_vtx), filterOptions) ;
        mu_pred_vtx = muMVfilt_re(:) + 1j * muMVfilt_im(:) ;
        mu_pred_vtx(nU*(nV-1)+1:nU*nV) = mu_pred_vtx(1:nU) ;
        
        % Difference between prediction and computed value
        diff_mu = mu_pred_vtx - reshape(mu_all(tidx, :), size(mu_pred_vtx)) ;
        % diff_mu(nU*(nV-1)+1:nU*nV) = diff_mu(1:nU) ;
        
        plotOpts.clim = [-0.5, 0.5] ;
        plotOpts.edgecolor = 'none' ;
        plotOpts.labels = {'$\Re \mu_{\textrm{pred}}$', ...
            '$\Re\mu$', ...
            '$\Re(\mu_{\textrm{pred}} - \mu)$', ...
            '$\Im\mu_{\textrm{pred}}$', ...
            '$\Im\mu$', ...
            '$\Im(\mu_{\textrm{pred}} - \mu)$', ...
            } ;
        plotOpts.makeCbar = [0,0,0,0,0,0] ;
        [~, ~, ~, plotOpts.xyzlims] = QS.getXYZLims() ;
        plotOpts.view = [0,0] ;
        plotOpts.axisOff = true ;
        plotOpts.masterCbar = true ;
        % 3d
        [axs, cbs, meshHandles] = nFieldsOnSurface(mesh, ...
            {real(mu_pred_vtx), ...
            real(squeeze(mu_all(tidx,:))), ...
            real(diff_mu), ...
            imag(mu_pred_vtx), ...
            imag(squeeze(mu_all(tidx,:))), ...
            imag(diff_mu)}, ...
            plotOpts) ;
        sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
            ' ', QS.timeUnits], 'Interpreter', 'latex')
        imfn = fullfile(predDir3d, [sprintf(QS.fileBase.name, tp) '.png']) ;
        disp(['Saving image: ' imfn])
        saveas(gcf, imfn)
        close all
        
        % 2d
        plotOpts.view = [0,90] ;
        plotOpts.makeCbar = [0,0,0,0,0,0] ;
        plotOpts.masterCbar = true ;
        plotOpts.clim = [-0.5, 0.5] ;
        plotOpts = rmfield(plotOpts, 'xyzlims') ;
        nFieldsOnSurface(refMesh2plot, ...
            {real(mu_pred_vtx), ...
            real(squeeze(mu_all(tidx,:))), ...
            real(diff_mu), ...
            imag(mu_pred_vtx), ...
            imag(squeeze(mu_all(tidx,:))), ...
            imag(diff_mu)}, ...
            plotOpts) ;
        sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
            ' ', QS.timeUnits], 'Interpreter', 'latex')
        imfn = fullfile(predDir2d, [sprintf(QS.fileBase.name, tp) '.png']) ;
        disp(['Saving image: ' imfn])
        saveas(gcf, imfn)
        close all

    end

    % Compute diff from each kymo 
    diff_apM = real(mu_apM) - indent_apM * 0.5 ;
    recompute = false ;
end

nU = QS.nU ;
nV = QS.nV ;
timePoints = QS.xp.fileMeta.timePoints ;


%% Plot kymograph of difference
if strcmpi(coordSys, 'ricci')
    fn = fullfile(sprintf(QS.dir.pathlines.quasiconformal, ...
        t0Pathlines), 'prediction_kymograph.png') ;
else
    fn = fullfile(sprintf(QS.dir.pathlines_uvprime.quasiconformal, ...
        t0Pathlines), 'prediction_kymograph.png') ;
end
if ~exist(fn, 'file')
    close all
    fields = {real(mu_apM), indent_apM * 0.5, real(mu_apM) - indent_apM * 0.5} ;
    labels = {'$\Re \mu_{\textrm{pred}}$', ...
            '$\Re\mu$', ...
            '$\Re(\mu_{\textrm{pred}} - \mu)$', ...
            '$\Im\mu_{\textrm{pred}}$', ...
            '$\Im\mu$', ...
            '$\Im(\mu_{\textrm{pred}} - \mu)$'} ; 
    for nn = 1:3
        subplot(1, 3, nn)
        if strcmpi(QS.timeUnits, 'min')
            imagesc(1:nU, (timePoints - QS.t0) / 60, fields{nn})
        else
            imagesc(1:nU, timePoints - QS.t0, fields{nn})
        end
        title(labels{nn}, 'interpreter', 'latex')
        axis tight
        caxis([-0.5, 0.5])
        xticks([])
        xlabel('ap position, $u$', 'interpreter', 'latex')
        if nn == 1
            ylabel('time [hr]', 'interpreter', 'latex')
        end
        colormap blueblackred
        cb = colorbar('location', 'southoutside') ;
    end
    disp(['saving fig: ' fn])
    saveas(gcf, fn) 
end

%% Obtain location of folds / features
featureOpts = struct() ;
featureOpts.field2 = 'radius' ;
if strcmpi(coordSys, 'ricci')
    featureIDs = QS.getPullbackPathlineFeatureIDs('vertices', featureOpts) ;
else
    featureIDs = QS.getUVPrimePathlineFeatureIDs('vertices', featureOpts) ;
end
QS.getFeatures()
fold_onset = QS.features.fold_onset ;
widths = [0.01, 0.02, 0.05] ;

%% Use features to estimate fold locations and plot constriction
close all
for width = widths
    close all

    if strcmpi(coordSys, 'ricci')
        fn = fullfile(sprintf(QS.dir.pathlines.quasiconformal, ...
            t0Pathlines), ['linear_approx_folds_width' num2str(round(width * nU)) '.png']) ;
    else
        fn = fullfile(sprintf(QS.dir.pathlines_uvprime.quasiconformal, ...
            t0Pathlines), ['linear_approx_folds_width' num2str(round(width * nU)) '.png']) ;    
    end

    nfolds = length(fold_onset) ;

    % Adjust time to hours if in minutes
    if strcmpi(QS.timeUnits, 'min')
        tt = timePoints / 60 ;
    else
        tt = timePoints ;
    end
    for qq = 1:length(fold_onset) 
        % tt = timePoints(fold_onset(qq):end) ;
        f1 = (featureIDs(qq) - round(width*nU)):(featureIDs(qq) + round(width*nU)) ;

        ax = subplot(nfolds, 1, qq) ;
        h1 = plot(tt, real(mean(mu_apM(:, f1, :), 2)), '-', 'color', QS.plotting.colors(1, :)) ;
        hold on;
        h2 = plot(tt, mean(indent_apM(:, f1, :), 2), '--', 'color', QS.plotting.colors(2, :)) ;
        ylabel(['fold ' num2str(qq)])
        if qq == length(fold_onset)
            if strcmpi(QS.timeUnits, 'min')
                xlabel('time [hr]', 'interpreter', 'latex')
            else
                xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
            end
        end
        
        % TODO: TAKE NMODES OF THE VARIATION, PLOT SCATTERPLOT OF FOLDS        
        
        if qq == 1
            legend({'$\mu$', '$\rho$'}, 'location', 'eastOutside', 'interpreter', 'latex')
            axpos = get(gca, 'position') ;
        else
            newpos = get(ax, 'position') ;
            set(gca, 'position', [newpos(1) newpos(2) axpos(3) axpos(4)]) ;
        end

    end
    saveas(gcf, fn)
end
