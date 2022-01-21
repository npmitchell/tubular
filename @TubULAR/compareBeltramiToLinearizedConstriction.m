function compareBeltramiToLinearizedConstriction(QS, options)
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
coordSys = 'ricci' ;
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
    outdir = sprintf(QS.dir.pathlines.quasiconformal, t0Pathlines) ; 
elseif strcmpi(coordSys, 'uvprime')
    v3d = load(sprintf(QS.fileName.pathlines_uvprime.v3d, t0Pathlines));
    load(sprintf(QS.fileName.pathlines_uvprime.refMesh, t0Pathlines), 'refMesh') ;
    outdir = sprintf(QS.dir.pathlines_uvprime.quasiconformal, t0Pathlines) ;
end
v3x = v3d.v3dPathlines.vXrs ;
v3y = v3d.v3dPathlines.vYrs ;
v3z = v3d.v3dPathlines.vZrs ;

% Load refMesh
u2d = [refMesh.u(:, 1), refMesh.u(:, 2), 0*refMesh.u(:, 2)] ;
mesh2 = {refMesh.f, u2d} ;

% Compare mu to indentation rho/2
outdir = fullfile(outdir, 'linear_approx') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Get xyzlim
[~,~,~,xyzlim] = QS.getXYZLims() ;

indent_computed = false ;
loaded_mu = false ;

recompute = true ;
if ~overwrite
    try
        % Load indentation kymograph
        if strcmpi(coordSys, 'ricci')
            indentAPfn = sprintf(...
                QS.fileName.pathlines.kymographs.indentation, t0Pathlines) ;
            muAPfn = sprintf(...
                QS.fileName.pathlines.kymographs.mu, t0Pathlines) ;
        elseif strcmpi(coordSys, 'uvprime')            
            indentAPfn = sprintf(...
                QS.fileName.pathlines_uvprime.kymographs.indentation, t0Pathlines) ;
            muAPfn = sprintf(...
                QS.fileName.pathlines_uvprime.kymographs.mu, t0Pathlines) ;
        end
        disp(['Attempting to load indent kymos from file: ' indentAPfn])
        load(indentAPfn, 'indent_apM', 'indent_dM', 'indent_vM', ...
            'indent_lM', 'indent_rM') ;
        disp('loaded')

        % Save mu Beltrami coefficient kymograph
        disp(['Loading mu kymos from file: ' muAPfn])
        load(muAPfn, 'mu_apM', 'mu_dM', 'mu_vM', 'mu_lM', 'mu_rM') ;

        % Compute diff from each kymo 
        diff_apM = real(mu_apM) - indent_apM * 0.5 ;
        recompute = false ;
    catch 
        disp('Could not load from disk, computing...')
    end
end

if recompute
    timePoints = QS.xp.fileMeta.timePoints ;
    mu_apM = zeros(length(timePoints), QS.nU) ;
    mu_dM = zeros(length(timePoints), QS.nU) ;
    mu_vM = zeros(length(timePoints), QS.nU) ;
    mu_lM = zeros(length(timePoints), QS.nU) ;
    mu_rM = zeros(length(timePoints), QS.nU) ;
    indent_apM = zeros(length(timePoints), QS.nU) ;
    indent_dM = zeros(length(timePoints), QS.nU) ;
    indent_vM = zeros(length(timePoints), QS.nU) ;
    indent_lM = zeros(length(timePoints), QS.nU) ;
    indent_rM = zeros(length(timePoints), QS.nU) ;
    diff_apM = zeros(length(timePoints), QS.nU) ;
    diff_dM = zeros(length(timePoints), QS.nU) ;
    diff_vM = zeros(length(timePoints), QS.nU) ;
    diff_lM = zeros(length(timePoints), QS.nU) ;
    diff_rM = zeros(length(timePoints), QS.nU) ;
    
    % Grab quarter indices for dv slices
    [dorsal, ventral, left, right] = QS.quarterIndicesDV(QS.nV) ;
    
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        disp(['t = ', num2str(tp)])

        %% GET DATA
        % Load quasiconformal results for this material frame (t0Pathlines)
        if ~loaded_mu

            if strcmpi(coordSys, 'ricci')
                fn = sprintf(QS.fileName.pathlines.quasiconformal, t0Pathlines) ;
            elseif strcmpi(coordSys, 'uvprime')
                fn = sprintf(QS.fileName.pathlines_uvprime.quasiconformal, t0Pathlines) ;
            end
            load(fn, 'mu_material_filtered')
            loaded_mu = true ;
        end
        muVtx = mu_material_filtered(tidx, :) ;
        
        % Load radial indentation in Lagrangian frame
        opts.coordSys = 'uvprime' ;
        opts.overwrite = false ;
        if ~indent_computed
            indentation = QS.measurePathlineIndentation(opts) ;
            indent_computed = true ;
        end
        indent = squeeze(indentation(tidx, :, :)) ;

        % Take difference between linear radial contraction description and
        % indentation
        xx = v3x(tidx, :, :) ;
        yy = v3y(tidx, :, :) ;
        zz = v3z(tidx, :, :) ; 
        v3d = [xx(:), yy(:), zz(:)] ;

        % refMesh2glue = refMesh ;
        % refMesh2glue.v = v3d ;
        % glueMesh = glueCylinderCutMeshSeam(refMesh2glue) ;
        % 
        % [V2F, F2V] = meshAveragingOperators(glueMesh.f, glueMesh.v) ;
        % % difference = real(mu(:)) - V2F * indent(:) * 0.5 ;
        % muVtx = F2V*mu(:) ;
        nU = refMesh.nU ;
        nV = refMesh.nV ;
        % muVtx(nU*(nV-1) + 1:nU*nV) = muVtx(1:nU) ;
        difference = real(muVtx(:)) - indent(:) * 0.5 ;


        %% MAKE FIGURE
        outfn = fullfile(outdir, sprintf('linear_quasiconformal_%06d.png', tp)) ;
        if ~exist(outfn, 'file') || overwrite
            mesh3 = {refMesh.f, v3d} ;
            meshes = {mesh3, mesh3, mesh3, mesh2, mesh2, mesh2} ;
            fields = {real(muVtx(:)), indent(:), difference, ...
                real(muVtx(:)), indent(:), difference} ;
            options = struct() ;
            options.makeCbar = [false, false, false, true, true, true] ;
            options.clim = climit ;
            close all
            [axs, cbs, meshHandles] = ...
                nFieldsOnSurface(meshes, fields, options) ;
            titles = {'$\mu$', '$\rho/2$', '$\mu-\rho/2$'} ;            
            for qq = 1:3
                set(gcf,'CurrentAxes', axs{qq})
                view(0,0)
                axis off
                xlim(xyzlim(1, :))
                ylim(xyzlim(2, :))
                zlim(xyzlim(3, :))
                title(titles{qq}, 'interpreter', 'latex')
            end
            for qq = 4:6
                set(gcf,'CurrentAxes', axs{qq})
                view(2)
                axis off
            end
            sgtitle(['quasiconformal constriction, $t=$', ...
                sprintf('%03d', (tp-t0)*QS.timeInterval), ' ', QS.timeUnits], ...
                'interpreter', 'latex')

            % Save outfn
            disp(['saving ' outfn])
            saveas(gcf, outfn)
        end

        %% ADD TO KYMOGRAPH
        muVtx = reshape(muVtx, [nU, nV]) ;
        mu_apM(tidx, :) = mean(muVtx(:, 1:nV-1), 2) ;
        mu_dM(tidx, :) = mean(muVtx(:, dorsal), 2) ;
        mu_vM(tidx, :) = mean(muVtx(:, ventral), 2) ;
        mu_lM(tidx, :) = mean(muVtx(:, left), 2) ;
        mu_rM(tidx, :) = mean(muVtx(:, right), 2) ;
        indent_apM(tidx, :) = mean(indent(:, 1:nV-1), 2) ;
        indent_dM(tidx, :) = mean(indent(:, dorsal), 2) ;
        indent_vM(tidx, :) = mean(indent(:, ventral), 2) ;
        indent_lM(tidx, :) = mean(indent(:, left), 2) ;
        indent_rM(tidx, :) = mean(indent(:, right), 2) ;

        % difference kymographs
        difference = reshape(difference, [nU, nV]) ;
        diff_apM(tidx, :) = mean(difference(:, 1:nV-1), 2) ;
        diff_dM(tidx, :) = mean(difference(:, dorsal), 2) ;
        diff_vM(tidx, :) = mean(difference(:, ventral), 2) ;
        diff_lM(tidx, :) = mean(difference(:, left), 2) ;
        diff_rM(tidx, :) = mean(difference(:, right), 2) ;
    end

    % Save indentation kymograph
    if strcmpi(coordSys, 'ricci')
        indentAPfn = sprintf(...
            QS.fileName.pathlines.kymographs.indentation, t0Pathlines) ;
        muAPfn = sprintf(...
            QS.fileName.pathlines.kymographs.mu, t0Pathlines) ;
    elseif strcmpi(coordSys, 'uvprime')
        indentAPfn = sprintf(...
            QS.fileName.pathlines_uvprime.kymographs.indentation, t0Pathlines) ;
        muAPfn = sprintf(...
            QS.fileName.pathlines_uvprime.kymographs.mu, t0Pathlines) ;
    end
        
    % Save indentation and mu Beltrami coefficient kymograph
    save(indentAPfn, 'indent_apM', 'indent_dM', 'indent_vM', ...
        'indent_lM', 'indent_rM') ;
    save(muAPfn, 'mu_apM', 'mu_dM', 'mu_vM', 'mu_lM', 'mu_rM') ;
end

% 
nU = QS.nU ;
nV = QS.nV ;
timePoints = QS.xp.fileMeta.timePoints ;

%% Plot kymograph of difference
if strcmpi(coordSys, 'ricci')
    fn = fullfile(sprintf(QS.dir.pathlines.quasiconformal, ...
        t0Pathlines), 'linear_approx_kymograph.png') ;
elseif strcmpi(coordSys, 'uvprime')
    fn = fullfile(sprintf(QS.dir.pathlines_uvprime.quasiconformal, ...
        t0Pathlines), 'linear_approx_kymograph.png') ;
end
if ~exist(fn, 'file')
    close all
    fields = {real(mu_apM), indent_apM * 0.5, real(mu_apM) - indent_apM * 0.5} ;
    labels = {'$\mu$', '$\rho/2$', '$\mu - \rho/2$'} ;
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

%% Plot without factor of 2 -- no good reason yet
if strcmpi(coordSys, 'ricci')
    fn = fullfile(sprintf(QS.dir.pathlines.quasiconformal, ...
        t0Pathlines), 'other_approx_kymograph.png') ;
elseif strcmpi(coordSys, 'uvprime')
    fn = fullfile(sprintf(QS.dir.pathlines_uvprime.quasiconformal, ...
        t0Pathlines), 'other_approx_kymograph.png') ;
end
if ~exist(fn, 'file')
    close all
    fields = {real(mu_apM), indent_apM, real(mu_apM) - indent_apM} ;
    labels = {'$\mu$', '$\rho$', '$\mu - \rho$'} ;
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
        colormap blueblackred
        xticks([])
        xlabel('ap position, $u$', 'interpreter', 'latex')
        if nn == 1
            ylabel('time [hr]', 'interpreter', 'latex')
        end
        cb = colorbar('location', 'southoutside') ;
    end
    disp(['saving fig: ' fn])
    saveas(gcf, fn) 
end

%% Obtain location of folds / features
featureOpts = struct() ;
featureOpts.field2 = 'radius' ;
if strcmpi(coordSys, 'ricci')
    featureIDs = QS.getPathlineFeatureIDs('vertices', featureOpts) ;
elseif strcmpi(coordSys, 'uvprime')
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
    elseif strcmpi(coordSys, 'uvprime')
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
