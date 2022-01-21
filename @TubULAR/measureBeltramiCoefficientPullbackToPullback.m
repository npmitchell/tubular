function measureBeltramiCoefficientPullbackToPullback(QS, options)
%measureBeltramiCoefficientPullbackToPullback(QS, options)
% 
% Measure the Beltrami using the Ricci flow pullback of a (advected) mesh 
%   at some timepoint t1 relative to the pullback of the same 
%   (non-advected) mesh at t0. Since Ricci flow is a conformal map, these
%   Beltramis should equal those of the 3D embedding to 2D pullback 
%   comparison performed in QS.measureBeltramiCoefficient().
%
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   t0Pathlines : numeric (default = QS.t0) 
%       timepoint used to define Lagrangian/material frame for mu
%       measurements to be made with respect to
%   coordSys : str specifier ('ricci', 'uvprime' or 'sp')
%       coordinate system in which pathlines are computed
%   maxIter : int
%       only used if coordSys == 'ricci', #ricci iterations
%   overwrite : bool, overwrite previous results on disk
%   climit : limit for caxis of scalar fields Re(mu) and Im(mu)
%
% Returns
% -------
%
% See also
% --------
% QS.measureBeltramiCoefficient()
%
% NPMitchell 2021

%% Default options
save_ims = true ;
coordSys = 'ricci' ;    % or 'uvprime'
zwidth = 3 ;            % spectral filter averaging in x (tripulse)
nmodes = 5 ;            % spectral filter #modes to keep in circumferential dimension
maxIter = 200 ;         % only used if coordSys == 'ricci', #ricci iterations
if isempty(QS.t0)
    t0 = QS.t0set() ;
else
    t0 = QS.t0 ;
end
t0Pathlines = t0 ;
overwrite = false ;
if nargin < 2 
    options = struct() ;
end
climit = 1 ;
[~, ~, ~, xyzlim ] = QS.getXYZLims() ; 

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'save_ims')
    save_ims = options.save_ims ;
end
if isfield(options, 't0Pathlines')
    t0Pathlines = options.t0Pathlines ;
else
    disp('Using default t0 for pathlines')
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end

%% Make sure pathlines are on disk & load them
if contains(coordSys, 'ricci')
    vXYfn = sprintf(QS.fileName.pathlines.vXY, t0Pathlines) ;
    v3dfn = sprintf(QS.fileName.pathlines.v3d, t0Pathlines) ;
    refMeshFn = sprintf(QS.fileName.pathlines.refMesh, t0Pathlines) ;
    if ~exist(vXYfn, 'file') || ~exist(v3dfn, 'file') || ~exist(refMeshFn, 'file')
        QS.measurePullbackPathlines(options)
    end
    fnBase = QS.fullFileBase.pathlines.ricciQuasiconformal ;
elseif contains(coordSys, 'uvprime')
    vXYfn = sprintf(QS.fileName.pathlines_uvprime.vXY, t0Pathlines) ;
    v3dfn = sprintf(QS.fileName.pathlines_uvprime.v3d, t0Pathlines) ;
    refMeshFn = sprintf(QS.fileName.pathlines_uvprime.refMesh, t0Pathlines) ;
    if ~exist(vXYfn, 'file') || ~exist(v3dfn, 'file') || ~exist(refMeshFn, 'file')
        QS.measureUVPrimePathlines(options)
    end
    fnBase = QS.fullFileBase.pathlines_uvprime.ricci.mu ;
else
    vXYfn = sprintf(QS.fileName.pathlines.vXY, t0Pathlines) ;
    v3dfn = sprintf(QS.fileName.pathlines.v3d, t0Pathlines) ;
    refMeshFn = sprintf(QS.fileName.pathlines.refMesh, t0Pathlines) ;
end

% Load pathlines in conformal/(u',v')/other coordinates and their embedding
tmp = load(vXYfn) ;
vP2d = tmp.vertexPathlines ;
tmp = load(v3dfn) ;
vP3d = tmp.v3dPathlines ;
load(refMeshFn, 'refMesh') ;
if contains(coordSys, 'ricci')
    u_material = refMesh.u_ricci ;
    mu0 = 0 + 1j * 0 ;
else
    u_material = [refMesh.u(:, 1), refMesh.u(:, 2)] ;
    mu0 = mean(real(bc_metric(refMesh.f, u_material, refMesh.v, 3))) ;
end
try
    assert(abs(mu0) < 1e-4)
catch
    disp('refMesh is not particularly conformal...')
end

%% 
nTimePoints = length(QS.xp.fileMeta.timePoints) ;

if contains(coordSys, 'ricci')
    imDir = fullfile(sprintf(QS.dir.pathlines.ricci.quasiconformal, t0Pathlines), 'images') ;
elseif contains(coordSys, 'uvprime')
    imDir = fullfile(sprintf(QS.dir.pathlines_uvprime.ricci.quasiconformal, t0Pathlines), 'images') ;
end
if ~exist(fullfile(imDir, '2d_pb2pb'), 'dir')
    mkdir(fullfile(imDir, '2d_pb2pb'))
    mkdir(fullfile(imDir, '3d_pb2pb'))
end

todo1 = [111, 1:10:nTimePoints] ;
tidx2do = [todo1, setdiff(1:nTimePoints, todo1)] ;

first = true ;
for tidx = tidx2do
    disp(['tidx = ' num2str(tidx)])

    %% Set current time
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    imFn2d_material = sprintf(fullfile(imDir, '2d_pb2pb', 'mu2d_material_%06d.png'), tp);
    imFn3d_material = sprintf(fullfile(imDir, '3d_pb2pb', 'mu3d_material_%06d.png'), tp);

    QS.setTime(tp)

    % Define data filename for Beltramis
    fn = sprintf(fnBase, t0Pathlines, maxIter, tp) ;
    
    if (~exist(fn, 'file') || overwrite) && false
        try
            %% Measure Beltrami Coefficient
            v3d = zeros(size(vP3d.vX, 2) * size(vP3d.vX, 3), 3) ;
            v3d(:, 1) = reshape(squeeze(vP3d.vXrs(tidx, :, :)), [], 1) ;
            v3d(:, 2) = reshape(squeeze(vP3d.vYrs(tidx, :, :)), [], 1) ;
            v3d(:, 3) = reshape(squeeze(vP3d.vZrs(tidx, :, :)), [], 1) ;
            
            % Generate Ricci flow result
            opts = struct() ;
            tpMesh = struct() ;
            tpMesh.nU = refMesh.nU ;
            tpMesh.nV = refMesh.nV ;
            tpMesh.pathPairs = refMesh.pathPairs ;
            tpMesh.v = v3d ;
            tpMesh.u = refMesh.u ;
            tpMesh.f = refMesh.f ;
            opts.maxIter = maxIter ;
            opts.cutMesh = tpMesh ;
            opts.save_ims = false ;
            opts.t0Pathlines = t0Pathlines ;
            ricciMesh = QS.generatePathlineRicciMeshTimePoint(tp, opts) ;

            %% Material Beltrami coefficient
            u_deformed = [ricciMesh.rectangle.u, 0*ricciMesh.rectangle.u(:, 1)] ;
            mu_material = bc_metric(refMesh.f, u_material, u_deformed, 3) ;

            %% Mode filter mu_material 
            refMesh2glue = refMesh ;
            refMesh2glue.v = v3d ;
            glueMesh = glueCylinderCutMeshSeam(refMesh2glue) ;
            nU = refMesh.nU ;
            nV = refMesh.nV ;

            [V2F, F2V] = meshAveragingOperators(glueMesh.f, glueMesh.v) ;
            mu_material_vtx = F2V * squeeze(mu_material(tidx, :))' ;
            mu_material_vtx(nU*(nV-1)+1:nU*nV) = mu_material_vtx(1:nU) ;
            mu_material_vtx = reshape(mu_material_vtx, [refMesh.nU, refMesh.nV]) ;

            %% Check filtered image
            if first 
                fnfilter = fullfile(imDir, sprintf('filt_test.png')) ;
                if ~exist(fnfilter, 'file')
                    nmodes2do = 1:5 ;
                    clf
                    nrows = nmodes2do(end) + 1 ;
                    for nmodes_ii = nmodes2do
                        options.widthX = 3 ;
                        options.nmodesY = nmodes_ii ;
                        muMVfilt_re = modeFilterQuasi1D(real(mu_material_vtx), options) ;

                        subplot(nrows, 3, 1 + 3*(nmodes_ii-1))
                        imagesc(1:nU, 1:nV, real(mu_material_vtx)') ;
                        axis equal; axis tight; axis off
                        caxis([-1,1])
                        colormap blueblackred
                        if nmodes_ii == 1
                            title('raw $\mu$', 'interpreter', 'latex')
                        elseif nmodes_ii == nmodes2do(end)
                            pos = get(gca, 'position') ;
                            cb = colorbar('location', 'southOutside');
                            set(gca, 'position', pos)
                        end
                        subplot(nrows, 3, 2 + 3*(nmodes_ii-1)) 
                        imagesc(1:nU, 1:nV, muMVfilt_re')
                        caxis([-1,1])
                        colormap blueblackred
                        axis equal; axis tight; axis off
                        if nmodes_ii == 1
                            title('filtered $\mu$', 'interpreter', 'latex')
                        elseif nmodes_ii == nmodes2do(end)
                            pos = get(gca, 'position') ;
                            cb = colorbar('location', 'southOutside');
                            set(gca, 'position', pos)
                        end
                        subplot(nrows, 3, 3 + 3*(nmodes_ii-1)) 
                        imagesc(1:nU, 1:nV, muMVfilt_re'-real(mu_material_vtx)') ;
                        caxis([-1,1])
                        colormap blueblackred
                        axis equal; axis tight; axis off
                        if nmodes_ii == 1
                            title('filtered - raw', 'interpreter', 'latex')
                        elseif nmodes_ii == nmodes2do(end)
                            pos = get(gca, 'position') ;
                            cb = colorbar('location', 'southOutside');
                            set(gca, 'position', pos)
                            % set(cb, 'position', get(cb, 'position'));
                        end
                    end
                    sgtitle(['Lowest ', num2str(options.nmodesY), ' modes, ', ...
                        '$\sigma= $' num2str(options.widthX)], 'interpreter', 'latex')  
                    disp(['saving ' fnfilter])
                    saveas(gcf, fnfilter)
                end
            end

            %% Do filter on data
            filterOptions.widthX = zwidth ;
            filterOptions.nmodesY = nmodes ;
            filterOptions.preview = false ;
            muMVfilt_re = modeFilterQuasi1D(real(mu_material_vtx), filterOptions) ;
            muMVfilt_im = modeFilterQuasi1D(imag(mu_material_vtx), filterOptions) ;
            mu_material_filtered = muMVfilt_re(:) + 1j * muMVfilt_im(:) ;

            if save_ims

                close all
                labels = {'$\Re \mu$', '$\Im \mu$'} ;

                %% Plot mu_material in 3d
                options.labels = labels ;
                options.clim = climit ;
                [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                    twoScalarFieldsOnSurface({refMesh.f, v3d}, ...
                    real(mu_material_filtered(tidx, :)), ...
                    imag(mu_material_filtered(tidx, :)), options) ;
                sgtitle(['$\mu($embedding, material frame$)$, $t = $', ...
                    sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                    'interpreter', 'latex') 
                set(gcf,'CurrentAxes', ax1)
                view(0, 0)
                xlim(xyzlim(1, :))
                ylim(xyzlim(2, :))
                zlim(xyzlim(3, :))
                axis off
                set(gcf,'CurrentAxes', ax2)
                view(0, 0)
                xlim(xyzlim(1, :))
                ylim(xyzlim(2, :))
                zlim(xyzlim(3, :))
                axis off
                disp(['saving ' imFn3d_material])
                saveas(gcf, imFn3d_material)
                close all

                %% Save image of mu_material in 2d 
                set(gcf, 'visible', 'off')
                [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                    twoScalarFieldsOnSurface({refMesh.f, ...
                    [refMesh.u(:, 1) / max(refMesh.u(:, 1)), ...
                    refMesh.u(:, 2), 0 * refMesh.u(:, 2)]}, ...
                    real(mu_material_filtered(tidx, :)), ...
                    imag(mu_material_filtered(tidx, :)), options) ;
                sgtitle(['$\mu($embedding, material frame$)$, $t = $', ...
                    sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                    'interpreter', 'latex') ;
                set(gcf,'CurrentAxes', ax1)
                view(2)
                set(gcf,'CurrentAxes', ax2)
                view(2)
                disp(['saving ' imFn2d_material])
                saveas(gcf, imFn2d_material)
                close all
            end

            % No longer the first pass
            first = false ;
            
            % Save data as output
            disp(['saving ' fn])
            save(fn, 'mu_material', 'mu_material_filtered', 'filterOptions') 
        catch
            disp('could not generate Ricci mesh -- self intersections or bad face quality?')
        end
    end
    
    % Check that this is the same as computing from embedding only instead
    % of comparing PB(t0) to PB(t1)
end
disp('done with measureBeltramiCoefficient')


%% Compare to 3d->2d beltrami computation
if contains(coordSys, 'ricci')
    ricci2d3dFn = sprintf(QS.fileName.pathlines.quasiconformal, t0Pathlines) ;
elseif contains(coordSys, 'uvprime')
    ricci2d3dFn = sprintf(QS.fileName.pathlines_uvprime.quasiconformal, t0Pathlines) ;
end
r2d3d = load(ricci2d3dFn, 'mu_material', 'mu_material_filtered', 'filterOptions') ;
nU = refMesh.nU ; 
nV = refMesh.nV ;
m2d = refMesh ;
m2d.v = refMesh.u ;
m2d.v(:, 1) = m2d.v(:, 1) / max(m2d.v(:, 1)) ;
m2d.v(:, 2) = m2d.v(:, 2) / max(m2d.v(:, 2)) ;
 
for tidx = tidx2do
    disp(['tidx = ' num2str(tidx)])

    %% Set current time
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    imFn2d_material = sprintf(fullfile(imDir, '2d_pb2pb', 'mu2d_material_%06d.png'), tp);
    imFn3d_material = sprintf(fullfile(imDir, '3d_pb2pb', 'mu3d_material_%06d.png'), tp);

    QS.setTime(tp)

    % Define data filename for Beltramis
    fn = sprintf(fnBase, t0Pathlines, maxIter, tp) ;
    
    if exist(fn, 'file') 
        % Compare to 2d3d result
        load(fn, 'mu_material', 'mu_material_filtered', 'filterOptions')
        muf = mu_material ;
        muv = mu_material_filtered ;
        f3d = r2d3d.mu_material(tidx, :) ;
        v3d = r2d3d.mu_material_filtered(tidx, :) ;
        
        clf
        opts.visible = 'on' ;
        opts.view = [0, 90] ;
        opts.clims = {[], [], [], [0,1]} ;
        opts.axisOff = 'true' ;
        opts.colormap = 'gist_earth' ;
        opts.style = 'positive' ;
        delta = abs(f3d(:)-muf(:)) ;
        dfrac = abs(delta) ./ abs(f3d(:)) ;
        opts.labels = {['2D, max=' num2str(max(abs(muf(:))))], ...
            ['3D, max=' num2str(max(abs(f3d(:))))], ...
            ['$\Delta$, max=' num2str(max(delta))], ...
            ['$\Delta$/3D, max=' num2str(max(dfrac))]} ;
        nFieldsOnSurface(m2d, ...
            {abs(muf(:)), abs(f3d(:)), delta, dfrac}, opts)
        [tpTrue, timestr ] = trueTime(QS, tp, false) ;
        sgtitle(['t = ' timestr ', unfiltered'])
        saveas(gcf, fullfile(...
            sprintf(QS.dir.pathlines.ricci.quasiconformal, t0Pathlines), ...
            sprintf('%03d_diff_raw.png', tp)))
        
        clf
        opts.visible = 'on' ;
        opts.view = [0, 90] ;
        opts.clims = {[], [], [], [0,1]} ;
        opts.axisOff = 'true' ;
        opts.colormap = 'gist_earth' ;
        opts.style = 'positive' ;
        delta = abs(v3d(:)-muv(:)) ;
        dfrac = abs(delta) ./ abs(v3d(:)) ;
        opts.labels = {['2D, max=' num2str(max(abs(muv(:))))], ...
            ['3D, max=' num2str(max(abs(v3d(:))))], ...
            ['$\Delta$, max=' num2str(max(delta))], ...
            ['$\Delta$/3D, max=' num2str(max(dfrac))]} ;
        nFieldsOnSurface(m2d, ...
            {abs(muv(:)), abs(v3d(:)), delta, dfrac}, opts)
        [tpTrue, timestr ] = trueTime(QS, tp, false) ;
        sgtitle(['t = ' timestr ', filtered'])
        saveas(gcf, fullfile(...
            sprintf(QS.dir.pathlines.ricci.quasiconformal, t0Pathlines), ...
            sprintf('%03d_diff_filtered.png', tp)))
    end
    
end

