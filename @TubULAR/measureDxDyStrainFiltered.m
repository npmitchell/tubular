function measureDxDyStrainFiltered(QS, options) 
% measureDxDyStrainFiltered(QS, options)
%   Compute strain along bonds that lie along dx and dy in pullback space
%   from integrated pathlines deforming mesh vertices. Filter the results
%   heavily in space. 
%   Measurements are taken with respect to fixed Lagrangian frame. 
%   Plot results in 2d and/or 3d for each timepoint.   
%
% Parameters
% ----------
% QS : QuapSlap class instance
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

%% Parameter options
lambda_mesh = 0.002 ;
lambda = 0.01 ; 
debug = false ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = "Lagrangian" ;
LagrangianFrame = 'zetaphi' ;
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
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
bbr256 = blueblackred(256) ;
close all

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
%% Load pathlines to build metric from advected mesh faces
QS.loadPullbackPathlines(t0Pathline, 'vertexPathlines')
vP = QS.pathlines.vertices ;

% Output directory is inside pathline dir
outdir = sprintf(QS.dir.pathlines.strain, t0) ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

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

%% INTEGRATE STRAINRATE INTO STRAIN ON PATHLINES
% Compute or load all timepoints
load(sprintf(QS.fileName.pathlines.v3d, t0), 'v3dPathlines')
load(sprintf(QS.fileName.pathlines.vXY, t0), 'vertexPathlines')
vX3rs = v3dPathlines.vXrs ;
vY3rs = v3dPathlines.vYrs ;
vZ3rs = v3dPathlines.vZrs ;
t0 = v3dPathlines.t0 ;
tIdx0 = v3dPathlines.tIdx0 ;

if ~all(isfinite(vX3rs(:)))
    error('Some vertex 3D positions are infinite. Check interpolation')
end

% Define reference mesh
refMeshFn = fullfile(sprintf(QS.dir.pathlines.data, t0Pathline), ...
        sprintf('referenceMeshMaterialFrame_%04d.mat', t0Pathline)) ;
if exist(refMeshFn, 'file') || overwrite
    load(refMeshFn, 'refMesh')
else
    refMesh = struct() ; 
    vX = vertexPathlines.vX(tIdx0, :) ; 
    vY = vertexPathlines.vY(tIdx0, :) ;
    vXY = [vX(:), vY(:)] ;
    Lx = vertexPathlines.Lx ;
    Ly = vertexPathlines.Ly ;
    refMesh.f = defineFacesRectilinearGrid(vXY, QS.nU, QS.nV) ;
    refMesh.u = QS.XY2uv([Lx(tIdx0), Ly(tIdx0)], vXY, 1, 1) ;
    x0 = vX3rs(tIdx0, :) ;
    y0 = vY3rs(tIdx0, :) ;
    z0 = vZ3rs(tIdx0, :) ;
    refMesh.v = [ x0(:), y0(:), z0(:) ] ; 
    pathPairs = [ (1:nU)', (nV-1)*nU + (1:nU)' ] ;
    refMesh.pathPairs = pathPairs ;
    refMesh.nU = QS.nU ;
    refMesh.nV = QS.nV ;
    assert(numel(refMesh.u(:, 1)) == QS.nU * QS.nV)

    % Save reference Mesh as lagrangian frame
    save(refMeshFn, 'refMesh') ;
end

%% Prep directories
outdir = fullfile(sprintf(QS.dir.pathlines.data, ...
            t0Pathline), 'DxDyStrainFiltered') ;
egImDir2d = fullfile(outdir, 'strainDxDy_2d') ;    
egImDir3d = fullfile(outdir, 'strainDxDy_3d') ;    
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
if ~exist(egImDir2d, 'dir')
    mkdir(egImDir2d)
end
if ~exist(egImDir3d, 'dir')
    mkdir(egImDir3d)
end

%% Load refmesh
load(sprintf(QS.fileName.pathlines.refMesh, t0), 'refMesh')
glueMesh = glueCylinderCutMeshSeam(refMesh) ;
FF = glueMesh.f ;

%% Load pathline vertices for all time
tmp = load(sprintf(QS.fileName.pathlines.v3d, t0)) ;
vXrs = tmp.v3dPathlines.vXrs ;
vYrs = tmp.v3dPathlines.vYrs ;
vZrs = tmp.v3dPathlines.vZrs ;
clearvars tmp 

%% Compute for each timepoint
ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
tidx2do = 1:40:ntps ;
tidx2do = [tidx2do setdiff(1:ntps, tidx2do)] ;
for tidx = tidx2do
    % Identify current timepoint
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    disp(['t = ', num2str(tp)])

    % Compute the result
    dxdyFn = fullfile(outdir, sprintf(['dxdy_' QS.fileBase.strain], tp)) ;
    if ~exist(dxdyFn, 'file') || overwrite
        ffn = sprintf(QS.fullFileBase.pathlines.strain, t0, tp) ;
        try
            load(ffn, 'bondDxDy') ;
            dx_strain = bondDxDy.dx_strain ;
            dy_strain = bondDxDy.dy_strain ;
        catch
            disp('bondDxDy not on disk, first run QS.measurePathlineStrain()')
        end

        % Shift onto vertices -- note these are advected vertices
        QS.setTime(tp) ;
        % remove duplicate periodic seam from unwrapping of XYZ(tp)
        xx = squeeze(vXrs(tidx, :, :)) ;
        yy = squeeze(vYrs(tidx, :, :)) ;
        zz = squeeze(vZrs(tidx, :, :)) ;
        xg = xx(:, 1:end-1) ;
        yg = yy(:, 1:end-1) ;
        zg = zz(:, 1:end-1) ;
        [~, F2V] = meshAveragingOperators(FF, [xg(:), yg(:), zg(:)]) ;
        dx = reshape(F2V * dx_strain, [nU, nV-1]) ;
        dy = reshape(F2V * dy_strain, [nU, nV-1]) ;

        %% Heavily filter
        opts.nmodes = 3 ;
        opts.zwidth = 3 ;
        dxs = modeFilterQuasi1D(dx, opts) ;
        dys = modeFilterQuasi1D(dy, opts) ;

        %% Expand to cutMesh
        dxs(:, nV) = dx(:, 1) ;
        dys(:, nV) = dy(:, 1) ;

        %% Save result
        save(dxdyFn, 'dxs', 'dys') 

        %% Plot the result
        time_in_units = (tp - t0) * QS.timeInterval ;
        tstr = ['$t=$', sprintf('%03d', time_in_units), ' ', QS.timeUnits ];
        fn_dxdy2d = fullfile(egImDir2d, ...
            sprintf([QS.fileBase.spcutMeshSm '_DxDyFiltered.png'], tp));
        fn_dxdy3d = fullfile(egImDir3d, ...
            sprintf([QS.fileBase.spcutMeshSm '_DxDyFiltered.png'], tp));
        set(gcf, 'visible', 'off') ;
        if ~exist(fn_dxdy2d, 'file') || overwriteImages || overwrite
            % Panel 1
            subplot(1, 2, 1)
            trisurf(refMesh.f, ...
                refMesh.u(:, 1) / max(refMesh.u(:, 1)), ...
                refMesh.u(:, 2), 0 * refMesh.u(:, 2), ...
                dxs(:), 'edgecolor', 'none')
            daspect([1,1,1])
            cb = colorbar('location', 'southOutside') ;
            caxis([-1, 1])
            title('$\varepsilon_\zeta$', 'Interpreter', 'Latex')   
            colormap(bbr256)
            axis off; view(2)

            % Panel 2 
            subplot(1, 2, 2)
            trisurf(refMesh.f, ...
                refMesh.u(:, 1) / max(refMesh.u(:, 1)), ...
                refMesh.u(:, 2), 0 * refMesh.u(:, 2), ...
                dys(:), 'edgecolor', 'none')
            daspect([1,1,1])
            cb = colorbar('location', 'southOutside') ;
            caxis([-1, 1])
            title('$\varepsilon_\phi$', 'Interpreter', 'Latex')   
            colormap(bbr256)
            axis off; view(2)

            % Save the image
            sgtitle(['strain, ', tstr], 'Interpreter', 'latex') 
            saveas(gcf, fn_dxdy2d) ;
            clf
        end    
        close all

        %% Now plot dx dy in 3d
        set(gcf, 'visible', 'off') ;
        if ~exist(fn_dxdy3d, 'file') || overwrite
            % Panel 1 --> dx
            subplot(1, 2, 1)
            trisurf(refMesh.f, ...
                xx(:), yy(:), zz(:), ...
                dxs(:), 'edgecolor', 'none')
            daspect([1,1,1])
            cb = colorbar('location', 'southOutside') ;
            caxis([-1, 1])
            title('$\varepsilon_\zeta$', 'Interpreter', 'Latex')   
            colormap(bbr256)

            xyzlim = QS.plotting.xyzlim_um_buff ;
            tposZ = xyzlim(3, 2) + 0.5 * (xyzlim(3, 2) - xyzlim(3, 1)) ;
            tposY = xyzlim(2, 2) + 0.5 * (xyzlim(2, 2) - xyzlim(2, 1)) ;
            tposYZ = max(tposY, tposZ) ;
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            view(0, 0)
            tpos = get(get(gca, 'title'), 'position') ; 
            set(get(gca,'title'), 'Position', [tpos(1) tpos(2) tposYZ])
            colormap(bbr256)
            axpos = get(gca, 'position') ;
            set(cb, 'position', [0.130    0.1735    0.3347    0.0508])
            set(gca, 'position', axpos)
            axis off

            % Panel 2 --> dy
            subplot(1, 2, 2)
            trisurf(refMesh.f, ...
                xx(:), yy(:), zz(:), ...
                dys(:), 'edgecolor', 'none')
            daspect([1,1,1])
            cb = colorbar('location', 'southOutside') ;
            caxis([-1, 1])
            title('$\varepsilon_\phi$', 'Interpreter', 'Latex')   
            colormap(bbr256)
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            view(0, 0)
            tpos = get(get(gca, 'title'), 'position') ; 
            set(get(gca,'title'), 'Position', [tpos(1) tpos(2) tposYZ])
            colormap(bbr256)
            axpos = get(gca, 'position') ;
            set(cb, 'position', [0.5703    0.1735     0.3347    0.0508])
            set(gca, 'position', axpos)
            axis off

            % Save the image
            sgtitle(['orthogonal strain components, ', tstr], 'Interpreter', 'latex') 
            saveas(gcf, fn_dxdy3d) ;
            clf
        end    
        close all
    end
end
disp('done with integrated pathline strain calculations')
