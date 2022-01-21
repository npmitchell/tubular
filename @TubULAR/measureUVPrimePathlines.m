function measureUVPrimePathlines(QS, options)
%measurePullbackPathlines(QS, options)
%   Measure pathlines of optical flow in pullback space. 
%   These pathlines can then be used to query velocities or other
%   properties that are parameterized by pullback coordinates. 
%   For example, we average velocites along Lagrangian pathlines.
%   For another example, to build a Lagrangian-frame measure of divergence 
%   of tangential velocity, div(v_t), we can query div(v_t) at the coords
%   of the PullbackPathlines (via interpolation of div(v_t) defined on
%   pullback mesh vertices). For now assumes all pullbacks are the same
%   size -- for ex, 2000 x 2000 pixels.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   smoothing_sigma : float 
%       Gaussian kernel standard deviation if >0 for smoothing in-plane
%       displacements from piv readout
%
% Saves to disk
% -------------
% sprintf(QS.fileName.pathlines_uvprime.XY, t0) ;
% sprintf(QS.fileName.pathlines_uvprime.vXY, t0) ;
% sprintf(QS.fileName.pathlines_uvprime.fXY, t0) ;
% sprintf(QS.fileName.pathlines_uvprime.XYZ, t0) ;
%
% See also
% --------
% timeAverageVelocities(QS, options)
% pullbackPathlines(QS, options)
%
% NPMitchell 2020

%% Default options
overwrite = false ;
overwriteImages = false ;
preview = false ;
debug = false ;
timePoints = QS.xp.fileMeta.timePoints ;
samplingResolution = '1x' ;
nY2plot = 30 ;
scatterSz = 2 ;
movieScatterSz = 2 ;
t0 = QS.t0set() ;
smoothing_sigma = 2 ;  % Smoothing parameter for piv

if nargin < 2
    options = struct() ;
end

%% Unpack options
% Default values for options are to use sphi smoothed extended coords
% as PIV reference coord sys
if isfield(options, 't0Pathlines')
    t0 = options.t0Pathlines ;
    disp(['Setting t0 for Pathlines to be: ' num2str(t0)])
elseif isfield(options, 't0')
    t0 = options.t0 ;
    disp(['Setting t0 for Pathlines to be: ' num2str(t0)])
else
    disp('Using default t0 for pathlines')
end
if isfield(options, 'smoothing_sigma')
    smoothing_sigma = options.smoothing_sigma ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'nY2plot')
    nY2plot = options.nY2plot ;
end
if isfield(options, 'doubleCovered')
    doubleCovered = options.doubleCovered;
else
    doubleCovered = true ;
end
if doubleCovered
    Yoffset = 0.25 ;
else
    Yoffset = 0.0 ;
end

%% Unpack QS
% [rot, ~] = QS.getRotTrans() ;
% resolution = QS.APDV.resolution ; 
[~, ~, ~, xyzlim] = QS.getXYZLims() ;
axis_order = QS.data.axisOrder ;
blue = QS.plotting.colors(1, :) ;
red = QS.plotting.colors(2, :) ;
green = QS.plotting.colors(4, :) ;

% Define t0 at which pathlines form grid
tIdx0 = QS.xp.tIdx(t0) ;

%% Create directory for pathlines
pathlineDir = QS.dir.pathlines_uvprime ;
pdir = sprintf(pathlineDir.data, t0) ;
XYdir = sprintf(pathlineDir.XY, t0) ;
XYZdir = sprintf(pathlineDir.XYZ, t0) ;
vXYdir = sprintf(pathlineDir.vXY, t0) ;
fXYdir = sprintf(pathlineDir.fXY, t0) ;
v3ddir = sprintf(pathlineDir.v3d, t0) ;
f3ddir = sprintf(pathlineDir.f3d, t0) ;
dirs2make = {pdir, XYdir, XYZdir, vXYdir, fXYdir, v3ddir, f3ddir} ;
for qq = 1:length(dirs2make)
    dir2make = dirs2make{qq} ;
    if ~exist(dir2make, 'dir')
        disp(['Making dir: ' dir2make])
        mkdir(dir2make)
    end
end

%% Perform/Load pathline calculations

% Check if the time smoothed velocities exist already
% 2d velocities (pulled back), scaled by dilation of metric are v2dum 
% 2D velocities (pulled back) are v2d
% normal velocities on fieldfaces are vn
% 3d velocities on fieldfaces are v3d
% vertex-based velocities are vv
% face-based velocities are vf

QS.clearTime() ;

%% CREATE REFERENCE MESH 
refMeshFn = sprintf(QS.fileName.pathlines_uvprime.refMesh, t0) ;
if ~exist(refMeshFn, 'file') || overwrite 
    refMesh = load(sprintf(QS.fullFileBase.uvpcutMesh, t0), 'uvpcutMesh') ;
    refMesh = refMesh.uvpcutMesh.resampled ;
    umax0 = max(refMesh.u(:, 1)) ;
    vmax0 = max(refMesh.u(:, 2)) ;
    % Get size Lx and Ly for XY space
    im0 = imread(sprintf(QS.fullFileBase.im_uvprime_e, t0)) ;
    % for now assume all images are the same size
    disp('for now assuming all images are the same size')
    Lx = size(im0, 1) * ones(length(timePoints), 1) ;
    Ly = size(im0, 2) * ones(length(timePoints), 1) ;
    m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], refMesh.u, doubleCovered, umax0, vmax0) ;
    refMesh.XY = m0XY ;
    refMesh.vrs = QS.xyz2APDV(refMesh.v) ;
    try
        refMesh.mu ;
    catch
        refMesh.mu = bc_metric(refMesh.f, refMesh.u, refMesh.v, 3) ;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute dzeta -- material frame distance between DV hoops
    % (dzeta along the surface from each hoop to the next)
    % The distance from one hoop to another is the
    % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
    refMesh.dzeta = reshape(vecnorm(...
        diff(reshape(refMesh.vrs, [refMesh.nU, refMesh.nV, 3])), 2, 3), ...
        [refMesh.nU-1, refMesh.nV]) ;
    refMesh.dzeta_mean = nanmean(refMesh.dzeta, 2) ;
    
    save(refMeshFn, 'refMesh')
else
    load(refMeshFn, 'refMesh')
    try
        refMesh.mu ;
    catch
        refMesh.mu = bc_metric(refMesh.f, refMesh.u, refMesh.v, 3) ;
    end
end
mu0 = mean(real(refMesh.mu)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIV evaluation coordinates in XY pixel space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineXY = sprintf(QS.fileName.pathlines_uvprime.XY, t0) ;
if ~exist(plineXY, 'file') || overwrite
    disp('Computing piv Pathlines for XY(t0)')
    % Load 'initial' positions for pathlines to intersect at t=t0
    uvpPIVfn = fullfile(QS.dir.piv_uvprime, 'piv_results.mat') ;
    piv = load(uvpPIVfn) ;
    x0 = piv.x{QS.xp.tIdx(t0)} ;
    y0 = piv.y{QS.xp.tIdx(t0)} ;
    
    % Create pathlines emanating from (x0,y0) at t=t0
    % Additionally smooth the piv output by sigma
    if smoothing_sigma > 0
        disp(['Smoothing piv output with sigma=' num2str(smoothing_sigma)])
        for tidx = 1:length(QS.xp.fileMeta.timePoints)-1
            velx = piv.u_filtered{tidx} ;
            vely = piv.v_filtered{tidx} ;
            piv.u_filtered{tidx} = imgaussfilt(velx, smoothing_sigma) ;
            piv.v_filtered{tidx} = imgaussfilt(vely, smoothing_sigma) ;
        end
        disp('done smoothing')
    end
    
    options.piv = piv ;
    im0 = imread(sprintf(QS.fullFileBase.im_uvprime_e, t0)) ;
    % for now assume all images are the same size
    disp('for now assuming all images are the same size')
    Lx = size(im0, 1) * ones(length(timePoints), 1) ;
    Ly = size(im0, 2) * ones(length(timePoints), 1) ;
    options.Lx = Lx ;
    options.Ly = Ly ;
    [XX, YY] = QS.pullbackPathlines(x0, y0, t0, options) ;
    
    % Save pathlines of PIV evaluation points XY
    pivPathlines = struct() ;
    pivPathlines.XX = XX ;
    pivPathlines.YY = YY ;
    pivPathlines.Lx = Lx ;
    pivPathlines.Ly = Ly ;
    pivPathlines.t0 = t0 ;
    pivPathlines.tIdx0 =  tIdx0 ;
    save(plineXY, 'pivPathlines')
    
    computed_XY = true ;
else
    disp(['PIV XY pathlines already on disk: ' plineXY])
    computed_XY = false ;
end

%% Save image of result -- single image
plineFig = [plineXY(1:end-4) '.png'] ;
if ~exist(plineFig, 'file') || overwrite || overwriteImages
    close all
    set(gcf, 'visible', 'off')
    disp(['Saving PIV XY pathline image to disk: ' plineFig])
    if ~computed_XY 
        load(plineXY, 'pivPathlines')
        XX = pivPathlines.XX ;
        YY = pivPathlines.YY ;
        t0 = pivPathlines.t0 ;
        tIdx0 = pivPathlines.tIdx0 ;
        Lx = pivPathlines.Lx ;
        Ly = pivPathlines.Ly ;
        computed_XY = true ;
    end
    
    % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
    qsubX = round(size(XX, 2) / nY2plot * QS.a_fixed) ;
    if doubleCovered
        qsubY = round(size(XX, 3) / nY2plot * 2);
    else
        qsubY = round(size(XX, 3) / nY2plot) ;
    end
    colors = parula(length(timePoints)) ;
    clf
    for qq = fliplr(1:length(timePoints))
        xx = XX(qq, 1:qsubX:end, 1:qsubY:end) ;
        yy = YY(qq, 1:qsubX:end, 1:qsubY:end) ;
        scatter(xx(:) / double(Lx(qq)), ...
            (1 - (yy(:) / double(Ly(qq))) - Yoffset) / (1 - 2*Yoffset), ...
            scatterSz, ...
            'markerfacecolor', colors(qq, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.3) ;
        hold on;
    end
    axis equal
    ylim([0, 1])
    xlim([0, 1])
    daspect([1 QS.a_fixed 1])
    title('pullback pathlines: PIV coordinates', 'Interpreter', 'Latex')
    ylabel('circumferential position, $\phi / 2\pi$', ...
        'Interpreter', 'Latex')
    xlabel('ap position, $\zeta$', 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([(timePoints(1) - t0) * QS.timeInterval, ...
           (timePoints(end) - t0) * QS.timeInterval]) ;
    ylabel(cb, ['time [' QS.timeUnits ']'], ...
        'Interpreter', 'Latex')
    saveas(gcf, plineFig)
else
    disp(["PIV (u',v') pathline image already on disk: " plineFig])
end

%% Make movie -- PIV evaluation grid
for tidx = 1:length(timePoints)
    set(gcf, 'visible', 'off')
    if mod(tidx, 10) == 0
        disp(['t = ', num2str(timePoints(tidx))])
    end
    fn = fullfile(sprintf(pathlineDir.XY, t0), ...
        [sprintf(QS.fileBase.name, timePoints(tidx)) '.png']) ;
    
    if ~exist(fn, 'file') || overwrite || overwriteImages
        % Load pathlines if not already in RAM
        if ~computed_XY 
            load(plineXY, 'pivPathlines')
            XX = pivPathlines.XX ;
            YY = pivPathlines.YY ;
            t0 = pivPathlines.t0 ;
            tIdx0 = pivPathlines.tIdx0 ;
            Lx = pivPathlines.Lx ;
            Ly = pivPathlines.Ly ;
            computed_XY = true ;
        end
    
        clf ;
        set(gcf, 'visible', 'off')
        movieScatterSz = 1 ;
        hsc = scatter(XX(tidx, :) / double(Lx(tidx)), ...
            (1 - (YY(tidx, :) / double(Ly(tidx))) - Yoffset) / (1 - 2*Yoffset), ...
            movieScatterSz, 'markerfacecolor', QS.plotting.colors(1, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 1.0) ;
        axis equal
        ylim([0, 1])
        xlim([0, 1])
        daspect([1 QS.a_fixed 1])
        time_in_units = (timePoints(tidx) - t0)* QS.timeInterval ;
        tstr = ['$t = $' num2str(time_in_units) ' ' QS.timeUnits] ;
        title(['pullback pathlines, ' tstr], 'Interpreter', 'Latex')
        ylabel('circumferential position, $\phi / 2\pi$', ...
            'Interpreter', 'Latex')
        xlabel('ap position, $\zeta$', 'Interpreter', 'Latex')
        saveas(gcf, fn)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh vertex coordinates in XY pixel space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plinevXY = sprintf(QS.fileName.pathlines_uvprime.vXY, t0) ;
if ~exist(plinevXY, 'file') || overwrite || overwriteImages 
    % Create pathlines emanating from vertex positions (vX,vY) at t=t0
    % Load pullback mesh vertex positions to get Lx and Ly (will be done
    % later in pullbackPathlines() if not now, so no extra cost to do it
    % now instead). 
    % Load 'initial' positions for pathlines to intersect at t=t0
    uvpPIVfn = fullfile(QS.dir.piv_uvprime, 'piv_results.mat') ;
    piv = load(uvpPIVfn) ;
    
    % Additionally smooth the piv output by sigma
    if smoothing_sigma > 0
        disp(['Smoothing piv output with sigma=' num2str(smoothing_sigma)])
        for tidx = 1:length(QS.xp.fileMeta.timePoints)-1
            velx = piv.u_filtered{tidx} ;
            vely = piv.v_filtered{tidx} ;
            piv.u_filtered{tidx} = imgaussfilt(velx, smoothing_sigma) ;
            piv.v_filtered{tidx} = imgaussfilt(vely, smoothing_sigma) ;
        end
        disp('done smoothing')
    end
    
    % Create pathlines emanating from (x0,y0) at t=t0
    im0 = imread(sprintf(QS.fullFileBase.im_uvprime_e, t0)) ;
    % for now assume all images are the same size
    disp('for now assuming all images are the same size')
    Lx = size(im0, 1) * ones(length(timePoints), 1) ;
    Ly = size(im0, 2) * ones(length(timePoints), 1) ;
    
    disp('Loading uvprime mesh at t0 to obtain vertices for advection')
    mesh0 = load(sprintf(QS.fullFileBase.uvpcutMesh, t0), 'uvpcutMesh') ;
    mesh0 = mesh0.uvpcutMesh.resampled ;
    umax0 = max(mesh0.u(:, 1)) ;
    vmax0 = max(mesh0.u(:, 2)) ;
    m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], mesh0.u, doubleCovered, ...
        umax0, vmax0) ;
    m0X = m0XY(:, 1) ;
    m0Y = m0XY(:, 2) ;
    
    % Create pathlines emanating from vertex positions (vX,vY) at t=t0
    options.piv = piv ;
    options.Lx = Lx ;
    options.Ly = Ly ;
    [vX, vY] = QS.pullbackPathlines(m0X, m0Y, t0, options) ;
    
    vUV = QS.XY2uv([Lx(1), Ly(1)], [vX(:), vY(:)], doubleCovered, 1.0, 1.0) ;
    
    % Reshape into grid
    vX = reshape(vX, [length(timePoints), QS.nU, QS.nV]) ;
    vY = reshape(vY, [length(timePoints), QS.nU, QS.nV]) ;
    vU = reshape(vUV(:, 1), [length(timePoints), QS.nU, QS.nV]) ;
    vV = reshape(vUV(:, 2), [length(timePoints), QS.nU, QS.nV]) ;
    
    % maximum u' should be 1.0
    assert(max(vU(:)) < 1.05)
    
    % UNIQUE TO UVPRIME: Also compute relaxed affine coordinates
    % (conformal)
    ars = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
    disp('Loading relaxation factors for affine stretch -> conformal maps')
    for tidx = 1:length(QS.xp.fileMeta.timePoints) 
        tp = QS.xp.fileMeta.timePoints(tidx) ;
        uvpcutMeshFn = sprintf(QS.fullFileBase.uvpcutMesh, tp) ;
        tmp = load(uvpcutMeshFn) ;
        ars(tidx) = tmp.uvpcutMesh.raw.ar ;
    end
    vUa = ars .* vU ;
    
    ars_fn = fullfile(sprintf(QS.dir.pathlines_uvprime.data, t0), ...
        'affineRelaxFactors.png') ;
    plot(QS.xp.fileMeta.timePoints, ars)
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex') ;
    ylabel('affine factor, $a$', 'interpreter', 'latex')
    saveas(gcf, ars_fn)
    
    % Save pathlines of mesh locations u in XY coords, vXY
    vertexPathlines = struct() ;
    vertexPathlines.vX = vX ;
    vertexPathlines.vY = vY ;
    vertexPathlines.vU = vU ;
    vertexPathlines.vV = vV ;
    vertexPathlines.vU_affine = vUa ;
    vertexPathlines.affineRelaxFactors = ars ;
    vertexPathlines.t0 = t0 ;
    vertexPathlines.tIdx0 = tIdx0 ;
    vertexPathlines.Lx = Lx ;
    vertexPathlines.Ly = Ly ;
    save(plinevXY, 'vertexPathlines')
    
    computed_vXY = true ;
else
    disp(['Mesh vertex XY pathlines already on disk: ' plinevXY])
    computed_vXY = false ;
end

% Save image of result
plineFig1 = [plinevXY(1:end-4) '.png'] ;
plineFig2 = [plinevXY(1:end-4) '_relaxed.png'] ;
if ~exist(plineFig1, 'file') || ~exist(plineFig2, 'file') || overwrite || ...
        overwriteImages
    disp(['Saving mesh vertex XY pathline image to disk: ' plineFig1])
    if ~computed_vXY 
        load(sprintf(QS.fileName.pathlines_uvprime.vXY, t0), 'vertexPathlines')
        vU = vertexPathlines.vU ;
        vUa = vertexPathlines.vU_affine ;
        vV = vertexPathlines.vV ;
        t0 = vertexPathlines.t0 ;
        tIdx0 = vertexPathlines.tIdx0 ;
        Lx = vertexPathlines.Lx ;
        Ly = vertexPathlines.Ly ;
    end
    
    % Plot both (u'v') and (u'_relaxed, v')
    for pp = 1:2
        
        % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
        qsubX = round(size(vU, 2) / nY2plot) ;
        qsubY = round(size(vU, 3) / nY2plot * QS.a_fixed) ;
        colors = parula(length(timePoints)) ;
        clf
        for qq = fliplr(1:length(timePoints))
            if pp == 1
                xx = vU(qq, 1:qsubX:end, 1:qsubY:end) ;
                yy = vV(qq, 1:qsubX:end, 1:qsubY:end) ;
            else
                xx = vUa(qq, 1:qsubX:end, 1:qsubY:end) ;
                yy = vV(qq, 1:qsubX:end, 1:qsubY:end) ;
            end
            scatter(xx(:), yy(:), ...
                scatterSz, ...
                'markerfacecolor', colors(qq, :), ...
                'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.5) ;
            hold on;
            pause(0.1) ;
        end
        
        if pp == 1
            axis equal
            ylim([0, 1])
            xlim([0, 1])
            daspect([1 QS.a_fixed 1])
        end
        title("$(u',v')$ pullback pathlines: mesh vertices", 'Interpreter', 'Latex')
        ylabel("circumferential position, $v'$", ...
            'Interpreter', 'Latex')
        xlabel("axial position, $u'$", 'Interpreter', 'Latex')
        cb = colorbar() ;
        caxis([(timePoints(1) - t0) * QS.timeInterval, ...
               (timePoints(end) - t0) * QS.timeInterval]) ;
        ylabel(cb, ['time [' QS.timeUnits ']'], ...
            'Interpreter', 'Latex')
        if pp == 1
            saveas(gcf, plineFig1)
        else
            saveas(gcf, plineFig2)
        end
    end
else
    disp(['Mesh vertex XY pathline image already on disk: ' plineFig])
end

% Make movie -- vertex pathlines
for tidx = 1:length(timePoints)    
    if mod(tidx, 10) == 0
        disp(['t = ', num2str(timePoints(tidx))])
    end
    fn = fullfile(sprintf(pathlineDir.vXY, t0), ...
        [sprintf(QS.fileBase.name, timePoints(tidx)) '.png']) ;
    
    if ~exist(fn, 'file') || overwrite || overwriteImages
        % Load pathlines if not already in RAM
        if ~computed_vXY 
            load(plinevXY, 'vertexPathlines')
            vUa = vertexPathlines.vU_affine ;
            vV = vertexPathlines.vV ;
            t0 = vertexPathlines.t0 ;
            tIdx0 = vertexPathlines.tIdx0 ;
            ars = vertexPathlines.affineRelaxFactors ;
            Lx = vertexPathlines.Lx ;
            Ly = vertexPathlines.Ly ;
            maxY = max(vV(:)) ;
            minY = min(vV(:)) ;
            computed_vXY = true ;
        elseif tidx == 1
            maxY = max(vV(:)) ;
            minY = min(vV(:)) ;            
        end
    
        clf ;
        movieScatterSz = 1 ;
        % scatter(vX(tidx, :) / double(Lx(tidx)), ...
        %     (1 - (vY(tidx, :) / double(Ly(tidx))) - Yoffset) / (1 - 2*Yoffset), ...
        %     movieScatterSz, 'markerfacecolor', QS.plotting.colors(1, :), ...
        %     'markeredgecolor', 'none', 'MarkerFaceAlpha', 1.0) ;
        scatter(vUa(tidx, :), vV(tidx, :), ...
            movieScatterSz, 'markerfacecolor', QS.plotting.colors(1, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 1.0) ;
        axis equal
        ylim([max(-0.5, minY), min(1.5, maxY)])
        xlim([0, max(ars)])
        % axis equal
        % daspect([1 QS.a_fixed 1])
        time_in_units = (timePoints(tidx) - t0)* QS.timeInterval ;
        tstr = ['$t = $' num2str(time_in_units) ' ' QS.timeUnits] ;
        title(['$(u'',v'')$ pullback pathlines, ' tstr], 'Interpreter', 'Latex')
        ylabel('circumferential position, $v''$', ...
            'Interpreter', 'Latex')
        yticks([0, 1])
        xlabel("axial position, $u'$", 'Interpreter', 'Latex')
        saveas(gcf, fn)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh face barycenters in XY pixel space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineFXY = sprintf(QS.fileName.pathlines_uvprime.fXY, t0) ;
if ~exist(plineFXY, 'file') || overwrite
    % Create pathlines emanating from barycenters (bcX,bcY) at t=t0
    % Load pullback mesh vertex positions to get Lx and Ly (will be done
    % later in pullbackPathlines() if not now, so no extra cost to do it
    % now instead).
    
    % Load 'initial' positions for pathlines to intersect at t=t0
    uvpPIVfn = fullfile(QS.dir.piv_uvprime, 'piv_results.mat') ;
    piv = load(uvpPIVfn) ;
    
    % Additionally smooth the piv output by sigma
    if smoothing_sigma > 0
        disp(['Smoothing piv output with sigma=' num2str(smoothing_sigma)])
        for tidx = 1:length(QS.xp.fileMeta.timePoints)-1
            velx = piv.u_filtered{tidx} ;
            vely = piv.v_filtered{tidx} ;
            piv.u_filtered{tidx} = imgaussfilt(velx, smoothing_sigma) ;
            piv.v_filtered{tidx} = imgaussfilt(vely, smoothing_sigma) ;
        end
        disp('done smoothing')
    end
    
    % Create pathlines emanating from (x0,y0) at t=t0
    im0 = imread(sprintf(QS.fullFileBase.im_uvprime_e, t0)) ;
    % for now assume all images are the same size
    disp('for now assuming all images are the same size')
    Lx = size(im0, 1) * ones(length(timePoints), 1) ;
    Ly = size(im0, 2) * ones(length(timePoints), 1) ;
    
    mesh0 = load(sprintf(QS.fullFileBase.uvpcutMesh, t0), 'uvpcutMesh') ;
    mesh0 = mesh0.uvpcutMesh.resampled ;
    umax0 = max(mesh0.u(:, 1)) ;
    vmax0 = max(mesh0.u(:, 2)) ;
    m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], mesh0.u, doubleCovered, umax0, vmax0) ;
    
    % Obtain barycenters of faces
    f0XY = barycenter(m0XY, mesh0.f) ;  % this is in gptoolbox
    f0X = f0XY(:, 1) ;
    f0Y = f0XY(:, 2) ;
    
    % Create pathlines emanating from barycenters (bcX,bcY) at t=t0
    options.piv = piv ;
    options.Lx = Lx ;
    options.Ly = Ly ;
    [fX, fY] = QS.pullbackPathlines(f0X, f0Y, t0, options) ;
        
    fUV = QS.XY2uv([Lx(1), Ly(1)], [fX(:), fY(:)], doubleCovered, 1.0, 1.0) ;
    
    % Reshape into grid
    fU = reshape(fUV(:, 1), [length(timePoints), size(f0XY, 1) ]) ;
    fV = reshape(fUV(:, 2), [length(timePoints), size(f0XY, 1) ]) ;
    
    % UNIQUE TO UVPRIME: Also compute relaxed affine coordinates
    % (conformal)
    ars = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
    disp('Loading relaxation factors for affine stretch -> conformal maps')
    for tidx = 1:length(QS.xp.fileMeta.timePoints) 
        tp = QS.xp.fileMeta.timePoints(tidx) ;
        uvpcutMeshFn = sprintf(QS.fullFileBase.uvpcutMesh, tp) ;
        tmp = load(uvpcutMeshFn) ;
        ars(tidx) = tmp.uvpcutMesh.raw.ar ;
    end
    fUa = ars .* fU ;
    
    ars_fn = fullfile(sprintf(QS.dir.pathlines_uvprime.data, t0), ...
        'affineRelaxFactors_faces.png') ;
    plot(QS.xp.fileMeta.timePoints, ars)
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex') ;
    ylabel('affine factor, $a$', 'interpreter', 'latex')
    saveas(gcf, ars_fn)
    
    % Save pathlines of mesh face barycenters f in XY coords, fXY
    facePathlines = struct() ;
    facePathlines.fX = fX ;
    facePathlines.fY = fY ;
    facePathlines.fU = fU ;
    facePathlines.fV = fV ;
    facePathlines.fU_affine = fUa ;
    facePathlines.t0 = t0 ;
    facePathlines.tIdx0 = tIdx0 ;
    facePathlines.affineRelaxFactors = ars ;
    facePathlines.Lx = Lx ;
    facePathlines.Ly = Ly ;
    save(plineFXY, 'facePathlines')
    
    computed_fXY = true ;
else
    computed_fXY = false ;
end

% Save image of result
plineFig = [plineFXY(1:end-4) '.png'] ;
if ~exist(plineFig, 'file') || overwrite || overwriteImages
    if ~computed_fXY 
        load(QS.fileName.pathlines_uvprime.fXY, 'facePathlines')
        fX = facePathlines.fX ;
        fY = facePathlines.fY ;
        t0 = facePathlines.t0 ;
        tIdx0 = facePathlines.tIdx0  ;
    end
    
    % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
    qsubX = round(sqrt(size(fX, 2)) / nY2plot * 8) ;
    colors = parula(length(timePoints)) ;
    clf
    for qq = fliplr(1:length(timePoints))
        xx = fX(qq, 1:qsubX:end) ;
        yy = fY(qq, 1:qsubX:end) ;
        scatter(xx(:) / double(Lx(qq)), ...
            (1 - (yy(:) / double(Ly(qq))) - Yoffset) / (1 - 2*Yoffset), ...
            scatterSz, ...
            'markerfacecolor', colors(qq, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.3) ;
        hold on;
    end
    axis equal
    ylim([0, 1])
    xlim([0, 1])
    daspect([1 QS.a_fixed 1])
    title('pullback pathlines: mesh vertices', 'Interpreter', 'Latex')
    ylabel('circumferential position, $\phi / 2\pi$', ...
        'Interpreter', 'Latex')
    xlabel('ap position, $\zeta$', 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([(timePoints(1) - t0) * QS.timeInterval, ...
           (timePoints(end) - t0) * QS.timeInterval]) ;
    ylabel(cb, ['time [' QS.timeUnits ']'], ...
        'Interpreter', 'Latex')
    saveas(gcf, plineFig)
else
    disp(['Mesh face barycenter XY pathline image already on disk: ' plineFig])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Push forward from pullback to embedding coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineXYZ = sprintf(QS.fileName.pathlines_uvprime.XYZ, t0) ;
if ~exist(plineXYZ, 'file') || overwrite
    disp('Compute 3d pathlines in embedding space for PIV eval coords')
    % Build 3d pathlines in embedding space for all vertex locations
    % Load 2d pathlines in pullback space if not already in RAM
    if ~computed_XY 
        load(plineXY, 'pivPathlines')
        XX = pivPathlines.XX ;
        YY = pivPathlines.YY ;
        t0 = pivPathlines.t0 ;
        tIdx0 = pivPathlines.tIdx0 ;
        Lx = pivPathlines.Lx ;
        Ly = pivPathlines.Ly ;
        computed_XY = true ;
    end
    
    % Preallocate v3d = (XX, YY, ZZ)
    X3 = zeros(size(XX)) ;
    Y3 = zeros(size(XX)) ;
    Z3 = zeros(size(XX)) ;
    fieldfacesCell = cell(length(timePoints), 1) ;
    
    % For each timepoint, push forward into embedding space
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        if mod(tidx, 10) == 0
            disp(['t = ', num2str(tp)])
        end
        
        % Load this timepoint's cutMesh
        mesh0 = load(sprintf(QS.fullFileBase.uvpcutMesh, tp), ...
                             'uvpcutMesh') ;
        % NOTE: use raw meshes not resampled
        mesh0 = mesh0.uvpcutMesh.raw ;
        umax0 = max(mesh0.u(:, 1)) ;
        vmax0 = max(mesh0.u(:, 2)) ;
        tileCount = [2 2] ;
        [tm0f, tm0v2d, tm0v3d] = tileAnnularCutMesh(mesh0, tileCount);
        tm0XY = QS.uv2XY([Lx(tidx) Ly(tidx)], tm0v2d, ...
                         doubleCovered, umax0, vmax0) ;
        Xeval = XX(tidx, :, :) ;
        Yeval = YY(tidx, :, :) ;
        [pt0, fieldfaces] = interpolate2Dpts_3Dmesh(tm0f, tm0XY,...
                                        tm0v3d, [Xeval(:), Yeval(:)]) ;
        % Store in grid 
        X3(tidx, :, :) = reshape(pt0(:, 1), [size(X3, 2), size(X3, 3)]) ;
        Y3(tidx, :, :) = reshape(pt0(:, 2), [size(Y3, 2), size(Y3, 3)]) ;
        Z3(tidx, :, :) = reshape(pt0(:, 3), [size(Z3, 2), size(Z3, 3)]) ;
        fieldfacesCell{tidx} = fieldfaces ;
    end
    
    % Save pathlines of PIV evaluation coords in embedding coords, XYZ
    piv3dPathlines = struct() ;
    piv3dPathlines.XX = X3 ;
    piv3dPathlines.YY = Y3 ;
    piv3dPathlines.ZZ = Z3 ;
    piv3dPathlines.t0 = t0 ;
    piv3dPathlines.tIdx0 = tIdx0 ;
    save(plineXYZ, 'piv3dPathlines')
end                             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Push forward vertex pathlines from pullback to embedding coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plinev3d = sprintf(QS.fileName.pathlines_uvprime.v3d, t0) ;
if ~exist(plinev3d, 'file') || overwrite || overwriteImages
    disp('Computing 3d pathlines in embedding space for Lagrangian/advected uvprimecutMesh vertex coords')
    % Build 3d pathlines in embedding space for all vertex locations
    % Load 2d pathlines in pullback space if not already in RAM
    if ~computed_vXY 
        load(sprintf(QS.fileName.pathlines_uvprime.vXY, t0), 'vertexPathlines')
        vX = vertexPathlines.vX ;
        vY = vertexPathlines.vY ;
        t0 = vertexPathlines.t0 ;
        tIdx0 = vertexPathlines.tIdx0 ;
        Lx = vertexPathlines.Lx ;
        Ly = vertexPathlines.Ly ;
        computed_vXY = true ;
    end
    
    % Preallocate v3d = (XX, YY, ZZ)
    vX3 = zeros(size(vX)) ;
    vY3 = zeros(size(vX)) ;
    vZ3 = zeros(size(vX)) ;
    vX3rs = zeros(size(vX)) ;
    vY3rs = zeros(size(vX)) ;
    vZ3rs = zeros(size(vX)) ;
    fieldfacesCell = cell(length(timePoints), 1) ;
    
    % For each timepoint, push advected vertices forward into embedding space
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        if mod(tidx, 10) == 0
            disp(['t = ', num2str(tp)])
        end
        
        % Load this timepoint's uvprime cutMesh
        mesh0 = load(sprintf(QS.fullFileBase.uvpcutMesh, tp), ...
                             'uvpcutMesh') ;
        mesh0 = mesh0.uvpcutMesh.raw ;
        umax0 = max(mesh0.u(:, 1)) ;
        vmax0 = max(mesh0.u(:, 2)) ;
        tileCount = [2 2] ;
        [tm0f, tm0v2d, tm0v3d] = tileAnnularCutMesh(mesh0, tileCount);
        tm0XY = QS.uv2XY([Lx(tidx) Ly(tidx)], tm0v2d, ...
                         doubleCovered, umax0, vmax0) ;
                     
        % Evaluate at advected Lagrangian-frame vertex positions
        Xeval = vX(tidx, :, :) ;
        Yeval = vY(tidx, :, :) ;
        
        % Adjust the XY positions so that none are out of bounds
        epsilon = 1e-7 ;
        minXallowed = 1 + epsilon ;
        Xeval = Xeval(:) ;
        Yeval = Yeval(:) ;
        Xeval(Xeval(:, 1) < minXallowed) = minXallowed ;
        Xeval(Xeval(:, 1) > Lx(tidx) - epsilon) = Lx(tidx) - epsilon ;
        
        [pt0, fieldfaces] = interpolate2Dpts_3Dmesh(tm0f, tm0XY,...
                                        tm0v3d, [Xeval(:), Yeval(:)]) ;
        assert(all(isfinite(pt0(:)))) 
        
        % Check that points are inside
        % plot(tm0XY(:, 1), tm0XY(:, 2), '-')
        % hold on ;
        % plot(Xeval(:), Yeval(:), '.')
        
        % Rotate and scale    
        vXYZrs = QS.xyz2APDV(pt0) ;
        
        % Store in grid 
        vX3(tidx, :, :) = reshape(pt0(:, 1), [size(vX3, 2), size(vX3, 3)]) ;
        vY3(tidx, :, :) = reshape(pt0(:, 2), [size(vY3, 2), size(vY3, 3)]) ;
        vZ3(tidx, :, :) = reshape(pt0(:, 3), [size(vZ3, 2), size(vZ3, 3)]) ;
        vX3rs(tidx, :, :) = reshape(vXYZrs(:, 1), [size(vX3, 2), size(vX3, 3)]) ;
        vY3rs(tidx, :, :) = reshape(vXYZrs(:, 2), [size(vY3, 2), size(vY3, 3)]) ;
        vZ3rs(tidx, :, :) = reshape(vXYZrs(:, 3), [size(vZ3, 2), size(vZ3, 3)]) ;
        fieldfacesCell{tidx} = fieldfaces ;
    end
    
    % Save pathlines of PIV evaluation coords in embedding coords, XYZ
    v3dPathlines = struct() ;
    v3dPathlines.vX = vX3 ;
    v3dPathlines.vY = vY3 ;
    v3dPathlines.vZ = vZ3 ;
    v3dPathlines.vXrs = vX3rs ;
    v3dPathlines.vYrs = vY3rs ;
    v3dPathlines.vZrs = vZ3rs ;
    v3dPathlines.t0 = t0 ;
    v3dPathlines.tIdx0 = tIdx0 ;
    save(plinev3d, 'v3dPathlines')
    computed_v3d = true ;
else
    computed_v3d = false ;
end

%% Plot pathlines in 3d
% Save image of result
[~,~,~,xyzlims] = QS.getXYZLims() ;
nTimePoints = 20 ;
first = true ;
slab = 20 ;
black_figs = true ;

if ~computed_v3d 
    load(sprintf(QS.fileName.pathlines_uvprime.v3d, t0), 'v3dPathlines')
    vX3rs = v3dPathlines.vXrs ;
    vY3rs = v3dPathlines.vYrs ;
    vZ3rs = v3dPathlines.vZrs ;
    t0 = v3dPathlines.t0 ;
    tIdx0 = v3dPathlines.tIdx0 ;
    computed_v3d = false ;
end

% Get subsampling
qsubX = round(size(vX3rs, 2) / nY2plot) - 1 ;

% Obtain indices for subsampling
idmat = false([size(vX3rs, 2), size(vX3rs, 3)]) ;
% This would make gridlike sampling
% ap2do = 1:qsubX:(size(vX3rs, 2)-1)

% Saggital sampling
halfV = floor(QS.nV*0.5) ;
ap2do = [1, halfV] ;  % 2, QS.nV-1,halfV-1, halfV, halfV + 1] ;
for qq = ap2do 
    idmat(1:qsubX:end, qq) = true; 
end
inds = find(idmat) ;
clearvars idmat
xx = vX3rs(tIdx0, inds) ;
yy = vY3rs(tIdx0, inds) ;
zz = vZ3rs(tIdx0, inds) ;

% % Define the lagrangian points that lie in saggital plane at t=0
% yslab = find(abs(yy)<slab) ;
% plot3(xx(yslab), yy(yslab), zz(yslab), '.')

for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    
    plineFig = fullfile(sprintf(pathlineDir.v3d, t0), ...
        [sprintf(QS.fileBase.name, tp) '.png']) ;
    
    if ~exist(plineFig, 'file') || overwrite
        % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
        close all
        set(gcf, 'visible', 'off')
        tp2plot = max(1, tidx-nTimePoints):tidx ;
        
        if black_figs
            colors = bone(length(tp2plot)) ;
            colormap(bone)
        else
            colors = flipud(bone(length(tp2plot))) ;
            colormap(flipud(bone))
        end
        if size(colors, 1) == 1
            colors = bone(2) ;
            colors = colors(2, :) ;
        end
        
        alphas = linspace(0, 1, length(tp2plot)) ;
        for qind = 1:length(tp2plot)
            qq = tp2plot(qind) ;
            xx = vX3rs(qq, inds) ;
            yy = vY3rs(qq, inds) ;
            zz = vZ3rs(qq, inds) ;
            % xx = xx(yslab) ;
            % yy = yy(yslab) ;
            % zz = zz(yslab) ;
            scatter3(xx(:), yy(:), zz(:), ...
                scatterSz * 2, ...
                'markerfacecolor', colors(qind, :), ...
                'markeredgecolor', 'none', 'markerfacealpha', alphas(qind)) ;
            % s.MarkerFaceAlpha = (max(1.0 - 0.1 * yy(:), 0)) ;
            hold on;
        end
        view(0,0)
        axis equal
        scatter([], [])
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))
        % title('pathlines: mesh vertices', 'Interpreter', 'Latex')
        xlabel(['ap position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        ylabel(['lateral position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        zlabel(['dv position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        
        axis off
        
        if black_figs
            set(gcf, 'color', 'k')
            set(gca, 'color', 'k', 'xcol', 'k', 'ycol', 'k', 'zcol', 'k')
        else
            set(gca, 'color', 'w')
        end

        if black_figs
            cb = colorbar('eastOutside', 'color', 'w') ;
        else
            cb = colorbar('eastOutside') ;
        end
        
        % Adjust size of colorbar
        cbheight = cb.Position(4) ;
        cb.Position(1) = cb.Position(1) + 0.02 ;
        cb.Position(4) = cb.Position(4) * 0.6 ;
        cb.Position(2) = cb.Position(2) + 0.15 * cbheight ;
        ax = gca ;
        ax.Position(1) = ax.Position(1) - 0.05 ; 
        
        caxis([(QS.xp.fileMeta.timePoints(max(tp2plot))-t0-nTimePoints) * QS.timeInterval, ...
           (QS.xp.fileMeta.timePoints(max(tp2plot))-t0) * QS.timeInterval]) ;
        ylabel(cb, ['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')

        % saveas(gcf, plineFig)
        export_fig(plineFig, '-nocrop', '-r150')
        close all
    else
        disp(['Mesh advected vertex XYZ pathline image already on disk: ' plineFig])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Push forward vertex pathlines from pullback to embedding coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pliner3d = sprintf(QS.fileName.pathlines_uvprime.radius, t0) ;
if ~exist(pliner3d, 'file') || overwrite || overwriteImages
    disp('Computing radii along pathlines in embedding space for Lagrangian/advected uvprimecutMesh vertex coords')
    % Build 3d pathlines in embedding space for all vertex locations
    % Load 2d pathlines in pullback space if not already in RAM
    if ~computed_vXY 
        load(sprintf(QS.fileName.pathlines_uvprime.vXY, t0), 'vertexPathlines')
        vX = vertexPathlines.vX ;
        vY = vertexPathlines.vY ;
        t0 = vertexPathlines.t0 ;
        tIdx0 = vertexPathlines.tIdx0 ;
        Lx = vertexPathlines.Lx ;
        Ly = vertexPathlines.Ly ;
        computed_vXY = true ;
    end
    
    % Preallocate v3d = (XX, YY, ZZ)
    vRad = zeros(size(vX)) ;
    
    % For each timepoint, push advected vertices forward into embedding space
    for tidx = 1:length(timePoints)
        tp = timePoints(tidx) ;
        if mod(tidx, 10) == 0
            disp(['t = ', num2str(tp)])
        end
        
        % Load this timepoint's uvprime cutMesh
        mesh0 = load(sprintf(QS.fullFileBase.uvpcutMesh, tp), ...
                             'uvpcutMesh') ;
        mesh0 = mesh0.uvpcutMesh.raw ;
        rad0 = mesh0.radius_um(:) ;
        umax0 = max(mesh0.u(:, 1)) ;
        vmax0 = max(mesh0.u(:, 2)) ;
        nU = mesh0.nU ;
        nV = mesh0.nV ;
        tileCount = [1 1] ;
        [tm0f, tm0v2d, tm0v3d] = tileAnnularCutMesh(mesh0, tileCount);
        rad0Tiled = [rad0(1:nU*(nV-1)); rad0; rad0(nU+1:end)] ;
        tm0XY = QS.uv2XY([Lx(tidx) Ly(tidx)], tm0v2d, ...
                         doubleCovered, umax0, vmax0) ;
                     
        % Evaluate at advected Lagrangian-frame vertex positions
        Xeval = vX(tidx, :, :) ;
        Yeval = vY(tidx, :, :) ;
        
        % Adjust the XY positions so that none are out of bounds
        epsilon = 1e-7 ;
        minXallowed = 1 + epsilon ;
        Xeval = Xeval(:) ;
        Yeval = Yeval(:) ;
        Xeval(Xeval(:, 1) < minXallowed) = minXallowed ;
        Xeval(Xeval(:, 1) > Lx(tidx) - epsilon) = Lx(tidx) - epsilon ;
        
        radi = scatteredInterpolant(tm0XY(:, 1), tm0XY(:, 2), rad0Tiled,...
            'linear', 'nearest') ;
        radv = radi(Xeval(:), Yeval(:)) ;
        assert(all(isfinite(radv(:)))) 
        
        % Check that points are inside
        % plot(tm0XY(:, 1), tm0XY(:, 2), '-')
        % hold on ;
        % plot(Xeval(:), Yeval(:), '.')
                
        % Store in grid 
        vRad(tidx, :, :) = reshape(radv, [size(vRad, 2), size(vRad, 3)]) ;
    end
    
    % Save pathlines of PIV evaluation coords in embedding coords, XYZ
    vRadiusPathlines = struct() ;
    vRadiusPathlines.radii = vRad ;
    vRadiusPathlines.t0 = t0 ;
    vRadiusPathlines.tIdx0 = tIdx0 ;
    disp(['Saving pathline radii to ', pliner3d])
    save(pliner3d, 'vRadiusPathlines')
    computed_rad = true ;
else
    computed_rad = false ;
end

if ~exist(sprintf(pathlineDir.radius, t0), 'dir')
    mkdir(sprintf(pathlineDir.radius, t0))
end

%% Compute radius kymograph along AP
kymoDir = sprintf(QS.dir.pathlines_uvprime.kymographs, t0) ;
if ~exist(kymoDir, 'dir')
    mkdir(kymoDir) 
end
plinerAP = sprintf(QS.fileName.pathlines_uvprime.kymographs.radius, t0) ;
if ~exist(plinerAP, 'file') || overwrite || overwriteImages
    if ~computed_rad
        load(pliner3d, 'vRadiusPathlines')
        computed_rad = true ;
    end
    timePoints = QS.xp.fileMeta.timePoints;
    nU = size(vRadiusPathlines.radii, 2) ;
    nV = size(vRadiusPathlines.radii, 3) ;
    radius_apM = zeros(length(timePoints), nU) ;
    radius_dM = zeros(length(timePoints), nU) ;
    radius_vM = zeros(length(timePoints), nU) ;
    radius_lM = zeros(length(timePoints), nU) ;
    radius_rM = zeros(length(timePoints), nU) ;
    for tidx = 1:length(timePoints)
        radi = squeeze(vRadiusPathlines.radii(tidx, :, :)) ;
        [dorsal, ventral, left, right] = QS.quarterIndicesDV(nV) ;
        radius_apM(tidx, :) = mean(radi(:, 1:nV-1), 2) ;
        radius_dM(tidx, :) = mean(radi(:, dorsal), 2) ;
        radius_vM(tidx, :) = mean(radi(:, ventral), 2) ;
        radius_lM(tidx, :) = mean(radi(:, left), 2) ;
        radius_rM(tidx, :) = mean(radi(:, right), 2) ;
    end
    
    save(plinerAP, 'radius_apM', 'radius_dM', 'radius_vM', ...
        'radius_lM', 'radius_rM')
    
    % Plot each kymo
    kymos = {radius_apM, radius_dM, radius_vM, radius_lM, radius_rM} ;
    titles = {'$(u'',v'')$ pathline radii', ...
        '$(u'',v'')$ pathline radii, dorsal side', ...
        '$(u'',v'')$ pathline radii, ventral side', ...
        '$(u'',v'')$ pathline radii, left side', ...
        '$(u'',v'')$ pathline radii, right side'} ;
    fns = {[plinerAP(1:end-4) '_ap.png'], ...
        [plinerAP(1:end-4) '_d.png'], ...
        [plinerAP(1:end-4) '_v.png'], ...
        [plinerAP(1:end-4) '_l.png'], ...
        [plinerAP(1:end-4) '_r.png']} ;
    uspace = linspace(0, 1, nU) ;
    for qq = 1:length(kymos)
        if strcmpi(QS.timeUnits, 'min')
            imagesc(uspace, (timePoints - vRadiusPathlines.t0)/60, kymos{qq})
        else
            imagesc(uspace, timePoints - vRadiusPathlines.t0, radius_apM)
        end
        xlabel('ap position, $u$', 'interpreter', 'latex')
        if strcmpi(QS.timeUnits, 'min')
            ylabel('time, [hr]', 'interpreter', 'latex')
        else
            ylabel(['time, [' QS.timeUnits ']'], 'interpreter', 'latex')
        end
        title(titles{qq}, 'interpreter', 'latex')
        cb = colorbar() ;
        ylabel(cb, ['radius, [' QS.spaceUnits ']'], 'interpreter', 'latex')
        saveas(gcf, fns{qq})
    end
end


% Plot radius in 3d
for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    
    plineRadFig = fullfile(sprintf(pathlineDir.radius, t0), ...
        [sprintf(QS.fileBase.name, tp) '.png']) ;
    
    if ~exist(plineRadFig, 'file') || overwrite
        if ~computed_rad
            load(pliner3d, 'vRadiusPathlines')
            computed_rad = true ;
        end
        
        % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
        close all
        set(gcf, 'visible', 'off')
        clf
        
        xx = vX3rs(tidx, :) ;
        yy = vY3rs(tidx, :) ;
        zz = vZ3rs(tidx, :) ;
        opts = struct() ;
        opts.axisOff = true ;
        opts.style = 'positive' ;
        opts.axPosition = [0.1 0.11 0.85 0.8] ;
        opts.fig = gcf ;
        opts.sscale = max(vRadiusPathlines.radii(:)) ;
        [~, cb] = scalarFieldOnSurface(refMesh.f, [xx(:), yy(:), zz(:)], ...
            vRadiusPathlines.radii(tidx, :), opts) ;
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))
        view(0, 0) 
        xlabel(['ap position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        ylabel(['lateral position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        zlabel(['dv position [' QS.spaceUnits ']'], 'Interpreter', 'Latex')
        
        % axis off
        
        if black_figs
            set(gcf, 'color', 'k')
            set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w', 'zcol', 'w')
        else
            set(gca, 'color', 'w')
        end

        if black_figs
            set(cb, 'color', 'w') ;
        end
        
        % Adjust size of colorbar
        cbheight = cb.Position(4) ;
        cb.Position(1) = cb.Position(1) - 0.02 ;
        cb.Position(4) = cb.Position(4) ;
        cb.Position(2) = cb.Position(2) + 0.15 * cbheight ;
        ax = gca ;
        ax.Position(1) = ax.Position(1) - 0.05 ; 
        ylabel(cb, ['radius [' QS.spaceUnits ']' ], 'Interpreter', 'Latex')

        saveas(gcf, plineRadFig)
        close all
    else
        disp(['Mesh radii pathline image already on disk: ' plineFig])
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Streamline version has issues but included here for reference:
% % Build pathlines starting at t=t0 
% % NOTE: If timePoints does not have a uniform dt, there is no problem:
% % we define the brick of worldlines with uniform spacing and stepsize=1
% % so that paths follow optical flow between frames and vz is 1: ie, the
% % velocities connect adjacent frames, whatever the dt between them.
% xG = repmat(x0, [1, 1, ntps]) ;
% yG = repmat(y0, [1, 1, ntps]) ;
% % make z grid increasing along third (z) dimension
% zG = ones(size(xG)) ;
% for qq = 1:ntps
%     zG(:, :, qq) = qq ;
% end
% % check each xyzgrid as image
% if preview
%     for qq = 1:ntps
%         imagesc(squeeze(xG(:, :, qq)))
%         title(['t=' num2str(qq)])
%         pause(0.1)
%     end
%     for qq = 1:ntps
%         imagesc(squeeze(yG(:, :, qq)))
%         title(['t=' num2str(qq)])
%         pause(0.1)
%     end
%     for qq = 1:ntps
%         imagesc(squeeze(zG(:, :, qq)))
%         title(['t=' num2str(qq)])
%         caxis([0, max(max(zG(:)), max(xG(:)))])
%         pause(0.1)
%     end
% end
%
% vzdat = ones(size(xG)) ;
% sLineOptions = [1] ;  % stepsize is dt 
% % z0 is the starting time array in the worldline loaf
% z0 = QS.xp.tIdx(t0) * ones(size(x0)) ; 
% 
% strLineM = streamline(xG, yG, zG, ...
%     squeeze(vPIV(:, :, :, 1)), ...
%     squeeze(vPIV(:, :, :, 2)), vzdat, ...
%     x0(:), y0(:), z0(:), sLineOptions) ;
% % reshape pathlines into grid format
% strLineM = reshape(strLineM, [ntps, size(x0, 1), size(x0, 2)], 2) ;

