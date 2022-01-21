function measurePullbackPathlines(QS, options)
%measurePullbackPathlines(QS, options)
%   Measure pathlines of optical flow in pullback space. 
%   These pathlines can then be used to query velocities or other
%   properties that are parameterized by pullback coordinates. 
%   For example, we average velocites along Lagrangian pathlines.
%   For another example, to build a Lagrangian-frame measure of divergence 
%   of tangential velocity, div(v_t), we can query div(v_t) at the coords
%   of the PullbackPathlines (via interpolation of div(v_t) defined on
%   pullback mesh vertices).
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   maxIter : int (default=100)
%       #Ricci steps for refMesh.u_ricci and beltrami calculation
%
% Saves to disk
% -------------
% refMesh: sprintf(QS.fileName.pathlines.refMesh, t0) 
%   reference mesh in 3d and its pullback that is conformal mapping to the 
%   plane, created using Ricci flow at time t0
% piv_pathlines_v3d:XY: sprintf(QS.fileName.pathlines.XY, t0) 
%   Advected PIV evaluation coordinates in XY pixel space make pathlines
%   in XY.
% piv_pathlines_v3d:vXY: sprintf(QS.fileName.pathlines.vXY, t0) 
%   Advect refMesh vertex coordinates in XY pixel space in optical flow to 
%   make pathlines in XY (pixel space).
% piv_pathlines_v3d:fXY: sprintf(QS.fileName.pathlines.fXY, t0) 
%   Advect refMesh face barycenters in XY pixel space to make pathlines
% piv_pathlines_v3d:XYZ: sprintf(QS.fileName.pathlines.XYZ, t0) 
%   Embedding coordinates of advected PIV evaluation coordinates, pushed
%   into APDV coordinate system
% piv_pathlines_v3d:v3d: sprintf(QS.fileName.pathlines.v3d, t0) 
%   Embedding coordinates of advected refMesh vertices, pushed into APDV
%   coordinate system
%   
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
pivimCoords = QS.piv.imCoords ;
timePoints = QS.xp.fileMeta.timePoints ;
samplingResolution = '1x' ;
nY2plot = 30 ;
scatterSz = 2 ;
movieScatterSz = 2 ;
maxIter = 200 ;
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
    t0 = QS.t0set() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
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
if strcmp(pivimCoords(end), 'e')
    doubleCovered = true ;
    Yoffset = 0.25 ;
else
    doubleCovered = false ;
    Yoffset = 0.0 ;
end
if isfield(options, 'maxIter')
    maxIter = options.maxIter ;
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
pathlineDir = QS.dir.pathlines ;
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
% Use Ricci flow for beltrami quasiconformal coefficient mu, etc
refMeshFn = sprintf(QS.fileName.pathlines.refMesh, t0) ;
if ~exist(refMeshFn, 'file') || overwrite 
    if strcmp(QS.piv.imCoords, 'sp_sme')
        refMesh = load(sprintf(QS.fullFileBase.spcutMeshSm, t0), 'spcutMeshSm') ;
        refMesh = refMesh.spcutMeshSm ;
    else
        error('handle this imCoords here')
    end
    umax0 = max(refMesh.u(:, 1)) ;
    vmax0 = max(refMesh.u(:, 2)) ;
    % Get size Lx and Ly for XY space
    im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, t0)) ;
    % for now assume all images are the same size
    disp('for now assuming all images are the same size')
    Lx = size(im0, 1) * ones(length(timePoints), 1) ;
    Ly = size(im0, 2) * ones(length(timePoints), 1) ;
    m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], refMesh.u, doubleCovered, umax0, vmax0) ;
    refMesh.XY = m0XY ;
    
    try 
        refMesh.vrs;
    catch
        refMesh.vrs = QS.xyz2APDV(refMesh.v) ;
    end
    
    % Dump ricci result for this reference timepoint into refMesh
    ricciMeshFn = sprintf(QS.fullFileBase.ricciMesh, maxIter, t0) ;
    if exist(ricciMeshFn, 'file')
        load(ricciMeshFn, 'ricciMesh')
        refMesh.u_ricci = ricciMesh.rectangle.u ;

        % beltrami for this t0
        ricciMuFn = sprintf(QS.fullFileBase.ricciMu, maxIter, t0) ;
        if exist(ricciMuFn, 'file')
            load(ricciMuFn, 'mu_rectangle')
            refMesh.mu = mu_rectangle ;
        else
            disp('could not find beltrami for ricciMesh on disk, computing...')
            refMesh.mu = bc_metric(refMesh.f, refMesh.u_ricci, refMesh.vrs, 3) ;    
        end
    else
        opts = struct() ;
        opts.maxIter = maxIter ;
        [ricciMesh, ricciMu] = QS.generateRicciMeshTimePoint(t0, opts) ;
        refMesh.u_ricci = ricciMesh.rectangle.u ;
        refMesh.mu = ricciMu.rectangle ;
    end
    refMesh.readme = 'mu is computed after Ricci flow generates refMesh.u_ricci' ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute dzeta -- material frame distance between DV hoops
    % (dzeta along the surface from each hoop to the next)
    % The distance from one hoop to another is the
    % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
    refMesh.dzeta = reshape(vecnorm(...
        diff(reshape(refMesh.vrs, [refMesh.nU, refMesh.nV, 3])), 2, 3), ...
        [refMesh.nU-1, refMesh.nV]) ;
    refMesh.dzeta_mean = nanmean(refMesh.dzeta, 2) ;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save refMesh
    save(refMeshFn, 'refMesh')
else
    load(refMeshFn, 'refMesh')
    try
        refMesh.mu ;
    catch
        refMesh.mu = bc_metric(refMesh.f, refMesh.u_ricci, refMesh.vrs, 3) ;
    end
end

%% Save metric images of refMesh
xx = refMesh.u ;
aspectShear = (1 - mean(real(refMesh.mu))) / (1 + mean(real(refMesh.mu))) ;
xx(:, 1) = refMesh.u(:, 1) / max(refMesh.u(:, 1)) * aspectShear ; 
[gcell, bcell] = constructFundamentalForms(refMesh.f, refMesh.vrs, xx) ;
gg = zeros(length(gcell), 4) ;
bb = zeros(length(gcell), 4) ;
for qq = 1:length(gcell)
    gg(qq, 1) = gcell{qq}(1, 1) ;
    gg(qq, 2) = gcell{qq}(1, 2) ;
    gg(qq, 3) = gcell{qq}(2, 1) ;
    gg(qq, 4) = gcell{qq}(2, 2) ;
    bb(qq, 1) = bcell{qq}(1, 1) ;
    bb(qq, 2) = bcell{qq}(1, 2) ;
    bb(qq, 3) = bcell{qq}(2, 1) ;
    bb(qq, 4) = bcell{qq}(2, 2) ;
end
strClims = {'climVariable', 'climUniform'} ;
for pp = 1:2
    close all
    opts = struct() ;
    if pp == 1
        opts.clims = {max(abs(gg(:, 1))) * [-1, 1], ...
            max(abs(gg(:, 2))) * [-1, 1], ...
            max(abs(gg(:, 3))) * [-1, 1], ...
            max(abs(gg(:, 4))) * [-1, 1]} ;
    else
        opts.clim = max(abs(gg(:))) * [-1, 1] ;
    end
    labels = {'$\mathbf{g}_{\zeta\zeta}$', ...
        '$\mathbf{g}_{\zeta\phi}$', ...
        '$\mathbf{g}_{\phi\zeta}$', ...
        '$\mathbf{g}_{\phi\phi}$'} ;
    m2view = refMesh ;
    m2view.v = refMesh.vrs ;
    opts.labels = labels ;
    opts.views = [0, 0] ;
    opts.axisOff = true ;
    [axs, cbs, meshHandles] = ...
        nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
        {gg(:, 1), gg(:, 2), gg(:, 3), gg(:, 4)}, opts) ;
    fn = sprintf(QS.fileName.pathlines.refMesh, t0) ;
    fn = [fn(1:end-4) '_dirichlet_g_fundForm' strClims{pp} '.png'] ;
    saveas(gcf, fn)
    % plot bb
    opts = struct() ;
    if pp == 1
        opts.clims = {max(abs(bb(:, 1))) * [-1, 1], ...
            max(abs(bb(:, 2))) * [-1, 1], ...
            max(abs(bb(:, 3))) * [-1, 1], ...
            max(abs(bb(:, 4))) * [-1, 1]} ;
    else
        opts.clim = max(abs(bb(:))) * [-1, 1] ;
    end
    labels = {'$\mathbf{b}_{\zeta\zeta}$', ...
        '$\mathbf{b}_{\zeta\phi}$', ...
        '$\mathbf{b}_{\phi\zeta}$', ...
        '$\mathbf{b}_{\phi\phi}$'} ;
    opts.labels = labels ;
    opts.views = [0, 0] ;
    opts.axisOff = true ;
    [axs, cbs, meshHandles] = ...
        nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
        {bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4)}, opts) ;
    fn = sprintf(QS.fileName.pathlines.refMesh, t0) ;
    fn = [fn(1:end-4) '_dirichlet_b_fundForm' strClims{pp} '.png'] ;
    saveas(gcf, fn)
end


%% Save metric images of refMesh -- ricci result
xx = refMesh.u_ricci ;
[gcell, bcell] = constructFundamentalForms(refMesh.f, refMesh.vrs, xx) ;
gg = zeros(length(gcell), 4) ;
bb = zeros(length(gcell), 4) ;
for qq = 1:length(gcell)
    gg(qq, 1) = gcell{qq}(1, 1) ;
    gg(qq, 2) = gcell{qq}(1, 2) ;
    gg(qq, 3) = gcell{qq}(2, 1) ;
    gg(qq, 4) = gcell{qq}(2, 2) ;
    bb(qq, 1) = bcell{qq}(1, 1) ;
    bb(qq, 2) = bcell{qq}(1, 2) ;
    bb(qq, 3) = bcell{qq}(2, 1) ;
    bb(qq, 4) = bcell{qq}(2, 2) ;
end
strClims = {'climVariable', 'climUniform'} ;
for pp = 1:2
    close all
    opts = struct() ;
    if pp == 1
        opts.clims = {max(abs(gg(:, 1))) * [-1, 1], ...
            max(abs(gg(:, 2))) * [-1, 1], ...
            max(abs(gg(:, 3))) * [-1, 1], ...
            max(abs(gg(:, 4))) * [-1, 1]} ;
    else
        opts.clim = max(abs(gg(:))) * [-1, 1] ;
    end
    labels = {'$\mathbf{g}_{\zeta\zeta}$', ...
        '$\mathbf{g}_{\zeta\phi}$', ...
        '$\mathbf{g}_{\phi\zeta}$', ...
        '$\mathbf{g}_{\phi\phi}$'} ;
    m2view = refMesh ;
    m2view.v = refMesh.vrs ;
    opts.labels = labels ;
    opts.views = [0, 0] ;
    opts.axisOff = true ;
    [axs, cbs, meshHandles] = ...
        nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
        {gg(:, 1), gg(:, 2), gg(:, 3), gg(:, 4)}, opts) ;
    fn = sprintf(QS.fileName.pathlines.refMesh, t0) ;
    fn = [fn(1:end-4) '_ricci_g_fundForm' strClims{pp} '.png'] ;
    saveas(gcf, fn)
    % plot bb
    opts = struct() ;
    if pp == 1
        opts.clims = {max(abs(bb(:, 1))) * [-1, 1], ...
            max(abs(bb(:, 2))) * [-1, 1], ...
            max(abs(bb(:, 3))) * [-1, 1], ...
            max(abs(bb(:, 4))) * [-1, 1]} ;
    else
        opts.clim = max(abs(bb(:))) * [-1, 1] ;
    end
    labels = {'$\mathbf{b}_{\zeta\zeta}$', ...
        '$\mathbf{b}_{\zeta\phi}$', ...
        '$\mathbf{b}_{\phi\zeta}$', ...
        '$\mathbf{b}_{\phi\phi}$'} ;
    opts.labels = labels ;
    opts.views = [0, 0] ;
    opts.axisOff = true ;
    [axs, cbs, meshHandles] = ...
        nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
        {bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4)}, opts) ;
    fn = sprintf(QS.fileName.pathlines.refMesh, t0) ;
    fn = [fn(1:end-4) '_ricci_b_fundForm' strClims{pp} '.png'] ;
    saveas(gcf, fn)
end

%% For reference -- Dirichlet approach to mu minimization
% refMeshFn = sprintf(QS.fileName.pathlines.refMesh, t0) ;
% if ~exist(refMeshFn, 'file') || overwrite 
%     if strcmp(QS.piv.imCoords, 'sp_sme')
%         refMesh = load(sprintf(QS.fullFileBase.spcutMeshSm, t0), 'spcutMeshSm') ;
%         refMesh = refMesh.spcutMeshSm ;
%     else
%         error('handle this imCoords here')
%     end
%     umax0 = max(refMesh.u(:, 1)) ;
%     vmax0 = max(refMesh.u(:, 2)) ;
%     % Get size Lx and Ly for XY space
%     im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, t0)) ;
%     % for now assume all images are the same size
%     disp('for now assuming all images are the same size')
%     Lx = size(im0, 1) * ones(length(timePoints), 1) ;
%     Ly = size(im0, 2) * ones(length(timePoints), 1) ;
%     m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], refMesh.u, doubleCovered, umax0, vmax0) ;
%     refMesh.XY = m0XY ;
%     uv = refMesh.u ;
%     uv(:, 1) = refMesh.u(:, 1) / max(refMesh.u(:, 1)) ;
%     refMesh.readme = 'mu is computed after rescaling refMesh.u(:, 1) to range from 0 to 1' ;
%     refMesh.mu = bc_metric(refMesh.f, uv, refMesh.v, 3) ;
%     
%     try 
%         refMesh.vrs;
%     catch
%         refMesh.vrs = QS.xyz2APDV(refMesh.v) ;
%     end
%     save(refMeshFn, 'refMesh')
% else
%     load(refMeshFn, 'refMesh')
%     try
%         refMesh.mu ;
%     catch
%         refMesh.mu = bc_metric(refMesh.f, refMesh.u, refMesh.v, 3) ;
%     end
% end
% 
% %% Save metric images of refMesh
% xx = refMesh.u ;
% aspectShear = (1 - mean(real(refMesh.mu))) / (1 + mean(real(refMesh.mu))) ;
% xx(:, 1) = refMesh.u(:, 1) / max(refMesh.u(:, 1)) * aspectShear ; 
% [gcell, bcell] = constructFundamentalForms(refMesh.f, refMesh.vrs, xx) ;
% gg = zeros(length(gcell), 4) ;
% bb = zeros(length(gcell), 4) ;
% for qq = 1:length(gcell)
%     gg(qq, 1) = gcell{qq}(1, 1) ;
%     gg(qq, 2) = gcell{qq}(1, 2) ;
%     gg(qq, 3) = gcell{qq}(2, 1) ;
%     gg(qq, 4) = gcell{qq}(2, 2) ;
%     bb(qq, 1) = bcell{qq}(1, 1) ;
%     bb(qq, 2) = bcell{qq}(1, 2) ;
%     bb(qq, 3) = bcell{qq}(2, 1) ;
%     bb(qq, 4) = bcell{qq}(2, 2) ;
% end
% strClims = {'climVariable', 'climUniform'} ;
% for pp = 1:2
%     opts = struct() ;
%     if pp == 1
%         opts.clims = {max(abs(gg(:, 1))) * [-1, 1], ...
%             max(abs(gg(:, 2))) * [-1, 1], ...
%             max(abs(gg(:, 3))) * [-1, 1], ...
%             max(abs(gg(:, 4))) * [-1, 1]} ;
%     else
%         opts.clim = max(abs(gg(:))) * [-1, 1] ;
%     end
%     labels = {'$\mathbf{g}_{\zeta\zeta}$', ...
%         '$\mathbf{g}_{\zeta\phi}$', ...
%         '$\mathbf{g}_{\phi\zeta}$', ...
%         '$\mathbf{g}_{\phi\phi}$'} ;
%     m2view = refMesh ;
%     m2view.v = refMesh.vrs ;
%     opts.labels = labels ;
%     opts.views = [0, 0] ;
%     opts.axisOff = true ;
%     [axs, cbs, meshHandles] = ...
%         nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
%         {gg(:, 1), gg(:, 2), gg(:, 3), gg(:, 4)}, opts) ;
%     fn = sprintf(QS.fileName.pathlines.refMesh, t0) ;
%     fn = [fn(1:end-4) '_g_fundForm' strClims{pp} '.png'] ;
%     saveas(gcf, fn)
%     % plot bb
%     opts = struct() ;
%     if pp == 1
%         opts.clims = {max(abs(bb(:, 1))) * [-1, 1], ...
%             max(abs(bb(:, 2))) * [-1, 1], ...
%             max(abs(bb(:, 3))) * [-1, 1], ...
%             max(abs(bb(:, 4))) * [-1, 1]} ;
%     else
%         opts.clim = max(abs(bb(:))) * [-1, 1] ;
%     end
%     labels = {'$\mathbf{b}_{\zeta\zeta}$', ...
%         '$\mathbf{b}_{\zeta\phi}$', ...
%         '$\mathbf{b}_{\phi\zeta}$', ...
%         '$\mathbf{b}_{\phi\phi}$'} ;
%     opts.labels = labels ;
%     opts.views = [0, 0] ;
%     opts.axisOff = true ;
%     [axs, cbs, meshHandles] = ...
%         nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
%         {bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4)}, opts) ;
%     fn = sprintf(QS.fileName.pathlines.refMesh, t0) ;
%     fn = [fn(1:end-4) '_b_fundForm' strClims{pp} '.png'] ;
%     saveas(gcf, fn)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advect PIV evaluation coordinates in XY pixel space to make pathlines
% in XY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineXY = sprintf(QS.fileName.pathlines.XY, t0) ;
if ~exist(plineXY, 'file') || overwrite
    disp('Computing piv Pathlines for XY(t0)')
    % Load 'initial' positions for pathlines to intersect at t=t0
    QS.getPIV() 
    
    % NOTE smoothed PIV results will be used if QS.piv.smoothing_sigma > 0
    piv = QS.piv.smoothed ;
    x0 = piv.x{QS.xp.tIdx(t0)} ;
    y0 = piv.y{QS.xp.tIdx(t0)} ;
    
    % Create pathlines emanating from (x0,y0) at t=t0
    [XX, YY] = QS.pullbackPathlines(x0, y0, t0, options) ;
    Lx = QS.piv.Lx ;
    Ly = QS.piv.Ly ;
    
    % Save pathlines of PIV evaluation points XY
    pivPathlines = struct() ;
    pivPathlines.XX = XX ;
    pivPathlines.YY = YY ;
    pivPathlines.Lx = Lx ;
    pivPathlines.Ly = Ly ;
    pivPathlines.t0 = t0 ;
    pivPathlines.tIdx0 =  tIdx0 ;
    
    smoothing_sigma = QS.piv.smoothing_sigma ;
    disp(['Saving ' plineXY])
    save(plineXY, 'pivPathlines', 'smoothing_sigma')
    
    computed_XY = true ;
else
    disp(['PIV XY pathlines already on disk: ' plineXY])
    computed_XY = false ;
end

% Save image of result
plineFig = [plineXY(1:end-4) '.png'] ;
if ~exist(plineFig, 'file') || overwrite
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
    disp(['PIV XY pathline image already on disk: ' plineFig])
end

% Make movie
for tidx = 1:length(timePoints)    
    if mod(tidx, 10) == 0
        disp(['t = ', num2str(timePoints(tidx))])
    end
    fn = fullfile(sprintf(pathlineDir.XY, t0), ...
        [sprintf(QS.fileBase.name, timePoints(tidx)) '.png']) ;
    
    if ~exist(fn, 'file') || overwrite
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
        movieScatterSz = 3 ;
        scatter(XX(tidx, :) / double(Lx(tidx)), ...
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
% Advect Mesh vertex coordinates in XY pixel space in optical flow to make 
% pathlines in XY (pixel space).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plinevXY = sprintf(QS.fileName.pathlines.vXY, t0) ;
if ~exist(plinevXY, 'file') || overwrite 
    % Create pathlines emanating from vertex positions (vX,vY) at t=t0
    % Load pullback mesh vertex positions to get Lx and Ly (will be done
    % later in pullbackPathlines() if not now, so no extra cost to do it
    % now instead). 
    QS.getPIV() 
    Lx = QS.piv.Lx ;
    Ly = QS.piv.Ly ;
    
    if strcmp(QS.piv.imCoords, 'sp_sme')
        mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, t0), 'spcutMeshSm') ;
    else
        error('handle this coord sys here')
    end
    mesh0 = mesh0.spcutMeshSm ;
    umax0 = max(mesh0.u(:, 1)) ;
    vmax0 = max(mesh0.u(:, 2)) ;
    m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], mesh0.u, doubleCovered, ...
        umax0, vmax0) ;
    m0X = m0XY(:, 1) ;
    m0Y = m0XY(:, 2) ;
    
    % Create pathlines emanating from vertex positions (vX,vY) at t=t0
    [vX, vY] = QS.pullbackPathlines(m0X, m0Y, t0, options) ;
    vX = reshape(vX, [length(timePoints), QS.nU, QS.nV]) ;
    vY = reshape(vY, [length(timePoints), QS.nU, QS.nV]) ;
    
    % Save pathlines of mesh locations u in XY coords, vXY
    vertexPathlines = struct() ;
    vertexPathlines.vX = vX ;
    vertexPathlines.vY = vY ;
    vertexPathlines.t0 = t0 ;
    vertexPathlines.tIdx0 = tIdx0 ;
    vertexPathlines.Lx = Lx ;
    vertexPathlines.Ly = Ly ;
    
    smoothing_sigma = QS.piv.smoothing_sigma ;
    save(plinevXY, 'vertexPathlines', 'smoothing_sigma')
    
    computed_vXY = true ;
else
    disp(['Mesh vertex XY pathlines already on disk: ' plinevXY])
    computed_vXY = false ;
end

% Save image of result
plineFig = [plinevXY(1:end-4) '.png'] ;
if ~exist([plinevXY(1:end-4) '.png'], 'file') || overwrite
    disp(['Saving mesh vertex XY pathline image to disk: ' plineFig])
    if ~computed_vXY 
        load(QS.fileName.pathlines.vXY, 'vertexPathlines')
        vX = vertexPathlines.vX ;
        vY = vertexPathlines.vY ;
        t0 = vertexPathlines.t0 ;
        tIdx0 = vertexPathlines.tIdx ;
        Lx = vertexPathlines.Lx ;
        Ly = vertexPathlines.Ly ;
    end
    
    % Plot the pathlines -- NOTE that we flip YY coordinates (MATLAB)
    qsubX = round(size(vX, 2) / nY2plot) ;
    qsubY = round(size(vX, 3) / nY2plot * QS.a_fixed) ;
    colors = parula(length(timePoints)) ;
    clf
    for qq = fliplr(1:length(timePoints))
        xx = vX(qq, 1:qsubX:end, 1:qsubY:end) ;
        yy = vY(qq, 1:qsubX:end, 1:qsubY:end) ;
        scatter(xx(:) / double(Lx(qq)), ...
            (1 - (yy(:) / double(Ly(qq))) - Yoffset) / (1 - 2*Yoffset), ...
            scatterSz, ...
            'markerfacecolor', colors(qq, :), ...
            'markeredgecolor', 'none', 'MarkerFaceAlpha', 0.3) ;
        hold on;
        pause(0.1) ;
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
    disp(['Mesh vertex XY pathline image already on disk: ' plineFig])
end

% Make movie
for tidx = 1:length(timePoints)    
    if mod(tidx, 10) == 0
        disp(['t = ', num2str(timePoints(tidx))])
    end
    fn = fullfile(sprintf(pathlineDir.vXY, t0), ...
        [sprintf(QS.fileBase.name, timePoints(tidx)) '.png']) ;
    
    if ~exist(fn, 'file') || overwrite
        % Load pathlines if not already in RAM
        if ~computed_vXY 
            load(plineXY, 'pivPathlines')
            vX = pivPathlines.vX ;
            vY = pivPathlines.vY ;
            t0 = pivPathlines.t0 ;
            tIdx0 = pivPathlines.tIdx0 ;
            Lx = pivPathlines.Lx ;
            Ly = pivPathlines.Ly ;
            computed_vXY = true ;
        end
    
        clf ;
        movieScatterSz = 3 ;
        scatter(vX(tidx, :) / double(Lx(tidx)), ...
            (1 - (vY(tidx, :) / double(Ly(tidx))) - Yoffset) / (1 - 2*Yoffset), ...
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
%% Advect Mesh face barycenters in XY pixel space to make pathlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineFXY = sprintf(QS.fileName.pathlines.fXY, t0) ;
if ~exist(plineFXY, 'file') || overwrite
    % Create pathlines emanating from barycenters (bcX,bcY) at t=t0
    % Load pullback mesh vertex positions to get Lx and Ly (will be done
    % later in pullbackPathlines() if not now, so no extra cost to do it
    % now instead).
    QS.loadPIV() 
    Lx = QS.piv.Lx ;
    Ly = QS.piv.Ly ;
    
    if strcmp(QS.piv.imCoords, 'sp_sme')
        mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, t0), 'spcutMeshSm') ;
    else
        error('handle coord sys here')
    end
    mesh0 = mesh0.spcutMeshSm ;
    umax0 = max(mesh0.u(:, 1)) ;
    vmax0 = max(mesh0.u(:, 2)) ;
    m0XY = QS.uv2XY([Lx(tIdx0), Ly(tIdx0)], mesh0.u, doubleCovered, umax0, vmax0) ;
    
    % Obtain barycenters of faces
    f0XY = barycenter(m0XY, mesh0.f) ;  % this is in gptoolbox
    f0X = f0XY(:, 1) ;
    f0Y = f0XY(:, 2) ;
    
    % Create pathlines emanating from barycenters (bcX,bcY) at t=t0
    [fX, fY] = QS.pullbackPathlines(f0X, f0Y, t0, options) ;
    
    % Save pathlines of mesh face barycenters f in XY coords, fXY
    facePathlines = struct() ;
    facePathlines.fX = fX ;
    facePathlines.fY = fY ;
    facePathlines.t0 = t0 ;
    facePathlines.tIdx0 = tIdx0 ;
    facePathlines.Lx = Lx ;
    facePathlines.Ly = Ly ;
    save(plineFXY, 'facePathlines')
    
    computed_fXY = true ;
else
    computed_fXY = false ;
end

% Save image of result
plineFig = [plineFXY(1:end-4) '.png'] ;
if ~exist(plineFig, 'file') || overwrite
    if ~computed_fXY 
        load(QS.fileName.pathlines.fXY, 'facePathlines')
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
%% Push forward from pullback XY pathlines to embedding coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plineXYZ = sprintf(QS.fileName.pathlines.XYZ, t0) ;
if ~exist(plineXYZ, 'file') || overwrite
    disp('Compute 3d pathlines in embedding space for PIV eval coords')
    % Build 3d pathlines in embedding space for all vertex locations
    % Load 2d pathlines in pullback space if not already in RAM
    if ~computed_XY 
        load(plineXY, 'pivPathlines', 'smoothing_sigma')
        XX = pivPathlines.XX ;
        YY = pivPathlines.YY ;
        t0 = pivPathlines.t0 ;
        tIdx0 = pivPathlines.tIdx0 ;
        Lx = pivPathlines.Lx ;
        Ly = pivPathlines.Ly ;
        computed_XY = true ;
        try
            assert(smoothing_sigma == QS.piv.smoothing_sigma)
        catch
            error('Smoothing_sigma of QS.piv does not match that stored on file in pivPathlines')
        end
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
        mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), ...
                             'spcutMeshSm') ;
        mesh0 = mesh0.spcutMeshSm ;
        umax0 = max(mesh0.u(:, 1)) ;
        vmax0 = max(mesh0.u(:, 2)) ;
        tileCount = [2 2] ;
        [tm0f, tm0v2d, tm0v3d, ~] = tileAnnularCutMesh(mesh0, tileCount);
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
    
    save(plineXYZ, 'piv3dPathlines', 'smoothing_sigma')
end                             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Push forward vertex pathlines from pullback to embedding coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plinev3d = sprintf(QS.fileName.pathlines.v3d, t0) ;
if ~exist(plinev3d, 'file') || overwrite 
    disp('Compute 3d pathlines in embedding space for Lagrangian/advected spcutMesh vertex coords')
    % Build 3d pathlines in embedding space for all vertex locations
    % Load 2d pathlines in pullback space if not already in RAM
    if ~computed_vXY 
        load(sprintf(QS.fileName.pathlines.vXY, t0), 'vertexPathlines', 'smoothing_sigma')
        vX = vertexPathlines.vX ;
        vY = vertexPathlines.vY ;
        t0 = vertexPathlines.t0 ;
        tIdx0 = vertexPathlines.tIdx0 ;
        Lx = vertexPathlines.Lx ;
        Ly = vertexPathlines.Ly ;
        computed_vXY = true ; % we have loaded the result, so we've already computed it
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
        
        % Load this timepoint's cutMesh
        mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), ...
                             'spcutMeshSm') ;
        mesh0 = mesh0.spcutMeshSm ;
        umax0 = max(mesh0.u(:, 1)) ;
        vmax0 = max(mesh0.u(:, 2)) ;
        tileCount = [2 2] ;
        [tm0f, tm0v2d, tm0v3d, ~] = tileAnnularCutMesh(mesh0, tileCount);
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
    save(plinev3d, 'v3dPathlines', 'smoothing_sigma')
    computed_v3d = true ;
else
    computed_v3d = false ;
end

%% Plot pathlines in 3d
% Save image of result
[~,~,~,xyzlims] = QS.getXYZLims() ;
nTimePoints = 45 ;
first = true ;
slab = 20 ;
black_figs = true ;

if ~computed_v3d 
    load(sprintf(QS.fileName.pathlines.v3d, t0), 'v3dPathlines')
    vX3rs = v3dPathlines.vXrs ;
    vY3rs = v3dPathlines.vYrs ;
    vZ3rs = v3dPathlines.vZrs ;
    t0 = v3dPathlines.t0 ;
    tIdx0 = v3dPathlines.tIdx0 ;
    computed_v3d = false ;
end

% Get subsampling
qsubX = round(size(vX3rs, 2) / nY2plot) - 1 ;
qsub = 1 ;

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
        daspect([1 QS.a_fixed 1])
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
        disp(['Saving ' plineFig])
        export_fig(plineFig, '-nocrop', '-r150')
        close all
    else
        disp(['Mesh advected vertex XYZ pathline image already on disk: ' plineFig])
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RADII -- Push forward vertex pathlines from pullback to embedding coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pliner3d = sprintf(QS.fileName.pathlines.radius, t0) ;
if ~exist(pliner3d, 'file') || overwrite || overwriteImages
    disp('Computing radii along pathlines in embedding space for Lagrangian/advected uvprimecutMesh vertex coords')
    % Build 3d pathlines in embedding space for all vertex locations
    % Load 2d pathlines in pullback space if not already in RAM
    if ~computed_vXY 
        load(sprintf(QS.fileName.pathlines.vXY, t0), 'vertexPathlines')
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
        
        % Load this timepoint's smoothed cutMesh
        if strcmpi(QS.piv.imCoords, 'sp_sme')
            mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), ...
                                 'spcutMeshSm') ;
        else
            error('handle these piv coords here')
        end
        mesh0 = mesh0.spcutMeshSm ;
        
        if ~isfield(mesh0, 'radius_um')
            error('no radius stored with mesh -- can add in post via code here')
            % 
            % for tidx = 1:length(QS.xp.fileMeta.timePoints)
            %     tp = QS.xp.fileMeta.timePoints(tidx) ;
            %     % Load this timepoint spcutMeshSm
            %     QS.setTime(tp)
            %     spcutMeshSm = QS.getCurrentSPCutMeshSm() ;
            %     spcutMeshSmRS = QS.getCurrentSPCutMeshSmRS() ;
            %     spcutMeshSmRSC = QS.getCurrentSPCutMeshSmRSC() ;
            % 
            %     % Make avgpts in pixel space (not RS)
            %     fprintf('Resampling uvgrid3d curves in pix...\n')
            %     nU = spcutMeshSm.nU ;
            %     nV = spcutMeshSm.nV ;
            %     curves3d_pix = reshape(spcutMeshSm.v, [nU, nV, 3]) ;
            %     c3d_dsv_pix = zeros(size(curves3d_pix)) ;  % in units of pix
            %     avgpts_pix = zeros(nU, 3) ;
            %     radius_pix = zeros(nU, nV) ;
            %     for i=1:nU
            %         % Note: no need to add the first point to the curve
            %         % since the endpoints already match exactly in 3d and
            %         % curvspace gives a curve with points on either
            %         % endpoint (corresponding to the same 3d location).
            %         c3d_dsv_pix(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d_pix(i, :, :)), nV, true) ;
            %         if vecnorm(squeeze(c3d_dsv_pix(i, 1, :)) - squeeze(c3d_dsv_pix(i, end, :))) > 1e-7
            %             error('endpoints do not join! Exiting')
            %         end
            %         % Drop the final endpoint in the mean pt determination
            %         avgpts_pix(i, :) = mean(squeeze(c3d_dsv_pix(i, 1:end-1, :)), 1) ; 
            %         radius_pix(i, :) = vecnorm(squeeze(curves3d_pix(i, :, :)) - avgpts_pix(i, :), 2, 2) ;
            %     end
            % 
            %     % uvpcutMesh.raw.avgpts_pix = avgpts_pix ;
            %     % uvpcutMesh.raw.radius_pix = radius_pix ;
            %     spcutMeshSm.avgpts_um = QS.xyz2APDV(avgpts_pix) ;
            %     spcutMeshSm.radius_um = radius_pix * QS.APDV.resolution ;
            % 
            %     % Add to RS and RSC versions
            %     spcutMeshSmRS.avgpts_um = spcutMeshSm.avgpts_um ;
            %     spcutMeshSmRS.radius_um = spcutMeshSm.radius_um ;
            %     % RSC
            %     spcutMeshSmRSC.avgpts_um = spcutMeshSm.avgpts_um ;
            %     spcutMeshSmRSC.radius_um = spcutMeshSm.radius_um(:, 1:end-1) ;
            % 
            %     % Save spcutMeshSm / RS / RSC
            %     save(sprintf(QS.fullFileBase.spcutMeshSm, QS.currentTime), ...
            %         'spcutMeshSm')
            %     save(sprintf(QS.fullFileBase.spcutMeshSmRS, QS.currentTime), ...
            %         'spcutMeshSmRS')
            %     save(sprintf(QS.fullFileBase.spcutMeshSmRSC, QS.currentTime), ...
            %         'spcutMeshSmRSC')
            % 
            %     if mod(tidx, 10) < 3
            %         clf; set(gcf, 'visible', 'on')
            %         trisurf(triangulation(spcutMeshSmRSC.f, spcutMeshSmRSC.v), 'edgecolor', 'none')
            %         axis equal 
            %         title(num2str(tidx))
            %         pause(1)
            %         cla
            %         trisurf(triangulation(spcutMeshSmRS.f, spcutMeshSmRS.v), 'edgecolor', 'k')
            %         axis equal
            %         title([num2str(tidx) ' open'])
            %         pause(0.0001)
            %     end
            % end
        end
        
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
kymoDir = sprintf(QS.dir.pathlines.kymographs, t0) ;
if ~exist(kymoDir, 'dir')
    mkdir(kymoDir) 
end
plinerAP = sprintf(QS.fileName.pathlines.kymographs.radius, t0) ;
if ~exist(plinerAP, 'file') || overwrite 
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
else
    load(plinerAP, 'radius_apM', 'radius_dM', 'radius_vM', ...
        'radius_lM', 'radius_rM')
end

% Image filenames
fns = {[plinerAP(1:end-4) '_ap.png'], ...
    [plinerAP(1:end-4) '_d.png'], ...
    [plinerAP(1:end-4) '_v.png'], ...
    [plinerAP(1:end-4) '_l.png'], ...
    [plinerAP(1:end-4) '_r.png']} ;
fns_exist = exist(fns{1}, 'file') ;
if ~fns_exist || overwrite || overwriteImages
    % Get nU and nV
    if ~computed_rad
        load(pliner3d, 'vRadiusPathlines')
        computed_rad = true ;
    end
    timePoints = QS.xp.fileMeta.timePoints;
    nU = size(vRadiusPathlines.radii, 2) ;
    nV = size(vRadiusPathlines.radii, 3) ;
    
    % Plot each kymo
    clf
    kymos = {radius_apM, radius_dM, radius_vM, radius_lM, radius_rM} ;
    titles = {'$(u'',v'')$ pathline radii', ...
        '$(u'',v'')$ pathline radii, dorsal side', ...
        '$(u'',v'')$ pathline radii, ventral side', ...
        '$(u'',v'')$ pathline radii, left side', ...
        '$(u'',v'')$ pathline radii, right side'} ;
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

%% Plot radius in 3d
close all
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
        cb.Position(4) = 0.6 * cb.Position(4) ;
        cb.Position(2) = cb.Position(2) + 0.15 * cbheight ;
        ax = gca ;
        ax.Position(1) = ax.Position(1) - 0.05 ; 
        ylabel(cb, ['radius [' QS.spaceUnits ']' ], 'Interpreter', 'Latex')
        
        title(['$t=$' sprintf('%03d', (tp-t0)*QS.timeInterval) ' ' QS.timeUnits ], ...
            'interpreter', 'latex')
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

