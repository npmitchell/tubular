function helmholtzHodge(tubi, options) 
%helmholtzHodge(tubi, options) 
%   - Compute DEC divergence and DEC "curl" on 2d evolving surface in 3d.
%   - Compute laplacian in covariant frame 
%   - Compute scalar fields representing rotational and harmonic components
%   - Store all on disk. 
%
% Parameters
% ----------
% tubi : TubULAR class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   averagingStyle : str ('Lagrangian', 'simple')
%       style in which velocities were averaged over time
%   alphaVal : float
%       the opacity of the heatmap to overlay
%   invertImage : bool
%       invert the data pullback under the velocity map
%   clipDiv : float (default=5.0)
%       max allowed divergence measurement  
%   clipRot : float (default=0.5)
%       max allowed vorticity measurement
%   computeLaplacian : bool (default = false)
%       compute the laplacian of the velocity field
%
% NPMitchell 2020/2021

%% Default options
overwrite = false ;
overwrite_images = false ;
qsubU = 5 ; 
qsubV = 10 ;
niter_smoothing = [1000, 1000, 1000] ;
plot_dec_pullback = true ;
plot_lap_pullback = true ;
plot_dec_texturepatch = false ;
preview = false ;
pivimCoords = tubi.piv.imCoords ;
lambda_smooth = tubi.smoothing.lambda ;
lambda_mesh = tubi.smoothing.lambda_mesh ;
nmodes = tubi.smoothing.nmodes ;
zwidth = tubi.smoothing.zwidth ;
samplingResolution = '1x' ;
averagingStyle = 'Lagrangian' ;
clipDiv = 5.0 ;                     % max allowed divergence measurement  
clipRot = 0.5 ;                     % max allowed vorticity measurement
sscaleDiv = 0.5 ;                   % climit (color limit) for divergence
sscaleRot = 0.15 ;                  % climit (color limit) for vorticity
sscaleLap = 0.05 ;                  % climit (color limit) for Laplacian of vel
computeLaplacian = false ;

%% Unpack options
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end
if isfield(options, 'qsubU')
    qsubU = options.qsubU ;
end
if isfield(options, 'qsubV')
    qsubV = options.qsubV ;
end
if isfield(options, 'clipDiv')
    clipDiv = options.clipDiv ;
end
if isfield(options, 'clipRot')
    clipRot = options.clipRot ;
end
if isfield(options, 'sscaleDiv')
    sscaleDiv = options.sscaleDiv ;
end
if isfield(options, 'sscaleRot')
    sscaleRot = options.sscaleRot ;
end
if isfield(options, niter_smoothing)
    niter_smoothing = options.niter_smoothing ;
end
if isfield(options, plot_dec_pullback)
    plot_dec_pullback = options.plot_dec_pullback ;
end
if isfield(options, plot_dec_texturepatch)
    plot_dec_texturepatch = options.plot_dec_texturepatch ;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
end
if isfield(options, 'lambda')
    lambda_smooth = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'computeLaplacian')
    computeLaplacian = options.computeLaplacian ;
end

% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack tubi
timePoints = tubi.xp.fileMeta.timePoints ;
fname = tubi.fileBase.name ;
[rot, trans] = tubi.getRotTrans() ;
resolution = tubi.APDV.resolution ;
if doubleResolution
    nU = 2 * tubi.nU - 1 ;
    nV = 2 * tubi.nV - 1 ;
    piv3dFileBase = tubi.fullFileBase.piv3d2x ;
else
    nU = tubi.nU ;
    nV = tubi.nV ;
    piv3dFileBase = tubi.fullFileBase.piv3d ;
end
if strcmp(averagingStyle, 'Lagrangian')
    if doubleResolution
        tubi.getVelocityAverage2x('vf', 'v2dum') ;
        vf = tubi.velocityAverage2x.vf ;
        v2dum = tubi.velocityAverage2x.v2dum ;
        decDirRoot = tubi.dir.piv.avgDEC2x ;
        decFnBase = tubi.fullFileBase.decAvg2x ;
    else
        tubi.getVelocityAverage('vf', 'v2dum') ;
        vf = tubi.velocityAverage.vf ;
        v2dum = tubi.velocityAverage.v2dum ;
        decDirRoot = tubi.dir.piv.avgDEC ;
        decFnBase = tubi.fullFileBase.decAvg ;    
    end
elseif strcmp(averagingStyle, 'simple') 
    if doubleResolution
        tubi.getVelocitySimpleAverage2x('vf', 'v2dum') ;
        vf = tubi.velocitySimpleAverage2x.vf ;
        v2dum = tubi.velocitySimpleAverage2x.v2dum ;
        decDirRoot = tubi.dir.pivSimAvgDEC2x ;
        nU = 2 * tubi.nU - 1 ;
        nV = 2 * tubi.nV - 1 ;
        piv3dFileBase = tubi.fullFileBase.piv3d2x ;
        decFnBase = tubi.fullFileBase.decSimAvg2x ;
    else
        tubi.getVelocitySimpleAverage('vf', 'v2dum') ;
        vf = tubi.velocitySimpleAverage.vf ;
        v2dum = tubi.velocitySimpleAverage.v2dum ;
        decDirRoot = tubi.dir.pivSimAvgDEC ;
        nU = tubi.nU ;
        nV = tubi.nV ;
        piv3dFileBase = tubi.fullFileBase.piv3d ;
        decFnBase = tubi.fullFileBase.decSimAvg ;
    end
elseif strcmp(averagingStyle, 'none') 
    if doubleResolution
        tubi.getVelocityRaw('vf', 'v2dum') ;
        vf = tubi.velocityRaw2x.vf ;
        v2dum = tubi.velocityRaw2x.v2dum ;
        decDirRoot = tubi.dir.pivRawDEC2x ;
        nU = 2 * tubi.nU - 1 ;
        nV = 2 * tubi.nV - 1 ;
        piv3dFileBase = tubi.fullFileBase.piv3d2x ;
        decFnBase = tubi.fullFileBase.dec2x ;
    else
        tubi.getVelocityRaw('vf', 'v2dum') ;
        vf = tubi.velocityRaw.vf ;
        v2dum = tubi.velocityRaw.v2dum ;
        decDirRoot = tubi.dir.pivRawDEC ;
        nU = tubi.nU ;
        nV = tubi.nV ;
        piv3dFileBase = tubi.fullFileBase.piv3d ;
        decFnBase = tubi.fullFileBase.dec ;
    end
else
    error('averagingStyle not recognized. Use Lagrangian or simple')
end
 
t0 = tubi.t0set() ;
[~, ~, ~, xyzlim] = tubi.getXYZLims() ;

%% define the output dirs
decDir = decDirRoot.data ;
dDir2d = decDirRoot.div2d ;
dDir3d = decDirRoot.div3d ;
dDir3dt = decDirRoot.div3dTexture ;
rDir2d = decDirRoot.rot2d ;
rDir3d = decDirRoot.rot3d ;
rDir3dt = decDirRoot.rot3dTexture ;
hDir2d = decDirRoot.harm2d ;
hDir3d = decDirRoot.harm3d ;
LDir2d = decDirRoot.lap2d ;
LDir3d = decDirRoot.lap3d ;
LDir3dt = decDirRoot.lap3dTexture ;
% create the output dirs
dirs2do = {decDir, dDir2d, dDir3d, dDir3dt, ...
    rDir2d, rDir3d, rDir3dt, ...
    LDir2d, LDir3d, LDir3dt, ...
    hDir2d, hDir3d} ;
for i = 1:length(dirs2do)
    ensureDir(dirs2do{i}) ;
end

%% Perform Decomposition
% Consider each timepoint and plot the div and curl
tidxAll = 1:length(timePoints(1:end-1)) ;
tidx2do = 1:50:length(timePoints(1:end-1)) ;
tidx10 = 1:10:length(timePoints(1:end-1)) ;
% add next sparsity level (every 10)
tidxNext = setdiff(tidx10, tidx2do) ;
tidx2do = [tidx2do, tidxNext] ;
% add final sparsity level (every frame)
tidxOther = setdiff(tidxAll, tidx2do) ;
tidx2do = [tidx2do, tidxOther] ;
for tidx = tidx2do
    tp = timePoints(tidx) ;
    disp(['t = ', num2str(tp)])
    
    % DEC data filename to save
    decDataFn = sprintfm(decFnBase, tp) ;
    
    % Prepare filenames
    % pullback view of divergence
    div2dfn = fullfile(dDir2d, [sprintf(fname, tp) '_div2d.png']) ; 
    % lateral view of surface plot
    div3dfn_lateral = fullfile(dDir3d, ['lateral_' sprintf(fname, tp) '_div3d.png']) ;
    % ventral view of surface plot
    div3dfn_ventral = fullfile(dDir3d, ['ventral_' sprintf(fname, tp) '_div3d.png']) ;
    % texture patch overlay
    div3dtfn = fullfile(dDir3dt, [sprintf(fname, tp) '_divt3d.png']) ;
    
    lap2dfn = fullfile(LDir2d, [sprintf(fname, tp) '_lap2d.png']) ; 
    lap3dfn_lateral = fullfile(LDir3d, ['lateral_' sprintf(fname, tp) '_lap3d.png']) ;
    lap3dfn_ventral = fullfile(LDir3d, ['ventral_' sprintf(fname, tp) '_lap3d.png']) ;
    lap3dtfn = fullfile(LDir3dt, [sprintf(fname, tp) '_lapt3d.png']) ;
    
    rot2dfn = fullfile(rDir2d, [sprintf(fname, tp) '_rot2d.png']) ;
    rot3dfn_lateral = fullfile(rDir3d, ['lateral_' sprintf(fname, tp) '_rot3d.png']) ;
    rot3dfn_ventral = fullfile(rDir3d, ['ventral_' sprintf(fname, tp) '_rot3d.png']) ;
    rot3dtfn = fullfile(rDir3dt, [sprintf(fname, tp) '_rot3d.png']) ;
    harm2dfn = fullfile(hDir2d, [sprintf(fname, tp) '_harm2d.png']) ; 
    harm3dfn = fullfile(hDir3d, [sprintf(fname, tp) '_harm3d.png']) ;

    % Discern whether to plot results in 2d and 3d or in texturepatch
    plot_dec_pullback_tidx = plot_dec_pullback && ...
        (~exist(div2dfn, 'file') || ...
        ~exist(div3dfn_lateral, 'file') || ...
        ~exist(div3dfn_ventral, 'file') || overwrite || overwrite_images) ;
    plot_lap_pullback_tidx = plot_lap_pullback && ...
        (~exist(lap2dfn, 'file') || ...
        ~exist(lap3dfn_lateral, 'file') || ...
        ~exist(lap3dfn_ventral, 'file') || overwrite || overwrite_images) ;
    plot_dec_texturepatch_tidx = plot_dec_texturepatch && ...
        (overwrite || ~exist(div2dfn, 'file') || ...
        ~exist(div3dfn_lateral, 'file') || ...
        ~exist(div3dfn_ventral, 'file') || overwrite || overwrite_images) ;
    
    if ~exist(decDataFn, 'file') || overwrite || true
        if ~exist(decDataFn, 'file')
            disp(['Creating DEC: ' decDataFn])
        else
            disp(['Overwriting DEC: ', decDataFn])
        end
        
        % Load piv3d
        piv3d = load(sprintfm(piv3dFileBase, tp)) ;
        piv3d = piv3d.piv3dstruct ;
    
        % Obtain smoothed velocities on all faces
        vfsm = squeeze(vf(tidx, :, :)) ;
        v2dsmum_ii = squeeze(v2dum(tidx, :, :)) ;
        
        % Use current time's tiled smoothed mesh
        % Note: vfsmM is in um/min rs
        FF = piv3d.m0f ;   % #facesx3 float: mesh connectivity list
        V2D = piv3d.m0XY ; % Px2 float: 2d mesh vertices in pullback image pixel space
        v3drs = tubi.xyz2APDV(piv3d.m0v3d) ;
        
        % Check it
        % figure;
        % set(gcf, 'visible', 'on')
        % plot(v3drs)
        % pause(1)
        % if tp > timePoints(1)
        %     plot(v3drs - oldv)
        %     disp('here')
        % end
        % oldv = v3drs ;
        % close all
        
        % Grab cutMesh from piv3d. Could grab from disk instead...
        cutM.f = FF ;
        cutM.u = V2D ;
        cutM.v = v3drs ;
        cutM.nU = nU ;
        cutM.nV = nV ;
        
        %% Compute divs and rots. Note that smoothing occurs inside func
        disp('Decomposing flow into div/rot...')
        Options = struct() ;
        Options.lambda = lambda_smooth ;   
        Options.lambda_mesh = lambda_mesh ;     
        Options.nSpectralFilterModes = nmodes ;
        Options.spectralFilterWidth = zwidth ;
        Options.outdir = tubi.dir.piv.avgDEC.data ;
        Options.do_calibration = (tidx == tidx2do(1)) ;
        Options.computeLaplacian = computeLaplacian ;
        [divs, rots, harms, lapvs, glueMesh] = ...
            helmHodgeDECRectGridPullback(cutM, vfsm, Options,...
            'niterSmoothing', niter_smoothing, ...
            'clipDiv', [-clipDiv, clipDiv], ...
            'clipRot', [-clipRot, clipRot], ...
            'preview', preview, 'method', 'both') ;
        
        [gV2F, ~] = meshAveragingOperators(glueMesh.f, glueMesh.v) ;
        if computeLaplacian
            [lapvn, lapvt, lapv2d, ~, ~, ~, dilation] = ...
                resolveTangentNormalVelocities(cutM.f, ...
                cutM.v, gV2F * lapvs.lapv(1:nU*(nV-1), :), ...
                1:length(cutM.f), V2D) ;

            % tmp = vecnorm(lapv0t, 2, 2) ;
            % plot(lapv0n, tmp, '.')

            lapv2dsc = lapv2d ./ dilation(:) * resolution ;
            lapvs.lapv2d = lapv2d ;
            lapvs.lapv2dsc = lapv2dsc ;
            lapvs.lapvn = lapvn ;
            lapvs.lapvt = lapvt ;
        else
            lapvs = [] ;
        end
        
        
        %% save divs, rots, and harms as structs in .mat file
        disp(['Saving DEC t=' num2str(tp) ': ' decDataFn])
        save(decDataFn, 'divs', 'rots', 'harms', 'lapvs', ...
            'lambda_smooth', 'lambda_mesh', 'nmodes', 'zwidth')
    elseif plot_dec_pullback_tidx || plot_lap_pullback_tidx || ...
            plot_dec_texturepatch_tidx
        disp(['Loading DEC t=' num2str(tp) ': ' decDataFn])
        load(decDataFn, 'divs', 'rots', 'lapvs')
    
        % Load piv3d
        piv3d = load(sprintfm(piv3dFileBase, tp)) ;
        piv3d = piv3d.piv3dstruct ;

        % Obtain smoothed velocities on all faces
        vfsm = squeeze(vf(tidx, :, :)) ;
        v2dsmum_ii = squeeze(v2dum(tidx, :, :)) ;

        % Use current time's tiled smoothed mesh
        % Note: vfsmM is in um/min rs
        FF = piv3d.m0f ;   % #facesx3 float: mesh connectivity list
        V2D = piv3d.m0XY ; % Px2 float: 2d mesh vertices in pullback image pixel space
        v3drs = tubi.xyz2APDV(piv3d.m0v3d) ;

        % Grab cutMesh from piv3d. Could grab from disk instead...
        cutM.f = FF ;
        cutM.u = V2D ;
        cutM.v = v3drs ;
        cutM.nU = nU ;
        cutM.nV = nV ;
    end

    %% Plot results
    if plot_dec_pullback_tidx
        disp('Plotting div/rot...')
        close all
        set(gcf, 'visible', 'off')
        if strcmp(pivimCoords, 'sp_sme')
            im = imread(sprintfm(tubi.fullFileBase.im_sp_sme, tp)) ;
            ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
        else
            error(['Have not coded for this pivimCoords option. Do so here: ' pivimCoords])
        end
        % convert im to rgb for no cmap change, if not already RGB
        if size(im, 3) ~= 3
            im = cat(3, im, im, im) ;
        end
        addTitleStr = [': $t=$', num2str((tp - t0)*tubi.timeInterval), ...
                       ' ', tubi.timeUnits] ;
        Options = struct() ;
        Options.addTitleStr = addTitleStr ;
        Options.div2dfn = div2dfn ;
        Options.div3dfn_lateral = div3dfn_lateral ;
        Options.div3dfn_ventral = div3dfn_ventral ;
        Options.rot2dfn = rot2dfn ;
        Options.rot3dfn_lateral = rot3dfn_lateral ;
        Options.rot3dfn_ventral = rot3dfn_ventral ;
        Options.qsubU = qsubU ; 
        Options.qsubV = qsubV ;
        Options.sscaleDiv = sscaleDiv ;
        Options.sscaleRot = sscaleRot ;
        Options.qscaleDiv3d = 0 ;
        Options.qscaleRot3d = 0 ;
        Options.qscaleDiv2d = 0 ;
        Options.qscaleRot2d = 0 ;
        Options.xyzlim = xyzlim ;
        opts2d.xlim = [0, size(im, 1)] ;
        opts2d.ylim = [0.25 * size(im, 2), 0.75 * size(im, 2) ] ;
        xy = {piv3d.x0, piv3d.y0} ;
        plotHelmHodgeDECPullback(im, cutM, vfsm, xy, v2dsmum_ii, ...
            divs, rots, Options, opts2d) ;
    end
    
    % Plot laplacian of tangential vector field (pushed onto vertices)
    if plot_lap_pullback_tidx && computeLaplacian
        disp('Plotting laplacian(v)')
        close all
        set(gcf, 'visible', 'off')
        if strcmp(pivimCoords, 'sp_sme')
            im = imread(sprintfm(tubi.fullFileBase.im_sp_sme, tp)) ;
            ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
        else
            error(['Have not coded for this pivimCoords option. Do so here: ' pivimCoords])
        end
        im = cat(3, im, im, im) ;  % convert to rgb for no cmap change
        addTitleStr = [': $t=$', num2str((tp - t0)*tubi.timeInterval), ...
                       ' ', tubi.timeUnits] ;
        Options = struct() ;
        Options.addTitleStr = addTitleStr ;
        Options.lapv2dfn = lap2dfn ;
        Options.lapv3dfn_lateral = lap3dfn_lateral ;
        Options.lapv3dfn_ventral = lap3dfn_ventral ;
        Options.qsubU = qsubU ; 
        Options.qsubV = qsubV ;
        Options.sscale = sscaleLap ;
        Options.qscale3d = 0 ;
        Options.qscale2d = 0 ;
        Options.xyzlim = xyzlim ;
        opts2d.xlim = [0, size(im, 1)] ;
        opts2d.ylim = [0.25 * size(im, 2), 0.75 * size(im, 2) ] ;
        xy = {piv3d.x0, piv3d.y0} ;
        plotDECLaplacianVPullback(im, cutM, vfsm, xy, lapvs, ...
            Options, opts2d) ;
    end
    
    if plot_dec_texturepatch_tidx
        % load current timepoint
        % (3D data for coloring mesh pullback)
        tubi.setTime(tp) ;
        tubi.getCurrentData() ;     
        IV = tubi.currentData.IV ;
        IV = imcomplement(IV) ; % does this work? debug 2020
        % IV = max(IV(:)) - IV ; % used to do this.

        addTitleStr = ['$t=$', num2str(tp - t0)] ;
        Options.addTitleStr = addTitleStr ;
        Options.div3dfn = div3dtfn ;
        Options.rot3dfn = rot3dtfn ;
        Options.sscaleDiv = 0.5 ;
        Options.sscaleRot = 0.2 ;
        Options.xyzlim = xyzlim ;

        % Load the cutmesh vertices and normals
        cutM.v = piv3d.m0v3d ;
        cutM.v3drs = ((rot * cutM.v')' + trans) * resolution ;
        plotHelmHodgeDECTexture3d(IV, cutM, divs, rots, rot, trans, Options) ;
    end
end
