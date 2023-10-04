function outStruct = measurePathlineVelocities(tubi, options)
%measurePullbackStreamlines(tubi, options)
%   Use pathlines of optical flow in pullback space to query velocities
%   and average along pathlines.
%
% Parameters
% ----------
% tubi : TubULAR class instance
% options : struct with fields 
%   overwrite : bool, default=false
%       overwrite previous results
%   preview : bool, default=false
%       view intermediate results
%   t0 : int, default=tubi.t0set()
%       timestamp at which the pathlines form a grid onto mesh vertices, 
%       mesh face barycenters, or PIV evaluation points
%
% Returns
% -------
% outStruct = struct('v2dum', v2dMum, ...
%     'v2d', v2dM, ...
%     'vn', vnM, ...
%     'v3d', vM, ...
%     'vf', vfM, ... 
%     'vv', vvM, ...
%     'v2dsmum', v2dsmMum, ... 
%     'v2dsm', v2dsmM, ...
%     'vnsm', vnsmM, ...
%     'v3dsm', vsmM, ...
%     'vfsm', vfsmM, ... 
%     'vvsm', vvsmM ) ;
% 
% Saves to disk
% -------------
% fvsmM   = sprintf(tubi.fileName.pathlines.velocities.v3dsm, t0) ;
% fvfsmM  = sprintf(tubi.fileName.pathlines.velocities.vfsm, t0) ;
% fvnsmM  = sprintf(tubi.fileName.pathlines.velocities.vnsm, t0) ;
% fvvsmM  = sprintf(tubi.fileName.pathlines.velocities.vvsm, t0) ;
% fv2dsmM = sprintf(tubi.fileName.pathlines.velocities.v2dsm, t0) ;
% fv2dsmMum = sprintf(tubi.fileName.pathlines.velocities.v2dsmum, t0) ;
%
% NPMitchell 2020

%% Default options
overwrite = false ;
overwriteImages = false ;
timePoints = tubi.xp.fileMeta.timePoints ;
samplingResolution = '1x' ;
imethod = 'linear' ;
if nargin < 2
    options = struct() ;
end

%% Unpack options
% Default values for options are to use sphi smoothed extended coords
% as PIV reference coord sys
pivimCoords = tubi.piv.imCoords ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if strcmp(pivimCoords(end), 'e')
    doubleCovered = true ;
else
    doubleCovered = false ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack tubi
t0 = tubi.t0set() ;
timePoints = tubi.xp.fileMeta.timePoints ;

%% Perform/Load Lagrangian averaging along pathline
disp('Performing/Loading Lagrangian averaged velocities')
% Create directories
fileNames = tubi.fileName.pathlines.velocities ;
% Apply t0 to fileNames
fieldnames = fields(fileNames) ;
for qq = 1:length(fieldnames)
    fileNames.(fieldnames{qq}) = sprintfm(fileNames.(fieldnames{qq}), t0) ;
end

% Check if the time smoothed velocities exist already
% 2d velocities (pulled back), scaled by dilation of metric are v2dum 
% 2D velocities (pulled back) are v2d
% normal velocities on fieldfaces are vn
% 3d velocities on fieldfaces are v3d
% vertex-based velocities are vv
% face-based velocities are vf

tubi.clearTime() ;

%% Build grids for averaging
if ~exist(fileNames.v2dum, 'file') || ~exist(fileNames.v2d, 'file') || ...
        ~exist(fileNames.vn, 'file') || ~exist(fileNames.v3d, 'file') || ...
        ~exist(fileNames.vf, 'file') || ~exist(fileNames.vv, 'file') || ...
        ~exist(fileNames.v2dsmum, 'file') || ~exist(fileNames.v2dsm, 'file') || ...
        ~exist(fileNames.vnsm, 'file') || ~exist(fileNames.v3dsm, 'file') || ...
        ~exist(fileNames.vfsm, 'file') || ~exist(fileNames.vvsm, 'file') || ...
        overwrite
    
    disp('Could not find time-smoothed velocities on disk')
    disp('Computing them...')
    
    % Load lagrangian pathlines
    disp('Loading pathlines for velocity computation: XY')
    load(sprintfm(tubi.fileName.pathlines.XY, t0), 'pivPathlines')
    load(sprintfm(tubi.fileName.pathlines.vXY, t0), 'vertexPathlines')
    load(sprintfm(tubi.fileName.pathlines.fXY, t0), 'facePathlines')
    
    ntps = length(timePoints)-1;
    
    %% Now load 3d piv and smooth in Lagrangian coords (along streamlines)
    first = true ; 
    for tidx = 1:ntps
        tp = timePoints(tidx) ;
        % dt = timePoints(i + 1) - tp ;
        
        disp(['Filling in velocity matrices, t=' num2str(tp)])
        tubi.setTime(tp) ;
        if doubleResolution
            tubi.getCurrentVelocity('piv3d2x') ;
            piv3d = tubi.currentVelocity.piv3d2x ;
        else
            tubi.getCurrentVelocity('piv3d') ;
            piv3d = tubi.currentVelocity.piv3d ;
        end
        
        % Allocate memory if this is the first timestep. Assume all grids
        % are equally sized.
        if first
            vM = zeros(ntps, size(piv3d.v0_rs, 1), size(piv3d.v0_rs, 2));
            vfM = zeros(ntps, size(piv3d.v3dfaces, 1), size(piv3d.v3dfaces, 2)); 
            vvM = zeros(ntps, ...
                size(piv3d.v3dvertices, 1), ...
                size(piv3d.v3dvertices, 2)); 
            vnM = zeros(ntps, size(piv3d.v0n_rs, 1), size(piv3d.v0n_rs, 2));
            v2dM = zeros(ntps, size(piv3d.v0t2d, 1), size(piv3d.v0t2d, 2));
            v2dMum = zeros(ntps, size(piv3d.v0t2d, 1), size(piv3d.v0t2d, 2));
            first = false ;
        end
        
        % Assert no NaNs    
        try
            assert(~any(isnan(piv3d.v0_rs(:))))
            assert(~any(isnan(piv3d.v3dfaces_rs(:))))
            assert(~any(isnan(piv3d.v0n_rs(:))))
            assert(~any(isnan(piv3d.v0t2d(:))))
        catch
           % disp('inpainting NaNs in pt0 & pt1')
           error(['There are NaNs in the velocity data. Could use ', ...
               'inpaint_nans, but why/how are they there?'])
           % pt0 = inpaint_nans(pt0) ;
           % pt1 = inpaint_nans(pt1) ;
           close all
           figure ;
           scatter(piv3d.x0(:), piv3d.y0(:), 10, piv3d.v0_rs(:, 1))
           bad = find(isnan(piv3d.v0_rs(:, 1))) ;
           hold on; 
           xx1 = piv3d.x0(:) ;
           yy1 = piv3d.y0(:) ;
           scatter(xx1(bad), yy1(bad), 30, 'k')
        end
        
        %% For all pathlines, query velocity at this current time
        % Query velocities at eval coords, mesh vertices, mesh faces via 
        % interpolation. Also query vn, v2d (ie vt), and v2dum (vt/|g|). 
               
        % Load streamline positions at this timePoint
        XX = pivPathlines.XX(tidx, :, :) ;
        YY = pivPathlines.YY(tidx, :, :) ;
        fX = facePathlines.fX(tidx, :, :) ;
        fY = facePathlines.fY(tidx, :, :) ;
        vX = vertexPathlines.vX(tidx, :, :) ;
        vY = vertexPathlines.vY(tidx, :, :) ;
        
        % Store x0,y0 PIV evaluation positions in pixels
        x0 = piv3d.x0 ;
        y0 = piv3d.y0 ;
        
        % 1. Interpolate v0_rs
        v0_rsX = reshape(piv3d.v0_rs(:, 1), size(x0)) ;
        v0_rsY = reshape(piv3d.v0_rs(:, 2), size(x0)) ;
        v0_rsZ = reshape(piv3d.v0_rs(:, 3), size(x0)) ;
        Fx = griddedInterpolant(x0', y0', v0_rsX', imethod, 'nearest') ;
        Fy = griddedInterpolant(x0', y0', v0_rsY', imethod, 'nearest') ;
        Fz = griddedInterpolant(x0', y0', v0_rsZ', imethod, 'nearest') ;
        % Query velocities at PIV evaluation pathlines
        v0_rs = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:)), Fz(XX(:), YY(:))] ;
        
        % 2. Interpolate velocities onto face barycenter pathlines 
        % Query velocities on face barycenters
        v3dfaces_rs = [Fx(fX(:), fY(:)), Fy(fX(:), fY(:)), Fz(fX(:), fY(:))] ;

        % 3. Interpolate v0n_rs onto PIV evaluation pt pathlines
        v0nrsN = reshape(piv3d.v0n_rs, size(x0)) ;
        Fn = griddedInterpolant(x0', y0', v0nrsN', imethod, 'nearest') ;
        % Query velocities on PIV evaluation pathlines
        v0n_rs = Fn(XX(:), YY(:)) ;

        % 4. Interpolate onto mesh vertices (v3dvertices)
        % Query velocities on vertices
        v3dvertices = [Fx(vX(:), vY(:)), Fy(vX(:), vY(:)), Fz(vX(:), vY(:))] ;
        
        % 4b: Could Query velocities on PIV evaluation pathlines
        % Chose to instead rely on metric kinematics measurement for this,
        % which resides in:
        %    fullfile(tubi.dir.metricKinematics.pathline.measurements, ...
        %       sprintfm('veln_pathline%04d_%06d.mat', t0Pathline, tp)))
        % v0nvertices_rs = Fn(vX(:), vY(:)) ;   --> see above     

        % 5. Interpolate v0t2d onto PIV evaluation pt pathlines
        v0t2dX = reshape(piv3d.v0t2d(:, 1), size(x0)) ;
        v0t2dY = reshape(piv3d.v0t2d(:, 2), size(x0)) ;
        Fx = griddedInterpolant(x0', y0', v0t2dX', imethod, 'nearest') ;
        Fy = griddedInterpolant(x0', y0', v0t2dY', imethod, 'nearest') ;
        % Query velocities at PIV evaluation pt pathlines
        v0t2d = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:))] ;

        % 6. Interpolate v0t2dum onto PIV evaluation pt pathlines
        v0t2dx = reshape(piv3d.v0t2d(:, 1) ./ piv3d.dilation, size(x0)) ;
        v0t2dy = reshape(piv3d.v0t2d(:, 2) ./ piv3d.dilation, size(x0)) ;
        Fx = griddedInterpolant(x0', y0', v0t2dx', imethod, 'nearest') ;
        Fy = griddedInterpolant(x0', y0', v0t2dy', imethod, 'nearest') ;
        % Query velocities at PIV evaluation pt pathlines
        v0t2dum = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:))] ;

        %% BUILD ARRAYS
        vM(tidx, :, :) = v0_rs ;             % in um/dt rs at PIV evaluation pathlines
        vfM(tidx, :, :) = v3dfaces_rs ;      % in um/dt rs at face barycenters
        vnM(tidx, :, :) = v0n_rs ;           % in um/dt rs at PIV evaluation pathlines
        vvM(tidx, :, :) = v3dvertices ;      % in um/min rs at vertex pathlines
        v2dM(tidx, :, :) = v0t2d ;           % in pixels/ min at PIV evaluation pathlines
        v2dMum(tidx, :, :) = v0t2dum ;       % in scaled pix/min, but proportional to um/min at PIV evaluation pathlines
        
    end
    
    clearvars first 
    disp('built v0 matrix')
    
    %% Save raw matrices
    outdir = sprintfm(tubi.dir.pathlines.velocities, t0) ;
    fvM   = sprintfm(tubi.fileName.pathlines.velocities.v3d, t0) ;
    fvfM  = sprintfm(tubi.fileName.pathlines.velocities.vf, t0) ;
    fvnM  = sprintfm(tubi.fileName.pathlines.velocities.vn, t0) ;
    fvvM  = sprintfm(tubi.fileName.pathlines.velocities.vv, t0) ;
    fv2dM = sprintfm(tubi.fileName.pathlines.velocities.v2d, t0) ;
    fv2dMum = sprintfm(tubi.fileName.pathlines.velocities.v2dum, t0) ; 
    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end
    disp(['Saving: ' fvM])
    save(fvM, 'vM')  
    disp(['Saving: ' fvfM])
    save(fvfM, 'vfM') 
    disp(['Saving: ' fvnM])
    save(fvnM, 'vnM') 
    disp(['Saving: ' fvvM])
    save(fvvM, 'vvM') 
    disp(['Saving: ' fv2dM])
    save(fv2dM, 'v2dM') 
    disp(['Saving: ' fv2dMum])
    save(fv2dMum, 'v2dMum')
    
    %% Filter in time axis -- light time smoothing tripulse of length 5
    % linfilt = 0.1 * ones(10, 1, 1) ;
    % ellipsoid = fspecial3('ellipsoid', [5, 1, 1]) ;
    disp('Building tripulse filter equivalent to tripuls()')
    tripulse3 = [ 0.3333; 0.6666; 1; 0.6666; 0.3333];
    tripulse3 = tripulse3 ./ sum(tripulse3(:)) ;
    tripulse3 = reshape(tripulse3, [length(tripulse3), 1]) ;
    
    % Check that no NaNs
    assert(~any(isnan(vM(:))))
    assert(~any(isnan(vvM(:))))
    assert(~any(isnan(vnM(:))))
    assert(~any(isnan(v2dM(:))))
    assert(~any(isnan(v2dMum(:))))
    assert(~any(isnan(vfM(:))))
    
    disp('filtering velocity matrices')
    vsmM = imfilter(vM, tripulse3, 'replicate');         % in um/min, rs
    vvsmM = imfilter(vvM, tripulse3, 'replicate');       % in um/min, rs
    vnsmM = imfilter(vnM, tripulse3, 'replicate');       % in um/min
    v2dsmM = imfilter(v2dM, tripulse3, 'replicate');     % in pix/min
    v2dsmMum = imfilter(v2dMum, tripulse3, 'replicate'); % in scaled pix/min, proportional to um/min  
    vfsmM = imfilter(vfM, tripulse3, 'replicate') ;      % in um/min, rs

    % Save the simpleminded averaging
    disp('Saving the time-smoothed velocities to disk')
    fvsmM   = sprintfm(tubi.fileName.pathlines.velocities.v3dsm, t0) ;
    fvfsmM  = sprintfm(tubi.fileName.pathlines.velocities.vfsm, t0) ;
    fvnsmM  = sprintfm(tubi.fileName.pathlines.velocities.vnsm, t0) ;
    fvvsmM  = sprintfm(tubi.fileName.pathlines.velocities.vvsm, t0) ;
    fv2dsmM = sprintfm(tubi.fileName.pathlines.velocities.v2dsm, t0) ;
    fv2dsmMum = sprintfm(tubi.fileName.pathlines.velocities.v2dsmum, t0) ;
    disp(['Saving: ' fvsmM])
    save(fvsmM, 'vsmM') ;          % in um/min
    disp(['Saving: ' fvvsmM])
    save(fvvsmM, 'vvsmM') ;        % in um/min, rs
    disp(['Saving: ' fvfsmM])
    save(fvfsmM, 'vfsmM') ;        % in um/min, rs
    disp(['Saving: ' fvnsmM])
    save(fvnsmM, 'vnsmM') ;        % in um/min
    disp(['Saving: ' fv2dsmM])
    save(fv2dsmM, 'v2dsmM') ;      % in pix/min
    disp(['Saving: ' fv2dsmMum])
    save(fv2dsmMum, 'v2dsmMum') ;  % in scaled pix/min, proportional to um/min 
    
    overwriteImages = true ;
    
    % package into struct for output
    if nargout > 0
        outStruct = struct('v2dum', v2dMum, ...
            'v2d', v2dM, ...
            'vn', vnM, ...
            'v3d', vM, ...
            'vf', vfM, ... 
            'vv', vvM, ...
            'v2dsmum', v2dsmMum, ... 
            'v2dsm', v2dsmM, ...
            'vnsm', vnsmM, ...
            'v3dsm', vsmM, ...
            'vfsm', vfsmM, ... 
            'vvsm', vvsmM ) ;
    end
    
    disp('done')
else
    disp('Lagrangian-pathline-averaged velocities already found on disk.')
    
    % Load these in case we plot them
    if nargout > 0 || overwriteImages
        load(fileNames.v2dsmum, 'v2dsmMum') ; 
        load(fileNames.vnsm, 'vnsmM') ; 
        load(fileNames.v3dsm, 'vsmM')  ;
    end
    if nargout > 0
        load(fileNames.v2dum, 'v2dMum') ;
        load(fileNames.v2d, 'v2dM') ;
        load(fileNames.vn, 'vnM') ;
        load(fileNames.v3d, 'vM') ;
        load(fileNames.vf, 'vfM') ; 
        load(fileNames.vv, 'vvM')  ;
        load(fileNames.v2dsm, 'v2dsmM')  ;
        load(fileNames.vfsm, 'vfsmM') ; 
        load(fileNames.vvsm, 'vvsmM') ; 
        outStruct = struct('v2dum', v2dMum, ...
            'v2d', v2dM, ...
            'vn', vnM, ...
            'v3d', vM, ...
            'vf', vfM, ... 
            'vv', vvM, ...
            'v2dsmum', v2dsmMum, ... 
            'v2dsm', v2dsmM, ...
            'vnsm', vnsmM, ...
            'v3dsm', vsmM, ...
            'vfsm', vfsmM, ... 
            'vvsm', vvsmM ) ;
    end
end


%% Plot each timepoint's pathline velocity
if overwriteImages
    disp('Loading pathlines for velocity computation: XY')
    load(sprintfm(tubi.fileName.pathlines.XY, t0), 'pivPathlines')
    % load(sprintfm(tubi.fileName.pathlines.vXY, t0), 'vertexPathlines')
    % load(sprintfm(tubi.fileName.pathlines.fXY, t0), 'facePathlines')
    
    sampleIDxFn = fullfile(sprintfm(tubi.dir.pathlines.velocities, t0), ...
        'sampleIDx_for_plotting_XY.mat') ;
    if ~exist(sampleIDxFn, 'file')
        xx = pivPathlines.XX(tubi.xp.tIdx(t0), :, :) ;
        yy = pivPathlines.YY(tubi.xp.tIdx(t0), :, :) ;
        disp('defining points at which to draw arrows')
        if length(xx(:)) > 200
            % seed a point in the center
            findpt = [mean(xx(:)), mean(yy(:))] ;
            seedIDx = pointMatch(findpt, [xx(:), yy(:)]) ;
            % Define faces for farthest point search
            ff = delaunay(xx(:), yy(:)) ;
            sampleIDx = farthestPointVertexSampling(200, ...
                [xx(:), yy(:)], ff, seedIDx, struct('preview', true)) ;
        else
            sampleIDx = 1:length(xx(:)) ;
        end
        save(sampleIDxFn, 'sampleIDx')
    else
        disp('Loading vertex IDs at which to plot vector field')
        load(sampleIDxFn, 'sampleIDx')
    end
    
    for tidx = 1:length(timePoints)-1
        tp = timePoints(tidx) ;
        tubi.setTime(tp) ;
        close all
        % Load streamline positions at this timePoint
        plotOptions.XX = pivPathlines.XX(tidx, :, :) ;
        plotOptions.YY = pivPathlines.YY(tidx, :, :) ;
        % Load velocities to plot
        plotOptions.vsm = squeeze(vsmM(tidx, :, :))  ;
        plotOptions.vnsm = squeeze(vnsmM(tidx, :, :))  ;
        plotOptions.v2dsmum = squeeze(v2dsmMum(tidx, :, :)) ;
        plotOptions.overwrite = overwriteImages ;
        plotOptions.sampleIDx = sampleIDx ;
        tubi.plotPathlineVelocitiesTimePoint(tp, plotOptions) ;
    end
end