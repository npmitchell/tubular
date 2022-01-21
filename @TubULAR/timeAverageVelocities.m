function timeAverageVelocities(QS, options)
%measurePullbackStreamlines(QS, options)
%   Use pathlines of optical flow in pullback space to query velocities
%   and average along Lagrangian pathlines. 
%   Default is to weight velocities contributions with a time width 5 
%   weighted by tripulse filter.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   plotOptions : optional struct with fields
%       vtscale : float
%       vnscale : float
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%
% Saves to disk
% -------------
% vsmM : (#timePoints-1) x (nX*nY) x 3 float array
%   3d velocities at PIV evaluation coordinates in spaceUnits/dt rs
% vfsmM : (#timePoints-1) x (2*nU*(nV-1)) x 3 float array
%   3d velocities at face barycenters in um/dt rs
% vnsmM : (#timePoints-1) x (nX*nY) float array
%   normal velocity at PIV evaluation coordinates in spaceUnits/dt rs
% vvsmM : (#timePoints-1) x (nU*nV) x 3 float array
%   3d velocities at (1x resolution) mesh vertices in spaceUnits/dt rs
% v2dsmM : (#timePoints-1) x (nX*nY) x 2 float array
%   2d velocities at PIV evaluation coordinates in pixels/dt
% v2dsmMum : (#timePoints-1) x (nX*nY) x 2 float array
%   2d velocities at PIV evaluation coordinates in scaled pix/dt, but 
%   proportional to spaceUnits/dt (scaled by dilation of map)
%
% NPMitchell 2020

%% Default options
overwrite = false ;
overwriteImages = false ;
preview = true ;
timePoints = QS.xp.fileMeta.timePoints ;
pivimCoords = QS.piv.imCoords ;
plotOptions = struct() ;
samplingResolution = '1x' ;
XYkernel = 0 ;        % sigma for light gaussian smoothing on full velocity 
                      % fields in XY pixel space before smoothing in time
imethod = 'linear' ;  % interpolation method for velocities onto pathlines
twidth = 2 ;          % average over (t-twidth, t+twidth) timepoints

%% Unpack options
% Default values for options are to use sphi smoothed extended coords
% as PIV reference coord sys
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'plotOptions')
    plotOptions = options.plotOptions ;
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
if isfield(options, 'XYkernel')
    XYkernel = options.XYkernel ;
end
if isfield(options, 'imethod')
    imethod = options.imethod ;
end
if isfield(options, 'twidth')
    twidth = options.twidth ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
pivDir = QS.dir.piv ;
piv3dfn = QS.fullFileBase.piv3d ;
ntps = length(timePoints) ;
% [rot, ~] = QS.getRotTrans() ;
% resolution = QS.APDV.resolution ; 
[~, ~, ~, xyzlim_APDV] = QS.getXYZLims() ;
axis_order = QS.data.axisOrder ;
blue = QS.plotting.colors(1, :) ;
red = QS.plotting.colors(2, :) ;
green = QS.plotting.colors(4, :) ;
t0 = QS.t0set() ;
timePoints = QS.xp.fileMeta.timePoints ;

%% Colormap
close all
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
close all

%% Perform/Load simple averaging
disp('Performing/Loading simple averaging')
% Create directories
if doubleResolution
    fileNames = QS.fileName.pivAvg2x ;
else
    fileNames = QS.fileName.pivAvg ;
end
% Check if the time smoothed velocities exist already
% 2d velocities (pulled back), scaled by dilation of metric are v2dum 
% 2D velocities (pulled back) are v2d
% normal velocities on fieldfaces are vn
% 3d velocities on fieldfaces are v3d
% vertex-based velocities are vv
% face-based velocities are vf

QS.clearTime() ;

%% Build grids for averaging
if ~exist(fileNames.v2dum, 'file') || ~exist(fileNames.v2d, 'file') || ...
        ~exist(fileNames.vn, 'file') || ~exist(fileNames.v3d, 'file') || ...
        ~exist(fileNames.vf, 'file') || ~exist(fileNames.vv, 'file') || ...
        overwrite
    
    disp('Could not find time-smoothed velocities on disk')
    disp('Computing them...')
                
    ntps = length(timePoints)-1;
    
    %% Now load 3d piv and smooth in Lagrangian coords (along streamlines)
    first = true ; 
    for tidx = 1:ntps
        tp = timePoints(tidx) ;
        % dt = timePoints(i + 1) - tp ;
        
        disp(['Filling in velocity matrices, t=' num2str(tp)])
        QS.setTime(tp) ;
        if doubleResolution
            QS.getCurrentVelocity('piv3d2x') ;
            piv3d = QS.currentVelocity.piv3d2x ;
        else
            QS.getCurrentVelocity('piv3d') ;
            piv3d = QS.currentVelocity.piv3d ;
        end
        
        % Allocate memory if this is the first timestep. Assume all grids
        % are equally sized.
        if first
            vsmM = zeros(ntps, size(piv3d.v0_rs, 1), size(piv3d.v0_rs, 2));
            vfsmM = zeros(ntps, size(piv3d.v3dfaces, 1), size(piv3d.v3dfaces, 2)); 
            vvsmM = zeros(ntps, ...
                size(piv3d.v3dvertices, 1), ...
                size(piv3d.v3dvertices, 2)); 
            vnsmM = zeros(ntps, size(piv3d.v0n_rs, 1), size(piv3d.v0n_rs, 2));
            v2dsmM = zeros(ntps, size(piv3d.v0t2d, 1), size(piv3d.v0t2d, 2));
            v2dsmMum = zeros(ntps, size(piv3d.v0t2d, 1), size(piv3d.v0t2d, 2));
            first = false ;
        end
        
        % Assert no NaNs    
        try
            assert(~any(isnan(piv3d.v0_rs(:))))
            assert(~any(isnan(piv3d.v3dfaces_rs(:))))
            assert(~any(isnan(piv3d.v0n_rs(:))))
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
        
        % Compute short lagrangian pathlines (t-twidth, t+twidth)
        % Note that we clip the timepoints to compute at 1 and ntps.
        disp(['Computing pathlines around t=', num2str(tp)])
        tp2do = (tp - twidth):(tp + twidth) ;
        tp2do = min(max(timePoints(1), tp2do), timePoints(end-1)) ;
        popts.timePoints = tp2do ;
        
        % PIV pathlines
        [XXpath, YYpath] = QS.pullbackPathlines(piv3d.x0, piv3d.y0, tp, popts) ;
        % face pathlines
        if strcmp(QS.piv.imCoords, 'sp_sme')
            im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
            mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp), 'spcutMeshSmRS') ;
            mesh0 = mesh0.spcutMeshSmRS ;
            umax = max(mesh0.u(:, 1)) ;
            vmax = max(mesh0.u(:, 2)) ;
            mXY = QS.uv2XY(im0, mesh0.u, doubleCovered, umax, vmax) ;
        end
        bc = barycenter(mXY, mesh0.f) ; % this is in gptoolbox
        [fXpath, fYpath] = QS.pullbackPathlines(bc(:, 1), bc(:, 2), tp, popts) ;
        % vertex pathlines
        [vXpath, vYpath] = QS.pullbackPathlines(mXY(:, 1), mXY(:, 2), tp, popts) ;
        
        % Check eval points
        % clf
        % plot(mXY(:, 1), mXY(:, 2), '.')
            
        %% Run through all timepoints in the window and compute a grid of
        % velocities to average over -- Grab PIV for each timepoint and 
        % interpolate onto pathline segment for all timepoints in tp2do]
        
        %% Filter in time: multiply each contribution by weight of tripulse
        % linfilt = 0.1 * ones(10, 1, 1) ;
        % ellipsoid = fspecial3('ellipsoid', [5, 1, 1]) ;
        disp('Building tripulse filter equivalent to tripuls()')
        if twidth == 2
            tripulse = [ 0.3333; 0.6666; 1; 0.6666; 0.3333];
            tripulse = tripulse ./ sum(tripulse(:)) ;
            % tripulse = reshape(tripulse, [length(tripulse), 1]) ;
        elseif twidth == 1
            tripulse = [ 0.5; 1; 0.5];
            tripulse = tripulse ./ sum(tripulse(:)) ;
        elseif twidth == 0
            tripulse = [ 1 ] ;
        else
            error(['build tripulse of twidth ' num2str(twidth) ' here'])
        end
        
        % Preallocate
        vMtp  = zeros(size(vsmM, 2), size(vsmM, 3)) ;    % in um/dt rs at PIV evaluation points
        vfMtp = zeros(size(vfsmM, 2), size(vfsmM, 3)) ;  % in um/dt rs at face barycenters
        vnMtp = zeros(size(vnsmM, 2), size(vnsmM, 3)) ;  % in um/dt rs at 
        vvMtp = zeros(size(vvsmM, 2), size(vvsmM, 3)) ;  % in um/min rs
        v2dMtp   = zeros(size(v2dsmM, 2), size(v2dsmM, 3)) ;  % in pixels/ min
        v2dMumtp = zeros(size(v2dsmMum, 2), size(v2dsmMum, 3)) ;  % in scaled pix/min, but proportional to um/min
        for qq = 1:length(tp2do)
            disp(['timeAverageVelocities: Interpolating t0=', ...
                num2str(tp) ' coords onto pathline at t=', ...
                num2str(tp2do(qq))])
            % Load streamline positions at this timePoint
            XX = XXpath(qq, :, :) ;
            YY = YYpath(qq, :, :) ;
            fX = fXpath(qq, :, :) ;
            fY = fYpath(qq, :, :) ;
            vX = vXpath(qq, :, :) ;
            vY = vYpath(qq, :, :) ;
            
            % Grab PIV for this timepoint and interpolate onto pathline
            % segment for all timepoints in tp2do          
            QS.setTime(tp2do(qq)) ;
            
            % Load piv if this is not the last timePoint
            if doubleResolution 
                QS.getCurrentVelocity('piv3d2x') ;
                piv3d = QS.currentVelocity.piv3d2x ;
            else
                QS.getCurrentVelocity('piv3d') ;
                piv3d = QS.currentVelocity.piv3d ;
            end

            % Store x0,y0 PIV evaluation positions in pixels
            x0 = piv3d.x0 ;
            y0 = piv3d.y0 ;

            % 1. Interpolate v0_rs
            v0_rsX = reshape(piv3d.v0_rs(:, 1), size(x0)) ;
            v0_rsY = reshape(piv3d.v0_rs(:, 2), size(x0)) ;
            v0_rsZ = reshape(piv3d.v0_rs(:, 3), size(x0)) ;
            
            % Optional: light smoothing with image kernel
            if XYkernel > 0
                v0_rsX = imgaussfilt(v0_rsX, XYkernel) ;
                v0_rsY = imgaussfilt(v0_rsY, XYkernel) ;
                v0_rsZ = imgaussfilt(v0_rsZ, XYkernel) ;
            end
            Fx = griddedInterpolant(x0', y0', v0_rsX', imethod, 'nearest') ;
            Fy = griddedInterpolant(x0', y0', v0_rsY', imethod, 'nearest') ;
            Fz = griddedInterpolant(x0', y0', v0_rsZ', imethod, 'nearest') ;
            % Query velocities
            v0_rs = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:)), Fz(XX(:), YY(:))] ;

            % 2. Interpolate velocities onto face barycenter pathlines 
            % Query velocities
            v3dfaces_rs = [Fx(fX(:), fY(:)), Fy(fX(:), fY(:)), Fz(fX(:), fY(:))] ;

            % 3. Interpolate v0n_rs
            v0nrsN = reshape(piv3d.v0n_rs, size(x0)) ;
            Fn = griddedInterpolant(x0', y0', v0nrsN', imethod, 'nearest') ;
            % Query velocities
            v0n_rs = Fn(XX(:), YY(:)) ;
            
            % Check it
            % imagesc(x0(:, 1), y0(1, :), reshape(piv3d.v0n_rs, size(x0)))
            % caxis([-2, 2]); colormap(bwr256)
            % imagesc(x0(:, 1), y0(1, :), reshape(v0n_rs, size(x0)))
            % caxis([-2, 2]); colormap(bwr256)
            % imagesc(x0(:, 1), y0(1, :), reshape(vnMtp, size(x0)))
            % caxis([-2, 2]); colormap(bwr256)

            % 4. Interpolate onto mesh vertices (v3dvertices)
            % Query velocities
            v3dvertices = [Fx(vX(:), vY(:)), Fy(vX(:), vY(:)), Fz(vX(:), vY(:))] ;

            % Check it
            if preview
                close all
                % mesh0.vrs = QS.xyz2APDV(mesh0.v);
                % vn_rs = QS.dx2APDV(mesh0.vn) ;
                % mesh0.vn_rs = vn_rs ./ vecnorm(vn_rs, 2, 2) ;
                subplot(3, 1, 1)
                vvn = dot(mesh0.vn, v3dvertices, 2) ;
                trisurf(triangulation(mesh0.f, mesh0.v), vvn, ...
                    'edgecolor', 'none') ;
                caxis([-2,2]); axis equal; colormap(bwr256); colorbar()
                title(['Computed result for ' num2str(tp2do(qq))])
                view(2)
                subplot(3, 1, 2)
                tmp = load(QS.fileName.pivAvg.vv, 'vvsmM') ;
                tmpv = squeeze(tmp.vvsmM(tidx, :, :)) ;
                vvn2 = dot(mesh0.vn, tmpv, 2) ;
                trisurf(triangulation(mesh0.f, mesh0.v), vvn2, ...
                    'edgecolor', 'none') ;
                caxis([-2,2]); axis equal; colormap(bwr256); colorbar()
                title('Loaded result of average from disk')
                view(2)
                subplot(3, 1, 3)
                trisurf(triangulation(mesh0.f, mesh0.v), vvn-vvn2, ...
                    'edgecolor', 'none') ;
                caxis([-2,2]); axis equal; colormap(bwr256); colorbar()
                title('difference from loaded result')
                view(2)
                pause(1) ;
            end
            
            % 5. Interpolate v0t2d
            v0t2dX = reshape(piv3d.v0t2d(:, 1), size(x0)) ;
            v0t2dY = reshape(piv3d.v0t2d(:, 2), size(x0)) ;
            Fx = griddedInterpolant(x0', y0', v0t2dX', imethod, 'nearest') ;
            Fy = griddedInterpolant(x0', y0', v0t2dY', imethod, 'nearest') ;
            % Query velocities
            v0t2d = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:))] ;
            
            % 6. Interpolate v0t2dum
            v0t2dx = reshape(piv3d.v0t2d(:, 1) ./ piv3d.dilation, size(x0)) ;
            v0t2dy = reshape(piv3d.v0t2d(:, 2) ./ piv3d.dilation, size(x0)) ;
            Fx = griddedInterpolant(x0', y0', v0t2dx', imethod, 'nearest') ;
            Fy = griddedInterpolant(x0', y0', v0t2dy', imethod, 'nearest') ;
            % Query velocities
            v0t2dum = [Fx(XX(:), YY(:)), Fy(XX(:), YY(:))] ;

            % Check it
            % if preview
            %     close all
            %     QS.getPIV(options) ;
            %     % get size of images to make
            %     gridsz = size(QS.piv.raw.x{1}) ;
            %     % Define Nx1 and Mx1 float arrays for xspace and yspace
            %     xx = QS.piv.raw.x{tidx}(1, :) ;
            %     yy = QS.piv.raw.y{tidx}(:, 1) ;
            %     % Load the image to put flow on top
            %     if strcmp(pivimCoords, 'sp_sme')
            %         im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
            %         ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
            %     else
            %         error(['Have not coded for this pivimCoords option. Do so here: ' pivimCoords])
            %     end
            %     im = cat(3, im, im, im) ;  % convert to rgb for no cmap change
            %     opts.fig = figure ;
            %     vectorFieldHeatPhaseOnImage(im, xx, yy, ...
            %         v0t2dum(:, 1), v0t2dum(:, 2), 5, opts) ;
            %     set(gcf, 'visible', 'on')
            %     opts.fig = figure ;
            %     vectorFieldHeatPhaseOnImage(im, xx, yy, ...
            %         v0t2dx, v0t2dY, 5, opts) ;
            %     set(gcf, 'visible', 'on')
            %     waitfor(gcf) ;
            % end
            
            %% BUILD ARRAYS -- AVERAGE ARRAY over time dimension w/ weight
            vMtp = vMtp + tripulse(qq) * v0_rs ;           % in um/dt rs at PIV evaluation points
            vfMtp = vfMtp + tripulse(qq) * v3dfaces_rs ;   % in um/dt rs at face barycenters
            vnMtp = vnMtp + tripulse(qq) * v0n_rs ;        % in um/dt rs at 
            vvMtp = vvMtp + tripulse(qq) * v3dvertices ;   % in um/min rs
            v2dMtp = v2dMtp + tripulse(qq) * v0t2d ;       % in pixels/ min
            v2dMumtp = v2dMumtp + tripulse(qq) * v0t2dum ; % in scaled pix/min, but proportional to um/min
        
            % Check it
        end
        
        %% BUILD ARRAYS by collating timepoints
        vsmM(tidx, :, :) = vMtp ;             % in um/dt rs, 3d velocities at PIV evaluation coordinates
        vfsmM(tidx, :, :) = vfMtp ;           % in um/dt rs, 3d velocities at face barycenters
        vnsmM(tidx, :, :) = vnMtp ;           % in um/dt rs, normal velocity at PIV evaluation coordinates
        vvsmM(tidx, :, :) = vvMtp ;           % in um/min rs, 3d velocities at mesh vertices
        v2dsmM(tidx, :, :) = v2dMtp ;         % in pixels/ min, 2d velocities at PIV evaluation coordinates
        v2dsmMum(tidx, :, :) = v2dMumtp ;     % in scaled pix/min, but proportional to um/min, 2d velocities at PIV evaluation coordinates

        % Check that no NaNs
        assert(~any(isnan(vsmM(:))))
        assert(~any(isnan(vvsmM(:))))
        assert(~any(isnan(vnsmM(:))))
        assert(~any(isnan(v2dsmM(:))))
        assert(~any(isnan(v2dsmMum(:))))
        assert(~any(isnan(vfsmM(:))))

        %% Plot this timepoint
        close all
        plotOptions.vsm = squeeze(vsmM(tidx, :, :))  ;
        plotOptions.vnsm = squeeze(vnsmM(tidx, :, :))  ;
        plotOptions.v2dsm = squeeze(v2dsmM(tidx, :, :)) ;
        plotOptions.v2dsmum = squeeze(v2dsmMum(tidx, :, :)) ;
        plotOptions.overwrite = overwriteImages ;
        QS.plotAverageVelocitiesTimePoint(tp, plotOptions) ;
    end
    
    clearvars first 
    disp('built v0 matrix')
    
    % Check that no NaNs
    assert(~any(isnan(vsmM(:))))
    assert(~any(isnan(vvsmM(:))))
    assert(~any(isnan(vnsmM(:))))
    assert(~any(isnan(v2dsmM(:))))
    assert(~any(isnan(v2dsmMum(:))))
    assert(~any(isnan(vfsmM(:))))
    
    %% Save raw matrices
    save(fileNames.v3d, 'vsmM') 
    save(fileNames.vf, 'vfsmM') 
    save(fileNames.vn, 'vnsmM') 
    save(fileNames.vv, 'vvsmM') 
    save(fileNames.v2d, 'v2dsmM') 
    save(fileNames.v2dum, 'v2dsmMum') 
        
    disp('done')

end