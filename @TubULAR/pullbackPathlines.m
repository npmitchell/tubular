function [XX, YY] = pullbackPathlines(QS, x0, y0, t0, options)
% pullbackPathlines(QS, x0, y0, t0, options)
%   Create paths in pullback space (in pixels, XY) by following optical
%   flow of PIV measured on QS.pivimCoords pullbacks.
%   For non-standard PIV pathline propagation (such as uvprime coords),
%   supply piv, Lx, and Ly in options.
%   Note: I've chosen to use spatial smoothing via Gaussian blur rather
%   than temporal smoothing of velocities (as is done in metric kinematics)
%   for this code. This can be changed by setting 
%   QS.piv.smoothing_sigma == 0. Ie, if QS.piv.smoothing_sigma > 0, then 
%   we advect in the smoothed piv results. 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% x0 : n*m float array 
%   x coordinates in pullback pixels to start pathlines at t0
% y0 : n*m float array 
%   y coordinates in pullback pixels to start pathlines at t0
% t0 : 
%   time at which to begin the pathlines, must be member of
%   QS.xp.fileMeta.timePoints
% options : struct with fields 
%   preview : bool
%       view intermediate results
%   timePoints : 1d int array
%       the timpepoints along which to map pathlines that intersect with
%       (x0, y0) at t0
% 
% Returns
% -------
% XX : #timePoints x size(x0, 1) x size(x0, 2) float array
%   the pullback pixel X coordinates of the pathline spanning timePoints
% YY : #timePoints x size(x0, 1) x size(x0, 2) float array
%   the pullback pixel Y coordinates of the pathline spanning timePoints
%
% NPMitchell 2020

%% Input checking
assert(all(size(x0) == size(y0))) 

%% Default options
preview = false ;
debug = false ;
pivimCoords = QS.piv.imCoords ;
timePoints = QS.xp.fileMeta.timePoints ;
samplingResolution = '1x' ;

%% Unpack options
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
end
% Doesn't actually matter whether pullback is a single or double cover.
% if strcmp(pivimCoords(end-1), 'e')
%     doubleCovered = true ;
% else
%     doubleCovered = false ;
% end


%% Set it up
first = true ;
ntps = length(timePoints)-1;

%% First load raw PIV and store (uu, vv) values in a grid
% Load raw PIV for (uu, vv) in pix/dt
% Load the positions of the velocity vectors in pixels

if isfield(options, 'piv') 
    disp('Using supplied PIV results')
    piv = options.piv ;
else
    disp('Loading raw PIV results')
    QS.getPIV()
    if QS.piv.smoothing_sigma > 0
        piv = QS.piv.smoothed ;
    else
        piv = QS.piv.raw ;
    end
end
if isfield(options, 'Lx')
    disp('QS.pullbackPathlines(): using supplied Lx')
    Lxs = options.Lx ;
else
    Lxs = QS.piv.Lx ;
end
if isfield(options, 'Ly')
    disp('QS.pullbackPathlines(): using supplied Ly')
    Lys = options.Ly ;
else
    Lys = QS.piv.Ly ;
end
disp('Building pathlines')

%% Collate velocities into 4d array
if debug
    % Template example for debugging
    % Debug this function: make fake PIV
    x0 = piv.x{1} ;
    y0 = piv.y{1} ;
    vPIV = zeros(ntps, size(x0, 1), size(x0, 2), 2);
    for ii = 1:ntps
        vPIV(ii, :, :, 1) = 100 * cos(2 * pi * ii / 20) + 5 ;
        vPIV(ii, :, :, 2) = 10 ;
    end
else
    % Load in true PIV
    for ii = 1:ntps
        tidx = QS.xp.tIdx(timePoints(ii)) ; 
        uu = piv.u_filtered{tidx} ;
        vv = piv.v_filtered{tidx} ; 

        % Ensure no NaNs in uu and vv
        if any(isnan(uu(:))) || any(isnan(vv(:)))
           disp('inpainting NaNs in uu & vv')
           uu = inpaint_nans(uu) ;
           vv = inpaint_nans(vv) ;
        end

        % Build PIV velocity grid (time, ugrid(1), ugrid(2), x/y)
        if first
            vPIV = zeros(ntps, size(uu, 1), size(uu, 2), 2);
            xpiv = piv.x{tidx} ;
            ypiv = piv.y{tidx} ;
            first = false ;
        else
            % Ensure that evaluation gridpts are same throughout
            assert(all(all(xpiv == piv.x{tidx}))) 
            assert(all(all(ypiv == piv.y{tidx}))) 
        end

        % Only give nonzero velocity if the next timepoint is different
        % than this one. Note that if we have timePoints=[1 2 3 3], then
        % allowing the fourth entry of vPIV(4, :, :) to be nonzero does no
        % harm, since we backward propagate STARTING at ii=3, not ii=4. 
        if ii < ntps
            if tidx < QS.xp.tIdx(timePoints(ii+1)) 
                vPIV(ii, :, :, 1) = uu ;             % in pix/dt
                vPIV(ii, :, :, 2) = vv ;             % in pix/dt
            end
        else
            vPIV(ii, :, :, 1) = uu ;             % in pix/dt
            vPIV(ii, :, :, 2) = vv ;             % in pix/dt            
        end
    end
end

%% Propagate along velocities forward and backward

% Could use streamline but there is an issue with moving out of the
% frame. Instead use griddedInterpolant with padded edges, clip at each
% step along the way to keep in frame, and use periodic BC for vertical
% direction
% Preallocate positions for all time
XX = zeros(length(timePoints), size(x0, 1), size(x0, 2)) ;
YY = zeros(length(timePoints), size(y0, 1), size(y0, 2)) ;
% Fill in position at starting time t0
idx0 = find(timePoints == t0) ;  % QS.xp.tIdx(t0) ;

% Permit duplicate timePoints for timeaverging at t~0 or t~end, but print 
% warning when doing so. Do not permit duplicate timePoints in the middle
% of the array
if length(idx0) > 1
    disp("WARNING: QS.pullbackPathlines(): multiple timepoints are identical and equal t0")
    for qq = 1:length(idx0)
        ttmp = idx0(qq) ;
        XX(ttmp, :, :) = x0 ;
        YY(ttmp, :, :) = y0 ;
    end
    if ~any(idx0 == 1) && ~any(idx0 == length(timePoints))
        error("WARNING: duplicate points in the middle of timePoint list")
    end
else
    XX(idx0, :, :) = x0 ;
    YY(idx0, :, :) = y0 ;
end

% Now we are done with input handling and timepoint t0.
if any(diff(timePoints) == 0)
    disp('pausing here')
end

% Propagate forward first: tIdx(t0)+2 onward
for qq = (max(idx0)+1):length(timePoints)
    disp(['tidx = ' num2str(qq)])
    % 1. Interpolate velocities at time qq-1
    uu = squeeze(vPIV(qq-1, :, :, 1)) ;
    vv = squeeze(vPIV(qq-1, :, :, 2)) ;
    % in transposed coords
    ui = griddedInterpolant(xpiv', ypiv', uu', 'linear', 'nearest') ; 
    vi = griddedInterpolant(xpiv', ypiv', vv', 'linear', 'nearest') ; 

    % 2. Evaluate at XY(qq-1) non-transposed coords
    xx = squeeze(XX(qq-1, :, :)) ;
    yy = squeeze(YY(qq-1, :, :)) ;
    % assert(all(abs(xx(:)) > 0))
    dx = reshape(ui(xx(:), yy(:)), size(xx)) ;
    dy = reshape(vi(xx(:), yy(:)), size(yy)) ;

    % 3. push XY
    Xqq = xx + dx ;
    Yqq = yy + dy ;

    % 4. Clip at x=0,Lx and wrap at y=0,2Ly
    % Load Lx, Ly from vector (match timepoint)
    tidx_match = find(QS.xp.fileMeta.timePoints == timePoints(qq)) ;
    try
        assert(length(tidx_match) == 1)
    catch
        error(['Multiple timepoints in the QS instance match ', ...
            'this timepoint. QS time must be monotonic'])
    end
    Lx = Lxs(tidx_match) ;
    Ly = Lys(tidx_match) ;
    [Xqq, Yqq] = QS.clipXY(Xqq, Yqq, Lx, Ly) ;
    XX(qq, :, :) = Xqq ;
    YY(qq, :, :) = Yqq ;

    if debug
        subplot(1, 2, 1)
        plot(qq, dx(1, 19, 19), 'b.')
        plot(qq, dy(1, 19, 19), 'r.')
        hold on;
        subplot(1, 2, 2) 
        plot(qq, XX(qq, 19, 19), 'b.')
        plot(qq, YY(qq, 19, 19), 'r.')
        hold on;

        %clf 
        %subplot(1, 2, 1)
        %imagesc(squeeze(dx))
        % caxis([0, 1])
        %title(['t = ' num2str(qq)])
        %colorbar()
        %pause(0.2)
    end
end

% Propagate backward in time if t0 > timePoints(2)
if min(idx0) > 1 
    backward_times = fliplr( 1:(min(idx0)-1) ) ;
    for qq = backward_times 
        disp(['tidx = ' num2str(qq)])
        % 1. Interpolate velocities of time qq at their advected 
        %    locations in qq+1.
        %
        %   .---->*     .<----*
        %     qq=.       qq+1=*
        % 
        %    Since the advected x0+u_qq,y0+v_qq are unstructured, we
        %    interpolate current velocity at the advected positions,
        %    ie what velocities got you there, evaluate the velocities 
        %    that pushed to the next positions at the next positions,
        %    and pull back along those velocities.
        %
        uu = squeeze(vPIV(qq, :, :, 1)) ;
        vv = squeeze(vPIV(qq, :, :, 2)) ;
        
        % Load Lx, Ly from vector (match timepoint)
        tidx_match = find(QS.xp.fileMeta.timePoints == timePoints(qq)) ;
        try
            assert(length(tidx_match) == 1)
        catch
            error(['Multiple timepoints in the QS instance match ', ...
                'this timepoint. QS time must be monotonic'])
        end
        Lx = Lxs(tidx_match) ;
        Ly = Lys(tidx_match) ;
        [xa, ya] = QS.clipXY(xpiv(:) + uu(:), ypiv(:) + vv(:), Lx, Ly) ;
        ui = scatteredInterpolant(xa, ya, uu(:), 'natural', 'nearest') ;
        vi = scatteredInterpolant(xa, ya, vv(:), 'natural', 'nearest') ;

        % 2. Evaluate at XY(qq+1) non-transposed coords 
        xx = squeeze(XX(qq+1, :, :)) ;
        yy = squeeze(YY(qq+1, :, :)) ;

        % 3. Pull XY back
        dx = reshape(ui(xx(:), yy(:)), size(xx)) ;
        dy = reshape(vi(xx(:), yy(:)), size(xx)) ;
        Xqq = xx - dx ;
        Yqq = yy - dy ;

        % 4. Clip at x=0,Lx and wrap at y=0,2Ly
        [Xqq, Yqq] = QS.clipXY(Xqq, Yqq, Lx, Ly) ;
        XX(qq, :, :) = Xqq ;
        YY(qq, :, :) = Yqq ;
    end
end

%% Check debugged pathlines
if debug
    if preview
        for qq = 1:50:size(x0, 1)
            for pp = 1:50:size(x0, 2)
                plot3(1:length(timePoints), XX(:, qq, pp), YY(:, qq, pp), '-')
                hold on;
            end
        end
        xlabel('t')
        ylabel('X [pix]')
        zlabel('Y [pix]')
        title('Artificial flow for streamline checking')
    end
    view(2)
    saveas(gcf, fullfile(QS.dir.pivAvg, 'streamline_test_xt.png'))
    view(90, 90)    
    saveas(gcf, fullfile(QS.dir.pivAvg, 'streamline_test_yt.png'))
    view(90, 0)
    saveas(gcf, fullfile(QS.dir.pivAvg, 'streamline_test_xy.png'))
    waitfor(gcf)
end