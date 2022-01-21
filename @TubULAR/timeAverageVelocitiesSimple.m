function timeAverageVelocitiesSimple(QS, samplingResolution, options)
%timeAverageVelocities(QS, options)
%   Average velocities in time with a tripulse fiter of full-width 5*dt 
%   (half width of 2.5*dt).
%   
%   Note that velocities are all saved as um/timeunits, 
%   except for v2dM, which is in full-resolution-data pixels / timeunits
%   and v2dMum, which is not truly in um, but in pixels/min/dilation, so it
%   is proportional to um/timeunits
%   
% Parameters
% ----------
% QS : QuapSlap class instance
% samplingResolution : str specifier ('1x' or '2x', 'single' or 'double)
%   whether to sample pullback velocities at nU x nV or 2*nU x 2*nV 
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%
%
% Returns
% -------
% Saved files on disk containing
%   vM  : , in um/min rs
%       3d velocities at PIV evaluation points in um/dt in APDV coords
%   vfM : ,  in um/min rs
%       3d velocities at face barycenters in um/dt in APDV coords
%   vnM : , in um/min rs
%       normal velocities at PIV evaluation points in um/dt in APDV coords
%   vvM : , in um
%   v2dM: , in pix/dt 
%       in pix/dt at PIV evaluation coords
%   v2dMum : , 
%       scaled pix/dt, but proportional to um/min at PIV evaluation coords
%
% NPMitchell 2020

%% Default options
overwrite = false ;
timePoints = QS.xp.fileMeta.timePoints ;

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
    if strcmp(pivimCoords(-1), 'e')
        doubleCovered = true ;
    else
        doubleCovered = false ;
    end
else
    pivimCoords = 'sp_sme' ;
    doubleCovered = true;
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

%% Perform/Load simple averaging
disp('Performing/Loading simple averaging')
% Create directories
if doubleResolution
    % pivSimAvgDir = QS.dir.pivSimAvg2x ;
    fileNames = QS.fileName.pivSimAvg2x ;
else
    % pivSimAvgDir = QS.dir.pivSimAvg ;
    fileNames = QS.fileName.pivSimAvg ;
end

% Check if the time smoothed velocities exist already
% 2d velocities (pulled back), scaled by dilation of metric are v2dum 
% 2D velocities (pulled back) are v2d
% normal velocities on fieldfaces are vn
% 3d velocities on fieldfaces are v3d
% vertex-based velocities are vv
% face-based velocities are vf

QS.clearTime() ;

if ~exist(fileNames.v2dum, 'file') || ~exist(fileNames.v2d, 'file') || ...
        ~exist(fileNames.vn, 'file') || ~exist(fileNames.v3d, 'file') || ...
        ~exist(fileNames.vf, 'file') || ~exist(fileNames.vv, 'file') || ...
        overwrite
    
    disp('Could not find time-smoothed velocities on disk')
    disp('Computing them...')
    first = true ;
    ntps = length(timePoints)-1;
    for i = 1:ntps
        tp = timePoints(i) ;
        dt = timePoints(i + 1) - tp ;
        
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
        
        % BUILD ARRAYS
        vM(i, :, :) = piv3d.v0_rs ;             % in um/min rs at PIV evaluation coords
        vfM(i, :, :) = piv3d.v3dfaces_rs ;      % in um/min rs at mesh face barycenters
        vnM(i, :, :) = piv3d.v0n_rs ;           % in um/min rs at PIV evaluation coords
        vvM(i, :, :) = QS.dx2APDV(piv3d.v3dvertices) ; % in um/min rs at mesh vertices
        v2dM(i, :, :) = piv3d.v0t2d ;           % in pix/dt at PIV evaluation coords
        v2dMum(i, :, 1) = piv3d.v0t2d(:, 1) ./ piv3d.dilation ; % in scaled pix/min, but proportional to um/min, at PIV evaluation coords
        v2dMum(i, :, 2) = piv3d.v0t2d(:, 2) ./ piv3d.dilation ; % in scaled pix/min, but proportional to um/min, at PIV evaluation coords
        
    end
    clearvars first 
    disp('built v0 matrix')
    % Filter in time axis
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
    save(fileNames.v2dum, 'v2dsmMum') ;  % in scaled pix/min, proportional to um/min 
    save(fileNames.v2d, 'v2dsmM') ;      % in pix/min
    save(fileNames.vn, 'vnsmM') ;        % in um/min
    save(fileNames.v3d, 'vsmM') ;          % in um/min
    save(fileNames.vv, 'vvsmM') ;          % in um/min, rs
    save(fileNames.vf, 'vfsmM') ;        % in um/min, rs

end

disp('done with timeAverageVelocitiesSimple()')

