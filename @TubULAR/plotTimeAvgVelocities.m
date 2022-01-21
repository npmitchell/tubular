function plotTimeAvgVelocities(QS, options) 
%plotTimeAvgVelSimple(QS, options) 
%   Save images of the velocity field over time that has been "simply"
%   averaged over time in-place in the surface-Lagrangian pullback. That
%   is, the velocity field at location (u,v) has been averaged in time with
%   previous and later timepoints at the same pullback location (u,v).
%   For samplingResolution, all fields plotted in this method depend only 
%   on the PIV sampling settings
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%   alphaVal : float
%       the opacity of the heatmap to overlay
%   invertImage : bool
%       invert the data pullback under the velocity map
%   vtscale = 5 ;      % um / min
%   vnscale = 2 ;      % um / min
%   vscale = 2 ;       % um / min
%   alphaVal = 0.7 ;   % alpha for normal velocity heatmap
%   washout2d = 0.5 ;  % lightening factor for data
%   qsubsample = 5 ;   % quiver subsampling in pullback space 
%
%
%
% NPMitchell 2020

%% Default options
% Declare plotting options for limits
plot_vxyz = false ;      
pivimCoords = QS.piv.imCoords ;  % coordinate system of the pullback images used in PIV
averagingStyle = 'Lagrangian' ;  % Lagrangian or simple, how velocities are averaged over time
samplingResolution = '1x' ;      % 1x or 2x, resolution of 
vtscale = 0 ;                    % if zero, default is used
vnscale = 0 ;                    % if zero, default is used
vscale = 0 ;                     % if zero, default is used
invertImage = false ;            % invert the data underneath velocity heatmaps
washout2d = 0.5 ;                % washout the data image under velocity heatmaps
%% Unpack options
if isfield(options, 'plot_vxyz')
    plot_vxyz = options.plot_vxyz ;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
end
if isfield(options, 'vtscale')
    vtscale = options.vtscale ;
end
if isfield(options, 'vnscale')
    vnscale = options.vnscale ;
end
if isfield(options, 'vscale')
    vscale = options.vscale ;
end
if isfield(options, 'invertImage')
    invertImage = options.invertImage ;
end
if isfield(options, 'washout2d')
    washout2d = options.washout2d ;
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

%% Unpack QS
if strcmp(averagingStyle, 'Lagrangian')
    if doubleResolution
        pivDir = QS.dir.piv.avg2x ;
    else
        pivDir = QS.dir.piv.avg ;
    end
elseif strcmp(averagingStyle, 'simple')
    if doubleResolution
        pivDir = QS.dir.pivSimAvg2x ;
    else
        pivDir = QS.dir.pivSimAvg ;
    end    
end
timePoints = QS.xp.fileMeta.timePoints ;

% Auxiliary function for plotting smoothed meshes

%% Define directories
pivImTDir = fullfile(pivDir, 'vtH') ;  % Heatmap
pivImGDir = fullfile(pivDir, 'vtG ') ;  % Gaussian smoothed in space
pivImSDir = fullfile(pivDir, 'vmag') ;  % speed |v_3D|
pivImNDir = fullfile(pivDir, 'vn') ;
% pivImQDir = fullfile(pivSimAvgDir, 'vtQ') ;  % Quiverplot
% pivImDvgDir = fullfile(pivSimAvgDir, 'dvg') ;
% pivImRotDir = fullfile(pivSimAvgDir, 'rot') ;
% pivImShearDir = fullfile(pivSimAvgDir, 'shear_dvphi_ds') ;

pivDir = QS.dir.piv.root ;
dilDir = fullfile(pivDir, 'dilation') ;
vxyorigDir = fullfile(pivDir, 'vxyorig') ;
if plot_vxyz
    ensureDir(pivAImXDir)
    ensureDir(pivAImYDir)
    ensureDir(pivAImZDir)
end
dirs2make = {pivImTDir, pivImGDir,...
    pivImNDir, dilDir, vxyorigDir, pivImSDir} ;
for pp = 1:length(dirs2make)
    ensureDir(dirs2make{pp}) ;
end

%% Load simple/Lagrangian average piv results
if strcmp(averagingStyle, 'Lagrangian')
    if doubleResolution
        disp('Loading double Resolution velocity sampling')
        QS.getVelocityAverage2x()
        velstruct = QS.velocityAverage2x ;
    else
        disp('Loading single Resolution velocity sampling')
        QS.getVelocityAverage()
        velstruct = QS.velocityAverage ;
    end
elseif strcmp(averagingStyle, 'Simple')
    if doubleResolution
        disp('Loading double Resolution simple velocity sampling')
        QS.getVelocitySimpleAverage2x()
        velstruct = QS.velocitySimpleAverage2x ;
    else
        disp('Loading single Resolution simple velocity sampling')
        QS.getVelocitySimpleAverage()
        velstruct = QS.velocitySimpleAverage ;
    end
end
vsmM = velstruct.v3d ;
v2dsmMum = velstruct.v2dum ;
vnsmM = velstruct.vn ;

%% Make plots
% Display the velocities
close all
fig = figure('visible', 'off') ;
t2do = 1:10:size(vsmM, 1) ;
t2do = [t2do, setdiff(1:size(vsmM, 1), t2do)] ;
for i = t2do
    tp = timePoints(i) ;
    disp(['t = ', num2str(tp)])
    
    % grab the tangential velocity for this timestep
    vsm_ii = squeeze(vsmM(i, :, :)) ;
    v2dsmum_ii = squeeze(v2dsmMum(i, :, :)) ;
    vnsm_ii = squeeze(vnsmM(i, :, :)) ;
    
    options.vsm = vsm_ii ;
    options.v2dsmum = v2dsmum_ii ;
    options.vnsm = vnsm_ii ;
    options.vtscale = vtscale ;
    options.vnscale = vnscale ;
    options.invertImage = invertImage ;
    QS.plotAverageVelocitiesTimePoint(tp, options)   
        
end

if strcmp(averagingStyle, 'simple')
    disp('Done plotting in-place-avgeraged velocitites')
else
    disp('Done plotting Lagrangian-avgeraged velocitites')
end


% RETIRED CODE: DIVERGENCE, CURL, and SHEAR (these are crude, replaced
% by DEC) 
%
%
% shearfn = fullfile(pivSimAvgImShearDir, [sprintf('%04d', tp) '.png']) ;
%
% % Plot divergence
% if ~exist(dvgfn, 'file') || overwrite 
% 
%     vxb = imgaussfilt(vx, 10) ;
%     vyb = imgaussfilt(vy, 10) ;
% 
%     % Interpolate dilation onto locations where curl is defined
%     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
%     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
%     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
%     % dilv = dilation_allfaces(subfieldfaces) ;
% 
%     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
%     dvg = divergence(piv.x{i}, piv.y{i}, vxb, vyb) .* dilum;
%     opts.label = '$\nabla \cdot v_t$ [min$^{-1}$]' ;
%     opts.title = [ '$t=$' num2str(tp) ' min' ] ;
%     opts.outfn = dvgfn ;
%     opts.qscale = 10 ;
%     opts.sscale = 1.0 ;
%     opts.alpha = 0.8 ;
%     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
%     scalarVectorFieldsOnImage(im, xx, yy, ...
%         dvg, xx, yy, vxb, vyb, opts) ;
%     clearvars opts
% end
% 
% % Plot curl
% if ~exist(curlfn, 'file') || overwrite || true
%     % xd = xx(1:10:end) ;
%     % yd = yy(1:10:end) ;
%     % xdgrid = xd .* ones(length(xd), length(yd)) ; 
%     % ydgrid = (xd .* ones(length(yd), length(xd)))' ; 
%     % check it
%     % scatter(xdgrid(:), ydgrid(:), 10, xdgrid(:))
% 
%     % Interpolate dilation onto locations where curl is defined
%     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
%     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
%     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
%     % dilv = dilation_allfaces(subfieldfaces) ;
% 
%     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
%     curlv = curl(piv.x{i}, piv.y{i}, vxb, vyb) ;
%     curlv = curlv .* dilum ;
%     opts.title = [ '$t=$' num2str(tp) ' min' ] ;
%     opts.label = '$\nabla \times v_t$ [min$^{-1}$]' ;
%     opts.outfn = curlfn ;
%     opts.qscale = 10 ;
%     opts.sscale = 0.5 ;
%     opts.alpha = 0.8 ;
%     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
%     scalarVectorFieldsOnImage(im, xx, yy, ...
%         curlv, xx, yy, vxb, vyb, opts) ;
%     clearvars opts
% end
% 
% % Plot strain dv_phi/ds on image
% if ~exist(shearfn, 'file') || overwrite || true
%     % xd = xx(1:10:end) ;
%     % yd = yy(1:10:end) ;
%     % xdgrid = xd .* ones(length(xd), length(yd)) ; 
%     % ydgrid = (xd .* ones(length(yd), length(xd)))' ; 
%     % check it
%     % scatter(xdgrid(:), ydgrid(:), 10, xdgrid(:))
% 
%     % Interpolate dilation onto locations where curl is defined
%     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
%     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
%     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
%     % dilv = dilation_allfaces(subfieldfaces) ;
% 
%     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
%     dvphidX = gradient(vyb, xx(2) - xx(1), yy(2) - yy(1)) ;
%     dvphids = dvphidX .* dilum ;
%     opts.title = [ '$t=$' num2str(tp) ' min' ] ;
%     opts.label = '$\nabla_s v_{\phi}$ [min$^{-1}$]' ;
%     opts.outfn = shearfn ;
%     opts.qscale = 10 ;
%     opts.sscale = 0.5 ;
%     opts.alpha = 0.8 ;
%     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
%     scalarVectorFieldsOnImage(im, xx, yy, ...
%         dvphids, xx, yy, vxb, vyb, opts) ;
%     clearvars opts
% end