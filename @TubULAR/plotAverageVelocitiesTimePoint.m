function plotAverageVelocitiesTimePoint(QS, tp, options)
% plotTimeAvgVelocitiesTimePoint(QS, tp, options)
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : optional struct with fields 
%   vnsm    : 
%   v2dsmum : 
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
%
%
% NPMitchell 2020

%% Unpack options
overwrite = false ;
vtscale = 5 ;                   % um / min
vnscale = 2 ;                   % um / min
vscale = 5 ;                    % um / min
alphaVal = 0.7 ;                % alpha for normal velocity heatmap
invertImage = false ;           % invert the data underneath velocity heatmaps
washout2d = 0.5 ;               % lightening factor for data if < 1
qsubsample = 10 ;               % quiver subsampling in pullback space 
pivimCoords = QS.piv.imCoords ; % coordinate system of the pullback images used in PIV
averagingStyle = 'Lagrangian' ; % Lagrangian or simple, how velocities are averaged over time
samplingResolution = '1x' ;     % 1x or 2x, resolution of 
v2dsmum_ii = options.v2dsmum ;
vnsm_ii = options.vnsm ;
vsm_ii = options.vsm ;

%% Unpack options
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'vtscale')
    if options.vtscale > 0
        vtscale = options.vtscale ;
    end
end
if isfield(options, 'vnscale')
    if options.vnscale > 0
        vnscale = options.vnscale ;
    end
end
if isfield(options, 'vscale')
    if options.vscale > 0
        vscale = options.vscale ;
    end
end
if isfield(options, 'invertImage')
    invertImage = options.invertImage ;
end
if isfield(options, 'washout2d')
    washout2d = options.washout2d ;
end
if isfield(options, 'alphaVal')
    alphaVal = options.alphaVal ;
end
if isfield(options, 'qsubsample')
    qsubsample = options.qsubsample ;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
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
tidx = QS.xp.tIdx(tp) ;

%% Load piv results
disp('Obtaining raw piv to get field size')
QS.getPIV(options) ;
piv = QS.piv.raw ;
% get size of images to make
gridsz = size(QS.piv.raw.x{1}) ;
% Define Nx1 and Mx1 float arrays for xspace and yspace
xx = QS.piv.raw.x{tidx}(1, :) ;
yy = QS.piv.raw.y{tidx}(:, 1) ;

%% Prepare for plots
QS.getFeatures('ssfold')
ssfold_frac = QS.features.ssfold / QS.nU ;
t0 = QS.t0set() ;
tunit = [' ' QS.timeUnits] ;

%% Figure out which averaging style directory to direct to
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

% Auxiliary function for plotting smoothed meshes

%% Define directories
pivImTDir = fullfile(pivDir, 'vtH') ;  % Heatmap
pivImGDir = fullfile(pivDir, 'vtG ') ;  % Gaussian smoothed in space
pivImSDir = fullfile(pivDir, 'vmag') ;  % speed |v_3D|
pivImNDir = fullfile(pivDir, 'vn') ;    % normal velocity
pivImGNDir = fullfile(pivDir, 'vtG_vn') ; % Combined smoothed tangential, plus normal vel
pivDir = QS.dir.piv.root ;
dilDir = fullfile(pivDir, 'dilation') ;
vxyorigDir = fullfile(pivDir, 'vxyorig') ;
% if plot_vxyz
%     ensureDir(pivAImXDir)
%     ensureDir(pivAImYDir)
%     ensureDir(pivAImZDir)
% end
dirs2make = {pivImTDir, pivImGDir,...
    pivImNDir, dilDir, vxyorigDir, pivImSDir, pivImGNDir} ;
for pp = 1:length(dirs2make)
    ensureDir(dirs2make{pp}) ;
end

%% Check if normal velocity plot exists
vnfn = fullfile(pivImNDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
vthfn = fullfile(pivImTDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
vtgfn = fullfile(pivImGDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
vtgvnfn = fullfile(pivImGNDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
vxyorigfn = fullfile(vxyorigDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
speedfn = fullfile(pivImSDir, [sprintf(QS.fileBase.name, tp) '.png']) ;

% Load the image to put flow on top
if strcmp(pivimCoords, 'sp_sme')
    im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
    ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
else
    error(['Have not coded for this pivimCoords option. Do so here: ' pivimCoords])
end
im = cat(3, im, im, im) ;  % convert to rgb for no cmap change

% Define the proper coordinates (approx)
% if ~exist(dvgfn, 'file') || ~exist(curlfn, 'file') || overwrite || true
%     dilation_allfaces = zeros(length(jac), 1) ;
%     for f = 1:length(jac)
%         qg = jac{f} * jac{f}' ;
%         dilation_allfaces(f) = sqrt(det(qg)) ;
%     end
% end
%     % load the spcutMesh
%     load(sprintf(spcutMeshBase, tp), 'spcutMesh')
%     ss = spcutMesh.sphi(:, 1) ;
%     pp = spcutMesh.sphi(:, 1) ;
%     [ss, pp] = gridDistancesInterpMesh3D() ;
% end
 
% if plot_vxyz
%     aux_plot_vxyz_simpleavg(im, vsm_ii, xx, yy, tp, vscale, ...
%         pivImXDir, pivImYDir, pivImZDir) ;
% end

% Plot the magnitude of the velocity 
if ~exist(speedfn, 'file') || overwrite
    disp(['Saving ' speedfn])
    close all
    fig = figure('units', 'normalized', ...
            'outerposition', [0 0 1 1], 'visible', 'off') ;
    colormap parula ;
    labelOpts.label = '$|v|$ [$\mu$m/min]' ;
    labelOpts.title = ['speed, $|v|$: $t=$' num2str(tp - t0) tunit] ;
    scalarFieldOnImage(im, [xx', yy], ...
        reshape(vecnorm(vsm_ii, 2, 2), gridsz),...
        alphaVal, vscale, labelOpts, 'Style', 'Positive') ;
    ylim([0.25 * size(im, 1), 0.75 * size(im, 1)])
    disp(['saving figure: ' speedfn])
    saveas(fig, speedfn) ;
    close all
end

% % Check normals
% quiver3(piv3d{i}.pt0(:, 1), piv3d{i}.pt0(:, 2), piv3d{i}.pt0(:, 3), ...
%     piv3d{i}.normals(:, 1), piv3d{i}.normals(:, 2), piv3d{i}.normals(:, 3)) ;
% 
% % Check normals rotated and scaled
% pt0_rs = ((rot * piv3d{i}.pt0')' + trans) * resolution  ;
% quiver3(pt0_rs(:, 1), pt0_rs(:, 2), pt0_rs(:, 3), ...
%     piv3d{i}.normals_rs(:, 1), piv3d{i}.normals_rs(:, 2), piv3d{i}.normals_rs(:, 3)) ;
%

% Look at smoothed 2d velocity fields
vn = reshape(vnsm_ii, gridsz) ;
vx = reshape(v2dsmum_ii(:, 1), gridsz) ;
vy = reshape(v2dsmum_ii(:, 2), gridsz) ;

% % Check normal velocity and rotations
% quiver3(piv.x{i}, piv.y{i}, 0*piv.x{i}, vx, vy, vn, 0)

% Get lobes for this timepoint
foldx = ssfold_frac(tidx, :) * size(im, 2) ;

% Plot the normal velocity on top
if ~exist(vnfn, 'file') || overwrite
    disp(['Saving ' vnfn])
    close all
    fig = figure('units', 'normalized', ...
            'outerposition', [0 0 1 1], 'visible', 'off') ;
    labelOpts.label = '$v_n$ [$\mu$m/min]' ;
    labelOpts.title = ['normal velocity, $v_n$: $t=$' num2str(tp - t0) tunit] ;
    labelOpts.cmap = twilight_shifted_mod(256) ;
    if invertImage
        imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
    else
        imw = im * washout2vnscaled + max(im(:)) * (1-washout2d) ;
    end
    scalarFieldOnImage(imw, [xx', yy], vn, alphaVal, vnscale, ...
        labelOpts) ;
    ylim(ylims)
    disp(['saving figure: ' vnfn])
    saveas(fig, vnfn) ;
    close all
end

% Plot the tangential velocity as heatmap on top of the image
if ~exist(vthfn, 'file') || overwrite
    disp(['Saving ' vthfn])
    if invertImage
        imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
    else
        imw = im * washout2d + max(im(:)) * (1-washout2d) ;
    end
    qopts.overlay_quiver = false ;
    qopts.qsubsample = qsubsample ;
    qopts.overlay_quiver = true ;
    qopts.qscale = 10 ;
    qopts.label = '$v_t$ [$\mu$m/min]' ;
    qopts.title = ['tangential velocity, $v_t$: $t=$' num2str(tp - t0) tunit] ;
    qopts.outfn = vthfn ;
    qopts.ylim = ylims ;
    xyf = struct() ;
    xyf.x = xx ;
    xyf.y = yy ;
    vectorFieldHeatPhaseOnImage(imw, xyf, vx, vy, vtscale, qopts) ;
    clear qopts 
end

% Gaussian smooth the velocities
if ~exist(vtgfn, 'file') || overwrite
    disp(['Saving ' vtgfn])
    vxb = imgaussfilt(vx, 4) ;
    vyb = imgaussfilt(vy, 4) ;
    if invertImage
        imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
    else
        imw = im * washout2d + max(im(:)) * (1-washout2d) ;
    end
    qopts.qsubsample = qsubsample ;
    qopts.overlay_quiver = true ;
    qopts.qscale = 10 ;
    qopts.label = '$v_t$ [$\mu$m/min]' ;
    qopts.title = ['Smoothed velocity, $v_t$: $t=$', ...
        num2str(tp - t0), tunit] ;
    qopts.outfn = vtgfn ;
    qopts.ylim = ylims ;
    xyf = struct() ;
    xyf.x = xx ;
    xyf.y = yy ;
    % Plot the coarse-grained tangential velocity as heatmap on top of im
    vectorFieldHeatPhaseOnImage(imw, xyf, vxb, vyb, vtscale, qopts) ;    
    clearvars qopts
end

    
% Combined with normal velocity -- Gaussian smooth the velocities
if ~exist(vtgvnfn, 'file') || overwrite || true 
    disp(['Saving ' vtgvnfn])
    
    close all
    clearvars qopts
    fig = figure('units', 'centimeters', ...
            'outerposition', [0 0 18 18], 'visible', 'off') ;
    ax1 = subplot(1, 2, 1) ;
       
    % VT
    vxb = imgaussfilt(vx, 4) ;
    vyb = imgaussfilt(vy, 4) ;
    if invertImage
        imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
    else
        imw = im * washout2d + max(im(:)) * (1-washout2d) ;
    end
    qopts.qsubsample = qsubsample ;
    qopts.overlay_quiver = true ;
    qopts.qscale = 10 ;
    qopts.label = '$v_\parallel$ [$\mu$m/min]' ;
    qopts.title = ['tangential tissue velocity, $v_\parallel$'] ;
    qopts.ylim = ylims ;
    qopts.axPosition = [0.1 0.55, 0.85, 0.4] ;
    qopts.cbPosition = [.9 .55 .02 .2] ;
    qopts.pbPosition = [0.9, 0.85, 0.07, 0.07] ;
    qopts.ax = ax1 ;
    xyf = struct() ;
    xyf.x = xx ;
    xyf.y = yy ;
    % Plot the coarse-grained tangential velocity as heatmap on top of im
    vectorFieldHeatPhaseOnImage(imw, xyf, vxb, vyb, vtscale, qopts) ;    
    clearvars qopts
    
    saveas(gcf, vtgvnfn) ;
    panel1 = imread(vtgvnfn) ; 
    
    % VN
    close all
    fig = figure('units', 'centimeters', ...
            'outerposition', [0 0 18 18], 'visible', 'off') ;
    axes('Position', [0.1 0.05, 0.85, 0.4])
    % qopts.axPosition = [0.2 0.55, 0.85, 0.4] ;
    labelOpts.cbPosition = [.9 .1 .02 .3] ;
    labelOpts.label = '$v_n$ [$\mu$m/min]' ;
    labelOpts.title = ['normal tissue velocity, $v_n$'] ;
    labelOpts.cmap = twilight_shifted_mod(256) ;
    if invertImage
        imw = max(im(:)) - im ;
        % imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
    else
        imw = im * washout2vnscaled + max(im(:)) * (1-washout2d) ;
    end
    scalarFieldOnImage(imw, [xx', yy], vn, alphaVal, vnscale, ...
        labelOpts) ;
    ylim(ylims)
    disp(['saving figure: ' vtgvnfn])
    saveas(gcf, vtgvnfn) ;
    close all
    panel2 = imread(vtgvnfn) ;
    sz = size(panel1) ;
    sz2 = size(panel2) ;
    assert(all(sz == sz2)) 
    outim = zeros(sz) ;
    szhalf = round(0.5 * sz(1));
    outim(1:szhalf,:, :) = mat2gray(panel1(1:szhalf,:, :) );
    outim(szhalf:end, :, :) = mat2gray(panel2(szhalf:end, :, :));
    % Write combined image now
    imwrite(outim, vtgvnfn)
end

% Find hyperbolic fixed points
%

% Plot original velocity -- note no smoothing done here ;)
if ~exist(vxyorigfn, 'file') || overwrite
    if invertImage
        imw = (max(im(:))-im) * washout2d + max(im(:)) * (1-washout2d) ;
    else
        imw = im * washout2d + max(im(:)) * (1-washout2d) ;
    end
    opts.label = '$\tilde{v}$ [pix/min]' ;
    opts.outfn = vxyorigfn ;
    opts.qscale = 15 ;
    opts.ylim = ylims ;
    opts.title = ['piv results, $t=$', num2str(tp - t0), tunit] ;
    xyf = struct() ;
    xyf.x = xx ;
    xyf.y = yy ;
    vectorFieldHeatPhaseOnImage(imw, xyf, ...
        piv.u_filtered{tidx}, piv.v_filtered{tidx}, 15, opts) ;
    clearvars opts
end