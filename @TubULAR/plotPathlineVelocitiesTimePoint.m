function plotPathlineVelocitiesTimePoint(QS, tp, options)
% plotPathlineVelocitiesTimePoint(QS, tp, options)
%   Plot the velocities defined on PIV evaluation points in Lagrangian
%   coordinates.
%   
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : optional struct with fields 
%   gridTopology : str
%       'rectilinear' or 'triangulated'
%       triangulate anew each timepoint or keep grid structure
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
vtscale = 5 ;      % um / min
vnscale = 2 ;      % um / min
vscale = 5 ;       % um / min
alphaVal = 0.7 ;   % alpha for normal velocity heatmap
washout2d = .5 ;   % lightening factor for data if < 1
qsubsample = 5 ;   % quiver subsampling in pullback space 
pivimCoords = QS.piv.imCoords ;  % coordinate system of the pullback images used in PIV
samplingResolution = '1x' ;      % 1x or 2x, resolution of 
gridTopology = 'triangulated' ;  % ('rectilinear' or 'triangulated') triangulate anew each timepoint or keep grid structure

%% Unpack options
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'vtscale')
    vtscale = options.vtscale ;
end
if isfield(options, 'vscale')
    vscale = options.vscale ;
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
if strcmp(pivimCoords(end), 'e')
    doubleCovered = true ;
else
    doubleCovered = false ;
end

%% Load or unpack data to plot
if isfield(options, 'v2dsmum')
    v2dsmum_ii = options.v2dsmum ;
end
if isfield(options, 'vnsm')
    vnsm_ii = options.vnsm ;
end
if isfield(options, 'vsm')
    vsm_ii = options.vsm ;
end
if isfield(options, 'XX')
    XX = options.XX ;
end
if isfield(options, 'YY')
    YY = options.YY ;
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

%% Prepare for plots
QS.getFeatures('ssfold')
ssfold_frac = QS.features.ssfold / QS.nU ;
t0 = QS.t0set() ;
tunit = [' ' QS.timeUnits] ;

%% Figure out which averaging style directory to direct to
if doubleResolution
    pivDir = sprintf(QS.dir.pathlines2x.velocities, t0) ;
else
    pivDir = sprintf(QS.dir.pathlines.velocities, t0) ;
end

% Auxiliary function for plotting smoothed meshes

%% Define directories
pivImTDir = fullfile(pivDir, 'vtH') ;  % Heatmap
pivImSDir = fullfile(pivDir, 'vmag') ;  % speed |v_3D|
pivImNDir = fullfile(pivDir, 'vn') ;
dirs2make = {pivImTDir,...
    pivImNDir, pivImSDir} ;
for pp = 1:length(dirs2make)
    ensureDir(dirs2make{pp}) ;
end

%% Check if normal velocity plot exists
vnfn = fullfile(pivImNDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
vthfn = fullfile(pivImTDir, [sprintf(QS.fileBase.name, tp) '.png']) ;
speedfn = fullfile(pivImSDir, [sprintf(QS.fileBase.name, tp) '.png']) ;

% Load the image to put flow on top
if strcmp(pivimCoords, 'sp_sme')
    im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
    ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
else
    error(['Have not coded for this pivimCoords option. Do so here: ' pivimCoords])
end
im = cat(3, im, im, im) ;  % convert to rgb for no cmap change

%% Plot the magnitude of the velocity at pathlines XX,YY
if ~exist(speedfn, 'file') || overwrite
    disp(['Saving ' speedfn])
    close all
    fig = figure('units', 'normalized', ...
            'outerposition', [0 0 1 1], 'visible', 'off') ;
    colormap parula ;
    labelOpts.label = '$|v|$ [$\mu$m/min]' ;
    labelOpts.title = ['speed, $|v|$: $t=$' num2str(tp - t0) tunit] ;
    xyf.v = [XX(:), YY(:)] ;
    if contains(lower(gridTopology), 'rect')
        % OPTION 1: keep rectilinear grid topology
        xyf.f = defineFacesRectilinearGrid([], size(XX, 1), size(XX, 2)) ;    
        out_of_bounds = YY(:) > 0.9 * size(im, 1) | YY(:) < 0.1 * size(im,1) ;
        out_of_bounds = find(out_of_bounds) ;
        % Remove out-of-bounds vertices from triangulation (pathlines are
        % periodic in Y, and these faces will be annoying to plot)
        xyf.f(any(ismember(xyf.f, out_of_bounds), 2), :) = [] ;
        [ xyf.f, xyf.v, oldVertexIDx ] = ...
            remove_vertex_from_mesh(xyf.f, xyf.v, out_of_bounds) ;
    else    
        % OPTION 2: re-triangulation
        xyf.f = delaunay(XX(:), YY(:)) ;
        oldVertexIDx = 1:size(vsm_ii, 1) ;
    end
    scalarFieldOnImage(im, xyf, vecnorm(vsm_ii(oldVertexIDx, :), 2, 2),...
        alphaVal, vscale, labelOpts, 'Style', 'Positive') ;
    
    ylim([0.25 * size(im, 1), 0.75 * size(im, 1)])
    
    % Save the figure
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
vx = v2dsmum_ii(:, 1) ;
vy = v2dsmum_ii(:, 2) ;

% % Check normal velocity and rotations
% quiver3(piv.x{i}, piv.y{i}, 0*piv.x{i}, vx, vy, vn, 0)

% Get lobes for this timepoint
foldx = ssfold_frac(tidx, :) * size(im, 2) ;

%% Plot the normal velocity on top at pathlines XX,YY
if ~exist(vnfn, 'file') || overwrite
    disp(['Saving ' vnfn])
    close all
    fig = figure('units', 'normalized', ...
            'outerposition', [0 0 1 1], 'visible', 'off') ;
    labelOpts.label = '$v_n$ [$\mu$m/min]' ;
    labelOpts.title = ['normal velocity, $v_n$: $t=$' num2str(tp - t0) tunit] ;
    
    xyf.v = [XX(:), YY(:)] ;
    if contains(lower(gridTopology), 'rect')
        % OPTION 1: keep rectilinear grid topology
        xyf.f = defineFacesRectilinearGrid([], size(XX, 1), size(XX, 2)) ;    
        out_of_bounds = YY(:) > 0.9 * size(im, 1) | YY(:) < 0.1 * size(im,1) ;
        out_of_bounds = find(out_of_bounds) ;
        % Remove out-of-bounds vertices from triangulation (pathlines are
        % periodic in Y, and these faces will be annoying to plot)
        xyf.f(any(ismember(xyf.f, out_of_bounds), 2), :) = [] ;
        [ xyf.f, xyf.v, oldVertexIDx ] = ...
            remove_vertex_from_mesh(xyf.f, xyf.v, out_of_bounds) ;
    else    
        % OPTION 2: re-triangulation
        xyf.f = delaunay(XX(:), YY(:)) ;
        oldVertexIDx = 1:numel(vnsm_ii) ;
    end
    scalarFieldOnImage(im, xyf, vnsm_ii(oldVertexIDx), alphaVal, ...
        vnscale, labelOpts) ;

    ylim([0.25 * size(im, 1), 0.75 * size(im, 1)])
    
    disp(['saving figure: ' vnfn])
    saveas(fig, vnfn) ;
    close all
end

%% Plot the tangential velocity as heatmap on top of the image at pathlines
if ~exist(vthfn, 'file') || overwrite
    disp(['Saving ' vthfn])
    imw = im * washout2d + max(im(:)) * (1-washout2d) ;
    qopts.overlay_quiver = false ;
    qopts.qsubsample = qsubsample ;
    qopts.overlay_quiver = true ;
    qopts.qscale = 10 ;
    qopts.label = '$v_t$ [$\mu$m/min]' ;
    qopts.title = ['tangential velocity, $v_t$: $t=$' num2str(tp - t0) tunit] ;
    qopts.outfn = vthfn ;
    qopts.ylim = ylims ;
    qopts.nPts = 200 ;
    
    % Use custom subsampling from supplied
    if isfield(options, 'sampleIDx')
        qopts.sampleIDx = options.sampleIDx ;
        qopts.subsamplingMethod = 'custom' ;
    else
        qopts.subsamplingMethod = 'farthestPoint' ;
    end
    xyf.v = [XX(:), YY(:)] ;
    
    if contains(lower(gridTopology), 'rect')
        % OPTION 1: keep rectilinear grid topology
        xyf.f = defineFacesRectilinearGrid([], size(XX, 1), size(XX, 2)) ;    
        out_of_bounds = YY(:) > 0.9 * size(im, 1) | YY(:) < 0.1 * size(im,1) ;
        out_of_bounds = find(out_of_bounds) ;
        % Remove out-of-bounds vertices from triangulation (pathlines are
        % periodic in Y, and these faces will be annoying to plot)
        xyf.f(any(ismember(xyf.f, out_of_bounds), 2), :) = [] ;
    else    
        % OPTION 2: re-triangulation
        xyf.f = delaunay(XX(:), YY(:)) ;
    end
    
    %% Now plot it
    xyf.x = xyf.v(:, 1) ;
    xyf.y = xyf.v(:, 2) ;
    qopts.subsamplingMethod = 'random'; 
    vectorFieldHeatPhaseOnImage(imw, xyf, vx, vy, vtscale, qopts) ;
    clear qopts 
end

% Find hyperbolic fixed points
%