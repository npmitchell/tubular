function plotPathlineVelocities(QS, options)
%plotPathlineVelocities(QS, options)
%   Use pathlines of optical flow in pullback space to query velocities
%   and average along pathlines.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool, default=false
%       overwrite previous results
%   preview : bool, default=false
%       view intermediate results
%   t0 : int, default=QS.t0set()
%       timestamp at which the pathlines form a grid onto mesh vertices, 
%       mesh face barycenters, or PIV evaluation points
%
%
% NPMitchell 2020

%% Default options
overwrite = false ;
timePoints = QS.xp.fileMeta.timePoints ;
samplingResolution = '1x' ;

%% Unpack options
% Default values for options are to use sphi smoothed extended coords
% as PIV reference coord sys
pivimCoords = QS.piv.imCoords ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
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
if isfield(options, 't0')
    t0 = options.t0 ;
else
    t0 = QS.t0set() ;
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
timePoints = QS.xp.fileMeta.timePoints ;

%% Perform/Load Lagrangian averaging along pathline
disp('Loading simple averaging')
% Create directories
fileNames = QS.fileName.pathlines.velocities ;
% Apply t0 to fileNames
fieldnames = fields(fileNames) ;
for qq = 1:length(fieldnames)
    fileNames.(fieldnames{qq}) = sprintf(fileNames.(fieldnames{qq}), t0) ;
end

% Check if the time smoothed velocities exist already
% 2d velocities (pulled back), scaled by dilation of metric are v2dum 
% 2D velocities (pulled back) are v2d
% normal velocities on fieldfaces are vn
% 3d velocities on fieldfaces are v3d
% vertex-based velocities are vv
% face-based velocities are vf

%% Build grids for averaging
% Load lagrangian pathlines
disp('Loading pathlines for plotting: XY')
load(sprintf(QS.fileName.pathlines.XY, t0), 'pivPathlines')
load(sprintf(QS.fileName.pathlines.vXY, t0), 'vertexPathlines')
load(sprintf(QS.fileName.pathlines.fXY, t0), 'facePathlines')

ntps = length(timePoints)-1;

%% Now load Lagrangian vels (along streamlines) & plot them
disp('Saving images of time-smoothed velocities to disk')
fvsmM   = sprintf(QS.fileName.pathlines.velocities.v3dsm, t0) ;
fvfsmM  = sprintf(QS.fileName.pathlines.velocities.vfsm, t0) ;
fvnsmM  = sprintf(QS.fileName.pathlines.velocities.vnsm, t0) ;
fvvsmM  = sprintf(QS.fileName.pathlines.velocities.vvsm, t0) ;
fv2dsmM = sprintf(QS.fileName.pathlines.velocities.v2dsm, t0) ;
fv2dsmMum = sprintf(QS.fileName.pathlines.velocities.v2dsmum, t0) ;
load(fvsmM, 'vsmM') ;          % in um/min
load(fvvsmM, 'vvsmM') ;        % in um/min, rs
load(fvfsmM, 'vfsmM') ;        % in um/min, rs
load(fvnsmM, 'vnsmM') ;        % in um/min
load(fv2dsmM, 'v2dsmM') ;      % in pix/min
load(fv2dsmMum, 'v2dsmMum') ;  % in scaled pix/min, proportional to um/min 

%% Plot each timepoint
tp2doA = 1:10:ntps ;
tp2doB = setdiff(1:ntps, tp2doA) ;
tp2do = [tp2doA, tp2doB] ;

% Predefine vector field subsampling for quiver
QS.getPIV() ;
xx = QS.piv.raw.x{1} ;
nYpts = 30 ;
if doubleCovered
    yspacing = round(size(xx, 2) / nYpts * 0.5 * QS.a_fixed) ;
else
    yspacing = round(size(xx, 2) / nYpts * QS.a_fixed) ;
end
if any(diff(xx(1, :)) > 0)
    sampleIDx = [] ;
    qsrow = 1:round(size(xx, 1) / nYpts):size(xx, 1) ;
    for qq = 0:yspacing:size(xx, 2) - 1
        nn = length(sampleIDx) ;
        sampleIDx(nn+1:nn + length(qsrow)) = qsrow + qq * size(xx, 1) ;
    end
elseif any(diff(QS.piv.raw.x{1}(:, 1)) < 0)
    error('handle the case of transposed xgrid here')
else
    error('Could not parse dimensions of piv xgrid')
end

for tidx = tp2do
    tp = timePoints(tidx) ;
    % dt = timePoints(i + 1) - tp ;
        
    close all
    % Load streamline positions at this timePoint
    plotOptions.XX = squeeze(pivPathlines.XX(tidx, :, :)) ;
    plotOptions.YY = squeeze(pivPathlines.YY(tidx, :, :)) ;
    plotOptions.fX = squeeze(facePathlines.fX(tidx, :, :)) ;
    plotOptions.fY = squeeze(facePathlines.fY(tidx, :, :)) ;
    plotOptions.vX = squeeze(vertexPathlines.vX(tidx, :, :)) ;
    plotOptions.vY = squeeze(vertexPathlines.vY(tidx, :, :)) ;
    % Load velocities to plot
    plotOptions.vsm = squeeze(vsmM(tidx, :, :))  ;
    plotOptions.vnsm = squeeze(vnsmM(tidx, :, :))  ;
    plotOptions.v2dsm = squeeze(v2dsmM(tidx, :, :)) ;
    plotOptions.v2dsmum = squeeze(v2dsmMum(tidx, :, :)) ;
    plotOptions.overwrite = overwrite ;
    plotOptions.sampleIDx = sampleIDx ; 
    QS.plotPathlineVelocitiesTimePoint(tp, plotOptions) ;
end

disp('done')


