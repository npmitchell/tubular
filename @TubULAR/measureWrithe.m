function [Wr, Wr_density, dWr, Length_t, clines_resampled] = measureWrithe(tubi, options)
%MEASUREWRITHE(TUBI, OPTIONS)
%   Compute the centerline length and writhe and save to disk
%
% Parameters
% ----------
% tubi : TubULAR class instance
%   Current TubULAR object
% options : struct with fields
%   overwrite : bool, default=false
%       overwrite results on disk
%   overwriteFigs : bool, default=false
%       overwrite figure results on disk
%   preivew : bool, default=QS.plotting.preview
%       display intermediate results
%   Wr_style : str specifier, default='Levitt' ;
%       which style of writhe calculation to plot
%   black_figs : bool
%       plot figures with black background rather than white
%   flipy_centerline : bool
%       flip the centerline data y --> -y (inverts sign of writhe)
%
% Returns
% -------
% <none>
%
% NPMitchell 2020

%% Unpack tubi
nU = tubi.nU ;
nV = tubi.nV ;
dvexten = sprintf('_nU%04d_nV%04d', nU, nV) ;
timePoints = tubi.xp.fileMeta.timePoints ;
timeInterval = tubi.timeInterval ;
timeUnits = tubi.timeUnits ;
spaceUnits = tubi.spaceUnits ;
preview = tubi.plotting.preview ;
[rot, trans] = tubi.getRotTrans() ;
resolution = tubi.APDV.resolution ;
clineDVhoopBase = tubi.fullFileBase.clineDVhoop ;
cylinderMeshCleanBase = tubi.fullFileBase.cylinderMeshClean ;
writheDir = tubi.dir.writhe ;
writheImDir = fullfile(writheDir, 'images') ;
if ~exist(writheImDir, 'dir')
    mkdir(writheImDir)
end
meshDir = tubi.dir.mesh ;
% get timestamps to indicate on writhe plot
try
    load(tubi.fileName.fold, 'fold_onset') ;
catch 
    disp('no features to load')
    fold_onset = [] ;
end
[~, ~, xyzlim] = getXYZLims(tubi) ;
flipy = tubi.flipy ;
flipy_centerline = false ;
t0 = tubi.t0set() ;
if isempty(t0)
    t0 = 0 ;
end
timePoints = timePoints - t0 ; 

%% Unpack options
if nargin < 2
    options = struct() ;
end
% overwrite QS preview option if present in local options
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
if isfield(options, 'overwriteFigs')
    overwriteFigs = options.overwriteFigs ;
else
    overwriteFigs = false ;
end
if isfield(options, 'Wr_style')
    Wr_style = options.Wr_style ;
else
    Wr_style = 'Levitt' ;
end
if isfield(options, 'flipy_centerline')
    flipy_centerline = options.flipy_centerline ;
end
if isfield(options, 'omit_endpts')
    omit_endpts = options.omit_endpts ;
else
    omit_endpts = 4 ;
end
if isfield(options, 'filter_curve')
    filter_curve = options.filter_curve ;
else
    filter_curve = 7 ;
end
if isfield(options, 'black_figs')
    black_figs = options.black_figs ;
else
    black_figs = false ;
end

%% First compute Writhe using the avgpts (DVhoop mean positions)
wrfn = tubi.fileName.writhe ;
if ~exist(wrfn, 'file') || overwrite
    disp('Computing length and writhe...')

    [Wr, Wr_density, dWr, Length_t, clines_resampled] = ...
        aux_compute_writhe(clineDVhoopBase, timePoints, ...
        filter_curve, omit_endpts, flipy_centerline, preview) ;
    
    % Save the fold locations as a mat file
    save(wrfn, 'Wr', 'Wr_density', 'dWr', 'Length_t', 'clines_resampled')
else
    disp('Loading length and writhe...')
    load(wrfn, 'Wr', 'Wr_density', 'dWr', 'Length_t', 'clines_resampled')
end

% Which writhe to plot is determined by Wr_style
tmpfn = fullfile(writheImDir, ['writhe_' Wr_style '_vs_time_comparison_DVhoop.png']) ;
if ~exist(tmpfn, 'file') || overwrite || overwriteFigs
    % Compute ringpath pathlength for results found using centerline
    area_volume_fn = fullfile(meshDir, 'surfacearea_volume_stab.mat') ;
    aux_plot_writhe(tubi, clines_resampled, ...
        Wr, Wr_density, dWr, Length_t, writheImDir, area_volume_fn, ...
        fold_onset, Wr_style, clineDVhoopBase, ...
        cylinderMeshCleanBase, rot, trans, resolution, flipy, ...
        omit_endpts, black_figs)
end

% Done with measuring writhe
disp('done measuring writhe')