function identifyFolds(QS, options)
%IDENTIFYFOLDS(QS, OPTIONS)
%   Identify fold locations & lobes in the spcutMesh of QS.
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   Current QuapSlap object
% options : struct with fields
%   overwrite : bool, default=false
%       overwrite results on disk
%   guess123 : length 3 increasing float array, default=[0.2, 0.5, 0.8]
%       guess for initial fractional AP positions of each fold
%       Values should all lie within (0, 1).
%   max_wander : bool, default=20
%       max amount that DV hoop id can wander per timepoint, in units of 
%       ap sampling steps
%   first_tp_allowed : int, default = -1
%       enforces that no folds are considered physical before
%       tp = first_tp_allowed
%   preivew : bool, default=QS.plotting.preview
%       display intermediate results
%
% Returns
% -------
% Saves fold locations to disk. 
%   Note that fold_onset is in units of QS.xp.fileMeta.timepoints, NOT 
%       indices into timepoints.
%   Data is saved here:
%       fullfile(lobeDir, ['fold_locations_sphi' dvexten '_avgpts.mat'])
%
% NPMitchell 2020

%% Unpack QS
lobeDir = QS.dir.lobe ;
nU = QS.nU ;
nV = QS.nV ;
t0 = QS.t0set() ;
timeInterval = QS.timeInterval ;
uvexten = QS.uvexten ;  % string like sprintf('_nU%04d_nV%04d', nU, nV) ;
timePoints = QS.xp.fileMeta.timePoints ;
preview = QS.plotting.preview ;
spcutMeshBase = QS.fullFileBase.spcutMesh ;

%% Unpack options
if nargin < 2
    options = struct() ;
end
if isfield(options, 'guess123')
    guess123 = options.guess123 ;
    if length(guess123) ~= 3 || ~all(diff(guess123) > 0)
        disp(guess123)
        disp(['length(guess123) = ', num2str(length(guess123))])
        disp('diff(guess123) = ')
        disp(diff(guess123))
        error('guess123 must be length 3 increasing float array')
    end
else
    guess123 = [0.2, 0.5, 0.8] ;
end
if isfield(options, 'max_wander')
    max_wander = options.max_wander ;
else
    max_wander = 20 ; 
end
if isfield(options, 'maxDistsFromGuess')
    maxDistsFromGuess = options.maxDistsFromGuess ;
else
    maxDistsFromGuess = 0.1 * [1,1,1] ; 
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
if isfield(options, 'first_tp_allowed')
    first_tp_allowed = options.first_tp_allowed ;
else
    first_tp_allowed = -1 ;
end
% overwrite QS preview option if present in local options
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'meshCoords')
    switch lower(options.meshCoords)
        case 'sp'
            spcutMeshBase = QS.fullFileBase.spcutMesh ;
        case 'sp_sm'
            spcutMeshBase = QS.fullFileBase.spcutMeshSm ;
        otherwise 
            error('Could not recognize options.meshCoords: sp, sm_sm, etc')
    end
end

%% First compute using the avgpts (DVhoop means)
disp('Identifying lobes...')
foldfn = fullfile(lobeDir, ['fold_locations_sphi' uvexten '_avgpts.mat']) ;
if exist(foldfn, 'file') && ~overwrite
    disp('Loading lobes')
    % Save the fold locations as a mat file
    load(foldfn, 'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset', ...
        'rssfold', 'rssfold_frac', 'rssmax', 'rmax')
    
else
    options.guess123 = guess123 ;
    options.max_wander = max_wander ;
    options.preview = false ;
    options.method = 'avgpts' ;  % whether to use avgpts_ss or ringpath_ss
    options.first_tp_allowed = first_tp_allowed ;
    options.maxDistsFromGuess = maxDistsFromGuess ;
    [folds, ssfold, ssfold_frac, ssmax, rmax, fold_onset] = ...
        identifyLobes(timePoints, spcutMeshBase, options) ;
    
    % Compute ringpath pathlength for results found using centerline
    disp('Converting folds to ringpath_ss locations...')
    [rssfold, rssfold_frac, rssmax] = rssFromFoldID(folds, ...
        timePoints, spcutMeshBase) ;
    
    readme = struct(...
        'timeInterval', 'float, dt for timestamps in QS.timeUnits', ...
        'fold_onset', '#folds x 1 float, timestamps (not indices) of fold onset', ...
        'folds', '#timepoints x #folds int, indices of nU sampling of folds', ...
        'ssmax', '#timepoints x 1 float, maximum length of the centerline', ...
        'rssmax', '#timepoints x 1 float, maximum proper length of the surface over time', ...
        'ssfold', '#timepoints x #folds float, positional pathlength along centerline of folds', ...
        'rssfold', '#timepoints x #folds float, positional proper length along surface of folds', ...
        'rmax', '#timepoints x 1 float, maximum radius of any part of the surface at each timepoint') ;

    % Save the fold locations as a mat file
    timeInterval = QS.timeInterval ;
    save(foldfn, 'rssfold', 'rssfold_frac', 'rssmax', 'rmax', ...
        'timeInterval', 'ssfold', 'folds', 'ssfold_frac', ...
        'ssmax', 'fold_onset', 'readme')
end
clearvars guess123 maxwander

% Plot results as both avgpts and ringpath distances
fold_ofn = dir(fullfile(lobeDir, ['radii_folds' uvexten '_avgpts*.png'])) ;
if (length(fold_ofn) < length(timePoints)) || overwrite || true
    disp('Plotting ss folds...')
    aux_plot_folds(folds, ssfold, ssfold_frac, ssmax, rmax, nU, ...
        timePoints, lobeDir, uvexten, spcutMeshBase, ...
        'avgpts', overwrite, timeInterval, t0)
end

% Plot the radii and the fold locations
fold_ofn = dir(fullfile(lobeDir, ['radii_folds' uvexten '_avgpts*.png'])) ;
if (length(fold_ofn) < length(timePoints)) || overwrite  || true
    disp('Plotting rss folds...')
    aux_plot_folds(folds, rssfold, rssfold_frac, rssmax, rmax, nU, ...
        timePoints, lobeDir, uvexten, spcutMeshBase, ...
        'ringpath', overwrite, timeInterval, t0)
end

% Done creating folds and images