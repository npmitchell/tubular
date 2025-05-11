function [folds, ssfold, ssfold_frac, ssmax, rmax, fold_onset] = ...
    identifyLobes(timePoints, sphiBase, options)
%idenitfyLobes(timePoints, sphiBase, guess1, guess2, guess3, preview) Find folds in meshes 
%   Load each spcutMesh, find the local minima in radius, and mark these as 
%   fold locations. Track the location of those local minima over time, and
%   mark fold locations before the appearance of a fold as the location
%   where it appears. The pathlength is computed from the mean centerline's
%   pathlength ('mcline'), not ringpath_ss, which uses the surface
%   pathlength.
%
% Parameters
% ----------
% timePoints : N x 1 float/int array
%   The timestamps of each file to load (1 per minute is assumed)
% sphiBase : string
% guess123 : array of #folds floats, each between 0 and 1
%   Guesses for fractional position along U of each fold
% maxDistsFromGuess : #folds x 1 numeric array
%   how far from the initial guess to consider a detected minimum in radius 
%   as a valid match, in units of total AP discretization length (ie in 
%   units of #pts along AP axis in gridded parameterization)
% max_wander : #folds x 1 numeric array
%   maximum distance a fold location can wander in a given timepoint once 
%   identified, in units of pathlength (ss)
% wander_units : str ('pcAP' or QS.spaceUnits ['$u$m', for ex])
%   The units of wandering 
% preview : bool
% method : ('ringpath' or 'avgpts')
%   which method to use to find local minima, 'avgpts' recommended
% first_tp_allowed : int, default=-1
%   First timepoint in which we allow detected local minumum to be 
%   considered a true, physical fold. There may be one value for each fold
%   passed to the function, for ex [-1 -1 10] would prevent a third fold
%   from being found before t=10. 
%
% Returns
% -------
% fold_onset : #folds x 1 float
%   timestamps (not indices) of fold onset
% folds : #timepoints x #folds int
%   indices of nU sampling of folds
% ssmax : #timepoints x 1 float
%   maximum length of the centerline for each timepoint
% ssfold : #timepoints x #folds float
%   positional pathlength along centerline of folds
% rmax : #timepoints x 1 float
%   maximum radius of any part of the surface at each timepoint
% 
% NPMitchell 2019

% Default options
guess123 = [0.3, 0.5, 0.8] ;
nfolds = length(guess123) ;
max_wander = 20 ;
preview = false ;
verbose = true ;
method = 'avgpts' ;  % whether to use avgpts_ss or ringpath_ss
first_tp_allowed = 0 ;
maxDistsFromGuess = 0.1 * ones(nfolds, 1) ;  % maximum allowed distance from initial guesses to allow minima
wander_units = 'pcAP' ;  % 'pcAP' for percent of ap length, or 'spaceUnits' ;

if isfield(options, 'guess123')
    guess123 = options.guess123 ;
end
if isfield(options, 'max_wander')
    max_wander = options.max_wander ;
end
if isfield(options, 'wander_units')
    wander_units = options.wander_units ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'verbose')
    preview = options.verbose ;
end
if isfield(options, 'method')
    method = options.method ;
end
if isfield(options, 'first_tp_allowed')
    first_tp_allowed = options.first_tp_allowed ;
end
if isfield(options, 'maxDistsFromGuess')
    maxDistsFromGuess = options.maxDistsFromGuess ;
end

% Default behavior: if single threshold time is enforced, enforce for all
% folds.
if length(first_tp_allowed) == 1
    first_tp_allowed = first_tp_allowed * ones(nfolds, 1) ;
end

% pathlength preallocations
ssfold = zeros(length(timePoints) + 1, nfolds) ;
ssfold_frac = zeros(length(timePoints) + 1, nfolds) ;
ssmax = zeros(length(timePoints), 1) ;
    
% Convert method to a boolean
if strcmp(method, 'avgpts')
    method = true ;
elseif strcmp(method, 'ringpath')
    method = false ;
else
    error('method for assigning ss must be avgpts or ringpath')
end

% Store maximum radius if desired as output
if nargout > 4
    rmax = zeros(length(timePoints), 1) ;
end

% Fold onset timepoint indices
firsts = true(nfolds, 1) ;
onset_kks = Inf * ones(nfolds, 1) ;

% Consider each timepoint
if length(timePoints) == 1
    disp('Static dataset (only one timepoint). Finding likely folds.')

    % Translate to which timestamp
    t = timePoints ;
    load(sprintf(sphiBase, t), 'spcutMesh') ;

    % Extract centerline pathlength from DVhoop mean curve
    if method
        ss = spcutMesh.avgpts_ss ;
    else
        ss = spcutMesh.ringpath_ss ;
    end

    maxss = max(ss) ;
    ssmax = maxss ;

    % Initialize the fold positions to be approximately .3, .5, .8
    folds = int16(ones(length(timePoints) + 1, nfolds)) ;
    foldguess = guess123 * length(spcutMesh.phi0s) ;
    for qq = 1:nfolds
        folds(:, qq) = foldguess(qq) .* folds(:, qq) ;            
    end
    maxDists = length(spcutMesh.phi0s) * maxDistsFromGuess ;

    % Find the minima in the radius. First make radius 1d
    rad = mean(spcutMesh.radii_from_mean_uniform_rs, 2) ;
    minidx = islocalmin(rad) ;
    if any(minidx) && any(t > (first_tp_allowed - 1))
        minidx = find(minidx) ;
        % Identify each minimum which fold it may be
        dd = zeros(length(minidx), 1) ;
        which_fold = zeros(length(minidx), 1) ;
        for jj = 1:length(minidx)
            [dd(jj), which_fold(jj)] = min(abs(minidx(jj) - folds(1, :))) ;
        end

        % Mark if this is the first appearance of a fold
        % Also make sure this is close enough to the guess to be allowed
        % Note: folds(kk, :) are the guess fold locations
        for qq = 1:length(firsts)
            if firsts(qq) && t > (first_tp_allowed(qq) - 1) && ...
                    any(which_fold == qq) && ...
                    any(dd(which_fold == qq) < maxDists(qq))
                firsts(qq) = false ;
                onset_kks(qq) = 1 ;
            end
        end
        % if any(t > first_tp_allowed)
        %     disp('pausing')
        % end

        % Sort indices to folds
        if length(minidx) > 3
            disp('there are more minima than folds. Sorting...')
            % Find closest fold candidate for each fold (1,2,3)
            for pp = 1:size(folds, 2)
                if length(find(which_fold == pp)) > 1
                    disp(['more than one fold #' num2str(pp) ' found. Choosing closest within range....'])
                    [~, closest_idx] = min(dd(which_fold == pp)) ;
                    options = find(which_fold ==pp) ;
                    choice1 = options(closest_idx) ;
                    to_remove = setdiff(options, choice1) ;
                    to_keep = setdiff(1:length(minidx), to_remove) ;
                    % Remove the further matches
                    minidx = minidx(to_keep) ;
                    which_fold = which_fold(to_keep) ;
                    dd = dd(to_keep) ;
                end
            end
        end
        folds_kk = folds(1, :) ;
        folds_kk(which_fold) = minidx ;

        % Ensure that no folds have wandered farther than max_wander unless
        % this is the first appearance of that fold
        for pp = 1:length(firsts)
            % If the fold has already appeared, consider its motion, and 
            % check if in bounds

            if ~firsts(pp)
                % fold has appeared previously, check that this timepoint
                % is after the onset of fold formation 
                if 1 > onset_kks(pp)
                    % not first appearance, so check distance
                    if strcmpi(wander_units, 'pcap')
                        if t == timePoints(end)
                            pause(2)
                        end
                        wander_too_far = abs(double(folds_kk(pp) - folds(kk, pp))) > max_wander ; 
                    elseif strcmpi(wander_units, 'spaceunits')
                        wander_too_far = abs(ss(folds_kk(pp)) - ssfold(kk, pp)) > max_wander ;
                    else
                        error(['Unknown unit for wander_units: should be pcAP (for %AP axis) or ' QS.spaceUnits])
                    end
                    if wander_too_far
                        % too far, mark as unchanged
                        disp('outside wandering range')
                        if ~any(firsts)
                            disp('pausing here')
                        end
                        folds_kk(pp) = folds(1, pp) ;
                    end
                end
            end
        end

        if verbose
            disp(['folds = ', num2str(folds(2, :))])
        end
    end

    % Append the updated fold positions if valid
    % First check if first fold position or not. 
    for qq = 1:nfolds
        if ~firsts(qq)
            folds(2, qq) = folds_kk(qq) ;
            ssfold(2, qq) = ss(folds_kk(qq)) ;
            ssfold_frac(2, qq) = ss(folds_kk(qq)) / maxss ; 
        end
    end

    if true
        clf
        set(gcf, 'visible', 'on')
        plot(rad)
        hold on;
        plot(folds(2, :), rad(folds(2, :)), 'o')
        xlabel('s/L')
        ylabel('radius')
        title(['fold id, t=', num2str(t)])
        pause(0.001)
    end

    % Store maximum radius if desired as output
    if nargout > 4
        rmax = max(spcutMesh.radii_from_mean_uniform_rs(:)) ;
    end

    % Plot the indices of the identified folds
    if preview
        clf
        plot(folds(:, 1)); hold on;
        plot(folds(:, 2)); hold on;
        plot(folds(:, 3)); hold on;
        title('Fold locations')
        xlabel('time [min]')
        pause(0.01)
    end

    % Truncate extra first entry which was the initial guess
    folds = folds(2:end, :) ;
    ssfold = ssfold(2:end, :) ;
    ssfold_frac = ssfold_frac(2:end, :) ;

    % Go back and find the ss (pathlength) value for pre-fold indices
    if method
        ss = spcutMesh.avgpts_ss ;
    else
        ss = spcutMesh.ringpath_ss ;
    end

    % Assign the correct ss value to each
    for pp = 1:length(onset_kks)
        if qq < onset_kks(pp)
            assert(ss(folds(qq, 1)) > 0)
            ssfold(qq, pp) = ss(folds(qq, pp)) ;
            assert(ssfold(qq, pp) ~= 0)
        end
    end
else
    for kk = 1:length(timePoints)
        % Translate to which timestamp
        t = timePoints(kk) ;
        load(sprintf(sphiBase, t), 'spcutMesh') ;

        % Extract centerline pathlength from DVhoop mean curve
        if method
            ss = spcutMesh.avgpts_ss ;
        else
            ss = spcutMesh.ringpath_ss ;
        end

        maxss = max(ss) ;
        ssmax(kk, :) = maxss ;

        if kk == 1
            % Initialize the fold positions to be approximately .3, .5, .8
            folds = int16(ones(length(timePoints) + 1, nfolds)) ;
            foldguess = guess123 * length(spcutMesh.phi0s) ;
            for qq = 1:nfolds
                folds(:, qq) = foldguess(qq) .* folds(:, qq) ;            
            end
        end
        maxDists = length(spcutMesh.phi0s) * maxDistsFromGuess ;

        % Find the minima in the radius. First make radius 1d
        rad = mean(spcutMesh.radii_from_mean_uniform_rs, 2) ;
        minidx = islocalmin(rad) ;
        if any(minidx) && any(t > (first_tp_allowed - 1))
            minidx = find(minidx) ;
            % Identify each minimum which fold it may be
            dd = zeros(length(minidx), 1) ;
            which_fold = zeros(length(minidx), 1) ;
            for jj = 1:length(minidx)
                [dd(jj), which_fold(jj)] = min(abs(minidx(jj) - folds(kk, :))) ;
            end

            % If there is an increase of two folds, name them separately
            % if (length(minidx) - nnz(folds)) > 1 && nnz(folds) < 2
            %    % There are now at least two more folds so two are assigned to 
            %    % the same fold index: assigned to nnz(folds) + 1.
            %    % If there are more than two more new assignments, this is no
            %    % good
            % ALTERNATIVE:
            % We instead say 0.3 is fold 1, 0.5 is fold 2, 0.8 is fold 3

            % Mark if this is the first appearance of a fold
            % Also make sure this is close enough to the guess to be allowed
            % Note: folds(kk, :) are the guess fold locations
            for qq = 1:length(firsts)
                if firsts(qq) && t > (first_tp_allowed(qq) - 1) && ...
                        any(which_fold == qq) && ...
                        any(dd(which_fold == qq) < maxDists(qq))
                    firsts(qq) = false ;
                    onset_kks(qq) = kk ;
                end
            end
            % if any(t > first_tp_allowed)
            %     disp('pausing')
            % end

            % Sort indices to folds
            if length(minidx) > 3
                disp('there are more minima than folds. Sorting...')
                % Find closest fold candidate for each fold (1,2,3)
                for pp = 1:size(folds, 2)
                    if length(find(which_fold == pp)) > 1
                        disp(['more than one fold #' num2str(pp) ' found. Choosing closest within range....'])
                        [~, closest_idx] = min(dd(which_fold == pp)) ;
                        options = find(which_fold ==pp) ;
                        choice1 = options(closest_idx) ;
                        to_remove = setdiff(options, choice1) ;
                        to_keep = setdiff(1:length(minidx), to_remove) ;
                        % Remove the further matches
                        minidx = minidx(to_keep) ;
                        which_fold = which_fold(to_keep) ;
                        dd = dd(to_keep) ;
                    end
                end
            end
            folds_kk = folds(kk, :) ;
            folds_kk(which_fold) = minidx ;

            % Ensure that no folds have wandered farther than max_wander unless
            % this is the first appearance of that fold
            for pp = 1:length(firsts)
                % If the fold has already appeared, consider its motion, and 
                % check if in bounds

                if ~firsts(pp)
                    % fold has appeared previously, check that this timepoint
                    % is after the onset of fold formation 
                    if kk > onset_kks(pp)
                        % not first appearance, so check distance
                        if strcmpi(wander_units, 'pcap')
                            if t == timePoints(end)
                                pause(2)
                            end
                            wander_too_far = abs(double(folds_kk(pp) - folds(kk, pp))) > max_wander ; 
                        elseif strcmpi(wander_units, 'spaceunits')
                            wander_too_far = abs(ss(folds_kk(pp)) - ssfold(kk, pp)) > max_wander ;
                        else
                            error(['Unknown unit for wander_units: should be pcAP (for %AP axis) or ' QS.spaceUnits])
                        end
                        if wander_too_far
                            % too far, mark as unchanged
                            disp('outside wandering range')
                            if ~any(firsts)
                                disp('pausing here')
                            end
                            folds_kk(pp) = folds(kk, pp) ;
                        end
                    end
                end
            end

            if verbose
                disp(['folds = ', num2str(folds(kk, :))])
            end
        end

        % Append the updated fold positions if valid
        % First check if first fold position or not. 
        for qq = 1:nfolds
            if ~firsts(qq)
                folds(kk + 1, qq) = folds_kk(qq) ;
                ssfold(kk + 1, qq) = ss(folds_kk(qq)) ;
                ssfold_frac(kk + 1, qq) = ss(folds_kk(qq)) / maxss ; 
            end
        end

        if true
            clf
            set(gcf, 'visible', 'on')
            plot(rad)
            hold on;
            plot(folds(kk+1, :), rad(folds(kk+1, :)), 'o')
            xlabel('s/L')
            ylabel('radius')
            title(['fold id, t=', num2str(t)])
            pause(0.001)
        end

        % Store maximum radius if desired as output
        if nargout > 4
            rmax(kk) = max(spcutMesh.radii_from_mean_uniform_rs(:)) ;
        end

        % Plot the indices of the identified folds
        if preview
            clf
            plot(folds(:, 1)); hold on;
            plot(folds(:, 2)); hold on;
            plot(folds(:, 3)); hold on;
            title('Fold locations')
            xlabel('time [min]')
            pause(0.0001)
        end
    end
    % Truncate extra first entry which was the initial guess
    folds = folds(2:end, :) ;
    ssfold = ssfold(2:end, :) ;
    ssfold_frac = ssfold_frac(2:end, :) ;

    % Identify location prior to appearance of each fold as first location
    for qq = 1:length(onset_kks)
        folds(1:(onset_kks(qq)-1), qq) = folds(onset_kks(qq), qq) ;
        ssfold_frac(1:(onset_kks(qq)-1), qq) = ssfold_frac(onset_kks(qq), qq) ;
    end

    % Go back and find the ss (pathlength) value for pre-fold indices
    for qq = 1:max(onset_kks)
        t = timePoints(qq) ;
        load(sprintf(sphiBase, t), 'spcutMesh') ;
        if method
            ss = spcutMesh.avgpts_ss ;
        else
            ss = spcutMesh.ringpath_ss ;
        end

        % Assign the correct ss value to each
        for pp = 1:length(onset_kks)
            if qq < onset_kks(pp)
                assert(ss(folds(qq, 1)) > 0)
                ssfold(qq, pp) = ss(folds(qq, pp)) ;
                assert(ssfold(qq, pp) ~= 0)
            end
        end

        % if qq < k2 
        %     assert(ss(folds(qq, 2)) > 0)
        %     ssfold(qq, 2) = ss(folds(qq, 2)) ;
        %     assert(ssfold(qq, 2) ~= 0)
        % end
        % if qq < k3 
        %     assert(ss(folds(qq, 3)) > 0)
        %     ssfold(qq, 3) = ss(folds(qq, 3)) ;
        %     assert(ssfold(qq, 3) ~= 0)
        % end
    end
end

fold_onset = timePoints(onset_kks) ;

end
