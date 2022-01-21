function featureIDs = measureUVPrimePathlineFeatureIDs(QS, pathlineType, options)
%measureUVPrimePathlineFeatureIDs(QS, options)
% Interactively identify one longitudinal (zeta) position per feature for 
% Lagrangian pathlines using field measurements along Lagrangian pathlines
% of the coordinate system (u',v') [conformally mapped to plane]
%
% One can use, for example, radius of the pathlines on the mesh from the 
% centerline, or the normal velocity, or dz -- the ratio of projected 
% (embedding space) proper distance of a unit vector along longitudinal 
% mesh direction to the pullback distance along the longitudinal direction,
% or dp -- similar for circumferential direction, divv -- the divergence of
% the velocity field, etc. Uses two fields to identify the feature
% positions along the longitudinal dimension.
%
% todo: generalize beyond vP pathline array
%
% See also
% --------
% measurePathlineFeatureIDs.m
% 
% Parameters
% ----------
% QS : QuapSlap class object instance
% pathlineType : str specifier ('vertices', 'piv', 'faces')
%   whether pathlines in question thread through mesh vertices, piv
%   evaluation coordinates, or face barycenters at t=t0Pathlines
% options : optional struct with fields
%   overwrite : bool
%   field1 : str specifier ('radius', 'dz', 'dp', 'divv', 'veln')
%   field2 : str specifier ('radius', 'dz', 'dp', 'divv', 'veln')
% 
% Returns
% -------
% featureIDs : #features x 1 int array
%   longitudinal pullback coordinate for features, as an index into the
%   pullback pathlines (which are grid-like).
% 
% NPMitchell 2020

% Default options
overwrite = false ;
field1 = 'radius' ;
field2 = 'veln' ;
starttime_for_detection = 20 ;
endtime_for_detection = 60 ; 
if nargin < 2
    pathlineType = 'vertices' ;
end
if isempty(QS.pathlines.t0)
    QS.pathlines.t0 = QS.t0set() ;
end
t0Pathline = QS.pathlines.t0 ;
nU = QS.nU ;
climit = 0 ;
bwr256 = bluewhitered(256) ;
guess_minmax = 'min' ;          % whether initial guess is based on scalar field 1's local max or min values
% todo: handle double resolution

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'guess_minmax')
    guess_minmax = options.guess_minmax ;
end
if isfield(options, 'starttime_for_detection')
    starttime_for_detection = options.starttime_for_detection ;
end
if isfield(options, 'endtime_for_detection')
    endtime_for_detection = options.endtime_for_detection ;
end
if isfield(options, 'field1')
    field1 = options.field1 ;
end
if isfield(options, 'field2')
    field2 = options.field2 ;
end

%% Identify pathline coodinates using desired scalar fields for inspection
if strcmpi(pathlineType, 'vertices')
    %% To grab fold location in Lagrangian coords robustly, find minima of 
    % radii from ap average as kymograph and grab folds and lobes indexed in 
    % Lagrangian coords
    % fIDfn is the feature ID filename for these Lagrangian data
    fIDfn = sprintf(QS.fileName.pathlines_uvprime.featureIDs, t0Pathline) ;
    if exist(fIDfn, 'file') && ~overwrite
        load(fIDfn, 'featureIDs')
    else
        %% Interactively choose the feature locations

        % Load scalar field data to choose featureIDs (ex, folds)
        if strcmpi(field1, 'radius') || strcmpi(field2, 'radius')
            apKymoFn = ...
                sprintf(QS.fileName.pathlines_uvprime.kymographs.radius, ...
                t0Pathline) ;
            try
                tmp = load(apKymoFn, 'radius_apM') ;
                if strcmpi(field1, 'radius') 
                    sf1 = tmp.radius_apM ;
                end
                if strcmpi(field2, 'radius')
                    sf2 = tmp.radius_apM ;
                end    
            catch
                error('Run QS.measureUVPrimePathlines() before QS.measureUVPrimePathlineFeatureIDs()')
            end
        end
        if strcmpi(field1, 'divv') || strcmpi(field2, 'divv')
            error('handle using divergence field in uvprime coordSys here')
            try
                tmp = load(apKymoFn, 'divv_apM') ;
                if strcmpi(field1, 'divv') 
                    sf1 = tmp.divv_apM ;
                end
                if strcmpi(field2, 'divv')
                    sf2 = tmp.divv_apM ;
                end
            catch
                error('Run QS.measureUVPrimePathlines() before QS.measureUVPrimePathlineFeatureIDs()')
            end
        end
        if strcmpi(field1, 'veln') || strcmpi(field2, 'veln')
            error('handle using divergence field in uvprime coordSys here')
            try
                tmp = load(apKymoFn, 'veln_apM') ;
                if strcmpi(field1, 'veln') 
                    sf1 = tmp.veln_apM ;
                end
                if strcmpi(field2, 'veln')
                    sf2 = tmp.veln_apM ;
                end
            catch
                error('Run QS.plotPathlineMetricKinematics() before QS.plotPathlineStrainRate()')
            end
        end

        %% Interactively adjust feature locations
        close all
        subplot(1, 2, 1)
        imagesc(sf1)
        title(field1)
        subplot(1, 2, 2)
        imagesc(sf2)
        title(field2)
        sgtitle('Scalar fields for detecting UVPrime Lagrangian featureIDs')
        
        nfeatureIDs = input('How many featureIDs (ex folds) to identify in Lagrangian data? [Default=3]') ;
        if isempty(nfeatureIDs)
            nfeatureIDs = 3 ;
        end
        close all

        % Make a guess as to the features using minima of field1
        tps = QS.xp.fileMeta.timePoints - QS.t0 ;
        % grab data of the scalar field from some time after onset set by
        % start and end times for detection
        sf1d = mean(sf1(tps > max(starttime_for_detection, min(tps)) & ...
            tps < min(max(tps), endtime_for_detection), :), 1) ;
        if any(isnan(sf1d)) || isempty(sf1d)
            sf1d = mean(sf1, 1) ;
        end
        % smooth the 1d curve to find peaks/valleys
        div1dsm = savgol(sf1d, 2, 11) ;
        
        if strcmpi(guess_minmax, 'min')
            [~, featureIDs] = maxk(abs(islocalmin(div1dsm) .* (div1dsm - mean(div1dsm))), nfeatureIDs) ;
        else
            [~, featureIDs] = maxk(abs(islocalmax(div1dsm) .* (div1dsm - mean(div1dsm))), nfeatureIDs) ;
        end 
        featureIDs = sort(featureIDs) ;

        % Show guess overlaying field1 data
        figure ;
        set(gcf, 'visible', 'on')
        subplot(1, 2, 1)
        imagesc((1:nU)/nU, tps, sf1)
        % If there are both positive and negative values, use diverging
        % colormap
        if any(sf1(:) > 0) && any(sf1(:) < 0) 
            try
                colormap(bwr256)
            catch
                debugMsg(1, 'Could not make colormap bwr256')
            end
        end
        if climit > 0
            caxis([-climit, climit])
        end
        hold on;
        for qq = 1:nfeatureIDs
            plot(featureIDs(qq)/nU * ones(size(tps)), tps) ;
        end
        title(['Guess for featureIDs on ' field1])
        subplot(1, 2, 2) ;
        imagesc((1:nU)/nU, tps, sf2)
        % If there are both positive and negative values, use diverging
        % colormap
        if any(sf2(:) > 0) && any(sf2(:) < 0) 
            try
                colormap(bwr256)
            catch
                debugMsg(1, 'Cound not set colormap to bwr256')
            end
        end
        caxis([-max(abs(sf2(:))), max(abs(sf2(:)))])
        hold on;
        for qq = 1:nfeatureIDs
            plot(featureIDs(qq)/nU * ones(size(tps)), tps) ;
        end
        title('Guess for featureIDs on normal velocity')
        disp(['Guessed automatic featureIDs to be: [' num2str(featureIDs) ']'])

        % Update the guess
        for qq = 1:nfeatureIDs
            qok = false ;
            while ~qok
                msg = 'What is the Lagrangian zeta ID of feature ' ;
                msg = [msg num2str(qq) '? '] ;
                msg = [msg '[Default=' num2str(featureIDs(qq)) ']'] ;
                newvalley = input(msg) ;
                if isa(newvalley, 'double')
                    if ~isempty(newvalley)
                        featureIDs(qq) = newvalley ;
                    end
                end

                % Show guess overlaying div(v) data
                clf;
                set(gcf, 'visible', 'on')
                ax1 = subplot(1, 2, 1) ;
                imagesc((1:nU)/nU, tps, sf1)
                colormap(bwr256)
                if climit > 0
                    caxis([-climit, climit])
                else
                    caxis([-max(abs(sf1(:))), max(abs(sf1(:)))])
                end
                hold on;
                title('Guess for featureIDs on divergence(v)')
                ax2 = subplot(1, 2, 2) ;
                imagesc((1:nU)/nU, tps, sf2)
                colormap(bwr256)
                caxis([-max(abs(sf2(:))), max(abs(sf2(:)))])
                hold on;
                for pp = 1:nfeatureIDs
                    axes(ax1)
                    plot(featureIDs(pp)/nU * ones(size(tps)), tps) ;
                    axes(ax2)
                    plot(featureIDs(pp)/nU * ones(size(tps)), tps) ;
                end
                axes(ax2)
                title('Guess for featureIDs on normal velocity')
                disp(['Guessed automatic featureIDs to be: [' num2str(featureIDs) ']'])

                % Check it -- is the new feature look good?
                qYN = input(['does feature ' num2str(qq) ' look ok? [y/n]'], 's') ;
                if ~contains(lower(qYN), 'n')
                    qok = true ;
                end
            end
        end
        %% Save featureIDs (valleys of div(v))
        disp(['Saving featureIDs for UVPrime pathline coordinates: ' fIDfn])
        save(fIDfn, 'featureIDs')
    end
    
    % Store in QS
    QS.pathlines_uvprime.featureIDs.vertices = featureIDs ;
else
    error('Code for this pathlineType here')
end





