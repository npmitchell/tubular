function cntrlines = cleanFastMarchingCenterlines(tubi, idOptions)
% cntrlines = cleanFastMarchingCenterlines(QS, idOptions)
% Identify anomalous centerlines in a time series, fix them and save.
%
% Parameters
% ----------
% tubi : TubULAR class instance
% idOptions : struct with fields
%   - ssr_thres : float
%       threshold sum of squared residuals between adjacent timepoint 
%       centerlines to consider the later timepoint centerline to be
%       anomalous & require cleaning
%   - overwrite : bool
%       overwrite the centerlines on disk
%
% Returns
% -------
% QS.cleanCntrlines : length(timePoints) x 1 cell
%   populated with cntrlines in XYZ coords (data space, pixels)
% QS.fileName.cleanCntrlines
%
% NPMitchell 2020

% Unpack options
overwrite = false ;
ssr_thres = 15 ;
if nargin < 2
    idOptions = struct() ;
end
if isfield(idOptions, 'overwrite')
    overwrite = idOptions.overwrite ;
end
if isfield(idOptions, 'ssr_thres')
    ssr_thres = idOptions.ssr_thres ;
end

% Unpack QS
disp('Searching for anomalous centerlines')
timePoints = tubi.xp.fileMeta.timePoints ;
[xyzlim_raw, ~, xyzlim_um] = tubi.getXYZLims();
cleanclinefn = tubi.fileName.cleanFMCenterlines ;
centerlineBase = tubi.fullFileBase.centerlineXYZ ;
[rot, trans] = tubi.getRotTrans() ;
resolution = tubi.APDV.resolution ;
centerlineDir = tubi.dir.cntrline ;

if exist(cleanclinefn, 'file') && ~overwrite
    disp('Anomalous centerline identification on disk. Loading...')
    load(cleanclinefn, 'cntrlines')
else
    ssrs = zeros(length(timePoints) - 1, 1) ;
    cntrlines = cell(length(timePoints), 1) ;
    cntrlines_rs = cell(length(timePoints), 1) ;
    anomalous = false(length(timePoints), 1) ;
    for kk = 1:length(timePoints)
        tt = timePoints(kk) ;
        tidx = tubi.xp.tIdx(tt) ;
        disp(['t = ' num2str(tt)])
        
        % load the centerline for this timepoint
        clfns = dir(replace(centerlineBase, ...
            tubi.timeStampStringSpec, ...
            num2str(tt, tubi.timeStampStringSpec))) ;
        cline = dlmread(fullfile(clfns(1).folder, clfns(1).name)) ;
        cntrlines{tidx} = cat(2, ss_from_xyz(cline), cline) ; 
        tmp = ((rot * cline')' + trans) * resolution ; 
        cntrlines_rs{tidx} = cat(2, ss_from_xyz(tmp), tmp) ;

        % compare to previous
        if tt > timePoints(1)
            idx = pointMatch(cline, prevcline) ;
            ssr = sum(vecnorm(cline - prevcline(idx, :), 2, 2)) / length(idx) ;
            anomal = ssr > ssr_thres ;
            % Declare that a bad centerline was found
            if anomal
                disp(['Found t= ' num2str(tt) ' to be anomalous'])
            end
        else
            ssr = 0 ;
            anomal = false ;
        end

        % rename as previous
        if ~anomal
            prevcline = cline ;
        end

        % Store for a recap
        ssrs(tidx) = ssr * resolution ;
        anomalous(tidx) = anomal ;
    end
    clearvars kk cline

    % Now adjust the anomalous ones to fix them up
    ssrfix = zeros(length(anomalous), 1) ;
    if any(anomalous)
        badID = find(anomalous) ;
        for qqID = 1:length(find(anomalous))
            qq = badID(qqID) ;
            % Find the last timepoint that was NOT anomalous and precedes 
            % the current one
            previd = find(~anomalous & (timePoints' < timePoints(qq)), 1, 'last') ;
            prevtp = timePoints(previd) ;
            % Find the next timepoint that was NOT anomalous and follows
            % the current one
            nextid = find(~anomalous & (timePoints' > timePoints(qq)), 1, 'first') ;
            cl1 = cntrlines{timePoints(previd)} ;
            cl2 = cntrlines{timePoints(nextid)} ;
            if size(cl1, 2) == 4
                cl1 = cl1(:, 2:4) ;
                cl2 = cl2(:, 2:4) ;
            end
            id = pointMatch(cl1, cl2) ;
            cline = (cl2(id, :) + cl1) * 0.5 ;
            cntrlines{timePoints(qq)} = cat(2, ss_from_xyz(cline), cline)  ;

            % Do the same for rs centerline
            cl1 = cntrlines_rs{timePoints(previd)} ;
            cl2 = cntrlines_rs{timePoints(nextid)} ;
            id = pointMatch(cl1(:, 2:4), cl2(:, 2:4)) ;
            cline = (cl2(id, 2:4) + cl1(:, 2:4)) * 0.5 ;
            cntrlines_rs{timePoints(qq)} = cat(2, ss_from_xyz(cline), cline) ;

            % Prepare to plot the fixed SSR    
            % load the centerline for prev timepoint
            clfns = dir(replace(centerlineBase, tubi.timeStampStringSpec, ...
                num2str(prevtp, tubi.timeStampStringSpec))) ;
            prevcline = dlmread(fullfile(clfns(1).folder, clfns(1).name)) ;
            % Convert to APDV [um]
            prevcline = ((rot * prevcline')' + trans) * resolution ; 
            % Find the SSR
            idx = pointMatch(cline, prevcline) ;
            ssrfix(qq) = sum(vecnorm(cline - prevcline(idx, :), 2, 2)) / length(idx) ;
        end
        clearvars cl1 cl2 qq
    else
        disp('all centerlines pass quality checks')
    end
    
    % Plot the anomalous determination
    close all; figure('visible', 'off')
    plot(timePoints, ssrs, '-')
    hold on;
    plot(timePoints(anomalous), ssrs(anomalous), 'o')
    plot(timePoints(anomalous), ssrfix(anomalous), 's')
    legend({'S.S.R.', 'anomalous', 'adjusted'})
    xlabel('time [min]')
    ylabel('Sum of squared residuals [\mum]')
    if any(anomalous)
        title('Identification of anomalous centerlines')
    else
        title('No anomalous centerlines found')
    end
    saveas(gcf, fullfile(centerlineDir, 'centerlines_anomalies_fixed.png'))
    close all

    % Save the corrected centerlines as a cell array
    save(cleanclinefn, 'cntrlines', 'cntrlines_rs') 

    % View the corrected centerlines
    for qid = 1:length(timePoints)
        qq = timePoints(qid) ;
        if qid > 1
            prevqq = timePoints(qid-1) ;
        else
            prevqq = qq ;
        end
        if qid < length(timePoints)
            nextqq = timePoints(qid+1) ;
        else
            nextqq = qq ;
        end
        first = tubi.xp.tIdx(max(min(timePoints), prevqq)) ;
        last = tubi.xp.tIdx(min(max(timePoints), nextqq)) ;
        for kk = first:last
            cl = cntrlines{kk} ;
            plot3(cl(:, 1), cl(:, 2), cl(:, 3))
            hold on;
        end
        hold off
        xlim(xyzlim_raw(1, :))
        ylim(xyzlim_raw(2, :))
        zlim(xyzlim_raw(3, :))
        title('Visualizing the corrected centerlines [pix]')
        pause(0.001)
    end
    clearvars cl
    pause(2)

    % View the corrected centerlines
    for qid = 1:length(timePoints)
        qq = timePoints(qid) ;
        if qid > 1
            prevqq = timePoints(qid-1) ;
        else
            prevqq = qq ;
        end
        if qid < length(timePoints)
            nextqq = timePoints(qid+1) ;
        else
            nextqq = qq ;
        end
        first = tubi.xp.tIdx(max(min(timePoints), prevqq)) ;
        last = tubi.xp.tIdx(min(max(timePoints), nextqq)) ;
        for kk = first:last
            cl = cntrlines_rs{kk} ;
            plot3(cl(:, 2), cl(:, 3), cl(:, 4))
            hold on;
        end
        hold off
        xlim(xyzlim_um(1, :))
        ylim(xyzlim_um(2, :))
        zlim(xyzlim_um(3, :))
        title('Visualizing the corrected RS centerlines [um]')
        pause(0.05)
    end
    clearvars cl
    close all
    clearvars ssrs ssr anomal anomalous
end

disp('Assignment to self')
tubi.cleanFMCenterlines = cntrlines ;

disp('done')
