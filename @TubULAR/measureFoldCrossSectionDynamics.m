function measureFoldCrossSectionDynamics(QS, options)
% measureFoldCrossSectionDynamics(QS, options)
% 
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct
%   
% Returns
% -------
%
% NPMitchell 2021

%% Default options
probChannel = 2 ;
endTime = 1.5 ;  % end time in hrs for plotting, or in QS.timeUnits if QS.timeUnits ~= min
thres = 0.5 ;
timePoints = QS.xp.fileMeta.timePoints ;
features = QS.getFeatures('folds') ;
folds = features.folds ;
overwrite = false ;
preview = false ;

%% Unpack options
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'probabilityChannel')
    probChannel = options.probabilityChannel ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'endTime')
    endTime = options.endTime ;
end

nFolds = size(folds, 2) ;

%% Load or compute the features of the fold from mips (area)
outfn = fullfile(QS.dir.lobe, 'constriction_MIPs', 'initialInPoints.mat') ;
if ~exist(outfn, 'file') || overwrite
   for foldID = 1:nFolds       
        imfn = fullfile(QS.dir.lobe, 'constriction_MIPs', ...
                        sprintf('fold%03d', foldID), ...
                        [sprintf(QS.fileBase.name, timePoints(1)) ...
                        sprintf('_fold%03d_Probabilities.h5', foldID)]) ;
        dat = h5read(imfn, '/exported_data') ;
        dd = squeeze(dat(probChannel, :, :)) ;
        
        imagesc(dd')
        title('Click point in fold aperture')
        cntr = ginput(1) ;
        initialInPoint(foldID, :) = cntr ;
    end
    save(outfn, 'initialInPoint') ;
else
    load(outfn, 'initialInPoint') ;
end

outfn = fullfile(QS.dir.lobe, 'constriction_MIPs', 'foldMeasurementsViaMIPs.mat') ;
if ~exist(outfn, 'file') || overwrite
    xareas = zeros(length(timePoints), nFolds) ;
    rmajor = zeros(length(timePoints), nFolds) ;
    rminor = zeros(length(timePoints), nFolds) ;
    radii = zeros(length(timePoints), nFolds) ;
    cntrd = zeros(length(timePoints), nFolds, 2) ; 
    for foldID = 1:nFolds
        tidx = 1 ;
        for tp = timePoints
            disp(['t = ' num2str(tp)])
            imfn = fullfile(QS.dir.lobe, 'constriction_MIPs', ...
                            sprintf('fold%03d', foldID), ...
                            [sprintf(QS.fileBase.name, tp) ...
                            sprintf('_fold%03d_Probabilities.h5', foldID)]) ;
            dat = h5read(imfn, '/exported_data') ;
            dd = squeeze(dat(probChannel, :, :)) ;

            bw = bwconncomp(dd > thres) ;

            if tp == timePoints(1)
                cntr = initialInPoint(foldID, :) ;
                idx = round(cntr) ;
                cntr_idx = sub2ind(size(dd), idx(2), idx(1)) ;
            else
                try
                    cntr_idx = sub2ind(size(dd), prev_cntrd(2), prev_cntrd(1)) ;
                catch
                    cntr_idx = 1 ;
                end
            end
            if preview
                disp(['cntr_idx = ', num2str(cntr_idx)])
            end
            
            match = 0 ;
            for qq = 1:length(bw.PixelIdxList)
                if match == 0
                    if ~isempty(intersect(bw.PixelIdxList{qq}, cntr_idx))
                        match = qq ;
                    end
                end
            end
            try
                assert(match > 0) 
                halt = false ;
            catch
                disp('Centroid was not inside a fold region')
                clf
                imshow(dd > thres) ;
                hold on;
                plot(prev_cntrd(1), prev_cntrd(2), 'o')
                title('Centroid was not inside a fold region')
                pause(0.1)
                halt = true ;
            end

            % Store the result if it worked
            if ~halt
                xareas(tidx, foldID) = length(bw.PixelIdxList{match}) * ...
                    QS.APDV.resolution.^2 ;

                bwim = false(size(dd)) ;
                bwim(bw.PixelIdxList{match}) = true ;
                bw2 = bwconncomp(bwim) ;
                stats = regionprops(bw2,'Centroid',...
                    'MajorAxisLength','MinorAxisLength') ;
                radii(tidx, foldID) = mean([stats.MajorAxisLength ...
                    stats.MinorAxisLength],2) * 0.5 * QS.APDV.resolution;
                rmajor(tidx, foldID) = stats.MajorAxisLength * 0.5 * QS.APDV.resolution;
                rminor(tidx, foldID) = stats.MinorAxisLength * 0.5 * QS.APDV.resolution;
                cntrd(tidx, foldID, :) = stats.Centroid ;
                
                minPrev = max(1, tidx-5) ;
                if tidx == 1
                    prev_cntrd= round(...
                        squeeze(cntrd(minPrev:tidx, foldID, :))) ;
                else    
                    prev_cntrd= round(median(...
                        squeeze(cntrd(minPrev:tidx, foldID, :)), 1)) ;
                end
                
                % Show it
                if preview
                    clf
                    imshow(bwim')
                    hold on; 
                    plot(stats.Centroid(2), stats.Centroid(1), 'ro');
                    title(['t=', num2str(tp)])
                    pause(0.0001) ;
                end
            end

            % Update time dummy index
            tidx = tidx + 1 ;
        end
    end
    save(outfn, 'xareas', 'radii', 'rmajor', 'rminor', 'cntrd')
else
    load(outfn, 'xareas', 'radii', 'rmajor', 'rminor', 'cntrd')
end

%% Find filtered line to subtract
xareas_filt = gaussFilter1D(xareas, 4, 7) ;
radii_filt = gaussFilter1D(radii, 4, 7) ;
rmajor_filt = gaussFilter1D(rmajor, 4, 7) ;
rminor_filt = gaussFilter1D(rminor, 4, 7) ;

%% Prepare output and labels
outdir = fullfile(QS.dir.lobe, 'constriction_MIPs') ;

t0 = QS.t0set() ;
if contains(lower(QS.timeUnits), 'min')
    time2plot = (timePoints - t0) * QS.timeInterval / 60 ;
    timeLabel = ['time [hr]'] ;
else
    time2plot = (timePoints - t0) * QS.timeInterval ;
    timeLabel = ['time [' QS.timeUnits ']'] ;
end

%% Plot areas
close all
outFigFn = fullfile(outdir, 'foldCrossSectionAreaMIPs.pdf') ;
plot(time2plot, xareas, '.-')
xlabel(timeLabel, 'interpreter', 'latex')
ylabel(['fold cross-sectional radii [' QS.spaceUnits ']'], ...
    'interpreter', 'latex')
legend({'anterior', 'middle', 'posterior'}, 'interpreter', 'latex')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)

%% Plot radii
close all
outFigFn = fullfile(outdir, 'foldCrossSectionRadiiMIPs.pdf') ;
plot(time2plot, radii, '.-')
xlabel(timeLabel, 'interpreter', 'latex')
ylabel(['fold cross-sectional radii [' QS.spaceUnits ']'], ...
    'interpreter', 'latex')
legend({'anterior', 'middle', 'posterior'}, 'interpreter', 'latex')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)

%% Plot rmajor / rminor
colors = define_colors() ;
close all
outFigFn = fullfile(outdir, 'foldCrossSectionMajorMinorMIPs.pdf') ;
for foldID = 1:nFolds
    plot(time2plot, rmajor(:, foldID), '-', 'linewidth', 2, 'color', colors(foldID, :))
    hold on;
    plot(time2plot, rminor(:, foldID), '-', 'linewidth', 1, 'color', colors(foldID, :))
end
xlabel(timeLabel, 'interpreter', 'latex')
ylabel(['fold cross-sectional major/minor radii [' QS.spaceUnits ']'], ...
    'interpreter', 'latex')
legend({'anterior major radius', ...
    'anterior minor radius', ...
    'middle major radius', ...
    'middle minor radius', ...
    'posterior major radis', ...
    'posterior minor radis'}, 'interpreter', 'latex', 'location', 'eastOutside')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot areas -- DERIVS
close all
outFigFn = fullfile(outdir, 'foldCrossSectionAreaMIPsSub.pdf') ;
plot(time2plot(1:end-1), diff(xareas, 1, 1) / QS.timeInterval, '.-')
xlabel(timeLabel, 'interpreter', 'latex')
ylabel(['fold cross-sectional radii change [' QS.spaceUnits '/' QS.timeUnits ']'], ...
    'interpreter', 'latex')
legend({'anterior', 'middle', 'posterior'}, 'interpreter', 'latex')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)

%% Plot radii -- MEAN SUBTRACTED
close all
outFigFn = fullfile(outdir, 'foldCrossSectionRadiiMIPsSub.pdf') ;
plot(time2plot(1:end-1), diff(radii, 1, 1) / QS.timeInterval, '.-')
xlabel(timeLabel, 'interpreter', 'latex')
ylabel(['fold cross-sectional radii change [' QS.spaceUnits '/' QS.timeUnits ']'], ...
    'interpreter', 'latex')
legend({'anterior', 'middle', 'posterior'}, 'interpreter', 'latex')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)

%% Plot rmajor / rminor -- MEAN SUBTRACTED
colors = define_colors() ;
close all
outFigFn = fullfile(outdir, 'foldCrossSectionMajorMinorMIPsSub.pdf') ;
for foldID = 1:nFolds
    subplot(nFolds, 1, foldID)
    plot(time2plot(1:end-1), diff(rmajor(:, foldID), 1, 1), ...
        '-', 'linewidth', 2, 'color', colors(foldID, :))
    hold on;
    plot(time2plot(1:end-1), diff(rminor(:, foldID), 1, 1), ...
        '-', 'linewidth', 1, 'color', colors(foldID, :))
    xlim([min(time2plot), endTime])
    ylabel(['fold ' num2str(foldID)], 'interpreter', 'latex')
end
xlabel(timeLabel, 'interpreter', 'latex')
sgtitle(['fold cross-sectional major/minor radii change [' QS.spaceUnits '/' QS.timeUnits ']'], ...
    'interpreter', 'latex')
% legend({'anterior major radius', ...
%     'anterior minor radius', ...
%     'middle major radius', ...
%     'middle minor radius', ...
%     'posterior major radis', ...
%     'posterior minor radis'}, 'interpreter', 'latex', 'location', 'eastOutside')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot areas -- MEAN SUBTRACTED
close all
outFigFn = fullfile(outdir, 'foldCrossSectionAreaMIPsSub.pdf') ;
plot(time2plot, xareas - xareas_filt, '.-')
xlabel(timeLabel, 'interpreter', 'latex')
ylabel(['fold cross-sectional radii variation [' QS.spaceUnits ']'], ...
    'interpreter', 'latex')
legend({'anterior', 'middle', 'posterior'}, 'interpreter', 'latex')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)

%% Plot radii -- MEAN SUBTRACTED
close all
outFigFn = fullfile(outdir, 'foldCrossSectionRadiiMIPsSub.pdf') ;
plot(time2plot, radii - radii_filt, '.-')
xlabel(timeLabel, 'interpreter', 'latex')
ylabel(['fold cross-sectional radii variation [' QS.spaceUnits ']'], ...
    'interpreter', 'latex')
legend({'anterior', 'middle', 'posterior'}, 'interpreter', 'latex')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)

%% Plot rmajor / rminor -- MEAN SUBTRACTED
colors = define_colors() ;
close all
outFigFn = fullfile(outdir, 'foldCrossSectionMajorMinorMIPsSub.pdf') ;
for foldID = 1:nFolds
    subplot(nFolds, 1, foldID)
    plot(time2plot, rmajor(:, foldID) - rmajor_filt(:, foldID), ...
        '-', 'linewidth', 2, 'color', colors(foldID, :))
    hold on;
    plot(time2plot, rminor(:, foldID) - rminor_filt(:, foldID), ...
        '-', 'linewidth', 1, 'color', colors(foldID, :))
    xlim([min(time2plot), endTime])
    ylabel(['fold ' num2str(foldID)], 'interpreter', 'latex')
end
xlabel(timeLabel, 'interpreter', 'latex')
sgtitle(['fold cross-sectional major/minor radii variation [' QS.spaceUnits ']'], ...
    'interpreter', 'latex')
% legend({'anterior major radius', ...
%     'anterior minor radius', ...
%     'middle major radius', ...
%     'middle minor radius', ...
%     'posterior major radis', ...
%     'posterior minor radis'}, 'interpreter', 'latex', 'location', 'eastOutside')
xlim([min(time2plot), endTime])
saveas(gcf, outFigFn)

% %% Load or compute the features of the fold from meshes
% outfn = fullfile(QS.dir.lobe, 'constriction_MIPs', 'foldMeasurementsViaMeshes.mat') ;
% if ~exist(outfn, 'file') || overwrite
%     xareas = zeros(length(timePoints), nFolds) ;
%     rmajor = zeros(length(timePoints), nFolds) ;
%     rminor = zeros(length(timePoints), nFolds) ;
%     radii = zeros(length(timePoints), nFolds) ;
%     cntrd = zeros(length(timePoints), nFolds, 2) ; 
%     for foldID = 1:nFolds
%         tidx = 1 ;
%         for tp = timePoints
%             QS.setTime(tp) ;
%             
%             % Load the mesh for this timepoint (pre-smoothed)
%             QS.getCurrentSPCutMesh()
%             halt = false ;
%             
%             % Store the result if it worked
%             if ~halt
%                 xareas(tidx, foldID) = length(bw.PixelIdxList{match}) * ...
%                     QS.APDV.resolution.^2 ;
% 
%                 bwim = false(size(dd)) ;
%                 bwim(bw.PixelIdxList{match}) = true ;
%                 bw2 = bwconncomp(bwim) ;
%                 stats = regionprops(bw2,'Centroid',...
%                     'MajorAxisLength','MinorAxisLength') ;
%                 radii(tidx, foldID) = mean([stats.MajorAxisLength ...
%                     stats.MinorAxisLength],2);
%                 rmajor(tidx, foldID) = stats.MajorAxisLength ;
%                 rminor(tidx, foldID) = stats.MinorAxisLength ;
%                 cntrd(tidx, foldID, :) = stats.Centroid ;
%                 prev_cntrd= round(squeeze(cntrd(tidx, foldID, :))) ;
%             end
% 
%             % Update time dummy index
%             tidx = tidx + 1 ;
%         end
%     end
% end

