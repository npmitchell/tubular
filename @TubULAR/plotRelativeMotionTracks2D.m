function plotRelativeMotionTracks2D(QS, Options)
%plotRelativeMotionTracks2D(QS, Options)
%
% Plot the relative motion between two sets of tracks (in different layers
% of the tissue/surface, for ex) on top of the pullbacks used to measure
% those tracks. Summarize results of geodesic distances over time in 1d
% plots. 
%
% Parameters
% ----------
% 
% 
% NPMitchell 2021




relMotionFn = fullfile(QS.dir.tracking, 'relative_motion_tracks.mat') ;
overwrite = false ;
coordSys = 'spsm' ;
subdir1 = Options.subdir1 ;
subdir2 = Options.subdir2 ;
doubleCovered = false ;
blueColor = [0,    0.4470,    0.7410] ;
orangeColor = [0.8500,    0.3250,    0.0980] ;

if isfield(Options, 'overwrite')
    overwrite = Options.overwrite ;
end
if isfield(Options, 'coordSys')
    coordSys = Options.coordSys ;
end

if isfield(Options, 'doubleCovered')
    doubleCovered = Options.doubleCovered ;
end

figX = 300 ;
figY = 200 ;


%% Load data to plot
disp(['Loading track pairs from disk: ' relMotionFn])
load(relMotionFn, 'dusEuclidean', 'dusGeodesic', ...
    'tracks1', 'tracks2', 'pairIDs', 'U0s', 'V0s', ...
    'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
    'euclideanDistanceTraveledU', 'euclideanDistanceTraveledV', ...
    'v3d_u', 'v3d_v', 'nSaved') ;

n1 = length(tracks1) ;
n2 = length(tracks2) ;
nTracks = n2 ;
nTimePoints = size(tracks2{1}, 1) ;
timePoints = QS.xp.fileMeta.timePoints ;
assert(nTimePoints == length(timePoints)) ;


%% Plot in 2D on data
close all
fig = figure('Position', [100 100 1600 1000], 'units', 'centimeters');
colors = distinguishable_colors(nTracks) ;
t0 = QS.t0set();
timestamps = (timePoints - t0) * QS.timeInterval ;
if contains(lower(QS.timeUnits), 'min')
    timeunits = 'hr' ;
    timestamps = timestamps / 60 ;
end

% Plot tracks on each timepoint image
zooms = {'full', 'inset'} ;
Xlims = [NaN, NaN; 620, 950  ] ;
Ylims = [NaN, NaN; 550, 780 ] ;
sizes1 = [50, 200] ;
sizes2 = [150, 400] ;

linewidths = [1, 5];
fs = 18 ;

overlayStyles = {'simple', 'vector', 'relativeVector'} ;

for overlayID = 1:length(overlayStyles)
    extenOverlay = overlayStyles{overlayID} ;
    for zoomID = 1:length(zooms)
        exten = zooms{zoomID} ;
        Xlim = Xlims(zoomID, :) ;
        Ylim = Ylims(zoomID, :) ;
        sz1 = sizes1(zoomID) ;
        sz2 = sizes2(zoomID) ;
        lw = linewidths(zoomID) ;

        % output directory for this zoom level
        outTrackImDir = fullfile(QS.dir.tracking, 'tracks', exten, extenOverlay) ;
        if ~exist(outTrackImDir, 'dir')
            mkdir(outTrackImDir)
        end


        % Plot each timepoint image
        for tidx = 1:nTimePoints
            % output file name
            outfn = fullfile(outTrackImDir, ...
                sprintf(['track_' exten '_' extenOverlay '_%06d.png'], timePoints(tidx))) ;

            if ~exist(outfn, 'file') || overwrite
                % obtain all paired tracks
                UU = zeros(nTracks, 2) ;
                VV = zeros(nTracks, 2) ;
                U0 = zeros(nTracks, 2) ;
                V0 = zeros(nTracks, 2) ;
                for ii = 1:nTracks
                    VV(ii, :) = tracks2{ii}(tidx, 1:2) ;
                    UU(ii, :) = tracks1{pairIDs(ii)}(tidx, 1:2) ;

                    % initial positions
                    V0(ii, :) = tracks2{ii}(1, 1:2) ;
                    U0(ii, :) = tracks1{pairIDs(ii)}(1, 1:2) ;
                end

                % obtain all other tracks
                otherTracks1 = zeros(n1, 2) ;
                for trackID = 1:n1
                    otherTracks1(trackID, :) = tracks1{trackID}(tidx, 1:2) ;
                end
                otherTracks2 = zeros(n2, 2) ;
                for trackID = 1:n2
                    otherTracks2(trackID, :) = tracks2{trackID}(tidx, 1:2) ;
                end

                % Obtain data image for this timepoint
                if strcmpi(coordSys, 'spsm') || strcmpi(coordSys, 'spsmrs')
                    im = imread(fullfile(QS.dir.im_sp_sm, subdir1, ...
                        sprintf(QS.fileBase.im_sp_sm, timePoints(tidx)))) ;
                    im2 = imread(fullfile(QS.dir.im_sp_sm, subdir2, ...
                        sprintf(QS.fileBase.im_sp_sm, timePoints(tidx)))) ;
                    assert(all(size(im) == size(im2)))
                end

                % Plot them for this timepoint
                clf
                % Axis 1
                blueGrayColor = [0.2, 0.4, 0.7] ;
                orangeColor = [0.85, 0.325, 0.098] ;
                ax0 = subtightplot(2, 2, 1) ;
                imshow(im);
                hold on;
                scatter(otherTracks1(:, 1), otherTracks1(:, 2), sz1, blueGrayColor, 'filled')
                scatter(UU(:, 1), UU(:, 2), sz1, colors, 'filled', 'markeredgecolor', 'k')
                title([subdir1 ' positions'], ...
                    'interpreter', 'Latex', 'fontsize', fs)
                if ~isnan(Xlim)
                    xlim(Xlim)
                    ylim(Ylim)
                end
                % Axis 2
                ax1 = subtightplot(2, 2, 2) ;
                imshow(im2);
                hold on;
                % scatter(otherTracks2(:, 1), otherTracks2(:, 2), sz*4, orangeColor, 's')
                scatter(VV(:, 1), VV(:, 2), sz2, colors, 's', 'linewidth', lw)
                title([subdir2 ' positions'], ...
                    'interpreter', 'Latex', 'fontsize', fs)
                if ~isnan(Xlim)
                    xlim(Xlim)
                    ylim(Ylim)
                end

                % OVERLAY AXIS ----------------------------------------
                % Axis 3
                ax2 = subtightplot(2, 1, 2) ;
                imshow(min(0.5*im+im2, 255));

                if zoomID > 1
                    set(ax2, 'Position', [0.05, 0.01, 0.9, 0.445])
                end

                timestamp = (timePoints(tidx)-t0) * QS.timeInterval ;
                hold on;
                if overlayID == 1
                    % scatter(otherTracks1(:, 1), otherTracks1(:, 2), sz, blueGrayColor, 'filled', 'markeredgecolor', 'k')
                    scatter(VV(:, 1), VV(:, 2), sz2, colors, 's', 'linewidth', lw)
                    scatter(UU(:, 1), UU(:, 2), sz1, colors, 'filled', 'markeredgecolor', 'k')
                    title(['Overlay, $t=' num2str(timestamp) '$ ' QS.timeUnits], ...
                        'interpreter', 'Latex', 'fontsize', fs)
                elseif overlayID == 2
                    dUV = VV - UU ;
                    % scatter(UU(:, 1), UU(:, 2), sz1, colors, 'filled', 'markeredgecolor', 'k')
                    quiver(UU(:, 1), UU(:, 2), dUV(:, 1), dUV(:, 2), 0, 'linewidth', 2, 'color', orangeColor, 'MaxHeadSize', 0.3)
                    title(['Relative position, $t=' num2str(timestamp) '$ ' QS.timeUnits], ...
                        'interpreter', 'Latex', 'fontsize', fs)
                else
                    dUV = (VV - V0s) - (UU - U0s) ;
                    % scatter(UU(:, 1), UU(:, 2), sz1, colors, 'filled', 'markeredgecolor', 'k')
                    quiver(UU(:, 1), UU(:, 2), dUV(:, 1), dUV(:, 2), 0, 'linewidth', 2, 'color', orangeColor, 'MaxHeadSize', 0.3)
                    title(['Relative displacement, $t=' num2str(timestamp) '$ ' QS.timeUnits], ...
                        'interpreter', 'Latex', 'fontsize', fs)
                end

                hold off;
                if ~isnan(Xlim)
                    xlim(Xlim)
                    ylim(Ylim)
                end

                disp(['Saving image of tracks: ' outfn])
                saveas(gcf, outfn)

            end
        end
    end
end


%% First plot raw relative distances, then displacements from
% starting positions

% Store in-plane positions -- Rough n Dirty method in UV as warmup for full Lagrangian treatment
Uall = nan(nTracks, nTimePoints, 2) ;
Vall = nan(nTracks, nTimePoints, 2) ;
phaseRaw = nan(nTracks, nTimePoints) ;
magUVRaw = nan(nTracks, nTimePoints) ;
phaseRelative = nan(nTracks, nTimePoints) ;
magUVRelative = nan(nTracks, nTimePoints) ;

if ~doubleCovered
    im = imread(sprintf(QS.fullFileBase.im_sp_sm, timePoints(1))) ;
else
    im = imread(sprintf(QS.fullFileBase.im_sp_sme, timePoints(1))) ;
end
umax = 1 ;
vmax = 1 ;
for ii = 1:nTracks
    % Get starting positions u0 and v0 for layer1 and 2
    Uall(ii, :, :) = tracks1{pairIDs(ii)}(:, 1:2) ;
    Vall(ii, :, :) = tracks2{ii}(:, 1:2) ;

    dY = Vall(ii, :, 2) - Uall(ii, :, 2) ;
    dX = Vall(ii, :, 1) - Uall(ii, :, 1) ;
    dXY = [dX(:), dY(:)] ;

    dUV = QS.dXY2duv(im, dXY, doubleCovered, umax, vmax) ;

    phaseRaw(ii, :) = atan2(dUV(:, 2), dUV(:, 1)) ;
    magUVRaw(ii, :) = sqrt(dUV(:, 1).^2 + dUV(:, 2).^2)  ;

    % Relative
    dY0 = Vall(ii, 1, 2) - Uall(ii, 1, 2) ;
    dX0 = Vall(ii, 1, 1) - Uall(ii, 1, 1) ;
    dXYm0 = [dX(:) - dX0(:) , dY(:) - dY0(:)] ;
    dUVm0 = QS.dXY2duv(im, dXYm0, doubleCovered, umax, vmax) ;
    phaseRelative(ii, :) = atan2(dUVm0(:, 2), dUVm0(:, 1)) ;
    magUVRelative(ii, :) = sqrt(dUVm0(:, 1).^2 + dUVm0(:, 2).^2) ;

end


figDir = QS.dir.tracking; 
for analysisMode = 1:2
    
    % Raw or relative
    if analysisMode == 1
        % raw displacement vectors
        duG = dusGeodesic ;
        duE = dusEuclidean ;
        phaseV = phaseRaw ;
        exten = '_raw' ;
    else
        % relative to starting points
        duG = dusGeodesic - dusGeodesic(:, 1) ;
        duE = dusEuclidean - dusEuclidean(:, 1) ;
        phaseV = phaseRelative ;
        exten = '' ;
    end


    %% 1D plot: curves for distance over time
    close all
    hf = figure('Position', [100 100 2*figX figY], 'units', 'centimeters');
    h1 = subplot(1, 3, 1) ;
    plot(timestamps, duG)
    ylabel(['relative motion of nuclei between layers [' QS.spaceUnits ']'], ...
        'interpreter', 'latex')
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
    if analysisMode == 1
        title('Geodesic separation', 'interpreter', 'latex')
    else
        title('Change in geodesic separation', 'interpreter', 'latex')
    end
    ylimsRelG = ylim ;

    h2 = subplot(1, 3, 2) ;
    plot(timestamps, duE)
    ylabel(['relative motion of nuclei between layers [' QS.spaceUnits ']'], ...
        'interpreter', 'latex')
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
    if analysisMode == 1
        title('Euclidean separation in 3D', 'interpreter', 'latex')
    else
        title('Change in Euclidean separation in 3D', 'interpreter', 'latex')
    end
    ylimsRelE = ylim ;

    h3 = subplot(1, 3, 3) ;
    plot(timestamps, euclideanDistanceTraveledU')
    ylabel(['distance traveled by endoderm nuclei [' QS.spaceUnits ']'], ...
        'interpreter', 'latex')
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
    title('Euclidean distance traveled by endoderm in 3D', ...
        'interpreter', 'latex')
    ylimsInp = ylim ;

    % set axis limits to match
    ylims = [min([ylimsRelG(1), ylimsRelE(1), ylimsInp(1)]), ...
        max([ylimsRelG(2), ylimsRelE(1), ylimsInp(2)])] ;
    ylim(ylims)
    axes(h1)
    ylim(ylims)
    axes(h2)
    ylim(ylims)
    saveas(gcf, fullfile(figDir, ['relative_motion_all_curves' exten '.pdf']))

    
    %% Highlight just two curves
    % 1D plot: curves for distance over time
    % otherBlue = [ 0.3010    0.7450    0.9330 ];
    % otherRed = [0.6350    0.0780    0.1840 ];
    % otherRed = [0.6350    0.3    0.32 ];
    colorCurv1 = [ 0.9290    0.6940    0.1250];
    colorCurv2 = [ 0.4660    0.6740    0.1880] + 0.1;
    
    ymin = min(min(euclideanDistanceTraveledU(:), min(duG(:)))) ;
    ymax = max(max(euclideanDistanceTraveledU(:), max(duG(:)))) ;
    
    close all
    hf = figure('Position', [100 100 figX figY], 'units', 'centimeters');
    h1 = subplot(1, 2, 1) ;
    indivID = [30, 33] ; % 15 fine, 25 too good, 
    others = setdiff(1:size(duG, 1), indivID) ;
    hold on;
    plot(timestamps, duG(others, :), 'color', [0.5, 0.5, 0.5]);
    plot(timestamps, duG(indivID(1), :), 'color', colorCurv1, ...
        'linewidth', 3); 
    plot(timestamps, duG(indivID(2), :), 'color', colorCurv2, ...
        'linewidth', 3); 
    ylabel(['relative motion of nuclei between layers [' QS.spaceUnits ']'], ...
        'interpreter', 'latex')
    xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
    % if analysisMode == 1
    %     title('Geodesic separation', 'interpreter', 'latex')
    % else
    %     title('Change in geodesic separation', 'interpreter', 'latex')
    % end
    ylim([ymin ymax])
    
    h2 = subplot(1, 2, 2) ;
    hold on;
    plot(timestamps, euclideanDistanceTraveledU(others, :)', ...
        'color', [0.5, 0.5, 0.5]);
    plot(timestamps, euclideanDistanceTraveledU(indivID(1), :), ...
        'color', colorCurv1, 'linewidth', 3); 
    plot(timestamps, euclideanDistanceTraveledU(indivID(2), :), ...
        'color', colorCurv2, 'linewidth', 3); 
    ylabel(['distance traveled by endoderm nuclei [' QS.spaceUnits ']'], ...
        'interpreter', 'latex')
    xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
    % title('Euclidean distance traveled by endoderm', 'interpreter', 'latex')
    ylim([ymin, ymax])
    h1pos = get(h1, 'position') ;
    h2pos = get(h2, 'position') ;
    % set(h2, 'position', [h2pos(1) h1pos(2) h1pos(3) h1pos(4)])

    saveas(gcf, fullfile(figDir, ['relative_motion_individual_curves' exten '.pdf']))
    save(fullfile(figDir, ['relative_motion_individual_curves' exten '.mat']), ...
        'timestamps', 'euclideanDistanceTraveledU', 'duG')






    % %% Histogram version
    % fewTidx = 1:15:nTimePoints ;
    % tColors = viridis(length(fewTidx)) ;
    % qq = 1 ;
    % clf
    % for tidx = fewTidx
    %     phase2plot = phaseV(:, tidx) ;
    %     phase2plot = phase2plot(~isnan(phase2plot)) ;
    %     polarhistogram(phase2plot, 15, ...
    %         'EdgeColor', 'none', 'FaceColor', tColors(qq, :), ...
    %         'FaceAlpha',.2)
    %     hold on;
    %     pause(0.2)
    %     qq = qq + 1 ;
    % end    
    % saveas(gcf, fullfile(figDir, ['relative_motion_histogram_phaseV_UV' exten '.pdf']))





    %% Polar/2d scatter plot over time -- Rough n Dirty method in UV
    close all
    hf = figure('Position', [100 100 figX figY], 'units', 'centimeters');
    if analysisMode == 1
        duV = magUVRaw ;
        phaseV = phaseRaw ;
    else
        duV = magUVRelative ;
        phaseV = phaseRelative ;
    end

    alphaVal = 0.5 ;
    for ii = 1:nTracks
        lastIdx = find(isnan(duG(ii, :))) ;
        if isempty(lastIdx)
            xx = abs(duV(ii, :)) .* cos(phaseV(ii, :)) ;
            yy = abs(duV(ii, :)) .* sin(phaseV(ii, :)) ;
            aColor = 1:nTimePoints ;
        else
            xx = abs(duV(ii, 1:lastIdx-1)) .* cos(phaseV(ii, 1:lastIdx-1)) ;
            yy = abs(duV(ii, 1:lastIdx-1)) .* sin(phaseV(ii, 1:lastIdx-1)) ;
            aColor = 1:lastIdx-1 ;
        end
        mfig = patch([xx NaN], [yy NaN], [aColor nTimePoints]);
        % alphamap = 0.2 * ones(length(xx), 1) ;
        % alphamap(end) = 0 ;
        set(mfig,'FaceColor','none','EdgeColor','flat',...
           'LineWidth',2, ...
           'FaceVertexAlphaData',alphaVal,...
           'EdgeAlpha',alphaVal)
    end
    xlabel('relative ap displacement [$\Delta x/L$]', ...
        'interpreter', 'latex')
    ylabel('relative dv displacement [$\Delta \phi/2\pi R]$', ...
        'interpreter', 'latex')
    title('Relative motion in pullback frame')
    cb = colorbar() ;
    ylims = get(cb, 'ylim') ;

    % colorbar label
    if contains(lower(QS.timeUnits), 'min') 
        ylabel(cb, 'time [hr]', 'interpreter', 'latex')
        set(cb, 'Ytick', [1, 16, 31, 46])
        set(cb, 'YtickLabel', ...
            {num2str(QS.timeInterval /60 * (timePoints(1) - min(timePoints))), ...
            num2str(QS.timeInterval /60 * (timePoints(16) - min(timePoints))), ...
            num2str(QS.timeInterval /60 * (timePoints(31) - min(timePoints))), ...
            num2str(QS.timeInterval /60 * (timePoints(46) - min(timePoints)))}) ;
    else
        ylabel(cg, ['time [' QS.timeUnits ']'], 'interpreter', 'latex')
    end
    saveas(gcf, fullfile(figDir, ['relative_motion_inPlaneUV' exten '.pdf']))

end

%% Parameterize by time -- mean+/- std
figDir = QS.dir.tracking ;
du_tp = dusGeodesic - dusGeodesic(:, 1) .* ones(size(dusGeodesic)) ;
means_du = mean(du_tp, 1, 'omitnan') ;
stds_du = std(du_tp, 1, 'omitnan') ;

means_dd = mean(euclideanDistanceTraveledU, 1, 'omitnan') ;
stds_dd = std(euclideanDistanceTraveledU, 1, 'omitnan') ;

% Handle minutes to hours conversion
xlabelString = ['time [' timeunits ']'] ;

% Plot it
close all
bigRed = [ 0.62    0.76    0.84 ];
bigBlue = [ 0.90    0.55    0.55] ;
big3 = [ 0.89    0.10    0.11] ;
big4 = [ 0.12    0.47    0.70] ;
lineProps1 = {'-','color', bigRed} ;
lineProps2 = {'-','color', bigBlue} ;

hf = figure('Position', [100 100 180 120], 'units', 'centimeters');
h1 = shadedErrorBar(timestamps, ...
    means_du, stds_du, 'lineProps', lineProps1) ;
hold on;
h2 = shadedErrorBar(timestamps, ...
    means_dd, stds_dd, 'lineProps', lineProps2) ;
legend('relative motion', 'motion')
ylabelString = ['displacement [' QS.spaceUnits ']'] ;
xlabel(xlabelString, 'interpreter', 'latex')
ylabel(ylabelString, 'interpreter', 'latex')
ylim([0, max(means_dd + stds_dd)])
saveas(gcf, fullfile(figDir, ['relative_motion.pdf']))
save( fullfile(figDir, 'relative_motion.mat'), 'timestamps', 'means_du', 'stds_du', 'means_dd', 'stds_dd')


%% Plot in 3d using initial du
% du0 = dusGeodesic(1, :) ;
% [du0, indx] = sort(du0) ;
% assert(length(du0) == nTimePoints)
% for qq = 1:nTimePoints
%     curv = dusGeodesic(indx(qq), :) ;
%     deltau = euclideanDistanceTraveledU(indx(qq), :) ;
%     plot3(du0(qq)* ones(size(curv)), timestamps, curv, '.-', 'color', blueColor)
%     plot3(du0(qq)* ones(size(curv)), timestamps, deltau, '.-', 'color', orangeColor)
%     hold on;
% end
% pause

disp('done!')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        