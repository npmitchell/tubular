function plotSegmentationStatisticsLobes(QS, options)
% estimateIntercalationRate(QS, options)
%
% Parameters
% ----------
% options: struct with fields
%   aspectLims : length 2 numeric array
%       aspect ratio limits for aspect*cos2theta
%
% NPMitchell 2021


% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;
maxCellSize = 200 ;  % maximum size allowed to include a cell
overwrite = false ;
corrected = false ;
overwriteImages = false ;
debug = false ;
[~, ~, ~, xyzlims] = QS.getXYZLims() ;
t0 = QS.t0set() ;
nDists = 5 ;
arlims = [-6, 6] ;
distlims = [0, 0.25] ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'correctedSegmentation')
    corrected = options.correctedSegmentation ;
elseif isfield(options, 'corrected')
    corrected = options.corrected ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'xyzlims')
    xyzlims = options.xyzlims ;
end
if isfield(options, 'aspectLims')
    arlims = options.aspectLims ;
end

% Directory preparation
if corrected
    segDir = fullfile(QS.dir.segmentation, 'seg3d_corrected') ;
    imDir = fullfile(segDir, 'images') ;
else
    segDir = fullfile(QS.dir.segmentation, 'seg3d') ;
    imDir = fullfile(QS.dir.segmentation, 'seg3d', 'images') ;
end
if ~exist(imDir, 'dir')
    mkdir(imDir)
end
if ~exist(fullfile(imDir, 'distFromFeatures'), 'dir')
    mkdir(fullfile(imDir, 'distFromFeatures'))
end

% Load indentation
pathlines = QS.getPullbackPathlines() ;
indent = QS.getPathlineIndentation() ;

%% Plot each lobe
features = QS.getFeatures() ;
folds = features.folds ;
try
    tidxs = QS.xp.tIdx(timePoints) ;
catch
    tidxs = zeros(size(timePoints), 'uint8') ;
    for qq = 1:length(timePoints)
        tidxs(qq) = QS.xp.tIdx(timePoints(qq)) ;
    end
end

% define colors
colors = define_colors() ;

% preallocate
cos2thetaM = cell(length(folds), 1) ;
sin2thetaM = cos2thetaM ;
aspectM = cos2thetaM ;
nbins = 10 ;
for qq = 1:length(folds) + 1
    cos2thetaM{qq} = zeros(nbins, length(tidxs)) ;
    sin2thetaM{qq} = zeros(nbins, length(tidxs)) ;
    aspectM{qq} = zeros(nbins, length(tidxs)) ;
    distFits{qq} = zeros(length(tidxs), 4) ;
    distCorrs{qq} = zeros(length(tidxs), 1) ;
    cellIndentation{qq} = cell(length(tidxs), 1) ;
    distFold{qq} = cell(length(tidxs), 1) ;
    anisotropy{qq} = cell(length(tidxs), 1) ;
end
dmy = 1 ;
for tidx = tidxs 
    tp = timePoints(dmy) ;
    QS.setTime(tp) ;
    nU = QS.nU ;
    fold0 = double(folds(tidx, :)) / double(nU) ;
    nLobes = length(fold0(:)) + 1 ;
    foldt = [0; fold0(:); 1] ;
    
    % Load segs
    if corrected
        seg2d = QS.getCurrentSegmentation2DCorrected() ;
        seg3d = QS.getCurrentSegmentation3DCorrected() ;
    else
        seg2d = QS.getCurrentSegmentation2D() ;
        seg3d = QS.getCurrentSegmentation3D() ;
    end
    coordSys = seg2d.coordSys ;
    segIm = seg2d.segIm ;
    seg2d = seg2d.seg2d ;
    Lx = size(segIm, 2) ;
    % ensure 3d and 2d segmentation are based on the same coordsys
    assert(strcmpi(seg3d.coordSys, coordSys))
    seg3d = seg3d.seg3d ;    
    
    % divide into chambers
    close all
    for lobe = 1:nLobes
        % find cells in this current chamber
        if strcmpi(coordSys, 'spsme') || strcmpi(coordSys, 'sprsme')
            inLobe = find(seg2d.cdat.centroid(:, 1) > Lx * foldt(lobe) & ...
                seg2d.cdat.centroid(:, 1) < Lx * foldt(lobe+1)) ;
        else
            error(['code for fold locations in this coordSys here: ' coordSys])
        end
        keep = intersect(inLobe, seg3d.statistics.keep) ;
        
        ars = sqrt(seg3d.qualities.moment2 ./ seg3d.qualities.moment1) ;
        strength = ars - 1 ;
        edges = linspace(-1, 1, nbins+1) ;
        edgesAR = linspace(1, 5, nbins+1) ;
        tmp = histcounts(cos(2* seg3d.qualities.ang1(keep)), edges) ;
        cos2thetaM{lobe}(:, dmy) = tmp / sum(tmp) ;
        tmp = histcounts(sin(2* seg3d.qualities.ang1(keep)), edges) ;
        sin2thetaM{lobe}(:, dmy) = tmp / sum(tmp) ;
        tmp = histcounts(ars(keep), edgesAR) ;
        aspectM{lobe}(:, dmy) = tmp / sum(tmp) ;

        % anisotropy versus distance from fold
        distFold{lobe}{dmy} = min(abs(seg2d.cdat.centroid(keep, 1) / Lx - fold0), [], 2) ;
        % check it
        % scatter(seg2d.cdat.centroid(:, 1), seg2d.cdat.centroid(:, 2), 5, distFold{lobe}{dmy})
        anisotropy{lobe}{dmy} = (ars(keep) -1).* cos(2* seg3d.qualities.ang1(keep)) ;
        
        % STORE ALL FOR TRAJECTORIES
        strengthsAll{lobe}{dmy} = strength(keep) ;
        thetasAll{lobe}{dmy} = seg3d.qualities.ang1(keep) ;
        weightsAll{lobe}{dmy} = seg3d.statistics.weights(keep) ;       
        
        % Interpolate centroids onto pathline triangulation
        XX = pathlines.vertices.vX(tidx, :, :) ;
        YY = pathlines.vertices.vY(tidx, :, :) ;
        indtV = indent(tidx,:, :) ;
        indtI = scatteredInterpolant(XX(:), YY(:), indtV(:), ...
            'natural', 'nearest') ;
        xV = seg2d.cdat.centroid(keep, 1) ;
        yV = seg2d.cdat.centroid(keep, 2) ;
        cellIndentation{lobe}{dmy} = indtI(xV, yV) ;
    end
    
    % Plot anisotropy vs distfold
    imfn = fullfile(imDir, 'distFromFeatures', sprintf('anisotropy_%06d.png', tp)) ;
    plot_fit = ~exist(imfn, 'file') || overwrite || overwriteImages ;
    for lobe = 1:nLobes
        
        % Fit to line
        [cc, SS] = polyfit(distFold{lobe}{dmy}, anisotropy{lobe}{dmy},1);
       
        cf = fit(distFold{lobe}{dmy}, anisotropy{lobe}{dmy},'poly1'); 
        cf_coeff = coeffvalues(cf);
        cf_confint = confint(cf);
        aa = cf_coeff(1);
        bb = cf_coeff(2);
        aa_unc = (cf_confint(2,1) - cf_confint(1,1))/2;
        bb_unc = (cf_confint(2,2) - cf_confint(1,2))/2;

        % Display evaluated equation y = m*x + b
        disp(['Equation is y = ' num2str(aa) '*x + ' num2str(bb)])
        % Evaluate fit equation using polyval
        
        xx = linspace(distlims(1), distlims(2), nDists) ;
        [yy, deltay] = polyval(cc, xx, SS) ;
        yy(all(xx > distFold{lobe}{dmy})) = NaN ;
        deltay(all(xx > distFold{lobe}{dmy})) = NaN ;
        
        if plot_fit
            plot(distFold{lobe}{dmy}, anisotropy{lobe}{dmy}, ...
                '.', 'color', colors(lobe, :)) ;
            hold on;
            ylim(arlims)
            xlim(distlims) 
            errorbar(xx, yy, deltay, 'color', colors(lobe, :)) ;
        end
        tmp = corrcoef(distFold{lobe}{dmy}, anisotropy{lobe}{dmy}) ;
        distFits{lobe}(dmy, :) = [aa, bb, aa_unc, bb_unc] ;
        distCorrs{lobe}(dmy) = tmp(1, 2) ;
        fitEvalDist{lobe}(dmy, :) = xx ;
        fitEvalElong{lobe}(dmy, :) = yy ;
        fitEvalUnc{lobe}(dmy, :) = deltay ;
    end
    
    if plot_fit
        xlabel(['distance from fold, [$\zeta/L$]'], ...
            'interpreter', 'latex')
        ylabel('anisotropy, $(a/b-1) \cos 2\theta $', ...
            'interpreter', 'latex')
        title(['$t = $' num2str(tp - t0) ' ' QS.timeUnits], 'interpreter', 'latex')
        disp(['Saving ' imfn])
        saveas(gcf, imfn)   
    end
    
    dmy = dmy + 1 ;
end

%% Colormaps
colormap(cividis)
try
    getPyPlot_cMap('cividis', length(tidxs)) ;
    cmap0 = cividis(length(tidxs)) ;
catch
    tmp = cividis ;
    cmap0 = tmp(round(linspace(1, length(tmp), length(tidxs))), :) ;
end

%% anisotropy change as kymograph
clf
mats = {} ;
xx = linspace(distlims(1), distlims(2), nDists) ;
for lobe = 1:nLobes
    dmy = 1 ;
    % Store each 2d array in a matrix
    mat = zeros(length(tidxs), nDists) ;
    for tidx = tidxs 
        % errorbar(fitEvalDist{lobe}(dmy, :), ...
        %     fitEvalElong{lobe}(dmy, :), ...
        %     fitEvalUnc{lobe}(dmy, :), 'Color', cmap0(dmy, :)) ;
        mat(dmy, :) = fitEvalElong{lobe}(dmy, :) ;
        dmy = dmy + 1 ;
    end
    subplot(2, 2, lobe)
    imagesc(xx, timePoints-t0, mat); 
    if mod(lobe, 2) == 1
        ylabel(['time, [' QS.timeUnits ']'], 'interpreter', 'latex')
    end
    if lobe > 2
        xlabel('distance, $\zeta/L$', 'interpreter', 'latex')
    end
    cb = colorbar() ;
    ylabel(cb, '$\sqrt{I_1/I_2}$', 'interpreter', 'latex')
    mats{lobe} = mat ;
end
saveas(gcf, fullfile(segDir, 'dynamics_shapeChanges.png'))

% Take derivative
for lobe = 1:nLobes
    subplot(2, 2, lobe)
    dmat = diff(mats{lobe}, 1) ;
    dmat(isnan(dmat(:))) = 0 ;
    imagesc(xx, timePoints-t0, dmat); 
    caxis([-2, 2])
    colormap(blueblackred)
    cb = colorbar() ;
end

%% Plot anisotropy versus indentation
allI = [] ;
allE = [] ;
clf
dmy = 1 ;
for tidx = tidxs 
    for lobe = 1:nLobes
        scatter(cellIndentation{lobe}{dmy}, anisotropy{lobe}{dmy}, ...
            3, 'MarkerEdgeColor', colors(lobe, :)) ;
        hold on;
        allI = [allI; cellIndentation{lobe}{dmy}(:)] ;
        allE = [allE; anisotropy{lobe}{dmy}(:)] ;
    end
    dmy = dmy + 1 ;
end
xlim([-0.25, 1])
ylim(arlims)
legend({'lobe 1', 'lobe 2', 'lobe 3', 'lobe 4'})
xlabel('indentation, $\delta r / r_0$', 'interpreter', 'latex')
ylabel('anisotropy, $(a/b-1) \cos 2\theta $', ...
    'interpreter', 'latex')
saveas(gcf, fullfile(segDir, 'anisotropy_indentation_lobes.png'))

%% Plot anisotropy versus indentation
clf
dmy = 1 ;
% cmap0 = cubehelix(length(tidxs), 0.69, -0.6, 1.8, 1.13, [0.5, 0.5], [0.23, 0.96]) ;
for tidx = tidxs
    for lobe = 1:nLobes
        scatter(cellIndentation{lobe}{dmy}, anisotropy{lobe}{dmy}, ...
            3, 'MarkerEdgeColor', cmap0(dmy, :)) ;
        hold on;
    end
    dmy = dmy + 1 ;
end
xlim([-0.25, 1])
ylim(arlims)
xlabel('indentation, $\delta r / r_0$', 'interpreter', 'latex')
ylabel('anisotropy, $(a/b-1) \cos 2\theta $', ...
    'interpreter', 'latex')
cb = colorbar() ;
%Create 8 ticks from zero to 1
nTicks = 8 ;
cb.Ticks = linspace(0, 1, nTicks) ; 
cb.TickLabels = num2cell(round( ...
    linspace(min(timePoints-t0), max(timePoints-t0), nTicks))) ; 
ylabel(cb, ['time [' QS.timeUnits ']'], 'interpreter', 'latex')
saveas(gcf, fullfile(segDir, 'anisotropy_indentation_time.png'))


%% Plot histograms
clim = 5 / nbins ;
for lobe = 1:nLobes
    imfn = fullfile(segDir, sprintf('cell_anisotropy_hist_lobe%d.png', lobe)) ;
    clf
    colormap(cividis)
    subplot(2, 2, 1)
    imagesc(timePoints - t0, edges, cos2thetaM{lobe})
    set(gca,'YDir','normal')
    caxis([0, clim])
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
    ylabel('nematic orientation $\cos2\theta$',   'interpreter', 'latex')
    subplot(2, 2, 2)
    imagesc(timePoints - t0, edges, sin2thetaM{lobe})
    set(gca,'YDir','normal')
    caxis([0, clim])
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
    ylabel('nematic orientation $\sin2\theta$',   'interpreter', 'latex')
    subplot(2, 2, 3)
    imagesc(timePoints - t0, edgesAR, aspectM{lobe})
    set(gca,'YDir','normal')
    caxis([0, clim])
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
    ylabel('aspect ratio $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
    subplot(4, 2, 6)
    cb = colorbar('location', 'south') ;
    ylabel(cb, 'probability', 'interpreter', 'latex')
    caxis([0, clim])
    axis off
    sgtitle(['nematic dynamics for chamber ' num2str(lobe)], ...
        'interpreter', 'latex')
    disp(['Saving fig: ' imfn])
    saveas(gcf, imfn) 
    
end

%% Store lobes 1 and 2 together for comparison between optogenetic and WT
imfn = fullfile(segDir, 'cell_anisotropy_traj_lobes12.png') ;    
clf
confidence = 0.5 ;
clipping_radius = inf ;
nTimes = length(strengthsAll{1}) ;
meanxs = zeros(nTimes, 1) ;
meanys = zeros(nTimes, 1) ;
stdmeanxs = zeros(nTimes, 1) ;
stdmeanys = zeros(nTimes, 1) ;
stdxs = zeros(nTimes, 1) ;
stdys = zeros(nTimes, 1) ;
covs = zeros(nTimes, 2, 2) ;

maxx = zeros(nTimes, 1) ;
maxy = zeros(nTimes, 1) ;
for dmy = 1:nTimes
    
    x1 = strengthsAll{1}{dmy} .* cos(2*thetasAll{1}{dmy}) ;
    y1 = strengthsAll{1}{dmy} .* sin(2*thetasAll{1}{dmy}) ;
    x2 = strengthsAll{2}{dmy} .* cos(2*thetasAll{2}{dmy}) ;
    y2 = strengthsAll{2}{dmy} .* sin(2*thetasAll{2}{dmy}) ;
    
    xx = [x1(:); x2(:)] ;
    yy = [y1(:); y2(:)] ;
    
    weights = [weightsAll{1}{dmy}; weightsAll{2}{dmy}] ;
    
    % meanx = sum(weightsAll{lobe}{dmy} .* xx) / sum(weightsAll{lobe}{dmy});
    % meany = sum(weightsAll{lobe}{dmy} .* yy) / sum(weightsAll{lobe}{dmy}) ;
    [meanxs(dmy), stdmeanxs(dmy), stdxs(dmy)] = ...
        weightedStats(xx, weights, 'w') ;
    [meanys(dmy), stdmeanys(dmy), stdys(dmy)] = ...
        weightedStats(yy, weights, 'w') ;

    h = scatter(xx, yy, 2, colors(1, :), ...
        'filled', 'markeredgecolor', 'none') ;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h.MarkerFaceAlpha = 0.2 ;

    cc = cov(xx, yy) ;
    covs(dmy, :, :) = cc ;
    % [h, ellip] = error_ellipse(cc, [meanxs(dmy), meanys(dmy)]) ;
    [xx,yy] = covarianceEllipse2d(cc, [meanxs(dmy), meanys(dmy)], ...
        confidence, clipping_radius) ;
    maxx(dmy) = max(xx(:)) ;
    minx(dmy) = min(xx(:)) ;
    maxy(dmy) = max(yy(:)) ;
    miny(dmy) = min(yy(:)) ;
    ellipses{dmy} = [xx, yy] ;
    % set(h, 'color', colorList(dmy, :))
    % set(h, 'linewidth', 2)
    % the following line skip the name of the previous plot from the legend
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % set(h, 'alpha', 0.2)
    hold on;

end

covariancesL12 = covs ;
maxysL12 = maxy ;
minysL12 = miny ;
maxxsL12 = maxx ;
minxsL12 = minx ;
meanxsL12 = meanxs;
meanysL12 = meanys ;
% hs = (meanxsL12, meanysL12, '.-', 'linewidth', 1, 'color', colors(1, :)) ;
lineProps = {'.-', 'color', colors(1, :)} ;
hs=shadedErrorBar(meanxsL12, meanysL12, stdys, 'lineProps', lineProps) ;

axis equal
% plot(meanxsL12, maxysL12, '-', 'color', colors(1, :)) ;
% plot(meanxsL12, minysL12, '-', 'color', colors(1, :)) ;
% plot(meanxsL12, meanysL12+stdys, '-', 'color', colors(1, :)) ;
% plot(meanxsL12, meanysL12-stdys, '-', 'color', colors(1, :)) ;
plot(meanxsL12, meanysL12+stdmeanys, '--', 'color', colors(1, :)) ;
plot(meanxsL12, meanysL12-stdmeanys, '--', 'color', colors(1, :)) ;

ylim([-3, 3])
legend(hs.patch, {'lobes 1 & 2'})
axis equal
ylim([-3, 3])
daspect([1 1 1]) 
xlabel('$2|Q| \cos 2 \theta$, lobes 1\&2', 'interpreter', 'latex')
ylabel('$2|Q| \sin 2 \theta$, lobes 1\&2', 'interpreter', 'latex')
saveas(gcf, imfn)

%% Save statistics summary Lobes 1&2
summaryFn = fullfile(segDir, 'stats_summary_L12.mat') ;
timeStamps = (timePoints - t0) * QS.timeInterval ;
timeUnits = QS.timeUnits ;
save(summaryFn, 'timePoints', 'timeStamps', 'timeUnits', ...
    'meanxsL12', 'meanysL12', 'covariancesL12', ...
    'stdmeanxs', 'stdmeanys', 'stdxs', 'stdys') ;

%% Plot trajectory in phase space
% Plot trajectory in phase space
imfn = fullfile(segDir, 'cell_anisotropy_traj_lobes.png') ;    
clf
confidence = 0.5 ;
clipping_radius = inf ;
maxysAll = {} ;
minysAll = {} ;
maxxsAll = {} ;
minxsAll = {} ;
stdmeanxsAll = {} ;
stdmeanysAll = {} ;
stdxsAll = {} ;
stdysAll = {} ;
ellipsesAll = {} ;
covariancesAll = {} ;
meanxsAll = {} ;
meanysAll = {} ;
for lobe = 1:nLobes
    nTimes = length(strengthsAll{lobe}) ;
    % colorList = cmap0(round(linspace(1, size(cmap0, 1), nTimes)), :) ;
    % colorList = cubehelix(nTimes, 2.54, 1, 1.3, 1.3, [.6, 0], [0.5, 1]) ;
    colorList = squeeze(outerProduct(colors(lobe, :), linspace(0.25, 1, nTimes)))' ;
    meanxs = zeros(nTimes, 1) ;
    meanys = zeros(nTimes, 1) ;
    stdmeanxs = zeros(nTimes, 1) ;
    stdmeanys = zeros(nTimes, 1) ;
    stdxs = zeros(nTimes, 1) ;
    stdys = zeros(nTimes, 1) ;
    covs = zeros(nTimes, 2, 2) ;

    maxx = zeros(nTimes, 1) ;
    maxy = zeros(nTimes, 1) ;
    for dmy = 1:nTimes
        xx = strengthsAll{lobe}{dmy} .* cos(2*thetasAll{lobe}{dmy}) ;
        yy = strengthsAll{lobe}{dmy} .* sin(2*thetasAll{lobe}{dmy}) ;
        % meanx = sum(weightsAll{lobe}{dmy} .* xx) / sum(weightsAll{lobe}{dmy});
        % meany = sum(weightsAll{lobe}{dmy} .* yy) / sum(weightsAll{lobe}{dmy}) ;
        [meanxs(dmy), stdmeanxs(dmy), stdxs(dmy)] = ...
            weightedStats(xx, weightsAll{lobe}{dmy}, 'w') ;
        [meanys(dmy), stdmeanys(dmy), stdys(dmy)] = ...
            weightedStats(yy, weightsAll{lobe}{dmy}, 'w') ;

        h = scatter(xx, yy, 2, colorList(dmy,:), ...
            'filled', 'markeredgecolor', 'none') ;
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.MarkerFaceAlpha = 0.2 ;
        
        cc = cov(xx, yy) ;
        covs(dmy, :, :) = cc ;
        % [h, ellip] = error_ellipse(cc, [meanxs(dmy), meanys(dmy)]) ;
        [xx,yy] = covarianceEllipse2d(cc, [meanxs(dmy), meanys(dmy)], ...
            confidence, clipping_radius) ;
        maxx(dmy) = max(xx(:)) ;
        minx(dmy) = min(xx(:)) ;
        maxy(dmy) = max(yy(:)) ;
        miny(dmy) = min(yy(:)) ;
        ellipses{dmy} = [xx, yy] ;
        % set(h, 'color', colorList(dmy, :))
        % set(h, 'linewidth', 2)
        % the following line skip the name of the previous plot from the legend
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        % set(h, 'alpha', 0.2)
        hold on;
        
    end
    
    covariancesAll{lobe} = covs ;
    ellipsesAll{lobe} = ellipses ;
    stdmeanxsAll{lobe} = stdmeanxs ;
    stdmeanysAll{lobe} = stdmeanys ;
    stdxsAll{lobe} = stdxs ;
    stdysAll{lobe} = stdys ;
    maxysAll{lobe} = maxy ;
    minysAll{lobe} = miny ;
    maxxsAll{lobe} = maxx ;
    minxsAll{lobe} = minx ;
    meanxsAll{lobe} = meanxs;
    meanysAll{lobe} = meanys ;
end
hs = [] ;
for lobe = 1:nLobes
    meanxs = meanxsAll{lobe} ;
    meanys = meanysAll{lobe} ;
    maxys = maxysAll{lobe} ;
    minys = minysAll{lobe} ;
    hs(lobe) = plot(meanxs, meanys, '.-', ...
        'linewidth', 1, 'color', colors(lobe, :)) ;
    axis equal
    plot(meanxs, maxys, '--', 'color', colors(lobe, :)) ;
    plot(meanxs, minys, '--', 'color', colors(lobe, :)) ;
end
ylim([-3, 3])
legend(hs, {'lobe 1', 'lobe 2', 'lobe 3', 'lobe 4'})
axis equal
ylim([-3, 3])
daspect([1 1 1]) 
xlabel('$|Q| \cos 2 \theta$', 'interpreter', 'latex')
ylabel('$|Q| \sin 2 \theta$', 'interpreter', 'latex')
saveas(gcf, imfn)

%% Save statistics summary
summaryFn = fullfile(segDir, 'stats_summary.mat') ;
timeStamps = (timePoints - t0) * QS.timeInterval ;
timeUnits = QS.timeUnits ;
save(summaryFn, 'timePoints', 'timeStamps', 'timeUnits', ...
    'meanxsAll', 'meanysAll', 'covariancesAll', ...
    'stdmeanxs', 'stdmeanys', 'stdxs', 'stdys')
