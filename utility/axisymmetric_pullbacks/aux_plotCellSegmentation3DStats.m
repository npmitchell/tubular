% aux_plotCellSegmentation3DStats

%% Define some colors
colors = define_colors() ;
bluecol = colors(1, :) ;
redcol = colors(2, :) ;
yelcol = colors(3, :) ;
greencol = colors(5, :) ;
CScolors = colors([9,5], :) ;
figW = 9 ;
figH = 4.5 ;

%% Plot meanQ_A AR +/- unc, meanQ Theta +/- unc over time
imfn = fullfile(tubi.dir.segmentation, segSubDir, 'cell_meanQAaspectRatio') ;
timestamps = timePoints - t0 ;
if contains(lower(tubi.timeUnits), 'min')
    timestamps = timestamps / 60 ;
    timeunits = 'hr' ;
else
    timeunits = tubi.timeUnits ;
end
xlims = [min(timestamps) - 0.25*min(diff(timestamps)), ...
        max(timestamps) + 0.25*min(diff(timestamps))] ;

close all 
fig = figure('units', 'centimeters', 'position', [0, 0, figW, figH]) ;

% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timestamps, fliplr(timestamps)] ;
subplot(1, 2, 1)
fill(x2, [mQAar - mQAarStd, ...
    fliplr(mQAar + mQAarStd)], ...
    bluecol, 'facealpha', 0.3, 'edgecolor', 'none');

hold on;
errorbar(timestamps, mQAar, mQAarSte, ...
    '-', 'color', bluecol)
xlim(xlims) 
xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
ylabel('cell aspect ratio, $a/b$',   'interpreter', 'latex')

% Theta
subplot(1, 2, 2)
fill(x2, [0.5*mod(mQAtheta*2 - mQAthetaStd*2, 2*pi), ...
    fliplr(0.5*mod(mQAtheta*2 + mQAthetaStd*2, 2*pi))],...
    greencol, 'facealpha', 0.3, 'edgecolor', 'none');

hold on;
errorbar(timestamps, 0.5*mod(2*mQAtheta, 2*pi), ...
    mod(mQAthetaSte, 2*pi), ... 
    '-', 'color', greencol)
ylim([0,pi])
yticks([0, pi*0.5, pi])
yticklabels({'0', '\pi/2', '\pi'})
ylabel('cell orientation, $\theta$',   'interpreter', 'latex')

xlim(xlims)
xlabel(['time [' timeunits ']'], 'interpreter', 'latex')

saveas(gcf, [imfn '.png'])
saveas(gcf, [imfn '.pdf'])


%% Plot meanQ_E AR +/- unc, meanQ Theta +/- unc over time
imfn = fullfile(tubi.dir.segmentation, segSubDir, 'cell_meanQEaspectRatio') ;

close all 
fig = figure('units', 'centimeters', 'position', [0, 0, figW, figH]) ;

% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timestamps, fliplr(timestamps)] ;
subplot(1, 2, 1)
fill(x2, [mQEar - mQEarStd, ...
    fliplr(mQEar + mQEarStd)], ...
    bluecol, 'facealpha', 0.3, 'edgecolor', 'none');

hold on;
errorbar(timestamps, mQEar, mQEarSte, ...
    '-', 'color', bluecol)
xlim(xlims) 
xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
ylabel('cell aspect ratio, $a/b$',   'interpreter', 'latex')

% Theta
subplot(1, 2, 2)
fill(x2, [0.5*mod(mQEtheta*2 - mQEthetaStd*2, 2*pi), ...
    fliplr(0.5*mod(mQEtheta*2 + mQEthetaStd*2, 2*pi))],...
    CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');

hold on;
errorbar(timestamps, 0.5*mod(2*mQEtheta, 2*pi), ...
    mod(mQEthetaSte, 2*pi), ... 
    '-', 'color', CScolors(1, :))
ylim([0,pi])
yticks([0, pi*0.5, pi])
yticklabels({'0', '\pi/2', '\pi'})
ylabel('cell orientation, $\theta$',   'interpreter', 'latex')

xlim(xlims)
xlabel(['time [' timeunits ']'], 'interpreter', 'latex')

saveas(gcf, [imfn '.png'])
saveas(gcf, [imfn '.pdf'])


%% Plot meanQ_E AR +/- unc, meanQ Theta +/- unc over time
imfn = fullfile(tubi.dir.segmentation, segSubDir, 'cell_meanQEeccentricity') ;

close all 
fig = figure('units', 'centimeters', 'position', [0, 0, figW, figH]) ;

% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timestamps, fliplr(timestamps)] ;
subplot(1, 2, 1)

% Convert to eccentricity
[mQEe, mQEeStd] = aspectRatioToEccentricity(mQEar, mQEarStd) ;
[~, mQEeSte] = aspectRatioToEccentricity(mQEar, mQEarSte) ;
fill(x2, [mQEe - mQEeStd, ...
    fliplr(mQEe + mQEeStd)], ...
    bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
errorbar(timestamps, mQEe, mQEeSte, ...
    '-', 'color', bluecol)
xlim(xlims) 
xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
ylabel('cell eccentricity, $e$',   'interpreter', 'latex')

% Theta
subplot(1, 2, 2)
fill(x2, [0.5*mod(mQEtheta*2 - mQEthetaStd*2, 2*pi), ...
    fliplr(0.5*mod(mQEtheta*2 + mQEthetaStd*2, 2*pi))],...
    greencol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
errorbar(timestamps, 0.5*mod(2*mQEtheta, 2*pi), ...
    mod(mQEthetaSte, 2*pi), ... 
    '-', 'color', greencol)
ylim([0,pi])
yticks([0, pi*0.5, pi])
yticklabels({'0', '\pi/2', '\pi'})
ylabel('cell orientation, $\theta$',   'interpreter', 'latex')

xlim(xlims)
xlabel(['time [' timeunits ']'], 'interpreter', 'latex')

saveas(gcf, [imfn '.png'])
saveas(gcf, [imfn '.pdf'])


%% Plot mean +/- pctile over time --- RAW GEOMETRY, NOT ALIGNED OR AVERAGED
imfn = fullfile(tubi.dir.segmentation, segSubDir, 'cell_anisotropyRaw') ;
timestamps = timePoints - t0 ;
if contains(lower(tubi.timeUnits), 'min')
    timestamps = timestamps / 60 ;
    timeunits = 'hr' ;
else
    timeunits = tubi.timeUnits ;
end
xlims = [min(timestamps) - 0.25*min(diff(timestamps)), ...
        max(timestamps) + 0.25*min(diff(timestamps))] ;

close all 
fig = figure('units', 'centimeters', 'position', [0, 0, figW, figH]) ;
stdStyles = {'_quartiles', '_std'} ;
thetaStyles = {'_c2ts2t', '_theta'} ;
for tStyle = 1:2
    for stdStyle = 1:2
        if ~(stdStyle == 1 && tStyle==2)
            disp('plotting raw aratios and thetas')
            % shade(timePoints - t0, bndlow, timePoints, bndhigh)
            x2 = [timestamps, fliplr(timestamps)] ;
            subplot(1, 2, 1)
            if stdStyle == 1
                fill(x2, [aratio_low25, fliplr(aratio_high75)], bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
            else
                fill(x2, [aratio_mean - aratio_std, ...
                    fliplr(aratio_mean + aratio_std)], ...
                    bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
            end
            hold on;
            errorbar(timestamps, aratio_mean, aratio_ste, ...
                '-', 'color', bluecol)
            xlim(xlims) 
            xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
            ylabel('cell aspect ratio, $a/b$',   'interpreter', 'latex')
            title('endoderm anisotropy over time', 'interpreter', 'latex')

            if tStyle == 1
                % Cos/Sin(2theta)
                set(fig,'defaultAxesColorOrder', CScolors);
                subplot(1, 2, 2)
                yyaxis left
                if stdStyle == 1
                    fill(x2, [c2t_low25, fliplr(c2t_high75)], CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
                else
                    fill(x2, [min(1, max(-1, c2t_mean - c2t_std)), ...
                        fliplr(min(1, max(-1,c2t_mean + c2t_std)))],...
                        CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
                end
                hold on;
                % shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
                errorbar(timestamps, c2t_mean, ...
                    min(c2t_ste, c2t_mean+1), ...  % bound from below at -1      
                    min(c2t_ste, 1-c2t_mean), ... % bound from above at 1
                    '-', 'color', CScolors(1, :))
                ylim([-1,1])
                yyaxis right ;
                if stdStyle == 1
                    fill(x2, [s2t_low25, fliplr(s2t_high75)], CScolors(2, :), 'facealpha', 0.3, 'edgecolor', 'none');
                else
                    fill(x2, [min(1, max(-1, s2t_mean - s2t_std)), ...
                        fliplr(min(1, max(-1,s2t_mean + s2t_std)))], ...
                        CScolors(2, :), 'facealpha', 0.3, 'edgecolor', 'none');
                end
                errorbar(timestamps, s2t_mean, ...
                    min(s2t_ste, s2t_mean+1), ... % bound from below at -1
                    min(s2t_ste, 1-s2t_mean), ... % bound from above at 1
                    '-', 'color', CScolors(2, :))
                ylim([-1,1])

                yyaxis left
                ylabel('cell orientation, $\cos 2\theta$',   'interpreter', 'latex')
                yyaxis right
                ylabel('cell orientation, $\sin 2\theta$',   'interpreter', 'latex')
            else

                % Theta
                subplot(1, 2, 2)
                if stdStyle == 2
                    fill(x2, [0.5*mod(theta_mean*2, 2*pi)-0.5*mod(theta_std*2, 2*pi), ...
                        fliplr(0.5*mod(theta_mean*2, 2*pi) + 0.5*mod(theta_std*2, 2*pi))],...
                        greencol, 'facealpha', 0.3, 'edgecolor', 'none');
                else
                    error('should not enter this case -- no quartile defined')
                end
                hold on;
                errorbar(timestamps, 0.5*mod(2*theta_mean, 2*pi), ...
                    theta_ste, ... 
                    '-', 'color', greencol)
                ylim([0,pi])
                yticks([0, pi*0.5, pi])
                yticklabels({'0', '\pi/2', '\pi'})
                ylabel('cell orientation, $\theta$',   'interpreter', 'latex')
            end

            xlim(xlims)
            xlabel(['time [' timeunits ']'], 'interpreter', 'latex')

            saveas(gcf, [imfn, stdStyles{stdStyle}, thetaStyles{tStyle}, '.png'])
            saveas(gcf, [imfn, stdStyles{stdStyle}, thetaStyles{tStyle}, '.pdf'])
            
        end
        clf
    end
end

%% Plot nematic strength and direction for each lobe 
imfn = fullfile(tubi.dir.segmentation, segSubDir, ...
    'cell_anisotropy_lobes') ;
AEstr = {'_aspectWeighted', '_eccentricityWeighted'} ;

for AorE = 1:2
    close all
    fig = figure('units','centimeters','position',[0,0,figW, figW], 'visible', 'off') ;
    set(fig,'defaultAxesColorOrder', [bluecol; greencol]);
    for lobe = 1:nLobes
        subplot(ceil(nLobes * 0.5), 2, lobe)
        if AorE ==1
            midline = squeeze(meanQALobeAspects(lobe, :)) ;
            uncs = squeeze(meanQALobeAspectStds(lobe, :)) ;
            errs = squeeze(meanQALobeAspectStes(lobe, :)) ;
        else
            midline = squeeze(meanQELobeAspects(lobe, :)) ;
            uncs = squeeze(meanQELobeAspectStds(lobe, :)) ;
            errs = squeeze(meanQELobeAspectStes(lobe, :)) ;
        end
        timestamps = timePoints - t0 ;
        if contains(tubi.timeUnits, 'min')
            timestamps = timestamps / 60 ;
        end
        x2 = [timestamps, fliplr(timestamps)] ;
        fill(x2,[midline-uncs, fliplr(midline+uncs)], ...
            bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
        hold on;
        errorbar(timestamps, midline, errs, 'color', bluecol)
        ylim([1, Inf])

        % THETA (as angle)
        yyaxis right
        
        if AorE ==1
            midline = 0.5*(mod(2*squeeze(meanQALobeThetas(lobe, :)), 2*pi)) ;
            uncs = 0.5*(mod(2*squeeze(meanQALobeThetaStds(lobe, :)), 2*pi)) ;
            errs = 0.5*(mod(2*squeeze(meanQALobeThetaStes(lobe, :)), 2*pi)) ;
        else
            midline = 0.5*(mod(2*squeeze(meanQELobeThetas(lobe, :)), 2*pi)) ;
            uncs = 0.5*(mod(2*squeeze(meanQELobeThetaStds(lobe, :)), 2*pi)) ;
            errs = 0.5*(mod(2*squeeze(meanQELobeThetaStes(lobe, :)), 2*pi)) ;
        end
        fill(x2,[midline-uncs, fliplr(midline+uncs)], ...
            greencol, 'facealpha', 0.3, 'edgecolor', 'none');
        hold on;
        errorbar(timestamps, midline, errs, 'color', greencol)
        % fill(x2, [c2t_low, fliplr(c2t_high)], redcol, 'facealpha', 0.3, 'edgecolor', 'none');
        ylim([0, pi])
        yticks([0, pi*0.5, pi])
        yticklabels({'0', '\pi/2', '\pi'})

        if mod(lobe, 2) == 1 
            yyaxis left
            ylabel(['lobe ' num2str(lobe) ' aspect ratio'],   'interpreter', 'latex')
        else
            yyaxis right
            ylabel('nematic orientation $\theta$',   'interpreter', 'latex')
        end

        % Time label
        if lobe > nLobes - 2
            if contains(tubi.timeUnits, 'min')
                xlabel('time [hr]', 'interpreter', 'latex')  
            else
                xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')  
            end
        end
    end
    saveas(gcf, [imfn AEstr{AorE} '.png'])
    saveas(gcf, [imfn AEstr{AorE} '.pdf'])
end




%% Plot nematic strength and direction for each lobe as shaded errorplot
imfn = fullfile(tubi.dir.segmentation, segSubDir, ...
    'cell_anisotropy_lobes_signed') ;

for AorE = 1:2
    close all
    fig = figure('units','centimeters','position',[0,0,figW,figW]) ;
    for lobe = fliplr(1:nLobes)
        % subplot(ceil(nLobes * 0.5), 2, lobe)
        if AorE == 1
            c2t = cos(2*meanQALobeThetas(lobe, :))  ;
            s2t = sin(2*meanQALobeThetas(lobe, :))  ;
            stdt = squeeze(meanQALobeThetaStds(lobe, :)) ;
            stet = squeeze(meanQALobeThetaStes(lobe, :)) ;
            arA = (squeeze(meanQALobeAspects(lobe, :)) - 1) ;
            midline = c2t .* (squeeze(meanQALobeAspects(lobe, :)) - 1) ;
            uncs = sqrt(...
                (c2t .* stdt).^2 + ...
                (2*s2t .* stdt .* arA).^2 ) ;
            errs = sqrt(...
                (c2t .* stet).^2 + ...
                (2*s2t .* stet .* arA).^2 ) ;
        else
            c2t = cos(2*meanQELobeThetas(lobe, :))  ;
            s2t = sin(2*meanQELobeThetas(lobe, :))  ;
            midline = c2t .* (squeeze(meanQELobeAspects(lobe, :)) - 1) ;
            stdt = squeeze(meanQELobeThetaStds(lobe, :)) ;
            stet = squeeze(meanQELobeThetaStes(lobe, :)) ;
            arE = (squeeze(meanQELobeAspects(lobe, :)) - 1) ;
            uncs = sqrt(...
                (c2t .* stdt).^2 + ...
                (2*s2t .* stdt .* arE).^2 ) ;
            errs = sqrt(...
                (c2t .* stet).^2 + ...
                (2*s2t .* stet .* arE).^2 ) ;
        end
        timestamps = timePoints - t0 ;
        if contains(tubi.timeUnits, 'min')
            timestamps = timestamps / 60 ;
        end
        x2 = [timestamps, fliplr(timestamps)] ;
        fill(x2,[midline-abs(uncs), fliplr(midline+abs(uncs))], ...
            colors(lobe, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
            'HandleVisibility', 'off');
        hold on;
        hs{lobe} = errorbar(timestamps, midline, errs, 'color', colors(lobe, :)) ;

        legendentries{lobe} = ['chamber ' num2str(lobe)] ;
    end
    % ylims = ylim() ;
    % ylim([-max(abs(ylims)), max(abs(ylims))])

    % Mark zero line
    plot(timestamps, 0*timestamps, 'k--', 'HandleVisibility','off')
    % Labels
    legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
    ylabel('cell anisotropy, $(a/b-1) \cos 2\theta$', ...
        'interpreter', 'latex')
    if contains(tubi.timeUnits, 'min')
        xlabel('time [hr]', 'interpreter', 'latex')  
    else
        xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')  
    end
    sgtitle('endoderm orientation over time', 'interpreter', 'latex')

    saveas(gcf, [imfn, AEstr{AorE} '.png'])
    saveas(gcf, [imfn, AEstr{AorE} '.pdf'])

end



%% Plot histograms
imfn = fullfile(tubi.dir.segmentation, segSubDir, 'cell_anisotropy_hist') ;

close all
fig = figure('units','centimeters','position',[0,0,7, 7]) ;
colormap(inferno)
subplot(2, 2, 1)
imagesc(timePoints - t0, edges, cos2thetaM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
ylabel('cell orientation, $\cos2\theta$',   'interpreter', 'latex')
subplot(2, 2, 2)
imagesc(timePoints - t0, edges, sin2thetaM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
ylabel('cell orientation, $\sin2\theta$',   'interpreter', 'latex')
subplot(2, 2, 3)
imagesc(timePoints - t0, edgesAR, aspectM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
ylabel('cell aspect ratio',   'interpreter', 'latex')
subplot(4, 2, 6)
cb = colorbar('location', 'south') ;
ylabel(cb, 'probability', 'interpreter', 'latex')
caxis([0, 0.05])
axis off
set(gcf, 'color', 'white')

% export_fig(gcf, [imfn '_alt.png'], '-nocrop', '-r600')
saveas(gcf, [imfn '.pdf']) 
subplot(2, 2, 2)
xlabel('')
ylabel('')
Fim = getFigureFrameCData(2000) ;
imwrite(Fim, [imfn '.png'])

% Raw thetas instead
subplot(2, 2, 1)
imagesc(timePoints - t0, edges, thetaM)
set(gca,'YDir','normal')
caxis([0, 0.05])
xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
ylabel('cell orientation, $\theta$',   'interpreter', 'latex')
subplot(2, 2, 2)
cla ; axis off
saveas(gcf, [imfn '_alt.pdf']) 
Fim = getFigureFrameCData(2000) ;
imwrite(Fim, [imfn '_alt.png'])


%% Aspect ratio distributions, variation of mean shape aspect
imfn = fullfile(tubi.dir.segmentation, segSubDir, 'cell_anisotropy_mratio') ;
clf
colors = define_colors() ;
bluecol = colors(1, :) ;
% shade(timePoints - t0, bndlow, timePoints, bndhigh)
x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
fill(x2, [aratio_low25, fliplr(aratio_high75)], ...
    bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
hold on;
% shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
plot(timePoints - t0, aratio_median, '.-')
plot(timePoints - t0, aratio_mean, '.-')
plot(timePoints - t0, mQAar, '.-')
legend({'25-75\%', 'median $a/b$', ...
    '$\sqrt{\lambda_1^{\langle I \rangle}/\lambda_2^{\langle I \rangle}}$', ...
    '$2||\langle Q\rangle|| + 1$'}, 'interpreter', 'latex')

xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
ylabel('aspect ratio',   'interpreter', 'latex')
title('endoderm orientation over time', 'interpreter', 'latex')

saveas(gcf, [imfn, '.png'])
saveas(gcf, [imfn, '.pdf'])



%% Plot as a function of AP position and time (kymograph)
imfn = fullfile(tubi.dir.segmentation, segSubDir, 'ap_kymographs_qc2t_qs2t') ;
clf
if ~exist(imfn, 'file') || overwriteImages
    subplot(1, 2, 1)
    % imagesc(mid_ap, timePoints-t0, medfilt2(mean_qc2ts, [3, 1])) ;
    imagesc(mid_ap, timestamps, mean_QAc2ts) ;
    caxis([-2, 2])
    colormap(blueblackred)
    % gray out NaNs
    whiteout = double(isnan(mean_QAc2ts)) ;
    hold on; 
    wh = imagesc(mid_ap,timestamps, cat(3, whiteout, whiteout, whiteout)) ;
    set(wh, 'alphaData', whiteout)
    xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
    ylabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
    cb = colorbar('location', 'southOutside') ;
    ylabel(cb, '$\langle(a/b-1)\cos 2\theta\rangle$', ...
        'interpreter', 'latex')
    subplot(1, 2, 2)
    % imagesc(mid_ap, timestamps, medfilt2(mean_qs2ts, [3, 1])) ;
    imagesc(mid_ap, timestamps, mean_QAs2ts) ;
    caxis([-2, 2])
    colormap(blueblackred)
    hold on; 
    wh = imagesc(mid_ap,timestamps, cat(3, whiteout, whiteout, whiteout)) ;
    set(wh, 'alphaData', whiteout)
    xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
    ylabel(['time [' timeunits ']'], 'interpreter', 'latex')
    cb = colorbar('location', 'southOutside') ;
    ylabel(cb, '$\langle(a/b-1)\sin 2\theta \rangle$', ...
        'interpreter', 'latex')
    sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
    saveas(gcf, [imfn, '.png'])
    saveas(gcf, [imfn, '.pdf'])
end

%% Time derivative of filtered image AP position
imfn = fullfile(tubi.dir.segmentation, segSubDir, 'ap_kymographs_dc2t_ds2t') ;
clf
clim = 0.05 ;
if ~exist(imfn, 'file') || overwriteImages
    % cfiltered = medfilt2(mean_qc2ts, [3, 1]) ;
    [~, dc2t] = gradient(mean_QAc2ts, diff(timeunits)) ;
    % sfiltered = medfilt2(mean_qs2ts, [3, 1]) ;
    [~, ds2t] = gradient(mean_QAs2ts, diff(timeunits)) ;
    subplot(1, 2, 1)
    imagesc(mid_ap, timestamps, dc2t) ;
    % gray out NaNs
    whiteout = double(isnan(ds2t)) ;
    hold on; 
    wh = imagesc(mid_ap,timestamps, cat(3, whiteout, whiteout, whiteout)) ;
    set(wh, 'alphaData', whiteout)
    caxis([-clim, clim])
    xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
    ylabel(['time [' timeunits ']'], 'interpreter', 'latex')
    colormap(blueblackred)
    cb = colorbar('location', 'southOutside') ;
    ylabel(cb, ...
        ['$\partial_t \langle(a/b-1) \cos 2\theta\rangle$ [' timeunits '$^{-1}$]'], ...
        'interpreter', 'latex')
    subplot(1, 2, 2)
    imagesc(mid_ap, timestamps, ds2t) ;
    % gray out NaNs
    hold on; 
    wh = imagesc(mid_ap,timestamps, cat(3, whiteout, whiteout, whiteout)) ;
    set(wh, 'alphaData', whiteout)
    caxis([-clim clim])
    xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
    ylabel(['time [' timeunits ']'], 'interpreter', 'latex')
    cb = colorbar('location', 'southOutside') ;
    ylabel(cb,  ...
        ['$\partial_t\langle (a/b-1) \sin 2\theta \rangle$ [' timeunits '$^{-1}$]'], ...
        'interpreter', 'latex')
    sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
    saveas(gcf, [imfn '.png'])
    saveas(gcf, [imfn '.pdf'])
end






% OLD CODE FROM generateCellSegmentationPathlines3D.m for plotting
% imfn = fullfile(tubi.dir.segmentation, 'pathlines', [segSubDir 'cell_anisotropy.png']) ;
% timestamps = timePoints - t0 ;
% if contains(lower(tubi.timeUnits), 'min')
%     timestamps = timestamps / 60 ;
%     timeunits = 'hr' ;
% else
%     timeunits = tubi.timeUnits ;
% end
% xlims = [min(timestamps) - 0.25*min(diff(timestamps)), ...
%         max(timestamps) + 0.25*min(diff(timestamps))] ;
% 
% close all 
% fig = figure('units', 'centimeters', 'position', [0, 0, figW, figH]) ;
% stdStyles = {'_quartiles', '_std'} ;
% for stdStyle = 1:2
%     % shade(timePoints - t0, bndlow, timePoints, bndhigh)
%     x2 = [timestamps, fliplr(timestamps)] ;
%     subplot(1, 2, 1)
%     if stdStyle == 1
%         fill(x2, [mratio_low25, fliplr(mratio_high75)], bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
%     else
%         fill(x2, [mean_mratio - mratio_std, ...
%             fliplr(mean_mratio + mratio_std)], ...
%             bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
%     end
%     hold on;
%     errorbar(timestamps, mean_mratio, mratio_ste, ...
%         '-', 'color', bluecol)
%     xlim(xlims) 
%     xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
%     ylabel('cell aspect ratio, $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
%     title('advected endoderm anisotropy over time', 'interpreter', 'latex')
%     
%     % Cos/Sin(2theta)
%     set(fig,'defaultAxesColorOrder', CScolors);
%     subplot(1, 2, 2)
%     yyaxis left
%     if stdStyle == 1
%         fill(x2, [c2t_low25, fliplr(c2t_high75)], CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
%     else
%         fill(x2, [min(1, max(-1, mc2t - c2t_std)), ...
%             fliplr(min(1, max(-1,mc2t + c2t_std)))],...
%             CScolors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
%     end
%     hold on;
%     % shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
%     errorbar(timestamps, mc2t, ...
%         min(c2t_ste, mc2t+1), ...  % bound from below at -1      
%         min(c2t_ste, 1-mc2t), ... % bound from above at 1
%         '-', 'color', CScolors(1, :))
%     ylim([-1,1])
%     yyaxis right ;
%     if stdStyle == 1
%         fill(x2, [s2t_low25, fliplr(s2t_high75)], CScolors(2, :), 'facealpha', 0.3, 'edgecolor', 'none');
%     else
%         fill(x2, [min(1, max(-1, ms2t - s2t_std)), ...
%             fliplr(min(1, max(-1,ms2t + s2t_std)))], ...
%             CScolors(2, :), 'facealpha', 0.3, 'edgecolor', 'none');
%     end
%     errorbar(timestamps, ms2t, ...
%         min(s2t_ste, ms2t+1), ... % bound from below at -1
%         min(s2t_ste, 1-ms2t), ... % bound from above at 1
%         '-', 'color', CScolors(2, :))
%     ylim([-1,1])
% 
%     yyaxis left
%     ylabel('cell orientation, $\cos 2\theta$',   'interpreter', 'latex')
%     yyaxis right
%     ylabel('cell orientation, $\sin 2\theta$',   'interpreter', 'latex')
%     
%     xlim(xlims)
%     xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
%     title('advected endoderm orientation over time', 'interpreter', 'latex')
% 
%     saveas(gcf, [imfn, stdStyles{stdStyle} '.png'])
%     saveas(gcf, [imfn, stdStyles{stdStyle} '.pdf'])
%     clf
% end
% 
% 
% %% Plot nematic strength and direction for each lobe 
% imfn = fullfile(tubi.dir.segmentation, 'pathlines', 'cell_anisotropy_lobes') ;
% close all
% fig = figure('units','centimeters','position',[0,0,figH, figH], 'visible', 'off') ;
% 
% for lobe = 1:nLobes
%     subplot(ceil(nLobes * 0.5), 2, lobe)
%     midline = squeeze(meanQLobeAspects(lobe, :)) ;
%     uncs = squeeze(meanQLobeAspectStds(lobe, :)) ;
%     timestamps = timePoints - t0 ;
%     if contains(tubi.timeUnits, 'min')
%         timestamps = timestamps / 60 ;
%     end
%     x2 = [timestamps, fliplr(timestamps)] ;
%     fill(x2,[midline-uncs, fliplr(midline+uncs)], ...
%         bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
%     hold on;
%     plot(timestamps, midline, '.-', 'color', bluecol)
%     ylim([1, Inf])
%     
%     yyaxis right
%     % fill(x2, [c2t_low, fliplr(c2t_high)], redcol, 'facealpha', 0.3, 'edgecolor', 'none');
%     hold on;
%     plot(timestamps, mod(meanQLobeThetas(lobe, :), pi)/pi, '.-')
%     % plot(timestamps, sin(2*meanQLobeThetas(lobe, :)), '.-')
%     % 'color', redcol)
%     ylim([0, 1])
%     
%     if mod(lobe, 2) == 1 
%         yyaxis left
%         ylabel(['lobe ' num2str(lobe) ' aspect ratio'],   'interpreter', 'latex')
%     else
%         yyaxis right
%         ylabel('nematic orientation $\theta/\pi$',   'interpreter', 'latex')
%     end
%     
%     % Time label
%     if lobe > nLobes - 2
%         if contains(tubi.timeUnits, 'min')
%             xlabel('time [hr]', 'interpreter', 'latex')  
%         else
%             xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')  
%         end
%     end
% end
% sgtitle('advected endoderm orientation over time', 'interpreter', 'latex')
% saveas(gcf, [imfn '.png'])
% saveas(gcf, [imfn '.pdf'])
% 
% %% Plot nematic strength and direction for each lobe 
% imfn = fullfile(tubi.dir.segmentation, 'pathlines', ...
%     'cell_anisotropy_lobes_signed') ;
% close all
% fig = figure('units','centimeters','position',[0,0,figH, figH]) ;
% 
% for lobe = 1:nLobes
%     % subplot(ceil(nLobes * 0.5), 2, lobe)
%     c2t = cos(2*meanQLobeThetas(lobe, :))  ;
%     midline = 0.5 *c2t .* (squeeze(meanQLobeAspects(lobe, :)) - 1) ;
%     uncs = 0.5 *c2t .* (squeeze(meanQLobeAspectStds(lobe, :)) - 1) ;
%     timestamps = timePoints - t0 ;
%     if contains(tubi.timeUnits, 'min')
%         timestamps = timestamps / 60 ;
%     end
%     x2 = [timestamps, fliplr(timestamps)] ;
%     fill(x2,[midline-abs(uncs), fliplr(midline+abs(uncs))], ...
%         colors(lobe, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
%         'HandleVisibility', 'off');
%     hold on;
%     hs{lobe} = plot(timestamps, midline, '.-', 'color', colors(lobe, :)) ;
%     
%     legendentries{lobe} = ['chamber ' num2str(lobe)] ;
% end
% % ylims = ylim() ;
% % ylim([-max(abs(ylims)), max(abs(ylims))])
% 
% % Mark zero line
% plot(timestamps, 0*timestamps, 'k--', 'HandleVisibility','off')
% % Labels
% legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
% ylabel('elongation, $2||Q|| \cos 2\theta$',   'interpreter', 'latex')
% if contains(tubi.timeUnits, 'min')
%     xlabel('time [hr]', 'interpreter', 'latex')  
% else
%     xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')  
% end
% sgtitle('endoderm orientation over time', 'interpreter', 'latex')
% 
% saveas(gcf, [imfn, '.png'])
% saveas(gcf, [imfn, '.pdf'])
% 
% %% Plot histograms
% imfn = fullfile(tubi.dir.segmentation, 'pathlines', 'cell_anisotropy_hist.png') ;
% clf
% colormap(cividis)
% subplot(2, 2, 1)
% imagesc(timePoints - t0, edges, cos2thetaM)
% set(gca,'YDir','normal')
% caxis([0, 0.05])
% xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
% ylabel('nematic orientation $\cos2\theta$',   'interpreter', 'latex')
% subplot(2, 2, 2)
% imagesc(timePoints - t0, edges, sin2thetaM)
% set(gca,'YDir','normal')
% caxis([0, 0.05])
% xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
% ylabel('nematic orientation $\sin2\theta$',   'interpreter', 'latex')
% subplot(2, 2, 3)
% imagesc(timePoints - t0, edgesAR, aspectM)
% set(gca,'YDir','normal')
% caxis([0, 0.05])
% xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
% ylabel('aspect ratio $\sqrt{I_{1}/I_{2}}$',   'interpreter', 'latex')
% subplot(4, 2, 6)
% cb = colorbar('location', 'south') ;
% ylabel(cb, 'probability', 'interpreter', 'latex')
% caxis([0, 0.05])
% axis off
% saveas(gcf, imfn) 
% 
% %% Aspect ratio distributions, variation of mean shape aspect
% imfn = fullfile(tubi.dir.segmentation, 'pathlines', 'cell_anisotropy_mratio.png') ;
% clf
% colors = define_colors() ;
% bluecol = colors(1, :) ;
% % shade(timePoints - t0, bndlow, timePoints, bndhigh)
% x2 = [timePoints - t0, fliplr(timePoints - t0)] ;
% fill(x2, [mratio_low, fliplr(mratio_high)], ...
%     bluecol, 'facealpha', 0.3, 'edgecolor', 'none');
% hold on;
% % shadedErrorBar(timePoints - t0, mean(y,1),std(y),'lineProps','g');
% plot(timePoints - t0, median_moiratio, '.-')
% plot(timePoints - t0, mean_moiratio, '.-')
% plot(timePoints - t0, mean_mratio, '.-')
% legend({'25-75\%', 'median $\sqrt{I_1/I_2}$', ...
%     '$\sqrt{\lambda_1^{\langle I \rangle}/\lambda_2^{\langle I \rangle}}$', ...
%     '$2||\langle Q\rangle|| + 1$'}, 'interpreter', 'latex')
% 
% xlabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
% ylabel('aspect ratio',   'interpreter', 'latex')
% title('endoderm orientation over time', 'interpreter', 'latex')
% saveas(gcf, imfn)
% 
% 
% 
% %% Plot as a function of AP position and time (kymograph)
% imfn = fullfile(tubi.dir.segmentation, 'pathlines', 'ap_kymographs_c2t_s2t.png') ;
% clf
% subplot(1, 2, 1)
% imagesc(mid_ap, timePoints-t0, medfilt2(mean_c2ts, [3, 1])) ;
% caxis([-10, 10])
% colormap(blueblackred)
% xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
% ylabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
% cb = colorbar('location', 'southOutside') ;
% ylabel(cb, '$\sqrt{\lambda_1/\lambda_2}\cos 2\theta$', ...
%     'interpreter', 'latex')
% subplot(1, 2, 2)
% imagesc(mid_ap, timePoints-t0, medfilt2(mean_s2ts, [3, 1])) ;
% caxis([-10, 10])
% colormap(blueblackred)
% xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
% ylabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
% cb = colorbar('location', 'southOutside') ;
% ylabel(cb, '$\sqrt{\lambda_1/\lambda_2}\sin 2\theta$', ...
%     'interpreter', 'latex')
% sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
% saveas(gcf, imfn)
% 
% %% Time derivative of filtered image AP position
% imfn = fullfile(tubi.dir.segmentation, 'pathlines', 'ap_kymographs_dc2t_ds2t.png') ;
% clf
% cfiltered = medfilt2(mean_c2ts, [3, 1]) ;
% [~, dc2t] = gradient(cfiltered) ;
% sfiltered = medfilt2(mean_s2ts, [3, 1]) ;
% [~, ds2t] = gradient(sfiltered) ;
% subplot(1, 2, 1)
% imagesc(mid_ap, timePoints-t0, imgaussfilt(medfilt2(dc2t, [3, 1]), 0.5)) ;
% caxis([-3, 3])
% xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
% ylabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
% colormap(blueblackred)
% cb = colorbar('location', 'southOutside') ;
% ylabel(cb, '$\partial_t\left(\sqrt{\lambda_1/\lambda_2} \cos 2\theta\right)$', ...
%     'interpreter', 'latex')
% subplot(1, 2, 2)
% imagesc(mid_ap, timePoints-t0, imgaussfilt(medfilt2(ds2t, [3, 1]), 0.5)) ;
% caxis([-3, 3])
% xlabel('ap position, $\zeta/L$', 'interpreter', 'latex')
% ylabel(['time [' tubi.timeUnits ']'], 'interpreter', 'latex')
% cb = colorbar('location', 'southOutside') ;
% ylabel(cb, '$\partial_t\left(\sqrt{\lambda_1/\lambda_2} \sin 2\theta\right)$', ...
%     'interpreter', 'latex')
% sgtitle('cell anisotropy kymographs', 'interpreter', 'latex')
% saveas(gcf, imfn)