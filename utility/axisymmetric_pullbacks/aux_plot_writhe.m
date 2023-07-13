function aux_plot_writhe(tubi, clines_resampled, ...
    Wr, Wr_density, dWr, Length_t, wrfigdir, area_volume_fn, fold_onset, Wr_style, ...
    clineDVhoopBase, cylinderMeshCleanBase, rot, trans, resolution, ...
    flipy, omit_endpts, black_figs)
%AUX_PLOT_WRITHE(Wr, Wr_density, dWr, Length_t)
% auxiliary function for plotting writhe over time along with geometric
% properties of length, area, and volume

%
timePoints = tubi.xp.fileMeta.timePoints ;
timeInterval = tubi.timeInterval ; 
timeUnits = tubi.timeUnits ; 
spaceUnits = tubi.spaceUnits ;
t0 = tubi.t0 ;
nU = tubi.nU ;
[~,~,~,xyzlim] = tubi.getXYZLims ;

% figure parameters
xwidth = 32 ; % cm
ywidth = 20 ; % cm

% Unpack input structs
if strcmp(Wr_style, 'polar')
    Wrc = Wr.polar ;
    wr_densities = Wr_density.polar ;
    dwr = dWr.polar ;
elseif strcmp(Wr_style, 'Levitt')
    Wrc = Wr.Levitt ;
    wr_densities = Wr_density.Levitt ;
    dwr = dWr.Levitt ;
elseif strcmp(Wr_style, 'Gauss')
    Wrc = Wr.Gauss ;
    wr_densities = Wr_density.Gauss ;
    dwr = dWr.Gauss ;
end

Wrp = Wr.polar ;
Wrpseg1 = Wr.polar_seg1 ;
WrL = Wr.Levitt ;
WrG = Wr.Gauss ;
lengths = Length_t.lengths ;
dl = Length_t.dl ;


%% Save simple figures
% Save chirality as a figure -- todo
% close all
% fig = figure('Visible', 'Off') ;
% plot(timepoints, Chirality) ;
% xlabel('time [' timeUnits ']', 'Interpreter', 'Latex')
% ylabel('Chirality, $\Delta \phi$', 'Interpreter', 'Latex')
% title('Chirality over time', 'Interpreter', 'Latex')
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);        
% saveas(fig, fullfile(figoutdir, 'chirality_vs_time_DVhoop.pdf'))
% saveas(fig, fullfile(figoutdir, 'chirality_vs_time_DVhoop.png'))

% Save writhe as a figure
close all
fig = figure('Visible', 'Off') ;
plot((timePoints-t0)*timeInterval, Wrp) ;
hold on 
plot((timePoints-t0)*timeInterval, Wrpseg1) ;
plot((timePoints-t0)*timeInterval, WrL, '--') ;
plot((timePoints-t0)*timeInterval, WrG, ':') ;
xlabel(['time [' timeUnits ']'], 'Interpreter', 'Latex')
ylabel('Writhe')
title('Writhe over time', 'Interpreter', 'Latex')
legend({'polar', 'polar (segment 1 only)', 'Levitt', 'Gauss'}, 'location', 'best')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);        

outfn = fullfile(wrfigdir, ['writhe_' Wr_style '_vs_time_comparison_DVhoop']) ;
disp(['Saving figure to: ' outfn ' as pdf/png'])
saveas(fig, [outfn '.pdf'])
saveas(fig, [outfn '.png'])
set(gcf, 'visible', 'on')

%% Save writhe as a figure with length, surfarea, volume
% Do it again in stages --> the last build isn't finished right now
% load 'aas', 'vvs', 'dt' 
load(area_volume_fn, 'aas', 'vvs')
% load the folding timepoints
try
    t2 = fold_onset(2) ;
    t1 = fold_onset(1) ;
    t3 = fold_onset(3) ;
    indicateFolds = true ;
catch
    indicateFolds = false ;
end

ind = find(timePoints == t0) ;
[colors, ~] = define_colors ;

for ii = 1:4
    disp(['ii = ' num2str(ii)])
    close all
    fig = figure('units', 'centimeters', ...
            'outerposition', [0 0 12 12], 'Visible', 'Off') ;
    axis square
    pbaspect([1.5 1 1])
    yyaxis left
    
    % fill in boxes
    % if ii == 5
    %     % To make fills, get xlimits from a handle that we clear
    %     vh = plot(timepoints - t0, vvs / vvs(ind)) ;
    %     xlims = xlim ;
    %     cla
    % 
    %     % Get the limits to shade in the background
    %     ymax = 3.2 ;
    %     twrithe = 40 ;
    %     tend = 70 ;
    %     xx1 = [xlims(1) 0, 0, xlims(1)] ;
    %     xx2 = [0, twrithe, twrithe, 0] ;
    %     xx3 = [twrithe, tend, tend, twrithe] ;
    %     yfill = [0, ymax, ymax, 0 ] ;
    %     hf1 = fill(xx1, yfill, color1, 'facealpha', 0.5);
    %     hold on;
    %     hf2 = fill(xx2, yfill, color2, 'facealpha', 0.5);
    %     hf3 = fill(xx3, yfill, color3, 'facealpha', 0.5);
    % end
    
    vh = plot((timePoints - t0)* timeInterval, vvs / vvs(ind),...
        'color', colors(1, :), 'linewidth', 2) ;
    ylabel('$V$', 'Interpreter', 'Latex')
    hold on ;
    if ii > 1
        ah = plot((timePoints - t0)* timeInterval, aas / aas(ind),...
        'color', colors(1, :), 'linewidth', 2) ;
        ylabel('$V$, $A$', 'Interpreter', 'Latex')
    end
    if ii > 2
        lh = plot((timePoints - t0)* timeInterval, lengths / lengths(ind), ...
            'color', colors(1, :), 'linewidth', 2) ;
        ylabel('$V$, $A$, $L$', 'Interpreter', 'Latex')
    end
    
    if ii == 1
        % writhe on right
        h2 = plot([1], [1]) ;
        h3 = plot([1], [1]) ;
        legend({'volume', ' ', ' '}, ...
            'location', 'northwest', 'AutoUpdate', 'off')
    elseif ii == 2
        h3 = plot([1], [1]) ;
        legend({'volume', 'area', ''}, ...
            'location', 'northwest', 'AutoUpdate', 'off')
    elseif ii == 3
        legend({'volume', 'area', 'length'}, ...
            'location', 'northwest', 'AutoUpdate', 'off')
    elseif ii > 3
        % writhe on right
        yyaxis right
        wh = plot((timePoints - t0)* timeInterval, Wrc, 'linewidth', 2) ;
        xlims = get(gca, 'xlim') ;
        ylims = get(gca, 'ylim') ;
        ylim([ylims(1), ylims(2) + 0.8])
        ylabel('Writhe, $Wr$', 'Interpreter', 'Latex')
        legend({'volume', 'area', 'length'}, ...
            'location', 'northwest', 'AutoUpdate', 'off')
    end
    % title('Gut dynamics', 'Interpreter', 'Latex')
    % Plot folding events
    if indicateFolds
        if ii > 3
            plot((t1 - t0)*timeInterval, Wrc(t1-t0+ind), 'ks') ;
            plot((t3 - t0)*timeInterval, Wrc(t3-t0+ind), 'k^') ;
            plot(0, Wrc(ind), 'ko') ;
        end

        yyaxis left
        plot((t1 - t0)*timeInterval, vvs(t1-t0+ind)/vvs(ind), 'ks') ;
        plot((t3 - t0)*timeInterval, vvs(t3-t0+ind)/vvs(ind), 'k^') ;
        plot(0, 1, 'ko') ;
        if ii > 1
            plot((t1 - t0)*timeInterval, aas(t1-t0+ind)/aas(ind), 'ks') ;
            plot((t3 - t0)*timeInterval, aas(t3-t0+ind)/aas(ind), 'k^') ;
            plot(0, 1, 'ko') ;
        end
        if ii > 2
            plot((t1 - t0)*timeInterval, lengths(t1-t0+ind)/lengths(ind), 'ks') ;
            plot((t3 - t0)*timeInterval, lengths(t3-t0+ind)/lengths(ind), 'k^') ;
            plot(0, 1, 'ko') ;
        end
    end
    
    if ii < 4
        yyaxis right
        yticks([])
    end
    
    yyaxis left
    ylimits = ylim ;
    if ii < 3
        ymax = 1.8 ;
    else
        ymax = 3.2 ;
    end
    if ylimits(2) < ymax
        ylimits(2) = ymax ;
    end
    ylim([0 ylimits(2)])
    
    % Set limits
    xlabel(['time [' timeUnits ']'], 'Interpreter', 'Latex')
    set(gca, 'Position', [0.25, 0.25, 0.6, 0.6])
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);   
    set(gca, 'color', 'w', 'xcol', 'w', 'ycol', 'w')
    set(fig, 'color', 'k')
    
    if ii > 3
        fn = fullfile(wrfigdir, ['writhe_dynamics_DVhoop_build' num2str(ii) '_' Wr_style '.png']) ;
    else
        fn = fullfile(wrfigdir, ['writhe_dynamics_DVhoop_build' num2str(ii) '.png']) ;        
    end
    disp(['Saving figure ' fn])
    % saveas(gcf, fn)
    export_fig(fn, '-r300', '-nocrop')
end

%% Save writhe as a figure with length, surfarea, volume
% load 'aas', 'vvs', 'dt'
load(area_volume_fn, 'aas', 'vvs', 'dv', 'da')
% load the folding timepoints
ind = find(timePoints == t0) ;
ylims_derivs = [-0.045, 0.045];
close all
fig = figure('Visible', 'On') ;
s1 = subplot(2, 1, 1) ;
hold on;
yyaxis left
vh = plot((timePoints - t0)*timeInterval, vvs / vvs(ind)) ;
ah = plot((timePoints - t0)*timeInterval, aas / aas(ind)) ;
lh = plot((timePoints - t0)*timeInterval, lengths / lengths(ind)) ;
ylabel('$V$, $A$, $L$', 'Interpreter', 'Latex')
% writhe on right
yyaxis right
wh = plot((timePoints - t0)*timeInterval, Wrc ) ;
xlims = get(gca, 'xlim') ; 
ylabel('Writhe, $Wr$', 'Interpreter', 'Latex')
legend({'volume', 'area', 'length'}, ...
    'location', 'northwest', 'AutoUpdate', 'off')
title('Gut dynamics', 'Interpreter', 'Latex')
% Plot folding events
if indicateFolds
    plot((t1 - t0)*timeInterval, Wrc(t1-t0+ind), 'ks') ;
    plot((t3 - t0)*timeInterval, Wrc(t3-t0+ind), 'k^') ;
end
plot(0, Wrc(ind), 'ko') ;
yyaxis left
if indicateFolds
    plot((t1 - t0)*timeInterval, aas(t1-t0+ind)/aas(ind), 'ks') ;
    plot((t3 - t0)*timeInterval, aas(t3-t0+ind)/aas(ind), 'k^') ;
    plot(0, 1, 'ko') ;
    plot((t1 - t0)*timeInterval, vvs(t1-t0+ind)/vvs(ind), 'ks') ;
    plot((t3 - t0)*timeInterval, vvs(t3-t0+ind)/vvs(ind), 'k^') ;
    plot(0, 1, 'ko') ;
    plot((t1 - t0)*timeInterval, lengths(t1-t0+ind)/lengths(ind), 'ks') ;
    plot((t3 - t0)*timeInterval, lengths(t3-t0+ind)/lengths(ind), 'k^') ;
    plot(0, 1, 'ko') ;
end


% Second plot below
s2 = subplot(2, 1, 2) ;
hold on;
% Plot derivatives
windowSize = 7; 
di = (windowSize:length(dwr)-windowSize) ;
% vcolor = get(vh, 'color') ;
% acolor = get(ah, 'color') ;
% lcolor = get(lh, 'color') ;
% wcolor = get(wh, 'color') ;
yyaxis left
plot((timePoints(di) - t0)*timeInterval, dv(di) / vvs(ind)) ; %, 'Color', vcolor) ;
plot((timePoints(di) - t0)*timeInterval, da(di) / aas(ind)) ; %, 'Color', acolor) ;
plot((timePoints(di) - t0)*timeInterval, dl(di) / lengths(ind)) ; %, 'Color', lcolor) ;
set(gca, 'ylim', ylims_derivs) ;
ylabel('$\partial_t V / V_0$, $\partial_t A / A_0$, $\partial_t L/ L_0$',...
    'Interpreter', 'Latex')

yyaxis right
plot((timePoints(di) - t0)*timeInterval, dwr(di)) ; %, 'Color', wcolor) ;
ylabel('$\partial_t Wr$', 'Interpreter', 'Latex')
set(gca, 'ylim', ylims_derivs) ;

% Set limits
set(s1, 'xlim', xlims)
set(s2, 'xlim', xlims)
xlabel(['time [' timeUnits ']'], 'Interpreter', 'Latex')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 xwidth 2*ywidth]);   
outfn = fullfile(wrfigdir, ['writhe_dynamics_' Wr_style '_DVhoop']) ;
disp(['Saving figure ' outfn])
saveas(fig, [outfn '.pdf'])
saveas(fig, [outfn '.png'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% similar figure  but one panel for everything %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
fig = figure('Visible', 'Off') ;
hold on;
toff = t0 - timePoints(1) ;
disp(['time offset is ' num2str(toff)])
vh = plot((timePoints - t0)*timeInterval, vvs / vvs(ind)) ;
ah = plot((timePoints - t0)*timeInterval, aas / aas(ind)) ;
lh = plot((timePoints - t0)*timeInterval, lengths / lengths(ind)) ;
wh = plot((timePoints - t0)*timeInterval, Wrc ) ;
vcolor = get(vh, 'color') ;
acolor = get(ah, 'color') ;
xlims = get(gca, 'xlim') ; 
ylims = get(gca, 'ylim') ; 
plot((timePoints - t0)*timeInterval, dv / vvs(ind) * 100, '--', 'Color', vcolor) ;
plot((timePoints - t0)*timeInterval, da / aas(ind) * 100, '--', 'Color', acolor) ;
plot(t0, Wrc(max(toff, 1)), 'k.') ;
if indicateFolds
    plot(t1 - t0, Wrc(max(t1 - t0 + toff, 1)), 'ko') ;
    plot(t3 - t0, Wrc(max(t3 - t0 + toff, 1)), 'ks') ;
end
set(gca, 'xlim', xlims)
set(gca, 'ylim', ylims)
xlabel(['time [' timeUnits ']'], 'Interpreter', 'Latex')
ylabel('Volume, Area, Length, \& Writhe', 'Interpreter', 'Latex')
title('Gut dynamics', 'Interpreter', 'Latex')
legend({'volume', 'area', 'length'}, 'location', 'northwest')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);  
outfn = fullfile(wrfigdir, ['writhe_dynamics2_DVhoop_' Wr_style]) ;
disp(['Saving figure to ' outfn])
saveas(fig, [outfn '.pdf'])
saveas(fig, [outfn '.png'])




%% Save figures of writhe density
fig = figure('Visible', 'Off') ;
for ii=1:length(wr_densities)
    t = timePoints(ii) ;
    xpt = clines_resampled{ii}.xyzp(:, 1) ;
    ypt = clines_resampled{ii}.xyzp(:, 2) ;
    zpt = clines_resampled{ii}.xyzp(:, 3) ;    
    
    if ~strcmp(Wr_style, 'polar')
        xpt = (xpt(1:end-1) + xpt(2:end)) * 0.5 ;
        ypt = (ypt(1:end-1) + ypt(2:end)) * 0.5 ;
        zpt = (zpt(1:end-1) + zpt(2:end)) * 0.5 ;
    end
    
    timepointstr = sprintf('%06d', timePoints(ii));
    if length(omit_endpts) == 1
        keep = omit_endpts:(length(xpt) - omit_endpts) ;
    else
        keep = omit_endpts(1):(length(xpt) - omit_endpts(2)) ;
    end
    % scatter3(xpt(keep), ypt(keep), zpt(keep), 10, wr_densities{ii}, 'fill')

    % Interpolate to have a thick curve
    tmpx = linspace(0, 1, length(xpt(keep))) ;
    ix = linspace(0, 1, 1000) ;
    xi = interp1(tmpx, xpt(keep), ix) ;
    yi = interp1(tmpx, ypt(keep), ix) ;
    zi = interp1(tmpx, zpt(keep), ix) ;
    wi = interp1(tmpx, wr_densities{ii}, ix) ;

    scatter3(xi, yi, zi, 10, wi, 'fill')
    hold on; 
    % Load the mesh
    if black_figs 
        tcolor = 'w' ;
        falph = 0.1 ;
    else
        tcolor = 'k' ;
        falph = 0.05 ;
    end
    mesh = read_ply_mod(sprintf(cylinderMeshCleanBase, t)) ;
    mesh.v = ((rot * mesh.v')' + trans) * resolution ;
    if flipy
        mesh.v(:, 2) = -mesh.v(:, 2) ;
    end
    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        'EdgeColor', 'none', 'FaceAlpha', falph, 'FaceColor', tcolor)
    
    
    % % Add other points of centerline
    % fn = sprintf(clineDVhoopBase, t) ;
    % disp(['t = ' num2str(t)])
    % disp(['Loading DVhoop centerline from ' fn])
    % load(fn, 'mss', 'mcline', 'avgpts')
    % hold on; plot3(mcline(:, 1), mcline(:, 2), mcline(:, 3), 'k--')
    
    
    % Note: must have bwr colormap defined
    colormap(bwr)
    
    if strcmp(Wr_style, 'polar')
        caxis([-.02, .02])
    else
        % caxis([-.002, .002])
        caxis([-0.003, 0.003])
    end
    if black_figs
        cb = colorbar('eastOutside', 'color', 'w') ;
    else
        cb = colorbar('eastOutside') ;
    end
    set(cb, 'units', 'normalized', 'position', [.85, .3, 0.05, 0.4])
    set(get(cb, 'label'), 'string', 'Writhe density')
    
    thandle = title(['$t=$' sprintf('%03d', (t - t0)*timeInterval) ' ' timeUnits], ...
        'Interpreter', 'latex') ; 
    set(get(gca,'title'), 'Position',[100, 0, 100 ])
    % set(thandle, 'position', get(thandle, 'position') - [0, 0, 10])
    if black_figs 
        set(thandle, 'color', 'w')
    end
    xlabel(['x [' spaceUnits ']']) ;
    ylabel(['y [' spaceUnits ']']) ;
    zlabel(['z [' spaceUnits ']']) ;
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); 
    axis equal
    xlim([xyzlim(1, 1), xyzlim(1, 2)])
    ylim([xyzlim(2, 1), xyzlim(2, 2)])
    zlim([xyzlim(3, 1), xyzlim(3, 2)])  
    % axis off ;
    grid off ;
    if black_figs
        set(gcf, 'color', 'k')
        set(gca, 'color', 'k', 'xcol', 'k', 'ycol', 'k', 'zcol', 'k')
    else
        set(gca, 'color', 'w')
        set(gcf, 'color', 'w')
    end
    outfn = fullfile(wrfigdir, ['writhe_densities_DVhoop_' Wr_style '_' timepointstr '.png']) ;
    disp(['Saving figure to ' outfn])
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);  
    export_fig(outfn, '-r150', '-nocrop')
    clf
end






%% Save figures of writhe density, with SURFACE colored by writhe density
close all
for ii=1:length(wr_densities)
    t = timePoints(ii) ;
    tubi.setTime(t) ;

    xpt = clines_resampled{ii}.xyzp(:, 1) ;
    timepointstr = sprintf('%06d', timePoints(ii));
    if length(omit_endpts) == 1
        keep = (omit_endpts+1):(length(xpt) - omit_endpts) ;
    else
        keep = (omit_endpts(1)+1):(length(xpt) - omit_endpts(2)) ;
    end

    % duplicate Wr_density values to cover mesh surface
    wrd = wr_densities{ii} ;

    % paint on wrd onto mesh
    mesh = tubi.getCurrentSPCutMeshSmRSC ;
    uv = mesh.u ;
    uv(:, 1) = uv(:, 1) / max(uv(:, 1)) ;
    wrd = interp1(uv(keep, 1), wrd, uv(:, 1), 'linear', 'extrap') ;

    % Load the mesh
    %subplot(2,1,1)
   %  %trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
   %      'EdgeColor', 'none', 'FaceVertexCData', wrd)
   % axis equal
   % view([0,0])
   % 
   %  subplot(2,1,2)
    trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
        'EdgeColor', 'none', 'FaceVertexCData', wrd, 'faceAlpha', 0.5)
   axis equal
   % view([180,0])

    
    % Note: must have bwr colormap defined
    colormap(bwr)
    
    if strcmp(Wr_style, 'polar')
        caxis([-.02, .02])
    else
        caxis([-.002, .002])
    end
    if black_figs
        cb = colorbar('eastOutside', 'color', 'w') ;
    else
        cb = colorbar('eastOutside') ;
    end
    set(cb, 'units', 'normalized', 'position', [.85, .3, 0.05, 0.4])
    set(get(cb, 'label'), 'string', 'Writhe density')
    
    thandle = title(['$t=$' sprintf('%03d', (t - t0)*timeInterval) ' ' timeUnits], ...
        'Interpreter', 'latex') ; 
    set(get(gca,'title'), 'Position',[100, 0, 100 ])
    % set(thandle, 'position', get(thandle, 'position') - [0, 0, 10])
    if black_figs 
        set(thandle, 'color', 'w')
    end
    xlabel(['x [' spaceUnits ']']) ;
    ylabel(['y [' spaceUnits ']']) ;
    zlabel(['z [' spaceUnits ']']) ;
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); 
    axis equal
    xlim([xyzlim(1, 1), xyzlim(1, 2)])
    ylim([xyzlim(2, 1), xyzlim(2, 2)])
    zlim([xyzlim(3, 1), xyzlim(3, 2)])  
    % axis off ;
    grid off ;
    if black_figs
        set(gcf, 'color', 'k')
        set(gca, 'color', 'k', 'xcol', 'k', 'ycol', 'k', 'zcol', 'k')
    else
        set(gca, 'color', 'w')
    end
    outfn = fullfile(wrfigdir, ['writhe_densities_surfaceDVhoop_' Wr_style '_' timepointstr '.png']) ;
    disp(['Saving figure to ' outfn])
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 xwidth ywidth]);  
    export_fig(outfn, '-r150', '-nocrop')
    clf
end

