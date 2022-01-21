function plotCellDensityKymograph(QS, options)
% plotCellDensityKymograph(QS, options)
%   Plot the cell density over time in reduced dimensionality
%
% Parameters
% ----------
% QS : QuapSlap class instance
% nuclei_or_membrane : str specifier ('nuclei' or 'membrane')
%   whether to use membrane or nuclear method to measure
%   density
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
% 
%
% NPMitchell 2020

eu = 1 ;
ev = 4 ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
if isfield(options, 'preview')
    preview = options.preview ;
else
    preview = false ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
else
    timePoints = QS.xp.fileMeta.timePoints ;
end

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
imDir = fullfile(QS.dir.cellID, 'kymographs') ;
if ~exist(imDir, 'dir')
    mkdir(imDir)
end
[~, ~, ~, xyzlim] = QS.getXYZLims() ;
% Add a bit of axis limit buffer for plotting
xyzlim = xyzlim + [-5, 5] ;
t0 = QS.t0set() ;
QS.getFeatures('folds', 'fold_onset')
colors_noyellow = QS.plotting.colors([1,2,4:7], :) ;

%% Declare figures to create
ufigfn_area = fullfile(imDir, 'area_u.png') ;
vfigfn_area = fullfile(imDir, 'area_v.png') ;
ufigfn_dens = fullfile(imDir, 'density_u.png') ;
vfigfn_dens = fullfile(imDir, 'density_v.png') ;
ufigfn_dens_sm = fullfile(imDir, 'density_u_sm.png') ;
vfigfn_dens_sm = fullfile(imDir, 'density_v_sm.png') ;
ufigfn_dDdt_sm = fullfile(imDir, 'dDdt_u_sm.png') ;
vfigfn_dDdt_sm = fullfile(imDir, 'dDdt_v_sm.png') ;
ufigfn_dDdt_sm_frac = fullfile(imDir, 'dDdt_u_sm_frac.png') ;
vfigfn_dDdt_sm_frac = fullfile(imDir, 'dDdt_v_sm_frac.png') ;
kymoU = zeros(length(timePoints), nU-eu*2) ;
kymoV = zeros(length(timePoints), nV-ev*2) ;
kymoU = inpaint_nans(kymoU) ;
kymoV = inpaint_nans(kymoV) ;

% if ~exist(figfn, 'file') || overwrite
for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    datfn = sprintf(QS.fullFileBase.cellID, tp) ;
    
    % Check that figures have been saved of the cell density
    disp(['Loading cells from disk:' datfn])
    load(datfn, 'cells')
    
    % % Smoothing of voronoi area in space
    % lambda = 0.01 ;
    % vorarea_sm = laplacian_smooth(cells.xyz, cells.faces, 'cotan', [],...
    %                         lambda, 'explicit', cells.vorarea) ;
    % % Check it
    % figure ;
    % subplot(1, 2, 1)
    % trisurf(cells.faces, cells.xyz(:, 1), cells.xyz(:, 2), ...
    %     cells.xyz(:, 3), cells.vorarea);
    % axis equal
    % caxis([0, 100])
    % subplot(1, 2, 2)
    % trisurf(cells.faces, cells.xyz(:, 1), cells.xyz(:, 2), ...
    %     cells.xyz(:, 3), vorarea_sm);
    % axis equal
    % caxis([0, 100])
    % pause(1)
    
    % Unpack single cover
    vorarea = cells.vorarea ;
    % xy = cells.uv ;
    % umax = cells.umax ;
    % vmax = cells.vmax ;
    
    % Binned data
    uBinID = cells.uBinID ;
    vBinID = cells.vBinID ;
    
    % Filter vorarea a bit
    ok_size = vorarea > 20 & vorarea < 175 ;
    inds = ((uBinID > eu & uBinID < nU-eu+1) & ...
        (vBinID > ev & vBinID < nV-ev+1)) & ok_size;
    vorarea = vorarea(inds) ;
    uBinID = uBinID(inds) ;
    vBinID = vBinID(inds) ;

    % count the number of points at each unique x bin
    counts = accumarray(uBinID, 1.0);  
    % average the voronoi areas that fall into each unique x bin
    sums = accumarray(uBinID, vorarea);
    umeans = sums ./ counts ;
    umeans = umeans(eu+1:end) ;
    if length(umeans) < nU-eu*2
        umeans(nU-eu*2) = NaN ;
    end

    % count the number of points at each unique x bin
    counts = accumarray(vBinID, 1.0);  
    % average the voronoi areas that fall into each unique x bin
    sums = accumarray(vBinID, vorarea);
    vmeans = sums ./ counts ;
    vmeans = vmeans(ev+1:end) ;
    if length(vmeans) < nV-ev*2
        umeans(nV-ev*2) = NaN ;
    end

    kymoU(tidx, :) = umeans ;
    kymoV(tidx, :) = vmeans ;
    
end 
uus = linspace(eu/nU, 1-eu/nU, nU-2*eu) ;
vvs = linspace(ev/nV, 1-ev/nV, nV-2*ev) ;

%% Plot the cell areas
for ii = 1:2
    close all
    figure('visible', 'off')
    tps = timePoints - t0 ;
    if ii == 1
        imagesc(uus, tps, kymoU) ;
        hold on;
        for qq=1:length(QS.features.fold_onset)
            fons = QS.features.fold_onset(qq) ;
            ntps = length(tps(fons:end)) ;
            plot(QS.features.folds(fons:fons+ntps-1, qq) / nU, ...
                tps(fons:end), '-', ...
                'color', colors_noyellow(qq, :))
            plot(QS.features.folds(fons, qq) / nU, ...
                tps(fons), QS.plotting.markers{qq}, ...
                'color', colors_noyellow(qq, :))
        end
    else
        imagesc(vvs, tps, kymoV) ;
    end
    
    % Labels
    if ii == 1
        xlabel('ap position, $\zeta/L$', 'Interpreter', 'Latex')
    else
        xlabel('circumferential position, $\phi$',...
            'Interpreter', 'Latex')
    end
    ylabel(['time [' QS.timeunits ']'], 'Interpreter', 'Latex')
    caxis([0, 100]) ;
    cb = colorbar() ;
    ylabel(cb, 'cell area [$\mu$m$^2$]', 'Interpreter', 'Latex') ;
    title('Cell areas, [$\mu$m$^2$]', 'Interpreter', 'Latex')
    if ii == 1
        saveas(gcf, ufigfn_area) ;
    else
        saveas(gcf, vfigfn_area) ;
    end
end


%% Plot the cell densities
for ii = 1:2
    close all
    figure('visible', 'off')
    if ii ==1 
        imagesc(uus, tps, 1 ./ kymoU) ;
        hold on;
        for qq=1:length(QS.features.fold_onset)
            fons = QS.features.fold_onset(qq) ;
            ntps = length(tps(fons:end)) ;
            plot(QS.features.folds(fons:fons+ntps-1, qq) / nU, ...
                tps(fons:end), '-', ...
                'color', colors_noyellow(qq, :))
            plot(QS.features.folds(fons, qq) / nU, ...
                tps(fons), QS.plotting.markers{qq}, ...
                'color', colors_noyellow(qq, :))
        end
    else
        imagesc(vvs, tps, 1 ./ kymoV) ;
    end
    % Labels
    if ii == 1
        xlabel('ap position, $\zeta/L$', 'Interpreter', 'Latex')
    else
        xlabel('circumferential position, $\phi$',...
            'Interpreter', 'Latex')
    end
    ylabel(['time, ' QS.timeunits], 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([0, 0.04])
    ylabel(cb, 'cell density [$\mu$m$^{-2}$]', 'Interpreter', 'Latex') ;
    title('Cell density, [$\mu$m$^{-2}$]', 'Interpreter', 'Latex')
    if ii == 1
        saveas(gcf, ufigfn_dens) ;
    else
        saveas(gcf, vfigfn_dens) ;        
    end
    close all
end

%% Smooth and save
% Prep filter
disp('Building tripulse filter equivalent to tripuls(-0.5:0.1:0.5)')
% Tripulse in time
windowsz_time = 9 ;
tripulse = 0:0.2:1 ;
tripulse = [tripulse, fliplr(tripulse(1:end-1))] ;
tripulse = tripulse ./ sum(tripulse(:)) ;
tripulse = reshape(tripulse, [length(tripulse), 1]) ;
kusm = nanconv(1.0 ./ kymoU, tripulse, 'edge') ;
kvsm = nanconv(1.0 ./ kymoV, tripulse, 'edge') ;
% Smaller tripulse in space
windowsz_space = 3 ;
tripulse = 0:0.5:1 ;
tripulse = [tripulse, fliplr(tripulse(1:end-1))] ;
tripulse = tripulse ./ sum(tripulse(:)) ;
tripulse = reshape(tripulse, [length(tripulse), 1]) ;
kusm = nanconv(kusm', tripulse, 'edge')' ;
kvsm = nanconv(kvsm', tripulse, 'edge')'  ;

% plot it
for ii = 1:2
    figure('visible', 'off')
    if ii == 1
        imagesc(uus, tps, kusm) ;
        hold on;
        for qq=1:length(QS.features.fold_onset)
            fons = QS.features.fold_onset(qq) ;
            ntps = length(tps(fons:end)) ;
            plot(QS.features.folds(fons:fons+ntps-1, qq) / nU, ...
                tps(fons:end), '-', ...
                'color', colors_noyellow(qq, :))
            plot(QS.features.folds(fons, qq) / nU, ...
                tps(fons), QS.plotting.markers{qq}, ...
                'color', colors_noyellow(qq, :))
        end
    else
        imagesc(vvs, tps, kvsm) ;        
    end
    
    % Labels
    if ii == 1
        xlabel('ap position, $\zeta/L$', 'Interpreter', 'Latex')
    else
        xlabel('circumferential position, $\phi$',...
            'Interpreter', 'Latex')
    end
    ylabel(['time, ' QS.timeunits], 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([0, 0.04])
    ylabel(cb, 'cell density [$\mu$m$^{-2}$]', 'Interpreter', 'Latex') ;
    title('Cell density, [$\mu$m$^{-2}$]', 'Interpreter', 'Latex')
    if ii == 1
        saveas(gcf, ufigfn_dens_sm) ;
    else
        saveas(gcf, vfigfn_dens_sm) ;
    end
    close all
end

%% Take time derivative
[~, uDdt] = gradient(kusm) ;
[~, vDdt] = gradient(kvsm) ;
for ii = 1:2
    close all
    figure('visible', 'off') ;
    if ii == 1
        imagesc(uus, tps, uDdt) 
        colormap bwr
        hold on;
        for qq=1:length(QS.features.fold_onset)
            fons = QS.features.fold_onset(qq) ;
            ntps = length(tps(fons:end)) ;
            plot(QS.features.folds(fons:fons+ntps-1, qq) / nU, ...
                tps(fons:end), '-', ...
                'color', colors_noyellow(qq, :))
            plot(QS.features.folds(fons, qq) / nU, ...
                tps(fons), QS.plotting.markers{qq}, ...
                'color', colors_noyellow(qq, :))
        end
    else
        imagesc(vvs, tps, vDdt) 
        colormap bwr
    end
    
    % Labels
    if ii == 1
        xlabel('ap position, $\zeta/L$', 'Interpreter', 'Latex')
    else
        xlabel('circumferential position, $\phi$',...
            'Interpreter', 'Latex')
    end
    ylabel(['time, ' QS.timeunits], 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([-0.002, 0.002])
    ylabel(cb, 'cell density change, $\partial_t n$ [$\mu$m$^{-2}$]', ...
        'Interpreter', 'Latex') ;
    title('Cell density change, $\partial_t n$ [$\mu$m$^{-2}$]',...
        'Interpreter', 'Latex')
    if ii == 1
        saveas(gcf, ufigfn_dDdt_sm) ;
    else
        saveas(gcf, vfigfn_dDdt_sm) ;
    end
    close all
end

%% Plot fractional change (like gdot)
for ii = 1:2
    close all
    figure('visible', 'off')
    if ii == 1
        imagesc(uus, tps, uDdt ./ kusm) 
        colormap bwr
        hold on;
        for qq=1:length(QS.features.fold_onset)
            fons = QS.features.fold_onset(qq) ;
            ntps = length(tps(fons:end)) ;
            plot(QS.features.folds(fons:fons+ntps-1, qq) / nU, ...
                tps(fons:end), '-', ...
                'color', colors_noyellow(qq, :))
            plot(QS.features.folds(fons, qq) / nU, ...
                tps(fons), QS.plotting.markers{qq}, ...
                'color', colors_noyellow(qq, :))
        end
    else
        imagesc(vvs, tps, vDdt ./ kvsm) 
        colormap bwr
    end
    
    % Labels
    if ii == 1
        xlabel('ap position, $\zeta/L$', 'Interpreter', 'Latex')
    else
        xlabel('circumferential position, $\phi$',...
            'Interpreter', 'Latex')
    end
    ylabel(['time, ' QS.timeunits], 'Interpreter', 'Latex')
    cb = colorbar() ;
    caxis([-0.05, 0.05])
    ylabel(cb, 'fractional cell density change, $\partial_t n /n$ [$\mu$m$^{-2}$]', ...
        'Interpreter', 'Latex') ;
    title('Fractional cell density change, $\partial_t n /n$ [$\mu$m$^{-2}$]',...
        'Interpreter', 'Latex')
    if ii == 1
        saveas(gcf, ufigfn_dDdt_sm_frac) ;
    else
        saveas(gcf, vfigfn_dDdt_sm_frac) ;
    end
    close all
end

%%%%%%%%%%%%%%%%%%%%%%
%% Save raw kymographs
%%%%%%%%%%%%%%%%%%%%%%
kymographs.area_u = kymoU ;
kymographs.area_v = kymoV ;
kymographs.density_u = kusm ;
kymographs.density_dt_u = uDdt ;
kymographs.windowsz_space = windowsz_space ;
kymographs.windowsz_time = windowsz_time ;
save(fullfile(imDir, 'kymographs.mat'), 'kymographs')

