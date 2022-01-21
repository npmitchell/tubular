function plotCellDensity(QS, options)
% measureCellDensity(QS, nuclei_or_membrane, options)
%   Measure the cell density over time using relaxed smoothed-mesh 
%   pullbacks.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%
% Returns
% -------
% <none>
%
% NPMitchell 2020

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
if isfield(options, 'doubleCovered')
    doubleCovered = options.doubleCovered ; 
else
    doubleCovered = true;
end

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
imDir1 = fullfile(QS.dir.cellID, 'images') ;
if ~exist(imDir1, 'dir')
    mkdir(imDir1)
end
imDir2 = fullfile(QS.dir.cellID, 'images2d') ;
if ~exist(imDir2, 'dir')
    mkdir(imDir2)
end

[~, ~, ~, xyzlim] = QS.getXYZLims() ;
% Add a bit of axis limit buffer for plotting
xyzlim = xyzlim + [-5, 5] ;
t0 = QS.t0set() ;
blue = QS.plotting.colors(1,:) ;
sky = QS.plotting.colors(6,:) ;

%% Plots
for tp=timePoints
    
    % Check that figures have been saved of the cell density
    datfn = sprintf(QS.fullFileBase.cellID, tp) ;
    figfn1 = fullfile(imDir1, [sprintf(QS.fileBase.cellID, tp) '.png']) ;
    figfn2 = fullfile(imDir2, [sprintf(QS.fileBase.cellID, tp) '.png']) ;
        
    figs_missing = ~exist(figfn1, 'file') || overwrite || ...
        ~exist(figfn2, 'file')  ;
    
    % Load data if needed
    if figs_missing 
        disp(['Loading cells from disk:' datfn])
        load(datfn, 'cells')

        % load image too
        fn = sprintf(QS.fullFileBase.im_r_sme, tp) ;
        im = imread(fn) ;
    end
    
    if ~exist(figfn1, 'file') || overwrite
        % Save image of density
        xx = cells.xyz ;
        vor = cells.vorarea ;
        ff = cells.faces ;
        close all ;
        set(gcf, 'visible', 'off')
        npanel = 6 ;
        subplot(1, npanel, 4:npanel) ;
        trisurf(ff, xx(:, 1), xx(:, 2), xx(:, 3), vor, ...
            'edgecolor', 'none')
        axis equal
        xlim(xyzlim(1, :)) ;
        ylim(xyzlim(2, :)) ;
        zlim(xyzlim(3, :)) ;
        xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
        title(['cell area, $t=$' sprintf('%03d', tp - t0) ' ' QS.timeunits], ...
             'Interpreter', 'Latex')
        % Colorbar
        cb = colorbar('location', 'southoutside') ;
        ylabel(cb, 'cell area [$\mu$m$^2$]', 'Interpreter', 'Latex') ;
        cpos = get(cb, 'position') ;
        cwidth = 0.25 ;
        cxpos = cpos(1) + 0.5*(cpos(3) - cwidth) - 0.015 ;
        set(cb,'position', [cxpos cpos(2) cwidth .02])
        caxis([0, 150])
        view(0, 0)

        % 2D cell density heatmap
        subplot(1, npanel, 1:2)
        imc = imcomplement(im) ;
        imshow(cat(3, imc, imc, imc)); hold on;
        XYtc = cells.tripleCover.XY ;
        disp('Plotting voronoi cells')
        [v, c] = voronoin(XYtc) ;
        for i=1:length(c)
            fill(v(c{i},1), v(c{i},2), cells.tripleCover.vorarea(i),...
                'FaceAlpha', 0.35)
        end
        caxis([0, 150])
        % Limit plot to the single cover
        xlim([0, size(im, 2)])
        if doubleCovered
            ylim([0.25, 0.75] * size(im, 1))
        else
            ylim([0, size(im, 1)])
        end

        % Save the figure
        disp(['Saving ' figfn1])
        saveas(gcf, figfn1)
    end
    
    % 2D image with histograms on sides
    if ~exist(figfn2, 'file') || overwrite
        % 2D cell density heatmap
        close all
        xw = 1000 ;
        yh = 500 ;
        fig = figure('Position', [10 10 xw yh], 'visible', 'off') ;
        imc = imcomplement(im) ;
        imshow(cat(3, imc, imc, imc)); hold on;
        XYtc = cells.tripleCover.XY ;
        disp('Plotting voronoi cells')
        [v, c] = voronoin(XYtc) ;
        for i=1:length(c)
            plot(v(c{i},1), v(c{i},2), '-', 'color', blue)
        end
        % Limit plot to the single cover
        set(fig, 'Position', [10 10 xw yh])
        xlim([0, size(im, 2)])
        if doubleCovered
            ylim([0.25, 0.75] * size(im, 1))
        else
            ylim([0, size(im, 1)])
        end
        
        % Now create density histograms
        uBinID = cells.uBinID ;
        vBinID = cells.vBinID ;
        vorarea = cells.vorarea ;
        % count the number of points at each unique x bin
        counts = accumarray(uBinID, 1.0);  
        % average the voronoi areas that fall into each unique x bin
        sums = accumarray(uBinID, vorarea);
        umeans = sums ./ counts ;
        du = 1.0 ./ umeans ;

        % count the number of points at each unique x bin
        counts = accumarray(vBinID, 1.0);  
        % average the voronoi areas that fall into each unique x bin
        sums = accumarray(vBinID, vorarea);
        vmeans = sums ./ counts ;
        dv = 1.0 ./ vmeans ;
        
        % Set position of ax0
        ax0 = gca ;
        pos = get(ax0, 'position') ;
        set(ax0, 'position', [pos(1) 0.1 pos(3) 0.65]) ;
        
        % Plot density along AP
        dmax = 0.04 ;
        pos = plotboxpos(ax0) ;
        subh = 0.1 ;
        ax2 = axes('position', [pos(1) 0.8 pos(3) subh]) ;
        plot(linspace(0, 1, nU), du) ;
        % xlabel('ap position, $\zeta/L$', 'Interpreter', 'Latex')
        xticks([])
        ylabel('cell density [$\mu$m$^{-2}$]', 'Interpreter', 'Latex') ;
        ylim([0, dmax])
        yticks([0, dmax])
        
        % Plot density along DV
        ax3 = axes('position', ...
            [pos(1)+pos(3)+(0.05*yh/xw), pos(2), yh/xw*subh, pos(4)]) ;
        plot(dv, linspace(0, 1, nV)) ;
        % ylabel('circumferential position, $\phi$',...
        %    'Interpreter', 'Latex', 'position', 'right')
        yticks([])
        xlabel('cell density [$\mu$m$^{-2}$]', 'Interpreter', 'Latex') ;
        set ( ax3, 'xdir', 'reverse' )
        xlim([0, dmax])
        xticks([0, dmax])
        set(fig, 'Position', [10 10 xw yh])
                
        % Save the figure
        % saveas(fig, figfn2)
        export_fig(figfn2, '-nocrop', '-r300', '-transparent')
    end
    
end