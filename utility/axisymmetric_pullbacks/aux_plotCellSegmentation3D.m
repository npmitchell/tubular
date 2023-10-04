function aux_plotCellSegmentation3D(QS, tp, seg3d, imdir, overwrite,...
    xyzlims, vertexBasedTiling, mesh4occlusion)
% Plot cell segmentation on tissue in 3D
%
%


% Parse options
occlude = true ;
occludeStr = '' ;
if nargin > 7 
    if ~isempty(mesh4occlusion)
        occlude = true ;
        occludeStr = 'Occluded';
        metafn = fullfile(QS.dir.texturePatchIm, 'overlays', 'metadat.mat') ;
        load(metafn, 'metadat', 'Options')
                    
        disp('Overriding default xyzlims with supplied metadata')
        xyzlims = metadat.xyzlim ;
    end
end

% imWpix = 1680 ; % 4*420 ;
% imHpix = 1238 ; % 4*570 ;
imW = 1260 ;
imH = 787 ;
xwidth = 16 ;
ywidth = 10 ;
areaMaxClim = 85;


%% Draw cells colored by area
t0 = QS.t0() ;
titlestr = ['$t=$' sprintf('%03d', tp-t0) ' ' QS.timeUnits] ;

% Easiest way is to triangulate all polygons using centroids
% This is fine if the cells are all convex
% see also: polygonsToPatchTriangles3D()
faces = seg3d.cdat.polygons ;
keep = seg3d.statistics.keep ;

areas = seg3d.qualities.areas ; 
ang1 = seg3d.qualities.ang1 ; 
% ang2 = seg3d.qualities.ang2 ; 
mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
boa2 = seg3d.qualities.moment1 ./ seg3d.qualities.moment2 ;
moinertia = seg3d.qualities.mInertia ;
c3d = seg3d.vdat.xyzrs ;
cellCntrd = seg3d.cdat.centroids_3d ;
ars = sqrt(mratio) ;

opts = struct() ;
opts.keep = keep ;
opts.centroids = cellCntrd ;
opts.vertexBasedTiling = vertexBasedTiling ;
if ~vertexBasedTiling
    % If cellIDs is passed as cell array, convert to linear
    if iscell(faces)
        cellIDs = zeros(length(faces), 1) ;
        dmyk = 1 ;
        for qq = 1:length(faces)
            for pp = 1:length(faces{qq})
                cellIDs(dmyk) = qq ;
                dmyk = dmyk + 1 ;
            end
        end
    end
    [ff, vv, faceMemberIDs] = ...
        polygonsToPatchTriangles3D(c3d, cellIDs, opts)  ;
else
    [ff, vv, faceMemberIDs] = ...
        polygonsToPatchTriangles3D(c3d, faces, opts)  ;
end


close all
areaV = areas(faceMemberIDs) ;
ang1V = ang1(faceMemberIDs) ;
% IxxV = moinertia(faceMemberIDs, 1) ;
% IyyV = moinertia(faceMemberIDs, 3) ;
aa = sqrt(mratio) ;
aaV = aa(faceMemberIDs) ;
QAxx = (sqrt(mratio) - 1) .* cos(2*ang1) ;
QAxxV = QAxx(faceMemberIDs) ;
ee =  sqrt(1-boa2) ;
eeV = ee(faceMemberIDs) ;
QExx = sqrt(1-boa2) .* cos(2*ang1) ;
QExxV = QExx(faceMemberIDs) ;

%% dirs
if ~exist(fullfile(imdir, ['areas' occludeStr]), 'dir')
    mkdir(fullfile(imdir, ['areas' occludeStr]))
end
if ~exist(fullfile(imdir, ['aspectRatio' occludeStr]), 'dir')
    mkdir(fullfile(imdir, ['aspectRatio' occludeStr]))
end
if ~exist(fullfile(imdir, ['anisotropyQ' occludeStr]), 'dir')
    mkdir(fullfile(imdir, ['anisotropyQ' occludeStr]))
end
if ~exist(fullfile(imdir, ['eccentricity' occludeStr]), 'dir')
    mkdir(fullfile(imdir, ['eccentricity' occludeStr]))
end
if ~exist(fullfile(imdir, ['eccentricityQ' occludeStr]), 'dir')
    mkdir(fullfile(imdir, ['eccentricityQ' occludeStr]))
end
if ~exist(fullfile(imdir, ['order' occludeStr]), 'dir')
    mkdir(fullfile(imdir, ['order' occludeStr]))
end
if ~exist(fullfile(imdir, ['anisotropyQr' occludeStr]), 'dir')
    mkdir(fullfile(imdir, ['anisotropyQr' occludeStr]))
end


%% Color segmentation by area
areaimfn = fullfile(imdir, ['areas' occludeStr], sprintf('cellseg3d_area_%06d', tp)) ;

%% Nematic anisotropy tensor WITH BONDS       
aimfn = fullfile(imdir, ['aspectRatio' occludeStr], sprintf('cellseg3d_aspectRatioBonds_%06d', tp)) ;
aQimfn = fullfile(imdir, ['anisotropyQ' occludeStr],sprintf('cellseg3d_anisotropyQBonds_%06d', tp)) ;
aQimfn_reverse = fullfile(imdir, ['anisotropyQr' occludeStr],sprintf('cellseg3d_anisotropyQreverseBonds_%06d', tp)) ;


%% Nematic eccentricity tensor WITH BONDS    
eimfn = fullfile(imdir, ['eccentricity' occludeStr], sprintf('cellseg3d_eccentricityBonds_%06d', tp)) ;
eQimfn = fullfile(imdir, ['eccentricityQ' occludeStr], sprintf('cellseg3d_eccentricityQBonds_%06d', tp)) ;

%% Nematic order  WITH BONDS    
oimfn = fullfile(imdir, ['order' occludeStr], ...
    sprintf('cellseg3d_orderBonds_%06d', tp)) ;

imfns = {areaimfn, aimfn, aQimfn, eimfn, eQimfn, oimfn, aQimfn_reverse} ;
colorVs = {areaV, aaV, QAxxV, eeV, QExxV, cos(2*ang1V(:)), -QAxxV} ;

for qq = 1:length(imfns)
    for viewID = 0:1
        if viewID == 0
            imfn = [imfns{qq} '.png'];
        else
            imfn = [imfns{qq} '_persp.png'] ;
        end
        close all
        fig = figure('Visible', 'Off',  'units', 'centimeters', ...
        'position', [0,0,xwidth,ywidth]) ;

        if ~exist(imfn, 'file') || overwrite
            colorV = colorVs{qq} ;
            if ~vertexBasedTiling

                % draw cells
                patch('Faces',ff,'Vertices',vv,...
                    'FaceVertexCData',colorV,'FaceColor','flat', ...
                    'Edgecolor', 'none');
                hold on;

                % draw bonds/contours
                for cid = 1:max(cellIDs)
                    plot3(c3d(cellIDs==cid, 1), c3d(cellIDs==cid, 2), ...
                        c3d(cellIDs==cid, 3), '-', 'color', [0.5, 0.5, 0.5]); 
                    hold on;
                end
                cb = colorbar ;

                if qq == 1                    
                    ylabel(cb, ['area [' QS.spaceUnits '$^2$]'],   'interpreter', 'latex')
                    colormap inferno
                    caxis([0, areaMaxClim])
                elseif qq == 2
                    ylabel(cb, 'cell aspect ratio, $a/b$',   'interpreter', 'latex')
                    colormap inferno
                    caxis([1., 3])
                elseif qq == 3
                    ylabel(cb, 'cell anisotropy, $\left(a/b-1 \right)\cos 2\theta$',   'interpreter', 'latex')
                    caxis([-1.5, 1.5])
                    % colormap(blueblackred)
                    colormap(brewermap(256, '*RdBu'))
                elseif qq == 4
                    ylabel(cb, 'cell eccentricity, $e$',   'interpreter', 'latex')
                    caxis([0, 1])
                    colormap inferno
                elseif qq == 5
                    ylabel(cb, 'cell anisotropy, $e \cos 2\theta$',   'interpreter', 'latex')
                    caxis([-1, 1])
                    % colormap(blueblackred)
                    colormap(brewermap(256, '*RdBu'))
                elseif qq == 6
                    ylabel(cb, 'cell orientation, $\cos(2\theta)$', ...
                        'interpreter', 'latex')
                    caxis([-1, 1])
                    % colormap(blueblackred)
                    colormap(brewermap(256, '*RdBu'))
                elseif qq == 7
                    ylabel(cb, 'cell anisotropy, $\left(1-a/b \right)\cos 2\theta$',   'interpreter', 'latex')
                    caxis([-1.5, 1.5])
                    % colormap(blueblackred)
                    colormap(brewermap(256, '*RdBu'))
                else
                    error('handle this field name here')
                end
                
                %% OCCLUDE THE FAR SIDE OF THE ORGAN
                if occlude
                    trisurf(triangulation(mesh4occlusion.f, mesh4occlusion.v), ...
                        'facecolor', 'k')
                end
                
                %% FORMAT VIEW
                axis equal
                if viewID == 0
                    view(0,0)
                else
                    % plot perspective view
                    view(-20, 20)
                end
                xlim(xyzlims(1, :))
                ylim(xyzlims(2, :))
                zlim(xyzlims(3, :))
                xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                title(titlestr, 'interpreter', 'latex')

                %% Match texturepatch export options
                % formatting and export
                set(fig, 'PaperUnits', 'centimeters');
                set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

                set(gcf, 'color', 'w')
                if occlude
                    set(gca, 'color', 'k')
                end
                drawnow ;
                set(gcf, 'position', [0, 0, xwidth, ywidth])
                
                disp(['Saving figure to: ' imfn])
                export_fig(imfn, '-nocrop', '-r200')
                Fim = imread(imfn) ;
                % resize output image
                if ~all(size(Fim) == [imH,   imW, 3])
                    Fim = imresize(Fim, [imH ,  imW]) ;
                    imwrite(Fim, imfn) ;
                end

                % Simpler version
                % Fim = getFigureFrameCData(imWpix, imHpix) ;
                %imwrite(Fim, imfn) ;
                
                clf
                
                %% Add patches on top of texturePatch image
                if occlude
                    

                    opts = metadat ;
                    opts.subdir = 'overlays'; 
                    opts.timePoints = tp ;
                    opts.plot_perspective = true ;
                    opts.plot_dorsal = false ;
                    opts.plot_ventral = false ;
                    opts.plot_left = true ;
                    opts.plot_right = false ;
                    opts.blackFigure = false ;
                    opts.makeColorbar = true ;
                    QS.plotSeriesOnSurfaceTexturePatch(opts, Options)

                    figoutdir = fullfile(QS.dir.texturePatchIm , 'overlays') ;
                    if viewID == 0
                        axisROIFn = fullfile(imdir, 'lateral_embeddingROI_for_overlay.mat') ;
                        figPerspDir = fullfile(figoutdir, 'lateral1') ;
                        timFn = fullfile(figPerspDir, sprintf('patch_lateral1_%06d.png', tp)) ;
                    
                    else
                        axisROIFn = fullfile(imdir, 'perspective_embeddingROI_for_overlay.mat') ;
                        figPerspDir = fullfile(figoutdir, 'perspective') ;
                        timFn = fullfile(figPerspDir, sprintf('patch_persp_%06d.png', tp)) ;
                    
                    end
                    tim = imread(timFn) ;
                    size(tim)

                    % Add images together -- area to add is in axisROI
                    % tim on outside, pim right of colorbarX
                    im2 = tim ;
                    if ~exist(axisROIFn, 'file')
                        axisROI = roipoly(0.5*Fim+0.5*tim) ;
                        colorbarX = 1010 ;
                        save(axisROIFn, 'axisROI', 'colorbarX')    
                    else
                        load(axisROIFn, 'axisROI', 'colorbarX')
                    end
                    % sum together inside axisROI
                    sim = tim + 0.9 * Fim ;
                    pR = squeeze(sim(:, :, 1)) ;
                    pG = squeeze(sim(:, :, 2)) ;
                    pB = squeeze(sim(:, :, 3)) ;
                    im2R = squeeze(im2(:,:,1)) ;
                    im2G = squeeze(im2(:,:,1)) ;
                    im2B = squeeze(im2(:,:,1)) ;
                    im2R(axisROI) = pR(axisROI) ;
                    im2G(axisROI) = pG(axisROI) ;
                    im2B(axisROI) = pB(axisROI) ;
                    im2 = cat(3, im2R, im2G, im2B) ;
                    % pim right of colorbarX
                    im2(:, colorbarX:end, :) = Fim(:, colorbarX:end, :) ;

                    disp(['Saving combined image to: ' imfn])
                    imwrite(im2, imfn)
                end
                
            else
                error('handle here')
            end
        else
            disp('image exists on file')
        end
    end
end

%% Draw bonds
imfn = fullfile(imdir, sprintf('cellseg3d_bonds_full_%06d.png', tp)) ;
if ~exist(imfn, 'file') || overwrite 
    close all

    fig = figure('Visible', 'Off',  'units', 'centimeters', ...
        'position', [0,0,xwidth,ywidth]) ;

    for fullID = [0, 1] 
        if vertexBasedTiling
            Xs = zeros(4*length(c3d(:, 1)), 1) ;
            Ys = Xs ;
            Zs = Xs ;
            Us = Xs ;
            Vs = Xs ;
            Ws = Xs ;
            dmyk = 1 ;
            for qq = 1:nVertices
                for id = seg3d.vdat.NL(qq, :)
                    if id > 0
                        Xs(dmyk) = c3d(qq, 1) ;
                        Ys(dmyk) = c3d(qq, 2) ; 
                        Zs(dmyk) = c3d(qq, 3) ;
                        Us(dmyk) = c3d(id, 1) - c3d(qq, 1) ;
                        Vs(dmyk) = c3d(id, 2) - c3d(qq, 2) ; 
                        Ws(dmyk) = c3d(id, 3) - c3d(qq, 3) ;
                        dmyk = dmyk + 1 ;
                    end
                end
            end
            plot3(c3d(:, 1), c3d(:, 2), c3d(:, 3), '.')
            hold on;
            q = quiver3(Xs,Ys,Zs, Us, Vs, Ws, 0, 'color', [ 0.8500    0.3250    0.0980]);
            axis equal
            q.ShowArrowHead = 'off';
            [~, ~, ~, xyzlims] = QS.getXYZLims() ;
            xlim(xyzlims(1, :))
            if fullID
                ylim(xyzlims(2, :))
                imfn = fullfile(imdir, sprintf('cellseg3d_bonds_full_%06d.png', tp)) ;
            else
                ylim([xyzlims(2, 1), 0])
                imfn = fullfile(imdir, sprintf('cellseg3d_bonds_%06d.png', tp)) ;    
            end
            zlim(xyzlims(3, :))
            view(0,0)
            xlabel(['ap position [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            ylabel(['lateral position [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            zlabel(['dv position [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            title(titlestr, 'interpreter', 'latex')

        else
            clf
            for cid = 1:max(cellIDs)
                plot3(c3d(cellIDs==cid, 1), c3d(cellIDs==cid, 2), ...
                    c3d(cellIDs==cid, 3), '-'); 
                hold on;
            end
            axis equal
            [~, ~, ~, xyzlims] = QS.getXYZLims() ;
            xlim(xyzlims(1, :))
            if fullID
                ylim(xyzlims(2, :))
                imfn = fullfile(imdir, sprintf('cellseg3d_bonds_full_%06d.png', tp)) ;
            else
                ylim([xyzlims(2, 1), 0])
                imfn = fullfile(imdir, sprintf('cellseg3d_bonds_%06d.png', tp)) ;    
            end
            zlim(xyzlims(3, :))
            view(0,0)
            xlabel(['ap position [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            ylabel(['lateral position [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            zlabel(['dv position [' QS.spaceUnits ']'], 'Interpreter', 'latex')
            title(titlestr, 'interpreter', 'latex')


        end

        % formatting and export
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'PaperPosition', [0 0 xwidth ywidth]);

        set(gcf, 'color', 'w')
        drawnow ;

        % Match texturepatch export options
        set(gcf, 'position', [0, 0, xwidth, ywidth])
        export_fig(imfn, '-nocrop', '-r200')
        Fim = imread(imfn) ;
        % resize output image
        if ~all(size(Fim) == [imH,   imW, 3])
            Fim = imresize(Fim, [imH ,  imW]) ;
            imwrite(Fim, imfn) ;
        end

        % Simpler version
        % Fim = getFigureFrameCData(imWpix, imHpix) ;
        clf
    end
end

%% Statistics
statsfn = fullfile(imdir, sprintf('stats_%06d.png', tp)) ;
if ~exist(statsfn, 'file') || overwrite
    clf
    % plot(sqrt(i11), areas(keep), '.') ; hold on;
    % plot(sqrt(i22), areas(keep), '.') ; hold on;
    subplot(2, 1, 1)
    plot(areas(keep), sqrt(seg3d.qualities.moment2(keep)), '.') ; hold on;
    plot(areas(keep), sqrt(seg3d.qualities.moment1(keep)), '.') ; 
    xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
    ylabel('$\sqrt{\lambda_2}, \sqrt{\lambda_1}$', 'interpreter', 'latex')
    subplot(2, 2, 3)
    ars_tmp = sqrt(seg3d.qualities.moment2(keep)./seg3d.qualities.moment1(keep)) ;
    plot(areas(keep), ars_tmp, '.') ; 
    xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
    ylabel('$\sqrt{\lambda_2 / \lambda_1}$', 'interpreter', 'latex')
    % Fit to line to see if there is variation
    [cc, SS] = polyfit(areas(keep), ars_tmp, 1) ;
    uncs = sqrt(abs(SS.R)) / SS.df ;
    title(['$\sqrt{\lambda_2 / \lambda_1} = ($' ...
        num2str(round(cc(1), 1, 'significant')) '$\pm$' ...
        num2str(round(uncs(1,1), 1, 'significant')) '$)A + $' ...
        num2str(round(cc(2), 3, 'significant')) '$\pm$' ...
        num2str(round(uncs(2,2), 1, 'significant'))], ...
        'interpreter', 'latex')

    subplot(2, 2, 4)
    plot(areas(keep), ars_tmp, '.') ;
    corrs = corrcoef(areas(keep), ars_tmp) ; 
    title(['correlation = ' ...
        num2str(round(corrs(1, 2), 2, 'significant'))], ...
        'interpreter', 'latex')

    ylim([0.5, 5])
    xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
    ylabel('$\sqrt{\lambda_2 / \lambda_1}$', 'interpreter', 'latex')
    sgtitle(titlestr, 'interpreter', 'latex')
    saveas(gcf, statsfn)
end
end




%% Original non-pathline code:
% 
% function plotCells(QS, tp, seg3d, imdir, overwrite, xyzlims, corrected)
% 
%     %% Draw cells colored by area
%     t0 = QS.t0() ;
%     titlestr = ['$t=$' sprintf('%03d', tp-t0) ' ' QS.timeUnits] ;
%     
%     % Easiest way is to triangulate all polygons using centroids
%     % This is fine if the cells are all convex
%     if ~corrected
%         faces = seg3d.cdat.polygons ;
%     else
%         faces = seg3d.cdat.polygonVertexID ;
%     end
%     keep = seg3d.statistics.keep ;
%     
%     areas = seg3d.qualities.areas ; 
%     ang1 = seg3d.qualities.ang1 ; 
%     ang2 = seg3d.qualities.ang2 ; 
%     mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
%     moinertia = seg3d.qualities.mInertia ;
%     c3d = seg3d.vdat.xyzrs ;
%     cellCntrd = seg3d.cdat.centroids_3d ;
%     keep = seg3d.statistics.keep ;
%     Qxx = seg3d.qualities.nematicStrength .* cos(2*ang1) ;
%     
%     nCells = length(faces) ;
%     nVertices = size(seg3d.vdat.uv, 1) ;
%     dmyk = 1 ;
%     ff = zeros(nCells * 7, 3) ;
%     areaV = NaN * zeros(nCells * 7, 1) ;
%     ang1V = areaV ;
%     mratioV = areaV ;
%     oparmV = areaV ;
%     IxxV = areaV ;
%     IxyV = areaV ;
%     IyyV = areaV ;
%     for cid = 1:nCells
%         if ismember(cid, keep)
%             face = faces{cid} ;
%             if numel(face) > 1
%                 for vid = 1:length(face)
%                     if vid < length(face)
%                         addface = [face(vid), face(vid+1), nVertices + cid] ;
%                     else
%                         addface = [face(vid), face(1), nVertices + cid] ;
%                     end
%                     ff(dmyk, :) = addface ;
%                     areaV(dmyk) = areas(cid) ;
%                     ang1V(dmyk) = ang1(cid) ;
%                     mratioV(dmyk) = mratio(cid) ;
%                     % qualities
%                     IxxV(dmyk) = moinertia(cid, 1) ;
%                     IxyV(dmyk) = moinertia(cid, 2) ;
%                     IyyV(dmyk) = moinertia(cid, 3) ;
%                     % order parameter
%                     oparmV(dmyk) = (mratio(cid) - 1) * cos(2*ang1(cid)) ;
%                     % nematic strength
%                     QxxV(dmyk) = Qxx(cid) ;
%                     
%                     dmyk = dmyk + 1 ;
%                 end
% 
%                 % check it 
%                 % clf
%                 % patch('Faces',ff(1:dmyk-1, :),'Vertices', [c3d; cellCntrd],...
%                 %     'FaceVertexCData',areaV(1:dmyk-1),'FaceColor','flat', ...
%                 %    'Edgecolor', 'k');
%                 % hold on;
%                 % patch('Faces',ff(dmyk-1:dmyk-1, :),'Vertices', [c3d; cellCntrd],...
%                 %     'FaceVertexCData',0,'FaceColor','flat', ...
%                 %    'Edgecolor', 'k');
%                 % scatter3(c3d(1:123, 1), c3d(1:123, 2), c3d(1:123, 3), 40, 1:123, 'filled')
%             end
%         end
%     end
%     
%     % [ff, vv] = poly2fv({x1, x2, x3}, {y1, y2, y3});
%     
%     ff = ff(1:dmyk-1, :) ;
%     areaV = areaV(1:dmyk-1) ;
%     ang1V = ang1V(1:dmyk-1) ;
%     IxxV = IxxV(1:dmyk-1) ;
%     IxyV = IxyV(1:dmyk-1) ;
%     IyyV = IyyV(1:dmyk-1) ;
%     mratioV = mratioV(1:dmyk-1) ;
%     oparmV = oparmV(1:dmyk-1) ;
%     QxxV = QxxV(1:dmyk-1) ;
%     
%     % Extend vertices to include centroids
%     vv = [c3d; cellCntrd] ;
%     
%     %% Color segmentation by Qxx
%     imfn = fullfile(imdir, sprintf('cellseg3d_Qxx_%06d.png', tp)) ;
%     if ~exist(imfn, 'file') || overwrite
%         clf
%         colorV = QxxV(:) ;
%         patch('Faces',ff,'Vertices',vv,...
%             'FaceVertexCData',colorV(:),'FaceColor','flat', ...
%             'Edgecolor', 'none');
%         cb = colorbar ;
%         ylabel(cb, ['cell elongation $Q_{xx}$'],   'interpreter', 'latex')
%         axis equal
%         view(0,0)
%         caxis([-1.5, 1.5])
%         colormap(blueblackred)
%         xlim(xyzlims(1, :))
%         ylim(xyzlims(2, :))
%         zlim(xyzlims(3, :))
%         xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         title(titlestr, 'interpreter', 'latex')
%         saveas(gcf, imfn)
%     end
%     
%     %% Color segmentation by area
%     imfn = fullfile(imdir, sprintf('cellseg3d_area_%06d.png', tp)) ;
%     if ~exist(imfn, 'file') || overwrite
%         clf
%         colorV = areaV(:) ;
%         patch('Faces',ff,'Vertices',vv,...
%             'FaceVertexCData',colorV(:),'FaceColor','flat', ...
%             'Edgecolor', 'none');
%         cb = colorbar ;
%         ylabel(cb, ['area [' QS.spaceUnits '$^2$]'],   'interpreter', 'latex')
%         caxis([0, min(max(areas), nanmean(areas) + 2.5*nanstd(areas))])
%         axis equal
%         view(0,0)
%         xlim(xyzlims(1, :))
%         ylim(xyzlims(2, :))
%         zlim(xyzlims(3, :))
%         colormap viridis
%         xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         title(titlestr, 'interpreter', 'latex')
%         saveas(gcf, imfn)
%     end
%     
%     %% Color segmentation by moi ratio -- sp coord sys and principal coords
%     imfns = {fullfile(imdir, sprintf('cellseg3d_mratioSP_log_%06d.png', tp)), ...
%         fullfile(imdir, sprintf('cellseg3d_mratioSP_%06d.png', tp))};
%     for tmp = 1:2
%                 
%         if ~exist(imfns{tmp}, 'file') || overwrite 
%             clf
%             bbr256 = blueblackred ;
%             if tmp == 1
%                 patch('Faces',ff,'Vertices',vv,...
%                     'FaceVertexCData',real(log10(sqrt(IyyV ./ IxxV))), ...
%                     'FaceColor','flat', ...
%                     'Edgecolor', 'none');
%                 cb = colorbar ;
%                 ylabel(cb, '$\log_{10} \sqrt{I_{\phi\phi}/I_{\zeta\zeta}}$',   'interpreter', 'latex')
%                 caxis([-1, 1])
%             else
%                 patch('Faces',ff,'Vertices',vv,...
%                     'FaceVertexCData',real(sqrt(IyyV ./ IxxV)), ...
%                     'FaceColor','flat', ...
%                     'Edgecolor', 'none');
%                 cb = colorbar ;
%                 ylabel(cb, '$\sqrt{I_{\phi\phi}/I_{\zeta\zeta}}$',   'interpreter', 'latex')
%                 caxis([0, 2])
%             end
%             colormap(bbr256)
%             axis equal
%             view(0,0)
%             xlim(xyzlims(1, :))
%             ylim(xyzlims(2, :))
%             zlim(xyzlims(3, :))
%             xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%             ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%             zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%             title(titlestr, 'interpreter', 'latex')
%             saveas(gcf, imfns{tmp})
%         end
%     end
%     
%     %% Order parameter    
%     imfn = fullfile(imdir, sprintf('cellseg3d_order_%06d.png', tp)) ;
%     if ~exist(imfn, 'file') || overwrite
%         patch('Faces',ff,'Vertices',vv,...
%             'FaceVertexCData',cos(2*ang1V(:)),'FaceColor','flat', ...
%             'Edgecolor', 'none');
%         cb = colorbar ;
%         ylabel(cb, '$\cos 2\theta$',   'interpreter', 'latex')
%         caxis([-1, 1])
%         colormap(blueblackred)
%         axis equal
%         view(0,0)
%         xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
%         title(titlestr, 'interpreter', 'latex')
%         saveas(gcf, imfn)
%     end
%     
%     
%     %% Draw bonds
%     imfn = fullfile(imdir, sprintf('cellseg3d_bonds_full_%06d.png', tp)) ;
%     if ~exist(imfn, 'file') || overwrite 
%         % Draw full and half of organ (up to midsaggital plane)
%         for fullID = [0, 1] 
%             if corrected
%                 clf
% 
%                 tmp = patch('Faces',ff,'Vertices',vv,...
%                     'FaceVertexCData',cos(2*ang1V(:)),'FaceColor','flat', ...
%                     'Edgecolor', 'none', 'FaceVertexAlphaData', 0.2);
%                 tmp.FaceAlpha = 'flat' ;
%                 hold on;
%                 for pid = 1:length(seg3d.cdat.polygons)
%                     pgon = seg3d.cdat.polygons{pid} ;
%                     plot3(pgon(:, 1), pgon(:, 2), pgon(:, 3), 'k-')
%                     hold on;
%                 end
%             else
%                 % Plot as network
%                 Xs = zeros(4*length(c3d(:, 1)), 1) ;
%                 Ys = Xs ;
%                 Zs = Xs ;
%                 Us = Xs ;
%                 Vs = Xs ;
%                 Ws = Xs ;
%                 dmyk = 1 ;
%                 for qq = 1:nVertices
%                     for id = seg3d.vdat.NL(qq, :)
%                         if id > 0
%                             Xs(dmyk) = c3d(qq, 1) ;
%                             Ys(dmyk) = c3d(qq, 2) ; 
%                             Zs(dmyk) = c3d(qq, 3) ;
%                             Us(dmyk) = c3d(id, 1) - c3d(qq, 1) ;
%                             Vs(dmyk) = c3d(id, 2) - c3d(qq, 2) ; 
%                             Ws(dmyk) = c3d(id, 3) - c3d(qq, 3) ;
%                             dmyk = dmyk + 1 ;
%                         end
%                     end
%                 end
%                 plot3(c3d(:, 1), c3d(:, 2), c3d(:, 3), '.')
%                 hold on;
%                 q = quiver3(Xs,Ys,Zs, Us, Vs, Ws, 0, 'color', [ 0.8500    0.3250    0.0980]);
%                 axis equal
%                 q.ShowArrowHead = 'off';
%             end
%             axis equal
%             [~, ~, ~, xyzlims] = QS.getXYZLims() ;
%             xlim(xyzlims(1, :))
%             if fullID
%                 ylim(xyzlims(2, :))
%                 imfn = fullfile(imdir, sprintf('cellseg3d_bonds_full_%06d.png', tp)) ;
%             else
%                 ylim([xyzlims(2, 1), 0])
%                 imfn = fullfile(imdir, sprintf('cellseg3d_bonds_%06d.png', tp)) ;    
%             end
%             zlim(xyzlims(3, :))
%             view(0,0)
%             xlabel(['ap position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
%             ylabel(['lateral position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
%             zlabel(['dv position, [' QS.spaceUnits ']'], 'Interpreter', 'latex')
%             title(titlestr, 'interpreter', 'latex')
%             saveas(gcf, imfn) 
%             clf
%         end
%     end
%     
%     %% Statistics
%     statsfn = fullfile(imdir, sprintf('stats_%06d.png', tp)) ;
%     if ~exist(statsfn, 'file') || overwrite || true
%         clf
%         % plot(sqrt(i11), areas(keep), '.') ; hold on;
%         % plot(sqrt(i22), areas(keep), '.') ; hold on;
%         subplot(2, 1, 1)
%         plot(areas(keep), sqrt(seg3d.qualities.moment2(keep)), '.') ; hold on;
%         plot(areas(keep), sqrt(seg3d.qualities.moment1(keep)), '.') ; 
%         xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
%         ylabel('$\sqrt{\lambda_2}, \sqrt{\lambda_1}$', 'interpreter', 'latex')
%         subplot(2, 2, 3)
%         ars_tmp = sqrt(seg3d.qualities.moment2(keep)./seg3d.qualities.moment1(keep)) ;
%         plot(areas(keep), ars_tmp, '.') ; 
%         xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
%         ylabel('$\sqrt{\lambda_2 / \lambda_1}$', 'interpreter', 'latex')
%         % Fit to line to see if there is variation
%         [cc, SS] = polyfit(areas(keep), ars_tmp, 1) ;
%         uncs = sqrt(abs(SS.R)) / SS.df ;
%         title(['$\sqrt{\lambda_2 / \lambda_1} = ($' ...
%             num2str(round(cc(1), 1, 'significant')) '$\pm$' ...
%             num2str(round(uncs(1,1), 1, 'significant')) '$)A + $' ...
%             num2str(round(cc(2), 3, 'significant')) '$\pm$' ...
%             num2str(round(uncs(2,2), 1, 'significant'))], ...
%             'interpreter', 'latex')
%         
%         subplot(2, 2, 4)
%         plot(areas(keep), ars_tmp, '.') ;
%         corrs = corrcoef(areas(keep), ars_tmp) ; 
%         title(['correlation = ' ...
%             num2str(round(corrs(1, 2), 2, 'significant'))], ...
%             'interpreter', 'latex')
%         
%         ylim([0.5, 5])
%         xlabel(['area [' QS.spaceUnits '$^2$]'], 'interpreter', 'latex')
%         ylabel('$\sqrt{\lambda_2 / \lambda_1}$', 'interpreter', 'latex')
%         sgtitle(titlestr, 'interpreter', 'latex')
%         saveas(gcf, statsfn)
%     end
% end
