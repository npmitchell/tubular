function plotClineXSections(QS, options)
% Plot fancy "cross-section" view of centerlines, saved to DVhoop
% centerline directory
%
% NPMitchell 2021

% Default options
overwrite = false ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end

% Define output dir
fancyClineDir = fullfile(QS.dir.clineDVhoop, 'fancyCenterlineFigures') ;
if ~exist(fancyClineDir, 'dir')
    mkdir(fancyClineDir)
end
% Define output filename base
figfn = fullfile(fancyClineDir, 'xsection_cline_dvhoop_%06d.png') ;
figfn2 = fullfile(fancyClineDir, 'outline_cline_dvhoop_%06d.png') ;
[~,~,~,xyzlim_um] = QS.getXYZLims() ;
xwidth = 12 ;
ywidth = 16 ;

colors = define_colors() ;
facecolor = colors(5, :) ; % green
linecolor = colors(9, :) ; % brown 
surfcolor = colors(4, :) ; % purple

timePoints =  QS.xp.fileMeta.timePoints ;
tidx2do = 1:50:length(timePoints) ;
tidx2do = [tidx2do, setdiff(1:10:length(timePoints), tidx2do) ] ; 
tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do) ] ; 

tp2do = timePoints(tidx2do) ;
% tp2do = [83, 123, 183, 253:263] ;
% tp2do = [tp2do, setdiff(timePoints(tidx2do), tp2do)] ;

for tt = tp2do
    QS.setTime(tt) ;

    fig_saved = ~isempty(dir(sprintf(figfn, tt))) ;
    if (~fig_saved || overwrite)
        % %% Plot and save
        % if resolution_matches
        %     disp(['Loading skelrs from disk: ' skeloutfn])
        %     skeloutfn = [skel_rs_outfn '.txt'] ;
        % else
        %     other_skelrs = dir([fullfile(outdir, name) '_centerline_scaled*']) ;
        %     skeloutfn = fullfile(other_skelrs(1).folder, other_skelrs(1).name) ;
        % end
        % ssskelrs = importdata(skeloutfn) ;
        % skelrs = ssskelrs(:, 2:4) ;

        [mss, mcline, avgpts] = QS.getCurrentClineDVhoop() ;
        
        % Save plot of rotated and translated mesh
        disp('Saving rotated & translated figure (xy)...')    
        close all
        % Load rotated mesh
        % if useSavedAPDVMeshes
        meshfn = sprintf(QS.fullFileBase.alignedMesh, tt) ;
        disp(['Loading ' meshfn])
        mesh = read_ply_mod(meshfn) ;
        xyzrs = mesh.v ;
        % If we are plotting the mesh, reverse the triangle ordering
        % for ambient occlusion to work properly, but the APDV mesh
        % should ALREADY have the correct y values (mirrored XZ)
        if QS.flipy
            tri = mesh.f(:, [1 3 2]) ;
        else
            tri = mesh.f ;
        end
        
        % Plot the result
        fig = figure('Visible', 'Off') ;
        tmp = trisurf(tri, xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), ...
            'edgecolor', 'none', 'facecolor', facecolor, ...
            'FaceAlpha', 0.4) ;
        
        hold on;
        % plot the skeleton
        plot3(avgpts(:,1), avgpts(:,2), avgpts(:,3), ...
            '-','Color', linecolor, 'LineWidth', 2);
        
        % Find spots where y component normal is nearly zero
        % mm = QS.getCurrentSPCutMeshSmRS()  ;
        % Get spine with normal nearly in plane
        % vn = reshape(mm.vn, [mm.nU, mm.nV, 3]) ;
        % vv = reshape(mm.v, [mm.nU, mm.nV, 3]) ;
        % 
        % for ring = 1:mm.nU 
        %     [~, idx] = find(min(abs(vn(ring, :, 2)))) ;
        %     if ring == 1
        %         spine = squeeze(vv(ring, idx, :))' ;
        %     else
        %         spine = [spine; squeeze(vv(ring, idx, :))'] ;
        %     end
        % end
        % plot3(spine(:, 1), spine(:, 2), spine(:, 3), 'color', linecolor) 
        
        % annotate figure with APDV
        xlabel('x [$\mu$m]', 'Interpreter', 'Latex');
        ylabel('y [$\mu$m]', 'Interpreter', 'Latex');
        zlabel('z [$\mu$m]', 'Interpreter', 'Latex');
        axis equal
        grid off
        % yz
        disp(['Saving rotated & translated figure (yz): ' sprintf(figfn, tt)])    
        view(0, 0);
        xlim(xyzlim_um(1, :)); 
        ylim(xyzlim_um(2, :)); 
        zlim(xyzlim_um(3, :)) ;
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperPosition', [0 0 xwidth ywidth]); %x_width=10cm y_width=15cm
        export_fig(sprintf(figfn, tt), '-nocrop', '-transparent', '-r200')
        
        % Load and resave with border
        clf
        ww = 3 ;
        im = imread(sprintf(figfn, tt)) ;
        im2 = 0*squeeze(im(:, :, 1)) ;
        im2(im(:,:,1) < 255) = 255 ;
        biggest = bwareafilt(logical(im2), 1) ;
        im2(biggest) = 0 ;
        se = offsetstrel('ball',ww,ww);
        im3 = imdilate(im2, se);
        im2 = imerode(im2, se) ;
        im2(im3(:) < 100) = 0 ;
        im2(im3(:) > 255) = 1 ;
        im3(logical(im2)) = 0 ;
        im3(im3(:) < 100) = 0 ;
        im3(im3(:) > 255) = 255 ;
        im3 = logical(im3) ;
        redC = im(:, :, 1) ; 
        bluC = im(:, :, 2) ; 
        grnC = im(:, :, 3) ; 
        redC(im3) = round(surfcolor(1) * 255) ;
        bluC(im3) = round(surfcolor(2) * 255) ;
        grnC(im3) = round(surfcolor(3) * 255) ;
        outim = (cat(3, redC, bluC, grnC)) ;
        imwrite(outim, sprintf(figfn2, tt)) ;
        
        % imshow(im)
        % hold on;
        % [C,h] = imcontour(im2, 1) ;
        % h.LineColor = surfcolor ;
        % h.LineWidth = 3 ;
        
    end
end
