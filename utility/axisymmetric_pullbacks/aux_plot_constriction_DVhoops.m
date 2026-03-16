function aux_plot_constriction_DVhoops(folds, fold_onset, outDir, dvexten, ...
    save_ims, overwrite_lobeims, timePoints, spcutMeshBase, ...
    alignedMeshBase, normal_shift, rot, trans, resolution, colors, ...
    xyzlim, flipy, t0, timeInterval, timeUnits, spaceUnits, tidxMap)
%AUX_PLOT_AVGPTCLINE_LOBES auxiliary function for plotting the motion of
%the constrictions between lobes and the centerlines over time
% 
% Parameters
% ----------
% timePoints : N x 1 int array
%   the timepoints in the experiment (xp.fileMeta.timePoints)
% tidxMap : #timePoints in xp x 1 array 
%   mapping from timePoint values to indices
%
% Returns
% -------
%
% NPMitchell 2020 


tp = (timePoints - t0) * timeInterval ;
c1 = colors(1, :) ;
c2 = colors(2, :) ;
c3 = colors(3, :) ;
shift_mag = (normal_shift * resolution) + 1 ; 
xlims = xyzlim(1, :) + [-shift_mag, shift_mag] ;
ylims = xyzlim(2, :) + [-shift_mag, shift_mag] ;
zlims = xyzlim(3, :) + [-shift_mag, shift_mag] ;

fold_ant_figfn = fullfile(outDir, 'constriction_anterior_%06d') ;
fold_lat_figfn = fullfile(outDir, 'constriction_lat_%06d') ;
fns = dir(strrep(fold_ant_figfn, '%06d', '*')) ;
if save_ims && (~(length(fns)==length(timePoints)) || overwrite_lobeims)
    disp('Creating constriction dynamics plots...')
    
    disp('Loading hoops for each fold...')
    for kk = 1:length(timePoints)
        % Translate to which timestamp
        t = timePoints(kk) ;
        tp4title = tp(kk) ;
        load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
        mesh = read_ply_mod(sprintf(alignedMeshBase, t)) ;
        % shift to match spcutMesh > note assume inward normal, but shift 
        % is outward, so subtract.
        mesh.v = mesh.v - normal_shift * resolution * mesh.vn ;
        nU = spcutMesh.nU ;
        nV = spcutMesh.nV ;

        % Load the centerline too
        avgpts = spcutMesh.ringpath_ss ;

        % rename the fold indices (in U)
        tidx = find(tidxMap == t) ;
        f1 = folds(tidx, 1) ;
        f2 = folds(tidx, 2) ;
        f3 = folds(tidx, 3) ;

        % store distance from x axis of folds
        % f1pt = avgpts(f1, :) ;
        % f2pt = avgpts(f2, :) ;
        % f3pt = avgpts(f3, :) ;
        
        % rotate vertices 
        vrs = ((rot * spcutMesh.v')' + trans) * resolution ;
        if flipy 
            vrs(:, 2) = - vrs(:, 2) ; 
        end
        
        % plot DVhoop
        f1 = double(f1) ;
        f2 = double(f2) ;
        f3 = double(f3) ;
        hoop1 = vrs(f1+nU*(0:(nV-1)), :) ;
        hoop2 = vrs(f2+nU*(0:(nV-1)), :) ;
        hoop3 = vrs(f3+nU*(0:(nV-1)), :) ;

        % store distance from x axis of folds
        f1pt = mean(hoop1(1:end-1, :)) ;
        f2pt = mean(hoop2(1:end-1, :)) ;
        f3pt = mean(hoop3(1:end-1, :)) ;
        
        % Plot it
        close all
        alph = 0.04 ;
        fig = figure('visible', 'off', 'position', [0,0,800,600]); 
        hold on;
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            mesh.v(:, 1), 'EdgeColor', 'none', 'FaceAlpha', alph)

        sz = 20 ;
        msz = 150 ;
        lw = 4 ;
        plot3(hoop1(:, 1), hoop1(:, 2), hoop1(:, 3), '-', 'color', c1, 'linewidth', lw); hold on;
        plot3(hoop2(:, 1), hoop2(:, 2), hoop2(:, 3), '-', 'color', c2, 'linewidth', lw);
        plot3(hoop3(:, 1), hoop3(:, 2), hoop3(:, 3), '-', 'color', c3, 'linewidth', lw);

        
        scatter3(f1pt(1), f1pt(2), f1pt(3), msz, ...
           'o', 'filled', 'markerfacecolor', c1); 
        scatter3(f2pt(1), f2pt(2), f2pt(3), msz, ...
           's', 'filled', 'markerfacecolor', c2);
        scatter3(f3pt(1), f3pt(2), f3pt(3), msz, ...
           '^', 'filled', 'markerfacecolor', c3);

        
        % if t < fold_onset(1)    
        %     plot3(f1pt(1), f1pt(2), f1pt(3), ...
        %         '.', 'color', c1, 'markersize', sz);
        % else
        %     scatter3(f1pt(1), f1pt(2), f1pt(3), msz, ...
        %         'o', 'filled', 'markerfacecolor', c1); 
        % end
        % if t < fold_onset(2)
        %     plot3(f2pt(1), f2pt(2), f2pt(3), ...
        %         '.', 'color', c2, 'markersize', sz);
        % else
        %     scatter3(f2pt(1), f2pt(2), f2pt(3), msz, ...
        %         's', 'filled', 'markerfacecolor', c2);
        % end
        % if t < fold_onset(3)
        %     plot3(f3pt(1), f3pt(2), f3pt(3),...
        %         '.', 'color', c3, 'markersize', sz);
        % else
        %     scatter3(f3pt(1), f3pt(2), f3pt(3), msz, ...
        %         '^', 'filled', 'markerfacecolor', c3);
        % end
        axis equal
        xlim(xlims)
        ylim(ylims)
        zlim(zlims)
        xlabel(['AP position [' spaceUnits ']'], 'interpreter', 'latex')
        ylabel(['lateral position [' spaceUnits ']'], 'interpreter', 'latex')
        zlabel(['DV position [' spaceUnits ']'], 'interpreter', 'latex')
        title(sprintf(['Constriction dynamics, t=%03d ' timeUnits], tp4title), ...
            'interpreter', 'latex')
        ofn = sprintf(fold_lat_figfn, t) ;
        disp(['Saving figure to ' ofn])
        view(0, 0) ;
        saveas(fig, [ofn '.png']) ;
        
        
        ylim(xlims - mean(xlims))
        ofn = sprintf(fold_ant_figfn, t) ;
        disp(['Saving figure to ' ofn])
        view(-90, 0) 
        saveas(fig, [ofn '.png']) ;
        
        % save axis limits
        if kk == 1
            clf
            axis equal
            xlim(xlims)
            ylim(ylims)
            zlim(zlims)
            view(0, 0) ;
            ofn = sprintf(fold_lat_figfn, t) ;
            saveas(gcf, [ofn '_axis.pdf']) ; % '-r300') ;
            view(-90, 0) ;
            ofn = sprintf(fold_ant_figfn, t) ;
            saveas(gcf, [ofn '_axis.pdf']) ; % '-r300') ;
            close all
        end
    end
else
    disp(['Skipping constriction hoop dynamics plots since they exist (' fold_ant_figfn ')...'])
end


