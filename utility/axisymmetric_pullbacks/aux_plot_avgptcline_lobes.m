function aux_plot_avgptcline_lobes(folds, fold_onset, lobeDir, dvexten, ...
    save_ims, overwrite_lobeims, timePoints, spcutMeshBase, clineDVhoopBase, ...
    t0, timeInterval, timeUnits, spaceUnits, tidxMap)
%AUX_PLOT_AVGPTCLINE_LOBES auxiliary function for plotting the motion of
%the constrictions between lobes and the centerlines over time
% 
% Parameters
% ----------
% tp : N x 1 int array
%   xp.fileMeta.timePoints 
% timePoints : N x 1 int array
%   the timepoints in the experiment (xp.fileMeta.timePoints)
% t0, min(fold_onset)
% 
% Returns
% -------
% <none>
%
% NPMitchell 2020 


tp = (timePoints - t0) * timeInterval ;
fold_dynamics_figfn = fullfile(lobeDir, ['constriction_dynamics' dvexten '.png']) ;
if save_ims && (~exist(fold_dynamics_figfn, 'file') || overwrite_lobeims)
    disp('Creating constriction dynamics plots...')
    f1pts = zeros(length(tp), 3) ;
    f2pts = zeros(length(tp), 3) ;
    f3pts = zeros(length(tp), 3) ;

    disp('Loading centerline position of each fold...')
    for kk = 1:length(timePoints)
        % Translate to which timestamp
        t = timePoints(kk) ;
        load(sprintf(spcutMeshBase, t), 'spcutMesh') ;

        % Load the centerline too
        avgpts = spcutMesh.avgpts ;

        % rename the fold indices (in U)
        f1 = folds(kk, 1) ;
        f2 = folds(kk, 2) ;
        f3 = folds(kk, 3) ;

        % store distance from x axis of folds
        f1pts(kk, :) = avgpts(f1, :) ;
        f2pts(kk, :) = avgpts(f2, :) ;
        f3pts(kk, :) = avgpts(f3, :) ;
    end
    close all
    alph = 0.2 ;
    cmap = colormap ;
    fig = figure('visible', 'off'); 
    hold on;
    for qq = 1:length(tp)
        t = timePoints(qq) ;
        % Load the centerline
        fn = sprintf(clineDVhoopBase, t) ;
        load(fn, 'avgpts')

        color = cat(2, cmap(uint8(max(1, qq * length(cmap)/length(tp))), :), alph);
        % plot([f1pts(qq, 2), f2pts(qq, 2), f3pts(qq, 2)], ...
        %     [f1pts(qq, 3), f2pts(qq, 3), f3pts(qq, 3)], '-', 'color', color)
        plot(avgpts(:, 2), avgpts(:, 3), '-', 'color', color)
    end
    sz = 40 ;
    msz = 20 ;
    scatter(f1pts(:, 2), f1pts(:, 3), sz, tp, 'o'); hold on;
    scatter(f2pts(:, 2), f2pts(:, 3), sz, tp, 's');
    scatter(f3pts(:, 2), f3pts(:, 3), sz, tp, '^');
    fons1 = max(fold_onset(1), 1) ;
    fons2 = max(fold_onset(2), 1) ;
    fons3 = max(fold_onset(3), 1) ;
    plot(f1pts(tidxMap(fons1), 2), f1pts(tidxMap(fons1), 3), 'ko', 'markersize', msz); 
    plot(f2pts(tidxMap(fons2), 2), f2pts(tidxMap(fons2), 3), 'ks', 'markersize', msz);
    plot(f3pts(tidxMap(fons3), 2), f3pts(tidxMap(fons3), 3), 'k^', 'markersize', msz);
    axis equal
    xlabel(['y [' spaceUnits ']'], 'interpreter', 'latex')
    ylabel(['z [' spaceUnits ']'], 'interpreter', 'latex')
    title('Constriction dynamics')
    cb = colorbar() ;
    cb.Label.String = ['time [' timeUnits ']'] ;
    cb.Label.Interpreter = 'Latex' ;
    disp(['Saving figure to ' fold_dynamics_figfn])
    saveas(fig, fold_dynamics_figfn) ;
    close all
else
    disp(['Skipping constriction dynamics plots since they exist (' fold_dynamics_figfn ')...'])
end


