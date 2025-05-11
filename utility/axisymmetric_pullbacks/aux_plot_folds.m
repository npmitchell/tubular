function aux_plot_folds(folds, ssfold, ssfold_frac, ssmax, rmax, nU, ...
    timePoints, lobeDir, dvexten, sphiBase, method, overwrite, ...
    timeInterval, t0) 
% aux_plot_folds(folds, ssfold, ssfold_frac, ssmax, nU, timePoints, lobeDir, dvexten)
%
% Parameters
% ----------
% folds :
% ssfold : N x 1 float array
% ssfold_frac : N x 1 float array
%   Values between 0 and 1 of x/L
% ssmax : N x 1 float
%   Maximum pathlength for this timepoint in pixels
% rmax : N x 1 float
%   Maximum radius for this timepoint in pixels
% nU : int or double
%   Number of sampled points along the U axis of the mesh 
% timePoints : N x 1 int array
% resolution : float
%   Conversion from mesh coordinates to microns
% lobeDir : str
%   path to the directory where plots are saved
% dvexten : str
%   extension to specify APDV sampling, of the form '_nU%03d_nV%03d'
% method : 'avgpts' or 'ringpath'
%   Whether distances (pathlengths) were determined via DVhoop means or by
%   the ringpath (mean distances along the surface between hoops)
% 
% 
% Returns
% -------
%
% NPMitchell 2019

foldImDir = fullfile(lobeDir, ['images_foldID_' method]) ;
if ~exist(foldImDir, 'dir')
    mkdir(foldImDir)
end
mexten = ['_' method] ;

% Plot kymograph
close all; figure('visible', 'off'); hold on
if length(timePoints) > 1
    plot(folds(:, 1), (timePoints - t0) * timeInterval)
    plot(folds(:, 2), (timePoints - t0) * timeInterval)
    plot(folds(:, 3), (timePoints - t0) * timeInterval)
    ylim([(timePoints(1)- t0) * timeInterval, (timePoints(end)- t0) * timeInterval])
else
    plot(folds(:, 1), (timePoints - t0) * timeInterval, 'o')
    plot(folds(:, 2), (timePoints - t0) * timeInterval, 's')
    plot(folds(:, 3), (timePoints - t0) * timeInterval, '^')
end
xlim([0, nU])
ylabel('time [min]')
xlabel('position index')
title('Fold positions')
set(gca,'YDir','reverse')
saveas(gcf, fullfile(lobeDir, ['fold_idx' dvexten mexten '.png']))

% Plot the locations as a fraction of the total length as a kymograph
close all; figure('visible', 'off'); hold on
if length(timePoints) > 1
    plot(ssfold_frac(:, 1), (timePoints - t0) * timeInterval)
    plot(ssfold_frac(:, 2), (timePoints - t0) * timeInterval)
    plot(ssfold_frac(:, 3), (timePoints - t0) * timeInterval)
    ylim([(timePoints(1)- t0) * timeInterval, (timePoints(end)- t0) * timeInterval])
else
    plot(ssfold_frac(:, 1), (timePoints - t0) * timeInterval, 'o')
    plot(ssfold_frac(:, 2), (timePoints - t0) * timeInterval, 's')
    plot(ssfold_frac(:, 3), (timePoints - t0) * timeInterval, '^')
    
end
xlim([0, 1])

ylabel('time [min]')
xlabel('position [x/L]')
title('Fold positions')
set(gca,'YDir','reverse')
saveas(gcf, fullfile(lobeDir, ['fold_ssfrac' dvexten mexten '.png']))

% Plot the locations (pathlength) as a kymograph
close all; figure('visible', 'off'); hold on
if length(timePoints) > 1
    plot(ssfold(:, 1), (timePoints - t0) * timeInterval)
    plot(ssfold(:, 2), (timePoints - t0) * timeInterval)
    plot(ssfold(:, 3), (timePoints - t0) * timeInterval)
    plot(ssmax, (timePoints - t0) * timeInterval, 'k--')
    ylim([(timePoints(1)- t0) * timeInterval, (timePoints(end)- t0) * timeInterval])
else
    plot(ssfold(:, 1), (timePoints - t0) * timeInterval, 'o')
    plot(ssfold(:, 2), (timePoints - t0) * timeInterval, 's')
    plot(ssfold(:, 3), (timePoints - t0) * timeInterval, '^')
    plot(ssmax, (timePoints - t0) * timeInterval, 'ks')
end
xlim([0, max(ssmax)])
ylabel('time [min]')
xlabel('position [\mum]')
title('Fold positions')
set(gca,'YDir','reverse')
saveas(gcf, fullfile(lobeDir, ['fold_ss' dvexten mexten '.png']))
close all

% Plot location of each fold superimposed on the mean radius for uniformly
% sampled DVhoops
blue = [0 0.4470 0.7410] ;
red = [0.8500    0.3250    0.0980] ;
maroon = [0.6350    0.0780    0.1840]; 
yellow = [0.9290    0.6940    0.1250 ]; 
maxrmax = max(rmax) ;
maxss = max(ssmax) ;
fig = figure('visible', 'off') ;
for kk = 1:length(timePoints)
    % Translate to which timestamp
    t = timePoints(kk) ;
    timestr = sprintf('_%04d', t) ;

    ofn = fullfile(foldImDir, ['radii_folds' dvexten mexten timestr '.png']) ;
    if ~exist(ofn, 'file') || overwrite        
        load(sprintf(sphiBase, t), 'spcutMesh') ;

        % Find the average radius over the uniformly sampled hoops. 
        % (make radius 1d)
        rad = mean(spcutMesh.radii_from_mean_uniform_rs, 2) ;
        minrad = min(spcutMesh.radii_from_mean_uniform_rs, [], 2) ;
        maxrad = max(spcutMesh.radii_from_mean_uniform_rs, [], 2) ;
        
        stdradlo = rad - std(spcutMesh.radii_from_mean_uniform_rs, [], 2) ;
        stdradhi = rad + std(spcutMesh.radii_from_mean_uniform_rs, [], 2) ;

        if strcmp(method, 'avgpts')
            ss = spcutMesh.avgpts_ss ;
        elseif strcmp(method, 'ringpath')
            ss = spcutMesh.ringpath_ss ;
        end

        % Save plot of radius with fold locations marked
        fill([ss; flipud(ss)], [minrad; flipud(maxrad)], blue, ...
            'facealpha', 0.3, 'edgecolor', 'none')
        hold on;
        fill([ss; flipud(ss)], [stdradlo; flipud(stdradhi)], maroon, ...
            'facealpha', 0.3, 'edgecolor', 'none')
        plot(ss, rad, 'Color', 'k')
        plot(ss(folds(kk, :)), rad(folds(kk, :)), 'o', 'Color', yellow)
        xlim([0, maxss])
        ylim([0, maxrmax])
        % axis equal
        xlabel('AP position [\mum]')
        ylabel('radius [\mum]')
        title('Fold locations')
        disp(['Saving radius vs ss plot for t=' num2str(t) ' to: ' ofn])
        saveas(fig, ofn)
        clf
    end
end

