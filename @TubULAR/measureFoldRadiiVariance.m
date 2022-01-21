function measureFoldRadiiVariance(QS, options)
%foldRadiiVariance(QS, OPTIONS)
%   Identify azimuthal variance in radius around AP feature locations/folds
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   Current QuapSlap object
% options : struct with fields
%   overwrite : bool, default=false
%       overwrite results on disk
%   preivew : bool, default=QS.plotting.preview
%       display intermediate results
%
% Returns
% -------
% Saves radii and vairances to disk. 
%   Data is saved here:
%       fullfile(lobeDir, ['fold_locations_sphi' dvexten '_avgpts.mat'])
%
% NPMitchell 2020

%% Unpack QS
lobeDir = QS.dir.lobe ;
nU = QS.nU ;
nV = QS.nV ;
uvexten = QS.uvexten ;  % string like sprintf('_nU%04d_nV%04d', nU, nV) ;
timePoints = QS.xp.fileMeta.timePoints ;
preview = QS.plotting.preview ;
spcutMeshBase = QS.fullFileBase.spcutMesh ;

%% Unpack options
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
% overwrite QS preview option if present in local options
if isfield(options, 'preview')
    preview = options.preview ;
end

%% First compute using the avgpts (DVhoop means)
disp('Identifying lobes...')

foldfn = fullfile(lobeDir, ['fold_locations_sphi' uvexten '_avgpts.mat']) ;
if exist(foldfn, 'file') && ~overwrite
    disp('Loading lobes')
    % Save the fold locations as a mat file
    load(foldfn, 'ssfold', 'folds', 'ssfold_frac', 'ssmax', 'fold_onset', ...
        'rssfold', 'rssfold_frac', 'rssmax', 'rmax')    
    
    error('here')
    
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
        avgpts = spcutMesh.avgpts ;

        % rename the fold indices (in U)
        f1 = folds(kk, 1) ;
        f2 = folds(kk, 2) ;
        f3 = folds(kk, 3) ;

        % store distance from x axis of folds
        f1pt = avgpts(f1, :) ;
        f2pt = avgpts(f2, :) ;
        f3pt = avgpts(f3, :) ;

        % rotate vertices 
        vrs = QS.xyz2APDV(spcutMesh.v) ;
        
        % plot DVhoop
        hoop1 = vrs(f1+nU*(0:(nV-1)), :) ;
        hoop2 = vrs(f2+nU*(0:(nV-1)), :) ;
        hoop3 = vrs(f3+nU*(0:(nV-1)), :) ;
        
        % Plot it
        close all
        alph = 0.1 ;
        fig = figure('visible', 'off'); 
        hold on;
        trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
            mesh.v(:, 1), 'EdgeColor', 'none', 'FaceAlpha', alph)

        sz = 20 ;
        msz = 50 ;
        plot3(hoop1(:, 1), hoop1(:, 2), hoop1(:, 3), '-', 'color', c1); hold on;
        plot3(hoop2(:, 1), hoop2(:, 2), hoop2(:, 3), '-', 'color', c2);
        plot3(hoop3(:, 1), hoop3(:, 2), hoop3(:, 3), '-', 'color', c3);
        if t < fold_onset(1)    
            plot3(f1pt(1), f1pt(2), f1pt(3), ...
                '.', 'color', c1, 'markersize', sz);
        else
            scatter3(f1pt(1), f1pt(2), f1pt(3), msz, ...
                'o', 'filled', 'markerfacecolor', c1); 
        end
        if t < fold_onset(2)
            plot3(f2pt(1), f2pt(2), f2pt(3), ...
                '.', 'color', c2, 'markersize', sz);
        else
            scatter3(f2pt(1), f2pt(2), f2pt(3), msz, ...
                's', 'filled', 'markerfacecolor', c2);
        end
        if t < fold_onset(3)
            plot3(f3pt(1), f3pt(2), f3pt(3),...
                '.', 'color', c3, 'markersize', sz);
        else
            scatter3(f3pt(1), f3pt(2), f3pt(3), msz, ...
                '^', 'filled', 'markerfacecolor', c3);
        end
        axis equal
        zlim(xlims)
        ylim(ylims)
        zlim(zlims)
        xlabel('AP position [\mum]')
        ylabel('lateral position [\mum]')
        zlabel('DV position [\mum]')
        title(sprintf('Constriction dynamics, t=%03d min', tp4title))
        ofn = sprintf(fold_lat_figfn, t) ;
        disp(['Saving figure to ' ofn])
        view(0, 0) ;
        saveas(fig, ofn) ;
        ofn = sprintf(fold_ant_figfn, t) ;
        disp(['Saving figure to ' ofn])
        view(-90, 0) 
        saveas(fig, ofn) ;
        close all
    end
end


