function plotRelativeMotionTracks(QS, options)
% plotRelativeMotionTracks(QS, options)
%   Load tracked nuclei/objects in two layers of the evolving surface and
%   plot their pathlines in 3D on the surface.
%
% Parameters
% ----------
% QS : QuapSlap object
% options : struct with fields
%   layerLabel : str
%   relMotionFn : str
%       default is fullfile(QS.dir.tracking, 'relative_motion_tracks.mat')
%   lwTrail : float (default=2)
%   uniqueCorrespondence : bool (default = false)
%       plot only one endoderm cell for each muscle cell AND only one
%       muscle cell for each endoderm cell
%   jitterCorrespondence : float >=0 (default = 0) 
%       If more than one match in object {A} list for a given object B in
%       {B}, jitter matching object {A} by this amount in 3d
%   specifyTracks :
%
%
% Saves to disk
% -------------
% pngs: fullfile(QS.dir.tracking, 'relativeMotionImages',
%                'tracks_%06d.png')
%
%
% NPMitchell 2021

overlayTexturePatch = true ;
subdir_layer1 = 'overlaysEndoderm' ;
subdir_layer2 = 'overlaysMuscle' ;
styles2do = 1:4 ;
layer1_adjustDataLim = [8, 4191] ;
layer2_adjustDataLim = [8, 19241] ;
textureFactor = 0.7 ;


%% Prepare plot settings
if isfield(options, 'layerLabel')
    subdir = options.layerLabel ;
    exten = [ '_' subdir ] ;
else
    exten = '' ;
end
if isfield(options, 'overlayTexturePatch')
    overlayTexturePatch = options.overlayTexturePatch ;
end
if isfield(options, 'subdir_layer1')
    subdir_layer1 = options.subdir_layer1 ;
end
if isfield(options, 'subdir_layer2')
    subdir_layer2 = options.subdir_layer2 ;
end
if isfield(options, 'styles2do')
    styles2do = options.styles2do ;
end
opts = load(fullfile(QS.dir.texturePatchIm, ['metadat' exten '.mat'])) ;
xyzlims = opts.metadat.xyzlim ;
smoothing_lambda = opts.metadat.smoothing_lambda ;
normal_shift = opts.metadat.normal_shift ;
backdrop_normal_shift = opts.metadat.normal_shift ;
% Rotation = opts.options.Rotation ;
% Translation = opts.options.Translation ;
% Dilation = opts.options.Dilation ;
flipy = QS.flipy ;
meshFileBase = QS.fullFileBase.mesh ;
initialDistanceThres = Inf ;
overwrite = false ;
lwTrail= 1 ;  % linewidth for the trailing track in 3d
uniqueCorrespondence = false ;
jitterCorrespondence = 0 ;

%% Unpack options
relMotionFn = fullfile(QS.dir.tracking, 'relative_motion_tracks.mat') ;
timePoints = QS.xp.fileMeta.timePoints ;
if isfield(options, 'relMotionFn')  
    relMotionFn = options.relMotionFn ;
end
if isfield(options, 'specifyTracks')
    tracks2plot = options.specifyTracks ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'normal_shift')
    normal_shift = options.normal_shift ;
elseif isfield(options, 'normalShift')
    normal_shift = options.normalShift ;
end
if isfield(options, 'initialDistanceThres')
    initialDistanceThres = options.initialDistanceThres ;
end
if isfield(options, 'overwrite')  
    overwrite = options.overwrite ;
end
if isfield(options, 'lwTrail')  
    lwTrail = options.lwTrail ;
end
if isfield(options, 'jitterCorrespondence')
    jitterCorrespondence = options.jitterCorrespondence ;
end
if isfield(options, 'uniqueCorrespondence')
    uniqueCorrespondence = options.uniqueCorrespondence ;
end

% Load relative tracks
load(relMotionFn, 'dusEuclidean', 'dusGeodesic', ...
        'tracks1', 'tracks2', 'pairIDs', 'U0s', 'V0s', ...
        'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
        'euclideanDistanceTraveledU', 'euclideanDistanceTraveledV', ...
        'v3d_u', 'v3d_v', 'nSaved') 
        
nTracks = size(v3d_u, 1) ;
assert(nSaved == nTracks) ;
if ~exist('tracks2plot', 'var')
    % No specified tracks, so just choose the ones who are closer than
    % threshold initial distance
    tracks2plot = find(dusEuclidean(:, 1) < initialDistanceThres) ;
    
    % Optionally, filter pairs that are not unique
    if uniqueCorrespondence
        disp('Removing duplicate correspondences in candidate tracks... (uniqueCorrespondence)')
        [~, w] = unique( pairIDs, 'stable' );
        duplicate_indices = setdiff( 1:numel(pairIDs), w ) ;
        duplicate_indices = intersect(tracks2plot, duplicate_indices) ;
        tracks2plot = setdiff(tracks2plot, duplicate_indices) ;
    end
else
    % We specified initial tracks2plot. Did we specify these as integer
    % indices or as coordinates in 3d?
    if class(tracks2plot)
        tracks2plot = pointMatch(tracks2plot, squeeze(v3d_u(:, 1, :))) ;
    end
end

% Prepare I/O
if isfield(options, 'specifyTracks')
    subdir = 'specifiedTracks' ;
elseif isfinite(initialDistanceThres)
    subdir = sprintf('initialDistThres%0.2f', initialDistanceThres) ;
else
    subdir = 'allTracks' ;
end
if jitterCorrespondence > 0
    subdir = [subdir '_jitterCorrespondences'] ;
elseif uniqueCorrespondence > 0
    try
        assert(~isfield(options, 'specifyTracks'))
    catch
        error('Cannot both specify tracks and demand unique Correspondences, for intelligibility')
    end
    subdir = [subdir '_uniqueCorrespondences'] ;
end
lateralDir = fullfile(QS.dir.tracking, 'relativeMotionImages', subdir, 'lateral1') ;
ventralDir = fullfile(QS.dir.tracking, 'relativeMotionImages', subdir, 'ventral') ;
lateralVentralDir = fullfile(QS.dir.tracking, 'relativeMotionImages', subdir, 'lateralAndVentral') ;
if ~exist(lateralDir, 'dir')
    mkdir(lateralDir)
end
if ~exist(ventralDir, 'dir')
    mkdir(ventralDir)
end
if ~exist(lateralVentralDir, 'dir')
    mkdir(lateralVentralDir)
end

% Title and timing
t0 = QS.t0set() ;

close all
% Pick colors maximally distinguishable from both black //and white
ntracks2plot = length(tracks2plot) ;
if ntracks2plot < 8
    trackColors = define_colors(ntracks2plot) ;
    trackColors = trackColors ./ max(trackColors, [], 2) ;
else
    trackColors = distinguishable_colors(ntracks2plot, [0,0,0]) ;
    trackColors = trackColors ./ max(trackColors, [], 2) ;
end
colormap(trackColors)

tidx2do = [1, 30, 31, 60] ;%1:50:length(timePoints) ;
tidx2do = [tidx2do, setdiff(1:30:length(timePoints), tidx2do)] ;
tidx2do = [tidx2do, setdiff(1:10:length(timePoints), tidx2do)] ;
tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do)] ;
for tidx = tidx2do
    % Set the current time ------------------------------------------------
    tp = timePoints(tidx) ;
    QS.setTime(tp) 
    
    % Obtain mesh for this timepoint as blackdrop (black backdrop)---------
    % Read in the mesh file 
    disp('Reading mesh...')
    % Specfiy the mesh file to load
    meshfn = sprintf( meshFileBase, tp );
    mesh = read_ply_mod( meshfn );

    % If we smooth before pushing along the normal
    if smoothing_lambda > 0 
        disp('smoothing mesh via laplacian filter')
        mesh.v = laplacian_smooth(...
            mesh.v, mesh.f, 'cotan', [], smoothing_lambda, 'implicit') ;

    end

    % Make sure vertex normals are normalized -----------------------------
    mesh.vn = mesh.vn ./ sqrt( sum( mesh.vn.^2, 2 ) );
    mesh0 = mesh ;
    % Normally evolve vertices
    mesh.v = mesh.v + backdrop_normal_shift .* mesh.vn;
    
    % transform into APDV coordsys
    % Allow for overall flip
    % --> apply rotation and translation and dilation BEFORE flipping
    VV = mesh.v ;
    V0 = mesh0.v ;
    
    opts = struct('Rotation', QS.APDV.rot, ...
        'Translation', QS.APDV.trans, ...
        'Dilation', QS.APDV.resolution) ;
    if isfield(opts, 'Rotation')
        disp('rotating...')
        VV = (opts.Rotation * VV')' ;
        V0 = (opts.Rotation * V0')' ;
    end
    if isfield(opts, 'Translation')
        disp('translating...')
        VV = VV + opts.Translation ;
        V0 = V0 + opts.Translation ;
    end
    if isfield(opts, 'Dilation')
        disp('dilating...')
        VV = VV * opts.Dilation ;
        V0 = V0 * opts.Dilation ;
        dilation = opts.Dilation ;
    else
        dilation = 1 ;
    end
    if flipy
        VV(:, 2) = -VV(:, 2) ;
        V0(:, 2) = -V0(:, 2) ;
    end
    
    
    % Push muscle layer tracks by normal_shift
    vnormals = normals(V0, mesh.f) ;
    vnormals = vnormals ./ vecnorm(vnormals, 2, 2) ;
    faceIndices = ptFaceLocations(tracks2plot, tidx, 1) ;
    nanID = find(isnan(faceIndices)) ;
    faceIndices(nanID) = 1 ;
    norms2add = vnormals(faceIndices, :) ;
    norms2add(nanID, :) = 0 ;
    
    uus = squeeze(v3d_u(tracks2plot, tidx, :)) ;
    vvs = squeeze(v3d_v(tracks2plot, tidx, :)) + ...
        (normal_shift * dilation) .* norms2add ;
    
    % If jitterCorrespondence is nonzero then jitter any points that have
    % multiple pairIDs
    if jitterCorrespondence > 0 
        [~, w] = unique( pairIDs, 'stable' );
        duplicate_indices = setdiff( 1:numel(pairIDs), w ) ;
        duplicate_indices = intersect(tracks2plot, duplicate_indices) ;
        uind2jitter = ismember(tracks2plot, duplicate_indices) ;
        uus(uind2jitter, :) = uus(uind2jitter, :) + ...
            jitterCorrespondence * ...
            rand(size( uus(uind2jitter, :), 1), size( uus(uind2jitter, :), 2)) ;
    elseif uniqueCorrespondence
        % We should have no duplicates by construction. Assert!
        [~, w] = unique( pairIDs, 'stable' );
        duplicate_indices = setdiff( 1:numel(pairIDs), w ) ;
        duplicate_indices = intersect(tracks2plot, duplicate_indices) ;
        assert(isempty(duplicate_indices))
    end
    
    bcs = barycenter(V0, mesh.f) ;
    
    % check it
    % quiver3(bcs(:, 1), bcs(:, 2), bcs(:, 3), ...
    %     vnormals(:, 1), vnormals(:, 2), vnormals(:, 3), 1)

    
    % Alternative mesh
    % mesh = QS.getCurrentSPCutMeshSmRSC() ;
    
    % Which tracks are initially on the left side (near from view(0,0))?
    nearSide = v3d_u(tracks2plot, 1, 2) < 0;
    
    % Plot the tracks at this time ----------------------------------------
    styles = {'textureOverlay', 'backdrop', 'surf', 'nearOnly'} ;
    for styleID = styles2do
        exten = styles{styleID} ;
        latDir = fullfile(lateralDir, styles{styleID}) ;
        venDir = fullfile(ventralDir, styles{styleID}) ;
        latVenDir = fullfile(lateralVentralDir, styles{styleID}) ;
        if ~exist(latDir, 'dir')
            mkdir(latDir)
        end
        if ~exist(venDir, 'dir')
            mkdir(venDir)
        end
        if ~exist(latVenDir, 'dir')
            mkdir(latVenDir)
        end
        
        %% Lateral and ventral together in one image
        outputfnBase = fullfile(latVenDir, 'tracks_LateralAndVentral_%s_%06d.png') ;
        if ~exist(sprintf(outputfnBase, exten, tp), 'file')   || overwrite
            yzlims = [min(xyzlims(2:3, 1)), max(xyzlims(2:3, 2))] ;
            
            clf; 
            sploth = cell(2, 1) ;
            for subplotID = 1:2
                sploth{subplotID} = subplot(2, 1, subplotID) ;
                hold on;

                % make a backdrop mesh for effect
                if styleID ~= 2
                    trisurf(triangulation(mesh.f, VV + [0, 200, 0]), 'facecolor', 'k')
                    trisurf(triangulation(mesh.f, VV + [0, 0, 200]), 'facecolor', 'k')
                else
                    % make a foreground mesh for opacity/depth perception
                    trisurf(triangulation(mesh.f, V0), 'facecolor', 'k', ...
                        'edgecolor', 'none', 'facealpha', 0.2)
                end

                for ii = 1:length(tracks2plot)
                    if styleID < 4 || nearSide(ii)
                        % trackii = tracks2plot(ii) ;
                        scatter3(uus(ii, 1), uus(ii, 2), uus(ii, 3), ...
                            50, 'filled', 'markerfacecolor', ...
                            trackColors(ii, :), ...
                            'markeredgecolor', 'none') ;
                        scatter3(vvs(ii, 1), vvs(ii, 2), vvs(ii, 3), ...
                            50, '^', 'filled', 'markerfacecolor',...
                            trackColors(ii, :), ...
                            'markeredgecolor', 'none') ;
                    end
                end

                % Plot trailing track
                if tidx > 1
                    for ii = 1:length(tracks2plot)
                        if styleID < 3 || nearSide(ii)
                            trackid = tracks2plot(ii) ;
                            colorii = trackColors(ii, :) ;
                            plot3(v3d_u(trackid, 1:tidx, 1), ...
                                v3d_u(trackid, 1:tidx, 2), ...
                                v3d_u(trackid, 1:tidx, 3), '-', ...
                                'color', colorii, 'linewidth', lwTrail)
                            plot3(v3d_v(trackid, 1:tidx, 1), ...
                                v3d_v(trackid, 1:tidx, 2), ...
                                v3d_v(trackid, 1:tidx, 3), '-', ...
                                'color', colorii, 'linewidth', lwTrail)
                        end
                    end
                end
            
                if subplotID == 1
                    % Lateral view
                    view(0, 0)

                    % Format the figure
                    axis equal
                    xlim(xyzlims(1, :))
                    ylim(yzlims + [0,  200])
                    zlim(yzlims)
                    ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    timeStamp = num2str((timePoints(tidx) - t0) * QS.timeInterval) ;
                    title(['nuclear tracks, $t=$' timeStamp ' ' QS.timeUnits], 'interpreter', 'latex')
                    grid off
                    
                    sh1 = scatter3(NaN, NaN, NaN, ...
                        50, 'filled', 'markerfacecolor', ...
                        'k', ...
                        'markeredgecolor', 'none') ;
                    sh2 = scatter3(NaN, NaN, NaN, ...
                        50, '^', 'filled', 'markerfacecolor','k', ...
                        'markeredgecolor', 'none') ;
                    lh = legend([sh1, sh2], {'endoderm', 'muscle'},...
                        'position', [0.80    0.7182    0.1742    0.0724], ...
                        'interpreter', 'latex') ;
                    % drawnow
                    % axpos = get(sploth, 'position') ;
                else

                    % Ventral view
                    view(0, 270)

                    % Format the figure
                    axis equal
                    xlim(xyzlims(1, :))
                    ylim(yzlims)
                    zlim(yzlims + [0, 200])
                    xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
                    grid off
                    % drawnow
                    % ax2pos = get(sploth, 'position') ;
                    % try
                    %     assert(axpos(1) == ax2pos(1))
                    %     assert(axpos(3) == ax2pos(3))
                    %     assert(axpos(4) == ax2pos(4))
                    % catch
                    %     disp('Axes positions unequal, adjusting...')
                    %     set(gca, 'position', [axpos(1),  ax2pos(2), axpos(3), ax2pos(4)])
                    % end
                    % drawnow
                end
                
            end
            
            % Align axes
            axFixPos = get(sploth{1}, 'position') ;
            for subplotIDj = 1:2
                axpos = get(sploth{subplotIDj}, 'position') ;
                set(sploth{subplotIDj}, 'position', ...
                    [axFixPos(1), axpos(2), axFixPos(3), axFixPos(4)])
            end
            
            % Save the figure
            set(gcf,'color','w');
            outfn = sprintf(outputfnBase, exten, tp) ;
            disp(['Saving figure: ' outfn])
            export_fig(outfn, '-nocrop', '-r150')
            
        end
        
        %% ventral and lateral separately, for overlay with fluor data.
        outputfnLateral = fullfile(latDir, 'tracks_lateral_%s_%06d.png') ;
        outputfnVentral = fullfile(venDir, 'tracks_ventral_%s_%06d.png') ;
        if ~exist(sprintf(outputfnLateral, exten, tp), 'file') || ...
            ~exist(sprintf(outputfnVentral, exten, tp), 'file') || overwrite

            close all
            
            fig = figure('visible', 'off', ...
                'units', 'centimeters', 'position', [0, 0, 16, 10]) ;
            clf; hold on;

            % make a backdrop mesh for effect
            if styleID ~= 3 && styleID ~= 1
                trisurf(triangulation(mesh.f, VV + [0, 200, 0]), 'facecolor', 'k')
                trisurf(triangulation(mesh.f, VV + [0, 0, 200]), 'facecolor', 'k')
            elseif styleID > 1
                % make a foreground mesh for opacity/depth perception
                trisurf(triangulation(mesh.f, V0), 'facecolor', 'k', ...
                    'edgecolor', 'none', 'facealpha', 0.2)
            end

            for ii = 1:length(tracks2plot)
                if styleID < 3 || nearSide(ii)
                    % trackii = tracks2plot(ii) ;
                    scatter3(uus(ii, 1), uus(ii, 2), uus(ii, 3), ...
                        50, 'filled', 'markerfacecolor', ...
                        trackColors(ii, :), ...
                        'markeredgecolor', 'none')
                    scatter3(vvs(ii, 1), vvs(ii, 2), vvs(ii, 3), ...
                        50, '^', 'filled', 'markerfacecolor', ...
                        trackColors(ii, :), ...
                        'markeredgecolor', 'none')
                end
            end

            % Plot trailing track
            if tidx > 1
                for ii = 1:length(tracks2plot)
                    if styleID < 3 || nearSide(ii)
                        trackid = tracks2plot(ii) ;
                        colorii = trackColors(ii, :) ;
                        plot3(v3d_u(trackid, 1:tidx, 1), ...
                            v3d_u(trackid, 1:tidx, 2), ...
                            v3d_u(trackid, 1:tidx, 3), '-', ...
                            'color', colorii, 'linewidth', lwTrail)
                        plot3(v3d_v(trackid, 1:tidx, 1), ...
                            v3d_v(trackid, 1:tidx, 2), ...
                            v3d_v(trackid, 1:tidx, 3), '-', ...
                            'color', colorii, 'linewidth', lwTrail)
                    end
                end
            end
            % Lateral view
            view(0, 0)

            % Format the figure
            set(gcf,'color','w');
            if styleID == 1
                set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w', 'zcol', 'w')
                set(gcf, 'InvertHardCopy', 'off');
                set(gca, 'color', 'k')
                set(gcf, 'color', 'k')
            end
            axis equal
            xlim(xyzlims(1, :))
            ylim([xyzlims(2, 1), xyzlims(2, 2) + 200])
            zlim(xyzlims(3, :))
            xlabel(['AP position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            zlabel(['DV position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            timeStamp = num2str((timePoints(tidx) - t0) * QS.timeInterval) ;
            
            titlestr = ['$t=$' timeStamp ' ' QS.timeUnits] ;
            if styleID == 1
                title(titlestr, 'Interpreter', 'Latex', 'Color', 'white') 
            else
                title(titlestr, 'Interpreter', 'Latex', 'Color', 'k') 
            end
            grid off
                    
            set(fig, 'PaperUnits', 'centimeters')
            set(fig, 'PaperPosition', [0,0,16,10])
            
            % Save the figure
            outfnL = sprintf(outputfnLateral, exten, tp) ;
            disp(['Saving figure: ' outfnL])
            export_fig(outfnL, '-nocrop', '-r200')

            % Ventral view
            view(0, 270)

            % Format the figure
            axis equal
            xlim(xyzlims(1, :))
            ylim(xyzlims(2, :))
            zlim(xyzlims(3, :) + [0, 200])
            xlabel(['ap position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            ylabel(['lateral position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            zlabel(['dv position [' QS.spaceUnits ']'], 'interpreter', 'latex')
            timeStamp = num2str((timePoints(tidx) - t0) * QS.timeInterval) ;
            title(['$t=$' timeStamp ' ' QS.timeUnits], 'interpreter', 'latex')
            grid off
            
            % Save the figure
            outfnV = sprintf(outputfnVentral, exten, tp) ;
            disp(['Saving figure: ' outfnV])
            export_fig(outfnV, '-nocrop', '-r200')
            
            
            % Overlay with texturepatch
            if styleID == 1 
                
                % Create texturepatch image in subdir
                disp('Creating texture patch image for layer 1...')
                metafn = fullfile(QS.dir.texturePatchIm, ...
                    subdir_layer1, 'metadat.mat') ;
                load(metafn, 'metadat', 'Options')
                opts = metadat ;
                opts.subdir = subdir_layer1 ; 
                opts.timePoints = tp ;
                opts.plot_perspective = false ;
                opts.plot_dorsal = false ;
                opts.plot_ventral = true ;
                opts.plot_left = true ;
                opts.plot_right = false ;
                opts.blackFigure = true ;
                opts.makeColorbar = false ;
                QS.clearTime()
                QS.setTime(tp) ;
                QS.data.adjustlow = layer1_adjustDataLim(1) ;
                QS.data.adjusthigh = layer1_adjustDataLim(2) ;
                % QS.setDataLimits(QS.xp.fileMeta.timePoints(1), 1.0, 99.0)
                QS.plotSeriesOnSurfaceTexturePatch(opts, Options)
                
                % Create texturepatch image in subdir
                disp('Creating texture patch image for layer 2...')
                metafn = fullfile(QS.dir.texturePatchIm, ...
                    subdir_layer2, 'metadat.mat') ;
                load(metafn, 'metadat', 'Options')
                opts = metadat ;
                opts.subdir = subdir_layer2 ; 
                opts.timePoints = tp ;
                opts.plot_perspective = false ;
                opts.plot_dorsal = false ;
                opts.plot_ventral = true ;
                opts.plot_left = true ;
                opts.plot_right = false ;
                opts.blackFigure = true ;
                opts.makeColorbar = false ;
                QS.clearTime()
                QS.setTime(tp) ;
                QS.data.adjustlow = layer2_adjustDataLim(1) ;
                QS.data.adjusthigh = layer2_adjustDataLim(2) ;
                % QS.setDataLimits(QS.xp.fileMeta.timePoints(1), 1.0, 99.995)
                QS.plotSeriesOnSurfaceTexturePatch(opts, Options)
                    
                lim = imread(outfnL) ;
                vim = imread(outfnV) ;
                lat1fn = fullfile(QS.dir.texturePatchIm, ...
                    subdir_layer1, 'lateral1', ...
                    sprintf('patch_lateral1_%06d.png', tp)) ;
                lat2fn = fullfile(QS.dir.texturePatchIm, ...
                    subdir_layer2, 'lateral1', ...
                    sprintf('patch_lateral1_%06d.png', tp)) ;
                ven1fn = fullfile(QS.dir.texturePatchIm, ...
                    subdir_layer1, 'ventral', ...
                    sprintf('patch_ventral_%06d.png', tp)) ;
                ven2fn = fullfile(QS.dir.texturePatchIm, ...
                    subdir_layer2, 'ventral', ...
                    sprintf('patch_ventral_%06d.png', tp)) ;
                tlim1 = imread(lat1fn) ;
                tlim2 = imread(lat2fn) ;
                tvim1 = imread(ven1fn) ;
                tvim2 = imread(ven2fn) ;
                % summed texture images
                stlim = textureFactor * (tlim1 + tlim2) ;
                stvim = textureFactor * (tvim1 + tvim2) ;
                
                %% Lateral correction
                % GET ROI for adding images inside canvas
                axisROIFn = fullfile(lateralDir, 'embeddingROI_for_overlay.mat') ;
                if ~exist(axisROIFn, 'file')
                    axisROI = roipoly(stlim) ;
                    save(axisROIFn, 'axisROI')    
                else
                    load(axisROIFn, 'axisROI')
                end
                lroi = axisROI ;
                
                % sum together inside axisROI
                im2L = tlim2 ;
                sRL = squeeze(stlim(:, :, 1)) ; % summed image
                sGL = squeeze(stlim(:, :, 2)) ;
                sBL = squeeze(stlim(:, :, 3)) ;
                im2R = squeeze(im2L(:,:,1)) ;
                im2G = squeeze(im2L(:,:,1)) ;
                im2B = squeeze(im2L(:,:,1)) ;
                im2R(lroi) = sRL(lroi) ;
                im2G(lroi) = sGL(lroi) ;
                im2B(lroi) = sBL(lroi) ;
                
                % Add track image inside ROI
                trackMeanRGB = mean(lim, 3) ;
                maskTrack = trackMeanRGB > 5 ;
                trackimR = squeeze(lim(:, :, 1)) ;
                trackimG = squeeze(lim(:, :, 2)) ;
                trackimB = squeeze(lim(:, :, 3)) ;
                im2R(maskTrack) = trackimR(maskTrack) ;
                im2G(maskTrack) = trackimG(maskTrack) ;
                im2B(maskTrack) = trackimB(maskTrack) ;
                im2L = cat(3, im2R, im2G, im2B) ;
                imwrite(im2L, outfnL)

                %% VENTRAL correction
                axisROIFn = fullfile(ventralDir, 'embeddingROI_for_overlay.mat') ;
                if ~exist(axisROIFn, 'file')
                    axisROI = roipoly(stvim) ;
                    save(axisROIFn, 'axisROI')    
                else
                    load(axisROIFn, 'axisROI')
                end
                vroi = axisROI ;

                % sum together inside axisROI
                im2V = tvim2 ;
                sRV = squeeze(stvim(:, :, 1)) ; % summed image
                sGV = squeeze(stvim(:, :, 2)) ;
                sBV = squeeze(stvim(:, :, 3)) ;
                im2R = squeeze(im2V(:,:,1)) ;
                im2G = squeeze(im2V(:,:,1)) ;
                im2B = squeeze(im2V(:,:,1)) ;
                im2R(vroi) = sRV(vroi) ;
                im2G(vroi) = sGV(vroi) ;
                im2B(vroi) = sBV(vroi) ;
                
                % Add track image inside ROI
                trackMeanRGB = mean(vim, 3) ;
                maskTrack = trackMeanRGB > 5 ;
                trackimR = squeeze(vim(:, :, 1)) ;
                trackimG = squeeze(vim(:, :, 2)) ;
                trackimB = squeeze(vim(:, :, 3)) ;
                im2R(maskTrack) = trackimR(maskTrack) ;
                im2G(maskTrack) = trackimG(maskTrack) ;
                im2B(maskTrack) = trackimB(maskTrack) ;
                im2V = cat(3, im2R, im2G, im2B) ;
                imwrite(im2V, outfnV)
            end
            
        end
    end
end
