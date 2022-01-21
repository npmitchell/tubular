function [dusEuclidean, dusGeodesic, ...
        tracks1, tracks2, pairIDs, U0s, V0s, ...
        geodesicPaths, ptBarycenters, ptFaceLocations, ...
        euclideanDistanceTraveledU, euclideanDistanceTraveledV, ...
        v3d_u, v3d_v, nSaved] = ...
        measureRelativeMotionTracks(QS, track1fn, track2fn, Options)
%measureRelativeMotionBetweenCellTracks(QS, track1fn, track2fn, Options)
%
% Measure the relative motion of a set of tracks2 (ex, muscle) against the
% tracks in tracks1 (ex, endoderm). 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% track1fn : str
%   full path to mat file with set of tracks to compare against
%   (endoder/lagrangian tracks) in 2d pullback coordinates
%   (Options.coordSys assigns the coordinate system used and whether single
%   or double cover of that system)
% track2fn : str
%   full path to other set of tracks whose motion we measure against
%   tracks1, in 2d pullback coordinates
%   (Options.coordSys assigns the coordinate system used and whether single
%   or double cover of that system)
% Options : struct with fields (optional)
%   timePoints : nTimePoints x 1 numeric array
%       timepoints of the tracks, so that track{ii}(tidx, :) is the
%       XY position of track ii at time timePoints(tidx). Note this cannot 
%       be some arbitrary subset of the experiment timepoints
%       QS.xp.fileMeta.timePoints but must instead match the true track 
%       info!
%   max_du : float in units of QS.spaceUnits (optional, default=Inf)
%       maximum possible physical displacement between nuclei, for 
%       filtering out unreasonable results
%   coordSys : str specifier ('spsme', 'spsm', 'uv', 'uve', 'ricci', etc)
%       coordinate system in which 2d tracks UV are embedded
%   
%
% Returns 
% -------
% dusEuclidean : Euclidean distance between 
% dusGeodesic, ...
% tracks1 : cell array of #timepoints x 2 float/numeric arrays
%   loaded tracks from disk for set 1/A
% tracks2 : cell array of #timepoints x 2 float/numeric arrays
%   loaded tracks from disk for set 2/B
% pairIDs :
% U0s :
% V0s :
% geodesicPaths
% ptBarycenters
% ptFaceLocations : 
% euclideanDistanceTraveledU 
% euclideanDistanceTraveledV 
% v3d_u :
% v3d_v :
% nSaved
% 
%
% NPMitchell 2021

relMotionFn = fullfile(QS.dir.tracking, 'relative_motion_tracks.mat') ;
relMotionStartFn = fullfile(QS.dir.tracking, 'relative_motion_trackStart.png') ;
overwrite = false ;
preview = false ;
subdir1 = Options.subdir1 ;
subdir2 = Options.subdir2 ;

if nargin < 4
    Options = struct() ;
end
if isfield(Options, 'overwrite')
    overwrite = Options.overwrite ;
end

tmp1 = load(track1fn) ;
tmp2 = load(track2fn) ;

if ~isfield(tmp1, 'tracks') && isfield(tmp1, 'GG')
    tracks1 = trackingGraph2Cell(tmp1.GG) ;
else
    tracks1 = tmp1.tracks ;
end
if ~isfield(tmp2, 'tracks') && isfield(tmp2, 'GG')
    tracks2 = trackingGraph2Cell(tmp2.GG) ;
else
    tracks2 = tmp2.tracks ;
end

n1 = length(tracks1) ;
n2 = length(tracks2) ;

% Options handling
doubleCovered = false ;
umax = 1.0 ;
vmax = 1.0 ;
coordSys = 'spsmrs' ;
max_du = Inf ;
timePoints = QS.xp.fileMeta.timePoints ;
if isfield(Options, 'doubleCovered')
    doubleCovered = Options.doubleCovered ;
end
if isfield(Options, 'umax')
    umax = Options.umax ;
end
if isfield(Options, 'vmax')
    vmax = Options.vmax ;
end
if isfield(Options, 'max_du')
    max_du = Options.max_du ;
end
if isfield(Options, 'coordSys')
    coordSys = Options.coordSys ;
end
if isfield(Options, 'timePoints')
    timePoints = Options.timePoints ;
end

% Colormap and shorthand
nTracks = length(tracks2) ;
colors = jetshuffle(nTracks) ;
nTimePoints = length(timePoints) ;

% Pre-allocate positions of starting points in tracks of tracks1 (endoderm)
nearby = zeros(n1, 2) ;
for trackID = 1:n1
    nearby(trackID, :) = tracks1{trackID}(1, 1:2) ;
end

overwrite = true ;
if exist(relMotionFn, 'file') && ~overwrite
    load(relMotionFn, 'dusEuclidean', 'dusGeodesic', ...
        'tracks1', 'tracks2', 'pairIDs', 'U0s', 'V0s', ...
        'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
        'euclideanDistanceTraveledU', 'euclideanDistanceTraveledV', ...
        'v3d_u', 'v3d_v', 'nSaved') 
    % track1 3 & track2 191: t = 17 geodesic distance
    if exist(relMotionFn, 'file')
        movefile(relMotionFn, [relMotionFn '_backup'])
    end

else
    % For each track in tracks2, find initially nearest track in tracks1
    dusGeodesic = nan(nTracks, nTimePoints) ;
    dusEuclidean = nan(nTracks, nTimePoints) ;
    pairIDs = zeros(nTracks, 1) ;
    U0s = zeros(nTracks, 2) ;
    V0s = zeros(nTracks, 2) ;
    geodesicPaths = cell(nTracks, 1) ;
    ptBarycenters = nan(nTracks, nTimePoints, 2, 3) ;
    ptFaceLocations = nan(nTracks, nTimePoints, 2) ;
    euclideanDistanceTraveledU = dusGeodesic ;
    euclideanDistanceTraveledV = dusGeodesic ;
    v3d_u = zeros(nTracks, nTimePoints, 3) ;
    v3d_v = zeros(nTracks, nTimePoints, 3) ;
    nSaved = 1 ;
end

for ii = nSaved:nTracks
    disp(['matching track ' num2str(ii)])

    % Get starting positions u0 and v0 for layer1 and 2
    V0 = tracks2{ii}(1, 1:2) ;
    V0s(ii, :) = V0 ;
    VV = tracks2{ii}(:, 1:2) ;
    
    % Look for initially nearby nuclei in pullback space
    dists = vecnorm(nearby - V0, 2, 2) ;
    [~, minID] = min(dists) ;
    UU = tracks1{minID}(:, 1:2) ;
    U0s(ii, :) = UU(1, 1:2) ;

    % Project into 3D
    if strcmpi(coordSys, 'spsm') || strcmpi(coordSys, 'spsmrs')
        im = imread(fullfile(QS.dir.im_sp_sm, subdir1, ...
            sprintf(QS.fileBase.im_sp_sm, timePoints(1)))) ;
        im2 = imread(fullfile(QS.dir.im_sp_sm, subdir2, ...
            sprintf(QS.fileBase.im_sp_sm, timePoints(1)))) ;
        assert(all(size(im) == size(im2)))
    end
    uu = QS.XY2uv(im, UU, doubleCovered, umax, vmax) ;
    vv = QS.XY2uv(im, VV, doubleCovered, umax, vmax) ;

    % For each timepoint, project into 3d for that mesh
    geodesics = cell(nTimePoints, 1) ;
    ptBarys = nan(nTimePoints, 2, 3) ;
    ptFaces = nan(nTimePoints, 2) ;
    duEuclidean = zeros(nTimePoints, 1) ;
    duGeodesic = zeros(nTimePoints, 1) ;
    u3ds = zeros(nTimePoints, 3) ;
    v3ds = zeros(nTimePoints, 3) ;
    for tidx = 1:nTimePoints
        tp = timePoints(tidx) ;
        QS.setTime(tp) 

        if ~any(isnan(uu(tidx, :))) && ~any(isnan(vv(tidx, :))) && ...
            all(uu(tidx, :) > 0) && all(vv(tidx, :) > 0)
            [u3d, fieldfacesU, ~, barycU] = QS.uv2APDV(uu(tidx, :), coordSys) ;
            [v3d, fieldfacesV, ~, barycV] = QS.uv2APDV(vv(tidx, :), coordSys) ;

            % Note: if we define mesh as:
            % mesh = QS.getCurrentSPCutMeshSmRS() ;
            %  with normalization such that
            % mesh.u(:, 1) = mesh.u(:, 1) * umax / max(mesh.u(:, 1)) ;
            %  then we can use baryc to recover the position in pullback
            % uu_check = [sum(barycU .* mesh.u(mesh.f(fieldfacesU, :), 1)', 2), ...
            %            sum(barycU .* mesh.u(mesh.f(fieldfacesU, :), 2)', 2)] ;
            % assert(all(uu_check == uu(tidx, :)))

            if tidx == 1 
                u3d0 = u3d ;
                v3d0 = v3d ;
                fieldfacesU0 = fieldfacesU ;
                fieldfacesV0 = fieldfacesV ;
            end

            progressString = ['trackA ' num2str(ii) ' & trackB ' num2str(minID)...
                    ': t = ' num2str(timePoints(tidx)) ...
                    ' | U=[' sprintf('%0.0f,%0.0f', UU(tidx, 1), UU(tidx, 2)) ...
                    '], V=[' sprintf('%0.0f,%0.0f', VV(tidx, 1), VV(tidx, 2)) ']'] ;

            % Euclidean distance in 3d     
            duEuclidean(tidx) = vecnorm(u3d-v3d, 2, 2) ;

            if fieldfacesU == fieldfacesV 
                % Points are very close --> on the same face: 
                % Euclidean distance in 3D is equal to the geodesic 
                % distance      
                disp(progressString)
                duGeodesic(tidx) = duEuclidean(tidx) ;
            else
                % Measure distance of these tracks over surface as geodesic
                disp([progressString ' geodesic distance'])
                mesh = QS.getCurrentSPCutMeshSmRSC() ;
                nvtx = size(mesh.v, 1) ;

                % Acquire the geodesics
                [geodesicP, pointLocations] = surfaceGeodesicPairs( mesh.f, ...
                        mesh.v, [nvtx+1, nvtx+2], [mesh.v; u3d; v3d] ) ;
                ptFace = [ pointLocations(nvtx+1).face, ...
                    pointLocations(nvtx+2).face ] ;
                assert(all(ptFace == [fieldfacesU, fieldfacesV] ))

                % note: ptBary definition here has transpositions. Use
                % barycU and barcV as defined above instead.
                % ptBary = [ pointLocations(nvtx+1).barycentricCoordinates; ...
                %     pointLocations(nvtx+2).barycentricCoordinates ] ;
                % assert(all(barycU == ptBary(1, :)))
                % assert(all(barycV == ptBary(2, :)))

                pp = geodesicP{1} ;
                duGeodesic(tidx) = sum(vecnorm(diff(pp), 2, 2)) ;

                geodesics{tidx} = pp ;                    
            end
            ptBarys(tidx, :, :) = [barycU; barycV] ;
            ptFaces(tidx, :) = [fieldfacesU, fieldfacesV] ;

            % Keep 3d positions
            u3ds(tidx, :) = u3d ;
            v3ds(tidx, :) = v3d ;
        else
            disp('NaN for one track!')
            duGeodesic(tidx) = NaN ;
            u3ds(tidx, :) = NaN ;
            v3ds(tidx, :) = NaN ;
        end
    end

    % Compare to distance traveled
    % smooth u3d a little
    u3d_sm = u3ds ;
    u3d_sm(:, 1) = smoothdata(u3d_sm(:, 1),'sgolay');
    u3d_sm(:, 2) = smoothdata(u3d_sm(:, 2),'sgolay');
    u3d_sm(:, 3) = smoothdata(u3d_sm(:, 3),'sgolay');
    v3d_sm = v3ds ;
    v3d_sm(:, 1) = smoothdata(v3d_sm(:, 1), 'sgolay');
    v3d_sm(:, 2) = smoothdata(v3d_sm(:, 2), 'sgolay');
    v3d_sm(:, 3) = smoothdata(v3d_sm(:, 3), 'sgolay');
    Deltau = cumsum(vecnorm(diff(u3d_sm), 2, 2)) ;
    Deltav = cumsum(vecnorm(diff(v3d_sm), 2, 2)) ;
    Deltau = [0; Deltau ] ;
    Deltav = [0; Deltav ] ;

    % Store all these vectors
    geodesicPaths{ii} = geodesics ;
    ptBarycenters(ii, :, :, :) = ptBarys ;
    ptFaceLocations(ii, :, :) = ptFaces ;
    euclideanDistanceTraveledU(ii, :) = Deltau ;
    euclideanDistanceTraveledV(ii, :) = Deltav ;
    dusEuclidean(ii, :) = duEuclidean ;
    dusGeodesic(ii, :) = duGeodesic ;
    pairIDs(ii) = minID ;
    v3d_u(ii, :, :) = u3ds ;
    v3d_v(ii, :, :) = v3ds ;
    nSaved = ii ;
    
    
    %% Filter bad points
    % Bad row, columns == [25, 26], [42, 42]
    if ~isnan(max_du) && isfinite(max_du)
        if any(dusGeodesic(:) > max_du) || any(dusEuclidean(:) > max_du) 
            response = input('Warning: sure you want to filter out large distances?', 's') ;
            if contains(lower(response), 'y')
                disp('filtering...')
                dusGeodesic(dusGeodesic(:) > max_du) = NaN ;
                dusEuclidean(dusEuclidean(:) > max_du) = NaN ;
            end
        end
    end
    
    %% Save it
    disp(['Saving track pairs [' num2str(ii) '/' num2str(nTracks) '] to ' relMotionFn])
    save(relMotionFn, 'dusEuclidean', 'dusGeodesic', ...
        'tracks1', 'tracks2', 'pairIDs', 'U0s', 'V0s', ...
        'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
        'euclideanDistanceTraveledU', 'euclideanDistanceTraveledV', ...
        'v3d_u', 'v3d_v', 'nSaved') 

    %% Preview the resulting initial match
    Xlim = U0s(ii, 1) + [-200, 200] ;
    Ylim = U0s(ii, 2) + [-200, 200] ;

    % Plot starting correspondences
    if preview 
        clf
        sz = 200 ;
        % Axis 1
        ax0 = subtightplot(2, 2, 1) ;
        imshow(0.25*im);
        hold on;
        scatter(nearby(:, 1), nearby(:, 2), sz, [0, 0.447, 0.741], 'filled')
        scatter(U0s(:, 1), U0s(:, 2), sz, colors, 'filled', 'markeredgecolor', 'k')
        title([subdir1 ' positions'])
        xlim(Xlim)
        ylim(Ylim)
        % Axis 2
        ax1 = subtightplot(2, 2, 2) ;
        imshow(im2);
        hold on;
        scatter(V0s(:, 1), V0s(:, 2), sz*4, colors, 's', 'linewidth', 5)
        title([subdir2 ' positions'])
        xlim(Xlim)
        ylim(Ylim)
        % Axis 3
        ax2 = subtightplot(2, 1, 2) ;
        imshow(min(0.25*im+im2, 255));
        hold on;
        scatter(nearby(:, 1), nearby(:, 2), sz, [0, 0.447, 0.741], 'filled')
        scatter(V0s(:, 1), V0s(:, 2), sz*4, colors, 's', 'linewidth', 5)
        scatter(U0s(:, 1), U0s(:, 2), sz, colors, 'filled', 'markeredgecolor', 'k')
        for jj = 1:nTracks
            plot([U0s(jj, 1), V0s(jj, 1)], [U0s(jj, 2), V0s(jj, 2)], '-', ...
                'color', colors(jj, :))
        end
        hold off;
        xlim(Xlim)
        ylim(Ylim)
        pause(1)
    end
end

%% Plot all starting correspondences
clf
sz = 20 ;
% Axis 1
ax0 = subtightplot(2, 2, 1) ;
imshow(0.25*im);
hold on;
scatter(nearby(:, 1), nearby(:, 2), sz*0.2, [0, 0.447, 0.741])
scatter(U0s(:, 1), U0s(:, 2), sz, colors, 'filled')
title([subdir1 ' positions'])
% Axis 2
ax1 = subtightplot(2, 2, 2) ;
imshow(im2);
hold on;
scatter(V0s(:, 1), V0s(:, 2), sz*2, colors, 's')
title([subdir2 ' positions'])
% Axis 3
ax2 = subtightplot(2, 1, 2) ;
imshow(min(0.25*im+im2, 255));
hold on;
scatter(U0s(:, 1), U0s(:, 2), sz, colors, 'filled')
scatter(V0s(:, 1), V0s(:, 2), sz*2, colors, 's')
for jj = 1:nTracks
    plot([U0s(jj, 1), V0s(jj, 1)], [U0s(jj, 2), V0s(jj, 2)], '-', ...
        'color', colors(jj, :))
end
hold off;
saveas(gcf, relMotionStartFn)

