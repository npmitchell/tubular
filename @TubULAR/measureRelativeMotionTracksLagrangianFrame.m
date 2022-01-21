function measureRelativeMotionTracksLagrangianFrame(QS, Options)
%
%
%
%


relMotionFn = fullfile(QS.dir.tracking, 'relative_motion_tracks.mat') ;
relMotionLagrangianFn = fullfile(QS.dir.tracking, 'relative_motion_tracks_Lagrangian.mat') ;
overwrite = false ;
tracksDoubleCovered = false ;
pathlinesDoubleCovered = true ;
preview = false ;
xlimit = 30 ;
tcut1 = 0.5 ;
tcut2 = 1 ;
tcut3 = 1.5 ;

if isfield(Options, 'overwrite')
    overwrite = Options.overwrite ;
end
if isfield(Options, 'preview')
    preview = Options.preview ;
end
if isfield(Options, 'tracksDoubleCovered')
    tracksDoubleCovered = Options.tracksDoubleCovered ;
end
if isfield(Options, 'pathlinesDoubleCovered')
    pathlinesDoubleCovered = Options.pathlinesDoubleCovered ;
end
if isfield(Options, 'xlimit')
    xlimit = Options.xlimit ;
end

%% Interpret options
singleCover = ~tracksDoubleCovered ;

%% Load data to plot
disp(['Loading track pairs from disk: ' relMotionFn])
load(relMotionFn, 'tracks1', 'tracks2', 'pairIDs') ;

n1 = length(tracks1) ;
n2 = length(tracks2) ;
nTimePoints = size(tracks1{1}, 1) ;
timePoints = QS.xp.fileMeta.timePoints ;
assert(length(timePoints) == nTimePoints) 
nTracks = n2 ;



%% Store in-plane positions as arrays -- careful Lagrangian frame method
pausetime = 0.01 ;

% Load pathlines and the ricci mesh (reference mesh)
pline = QS.getPullbackPathlines() ;
refMesh = pline.refMesh ;
% Convert to APDV coordinates
refMesh.v = QS.xyz2APDV(refMesh.v) ;

if ~exist(relMotionLagrangianFn, 'file') || overwrite
    % Ricci 2d positions pre-allocate
    UallRicci = nan(nTracks, nTimePoints, 2) ;
    VallRicci = nan(nTracks, nTimePoints, 2) ;
    U0Riccis = zeros(nTracks, 2) ;
    V0Riccis = zeros(nTracks, 2) ;
    UallOriginal = UallRicci ;
    VallOriginal = VallRicci ;

    % geodesic info pre-allocate
    dusGeodesic = nan(nTracks, nTimePoints) ;
    geodesicPaths = cell(nTracks, 1) ;
    ptBarycenters = nan(nTracks, nTimePoints, 2, 3) ;
    ptFaceLocations = nan(nTracks, nTimePoints, 2) ;
    lagrangianDistanceTraveledU = dusGeodesic ;
    lagrangianDistanceTraveledV = dusGeodesic ;
    v3d_uLagr = zeros(nTracks, nTimePoints, 3) ;
    v3d_vLagr = zeros(nTracks, nTimePoints, 3) ;
    v3d_uLagr_smoothed = zeros(nTracks, nTimePoints, 3) ;
    v3d_vLagr_smoothed = zeros(nTracks, nTimePoints, 3) ;

    phaseRaw = nan(nTracks, nTimePoints) ;
    magXYRaw = nan(nTracks, nTimePoints) ;
    phaseRelative = nan(nTracks, nTimePoints) ;

    for ii = 1:nTracks
        disp(['following track pair ' num2str(ii) ' in pullback pathlines'])

        % Get starting positions u0 and v0 for layer1 and 2
        UU = tracks1{pairIDs(ii)}(:, 1:2) ;
        VV = tracks2{ii}(:, 1:2) ;
        UallOriginal(ii, :, :) = UU ; 
        VallOriginal(ii, :, :) = VV ; 

        % Compute U_ricci and V_ricci for each timepoint
        for tidx = 1:nTimePoints
            QS.setTime(timePoints(tidx)) ;
            % mesh = QS.getCurrentSPCutMeshSmRSC() ;
            % [gcell, bcell] = constructFundamentalForms(mesh.f, mesh.v, mesh.u) ;
            % [~, v0t, v0t2d] = ...
            %     resolveTangentNormalVelocities(faces, vertices, v0s, fieldfaces, ...
            %         mesh.u) ;

            vX = pline.vertices.vX(tidx, :, :) ;
            vY = pline.vertices.vY(tidx, :, :) ;
            uu = UU(tidx, :) ;
            vv = VV(tidx, :) ;

            if preview
                figure(1)
                clf
                plot([0,0], [0, 1000], 'k') ; hold on;
                plot([0,2000], [1000, 1000], 'k')
                plot([2000,2000], [0, 1000], 'k')
                plot([0,2000], [0, 0], 'k')
                scatter(uu(:, 1), uu(:, 2), 'filled', 'markerfacecolor', 'r')
                scatter(vv(:, 1), vv(:, 2), 'filled', 'markerfacecolor', 'g')
            end

            if singleCover 
                uu(:, 2) = uu(:, 2) + min(pline.refMesh.XY(:, 2)) - 1 ;
                vv(:, 2) = vv(:, 2) + min(pline.refMesh.XY(:, 2)) - 1 ;
            end

            tridx = triangulation(pline.refMesh.f, [vX(:), vY(:)]) ;

            if any(isnan(uu)) || any(isnan(vv))
                fieldfacesU = NaN ;
                barycU = [NaN, NaN, NaN] ;
                fieldfacesV = NaN ;
                barycV = [NaN, NaN, NaN] ;
                U_ricci = [NaN, NaN ] ;
                V_ricci = [NaN, NaN ] ;
            else
                [fieldfacesU, barycU] = pointLocation(tridx, uu) ;
                [fieldfacesV, barycV] = pointLocation(tridx, vv) ;

                % Inspect that we are sampling the mesh in the right spot        
                if preview
                    figure(2)
                    clf
                    triplot(tridx, 'color', 'k') ; hold on;
                    scatter(uu(:, 1), uu(:, 2), 'filled', 'markerfacecolor', 'r')
                    scatter(vv(:, 1), vv(:, 2), 'filled', 'markerfacecolor', 'g')
                    pause(pausetime)
                end

                % Find the right branchcut in which the point resides and assign
                % fieldfaces and barycenters in lagrangian frame (advected pathline
                % mesh)
                if isnan(fieldfacesU) || isnan(fieldfacesV)
                    if pathlinesDoubleCovered && singleCover
                        Ly2add = pline.vertices.Ly(tidx) * 0.5 ;
                    elseif pathlinesDoubleCovered && ~singleCover
                        Ly2add = pline.vertices.Ly(tidx) ;
                    elseif ~pathlinesDoubleCovered && singleCover
                        Ly2add = pline.vertices.Ly(tidx) ;
                    elseif ~pathlinesDoubleCovered && ~singleCover
                        Ly2add = pline.vertices.Ly(tidx) * 2 ;
                    else
                        error('did not understand what is single/double covered')
                    end

                    % try adding or subtracting Ly to uu to put into correct branch
                    % cut
                    if isnan(fieldfacesU)
                        disp('Looking on different branch cuts for match of this (ux,uy) in pathline mesh')
                        uuNext = uu + [0, Ly2add] ;
                        [fieldfacesU, barycU] = pointLocation(tridx, uuNext) ;   
                        if isnan(fieldfacesU)
                            uuNext = uu - [0, Ly2add] ;
                            [fieldfacesU, barycU] = pointLocation(tridx, uuNext) ; 
                            if isnan(fieldfacesU)
                                error(['Could not find branch cut on which UU is ', ...
                                    'in pathline mesh -- probably UU is not ', ...
                                    'in Lagrangian frame (appears from outside ', ...
                                    'boundary of material coords at t0)'])
                            end
                        end     
                    end
                    if isnan(fieldfacesV)
                        disp('Looking on different branch cuts for match of this (vx,vy) in pathline mesh')
                        vvNext = vv+ [0, Ly2add] ;
                        [fieldfacesV, barycV] = pointLocation(tridx, vvNext) ;   
                        if isnan(fieldfacesV)
                            vvNext = vv - [0, Ly2add] ;
                            [fieldfacesV, barycV] = pointLocation(tridx, vvNext) ; 
                            if isnan(fieldfacesV)
                                error(['Could not find branch cut on which VV is ', ...
                                    'in pathline mesh -- probably VV is not ', ...
                                    'in Lagrangian frame (appears from outside ', ...
                                    'boundary of material coords at t0)'])
                            end
                        end     
                    end
                end

                %% Now what are these coordinates in the Lagrangian frame?
                % Note: if we define mesh as:
                % mesh = QS.getCurrentSPCutMeshSmRS() ;
                %  with normalization such that
                % mesh.u(:, 1) = mesh.u(:, 1) * umax / max(mesh.u(:, 1)) ;
                %  then we can use baryc to recover the position in pullback
                % uu_check = [sum(barycU .* mesh.u(mesh.f(fieldfacesU, :), 1)', 2), ...
                %            sum(barycU .* mesh.u(mesh.f(fieldfacesU, :), 2)', 2)] ;
                % assert(all(uu_check == uu(tidx, :)))
                U_ricci = [sum(barycU .* refMesh.u_ricci(refMesh.f(fieldfacesU, :), 1)', 2), ...
                         sum(barycU .* refMesh.u_ricci(refMesh.f(fieldfacesU, :), 2)', 2)]  ;
                V_ricci = [sum(barycV .* refMesh.u_ricci(refMesh.f(fieldfacesV, :), 1)', 2), ...
                         sum(barycV .* refMesh.u_ricci(refMesh.f(fieldfacesV, :), 2)', 2)]  ;

                
            end

            % Store location in t0 pullback via pathlines
            UallRicci(ii, tidx, :) = U_ricci ;
            VallRicci(ii, tidx, :) = V_ricci ;
            if tidx == 1
                U0Riccis(ii, :) = U_ricci ;
                V0Riccis(ii, :) = V_ricci ;
            end
            ptFaceLocations(ii, tidx, :) = [fieldfacesU, fieldfacesV] ;
            ptBarycenters(ii, tidx, :, :) = [barycU; barycV] ;

            % Check it
            if preview
                figure(3)
                clf
                tridx_ricci = triangulation(refMesh.f, refMesh.u_ricci) ;
                triplot(tridx_ricci, 'color', 'k') ; hold on;
                scatter(U_ricci(:, 1), U_ricci(:, 2), 'filled', 'markerfacecolor', 'r')
                scatter(V_ricci(:, 1), V_ricci(:, 2), 'filled', 'markerfacecolor', 'g')
                pause(pausetime)
            end
        end

        % Convert to 3d locations
        validIdU = find(~isnan(squeeze(UallRicci(ii, :, 1)))) ;
        validIdV = find(~isnan(squeeze(UallRicci(ii, :, 1)))) ;
        [u3d, ffacesCheckU, ~, bcCheckU] = interpolate2Dpts_3Dmesh(...
            refMesh.f, refMesh.u_ricci, refMesh.v, squeeze(UallRicci(ii, validIdU, :))) ;
        [v3d, ffacesCheckV, ~, bcCheckV] = interpolate2Dpts_3Dmesh(...
            refMesh.f, refMesh.u_ricci, refMesh.v, squeeze(VallRicci(ii, validIdV, :))) ;

        % Stuff valid points into full array
        u3dnew = nan(nTimePoints, 3) ;
        v3dnew = nan(nTimePoints, 3) ;
        u3dnew(validIdU, :) = u3d ;
        v3dnew(validIdV, :) = v3d ;
        % keep only the valid (nonNaN) component in store for smoothing
        u3dvalid = u3d ;
        v3dvalid = v3d ;
        u3d = u3dnew ;
        v3d = v3dnew ;

        % inspect the assignments
        % triplot(triangulation(refMesh.f, refMesh.u_ricci)); hold on;
        % plot(squeeze(UallRicci(ii, :, 1)), squeeze(UallRicci(ii, :, 2)), 'o')

        % Check against ptFaceLocations(ii, :, :) and ptBarycenters(ii, :, :, :)
        assert(all(ffacesCheckU == squeeze(ptFaceLocations(ii, validIdU, 1))'))
        assert(all(ffacesCheckV == squeeze(ptFaceLocations(ii, validIdV, 2))'))
        assert(all(all(abs(bcCheckU - squeeze(ptBarycenters(ii, validIdU, 1, :))) < 1e-11)))
        assert(all(all(abs(bcCheckV - squeeze(ptBarycenters(ii, validIdV, 2, :))) < 1e-11)))

        % Acquire the geodesics
        assert(length(validIdU) == length(validIdV)) ;
        nValid = length(validIdU) ;

        disp(['track pair ' num2str(ii) ' > finding geodesics...'])
        glueMesh = glueRectCylinderCutMeshSeam(refMesh) ;
        nvtx = size(glueMesh.v, 1) ;
        uPtIndices = (nvtx + (1:nValid))' ;
        vPtIndices = (nvtx + nValid + (1:nValid))' ;
        pointPairs = [ uPtIndices, vPtIndices] ;
        [geodesicP, pointLocations] = surfaceGeodesicPairs( glueMesh.f, ...
                glueMesh.v, pointPairs, [glueMesh.v; u3d(validIdU, :); v3d(validIdV, :)] ) ;

        % Consistency check: faces are the same as we assigned before
        ffacesCheckU = zeros(length(uPtIndices), 1) ;
        ffacesCheckV = zeros(length(vPtIndices), 1) ;
        for qq = 1:length(uPtIndices)
            ffacesCheckU(qq) = pointLocations(uPtIndices(qq)).face ; 
            ffacesCheckV(qq) = pointLocations(vPtIndices(qq)).face ;
        end    
        assert(all(ffacesCheckU == squeeze(ptFaceLocations(ii, validIdU, 1))'))
        assert(all(ffacesCheckV == squeeze(ptFaceLocations(ii, validIdV, 2))'))

        geodesicPnew = cell(nTimePoints, 1) ;
        for qq = 1:nValid
            geodesicPnew{validIdU(qq)} = geodesicP{qq} ;
        end
        geodesicP = geodesicPnew ;

        % note: ptBary definition here has transpositions. Use
        % barycU and barcV as defined above instead.
        % ptBary = [ pointLocations(nvtx+1).barycentricCoordinates; ...
        %     pointLocations(nvtx+2).barycentricCoordinates ] ;
        % assert(all(barycU == ptBary(1, :)))
        % assert(all(barycV == ptBary(2, :)))

        duGeodesic = zeros(nTimePoints, 1) ;
        for tidx = 1:nTimePoints
            pp = geodesicP{tidx} ;
            duGeodesic(tidx) = sum(vecnorm(diff(pp), 2, 2)) ;
        end

        % Distance traveled in the lagrangian frame
        % smooth u3d a bit
        u3d_sm = u3dvalid ;
        u3d_sm(:, 1) = smoothdata(u3d_sm(:, 1),'sgolay');
        u3d_sm(:, 2) = smoothdata(u3d_sm(:, 2),'sgolay');
        u3d_sm(:, 3) = smoothdata(u3d_sm(:, 3),'sgolay');
        v3d_sm = v3dvalid ;
        v3d_sm(:, 1) = smoothdata(v3d_sm(:, 1),'sgolay');
        v3d_sm(:, 2) = smoothdata(v3d_sm(:, 2),'sgolay');
        v3d_sm(:, 3) = smoothdata(v3d_sm(:, 3),'sgolay');
        Deltau = cumsum(vecnorm(diff(u3d_sm), 2, 2)) ;
        Deltav = cumsum(vecnorm(diff(v3d_sm), 2, 2)) ;
        Deltau = [0; Deltau ] ;
        Deltav = [0; Deltav ] ;

        % Full vectors have nan for non-valid timepoints where track is lost
        u3d_sm_new = nan(nTimePoints, 3) ;
        v3d_sm_new = nan(nTimePoints, 3) ;
        Deltau_new = nan(nTimePoints, 1) ;
        Deltav_new = nan(nTimePoints, 1) ;
        u3d_sm_new(validIdU, :) = u3d_sm ;
        v3d_sm_new(validIdV, :) = v3d_sm ;
        Deltav_new(validIdU) = Deltau ;
        Deltav_new(validIdV) = Deltav ;
        u3d_sm = u3d_sm_new ;
        v3d_sm = v3d_sm_new ;
        Deltau = Deltau_new ;
        Deltav = Deltav_new ;

        % Check smoothing
        if preview
            figure(4)
            clf
            plot3(u3d(:, 1), u3d(:, 2), u3d(:, 3), '.-') ; hold on;
            plot3(v3d(:, 1), v3d(:, 2), v3d(:, 3), '.-') ;
            plot3(u3d_sm(:, 1), u3d_sm(:, 2), u3d_sm(:, 3), '.-') ;
            plot3(v3d_sm(:, 1), v3d_sm(:, 2), v3d_sm(:, 3), '.-') ; 
            axis equal
            xlabel('x'); ylabel('y'); zlabel('z')
            title('distance traveled in Lagrangian frame')
        end

        geodesicPaths{ii} = geodesicP ;
        % ptBarycenters(ii, :, :, :) = ptBarys ; <--- already assigned
        % ptFaceLocations(ii, :, :) = ptFaces ;  <--- already assigned
        lagrangianDistanceTraveledU(ii, :) = Deltau ;
        lagrangianDistanceTraveledV(ii, :) = Deltav ;
        dusGeodesic(ii, :) = duGeodesic ;

        v3d_uLagr(ii, :, :) = u3d ;
        v3d_vLagr(ii, :, :) = v3d ;
        v3d_uLagr_smoothed(ii, :, :) = u3d_sm ;
        v3d_vLagr_smoothed(ii, :, :) = v3d_sm ;
        nSaved = ii ;


        % Get phase from conformal coordinates (Lagrangian frame)
        dY = VallRicci(ii, :, 2) - UallRicci(ii, :, 2) ;
        dX = VallRicci(ii, :, 1) - UallRicci(ii, :, 1) ;

        phaseRaw(ii, :) = atan2(dY, dX) ;
        magXYRaw(ii, :) = vecnorm(VallRicci(ii, :, :) - UallRicci(ii, :, :), 2, 3) ;
        
        dY0 = VallRicci(ii, 1, 2) - UallRicci(ii, 1, 2)  ;
        dX0 = VallRicci(ii, 1, 1) - UallRicci(ii, 1, 1)  ;

        phaseRelative(ii, :) = atan2(dY - dY0, dX - dX0) ;

    end

    disp(['Saving Lagrangian result to : ' relMotionLagrangianFn])
    save(relMotionLagrangianFn, 'dusGeodesic', ...
        'tracks1', 'tracks2', 'pairIDs', ...
        'U0Riccis', 'V0Riccis', ...
        'UallRicci', 'VallRicci', ...
        'magXYRaw', 'phaseRaw', 'phaseRelative', ...
        'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
        'lagrangianDistanceTraveledU', 'lagrangianDistanceTraveledV', ...
        'v3d_uLagr', 'v3d_vLagr', ...
        'v3d_uLagr_smoothed', 'v3d_vLagr_smoothed', 'nSaved') ;
else

    load(relMotionLagrangianFn, 'dusGeodesic', ...
        'tracks1', 'tracks2', 'pairIDs', ...
        'U0Riccis', 'V0Riccis', ...
        'UallRicci', 'VallRicci', ...
        'magXYRaw', 'phaseRaw', 'phaseRelative', ...
        'geodesicPaths', 'ptBarycenters', 'ptFaceLocations', ...
        'lagrangianDistanceTraveledU', 'lagrangianDistanceTraveledV', ...
        'v3d_uLagr', 'v3d_vLagr', ...
        'v3d_uLagr_smoothed', 'v3d_vLagr_smoothed', 'nSaved') ;
    
end

%% Plot results
eulerianResult = load(relMotionFn, 'dusGeodesic') ;

%% Polar/2d scatter plot over time 
% Could consider both raw and relative motion from initial condition. This
% gets a bit tricky for distances but could approximate using isothermal
% metric and fieldfaces
close all

figDir = QS.dir.tracking ;
timestamps = (QS.xp.fileMeta.timePoints - QS.t0set()) * QS.timeInterval ;

if contains(lower(QS.timeUnits), 'min')
    timestamps = timestamps / 60 ;
    timeunits = 'hr' ;
else
    timeunits = QS.timeUnits ;
end

for absolute_relative = 1:2
    
    fig = figure('Position', [100 100 180 180], 'units', 'centimeters');
    if absolute_relative == 1
        exten = '_raw' ;
    else
        exten = '_relative' ;
    end
    
    % consider geodesic distances in eulerian frame (deformed
    % configuration), then Lagrangian. Use phase from Lagrangian
    % in both cases.
    for analysisMode = 1:2

        %% Polar/2d scatter plot over time -- careful method in Lagrangian
        clf
        
        fig2 = figure('Position', [100 100 240 240], 'units', 'centimeters');
        if analysisMode == 1
            duV = eulerianResult.dusGeodesic ;
            phaseV = phaseRaw ;
            distMode = 'rhoEulerian' ;
        else
            duV = dusGeodesic ;
            phaseV = phaseRaw ;
            distMode = 'rhoLagrangian' ;
        end

        % Scatterhist first, then colored paths
        if absolute_relative == 1
            allduX = duV(:) .* cos(phaseV(:)) ;
            allduY = duV(:) .* sin(phaseV(:)) ;
        else
            relduX = duV .* cos(phaseV) - duV(:, 1) .* cos(phaseV(:, 1)) ;
            relduY = duV .* sin(phaseV) - duV(:, 1) .* sin(phaseV(:, 1)) ;
            allduX = relduX(:) ;
            allduY = relduY(:) ;
        end
        timeAll = ones(nTracks, nTimePoints) .* timestamps ;
        timeAll = timeAll(:) ;
        nonnan = ~isnan(allduX) & timeAll < tcut3 ;
        allduX = allduX(nonnan) ;
        allduY = allduY(nonnan) ;
        timeAll = timeAll(nonnan) ;
        
        early = timeAll < tcut1 ;
        mid = timeAll > tcut1 & timeAll < tcut2 ;
        late = timeAll > tcut2 & timeAll < tcut3 ;
        earlyMidLate = zeros(size(timeAll)) ;
        earlyMidLate(early) = 1 ;
        earlyMidLate(mid) = 2 ;
        earlyMidLate(late) = 3 ;
        earlyMidLate(earlyMidLate == 0 ) = NaN ;
        % earlyMidLate = mat2cell(earlyMidLate, length(timeAll), 1) ;
        sh = scatterhist(allduX, allduY, 'group', earlyMidLate, ...
            'kernel', 'on', 'Direction', 'out', ...
            'Location', 'southwest', 'LineStyle',{'-','-','-'}, ...
            'Marker', 'nnnn', 'Legend', 'off');
        hold on;
        alphaVal = 0.5 ;
        for ii = 1:nTracks
            lastIdx = find(isnan(duV(ii, :)) | timestamps > tcut3) ;
            if isempty(lastIdx)
                if absolute_relative == 1
                    xx = abs(duV(ii, :)) .* cos(phaseV(ii, :)) ;
                    yy = abs(duV(ii, :)) .* sin(phaseV(ii, :)) ;
                else
                    xx = relduX(ii, :) ;
                    yy = relduY(ii, :) ;
                end
                aColor = 1:length(find(timestamps < tcut3)) ;
            else
                
                if absolute_relative == 1
                    xx = abs(duV(ii, 1:lastIdx-1)) .* cos(phaseV(ii, 1:lastIdx-1)) ;
                    yy = abs(duV(ii, 1:lastIdx-1)) .* sin(phaseV(ii, 1:lastIdx-1)) ;
                else
                    xx = relduX(ii, 1:lastIdx-1) ;
                    yy = relduY(ii, 1:lastIdx-1) ;
                end
                aColor = 1:lastIdx-1 ;
            end
            mfig = patch([xx NaN], [yy NaN], [aColor length(find(timestamps < tcut3))]);
            % alphamap = 0.2 * ones(length(xx), 1) ;
            % alphamap(end) = 0 ;
            set(mfig,'FaceColor','none','EdgeColor','flat',...
               'LineWidth',1, ...
               'FaceVertexAlphaData',alphaVal,...
               'EdgeAlpha',alphaVal)
        end
        xlabel(['ap displacement [' QS.spaceUnits ']'], ...
            'interpreter', 'latex')
        ylabel(['dv displacement [' QS.spaceUnits ']'], ...
            'interpreter', 'latex')
        
        if analysisMode == 1
            sgtitle('Relative motion in Eulerian frame', 'interpreter', 'latex')
        else
            sgtitle('Relative motion in Lagrangian frame', 'interpreter', 'latex')
        end
        
        colormap(isolum)
        cb = colorbar() ;
        ylims = get(cb, 'ylim') ;
        
        % colorbar label
        if contains(lower(QS.timeUnits), 'min') 
            ylabel(cb, 'time [hr]', 'interpreter', 'latex')
            set(cb, 'Ytick', [1, 16, 31, 46])
            set(cb, 'YtickLabel', ...
                {num2str(QS.timeInterval /60 * (timePoints(1) - min(timePoints))), ...
                num2str(QS.timeInterval /60 * (timePoints(16) - min(timePoints))), ...
                num2str(QS.timeInterval /60 * (timePoints(31) - min(timePoints))), ...
                num2str(QS.timeInterval /60 * (timePoints(46) - min(timePoints)))}) ;
        else
            ylabel(cg, ['time [' QS.timeUnits ']'], 'interpreter', 'latex')
        end
        
        axis equal
        drawnow
        if ~xlimit
            lims = [abs(xlim) abs(ylim)] ;
            xlimit = max(lims) ;
        end

        % Axes are main, above, left
        for axID = 1:length(sh)
            axes(sh(axID))
                xlim([-xlimit, xlimit])
            if axID <2 
                ylim([-xlimit, xlimit])
            end
        end
        s2pos = get(sh(2), 'position') ;
        s3pos = get(sh(3), 'position') ;
        subheight = min([s2pos(3:4) s3pos(3:4)]) * 2 ;
        s2ylim = get(sh(2), 'ylim') ;
        s3ylim = get(sh(3), 'ylim') ;
        set(sh(2), 'position', [s2pos(1), s2pos(2), s2pos(3), subheight])
        set(sh(3), 'position', [s3pos(1)-0.1, s3pos(2), subheight, s3pos(4)])
        set(sh(2), 'ylim', [0, max([s2ylim s3ylim])])
        set(sh(3), 'ylim', [0, max([s2ylim s3ylim])])
        
        saveas(gcf, fullfile(figDir, ['relative_motion_Lagrangian_' distMode '_inPlaneLagrangian' exten '.pdf']))
        close 


        %% Take means for each timepoint -- Lagrangian frame
        clf
        if absolute_relative == 1
            % absolute positions
            meanX = mean(duV .* cos(phaseV), 1, 'omitnan') ;
            meanY = mean(duV .* sin(phaseV), 1, 'omitnan') ;
            stdX = std(duV .* cos(phaseV), 1, 'omitnan') ;
            stdY = std(duV .* sin(phaseV), 1, 'omitnan') ;
        else
            % relative positions
            meanX = mean(relduX, 1, 'omitnan') ;
            meanY = mean(relduY, 1, 'omitnan') ;
            stdX = std(relduX, 1, 'omitnan') ;
            stdY = std(relduY, 1, 'omitnan') ;
        end
        timestamps = (timePoints - QS.t0set()) * QS.timeInterval ;
        
        if contains(lower(QS.timeUnits), 'min')
            timestamps = timestamps / 60 ;
        end
        % subplot(1, 2, 1)
        % plot(timestamps, meanX, '-'); hold on;
        % plot(timestamps, meanY, '-')
        % subplot(1, 2, 2)
        blueColor = [0    0.4470    0.7410] ;
        orangeColor = [0.8500    0.3250    0.0980] ;
        lineProps = {'-','color', orangeColor} ;
        h2 = shadedErrorBar(timestamps, meanY, stdY, 'lineprops', lineProps);
        lineProps = {'-','color', blueColor} ;
        h1 = shadedErrorBar(timestamps, meanX, stdX, 'lineprops', lineProps);
        hold on;
        legend({'$\delta s$', '$\delta \phi$'}, 'interpreter', 'latex')
        ylabel(['relative displacement [' QS.spaceUnits ']'], 'interpreter', 'latex')
        if contains(lower(QS.timeUnits), 'min')
            xlabel(['time [hr]'], 'interpreter', 'latex')
        else
            xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
        end
        
        saveas(gcf, fullfile(figDir, ['relative_motion_Lagrangian_' distMode '_inPlaneLagrangianErr' exten '.pdf'])) ;
        
        
        if analysisMode == 1
            stdX_eulerian = stdX ;
            stdY_eulerian = stdY ;
        else
            stdX_lagrangian = stdX ;
            stdY_lagrangian = stdY ;
        end
        
        % %% Plot errorelipse over time
        % % alphamap(end) = 0 ;
        % set(mfig,'FaceColor','none','EdgeColor','flat',...
        %    'LineWidth',2, ...
        %    'FaceVertexAlphaData',alphaVal,...
        %    'EdgeAlpha',alphaVal)
        % 
        % xlabel(['relative ap displacement [' QS.spaceUnits ']'], ...
        %     'interpreter', 'latex')
        % ylabel(['relative dv displacement [' QS.spaceUnits ']'], ...
        %     'interpreter', 'latex')
        % title('Relative motion in pullback frame')
        % cb = colorbar() ;
        % ylims = get(cb, 'ylim') ;
        % 
        % % colorbar label
        % if contains(lower(QS.timeUnits), 'min') 
        %     ylabel(cb, 'time [hr]', 'interpreter', 'latex')
        %     set(cb, 'Ytick', [1, 16, 31, 46])
        %     set(cb, 'YtickLabel', ...
        %         {num2str(QS.timeInterval /60 * (timePoints(1) - min(timePoints))), ...
        %         num2str(QS.timeInterval /60 * (timePoints(16) - min(timePoints))), ...
        %         num2str(QS.timeInterval /60 * (timePoints(31) - min(timePoints))), ...
        %         num2str(QS.timeInterval /60 * (timePoints(46) - min(timePoints)))}) ;
        % else
        %     ylabel(cg, ['time [' QS.timeUnits ']'], 'interpreter', 'latex')
        % end
        % saveas(gcf, fullfile(figDir, ['relative_motion_Lagrangian_' distMode '_inPlaneLagrangianErr' exten '.pdf'])) ;
    end
    
    %% Compare muscle prediction based on length
    close all
    fig = figure('Position', [100 100 820 180], 'units', 'centimeters');
    tissue_length = zeros(nTimePoints, 1) ;
    for tidx = 1:nTimePoints
        % load spcutMesh, which has a measure of length
        QS.setTime(timePoints(tidx))
        mesh = QS.getCurrentSPCutMeshSmRS() ;
        maxUvals = mesh.u(QS.nU:QS.nU:end, 1) ;
        assert(all(maxUvals == maxUvals(1)))
        tissue_length(tidx) = maxUvals(1) ;
    end
    % X variations, raw
    h1 = subplot(1, 5, 1) ;
    scatter( stdX_lagrangian, stdX_eulerian, 20, timestamps, 'filled' )
    xlabel({'transverse relative motion'; ...
        'in material coordinates, $\sigma^L_\perp$'}, 'interpreter', 'latex')
    ylabel({'parallel relative'; 'motion $\sigma^E_{\perp}$'}, ...
        'interpreter', 'latex')
    
    % Title with corr coeff and linear fit
    % correlation coefficient
    RR = corrcoef(stdX_lagrangian, stdX_eulerian) ;
    rho = RR(1, 2) ;
    % fit to line
    pp = polyfit(stdX_lagrangian, stdX_eulerian, 1) ;
    corrString = ['$\rho=$' sprintf('%0.3f', rho) ];
    fitString = [' $\sigma_{\perp}^{E}=$' sprintf('%0.3f', pp(1)) ...
        '$\sigma_{\perp}^{L} +$' sprintf('%0.3f', pp(2)) ] ;
    % Title
    titleString = [corrString fitString] ;
    th1 = title(titleString, 'interpreter', 'latex') ;    
    t1pos = get(th1, 'Position') ;
    t1pos(2) = t1pos(2) ;
    set(th1, 'position', t1pos)
    
    % format h1
    axis equal
    drawnow
    h1pos = get(h1, 'position') ;
    
    % X variations, rescaled
    h2 = subplot(1, 5, 3) ;
    foldIncreaseLength = tissue_length / tissue_length(QS.xp.tIdx(QS.t0set())) ;
    stdX_LagrNorm = stdX_lagrangian'.* foldIncreaseLength ;
    scatter( stdX_LagrNorm, stdX_eulerian, 20, timestamps, 'filled' )
    hold on;
    identityCurv = [min(max(stdX_LagrNorm), min(stdX_eulerian)), ...
                 max(max(stdX_LagrNorm), max(stdX_eulerian))] ;
    plot(identityCurv, identityCurv, 'k--')
    ylabel({'transverse relative'; ...
        ['motion, $\sigma^E_\perp$ [' QS.spaceUnits ']']}, 'interpreter', 'latex')
    xlabel({'rescaled transverse relative motion'; ['in material coordinates $\sigma_{\perp}^{L} L_s / L_s^0$ [' QS.spaceUnits ']']}, ...
        'interpreter', 'latex')
    
    % Title with corr coeff and linear fit
    % correlation coefficient
    RR = corrcoef(stdX_LagrNorm, stdX_eulerian) ;
    rho = RR(1, 2) ;
    % fit to line
    pp = polyfit(stdX_LagrNorm, stdX_eulerian, 1) ;
    corrString = ['$\rho=$' sprintf('%0.3f', rho) ];
    fitString = [' $\sigma_{\perp}^{E}=$' sprintf('%0.3f', pp(1)) ...
        '$(\sigma_{\perp}^{L} L_s / L_s^0) +$' sprintf('%0.3f', pp(2)) ] ;
    % Title
    titleString = [corrString fitString] ;
    th2 = title(titleString, 'interpreter', 'latex') ;
    t2pos = get(th2, 'Position') ;
    t2pos(2) = t2pos(2) + 8 ;
    set(th2, 'position', t2pos)
    
    % format h2
    axis equal
    drawnow
    h2pos = get(h2, 'position') ;
    set(h2, 'position', [h2pos(1) h2pos(2) h1pos(3) h1pos(4)])
    Xlim2 = xlim;
    Ylim2 = ylim ;
    
    
    % Y variations
    h3 = subplot(1, 5, 5) ;
    scatter(stdY_lagrangian, stdY_eulerian, 20, timestamps, 'filled');
    hold on;
    identityCurv = [min(min(stdY_lagrangian), min(stdY_eulerian)), ...
        max(max(stdY_lagrangian), max(stdY_eulerian))] ;
    plot(identityCurv, identityCurv, 'k--')
    ylabel({'parallel relative'; ...
        ['motion, $\sigma^E_\parallel$ [' QS.spaceUnits ']']}, 'interpreter', 'latex')
    xlabel({'parallel relative motion'; ['in material coordinates $\sigma_{\parallel}^{L}$ [' QS.spaceUnits ']']}, ...
        'interpreter', 'latex')
    
    % Title with corr coeff and linear fit
    % correlation coefficient
    RR = corrcoef(stdY_lagrangian, stdY_eulerian) ;
    rho = RR(1, 2) ;
    % fit to line
    pp = polyfit(stdY_lagrangian, stdY_eulerian, 1) ;
    corrString = ['$\rho=$' sprintf('%0.3f', rho) ];
    fitString = [' $\sigma_{\parallel}^{E} =$' sprintf('%0.3f', pp(1)) ...
        '$\sigma_{\parallel}^{L} +$' sprintf('%0.3f', pp(2)) ] ;
    % Title
    titleString = [corrString fitString] ;
    th3 = title(h3, titleString, 'interpreter', 'latex') ;
    t3pos = get(th3, 'Position') ;
    t3pos(2) = t3pos(2) + 5 ;
    set(th3, 'position', t3pos)
    
    % Colorbar and formatting
    axis equal
    colorbar
    colormap(isolum)
    drawnow
    h3pos = get(h3, 'position') ;
    set(h3, 'position', [h3pos(1) h3pos(2) h1pos(3) h1pos(4)])
    Xlim3 = xlim;
    Ylim3 = ylim ;
    
    lims = [Xlim2, Xlim3, Ylim2, Ylim3];
    commonLim = [min(lims), max(lims)] ;
    %set(h2, 'xlim', commonXlim)
    %set(h3, 'xlim', commonXlim)
    xlim(h2, commonLim)
    ylim(h2, commonLim)
    xlim(h3, commonLim)
    ylim(h3, commonLim)
    
    
    axes(h1)
    Xlim1 = xlim ;
    Ylim1 = ylim ;
    lims1 = [Xlim1, Ylim1] ;
    xlim(h1, [min(lims1), max(lims1)])
    ylim(h1, [min(lims1), max(lims1)])
    
    h3pos = get(h3, 'position') ;
    set(h1, 'position', [h1pos(1) h1pos(2) h3pos(3) h3pos(4)])
    
     
    
    % SAVE THE FIGURE
    saveas(gcf, fullfile(figDir, ['relative_motion_Lagrangian_' distMode '_anisotropicPrediction' exten '.pdf'])) ;
    save(fullfile(figDir, ['relative_motion_Lagrangian_' distMode '_anisotropicPrediction' exten '.mat']), ...
        'stdX_LagrNorm', 'stdX_eulerian', 'stdX_lagrangian', 'tissue_length', 'timestamps', ...
        'stdY_eulerian', 'stdY_lagrangian') ;
    
    
    % Organ length contribution
    close all
    fig = figure('Position', [100 100 300 300], 'units', 'centimeters');
    tissue_length = zeros(nTimePoints, 1) ;
    for tidx = 1:nTimePoints
        % load spcutMesh, which has a measure of length
        QS.setTime(timePoints(tidx))
        mesh = QS.getCurrentSPCutMeshSmRS() ;
        maxUvals = mesh.u(QS.nU:QS.nU:end, 1) ;
        assert(all(maxUvals == maxUvals(1)))
        tissue_length(tidx) = maxUvals(1) ;
    end
    plot(timestamps, tissue_length)
    xlabel('time [hr]', 'interpreter', 'latex')
    ylabel(['organ length [' QS.spaceUnits ']'], 'interpreter', 'latex')

    saveas(gcf, fullfile(figDir, ['organLength.pdf']))
    
end
