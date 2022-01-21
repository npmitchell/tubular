function generateCellSegmentation3D(QS, options)
% generateCellSegmentation3D(QS, options)
% In the code below, "aspect" typically means (a/b-1) = ars - 1
%
% Define Q tensor for each cell by
%   Q =  q (n^T n - II/2), where q = (sqrt(I_1/I_2) - 1) is the
%   magnitude of the anisotropy
% 
%
% The Q tensor is related to the aspect ratio of cells via:
%   abs(norm(Q)) * 2 = |ar| - 1
% 
% The mratio 
% 
% The factor of two arises on the LHS since the 
% norm(n'*n - Identity*0.5) = 0.5 for any unit vector n, and an isotropic
% cell will have norm(Q) = 0.  
%
% NPMitchell 2021

% Colors for cosine, sine:
tmp = define_colors ;
% CScolors = [0.90    0.55    0.55 ;
%         0.62    0.76    0.84 ] ;
CScolors = tmp([9,5], :) ;


% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;
maxCellSize = Inf ;  % maximum size allowed to include a cell
overwrite = false ;
corrected = false ;
coordSys = 'spsme' ;
overwriteImages = false ;
debug = false ;
[~, ~, ~, xyzlims] = QS.getXYZLims() ;
xyzlims = xyzlims + [-10, 10 ] ;
t0 = QS.t0set() ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
    coordSys = lower(erase(coordSys, '_')) ;
    % If the coordinate system is (s,phi) relaxed aspect ratio smoothed
    % extended, correct for variants of naming convention.
    if strcmpi(coordSys, 'rsme') || strcmpi(coordSys, 'rspsme')
        coordSys = 'sprsme' ;
    end
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'xyzlims')
    xyzlims = options.xyzlims ;
end
if isfield(options, 'correctedSegmentation')
    corrected = options.correctedSegmentation ;
end

% Directory preparation
if corrected
    segDir = fullfile(QS.dir.segmentation, 'seg3d_corrected') ;
    if ~exist(segDir, 'dir')
        mkdir(segDir)
    end
    imdir = fullfile(QS.dir.segmentation, 'seg3d_corrected', 'images') ;
else
    imdir = fullfile(QS.dir.segmentation, 'seg3d', 'images') ;
end
if ~exist(imdir, 'dir')
    mkdir(imdir)
end

% Prep for dicing data by lobes
features = QS.getFeatures() ;
folds = features.folds ;
nLobes = size(folds, 2) + 1 ;

% Setting the current timepoint clears non-timepoint segmentations
close all
aratio_median = [] ;
aratio_mean = [] ;
aratio_low25 = [] ;
aratio_high75 = [] ;
aratio_std = [] ;
aratio_ste = [] ;

theta_mean = [] ;
theta_std = [] ;
theta_ste = [] ;
c2t_low25 = [] ;
c2t_high75 = [] ;
s2t_low25 = [] ;
s2t_high75 = [] ;
c2t_mean = [] ;
s2t_mean = [] ;
c2t_std = [] ;
s2t_std = [] ;
c2t_ste = [] ;
s2t_ste = [] ;
% aspect-weighted
mQAar       = [] ;
mQAarStd    = [] ;
mQAarSte    = [] ;
mQAtheta    = [] ;
mQAthetaStd = [] ;
mQAthetaSte = [] ;
mQAc2t      = [] ;
mQAs2t      = [] ;
mQAc2tStd   = [] ;
mQAs2tStd   = [] ;
mQAc2tSte   = [] ;
mQAs2tSte   = [] ;
% eccentricity-weighted
mQEar       = [] ;
mQEarStd    = [] ;
mQEarSte    = [] ;
mQEtheta    = [] ;
mQEthetaStd = [] ;
mQEthetaSte = [] ;
mQEc2t      = [] ;
mQEs2t      = [] ;
mQEc2tStd   = [] ;
mQEs2tStd   = [] ;
mQEc2tSte   = [] ;
mQEs2tSte   = [] ;

    
dmy = 1 ;
cos2thetaM = zeros(99, 1) ;
sin2thetaM = cos2thetaM ;
aspectM = cos2thetaM ;
thetaM = cos2thetaM ;
nAPBins = 20 ;

% (aspect-1)-weighted Q tensor stats
mean_QAc2ts = zeros(length(timePoints), nAPBins) ;
mean_QAs2ts = zeros(length(timePoints), nAPBins) ;
std_QAc2ts = zeros(length(timePoints), nAPBins) ;
std_QAs2ts = zeros(length(timePoints), nAPBins) ;
ste_QAc2ts = zeros(length(timePoints), nAPBins) ;
ste_QAs2ts = zeros(length(timePoints), nAPBins) ;
meanQALobeAspects = zeros(nLobes, length(timePoints)) ;
meanQALobeAspectStds = zeros(nLobes, length(timePoints)) ;
meanQALobeAspectStes = zeros(nLobes, length(timePoints)) ;
meanQALobeThetas = zeros(nLobes, length(timePoints)) ;
meanQALobeThetaStds = zeros(nLobes, length(timePoints)) ;
meanQALobeThetaStes = zeros(nLobes, length(timePoints)) ;

% eccentricity-weighted Q tensor stats
mean_QEc2ts = zeros(length(timePoints), nAPBins) ;
mean_QEs2ts = zeros(length(timePoints), nAPBins) ;
std_QEc2ts = zeros(length(timePoints), nAPBins) ;
std_QEs2ts = zeros(length(timePoints), nAPBins) ;
ste_QEc2ts = zeros(length(timePoints), nAPBins) ;
ste_QEs2ts = zeros(length(timePoints), nAPBins) ;
meanQELobeAspects = zeros(nLobes, length(timePoints)) ;
meanQELobeAspectStds = zeros(nLobes, length(timePoints)) ;
meanQELobeAspectStes = zeros(nLobes, length(timePoints)) ;
meanQELobeThetas = zeros(nLobes, length(timePoints)) ;
meanQELobeThetaStds = zeros(nLobes, length(timePoints)) ;
meanQELobeThetaStes = zeros(nLobes, length(timePoints)) ;

for tp = timePoints
    sprintf(['t = ' num2str(tp)])
    QS.setTime(tp)
    tidx = QS.xp.tIdx(tp) ;
    
    if corrected
        outfn = sprintf(QS.fullFileBase.segmentation3dCorrected, QS.currentTime) ;
    else
        outfn = sprintf(QS.fullFileBase.segmentation3d, QS.currentTime) ;
    end
    if ~exist(outfn, 'file') || overwrite

        % Obtain the segmentation in 2d 
        if corrected
            seg2d = QS.getCurrentSegmentation2DCorrected(options) ;
        else
            seg2d = QS.getCurrentSegmentation2D(options) ;
        end
        
        % obtain the current cut mesh in APDV coordinates -->
        % 2d coordinates in plane of MESH
        if strcmp(coordSys, 'spsme') || strcmp(coordSys, 'sprsme') 
            cutMesh = QS.getCurrentSPCutMeshSmRS() ;

            % Normalize the u coordinate
            cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
        else
            error(['Code for this coordSys: ' coordSys])
        end

        % tile annular cut mesh
        tileCount = [1,  1] ;
        [ faces, v2D, v3D ] = tileAnnularCutMesh(cutMesh, tileCount) ;

        % Collate cell vertices as XY
        % XY = zeros(nVertices, 2) ;
        % for qq = 1:nVertices
        %     XY(qq, :) = [seg2d.seg2d.Vdat(qq).vertxcoord, seg2d.seg2d.Vdat(qq).vertycoord] ;
        % end
        if corrected
            % corrected segmentation stores polygons as vertices
            polygons = seg2d.seg2d.cdat.polygons ;
            XY = [] ;
            pgonIDs = [] ;
            for pp = 1:length(polygons)
                pgon = polygons{pp} ;
                if ~isempty(pgon)
                    pgon = pgon(:, [2, 1]) ;
                    if isempty(XY)
                        XY = pgon ;
                    else
                        XY = [XY; pgon] ;
                    end
                    nn = length(pgonIDs) ;
                    pgonIDs(nn+1:nn+size(pgon, 1)) = pp ;
                end
            end
        else
            % seg2d stores a network of vertices and bonds
            XY = seg2d.seg2d.vdat.v ;
        end
        
        % Collate cell centroids
        if corrected
            centroids = seg2d.seg2d.cdat.centroid ;
            nCells = size(centroids, 1) ;
        else
            nCells = length(seg2d.seg2d.Cdat) ;
            centroids = zeros(nCells, 2) ;
            for qq = 1:nCells
                centroids(qq, :) = seg2d.seg2d.Cdat(qq).centroid.coord ;
            end    
        end
        
        % Check that we are interpolating over the same domain of
        % parameterization
        if strcmpi(seg2d.coordSys, 'spsme') || ...
            strcmpi(seg2d.coordSys, 'sprsme')
            Ly = size(seg2d.segIm, 1) ;
            Lx = size(seg2d.segIm, 2) ;
            doubleCovered = true ;
        else
            error(['Did not code for coordSys ' coordSys 'yet'])
        end
        % uv coordinates of the cell vertices
        umax = max(cutMesh.u(:, 1)) ;
        vmax = max(cutMesh.u(:, 2)) ;
        uv = QS.XY2uv([Lx, Ly], XY, doubleCovered, umax, vmax) ;

        % Centroids in uv coordinates
        cntrds = QS.XY2uv([Lx, Ly], centroids, doubleCovered, umax, vmax) ;

        % cell vertices in 3d
        [c3d, vertexMeshFaces] = ...
            interpolate2Dpts_3Dmesh(faces, v2D, v3D, uv) ;

        if corrected
            poly3d = cell(length(polygons), 1) ;
            polygonVtxIds = cell(length(polygons), 1) ;
            cellVtx2D = cell(length(polygons), 1) ;
            for pp = 1:length(polygons)
                poly3d{pp} = c3d(pgonIDs == pp, :) ;
                polygonVtxIds{pp} = find(pgonIDs == pp) ;
                cellVtx2D{pp} = uv(polygonVtxIds{pp}, :) ;
            end

            % Obtain cell vertices in 3d
            % cell2d0 = seg2d.seg2d.cdat.polygons{cid} ;
            % cellVtx2D = uv(pgonIDs == cid, :) ;
            % cellVtx0 = c3d(pgonIDs == cid, :) ;
            [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
                moment2, ang2, moinertia, cellQ2d, ...
                cellMeshFaces, vertexMeshFaces, cellNormals ] = ...
                polygonNetwork3dMeasurements(faces, v3D, v2D, uv, polygonVtxIds, cntrds) ;
        else
            % Obtain cell vertices in 3d
            % original segmentation stores indices into vertices for
            % polygon shapes
            % cell2d0 = seg2d.seg2d.vdat.v(seg2d.seg2d.cdat.polygons{cid}, :) ;
            polygonVtxIds = seg2d.seg2d.cdat.polygons ;
            for pp = 1:length(polygonVtxIds)
                cellVtx2D{pp} = uv(polygonVtxIds{cid}, :) ;
                % cellVtx0 = c3d(seg2d.seg2d.cdat.polygons{cid}, :) ;
            end
            [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
                moment2, ang2, moinertia, cellQ2d,  cellMeshFaces, vertexMeshFaces] = ...
                polygonNetwork3dMeasurements(faces, v3D, v2D, uv, polygonVtxIds, cntrds) ;
        end
        
        %% Save results stored in struct
        if corrected
            seg2d = QS.currentSegmentation.seg2dCorrected ;
        else
            seg2d = QS.currentSegmentation.seg2d ;
        end
        seg3d = struct('vdat', struct(), 'cdat', struct(), ...
            'qualities', struct(), 'map', struct()) ;
        
        % Mesh information for mapping
        seg3d.map.f = faces ;
        seg3d.map.v2d = v2D ;
        seg3d.map.v3d = v3D ;
        
        % vertex data in 3d
        seg3d.vdat.xyzrs = c3d ;
        seg3d.vdat.uv = uv ;
        seg3d.vdat.meshFaces = vertexMeshFaces ;
        if ~corrected
            seg3d.vdat.NL = seg2d.seg2d.vdat.NL ;
            seg3d.vdat.BL = seg2d.seg2d.vdat.BL ;
            seg3d.vdat.fourfold = seg2d.seg2d.vdat.fourfold ;
        end
        
        % cell data in 3d
        seg3d.cdat.centroids_uv = cntrds ;
        seg3d.cdat.centroids_3d = cellCntrd3d ;
        seg3d.cdat.meshFaces = cellMeshFaces ;
        seg3d.cdat.normals = cellNormals ;
        if corrected
            seg3d.cdat.polygons = poly3d ;
            seg3d.vdat.vertexPolygonMemberID = pgonIDs ;
            seg3d.cdat.polygonVertexID = cell(length(polygons), 1) ;
            dmyk = 1 ;
            for pp = 1:length(polygons)
                if length(polygons{pp}) > 1
                    seg3d.cdat.polygonVertexID{pp} = dmyk:(dmyk+length(polygons{pp})-1) ;
                    dmyk = dmyk + length(polygons{pp}) ;
                end
            end
        else
            seg3d.cdat.polygons = seg2d.seg2d.cdat.polygons ;
        end
        
        % Cell nematic tensor Q and strength of elongation qstrength
        % --> add to struct after filtering
        % mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
        % strength = zeros(nCells, 1) ;
        % QQ = zeros(nCells, 2, 2) ;
        % for qq = 1:nCells
        %     if ~isempty(intersect(keep, qq))
        %         tt = mod(seg3d.qualities.ang1(qq), pi) ;
        %         nn = [cos(tt), sin(tt)] ;
        %         % Create traceless symmetric matrix using unit vec
        %         strength(qq) = abs(sqrt(mratio(qq))) - 1 ;
        %         QQ(qq, :, :) = nn' * nn - 0.5 * [1, 0; 0, 1] ;
        %     end
        % end
        
        % cell qualities, including Q= q*(n^T n - II/2)
        seg3d.qualities = struct() ;
        seg3d.qualities.areas = areas ;
        seg3d.qualities.perim = perim ;
        seg3d.qualities.moment1 = moment1 ;
        seg3d.qualities.moment2 = moment2 ;
        seg3d.qualities.ang1 = ang1 ;
        seg3d.qualities.ang2 = ang2 ;
        seg3d.qualities.mInertia = moinertia ;
        seg3d.qualities.cellQ2d = cellQ2d ;
        % seg3d.qualities.nematicTensor = QQ ;   % added after filtering
        % seg3d.qualities.nematicStrength = strength ; % added after filtering
        seg3d.qualities.readme = struct(...
            'areas', ['area of each cell in squared ' QS.spaceUnits], ...
            'perim', ['perimeter of each cell in ' QS.spaceUnits], ...
            'moment1', ['smaller moment of inertia (about long axis) ', ...
                'in density * squared ' QS.spaceUnits], ...
            'moment2', ['larger moment of inertia (about short axis) ', ...
                'in density * squared ' QS.spaceUnits], ...
            'mInertia', ['centroidal moment of inertia in coordSys ', ...
                'coords Iuu Iuv Ivv, with uv parallel with XY but ', ...
                'offset to centroid, in density * squared ', ...
                QS.spaceUnits], ...
            'ang1', 'angle of long axis in coordSys, in radians', ...
            'ang2', 'angle of short axis in coordSys, in radians', ...
            'cellQ2d', ['quasi-2d cell polygon from embedding space, ', ...
                    'but with cell centroid surface normal rotated ', ...
                    'to be along x axis'], ...
            'nematicTensor', 'n^T * n - 0.5 * [1, 0; 0, 1], where n is along long axis', ...
            'nematicStrength', 'abs(sqrt(MOIEigenvalueRatio)) - 1, strength of elongation' ) ;
            
        % cell statistics 
        % find which are "good" cells to consider
        keep = find(~isnan(ang1) & (areas < maxCellSize) & ...
            moment1 > 0 & moment2 > 0) ;
        seg3d.statistics.keep = keep ;
        seg3d.statistics.maxCellSize = maxCellSize ;
         
        %% SAVE BEFORE STATS
        % which coordinate system has been used for segmentation
        coordSys = seg2d.coordSys ;
        save(outfn, 'seg3d', 'coordSys')
        
        
        
        %% -----------------------------------------------------------------
        % NOW DO STATS
        %%-----------------------------------------------------------------
        
        %% Medians of orientation and moment ratio over TIME    
        % Compute cell statistics
        nCells = length(seg3d.qualities.areas) ;
        keep = seg3d.statistics.keep ;
        areas = seg3d.qualities.areas ;

        % New code here
        % ----------------------------------------
        % % Compute mean MOI by first making the tensor unit 
        iuu = nanmean(seg3d.qualities.mInertia(keep, 1, 1)) ;
        iuv = - nanmean(seg3d.qualities.mInertia(keep, 1, 2)) ;
        ivv = nanmean(seg3d.qualities.mInertia(keep, 2, 2)) ;
        iii = [ iuu  -iuv ;
             -iuv   ivv ];
        [ eig_vec, eig_val ] = eig(iii);
        meanMOI = iii ;
        meanMOIMoment1 = eig_val(1,1);
        meanMOIMoment2 = eig_val(2,2);
        meanMOIAng1 = atan2( eig_vec(2,1), eig_vec(1,1) );

        % Compute mean "MOI" by first computing Q tensor then area-weight
        % --> is MOI already area-weighted? No, it's weighted oddly. 

        % Could rescale MOIs by sqrt(determinant) so each is nematic tensor 
        %   with unit 'size'. How to do this properly? Define Q tensor for each
        %   by Q =  q (n^T n - II/2), where q = (sqrt(I_1/I_2) - 1) is the
        %   magnitude of the anisotropy

        % Compute nematic tensor for each
        mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
        mDM2 = seg3d.qualities.moment1 ./ seg3d.qualities.moment2 ;
        strength = zeros(nCells, 1) ;
        eccentricity = zeros(nCells, 1) ;
        QQ = zeros(nCells, 2, 2) ;
        for qq = 1:nCells
            if ~isempty(intersect(keep, qq))
                tt = mod(seg3d.qualities.ang1(qq), pi) ;
                nn = [cos(tt), sin(tt)] ;
                % Create traceless symmetric matrix using unit vec
                strength(qq) = abs(sqrt(mratio(qq))) - 1 ;
                eccentricity(qq) = sqrt(1 - mDM2(qq)) ;
                QQ(qq, :, :) = nn' * nn - 0.5 * [1, 0; 0, 1] ;
            end
        end

        % Statistics on Q tensors
        nU = QS.nU ;
        fold0 = double(folds(tidx, :)) / double(nU) ;
        nLobes = length(fold0(:)) + 1 ;
        foldt = [0; fold0(:); 1] ;
        ap_pos = seg3d.cdat.centroids_uv(keep, 1) ;
        xedges = linspace(0, 1, nAPBins + 1) ;
        % QtensorStats compute here
        outstats = aux_computeQstats(QQ, ang1, mratio, ...
            keep, areas, strength, eccentricity, ...
            foldt, nLobes, ap_pos, xedges, maxCellSize) ;

        % Dump stats into seg3d
        % ensure that no fields are missing from outstats not already in
        % seg3d.statistics as field
        seg3d.statistics = outstats ;
        seg3d.statistics.meanMoI = meanMOI ;
        seg3d.statistics.meanMoImoment1 = meanMOIMoment1 ;
        seg3d.statistics.meanMoImoment2 = meanMOIMoment2 ;
        seg3d.statistics.meanMoITheta = meanMOIAng1 ;

        %% SAVE seg3d with stats
        save(outfn, 'seg3d', 'coordSys')

    else
        if corrected
            seg3d = QS.loadCurrentSegmentation3DCorrected() ;
        else
            seg3d = QS.loadCurrentSegmentation3D() ;
        end
        coordSys = seg3d.coordSys ;
        seg3d = seg3d.seg3d ;

        keep = seg3d.statistics.keep ;
        mratio = seg3d.qualities.moment2 ./ seg3d.qualities.moment1 ;
    end
    
    
    %% Store for later plotting
    % Collate for later plotting -- raw stats and meanQs (global)
    
    % shorthands for some vars
    % meanQAWB = seg3d.statistics.meanQ.aspectWeighted.meanQWeightBounded ;
    meanQAAspectWB = seg3d.statistics.meanQ.aspectWeighted.meanQAspectWeightBounded ;
    meanQAAspectStdWB = seg3d.statistics.meanQ.aspectWeighted.meanQAspectStdWeightBounded ;
    meanQAAspectSteWB = seg3d.statistics.meanQ.aspectWeighted.meanQAspectSteWeightBounded ;
    meanQAThetaWB = seg3d.statistics.meanQ.aspectWeighted.meanQThetaWeightBounded ;
    meanQAThetaStdWB = seg3d.statistics.meanQ.aspectWeighted.meanQThetaStdWeightBounded ;
    meanQAThetaSteWB = seg3d.statistics.meanQ.aspectWeighted.meanQThetaSteWeightBounded ;
        
    % meanQEWB = seg3d.statistics.meanQ.eccentricityWeighted.meanQWeightBounded ;
    meanQEAspectWB = seg3d.statistics.meanQ.eccentricityWeighted.meanQAspectWeightBounded ;
    meanQEAspectStdWB = seg3d.statistics.meanQ.eccentricityWeighted.meanQAspectStdWeightBounded ;
    meanQEAspectSteWB = seg3d.statistics.meanQ.eccentricityWeighted.meanQAspectSteWeightBounded ;
    meanQEThetaWB = seg3d.statistics.meanQ.eccentricityWeighted.meanQThetaWeightBounded ;
    meanQEThetaStdWB = seg3d.statistics.meanQ.eccentricityWeighted.meanQThetaStdWeightBounded ;
    meanQEThetaSteWB = seg3d.statistics.meanQ.eccentricityWeighted.meanQThetaSteWeightBounded ;

    mratio_principal = mratio(keep) ;
    ars = sqrt(mratio_principal(:)) ;
    
    % Multiplying by two here since norm gives 1/2*(strength_Q)
    aratio_median = [aratio_median, median(ars)] ;
    aratio_mean = [aratio_mean, mean(ars) ] ;
    aratio_low25 = [aratio_low25, seg3d.statistics.aspect25 ] ;
    aratio_high75 = [aratio_high75, seg3d.statistics.aspect75 ] ;
    aratio_std = [aratio_std, seg3d.statistics.aspectStd ] ;
    aratio_ste = [aratio_ste, seg3d.statistics.aspectSte ] ;

    theta_mean = [theta_mean, seg3d.statistics.thetaMean ] ;
    theta_std = [theta_std, seg3d.statistics.thetaStd ] ;
    theta_ste = [theta_ste, seg3d.statistics.thetaSte ] ;
    c2t_low25 = [c2t_low25, seg3d.statistics.cos2theta25] ;
    c2t_high75 = [c2t_high75, seg3d.statistics.cos2theta75 ] ;
    s2t_low25 = [s2t_low25, seg3d.statistics.sin2theta25 ] ;
    s2t_high75 = [s2t_high75, seg3d.statistics.sin2theta75 ] ;
    c2t_mean = [c2t_mean, seg3d.statistics.cos2thetaMean ] ;
    s2t_mean = [s2t_mean, seg3d.statistics.sin2thetaMean ] ;
    c2t_std = [c2t_std, seg3d.statistics.cos2thetaStd ] ;
    s2t_std = [s2t_std, seg3d.statistics.sin2thetaStd ] ;
    c2t_ste = [c2t_ste, seg3d.statistics.cos2thetaSte ] ;
    s2t_ste = [s2t_ste, seg3d.statistics.sin2thetaSte ] ;
    % aspect-weighted
    mQAar       = [mQAar, meanQAAspectWB] ;
    mQAarStd    = [mQAarStd, meanQAAspectStdWB] ;
    mQAarSte    = [mQAarSte, meanQAAspectSteWB] ;
    mQAtheta    = [mQAtheta, meanQAThetaWB] ;
    mQAthetaStd = [mQAthetaStd, meanQAThetaStdWB] ;
    mQAthetaSte = [mQAthetaSte, meanQAThetaSteWB] ;
    mQAc2t      = [mQAc2t, cos(2*meanQAThetaWB)] ;
    mQAs2t      = [mQAs2t, sin(2*meanQAThetaWB)] ;
    mQAc2tStd   = [mQAc2tStd, abs(sin(2*meanQAThetaWB))*2*meanQAThetaStdWB ] ;
    mQAs2tStd   = [mQAs2tStd, abs(cos(2*meanQAThetaWB))*2*meanQAThetaStdWB ] ;
    mQAc2tSte   = [mQAc2tSte, abs(sin(2*meanQAThetaWB))*2*meanQAThetaSteWB ] ;
    mQAs2tSte   = [mQAs2tSte, abs(cos(2*meanQAThetaWB))*2*meanQAThetaSteWB ] ;
    % eccentricity-weighted
    mQEar       = [mQEar, meanQEAspectWB] ;
    mQEarStd    = [mQEarStd, meanQEAspectStdWB] ;
    mQEarSte    = [mQEarSte, meanQEAspectSteWB] ;
    mQEtheta    = [mQEtheta, meanQEThetaWB] ;
    mQEthetaStd = [mQEthetaStd, meanQEThetaStdWB] ;
    mQEthetaSte = [mQEthetaSte, meanQEThetaSteWB] ;
    mQEc2t      = [mQEc2t, cos(2*meanQEThetaWB)] ;
    mQEs2t      = [mQEs2t, sin(2*meanQEThetaWB)] ;
    mQEc2tStd   = [mQEc2tStd, abs(sin(2*meanQEThetaWB))*2*meanQEThetaStdWB ] ;
    mQEs2tStd   = [mQEs2tStd, abs(cos(2*meanQEThetaWB))*2*meanQEThetaStdWB ] ;
    mQEc2tSte   = [mQEc2tSte, abs(sin(2*meanQEThetaWB))*2*meanQEThetaSteWB ] ;
    mQEs2tSte   = [mQEs2tSte, abs(cos(2*meanQEThetaWB))*2*meanQEThetaSteWB ] ;
   
    %% Collate for plotting -- lobes
    % anisotropy-weighted (aspect - 1) -- LOBES 
    meanQALobeAspects(:, dmy) = ...
        seg3d.statistics.lobes.aspectWeighted.meanQLobeAspect ;
    meanQALobeAspectStds(:, dmy) = ...
        seg3d.statistics.lobes.aspectWeighted.meanQLobeAspectStd ;
    meanQALobeAspectStes(:, dmy) = ...
        seg3d.statistics.lobes.aspectWeighted.meanQLobeAspectSte ;
    meanQALobeThetas(:, dmy) = ...
        seg3d.statistics.lobes.aspectWeighted.meanQLobeTheta ;
    meanQALobeThetaStds(:, dmy) = ...
        seg3d.statistics.lobes.aspectWeighted.meanQLobeThetaStd ;
    meanQALobeThetaStes(:, dmy) = ...
        seg3d.statistics.lobes.aspectWeighted.meanQLobeThetaSte ;
    % eccentricity-weighted (aspect - 1)    
    meanQELobeAspects(:, dmy) = ...
        seg3d.statistics.lobes.eccentricityWeighted.meanQLobeAspect ;
    meanQELobeAspectStds(:, dmy) = ...
        seg3d.statistics.lobes.eccentricityWeighted.meanQLobeAspectStd ;
    meanQELobeAspectStes(:, dmy) = ...
        seg3d.statistics.lobes.eccentricityWeighted.meanQLobeAspectSte ;
    meanQELobeThetas(:, dmy) = ...
        seg3d.statistics.lobes.eccentricityWeighted.meanQLobeTheta ;
    meanQELobeThetaStds(:, dmy) = ...
        seg3d.statistics.lobes.eccentricityWeighted.meanQLobeThetaStd ;
    meanQELobeThetaStes(:, dmy) = ...
        seg3d.statistics.lobes.eccentricityWeighted.meanQLobeThetaSte ;

    %% collate for plotting -- AP stats
    % aspect-weighted
    mean_QAc2ts(dmy, :) = ...
        seg3d.statistics.apStats.aspectWeighted.apCos2Theta ;
    mean_QAs2ts(dmy, :) = ...
        seg3d.statistics.apStats.aspectWeighted.apSin2Theta ;
    std_QAc2ts(dmy, :) = ...
        seg3d.statistics.apStats.aspectWeighted.apCos2ThetaStd ;
    std_QAs2ts(dmy, :) = ...
        seg3d.statistics.apStats.aspectWeighted.apSin2ThetaStd ;
    ste_QAc2ts(dmy, :) = ...
        seg3d.statistics.apStats.aspectWeighted.apCos2ThetaSte ;
    ste_QAs2ts(dmy, :) = ...
        seg3d.statistics.apStats.aspectWeighted.apSin2ThetaSte ;
    % eccentricity-weighted
    mean_QEc2ts(dmy, :) = ...
        seg3d.statistics.apStats.eccentricityWeighted.apCos2Theta ;
    mean_QEs2ts(dmy, :) = ...
        seg3d.statistics.apStats.eccentricityWeighted.apSin2Theta ;
    std_QEc2ts(dmy, :) = ...
        seg3d.statistics.apStats.eccentricityWeighted.apCos2ThetaStd ;
    std_QEs2ts(dmy, :) = ...
        seg3d.statistics.apStats.eccentricityWeighted.apSin2ThetaStd ;
    ste_QEc2ts(dmy, :) = ...
        seg3d.statistics.apStats.eccentricityWeighted.apCos2ThetaSte ;
    ste_QEs2ts(dmy, :) = ...
        seg3d.statistics.apStats.eccentricityWeighted.apSin2ThetaSte ;
    
    mid_ap = seg3d.statistics.apStats.aspectWeighted.apBins ;
    
    %% Plot this timepoint's segmentation in 3d
    % aux_plotCellSegmentation3D(QS, tp, seg3d, imdir, ...
    %    overwriteImages, xyzlims, ~corrected)
    
    glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
    if ~QS.flipy
        glueMesh.vn = - glueMesh.vn ;
    end
    glueMesh.v = glueMesh.v + glueMesh.vn .* 6*QS.APDV.resolution ;
        
    aux_plotCellSegmentation3D(QS, tp, seg3d, imdir, ...
        overwriteImages, xyzlims, ~corrected, glueMesh)
    
    %% Plot as histogram
    edges = linspace(-1, 1, 100) ;
    edgesAR = linspace(1, 3, 100) ;
    edgesTheta = linspace(0, pi, 100) ;
    tmp = histcounts(cos(2* seg3d.qualities.ang1), edges) ;
    cos2thetaM(:, dmy) = tmp / sum(tmp) ;
    tmp = histcounts(sin(2* seg3d.qualities.ang1), edges) ;
    sin2thetaM(:, dmy) = tmp / sum(tmp) ;
    tmp = histcounts(ars, edgesAR) ;
    aspectM(:, dmy) = tmp / sum(tmp) ;
    tmp = histcounts(mod(seg3d.qualities.ang1, pi), edgesTheta) ;
    thetaM(:, dmy) = tmp / sum(tmp) ;
    
    dmy = dmy + 1 ;
end


%% Plot mean +/- pctile over time
if corrected
    segSubDir = 'seg3d_corrected' ;
else
    segSubDir = 'seg3d' ;
end

aux_plotCellSegmentation3DStats


%% Testing for shape characterization
% % Is the result dependent on the distribution of vertices? No! yay.
% xx = 2 * [0, 1, 1, 1, 1, 0];
% yy = [0, 0, 0.5, 2, 3, 1];
% [ geom, iner, cpmo ] = polygeom( xx, yy ) ;
% plot([xx xx(1)], [yy yy(1)], '.-')
% axis equal; hold on;
% 
% 
% xx2 = 2 * [0, 1, 1, 0];
% yy2 = [0, 0, 3, 1];
% [ geom2, iner2, cpmo2 ] = polygeom( xx2, yy2 ) ;
% plot([xx2 xx2(1)], [yy2 yy2(1)], 'o--')
% axis equal
% 
% moi = [iner(4) -iner(6); -iner(6) iner(5)]
% [tmp, tmp2] = eig(moi)
%
% areas(cid) = geom(1) ;
% perim(cid) = geom(4) ;
% pcentroid(cid) = geom(2:3) ;
% moment1(cid) = cpmo(1) ;
% ang1(cid) = cpmo(2) ;
% moment2(cid) = cpmo(3) ;
% ang2(cid) = cpmo(4) ;
% mratio(cid) = moment2(cid) / moment1(cid) ;
% moinertia(cid, :) = [iner(4) iner(6) iner(5)] ;

end

