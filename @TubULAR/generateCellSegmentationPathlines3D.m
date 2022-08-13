function generateCellSegmentationPathlines3D(QS, options)
% generateCellSegmentation3D(QS, options)
% 
%  Load segmentation results from a timepoint t0Pathlines, advect the cell
%  polygons along pathlines from PIV, measure the segmentation properties
%  of the advected "tissue" pattern frozen in the Lagrangian frame as the
%  material deforms.
%
% NPMitchell 2021


% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;
maxCellSize = Inf ;  % maximum size allowed to include a cell
overwrite = false ;
overwriteImages = true ;
useCorrected = true ;
debug = false ;
preview = false ;

[~, ~, ~, xyzlims] = QS.getXYZLims() ;
t0 = QS.t0set() ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'overwritePathlines')
    overwritePathlines = options.overwritePathlines ;
end
if isfield(options, 'preview')
    preview = options.preview ;
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
if isfield(options, 'useCorrectedSegmentation')
    useCorrected = options.useCorrectedSegmentation ;
elseif isfield(options, 'useCorrected')
    useCorrected = options.useCorrected ;
end
if isfield(options, 'segmentationPathlines')
    tmp = options.segmentationPathlines ;
    segVertexPathlines2D = tmp.segVertexPathlines2D ;
    segVertexPathlines3D = tmp.segVertexPathlines3D ;
    cellIDs = tmp.cellIDs ;
    segPathlinesPassedAsOption = true ;
else
    segPathlinesPassedAsOption = false ;
end

% Directory preparation
imdir = fullfile(QS.dir.segmentation, 'pathlines', 'images') ;
if ~exist(imdir, 'dir')
    mkdir(imdir)
end

% Prep for dicing data by lobes
features = QS.getFeatures() ;
folds = features.folds ;
nLobes = size(folds, 2) + 1 ;

%% Create synthetic pathlines for cells from t=0
% synethetic cell vertex pathlines in pullback space
if ~exist(fullfile(QS.dir.segmentation, 'pathlines'), 'dir')
    mkdir(fullfile(QS.dir.segmentation, 'pathlines'))
end
cellVertexPathlineFn = fullfile(QS.dir.segmentation, 'pathlines', ...
    sprintf('cellVertexPathlines_%06dt0.mat', t0)) ;
if ~exist(cellVertexPathlineFn, 'file') || overwritePathlines 
    disp('Creating vertex pathlines from scratch...')
    QS.setTime(t0) ;
    if useCorrected
        seg2d = getCurrentSegmentation2DCorrected(QS) ;
        cellV0 = [] ;
        cellIDs = [] ;
        for cellID = 1:length(seg2d.seg2d.cdat.polygons)
            cVtx = seg2d.seg2d.cdat.polygons{cellID} ;
            if ~isempty(cVtx)
                cellV0 = [cellV0; [cVtx(:, 2), cVtx(:, 1)] ] ;            
                cellIDs = [cellIDs; cellID * ones(size(cVtx, 1), 1)] ;
            end
        end
    else
        seg2d = getCurrentSegmentation2D(QS) ;
        cellV0 = seg2d.seg2d.vdat.v ;
        cellIDs = 1:length(seg2d.seg2d.cdat.polygons) ;
    end
    opts = struct('preview', true) ;
    [segVertexPathlines2D, segVertexPathlines3D] = ...
        QS.samplePullbackPathlines(cellV0, opts) ;
    save(cellVertexPathlineFn, 'segVertexPathlines2D', ...
        'segVertexPathlines3D', 'cellIDs')
else
    disp('Loading vertex pathlines from segmentation advection...')
    if ~segPathlinesPassedAsOption
        load(cellVertexPathlineFn, 'segVertexPathlines2D', ...
            'segVertexPathlines3D', 'cellIDs')
    end
end

% Load cell segmentation 2d at t0
QS.setTime(t0) ;
if useCorrected
    seg02d = QS.getCurrentSegmentation2DCorrected() ;
else
    seg02d = QS.getCurrentSegmentation2D() ;
end

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
    
    outfn = fullfile(QS.dir.segmentation, 'pathlines', ...
        sprintf([QS.fileBase.segmentation3d '.mat'], QS.currentTime)) ;
    if ~exist(outfn, 'file') || overwrite 

        % Obtain the segmentation in 2d --> use reference segmentation
        seg2d = seg02d ;
        % replace vertices with advected ones from t0
        seg2d.seg2d.vdat.v = squeeze(segVertexPathlines2D(tidx, :, :)) ;
        
        % replace centroids with correct ones (ie advected in XY plane)
        if ~useCorrected
            % Polygons are indices into vertices
            nCells = length(seg2d.seg2d.Cdat) ;
            for pp = 1:nCells 
                poly = seg2d.seg2d.cdat.polygons{pp} ;
                if ~isempty(poly)
                    try
                        geom = polygeom( seg2d.seg2d.vdat.v(poly, 1), ...
                            seg2d.seg2d.vdat.v(poly, 2) ) ;
                        seg2d.seg2d.cdat.centroid(pp, :) = geom(2:3) ; 
                    catch
                        error('Could not get polygon information in 3D')
                    end
                end
            end
            
            % Flush superfluous info that was loaded from memory
            seg2d.seg2d.Cdat = rmfield(seg2d.seg2d.Cdat, 'centroid') ;
            seg2d.seg2d.Vdat = rmfield(seg2d.seg2d.Vdat, 'vertxcoord') ;
            seg2d.seg2d.Vdat = rmfield(seg2d.seg2d.Vdat, 'vertycoord') ;
        else
            % Polygons are computed from vertices indexed by cid
            nCells = length(seg2d.seg2d.cdat.polygons) ;
            
            if preview
                tmpCOMs = seg2d.seg2d.cdat.centroid ;
            end
            
            for cid = 1:nCells
                poly = seg2d.seg2d.vdat.v(cellIDs == cid, :) ;
                if ~isempty(poly)
                    geom = polygeom( poly(:, 1), poly(:, 2) ) ;
                    seg2d.seg2d.cdat.centroid(cid, :) = geom(2:3) ; 
                end
            end
            
            % THe COMs should appear suitably advected
            if preview
                plot(tmpCOMs(:, 1), tmpCOMs(:, 2), '.')
                hold on; plot(seg2d.seg2d.cdat.centroid(:, 1), ...
                    seg2d.seg2d.cdat.centroid(:, 2), 'o')
            end
        end
        
        % obtain the current cut mesh in APDV coordinates -->
        % 2d coordinates in plane of MESH
        cutMesh = QS.getCurrentSPCutMeshSmRS() ;

        % Normalize the u coordinate
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;

        % tile annular cut mesh
        tileCount = [1,  1] ;
        [ faces, v2D, v3D ] = tileAnnularCutMesh(cutMesh, tileCount) ;

        % Collate cell vertices as XY
        XY = seg2d.seg2d.vdat.v ;

        % Collate cell centroids
        centroids = seg2d.seg2d.cdat.centroid ;
        
        % Check that we are interpolating over the same domain of
        % parameterization
        if strcmpi(seg2d.coordSys, 'spsme')
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

        
        %% COULD USE THIS LINE INSTEAD OF BELOW!
        % cell vertices in 3d
        % [c3d, vertexMeshFaces] = ...
        %     interpolate2Dpts_3Dmesh(faces, v2D, v3D, uv) ;
        % 
        % if useCorrected
        % 
        %     polygons = seg2d.seg2d.cdat.polygons ;
        %     % prepare structure for use with polygonNetwork3dMeasurements()
        %     poly3d = cell(length(polygons), 1) ;
        %     polygonVtxIds = cell(length(polygons), 1) ;
        %     cellVtx2D = cell(length(polygons), 1) ;
        %     for pp = 1:length(polygons)
        %         poly3d{pp} = c3d(pgonIDs == pp, :) ;
        %         polygonVtxIds{pp} = find(pgonIDs == pp) ;
        %         cellVtx2D{pp} = uv(polygonVtxIds{pp}, :) ;
        %     end
        % 
        %     % Obtain cell vertices in 3d
        %     % cell2d0 = seg2d.seg2d.cdat.polygons{cid} ;
        %     % cellVtx2D = uv(pgonIDs == cid, :) ;
        %     % cellVtx0 = c3d(pgonIDs == cid, :) ;
        %     [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
        %         moment2, ang2, moinertia, cellQ2d, ...
        %         cellMeshFaces, vertexMeshFaces, cellNormals ] = ...
        %         polygonNetwork3dMeasurements(faces, v3D, v2D, uv, polygonVtxIds, cntrds) ;
        % else
        %     % Obtain cell vertices in 3d
        %     % original segmentation stores indices into vertices for
        %     % polygon shapes
        %     % cell2d0 = seg2d.seg2d.vdat.v(seg2d.seg2d.cdat.polygons{cid}, :) ;
        %     polygonVtxIds = seg2d.seg2d.cdat.polygons ;
        %     for pp = 1:length(polygonVtxIds)
        %         cellVtx2D{pp} = uv(polygonVtxIds{cid}, :) ;
        %         % cellVtx0 = c3d(seg2d.seg2d.cdat.polygons{cid}, :) ;
        %     end
        %     [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
        %         moment2, ang2, moinertia, cellQ2d,  cellMeshFaces, vertexMeshFaces] = ...
        %         polygonNetwork3dMeasurements(faces, v3D, v2D, uv, polygonVtxIds, cntrds) ;
        % end
        
        %%
        % cell vertices in 3d
        [c3d, vertexMeshFaces] = ...
            interpolate2Dpts_3Dmesh(faces, v2D, v3D, uv) ;
        
        % Get fieldfaces for centroids
        [cellCntrd3d, cellMeshFaces] = interpolate2Dpts_3Dmesh(faces, v2D, v3D, cntrds) ;
        fN = faceNormal(triangulation(faces, v3D)) ;
        cellNormals = fN(cellMeshFaces, :) ; 
        % Not needed: already normalized by construction
        % cellNormals = cellNormals ./ vecnorm(cellNormals, 2, 2) ;
        
        % Also get vectors point along zeta from cell centroid to lie along y
        jac2d3d = jacobian2Dto3DMesh(v2D, v3D, faces) ;
        % jac3d2d = jacobian3Dto2DMesh(v2D, v3D, faces) ;
        
        %% Get aspect ratios of polygons 
        % rotate to 2d plane, with x along zeta, and y along phi
        areas = nan(nCells, 1) ;
        perim = nan(nCells, 1) ;
        moment1 = nan(nCells, 1) ;
        ang1 = nan(nCells, 1) ;
        moment2 = nan(nCells, 1) ;
        ang2 = nan(nCells, 1) ;
        moinertia = nan(nCells, 2, 2) ;
        cellQ2d = {} ;
        
        for cid = 1:nCells
            % Obtain cell vertices in 3d
            if useCorrected
                cell2d0 = seg2d.seg2d.vdat.v(cellIDs == cid, :) ;
                cellVtx0 = c3d(cellIDs == cid, :) ;
            else
                cell2d0 = seg2d.seg2d.vdat.v(seg2d.seg2d.cdat.polygons{cid}, :) ;
                cellVtx0 = c3d(seg2d.seg2d.cdat.polygons{cid}, :) ;
            end
            cellVtx = cellVtx0 - cellCntrd3d(cid, :) ;
        
            if ~isempty(cellVtx)
        
                % Note: this approach is tricky since we have to map to
                % 3d and then back to 2d.
                % To figure out which direction to take to z, map vec to 3d
        
                % dzeta3d points towards the mapped z axis but in original 3d
                % embedding space
                dzeta3d = (jac2d3d{cellMeshFaces(cid)} * [0, 1]')' ;
                dzeta3d = dzeta3d / vecnorm(dzeta3d, 2, 2) ;
        
                % rotate cell to nearly yz plane
                rot2xy = rotate3dToAlignAxis(cellNormals(cid, :), dzeta3d) ;
                % Note: this one doesn't do the secondary rotation
                % rot = createRotationVector3d(cellNormals(cid, :), [0, 0, 1]) ;
                cell_quasi2d = (rot2xy * cellVtx')' ;
        
                % % Alternative method based on jacobian --> actually this 
                % % is a bit limited since it assumes small cells 
                % % (same size or smaller than the face) to get the 
                % % stretching in each dimension correct
                % jac = jac3d2d{cellMeshFaces(cid)} ;
                % dilation = jacobian2dilation(jac) ;
                % cell2d = (jac * cellVtx0')' ;
                % mapped centroid should be the same as pre-mapped centroid
                % up to rescaling in each dimension separately
                % cntrd2d = jac * cellCntrd(cid, :)' ;
                % cntrds(cid, :)
        
                % Check 2d cell polygon
                if debug
                    disp('debugging cell polygon in 3d and quasi2d...')
                    clf
                    % dchi3d points towards the mapped x axis in embedding space
                    dchi3d = (jac2d3d{cellMeshFaces(cid)} * [1, 0]')' ;
                    dchi3d = dchi3d / vecnorm(dchi3d, 2, 2) ;
        
                    % Plot the cell in 3d and 2d
                    subplot(2, 2, 1)
                    plot(cell2d0(:, 1), cell2d0(:, 2), '.-');
                    axis equal
                    title('cell in 2d pullback XY', 'interpreter', 'latex')
                    subplot(2, 2, 2)
        
                    % Vector to transform = dzeta since this emanates from
                    % centroid just as the normal does.
                    zeta2d = (rot2xy * dzeta3d')' ;
                    try
                        assert(all(abs(zeta2d - [0, 0, 1]) < 1e-7))
                    catch
                        error('Rotation was not successful!')
                    end
        
                    % plot it
                    plot3(cell_quasi2d(:, 1), cell_quasi2d(:, 2), ...
                        cell_quasi2d(:, 3), '.-'); 
                    axis equal
                    title('cell projected onto tangent plane and oriented', ...
                        'interpreter', 'latex')
                    hold on; 
        
                    subplot(2, 2, 3)
                    plot3(cellVtx0(:, 1), cellVtx0(:, 2), cellVtx0(:, 3), ...
                        '.-')
                    hold on;
                    zplus = cellCntrd3d(cid, :) + dzeta3d * mean(var(cellVtx0)); 
                    xplus = cellCntrd3d(cid, :) + dchi3d * mean(var(cellVtx0));
                    plot3dpts([cellCntrd3d(cid, :); zplus], 'r-')
                    plot3dpts([cellCntrd3d(cid, :); xplus], 'g-')
                    title('cell in 3d with computed $(\zeta, \chi)$ coordinates', ...
                        'interpreter', 'latex')
                    axis equal
                end
        
                % Look at yz plane in which cell lives (very nearly, tangent plane)
                %   POLYGEOM( X, Y ) returns area, X centroid,
                %   Y centroid and perimeter for the planar polygon
                %   specified by vertices in vectors X and Y.
                %
                %   [ GEOM, INER, CPMO ] = POLYGEOM( X, Y ) returns
                %   area, centroid, perimeter and area moments of 
                %   inertia for the polygon.
                %   GEOM = [ area   X_cen  Y_cen  perimeter ]
                %   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
                %     u,v are centroidal axes parallel to x,y axes.
                %   CPMO = [ I1     ang1   I2     ang2   J ]
                %     I1,I2 are centroidal principal moments about axes
                %         at angles ang1,ang2.
                %     ang1 and ang2 are in radians.
                %     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv
        
                % Discard 3d info and compute
                if size(cellVtx0, 1) > 2
                    try
                        [ geom, iner, cpmo ] = polygeom( cell_quasi2d(:, 2), ...
                            cell_quasi2d(:, 3) ) ;
                    catch
                        error('here')
                    end
                    areas(cid) = geom(1) ;
                    perim(cid) = geom(4) ;
                    cellQ2d{cid} = cell_quasi2d ;
                    moment1(cid) = cpmo(1) ;
                    ang1(cid) = cpmo(2) ;
                    moment2(cid) = cpmo(3) ;
                    ang2(cid) = cpmo(4) ;
                    mratio(cid) = moment2(cid) / moment1(cid) ;
                    moinertia(cid, :, :) = [iner(4) -iner(6); -iner(6) iner(5)] ;
        
                end
            else
                disp(['bad cell: ' num2str(cid)])
            end
        end
        
        %% COULD USE THIS LINE INSTEAD OF THE ABOVE!
        %[c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
        % moment2, ang2, moinertia, cellQ2d, cellMeshFaces, vertexMeshFaces, cellNormals] = ...
        % polygonNetwork3dMeasurements(faces, v3D, v2D, cellVtx2D, polygons, cntrds)
        %%
        
        
        %% Save results stored in struct
        disp('building advected segmentation in struct...')
        
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
        if useCorrected
            seg3d.cellIDs = cellIDs ;
        else
            seg3d.vdat.NL = seg2d.seg2d.vdat.NL ;
            seg3d.vdat.BL = seg2d.seg2d.vdat.BL ;
            seg3d.vdat.fourfold = seg2d.seg2d.vdat.fourfold ;
        end
        
        % cell data in 3d
        seg3d.cdat.centroids_uv = cntrds ;
        seg3d.cdat.centroids_3d = cellCntrd3d ;
        seg3d.cdat.meshFaces = cellMeshFaces ;
        seg3d.cdat.normals = cellNormals ;
        seg3d.cdat.polygons = seg2d.seg2d.cdat.polygons ;
        
        % cell qualities
        seg3d.qualities = struct() ;
        seg3d.qualities.areas = areas ;
        seg3d.qualities.perim = perim ;
        seg3d.qualities.moment1 = moment1 ;
        seg3d.qualities.moment2 = moment2 ;
        seg3d.qualities.ang1 = ang1 ;
        seg3d.qualities.ang2 = ang2 ;
        seg3d.qualities.mInertia = moinertia ;
        seg3d.qualities.cellQ2d = cellQ2d ;
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
            'ang2', 'angle of long axis in coordSys, in radians', ...
            'cellQ2d', ['quasi-2d cell polygon from embedding space, ', ...
                    'but with cell centroid surface normal rotated ', ...
                    'to be along x axis'], ...
            'nematicTensor', ['n^T * n - 0.5 * [1, 0; 0, 1], ', ...
                    'where n is along long axis'], ...
            'nematicStrength', 'abs(sqrt(MOIEigenvalueRatio)) - 1, strength of elongation' ) ;
            
        % cell statistics 
        % find which are "good" cells to consider
        keep = find(~isnan(ang1) & (areas < maxCellSize) & ...
            moment1 > 0 & moment2 > 0) ;
        seg3d.statistics.keep = keep ;
        seg3d.statistics.maxCellSize = maxCellSize ;
        
        % which coordinate system has been used for segmentation
        coordSys = seg2d.coordSys ;
        disp(['saving segmentation in 3d to: ' outfn])
        save(outfn, 'seg3d', 'coordSys')

        %% Medians of orientation and moment ratio over TIME    
        % Compute cell statistics
        nCells = length(seg3d.qualities.areas) ;
        keep = seg3d.statistics.keep ;
        areas = seg3d.qualities.areas ;

        % Compute mean MOI
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

        % Compute mean MOI, area-weighted 
        % --> isn't MOI already area-weighted? No, it's weighted oddly. 
        weight = areas(keep) ;
        weights = weight ./ nansum(weight) ;

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
        disp(['Saving seg3d now with statistics to: ' outfn])
        save(outfn, 'seg3d', 'coordSys')
    else
        % seg3d = QS.loadCurrentSegmentation3D() ;
        disp(['Loading seg3d with statistics from: ' outfn])
        seg3d = load(outfn) ;
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
    glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
    if ~QS.flipy
        glueMesh.vn = -glueMesh.vn ;
    end
    glueMesh.v = glueMesh.v + glueMesh.vn .* 6*QS.APDV.resolution ;
    
    aux_plotCellSegmentation3D(QS, tp, seg3d, imdir, ...
        overwriteImages, xyzlims, ~useCorrected, glueMesh)
    
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

segSubDir = 'pathlines' ;
mid_ap = seg3d.statistics.apStats.aspectWeighted.apBins ;
aux_plotCellSegmentation3DStats



%% Compare to true segmentation GLOBALLY -- nematic strength and direction 
timeList = {timePoints, timePoints(timePoints < t0 + 75 & timePoints > -30)} ;
timeStr = {'', '_tlimit'} ;
stdsteStr = {'_std', '_ste'};
for std_ste = 1:2
    for pp = 1:2
        close all 
        fig = figure('units', 'centimeters', 'position', [0,0,figH,figH]) ;
        time2do = timeList{pp} ;
        pInds = ismember(timePoints, time2do) ;

        imfn = fullfile(QS.dir.segmentation, 'pathlines', ...
            ['cell_anisotropy_global_signed_COMPARE', timeStr{pp}]) ;
        % Collate results
        kk = 1;
        meanQAspectsTrue = zeros(length(time2do), 1) ;
        meanQAspectStdsTrue = zeros(length(time2do), 1) ;
        meanQAspectStesTrue = zeros(length(time2do), 1) ;
        meanQThetasTrue = meanQAspectsTrue ;
        for tp = time2do
            QS.setTime(tp) ;
            if useCorrected     
                seg3d = QS.getCurrentSegmentation3DCorrected() ; 
            else
                seg3d = QS.getCurrentSegmentation3D() ; 
            end
            meanQAspectsTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQAspect ;
            meanQAspectStdsTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQAspectStd ;
            meanQAspectStesTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQAspectSte ;
            meanQThetasTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQTheta ;
            meanQThetaStdsTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQThetaStd ;
            meanQThetaStesTrue(kk) = seg3d.seg3d.statistics.meanQ.aspectWeighted.meanQThetaSte ;
            kk = kk + 1;
        end


        c2tTrue = cos(2*meanQThetasTrue(pInds))  ;
        s2tTrue = sin(2*meanQThetasTrue(pInds))  ;
        trueAn = meanQAspectsTrue(pInds) - 1 ;
        trueAnStd = meanQAspectStdsTrue(pInds)  ;
        trueAnSte = meanQAspectStesTrue(pInds)  ;
        trueThetaStd = meanQThetaStdsTrue(pInds) ;
        trueThetaSte = meanQThetaStesTrue(pInds) ;
        trueLine = c2tTrue .* trueAn ;
        trueStds = sqrt( (c2tTrue .* trueAnStd).^2 + ...
            (trueAn .* s2tTrue * 2 .* trueThetaStd(:)).^2) ;
        trueStes = sqrt( (c2tTrue .* trueAnSte).^2 + ...
            (trueAn .* s2tTrue * 2 .* trueThetaSte(:)).^2) ;
        % ratios
        trueAR = meanQAspectsTrue(pInds) ;
        trueLineR = c2tTrue .* trueAR ;
        trueStdsR = sqrt( (c2tTrue .* trueAnStd).^2 + ...
            (trueAR .* s2tTrue * 2 .* trueThetaStd(:)).^2) ;
        trueStesR = sqrt( (c2tTrue .* trueAnSte).^2 + ...
            (trueAR .* s2tTrue * 2 .* trueThetaSte(:)).^2) ;

        c2t = cos(2*mQAtheta(pInds))  ;
        s2t = sin(2*mQAtheta(pInds))  ;
        An = squeeze(mQAar(pInds)) - 1 ;
        AnStd = squeeze(mQAarStd(pInds)) ;
        AnSte = squeeze(mQAarSte(pInds)) ;
        ThetaStd = squeeze(mQAthetaStd(pInds)) ;
        ThetaSte = squeeze(mQAthetaSte(pInds)) ;
        midline = c2t .* An ;
        stds = sqrt( (c2t .* AnStd).^2 + (An .* s2t *2 .* ThetaStd).^2) ;
        stes = sqrt( (c2t .* AnSte).^2 + (An .* s2t *2 .* ThetaSte).^2) ;
        % ratios
        AR = squeeze(mQAar(pInds)) ;
        midlineR = c2t .* AR ;
        stdsR = sqrt( (c2t .* AnStd).^2 + (AR .* s2t *2 .* ThetaStd).^2) ;
        stesR = sqrt( (c2t .* AnSte).^2 + (AR .* s2t *2 .* ThetaSte).^2) ;


        timestamps = timePoints - t0 ;
        if contains(QS.timeUnits, 'min')
            timestamps = timestamps / 60 ;
            timeunits = 'hr';
        else
            timeunits = QS.timeUnits ;
        end

        % Transformation (flip sign or subtract initial value)
        midline0 = midline ; 
        trueLine0 = trueLine ;
        midline = midline - midline(1) ;
        trueLine = trueLine - trueLine(1) ;


        if std_ste == 1
            hs = errorbar(midline, trueLine, trueStds, trueStds, stds, stds) ;
        else
            [midline0, inds] = sort(midline) ;
            trueLine0 = trueLine(inds) ;
            trueStds0 = trueStds(inds) ;
            trueStes0 = trueStes(inds) ;
            bounds = [trueLine0-abs(trueStds0), trueLine0+abs(trueStds0)] ;
            lowerB = bounds(:, 1)' ;
            upperB = bounds(:, 2)' ;
            fill([midline0, fliplr(midline0)], ...
                [lowerB, fliplr(upperB)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs = errorbar(midline0, trueLine0, ...
                trueStes0, trueStes0, 'color', colors(1, :)) ;
        end
        hold on;
        ylims = [ylim() xlim()] ;

        % Mark zero line
        plot([0,2], [0,2], 'k--', 'HandleVisibility','off')
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        ylabel('cell shear, $\Delta(a/b)$',   'interpreter', 'latex')
        xlabel('tissue shear, $\Delta(a/b)$',   'interpreter', 'latex')
        axis equal
        % ylim([-max(abs(ylims)), max(abs(ylims))])
        % xlim([-max(abs(ylims)), max(abs(ylims))])
        
        ylim([-Inf, max(abs(ylims))])
        xlim([-Inf, max(abs(ylims))])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '.pdf'])
        


        %% Global Difference over time
        clf
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        if std_ste == 1
            hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            hold on;
            hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
        else
            fill(x2, [midline - stds, fliplr(midline + stds)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            fill(x2, [trueLine - trueStds; flipud(trueLine + trueStds)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs2 = errorbar(timestamps(pInds), trueLine(:),trueStes(:), 'color', colors(1, :)) ;
            hold on;
            hs = errorbar(timestamps(pInds), midline(:), stes(:), 'color', colors(2, :)) ;
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear',   'interpreter', 'latex')
        legend([hs, hs2], {'tissue shear', 'cell shape change'})
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff.pdf'])
        
        
        
        %% Global curves and difference over time -- shaded ste or shaded std
        clf
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        diffLine = midline(:) - trueLine(:) ;
        diffStds = sqrt(stds(:).^2 + trueStds(:).^2) ;
        diffStes = sqrt(stes(:).^2 + trueStes(:).^2) ;
        if std_ste == 1
            % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            % hold on;
            % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
            hs = fill(x2, [midline - stds, fliplr(midline + stds)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [trueLine - trueStds; flipud(trueLine + trueStds)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [diffLine - diffStds; flipud(diffLine + diffStds)], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            errorbar(timestamps(pInds), trueLine(:),trueStes(:), 'color', colors(1, :)) ;
            hold on;
            errorbar(timestamps(pInds), midline(:), stes(:), 'color', colors(2, :)) ;
            errorbar(timestamps(pInds), diffLine(:), diffStes(:), 'color', colors(3, :)) ;
        else
            hs = fill(x2, [midline - stes, fliplr(midline + stes)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [trueLine - trueStes; flipud(trueLine + trueStes)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [diffLine - diffStes; flipud(diffLine + diffStes)], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
         legend([hs, hs2, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2.pdf'])
        set(gca, 'YScale', 'log')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_semilog.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_semilog.pdf'])
        
        
        %% Global curves and difference over time, relative to start of plot, moving mean -- shaded ste or shaded std
        clf
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        diffLine = (midline(:) - trueLine(:)) ;
        diffStds = sqrt(stds(:).^2 + trueStds(:).^2) ;
        diffStes = sqrt(stes(:).^2 + trueStes(:).^2) ;
        if std_ste == 1
            % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            % hold on;
            % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
            hs = fill(x2, [movmean(midline(:) - stds(:), 3); ...
                flipud(movmean(midline(:) + stds(:), 3))], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [movmean(trueLine(:) - trueStds(:), 3); ...
                flipud(movmean(trueLine(:) + trueStds(:), 3))], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [movmean(diffLine - diffStds, 3); ...
                flipud(movmean(diffLine + diffStds, 3))], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            errorbar(timestamps(pInds), movmean(trueLine, 3),...
                movmean(trueStes(:), 3), 'color', colors(1, :)) ;
            hold on;
            errorbar(timestamps(pInds), movmean(midline, 3), ...
                movmean(stes(:), 3), 'color', colors(2, :)) ;
            errorbar(timestamps(pInds), movmean(diffLine(:), 3), ...
                movmean(diffStes(:), 3), 'color', colors(3, :)) ;
        else
            hs = fill(x2, [movmean(midline - stes, 3), ...
                fliplr(movmean(midline + stes, 3))], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [movmean(trueLine - trueStes, 3); ...
                flipud(movmean(trueLine + trueStes, 3))], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [movmean(diffLine - diffStes, 3); ...
                flipud(movmean(diffLine + diffStes, 3))], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
        try
         legend([hs, hs2, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        catch
         legend([hs, hs2(1), hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        end
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean.pdf'])
        set(gca, 'YScale', 'log')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog.pdf'])
        xlim([-0.25, Inf])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_tx.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_tx.pdf'])
        
        %% Global curves and difference over time, relative to t0, with moving mean -- shaded ste or shaded std
        clf
        idx = find(timestamps < 0, 1, 'last') ;
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        midlineT0 = midline0(:) - midline0(idx) ;
        trueLineT0 = trueLine0(:)-trueLine0(idx) ;
        diffLine = (midlineT0 - trueLineT0) ;
        diffStds = sqrt(stds(:).^2 + trueStds(:).^2) ;
        diffStes = sqrt(stes(:).^2 + trueStes(:).^2) ;
        dataPlotted = struct('tissueShear', movmean(midlineT0(:), 3) , ...
            'cellShear', movmean(trueLineT0(:), 3), ...
            'intercalations', movmean(diffLine(:), 3), ...
            'std_tissueShear', movmean(stds(:), 3) , ...
            'std_cellShear', movmean(trueStds(:), 3), ...
            'std_intercalations', movmean(diffStds(:), 3), ...
            'ste_tissueShear', movmean(stes(:), 3) , ...
            'ste_cellShear', movmean(trueStes(:), 3), ...
            'ste_intercalations', movmean(diffStes(:), 3)) ;
        if std_ste == 1
            % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            % hold on;
            % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
            hs2 = fill(x2, [movmean(midlineT0(:) - stds(:), 3); ...
                flipud(movmean(midlineT0(:) + stds(:), 3))], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs = fill(x2, [movmean(trueLineT0(:) - trueStds(:), 3); ...
                flipud(movmean(trueLineT0(:) + trueStds(:), 3))], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [movmean(diffLine - diffStds, 3); ...
                flipud(movmean(diffLine + diffStds, 3))], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            errorbar(timestamps(pInds), movmean(trueLineT0, 3),...
                movmean(trueStes(:), 3), 'color', colors(1, :)) ;
            hold on;
            errorbar(timestamps(pInds), movmean(midlineT0, 3), ...
                movmean(stes(:), 3), 'color', colors(2, :)) ;
            errorbar(timestamps(pInds), movmean(diffLine(:), 3), ...
                movmean(diffStes(:), 3), 'color', colors(3, :)) ;
        else
            hs2 = fill(x2, [movmean(midlineT0 - stes(:), 3), ...
                flipud(movmean(midlineT0 + stes(:), 3))], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs = fill(x2, [movmean(trueLineT0 - trueStes, 3); ...
                flipud(movmean(trueLineT0 + trueStes, 3))], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [movmean(diffLine - diffStes, 3); ...
                flipud(movmean(diffLine + diffStes, 3))], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
        try
         legend([hs2, hs, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        catch
         legend([hs2(1), hs, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        end
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_t0.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_t0.pdf'])
        set(gca, 'YScale', 'log')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_t0.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_t0.pdf'])
        xlim([0, Inf])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_t0x.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_movmean_semilog_t0x.pdf'])
        save([imfn, stdsteStr{std_ste}, '_diff2_movmean_DATA.mat'], 'dataPlotted')
        
        
        %% Global curves and difference over time divided by first entry RATIO -- shaded ste or shaded std
        clf
        x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        diffLine1 = midline(:) - trueLine(:) + 1;
        diffStds = sqrt(stds(:).^2 + trueStds(:).^2) ;
        diffStes = sqrt(stes(:).^2 + trueStes(:).^2) ;
        if std_ste == 1
            % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
            % hold on;
            % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
            hs = fill(x2, [midlineR - stdsR, fliplr(midlineR + stdsR)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [trueLineR - trueStdsR; flipud(trueLineR + trueStdsR)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [diffLine1 - diffStds; flipud(diffLine1 + diffStds)], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            errorbar(timestamps(pInds), trueLineR(:),trueStesR(:), 'color', colors(1, :)) ;
            hold on;
            errorbar(timestamps(pInds), midlineR(:), stesR(:), 'color', colors(2, :)) ;
            errorbar(timestamps(pInds), diffLine1(:), diffStes(:), 'color', colors(3, :)) ;
        else
            hs = fill(x2, [midlineR - stesR, fliplr(midlineR + stesR)], ...
                 colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
             hold on;
            hs2 = fill(x2, [trueLineR - trueStesR; flipud(trueLineR + trueStesR)], ...
                 colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
            hs3 = fill(x2, [diffLine1 - diffStes; flipud(diffLine1 + diffStes)], ...
                 colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
                 'HandleVisibility', 'off');
        end
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
         legend([hs, hs2, hs3], {'tissue stretch ratio', 'cell stretch ratio', 'inferred intercalations'})
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diffRatio.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diffRatio.pdf'])
        set(gca, 'YScale', 'log')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diffRatio_semilog.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_diffRatio_semilog.pdf'])
        
        % %% plot on semilog scale like in tissue tectonics
        % clf
        % x2 = [timestamps(pInds), fliplr(timestamps(pInds))] ;
        % if std_ste == 1
        %     % hs2 = errorbar(timestamps(pInds), trueLine(:),trueStds(:)) ;
        %     % hold on;
        %     % hs = errorbar(timestamps(pInds), midline(:), stds(:)) ;
        %     hs2 = fill(x2, [midline0 - stds, fliplr(midline + stds)], ...
        %          colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        %          'HandleVisibility', 'off');
        %      hold on;
        %     hs = fill(x2, [trueLine0 - trueStds; flipud(trueLine + trueStds)], ...
        %          colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        %          'HandleVisibility', 'off');
        %     hs3 = fill(x2, [diffLine - diffStds; flipud(diffLine + diffStds)], ...
        %          colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        %          'HandleVisibility', 'off');
        %     errorbar(timestamps(pInds), trueLine(:),trueStes(:), 'color', colors(1, :)) ;
        %     hold on;
        %     errorbar(timestamps(pInds), midline(:), stes(:), 'color', colors(2, :)) ;
        %     errorbar(timestamps(pInds), diffLine(:), diffStes(:), 'color', colors(3, :)) ;
        % else
        %     hs2 = fill(x2, [midline0 - stes, fliplr(midline + stes)], ...
        %          colors(2, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        %          'HandleVisibility', 'off');
        %      hold on;
        %     hs = fill(x2, [trueLine0 - trueStes; flipud(trueLine + trueStes)], ...
        %          colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        %          'HandleVisibility', 'off');
        %     hs3 = fill(x2, [diffLine - diffStes; flipud(diffLine + diffStes)], ...
        %          colors(3, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        %          'HandleVisibility', 'off');
        % end
        % % Labels
        % % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        % xlabel(['time [' timeunits ']'],   'interpreter', 'latex')
        % ylabel('tissue shear, cell shear, intercalations',   'interpreter', 'latex')
        %  legend([hs, hs2, hs3], {'tissue shear', 'cell shape change', 'inferred intercalations'})
        % saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_semilog.png'])
        % saveas(gcf, [imfn, stdsteStr{std_ste}, '_diff2_semilog.pdf'])
    end
end





%% Compare to true segmentation  -- nematic strength and direction for each lobe 
timeList = {timePoints, timePoints(timePoints < t0 + 75 & timePoints > -30)} ;
timeStr = {'', '_tlimit'} ;
stdsteStr = {'_std', '_ste'};
for std_ste = 1:2
    % for each timespan
    for pp = 1:2
        close all 
        fig = figure('units', 'centimeters', 'position', [0,0,figH,figH]) ;
        time2do = timeList{pp} ;
        pInds = ismember(timePoints, time2do) ;

        imfn = fullfile(QS.dir.segmentation, 'pathlines', ...
            ['cell_anisotropy_lobes_signed_COMPARE', timeStr{pp}]) ;
        % Collate results
        kk = 1;
        meanQLobeAspectsTrue = zeros(nLobes, length(time2do)) ;
        meanQLobeAspectStdsTrue = zeros(nLobes, length(time2do)) ;
        meanQLobeAspectStesTrue = zeros(nLobes, length(time2do)) ;
        meanQLobeThetasTrue = meanQLobeAspectsTrue ;
        for tp = time2do
            QS.setTime(tp) ;
            if useCorrected     
                seg3d = QS.getCurrentSegmentation3DCorrected() ; 
            else
                seg3d = QS.getCurrentSegmentation3D() ; 
            end
            meanQLobeAspectsTrue(:, kk) = seg3d.seg3d.statistics.lobes.aspectWeighted.meanQLobeAspect ;
            meanQLobeAspectStdsTrue(:, kk) = seg3d.seg3d.statistics.lobes.aspectWeighted.meanQLobeAspectStd ;
            meanQLobeAspectStesTrue(:, kk) = seg3d.seg3d.statistics.lobes.aspectWeighted.meanQLobeAspectSte ;
            meanQLobeThetasTrue(:, kk) = seg3d.seg3d.statistics.lobes.aspectWeighted.meanQLobeTheta ;
            meanQLobeThetaStdsTrue(:, kk) = seg3d.seg3d.statistics.lobes.aspectWeighted.meanQLobeThetaStd ;
            meanQLobeThetaStesTrue(:, kk) = seg3d.seg3d.statistics.lobes.aspectWeighted.meanQLobeThetaSte ;
            kk = kk + 1;
        end

        hs = cell(nLobes, 1) ;
        for lobe = 1:nLobes

            c2tTrue = cos(2*meanQLobeThetasTrue(lobe, pInds))  ;
            s2tTrue = sin(2*meanQLobeThetasTrue(lobe, pInds))  ;
            trueAn = squeeze(meanQLobeAspectsTrue(lobe, pInds)) - 1 ;
            trueAnStd = squeeze(meanQLobeAspectStdsTrue(lobe, pInds))  ;
            trueAnSte = squeeze(meanQLobeAspectStesTrue(lobe, pInds))  ;
            trueThetaStd = squeeze(meanQLobeThetaStdsTrue(lobe, pInds)) ;
            trueThetaSte = squeeze(meanQLobeThetaStesTrue(lobe, pInds)) ;
            trueLine = c2tTrue .* trueAn ;
            trueStds = sqrt( (c2tTrue .* trueAnStd).^2 + ...
                (trueAn .* s2tTrue * 2 .* trueThetaStd).^2) ;
            trueStes = sqrt( (c2tTrue .* trueAnSte).^2 + ...
                (trueAn .* s2tTrue * 2 .* trueThetaSte).^2) ;
            % ratios
            trueAR = squeeze(meanQLobeAspectsTrue(lobe, pInds))  ;
            trueLineR = c2tTrue .* trueAR ;
            trueStdsR = sqrt( (c2tTrue .* trueAnStd).^2 + ...
                (trueAR .* s2tTrue * 2 .* trueThetaStd).^2) ;
            trueStesR = sqrt( (c2tTrue .* trueAnSte).^2 + ...
                (trueAR .* s2tTrue * 2 .* trueThetaSte).^2) ;            

            % subplot(ceil(nLobes * 0.5), 2, lobe)
            c2t = cos(2*meanQALobeThetas(lobe, pInds))  ;
            s2t = sin(2*meanQALobeThetas(lobe, pInds))  ;
            An = squeeze(meanQALobeAspects(lobe, pInds)) - 1 ;
            AnStd = squeeze(meanQALobeAspectStds(lobe, pInds)) ;
            AnSte = squeeze(meanQALobeAspectStes(lobe, pInds)) ;
            ThetaStd = squeeze(meanQALobeThetaStds(lobe, pInds)) ;
            ThetaSte = squeeze(meanQALobeThetaStes(lobe, pInds)) ;
            midline = c2t .* An ;
            stds = sqrt( (c2t .* AnStd).^2 + (An .* s2t *2 .* ThetaStd).^2) ;
            stes = sqrt( (c2t .* AnSte).^2 + (An .* s2t *2 .* ThetaSte).^2) ;
            % ratios
            AR = squeeze(meanQALobeAspects(lobe, pInds)) ;
            midlineR = c2t .* AR ;
            stdsR = sqrt( (c2t .* AnStd).^2 + (AR .* s2t *2 .* ThetaStd).^2) ;
            stesR = sqrt( (c2t .* AnSte).^2 + (AR .* s2t *2 .* ThetaSte).^2) ;
                       
            
            timestamps = timePoints - t0 ;
            if contains(QS.timeUnits, 'min')
                timestamps = timestamps / 60 ;
                timeunits = 'hr';
            else
                timeunits = QS.timeUnits ;
            end

            % Transformation (flip sign or subtract initial value)
            midline = midline - midline(1) ;
            trueLine = trueLine - trueLine(1) ;
            
            
            if std_ste == 1
                hs{lobe} = errorbar(midline, trueLine, trueStds, trueStds, stds, stds) ;
            else
                [midline, inds] = sort(midline) ;
                trueLine = trueLine(inds) ;
                trueStds = trueStds(inds) ;
                trueStes = trueStes(inds) ;
                bounds = [trueLine-abs(trueStds); trueLine+abs(trueStds)] ;
                lowerB = min(bounds, [], 1) ;
                upperB = max(bounds, [], 1) ;
                fill([midline, fliplr(midline)], ...
                    [lowerB, fliplr(upperB)], ...
                     colors(lobe, :), 'facealpha', 0.2, 'edgecolor', 'none', ...
                     'HandleVisibility', 'off');
                 hold on;
                hs{lobe} = errorbar(midline, trueLine, ...
                    trueStes, trueStes, 'color', colors(lobe, :)) ;
            end
            hold on;
            %hs{lobe} = plot(midline, tr, '.-', 'color', colors(lobe, :)) ;

            legendentries{lobe} = ['chamber ' num2str(lobe)] ;
        end
        ylims = [ylim() xlim()] ;

        % Mark zero line
        plot([0,2], [0,2], 'k--', 'HandleVisibility','off')
        % Labels
        % legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        ylabel('cell shear, $\Delta(a/b)$',   'interpreter', 'latex')
        xlabel('tissue shear, $\Delta(a/b)$',   'interpreter', 'latex')
        axis equal
        % ylim([-max(abs(ylims)), max(abs(ylims))])
        % xlim([-max(abs(ylims)), max(abs(ylims))])
        
        ylim([-Inf, max(abs(ylims))])
        xlim([-Inf, max(abs(ylims))])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '.png'])
        saveas(gcf, [imfn, stdsteStr{std_ste}, '.pdf'])
        
        % save with legend
        legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
        saveas(gcf, [imfn, stdsteStr{std_ste}, '_legend.pdf'])
        
    end
end



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


