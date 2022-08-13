function visualizeTracking3D(QS, options)
% visualizeTracking3D(QS, options)
% Draw cell contours (if method == segmentation) or nuclear track points
% (if method == nuclei) and possibly connections between pairs of cells (if
% drawGeodesics == true) in 3D on top of texturePatch image.
%
% Parameters
% ----------
% options : struct with fields
%   timePoints : timepoints to visualize
%   method : 'segmentation' or 'nuclei'
%   subdir : 'labeled_groundTruth'
%   selectPairs : int (default=0 for no pair selection)
%       number of cell pairs to select from segmentation for tracking.
%       Selection occurs at t=t0.
%   drawGeodesics : bool
%   geodesicStyle : 'pullback' or 'geodesic'
%       'pullback' draws a curve in 3d that connects the 2d segmentation pts
%       as straight line in pullback space projected into 3d
%       'geodesic' computes geodesic pairs on the surface using CGAL
%   pairColors : selectPairs x 1 cell array
%       selectPairs x 1 cell array of colors (either string
%       specifiers or 3x1 float array colors with values between 0-1.0)
%       pairColors{ii} is the color for the pair path connecting tracked
%       cells in pairIDs(ii, :) 
%   viewAngles : [-20,20] is default
%       must match texturePatch angles for overlays to be accurate
%   t0forPairs : timepoint at which we select which pairs to visualize
%   overlayStyle : 'mask', 'add'
%       whether to add or mask the segmentation/tracking/geodesics over the
%       texturePatch image
%   lw_geoP : float or int 
%       linewidth for the geodesic paths
%   roiColor : string specifier for color or 3x1 float between 0-1.0
%       color for the rectangular ROI boundary
%   overwriteROI : bool
%       if the ROI polygon (projected rectangle at t0) exists on disk,
%       overwrite it with a newly calculated ROI
%

xwidth = 16 ; % cm
ywidth = 10 ; % cm
lw = 0.01 ;
lw_geoP = 2 ;
roiColor = [226, 212, 200]/255. ;
viewAngles = [-20, 20] ;


% Glean options from texturePatch run
metafn = fullfile(QS.dir.texturePatchIm, 'overlays', 'metadat.mat') ;
load(metafn, 'metadat', 'Options')
metadat.subdir = 'overlays';

% Default options
t0 = QS.t0set() ;
timePoints = QS.xp.fileMeta.timePoints ;
maxNtracks = 300 ;
subdir = 'labeled_groundTruth';
method = 'segmentation' ;
overwrite = false ;
preview = false ;
coordSys = 'spsme' ;
blackFigure = false ;
selectPairs = 0 ;
t0forPairs = QS.t0set() ;
drawGeodesics = true ;
drawPairRectangle = false ;
overwriteROI = false ;
geodesicStyle = 'pullback' ; % 'geodesic', 'pullback', 'Lagrangian'
overlayStyle = 'add' ; % 'mask', 'add'
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'blackFigure')
    blackFigure = options.blackFigure ;
end
if isfield(options, 'method')
    method = options.method ;
end
if isfield(options, 'viewAngles')
    viewAngles = options.viewAngles ;
end
if isfield(options, 'selectPairs')
    selectPairs = options.selectPairs ;
end
if isfield(options, 'pairIDs')
    pairIDs = options.pairIDs ;
end
if isfield(options, 't0forPairs')
    t0forPairs = options.t0forPairs ;
end
if isfield(options, 'drawPairRectangle')
    drawPairRectangle = options.drawPairRectangle ;
end
if isfield(options, 'overwriteROI')
    overwriteROI = options.overwriteROI ;
end
if isfield(options, 'drawGeodesics')
    drawGeodesics = options.drawGeodesics ;
end
if isfield(options, 'geodesicStyle')
    geodesicStyle = options.geodesicStyle ;
end
if isfield(options, 'overlayStyle')
    overlayStyle = options.overlayStyle ;
end
if isfield(options, 'pairColors')
    pairColors = options.pairColors ;
else
    pairColors = cell(selectPairs, 1) ;
    for ii = 1:selectPairs 
        pairColors{ii} = 'w' ;
    end
end
if isfield(options, 'lw_geoP')
    lw_geoP = options.lw_geoP ;
end

% Obtain tracks' filename on disk -- segmentation or nuclei/points
if strcmpi(method, 'segmentation')
    trackOutfn = fullfile(QS.dir.tracking, subdir, 'seg2d_%06d.mat') ;
else
    trackOutfn = fullfile(QS.dir.tracking,  'manualTracks.mat') ;
end
if isfield(options, 'trackOutfn')
    trackOutfn = options.trackOutfn ;
end
% filename for roi rectangle as pathlines for all timepoints
roiPathlineRectangleFn = fullfile(QS.dir.tracking, subdir, ...
    sprintf('roiPathlineRectangle_t0_%06d.mat', t0forPairs)) ;

% More options
if isfield(options, 'buffer') || isfield(options, 'bufferXYZ')
    if isfield(options, 'buffer')
        bufferXYZ = options.buffer ;
    elseif isfield(options, 'bufferXYZ')
        bufferXYZ = options.bufferXYZ ;
    end
    % interpret bufferXYZ if length-1 float or length(3) array
    if numel(bufferXYZ) == 1
        bufferXYZ = [-bufferXYZ, bufferXYZ] ;
    elseif numel(bufferXYZ) == 3
        bufferXYZ = [-bufferXYZ(1), bufferXYZ(1); ...
            -bufferXYZ(2), bufferXYZ(2); ...
            -bufferXYZ(3), bufferXYZ(3)] ;
    else
        assert(all(size(bufferXYZ) == [3, 2]))
    end
    % Get XYZlimits
    [~, ~, ~, xyzlims] = QS.getXYZLims() ;
    xyzlims = xyzlims + bufferXYZ ;
    metadat.xyzlim = xyzlims ;
else
    xyzlims = metadat.xyzlim ;
end
if isfield(options, 'subdir')
    subdir = options.subdir ;
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end

%% Select which pairs to visualize
if selectPairs
    if isempty(pairIDs) || any(isnan(pairIDs(:)))
        tidx = QS.xp.tIdx(t0forPairs) ;

        % Load image for t0
        if strcmpi(coordSys, 'spsme')
            % Get the image size and mesh
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, t0forPairs)) ;
        elseif strcmpi(coordSys, 'spsm')
            % Get the image size and mesh
            im = imread(sprintf(QS.fullFileBase.im_sp_sm, t0forPairs)) ;
        else
            error('handle coordSys here')
        end

        if strcmpi(method, 'segmentation')
            error('handle case here')
        else
            seg2dFn = trackOutfn ;
            if exist(seg2dFn, 'file') 
                load(seg2dFn, 'tracks')       
                centroids = nan(length(tracks), 2) ;
                for trackID = 1:length(tracks)
                    centroids(trackID, :) = tracks{trackID}(tidx, 1:2) ;     
                end
            else
                error('Could not find tracks on disk')
            end
        end
        close all
        imshow(im)
        hold on;
        plot(centroids(:, 1), centroids(:, 2), '.')
        
        if isempty(pairIDs)
            pxy = ginput(2 * selectPairs) ;
            pidx = pointMatch(pxy, centroids) ;
        else
            id2plot = pairIDs(~isnan(pairIDs)) ;
            plot(centroids(id2plot, 1), ...
                centroids(id2plot, 2), 's')
            pidx = pairIDs ;
            id2fill = find(isnan(pairIDs)) ;
            xypts = ginput(length(id2fill)) ;
            pidx(id2fill) = pointMatch(xypts, centroids) ;
            pidx = pidx(:) ;
        end
        
        % break into pairs
        pairIDs = reshape(pidx, [length(pidx)*0.5, 2])  ;
    end
    disp(['PairIDs = '])
    pairIDs    
end
            
%% Plot each timepoint in 3d      
allAreas = nan(length(timePoints), maxNtracks) ;
timePoints = sort(timePoints) ;
tidx2do = find( ismember(QS.xp.fileMeta.timePoints, timePoints)) ;
dmyii = 1  ;
for tidx = tidx2do 
    tp = timePoints(dmyii) ;
    QS.setTime(tp) ;

    %% PREPARE MESH
    % load mesh and image
    if strcmpi(coordSys, 'spsme')
        % Get the image size and mesh
        im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
        cutMesh = QS.getCurrentSPCutMeshSm() ;
        glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
        glueMesh.u(:, 1) = glueMesh.u(:, 1) / max(glueMesh.u(:, 1)) ;
        shiftY = size(im, 1) * 0.5 ;
        doubleCovered = true ;
    elseif strcmpi(coordSys, 'spsm')
        % Get the image size and mesh
        im = imread(sprintf(QS.fullFileBase.im_sp_sm, tp)) ;
        cutMesh = QS.getCurrentSPCutMeshSm() ;
        glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
        glueMesh.u(:, 1) = glueMesh.u(:, 1) / max(glueMesh.u(:, 1)) ;
        shiftY = 0 ;
        doubleCovered = false ;
    else
        error('did not recognize coordSys')
    end
    
    % Load/push/tile annular cutmesh
    disp('Loading meshes')
    cutMesh = QS.getCurrentSPCutMeshSmRS() ;
    glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
    cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1) ) ;
    glueMesh.u(:, 1) = glueMesh.u(:, 1) / max(glueMesh.u(:, 1) ) ;


    % Make sure vertex normals are normalized and outward facing for
    % cutMesh, but not for glueMesh
    if isfield(options, 'normal_shift')
        normal_shift = options.normal_shift ;
    else
        normal_shift = metadat.normal_shift * QS.APDV.resolution ;
    end
    if QS.flipy
        % glueMesh vertex normals point IN, so we can shrink it a bit
        cutMesh.vn = - cutMesh.vn ;
    else
        % glueMesh vertex normals point IN, so we can shrink it a bit
        glueMesh.vn = - glueMesh.vn ;
    end
    cutMesh.vn = cutMesh.vn ./ sqrt( sum( cutMesh.vn.^2, 2 ) );
    % Normally evolve vertices
    cutMesh.v = cutMesh.v + normal_shift .* cutMesh.vn;
    glueMesh.v = glueMesh.v + normal_shift .* glueMesh.vn ;
    % Check that faces are outward facing
    % trisurf(triangulation(cutMesh.f, cutMesh.v), cutMesh.vn(:, 1), ...
    %     'edgecolor', 'none')


    %% PREPARE CELLS OR NUCLEI
    if strcmpi( method, 'segmentation')
        seg2dFn = sprintf(trackOutfn, tp) ;
        outFigFn = fullfile(QS.dir.tracking, subdir, 'embeddingFigure', sprintf('trackSeg3d_%06d.png', tp)) ;
        outImFn = fullfile(QS.dir.tracking, subdir, 'embeddingTexture', sprintf('trackSeg3d_%06d.png', tp)) ;
        if ~exist(fullfile(QS.dir.tracking, subdir, 'embeddingFigure'), 'dir')
            mkdir(fullfile(QS.dir.tracking, subdir, 'embeddingFigure'))
        end
        if ~exist(fullfile(QS.dir.tracking, subdir, 'embeddingTexture'), 'dir')
            mkdir(fullfile(QS.dir.tracking, subdir, 'embeddingTexture'))
        end
        if exist(seg2dFn, 'file') && ~overwrite
            load(seg2dFn, 'seg2d', 'segIm')
            centroids = seg2d.cdat.centroid ;
            polygons = seg2d.cdat.polygons ;
            cellVtx2D = seg2d.vdat.v ;
            cellIDs = seg2d.cdat.vertexCellIDs ;
        else
            % load segmentation tracks
            segIm = load(fullfile(QS.dir.tracking, subdir, ...
                    sprintf('tracks_label_%06d.mat', tp))) ;
            segIm = segIm.imlabel ;
            nCells = max(segIm(:)) ;

            % check it
            if preview 
                % labeledImage = bwlabel(segIm);
                tmp = label2rgb(segIm+1, 'jet', 'k', 'shuffle') ;
                imshow(tmp)
            end

            % polygons and centroids in 2d 
            props = regionprops(segIm, 'centroid') ;
            seg2d.cdat = struct() ;
            seg2d.cdat.centroid = zeros(max(segIm(:)), 2) ;
            for cid = 1:nCells
                seg2d.cdat.centroid(cid, :) = props(cid).Centroid ;
            end
            polygons = polygonsFromSegmentation(segIm) ;        
            seg2d.cdat.polygons = polygons ;
            centroids = seg2d.cdat.centroid ;

            % Unravel polygons into cellIDs and cellVtx2D
            cellVtx2D = nan(100*length(polygons), 2) ;
            cellIDs = cell(nCells, 1) ;
            dmyk = 1 ;
            for cid = 1:nCells
                poly = polygons{cid} ;
                nAdd = size(poly, 1) ;
                cellVtx2D(dmyk:dmyk+nAdd-1, :) = poly ;
                cellIDs{cid} = dmyk:(dmyk+nAdd-1) ;
                dmyk = dmyk + nAdd ;
            end
            cellVtx2D = cellVtx2D(1:dmyk-1, :) ;
            seg2d.vdat = struct() ;
            seg2d.vdat.v = cellVtx2D ;
            seg2d.cdat.vertexCellIDs = cellIDs ;

            % save processed segmentation2d
            save(seg2dFn, 'seg2d', 'segIm')
        end


        %% PUSH INTO 3D
        seg3dFn = fullfile(QS.dir.tracking, subdir, sprintf('seg3d_%06d.mat', tp)) ;
        if ~exist(seg3dFn, 'file') || overwrite 
            % tile annular cutmesh for triple covering 
            [ TF, TV2D, TV3D, TVN3D ] = tileAnnularCutMesh(cutMesh, [1,1]) ;

            % Allow tracks to disappear for this frame
            disp('Subset of the tracks that are active')
            keepBinary = ~isnan(centroids(:, 1)) ;
            keep = find(keepBinary) ;
            % check that no cells have missing x position but not missing y pos
            assert(all(keepBinary == ~isnan(centroids(:, 2)))) ;

            % Filter
            nCells = max(segIm(:)) ;
            centroids = centroids(keep, :) ;
            cellIDs = cellIDs(keep) ;

            % measurements in 3D of tracked cells
            disp('Perform measurements in 3D')
            cellVtxU = QS.XY2uv(im, cellVtx2D(:, [2, 1]), imDoubleCovered, 1, 1) ;
            cntrds_uv = QS.XY2uv( im, centroids, imDoubleCovered, 1, 1) ;

            % embed with measurements
            [c3d, cellCntrd3d, areas, perim, moment1, ang1, ...
                moment2, ang2, moinertia, cellQ2d, cellMeshFaces, vertexMeshFaces, cellNormals] = ...
                polygonNetwork3dMeasurements(TF, TV3D, TV2D, cellVtxU, cellIDs, cntrds_uv) ;

            % Uncertainty in areas is perim * resolution
            disp('Error estimation/propagation')
            TVXY = QS.uv2XY(im, TV2D, doubleCovered, 1, 1) ;
            gg = constructFundamentalForms(TF, TV3D, TVXY) ;
            cntrdSqrtDetg = zeros(length(keep), 1) ;
            perim2d = cntrdSqrtDetg ;
            for cid = 1:length(keep)
                faceID = cellMeshFaces(cid) ;
                cntrdSqrtDetg(cid) = sqrt(det(gg{faceID})) ;
                tmp = regionprops(segIm == keep(cid), 'perimeter') ;
                perim2d(cid) = tmp.Perimeter ;
            end
            unc_areas = cntrdSqrtDetg .* perim2d ;

            save(seg3dFn, 'keep', 'keepBinary', 'cellIDs', ...
                 'centroids', 'cntrds_uv', 'cellCntrd3d', ...
                'c3d', 'cellVtxU', 'areas', 'unc_areas', ...
                'moment1', 'ang1', 'moment2', 'moinertia', 'cellMeshFaces')
        else
            load(seg3dFn, 'keep', 'keepBinary', 'cellIDs', ...
                 'centroids', 'cntrds_uv', 'cellCntrd3d', ...
                'c3d', 'cellVtxU', 'areas', 'unc_areas', ...
                'moment1', 'ang1', 'moment2', 'moinertia', 'cellMeshFaces')
        end
        
        % store areas
        disp('Storing areas in large array')
        allAreas(dmyii, keepBinary) = areas ;
        
    else
        seg2dFn = trackOutfn ;
        outFigFn = fullfile(QS.dir.tracking, subdir, 'embeddingFigure', sprintf('trackSeg3d_%06d.png', tp)) ;
        outImFn = fullfile(QS.dir.tracking, subdir, 'embeddingTexture', sprintf('trackSeg3d_%06d.png', tp)) ;
        if ~exist(fullfile(QS.dir.tracking, subdir, 'embeddingFigure'), 'dir')
            mkdir(fullfile(QS.dir.tracking, subdir, 'embeddingFigure'))
        end
        if ~exist(fullfile(QS.dir.tracking, subdir, 'embeddingTexture'), 'dir')
            mkdir(fullfile(QS.dir.tracking, subdir, 'embeddingTexture'))
        end
        if exist(seg2dFn, 'file') 
            load(seg2dFn, 'tracks')       
            cellIDs = 1:length(tracks) ;
            centroids = nan(length(tracks), 2) ;
            for trackID = 1:length(tracks)
                centroids(trackID, :) = tracks{trackID}(tidx, 1:2) ;     
            end
            
            if selectPairs
                centroids = centroids(pairIDs(:), :) ;
            end
        else
            error('Could not find tracks on disk')
        end
        
         %% PUSH INTO 3D
        seg3dFn = fullfile(QS.dir.tracking, subdir, sprintf('seg3d_%06d.mat', tp)) ;
        if ~exist(seg3dFn, 'file') || overwrite 
            % tile annular cutmesh for triple covering 
            [ TF, TV2D, TV3D, TVN3D ] = tileAnnularCutMesh(cutMesh, [1,1]) ;

            % Allow tracks to disappear for this frame
            disp('Subset of the tracks that are active')
            keepBinary = ~isnan(centroids(:, 1)) ;
            keep = find(keepBinary) ;
            badidx = find(~keepBinary) ;
            % check that no cells have missing x position but not missing y pos
            assert(all(keepBinary == ~isnan(centroids(:, 2)))) ;

            % Filter
            nCells = max(cellIDs(:)) ;
            centroids = centroids(keep, :) ;
            cellIDs = cellIDs(keep) ;

            % measurements in 3D of tracked cells
            disp('Push tracks to 3D')
            cntrds_uv = QS.XY2uv( im, centroids, doubleCovered, 1, 1) ;

            % embed with measurements
            if strcmpi(coordSys, 'spsme') || strcmpi(coordSys, 'spsm')
                coordSysTemp = 'spsmrs' ;
            else
                error('handle this coordSys here')
            end
            [cellCntrd3d, cellMeshFaces] = QS.uv2APDV(cntrds_uv, ...
                coordSysTemp, 1.0, 1.0) ;
            
            save(seg3dFn, 'keep', 'keepBinary', 'cellIDs', ...
                 'centroids', 'cntrds_uv', 'cellCntrd3d', 'cellMeshFaces')
        else
            load(seg3dFn, 'keep', 'keepBinary', 'cellIDs', ...
                 'centroids', 'cntrds_uv', 'cellCntrd3d', 'cellMeshFaces')
        end

    end
    
    %% Geodesic and ROI handling
    if drawGeodesics
        if strcmpi(geodesicStyle, 'geodesic')
            error('I think the following does not quite work -- check that pt is on mesh')
            geoP = surfaceGeodesicPairs(glueMesh.f, glueMesh.v, [1,2], ...
                cellCntrd3d) ;
        elseif strcmpi(geodesicStyle, 'pullback')
            % curve connecting the two cells in 3D based on sphi
            disp('Push connecting curve to 3D')

            % embed with measurements
            if strcmpi(coordSys, 'spsme') || strcmpi(coordSys, 'spsm')
                coordSysTemp = 'spsmrs' ;
            else
                error('handle this coordSys here')
            end
            lspace = linspace(0, 1, 101) ;
            for pidx = 1:size(pairIDs, 1)
                try
                    pbCurvX = centroids(pidx*2-1, 1) * lspace + centroids(2*pidx, 1) * (1-lspace) ;
                    pbCurvY = centroids(pidx*2-1, 2) * lspace + centroids(2*pidx, 2) * (1-lspace) ;
                catch
                    error(['Indices do not match centroid size -- there are nans:' num2str(badidx)])
                end
                pbCurv_uv = QS.XY2uv( im, [pbCurvX(:), pbCurvY(:)], doubleCovered, 1, 1) ;

                geoP{pidx} = QS.uv2APDV(pbCurv_uv, ...
                    coordSysTemp, 1.0, 1.0) ;
            end
        end
    end

    %% Compute the rectangular patch that is formed at t0 for advection
    % to other timepoints
    if drawPairRectangle
            
        % if this is the first time that we are finding cells
        if dmyii == 1
            % Load the pathlines for the Rectangular ROI or compute them
            if ~exist(roiPathlineRectangleFn, 'file') || overwriteROI 

                seg3dFn0 = fullfile(QS.dir.tracking, subdir,...
                    sprintf('seg3d_%06d.mat', t0forPairs)) ;
                tmp = load(seg3dFn0, 'centroids') ;
                vtxC = tmp.centroids ;
                minU = min(vtxC(:, 1)) ;
                minV = min(vtxC(:, 2)) ;
                maxU = max(vtxC(:, 1)) ;
                maxV = max(vtxC(:, 2)) ;
                corners = [minU, minV; minU, maxV; maxU, maxV; maxU, minV];
                
                if ~contains(coordSys, 'e') && contains(QS.piv.imCoords, 'e')
                    disp('adding 1/4 of image for extended/nonextended mismatch')
                    corners = corners + [0, size(im, 2)*0.25 ] ;
                end
                XY0ROI = corners ;
                pbOpts.t0 = t0forPairs ;
                [roi2d, roi3d] = QS.samplePullbackPathlines(XY0ROI, pbOpts) ;

                % Save the ROI pullback pathlines
                save(roiPathlineRectangleFn, 'roi2d', 'roi3d')
            else
                load(roiPathlineRectangleFn, 'roi2d', 'roi3d') ;
            end
        end
    end

    %% Plot in 3d
    if ~exist(outImFn, 'file') || overwrite || true 
        close all
        fig = figure('Visible', 'Off', 'units', 'centimeters', ...
            'position', [0,0,xwidth,ywidth]) ;

        % Draw cell faces or NUCLEI
        if strcmpi(method, 'segmentation')
            % Draw cell faces
            opts = struct() ;
            opts.centroids = cellCntrd3d ;
            opts.vertexBasedTiling = true ;
            [ff, vv, faceMemberIDs] = ...
                polygonsToPatchTriangles3D(c3d, cellIDs, opts) ;
            patch('Faces',ff,'Vertices',vv,...
                'FaceVertexCData',areas(faceMemberIDs),'FaceColor','flat', ...
                'Edgecolor', 'none', 'linewidth', 0.01);
            hold on;

            % Draw contours of cells
            for cid = 1:length(cellIDs)
                poly = c3d(cellIDs{cid}, :) ;
                plot3(poly(:, 1), poly(:, 2), poly(:, 3), 'k-', ...
                    'linewidth', lw)
            end
        else
            % Draw nuclei
            plot3(cellCntrd3d(:, 1), cellCntrd3d(:, 2), cellCntrd3d(:, 3), ....
                '.w')
                hold on
        end

        % draw Geodesic paths between pairs
        if drawGeodesics
            for pidx = 1:size(pairIDs, 1)
                plot3(geoP{pidx}(:, 1), geoP{pidx}(:, 2),...
                    geoP{pidx}(:, 3), 'color', pairColors{pidx}, ...
                    'linewidth', lw_geoP)
            end
        end
        
        % draw patch for ROI to show deformation of rectangle
        if drawPairRectangle

            corners = squeeze(roi2d(tidx, :, :)) ;
            if ~contains(coordSys, 'e') && contains(QS.piv.imCoords, 'e')
                disp('adding 1/4 of image for extended/nonextended mismatch')
                corners = corners - [0, size(im, 2)*0.25 ] ;
            end
            lspace = linspace(0, 1, 101) ;
            X0ROI = [] ;
            Y0ROI = [] ;
            for pathPart = 1:4
                nextID = pathPart + 1 ;
                if nextID > 4
                    nextID = 1 ;
                end
                X0ROI = [X0ROI, corners(pathPart, 1) * (1-lspace) + ...
                    corners(nextID, 1) * lspace ] ;
                Y0ROI = [Y0ROI, corners(pathPart, 2) * (1-lspace) + ...
                    corners(nextID, 2) * lspace ];
            end
            XY0ROI = [X0ROI(:), Y0ROI(:)] ;
            uvroi = QS.XY2uv(im, XY0ROI, doubleCovered, 1, 1) ;
            
            % check it
            tmp = figure; 
            plot(XY0ROI(:, 1), XY0ROI(:, 2), '.')
            close(tmp) 
            
            if strcmpi(coordSys, 'spsme') || strcmpi(coordSys, 'spsm')
                coordSysTemp = 'spsmrs' ;
            else
                error('handle this coordSys here')
            end
            QS.setTime(tp)
            roi3d = QS.uv2APDV( uvroi, coordSysTemp, 1, 1) ;
            plot3(roi3d(:, 1), roi3d(:, 2), ...
                roi3d(:, 3), 'color', roiColor, 'linewidth', lw_geoP) ;
        end
        
        axis equal
        xlim(xyzlims(1, :))
        ylim(xyzlims(2, :))
        zlim(xyzlims(3, :))

        colormap inferno
        caxis([0, 80])
        cb = colorbar() ;
        ylabel(cb, ['area, ' QS.spaceUnits '$^{2}$]'], ...
             'interpreter', 'latex')

        timeinterval = QS.timeInterval ;
        timeunits = QS.timeUnits ;
        t0 = QS.t0set() ;
        titlestr = ['$t = $' num2str(tp*timeinterval-t0) ' ' timeunits] ;
        if blackFigure
            title(titlestr, 'Interpreter', 'Latex', 'Color', 'white') 
        else
            title(titlestr, 'Interpreter', 'Latex', 'Color', 'k') 
        end
        xlabel('AP position [$\mu$m]', 'Interpreter', 'Latex')
        ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex')
        zlabel('DV position [$\mu$m]', 'Interpreter', 'Latex')

        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
        view(viewAngles) ;

        % Plot mesh for occlusion
        trisurf(triangulation(glueMesh.f, glueMesh.v), 'faceColor', 'k') ;


        % Make background black & Make tick labels white
        if blackFigure
            set(gca, 'color', 'k', 'xcol', 'w', 'ycol', 'w', 'zcol', 'w')
            set(gcf, 'InvertHardCopy', 'off');
            set(gcf, 'Color', 'k')
            set(gcf, 'color', 'k')
        else
            set(gcf, 'color', 'w')
            set(gca, 'color', 'k')
        end

        % Use export_fig instead, from plotting/export_fig/
        % saveas(fig, fullfile(figvdir, fnv))
        export_fig(outFigFn, '-nocrop', '-r200')

        pim = imread(outFigFn) ;
        size(pim)

        % Texturepatch counterpart  
        opts = metadat ;
        opts.timePoints = timePoints(dmyii) ;
        opts.plot_perspective = true ;
        opts.plot_dorsal = false ;
        opts.plot_ventral = false ;
        opts.plot_left = false ;
        opts.plot_right = false ;
        opts.blackFigure = false ;
        opts.makeColorbar = true ;
        QS.plotSeriesOnSurfaceTexturePatch(opts, Options)


        figoutdir = fullfile(QS.dir.texturePatchIm , opts.subdir) ;
        figPerspDir = fullfile(figoutdir, 'perspective') ;
        timFn = fullfile(figPerspDir, sprintf('patch_persp_%06d.png', tp)) ;
        tim = imread(timFn) ;
        size(tim)

        % Add images together -- area to add is in axisROI
        % tim on outside, pim right of colorbarX
        im2 = tim ;
        axisROIFn = fullfile(QS.dir.tracking, subdir, 'embeddingROI_for_overlay.mat') ;
        if ~exist(axisROIFn, 'file') 
            if isa(pim, 'uint8')
                pim = mat2gray(pim) ;
                tim = mat2gray(tim) ;
            end
            axisROI = roipoly(0.5*pim+0.5*tim) ;
            colorbarX = 1010 ;
            save(axisROIFn, 'axisROI', 'colorbarX')    
        else
            load(axisROIFn, 'axisROI', 'colorbarX')
        end
        % sum together inside axisROI
        if strcmpi(overlayStyle, 'add')
            sim = tim + 0.9 * pim ;
        else
            simR = tim(:, :, 1) ;
            simG = tim(:, :, 2) ;
            simB = tim(:, :, 3) ;
            inds = pim(:,:,1) > 0 | pim(:,:,2) > 0 | ...
                pim(:,:,3) > 0 ;
            pimR = pim(:, :, 1) ;
            pimG = pim(:, :, 2) ;
            pimB = pim(:, :, 3) ;
            simR(inds) = pimR(inds) ;
            simG(inds) = pimG(inds) ;
            simB(inds) = pimB(inds) ;
            sim = cat(3, simR, simG, simB) ;
        end
        pR = squeeze(sim(:, :, 1)) ;
        pG = squeeze(sim(:, :, 2)) ;
        pB = squeeze(sim(:, :, 3)) ;
        im2R = squeeze(im2(:,:,1)) ;
        im2G = squeeze(im2(:,:,1)) ;
        im2B = squeeze(im2(:,:,1)) ;
        im2R(axisROI) = pR(axisROI) ;
        im2G(axisROI) = pG(axisROI) ;
        im2B(axisROI) = pB(axisROI) ;
        im2 = cat(3, im2R, im2G, im2B) ;
        % pim right of colorbarX
        im2(:, colorbarX:end, :) = pim(:, colorbarX:end, :) ;

        disp(['Saving combined image to: ' outImFn])
        imwrite(im2, outImFn)

        if preview
            set(gcf, 'visible', 'on')
            imshow(im2)
            pause(3)
        end
    end
    dmyii = dmyii + 1 ;
end

% Plot areas over time
nCells = size(allAreas, 2) ;
timestamps = (timePoints-t0)*QS.timeInterval ;
if contains(lower(QS.timeUnits), 'min') && (max(timestamps)-min(timestamps)) > 60
    timestamps = timestamps / 60 ;
    timeunits = 'hr';
else
    timeunits = QS.timeUnits ;
end

close all
figure('units', 'centimeters', 'position', [0, 0,7,7]) 
newA = nan(size(allAreas)) ;
normA = nan(size(allAreas)) ;
allAreas(allAreas == 0) = NaN ;
for cid = 1:nCells
    ma = medfilt1m(allAreas(:, cid), 2) ;
    ma = movmean(ma, 3, 'omitnan') ;
    ma(ma  == 0) = NaN;
    newA(:, cid) = ma ;
    normA(:, cid) = ma ./ ma(find(ma, 1)) ;
  
end

%% Stats on true values 
figure('units', 'centimeters', 'position', [0, 0,7,7]) ;
colors = define_colors ;
% check it
plot(timestamps, newA, 'color', 0.7 * [1,1,1])
hold on;
meansA = nanmean(newA') ;
stdsA = nanstd(newA') ;
stesA = stdsA ./ sum(~isnan(newA'), 1) ;
lowerA = meansA - stdsA ;
upperA = meansA + stdsA ;
x2 = [timestamps, fliplr(timestamps)] ;
fill(x2, [lowerA, fliplr(upperA)], colors(1, :), ...
    'facealpha', 0.3, 'edgecolor', 'none')   
hold on;
plot(timestamps, meansA, 'color', colors(1, :))
xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
ylabel(['cell areas, $A$'], 'interpreter', 'latex')
saveas(gcf, fullfile(QS.dir.tracking, subdir, 'areas_over_time_true_mean.pdf'))

errorbar(timestamps, meansA, stesA, 'color', colors(1, :))
saveas(gcf, fullfile(QS.dir.tracking, subdir, 'areas_over_time_true_err.pdf'))


%% Stats on norm changes 
close all
for ii = 1:2
    fig = figure('units', 'centimeters', 'position', [0, 0,7,7]) ;
    if ii == 1
        plot(timestamps, normA, 'color', 0.7 * [1,1,1])
    end
    hold on;
    colors = define_colors ;
    meansN = nanmean(normA') ;
    stdsN = nanstd(normA') ;
    stesN = stdsN ./ sum(~isnan(normA'), 1) ;
    lowerN = meansN - stdsN ;
    upperN = meansN + stdsN ;
    x2 = [timestamps, fliplr(timestamps)] ;
    fill(x2, [lowerN, fliplr(upperN)], colors(1, :), ...
        'facealpha', 0.3, 'edgecolor', 'none')   
    hold on;
    if ii == 1
        plot(timestamps, meansN, 'color', colors(1, :))
    else
        errorbar(timestamps, meansN, stesN, 'color', colors(1, :))
    end
    xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
    ylabel(['cell areas, $A/A_0$'], 'interpreter', 'latex')
    if ii == 1
        saveas(gcf, fullfile(QS.dir.tracking, subdir, 'areas_over_time_stats_Curves.pdf'))
    else
        saveas(gcf, fullfile(QS.dir.tracking, subdir, 'areas_over_time_stats.pdf'))
    end
end

%% Stats on mean values 
close all
fig = figure('units', 'centimeters', 'position', [0, 0,7,7]) ;
colors = define_colors ;
hold on;
x2 = [timestamps, fliplr(timestamps)] ;
fill(x2, [lowerA / meansA(1), fliplr(upperA / meansA(1))], colors(1, :), ...
    'facealpha', 0.3, 'edgecolor', 'none')   
hold on;
plot(timestamps, meansA / meansA(1), 'color', colors(1, :))
xlabel(['time [' timeunits ']'], 'interpreter', 'latex')
ylabel(['cell areas, $A/\langle A_0 \rangle$'], 'interpreter', 'latex')
saveas(gcf, fullfile(QS.dir.tracking, subdir, 'areas_over_time_divMean0_mean.pdf'))

errorbar(timestamps, meansA / meansA(1), stesA / meansA(1), 'color', colors(1, :))
saveas(gcf, fullfile(QS.dir.tracking, subdir, 'areas_over_time_divMean0_err.pdf'))


