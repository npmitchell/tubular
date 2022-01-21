function visualizeDemoTracks(QS, Options) 
% Load demo tracks (segmented or tracked objects) and plot over pullback
% data.
%
% Parameters
% ----------
% QS : 
% Options : optional struct with fields
%
% Returns
% -------
% <none>
%
% Saves to disk
% -------------
% image sequence of RGB overlays
%
% NPMitchell 2021

close all
coordSys = 'spsme' ;
timePoints = QS.xp.fileMeta.timePoints ;
floatingFrame = true ; % true for plotting regions in fixed coordinate systems
t0 = QS.t0set() ;
alphaVal = 0.3 ;
preview = false ;
imW = 20 ;  % image halfwidth in microns
imH = imW ;
overwrite = false ;
% Look for what tracks to make movies of
tracks2demo = dir(fullfile(QS.dir.tracking, 'demoTracks', 'demoTracks*.mat')) ;
outputResolution = 4 / QS.APDV.resolution ;
scaleByMetric = false ;
scaleByMetricComponents = true ;
detClim = 1.3 ;
g11g22Clim = 1.3 ;

% Unpack Options
if nargin < 2
    Options = struct() ;
end

if isfield(Options, 'coordSys')
    coordSys = Options.coordSys ;
end
if isfield(Options, 'overwrite')
    overwrite = Options.overwrite ;
end
if isfield(Options, 'floatingFrame')
    floatingFrame = Options.floatingFrame ;
end
if ~floatingFrame 
    buffer = 50 ;
else
    buffer = 0.08 ;
end
if isfield(Options, 'buffer') 
    buffer = Options.buffer ; 
end
if isfield(Options, 'preview') 
    preview = Options.preview ; 
end
if isfield(Options, 'imW')
    imW = Options.imW ;
    if ~isfield(Options, 'imH')
        imH = imW ;
    end
end
if isfield(Options, 'imH')
    imH = Options.imH ;
    if ~isfield(Options, 'imW')
        imW = imH ;
    end
end
if isfield(Options, 'scaleByMetric')
    scaleByMetric = Options.scaleByMetric ;
elseif isfield(Options, 'scaledByMetric')
    % allow typo
    scaleByMetric = Options.scaledByMetric ;
end
if isfield(Options, 'scaleByMetricComponents')
    scaleByMetricComponents = Options.scaleByMetricComponents ;
elseif isfield(Options, 'scaledByMetricComponents')
    % allow typo
    scaleByMetricComponents = Options.scaledByMetricComponents ;
end
if isfield(Options, 'outputResolution')
    outputResolution = Options.outputResolution ;
end
if isfield(Options, 'tracks2demo')
    tracks2demo = Options.tracks2demo ;
end

% Mesh is not moving for floatingFrame. For floatingFrame, we compute a
% paratmeterization of the patch containing all segmented cells in the
% track
if ~floatingFrame
    for qq = 1:length(tracks2demo)
        tname = strsplit(tracks2demo(qq).name, '.mat') ;
        tname = tname{1} ;
        outDir = fullfile(QS.dir.tracking, 'demoTracks', tname) ;
        if ~exist(outDir, 'dir')
            mkdir(outDir)
        end
        tmp = load(fullfile(tracks2demo(qq).folder, tracks2demo(qq).name)) ;
        pos = tmp.positions ;
        timePointsUsed = tmp.timePointsUsed ;
        Ncells = size(pos, 2) ;    

        colors = distinguishable_colors(Ncells, [0,0,0;1,1,1]) ;

        % Get XY limits for regions
        if strcmpi(coordSys, 'spsme')
            % Get the image size just once
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, timePointsUsed(1))) ;
            fixedImageSize = true ;
            doubleCovered = true ;
        else
            error('did not recognize coordSys')
        end

        % XYlims for plots
        xmin = Inf ;
        xmax = -Inf ;
        ymin = Inf ;
        ymax = -Inf ;
        for tidx = 1:length(timePointsUsed)

            % if no fixed frame, every timepoint has different size image
            if ~fixedImageSize
                im = imread(sprintf(QS.fullFileBase.im_sp_sme, timePointsUsed(1))) ;
            end

            % Get XYlim for this all cells at this timepoint
            for cellID = 1:Ncells
                cellPos = pos{tidx, cellID} ;
                [cx, cy] = ind2sub(size(im), cellPos) ;
                xmin = min(xmin, min(cx)) ;
                xmax = max(xmax, max(cx)) ;
                ymin = min(ymin, min(cy)) ;
                ymax = max(ymax, max(cy)) ;
            end
        end
        xmin = max(1, xmin - buffer) ;
        ymin = max(1, ymin - buffer) ;
        xmax = min(size(im, 1), xmax + buffer) ;
        ymax = min(size(im, 2), ymax + buffer) ;

        % Make video for all timepoints
        for tidx = 1:length(timePointsUsed)
            tp = timePointsUsed(tidx) ;
            outFigFn = fullfile(outDir, sprintf([tname '_%06d.png'], tp)) ;
            outImFn = fullfile(outDir, sprintf(['im_' tname '_%06d.png'], tp)) ;        
                          
            if ~exist(outFigFn, 'file') || ~exist(outImFn, 'file') || ...
                    overwrite
                % Load pullback
                if strcmpi(coordSys, 'spsme')
                    im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
                end

                segR = zeros(size(im)) ;
                segG = zeros(size(im)) ;
                segB = zeros(size(im)) ;
                segA = zeros(size(im)) ; % alpha channel
                for cellID = 1:Ncells
                    cellPos = pos{tidx, cellID} ;
                    if size(cellPos, 2) == 1
                        % bw = false(size(im)) ;
                        % bw(cellPos) = true ;
                        segR(cellPos) = colors(cellID, 1) ;
                        segG(cellPos) = colors(cellID, 2) ;
                        segB(cellPos) = colors(cellID, 3) ;
                        segA(cellPos) = 1 ;
                    else
                        error("handle XY coordinates here")
                        seg(cellPos(:, 1), cellPos(:, 2), 1) = colors(cellID, 1) ;
                        seg(cellPos(:, 1), cellPos(:, 2), 2) = colors(cellID, 2) ;
                        seg(cellPos(:, 1), cellPos(:, 2), 3) = colors(cellID, 3) ;
                    end
                end
                clf
                imRGB = cat(3, im,im,im) ;
                imRGB = imRGB(xmin:xmax, ymin:ymax, :) ;
                seg = cat(3, segR, segG, segB) ;
                seg = seg(xmin:xmax, ymin:ymax, :) ;

                % Plot it
                imshow(imRGB) ;
                hold on;
                sh = imshow(seg) ;
                set(sh, 'alphaData', alphaVal * segA(xmin:xmax, ymin:ymax))
                title(['$t = $' num2str((timePoints(tii)- t0) * QS.timeInterval ) ' ' QS.timeUnits], ...
                    'interpreter', 'latex')
                set(gcf, 'Color', 'w')

                drawnow;
                
                % Save as figure
                FF = getframe(gcf) ;
                if tidx == 1
                    Fsz = size(Fsz) ;
                else
                    
                end
                imwrite(FF.cdata, outFigFn)
                % export_fig(outFigFn, '-r150', '-nocrop')

                % Save as image
                FF = getframe(gca) ;
                imwrite(FF.cdata, outImFn)
            end
        end
    end
else
    % Floating Frame! We use just a patch of the image -- a submanifold
    % "subm" that includes all the faces involved in the track, generate an
    % as-rigid-as-possible transformation of that patch over time with 
    % free boundaries, align these to each other, plot the segmentation as
    % a patch object over the image (using the boundary of the
    % segmentation).
    % See Dillon's rigidParameterizationAlignment.m and 
    %              rigidParameterizationAlignment_test.m
    %
    for qq = 1:length(tracks2demo)
        tname = strsplit(tracks2demo(qq).name, '.mat') ;
        tname = tname{1} ;
        if scaleByMetric
            scalingStr = '_ARAPscaledByMetric' ;
        elseif scaleByMetricComponents
            scalingStr = '_ARAPscaledByMetricComponents' ;
        else
            scalingStr = '_ARAP' ;
        end
        outDir = fullfile(QS.dir.tracking, 'demoTracks' , ...
            [tname scalingStr]) ;
        if ~exist(outDir, 'dir')
            mkdir(outDir)
        end
        tmp = load(fullfile(tracks2demo(qq).folder, tracks2demo(qq).name)) ;
        pos = tmp.positions ;
        timePointsUsed = tmp.timePointsUsed ;
        Ncells = size(pos, 2) ;    

        colors = distinguishable_colors(Ncells, [0,0,0;1,1,1]) ;

        % Load pullback for fieldfaces determination
        if strcmpi(coordSys, 'spsme')
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, timePointsUsed(1))) ;
            bw = false(size(im)) ;
            fixedImageSize = true ;
            doubleCovered = true ;
        else
            error('handle this coordSys here')
        end

        % Remember the scalings used
        dXs = zeros(size(timePointsUsed)) ;
        dYs = zeros(size(timePointsUsed)) ;
        imsizes = zeros(size(timePointsUsed, 1), 2) ;
        
        % First pass quickly through timepoints, then fill in interstiftial
        % frames
        tidx2do = [1, length(timePointsUsed)] ; 
        tidx2do = [tidx2do, setdiff(1:30:length(timePointsUsed), tidx2do)]   ;
        tidx2do = [tidx2do, setdiff(1:10:length(timePointsUsed), tidx2do)] ;
        tidx2do = [tidx2do, setdiff(1:length(timePointsUsed), tidx2do)] ;
        for tidx = tidx2do
            tp = timePointsUsed(tidx) ;

            outTextImDir = fullfile(outDir, 'texturePatches') ;
            if ~exist(outTextImDir, 'dir')
                mkdir(outTextImDir)
            end
            outTextFullImDir = fullfile(outDir, 'texturePatchesFullSize') ;
            if ~exist(outTextFullImDir, 'dir')
                mkdir(outTextFullImDir)
            end
            outOverlayImDir = fullfile(outDir, 'segOverlayTexturePatches') ;
            if ~exist(outOverlayImDir, 'dir')
                mkdir(outOverlayImDir)
            end
            outSegImDir = fullfile(outDir, 'segmentationPatches') ;
            if ~exist(outSegImDir, 'dir')
                mkdir(outSegImDir)
            end
            scaleBarFn = fullfile(outDir, ...
                sprintf('scaleBar_%06d.txt', tp)) ;
            
            segImOutFn = fullfile(outSegImDir, ...
                sprintf('seg_%06d.png', tp)) ;
            overlayImOutFn = fullfile(outOverlayImDir, ...
                sprintf('overlay_%06d.png', tp)) ;
            textureImOutFn = fullfile(outTextImDir, ...
                sprintf('rigidTextureMap_%06d.png', tp)) ;
            textureFullImOutFn = fullfile(outTextFullImDir, ...
                sprintf('rigidTextureMapFullSize_%06d.png', tp)) ;
            textureFullBoxImOutFn = fullfile(outTextFullImDir, ...
                sprintf('rigidTextureMapFullSizeWithBox_%06d.png', tp)) ;
            outdatFn = fullfile(outDir, ...
                sprintf('segmentationFloatingFrame_%06d.mat', tp)) ;
            
            if ~exist(segImOutFn, 'file') || ~exist(outdatFn, 'file') || ...
                    ~exist(textureImOutFn, 'file') || overwrite
                QS.setTime(timePointsUsed(tidx)) ;

                cutMesh = QS.getCurrentSPCutMeshSm() ;
                glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
                cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
                glueMesh.u(:, 1) = glueMesh.u(:, 1) / max(glueMesh.u(:, 1)) ;
                bc = barycenter(cutMesh.u, cutMesh.f) ;

                % if no fixed image size, every timepoint has different size
                if ~fixedImageSize
                    im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
                    bw = false(size(im)) ;
                end

                xmin = Inf ;
                xmax = -Inf ;
                ymin = Inf ;
                ymax = -Inf ;
                % Get XYlim for  all cells at this timepoint
                cellU = cell(Ncells, 1) ;
                segBU = cell(Ncells, 1) ;
                segCOMU_eachCell = zeros(Ncells, 2) ;
                for cellID = 1:Ncells
                    cellPos = pos{tidx, cellID} ;

                    % Get bounding box for all cells
                    [cx, cy] = ind2sub(size(im), cellPos) ;
                    uv = QS.XY2uv(im, [cy, cx], doubleCovered, 1, 1) ;
                    xmin = min(xmin, min(uv(:, 1))) ;
                    xmax = max(xmax, max(uv(:, 1))) ;
                    ymin = min(ymin, min(uv(:, 2))) ;
                    ymax = max(ymax, max(uv(:, 2))) ;

                    % Get boundary of cell segmentation for this cell
                    bwI = bw ;
                    bwI(cellPos) = true ;
                    segbnds = bwboundaries(bwI) ;
                    assert(length(segbnds) == 1)
                    bxy = segbnds{1} ;
                    cellU{cellID} = uv ;
                    segBU{cellID} = QS.XY2uv(im, [bxy(:, 2), bxy(:, 1)], ...
                        doubleCovered, 1, 1) ;
                    segCOMU_eachCell(cellID, :) = [mean(uv(:, 1)), mean(uv(:, 2))] ;
                end
                xmin = max(0, xmin - buffer) ;
                ymin = max(0, ymin - buffer) ;
                xmax = min(1, xmax + buffer) ;
                ymax = min(1, ymax + buffer) ;

                % Check it
                % disp('segBU = ')
                % segBU

                % Incude all faces in [xmin,xmax], [ymin,ymax]
                faces2rm_initial = bc(:, 1) < xmin | bc(:, 1) > xmax |...
                    bc(:, 2) < ymin | bc(:, 2) > ymax ;

                % Figure out which faces to remove, including ears
                [ submV, submF, facesKept ] = removeMeshFaces(glueMesh.v, glueMesh.f,...
                    faces2rm_initial ) ;
                [submVnew, submFnew, vtxKept, facesKept2] = removeMeshEars(submV, submF) ;
                % [newFace, newFaceIdx, newVertex, newVertexIdx] = clip_mesh_mod(submF, submV) ;
                others2rm = setdiff(1:size(submF, 1), facesKept2) ;
                others2rm = facesKept(others2rm) ;
                faces2rm = unique([find(faces2rm_initial); others2rm]) ;
                assert(length(faces2rm) == length(find(faces2rm_initial)) + length(others2rm))

                % Now create ear-free submesh
                [ submTV, submTF ] = removeMeshFaces(cutMesh.v, cutMesh.f,...
                    faces2rm ) ;
                [ submV, submF ] = removeMeshFaces(glueMesh.v, glueMesh.f,...
                    faces2rm ) ;
                [ submU, ~ ] = removeMeshFaces(glueMesh.u, glueMesh.f,...
                    faces2rm ) ;
                subm = struct('f', submF, 'v', submV, 'u', submU, ...
                    'TV', submTV, 'submTF', submTF) ;

                % Check the submesh
                % trisurf(triangulation(subm.f, subm.v))


                % SURFACE PARAMETERIZATION
                param = struct();
                param.borderType = 2; % FIXED = 1, FREE = 2
                param.fixedType = 1; % ARC_LENGTH = 1, UNIFORM = 2
                param.fixedShape = 1; % CIRCLE = 1, SQUARE = 2
                % FIXED: BARYCENTRIC = 1, AUTHALIC = 2, CONFORMAL = 3, MEAN = 4
                % FREE: LSCM = 1, ARAP = 2
                param.paramMethod =  2;
                % param.corners = [];
                % param.fixedPoints = [];
                fprintf('Generating surface parameterization... ');
                [ V2D, outFaces ] = surface_parameterization( subm.f, subm.v, param );

                % Check for spuriously added points -- sometimes a point at the
                % origin is added at the end it seems(?). Not sure why
                if size(V2D, 1) ~= size(subm.v, 1)
                    error('Handle this case here. This is abnormal')
                    if size(V2D, 1) == size(subm.v, 1) + 1
                        assert(all(V2D( size(subm.v, 1) + 1, :) == 0))
                        [ outFaces, V2D] = ...
                            remove_vertex_from_mesh(outFaces, V2D,  size(subm.v, 1) + 1)
                    end
                end

                if param.borderType == 1
                    if param.fixedShape == 1
                        V2D = 2 .* ( V2D - 0.5 ); % Map disk to unit disk
                    end
                end

                % Register to orginal mesh
                [phi, u, v] = rigidParameterizationAlignment( ...
                    submF, submV, submU, submF, submV, V2D ) ;
                rotm = [cos(phi), sin(phi); -sin(phi), cos(phi)] ;
                V2Dr = (V2D + [u, v]) * rotm  ;
                subm.V2D = V2Dr ;
                
                % Careful checking of the order of operations
                % transV2D = [ V2D(:,1) + u, V2D(:,2) + v ];
                % newV2D = [ cos(phi).*transV2D(:,1)-sin(phi).*transV2D(:,2), ...
                %            sin(phi).*transV2D(:,1)+cos(phi).*transV2D(:,2) ];

                
                % Account for possible total dilation/contraction near the
                % tracked cells -- this is typically small, a few percent.
                if scaleByMetric || scaleByMetricComponents
                    [gg, ~] = constructFundamentalForms(submF, submV, V2Dr) ;
                    cellu = barycentricMap2d(submF, submU, V2Dr, cellU{1});
                    for cellID = 2:Ncells
                        cellu = [cellu; ...
                            barycentricMap2d(submF, submU, V2Dr, cellU{cellID})] ;
                    end
                    % comu = barycentricMap2d(submF, submU, V2Dr, segCOMU_eachCell) ;
                    cellFaces = unique(pointLocation(...
                        triangulation(submF, V2Dr), cellu)) ;
                    gF = zeros(length(cellFaces), 2, 2) ;
                    for cfid = 1:length(cellFaces)
                        gF(cfid, :, :) = gg{cellFaces(cfid)} ;
                    end
                    
                    if scaleByMetric
                        % Scale by det(g)
                        disp('scaling by metric dilation det(g)')
                        cellg = [mean(gF(:, 1, 1)), mean(gF(:, 1, 2));...
                            mean(gF(:, 2, 1)), mean(gF(:, 2, 2))] ;
                        celldetg = sqrt(det(cellg)) ;
                        V2Dr = V2Dr * sqrt(celldetg) ;
                        subm.V2Dr_scaled = V2Dr ;
                    elseif scaleByMetricComponents
                        % Scale each dim separately
                        disp('scaling by metric components g11, g22')
                        dilX = mean(gF(:, 1, 1)) ;
                        dilY = mean(gF(:, 2, 2)) ;
                        V2Dr(:, 1) = V2Dr(:, 1) * sqrt(dilX) ;
                        V2Dr(:, 2) = V2Dr(:, 2) * sqrt(dilY) ;
                        subm.V2Dr_scaled = V2Dr ;
                    else
                        error('should not end up here')
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Save g result: scaling and its variation AFTER SCALING
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [gg, ~] = constructFundamentalForms(submF, submV, V2Dr) ;
                gF = zeros(size(submF, 1), 2, 2) ;
                for cfid = 1:size(submF, 1)
                    gF(cfid, :, :) = gg{cfid} ;
                end

                % Check it
                if preview
                    clf
                    % trisurf(triangulation(submF, submV)) ;
                    hold on;
                    triplot(submF, submV(:, 1), submV(:, 2))
                    triplot(submF, V2D(:, 1), V2D(:, 2))
                    triplot(submF, V2Dr(:, 1), V2Dr(:, 2)) ;
                    axis equal
                    pause(5)

                end

                % Texture patch the submesh
                IV = QS.getCurrentData() ;
                [tim, imref] = texture_patch_to_image(submF, V2Dr, ...
                    submTF, submTV, IV) ;

                close all;
                fig = figure('visible', 'off') ;
                imshow(tim)


                % Get barycentric coordinates for segmentation boundaries
                % (patch edges)
                dX = imref.ImageSize(2) / imref.ImageExtentInWorldX ;
                dY = imref.ImageSize(1) / imref.ImageExtentInWorldY ;
                segV2D = cell(Ncells, 1) ;
                for cellID = 1:Ncells
                    segu = barycentricMap2d(submF, submU, V2Dr, segBU{cellID}) ;
                    % use imref to transform to image coordinates
                    segX = dX * (segu(:, 1) - imref.XWorldLimits(1)) ;
                    segY = dY * (segu(:, 2) - imref.YWorldLimits(1)) ;
                    segV2D{cellID} = [segX, segY] ;

                    % Check it       
                    if preview
                        if cellID == 1
                            figure(2); clf;
                            subplot(1, 3, 1)
                            triplot(triangulation(submF, V2Dr), 'color', 'k'); hold on;
                            subplot(1, 3, 2) ;
                            triplot(triangulation(submF, submU), 'color', 'k'); hold on;
                            subplot(1, 3, 3)
                            V2DrX = dX * (V2Dr(:, 1) - imref.XWorldLimits(1)) ;
                            V2DrY = dY * (V2Dr(:, 2) - imref.YWorldLimits(1)) ;
                            triplot(triangulation(submF, [V2DrX, V2DrY]), 'color', 'k'); hold on;
                        end
                        figure(2) ;
                        subplot(1, 3, 1) ;
                        plot(segu(:, 1), segu(:, 2), '.-'); hold on;
                        axis equal
                        subplot(1, 3, 2) ;
                        plot(segBU{cellID}(:, 1), segBU{cellID}(:, 2), '.-'); hold on;
                        axis equal
                        subplot(1, 3, 3)
                        plot(segX, segY, '.-'); hold on; axis equal
                        % Return attention to the first figure
                        figure(fig) ;
                    end      
                end

                comu = barycentricMap2d(submF, submU, V2Dr, segCOMU_eachCell) ;
                % use imref to transform COM of cells to image coordinates
                segCOMX = dX * (comu(:, 1) - imref.XWorldLimits(1)) ;
                segCOMY = dY * (comu(:, 2) - imref.YWorldLimits(1)) ;
                comV2D_eachCell = [segCOMX, segCOMY] ;
                comV2D = mean(comV2D_eachCell, 1) ;

                % Place the cells in the center of the image and crop  
                comV2D_p = comu ;
                XlimU = [comV2D_p(1) - imW, comV2D_p(1) + imW] ;
                YlimU = [comV2D_p(2) - imH, comV2D_p(2) + imH] ;
                Xlim = [comV2D(1)-imW * dX, comV2D(1)+imW * dX] ;
                Ylim = [comV2D(2)-imH * dY, comV2D(2)+imH * dY] ;
                xlim(Xlim)
                ylim(Ylim)
                FF = getframe(gca) ;
                textureIm = FF.cdata ;            
                textureIm = imresize(textureIm, ...
                    [2*imH* outputResolution, 2*imW*outputResolution]) ;

                % Plot the segmentation as overlay
                for cellID = 1:Ncells
                    hold on;
                    patch(segV2D{cellID}(:, 1), segV2D{cellID}(:, 2), ...
                        colors(cellID, :), ...
                        'FaceAlpha', alphaVal, 'linestyle', 'none')
                end
                
                FF = getframe(gca) ;
                overlayIm = FF.cdata ;
                overlayIm = imresize(overlayIm, ...
                    [2*imH* outputResolution, 2*imW*outputResolution]) ;

                % Save texture image only with cells in center of image
                imwrite(textureIm, textureImOutFn)

                % Save texture image with segmentation patch overlay
                imwrite(overlayIm, overlayImOutFn)

                % Save full texture image (not cropped)
                imwrite(tim, textureFullImOutFn) ;
                
                % Zoom out in figure and place box over FOV, save as
                % fullSize image withBox
                axis equal
                xlim([0, size(tim, 2)])
                ylim([0, size(tim, 1)])
                plot([Xlim(1), Xlim(1), Xlim(2), Xlim(2), Xlim(1)], ...
                    [Ylim(1), Ylim(2), Ylim(2), Ylim(1), Ylim(1)], '-')
                FF = getframe(gca) ;
                boxIm = FF.cdata ;
                imwrite(boxIm, textureFullBoxImOutFn)

                % Save image sizes for reporting
                dXs(tidx) = dX ;
                dYs(tidx) = dY ;
                imsizes(tidx, :) = size(tim) ;
                xlims(tidx, :) = Xlim ;
                ylims(tidx, :) = Ylim ;

                % Save segmentation only, for manipulation of LUT etc in post
                clf
                for cellID = 1:Ncells
                    hold on;
                    patch(segV2D{cellID}(:, 1), segV2D{cellID}(:, 2), ...
                        colors(cellID, :), ...
                        'FaceAlpha', alphaVal, 'linestyle', 'none')
                end
                axis off
                axis equal
                xlim(Xlim)
                ylim(Ylim)
                FF = getframe(gca) ;
                segIm = FF.cdata ;
                segIm = imresize(segIm, ...
                    [2*imH* outputResolution, 2*imW*outputResolution]) ;
                imwrite(segIm, segImOutFn)
                saveas(gcf, [segImOutFn(1:end-4) '.pdf'])

                % Save (approximate) scale in text file
                scaleBarX = diff(XlimU) / (2*imW* outputResolution) ;
                scaleBarY = diff(YlimU) / (2*imH* outputResolution) ;
                scaleBar = [scaleBarX, scaleBarY] ;
                spuHeader = strrep(QS.spaceUnits, '\m', 'm') ;
                header = ['scale of (cropped) image in ' spuHeader ' / pix'] ;
                write_txt_with_header(scaleBarFn, scaleBar, header)

                
                %% Save the frame info as mat
                save(outdatFn, 'segCOMU_eachCell', ...
                    'comV2D_eachCell', 'cellU', 'segV2D', ...
                    'subm', 'dX', 'dY', 'imref', 'gg', 'scaleBar')

                %% Plot scaling by metric if this was done
                close all
                fig = figure('units','centimeters','position',[0,0,18,9], ...
                    'visible', 'off') ;
                subplot(1, 3, 1)

                % Determinant panel
                dets = zeros(size(submF, 1), 1) ;
                for cfid = 1:size(submF, 1)
                    dets(cfid) = sqrt(det(gg{cfid})) ;
                end
                trisurf(submF, V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                    dets, 'edgecolor', 'none')
                view(2); grid off
                axis equal
                xlim(XlimU)
                ylim(YlimU)
                climit = detClim-1 ;
                caxis([-1,1])
                colormap(bluewhitered)
                caxis([1-climit, 1+climit])
                cb = colorbar ;
                set(gca,'XTick',[])
                set(gca,'YTick',[])
                ylabel(cb, '$\sqrt(\det g)$', 'interpreter', 'latex')

                % g11 panel
                subplot(1, 3, 2)
                trisurf(submF, V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                    gF(:, 1, 1), 'edgecolor', 'none')
                view(2) ; grid off ;
                axis equal
                xlim(XlimU)
                ylim(YlimU)
                climit1 = round(10*max(abs(gF(:, 1, 1)-1)))*0.1 ;
                climit2 = round(10*max(abs(gF(:, 2, 2)-1)))*0.1 ;
                caxis([-1, 1])
                colormap(bluewhitered)
                climit = g11g22Clim-1 ;
                caxis([1-climit, 1+climit])
                cb = colorbar ;
                set(gca,'XTick',[])
                set(gca,'YTick',[])
                ylabel(cb, '$g_{11}$', 'interpreter', 'latex')

                % g22 panel
                subplot(1, 3, 3)
                trisurf(submF, V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                    gF(:, 2, 2), 'edgecolor', 'none')
                view(2) ; grid off ;
                axis equal
                xlim(XlimU)
                ylim(YlimU)
                caxis([-1, 1])
                colormap(bluewhitered)
                caxis([1-climit, 1+climit])
                cb = colorbar ;
                set(gca,'XTick',[])
                set(gca,'YTick',[])
                ylabel(cb, '$g_{22}$', 'interpreter', 'latex')


                % Output figure
                set(gcf, 'color', 'white')
                sgtitle('map distortion in ARAP mapping', 'interpreter', 'latex')
                saveas(gcf, fullfile(outTextFullImDir, sprintf('scaling_detg11g22_%06d.pdf', tp)))

            end
        end  
        disp(['dXs = ', num2str( dXs)])
        disp(['dYs = ', num2str( dYs)])
        disp(['imXs = ', num2str( imsizes(:, 2)')])
        disp(['imYs = ', num2str( imsizes(:, 1)')])
        
    end    
end
disp('done with demoTrack visualization!')
