function visualizeSegmentationPatch(QS, Options) 
% Load demo tracks (segmented or tracked objects) and plot over pullback
% data from small patches of the surface mapped to the plane and
% registered.
%
% Todo: allow for pathline points near periodic boundary.
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : optional struct with fields
%   coordSys : str specifier
%   overwrite : bool
%   buffer : float or int or numeric
%   bufferX : float or int or numeric, trumps buffer for X dim
%   bufferY : float or int or numeric, trumps buffer for Y dim
%   preview : preview (default=false)
%   scaleByMetric : bool (default=false)
%       rescale parameterization to have det(g) approx 1 on either:
%           - faces queried by pts (if both imW and imH are not supplied)
%           - in window com(pts) +/- [imW*0.5, imH*0.5] 
%   scaleByMetricComponents : bool (default=true)
%       rescale parameterization to approximate isothermal coordinates
%   imW : int or numeric (default=25)
%       halfwidth of output image in QS.spaceUnits (approx)
%       (in units of quasi-embedding space --> ie registered/flattened
%       embedding submesh units, should be close to units of QS.spaceUnits if
%       submesh is not too strongly curved in embedding space).
%   imH : int or numeric (default=25)
%       halfheight of output image in QS.spaceUnits (approx)
%       (in units of quasi-embedding space --> ie registered/flattened
%       embedding submesh units, should be close to units of QS.spaceUnits if
%       submesh is not too strongly curved in embedding space).
%
%
% Returns
% -------
% <none>
%
% Saves to disk
% -------------
% image sequence of RGB overlays
%
% See also
% --------
% computeLocalSurfacePatch(QS, pts, options)
% visualizeDemoTracks(QS, options)
%
% NPMitchell 2021

close all
coordSys = 'spsme' ;
allTimePoints = QS.xp.fileMeta.timePoints ;
timePoints = allTimePoints ;
floatingFrame = true ; % true for plotting regions in fixed coordinate systems
t0 = QS.t0set() ;
t0Pathlines = t0 ;
alphaVal = 0.3 ;
preview = false ;
imW = 25 ;  % image halfwidth in microns
imH = imW ;
overwrite = false ;
% Look for what tracks to make movies of
tracks2demo = dir(fullfile(QS.dir.tracking, 'demoTracks', 'demoTracks*.mat')) ;
outputResolution = 4 / QS.APDV.resolution ;
scaleByMetric = false ;
scaleByMetricComponents = true ;
bufferX = 0.12 ;
bufferY = 0.08 ;
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
if isfield(Options, 'timePoints')
    timePoints = Options.timePoints ;
end
if isfield(Options, 't0Pathlines')
    t0Pathlines = Options.t0Pathlines ;
end
if isfield(Options, 'floatingFrame')
    floatingFrame = Options.floatingFrame ;
end
if isfield(Options, 'buffer') 
    bufferX = Options.buffer ; 
    bufferY = Options.buffer ; 
end
if isfield(Options, 'bufferX') 
    bufferX = Options.bufferX ; 
end
if isfield(Options, 'bufferY') 
    bufferY = Options.bufferY ; 
end
if isfield(Options, 'preview') 
    preview = Options.preview ; 
end
if isfield(Options, 'detClim') 
    detClim = Options.detClim ; 
end
if isfield(Options, 'g11g22Clim') 
    g11g22Clim = Options.g11g22Clim ; 
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
end

if isfield(Options, 'outputResolution')
    outputResolution = Options.outputResolution ;
end

% For each selected patch location, create both fixedFrame images AND
% floatingFrame registered patches.
if isfield(Options, 'demoPatchName')
    dpname = Options.demoPatchName ;
else
    outDirs = dir(fullfile(QS.dir.segmentation, 'demoPatch*')) ;
    nDirsThatExistAleady = 0 ;
    for qq = 1:length(outDirs)
        if exist(fullfile(outDirs.folder, outDirs.name), 'dir')
            nDirsThatExistAleady = nDirsThatExistAleady + 1 ;
        end
    end
    dpname = sprintf('demoPatch%03d', nDirsThatExistAleady) ;
end

% Now make the output directory
outDir = fullfile(QS.dir.segmentation, dpname) ;
if ~exist(outDir, 'dir')
    mkdir(outDir)
end
outFixedDir = fullfile(outDir, 'fixedFrame') ;
if ~exist(outFixedDir, 'dir')
    mkdir(outFixedDir)
end

% Get XY limits for regions
if strcmpi(coordSys, 'spsme')
    % Get the image size just once
    im = imread(sprintf(QS.fullFileBase.im_sp_sme, timePoints(1))) ;
    shiftY = size(im, 1) * 0.5 ;
    fixedImageSize = true ;
    doubleCovered = true ;
else
    error('did not recognize coordSys')
end


%% Find tissue patch position from pathline emanating from original point
startingPtFn = fullfile(outDir, 'startingPointPathlines.mat') ;
if exist(startingPtFn, 'file')
    disp(['Loading pathline to follow for tissue patch: ', startingPtFn])
    load(startingPtFn, 'XY0', 'XYpathline', 'v3dpathline') ;
else
    disp('choose starting point for folling tissue patch segmentation')
    im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, t0Pathlines)) ;
    imshow(im0) ;
    XY0 = ginput(1) ;
    opts = struct() ;
    opts.coordSys = coordSys ;
    opts.doubleCovered = doubleCovered ;
    opts.t0 = t0Pathlines ;
    [XYpathline, v3dpathline] = QS.samplePullbackPathlines(XY0, opts) ;

    save(startingPtFn, 'XY0', 'XYpathline', 'v3dpathline') ;
end


%% FIXED FRAME VIDEO 
% Get XY limits for fixed frame
Xbuffer = bufferX * size(im, 1) ;
Ybuffer = bufferY * size(im, 1) ;

% Make video for all timepoints
for tii = 1:length(timePoints)
    tp = timePoints(tii) ;
    tidx = QS.xp.tIdx(tp) ;
    
    outFigFn = fullfile(outFixedDir, sprintf([dpname '_%06d.png'], tp)) ;
    outImFn = fullfile(outFixedDir, sprintf(['im_' dpname '_%06d.png'], tp)) ;        

    if ~exist(outFigFn, 'file') || ~exist(outImFn, 'file') || ...
            overwrite
        % Load pullback
        if strcmpi(coordSys, 'spsme')
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
        end

        % Load segmentation of this timepoint
        seg = load(sprintf(QS.fullFileBase.segmentation2dCorrected, coordSys, tp)) ;
        colors = distinguishable_colors(max(seg.segIm(:))-1) ;
        % colors = shuffleArray(colors) ;
        overlayIm = seg.segIm - 1; 
        if doubleCovered
            overlayIm(round(size(im,1)*0.75+1):end, :) = ...
                max(...
                overlayIm(round(size(im,1)*0.75+1):end, :), ...
                overlayIm(round(size(im,1)*0.25+1):round(size(im,1)*0.5), :)) ;
            overlayIm(1:round(size(im,1)*0.25), :) = ...
                max( ...
                overlayIm(1:round(size(im,1)*0.25), :), ...
                overlayIm(round(size(im,1)*0.5+1):1500, :)) ;
        end
        segRGB = label2rgb(overlayIm, colors) ;
        alphaIm = alphaVal * (overlayIm  > 0) ;
        
        close all
        fig = figure('units','pixels','position',[0,0,640,512]) ;
        xmin = min(XYpathline(tidx, 1)) - Xbuffer ;
        xmax = max(XYpathline(tidx, 1)) + Xbuffer ;
        ymin = min(XYpathline(tidx, 2)) - Ybuffer ;
        ymax = max(XYpathline(tidx, 2)) + Ybuffer ;
        xminIm = max(1, round(xmin - 5)) ;
        yminIm = max(1, round(ymin - 5)) ;
        xmaxIm = min(size(im, 2), round(xmax + 5)) ;
        ymaxIm = min(size(im, 2), round(ymax + 5)) ;
        imc = im(yminIm:ymaxIm, xminIm:xmaxIm) ;
        imRGB_c = cat(3, imc,imc,imc) ;
        segRGB_c = segRGB(yminIm:ymaxIm, xminIm:xmaxIm, :) ;
        alphaIm_c = alphaIm(yminIm:ymaxIm, xminIm:xmaxIm) ;

        % Plot it
        imshow(imRGB_c) ;
        hold on;
        sh = imshow(segRGB_c) ;
        set(sh, 'alphaData', alphaIm_c)
        title(['$t = $' num2str((timePoints(tii)- t0) * QS.timeInterval ) ' ' QS.timeUnits], ...
            'interpreter', 'latex')
        set(gcf, 'Color', 'w')

        % Set XYlim
        xlim([xmin - xminIm, xmin - xminIm + (xmax-xmin)])
        ylim([ymin - yminIm, ymin - yminIm + (ymax-ymin)])
        drawnow;

        % Save as figure
        FF = getframe(gcf) ;
        figIm = FF.cdata ;
        if tii == 1
            Fsz = size(figIm) ;
        else
            if ~exist('Fsz', 'var')
                outFigFn0 = fullfile(outFixedDir, sprintf([dpname '_%06d.png'], timePoints(1))) ;
                Fsz = size(imread(outFigFn0)) ;
            end
            if any(size(figIm) ~= Fsz)
                figIm = imresize(figIm, [Fsz(1), Fsz(2)]) ;
            end
        end
        imwrite(figIm, outFigFn)
        % export_fig(outFigFn, '-r150', '-nocrop')

        % Save as image
        FF = getframe(gca) ;
        axIm = FF.cdata ;
        if tii == 1
            Asz = size(FF.cdata) ;
        else
            if ~exist('Asz', 'var')
                outImFn0 = fullfile(outFixedDir, sprintf(['im_' dpname '_%06d.png'], timePoints(1))) ;
                Asz = size(imread(outImFn0)) ;
            end
            if any(size(axIm) ~= Asz)
                axIm = imresize(axIm, [Asz(1), Asz(2)]) ;
            end
        end
        imwrite(axIm, outImFn)
    end
end


%% Now Floating Frame! 
% We use just a patch of the image -- a submanifold
% "subm" that includes all the faces involved in the track, generate an
% as-rigid-as-possible transformation of that patch over time with 
% free boundaries, align these to each other, plot the segmentation as
% a patch object over the image (using the boundary of the
% segmentation).
% See Dillon's rigidParameterizationAlignment.m and 
%              rigidParameterizationAlignment_test.m
%
if scaleByMetricComponents
    scalingStr = '_ARAPscaledByMetricComponents' ;
elseif scaleByMetric
    scalingStr = '_ARAPscaledByMetric' ;
else
    scalingStr = '_ARAP' ;
end
if ~exist(outDir, 'dir')
    mkdir(outDir)
end


% Remember the scalings used
dXs = zeros(size(timePoints)) ;
dYs = zeros(size(timePoints)) ;
imsizes = zeros(size(timePoints, 1), 2) ;

% First pass quickly through timepoints, then fill in interstiftial
% frames
tidx2do = [1, length(timePoints)] ; 
tidx2do = [tidx2do, setdiff(1:30:length(timePoints), tidx2do)]   ;
tidx2do = [tidx2do, setdiff(1:10:length(timePoints), tidx2do)] ;
tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do)] ;
for tii = tidx2do
    tp = timePoints(tii) ;
    tidx = QS.xp.tIdx(tp) ;

    outTextImDir = fullfile(outDir, ['texturePatches' scalingStr]) ;
    if ~exist(outTextImDir, 'dir')
        mkdir(outTextImDir)
    end
    outTextFullImDir = fullfile(outDir, ['texturePatchesFullSize' scalingStr]) ;
    if ~exist(outTextFullImDir, 'dir')
        mkdir(outTextFullImDir)
    end
    outSegImDir = fullfile(outDir, ['segPatches' scalingStr]) ;
    if ~exist(outSegImDir, 'dir')
        mkdir(outSegImDir)
    end
    outOverlayImDir = fullfile(outDir, ['segOverlayTexturePatches' scalingStr]) ;
    if ~exist(outOverlayImDir, 'dir')
        mkdir(outOverlayImDir)
    end
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
    scaleBarFn = fullfile(outDir, ...
        sprintf('scaleBar_%06d.txt', tp)) ;
    
    if ~exist(segImOutFn, 'file') || ~exist(outdatFn, 'file') || ...
            ~exist(textureImOutFn, 'file') || overwrite
        
        % Load this timepoint
        disp(['setting QS time to ' num2str(timePoints(tii))])
        QS.setTime(timePoints(tii)) ;

        cutMesh = QS.getCurrentSPCutMeshSm() ;
        glueMesh = QS.getCurrentSPCutMeshSmRSC() ;
        cutMesh.u(:, 1) = cutMesh.u(:, 1) / max(cutMesh.u(:, 1)) ;
        glueMesh.u(:, 1) = glueMesh.u(:, 1) / max(glueMesh.u(:, 1)) ;
        bc = barycenter(cutMesh.u, cutMesh.f) ;

        % if no fixed image size, every timepoint has different size
        if ~fixedImageSize
            im = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
        end

        % segmentation for this timepoint
        seg = load(sprintf(QS.fullFileBase.segmentation2dCorrected, coordSys, tp)) ;
        centroids = seg.seg2d.cdat.centroid ;
        pp = seg.seg2d.cdat.polygons ; 
        
        % Rough filtering
        % filter: keep only polygons that might be "somewhat nearby", here
        % assuming that the cells are not much larger than the FOV.
        nearby = max(Xbuffer, Ybuffer) * 2 ;
        % Keep cells nearby to the pathline point on some cover of the
        % manifold.
        if doubleCovered
            keep = abs(centroids(:, 1) - XYpathline(tidx, 1)) < nearby & ...
                (abs(centroids(:, 2) - XYpathline(tidx, 2)) < nearby | ...
                abs(centroids(:, 2) - XYpathline(tidx, 2) + shiftY) < nearby | ...
                abs(centroids(:, 2) - XYpathline(tidx, 2) - shiftY) < nearby) ;
        else
            keep = abs(centroids(:, 1) - XYpathline(tidx, 1)) < nearby & ...
                abs(centroids(:, 2) - XYpathline(tidx, 2)) < nearby ;
        end
        pp = pp(keep) ;
        centroids = centroids(keep, :) ;
        
        % Careful filtering
        ind2keep = false(length(pp), 1) ;
        for cellID = 1:length(pp)
            if ~isempty(pp{cellID})
                if doubleCovered
                    if any(abs(pp{cellID}(:,2) - XYpathline(tidx, 1)) < Xbuffer) && ...
                            (any(abs(pp{cellID}(:,1) - XYpathline(tidx, 2)) < Ybuffer) || ...
                            any(abs(pp{cellID}(:,1) - XYpathline(tidx, 2) - shiftY) < Ybuffer) || ...
                            any(abs(pp{cellID}(:,1) - XYpathline(tidx, 2) + shiftY) < Ybuffer)) 
                        ind2keep(cellID) = true ;
                    end
                else
                    if any(abs(pp{cellID}(:,2) - XYpathline(tidx, 1)) < Xbuffer) && ...
                            any(abs(pp{cellID}(:,1) - XYpathline(tidx, 2)) < Ybuffer)                        
                        ind2keep(cellID) = true ;
                    end
                end  
            end
        end
        % By convention, polygons are stored as (Y,X) so swap them
        pp = pp(ind2keep) ;
        segBXY = cell(size(pp)) ;
        segBU = cell(size(pp)) ;
        for cellID = 1:length(pp) 
            segBXY{cellID} = pp{cellID}(:, [2, 1]) ;
            segBU{cellID} = QS.XY2uv(im, segBXY{cellID}, doubleCovered, 1., 1.) ;
        end
        centroidsXY = centroids(ind2keep, :) ;
        Ncells = length(pp) ;
        colors = distinguishable_colors(Ncells) ;
        
        % Get XYlim for  all cells at this timepoint
        xmin = XYpathline(tidx, 1) - Xbuffer ;
        xmax = XYpathline(tidx, 1) + Xbuffer ;
        ymin = XYpathline(tidx, 2) - Ybuffer ;
        ymax = XYpathline(tidx, 2) + Ybuffer ;
        
        % Incude all faces in [xmin,xmax], [ymin,ymax]
        bcXY = QS.uv2XY(im, bc, doubleCovered, 1., 1.) ;
        faces2rm_initial = bcXY(:, 1) < xmin | bcXY(:, 1) > xmax |...
            (abs(bcXY(:, 2) - XYpathline(tidx, 2)) > Ybuffer & ...
            abs(bcXY(:, 2) - XYpathline(tidx, 2) - shiftY) > Ybuffer & ...
            abs(bcXY(:, 2) - XYpathline(tidx, 2) + shiftY) > Ybuffer);

        % Figure out which faces to remove, including ears
        [ submV, submF, facesKept ] = removeMeshFaces(glueMesh.v, glueMesh.f,...
            faces2rm_initial ) ;
        [~, ~, ~, facesKept2] = removeMeshEars(submV, submF) ;
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
        [ submUcut, submFcut ] = removeMeshFaces(cutMesh.u, cutMesh.f,...
            faces2rm ) ;

        % There are a few options for how this looks:
        % 1. submesh is within single Cover (most common, disk-like)
        %       --> no problems in sight, all triangles are oriented
        %       correctly
        % 2. submesh has one edge on another cover, with some long faces,
        % but is still a topological disk
        %   
        % 3. submesh has faces on multiple covers (topological disk)
        % 
        % 4. submesh spans the whole single Cover (topological annulus)
        bcdiff = barycenter(submUcut, submFcut) - barycenter(submU, submF) ;
        assert(~any(abs(bcdiff(:, 1)) > 1e-13))
        if ~any(abs(bcdiff(:, 2) > 1e-13 ))
            connectivityCase = 1;
        elseif eulerCharacteristic(subm) == 1
            % must be either case 2 or 3
            if all(size(submUcut) == size(submU))
                % Same number of vertices, so we haven't actually crossed
                % the branch cut, just are hugging it
                connectivityCase = 2 ;
                
                if median(submUcut(:, 2)) > 0.5 
                    submU(submU(:, 2) == 0, 2) = 1 ;
                elseif median(submUcut(:, 2)) < 0.5
                    submU(submU(:, 2) == 1, 2) = 0 ;
                else
                    error('handle how to cut the submesh here')
                end
            else
                connectivityCase = 3 ;
            end
        elseif eulerCharacteristic(subm) == 0
            connectivityCase = 4; 
        end
        % Package into struct
        subm = struct('f', submF, 'v', submV, 'u', submU, ...
            'connectivityCase', connectivityCase) ;
        if connectivityCase > 1
            subm.Vcut = submTV ;
            subm.Fcut = submTF ;
            subm.Ucut = submUcut ;
        end
        
        % Check the submesh
        if preview
            trisurf(triangulation(subm.f, subm.v))
        end
        
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
        disp('Generating surface parameterization... ');
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
        % Note: there could be some issues with the periodic BC felt by
        % submU's influence here!! todo
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
        % tracked cells -- this can be significant for large patches
        if scaleByMetric || scaleByMetricComponents
            % Compute metric in a small window near the origin
            [gg, ~] = constructFundamentalForms(submF, submV, V2Dr) ;
            bcV2D = barycenter(V2Dr, submF) ;
            Upathline = QS.XY2uv(im, XYpathline(tidx,:), doubleCovered, 1., 1.) ;
            comV2D0 = barycentricMap2d(submF, submU, V2Dr, Upathline) ;
            inBox = find(abs(bcV2D(:, 1) - comV2D0(1)) <imW*0.5 &  ...
                abs(bcV2D(:, 2) - comV2D0(2)) <imH*0.5);
            gF = zeros(size(inBox, 1), 2, 2) ;
            for cfid = 1:size(inBox, 1)
                gF(cfid, :, :) = gg{inBox(cfid)} ;
            end
            
            % Check that we are looking in the box
            if preview
                figure(3);
                triplot(triangulation(submF, V2Dr))
                hold on; 
                plot(bcV2D(inBox, 1), bcV2D(inBox, 2), 'o')
                if scaleByMetricComponents
                    title('Triangles whose metric we force to be approx isothermal')
                else
                    title('Triangles whose metric we scale to be approx det=1')
                end
                waitfor(gcf)
            end
            
            % Scale by det(g) -- either mean dilation 
            if scaleByMetricComponents
                % Scale each dim separately
                dilX = mean(gF(:, 1, 1)) ;
                dilY = mean(gF(:, 2, 2)) ;
                V2Dr(:, 1) = V2Dr(:, 1) * sqrt(dilX) ;
                V2Dr(:, 2) = V2Dr(:, 2) * sqrt(dilY) ;
                subm.V2D_scaled_g11g22 = V2Dr ;
                
                % plot preview of the scaling and its variation
                if preview
                    close all
                    subplot(1, 2, 1)
                    trisurf(submF(inBox, :), V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                        gF(:, 1, 1), 'edgecolor', 'none')
                    view(2) ; grid off ;
                    axis equal
                    cb = colorbar ;
                    ylabel(cb, '$g_{11}$', 'interpreter', 'latex')
                    
                    subplot(1, 2, 2)
                    trisurf(submF(inBox, :), V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                        gF(:, 2, 2), 'edgecolor', 'none')
                    view(2) ; grid off ;
                    axis equal
                    cb = colorbar ;
                    ylabel(cb, '$g_{22}$', 'interpreter', 'latex')
                    sgtitle('scaling in ARAP mapping', 'interpreter', 'latex')
                    set(gcf, 'color', 'white')
                    pause(5) ;
                end
            else
                cellg = [mean(gF(:, 1, 1)), mean(gF(:, 1, 2));...
                    mean(gF(:, 2, 1)), mean(gF(:, 2, 2))] ;
                celldetg = sqrt(det(cellg)) ;
                V2Dr = V2Dr * sqrt(celldetg) ;
                subm.V2D_scaled_detg = V2Dr ;

                % plot preview of the scaling and its variation
                if preview
                    clf
                    dets = zeros(size(submF, 1), 1) ;
                    for cfid = 1:size(submF, 1)
                        dets(cfid) = sqrt(det(gg{cfid})) ;
                    end
                    trisurf(submF, V2Dr(:, 1), V2Dr(:, 2), 0*V2Dr(:, 1), ...
                        dets, 'edgecolor', 'none')
                    view(2)
                    axis equal
                    cb = colorbar ;
                    ylabel(cb, '$\sqrt(\det g)$', 'interpreter', 'latex')
                    title('scaling in ARAP mapping', 'interpreter', 'latex')
                    pause(5) ;
                end
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
            subplot(1, 2, 1)
            triplot(submF, submU(:, 1), submU(:, 2))
            subplot(1, 2, 2)
            triplot(submF, V2D(:, 1), V2D(:, 2))
            triplot(submF, V2Dr(:, 1), V2Dr(:, 2)) ;
            axis equal
            pause(5)
        end

        % Texture patch the submesh
        IV = QS.getCurrentData() ;
        [tim, imref] = texture_patch_to_image(submF, V2Dr, ...
            submTF, submTV, IV) ;

        % Create figure
        close all
        fig = figure('units','pixels','position',[0,0,640,512], ...
            'visible', 'off') ;
        imshow(tim)
        
        % Get barycentric coordinates for segmentation boundaries
        % (patch edges)
        dX = imref.ImageSize(2) / imref.ImageExtentInWorldX ;
        dY = imref.ImageSize(1) / imref.ImageExtentInWorldY ;
        segV2D = cell(Ncells, 1) ;
        for cellID = 1:Ncells
            % Todo:
            % To allow periodic BCs exponentiate the pullback
            if connectivityCase == 1
                segu = barycentricMap2d(submF, submU, V2Dr, segBU{cellID}) ;
            elseif connectivityCase == 2
                % submUcut is same size as V2Dr, so correspondence is
                % rescued by exponentiating
                % zCell = exp(1j * 2*pi* segBU{cellID}(:, 2)) ;
                % magCell = log(1 + 4*segBU{cellID}(:, 1)) ;
                % cellRPhi = magCell .* [real(zCell), imag(zCell)] ;
                % zU = exp(1j * 2*pi * submUcut(:, 2)) ;
                % magU = log(1+4*submUcut(:, 1)) ;
                % suRPhi = magU .* [real(zU), imag(zU)] ;
                % segu = barycentricMap2d(submFcut, suRPhi, V2Dr, cellRPhi) ;
                
                % THE ABOVE ONLY WORKS IF WE HAVE CORRESPONDENCE BETWEEN
                % THE VERTICES V2DR AND SUBMU, WHICH WE DO NOT 
                
                % For now just cut out cells that are near the edge         
                segu = barycentricMap2d(submF, submU, V2Dr, segBU{cellID}) ;
            elseif connectivityCase == 3
                error('handle this case here: patch crosses periodic boundary')
            else
                error('handle this case here: patch is an annulus')
            end
            
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

        Upathline = QS.XY2uv(im, XYpathline(tidx,:), doubleCovered, 1., 1.) ;
        comu = barycentricMap2d(submF, submU, V2Dr, Upathline) ;
        % use imref to transform COM of cells to image coordinates
        segCOMX = dX * (comu(:, 1) - imref.XWorldLimits(1)) ;
        segCOMY = dY * (comu(:, 2) - imref.YWorldLimits(1)) ;
        comV2D = [segCOMX, segCOMY] ;

        % Place the cells in the center of the image and crop  
        comV2D_p = barycentricMap2d(submF, submU, V2Dr, Upathline) ;
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

        % Plot the segmentation as patch
        for cellID = 1:Ncells
            hold on;
            patch(segV2D{cellID}(:, 1), segV2D{cellID}(:, 2), ...
                colors(cellID, :), ...
                'FaceAlpha', alphaVal, 'linestyle', 'none')
            patch(segV2D{cellID}(:, 1), segV2D{cellID}(:, 2), ...
                'w', ...
                'FaceAlpha', 0, 'linestyle', '-', 'edgecolor', 'w', ...
                'edgealpha', 0.5*alphaVal)
        end
        axis equal
        xlim(Xlim)
        ylim(Ylim)

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
        dXs(tii) = dX ;
        dYs(tii) = dY ;
        imsizes(tii, :) = size(tim) ;
        xlims(tii, :) = Xlim ;
        ylims(tii, :) = Ylim ;
        
        % Save segmentation only, for manipulation of LUT etc in post
        clf
        for cellID = 1:Ncells
            hold on;
            patch(segV2D{cellID}(:, 1), segV2D{cellID}(:, 2), ...
                colors(cellID, :), ...
                'FaceAlpha', 1, 'linestyle', '-', 'edgecolor', 'w')
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
        scaleBarX = diff(Xlim) / (2*imW* outputResolution) ;
        scaleBarY = diff(Ylim) / (2*imH* outputResolution) ;
        scaleBar = [scaleBarX, scaleBarY] ;
        header = ['scale of (cropped) image in ' QS.spaceUnits ' / pix'] ;
        write_txt_with_header(scaleBarFn, scaleBar, header)


        % Save the frame info as mat
        save(outdatFn, 'comV2D', 'segBU', 'segV2D', 'centroidsXY', ...
            'subm', 'dX', 'dY', 'imref', 'gg', 'scaleBar')
        
        % Plot scaling by metric if this was done
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

disp('done with demoSegmentationPatch visualization!')
