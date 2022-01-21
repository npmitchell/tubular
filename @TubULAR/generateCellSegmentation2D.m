function generateCellSegmentation2D(QS, options)
% Unfinished code -- Lin working on it.
% ToDo: 
%
% Parameters
% ----------
%   QS : 
%   options : struct with fields
%       coordSys : coordinate system on which to compute the segmentation
%
% tissueAnalysisSuite fields are different than QuapSlap's. 
% Tor tissueAnalysisSuite, we have:
% vdat: 
%     nverts : neighbor list for vertices
%     ncells : 
%     vertxcoord : column of data in which vertex lives
%     vertycoord : row of data in which each vertex lives
% cdat : cell data    
%     ncells : indices of neighboring cells
%     nverts : vertices that define the cell
%     centroid.coord(1) x position of cell centroid
%     centroid.coord(2) y position of cell centroid
% bdat : 
%     nverts
%     ncells
%     pix : linear indices of the pixels associated with that bond
%
% Note on exterior calculus objects
% d0 is an e x c matrix of exterior derivatives with +1 and -1s
% at the endpts of each bond
% d1 is a v x e matrix of exterior derivatives. Upstream is +1,
% downstream is -1 when moving counterclockwise around a
% tension plaquette.
% d0 and d1 are matrices that take derivatives 

%% Parameters
overwrite = false ;
skipPolygons = false ;
overwriteImages = false ;
% timepoints to process
timePoints = QS.xp.fileMeta.timePoints ;
% how far in pixels is too far to connect two cell vertices
very_far = 250 ;
% which coordinate system to use for segmentation
coordSys = QS.currentSegmentation.coordSys ; 
% Toggle for iLastik version control -- zero for newer version
iLastikVersion = 0;
cellSize = 50 ;
strelRadius = 0 ;
gaussKernel = 1 ;
heightMinimum = 3.5 ;


%% unpack options
if nargin < 2 
    options = struct() ;
end

if isfield(options, 'overwrite') 
    overwrite = options.overwrite ;
end
if isfield(options, 'skipPolygons') 
    skipPolygons = options.skipPolygons ;
end
if isfield(options, 'overwriteImages') 
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'timePoints') 
    timePoints = options.timePoints ;
end
if isfield(options, 'very_far') 
    very_far = options.very_far ;
end
if isfield(options, 'coordSys') 
    coordSys = options.coordSys ;
end
if isfield(options, 'iLastikVersion') 
    iLastikVersion = options.iLastikVersion ;
end
if isfield(options, 'cellSize') 
    cellSize = options.cellSize ;
end
if isfield(options, 'strelRadius') 
    strelRadius = options.strelRadius ;
end
if isfield(options, 'gaussKernel') 
    gaussKernel = options.gaussKernel ;
end
if isfield(options, 'heightMinimum') 
    heightMinimum = options.heightMinimum ;
end

%% Load in h5 from ilastik.
if strcmpi(erase(coordSys, '_'), 'spsme') 
    Folder = [QS.dir.im_sp_sme, '_pixelClassification'] ;
    if ~exist(Folder, 'dir')
        mkdir(Folder)
        error(['Populate ' Folder ' with pixelClassification on pullbacks with coordSys ' coordSys])
    end
    filebase = [QS.fileBase.im_sp_sme(1:end-4) '_Probabilities.h5'] ;
elseif strcmpi(erase(coordSys, '_'), 'sprsme') || ...
        strcmpi(erase(coordSys, '_'), 'rspsme') || ...
        strcmpi(erase(coordSys, '_'), 'rsme')
    coordSys = 'sprsme' ;
    Folder = [QS.dir.im_r_sme, '_pixelClassification'] ;
    if ~exist(Folder, 'dir')
        mkdir(Folder)
        error(['Populate ' Folder ' with pixelClassification on pullbacks with coordSys ' coordSys])
    end
    filebase = [QS.fileBase.im_r_sme(1:end-4) '_Probabilities.h5'] ;
else
    error('Have not coded for this coordinate system yet. Do so here')
end

for tp = timePoints
    
    outfn = sprintf(QS.fullFileBase.segmentation2d, tp) ;
    if ~exist(outfn, 'file') || overwrite

        % Define path to this timePoint's hdf5 probabilities file
        h5fn = fullfile(Folder, sprintf(filebase, tp)) ;
        [ mem ] = load.ilastikh5Single( h5fn, iLastikVersion );

        %% Segment the membrane.
        % cellSize : int (default=200)
        %   kernel size for laplacian of Gaussian : set to scale of curvature
        %   picking out, around a cell size or higher, in units of area (pix^2)
        % strelRadius : int (default=1)
        %   strel disk radius for dilation of segmented image
        % gaussKernel : float (default=2)
        %   kernel size for Gaussian filter. Set to a couple pixels. 
        %   Has units of length (pix).
        % heightMinimum : float (default=3.5)
        %   height of any local minima to merge, to reduce noise at rugged
        %   minima
        L = seg.memWS(mem, cellSize, strelRadius, gaussKernel, heightMinimum);
        % Set bond=0 and clear_border = 1
        [L, Struct] = seg.generate_structs(L, 0, 1, 0, very_far);
        % Bad cells are bubble cells, which is a segmentation that forked and
        % reconnected.
        % ToDo: Should we do this? Can we skip it or does that lead to issues?
        disp('removing bad cells')
        L = seg.removeBadCells(Struct, L);
        % Now change label matrix after removing bad cells
        disp('relabelling cells')
        L = seg.relabelL(L);
        % Now also synchronize Struct after removing bad cells
        [segIm, seg2d] = seg.generate_structs(L, 0, 0, 0, very_far);
        disp('done with segmentation')

        %% Prepare data structure for inverse (optional? Does this improve segmentation?)
        % % put a parameter in the cdat of Struct, a boolean of whether every vertex
        % % is 3-fold.
        % Struct = seg.threefold_cell(Struct);
        % % generate the bdat structure in Struct
        % Struct = seg.recordBonds(Struct, L);
        % disp('generated the bond structure')
        % % Segment the curvature of each bond
        % Struct = seg.curvature(Struct, size(L));
        % disp('segmented the curvature of each bond')
        % % Remove all fourfold vertices, recursively if there are z>4
        % Struct = seg.removeFourFold(Struct);
        % disp('removed fourfold vertices')
        % % The inverse is ill-posed if we have convex cells, so hack those to be
        % % convex
        % Struct = seg.makeConvexArray(Struct);
        % disp('done with data preparation')

        %% Convert to simpler format
        disp('Constructing vdat')
        vdat = struct() ;
        vdat.v = zeros(length(seg2d.Vdat), 2) ;
        vdat.NL = zeros(length(seg2d.Vdat), 4) ;
        vdat.fourfold = false(length(seg2d.Vdat), 1) ;
        for qq = 1:length(seg2d.Vdat)
            vdat.v(qq, :) = [seg2d.Vdat(qq).vertxcoord, seg2d.Vdat(qq).vertycoord] ;
            nv = length(seg2d.Vdat(qq).nverts) ;
            try
                vdat.NL(qq, 1:nv) = seg2d.Vdat(qq).nverts ;
            catch
                % Increase the size of NL to accomodate more neighbors
                disp('Increasing NL size (dim 2)')
                swap = vdat.NL ;
                vdat.NL = zeros(length(seg2d.Vdat), nv) ;
                vdat.NL(1:qq, 1:size(swap, 2)) = swap(1:qq, :) ;
                vdat.NL(qq, 1:nv) = seg2d.Vdat(qq).nverts ;
            end
            vdat.fourfold(qq) = ~isempty(seg2d.Vdat(qq).fourfold) ;
        end    
        
        disp('generating bond list')
        BL = Vdat2BL(seg2d.Vdat) ;
        vdat.BL = BL ;
        
        % Assign polygons to vdat
        NL = vdat.NL ;
        coordSys = lower(erase(coordSys, '_')) ;
        if strcmpi(erase(coordSys, '_'), 'spsme') || ...
                strcmpi(erase(coordSys, '_'), 'rspsme') || ...
                strcmpi(erase(coordSys, '_'), 'sprsme')
            % provide ROI in [minx, maxx; miny, maxy]
            ROI = [-eps, size(segIm, 2) + eps; round([0.25, 0.75] * size(segIm, 1))] ;
        else
            error('Handle this coordSys here')
        end
        opts = struct() ;
        opts.roi = ROI ;
        
        seg2d.vdat = vdat ;
        seg2d.cdat = struct() ;
        if ~skipPolygons
            polygons = Cdat2polygons(seg2d.Cdat, vdat.v, BL, NL, opts) ;
            seg2d.cdat.polygons = polygons ;
        end
        seg2d.cdat.centroid = zeros(length(seg2d.Cdat), 2) ;
        for cid = 1:length(seg2d.Cdat)
            seg2d.cdat.centroid(cid, :) = seg2d.Cdat(cid).centroid.coord ;
        end
        
        %% Save the segmentation to disk
        if ~exist(fullfile(QS.dir.segmentation, 'seg2d'), 'dir')
            mkdir(fullfile(QS.dir.segmentation, 'seg2d'))
        end

        save(outfn, 'seg2d', 'segIm', 'coordSys')
    else
        disp(['already on disk: ' outfn])
        
        % %% Convert to simpler format
        load(outfn, 'seg2d', 'segIm', 'coordSys')
        
        
        
        % if ~isfield(seg2d, 'vdat')
        %     disp('Constructing vdat')
        %     vdat = struct() ;
        %     vdat.v = zeros(length(seg2d.Vdat), 2) ;
        %     vdat.NL = zeros(length(seg2d.Vdat), 4) ;
        %     vdat.fourfold = false(length(seg2d.Vdat), 1) ;
        %     for qq = 1:length(seg2d.Vdat)
        %         vdat.v(qq, :) = [seg2d.Vdat(qq).vertxcoord, seg2d.Vdat(qq).vertycoord] ;
        %         nv = length(seg2d.Vdat(qq).nverts) ;
        %         try
        %             vdat.NL(qq, 1:nv) = seg2d.Vdat(qq).nverts ;
        %         catch
        %             % Increase the size of NL to accomodate more neighbors
        %             disp('Increasing NL size (dim 2)')
        %             swap = vdat.NL ;
        %             vdat.NL = zeros(length(seg2d.Vdat), nv) ;
        %             vdat.NL(1:qq, 1:size(swap, 2)) = swap(1:qq, :) ;
        %             vdat.NL(qq, 1:nv) = seg2d.Vdat(qq).nverts ;
        %         end
        %         vdat.fourfold(qq) = ~isempty(seg2d.Vdat(qq).fourfold) ;
        %     end    
        % 
        %     disp('generating bond list')
        %     BL = Vdat2BL(seg2d.Vdat) ;
        %     vdat.BL = BL ;
        % 
        %     % Assign polygons to vdat
        %     NL = vdat.NL ;
        %     polygons = Cdat2polygons(seg2d.Cdat, vdat.v, BL, NL) ; 
        % 
        %     seg2d.cdat = struct() ;
        %     seg2d.cdat.polygons = polygons ;
        % 
        %     % Assign vdat to segmentation
        %     seg2d.vdat = vdat ;
        %     save(outfn, 'seg2d', 'segIm', 'coordSys')
        % end
        % 
        % % CLEANUP
        % if isfield(seg2d.vdat, 'polygons') || isfield(seg2d, 'polygons')
        %     if isfield(seg2d.vdat, 'polygons')
        %         pgs = seg2d.vdat.polygons ;
        %         seg2d.cdat = struct() ;
        %         seg2d.cdat.polygons = pgs ;            
        %         seg2d.vdat = rmfield(seg2d.vdat, 'polygons') ;
        %     end
        %     if isfield(seg2d, 'polygons')
        %         seg2d = rmfield(seg2d, 'polygons') ;
        %     end
        %     save(outfn, 'seg2d', 'segIm', 'coordSys')
        % end
        seg2d.cdat.centroid = zeros(length(seg2d.Cdat), 2) ;
        for cid = 1:length(seg2d.Cdat)
            seg2d.cdat.centroid(cid, :) = seg2d.Cdat(cid).centroid.coord ;
        end
        save(outfn, 'seg2d', 'segIm', 'coordSys')
    end
    
    %% Save image of the segmentation
    imfn = [outfn(1:end-3) 'png'] ;
    if ~exist(imfn, 'file') || overwrite || overwriteImages
        if strcmpi( coordSys, 'spsme')
            imageFn = sprintf(QS.fullFileBase.im_sp_sme, tp) ;
        elseif strcmpi( coordSys, 'sprsme')
            imageFn = sprintf(QS.fullFileBase.im_r_sme, tp) ;
        else
            error('Have not yet coded for this CoordSys');
        end
        im = imread(imageFn) ;
        clf
        Xs = zeros(size(seg2d.vdat.BL, 1), 1) ;
        Ys = Xs ;
        Us = Xs ;
        Vs = Xs ;
        dmyk = 1 ;
        for qq = 1:length(seg2d.Vdat)
            for id = seg2d.Vdat(qq).nverts
                Xs(dmyk) = seg2d.vdat.v(qq, 1) ;
                Ys(dmyk) = seg2d.vdat.v(qq, 2) ; 
                Us(dmyk) = seg2d.vdat.v(id, 1) - seg2d.vdat.v(qq, 1) ;
                Vs(dmyk) = seg2d.vdat.v(id, 2) - seg2d.vdat.v(qq, 2) ; 
                dmyk = dmyk + 1 ;
            end
        end
        % plot(seg2d.vdat.v(:, 1), seg2d.vdat.v(:, 2), '.')
        imshow(im) ;
        hold on;
        q = quiver(Xs,Ys, Us, Vs, 0, 'color', [ 0.8500    0.3250    0.0980]);
        axis equal
        t0 = QS.t0set() ;
        title(['t = ' sprintf('%03d', tp - t0) ' ' QS.timeUnits])
        q.ShowArrowHead = 'off';
        saveas(gcf, imfn)
    end
    
    
    %% Save image of the polygons
    imfn = [outfn(1:end-4) '_polygons.png'] ;
    if ~exist(imfn, 'file') || overwrite || overwriteImages || ~skipPolygons
        
        
        if strcmpi( coordSys, 'spsme')
            imageFn = sprintf(QS.fullFileBase.im_sp_sme, tp) ;
        elseif strcmpi( coordSys, 'sprsme')
            imageFn = sprintf(QS.fullFileBase.im_r_sme, tp) ;
        else
            error('Have not yet coded for this CoordSys');
        end
        im = imread(imageFn) ;
        
        opts = struct() ;
        [ff, vv] = polygons2triangulation(seg2d.cdat.polygons, ...
            seg2d.vdat.v, seg2d.cdat.centroid, opts) ;

        % vertex is part of input vertices in polygons or not
        colorV = zeros(size(vv, 1), 1) ;
        colorV(length(seg2d.vdat.v):end) = 255 ;
        
        clf
        imshow(cat(3, im, im, im)) ;
        hold on;
        Xs = zeros(size(seg2d.vdat.BL, 1), 1) ;
        Ys = Xs ;
        Us = Xs ;
        Vs = Xs ;
        dmyk = 1 ;
        for qq = 1:length(seg2d.Vdat)
            for id = seg2d.Vdat(qq).nverts
                Xs(dmyk) = seg2d.vdat.v(qq, 1) ;
                Ys(dmyk) = seg2d.vdat.v(qq, 2) ; 
                Us(dmyk) = seg2d.vdat.v(id, 1) - seg2d.vdat.v(qq, 1) ;
                Vs(dmyk) = seg2d.vdat.v(id, 2) - seg2d.vdat.v(qq, 2) ; 
                dmyk = dmyk + 1 ;
            end
        end
        q = quiver(Xs,Ys, Us, Vs, 0, 'color', [ 0.8500    0.3250    0.0980]);
        q.ShowArrowHead = 'off';
        axis equal
        hh = trisurf(ff, vv(:, 1), vv(:, 2), 0*vv(:, 2), ...
            'edgeColor', 'none', 'facealpha', 0.5) ;
        set(hh,'FaceColor','interp',...
           'FaceVertexCData',colorV,...
           'CDataMapping','scaled');
        t0 = QS.t0set() ;
        title(['t = ' sprintf('%03d', tp - t0) ' ' QS.timeUnits])
        colormap parula
        % saveas(gcf, imfn)
        disp(['saving image: ' imfn])
        current_ax = getframe(gca) ;
        imwrite(current_ax.cdata, imfn)
    end
end
