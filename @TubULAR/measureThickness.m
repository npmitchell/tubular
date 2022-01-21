function measureThickness(QS, options)
%measureThickness(options)
%
% Parameters
% ----------
% QS 
% options : optional struct with fields
%   coordSys : str (default='sp_r_sm')
%       coordinate system in which thickness stacks were trained
%
% Returns
% -------
% 
%
% NPMitchell 2020


%% Default options
overwrite = false ;
overwriteImages = false ;
coordSys = 'sp_r_sm' ;
preview = false ;
layer_spacing = 0.5 / QS.APDV.resolution ; % pixel resolution roughly matches xy
n_outward = 20 ;
n_inward = 40 ;
thres = 0.5 ;
tmin = 0 ;
tmax = 15 ;
foregroundChannel = 1 ;
t0 = QS.t0set() ;
if isempty(t0) 
    t0 = QS.xp.fileMeta.timePoints(1) ;
end

if nargin < 1
    options = struct() ;
end

%% Unpack options
if isfield(options, 'coordSys')
    coordSys = options.coordSys;
end
if isfield(options, 'thres')
    thres = options.thres ;
end
if isfield(options, 'foregroundChannel')
    foregroundChannel = options.foregroundChannel;
end



%% Obtain probability filenames as output from iLastik
res = layer_spacing * QS.APDV.resolution ;
if strcmpi(coordSys, 'sp_r_sm')
    fileNameBase = QS.fileBase.name ;
    spacingstr = strrep(sprintf('%0.2fum', res), '.', 'p') ;
    imfnBase = fullfile( QS.dir.im_r_sme_stack, ...
        [fileNameBase, '_%02d_%02d_' spacingstr '_Probabilities_thickness.h5']) ;
    extended = true ;
else
    error('Have not coded for this coordSys')
end

%% Extract thickness for each timepoint
outdir = sprintf(QS.dir.thickness, coordSys, n_outward, n_inward, res) ;
if ~exist(outdir, 'dir') 
    mkdir(outdir)
end
if ~exist(fullfile(outdir, 'images2d'), 'dir')
    mkdir(fullfile(outdir, 'images2d'))
end
if ~exist(fullfile(outdir, 'images3d'), 'dir')
    mkdir(fullfile(outdir, 'images3d'))
end

timePoints = QS.xp.fileMeta.timePoints ;
for qq = 1:length(timePoints)
    tp = timePoints(qq) ;
    outImFn2d = sprintf(QS.fullFileBase.thickness.ims2d, ...
        coordSys, n_outward, n_inward, res, tp) ; 
    outImFn3d = sprintf(QS.fullFileBase.thickness.ims3d, ...
        coordSys, n_outward, n_inward, res, tp) ; 
    outFn = sprintf(QS.fullFileBase.thickness.data, ...
        coordSys, n_outward, n_inward, res, tp) ; 
    
    if ~exist(outFn, 'file') || overwrite
        %% Load the segmentation
        imfn = sprintf(imfnBase, tp, n_outward, n_inward) ;
        disp(['Loading ' imfn])
        prob = h5read(imfn, '/exported_data') ;
        fg = squeeze(prob(foregroundChannel, :, :, :)) ;

        % Notes -- trash
        % bw = bwconncomp(fg > thres) ;
        % Only works in 2d
        % biggest = bwareafilt(fg > thres, 1);

        % Segment the data
        props = regionprops3(fg > thres, 'Volume') ;
        sortedVolumes = sort([props.Volume], 'descend') ;
        biggest = sortedVolumes(1);
        % Pull out 3 biggest into another 3-D logical image.
        isolated = bwareaopen(fg > thres, biggest);

        % Check the volume
        if preview
            for page = 1:5:size(isolated, 2)
                imshow(squeeze(isolated(:, page, :))) ;
                pause(0.01)
            end
        end

        % Compute the thickness from segmentation
        thick = zeros(size(isolated, 1), size(isolated, 2)) ;
        for ii = 1:size(isolated, 1)
            for jj = 1:size(isolated, 2)
                maxz = find(isolated(ii,jj,:),1,'last') ;
                minz = find(isolated(ii,jj,:),1,'first') ;
                if maxz
                    thick(ii, jj) = res * (maxz - minz) ;
                end
            end
        end

        %% Convert to sp/uv/etc coordinates
        [xx, yy] = ndgrid(1:size(thick, 1), 1:size(thick, 2)) ;
        xy = [xx(:), yy(:)] ;
        uv = QS.XY2uv(size(thick), xy, extended, 1.0, 1.0) ;

        %% Interpolate onto spcutMesh or other mesh coordinates
        uspace = linspace(min(uv(:, 1)), max(uv(:, 1)), nU) ;
        vspace = linspace(min(uv(:, 1)), max(uv(:, 1)), nV) ;
        thi = griddedInterpolant(uspace, vspace, thick) ;
        
        

        %% Save thickness to disk
    else
        % load from disk
        thickMat = load(outFn) ;
        error('here')
    end
    
    %% Save image in 2d
    if ~exist(outImFn2d, 'file') || overwrite || overwriteImages || true
        clf
        uspace = linspace(min(uv(:, 1)), QS.a_fixed * max(uv(:, 1)), size(thick, 1)) ;
        vspace = linspace(min(uv(:, 2)), max(uv(:, 2)), size(thick, 2)) ;
        imagesc(uspace, vspace, thick')
        cb = colorbar() ;
        ylabel(cb, ['thickness [' QS.spaceUnits ']'], 'interpreter', 'latex')
        axis equal
        axis off
        xlim([0, 1])
        ylim([0, 1])
        caxis([tmin, tmax])
        sgtitle(['$t=$' sprintf('%03d', (tp-t0) * QS.timeInterval) ' ' QS.timeUnits], ...
            'interpreter', 'latex')
        saveas(gcf, outImFn2d) ;
    end
    
    %% Save image in 3d
    if false && (~exist(outImFn3d, 'file') || overwrite || overwriteImages)
        clf
        trisurf(triangulation(mesh.f, mesh.v), thickv)
        cb = colorbar() ;
        ylabel(cb, ['thickness [' QS.spaceUnits ']'], 'interpreter', 'latex')
        axis equal
        axis off
        caxis([tmin, tmax])
        sgtitle(['$t=$' sprintf('%03d', (tp -t0)* QS.timeInterval) ' ' QS.timeUnits], ...
            'interpreter', 'latex')
        saveas(gcf, outImFn3d) ;
    end
    
    
end

end