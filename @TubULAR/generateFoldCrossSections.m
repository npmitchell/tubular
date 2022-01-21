function generateFoldCrossSections(QS, options)
% generateFoldCrossSections(QS, options)
% 
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct
%   
% Returns
% -------
%
% NPMitchell 2021

%% Default options
buff = 30 ; % extra space around mesh to visualize in mips
thickness = 10 ; % thickness in QS.spaceUnits to take mips of along AP axis
preview = false ;
filterWindowSize = 10 ;

%% Unpack options
if isfield(options, 'buff')
    buff = options.buff ;
end
if isfield(options, 'thickness')
    thickness = options.thickness ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'filterWindowSize')
    filterWindowSize = options.filterWindowSize ;
end

%% Make mips
folds = QS.getFeatures('folds') ;
folds = folds.folds ;

[rot, trans] = QS.getRotTrans() ;
pix2um = QS.APDV.resolution ;

% Colors
colors = define_colors(7) ; 
blue   = colors(1, :) ;
red    = colors(2, :) ;
yellow = colors(3, :) ;
purple = colors(4, :) ;
green  = colors(5, :) ;
sky    = colors(6, :) ;
maroon = colors(6, :) ;
lscolor = sky ;

% Get locations of folds in 3d and smooth
timePoints = QS.xp.fileMeta.timePoints ;
[~, ~, tmp] = QS.getXYZLims() ;
xyzlims = [tmp(:, 1) - buff, tmp(:, 2) + buff] ;

foldap = zeros(size(folds)) ;
for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    QS.setTime(tp) ;
    
    mesh = QS.getCurrentSPCutMesh() ;
    
    for foldID = 1:size(folds, 2)
        foldap(tidx, foldID) = mesh.avgpts(folds(tidx, foldID), 1) ;
        
        if preview
            vv = reshape(mesh.v, [mesh.nU, mesh.nV, 3]) ;
            foldRing = squeeze(vv(folds(tidx, foldID), :, :)) ;
            clf
            trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none')
            hold on;
            plot3(foldRing(:, 1), foldRing(:, 2), foldRing(:, 3), '.-')
        end
    end
end

% Smooth fold positions in time
sigma = 5;
foldapF = gaussFilter1D(foldap, sigma, filterWindowSize) ;
close all
plot(foldapF); hold on; plot(foldap)
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel(['ap position of folds [' QS.spaceUnits ']'], ...
    'interpreter', 'latex')
saveas(gcf, fullfile(QS.dir.lobe, 'fold_positions.png'))

% Now sample the image in the vecinity of this point
first = true ;
tidx2do = [1:20:length(timePoints)] ;
tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do)] ;
for tidx = tidx2do

    tp = timePoints(tidx) ;
    QS.setTime(tp) ;
    timestr = sprintf('%06d', tp) ;
    disp(['Considering time ' timestr])

    % Get which slice of data for MIP
    % mesh = QS.getCurrentSPCutMeshSmRS() ;
    % vv = reshape(mesh.v, [mesh.nU, mesh.nV, 3]) ;

    % tiffn_i = fullfile(projectDir, sprintf(tiffn, str2double(timestr))) ;
    raw = QS.getCurrentData() ;
    IV = raw{1} ;

    % For plotting
    xslice = round(size(raw, 1) * 0.5) ;
    yslice = round(size(raw, 2) * 0.5) ;
    zslice = round(size(raw, 3) * 0.5) ;

    % Interpolate the data and evaluate at rotated orthogonal planes     
    active_rotation = false ; % too slow for large data
    if active_rotation
        % ACTIVE ROTATION: TRANSLATE & ROTATE intensity XYZ points
        disp('trans & rotate...')
        xi = 1:size(rawclip, 1) ;
        yi = 1:size(rawclip, 2) ;
        zi = 1:size(rawclip, 3) ;
        [xr, yr, zr] = meshgrid(xi, yi, zi) ;
        xyzr = (rot * [xr(:), yr(:), zr(:)]')' + trans ;
        xxx = reshape(xyzr(:, 1), size(rawclip)) ;
        yyy = reshape(xyzr(:, 2), size(rawclip)) ;
        zzz = reshape(xyzr(:, 3), size(rawclip)) ;

        % Interpolate the data and evaluate at rotated orthogonal planes
        disp('defining interpolation...')
        xspace = linspace(min(xxx(:)), max(xxx(:)), 100) ;
        yspace = linspace(min(yyy(:)), max(yyy(:)), 100) ;
        zspace = linspace(min(zzz(:)), max(zzz(:)), 100) ;
        [yq, zq ] = meshgrid(xspace, yspace, zspace) ;
        interpI = scatteredInterpolant(xxx(:), yyy(:), zzz(:), double(rawclip(:))) ;
        newIx = interpI(0*yq, yq, zq) ;
        error('Complete active rotation here')
    else
        % PASSIVE ROTATION            
        % Obtain rotation matrix that undoes the APDV rotation 
        invRot = QS.invertRotation(rot) ;

        disp('Creating data space for interpolation')
        % Plane to image
        xMin = xyzlims(1, 1) ;
        xMax = xyzlims(1, 2) ;
        yMin = xyzlims(2, 1) ;
        yMax = xyzlims(2, 2) ;
        zMin = -max(xyzlims(3, :)) ;
        zMax = max(xyzlims(3, :)) ;
        oversample_factor = 1 ;
        nX = round((xMax - xMin) / pix2um * oversample_factor) ;
        nY = round((yMax - yMin) / pix2um * oversample_factor) ; 
        nZ = round((zMax - zMin) / pix2um * oversample_factor) ;
        
        xspace = linspace(xMin, xMax, nX) ;
        yspace = linspace(yMin, yMax, nY) ;
        zspace = linspace(zMin, zMax, nZ) ;
        [x1, y1, z1] = meshgrid(xspace, yspace, zspace) ; 
        
        xyzI = (invRot * (([x1(:), y1(:), z1(:)]/pix2um) - trans)')' ;
        
        xpx = reshape(xyzI(:, 1), [nX, nY, nZ]) ;
        xpy = reshape(xyzI(:, 2), [nX, nY, nZ]) ; 
        xpz = reshape(xyzI(:, 3), [nX, nY, nZ]) ;

        % preview with mesh
        if preview
            clf
            h = trisurf(triangulation(mesh.f, mesh.v), ...
                'EdgeColor', 'none', 'facecolor', lscolor) ;
            axis equal
            hold on;
            yslice = 100 ;
            scatter3(reshape(xpx(:, yslice, :), [nX*nZ, 1]), ...
                reshape(xpy(:, yslice, :), [nX*nZ, 1]), ...
                reshape(xpz(:, yslice, :), [nX*nZ, 1]), 10, 'markeredgecolor', 'r')
            waitfor(gcf)
            imagesc(yspace, zspace, ix')
        end
        disp('performing interpolation in Raw Data space')
        % Putting in transposition by hand for XY --> YX
        ix = interp3(double(IV), xyzI(:, 2), xyzI(:, 1), xyzI(:, 3)) ;
        ix = reshape(ix, size(x1)) ;
        
        % Preview the stack along longitudinal axis
        if preview
            clf
            for slice = 1:10:size(ix, 2)
                imagesc(squeeze(ix(:, slice, :)))
                pause(0.1)
            end
        end
        
        for foldID = 1:size(folds, 2)
            
            slices = find(abs(xspace - foldapF(tidx, foldID)) ...
                < thickness * 0.5) ;
            im = squeeze(max(ix(:, slices, :), [], 2)) ;
            imfn = fullfile(QS.dir.lobe, 'constriction_MIPs', ...
                [sprintf(QS.fileBase.name, tp) sprintf('_fold%03d.png', foldID)]) ;
            disp(['Writing ' imfn])
            im(isnan(im)) = 0 ;
            im = uint16(im) ;
            imwrite(flipud(im'), imfn, 'png')
        end
    end
end