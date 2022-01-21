function measurePolarity(QS, options)
% measurePolarity(QS, options)
% 
%
% Parameters
% ----------
% options : struct with fields
%   method : farthestPoint or pullback
%   todo: implement farthestPoint radon
%
% Returns
% -------
%
% NPMitchell 2021

% Default options
wd = 40 ;
coordSys = 'spsme' ;
method = 'pullback'; % 'farthestPoint' or 'pullback'
lambda_mesh_Hopf = 0.002 ;

% Unpack options
if isfield(options, 'width')
    wd = options.width ;
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end
if isfield(options, 'lambda_mesh_Hopf')
    lambda_mesh_Hopf = options.lambda_mesh_Hopf ;
end

if strcmpi(coordSys, 'spsme')
    imDir = fullfile(QS.dir.im_sp_sme, 'radon', 'images') ;
    outDir = fullfile(QS.dir.im_sp_sme, 'radon') ;
    pbDir = QS.fullFileBase.im_sp_sme ;
elseif strcmpi(coordSys, 'spsmre') || strcmpi(coordSys, 'rsme')
    imDir = fullfile(QS.dir.im_r_sme, 'radon', 'images') ;
    outDir = fullfile(QS.dir.im_r_sme, 'radon') ;
    pbDir = QS.fullFileBase.im_r_sme ;
else
    error('Code for this coordSys')
end
if ~exist(imDir, 'dir')
    mkdir(imDir)
end

for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    QS.setTime(tp) 
    
    % Compute or load polarity measurement in this pullback frame
    
    outfn = fullfile(outDir, sprintf('%03d_radon.mat', tp)) ;    
    if ~exist(outfn, 'file') || overwrite
        % Load pullback image
        fn = sprintf(pbDir, tp) ;
        
        % full image is "imF" for "image Full"
        imF = imread(fn) ;
        imF = squeeze(imF(:, :, 1)) ;

        % Divide into rectangles
        xv = wd:wd:size(imF, 1) - wd ;
        yv = wd:wd:size(imF, 2) - wd ;
        [xgrid, ygrid] = meshgrid(xv, yv) ;

        angleRaw = zeros(length(xv), length(yv)) ;
        magRaw = zeros(length(xv), length(yv)) ;
        for pp = 1:length(xv)
            disp(['p = ' num2str(pp)])
            for qq = 1:length(yv)

                % Get chunk of image
                minx = max(1, xv(pp) - wd) ;
                maxx = min(size(imF, 1), xv(pp) + wd) ;
                miny = max(1, yv(qq) - wd) ;
                maxy = min(size(imF, 2), yv(qq) + wd) ;
                im = imF(minx:maxx, miny:maxy) ;

                % Mask out a circle from the patch
                % create a xygrid
                [xsz, ysz] = size(im) ;
                xcenter = xsz * 0.5 ;
                ycenter = ysz * 0.5 ;
                [xc, yc] = meshgrid(1:xsz, 1:ysz) ;
                minsz = min(xsz, ysz) ;
                dist = (xc - xcenter) .^2 + (yc - ycenter) .^2 ;
                mask = dist' < (minsz*0.5)^2 ;
                chunk = im .* uint8(mask) ;

                % [R] = radon(chunk, theta);

                % Compute radon as a function of angle
                options = struct() ;
                options.res = 2 ;
                [angleRaw(pp, qq), magRaw(pp,qq)] = extractRadonNematic(chunk, options) ;
            end
        end

        % Smooth resulting nematic
        ct = magRaw .* cos(2*angleRaw) ;
        st = magRaw .* sin(2*angleRaw) ;
        ct = imgaussfilt(ct, 1) ;
        st = imgaussfilt(st, 1) ;
        mag = sqrt(ct.^2 + st.^2) ;
        angle = reshape(0.5* atan2(st(:), ct(:)), [length(xv), length(yv)]) ;

        % Save as 2d grid
        close all
        mesh2d = struct() ;
        [YY, XX] = meshgrid(yv, xv) ;
        mesh2d.v = [XX(:), YY(:), 0*XX(:)] ;
        mesh2d.f = defineFacesRectilinearGrid([XX(:), YY(:)], ...
            length(xv), length(yv)) ;
        opts = struct() ;
        opts.clim_mag = rms1d(mag(:)) ;
        opts.axisOff = true ;
        plotNematicField(mag, angle, opts)
        figfn = fullfile(QS.dir.im_sp_sme, 'radon', 'images', sprintf('raw_%03d_radon.png', tp)) ;
        saveas(gcf, figfn)

        % Interpolate onto spcutMesh
        magT = mag' ;
        angleT = mod(angle', pi); 
        spMesh = QS.getCurrentSPCutMeshSmRSC() ;
        spMesh.u(:, 1) = spMesh.u(:, 1) / max(spMesh.u(:, 1)) ;
        uv = QS.XY2uv(imF, [XX(:), YY(:)], true, 1.0, 1.0) ;
        mI = scatteredInterpolant(uv(:, 1), uv(:, 2), magT(:), 'natural', 'nearest') ;
        magMesh = mI(spMesh.u(:, 1), spMesh.u(:, 2)) ;
        cI = scatteredInterpolant(uv(:, 1), uv(:, 2), cos(2*angleT(:)), 'natural', 'nearest') ;
        sI = scatteredInterpolant(uv(:, 1), uv(:, 2), sin(2*angleT(:)), 'natural', 'nearest') ;
        cuv = cI(spMesh.u(:, 1), spMesh.u(:, 2)) ;
        suv = sI(spMesh.u(:, 1), spMesh.u(:, 2)) ;
        angleMesh = mod(0.5 * atan2(suv, cuv), pi) ;

        % % Check it
        % subplot(2, 1, 1)
        % scatter(uv(:, 1), uv(:, 2), 10, cos(2*angle(:))); 
        % hold on;
        % scatter(spMesh.u(:, 1), spMesh.u(:, 2), 5, cuv, 'filled')
        % caxis([-1,1])
        % colormap bwr ; colorbar
        % 
        % subplot(2, 1, 2)
        % scatter(uv(:, 1), uv(:, 2), 10, sin(2*angle(:))); 
        % hold on;
        % scatter(spMesh.u(:, 1), spMesh.u(:, 2), 5, suv, 'filled')
        % caxis([-1,1])
        % colormap bwr ; colorbar
        % 
        % figure; 
        % subplot(2, 1, 1)
        % imagesc(mod(angleRaw, pi))
        % caxis([0, pi])
        % subplot(2, 1, 2)
        % imagesc(mod(angleT', pi))
        % phasemap; phasebar
        % caxis([0, pi])
        % 
        % imagesc(reshape(2*angleMesh, [nU, nV-1])')
        % phasemap; phasebar
        % caxis([0, 2*pi])

        % save as image
        close all
        figfn = fullfile(QS.dir.im_sp_sme, 'radon', 'images', ...
            sprintf('mesh_%03d_radon.png', tp)) ;
        m2d = struct() ;
        tmp = QS.getCurrentSPCutMeshSmRS() ;
        m2d.f = tmp.f ;
        m2d.v = [spMesh.u(:, 1), spMesh.u(:, 2), 0*spMesh.u(:, 2)] ;
        m2d.v(nU*(nV-1)+1:nU*nV, 1) = m2d.v(1:nU, 1) ;
        m2d.v(nU*(nV-1)+1:nU*nV, 2) = 1.0 ; 
        m2d.v(nU*(nV-1)+1:nU*nV, 3) = 0 ; 
        magMesh2d = magMesh ;
        angleMesh2d = angleMesh ;
        magMesh2d(nU*(nV-1)+1:nU*nV) = magMesh2d(1:nU) ;
        angleMesh2d(nU*(nV-1)+1:nU*nV) = angleMesh2d(1:nU) ;

        % Save as figure
        opts = struct() ;
        opts.clim = rms1d(magMesh(:)) ;
        opts.view = {[0, 0], [0, 90]} ;
        opts.axisOff = true ;
        opts.cbarlabels = {"$|R_{peak}|$", "$|R_{peak}|$"};
        opts.label = {'radon polarity', ''} ;
        nFieldsOnSurface({spMesh, m2d}, ...
            {{magMesh, angleMesh}, {magMesh2d, angleMesh2d}}, opts)
        saveas(gcf, figfn)

        % Save the angle and magnitude
        if ~exist(fullfile(QS.dir.im_sp_sme, 'radon'), 'dir')
            mkdir(fullfile(QS.dir.im_sp_sme, 'radon'))
        end
        save(outfn, 'angle', 'magnitude', 'mesh2d', 'angleMesh', 'magMesh')
    else
        load(outfn, 'angle', 'magnitude', 'mesh2d', 'angleMesh', 'magMesh')
    end

    pol = dvAverageNematic(magMesh, angleMesh) ;
    polarity_ap(tidx, :) = pol ;
    
    %% Compare to Hopf differential
    hopfFn = sprintf(QS.fullFileBase.metric, coordSys, lambda_mesh_Hopf, tp) ;
    load(hopfFn)
    
end

% Kymograph
