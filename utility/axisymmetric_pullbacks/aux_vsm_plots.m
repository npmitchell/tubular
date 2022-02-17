% Auxiliary function for plotting smoothed meshes

%% Define directories
pivSimAvgImXDir = fullfile(pivSimAvgDir, 'vx') ;
pivSimAvgImYDir = fullfile(pivSimAvgDir, 'vy') ;
pivSimAvgImZDir = fullfile(pivSimAvgDir, 'vz') ;
pivSimAvgImTDir = fullfile(pivSimAvgDir, 'vtH') ;  % Heatmap
pivSimAvgImGDir = fullfile(pivSimAvgDir, 'vtG ') ;  % Gaussian smoothed in space
pivSimAvgImSDir = fullfile(pivSimAvgDir, 'vmag ') ;  % speed |v_3D|
pivSimAvgImQDir = fullfile(pivSimAvgDir, 'vtQ') ;  % Quiverplot
pivSimAvgImNDir = fullfile(pivSimAvgDir, 'vn') ;
% pivSimAvgImDvgDir = fullfile(pivSimAvgDir, 'dvg') ;
% pivSimAvgImCurlDir = fullfile(pivSimAvgDir, 'curl') ;
% pivSimAvgImShearDir = fullfile(pivSimAvgDir, 'shear_dvphi_ds') ;
dilDir = fullfile(pivDir, 'dilation') ;
vxyorigDir = fullfile(pivDir, 'vxyorig') ;
if plot_vxyz
    ensureDir(pivSimAvgImXDir)
    ensureDir(pivSimAvgImYDir)
    ensureDir(pivSimAvgImZDir)
end
dirs2make = {pivSimAvgImTDir, pivSimAvgImGDir,...
    pivSimAvgImNDir, dilDir, vxyorigDir, pivSimAvgImSDir, ...
    pivSimAvgImShearDir, pivSimAvgImDvgDir, pivSimAvgImCurlDir} ;
for pp = 1:length(dirs2make)
    ensureDir(dirs2make{pp}) ;
end

%% Make plots
% get size of images to make
gridsz = size(piv.x{1}) ;
bottom = round(gridsz(1) * 0.25) ;
top = round(gridsz(1) * 0.75) ;
% Display the velocities
close all
fig = figure('visible', 'off') ;
for i = 1:size(vsmM, 1)
    % Check if normal velocity plot exists
    vnfn = fullfile(pivSimAvgImNDir, [sprintf('%04d', time(i)) '.png']) ;
    vthfn = fullfile(pivSimAvgImTDir, [sprintf('%04d', time(i)) '.png']) ;
    vtgfn = fullfile(pivSimAvgImGDir, [sprintf('%04d', time(i)) '.png']) ;
    dilfn = fullfile(dilDir, [sprintf('%04d', time(i)) '.png']) ;
    vxyorigfn = fullfile(vxyorigDir, [sprintf('%04d', time(i)) '.png']) ;
    speedfn = fullfile(pivSimAvgImSDir, [sprintf('%04d', time(i)) '.png']) ;
    shearfn = fullfile(pivSimAvgImShearDir, [sprintf('%04d', time(i)) '.png']) ;

    % grab the tangential velocity for this timestep
    vsm_ii = squeeze(vsmM(i, :, :)) ;
    v2dsmum_ii = squeeze(v2dsmMum(i, :, :)) ;
    vnsm_ii = squeeze(vnsmM(i, :, :)) ;

    % Load the image to put flow on top
    fileName = split(fns(i).name, '.tif') ;
    fileName = fileName{1} ;
    im = imread(fullfile(fns(i).folder, fns(i).name)) ;
    im = cat(3, im, im, im) ;  % convert to rgb for no cmap change

    % Define Nx1 and Mx1 float arrays for xspace and yspace
    xx = piv.x{i}(1, :) ;
    yy = piv.y{i}(:, 1) ;

    % Define the proper coordinates (approx)
    % if ~exist(dvgfn, 'file') || ~exist(curlfn, 'file') || overwrite_vsm_plots || true
    %     dilation_allfaces = zeros(length(jac), 1) ;
    %     for f = 1:length(jac)
    %         qg = jac{f} * jac{f}' ;
    %         dilation_allfaces(f) = sqrt(det(qg)) ;
    %     end
    % end
    %     % load the spcutMesh
    %     load(sprintf(spcutMeshBase, time(i)), 'spcutMesh')
    %     ss = spcutMesh.sphi(:, 1) ;
    %     pp = spcutMesh.sphi(:, 1) ;
    %     [ss, pp] = gridDistancesInterpMesh3D() ;
    % end

    if plot_vxyz
        aux_plot_vxyz_simpleavg(im, vsm_ii, xx, yy, time(i), vscale, ...
            pivSimAvgImXDir, pivSimAvgImYDir, pivSimAvgImZDir) ;
    end

    % Plot the magnitude of the velocity 
    if ~exist(speedfn, 'file') || overwrite_vsm_plots
        disp(['Saving ' speedfn])
        close all
        fig = figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off') ;
        scalarFieldOnImage(im, xx, yy, vecnorm(vsm_ii, 2, 3), alphaVal, vtscale, '$|v|$ [$\mu$m/min]') ;
        ylim([0.25 * size(im, 1), 0.75 * size(im, 1)])
        saveas(fig, speedfn) ;
        close all
    end

    % % Check normals
    % quiver3(piv3d{i}.pt0(:, 1), piv3d{i}.pt0(:, 2), piv3d{i}.pt0(:, 3), ...
    %     piv3d{i}.normals(:, 1), piv3d{i}.normals(:, 2), piv3d{i}.normals(:, 3)) ;
    % 
    % % Check normals rotated and scaled
    % pt0_rs = ((rot * piv3d{i}.pt0')' + trans) * resolution  ;
    % quiver3(pt0_rs(:, 1), pt0_rs(:, 2), pt0_rs(:, 3), ...
    %     piv3d{i}.normals_rs(:, 1), piv3d{i}.normals_rs(:, 2), piv3d{i}.normals_rs(:, 3)) ;
    %

    % Look at smoothed 2d velocity fields
    vn = reshape(vnsm_ii, gridsz) ;
    vx = reshape(v2dsmum_ii(:, 1), gridsz) ;
    vy = reshape(v2dsmum_ii(:, 2), gridsz) ;

    % % Check normal velocity and rotations
    % quiver3(piv.x{i}, piv.y{i}, 0*piv.x{i}, vx, vy, vn, 0)

    % Get lobes for this timepoint
    foldx = ssfold_frac(i, :) * xesz ;

    % Plot the normal velocity on top
    if ~exist(vnfn, 'file') || overwrite_vsm_plots
        disp(['Saving ' vnfn])
        close all
        fig = figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off') ;
        scalarFieldOnImage(im, xx, yy, vn, alphaVal, vnscale, '$v_n$ [$\mu$m/min]') ;
        ylim([0.25 * size(im, 1), 0.75 * size(im, 1)])
        saveas(fig, vnfn) ;
        close all
    end

    % Plot the tangential velocity as heatmap on top of the image
    if ~exist(vthfn, 'file') || overwrite_vsm_plots
        disp(['Saving ' vthfn])
        imw = im * washout2d + max(im) * (1-washout2d) ;
        qopts.overlay_quiver = false ;
        qopts.label = '$v_t$ [$\mu$m/min]' ;
        qopts.outfn = vthfn ;
        vectorFieldHeatPhaseOnImage(imw, xx, yy, vx, vy, vtscale, qopts) ;
        clear qopts 
    end

    % Gaussian smooth the velocities
    if ~exist(vtgfn, 'file') || overwrite_vsm_plots
        disp(['Saving ' vtgfn])
        vxb = imgaussfilt(vx, 4) ;
        vyb = imgaussfilt(vy, 4) ;
        imw = im * washout2d + max(im) * (1-washout2d) ;
        qopts.qsubsample = qsubsample ;
        qopts.overlay_quiver = true ;
        qopts.qscale = 10 ;
        qopts.label = '$v_t$ [$\mu$m/min]' ;
        qopts.outfn = vtgfn ;
        % Plot the coarse-grained tang velocity as heatmap on top of the image
        vectorFieldHeatPhaseOnImage(imw, xx, yy, vxb, vyb, vtscale, qopts) ;    
        clearvars qopts
    end

    % Check dilation field
    if ~exist(dilfn, 'file') || overwrite_vsm_plots
        close all
        fig = figure('units', 'normalized', ...
                'outerposition', [0 0 1 1], 'visible', 'off') ;
        %imagesc(piv.x{i}(:), piv.y{i}(:), piv3d{i}.dilation)
        scalarFieldOnImage(im, xx, yy, ...
            reshape(log10(piv3d{i}.dilation), [length(xx), length(yy)]),...
            alphaVal, 0.5,...
            'dilation, $\log_{10}||J||$', 'style', 'diverging') ;
        ylim([size(im, 2) * 0.25, size(im, 2) * 0.75])
        saveas(gcf, dilfn)
        close all
    end

    % Find hyperbolic fixed points
    %

    % Plot original velocity -- note no smoothing done here ;)
    if ~exist(vxyorigfn, 'file') || overwrite_vsm_plots
        imw = im * washout2d + max(im) * (1-washout2d) ;
        opts.label = '$\tilde{v}$ [pix/min]' ;
        opts.outfn = vxyorigfn ;
        opts.qscale = 15 ;
        vectorFieldHeatPhaseOnImage(imw, xx, yy, ...
            piv.u_filtered{i}, piv.v_filtered{i}, 15, opts) ;
        clearvars opts
    end

    % % Plot divergence
    % if ~exist(dvgfn, 'file') || overwrite_vsm_plots 
    % 
    %     vxb = imgaussfilt(vx, 10) ;
    %     vyb = imgaussfilt(vy, 10) ;
    % 
    %     % Interpolate dilation onto locations where curl is defined
    %     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
    %     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
    %     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
    %     % dilv = dilation_allfaces(subfieldfaces) ;
    % 
    %     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
    %     dvg = divergence(piv.x{i}, piv.y{i}, vxb, vyb) .* dilum;
    %     opts.label = '$\nabla \cdot v_t$ [min$^{-1}$]' ;
    %     opts.title = [ '$t=$' num2str(time(i)) ' min' ] ;
    %     opts.outfn = dvgfn ;
    %     opts.qscale = 10 ;
    %     opts.sscale = 1.0 ;
    %     opts.alpha = 0.8 ;
    %     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
    %     scalarVectorFieldsOnImage(im, xx, yy, ...
    %         dvg, xx, yy, vxb, vyb, opts) ;
    %     clearvars opts
    % end
    % 
    % % Plot curl
    % if ~exist(curlfn, 'file') || overwrite_vsm_plots || true
    %     % xd = xx(1:10:end) ;
    %     % yd = yy(1:10:end) ;
    %     % xdgrid = xd .* ones(length(xd), length(yd)) ; 
    %     % ydgrid = (xd .* ones(length(yd), length(xd)))' ; 
    %     % check it
    %     % scatter(xdgrid(:), ydgrid(:), 10, xdgrid(:))
    % 
    %     % Interpolate dilation onto locations where curl is defined
    %     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
    %     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
    %     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
    %     % dilv = dilation_allfaces(subfieldfaces) ;
    % 
    %     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
    %     curlv = curl(piv.x{i}, piv.y{i}, vxb, vyb) ;
    %     curlv = curlv .* dilum ;
    %     opts.title = [ '$t=$' num2str(time(i)) ' min' ] ;
    %     opts.label = '$\nabla \times v_t$ [min$^{-1}$]' ;
    %     opts.outfn = curlfn ;
    %     opts.qscale = 10 ;
    %     opts.sscale = 0.5 ;
    %     opts.alpha = 0.8 ;
    %     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
    %     scalarVectorFieldsOnImage(im, xx, yy, ...
    %         curlv, xx, yy, vxb, vyb, opts) ;
    %     clearvars opts
    % end
    % 
    % % Plot strain dv_phi/ds on image
    % if ~exist(shearfn, 'file') || overwrite_vsm_plots || true
    %     % xd = xx(1:10:end) ;
    %     % yd = yy(1:10:end) ;
    %     % xdgrid = xd .* ones(length(xd), length(yd)) ; 
    %     % ydgrid = (xd .* ones(length(yd), length(xd)))' ; 
    %     % check it
    %     % scatter(xdgrid(:), ydgrid(:), 10, xdgrid(:))
    % 
    %     % Interpolate dilation onto locations where curl is defined
    %     % Di = scatteredInterpolant(tm0X, tm0Y, dilation_allfaces) ;
    %     % tr0 = triangulation(tm0f, [tm0X, tm0Y]) ;
    %     % [subfieldfaces, ~] = pointLocation(tr0, [xx(:), yy(:)]) ;
    %     % dilv = dilation_allfaces(subfieldfaces) ;
    % 
    %     dilum = reshape(piv3d{i}.dilation / resolution, size(vxb)) ;
    %     dvphidX = gradient(vyb, xx(2) - xx(1), yy(2) - yy(1)) ;
    %     dvphids = dvphidX .* dilum ;
    %     opts.title = [ '$t=$' num2str(time(i)) ' min' ] ;
    %     opts.label = '$\nabla_s v_{\phi}$ [min$^{-1}$]' ;
    %     opts.outfn = shearfn ;
    %     opts.qscale = 10 ;
    %     opts.sscale = 0.5 ;
    %     opts.alpha = 0.8 ;
    %     opts.ylim = [size(im, 2) * 0.25, size(im, 2) * 0.75] ;
    %     scalarVectorFieldsOnImage(im, xx, yy, ...
    %         dvphids, xx, yy, vxb, vyb, opts) ;
    %     clearvars opts
    % end

        
end