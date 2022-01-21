function generateCurrentSPCutMesh(QS, cutMesh, spcutMeshOptions)
% generateCurrentSPCutMesh(QS, cutMesh, spcutMeshOptions)
%
% Note that the only output in APDV (spaceUnits) coordinates are
%   mss, mcline, avgpts, avgpts_ss
%
%     spcutMesh.sphi 
%     spcutMesh.v 
%     spcutMesh.vn  
%     spcutMesh.ringpath_ss 
%     spcutMesh.radii_from_mean_uniform_rs  % from uniform DV sampling
%     spcutMesh.radii_from_avgpts 
%     spcutMesh.mss         % from uniform DV sampling, also stored in centerline
%     spcutMesh.mcline      % from uniform DV sampling, also stored in centerline
%     spcutMesh.avgpts      % from uniform DV sampling, also stored in centerline
%     spcutMesh.avgpts_ss   % pathlength from uniform sampling, also
%           stored in centerline, in QS.spaceUnits
%     spcutMesh.ar          % affine scaling in X that minimizes isoareal energy
%     spcutMesh.phi0s       % 
%     spcutMesh.phi0_fit    %
%
% Note that depending on the relative axisorder between data frame (TIFFs)
% and the APDV frame, you might need to set phi0_sign to true which would
% flip the sign of phi0 (added rather than subtracted from v in (u,v)).
% This can be true even if QS.flipy is false for some nonstandard 
% circumstances.
%
%
% Parameters
% ----------
% QS : QuapSlap class instance
%   Note that the following properties are used:
%       QS.phiMethod = ('3dcurves', 'texture', 'combined') 
%       QS.a_fixed = 2.0    
% cutMesh : cutMesh struct, optional
%   cutMesh with fields
% spcutMesh : spcutMesh struct, optional
% spcutMeshOptions : struct with fields, optional
%   overwrite : bool
%       overwrite previous results
%   save_phi0patch : bool
%       show the relaxation steps of phi0 determination
%   iterative_phi0 : bool
%       iteratively determine phi0 until convergence
%   smoothingMethod : str specifier (default='none')
%       method for smoothing phi0 wrt AP axis coordinate (ss)
%   textureAxisOrder : 'xyz', 'yxz', etc
%       
%
% Returns
% -------
% spcutMesh : struct with fields
%   sphi: equal dv, resampled 3d points based on spcutMesh.sphi0
%
% NPMitchell 2020

%% Default options
overwrite = false ;
save_phi0patch = false ;
iterative_phi0 = false ;
smoothingMethod = 'none' ;
maxPhiIterations = 10 ;     % maximum # iterations of minimizing phi via phi -> phi-phi0(s) 
textureAxisOrder = QS.data.axisOrder ;
save_ims = QS.plotting.save_ims ;
phi0_sign = QS.flipy ;  % NOTE: our convention is that new phi = v - phi_0 but for
                        % some reason when flipy is true even though we have
                        % flipped cutMesh, we must add a minus sign so that phiv_kk
                        % becomes phi_kk = v + phi0. 

%% Unpack options
if nargin < 2 || isempty(cutMesh)
    if isempty(QS.currentMesh.cutMesh)
        QS.loadCurrentCutMesh()
    end
    cutMesh = QS.currentMesh.cutMesh ;
end
if nargin > 2
    disp('Unpacking options')
    if isfield(spcutMeshOptions, 'overwrite')
        overwrite = spcutMeshOptions.overwrite ;
    end
    if isfield(spcutMeshOptions, 'save_phi0patch')
        save_phi0patch = spcutMeshOptions.save_phi0patch ;
    end
    if isfield(spcutMeshOptions, 'iterative_phi0')
        iterative_phi0 = spcutMeshOptions.iterative_phi0 ;
    end
    if isfield(spcutMeshOptions, 'smoothingMethod')
        smoothingMethod = spcutMeshOptions.smoothingMethod ;
    end
    if isfield(spcutMeshOptions, 't0_for_phi0')
        t0_for_phi0 = spcutMeshOptions.t0_for_phi0 ;
        disp(['Set t0 from options:' num2str(t0_for_phi0)])
    else
        % The default timepoint at which we set phi0=0 is either the first feature
        %   onset time or else the first timepoint
        try
            t0_for_phi0 = QS.t0set() ;
            disp(['Setting t0 = ' num2str(t0_for_phi0)])
        catch
            disp('Setting t0 = first timepoint')
            t0_for_phi0 = QS.xp.fileMeta.timePoints(1) ; 
        end
    end
    if isfield(spcutMeshOptions, 'maxPhiIterations')
        maxPhiIterations = spcutMeshOptions.maxPhiIterations ;
    end
    if isfield(spcutMeshOptions, 'textureAxisOrder')
        textureAxisOrder = spcutMeshOptions.textureAxisOrder ;
    end
    if isfield(spcutMeshOptions, 'save_ims')
        save_ims = spcutMeshOptions.save_ims ;
    end
    if isfield(spcutMeshOptions, 'phi0_sign')
        phi0_sign = spcutMeshOptions.phi0_sign ;
    end
end

% populate patchOpts from pbOptions
pbOptions = struct() ;
pbOptions.resolution = QS.xp.fileMeta.stackResolution(1) ;
pbOptions.axisorder = textureAxisOrder ;

%% Unpack QS
tt = QS.currentTime ;
nU = QS.nU ;
nV = QS.nV ;
a_fixed = QS.a_fixed ;
phi_method = QS.phiMethod ;
spcutMeshfn = sprintf(QS.fullFileBase.spcutMesh, tt) ;
fileNameBase = QS.fileBase.name ; 
phi0fitBase = QS.fullFileBase.phi0fit ;
spcutMeshBase = QS.fullFileBase.spcutMesh ;
[rot, trans] = getRotTrans(QS) ;
resolution = QS.APDV.resolution ;
QS.getCleanCntrlines ;
cleanCntrlines = QS.cleanCntrlines ;
cleanCntrline = cleanCntrlines{QS.xp.tIdx(tt)} ;
preview = QS.plotting.preview ;
[~, ~, xyzlim_um] = QS.getXYZLims() ;
clineDVhoopDir = QS.dir.clineDVhoop ;
clineDVhoopBase = QS.fullFileBase.clineDVhoop ;
sphiDir = QS.dir.spcutMesh ;

% Check that options are consistent
if strcmp(phi_method, '3dcurves') && strcmp(smoothingMethod, 'none') && iterative_phi0 
    % error(['smoothingMethod for phi0 is none, so iterating phi0 will ', ...
    %    'do no good when phi_method is geometric (not texture)'])
end

% Expand dirs for images
clineDVhoopImDir = fullfile(clineDVhoopDir, 'images') ;
clineDVhoopFigBase = fullfile(clineDVhoopImDir, 'clineDVhoop_%06d.png') ;
if ~exist(clineDVhoopImDir, 'dir')
    mkdir(clineDVhoopImDir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate s,phi coord system for rotated, scaled mesh (rs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Establishing s,phi coord system\n');
if ~exist(spcutMeshfn, 'file') || overwrite
    if overwrite
        disp('Overwriting spcutMesh...')
    else
        disp('spcutMesh not on disk. Generating ...')
    end

    % Transform from u,v coordinates to s, phi coordinates
    %----------------------------------------------------------------------
    % Generate tiled orbifold triangulation
    %----------------------------------------------------------------------
    tileCount = [1 1];  % how many above, how many below
    cutMeshrs = cutMesh;
    
    % DEBUG 12-03-2020
    % Rotate and translate TV3D
    cutMeshrs.v = QS.xyz2APDV(cutMesh.v) ;
    % USED TO DO THIS but does not flip in y
    % cutMeshrs.v = ((rot * cutMesh.v')' + trans) * resolution ;
    
    % DEBUG 12-03-2020 -- added flipping vn here
    cutMeshrs.vn = (rot * cutMesh.vn')' ;
    if QS.flipy
        cutMeshrs.vn(:, 2) = - cutMeshrs.vn(:, 2) ;
    end
    
    [ ~, ~, TV3D, TVN3D ] = tileAnnularCutMesh( cutMesh, tileCount );
    [ TF, TV2D, TV3Drs ] = tileAnnularCutMesh( cutMeshrs, tileCount );

    % check tiled mesh coords in uv space
    % plot(TV2D(:, 1), TV2D(:, 2), '.')
    % title('tiled cutMesh')
    % pause(2)
    % close all
    
    %----------------------------------------------------------------------
    % Calculate abbreviated centerline from cutMesh boundaries
    %----------------------------------------------------------------------
    % Load centerline from cleaned list
    cline = QS.xyz2APDV(cleanCntrline(:, 2:4)) ; 
    ss = cleanCntrline(:, 1) * QS.APDV.resolution ;
    disp('Finding relevant segment of centerline')
    [cseg, acID, pcID, ~, ~] = ...
        centerlineSegmentFromCutMesh(cline, TF, TV2D, TV3Drs) ;

    %----------------------------------------------------------------------
    % Generate surface curves of constant s
    %----------------------------------------------------------------------
    % For lines of constant phi
    disp('Creating crude uv curves with du=const to define uspace by ds(u)')
    % Make grid
    eps = 1e-14 ;
    uspace0 = linspace( eps, cutMesh.umax - eps, nU )' ;
    vspace = linspace( eps, 1-eps, nV )' ;

    disp('Casting crude (equal dU) points into 3D...')
    crude_ringpath_ss = ringpathsGridSampling(uspace0, vspace, TF, TV2D, TV3Drs) ;

    % Resample crude_ringpath_ds made from uspace0 (uspace0 had equal du, not equal ds_3D in u direction)
    [uspace_ds, eq_ringpath_ss] = equidistantSampling1D(linspace(0, 1, nU)', crude_ringpath_ss, nU, 'linear') ;
    % ensure that uspace is nU x 1, not 1 x nU
    uspace_ds = reshape(uspace_ds, [nU, 1]) ; 
    % hedge the first and last point to avoid NaNs
    eps = 1e-13 ;
    uspace_ds(1) = uspace_ds(1) + eps ;
    uspace_ds(end) = uspace_ds(end) - eps ;
    clearvars dsuphi curves3d 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Casting resampled points into 3D (approx equal ', ...
        'ds_3D in u dir, but variable ds_3D in v dir)...'])
    % NOTE: first dimension indexes u, second indexes v
    curves3d = zeros(nU, nV, 3) ;  % in units of um
    for kk = 1:nU
        if mod(kk, 20) == 0
            disp(['u = ' num2str(kk / nU)])
        end
        uv = [cutMesh.umax * uspace_ds(kk) * ones(size(vspace)), vspace] ;
        curves3d(kk, :, :) = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, uv) ;
    end 

    % Check the 3d curves 
    if preview
        figure ; hold on;
        for kk = 1:nU
            plot3(curves3d(kk, :, 1), curves3d(kk, :, 2), curves3d(kk, :, 3), '.') 
        end
        title('curves3d')
        axis equal
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Compute s(u) and radius(u) for "uniform"--> ', ...
            'evenly sample each DV hoop (0,1) so ds_3D=const \n']);
    % Note: c3d_dsv has equal spacing in both ds and dphi, units of um in
    % rotated & scaled (RS) coordinates --> is it flipped in Y?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Resample at evenly spaced dphi in embedding space (rs, in um)
    fprintf('Resampling curves...\n')
    c3d_dsv = zeros(size(curves3d)) ;  % in units of um
    for i=1:nU
        % Note: no need to add the first point to the curve
        % since the endpoints already match exactly in 3d and
        % curvspace gives a curve with points on either
        % endpoint (corresponding to the same 3d location).
        c3d_dsv(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d(i, :, :)), nV, true) ;
        if vecnorm(squeeze(c3d_dsv(i, 1, :)) - squeeze(c3d_dsv(i, end, :))) > 1e-7
            error('endpoints do not join! Exiting')
        end

        % Visualization for Troubleshooting:
        % triplot(TF, TV2D(:, 1), TV2D(:, 2))
        % hold on;
        % plot(uv(:, 1), uv(:, 2), '.')
    end

    % Check the 3d curves 
    if preview
        figure ; hold on;
        for kk = 1:nU
            plot3(c3d_dsv(kk, :, 1), c3d_dsv(kk, :, 2), c3d_dsv(kk, :, 3), '.') 
        end
        title('c3d_dsv')
        axis equal
    end

    fprintf('Finding s(u) and r(u) of resampled "uniform" c3ds [uniform ds in V dir]...\n')
    % mcline is the resampled centerline, with mss
    % avgpts is the Nx3 averaged hoops (hoops are resampled already), with 
    % avgpts_ss being the associated pathlength along the centerline
    [mss, mcline, radii_from_mean_uniform_rs, avgpts_ss, avgpts] = ...
        srFromDVCurves(c3d_dsv) ;

    % Used to find radius using original centerline
    % [ssv, radii, avgpts, cids] = srFromDVCurvesGivenCenterline(ss, cline, c3ds) ;
    % Could operate just on the centerline segment
    cseg_ss = ss(acID:pcID) ;
    % [ssv, radii, avgpts, cids] = srFromDVCurves(cseg_ss, cseg, c3ds) ;
    % 
    % Adjust the centerline indices to index into the full
    % centerline. Note that cseg_ss already does this for ss.
    % cids = cids + acID ;

    % Plot new centerline
    aux_plot_clineDVhoop(QS, avgpts, avgpts_ss, cseg, cline, cseg_ss, ...
        curves3d, xyzlim_um, clineDVhoopFigBase, tt)
    
    % Optional: clean curve with polynomial and point match
    % avgpts onto cleaned curve. Skipping for later.

    % Compute ringpath_ss, the mean distance traveled from one
    % line of constant u to the next
    disp('Computing ringpath_ss in "uniform" resampling (equal ds along DV)...')
    % The distance from one hoop to another is the
    % difference in position from (u_i, v_i) to (u_{i+1}, v_i).
    dsuphi = reshape(vecnorm(diff(c3d_dsv), 2, 3), [nU-1, nV]) ;
    ringpath_ds = nanmean(dsuphi, 2) ;
    ringpath_ss = cumsum([0; ringpath_ds]) ;
    clearvars dsuphi ringpath_ds

    % Save new centerline in rotated translated units
    fn = sprintf(clineDVhoopBase, tt) ;
    disp(['Saving new centerline to ' fn])
    save(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')

    % Note: radii_from_mean_uniform_rs is the radius of 
    % interpolated hoops, not the actual points

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Done making new centerline using uniformly sampled hoops\n') ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Create new3d, the grid of UV coords')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    onesUV = ones(nU, nV) ;
    uspace_ds_umax = uspace_ds * cutMesh.umax ;
    uu = uspace_ds_umax .* onesUV ;
    vv = (vspace .* onesUV')' ;
    % uv is equally spaced in ss, UNEQUALLY spaced in uu
    uv = [uu(:), vv(:)] ;
    % Note: here interpolate uv in the TV2D coord system, then
    % use uphi as the actual 2D coordinates for these vertices
    % NOTE: unlike curves3d, uvgrid3d is NOT rotated/translated/scaled
    uvgrid3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, uv) ;
    
    %% Make avgpts in pixel space (not RS)
    fprintf('Resampling uvgrid3d curves in pix...\n')
    curves3d_pix = reshape(uvgrid3d, [nU, nV, 3]) ;
    c3d_dsv_pix = zeros(size(curves3d_pix)) ;  % in units of pix
    avgpts_pix = zeros(nU, 3) ;
    for i=1:nU
        % Note: no need to add the first point to the curve
        % since the endpoints already match exactly in 3d and
        % curvspace gives a curve with points on either
        % endpoint (corresponding to the same 3d location).
        c3d_dsv_pix(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d_pix(i, :, :)), nV, true) ;
        if vecnorm(squeeze(c3d_dsv_pix(i, 1, :)) - squeeze(c3d_dsv_pix(i, end, :))) > 1e-7
            error('endpoints do not join! Exiting')
        end
        avgpts_pix(i, :) = mean(squeeze(c3d_dsv_pix(i, :, :)), 1) ; 
    end
    
    %%
    if tt == t0_for_phi0
        % Store for next timepoint
        phiv = (vspace .* ones(nU, nV))' ;
        phi0s = zeros(size(uspace_ds)) ;
        phi0_fit = phi0s ;
        compute_phi0 = false ;
    elseif tt > t0_for_phi0
        tidx = QS.xp.tIdx(tt) ;
        tp_for_comparison = QS.xp.fileMeta.timePoints(tidx-1) ;
        compute_phi0 = true ;
    elseif tt < t0_for_phi0 
        tidx = QS.xp.tIdx(tt) ;
        tp_for_comparison = QS.xp.fileMeta.timePoints(tidx+1) ;
        compute_phi0 = true ;
    end
    
    if compute_phi0
        % Load previous sphi vertices in 3d 
        if strcmp(phi_method, '3dcurves') || strcmp(phi_method, 'combined') 
            disp('Computing phi(v) via 3dcurve matching (geometric method)')
            % Load the previous spcutMesh and call it prev3d_sphi
            % Also note the previous spcutMesh pullback image's fn
            tmp = load(sprintf(spcutMeshBase, ...
                tp_for_comparison), 'spcutMesh') ;
            prevf = tmp.spcutMesh.f ;
            prev3d_sphi = reshape(tmp.spcutMesh.v, [nU, nV, 3]) ; 
            
            %% Resample the vertices evenly around each hoop!
            prev3d_sphi_dsdv = zeros(size(prev3d_sphi)) ;  % in units of pix
            for i=1:nU
                % Note: no need to add the first point to the curve
                % since the endpoints already match exactly in 3d and
                % curvspace gives a curve with points on either
                % endpoint (corresponding to the same 3d location).
                prev3d_sphi_dsdv(i, :, :) = resampleCurvReplaceNaNs(squeeze(prev3d_sphi(i, :, :)), nV, true) ;
                if vecnorm(squeeze(prev3d_sphi_dsdv(i, 1, :)) - squeeze(prev3d_sphi_dsdv(i, end, :))) > 1e-7
                    error('endpoints do not join! Exiting')
                end

                % Visualization for Troubleshooting:
                % triplot(TF, TV2D(:, 1), TV2D(:, 2))
                % hold on;
                % plot(uv(:, 1), uv(:, 2), '.')
            end
            
            %% Obtain previous avgpts to use for hoop matching
            % fn_prev_mcline = sprintf(clineDVhoopBase, tt - 1) ;
            disp('Grabbing previous centerline')
            % prev_avgpts = load(fn_prev_mcline, 'avgpts') ;
            % prev_avgpts = prev_avgpts.avgpts ;
            [~, ~, ~, prev_avgpts_ss_pix, prev_avgpts_pix] = ...
                srFromDVCurves(prev3d_sphi_dsdv) ;
            
            %% Obtain previous pullback image
            % prev2d_uphi = reshape(tmp.spcutMesh.uphi, [nU, nV, 2]) ;
            imfn_sp_prev = sprintf( QS.fullFileBase.im_sp, tp_for_comparison) ;

            % fit the shifts in the y direction
            dmyk = 0 ;
            phi0_fit = zeros(size(uspace_ds)) ;
            phi0s = zeros(size(uspace_ds)) ;
            phi0_fit_kk = 1 ; % for first pass, so that we enter loop              
            phiv_kk = (vspace .* ones(nU, nV))' ;
            ensureDir(fullfile(sphiDir, 'phi0_correction'))
            
            % Note: the last boolean ensures that if iterative_phi0 is
            % false, we only do one pass. 
            % If iterative_phi0 is true, then we do a second/third/etc pass
            % if the adjustment from the fit is large.
            do_iteration = true ;
            while any(phi0_fit_kk > 0.005) && do_iteration
                % Make sure we do only one pass if not iterative
                if ~iterative_phi0 || dmyk > maxPhiIterations
                    do_iteration = false ;
                end
                disp(['Iteration ' num2str(dmyk)])
                plotfn = sprintf(phi0fitBase, tt, dmyk);

                % Pass phiOptions to allow sliding along AP
                % passed to phiOffsetsFromPrevMesh
                phiOpts = struct(); 
                phiOpts.preview = preview ;
                phiOpts.avgpts = avgpts_pix ;
                phiOpts.c3d_dsv = c3d_dsv_pix ;
                phiOpts.prev_avgpts_ss = prev_avgpts_ss_pix ;
                phiOpts.prev_avgpts = prev_avgpts_pix ;
                
                % Will we save check pullbacks to preview the algo?
                % If so, create a struct to pass visualization options
                if save_phi0patch
                    if dmyk == 0 
                        % first pass --> vvals4plot is same as vvals
                        phi4plot = (vspace .* ones(nU, nV))' ;
                    end
                    patchImFn = sprintf( ...
                        fullfile(sphiDir, 'phi0_correction',...
                        [fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
                        tp_for_comparison)  ;
                    patchImFnRes = sprintf( ...
                        fullfile(sphiDir, 'phi0_correction',...
                        [fileNameBase, '_phi0residual_' num2str(dmyk) '.tif']), ...
                        tp_for_comparison)  ;
                    geomImFn = sprintf( ...
                        fullfile(sphiDir, 'phi0_correction', ...
                        ['3d' fileNameBase '_prephi0_' num2str(dmyk) '.tif']), ...
                        tp_for_comparison)  ;

                    % Load the intensity data for this timepoint
                    if isempty(QS.currentData.IV)
                        QS.getCurrentData()
                    end
                    IV = QS.currentData.IV ;
                    
                    % Texture patch options
                    Options.PSize = 5;
                    Options.EdgeColor = 'none';
                    % Texture image options
                    Options.imSize = ceil( 1000 .* [ 1 a_fixed ] );
                    Options.yLim = [0 1];

                    % Roll options into a struct
                    patchOpts = pbOptions ;
                    patchOpts.patchImFn = patchImFn ;
                    patchOpts.patchImFnRes = patchImFnRes ;
                    patchOpts.imfn_sp_prev = imfn_sp_prev ;
                    patchOpts.IV = IV ;
                    patchOpts.ringpath_ss = ringpath_ss ;
                    patchOpts.v3d = uvgrid3d ;
                    patchOpts.vvals4plot = phi4plot ;
                    % passed to texturePatchToImage
                    patchOpts.Options = Options ;
                    patchOpts.phi0_sign = phi0_sign ;
                else                
                    patchOpts = struct() ;
                end

                % Minimize difference in DV hoop positions wrt
                % previous pullback mesh   
                % Outputs are raw phi0 and smoothed phi0s
                %% THIS WORKS
                [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevMesh(TF,...
                    TV2D, TV3D, uspace_ds_umax, phiv_kk, ...
                    prev3d_sphi_dsdv, -0.45, 0.45, ...
                    save_ims, plotfn, smoothingMethod, ...
                    phiOpts, patchOpts) ;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % USED TO USE: uspace * cutMesh.umax
                % % BREAK IT DOWN
                % % Using simple offset
                % phi0s_kk = phiOffsetsFromPrevMesh(TF, TV2D, TV3D, ...
                %     uspace * cutMesh.umax, phiv_kk, prev3d_sphi, ...
                %     -0.45, 0.45, {preview}) ;
                % phi0_fit_kk = phi0s_kk;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Update the result
                dmyk = dmyk + 1;
                phi0_fit = phi0_fit + phi0_fit_kk ;
                phi0s = phi0s + phi0s_kk ;
                % phiv_kk is this iteration's updated correction 
                % NOTE: our convention is that new phi = v - phi_0 but for
                % some reason when flipy is true even though we have
                % flipped cutMesh, we must add a minus sign so that phiv_kk
                % becomes phi_kk = v + phi0. 
                if phi0_sign > 0
                    phiv_kk = (vspace .* ones(nU, nV))' + phi0_fit .* ones(nU, nV) ;
                    phi4plot = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
                    phi4plot = mod(phi4plot, 1) ;
                else
                    phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
                    % todo: check that this is right
                    phi4plot = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
                end
                
                %% Preview result
                if save_phi0patch
                    % Colored mesh with phi=0 contour
                    close all
                    figure('visible', 'off')
                    % plot mesh colored by the phase phi 
                    % previous timepoint
                    xprev = prev3d_sphi(:, :, 1) ;
                    yprev = prev3d_sphi(:, :, 2) ;
                    zprev = prev3d_sphi(:, :, 3) ;
                    phiprev = (vspace .* ones(nU, nV))' ;
                    colormap parula ;
                    % cmap = parula ;
                    % colors = cmap(max(1, uint8(colortmp(:) * length(parula))), :) ;
                    prevh = trisurf(prevf, xprev(:), yprev(:), zprev(:), phiprev(:), ...
                        'FaceColor', 'interp',...
                        'EdgeColor', 'none', 'FaceAlpha', 0.25) ;
                    axis equal
                    % before fitting
                    hold on;
                    pe0 = find(phiprev(:) < 1e-4 | phiprev(:) > 0.99) ;
                    pe0h = plot3(uvgrid3d(pe0, 1), uvgrid3d(pe0, 2), ...
                        uvgrid3d(pe0, 3), '.') ;
                    % after fitting
                    hold on;                   
                    pekk = find(phi4plot(:) < 1e-4 | phi4plot(:) > 0.99) ;
                    pekkh = plot3(uvgrid3d(pekk, 1), uvgrid3d(pekk, 2), ...
                        uvgrid3d(pekk, 3), '.') ;         
                    % format plot
                    legend({'prev mesh', 'before fitting', 'after fitting'},...
                        'location', 'eastoutside')
                    title(['phi0 preview, iteration ' num2str(dmyk)])
                    xlabel('x [\mum]')
                    ylabel('y [\mum]')
                    zlabel('z [\mum]')
                    view(2)
                    disp(['Saving phi0 preview to ' geomImFn])
                    saveas(gcf, geomImFn)
                    close all
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Plot vertices in 3d
                    % TF, TV2D, TV3D, uspace, vvals, prev3d_sphi
                    % for qq = 1:nU
                    %     % where are the displacement vectors
                    %     interpolate2Dpts_3Dmesh(TF, TV2D, ...
                    %         TV3Drs, [uspace(qq) * ones(nV, 1), ...
                    %         mod(vqq + phi_kk, 1)] ) ...
                    %         - prev3dvals, 2, 2) .^ 2)
                    % end
                    
                end
            end
        elseif strcmp(phi_method, 'texture') 
            disp('Computing phi(v) via texture matching (physical method)')
            error('adjust this method to use later timepoint if tt>t0_for_comparison')
            phi0_fit = QS.fitPhiOffsetsViaTexture(uspace_ds_umax, vspace) ;
        else
            error(["Could not recognize phi_method: ", ...
                "must be 'texture' or '3dcurves' or 'combined'"])
        end
        
        % If we use a combined method, use curves3d as initial guess for
        % texture method
        if strcmp(phi_method, 'combined') 
            error('adjust this method to use later timepoint if tt>t0_for_comparison')
            disp('Refining phi(v) via texture matching (physical method)')
            phi0_fit = QS.fitPhiOffsetsViaTexture(uspace_ds_umax, vspace,...
                        phi0_fit) ;
        end
        close all

        % Store to save at this timepoint
        % phiv = (vspace .* ones(nU, nV))' + phi0_fit .* ones(nU, nV) ;
        % phiv = phiv_kk ;
        phiv = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;    
                
    end

    % NOTE: We have coordinates u,phiv that we associate with
    % the 3d coordinates already mapped to uv
    uphi = [uu(:), phiv(:)] ;

    % Recompute radii_from_mean_uniform_rs as radii_from_avgpts 
    % NOTE: all radius calculations done in microns, not pixels
    sphi3d_rs = ((rot * uvgrid3d')' + trans) * resolution ;
    if QS.flipy
        sphi3d_rs(:, 2) = - sphi3d_rs(:, 2) ;
    end
    radii_from_avgpts = zeros(size(sphi3d_rs, 1), size(sphi3d_rs, 2)) ;
    for jj = 1:nU
        % Consider this hoop
        hoop = squeeze(sphi3d_rs(jj, :, :)) ;
        radii_from_avgpts(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
    end

    % Triangulate the sphigrid and store as its own cutMesh
    sv = ringpath_ss .* onesUV ;

    % Define path pairs for tiling the (s,phi) cut mesh
    spcutP1 = 1:nU;
    spcutP2 = nU*nV - fliplr(0:(nU-1)) ;
    spcutMesh.pathPairs = [ spcutP1', spcutP2' ];
    spcutMesh.f = defineFacesRectilinearGrid(uv, nU, nV) ;
    spcutMesh.nU = nU ;
    spcutMesh.nV = nV ;
    
    % First resampling
    spcutMesh.v0 = uvgrid3d ;
    spcutMesh.v0rs_equal_dsdv = c3d_dsv ;
    % spcutMesh.vrs0 = ((rot * new3d')' + trans) * resolution ;
    % Define normals based on the original mesh normals
    spvn03d = interpolate2Dpts_3Dmesh(TF, TV2D, TVN3D, uphi) ;
    spvn03d = spvn03d ./ vecnorm(spvn03d, 2, 2) ;
    spcutMesh.vn0 = spvn03d ;
    spcutMesh.sphi0 = [sv(:), phiv(:)] ;
    spcutMesh.uphi0 = uphi ;
    % Note: uv has no direct relation with cutMesh, just a grid
    % for utility and reference, but it does have unequal 
    % spacing in u in anticipation of building sphi0 as a 
    % near perfect grid.
    spcutMesh.uv = uv ;  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECOND RESAMPLING
    % Make a new grid
    slin = linspace(0, max(spcutMesh.sphi0(:, 1)), nU) ;
    plin = linspace(0, 1, nV) ;
    [ss, pp] = meshgrid(slin, plin) ;
    % Push the endpoints on each boundary in by epsilon to
    % avoid NaNs
    eps = 1e-14 ;
    ss(:, 1) = eps ;
    ss(:, end) = ss(:, end) - eps ;
    % Transpose so that x increases with increasing index first
    ss = ss' ;
    pp = pp' ;
    sp = [ss(:), pp(:)] ;

    % Tile the spcutMesh
    tileCount = [2, 2] ;
    spcutMesh.u = spcutMesh.sphi0 ;
    spcutMesh.v = spcutMesh.v0 ;
    spcutMesh.vn = spcutMesh.vn0 ;
    [ faces, v2d, v3d, vn3d ] = tileAnnularCutMesh( spcutMesh, tileCount );
    spcutMesh = rmfield(spcutMesh, 'u') ;
    spcutMesh = rmfield(spcutMesh, 'v') ;
    spcutMesh = rmfield(spcutMesh, 'vn') ;
    % Create the 3d vertices for sphi 2D vertices' correspondence
    spv3d = interpolate2Dpts_3Dmesh(faces, v2d, v3d, sp) ;
    % also interpolate the normals
    spvn3d = interpolate2Dpts_3Dmesh(faces, v2d, vn3d, sp) ;
    spvn3d = spvn3d ./ vecnorm(spvn3d, 2, 2) ;

    % Define new faces for second rectilinear resampling
    % NOTE: not necessary since we already defined the topology
    % from the guess [sv(:), phiv(:)] stored as spcutMesh.sphi0
    % spcutMesh.f = defineFacesRectilinearGrid(sp, nU, nV) ;
    spcutMesh.sphi = sp ;
    spcutMesh.v = spv3d ;
    spcutMesh.vn = spvn3d ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spcutMesh.ringpath_ss = ringpath_ss ;
    spcutMesh.radii_from_mean_uniform_rs = radii_from_mean_uniform_rs ;  % from uniform DV sampling
    spcutMesh.radii_from_avgpts = radii_from_avgpts ;
    spcutMesh.mss = mss ;       % from uniform DV sampling, also stored in centerline
    spcutMesh.mcline = mcline ; % from uniform DV sampling, also stored in centerline
    spcutMesh.avgpts = avgpts ; % from uniform DV sampling, also stored in centerline
    spcutMesh.avgpts_ss = avgpts_ss ; % from uniform sampling, also stored in centerline

    % Define optimal isoareal Affine dilation factor in s
    % USED TO SIMPLY TAKE ar FROM CUTMESH, NOW RECOMPUTE?
    tmp = spcutMesh.sphi ;
    tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
    spcutMesh.ar = minimizeIsoarealAffineEnergy( spcutMesh.f, spcutMesh.v, tmp );
    % clearvars tmp
    % spcutMesh.ar = cutMesh.ar ;
    disp(['old: ' num2str(cutMesh.ar)])
    disp(['new: ' num2str(spcutMesh.ar)])
    %error('check here')

    % todo: check that u coords have not shifted upon
    % redefinition of sphi0 -> sphi

    % Save s,phi and their 3D embedding
    spcutMesh.phi0s = phi0s ;
    spcutMesh.phi0_fit = phi0_fit ;
    save(spcutMeshfn, 'spcutMesh') ;
else
    disp('Loading spcutMesh from disk...')
    load(spcutMeshfn, 'spcutMesh') ;
    % QS.currentData.IVloaded = false ;

    % Load new centerline
    fn = sprintf(clineDVhoopBase, tt) ;
    disp(['Loading new centerline from ' fn])
    load(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')
end
fprintf('Done with generating S,Phi coords \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
