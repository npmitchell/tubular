function spcutMesh = generateCurrentSPCutMesh(tubi, cutMesh, spcutMeshOptions)
% generateCurrentSPCutMesh(QS, cutMesh, spcutMeshOptions)
%
% Compute the parameterization (s(u'),phi) where s is the
% circumferentially-averaged pathlenth along the surface from u=0 to u' and
% phi = v + phi0(u) is the cicumferential coordinate twisted by an amount
% phi0 relative to the uv conformal parameterization. phi0 is chosen for
% each discrete value of u based on either:
%   (1) [spcutMeshOptions.phiMethod=='3dcurves']
%       the geometric position of the circumferential curve in 3d space
%       relative to the previous timepoint (or next timepoint if tubi.t0 >
%       t, where t is the current timepoint in question). In other words,
%       we rotate each slice of the tube to more closely match the
%       geometric position of the analogous slice in a timepoint closer to
%       t0. 
%   (2) [spcutMeshOptions.phiMethod=='texture']
%       the optical correspondence of each slice in pullback space with a
%       'nearby' slice in pullback space of a previous timepoint (later
%       timepoint if t<tubi.t0). In other words, we match the material
%       position in the circumferential coordinate of a timepoint closer to
%       t0
%   (2) [spcutMeshOptions.phiMethod=='combined']
%       first based on the geometric position of the circumferential curve 
%       in 3d space relative to the previous timepoint 
%       (or next timepoint if t0 >t, where t is the current timepoint in 
%       question), THEN additionally apply an optical matching. This is
%       useful if there is a lot of jittery motion of the tissue. Note that
%       you can pass 
%       In other words, we rotate each slice of the tube to more closely 
%       match the geometric position of the analogous slice in a timepoint
%       closer to t0. 
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
% tubi : TubULAR class instance
%   Note that the following properties are used:
%       tubi.phiMethod = ('3dcurves', 'texture', 'combined') 
%       tubi.a_fixed = 2.0    
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
%   smoothingWidth : int 
%       width of kernel for smoothing of phi0 that takes v->phi=v-phi0.
%       Must be odd if smoothingMethod=='savgol', does not matter 
%       if smoothingMethod=='none'
%   phi0TextureOpts : struct with fields
%        lowerboundy : float
%        lowerboundy : float
%        step_phi0tile : float
%        width_phi0tile : float 
%        potential_sigmay : 
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
textureAxisOrder = tubi.data.axisOrder ;
save_ims = tubi.plotting.save_ims ;
phi0_sign = tubi.flipy ;  % NOTE: our convention is that new phi = v - phi_0 but for
                        % some reason when flipy is true even though we have
                        % flipped cutMesh, we must add a minus sign so that phiv_kk
                        % becomes phi_kk = v + phi0. 
smoothingWidth = 11 ;  % must be odd if smoothingMethod=='savgol', does not matter if smoothingMethod=='none'
smoothingOrder = 2 ;  % used if smoothingMethod=='savgol', polynomial order of filter

%% Unpack options
cutMesh = tubi.getCurrentCutMesh() ;
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
            t0_for_phi0 = tubi.t0set() ;
            disp(['Setting t0 = ' num2str(t0_for_phi0)])
        catch
            disp('Setting t0 = first timepoint')
            t0_for_phi0 = tubi.xp.fileMeta.timePoints(1) ; 
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
    if isfield(spcutMeshOptions, 'smoothingWidth')
        smoothingWidth = spcutMeshOptions.smoothingWidth ;
    end
    if isfield(spcutMeshOptions, 'smoothingOrder')
        smoothingOrder = spcutMeshOptions.smoothingOrder ;
    end
    if isfield(spcutMeshOptions, 'phi0TextureOpts')
         phi0TextureOpts = spcutMeshOptions.phi0TextureOpts ;
    end
end

% populate patchOpts from pbOptions
pbOptions = struct() ;
pbOptions.resolution = tubi.xp.fileMeta.stackResolution(1) ;
pbOptions.axisorder = textureAxisOrder ;

%% Unpack tubi
tt = tubi.currentTime ;
nU = tubi.nU ;
nV = tubi.nV ;
a_fixed = tubi.a_fixed ;
phi_method = tubi.phiMethod ;
spcutMeshfn = sprintfm(tubi.fullFileBase.spcutMesh, tt) ;
fileNameBase = tubi.fileBase.name ; 
phi0fitBase = tubi.fullFileBase.phi0fit ;
spcutMeshBase = tubi.fullFileBase.spcutMesh ;
[rot, trans] = getRotTrans(tubi) ;
resolution = tubi.APDV.resolution ;
tubi.getCleanFastMarchingCenterlines ;
cleanCntrlines = tubi.cleanFMCenterlines ;
cleanCntrline = cleanCntrlines{tubi.xp.tIdx(tt)} ;
preview = tubi.plotting.preview ;
[~, ~, xyzlim_um] = tubi.getXYZLims() ;
clineDVhoopDir = tubi.dir.clineDVhoop ;
clineDVhoopBase = tubi.fullFileBase.clineDVhoop ;
sphiDir = tubi.dir.spcutMesh ;

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
    % First initialize cutMeshrs as a full copy, but we will rotate and
    % scale the vertices and normals (rs = rotated and scaled)
    cutMeshrs = cutMesh;
    
    % DEBUG 12-03-2020
    % Rotate and translate TV3D
    cutMeshrs.v = tubi.xyz2APDV(cutMesh.v) ;
    
    % DEBUG 12-03-2020 -- added flipping vn here
    cutMeshrs.vn = (rot * cutMesh.vn')' ;
    if tubi.flipy
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
    cline = tubi.xyz2APDV(cleanCntrline(:, 2:4)) ; 
    ss = cleanCntrline(:, 1) * tubi.APDV.resolution ;
    disp('Finding relevant segment of centerline from the raw centerline')
    % Find the start and endpoints of a centerline closest to the 3D
    %   positions of the boundaries corresponding to u=0 and u=1 of the 2D
    %   mesh pullback representation.
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
    [crude_ringpath_ss, ~, uvgrid3d, uvgrid2d] = ...
        ringpathsGridSampling(uspace0, vspace, TF, TV2D, TV3Drs) ;

    % Resample crude_ringpath_ds made from uspace0 (uspace0 had equal du, not equal ds_3D in u direction)
    [uspace_ds, eq_ringpath_ss] = equidistantSampling1D(linspace(0, 1, nU)', crude_ringpath_ss, nU, 'linear') ;
    % ensure that uspace is nU x 1, not 1 x nU
    uspace_ds = reshape(uspace_ds, [nU, 1]) ; 
    % if we wish to squish the posterior end of the image as we go from uv 
    % to sphi then we might have:
    %   ^           / uspace_ds
    %   |          /
    %   |         /
    %   |      __/
    %   |   __/
    %   |__/
    %   |---------------> index
    % and crude_ringpath_ss would be a bit like its inverse but with units
    % of tubi.spaceUnits:
    %   ^         __/ crude_ringpath_ss. Note this is 'crude' since the
    %   |      __/         hoops are not evenly sampled in the computation.
    %   |   __/   
    %   |  / 
    %   | /
    %   |/
    %   |---------------> index
    % 
    
    
    % hedge the first and last point to avoid NaNs
    eps = 1e-13 ;
    uspace_ds(1) = uspace_ds(1) + eps ;
    uspace_ds(end) = uspace_ds(end) - eps ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Casting resampled points into 3D (approx equal ', ...
        'ds_3D in u dir, but variable ds_3D in v dir)...'])
    % NOTE: first dimension indexes u, second indexes v
    curves3d = zeros(nU, nV, 3) ;  % in units of um
    for kk = 1:nU
        if mod(kk, 20) == 0
            disp(['u = ' num2str(kk / nU)])
        end
        % here sample the mesh with equally spaced ds
        uv_ds = [cutMesh.umax * uspace_ds(kk) * ones(size(vspace)), vspace] ;
        curves3d(kk, :, :) = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, uv_ds) ;
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

    % The original centerline with endcaps truncated is cseg_ss for
    % comparison. cseg = centerline segment, ss = integrated pathlength
    % along original centerline, bu
    cseg_ss = ss(acID:pcID) ;
    % [ssv, radii, avgpts, cids] = srFromDVCurves(cseg_ss, cseg, c3ds) ;
    % 
    % Adjust the centerline indices to index into the full
    % centerline. Note that cseg_ss already does this for ss.
    % cids = cids + acID ;

    % Plot new centerline
    aux_plot_clineDVhoop(tubi, avgpts, avgpts_ss, cseg, cline, cseg_ss, ...
        curves3d, xyzlim_um, clineDVhoopFigBase, tt)
    
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
    fn = sprintfm(clineDVhoopBase, tt) ;
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
    uu_ds_max = uspace_ds_umax .* onesUV ;
    vv = (vspace .* onesUV')' ;
    % uv is equally spaced in ss, UNEQUALLY spaced in uu
    sv = [uu_ds_max(:), vv(:)] ;
    % Note: here interpolate uv in the TV2D coord system, then
    % use uphi as the actual 2D coordinates for these vertices
    % NOTE: unlike curves3d, uvgrid3d is NOT rotated/translated/scaled
    sv0grid3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, sv) ;
    
    % Here, sv0grid3d should be identical to curves3d but in pixel space
    assert(all(vecnorm(sv0grid3d - ...
        tubi.APDV2xyz( reshape(curves3d, [nU*nV, 3]) ), 2, 2) < 1e-10))
    
    %% Make avgpts in pixel space (not RS)
    fprintf('Resampling uvgrid3d curves in pix...\n')
    curves3d_pix = reshape(sv0grid3d, [nU, nV, 3]) ;
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
        tidx = tubi.xp.tIdx(tt) ;
        tp_for_comparison = tubi.xp.fileMeta.timePoints(tidx-1) ;
        compute_phi0 = true ;
    elseif tt < t0_for_phi0 
        tidx = tubi.xp.tIdx(tt) ;
        tp_for_comparison = tubi.xp.fileMeta.timePoints(tidx+1) ;
        compute_phi0 = true ;
    end
    
    if compute_phi0
        % Load previous sphi vertices in 3d 
        if strcmp(phi_method, '3dcurves') || strcmp(phi_method, 'combined') 
            disp('Computing phi(v) via 3dcurve matching (geometric method)')
            % Load the previous spcutMesh and call it prev3d_sphi
            % Also note the previous spcutMesh pullback image's fn
            tmp = load(sprintfm(spcutMeshBase, ...
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
            % fn_prev_mcline = sprintfm(clineDVhoopBase, tt - 1) ;
            disp('Grabbing previous centerline')
            % prev_avgpts = load(fn_prev_mcline, 'avgpts') ;
            % prev_avgpts = prev_avgpts.avgpts ;
            [~, ~, ~, prev_avgpts_ss_pix, prev_avgpts_pix] = ...
                srFromDVCurves(prev3d_sphi_dsdv) ;
            
            %% Obtain previous pullback image
            % prev2d_uphi = reshape(tmp.spcutMesh.uphi, [nU, nV, 2]) ;
            imfn_sp_prev = sprintfm( tubi.fullFileBase.im_sp, tp_for_comparison) ;

            % fit the shifts in the v direction as a curv along the u
            % direction.
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
                disp(['Fitting phi0 to minimize motion in 3D of parameterization: Iteration ' num2str(dmyk)])
                plotfn = sprintfm(phi0fitBase, tt, dmyk);

                % Pass phiOptions to allow sliding along AP
                % passed to phiOffsetsFromPrevMesh
                phiOpts = struct(); 
                phiOpts.preview = preview ;
                phiOpts.avgpts = avgpts_pix ;
                phiOpts.c3d_dsv = c3d_dsv_pix ;
                phiOpts.prev_avgpts_ss = prev_avgpts_ss_pix ;
                phiOpts.prev_avgpts = prev_avgpts_pix ;
                phiOpts.smoothingWidth = smoothingWidth ;
                phiOpts.smoothingOrder = smoothingOrder ;
                
                % Will we save check pullbacks to preview the algo?
                % If so, create a struct to pass visualization options
                if save_phi0patch
                    if dmyk == 0 
                        % first pass --> vvals4plot is same as vvals
                        phi4plot = (vspace .* ones(nU, nV))' ;
                    end
                    patchImFn = sprintfm( ...
                        fullfile(sphiDir, 'phi0_correction',...
                        [fileNameBase, '_prephi0_' num2str(dmyk) '.tif']), ...
                        tp_for_comparison)  ;
                    patchImFnRes = sprintfm( ...
                        fullfile(sphiDir, 'phi0_correction',...
                        [fileNameBase, '_phi0residual_' num2str(dmyk) '.tif']), ...
                        tp_for_comparison)  ;
                    geomImFn = sprintfm( ...
                        fullfile(sphiDir, 'phi0_correction', ...
                        ['3d' fileNameBase '_prephi0_' num2str(dmyk) '.tif']), ...
                        tp_for_comparison)  ;

                    % Load the intensity data for this timepoint
                    if isempty(tubi.currentData.IV)
                        tubi.getCurrentData()
                    end
                    IV = tubi.currentData.IV ;
                    
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
                    patchOpts.v3d = sv0grid3d ;
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
                    pe0h = plot3(sv0grid3d(pe0, 1), sv0grid3d(pe0, 2), ...
                        sv0grid3d(pe0, 3), '.') ;
                    % after fitting
                    hold on;                   
                    pekk = find(phi4plot(:) < 1e-4 | phi4plot(:) > 0.99) ;
                    pekkh = plot3(sv0grid3d(pekk, 1), sv0grid3d(pekk, 2), ...
                        sv0grid3d(pekk, 3), '.') ;         
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
            phi0_fit = tubi.fitPhiOffsetsViaTexture(uspace_ds_umax, vspace, [], ...
                phi0TextureOpts) ;
        else
            error(["Could not recognize phi_method: ", ...
                "must be 'texture' or '3dcurves' or 'combined'"])
        end
        
        % If we use a combined method, use curves3d as initial guess for
        % texture method
        if strcmp(phi_method, 'combined') 
            disp('Refining phi(v) via texture matching (physical method)')
            phi0_fit = tubi.fitPhiOffsetsViaTexture(uspace_ds_umax, vspace,...
                        phi0_fit, phi0TextureOpts) ;
        end
        close all

        % Store to save at this timepoint
        % phiv = (vspace .* ones(nU, nV))' + phi0_fit .* ones(nU, nV) ;
        % phiv = phiv_kk ;
        phiv = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;    
                
    end

    % NOTE: We have coordinates u,phiv that we associate with
    % the 3d coordinates previously mapped to uv
    uphi = [uu_ds_max(:), phiv(:)] ;

    % Recompute radii_from_mean_uniform_rs as radii_from_avgpts 
    % NOTE: all radius calculations done in microns, not pixels
    sphi3d_rs = ((rot * sv0grid3d')' + trans) * resolution ;
    if tubi.flipy
        sphi3d_rs(:, 2) = - sphi3d_rs(:, 2) ;
    end
    radii_from_avgpts = zeros(size(sphi3d_rs, 1), size(sphi3d_rs, 2)) ;
    for jj = 1:nU
        % Consider this hoop
        hoop = squeeze(sphi3d_rs(jj, :, :)) ;
        radii_from_avgpts(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
    end

    % Triangulate the sphigrid and store as its own cutMesh
    ssv = ringpath_ss .* onesUV ;

    % Define path pairs for tiling the (s,phi) cut mesh
    spcutP1 = 1:nU;
    spcutP2 = nU*nV - fliplr(0:(nU-1)) ;
    spcutMesh.pathPairs = [ spcutP1', spcutP2' ];
    spcutMesh.f = defineFacesRectilinearGrid(sv, nU, nV) ;
    spcutMesh.nU = nU ;
    spcutMesh.nV = nV ;
    
    % First resampling
    spcutMesh.v0 = sv0grid3d ;
    spcutMesh.v0rs_equal_dsdv = c3d_dsv ;
    % spcutMesh.vrs0 = ((rot * new3d')' + trans) * resolution ;
    % Define normals based on the original mesh normals
    spvn03d = interpolate2Dpts_3Dmesh(TF, TV2D, TVN3D, uphi) ;
    spvn03d = spvn03d ./ vecnorm(spvn03d, 2, 2) ;
    spcutMesh.vn0 = spvn03d ;
    spcutMesh.sphi0 = [ssv(:), phiv(:)] ;
    spcutMesh.uphi0 = uphi ;
    % Note: uv0 has no direct relation with cutMesh, just a grid
    % for utility and reference. It has unequal 
    % spacing in u (approx same as s) in anticipation of building sphi0 
    % Note in particular that the uv
    spcutMesh.uv0 = sv ;  

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
    % Define fields u,v, and vn just for tiling purposes
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
    spcutMesh.mss = mss ;       % 1000x1 float, pathlength of centerline from avgpts but finely sampled, from uniform DV sampling, also stored in centerline
    spcutMesh.mcline = mcline ; % 1000x3 float, centerline from avgpts but finely sampled, from uniform DV sampling, also stored in centerline
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
    spcutMesh.readme = struct('readme', ['principal fields are f, v, sphi for the mesh in 3d and 2d. ',...
        'vn are vertex normals. nU and nV are dimensions of the gridded pullback vertices.', ...
        'ringpath_ss, mss, avgpts, avgpts_ss'], ...
        'f', '#faces x 3 int, face connectivity list indexing in to vertices v0 (for uv mesh) or v (for sphi mesh)', ...
        'phi0s', 'nU x 1 float, values of phi0 which offset the pullback vertices along the second dimension', ...
        'phi0_fit', 'nU x 1 float, fitted values of phi0, which may be equal to phi0s if there is no smoothing applied',...
        'ar', 'float, aspect ratio of the pullback which minimizes IsoarealAffineEnergy. An aspect ratio of 1 is most conformal, but an aspect ratio of ar sometimes is helpful for visualization', ...
        'pathPairs', 'nUx2 int, indices of vertices which are identical across the periodic direciton', ...
        'nU', 'int, number of points along the longitudinal axis in the resampled grid uv or sphi', ...
        'nV', 'int, number of points along the circumferential axis in the resampled grid uv or sphi', ...
        'v0', 'nU*nV x 3 float, mesh embedding of uv vertices', ...
        'v0rs_equal_dsdv', 'nU*nV x 3 float, mesh embedding with vertices equidistant along the circumferential direction. Note there is no corresponding pullback grid of vertices for this mesh!', ...
        'vn0', 'nU*nV x 3 float, normal vectors at the vertices of the sphi0 vertices', ...
        'sphi0', 'advected nU*nV x 2, (s,phi) coordinates of v0 -- [ssv(:), phiv(:)] -- ie before second resampling into a grid in pullback space ',...
        'uphi0', 'nU*nV x 2 float, [uu_ds_max(:), phiv(:)] -- ie positions in (s,phi) space of uv grid', ...
        'uv0', 'nU*nV x 2 float, [uu_ds_max(:), v(:)] -- ie positions in (s,v) space of uv grid', ...
        'sphi', 'nU*nV x 2 float, 2D pullback positions of sphi coordinates (resampled as grid in s,phi space), corresponds to v and vn',...
        'v', 'nU*nV x 3 float, 3D vertex embedding positions',...
        'vn', 'nU*nV x 3 float, vertex normals of the sphicutMesh embedding', ...
        'ringpath_ss', 'nU x 1 float, pathlength along the centerline computed from average position of uniform-circumferential-pathlength-sampling of DV hoops', ...
        'radii_from_mean_uniform_rs', 'nU*nVx1 float, distance from a uniform-circumferential-pathlength-sampling of DV hoops to the avgpts centerline', ...
        'radii_from_avgpts', 'nU*nVx1 float, distance from the avgpts centerline to each vertex on the surface', ...
        'mss', '1000x3 float, pathlength along centerline from avgpts but finely sampled, from uniform DV sampling, also stored in centerline', ...
        'mcline', '1000x3 float, centerline from avgpts but finely sampled, from uniform DV sampling, also stored in centerline', ...
        'avgpts', 'from uniform DV sampling, also stored in centerline', ...
        'avgpts_ss', 'nU x 1 float from uniform sampling, also stored in centerline') ;
    
    save(spcutMeshfn, 'spcutMesh') ;
else
    disp('Loading spcutMesh from disk...')
    load(spcutMeshfn, 'spcutMesh') ;
    % QS.currentData.IVloaded = false ;

    % Load new centerline
    fn = sprintfm(clineDVhoopBase, tt) ;
    disp(['Loading new centerline from ' fn])
    load(fn, 'mss', 'mcline', 'avgpts', 'avgpts_ss')
end
fprintf('Done with generating S,Phi coords \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
