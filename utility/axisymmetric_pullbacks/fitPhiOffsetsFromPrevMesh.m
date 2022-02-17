function [phi0_fit, phi0s] = fitPhiOffsetsFromPrevMesh(TF, TV2D, TV3D,...
    uspace, vvals, prev3d_sphi_dsdv, lowerbound_phi0, upperbound_phi0, ...
    save_fit, plotfn, smoothingMethod, phiOpts, patchOpts)
%FITPHIOFFSETSFROMPREVMESH(TF, TV2D, TV3D, uspace, vspace, prev3d_sphi, lowerbound, upperbound, save_im, plotfn) 
%   Fit the offset phi values to add to V in UV coords to minimize
%   difference in 3D between current embedding mesh and previous one. This
%   rotates the hoops of the sphicutMesh.
%   - called by QuapSlap.generateCurrentSPCutMesh()
%   - NEW COORDINATES: phi = v - phi0, where v is input vvals, phi is
%   output coordinates, and phi0 is the result of the optimization.
%
% Parameters
% ----------
% TF : #faces x 3 int array
%   connectivity list of the mesh vertices into triangulated faces
% TV2D : #vertices x 2 float array
%   the UV mapped coordinates of mesh vertices
% TV3D : #vertices x 3 float array
%   the embedding coordinates of mesh vertices
% uspace : nU float array
%   The values of u for each line of constant v in pullback space
% vvals : nV float array OR nU x nV float array as grid
%   If nV x 1 float array, the values of v for each line of constant u in 
%   pullback space, otherwise the values for the whole grid. 
%   These are the values that will be shifted around by the fitting.
%   This allows
%   you to pass either a linspace for v (independent of u) or a series of
%   v values, one array for each u value
% prev3d_sphi_dsdv : nU x nV x 3 float array
%   The 3D coordinates of the embedding for the reference timepoint
%   (previous timepoint, for ex) at the 2D locations given by uspace and 
%   vspace. Note that uspace is not used explicitly, only nU is used to 
%   extract the strips over which we iterate, minimizing for phi0 for each
%   strip. Equally spaced in ds and in dv.
% lowerbound_phi0 : float 
%   lower bound for the fit of phi (offset to v). Must be > -1 and < 1
% upperbound_phi0 : float
%   upper bound for the fit of phi (offset to v). Must be > -1 and < 1.
%   Also must be <= lowerbound_phi0 + 1 to avoid phase ambiguity
% save_im : bool
% plotfn : str
% save_phi0patch : bool
%   save a patch colored by the phi0 motion deduced from the difference
%   between previous and current DVhoop coordinates
% smoothingMethod : str specifier
% phiOpts : struct passed to phiOffsetsFromPrevMesh, with fields
%   preview = preview ;
%   avgpts = avgpts_pix ;
%   c3d_dsv = c3d_dsv_pix ;
%       
%   prev_avgpts_ss = prev_avgpts_ss_pix ;
%       centerline pathlength for each average hoop position for previous
%       timepoint
%   prev_avgpts = prev_avgpts_pix ; 
%       average position of each hoop for previous timepoint
% patchOpts : optional options struct
%   preview   : bool (optional) visualize the progress of phi0
%   patchImFn : str  (optional) save the progress of phi0 as texture 
%                                 image saved to this path
%   imfn_sp_prev : str (optional) path to previous timepoint's sp pullback
%   IV        : MxNxP array (optional) intensity data to use for patch 
%   axisorder : 3x1 int array between 1-3 (optional) axis order for IV wrt
%               the mesh coordinate frame
%   v3d : 
%   ringpath_ss : 
%   preview : bool
%       whether to show intermediate results
%   vvals4plot : same as vvals but for plotting. The purpose of including
%       this is that there is a sign convention with what to do with phi0
%       after it has been found. If finding optimal phi0_kk iteratively by
%       running this code several times, then need to update 
%       vvals4plot as vvals --> vvals - phi0, even though phi0 is updated
%       as phi0_kk --> phi0_kk + phi0.
%   Options : struct, passed to texture_patch_to_image
%   
% Returns
% -------
% phi0_fit : nU x 1 float array
%   the smoothed best-fit rotation angles, in units of the v coordinate,
%   bounded by (lowerbound, upperbound) 
% phi0s : nU x 1 float array 
%   the best-fit rotation angles, bounded by (lowerbound, upperbound) 
% 
% NPMitchell 2019

try
    assert(strcmp(smoothingMethod, 'none') ||...
        strcmp(smoothingMethod, 'savgol'))
catch
    error('smoothingMethod must be none or savgol')
end

% If a patchOpts struct is passed, use it to save patch image preview
if ~isempty(fieldnames(patchOpts))
    try
        resolution = patchOpts.resolution ;
        patchImFn = patchOpts.patchImFn ;
        patchImFnRes = patchOpts.patchImFnRes ;
        imfn_sp_prev = patchOpts.imfn_sp_prev ;
        IV = patchOpts.IV ;
        axisorder = patchOpts.axisorder ;
        ringpath_ss = patchOpts.ringpath_ss ;
        v3d = patchOpts.v3d ;
        vvals4plot = patchOpts.vvals4plot ;
        Options = patchOpts.Options ;
        if isfield(patchOpts, 'save_phi0patch')
            save_phi0patch = patchOpts.save_phi0patch ;
        else
            save_phi0patch = true ;
        end
    catch
        error(['save_phi0patch is presumed true given that patchOpts ',...
            'were passed, but some patchOpts fields are missing!'])
    end
else
    disp('No patchOpts supplied, so save_phi0patch is false')
    save_phi0patch = false ;
end

% if ~isempty(varargin)
%     if isa(varargin{i},'logical')
%         continue;
%     end
% 
%     if ~isempty(regexp(varargin{i},'^[]ace[Nn]ormals','match'))
%         facenormals = varargin{i+1} ;
%     end
% end


% Consider each value of u in turn
% Fit for phi0 such that v = phi + phi0, meaning new
% phi = v - phi0.
disp('Minimizing phi0s...')

% Using simple offset
[phi0s, residuals] = phiOffsetsFromPrevMesh(TF, TV2D, TV3D, ...
    uspace, vvals, prev3d_sphi_dsdv, lowerbound_phi0, upperbound_phi0, phiOpts) ;

% Using dilation and offset
% [phi0s, ccoeffs] = phiOffsetsFromPrevMeshWithDilation(TF, TV2D, TV3D, ...
%      uspace, vspace, prev3d_sphi, preview) ;
% phi0s = mod(phi0s, 2 * pi) ;

% Convert phi0s to a smooth polynomial phi(u)
if strcmp(smoothingMethod, 'savgol')
    % Smoothing parameters
    framelen = 11 ;  % must be odd
    polyorder = 2 ;
    % Low pass filter (Savitsky-Golay)
    % Note we ignore the variations in ds -->
    % (instead use du=constant) to do this fit
    phi0_fit = savgol(phi0s, polyorder, framelen)' ;
elseif strcmp(smoothingMethod, 'none')    
    phi0_fit = phi0s ;
end
% Fit the smoothed curve
% phicoeffs = polyfit(uspace, phi0_fit, 14) ;
% phi0_fit = polyval(phicoeffs, uspace);

% Plot the fit for reference
if save_fit 
    close all
    fig = figure('visible', 'off') ;
    plot(uspace, phi0s, '.'); hold on;
    % plot(uspace, phix, '--')
    plot(uspace, phi0_fit, '-')
    if strcmp(smoothingMethod, 'savgol')
        legend({'measured shifts', 'SG filter'})
    elseif strcmp(smoothingMethod, 'none')
        legend({'measured shifts', 'interpolation'})
    end
    xlabel('u')
    ylabel('\phi_0')
    title('Shift \phi(u)')
    saveas(fig, plotfn)
    close all
end

%% OPTIONAL: PLOT CHANGE IN PHI OVER TEXTURE PATCH IMAGES OF PREVIOUS AND CURRENT IMAGE
if save_phi0patch
    % First check that we have vgrid, not vspace 
    nU = length(uspace) ;
    if any(size(vvals4plot) == 1)
        nV = length(vvals4plot) ;
        vgrid = (vvals4plot .* ones(nV, nU))' ;
    else
        vgrid = vvals4plot ;
        nV = size(vgrid, 1) ;
    end
    onesUV = ones(nU, nV) ;
    uu = uspace .* onesUV ; % only used for faces definition
    
    % Load previous sp pullback image
    if exist(imfn_sp_prev, 'file')
        im0 = double(imread(imfn_sp_prev)) / 255.0 ;
    else
        try
            im0 = double(imread(sprintf(QS.fullFileBase.im_sp, QS.currentTime))) ;
        catch
            error('Previous timepoint not available')
        end
    end
    tmp = ringpath_ss .* onesUV ;
    svcutMesh.u(:, 1) = tmp(:) ;
    svcutMesh.u(:, 2) = vgrid(:) ;
    svcutMesh.v = v3d ;
    
    % Generate cutMesh face triangulation
    vvtmp = linspace(0, 1, nV) .* ones(nV, nU) ; % only used for faces definition
    svcutMesh.f = defineFacesRectilinearGrid([uu(:), vvtmp(:)], nU, nV) ;
    
    % Generate cutpath pairs
    svcutP1 = 1:nU ;
    svcutP2 = nU*nV - fliplr(0:(nU-1)) ;
    svcutMesh.pathPairs = [ svcutP1', svcutP2' ] ;

    % [phi0_fit, phi0s] = fitPhiOffsetsFromPrevPullback() ;

    %% Generate Tiled Orbifold Triangulation -------------------
    disp('fitPhiOffsetsFromPrevMesh: generating temporary pullback')
    tileCount = [1 1];  % how many above, how many below
    [ TFtmp, TV2Dtmp, TV3Dtmp ] = tileAnnularCutMesh( svcutMesh, tileCount );

    % Create texture image
    if any(isnan(TV2Dtmp))
        error('here -- check for NaNs in TV2D ')
    end
    if any(isnan(TV3Dtmp))
        error('here -- check for NaNs in TV3D ')
    end
    disp('fitPhiOffsetsFromPrevMesh: Options = ')
    disp(Options)
    disp('fitPhiOffsetsFromPrevMesh: patchOpts = ')
    disp(patchOpts)
    disp('fitPhiOffsetsFromPrevMesh: patchOpts.Options = ')
    disp(patchOpts.Options)
    size(patchOpts.IV)
    patchIm = texture_patch_to_image( TFtmp, TV2Dtmp, TFtmp, TV3Dtmp(:, axisorder), ...
        IV, Options );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save a preview of the intermediate frame with phi0pix as heatmap,
    % overlaying current and previous patchIm
    close all
    figure('visible', 'off')
    disp(['Saving patchIm to ' patchImFn])
    % convert (ss, phi) to (x, y)
    xx = ringpath_ss * size(patchIm, 2) / max(ringpath_ss) ;
    yy = 1:100:size(patchIm, 1) ;
    phi0grid = (phi0_fit .* ones(length(uspace), length(yy)))' ;
    tmp = cat(3, patchIm, 0.5 * (patchIm + im0), im0) ;
    opts.label = '$\phi_0$' ;
    opts.qsubsample = 2 ;
    opts.qscale = 1000 ;
    xyfstruct = struct('x', xx, 'y', yy) ;
    [~, ~, ~, ax, ~] = vectorFieldHeatPhaseOnImage(tmp, xyfstruct, ...
        0*phi0grid, phi0grid, max(abs(phi0_fit))*2, opts) ;
    titlestr =  ['Blue = $t - 1$, yellow is current. ', ...
        '$\phi_0$ is to be subtracted from $v$.' ] ;
    title(ax, titlestr, 'Interpreter', 'Latex')
    % F = getframe(gca);
    % Image = frame2im(F);
    % imwrite(Image, patchImFn)
    disp(['Saving phi0patch: ' patchImFn])
    saveas(gcf, patchImFn) 
    close all
    
    %% Plot heatmap of residual over the image
    cfield = zeros(nU, nV) ;
    new3d = zeros(size(prev3d_sphi_dsdv)) ;
    for qq = 1:nU
        % 
        if any(size(vvals) == 1)
            vqq = vvals(qq, :)' ;
        else
            vqq = vvals(qq, :)' ;
        end
        % Old hoop vertex locations
        prev3dvals = squeeze(prev3d_sphi_dsdv(qq, :, :)) ;
        
        % Interpolated vertex positions
        new3d(qq, :, :) = interpolate2Dpts_3Dmesh(TF, TV2D, ...
                TV3D, [uspace(qq) * ones(nV, 1), mod(vqq + phi0s(qq), 1)]);
            
        % Distance field (should match residuals, but doesn't...)
        cfield(qq, :) = vecnorm(squeeze(new3d(qq, :, :)) - prev3dvals, 2, 2) ;
    end
    close all
    set(gcf, 'visible', 'off')
    % Plot residuals in 3d
    residual3d = residuals ;
    residual3d(:, nV) = residuals(:, 1) ;
    newx = new3d(:, :, 1) ;
    newx = newx(:) ;
    newy = new3d(:, :, 2) ;
    newy = newy(:) ;
    newz = new3d(:, :, 3) ;
    newz = newz(:) ;
    ox = prev3d_sphi_dsdv(:, :, 1) ;
    ox = ox(:) ;
    oy = prev3d_sphi_dsdv(:, :, 2) ;
    oy = oy(:) ;
    oz = prev3d_sphi_dsdv(:, :, 3) ;
    oz = oz(:) ;
    
    % Check hoops for alignment
    % plot3(newx, newy, newz, '.'); hold on;
    % plot3(ox, oy, oz, '.')
    
    trisurf(svcutMesh.f, ox, oy, oz, residual3d(:) * resolution,...
        'EdgeColor', 'none') ;
    hold on;
    scatter3(newx, newy, newz, 2, residual3d(:) * resolution)
    axis equal
    xlabel('AP position [\mum]'); 
    ylabel('lateral position [\mum]'); 
    zlabel('DV position [\mum]')
    title('Residual from geometric hoop alignment')
    cb = colorbar() ;
    ylabel(cb, 'residual [\mum]')
    caxis([0, 10])
    % Plot residuals in 2d
    % options_residual.flipy = false ;
    % options_residual.alpha = 0.3 ;
    % [fig, ax1, ax2, imhandle, shandle] = ...
    %     heatmap_on_image(tmp, xx, yy', cfield, options_residual)
    % title('Residual from hoop matching')
    disp(['Saving phi0patch: ' patchImFnRes])
    saveas(gcf, patchImFnRes) 
    close all
    
end