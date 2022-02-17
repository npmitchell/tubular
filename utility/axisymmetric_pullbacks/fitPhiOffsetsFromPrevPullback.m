function [phi0_fit, phi0s] = fitPhiOffsetsFromPrevPullback(IV,...
    new3d, umax, uspace, vgrid, ringpath_ss, imfn_sp_prev, ...
    lowerboundy, upperboundy, save_ims, plotfn, Options, step_phi0tile,...
    width_phi0tile, potential_sigmay, method, patchImFn)
%FITPHIOFFSETSFROMPREVPULLBACK(TF, TV2D, TV3D, uspace, vspace, prev3d_sphi, lowerbound, upperbound, save_im, plotfn, save_patchIm) 
%   Fit the offset phi values to add to V in UV coords to minimize
%   difference in 3D between current embedding mesh and previous one. This
%   rotates the hoops of the sphicutMesh. Texture method.
%   Smoothing is built-in to the function (hard-coded).
%
% Parameters
% ----------
% IV : 3d float values
%   intensity values of the raw data, in which the 3d mesh vertices reside
% new3d : nU*nV x 3 float array
%   the embedding coordinates of mesh vertices
% uspace : nU float array
%   The values of u for each line of constant v in pullback space
% vspace : nV float array OR nU x nV float array as grid
%   If nV x 1 float array, the values of v for each line of constant u in 
% imfn_sp_prev : str
%   full path to previous timepoint's pullback
% lowerbound : float 
%   lower bound for the fit of phi (offset to v)
% upperbound : float
%   upper bound for the fit of phi (offs
% save_im 
% plotfn
% Options : struct 
%   Options for the pullback generation to compare
% step : int
%   How far apart each column to correlate, distance in pixels in x dir
% width : int
%   How wide in pixels for each column to correlate
% method : str ('subpixel', 'integer')
%   Whether to use subpixel or integer resolution in determining phase
%   correlation peaks
% patchImFn : str
%   If not 'None', save the intermediate patchIm before optimization
%
% Returns
% -------
% phi0_fit
% phi0s : nU x 1 float array 
%   the additional rotation angles, bounded by (lowerbound, upperbound) 
% 
% NPMitchell 2020

% build svcutMesh
nU = length(uspace) ;
nV = size(vgrid, 1) ;
onesUV = ones(nU, nV) ;
uu = uspace .* onesUV ; % only used for faces generation

% Generate the connectivity list
% vv = (vspace .* ones(nV, nU))' ;
vvtmp = linspace(0, 1, nV) .* ones(nV, nU) ; % only used for faces generation
svcutMesh.f = defineFacesRectilinearGrid([uu(:), vvtmp(:)], nU, nV) ;

% Generate cutpath pairs
svcutP1 = 1:nU ;
svcutP2 = nU*nV - fliplr(0:(nU-1)) ;
svcutMesh.pathPairs = [ svcutP1', svcutP2' ] ;

% Load previous sp pullback image
if exist(imfn_sp_prev, 'file')
    im0 = double(imread(imfn_sp_prev)) / 255.0 ;
else
    error('Previous timepoint not available')
end
tmp = ringpath_ss .* onesUV ;
svcutMesh.u(:, 1) = tmp(:) ;
svcutMesh.u(:, 2) = vgrid(:) ;
svcutMesh.v = new3d ;
% [phi0_fit, phi0s] = fitPhiOffsetsFromPrevPullback() ;

%% Generate Tiled Orbifold Triangulation -------------------
disp('fitPhiOffsetsFromPrevPullback: generating temporary pullback')
tileCount = [1 1];  % how many above, how many below
[ TF, TV2D, TV3D ] = tileAnnularCutMesh( svcutMesh, tileCount );

% Create texture image
if any(isnan(TV2D))
    error('here -- check for NaNs case')
end
patchIm = texture_patch_to_image( TF, TV2D, TF, TV3D(:, [2 1 3]), ...
    IV, Options );
ysize0 = size(patchIm, 1) ;
xsize0 = size(patchIm, 2) ;

% Blue the image to increase correlation smoothness
patchIm = imgaussfilt(patchIm, 2) ;


%% Compare patchIm to im0, chopped into columns, to find phi0(u)
% Consider each vertical slice of the image
% Consider a dense sampling in the X/u direction
disp('fitPhiOffsetsFromPrevPullback: comparing to previous pullback')
if step_phi0tile < width_phi0tile
    Xtodo = 0:step_phi0tile:(xsize0+step_phi0tile) ;
else
    Xtodo = 0:step_phi0tile:(xsize0+width_phi0tile) ;
end    
phi0pix = zeros(length(Xtodo), 1) ;
dXpix = zeros(length(Xtodo), 1) ;
dYpix = dXpix ;
dxint = dXpix ;
dyint = dXpix ;
potential = 'none' ;
dim = 2 ;
for qq=1:length(Xtodo)
    xqq = Xtodo(qq) ;
    if mod(xqq, 50) == 0
        msg = ['u pixel index ' num2str(xqq) '/' num2str(xsize0) ] ;
        msg = [msg ' (' num2str(qq) '/' num2str(length(Xtodo)) ')'] ;
        disp(msg)
    end
    % Grab extra chunk of image if we are examining the image edge
    if qq == 1 || qq == length(Xtodo)
        xwidth = width_phi0tile * 2 ;
    end
    
    % grab the slice of current image
    a = patchIm(:, max(1, xqq-xwidth):min(xsize0, xqq+xwidth)) ;

    % grab the same slice of the previous image
    b = im0(:, max(1, xqq-xwidth):min(xsize0, xqq+xwidth)) ;

    % Brute force correlation
    % -----------------------
    % % Consider a range of shifts
    % yvals = -yrange:yrange ;
    % % preallocate the corrs
    % corrs = zeros(length(yvals), 1) ;
    % for rr = 1:length(yvals) 
    %     ashift = circshift(a, yvals(rr), 1) ;
    %     corrs(rr) = corr2(ashift, b);
    % end
    % corrs = smoothdata(corrs, 'gaussian', 100) ; 
    % [~, id] = max(corrs) ;
    % phi0pix(qq) = yvals(id) ;
    %
    % % preview results
    % if preview || true
    %     plot(yvals, corrs)
    %     xlabel('yshift [pix]')
    %     ylabel('correlation')
    %     hold on;
    %     plot([dY, dY], ylim(), '--')
    %     pause(0.001)
    % end

    % Use phase correlation
    if strcmp('method', 'subpixel')
        % USE SUBPIXEL RESOLUTION
        [dXpix(qq), dYpix(qq), dxint(qq), dyint(qq)] = ExtPhaseCorrelationPotentialBound(a, b, ...
            potential, dim, ...
            'lowerboundy', lowerboundy, 'upperboundy', upperboundy, ...
            'lowerboundx', -width_phi0tile*0.5, 'upperboundx', width_phi0tile*0.5, ...
            'PotentialSigmaY', potential_sigmay, 'SigmaFilter', 3) ;
        phi0pix(qq) = dYpix(qq) ;
    else
        % USE INTEGER RESOLUTION
        [dxint(qq), dyint(qq)] = PhaseCorrelationPotentialBound(a, b, ...
            potential, dim, ...
            'lowerboundy', lowerboundy, 'upperboundy', upperboundy, ...
            'lowerboundx', -width_phi0tile*0.5, 'upperboundx', width_phi0tile*0.5, ...
            'PotentialSigmaY', potential_sigmay, 'SigmaFilter', 3) ; %, 'verbose', true) ;
        phi0pix(qq) = dyint(qq) ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save a preview of the intermediate frame with phi0pix as heatmap,
% overlaying current and previous patchIm
if isempty(regexp(patchImFn, '^[Nn]one','match'))
    close all
    figure('visible', 'off')
    disp(['Saving patchIm to ' patchImFn])
    % imwrite(patchIm, patchImFn)
    yy = 1:100:size(patchIm, 1) ;
    phi0pixgrid = (phi0pix .* ones(length(Xtodo), length(yy)))' ;
    dXpixgrid = (dXpix .* ones(length(Xtodo), length(yy)))' ;
    tmp = cat(3, patchIm, patchIm, im0) ;
    opts.label = '$\phi_0$' ;
    opts.qsubsample = 2 ;
    opts.qscale = 1 ;
    [~, ~, ~, ax, ~] = vectorFieldHeatPhaseOnImage(tmp, Xtodo', yy', dXpixgrid, phi0pixgrid, max(abs(phi0pix))*2, opts) ;
    title(ax, 'Blue is prev timept, yellow is current')
    % F = getframe(gca);
    % Image = frame2im(F);
    % imwrite(Image, patchImFn)
    saveas(gcf, patchImFn) 
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Replace NaNs with nearby values
phi0pix = fillmissing(phi0pix, 'nearest') ;
dXpix = fillmissing(dXpix, 'nearest') ;
dYpix = fillmissing(dYpix, 'nearest') ;
dxint = fillmissing(dxint, 'nearest') ;
dyint = fillmissing(dyint, 'nearest') ;

% Relaxing assertions for now to allow for NaNs
if ~all(phi0pix > lowerboundy)
    disp(phi0pix)
    error('phi0pix out of bounds: below lower bound')
end
if ~all(phi0pix < upperboundy)
    disp(phi0pix)
    error('phi0pix out of bounds: above upper bound')
end

% Take moving median of the peaks
phi0pix2 = smoothdata(phi0pix, 'movmedian', 5) ;
% fit a curve to the result
% Convert phi0s to a smooth polynomial phi(u)
% Smoothing parameters
framelen = 7 ;  % must be odd
polyorder = 2 ;
% Low pass filter (Savitsky-Golay)
% Note we ignore the variations in ds -->
% (instead use du=constant) to do this fit
phi0s = savgol(phi0pix2, polyorder, framelen) ;

% resample phi0_fit so that it is the same length as
% ringpath_ss
samplex = Xtodo / xsize0 ;


% Could do this as a polynomial....
p = polyfit(samplex - 0.5, phi0s, 10) ;
% sample at ringpath_ss, note resize domain to (0,1)
phi0_fit = polyval(p, (ringpath_ss / max(ringpath_ss) - 0.5)) / ysize0 ;

% Or as interp1d on savgol curve and normalize to v = (0, 1)
% phi0_fit = interp1(samplex, phi0pix2, ringpath_ss / max(ringpath_ss), 'cubic') ;
% phi0_fit = phi0_fit / ysize0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE A PLOT OF THE PHI DETERMINATION AND FITTING
if save_ims
    if strcmp(method, 'subpixel')
        % view results (subpixel resolution)
        close all
        fig = figure('visible', 'off') ;
        subplot(2, 1, 1)
        plot(Xtodo, dYpix, 'o');
        hold on;
        plot(Xtodo, dyint, 's') ;
        plot(Xtodo, phi0pix2) ;
        plot(Xtodo, polyval(p, samplex- 0.5)) ;
        plot(Xtodo, dXpix, '.--') ;
        plot(Xtodo, dxint, '--') ;
        ylims = ylim() ;
        ylim0 = max(ylims(1), lowerboundy) ;
        ylim1 = min(ylims(2), upperboundy) ;
        ylim([ylim0, ylim1])
        xlabel('AP position [pullback pixels]')
        ylabel('Shift [pixels]')
        legend({'subpixel', 'integer', 'smoothed int',  'fit', ...
            'dX subpix', 'dX int'}, ...
            'location', 'northeastoutside')
        title(['Texture matching shifts in image: ' method ' method'])
        subplot(2, 1, 2)
        plot(ringpath_ss, phi0_fit) ;
        xlabel('AP position [\mum]')
        ylabel('Shift, \Delta\phi/\pi')
        disp(['Saving fig: ' plotfn])
        saveas(fig, plotfn)
    else
        % view results (integer resolution)
        close all
        fig = figure('visible', 'off') ;
        subplot(2, 1, 1)
        hold on;
        plot(Xtodo, dyint, 's') ;
        plot(Xtodo, phi0pix2, '.') ;
        plot(Xtodo, phi0s) ;
        % plot(Xtodo, polyval(p, samplex - 0.5), '--') ;
        plot(Xtodo, dxint, '--') ;
        ylims = ylim() ;
        ylim0 = max(ylims(1), lowerboundy) ;
        ylim1 = min(ylims(2), upperboundy) ;
        ylim([ylim0, ylim1])
        xlabel('AP position [pullback pixels]')
        ylabel('Shift [pixels]')
        legend({'integer corr', 'smoothed int', 'Sav-Gol', 'fit', ...
            'dX int'}, ...
            'location', 'northeastoutside')
        title(['Texture matching shifts in image: ' method ' method'])
        subplot(2, 1, 2)
        plot(ringpath_ss, phi0_fit) ;
        xlabel('AP position [\mum]')
        ylabel('Shift, \Delta\phi/\pi')
        disp(['Saving fig: ' plotfn])
        saveas(fig, plotfn)
    end
end
