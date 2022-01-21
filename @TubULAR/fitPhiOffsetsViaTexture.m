function phi0_fit = fitPhiOffsetsViaTexture(QS, ...
    uspace_ds_umax, vspace, phi0_init, phi0TextureOpts)
% fitPhiOffsetsViaTexture(QS)
%   fit offsets in the DV direction of pullback mesh by measuring phase
%   correlation between current and previous frame
%
% NPMitchell 2020

% Default parameters for texturepatch phase correlation computation of phi0
lowerboundy = -350 ;
upperboundy = 350 ;
step_phi0tile = 25 ;
width_phi0tile = 150 ;
potential_sigmay = 350 ;
if nargin > 4
    lowerboundy = phi0TextureOpts.lowerboundy ;
    upperboundy = phi0TextureOpts.upperboundy ;
    step_phi0tile = phi0TextureOpts.step_phi0tile ;
    width_phi0tile = phi0TextureOpts.width_phi0tile ;
    potential_sigmay = phi0TextureOpts.potential_sigmay ;
end

% Identify previous sphi pullback
imfn_sp_prev = sprintf( QS.fullFileBase.im_sp, QS.currentTime - 1 ) ;

% Load 3D data for coloring mesh pullback
if isempty(QS.currentData.IV)
    QS.xp.loadTime(tt);
    QS.xp.rescaleStackToUnitAspect();

    % Raw stack data
    IV = QS.xp.stack.image.apply();
    IV = imadjustn(IV{1});         
    QS.currentData.IV = IV ;
else
    % Raw stack data
    IV = QS.currentData.IV ;
end

% Unpack QS
nU = QS.nU ;
nV = QS.nV ;
sphiDir = QS.dir.spcutMesh ;
phi0fitBase = QS.fullFileBase.phi0fit ;
fileNameBase = QS.fileBase.name ;

% Texture patch options
Options.PSize = 5;
Options.EdgeColor = 'none';
% Texture image options
Options.imSize = ceil( 1000 .* [ 1 QS.a_fixed ] );
Options.yLim = [0 1];

% fit the shifts in the y direction
% todo: save uncorrected patchIms,
% Could try tiling twice, but since phase correlation is periodic (Fourier)
% it shouldn't make any difference at all.
dmyk = 0 ;
if nargin < 4
    phi0_fit = zeros(size(uspace_ds_umax)) ;
else
    phi0_fit = phi0_init ;
end

% phi0s = zeros(size(uspace)) ;  % no longer used
phi0_fit_kk = 1 ; % for first pass                
phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
ensureDir([sphiDir, '/phi0_correction/'])
while any(phi0_fit_kk > 0.002) && dmyk < 6
    disp(['Iteration ' num2str(dmyk)])
    plotfn = sprintf(phi0fitBase, QS.currentTime, dmyk);
    
    % Load previous pullback
    debugFnPattern = fullfile([sphiDir, '/phi0_correction/', ...
                fileNameBase, '_prephi0_' num2str(dmyk) '.tif']) ;
    patchImFn = sprintf( debugFnPattern, QS.currentTime - 1 )  ;
    [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevPullback(IV, ...
        uvgrid3d, cutMesh.umax, uspace_ds_umax, phiv_kk, ...
        ringpath_ss, imfn_sp_prev, lowerboundy, upperboundy, ...
        save_ims, plotfn, Options, ...
        step_phi0tile, width_phi0tile, potential_sigmay, 'integer', ...
        patchImFn) ;

    % Update the result: phi0_fit is what MATTERS
    dmyk = dmyk + 1;
    phi0_fit = phi0_fit + phi0_fit_kk ;
    % phi0s = phi0s + phi0s_kk ;  % no longer used
    
    % phiv_kk is this iteration's correction only
    phiv_kk = (vspace .* ones(nU, nV))' - phi0_fit .* ones(nU, nV) ;
end

