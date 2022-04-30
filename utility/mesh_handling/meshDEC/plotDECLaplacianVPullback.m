function plotDECLaplacianVPullback(im, cutMesh, vf, xy, lapvs, Options, opts2d)
% PLOTHELMHODGEDECPULLBACK(im, Options)
%   Plot the helmholtz-hodge decomposition of divergence and rotational
%   component of tangential velocities onto a pullback image
%
% Parameters
% ----------
% im : 2D image pullback
% cutMesh : struct with fields 
%   u : #vertices x 2 float array
%       2D vertices of mesh pullback
%   v : #vertices x 3 float array
%       3D vertices of mesh embedding
%   f : #faces x 3 int
%       face connectivity list
% vf : #faces x 3 float array 
%   original vector field before decomposition, in 3d
% vf2d : #faces x 3 float array 
%   original tangential vector field before decomposition, in 2d, properly
%   scaled by dilation (or not) for visualization
% lapvs : struct with field 
%   lapv 
%       scalar rotational field values (curl) defined on faces
% Options : struct with fields
%   lap2dfn : str
%       output filename for rotational flow in 2d
%   lap3dfn_lateral : str
%       output filename for rotational flow in 3d
%   lap3dfn_ventral : str
%       output filename for rotational flow in 3d
%   xyzlim : 3 x 2 float array
%       limits in each dimension for 3d plot
%   qsubU : int (optional)
%       subsampling factor for quiver overlay in pullback x coord
%   qsubV : int (optional)
%       subsampling factor for quiver overlay in pullback y coord
%   qscale2d : float (optional)
%       quiver arrow scale for 2d images for Lap images
%   qscale3d : float (optional)
%       quiver arrow scale for 3d images for Lap images
%   sscale : float (optional)
%       scalar field scale, so color limit is (0, +sscale)
%   addTitleStr : str (optional)
%       addition to title (for ex, timestamp specifier)
%   alpha : float (default=0.6)
%       opacity for heatmap to overlay image
%   
% 
% Returns
% -------
%
% Saves to disk
% -------------
%
% NPMitchell 2021

% Unpack the cutMesh
FF = cutMesh.f ;
V2D = cutMesh.u ;
v3drs = cutMesh.v ;
nU = cutMesh.nU ;
nV = cutMesh.nV ;

lapv3d = lapvs.lapv ;
lapv2d = lapvs.lapv2dsc ;
lapv2dang = atan2(lapv2d(:, 2), lapv2d(:, 1)) ;
lapv2dmag = vecnorm(lapvs.lapvt, 2, 2) ;
lapvn = lapvs.lapvn ;

% compute COM for each triangle in 2d and 3d --> note that we do this
% before gluing so that the 2d case is not messed up. The 3d case is
% the same either way. bc are the face barycenters
bc = cat( 3, v3drs(FF(:,1), :), v3drs(FF(:,2), :), v3drs(FF(:,3), :) );
bc = mean( bc, 3 ) ;
% bc2d = cat( 3, V2D(FF(:,1), :), V2D(FF(:,2), :), V2D(FF(:,3), :) );
% bc2d = mean( bc2d, 3 ) ;

% Unpack required Options
lapv2dfn = Options.lapv2dfn ;
lapv3dfnL = Options.lapv3dfn_lateral ;
lapv3dfnV = Options.lapv3dfn_ventral ;

% Unpack other options
qsubU = 2 ;
qsubV = 2 ;
sscale = 0.6 ;
qscale2d = 0 ;
qscale3d = 0 ;
addTitleStr = '' ;
labelUnit = '[$(\mu$m min$)^{-1}$]' ;
xyzlim = [] ;
doubleCovered = true ;
alphaVal = 0.6 ;
% Unpack options from struct
if isfield(Options, 'sscale')
    sscale = Options.sscale ;
end
% Basic default qscale set without 2d or 3d specification
if isfield(Options, 'qscale')
    qscale2d = Options.qscale ;
    qscale3d = Options.qscale ;
end
% Custom 2d qscale set with 2d specification
if isfield(Options, 'qscale2d')
    qscale2d = Options.qscale2d ;
end
% Custom 3d qscale set with 3d specification
if isfield(Options, 'qscale3d')
    qscale3d = Options.qscale3d ;
end
% Title and units
if isfield(Options, 'addTitleStr')
    addTitleStr = Options.addTitleStr ;
end
if isfield(Options, 'labelUnit')
    labelUnit = Options.labelUnit ;
end
% Quiver field subsampling
if isfield(Options, 'qsubU')
    qsubU = Options.qsubU ;
end
if isfield(Options, 'qsubV')
    qsubV = Options.qsubV ;
end
% Axis limits
if isfield(Options, 'xyzlim')
    xyzlim = Options.xyzlim ;
end
if isfield(Options, 'doubleCovered')
    doubleCovered = Options.doubleCovered;
end
if isfield(Options, 'alpha')
    alphaVal = Options.alpha ;
end

% Pack scales for scalar magnitude and quiver length
qscales2d = qscale2d ;
qscales3d = qscale3d ;
if doubleCovered
    ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
else
    ylims = [0, size(im, 1)] ;
end


% Save divergence and curl images
titlestr = {[ 'Laplacian of flow, $\nabla^2 v_t$' addTitleStr]} ;
% dd1 * hd1 * d0 * v :  v-> primal 1form -> dual 1form
% DEC Laplace = d*d(v)
labelstr = {['$\big($d$\star$d$v_\parallel^\flat\big)_\perp$, ' labelUnit], '', ...
    ['$\big($d$\star$d$v_\parallel^\flat\big)_\parallel$, ' labelUnit], ''} ;
fnstr2d = lapv2dfn ;
fnstr3dL = lapv3dfnL ;

%% Plot the laplacian(v) as 3d surface 
% 3D mesh plot with scalar and vector fields
opts = struct() ;
opts.style = 'diverging' ;
opts.xlabel = 'AP position [\mum]' ;
opts.ylabel = 'lateral position [\mum]' ;
opts.zlabel = 'DV position [\mum]' ;

opts.xyzlims = {xyzlim, [0, 1;0,1;-eps,eps], xyzlim, [0, 1;0,1;-eps,eps]} ;
opts.title = titlestr ;
opts.labels = labelstr ;
opts.clim = sscale ;
opts.cticks = [0, sscale] ;
opts.qscale = qscale3d ;
opts.linewidth = 0.8 ;
opts.axisOff = true ;
opts.polarStyle = 'polar'; 
opts.masterCbar = true ;
opts.view = {[0, 0], [0,90], [0,0],[0,90]};
inds = (0:qsubU:(nU-1))' * 2 * (nU-1) * ones(length(1:qsubV:(nV-1)), 1)' +...
     ones(length(1:qsubU:(nU-1)), 1) * 2 * (1:qsubV:(nV-1)) ;
% Plot the fields in 3d
close all
m2d = cutMesh ;
m2d.v = cutMesh.u ;
m2d.v(:, 1) = m2d.v(:, 1) ./ max(m2d.v(:, 1)) ;
m2d.v(:, 2) = m2d.v(:, 2) ./ max(m2d.v(:, 2)) ;
m2d = rmfield(m2d, 'u') ;

[axs, cbs, meshHandles] = nFieldsOnSurface({cutMesh, m2d, cutMesh, m2d}, ...
    {lapvn, lapvn, ...
     {lapv2dmag, lapv2dang}, {lapv2dmag, lapv2dang}}, opts) ;

% Save multiple views
% figure parameters
figWidth = 16 ;  % cm
figHeight = 10 ; % cm
outfn3d = fnstr3dL ; 
disp(['scalarVectorFieldsOnImage: saving ' outfn3d])
fig = gcf ;
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperPosition', [0 0 figWidth figHeight]);  
saveas(fig, outfn3d) ;  

% Change view
outfn3d = lapv3dfnV ; 
opts.view = {[0, 0], [0,270], [0,0],[0,270]};
opts.clim = sscale * 2 ;
[axs, cbs, meshHandles] = nFieldsOnSurface({cutMesh, m2d, cutMesh, m2d}, ...
    {lapvn, lapvn, ...
     {lapv2dmag, lapv2dang}, {lapv2dmag, lapv2dang}}, opts) ;
saveas(fig, outfn3d) ;  

close all    
% 
% %% Plot the 2d laplacian as quiver/polar field on pullback image
% sf = lapv2dmag ;
% xxs = V2D(1:nU, 1) ;
% yys = V2D(1:nU:nU*nV, 2) ;
% 
% % Handle cases of on-vertex and on-face separately
% set(gcf, 'visible', 'off')
% if qscales2d > 0 
%     % opts2d.outfn = fnstrs2d{divcurl} ;
%     % opts2d.qsubsample = 1 ;
%     % opts2d.faces = FF ;
%     % opts2d.title = titlestr ;
%     % opts2d.label = labelstr  ;
%     % opts2d.sscale = sscale ;
%     % opts2d.qscale = qscales2d ;
%     % opts2d.style = 'diverging' ;
%     % vfq = vf2ds{divcurl} ;
%     error('have not implemented quiver on image yet, use scalarVectorFieldOnImage()')
% else
%     labelOptions = struct() ;
%     labelOptions.label = labelstr ;
%     labelOptions.xlabel = 'ap position' ;
%     labelOptions.ylabel = 'circumferential position' ;
%     labelOptions.title = titlestr ;
%     if any(size(sf) == size(xxs, 1))
%         vectorFieldHeatPhaseOnImage(im, [xxs, yys], lapv2d(:, 1), lapv2d(:,2), ...
%             sscales, labelOptions) ;
%     elseif size(FF, 1) == size(sf, 1)
%         fxy = struct() ;
%         fxy.faces = FF ;
%         fxy.vertices = V2D ;
%         vectorFieldHeatPhaseOnImage(im, fxy, lapv2d(:, 1), lapv2d(:,2), sscale, ...
%             labelOptions) ;
%     else
%         error('sf is wrong size. Handle here')
%     end
%     ylim(ylims) 
%     saveas(gcf, fnstr2d)
%     close all
% end
    
end
