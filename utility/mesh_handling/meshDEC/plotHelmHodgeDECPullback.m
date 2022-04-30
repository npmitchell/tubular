function plotHelmHodgeDECPullback(im, cutMesh, vf, xy, vf2d, ...
    divs, rots, Options, opts2d)
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
% divs : struct with fields
%   divv 
%       scalar divergence field values (div) defined on vertices
%   divU
%       vector dilatational field flow (curl-free) defined on faces
% rots : struct with fields 
%   rotv 
%       scalar rotational field values (curl) defined on faces
%   rotU
%       vector rotational field flow (div-free) defined on faces
% Options : struct with fields
%   div2dfn : str
%       output filename for dilatational flow in 2d
%   div3dfn_lateral : str
%       output filename for dilatational flow in 3d
%   div3dfn_ventral : str
%       output filename for dilatational flow in 3d
%   rot2dfn : str
%       output filename for rotational flow in 2d
%   rot3dfn_lateral : str
%       output filename for rotational flow in 3d
%   rot3dfn_ventral : str
%       output filename for rotational flow in 3d
%   xyzlim : 3 x 2 float array
%       limits in each dimension for 3d plot
%   qsubU : int (optional)
%       subsampling factor for quiver overlay in pullback x coord
%   qsubV : int (optional)
%       subsampling factor for quiver overlay in pullback y coord
%   qscaleDiv2d : float (optional)
%       quiver arrow scale for 2d images for Divergence images
%   qscaleRot2d : float (optional)
%       quiver arrow scale for 2d images for Rotational images
%   qscaleDiv3d : float (optional)
%       quiver arrow scale for 3d images for Divergence images
%   qscaleRot3d : float (optional)
%       quiver arrow scale for 3d images for Rotational images
%   sscaleDiv : float (optional)
%       scalar field scale, so color limit is (-sscaleDiv, +sscaleDiv)
%   sscaleRot : float (optional)
%       scalar field scale, so color limit is (-sscaleRot, +sscaleRot)
%   addTitleStr : str (optional)
%       addition to title (for ex, timestamp specifier)
%   alpha : float (default=0.6)
%       opacity for heatmap to overlay image
%   
% 
% Returns
% -------
%
% NPMitchell 2020

% Unpack the cutMesh
FF = cutMesh.f ;
V2D = cutMesh.u ;
v3drs = cutMesh.v ;
nU = cutMesh.nU ;
nV = cutMesh.nV ;

% compute COM for each triangle in 2d and 3d --> note that we do this
% before gluing so that the 2d case is not messed up. The 3d case is
% the same either way. bc are the face barycenters
bc = cat( 3, v3drs(FF(:,1), :), v3drs(FF(:,2), :), v3drs(FF(:,3), :) );
bc = mean( bc, 3 ) ;
% bc2d = cat( 3, V2D(FF(:,1), :), V2D(FF(:,2), :), V2D(FF(:,3), :) );
% bc2d = mean( bc2d, 3 ) ;

% Unpack divs and rots
divv = divs.divv ;
% divU = divs.divU ;
% divU2d = divs.divU2d ;
rotv = rots.rotv ;
% rotU = rots.rotU ;
% rotU2d = rots.rotU2d ;

% Unpack required Options
div2dfn = Options.div2dfn ;
div3dfnL = Options.div3dfn_lateral ;
div3dfnV = Options.div3dfn_ventral ;
rot2dfn = Options.rot2dfn ;
rot3dfnL = Options.rot3dfn_lateral ;
rot3dfnV = Options.rot3dfn_ventral ;

% Unpack other options
qsubU = 2 ;
qsubV = 2 ;
sscaleDiv = 0.6 ;
sscaleRot = 0.15 ;
qscaleDiv2d = 0 ;
qscaleRot2d = 0 ;
qscaleDiv3d = 0 ;
qscaleRot3d = 0 ;
addTitleStr = '' ;
labelUnit = '[1/min]' ;
harmLabelUnit = '[$\mu$m/min]' ;
xyzlim = [] ;
doubleCovered = true ;
alphaVal = 0.6 ;
% Unpack options from struct
if isfield(Options, 'sscaleDiv')
    sscaleDiv = Options.sscaleDiv ;
end
if isfield(Options, 'sscaleRot')
    sscaleRot = Options.sscaleRot ;
end
% Basic default qscaleDiv/Rot set without 2d or 3d specification
if isfield(Options, 'qscaleDiv')
    qscaleDiv2d = Options.qscaleDiv ;
end
if isfield(Options, 'qscaleRot')
    qscaleRot2d = Options.qscaleRot ;
end
% Custom 2d qscaleDiv/Rot set with 2d specification
if isfield(Options, 'qscaleDiv2d')
    qscaleDiv2d = Options.qscaleDiv2d ;
end
if isfield(Options, 'qscaleRot2d')
    qscaleRot2d = Options.qscaleRot2d ;
end
% Custom 3d qscaleDiv/Rot set with 3d specification
if isfield(Options, 'qscaleDiv3d')
    qscaleDiv3d = Options.qscaleDiv3d ;
end
if isfield(Options, 'qscaleRot3d')
    qscaleRot3d = Options.qscaleRot3d ;
end
% Title and units
if isfield(Options, 'addTitleStr')
    addTitleStr = Options.addTitleStr ;
end
if isfield(Options, 'labelUnit')
    labelUnit = Options.labelUnit ;
end
if isfield(Options, 'harmLabelUnit')
    harmLabelUnit = Options.harmLabelUnit ;
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
sscales = [ sscaleDiv, sscaleRot ] ;
qscales2d = [ qscaleDiv2d, qscaleRot2d ] ;
qscales3d = [ qscaleDiv3d, qscaleRot3d ] ;
if doubleCovered
    ylims = [0.25 * size(im, 1), 0.75 * size(im, 1)] ;
else
    ylims = [0, size(im, 1)] ;
end


% Save divergence and curl images
titlestrs = {[ 'dilatational flow, $\nabla \cdot v_t$' addTitleStr], ...
    [ 'rotational flow, $\star \mathrm{d} v_t^\flat$' addTitleStr], ...
    [ 'harmonic component of $v_t$:  ' addTitleStr]} ;
labelstrs = {['$\star$d$\star\left(v_t^\flat\right)$, ' labelUnit], ...  % this is div(v)
    ['$\star$d$v_t^\flat$, ' labelUnit], ... 
    ... %Note that here the vector field curl v = \left(\star$d$v_t^\flat\right)^\sharp  
    ... % is NOT what is being plotted, since we are visualizing the 1form \star$d$v_t^\flat 
    ['harm$(v_t)$ ' harmLabelUnit]} ;
fnstrs2d = {div2dfn, rot2dfn} ;
fnstrs3dL = {div3dfnL, rot3dfnL} ;
fnstrs3dV = {div3dfnV, rot3dfnV} ;
sfs = {divv, rotv} ;
% vf3ds = {divU, rotU} ;
% vf2ds = {divU2d, rotU2d} ;

%% Consider both div, then rot
vf2ds = {vf2d, vf2d } ;
for divcurl = 1:2 
    %% Plot the div/rot as 3d surface 
    % 3D mesh plot with scalar and vector fields
    opts = struct() ;
    opts.style = 'diverging' ;
    opts.xlabel = 'AP position [\mum]' ;
    opts.ylabel = 'lateral position [\mum]' ;
    opts.zlabel = 'DV position [\mum]' ;
    if ~isempty(xyzlim)
        opts.xlim = xyzlim(1, :) ;
        opts.ylim = xyzlim(2, :) ;
        opts.zlim = xyzlim(3, :) ;
    end
    opts.view = [0, 0] ;
    opts.title = titlestrs{divcurl} ;
    opts.label = labelstrs{divcurl} ;
    opts.sscale = sscales(divcurl) ;
    opts.cticks = [-sscales(divcurl), 0, sscales(divcurl)] ;
    opts.qscale = qscales3d(divcurl) ;
    opts.linewidth = 0.8 ;
    inds = (0:qsubU:(nU-1))' * 2 * (nU-1) * ones(length(1:qsubV:(nV-1)), 1)' +...
         ones(length(1:qsubU:(nU-1)), 1) * 2 * (1:qsubV:(nV-1)) ;
    % vfq = vf3ds{divcurl} ;
    % Plot the fields in 3d
    close all
    fig = scalarVectorFieldsOnSurface(FF, v3drs, sfs{divcurl}(:), ...
                bc(inds,1), bc(inds,2), bc(inds,3), ...
                vf(inds, 1), vf(inds, 2), vf(inds, 3), opts) ;
    
    % Save multiple views
    % figure parameters
    figWidth = 16 ;  % cm
    figHeight = 10 ; % cm
    outfn3d = fnstrs3dL{divcurl} ; 
    disp(['scalarVectorFieldsOnImage: saving ' outfn3d])
    set(fig, 'PaperUnits', 'centimeters');
    set(fig, 'PaperPosition', [0 0 figWidth figHeight]);  
    saveas(fig, outfn3d) ;  
    
    % Change view
    outfn3d = fnstrs3dV{divcurl} ; 
    view(0, 270)
    saveas(fig, outfn3d) ;  
    
         
    close all    

    %% Plot the 2d div/rot as heatmap on pullback image
    sf = reshape(sfs{divcurl}, [nU, nV])' ;
    xxs = V2D(1:nU, 1) ;
    yys = V2D(1:nU:nU*nV, 2) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rot used to be on faces, so reshape each differently here
    % if divcurl == 1
    %     % DIVERGENCE COMPONENT
    %     sf = reshape(sfs{divcurl}, [nU, nV])' ;
    %     xxs = V2D(1:nU, 1) ;
    %     yys = V2D(1:nU:nU*nV, 2) ;
    % elseif divcurl == 2
    %     % ROTATIONAL COMPONENT        
    %     sf = sfs{divcurl} ;
    %     xxs = V2D(:, 1) ;
    %     yys = V2D(:, 2) ;
    % else
    %     error('Bad divcurl index. Should be 1 or 2')
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Handle cases of on-vertex and on-face separately
    set(gcf, 'visible', 'off')
    if qscales2d(divcurl) > 0 
        % opts2d.outfn = fnstrs2d{divcurl} ;
        % opts2d.qsubsample = 1 ;
        % opts2d.faces = FF ;
        % opts2d.title = titlestrs{divcurl} ;
        % opts2d.label = labelstrs{divcurl}  ;
        % opts2d.sscale = sscales(divcurl) ;
        % opts2d.qscale = qscales2d(divcurl) ;
        % opts2d.style = 'diverging' ;
        % vfq = vf2ds{divcurl} ;
        error('have not implemented quiver on image yet, use scalarVectorFieldOnImage()')
    else
        labelOptions = struct() ;
        labelOptions.label = labelstrs{divcurl} ;
        labelOptions.xlabel = 'ap position' ;
        labelOptions.ylabel = 'circumferential position' ;
        labelOptions.title = titlestrs{divcurl} ;
        if any(size(sf) == size(xxs, 1))
            scalarFieldOnImage(im, [xxs, yys], sf, alphaVal, ...
                sscales(divcurl), labelOptions) ;
        elseif size(FF, 1) == size(sf, 1)
            fxy.faces = FF ;
            fxy.vertices = V2D ;
            scalarFieldOnImage(im, fxy, sf, alphaVal, sscales(divcurl), ...
                labelOptions) ;
        else
            error('sf is wrong size. Handle here')
        end
        ylim(ylims) 
        saveas(gcf, fnstrs2d{divcurl})
        close all
    end
    
end

% % Plot the harmonic portion of the Helmholtz-Hodge decomposition
% % 3D harmonic field
% opts.style = 'phase' ;
% opts.xlabel = 'AP position [\mum]' ;
% opts.ylabel = 'lateral position [\mum]' ;
% opts.zlabel = 'DV position [\mum]' ;
% opts.xlim = xyzlim(1, :) ;
% opts.ylim = xyzlim(2, :) ;
% opts.zlim = xyzlim(3, :) ;
% opts.titlestr = titlestrs{3} ;
% opts.view = [0, 0] ;
% opts.label = labelstrs{3} ;
% opts.outfn = harm3dfn ;
% opts.sscale = 5 ;
% inds = (0:qsubU:(nU-1))' * 2 * (nU-1) * ones(length(1:qsubV:(nV-1)), 1)' +...
%      ones(length(1:qsubU:(nU-1)), 1) * 2 * (1:qsubV:(nV-1)) ;
% harmvphase = atan2(harmU2d(:, 2), harmU2d(:, 1)) ;
% scalarVectorFieldsOnSurface(FF, v3drs, harmvphase, ...
%             bc(inds,1), bc(inds,2), bc(inds,3), ...
%             harmU(inds, 1), harmU(inds, 2), harmU(inds, 3), opts) ;
% close all    
% 
% % 2D harmonic field
% options.outfn = harm2dfn ;
% vectorFieldHeatPhaseOnImage(im, V2D(:, 1), V2D(:, 2),...
%     harmU2d(:, 1), harmU2d(:, 2), 5, options)


% GPtoolbox for laplacian smoothing on VERTICES for div
%  laplacian_smooth(V,F,L_method,b,lambda,method,S,max_iter)