function plotHelmHodgeDECSurface3D(mesh, divs, rots, Options)
% PLOTHELMHODGEDECSURFACE3D(mesh, divs, rots, Options)
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
%   div3dfn : str
%       output filename for dilatational flow in 3d
%   rot2dfn : str
%       output filename for rotational flow in 2d
%   rot3dfn : str
%       output filename for rotational flow in 3d
%   xyzlim : 3 x 2 float array
%       limits in each dimension for 3d plot
%   qsubU : int (optional)
%       subsampling factor for quiver overlay in pullback x coord
%   qsubV : int (optional)
%       subsampling factor for quiver overlay in pullback y coord
%   sscaleDiv : float (optional)
%       scalar field scale, so color limit is (-sscaleDiv, +sscaleDiv)
%   sscaleRot : float (optional)
%       scalar field scale, so color limit is (-sscaleRot, +sscaleRot)
%   addTitleStr : str (optional)
%       addition to title (for ex, timestamp specifier)
%   
% 
% Returns
% -------
%
% NPMitchell 2020

% Unpack the cutMesh
FF = mesh.f ;
V2D = mesh.u ;
v3d = mesh.v ;
v3drs = mesh.v3drs ;
nU = mesh.nU ;
nV = mesh.nV ;

% compute COM for each triangle in 2d and 3d --> note that we do this
% % before gluing so that the 2d case is not messed up. The 3d case is
% % the same either way.
% bc = cat( 3, v3drs(FF(:,1), :), v3drs(FF(:,2), :), v3drs(FF(:,3), :) );
% bc = mean( bc, 3 ) ;
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
div3dfn = Options.div3dfn ;
rot3dfn = Options.rot3dfn ;

% Unpack other options
qsubU = 5 ;
qsubV = 10 ;
sscaleDiv = 0.5 ;
sscaleRot = 0.2 ;
qscaleDiv = 10 ;
qscaleRot = 25 ;
addTitleStr = '' ;
labelUnit = '[1/min]' ;
harmLabelUnit = '[$\mu$m/min]' ;
xyzlim = [] ;
% Unpack options from struct
if isfield(Options, 'sscaleDiv')
    sscaleDiv = Options.sscaleDiv ;
end
if isfield(Options, 'sscaleRot')
    sscaleRot = Options.sscaleRot ;
end
if isfield(Options, 'qscaleDiv')
    qscaleDiv = Options.qscaleDiv ;
end
if isfield(Options, 'qscaleRot')
    qscaleRot = Options.qscaleRot ;
end
if isfield(Options, 'addTitleStr')
    addTitleStr = Options.addTitleStr ;
end
if isfield(Options, 'labelUnit')
    labelUnit = Options.labelUnit ;
end
if isfield(Options, 'harmLabelUnit')
    harmLabelUnit = Options.harmLabelUnit ;
end
if isfield(Options, 'qsubU')
    qsubU = Options.qsubU ;
end
if isfield(Options, 'qsubV')
    qsubV = Options.qsubV ;
end
if isfield(Options, 'xyzlim')
    xyzlim = Options.xyzlim ;
end

% Pack scales for scalar magnitude and quiver length
sscales = [sscaleDiv, sscaleRot ] ;
qscales = [qscaleDiv, qscaleRot ] ;

% Save divergence and curl images
titlestrs = {[ 'dilatational flow, $\nabla \cdot v_t$:  ' addTitleStr], ...
    [ 'rotational flow, $\star \mathrm{d} v_t^\flat$:  ' addTitleStr], ...
    [ 'harmonic component of $v_t$:  ' addTitleStr]} ;
labelstrs = {['$\nabla \cdot v_t$, ' labelUnit], ...
    ['$\star$d$v_t^\flat$, ' labelUnit], ...
    ['harm$(v_t)$ ' harmLabelUnit]} ;
fnstrs3dt = {div3dfn, rot3dfn} ;
sfs = {divv, rotv} ;

for divcurl = 1:2 
    % 3D mesh plot with scalar and vector fields
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
    opts.outfn = fnstrs3dt{divcurl} ;
    opts.sscale = sscales(divcurl) ;
    
    scalarVectorFieldsOnSurface(faces, vertices, sfs{divcurl},...
        xxv, yyv, zzv, vx, vy, vz, opts)
    
end
