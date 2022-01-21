function generateCurrentPullbacks(QS, cutMesh, spcutMesh, ...
    spcutMeshSm, pbOptions)
%GENERATECURRENTPULLBACKS(QS, cutMesh, spcutMesh, spcutMeshSm, pbOptions)
%   Generate 2D images of tissue mapped from 3D
%
% Parameters
% ----------
% pbOptions : struct
%   overwrite : bool (default=false)
%       overwrite existing images on disk
%   generate_sphi : bool (default=false)
%       create pullbacks in uphi coords
%   generate_relaxed : bool (default=false)
%       create pullbacks in uphi coords
%   generate_uv : bool (default=false)
%       create pullbacks in uphi coords
%   generate_uphi : bool (default=false)
%       create pullbacks in uphi coords
%   generate_spsm : bool (default=false)
%       create pullbacks in uphi coords
%   generate_rsm : bool (default=false)
%       create pullbacks in uphi coords
%   generate_uvprime : bool (default=false)
%       create pullbacks in as-conformal-as-possible with minimally twisted
%       boundaries, in smoothed meshes
%   PSize : int (default=5)
%       how many interpolation points along each side of a patch
%   EdgeColor : colorspec (default='none')
%       color of mesh edges in pullback image
%   YLim : 1x2 float (default=[0 1])
%       Y extent of pullback image in v or phi coords
%   axisorder : length 3 ints (default = [1 2 3 ])
%       axis permutation for the texture mesh
%   preTextureLambda : float (default = 0)
%       If nonzero, apply laplacian smooth to mesh before rendering and
%       before moving along normal_shift (which occurs before texture
%       mapping)
%
%   Additional options as fields passed to texturePatch are
%       - pbOptions.imSize:       The size of the output image
%       - pbOptions.baseSize:     The side length in pixels of the smallest
%                               side of the output image when using the
%                               tight mesh bounding box
%       - pbOptions.xLim:         The x-bounds of the output image
%       - pbOptions.yLim:         The y-bounds of the output image
%       - pbOptions.pixelSearch:  The method of searching for the faces
%                               containing pixel centers
%                                   - 'AABB' (requires GPToolBox)
%                                   - 'Default' (MATLAB built-ins, faster than AABB)
%       - pbOptions.numLayers:    The number of onion layers to create
%                               Format is [ (num +), (num -) ]
%       - pbOptions.layerSpacing: The spacing between adjacent onion layers
%                               in units of pixels
%       - pbOptions.smoothIter:   Number of iterations of Laplacian mesh
%                               smoothing to run on the mesh prior to
%                               vertex normal displacement (requires
%                               GPToolBox) (Default is 0)
%       - pbOptions.vertexNormal: User supplied vertex unit normals to the
%                               texture triangulation
%       - pbOptions.Interpolant:  A pre-made texture image volume interpolant
%
% NPMitchell 2020

%% Unpack pbOptions
overwrite = false ;         % generate pullbacks even if they exist on disk
generate_sphi    = true  ;  % generate an (s, phi) coord system pullback
generate_relaxed = false ;  % generate a relaxed (s,phi) coord system pullback
generate_uv      = false ;  % generate a (u,v) coord system pullback
generate_uphi    = false ;  % generate a (u, phi) coord system pullback
generate_spsm    = false ;  % generate an (s, phi) coord system smoothed mesh pullback
generate_rsm     = false ;  % generate a relaxed (s, phi) coord system smoothed mesh pullback
generate_ricci   = false ;  % generate a ricci-flowed coord system pullback
generate_uvprime = false ;  % generate a (u',v') coord sys pullback, where u',v' are found by conformal map with minimally twisted boundaries
generate_ruvprime = false ;  % generate a (u',v') coord sys pullback, where u',v' are found by relaxed conformal map with minimally twisted boundaries
% Other options
save_as_stack    = false ;  % save data as stack for each timepoint, not MIP
channels = [] ;             % default is to image all channels (empty list)
axisorder = [1 2 3 ];   % Note: we should NOT use QS.data.axisOrder here, 
                        % since that is invoked upon loading IV instead of applying in post
normal_shift = 0 ;          % how much to shift mesh along vertex normals before rendering
preTextureLambda = 0 ;      % how much to smooth mesh before rendering

% Replace defaults
if nargin > 4
    disp('Unpacking options for which pullbacks to generate/overwrite')
    if isfield(pbOptions, 'overwrite')
        overwrite = pbOptions.overwrite ;
        pbOptions = rmfield(pbOptions, 'overwrite') ;
    end
    if isfield(pbOptions, 'save_as_stack')
        save_as_stack = pbOptions.save_as_stack ;
        pbOptions = rmfield(pbOptions, 'save_as_stack') ;
    end
    if isfield(pbOptions, 'generate_sphi')
        generate_sphi = pbOptions.generate_sphi ;
        pbOptions = rmfield(pbOptions, 'generate_sphi') ;
    end
    if isfield(pbOptions, 'generate_relaxed')
        generate_relaxed = pbOptions.generate_relaxed ;
        pbOptions = rmfield(pbOptions, 'generate_relaxed') ;
    end
    if isfield(pbOptions, 'generate_uv')
        generate_uv = pbOptions.generate_uv ;
        pbOptions = rmfield(pbOptions, 'generate_uv') ;
    end
    if isfield(pbOptions, 'generate_uphi_coord')
        generate_uphi = pbOptions.generate_uphi_coord ;
        pbOptions = rmfield(pbOptions, 'generate_uphi_coord') ;
    end
    if isfield(pbOptions, 'generate_spsm')
        generate_spsm = pbOptions.generate_spsm ;
        pbOptions = rmfield(pbOptions, 'generate_spsm') ;
    end
    if isfield(pbOptions, 'generate_rsm')
        generate_rsm = pbOptions.generate_rsm ;
        pbOptions = rmfield(pbOptions, 'generate_rsm') ;
    end
    if isfield(pbOptions, 'generate_ricci')
        generate_ricci = pbOptions.generate_ricci ;
        pbOptions = rmfield(pbOptions, 'generate_ricci') ;
    end
    if isfield(pbOptions, 'generate_uvprime')
        generate_uvprime = pbOptions.generate_uvprime ;
        pbOptions = rmfield(pbOptions, 'generate_uvprime') ;
    end
    if isfield(pbOptions, 'generate_ruvprime')
        generate_ruvprime = pbOptions.generate_ruvprime ;
        pbOptions = rmfield(pbOptions, 'generate_ruvprime') ;
    end
    if isfield(pbOptions, 'channels')
        channels = pbOptions.channels ;
        pbOptions = rmfield(pbOptions, 'channels') ;
    end
    if isfield(pbOptions, 'axisorder')
        axisorder = pbOptions.axisorder ;
        pbOptions = rmfield(pbOptions, 'axisorder') ;
    end
    if isfield(pbOptions, 'preTextureLambda')
        preTextureLambda = pbOptions.preTextureLambda ;
        pbOptions = rmfield(pbOptions, 'preTextureLambda') ;
    end
    if isfield(pbOptions, 'normal_shift')
        normal_shift = pbOptions.normal_shift ;
        pbOptions = rmfield(pbOptions, 'normal_shift') ;
    end
end

%% Unpack options
if nargin < 2 || isempty(cutMesh)
    if isempty(QS.currentMesh.cutMesh)
        QS.loadCurrentCutMesh()
    end
    cutMesh = QS.currentMesh.cutMesh ;
end

if (nargin < 3 || isempty(spcutMesh)) && (generate_uv || generate_relaxed || generate_sphi)
    if isempty(QS.currentMesh.spcutMesh)
        QS.loadCurrentSPCutMesh()
    end
    spcutMesh = QS.currentMesh.spcutMesh ;
    
    % Smooth embedding coordinates (optional)
    if preTextureLambda > 0
        tmpCutMesh = struct('f', spcutMesh.f, ...
            'v', spcutMesh.v, 'u', spcutMesh.sphi, ...
            'nU', spcutMesh.nU, 'nV', spcutMesh.nV, 'vn', spcutMesh.vn) ;
        glueMesh = glueRectCylinderCutMeshSeam(tmpCutMesh) ;
        glueMesh.v = laplacian_smooth(glueMesh.v, glueMesh.f, 'cotan', [], preTextureLambda, 'implicit') ;
        spcutMesh2 = cutRectilinearCylMesh(glueMesh) ;
        spcutMesh2.uv = spcutMesh.uv ;
        spcutMesh2.sphi = spcutMesh.sphi ;
        spcutMesh2.uphi = spcutMesh.uphi ;
        spcutMesh2.ar = spcutMesh.ar ;
        spcutMesh2.pathPairs = spcutMesh.pathPairs ;
        spcutMesh2.vn = per_vertex_normals(spcutMesh2.v, spcutMesh2.f, 'Weighting', 'angle') ;
        spcutMesh = spcutMesh2 ;
    end
    % Shift embedding coordinates along normal (optional)
    if abs(normal_shift) > 0
        spcutMesh.v = spcutMesh.v + normal_shift * spcutMesh.vn ;
    end
end

if (nargin < 4 || isempty(spcutMeshSm)) && (generate_rsm || generate_spsm)
    if isempty(QS.currentMesh.spcutMeshSm)
        QS.loadCurrentSPCutMeshSm()
    end
    spcutMeshSm = QS.currentMesh.spcutMeshSm ;
    
    % Smooth embedding coordinates (optional)
    if preTextureLambda > 0
        tmpCutMesh = struct('f', spcutMeshSm.f, ...
            'v', spcutMeshSm.v, 'u', spcutMeshSm.u, ...
            'nU', spcutMeshSm.nU, 'nV', spcutMeshSm.nV, 'vn', spcutMeshSm.vn) ;
        glueMesh = glueRectCylinderCutMeshSeam(tmpCutMesh) ;
        glueMesh.v = laplacian_smooth(glueMesh.v, glueMesh.f, 'cotan', [], preTextureLambda, 'implicit') ;
        spcutMeshSm2 = cutRectilinearCylMesh(glueMesh) ;
        spcutMeshSm2.ar = spcutMeshSm.ar ;
        spcutMeshSm2.pathPairs = spcutMeshSm.pathPairs ;
        spcutMeshSm2.vn = per_vertex_normals(spcutMeshSm2.v, spcutMeshSm2.f, 'Weighting', 'angle') ;
        spcutMeshSm = spcutMeshSm2 ;
    end
    % Shift embedding coordinates along normal (optional)
    if abs(normal_shift) > 0
        spcutMeshSm.v = spcutMeshSm.v + normal_shift * spcutMeshSm.vn ;
    end
end

if generate_uvprime || generate_ruvprime
    if isfield(pbOptions, 'uvpcutMesh') 
        uvpcutMesh = pbOptions.uvpcutMesh ;
    elseif isfield(pbOptions, 'uvpcutMeshSm') 
        uvpcutMesh = pbOptions.uvpcutMeshSm ;
    elseif isfield(pbOptions, 'uvprimecutMeshSm')
        uvpcutMesh = pbOptions.uvprimecutMeshSm ;
    else
        if isempty(QS.currentMesh.uvpcutMesh)
            QS.loadCurrentUVPrimeCutMesh()
        end
        uvpcutMesh = QS.currentMesh.uvpcutMesh ;
    end
    
    % Smooth embedding coordinates (optional)
    if preTextureLambda > 0
        uvpcutMesh.v = laplacian_smooth(uvpcutMesh.v, uvpcutMesh.f, 'cotan', [], preTextureLambda, 'implicit') ;
        uvpcutMesh.vn = per_vertex_normals(uvpcutMesh.v, uvpcutMesh.f, 'Weighting', 'angle') ;
    end
    % Shift embedding coordinates along normal (optional)
    if abs(normal_shift) > 0
        uvpcutMesh.v = uvpcutMesh.v + normal_shift * uvpcutMesh.vn ;
    end
end

if generate_ricci
    if isfield(pbOptions, 'ricciMesh')
        ricciMesh = pbOptions.ricciMesh ;
    end
end

%% Unpack QS
tt = QS.currentTime ;
a_fixed = QS.a_fixed ;
% fileNameBase = QS.fileBase.name ;
% imFolder = QS.dir.im ;
% imFolder_r = QS.dir.im_r ;
% imFolder_sp = QS.dir.im_sp ;
% imFolder_up = QS.dir.im_up ;
% imFolder_spsm = QS.dir.im_sp_sm ;
% imFolder_rsm = QS.dir.im_r_sm ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Checking whether to create pullback \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------
% Generate Output Image Files
%--------------------------------------------------------------
imfn_uv = sprintf( QS.fullFileBase.im_uv, tt); 
imfn_r = sprintf( QS.fullFileBase.im_r, tt) ;
imfn_sp = sprintf( QS.fullFileBase.im_sp, tt) ;
imfn_up = sprintf( QS.fullFileBase.im_up, tt) ;
if QS.dynamic
    imfn_spsm = sprintf( QS.fullFileBase.im_sp_sm, tt) ;
    imfn_rsm = sprintf( QS.fullFileBase.im_r_sm, tt) ;
end
imfn_uvprime = sprintf( QS.fullFileBase.im_uvprime, tt) ;
imfn_ruvprime = sprintf( QS.fullFileBase.im_r_uvprime, tt) ;
imfn_ricci = sprintf( QS.fullFileBase.im_ricci, tt) ;
do_pb1 = ~exist(imfn_uv, 'file') && generate_uv ;
do_pb2 = ~exist(imfn_r, 'file') && generate_relaxed ;
do_pb3 = ~exist(imfn_sp, 'file') && generate_sphi ;
do_pb4 = ~exist(imfn_up, 'file') && generate_uphi ;
if QS.dynamic
    do_pb5 = ~exist(imfn_spsm, 'file') && generate_spsm ;
    do_pb6 = ~exist(imfn_rsm, 'file') && generate_rsm ;
else
    do_pb5 = false ;
    do_pb6 = false ;
end
do_pb7 = ~exist(imfn_uvprime, 'file') && generate_uvprime ;
do_pb8 = ~exist(imfn_ruvprime, 'file') && generate_ruvprime ;
do_pb9 = ~exist(imfn_ricci, 'file') && generate_ricci;

do_pb = [do_pb1, do_pb2, do_pb3, do_pb4, do_pb5, do_pb6, do_pb7, do_pb8, do_pb9] ;
do_pullbacks = (any(do_pb) || overwrite) ;

if do_pullbacks
    % Declare what needs to be redone
    if overwrite
        disp('All pullback images will be recomputed & saved')
    else
        if do_pb1
            disp(['(u,v) PB will be generated: ', imfn_uv])
        end 
        if do_pb2
            disp(['Relaxed (s,phi) PB will be generated: ', imfn_r])
        end
        if do_pb3
            disp(['(s,phi) PB will be generated: ', imfn_sp])
        end
        if do_pb4
            disp(['(u,phi) PB will be generated: ', imfn_up])
        end
        if do_pb5
            disp(['Smooth (s,phi) PB will be generated: ', imfn_spsm])
        end
        if do_pb6
            disp(['Smooth relaxed (s,phi) PB will be generated: ', imfn_rsm])
        end
        if do_pb7
            disp(['Smooth uvprime PB will be generated: ', imfn_uvprime])
        end
        if do_pb8
            disp(['Smooth relaxed uvprime PB will be generated: ', imfn_ruvprime])
        end
        if do_pb8
            disp(['Ricci PB will be generated: ', imfn_ricci])
        end
    end     
    
    % Load 3D data for coloring mesh pullback
    QS.getCurrentData()
    % grab raw stack data
    IV = QS.currentData.IV ;
    
    % select channels
    if ~isempty(channels)
        IV = IV(channels) ;
    end
end

if (~exist(imfn_sp, 'file') || overwrite) && generate_sphi
    fprintf(['Generating SP output image: ' imfn_sp]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_sp,...
        pbOptions, axisorder, save_as_stack)
    spcutMesh = rmfield(spcutMesh, 'u') ;    
else
    disp('Skipping SP pullback image generation ')
end

if (~exist(imfn_up, 'file') || overwrite) && generate_uphi
    fprintf(['Generating uphi output image: ' imfn_up]);
    % Assigning field spcutMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = spcutMesh.uphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_up, ...
        pbOptions, axisorder, save_as_stack)
    spcutMesh = rmfield(spcutMesh, 'u') ;
else
    disp('Skipping UP pullback image generation ')
end

if (~exist(imfn_ricci, 'file') || overwrite) && generate_ricci
    fprintf(['Loading mesh for generating ricci output image: ' imfn_ricci]);
    % Assigning field ricciMesh.u to be [s, phi] (ringpath
    % and azimuthal angle)
    spcutMesh.u = ricciMesh.uphi ;
    aux_generate_orbifold( spcutMesh, a_fixed, IV, imfn_up, ...
        pbOptions, axisorder, save_as_stack)
    spcutMesh = rmfield(spcutMesh, 'u') ;
else
    disp('Skipping UP pullback image generation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Output Image File -- regular UV coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist(imfn_uv, 'file') || overwrite) && generate_uv
    % Generate output image in uv
    fprintf(['Generating UV output image: ' imfn_uv]);
    aux_generate_orbifold(cutMesh, a_fixed, IV, imfn_uv, ...
        pbOptions, axisorder, save_as_stack)
else
    disp('Skipping UV pullback image generation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save relaxed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if (~exist(imfn_r, 'file') || overwrite) && generate_relaxed
    disp('Generating relaxed image for sphi coords...')
    spcutMesh.u = spcutMesh.sphi ;
    aux_generate_orbifold(spcutMesh, spcutMesh.ar, IV, imfn_r, ...
        pbOptions, axisorder, save_as_stack)
    spcutMesh = rmfield(spcutMesh, 'u') ;
else
    disp('Skipping relaxed SP pullback image generation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save smoothed sp image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if QS.dynamic
    if generate_spsm && (~exist(imfn_spsm, 'file') || overwrite) 
        disp('Generating image for smoothed sphi coords...')
        aux_generate_orbifold(spcutMeshSm, a_fixed, IV, imfn_spsm, ...
            pbOptions, axisorder, save_as_stack)
    else
        disp('Skipping SPSm pullback image generation ')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save smoothed relaxed image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if QS.dynamic
    if (~exist(imfn_rsm, 'file') || overwrite) && generate_rsm
        disp('Generating relaxed image for sphi coords...')
        if ~isfield(spcutMeshSm, 'ar')
            % Compute relaxed aspect ratio
            tmp = spcutMeshSm.u ;
            tmp(:, 1) = tmp(:, 1) / max(tmp(:, 1)) ;
            arspsm = minimizeIsoarealAffineEnergy( spcutMeshSm.f, spcutMeshSm.v, tmp );
            spcutMeshSm.ar = arspsm ;
        end
        aux_generate_orbifold(spcutMeshSm, spcutMeshSm.ar, IV, imfn_rsm, ...
            pbOptions, axisorder, save_as_stack)
    else
        disp('Skipping relaxed SPSm pullback image generation ')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save smoothed u'v' image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if (~exist(imfn_uvprime, 'file') || overwrite) && generate_uvprime
    disp(['Generating image for smoothed uvprime coords: ' imfn_uvprime])
    aux_generate_orbifold(uvpcutMesh.raw, a_fixed, IV, imfn_uvprime, ...
        pbOptions, axisorder, save_as_stack)
else
    disp('Skipping UVPrimeSm pullback image generation ')
end

if (~exist(imfn_ruvprime, 'file') || overwrite) && generate_ruvprime
    disp(['Generating image for relaxed, smoothed uvprime coords: ' imfn_ruvprime])
    aux_generate_orbifold(uvpcutMesh.raw, uvpcutMesh.raw.ar_mu, IV, imfn_ruvprime, ...
        pbOptions, axisorder, save_as_stack)
else
    disp('Skipping relaxed UVPrimeSm pullback image generation ')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save ricci pullback image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if (~exist(imfn_ricci, 'file') || overwrite) && generate_ricci
    disp(['Generating image for ricci coords: ' imfn_ricci])
    aux_generate_orbifold(ricciMesh.rectangle, IV, imfn_ricci, ...
        pbOptions, axisorder, save_as_stack)
else
    disp('Skipping Ricci pullback image generation ')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save submesh array. Each cell element contains all the 
% submeshes for that TP, which in this case is just one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% meshStack{tidx} = cutMesh ;
% if generate_sphi_coord
%     spmeshStack{tidx} = spcutMesh ;
% end
fprintf('Done\n');