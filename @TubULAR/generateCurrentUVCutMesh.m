function uvcutMesh = generateCurrentUVCutMesh(tubi, options)
% generateCurrentUVCutMesh(tubi, cutMesh, spcutMeshOptions)
% Create a rectilinear (but periodic in 2nd dim) parameterization of the
% cutMesh with a nearly-conformal pullback projection uv. 
%     uvcutMesh.f 
%     uvcutMesh.uv
%     uvcutMesh.v 
%     uvcutMesh.vn  
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
% options : struct with fields, optional
%   overwrite : bool
%       overwrite previous results
%       
%
% Returns
% -------
% uvcutMesh : struct with fields
%   uv: 
%
% NPMitchell 2020

%% Default options
overwrite = false ;

%% Unpack options
cutMesh = tubi.getCurrentCutMesh() ;
if nargin > 1
    disp('Unpacking options')
    if isfield(options, 'overwrite')
        overwrite = options.overwrite ;
    end
end

%% Unpack tubi
tt = tubi.currentTime ;
nU = tubi.nU ;
nV = tubi.nV ;
phi_method = tubi.phiMethod ;
uvcutMeshfn = sprintf(tubi.fullFileBase.uvcutMesh, tt) ;
[rot, trans] = getRotTrans(tubi) ;
resolution = tubi.APDV.resolution ;
preview = tubi.plotting.preview ;
[~, ~, xyzlim_um] = tubi.getXYZLims() ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate s,phi coord system for rotated, scaled mesh (rs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Establishing s,phi coord system\n');
if ~exist(uvcutMeshfn, 'file') || overwrite
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
    [ TF, TV2D, TV3Drs, TVN3Drs ] = tileAnnularCutMesh( cutMeshrs, tileCount );

    %----------------------------------------------------------------------
    % Generate surface curves of constant s
    %----------------------------------------------------------------------
    % For lines of constant phi
    disp('Creating crude uv curves with du=const to define uspace by ds(u)')
    % Make grid
    uspace = linspace( 0, cutMesh.umax, nU )' ;
    vspace = linspace( 0, 1, nV )' ;

    disp('Casting grid of uv points into 3D...')
    
    % NOTE: first dimension indexes u, second indexes v
    for kk = 1:nU
        if mod(kk, 20) == 0
            disp(['u = ' num2str(kk / nU)])
        end
        uv_tmp = [uspace(kk) * ones(size(vspace)), vspace] ;
        uvgrid2d(kk, :, :) = uv_tmp ;
    end 
    uu = uvgrid2d(:, :, 1) ;
    vv = uvgrid2d(:, :, 2) ;
    uv2d = [uu(:), vv(:)] ;
    assert(all(size(uv2d)==[nU*nV, 2]))

    Xirs = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TV3Drs(:,1), 'linear', 'nearest') ;
    Yirs = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TV3Drs(:,2), 'linear', 'nearest') ;
    Zirs = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TV3Drs(:,3), 'linear', 'nearest') ;
    uv3dXrs = Xirs(uu(:), vv(:)) ;
    uv3dYrs = Yirs(uu(:), vv(:)) ;
    uv3dZrs = Zirs(uu(:), vv(:)) ;
    uv3drs = [uv3dXrs, uv3dYrs, uv3dZrs] ;
    % uvgrid3drs = reshape(uv3drs, [nU, nV, 3]) ;
    nXirs = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3Drs(:,1), 'linear', 'nearest') ;
    nYirs = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3Drs(:,2), 'linear', 'nearest') ;
    nZirs = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3Drs(:,3), 'linear', 'nearest') ;
    uv3dXnrs = nXirs(uu(:), vv(:)) ;
    uv3dYnrs = nYirs(uu(:), vv(:)) ;
    uv3dZnrs = nZirs(uu(:), vv(:)) ;
    uv3dnormals_rs = [uv3dXnrs, uv3dYnrs, uv3dZnrs] ;
    uv3dnormals_rs = uv3dnormals_rs ./ vecnorm(uv3dnormals_rs, 2, 2) ;
    assert(all(abs(vecnorm(uv3dnormals_rs, 2, 2) - 1) < 1e-14))
    % uv3dgridnormals_rs = reshape(uv3dnormals_rs, [nU, nV, 3]) ;
    
    Xi = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TV3D(:,1), 'linear', 'nearest') ;
    Yi = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TV3D(:,2), 'linear', 'nearest') ;
    Zi = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TV3D(:,3), 'linear', 'nearest') ;
    uv3dX = Xi(uu(:), vv(:)) ;
    uv3dY = Yi(uu(:), vv(:)) ;
    uv3dZ = Zi(uu(:), vv(:)) ;
    uv3d = [uv3dX, uv3dY, uv3dZ] ;
    % uvgrid3d = reshape(uv3d, [nU, nV, 3]) ;
    nXi = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3D(:,1), 'linear', 'nearest') ;
    nYi = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3D(:,2), 'linear', 'nearest') ;
    nZi = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2), TVN3D(:,3), 'linear', 'nearest') ;
    uv3dXn = nXi(uu(:), vv(:)) ;
    uv3dYn = nYi(uu(:), vv(:)) ;
    uv3dZn = nZi(uu(:), vv(:)) ;
    uv3dnormals = [uv3dXn, uv3dYn, uv3dZn] ;
    uv3dnormals = uv3dnormals ./ vecnorm(uv3dnormals, 2, 2) ;
    assert(all(abs(vecnorm(uv3dnormals, 2, 2) - 1) < 1e-14))
    % uv3dgridnormals = reshape(uv3dnormals, [nU, nV, 3]) ;
    
    % Equivalent procedure done in generateCurrentSPCutMesh:
    % I decided against this because it requires squeezing v values ever so
    % slightly into the domain (0,1) instead of leaving them [0, 1] with 0
    % and 1 inclusive.
    % [~, ~, uvgrid3d, uvgrid2d] = ...
    %    ringpathsGridSampling(uspace0, vspace, TF, TV2D, TV3Drs) ;

    
    % Now ensure that the pullback lies in ([0,1],[0,1]) unit square
    uv2d(:, 1) = uv2d(:, 1) ./ cutMesh.umax ;
    
    uvcutMesh.f = defineFacesRectilinearGrid(uvgrid2d, nU, nV) ;
    uvcutMesh.uv = uv2d ;
    uvcutMesh.v = uv3d ;
    uvcutMesh.vrs = uv3drs ;
    uvcutMesh.vn = uv3dnormals ;
    uvcutMesh.vnrs = uv3dnormals_rs ;
    uvcutMesh.pathPairs = [1:nU; (nV-1)*nU + 1:nV*nU]' ;
    uvcutMesh.nU = nU ;
    uvcutMesh.nV = nV ;
    uvcutMesh.readme = struct('readme', 'A nearly conformal mesh found by Dirichlet energy minimization. The principal fields are f, v, vrs, uv, as well as pathPairs, nU, and nV', ...
        'f', '#faces x 3 int, face connectivity list indexing in to vertices v0 (for uv mesh) or v (for sphi mesh)', ...
        'pathPairs', 'nUx2 int, indices of vertices which are identical across the periodic direciton', ...
        'nU', 'int, number of points along the longitudinal axis in the resampled grid uv or sphi', ...
        'nV', 'int, number of points along the circumferential axis in the resampled grid uv or sphi', ...
        'uv', 'nU x nV x 2, vertices at uv mesh embedding coordinates',...
        'v', 'nU x nV x 3, vertices at uv mesh embedding coordinates',...
        'vn', 'nU x nV x 3, vertex normals at sphi mesh embedding vertex coordinates in pixel space', ...
        'vnrs', 'nU x nV x 3, vertex normals at sphi mesh embedding vertex coordinates in APDV frame') ;
    disp(['Saving uvcutMesh to ' uvcutMeshfn])
    save(uvcutMeshfn, 'uvcutMesh') ;
else
    disp('Loading uvcutMesh from disk...')
    load(uvcutMeshfn, 'uvcutMesh') ;
end
fprintf('Done with generating U,V gridded coords \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
