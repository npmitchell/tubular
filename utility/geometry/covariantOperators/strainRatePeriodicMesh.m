function [strainrate, tre, dev, theta, outputStruct] = ...
    strainRatePeriodicMesh(cutMesh, vf, options)
%strainRateMesh(cutMesh, vf, options)
%   Compute strain rate tensor on faces of a mesh given the velocities on
%   faces, vf. Tiles the mesh before making computations. If tiled mesh is
%   not supplied in options, then it is computed. Assumes periodic in Y
%   direction only. Returns strain rate tensors on each face, the trace
%   (note this is the full trace, not 1/2 trace) of the strain rate in the
%   embedding Tr[inv(g) epsilon], the deviatoric component magnitude
%   (Frobenius norm of the deviator, which is
%   A = epsilon - (1/2)Tr[inv(g) epsilon] g, the angle of elongation of the
%   deviator in the embedding coordinates relative to the projection of the
%   zeta direction in the embedding coordinates, and additional output. 
%   Additional output returns angles theta_pb in the coordinates
%   of the pullback space (instead of in the embedding space relative
%   to the embedding-space projection of zeta_hat), the scaling factors
%   bondDxDy used to scale the eigenvectors in the theta determination, the
%   fundamental forms fundForms, and the symmetrized gradient of the
%   velocities.
%
% todo: generalize to tauri and other topologies
% 
% Parameters
% ----------
% cutMesh : struct with fields 
%   Rectilinear cutMesh with pullback space u and embedding space v. 
%   f : #faces x 3 int array
%       mesh face connectivity list 
%   v : #vertices x 3 float array
%       embedding space mesh vertices of the cut mesh
%   u : #vertices x 2 float array
%       pullback mesh vertices of the cut mesh
%   pathPairs : P x 2 int array
%       path along boundary that is glued together (identified) in the
%       closed form of the mesh
% vf : #faces x 3 float array
%   velocity vectors defined on faces
% options : optional struct with fields
%   mesh : struct with fields f, u, v, like cutMesh but without seam
%   debug : bool (default=false)
%       view intermediate results for debugging
%   pullbackTheta : bool (default=false)
%       return strain deviator elongation angles theta in pullback
%       coordinates rather than in embedding coordinates (angles relative
%       to the embedded projections (ie push forwards) of unit vectors 
%       vec{du} and vec{dv} from pullback space into 3D). Default behavior
%       is to return theta s as angles from the du vector pushed into 3D.
%   tiledMesh
%
% Returns
% -------
% strainrate : #faces x 1 cell array of 2x2 float matrices
% tre : #faces x 1 float array
% dev : #faces x 1 float array
% theta : #faces x 1 float array
% outputStruct
%   fundForms : struct with fields
%       gg : first fundamental form on each face
%       bb : second fundamental form on each face
%   bondDxDy : struct with fields
%       dx : length of dx in embedding space / length of dx in pullback space
%       dy : length of dy in embedding space / length of dy in pullback space
%   theta_pb : #faces x 1 float array
%       elongation angle, theta, in pullback space (differs from
%       theta by scaling of the eigenvectors by (dx, dy) returned in
%       bondDxDy
%   dvij : #faces x 1 cell array of 2x2 float matrix
%       symmetrized spatial velocity gradient 
%
% See Also
% --------
% strainRateMesh -- for non-periodic boundary conditions / open meshes
%
% NPMitchell 2020

%% Default options
debug = false ;
mesh = [] ;

%% Unpack options
if nargin > 2
    if isfield(options, 'debug')
        debug = options.debug ;
    end
    if isfield(options, 'mesh')
        mesh = options.mesh ;
    end
end

%% Input checks
try
    assert(size(vf, 1) == size(cutMesh.f, 1)) 
catch
    error(['velocity field must be supplied on faces, but ', ...
        'size(vf) = ', num2str(size(vf)), ' while size(cutMesh.f) = ', ...
        num2str(size(cutMesh.f))])
end

%% Tile the mesh and the associated velocity vectors to triple cover
tileCount = [1, 1] ;
[ TF, TV2D, TV3D, ~ ] = tileAnnularCutMesh( cutMesh, tileCount ) ;
% The face list should now be 3x the original size
assert(size(TF, 1) == 3 * size(cutMesh.f, 1)) ;
% Triple the face velocities
nF = size(vf, 1) ;         % number of faces == number of velocities
Tvf = zeros(3 * nF, 3) ;
Tvf(1:nF, :) = vf ;
Tvf((nF + 1):(2 * nF), :) = vf ;
Tvf((2 * nF + 1):(3 * nF), :) = vf ;

%% Decompose velocity into components and compute cov. derivative
[v0n, v0t] = resolveTangentNormalVector(TF, TV3D, Tvf) ;
[~, dvi] = vectorCovariantDerivative(v0t, TF, TV3D, TV2D) ;
dvij = cell(nF, 1) ;
for qq = 1:nF
    idx = nF + qq ;
    dvij{qq} = 0.5 * ( dvi{idx} + dvi{idx}' ) ;
end

%% Compute the fundamental forms
[gg, ~] = constructFundamentalForms(TF, TV3D, TV2D) ;
% Use tiled, open mesh (glued seam) to compute the second fundamental form
[~, bb] = constructFundamentalForms(TF, TV3D, TV2D) ;

%% Strain rate tensor
strainrate = cell(nF, 1) ;
for qq = 1:nF
    strainrate{qq} = dvij{qq} - v0n(qq + nF) .* bb{qq + nF} ;
end

%% Find dx and dy -- lengths of projected / lengths of pullback
[~, dbonds, ~] = labelRectilinearMeshBonds(cutMesh, options) ;
dx = vecnorm(dbonds.realSpace.u, 2, 2) ./ vecnorm(dbonds.baseSpace.u, 2, 2) ;
dy = vecnorm(dbonds.realSpace.v, 2, 2) ./ vecnorm(dbonds.baseSpace.v, 2, 2) ;

%% Metric strain -- separate trace and deviatoric strain comp, angle
tre = zeros(size(strainrate, 1), 1) ;  % traceful dilation
dev = zeros(size(strainrate, 1), 1) ;  % deviatoric magnitude
theta = zeros(size(strainrate, 1), 1) ;  % angle of elongation
for qq = 1:size(strainrate, 1)
    eq = strainrate{qq} ;
    gq = gg{qq} ;
    
    %% Trace / deviator / theta
    % take angle of deviator to be in embedding space, angle relative
    % to the embedding-space projection of zeta_hat. 
    [tre(qq), dev(qq), theta(qq), theta_pb(qq)] = ...
        traceDeviatorPullback(eq, gq, dx(qq), dy(qq)) ;
end
% Modulo is not necessary
% theta = mod(theta, pi) ;
% theta_pb = mod(theta_pb, pi) ;

%% Output the fundamental forms, embedding bond lengths, and theta_pb
if nargout > 4
    % fundamental forms
    fundForms = struct() ;
    fundForms.gg = gg ;
    fundForms.bb = bb ;
    
    % ratio of embedding bond lengths to pullback bond lengths
    bondDxDy = struct() ;
    bondDxDy.dx = dx ;
    bondDxDy.dy = dy ;

    % theta in pullback space -- theta_pb
    % symmetrized velocity gradient -- dvij
    
    % Pack them all into output struct
    outputStruct.fundForms = fundForms ;
    outputStruct.bondDxDy = bondDxDy ;
    outputStruct.theta_pb = theta_pb ;
    outputStruct.dvij = dvij ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debug -- check results against DEC
if debug
    for qq = 1:size(dvi,1)
        tre(qq) = trace(inv(gg{qq}) * dvij{qq}) - v0n(qq) .* Hf(qq) ;
        checkH(qq) = trace(inv(gg{qq}) * bb{qq}) ;        % == 2 * Hf(qq) ; 
        checkdiv(qq) = trace(inv(gg{qq}) * dvij{qq}) ;    % == divf(qq) ; 
    end

    % Compare directly
    clf;
    subplot(2, 1, 1)
    plot(checkdiv, divf, '.')
    hold on;
    plot(checkdiv, checkdiv, 'k--')
    axis equal
    xlabel('Tr$\nabla_i v_j$', 'Interpreter', 'Latex')
    ylabel('$\nabla \cdot \mathbf{v}$', 'Interpreter', 'Latex')
    subplot(2, 1, 2)
    plot(checkH, 2* Hf, 'o') 
    hold on;
    plot(checkH, checkH, 'k--') 
    axis equal
    xlabel('Tr$b_{ij}$', 'interpreter', 'latex')
    ylabel('$2H$', 'interpreter', 'latex')

    if preview
        %% Check Mean curvature
        subplot(2, 2, 1)
        trisurf(triangulation(mesh.f, mesh.v), checkH, 'edgecolor', 'none')
        axis equal; caxis([-.1, .1]); colorbar() ;
        title('Tr$\left[g^{-1} b\right]$', 'interpreter', 'latex')
        subplot(2, 2, 2)
        trisurf(triangulation(mesh.f, mesh.v), 2 * Hf, 'edgecolor', 'none')
        axis equal; caxis([-.1, .1]); colorbar() ;
        title('$2H$', 'interpreter', 'latex')
        subplot(2, 1, 2)
        trisurf(triangulation(mesh.f, mesh.v), checkH - 2 * Hf, 'edgecolor', 'none')
        axis equal; caxis([-.1, .1]); colorbar() ;
        title('Tr$[g^{-1}b] - 2H$', 'interpreter', 'latex')
        colormap(bwr)

        %% Check div(v)
        clf
        % set color limit, clim
        clim = max(2*std(abs(checkdiv)), 2*std(abs(divf))) ;
        subplot(2, 2, 1)
        trisurf(triangulation(mesh.f, mesh.v), checkdiv, ...
            'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        checkDivStr = '$\frac{1}{2}$Tr$\left[g^{-1} \left(\nabla_i v_j + \nabla_j v_i \right)\right]$' ;
        title(checkDivStr, 'interpreter', 'latex')
        subplot(2, 2, 2)
        trisurf(triangulation(mesh.f, mesh.v), divv, 'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title('$\nabla \cdot \mathbf{v}_{\parallel}$', 'interpreter', 'latex')
        subplot(2, 1, 2)
        trisurf(triangulation(mesh.f, mesh.v), checkdiv - divf, ...
            'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title([checkDivStr ' $-\nabla \cdot \mathbf{v}_{\parallel}$'], ...
            'interpreter', 'latex')
        colormap(bwr)

        %% Check velocities (raw vs averaged)
        clf
        clim = 0.1 ;
        subplot(2, 2, 1)
        trisurf(triangulation(mesh.f, mesh.v), vss(:,1)-vs(:, 1), ...
            'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        checkDivStr = '$\frac{1}{2}$Tr$\left[g^{-1} \left(\nabla_i v_j + \nabla_j v_i \right)\right]$' ;
        title(checkDivStr, 'interpreter', 'latex')
        xlabel('x'); ylabel('y'); zlabel('z')
        title('$v_x$')
        subplot(2, 2, 2)
        trisurf(triangulation(mesh.f, mesh.v), vss(:,2)-vs(:, 2), 'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title('$\nabla \cdot \mathbf{v}_{\parallel}$', 'interpreter', 'latex')
        xlabel('x'); ylabel('y'); zlabel('z')
        title('$v_y$')
        subplot(2, 1, 2)
        trisurf(triangulation(mesh.f, mesh.v), vss(:,3)-vs(:, 3), ...
            'edgecolor', 'none')
        axis equal; caxis([-clim, clim]); colorbar() ;
        title([checkDivStr ' $-\nabla \cdot \mathbf{v}_{\parallel}$'], ...
            'interpreter', 'latex')
        xlabel('x'); ylabel('y'); zlabel('z')
        title('$v_z$')
        colormap(bwr)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
