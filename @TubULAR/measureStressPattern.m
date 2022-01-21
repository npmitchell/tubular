function measureStressPattern(QS, options)
%measureStressPattern(QS, options)
%   Measure required stress pattern from Helm-Hodge decomp potential fields
%
% Parameters
% ----------
%
%
% Returns
% -------
%
% NPMitchell 2020

%% IN PLANE STRESS PATTERN
% sigma = (d^2 phi / da db)  
%          + 1/2 [ (d^2 psi / da* db) + (d^2 psi / da* db)] 
%          + harmonic derivs
% 
% active force:
%  f= 2 mu { 1/2 [laplace(v) + grad^j (\grad_k v^k)
%     + (b^{ij} b_i^m - 2 H b^{mj}) v_m]
%     - (grad_i vN) b^{ij} - 2 vN grad^j H }
% 
% b_ij is from fundamentalForm.m
% 
% need b^{ij} for b^{ij}b_i^m*v_m
%

%% Unpack QS
nU = QS.nU ;
nV = QS.nV ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
close all


%% Load velocities
% vvt : tangential velocity on vertices
% vft : tangential velocity on faces
% vn  : normal velocity on vertices
QS.getVelocityAverage('vv') ;
vf = QS.velocityAverage.vf ;
vv = QS.velocityAverage.vv ;

%% For this timepoint
%% Load mesh 
tp = 60 ;
tidx = QS.xp.tIdx(tp) ;
meshfn = sprintf(QS.fullFileBase.spcutMeshSmRSC, tp) ;
load(meshfn, 'spcutMeshSmRSC') ;
mesh = spcutMeshSmRSC ;
mesh.nU = nU ;
cutMesh = cutRectilinearCylMesh(mesh) ;

%% DEC
dec3d = DiscreteExteriorCalculus(mesh.f, mesh.v) ;

[V2F_3d, F2V_3d] = meshAveragingOperators(mesh.f, mesh.v) ;
H3d = sum(mesh.vn .* dec3d.laplacian(mesh.v), 2) * 0.5 ;
H2d = H3d ;
H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;

[V2F_2d, F2V_2d] = meshAveragingOperators(cutMesh.f, cutMesh.v) ;
% metric and 2ndFundForm are #Fx1 cell arrays
[gg, bb] = constructFundamentalForms(cutMesh.f, cutMesh.v, cutMesh.u) ;
[g3d, b3d] = constructFundamentalForms(mesh.f, mesh.v, mesh.u) ;


%% Velocity for this TP
vf_ii = squeeze(vf(tidx, :, :)) ;
vf_ii(:, 2) = - vf_ii(:, 2)
[vfn, vft] = resolveTangentNormalVector(mesh.f, mesh.v, vf_ii) ;
vvt = F2V_3d * vft ;
% Normal velocity is vertex-normal component of vertex vels
vv_ii = squeeze(vv(tidx, :, :)) ;
vv_ii(:, 2) = - vv_ii(:, 2) ;
vv3d = vv_ii(1:end-nU, :) ;
vvn = dot(mesh.vn, vv_ii(1:nU*(nV-1), :), 2) ;

%% Check it
% imagesc(reshape(vvn, [nU, nV-1]))
subplot(2, 2, 1)
trisurf(triangulation(mesh.f, mesh.v), vv3d(:, 1), 'edgecolor','none')
caxis([-1, 1]); xlabel('x') ; ylabel('y'); zlabel('z'); axis equal ;
subplot(2, 2, 2)
trisurf(triangulation(mesh.f, mesh.v), vv3d(:, 2), 'edgecolor','none')
caxis([-1, 1]); xlabel('x') ; ylabel('y'); zlabel('z'); axis equal ;
subplot(2, 2, 3)
trisurf(triangulation(mesh.f, mesh.v), vv3d(:, 3), 'edgecolor','none')
caxis([-1, 1]); xlabel('x') ; ylabel('y'); zlabel('z'); axis equal ;
subplot(2, 2, 4)
trisurf(triangulation(mesh.f, mesh.v), vvn, 'edgecolor','none')
caxis([-1, 1]); xlabel('x') ; ylabel('y'); zlabel('z'); axis equal ;
quiver3(mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), vv3d

figure;
subplot(2, 2, 1)
trisurf(triangulation(mesh.f, mesh.v), mesh.vn(:, 1), 'edgecolor','none')
caxis([-1, 1]); xlabel('x') ; ylabel('y'); zlabel('z'); axis equal ;
subplot(2, 2, 2)
trisurf(triangulation(mesh.f, mesh.v), mesh.vn(:, 2), 'edgecolor','none')
caxis([-1, 1]); xlabel('x') ; ylabel('y'); zlabel('z'); axis equal ;
subplot(2, 2, 3)
trisurf(triangulation(mesh.f, mesh.v), mesh.vn(:, 3), 'edgecolor','none')
caxis([-1, 1]); xlabel('x') ; ylabel('y'); zlabel('z'); axis equal ;


%% CHECKS
% face barycenters
bc = barycenter(cutMesh.u, cutMesh.f) ;
% Check that all g3d == g2d up until last faces
difg = zeros(length(gg), 2, 2) ;
for qq = 1:length(gg)
    difg(qq, :, :) = gg{qq} - g3d{qq} ;
end
scatter(bc(:, 1), bc(:, 2), 5, difg(:, 1, 1))
scatter(bc(:, 1), bc(:, 2), 5, difg(:, 1, 2))
caxis([-1, 1])
colormap(bwr256)

%% Take symmetrized gradient of tangent velocity field
% Gradient eats vertex vector field
dv = dec3d.gradient(vvt) ;
% Laplacian eats face vector field or primal 0-form (or 1 form, but let's 
% not go there)
%  --> returns 
lapv = dec3d.laplacian(vft) ;
% Divergence eats #FxD dual vector field or #Ex1 primal 1-form
%  --> returns #Vx1 primal 0-form
divv = dec3d.divergence(vft) ;

% Sum together
%  f= mu { laplace(v) + grad^j (\grad_k v^k)
%     + (b^{ij} b_i^m - 2 H b^{mj}) v_m
%     - 2 (grad_i vN) b^{ij} - 4 vN grad^j H }
%  term1 = laplace(v) + grad^j (\grad_k v^k)
%  term2 = (b^{ij} b_i^m v_m
%  term3 = - 2 H b^{mj} v_m
%  term4 = - 2 (grad_i vN) b^{ij}
%  term5 = - 4 vN grad^j H 
term1 = lapv + DEC.gradient(divv) ;
term2 = bb * dot(inv(gg) * bb, vft) ; 
term3 = -2 * H3d * inv(gg) * (inv(gg) * bb) * vft ;
term4 = -2 * dot(dec3d.gradient(vn), inv(gg) * (inv(gg) * bb)) ;
term5 = -4 * vn * dec3d.grad(H3d) ;

af = term1 + term2 + term3 + term4 + term5 ;




