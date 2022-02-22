%% COVARIANT_DERIVATIVE_TEST ==============================================
%
%   This script tests the validity of a construction for calculating the
%   covariant derivative of vector fields on surfaces
%
%   by Dillon Cislo 03/10/2020
%   cylinder check: NPMitchell 2020
%
%==========================================================================
clear; close all; clc;

% Add necessary directories to the path -----------------------------------
addpath(genpath('/home/dillon/Documents/MATLAB/GrowSurface/MeshCreation'));
addpath(genpath('/home/dillon/Documents/MATLAB/gptoolbox'));
addpath(genpath('/home/dillon/Documents/MATLAB/coordinate_transformations'));

% NPMitchell paths
addpath_recurse('/mnt/data/code/gut_matlab/mesh_handling/')


%% ************************************************************************
% *************************************************************************
%           PART ZERO: TANGENT VECTOR FIELD ON CYLINDER
% *************************************************************************
% *************************************************************************
nU = 100 ;
nV = 100 ;
zz = linspace(0, 10, nU) ;
phi = linspace(0, 2*pi, nV) ;
[zz, phi] = meshgrid(zz, phi) ;
zz = zz' ;
phi = phi' ;
uv = [ zz(:), phi(:) ] ;
vv = [ zz(:), cos(phi(:)), sin(phi(:)) ] ; 
faces = defineFacesRectilinearGrid(uv, 100, 100) ;
cutMesh.f = faces ;
cutMesh.v = vv ;
cutMesh.u = uv ;
cutMesh.nU = nU ;
cutMesh.nV = nV ;
mesh = glueCylinderCutMeshSeam(cutMesh) ;

% Obtain mean curvature H for checking against trace(b_ij)
DEC = DiscreteExteriorCalculus(mesh.f, mesh.v) ;
H3d = sum(mesh.vn .* DEC.laplacian(mesh.v), 2) * 0.5 ;
H2d = H3d ;
H2d(nU*(nV-1)+1:(nU*nV)) = H3d(1:nU) ;

% Define velocities on the vertices 
vx = (mean(mesh.v(:, 1)) - mesh.v(:, 1)).^2 ;
vs = [vx, 0 * ones(size(vx)), 0*vx] ;

% Push vectors onto faces
[V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
vf = V2F * vs ;

% Compute divergence and face-based mean curvature Hf
Hf = V2F * H3d ;
divv = DEC.divergence(vf) ;
divf = V2F * divv ;

% Compute the strain rate tensor
disp('Computing covariantDerivative')
[~, dvi] = vectorCovariantDerivative(vf, cutMesh.f, cutMesh.v, cutMesh.u) ;
dvij = cell(size(dvi)) ;
for qq = 1:length(dvi)
    dvij{qq} = 0.5 * ( dvi{qq} + dvi{qq} ) ;
end

% Compute the second fundamental form
[gg, ~] = constructFundamentalForms(cutMesh.f, cutMesh.v, cutMesh.u) ;
[~, bb] = constructFundamentalForms(mesh.f, mesh.v, mesh.u) ;

% Decompose velocity into components
[v0n, ~] = resolveTangentNormalVector(cutMesh.f, cutMesh.v, vf) ;

% Strain rate tensor
strainrate = cell(size(dvi)) ;
tre = zeros(size(dvi)) ;
checkH = zeros(size(dvi)) ;
checkdiv = zeros(size(dvi)) ;
for qq = 1:length(dvi)
    strainrate{qq} = dvij{qq} - v0n(qq) .* bb{qq} ;
    tre(qq) = trace(inv(gg{qq}) * dvij{qq}) - v0n(qq) .* Hf(qq) ;
    checkH(qq) = trace(inv(gg{qq}) * bb{qq}) ;        % == 2 * Hf(qq) ; 
    checkdiv(qq) = trace(inv(gg{qq}) * dvij{qq}) ;    % == divf(qq) ; 
end

clf;
subplot(2, 1, 1)
plot(checkdiv, divf * pi, '.')
hold on;
plot(checkdiv, checkdiv, 'k--')
axis equal
xlabel('Tr$b_{ij}$', 'interpreter', 'latex')
ylabel('$2H$', 'interpreter', 'latex')
subplot(2, 1, 2)
plot(checkH, 2* Hf, 'o') 
hold on;
plot(checkH, checkH, 'k--') 
axis equal
xlabel('Tr$\nabla_i v_j$', 'Interpreter', 'Latex')
ylabel('$\nabla \cdot \mathbf{v}$', 'Interpreter', 'Latex')

%% Check Mean curvature
clim = max(max(abs(checkH)), 2 * max(abs(Hf))) ;
subplot(2, 2, 1)
trisurf(triangulation(mesh.f, mesh.v), checkH, 'edgecolor', 'none')
axis equal; caxis([-clim, clim]); colorbar() ;
title('Tr$\left[g^{-1} b\right]$', 'interpreter', 'latex')
subplot(2, 2, 2)
trisurf(triangulation(mesh.f, mesh.v), 2 * Hf, 'edgecolor', 'none')
axis equal; caxis([-clim, clim]); colorbar() ;
title('$2H$', 'interpreter', 'latex')
subplot(2, 1, 2)
trisurf(triangulation(mesh.f, mesh.v), checkH - 2 * Hf, 'edgecolor', 'none')
axis equal; caxis([-clim, clim]); colorbar() ;
title('Tr$[g^{-1}b] - 2H$', 'interpreter', 'latex')
colormap(bwr)

%% Check div(v)
clf
clim = max(mean(abs(checkdiv)), mean(abs(divf))) ;
subplot(2, 2, 1)
trisurf(triangulation(mesh.f, mesh.v), checkdiv, ...
    'edgecolor', 'none')
axis equal; caxis([-clim, clim]); colorbar() ;
checkDivStr = '$\frac{1}{2}$Tr$\left[g^{-1} \left(\nabla_i v_j + \nabla_j v_i \right)\right]$' ;
title(checkDivStr, 'interpreter', 'latex')
subplot(2, 2, 2)
trisurf(triangulation(mesh.f, mesh.v), divf, 'edgecolor', 'none')
axis equal; caxis([-clim, clim]); colorbar() ;
title('$\nabla \cdot \mathbf{v}_{\parallel}$', 'interpreter', 'latex')
subplot(2, 1, 2)
trisurf(triangulation(mesh.f, mesh.v), checkdiv - divf, ...
    'edgecolor', 'none')
axis equal; caxis([-clim, clim]); colorbar() ;
title([checkDivStr ' $-\nabla \cdot \mathbf{v}_{\parallel}$'], ...
    'interpreter', 'latex')
colormap(bwr)

%% ************************************************************************
% *************************************************************************
%           PART ONE: TANGENT VECTOR FIELD ON THE UNIT SPHERE
% *************************************************************************
% *************************************************************************

clear; close all; clc;

%--------------------------------------------------------------------------
% Mesh geometry handling
%--------------------------------------------------------------------------

% Create a mesh of the unit sphere ----------------------------------------
sphereTri = sphereTriangulationVogel(5000);

F = sphereTri.ConnectivityList;
V = sphereTri.Points;

% Calculate the centroids of each face ------------------------------------
COM = mean( cat(3, V(F(:,1), :), V(F(:,2), :), V(F(:,3), :) ), 3 );

% Create face-to-vertex and vertex-to-face operators ----------------------
V2F = sparse( repmat(1:size(F,1), 3, 1), F.', 1/3, size(F,1), size(V,1) );

F2V = internalangles(V, F); % The internal angles of the mesh
F2V = sparse( F(:), repmat(1:size(F,1),1,3), F2V, size(V,1), size(F,1) );

% Create face-projection operator -----------------------------------------

% Face unit normals
FN = faceNormal(sphereTri);

P = cell( size(F,1), 1 );
for i = 1:size(F,1)
    P{i} = eye(3) - FN(i,:).' * FN(i,:);
end

%--------------------------------------------------------------------------
% Construct spherical tangent vectors
%--------------------------------------------------------------------------
X = V(:,1); Y = V(:,2); Z = V(:,3);

cosTheta = Z;
sinTheta = sin(acos(cosTheta));

phi = atan2(Y,X);
sinPhi = sin(phi);
cosPhi = cos(phi);

eTV = [ cosTheta .* cosPhi, cosTheta .* sinPhi, -sinTheta ];
ePV = [ -sinTheta .* sinPhi, sinTheta .* cosPhi, zeros(size(Z)) ];

eTF = V2F * eTV;
ePF = V2F * ePV;

%--------------------------------------------------------------------------
% Construct a primal (vertex-pased) vector field on the mesh
%--------------------------------------------------------------------------

UV = [ X .* Z .* ( Z.^2 - 1/4 ) - Y, ...
    Y .* Z .* ( Z.^2 - 1/4 ) + X, ...
    -( X.^2 + Y.^2 ) .* ( Z.^2 - 1/4 ) ];

% View results ------------------------------------------------------------
ssf = 15; % Sub-sampling factor for visualization
trisurf( sphereTri, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' );
hold on
quiver3( V(1:ssf:end, 1), V(1:ssf:end, 2), V(1:ssf:end, 3), ...
    UV(1:ssf:end, 1), UV(1:ssf:end, 2), UV(1:ssf:end, 3), ...
    1, 'LineWidth', 2, 'Color', 'm' );
hold off
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');

%--------------------------------------------------------------------------
% Construct a dual (face-based) vector field on the mesh
%--------------------------------------------------------------------------

% Average primal vector field onto faces
UF = V2F * UV;

% View results ------------------------------------------------------------
% ssf = 15; % Sub-sampling factor for visualization trisurf( sphereTri,
% 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' ); hold on quiver3(
% COM(1:ssf:end, 1), COM(1:ssf:end, 2), COM(1:ssf:end, 3), ...
%     UF(1:ssf:end, 1), UF(1:ssf:end, 2), UF(1:ssf:end, 3), ... 1,
%     'LineWidth', 2, 'Color', 'm' );
% hold off axis equal xlabel('X'); ylabel('Y'); zlabel('Z');

clear X Y Z ssf cosTheta sinTheta cosPhi sinPhi phi


%% Calculate the Covariant Derivative of the Vector Field =================

DU3D = vectorCovariantDerivative( UV, F, V );

%% Test Validity ==========================================================
clc; close all

% Construct analytic results ----------------------------------------------
Z = V(:,3);

% Covariant derivative along the theta direction
CTT = 3 .* Z.^3 - 9 .* Z ./ 4;
CTP = Z ./ sqrt(1-Z.^2);

CDtU = CTT .* eTV + CTP .* ePV;
CDtU = V2F * CDtU;
for i = 1:size(F,1)
    CDtU(i,:) = (P{i} * CDtU(i,:).').';
end

% Covariant derivative along the phi direction
CPT = -Z .* sqrt(1-Z.^2);
CPP = Z.^3 - Z/4;

CDpU = CPT .* eTV + CPP .* ePV;
CDpU = V2F * CDpU;
for i = 1:size(F,1)
    CDpU(i,:) = (P{i} * CDpU(i,:).').';
end

clear Z CTT CTP CPT CPP

% Calculate numerical results ---------------------------------------------
NDtU = zeros(size(F));
NDpU = zeros(size(F));
for i = 1:size(F,1)
    
    NDtU(i,:) = ( DU3D{i} * eTF(i,:).' ).';
    NDpU(i,:) = ( DU3D{i} * ePF(i,:).' ).';
    
end

% View results ------------------------------------------------------------

thetaErr = sqrt( sum( (CDtU-NDtU).^2, 2 ) ) ./ sqrt( sum( CDtU.^2, 2 ) );
maxThetaErr = max(thetaErr);
rmsThetaErr = sqrt( mean( thetaErr.^2 ) );
medThetaErr = median(thetaErr);

phiErr = sqrt( sum( (CDpU-NDpU).^2, 2 ) ) ./ sqrt( sum( CDpU.^2, 2 ) );
maxPhiErr = max(phiErr);
rmsPhiErr = sqrt( mean( phiErr.^2 ) );
medPhiErr = median(phiErr);

fprintf('Maximum relative error in theta direction = %0.5f\n', maxThetaErr);
fprintf('RMS relative error in theta direction = %0.5f\n', rmsThetaErr);
fprintf('Median relative error in theta direction = %0.5f\n', medThetaErr);
fprintf('Maximum relative error in phi direction = %0.5f\n', maxPhiErr);
fprintf('RMS relative error in phi direction = %0.5f\n', rmsPhiErr);
fprintf('Median relative error in phi direction = %0.5f\n', medPhiErr);

ssf = 15; % Sub-sampling factor for visualization

figure

subplot(1,2,1)
trisurf( sphereTri, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' );
hold on
quiver3( COM(1:ssf:end, 1), COM(1:ssf:end, 2), COM(1:ssf:end, 3), ...
    CDtU(1:ssf:end, 1), CDtU(1:ssf:end, 2), CDtU(1:ssf:end, 3), ...
    1, 'LineWidth', 2, 'Color', 'm' );
hold off
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Analytic Derivative in Theta');

subplot(1,2,2)
trisurf( sphereTri, 'FaceColor', 'flat', 'FaceVertexCData', thetaErr, ...
    'EdgeColor', 'none' );
axis equal
colorbar
set(gca, 'Clim', [0 1]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Relative Error in Theta Derivative');

figure

subplot(1,2,1)
trisurf( sphereTri, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' );
hold on
quiver3( COM(1:ssf:end, 1), COM(1:ssf:end, 2), COM(1:ssf:end, 3), ...
    CDpU(1:ssf:end, 1), CDpU(1:ssf:end, 2), CDpU(1:ssf:end, 3), ...
    1, 'LineWidth', 2, 'Color', 'm' );
hold off
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Analytic Derivative in Phi');

subplot(1,2,2)
trisurf( sphereTri, 'FaceColor', 'flat', 'FaceVertexCData', phiErr, ...
    'EdgeColor', 'none' );
axis equal
colorbar
set(gca, 'Clim', [0 1]);
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Relative Error in Phi Derivative');

clear maxThetaErr rmsThetaErr maxPhiErr rmsPhiErr ssf

%% ************************************************************************
% *************************************************************************
%           PART TWO: PULLING BACK THE COVARIANT DERIVATIVE
% *************************************************************************
% *************************************************************************

clear; close all; clc;

%--------------------------------------------------------------------------
% Mesh geometry handling
%--------------------------------------------------------------------------

% Create mesh of Gaussian surface -----------------------------------------
diskTri = diskTriangulation(100);

F = diskTri.ConnectivityList;
V2D = diskTri.Points;

AA = 1; % Maximum amplitude of the Gaussian
SS = 2; % Controls the width of the Gaussian
RR = 1; % Maximum radial distance in domain of parameterization

gauss3D = @(x,y) AA .* ( exp( -SS .* (x.^2 + y.^2) ) - ...
    exp( -SS .* RR.^2 ) ) ./ ( 1 - exp( -SS .* RR.^2 ) );

V3D = [V2D gauss3D(V2D(:,1), V2D(:,2))];

gaussTri = triangulation(F, V3D);

% Face centroids in 2D
COM2D = cat(3, V2D(F(:,1),:), V2D(F(:,2),:), V2D(F(:,3),:));
COM2D = mean(COM2D, 3);

% Face centroids in 3D
COM3D = cat(3, V3D(F(:,1),:), V3D(F(:,2),:), V3D(F(:,3),:));
COM3D = mean(COM3D, 3);

% Create face-to-vertex and vertex-to-face operators ----------------------
V2F = sparse(repmat(1:size(F,1), 3, 1), F.', 1/3, size(F,1), size(V3D,1));

F2V = internalangles(V3D, F); % The internal angles of the mesh
F2V = sparse(F(:), repmat(1:size(F,1),1,3), F2V, size(V3D,1), size(F,1));

% Create face-projection operator -----------------------------------------

% Face unit normals
FN = faceNormal(gaussTri);

P = cell( size(F,1), 1 );
for i = 1:size(F,1)
    P{i} = eye(3) - FN(i,:).' * FN(i,:);
end

%--------------------------------------------------------------------------
% Construct Gaussian surface tangent vectors
%--------------------------------------------------------------------------

% The derivative of the Z-coordinate with respect to x
gaussDx = @(x,y) 2 .* AA .* SS .* x .* exp(-SS .* (x.^2 + y.^2)) ./ ...
    (exp( -SS .* RR.^2) - 1);

% The derivative of the Z-coordinate with respect to y
gaussDy = @(x,y) 2 .* AA .* SS .* y .* exp(-SS .* (x.^2 + y.^2)) ./ ...
    (exp(-SS .* RR.^2) - 1);

eXV = [ ones(size(V3D,1),1) zeros(size(V3D,1),1), ...
    gaussDx(V3D(:,1), V3D(:,2)) ];
eYV = [ zeros(size(V3D,1),1) ones(size(V3D,1),1), ...
    gaussDy(V3D(:,1), V3D(:,2)) ];

eXF = V2F * eXV; 
eYF = V2F * eYV;

%--------------------------------------------------------------------------
% Construct a primal (vertex-based) vector field on the mesh
%--------------------------------------------------------------------------

UV = eXV;

% View results ------------------------------------------------------------
ssf = 25; % Sub-sampling factor for visualization
trisurf( gaussTri, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' );
hold on
quiver3( V3D(1:ssf:end, 1), V3D(1:ssf:end, 2), V3D(1:ssf:end, 3), ...
    UV(1:ssf:end, 1), UV(1:ssf:end, 2), UV(1:ssf:end, 3), ...
    1, 'LineWidth', 2, 'Color', 'm' );
hold off
axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');

%--------------------------------------------------------------------------
% Construct a dual (face-based) vector field on the mesh
%--------------------------------------------------------------------------

% Average primal vector field onto faces
UF = V2F * UV;
for i = 1:size(F,1)
    UF(i,:) = (P{i} * UF(i,:).').';
end

% View results ------------------------------------------------------------
% ssf = 15; % Sub-sampling factor for visualization trisurf( gaussTri,
% 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' ); hold on quiver3( ...
%     COM3D(1:ssf:end, 1), COM3D(1:ssf:end, 2), COM3D(1:ssf:end, 3), ...
%     UF(1:ssf:end, 1), UF(1:ssf:end, 2), UF(1:ssf:end, 3), ... 1,
%     'LineWidth', 2, 'Color', 'm' );
% hold off axis equal xlabel('X'); ylabel('Y'); zlabel('Z');

clear X Y Z ssf

%% Calculate the Covariant Derivative of the Vector Field =================

[DU3D, DU2D] = vectorCovariantDerivative(UV, F, V3D, V2D);

% Transform the Covariant Derivative as a (0,2)-Tensor ====================
% 
% J32 = jacobian3Dto2DMesh( V2D, V3D, F ); J23 = jacobian2Dto3DMesh( V2D,
% V3D, F );
% 
% for i = 1:size(F,1)
%     
%     % Type-(0,2) Tensor % DU{i} = J23{i}.' * DU{i}.' * J23{i};
%     
%     % Type-(1,1) Tensor % DU{i} = J23{i}.' * DU{i}.' * J32{i}.';
%     
%     % Type-(2,0) Tensor % DU{i} = J32{i} * DU{i}.' * J32{i}.';
%     
% end


%% Test Validity ==========================================================
clc;

%--------------------------------------------------------------------------
% Construct Analytic Results
%--------------------------------------------------------------------------

% Type-(0,2) Tensor -------------------------------------------------------

gaussCD = @(x,y) -4 .* AA.^2 .* SS.^2 .* ...
    exp( -2 .* SS .* (x.^2 + y.^2 - RR.^2) ) ./ ...
    (1 - exp(SS .* RR.^2)).^2;

X = V3D(:,1); Y = V3D(:,2);
CDUV = zeros(size(V3D,1), 4);
CDUV(:,1) = X .* (2 .* SS .* X.^2 - 1);
CDUV(:,2) = Y .* (2 .* SS .* X.^2 - 1);
CDUV(:,3) = 2 .* SS .* X.^2 .* Y;
CDUV(:,4) = 2 .* SS .* X .* Y.^2;
CDUV = gaussCD(X,Y) .* CDUV;

CDUF = V2F * CDUV;

CDU = cell(size(F,1), 1);
for i = 1:size(F,1)
    CDU{i} = [ CDUF(i,1) CDUF(i,2); CDUF(i,3) CDUF(i,4) ];
end

clear gaussCD X Y CDUV CDUF

%--------------------------------------------------------------------------
% Calculate Errors
% -------------------------------------------------------------------------

covErr = zeros(size(F,1), 1);
for i = 1:size(F,1)
    covErr(i) = sqrt(trace((CDU{i} - DU2D{i})*(CDU{i} - DU2D{i}).')) ./ ...
        sqrt(trace(CDU{i} * CDU{i}.'));
end

fprintf('Maximum relative error = %0.5f\n', max(covErr));
fprintf('Median relative error = %0.5f\n', median(covErr));
fprintf('RMS relative error = %0.5f\n', sqrt(mean(covErr.^2)));

%--------------------------------------------------------------------------
% View Results
%--------------------------------------------------------------------------
trisurf( gaussTri, 'FaceColor', 'flat', 'FaceVertexCData', covErr, ...
    'EdgeColor', 'none' );
axis equal
colorbar
set(gca, 'Clim', [0 1]);

%% ************************************************************************
% *************************************************************************
%       PART THREE: PULLING BACK THE COVARIANT DERIVATIVE (AGAIN)
% *************************************************************************
% *************************************************************************
% The vector field in part 2 is too simple - I'm not convinced

clear; close all; clc;

%--------------------------------------------------------------------------
% Mesh geometry handling
%--------------------------------------------------------------------------

% Create conformal mesh of a hemisphere -----------------------------------
diskTri = diskTriangulation(100);

F = diskTri.ConnectivityList;
F = bfs_orient(F);
V2D = diskTri.Points;

% Split up the 2D coordinates for easier use
x = V2D(:,1);
y = V2D(:,2);

r = sqrt(x.^2 + y.^2);

X = 2 .* x ./ (1 + r.^2);
Y = 2 .* y ./ (1 + r.^2);
Z = (1 - r.^2) ./ (1 + r.^2);

V3D = [ X, Y, Z ];

sphereTri = triangulation(F, V3D);

% Face centroids in 2D
COM2D = cat(3, V2D(F(:,1),:), V2D(F(:,2),:), V2D(F(:,3),:));
COM2D = mean(COM2D, 3);

% Face centroids in 3D
COM3D = cat(3, V3D(F(:,1),:), V3D(F(:,2),:), V3D(F(:,3),:));
COM3D = mean(COM3D, 3);

% Create face-to-vertex and vertex-to-face operators ----------------------
V2F = sparse(repmat(1:size(F,1), 3, 1), F.', 1/3, size(F,1), size(V3D,1));

F2V = internalangles(V3D, F); % The internal angles of the mesh
F2V = sparse(F(:), repmat(1:size(F,1),1,3), F2V, size(V3D,1), size(F,1));

% Create face-projection operator -----------------------------------------

% Face unit normals
FN = faceNormal(sphereTri);

P = cell( size(F,1), 1 );
for i = 1:size(F,1)
    P{i} = eye(3) - FN(i,:).' * FN(i,:);
end

%--------------------------------------------------------------------------
% Construct Hemispherical Surface Tangent Vectors
%--------------------------------------------------------------------------

eXV = [ 2 .* (1 - x.^2 + y.^2), -4 .* x .* y, -4 .* x ] ./ (1 + r.^2).^2;
eYV = [ -4 .* x .* y, 2 .* (1 + x.^2 - y.^2), -4 .* y ] ./ (1 + r.^2).^2;

eXF = V2F * eXV;
eYF = V2F * eYV;

for i = 1:size(F,1)
    eXF(i,:) = (P{i} * eXF(i,:).').';
    eYF(i,:) = (P{i} * eYF(i,:).').';
end

% plotVec = eYV;
% ssf = 1;
% trisurf(sphereTri, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' );
% hold on
% quiver3(X(1:ssf:end), Y(1:ssf:end), Z(1:ssf:end), ...
%     plotVec(1:ssf:end,1), plotVec(1:ssf:end,2), plotVec(1:ssf:end,3), ...
%     1, 'LineWidth', 2, 'Color', 'm' );
% hold off
% axis equal
% xlabel('X');
% ylabel('Y');
% zlabel('Z');

%--------------------------------------------------------------------------
% Construct a primal (vertex-based) vector field on the mesh
%--------------------------------------------------------------------------

UV = [ X .* Z .* ( Z.^2 - 1/4 ) - Y, ...
    Y .* Z .* ( Z.^2 - 1/4 ) + X, ...
    -( X.^2 + Y.^2 ) .* ( Z.^2 - 1/4 ) ];

% View results ------------------------------------------------------------
% ssf = 2; % Sub-sampling factor for visualization
% trisurf( sphereTri, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' );
% hold on
% quiver3( X(1:ssf:end), Y(1:ssf:end), Z(1:ssf:end), ...
%     UV(1:ssf:end, 1), UV(1:ssf:end, 2), UV(1:ssf:end, 3), ...
%     1, 'LineWidth', 2, 'Color', 'm' );
% hold off
% axis equal
% xlabel('X'); ylabel('Y'); zlabel('Z');

%--------------------------------------------------------------------------
% Construct a dual (face-based) vector field on the mesh
%--------------------------------------------------------------------------

% Average primal vector field onto faces
UF = V2F * UV;
for i = 1:size(F,1)
    UF(i,:) = (P{i} * UF(i,:).').';
end

% View results ------------------------------------------------------------
% ssf = 15; % Sub-sampling factor for visualization trisurf( sphereTri,
% 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none' ); hold on quiver3(
% COM(1:ssf:end, 1), COM(1:ssf:end, 2), COM(1:ssf:end, 3), ...
%     UF(1:ssf:end, 1), UF(1:ssf:end, 2), UF(1:ssf:end, 3), ... 1,
%     'LineWidth', 2, 'Color', 'm' );
% hold off axis equal xlabel('X'); ylabel('Y'); zlabel('Z');

%% Calculate the Covariant Derivative of the Vector Field =================

[DU3D, DU2D] = vectorCovariantDerivative(UV, F, V3D, V2D);

% Transform the Covariant Derivative as a (0,2)-Tensor ====================
% 
% J32 = jacobian3Dto2DMesh( V2D, V3D, F );
% J23 = jacobian2Dto3DMesh( V2D, V3D, F );
%
% for i = 1:size(F,1)
%     
%     % Type-(0,2) Tensor % DU{i} = J23{i}.' * DU{i}.' * J23{i};
%     
%     % Type-(1,1) Tensor % DU{i} = J23{i}.' * DU{i}.' * J32{i}.';
%     
%     % Type-(2,0) Tensor % DU{i} = J32{i} * DU{i}.' * J32{i}.';
%     
% end

%% Test Validity ==========================================================
clc;

%--------------------------------------------------------------------------
% Construct Analytic Results
%--------------------------------------------------------------------------

% Type-(0,2) Tensor -------------------------------------------------------

% [ DUx/Dx DUx/Dy DUy/Dx DUy/Dy ]
% CDUV = zeros(size(V2D,1), 4);

CDUV = zeros(size(F,1), 4);
x = COM2D(:,1); y = COM2D(:,2);

CDUV(:,1) = 3 - 3 .* x.^6 - 13 .* y.^2 + 13 .* y.^4 - 3 .* y.^6 - ...
    9 .* x.^4 .* (-5 + y.^2) + x.^2 .* (-45 + 58 .* y.^2 - 9 .* y.^4);

CDUV(:,2) = 4 .* (x.^6 + 8 .* x.^3 .* y + 8 .* x .* y .* (-1 + y.^2) + ...
    (-1 + y.^2) .* (1 + y.^2).^2 + x.^4 .* (1 + 3 .* y.^2) + ...
    x.^2 .* (-1 + 2 .* y.^2 + 3 .* y.^4));

CDUV(:,3) = -4 .* (x.^6 - 8 .* x.^3 .* y - 8 .* x .* y .* (-1 + y.^2) + ...
    (-1 + y.^2) .* (1 + y.^2).^2 + x.^4 .* (1 + 3 .* y.^2) + ...
    x.^2 .* (-1 + 2 .* y.^2 + 3 .* y.^4));

CDUV(:,4) = 3 - 3 .* x.^6 - 45 .* y.^2 + 45 .* y.^4 - 3 .* y.^6 + ...
    x.^4 .* (13 - 9 .* y.^2) + x.^2 .* (-13 + 58 .* y.^2 - 9 .* y.^4);

CDUV = CDUV ./ repmat( (1 + x.^2 + y.^2).^5, 1, 4 );

% CDUF = V2F * CDUV;
CDUF = CDUV;

CDU = cell(size(F,1), 1);
for i = 1:size(F,1)
    CDU{i} = [ CDUF(i,1) CDUF(i,3); CDUF(i,2) CDUF(i,4) ];
end

clear CDUV CDUF

%--------------------------------------------------------------------------
% Calculate Errors
%--------------------------------------------------------------------------

covErr = zeros(size(F,1), 1);
for i = 1:size(F,1)
    covErr(i) = norm( CDU{i} - DU2D{i}, 'fro' ) ./ norm( CDU{i}, 'fro' );
end

fprintf('Maximum relative error = %0.5f\n', max(covErr));
fprintf('Median relative error = %0.5f\n', median(covErr));
fprintf('RMS relative error = %0.5f\n', sqrt(mean(covErr.^2)));

%--------------------------------------------------------------------------
% View Results
%--------------------------------------------------------------------------
trisurf( sphereTri, 'FaceColor', 'flat', 'FaceVertexCData', covErr, ...
    'EdgeColor', 'none' );
axis equal
colorbar
set(gca, 'Clim', [0 1]);










