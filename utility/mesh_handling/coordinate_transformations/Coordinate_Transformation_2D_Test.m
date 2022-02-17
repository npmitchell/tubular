%% 2D Triangulation Coordinate Transformation Tests  ======================
% This is a script to test the functionality of the 'jacobian2Dto2DMesh.m'
% function in transforming piecewise-constant tensors defined on mesh faces
% between coordinate bases
% 
% Transforms a tensor defined as
% T = [x^2 + y^2, x^2 * y^3; ...
%      x^3 * y^2, x^2 - y^2 ]
% from a square domain X to a disk image U.
%
% by Dillon Cislo 2020/08/17
%==========================================================================

clear; close all; clc;

%--------------------------------------------------------------------------
% Construct a Basic 2D Triangulation
%--------------------------------------------------------------------------

% A triangulation of the unit disk ----------------------------------------
% TR = diskTriangulation(30);
% 
% F = TR.ConnectivityList; % Face connectivity list
% X = TR.Points; % Vertex coordinate list
% E = TR.edges; % Edge connectivity list

% A rectilinear grid ------------------------------------------------------

% Create basic grid points in 2D
L = 2; % The length of the domain
W = 2; % The width of the domain
NL = 50; % The number of vertices along the length of the cylinder
NW = round( W * NL ./ L); % The number of vertices along the width
[X, Y] = meshgrid( linspace(0, L, NL), linspace(0, W, NW) );

% Center on (0,0)
X = X-L/2;
Y = Y-W/2;

% Triangulate grid points in 2D
TR = delaunayTriangulation([X(:), Y(:)]);

F = bfs_orient(TR.ConnectivityList);
X = TR.Points;

TR = triangulation(F, X);

E = TR.edges;

clear L NL NW Y

% Calculate face centroids
COMX = cat( 3, X(F(:,1), :), X(F(:,2), :), X(F(:,3), :) );
COMX = mean( COMX, 3 );

%% Generate Analytic Results ==============================================

syms x y
assume(x, 'real');
assume(y, 'real');

% Construct a 2D mapping
u = x * sqrt(x^2 + y^2 - x^2 * y^2) ./ sqrt(x^2 + y^2);
v = y * sqrt(x^2 + y^2 - x^2 * y^2) ./ sqrt(x^2 + y^2);

% Construct the components of the Jacobian matrix of the transformation
% f(x): x --> u
JX2U_11 = gradient(u,x);
JX2U_12 = gradient(u,y);
JX2U_21 = gradient(v,x);
JX2U_22 = gradient(v,y);

JX2U_Mat = [JX2U_11 JX2U_12; JX2U_21 JX2U_22];

% Construct the components of the Jacobian matrix of the inverse
% transformation f^{-1}(u): u --> x
% JU2X_11 = gradient(x,u);
% JU2X_12 = gradient(x,v);
% JU2X_21 = gradient(y,u);
% JU2X_22 = gradient(y,v);
% 
% JU2X = repmat({ [JU2X_11 JU2X_12; JU2X_21 JU2X_22] }, size(F,1), 1);

% Construct the components of a tensor that will be transformed under the
% mapping f
T_11 = x^2 + y^2;
T_12 = x^2 * y^3;
T_21 = x^3 * y^2;
T_22 = x^2 - y^2;

T_Mat = [T_11 T_12; T_21 T_22];

%--------------------------------------------------------------------------
% Transform the Tensor T
%--------------------------------------------------------------------------

tensorType = [0 2];

if isequal(tensorType, [2 0])
    
    % Transform as a (2,0)-tensor
    TP_Mat = JX2U_Mat * T_Mat * (JX2U_Mat.');
    
elseif isequal(tensorType, [0 2])
    
    % Transform as a (0,2)-tensor (NOTICE THE MATRIX INVERSES)
    TP_Mat = inv(JX2U_Mat) * T_Mat * (inv(JX2U_Mat).');
    
else
    
    error('Supplied tensor type is not yet supported by this example');
    
end


clear JX2U_11 JX2U_12 JX2U_21 JX2U_22 T_11 T_12 T_21 T_22

%% Convert Symbolic Quantities to Numerical Quantities ====================

fprintf('Substituting numerical values for symbolic variables... ');

VX = X(:,1); VY = X(:,2);
FX = COMX(:,1); FY = COMX(:,2);

% The transformed mesh vertices
U = [ double(vpa(subs(u, {x,y}, {VX, VY}))), ...
    double(vpa(subs(v, {x,y}, {VX, VY}))) ];

% The Jacobian matrix
JX2U_Mat = double(vpa(subs(JX2U_Mat(:).', {x,y}, {FX, FY})));

% The original tensor
T_Mat = double(vpa(subs(T_Mat(:).', {x,y}, {FX, FY})));

% The transformed tensor
TP_Mat = double(vpa(subs(TP_Mat(:).', {x,y}, {FX, FY})));

fprintf('Done\n');

% Convert tensors on faces to cell arrays
JX2U = cell(size(F,1), 1);
T = cell(size(F,1), 1);
TP = cell(size(F,1), 1);

for f = 1:size(F,1)
    
    JX2U{f} = reshape(JX2U_Mat(f,:), 2, 2);
    T{f} = reshape(T_Mat(f,:), 2, 2);
    TP{f} = reshape(TP_Mat(f,:), 2, 2);
    
end

clear VX VY FX FY JX2U_Mat T_Mat TP_Mat

% Show transformed mesh
figure

subplot(1,2,1)
triplot(TR);
axis equal 
title('The Original Mesh (X)');

subplot(1,2,2);
triplot(triangulation(F,U));
axis equal
title('The Transformed Mesh (U)');

%% Construct Numerical Results ============================================

% Construct the Jacobian matrix on each mesh face
NJX2U = jacobian2Dto2DMesh(U, X, F);
NJU2X = jacobian2Dto2DMesh(X, U, F);

%--------------------------------------------------------------------------
% Transform the Tensor
%--------------------------------------------------------------------------

% Transform as a (2,0)-tensor ---------------------------------------------

NTP = cell(size(TP));

if isequal(tensorType, [2 0])
    
    % Transform as a (2,0)-tensor
    for f = 1:size(F,1)
        NTP{f} = NJX2U{f} * T{f} * (NJX2U{f}.');
    end
    
elseif isequal(tensorType, [0 2])
    
    % Transform as a (0,2)-tensor (NOTICE THE MATRIX INVERSES)
    for f = 1:size(F,1)
        NTP{f} = inv(NJX2U{f}) * T{f} * (inv(NJX2U{f}).');
    end
    
else
    
    error('Supplied tensor type is not yet supported by this example');
    
end

%% CHECK 1: Components of the Jacobian Matrix Against Analytic Results ===
clc; close all;

JErr = zeros(size(F,1), 1);
for f = 1:size(F,1)
    
    JErr(f) = norm(JX2U{f} - NJX2U{f}, 'fro') ./ norm(JX2U{f}, 'fro');
    
end

fprintf('RMS Relative Error in Jacobian = %0.5e\n', sqrt(mean(JErr.^2)));
fprintf('Max Relative Error in Jacobian = %0.5e\n', max(JErr));
fprintf('Median Relative Error in Jacobian = %0.5e\n', median(JErr));

patch('Faces', F, 'Vertices', U, ...
    'FaceVertexCData', JErr, 'FaceColor', 'flat', 'EdgeColor', 'none');
colorbar
title(['Numerical vs analytic Jacobian: ', ...
    '$\frac{J_{X\rightarrow U}-N[J_{X\rightarrow U}]}{|J_{X\rightarrow U}|}$'], ...
    'interpreter', 'latex')
set(gca, 'Clim', [0 0.1]);
axis equal;

%% CHECK 2: Components of the Inverse Jacobian ============================
% The inverse matrix on each face of the FORWARD Jacobian should match the
% BACKWARDS Jacobian of the inverse mapping
clc; close all;

JErr = zeros(size(F,1), 1);
for f = 1:size(F,1)
    
    JErr(f) = norm(inv(NJX2U{f})- NJU2X{f}, 'fro') ./ ...
        norm(inv(NJX2U{f}), 'fro');
    
end

fprintf('RMS Relative Error in Jacobian = %0.5e\n', sqrt(mean(JErr.^2)));
fprintf('Max Relative Error in Jacobian = %0.5e\n', max(JErr));
fprintf('Median Relative Error in Jacobian = %0.5e\n', median(JErr));

patch('Faces', F, 'Vertices', U, ...
    'FaceVertexCData', JErr, 'FaceColor', 'flat', 'EdgeColor', 'none');
colorbar
title(['Numerical vs analytic Jacobian: ', ...
    '$\frac{\left(J^{-1}_{X\rightarrow U}\right)^{-1}-J_{X\rightarrow U}}{|J_{X\rightarrow U}|}$'], ...
    'interpreter', 'latex')
set(gca, 'Clim', [0 0.1]);
axis equal;


%% CHECK 1: Components of the Transformed Tensor ==========================
clc; close all;

TErr = zeros(size(F,1), 1);
for f = 1:size(F,1)
    
    TErr(f) = norm(TP{f} - NTP{f}, 'fro') ./ norm(TP{f}, 'fro');
    
end

fprintf('RMS Relative Error = %0.5e\n', sqrt(mean(TErr.^2)));
fprintf('Max Relative Error = %0.5e\n', max(TErr));
fprintf('Median Relative Error= %0.5e\n', median(TErr));

patch('Faces', F, 'Vertices', U, ...
    'FaceVertexCData', TErr, 'FaceColor', 'flat', 'EdgeColor', 'none');
colorbar
tstr = '$T = \left[\begin{array}{cc}x^2 + y^2 & x^2 y^3 \\ x^3 y^2 & x^2 - y^2 \end{array}\right]$' ;
title(['tensor error, with ' tstr], 'interpreter', 'latex') 
set(gca, 'Clim', [0 0.1]);
axis equal;


    
