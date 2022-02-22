%% Induced Metric Functionality Test ======================================

clear; close all; clc;

% A basic triangulation of the unit disk
diskTri = diskTriangulation(30);

F = diskTri.ConnectivityList;
V2D = diskTri.Points;

%--------------------------------------------------------------------------
% Construct a Gaussian Surface
%--------------------------------------------------------------------------

AA = 1; % Maximum amplitude of the Gaussian
SS = 0.1; % Controls the width of the Gaussian
RR = 1; % Maximum radial distance in domain of parameterization

% The Z-coordinate of the surface
gauss3D = @(x,y) AA .* ( exp( -SS .* (x.^2 + y.^2) ) - ...
    exp( -SS .* RR.^2 ) ) ./ ( 1 - exp( -SS .* RR.^2 ) );

% The derivative of the Z-coordinate with respect to x
gaussDx = @(x,y) 2 .* AA .* SS .* x .* exp(-SS .* (x.^2 + y.^2)) ./ ...
    (exp( -SS .* RR.^2) - 1);

% The derivative of the Z-coordinate with respect to y
gaussDy = @(x,y) 2 .* AA .* SS .* y .* exp(-SS .* (x.^2 + y.^2)) ./ ...
    (exp(-SS .* RR.^2) - 1);

V3D = [ V2D, gauss3D( V2D(:,1), V2D(:,2) ) ];

%--------------------------------------------------------------------------
% Generate a Discrete Induced Metric
%--------------------------------------------------------------------------

g = inducedMetric(F, V3D, V2D);

%% CHECK 1: Does the Discrete Metric Return 3D Edge Lengths?

% The vertex IDs definig edges
E = diskTri.edges;

% The 3D directed edge vectors
e13D = V3D(F(:,3), :) - V3D(F(:,2), :);
e23D = V3D(F(:,1), :) - V3D(F(:,3), :);
e33D = V3D(F(:,2), :) - V3D(F(:,1), :);

% The 3D edge lengths
L1 = sqrt( sum( e13D.^2, 2 ) );
L2 = sqrt( sum( e23D.^2, 2 ) );
L3 = sqrt( sum( e33D.^2, 2 ) );

% The squared edge lengths
L = [ L1 L2 L3 ];
L = L.^2;

% The 2D directed edge vectors
e12D = V2D(F(:,3), :) - V2D(F(:,2), :);
e22D = V2D(F(:,1), :) - V2D(F(:,3), :);
e32D = V2D(F(:,2), :) - V2D(F(:,1), :);

LTest = zeros(size(L));
for i = 1:size(F,1)
    
    LTest(i,1) = e12D(i,:) * g{i} * e12D(i,:).';
    LTest(i,2) = e22D(i,:) * g{i} * e22D(i,:).';
    LTest(i,3) = e32D(i,:) * g{i} * e32D(i,:).';
    
end

LErr = abs(L(:)-LTest(:));
fprintf('Maximum edge length error = %0.5e\n', max(LErr));

%% CHECK 2: Does the Discrete Metric Match Analytic Results?

% The centroids of the pullback faces
COM = mean( cat(3, V2D(F(:,1), :), V2D(F(:,2), :), V2D(F(:,3), :)), 3);

% Construct analytic surface derivatices at face centroids
DFx = [ ones(size(COM(:,1))), zeros(size(COM(:,1))), ...
    gaussDx( COM(:,1), COM(:,2) ) ];

DFy = [ zeros(size(COM(:,1))), ones(size(COM(:,1))), ...
    gaussDy( COM(:,1), COM(:,2) ) ];

% Construct an analytic first fundamental form and calculate error
gTrue = cell(size(F,1), 1);
gErr = zeros(size(F,1), 1);
for i = 1:size(F,1)
    
    gTrue{i} = [ dot(DFx(i,:), DFx(i,:), 2), ...
        dot(DFx(i,:), DFy(i,:), 2); ...
        dot(DFy(i,:), DFx(i,:), 2), ...
        dot(DFy(i,:), DFy(i,:), 2) ];
    
    gErr(i) = norm( gTrue{i} - g{i}, 'fro' ) ./ norm(gTrue{i}, 'fro');
    
end

fprintf('\nMaximum error = %0.5e\n', max(gErr));
fprintf('Median error = %0.5e\n', median(gErr));
fprintf('RMS error = %0.5e\n', sqrt(mean(gErr.^2)));

trisurf( triangulation(F, V3D), 'FaceColor', 'flat', ...
    'EdgeColor', 'none', 'FaceVertexCData', gErr );
axis equal
colorbar
set(gca, 'Clim', [0 1]);





