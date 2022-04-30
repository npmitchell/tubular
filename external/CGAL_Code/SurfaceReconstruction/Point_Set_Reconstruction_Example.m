%% Point Set Surface Reconstruction Example ===============================
%
%   This is a script to test the functionality of a pipeline for generating
%   mesh triangulations of noise sparse point sets
%

addpath(genpath('/data/code/imsane_for_git_dillon_20200316/imsane/'));
addpath(genpath('/data/code/RicciFlow_MATLAB'));

%% Generate Initial Sparse Point Set ======================================
clear; close all; clc;

% % Basic points are part of a spherical cap
% [X, Y, Z] = sphere(25);
% P = [X(:), Y(:), Z(:)];
% P(P(:,3) < 0.2, :) = [];

% Basic parts are the surface of a cylinder
[X, Y, Z] = cylinder(0.25*ones(1, 25), 60);
P = [X(:), Y(:), Z(:)];

% Add noise to point positions
% rng(1,'twister');
% P = awgn(P, 30);

scatter3(P(:,1), P(:,2), P(:,3), 'filled', 'b');
axis equal

clear X Y Z

%% Infer Point Set Normals ================================================
close all; clc;

param = struct();
param.estimation_procedure = 2; % Uses PCA normal estimation
param.number_of_neighbors = 18;
orient_neighbors = 12;

P0 = P;
[ PN, P, oriented ] = point_set_normals( P, param, ...
    orient_neighbors );

% Check to see if any points were removed during normal orientation
if oriented
    warning(['Point set normal orientation procedure', ...
        'removed points from the input set!']);
end

% View results ------------------------------------------------------------
ssf = false( size(P,1), 1 );
sst = true( size(P, 1), 1 );
% sst = false( size(P,1), 1 );
% sst(1:50:end) = true;

scatter3(  P([sst; ssf; ssf]), P([ssf; sst; ssf]), ...
    P([ssf; ssf; sst]), 'filled' );
hold on
quiver3( P([sst; ssf; ssf]), P([ssf; sst; ssf]), ...
    P([ssf; ssf; sst]), ...
    PN([sst; ssf; ssf]), PN([ssf; sst; ssf]), ...
    PN([ssf; ssf; sst]), ...
    1, 'LineWidth', 2 );
axis equal
hold off

clear param orient_neighbors oriented ssf sst

%% Upsample Point Set =====================================================
close all; clc;

numPoints = 2000; % The number of output points
sharpnessAngle = 90; % Controls the sharpness of the results
edgeSensitvity = 0; % Controls sampling density near inferred edges
neighborRadius = eps; % Initial neighbor radius

[PU, PNU] = upsample_point_set( P, PN, numPoints, sharpnessAngle, ...
    edgeSensitvity, neighborRadius );

% View results ------------------------------------------------------------
ssf = false( size(PU,1), 1 );
sst = true( size(PU, 1), 1 );
% sst = false( size(P,1), 1 );
% sst(1:50:end) = true;

scatter3(  PU([sst; ssf; ssf]), PU([ssf; sst; ssf]), ...
    PU([ssf; ssf; sst]), 'filled' );
hold on
quiver3( PU([sst; ssf; ssf]), PU([ssf; sst; ssf]), ...
    PU([ssf; ssf; sst]), ...
    PNU([sst; ssf; ssf]), PNU([ssf; sst; ssf]), ...
    PNU([ssf; ssf; sst]), ...
    1, 'LineWidth', 2 );
axis equal
hold off

clear ssf sst
clear numPoints sharpnessAngle edgeSensitvity neighborRadius

%% Generate Surface Mesh ==================================================
close all; clc;

surfOptions = struct();
surfOptions.reconstructionMethod = 'advancingFront';

if strcmpi(surfOptions.reconstructionMethod, 'advancingFront')
    surfOptions.reconstructionMethod = 0;
elseif strcmp(surfOptions.reconstructionMethod, 'scaleSpace')
    surfOptions.reconstructionMethod = 1;
else
    error('Invalid surface reconstruction method supplied');
end

surfOptions.removeOutliers = false;
surfOptions.numNeighborsRO = 24;
surfOptions.thresholdPercentRO = 5;

surfOptions.simplifyPointSet = true;
surfOptions.numNeighborsSpacing = 6;
surfOptions.spacingScaleSIPS = 10;

surfOptions.smoothPointSet = true;
surfOptions.numNeighborsSMPS = 50; % 24;

surfOptions.smoothIterations = 10;
surfOptions.maxFacetLength = Inf;

[F, V] = surface_reconstruction(PU, surfOptions);
[F, V, ~, ~, ~] = remove_unreferenced_vertices_from_mesh(F, V);

E = edges(triangulation(F, V));

fprintf('Number of boundary components = %d\n', ...
    numel(DiscreteRicciFlow.compute_boundaries(F)));

eulerChi = size(F,1) + size(V,1) - size(E,1);
fprintf('Euler characteristic = %d\n', eulerChi);

trisurf(triangulation(F,V));
axis equal;


