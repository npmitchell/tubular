%% Point Set Surface Reconstruction Example ===============================
%
%   This is a script to test the functionality of a pipeline for generating
%   mesh triangulations of noise sparse point sets
%

addpath(genpath('/home/dillon/Documents/MATLAB/GrowSurface'));
addpath(genpath('/home/dillon/Documents/MATLAB/CGAL_Code'));
addpath(genpath('/home/dillon/Documents/MATLAB/RicciFlow_MATLAB'));
addpath(genpath('/home/dillon/Documents/MATLAB/gptoolbox'));

%% Load Inital Point Set from File ========================================
clear; close all; clc;

P = readOBJ('pointCloud_T295.obj');
% P = readOBJ('pointCloud_T400.obj');

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

numPoints = 3000; % The number of output points
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
surfOptions.reconstructionMethod = 'scaleSpace';

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

surfOptions.simplifyPointSet = false;
surfOptions.numNeighborsSpacing = 6;
surfOptions.spacingScaleSIPS = 2;

surfOptions.smoothPointSet = true;
surfOptions.numNeighborsSMPS = 24;

surfOptions.smoothIterations = 4;
surfOptions.maxFacetLength = 40;

[F, V] = surface_reconstruction(PU, surfOptions);
[F, V, ~, ~] = remove_isolated_mesh_components(F, V);

% F = fill_holes(V, F);

tar_length = 15;
num_iter = 5;
[F, V, ~, ~] = isotropic_remeshing(F, V, tar_length, num_iter);

E = edges(triangulation(F, V));

fprintf('Number of boundary components = %d\n', ...
    numel(DiscreteRicciFlow.compute_boundaries(F)));

eulerChi = size(F,1) + size(V,1) - size(E,1);
fprintf('Euler characteristic = %d\n', eulerChi);

trisurf(triangulation(F,V));
axis equal;

clear tar_length num_iter

%% Generate Figure -------------------------------------------------
close all; clc;

[x, y, z] = sphere(10);
cR = 15;
cColors = repmat([1 0 1], size(P0, 1), 1);

patch('Faces', F, 'Vertices', V, 'FaceVertexCData', 0.8 * ones(size(F)), ...
    'FaceColor', 'flat', 'FaceLighting', 'gouraud', 'EdgeColor', 'none');

hold on
for i = 1:size(P0,1)
    
    surf(cR*x+P0(i,1), cR*y+P0(i,2), cR*z+P0(i,3), ...
        'FaceColor', cColors(i,:), 'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud');
    
end
hold off

camlight
xlabel('x'); ylabel('y'); zlabel('z');
axis equal

view([-16 73]);
xlim([590.0 1233.3]);
ylim([501.1 1318.3]);
zlim([621.4 944.5]);


