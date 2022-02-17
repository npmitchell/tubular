%%
clear; close all; clc;

%%
%--------------------------------------------------------------------------
% Lets construct two 'rectilinear' meshes understood to live in their
% respective frames at timepoint 0 and time point 1
%--------------------------------------------------------------------------

% Create basic grid points in 2D
umax = 2; % The length of the domain
vmax = 2; % The width of the domain
NL = 50; % The number of vertices along the length of the cylinder
NW = round( vmax * NL ./ umax); % The number of vertices along the width
[X, Y] = meshgrid( linspace(0, umax, NL), linspace(0, vmax, NW) );

% Triangulate grid points in 2D
TR = delaunayTriangulation([X(:), Y(:)]);

F = bfs_orient(TR.ConnectivityList);
X = TR.Points;

% The mesh at time 0
F0 = F;
V0 = X;

% The mesh at time 1
F1 = F;
V1 = X;

clear X Y TR F X

figure

subplot(1,2,1)
triplot(triangulation(F0, V0));
axis equal
title('rectilinear mesh, t=0')

subplot(1,2,2)
triplot(triangulation(F1,V1));
axis equal
title('rectilinear mesh, t=1')

%%

%--------------------------------------------------------------------------
% While these two meshes are both valid partitions of their respective
% domains, they do NOT represent the actual transformation between the
% domains.  Suppose instead we have a set of material points with known
% locations at time point 0 and time point 1
%--------------------------------------------------------------------------

numPoints = size(V0, 1);

% The locations of material points at time 0
P0 = [ min(V0(:,1)) + (max(V0(:,1))-min(V0(:,1))) .* rand(numPoints,1), ...
    min(V0(:,2)) + (max(V0(:,2))-min(V0(:,2))) .* rand(numPoints,1) ];

% The locations of material points at time 1
maxChange = 0.01;
P1 = P0 + maxChange .* rand(numPoints, 2);

hold on
scatter(P0(:,1), P0(:,2), 'filled', 'c');
scatter(P1(:,1), P1(:,2), [], 'm', 's');
quiver(P0(:,1), P0(:,2), (P1(:,1)-P0(:,1)), (P1(:,2)-P0(:,2)), ...
    0, 'c');
hold off

axis equal

%--------------------------------------------------------------------------
% It is the motion of THESE points that defines the Lagrangian update of
% the domains of parameterizations, NOT the motion (or lack thereof) of the
% rectilinear meshes)
%--------------------------------------------------------------------------

%% 

%--------------------------------------------------------------------------
% If we want to define a transformation between the frame of time point 1
% and into the frame of time point 1 that ends up recapitulating the mesh
% at time point 1, we must ask "Where were the vertices of the mesh at time
% point 1 back during time point 0?".  We can extract this information
% using the motion of the material points
%--------------------------------------------------------------------------

% The 'inverse map'
invMapX = scatteredInterpolant( P1(:,1), P1(:,2), P0(:,1), 'natural' );
invMapY = scatteredInterpolant( P1(:,1), P1(:,2), P0(:,2), 'natural' );

% The locations of the rectilinear vertices at time 0
% 'Deformed Vertices'
DV1 = [ invMapX(V1(:,1), V1(:,2)), invMapY(V1(:,1), V1(:,2)) ];

triplot(triangulation(F1, DV1));
axis equal

%%

%--------------------------------------------------------------------------
% This 'Deformed Mesh' is the one that will end up at the rectilinear grid
% of time point 1 when it rides along the Lagrangian flow from time point 0
% to time point 1.  It is the transformation between THIS mesh and the
% rectilinear grid at time point 1 that should be used to calculate the
% forward transformation Jacobian
%--------------------------------------------------------------------------

J01 = jacobian2Dto2DMesh( V1, DV1, F1 );



