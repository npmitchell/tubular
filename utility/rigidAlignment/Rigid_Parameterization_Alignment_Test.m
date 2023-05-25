%% Rigid Parameterization Alignment Test ==================================
clear; close all; clc;

% Generate a triangulation of the unit disk
diskTri = diskTriangulation(30);

tarF = diskTri.ConnectivityList;
tarV2D = diskTri.Points;

% Generate the 3D vertex positions
r = sqrt( sum( tarV2D.^2, 2 ) ); % Radial position of each vertex
phi = atan2( tarV2D(:,2), tarV2D(:,1) ); % Angular argument of each vertex

thetaCap = pi/2;
tarV3D = [ sin( thetaCap .* r ) .* cos(phi) ./ sin(thetaCap), ...
    sin( thetaCap .* r ) .* sin(phi) ./ sin(thetaCap), ...
    ( cos(thetaCap.*r)-cos(thetaCap) ) ./ sin(thetaCap) ];

movF = tarF;
movV3D = tarV3D;

% Generate a random rigid transformation of the domain of parameterization
phi = 2 * pi * rand(1);
u = 10 * rand(1) - 5;
v = 10 * rand(1) - 5;

% movV2D = tarV2D + repmat([u v], size(tarV2D,1), 1);
% movV2D = [ cos(phi) .* movV2D(:,1) - sin(phi) .* movV2D(:,2), ...
%     sin(phi) .* movV2D(:,1) + cos(phi) .* movV2D(:,2) ];

movV2D = tarV2D;
movV2D = [ cos(phi) .* movV2D(:,1) - sin(phi) .* movV2D(:,2), ...
    sin(phi) .* movV2D(:,1) + cos(phi) .* movV2D(:,2) ];
movV2D = movV2D + repmat([u v], size(tarV2D,1), 1);


[phiTest, uTest, vTest] = rigidParameterizationAlignment( ...
    tarF, tarV3D, tarV2D, movF, movV3D, movV2D );
phiTest = wrapTo2Pi(-phiTest);
       
        
        
