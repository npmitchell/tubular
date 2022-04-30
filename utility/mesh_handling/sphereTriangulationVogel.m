function sphereTri = sphereTriangulationVogel( N, R )
%SPHERETRIANGULATIONVOGEL Creates a triangulation of a sphere using a
% spherical Fibonacci grid
%   INPUT PARAMETERS:
%       - N:            The number of vertices in the triangulation will
%                       equal (2*N+1)
%       - R:            The radius of the sphere.  Default is R = 1
%
%   OUTPUT PARAMETERS:
%       - sphereTri:    A MATLAB-style triangulation of the sphere
%
%   by Dillon Cislo

if nargin < 2
    R = 1;
end

% A shifted list of vertex indices
NVec = (-N:N)';

% A constant related to the Golden Ratio
GRC = ( sqrt(5) - 1 ) / 2;

% The azimuthal coordinates of each vertex
phi = 2 .* pi .* GRC .* NVec;

% Trigonometric functions of the polar coordinates of each vertex
cosTheta = 2 .* NVec ./ ( 2 .* N + 1 );
sinTheta = sin(acos(cosTheta));

% The vertex coordinates
vertex = R .* [ sinTheta .* cos(phi), sinTheta .* sin(phi), cosTheta ];

% Construct a Delaunay tetrahedralization of the spherical volume
sphereTri = delaunayTriangulation( vertex );

% The output triangulation is the free boundary of the previous result
sphereTri = triangulation( freeBoundary( sphereTri ), vertex );

end

