function [V, gradV, hessV] = stereographicSphere(U, k)
%STEREOGRAPHICSPHERE Stereographic parameterization of a spherical cap
%whose boundary is the unit disk
%
%   INPUT PARAMETERS:
%
%       - U:        #Vx2 2D vertex coordinate list
%
%       - k:        The curvature of the cap (k = 1/R)
%
%   OUTPUT PARAMETERS:
%
%       - V:        #Vx3 3D vertex coordinate list
%
%       - gradV:    #Vx3x2 surface gradients
%
%       - hessV:    #Vx3x3 surface Hessian values
%                   
%
%   by Dillon Cislo 03/21/2020

% Input Processing --------------------------------------------------------
if (nargin < 2), k = 1; end

validateattributes(U, {'numeric'}, {'2d', 'ncols', 2, 'finite', 'real'});
validateattributes(k, {'numeric'}, {'scalar', 'finite', 'real', '>=', 0});

% Construct Mapping -------------------------------------------------------

% Cartesian coordinates of the 2D vertices
x = U(:,1); y = U(:,2);

% Plane polar coordinates of the 2D vertices
rho = sqrt( x.^2 + y.^2 );

% The scale factor for the disk domain
r0 = ( 1 - sqrt(1-k.^2) ) ./ k.^2;

% The z-shift to match the boundary of the surface to the unit disk
z0 = ( 1 - k.^2 - sqrt(1-k.^2) ) ./ ( k .* sqrt(1-k.^2) - k );

% The stereographic surface coordinates
V = [ ...
    ( 2 .* r0 .* x ) ./ ( 1 + r0.^2 .* k.^2 .* rho.^2 ), ...
    ( 2 .* r0 .* y ) ./ ( 1 + r0.^2 .* k.^2 .* rho.^2 ), ...
    ( r0.^2 .* k.^2 .* rho.^2  - 1 ) ./ ...
    ( r0.^2 .* k.^3 .* rho.^2 + k ) + z0 ];

% Construct Gradients -----------------------------------------------------

gradV = [];
if (nargout > 1)
    
    % The denominator of all terms
    denom = ( 1 + r0.^2 .* k.^2 .* rho.^2 ).^2;
    
    % The partial derivatives dVdx
    DVx = [ 2 .* r0 + 2 .* k.^2 .* r0.^3 .* (y.^2 - x.^2), ...
        -4 .* k.^2 .* r0.^3 .* x .* y, ...
        4 .* k .* r0.^2 .* x ] ./ denom;
    
    % The partial derivatives dVdy
    DVy = [ -4 .* k.^2 .* r0.^3 .* x .* y, ...
        2 .* r0 + 2 .* k.^2 .* r0.^3 .* (x.^2 - y.^2), ...
        4 .* k .* r0.^2 .* y ] ./ denom;
    
    gradV = cat(3, DVx, DVy);
    
end

% Construct Hessians ------------------------------------------------------

hessV = [];
if (nargout > 2)
    
    % The denominator of all terms
    denom = ( 1 + r0.^2 .* k.^2 .* rho.^2 ).^3;
    
    % The partial derivatives d^2Vdx^2
    DVxx = [ 4 .* k.^2 .* r0.^3 .* x .* (-3 + k.^2 .* r0.^2 .* (x.^2 - 3 .* y.^2)), ...
        -4 .* y .* (k.^2 .* r0.^3 + k.^4 .* r0.^5 .* (y.^2 - 3 .* x.^2)), ...
        4 .* k.^2 .* r0.^2 .* (1 + k.^2 .* r0.^2 .* (y.^2 - 3 .* x.^2)) ] ./ denom;
    
    % The partial derivatives d^2Vdxdy
    DVxy = [ -4 .* y .* (k.^2 .* r0.^3 + k.^4 .* r0.^5 .* (y.^2 - 3 .* x.^2)), ...
        -4 .* x .* (k.^2 .* r0.^3 + k.^4 .* r0.^5 .* (x.^2 - 3 .* y.^2)), ...
        -16 .* k.^4 .* r0.^4 .* x .* y ] ./ denom;
    
    % The partial derivatives d^2Vdy^2
    DVyy = [ -4 .* x .* (k.^2 .* r0.^3 + k.^4 .* r0.^5 .* (x.^2 - 3 .* y.^2)), ...
        4 .* k.^2 .* r0.^3 .* y .* (-3 + k.^2 .* r0.^2 .* (y.^2 - 3 .* x.^2)), ...
        4 .* (k.^2 .* r0.^2 + k.^4 .* r0.^4 .* (x.^2 - 3 .* y.^2)) ] ./ denom;
    
    hessV = cat(3, DVxx, DVxy, DVyy);
    
end

end


