function X3D = stereographicParaboloid(X2D, h)
%STEREOGRAPHICPARABOLOID Construct a parabaloid of revolution using a
%stereographic parameterization of the unit disk
%
%   INTPUT PARAMETERS:
%
%       - X2D:      #Vx2 list of 2D input point coordinates
%
%       - h:        The height of the paraboloid
%
%   OUTPUT PARAMETERS:
%
%       - X3D:      #Vx3 list of 3D output point coordinates

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

if (nargin < 1), error('Please supply input point coordinates!'); end
if ( (nargin < 2) || isempty(h) ), h = 1; end

validateattributes(X2D, {'numeric'}, ...
    {'2d', 'real', 'finite', 'ncols', 2});
validateattributes(h, {'numeric'}, ...
    {'scalar', '>=', 0 'finite', 'real'});

%--------------------------------------------------------------------------
% Construct Surface Mapping
%--------------------------------------------------------------------------

% The polar radial coordinate of the input points
r = sqrt( sum( X2D.^2, 2 ) );

% % The polar angle of the input points
phi = atan2(X2D(:,2), X2D(:,1));

% % The height of each output point
% Z = ( h ./ (2 .* r.^2) ) .* ...
%     ( sqrt( 1 + 8 .* r.^2 ) - 1 - 2 .* r.^2 );
% 
% % Deal with NaN at r = 0
% Z(r <= 10*eps) = h;
% 
% % The polar radial coordinate of the output points
% rho =  ( sqrt( 1 + 8 .* r.^2 ) - 1 ) ./ ( 2 .* r );
% 
% X3D = [ rho .* cos(phi), rho .* sin(phi), Z ];

C = 1 + sqrt(1 + 4 * h^2);

rho = C .* r ./ ( 1 + sqrt( 1 + 4 .* h.^2 .* r.^2 ) );

% rho = sinh(r) ./ (2 * h);

U2D = rho .* exp( 1i .* phi );

U2D = [real(U2D), imag(U2D)];

Z = -h .* (rho.^2-1);

X3D = [ U2D, Z ];

end

