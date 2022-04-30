function vrot = rotate3D_rodrigues(v, k, phi)
%ROTATE3D_RODRIGUES Perform a 3D rotation of a set of points around a given
%axis of rotation using the Rodrigues formula
%
%   INPUT PARAMETERS:
%   
%       - v:    Nx3 input point coordinate list
%       - k:    3D vector specifying the axis of rotation
%       - phi:  The magnitude of the rotation in radians
%
%   OUTPUT PARAMETERS:
%
%       - vrot: Nx3 rotated output point coordinate list
%
% by Dillon Cislo 02/20/2020

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

validateattributes( v, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan', 'real'} );
validateattributes( k, {'numeric'}, ...
    {'vector', 'numel', 3, 'finite', 'nonnan', 'real'} );
validateattributes( phi, {'numeric'}, ...
    {'scalar', 'finite', 'nonnan', 'real'} );

% Turn k into a row vector if necessary
if ~(size(k,1) == 1), k = k.'; end

% Normalize k
k = k ./ sqrt(sum(k.^2));

% Match the dimensions of k to the dimensions of V
k = repmat(k, size(v,1), 1);

%--------------------------------------------------------------------------
% Perform Rotation
%--------------------------------------------------------------------------

crossKV = cross(k, v, 2);
dotKV = repmat(dot(k, v, 2), 1, 3);

vrot = cos(phi) .* v + sin(phi) .* crossKV + ...
    (1 - cos(phi)) .* dotKV .* k;


end
