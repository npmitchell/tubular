function colorsM = vectorHeatMap(mag, theta, options)
% polarHeatMap(mag, theta, options)
%   Transform vector field into colored grid
% 
% 
% theta : NxM float 
%   angle of vector field in radians
%
%
% NPMitchell

if nargin < 3
    options = struct() ;
end
if isfield(options, 'xx')
    xx = options.xx ;
else
    xx = 1:size(mag, 2) ;
end
if isfield(options, 'yy')
    yy = options.yy ;
else
    yy = 1:size(mag, 1) ;
end
if isfield(options, 'climit')
    climit = options.climit ;
else
    climit = max(mag(:)) ;
end
if isfield(options, 'colormap')
    cmapPolar = options.colormap ;
else
    cmapPolar = pm256 ;
end

% Map intensity from dev and color from the theta
indx = max(1, round(mod(theta(:), 2*pi)*size(cmapPolar, 1)/(2 * pi))) ;
colors = cmapPolar(indx, :) ;
devKclipped = min(mag(:) / climit, 1) ;
colorsM = devKclipped(:) .* colors ;
colorsM = reshape(colorsM, [size(mag, 1), size(mag, 2), 3]) ;
imagesc(xx, yy, colorsM)