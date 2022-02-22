function [h1, h2] = vectorFieldQuiverOnImage(im, xx, yy, vx, vy, vscale, ...
    options)
%VECTORFIELDQUIVERONIMAGE(im, xx, yy, vx, vy, vscale, options)
%   Plot a vector field (vx,vy) evaluated at grid[xx, yy] on an image im as
%   quiverplot (subsampled quiver)
%
% xx : N x 1 float array
%   x values of PIV grid evaluation points
% yy : M x 1 float array
%   y values of PIV grid evaluation points
% vx : N*M x 1 float array
%   velocity in x direction
% vy : N*M x 1 float array
%   velocity in y direction
% vscale : float
%   magnitude associated with maximum color/intensity in velocity image
% qopts : struct with fields
%   outfn : str
%       path to save image if given
%   label : str
%       colorbar label. Default is '$|v|$ [$\mu$m / min]' 
%   qsubsample : int
%       subsampling factor of the quiver field
%   qscale : float
%       overall scale of the quivers
%   outfn : str
%       output filename for figure as png 
%
% Returns
% -------
% h1 : handle for imshow
% h3 : handle for quiverplot
%
% NPMitchell 2020

% Default options
labelstr = '$|v|$ [$\mu$m / min]' ;
overlay_quiver = true ;
qsubsample = 10 ;
qscale = 10 ;

% Unpack options
if isfield(options, 'label')
    labelstr = options.label ;
end
if isfield(options, 'overlay_quiver')
    overlay_quiver = options.overlay_quiver ;
end
if isfield(options, 'qsubsample')
    qsubsample = options.qsubsample ;
end
if isfield(options, 'qscale') 
    qscale = options.qscale ;
end


% 
% vangle = reshape(mod(atan2(vy, -vx), 2* pi), gridsz) ;
% speed = reshape(vecnorm([v2dsm_ii(:, 1), v2dsm_ii(:, 2)], 2, 2), gridsz);
ww = length(xx) ;
hh = length(yy) ;
vangle = mod(atan2(vy, -vx), 2* pi) ;
speed = reshape(vecnorm([vx(:), vy(:)], 2, 2), [ww, hh]);

% Compute angle of the velocity vector
if ~all(size(vangle) == [ww, hh])
    vangle = reshape(vangle, [ww, hh]) ;
end

% Set up the figure
close all
fig = figure('units', 'normalized', ...
    'outerposition', [0 0 1 1], 'visible', 'off') ;
h1 = imshow(im) ;
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUIVER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vx = reshape(vx, [ww, hh]) ;
vy = reshape(vy, [ww, hh]) ;
QX = imresize(vx, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
QY = imresize(vy, [ww / qsubsample, hh / qsubsample], 'bicubic') ;
xq = 1:qsubsample:ww ;
yq = 1:qsubsample:hh ;
[xg, yg] = meshgrid(xx(xq), yy(yq)) ;

h2 = quiver(xg(:), yg(:), qscale * QX(:), qscale * QY(:), 0, 'k', 'LineWidth', 1.2) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phasemap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error('finish this function: make arrow to show scale')

if isfield(options, 'outfn')
    saveas(fig, options.outfn) ;   
    close all
end
