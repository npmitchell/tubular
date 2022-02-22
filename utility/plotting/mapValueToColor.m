function colors = mapValueToColor(vals, vlims, cmap)
% MAPVALUETOCOLOR(vals, vlims, cmap)
%   map values in vals to colors in colormap cmap given bounds vlims
%
% Parameters
% ----------
% vals : N x 1 float/int array
%   the values to map to colors
% vlims : 1 x 2 float/int array
%   the limits to map to the first and last color in cmap
% cmap : Mx3 or Mx4 float/int array
%   the colormap to map into
% 
% Returns
% -------
% colors : N x 3 float/int array 
%   the colors in the colormap mapped to vals
%
% NPMitchell 2020

% Solve linear mapping coeffs (m,b)
% y = mx+b, where 1 = m * vlims(1) + b and length(cmap) = m * vlims(2) + b
% 1 - m*vlims(1) = length(cmap) - m * vlims(2) 
% m* (vlims(2) - vlims(1)) = length(cmap) - 1
% m = [length(cmap) - 1] / (vlims(2) - vlims(1))
% b = 1 - m* vlims(1)
slope = double(length(cmap) - 1) / (vlims(2) - vlims(1)) ;
intcp = 1 - vlims(1) * slope ;
colors = cmap(round(max(1, min(slope * vals + intcp, length(cmap)))), :) ;

