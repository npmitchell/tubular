function [totalArea] = meshArea(v, f)
% MESHAREA(pts, tri) compute area of a triangulated mesh
%   Given a surface triangulation compute area
% 
% Parameters
% ----------
% v: N x 3 float
%   3D coordinates of the mesh vertices
% f: M x 3 integer 
%   face connectivity list indexing into pts 
% 
% Returns
% -------
% totalVolume: float
%   total volume enclosed
%
% See also
% --------
% meshVolumeArea.m
% 
% NPMitchell 2019

% Compute the vectors d13 and d12
d13= [(v(f(:,2), 1)-v(f(:,3), 1)), (v(f(:,2), 2)-v(f(:,3), 2)), (v(f(:,2), 3) - v(f(:,3), 3))];
d12= [(v(f(:,1), 1)-v(f(:,2), 1)), (v(f(:,1), 2)-v(f(:,2), 2)), (v(f(:,1), 3) - v(f(:,2), 3))];
cr = cross(d13, d12, 2);  %cross-product (vectorized)
area = 0.5 * sqrt(cr(:,1).^2+cr(:,2).^2+cr(:,3).^2);  % Area of each triangle
totalArea = nansum(area);
