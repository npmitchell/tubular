function [totalVolume, totalArea] = meshVolumeArea(pts, tri)
% MESHVOLUMEAREA(pts, tri) compute area and volume of a triangulated mesh
%   Given a surface triangulation, compute the volume enclosed using
%   divergence theorem.
%   Assumption: Triangle nodes are ordered correctly, i.e.,computed normal is outwards
% 
% Parameters
% ----------
% pts: Nx3 or 3xN float
%   3D coordinates of the mesh vertices
% tri: Mx3 or 3xM integer 
%   face connectivity list indexing into pts 
% 
% Returns
% -------
% totalVolume: float
%   total volume enclosed
% totalArea : float 
%   total area of surface  
%
% Authors: K. Suresh (suresh@engr.wisc.edu) 
% cleanup, speedup, error checking, and documentation by NPMitchell 2019

% Ensure pts are the right shape
if size(pts, 2) == 3
    p = pts' ;
elseif size(pts, 1) == 3 
    p = pts ;
else
    error('pts must be 3D')
end

% Ensure tri is the right shape
if size(pts, 2) == 3
    t = tri' ;
elseif size(pts, 1) == 3 
    t = tri ;
else
    error('pts must be 3D')
end

% Compute the vectors d13 and d12
d13= [(p(1,t(2,:))-p(1,t(3,:))); (p(2,t(2,:))-p(2,t(3,:)));  (p(3,t(2,:))-p(3,t(3,:)))];
d12= [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); (p(3,t(1,:))-p(3,t(2,:)))];
cr = cross(d13,d12,1);  %cross-product (vectorized)
area = 0.5 * sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);  % Area of each triangle
totalArea = nansum(area);
crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
zMean = (p(3,t(1,:))+p(3,t(2,:))+p(3,t(3,:)))/3;
nz = -cr(3,:) ./ crNorm;  % z component of normal for each triangle
volume = area .* zMean .* nz;  % contribution of each triangle
totalVolume = nansum(volume);  %divergence theorem
