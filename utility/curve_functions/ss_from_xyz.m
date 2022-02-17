function ss = ss_from_xyz(xyz)
% SS_FROM_XYZ(xyz) Obtain pathlength parameterization of curve xyz
% 
% Parameters
% ----------
% xyz : N x d array
%   d-Dimensional curve
%
% Returns
% -------
% ss
%
% NPMitchell 2019

% get distance increment
ds = vecnorm(diff(xyz), 2, 2) ;
% get pathlength at each skeleton point
ss = [0; cumsum(ds)] ;
return