function rcurve = resampleCurve( curve, ds )
%RESAMPLECURVE Resample a D-dimensional curve to be uniformly spaced
%   Supply an N x D curve and a step size
%
%   Input Parameters:
%       - curve:        NxD float array
%       - ds:           float, step size for resampling
%
%   Output Parameters:
%       - rcurve:       MxD float array, the resampled curve
%
%   NPMitchell 2019

% get distance increment
ds_curve = vecnorm(diff(xyz), 2, 2) ;
total_length = sum(ds_curve) ; 
rcurve = curvspace(curve, total_length / ds) ;
return 