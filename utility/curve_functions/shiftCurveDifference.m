function [ssdiff] = shiftCurveDifference(t0, curv1, curv2)
% Given two curves, interpolate the first in reference to the second and
% find the sum squared differences betweeen the two curves as a function of
% a time offset. This can be used to minimize difference between two time
% offset curves. 
%
% parameters
% ----------
% t0 : float
%   time to shift curv1 relative to curv2's t axis to match the two curves
% curv1 : N x 1 float array
%   first column is dependent variable t, second is a value f(t)
% curv1 : N x 1 float array
%   first column is dependent variable t, second is a value f(t)
%
% Returns 
% -------
% ssdiff : float
%   the sum of the squared differences between curv1(t2) and curv2(t2),
%   where t2 are the timestamps given in curv2

% Unpack 
t1 = curv1(:, 1) ;
t2 = curv2(:, 1) ;
p1 = curv1(:, 2) ;
p2 = curv2(:, 2) ;
toff = t1 + t0 ;

% Interpolate curv1 at times t2
curv1intrp = interp1(toff, p1, t2, 'pchip');
ssdiff = sum((curv1intrp(t2) - p2).^2) ;

end

