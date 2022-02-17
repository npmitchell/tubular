function [ssx, xp, yp, zp, coeffs] = ...
    smooth_curve_via_fit_3d(ss, skelrs, polyorder, framelen, fitorder)
%SMOOTH_CURVE_VIA_FIT Smooth a curve by fitting it to a ploynomial. 
%   Assumes reasonably uniform spacing in s for the fit (so equal weight 
%   for savgol filter).
% 
% Parameters
% ----------
% ss : the pathlength parameterization of the input curve
% skelrs : N x 3 input curve coordinates
% polyorder : int
% framelen : odd int
%
% Returns
% -------
% ssx : pathlength parameterization
% xp, yp, zp : N x 1 float arrays of the smoothed curve coords in 3d
% 
% NPMitchell 2019 

if nargin < 5
    fitorder = 10 ;
end

xskel = skelrs(:, 1) ;
yskel = skelrs(:, 2) ;
zskel = skelrs(:, 3) ;

% Note we ignore the variations in ds to do this fit
scx = savgol(xskel, polyorder, framelen)' ;
scy = savgol(yskel, polyorder, framelen)' ;
scz = savgol(zskel, polyorder, framelen)' ;
% sc = [scx, scy, scz] ;

% Fit the smoothed curve
[xcoeffs, ~, Mux] = polyfit(ss, scx, fitorder) ;
[ycoeffs, ~, Muy] = polyfit(ss, scy, fitorder) ;
[zcoeffs, ~, Muz] = polyfit(ss, scz, fitorder) ;
ssx = linspace(min(ss), max(ss), 100) ;
xp = polyval(xcoeffs, (ssx - Mux(1)) / Mux(2));
yp = polyval(ycoeffs, (ssx - Muy(1)) / Muy(2));
zp = polyval(zcoeffs, (ssx - Muz(1)) / Muz(2));

% If requested, return coefficients in struct
if nargout > 4
    coeffs.xcoeffs = xcoeffs ;
    coeffs.ycoeffs = ycoeffs ;
    coeffs.zcoeffs = zcoeffs ;
    coeffs.Mux = Mux ;
    coeffs.Muy = Muy ;
    coeffs.Muz = Muz ;
end

end

