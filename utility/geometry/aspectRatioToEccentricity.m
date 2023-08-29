function [ecc, unce] = aspectRatioToEccentricity(ar, unc)
% Convert aspect ratio(s) to eccentricity.
%
% Uncertainty is given by error propagation:
% \sqrt{\frac{{\delta a}^2}{a^4 \left(a^2-1\right)}}
%
% Parameters
% ----------
% ar : Nx1 float
%   aspect ratios
% unc : Nx1 float
%   uncertainties in aspect ratio
% 
% Returns
% -------
% ecc : Nx1 float
%   eccentricities
% unce : Nx1 float
%   uncertainties in eccentricity
%
% NPmitchell 2021

ecc = sqrt(1 - 1 ./(ar).^2) ;
unce = sqrt(unc.^2 ./ (ar.^4 .* (ar.^2 - 1) ) ) ;