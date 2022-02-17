function [Tw, tw] = twist(curv1, curv2, ss1, ss2)
%TWIST Compute the twist of one curve about another
%   Computes twist of curv1 about curv2 using TxT.r/|r|^3 integral.
% 
% Parameters
% ----------
% curv1 : N x 3 float array
%   xyz coordinates of curve 1, in order
% curv2 : M x 3 float array
%   xyz coordinates of curve 2, in order
% ss1 :  N x 3 float array (optional)
%   pathlength parameterization of curve 1
% ss2 :  M x 3 float array (optional)
%   pathlength parameterization of curve 2
% 
% Returns
% -------
% Tw : float
%   total twist of curv1 about curv2
% tw : N x 1 float array
%   twist density of curv1 about curv2 (not multiplied by ds1)
%
% NPMitchell 2019


% If the pathlength parameterization was not supplied, compute it
if nargin > 2
    if length(ss1) < 1
        ss1 = ss_from_xyz(curv1) ;
    end
    if nargin > 3
        if length(ss2) < 1
            ss2 = ss_from_xyz(curv2) ;
        end
    else
        ss2 = ss_from_xyz(curv2) ;
    end
else
    ss1 = ss_from_xyz(curv1) ;
    ss2 = ss_from_xyz(curv2) ;
end

% Compute tangent
tang1 = tangent_from_curve(ss1, curv1(:, 1), curv1(:, 2), curv1(:, 3)) ;
tang2 = tangent_from_curve(ss2, curv2(:, 1), curv2(:, 2), curv2(:, 3)) ;
ds1 = gradient(ss1) ;
ds2 = gradient(ss2) ;
tw = zeros(length(ss1), 1) ;
for jj=1:length(ss1)
    rmr = curv1(jj, :) - curv2 ;
    rmrmag = vecnorm(rmr')' ;
    txt = cross(tang1(jj,:) .* ones(length(ss2), 3), tang2);
    % Take row-wise inner product
    integrand = sum(sum(txt .* rmr, 2) ./ (rmrmag.^3 .* ones(size(txt))), 2) ;
    % Writhe per unit length is wr
    tw(jj) = sum(integrand .* ds2) ;
end
Tw = nansum(tw .* ds1) / (4 * pi) ;
tw = tw / (4 * pi) ;

end

