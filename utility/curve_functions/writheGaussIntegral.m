function [Wr, wr] = writheGaussIntegral(xyz, ss)
%WRITHEGAUSSINTEGRAL(xyz, s) Compute the writhe by directly evaluating the
%   Gauss integral. Also returns wr, the writhe density taken by
%   integrating only over ds', not ds.
%   
% Parameters
% ----------
% xyz : N x 3 float array
%   The 3D curve whose writhe we compute
% ss : N x 1 float array, optional
%   The cumulative pathlength from the endpoint of the curve xyz(1, :) to a
%   given curve point
% 
% Returns
% -------
% Wr : float 
%   The writhe of the curve as computed using the Gauss Integral
% wr : N x 1 float array 
%   The writhe "density" obtained by integrating only over ds', not over ds
% 
% NPMitchell 2019

% If the pathlength parameterization was not supplied, compute it
if nargin > 1
    if length(ss) < 1
        ss = ss_from_xyz(xyz) ;
    end
else
    ss = ss_from_xyz(xyz) ;
end

% Compute tangent
tangent = tangent_from_curve(ss, xyz(:, 1), xyz(:, 2), xyz(:, 3)) ;
ds = gradient(ss) ;
wr = zeros(length(ss), 1) ;
for jj=1:length(ss)
    oind = setdiff(1:length(ss), jj) ;
    rmr = xyz(jj, :) - xyz(oind, :) ;
    rmrmag = vecnorm(rmr')' ;
    txt = cross(tangent(jj,:) .* ones(length(oind), 3), tangent(oind, :));
    % Take row-wise inner product
    integrand = sum(sum(txt .* rmr, 2) ./ (rmrmag.^3 .* ones(size(txt))), 2) ;
    % Writhe per unit length is wr
    wr(jj) = sum(integrand.* ds(oind)) ;
    % wr(jj) = sum(integrand) ;
end
% Wr = nansum(wr .* ds) / (4 * pi) ;
Wr = nansum(wr .* ds) / (4 * pi) ;
wr = wr / (4 * pi) ;

% For some reason divide by 2pi
Wr = Wr / (2 * pi) ;
wr = wr / (2 * pi) ;

end

