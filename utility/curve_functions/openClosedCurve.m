function outcurv = openClosedCurve(curv, eps)
% openClosedCurve
%   very simple function that checks if curve is closed (endpt=startpt)
%   and removes endpoint if so. This is useful for optimization functions
%   on periodic data
%
% NPMitchell 2020

if nargin < 2 
    eps = 1e-6 ;
end

% Check if the endpt and startpt are within eps * mean spacing of each
% other
approx_ds = mean(vecnorm(diff(curv, 1), 2, 2)) ;
ept_is_spt = all(abs(curv(1, :) - curv(end, :)) < eps * approx_ds) ;

if ept_is_spt
    % input curv is indeed closed
    outcurv = curv(1:end-1, :) ;
else
    % input curv is not closed. Warn user
    disp('openClosedCurve: Input curv is not closed! Returning input curv')
    outcurv = curv;
end