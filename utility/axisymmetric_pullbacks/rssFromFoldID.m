function [rssfold, rssfold_frac, rssmax] = rssFromFoldID(folds, timePoints, sphiBase)
% RSSFROMFOLDID(folds, timePoints, sphiBase) 
%   Evaluate the ringpath pathlengths at given fold locations
%
% Parameters
% ----------
% folds : N x 3 int array
%   index locations (of U coordinate list) of fold locations
% timePoints : N x 1 int array
%   timestamps associated with each mesh
% sphiBase : string
%   base string to load spcutMesh
%
% Returns
% -------
% rssfold : N x 1 float array
%   ringpath proper pathlength (ringpath_ss) evaluated at fold locations
% rssfold : N x 1 float array
%   ringpath proper pathlength (ringpath_ss) evaluated at fold locations as
%   fraction of total ringpath length
% rssmax : N x 1 float array
%   Maximum proper ringpath pathlength at each timepoint
%
%
% NPMitchell 2019

% pathlength preallocations
rssfold = zeros(length(timePoints), 3) ;
rssfold_frac = zeros(length(timePoints), 3) ;
rssmax = zeros(length(timePoints), 1) ;

% Consider each timepoint
for kk = 1:length(timePoints)
    % Translate to which timestamp
    t = timePoints(kk) ;
    load(sprintf(sphiBase, t), 'spcutMesh') ;
    
    % Extract centerline pathlength from DVhoop distance
    ss = spcutMesh.ringpath_ss ;
    maxss = max(ss) ;
    rssmax(kk, :) = maxss ;
    rssfold(kk, :) = ss(folds(kk, :)) ;
    rssfold_frac(kk, :) = rssfold(kk, :) / maxss ;
end