function [ssv, radii, avgpts, cids] = srFromDVCurvesGivenCenterline(cseg_ss, cseg, curvesDV)
%srFromDVCurves Find pathlength match and radius for hoops in 3d 
%   Given the centerline pathlength cseg_ss, the centerline cseg, and the
%   array of 1d DV curves, compute the matching centerline pathlength for
%   each 
%
% 
% Returns
% -------
% ssv : N*M x 1 float
%   pathlength values of supplied centerline
% radii : N*M x 1 float
%   distance of each gridpt (on hoops) to pointmatch on supplied
%   centerline
% avgpts : N x 3 float
%   positions of the average of each hoop
% cids : N x 1 int
%   indices into supplied centerline that are closest on average to
%   points in each hoop

radii = zeros(size(curvesDV, 1), size(curvesDV, 2)) ;
cids = zeros(size(curvesDV, 1), 1) ;
avgpts = zeros(size(curvesDV, 1), 3) ;
% radii_from_mean = zeros(size(curvesDV, 1), size(curvesDV, 2)) ;

% Measure radius to nearest common centerline point for each bin
for jj = 1:size(curvesDV, 1)
    % Consider this hoop
    hoop = squeeze(curvesDV(jj, :, :)) ;

    % Compute radii for this bin
    % First find which centerline segment from which to compute radius
    rdist = zeros(length(cseg), 1) ;
    for kk = 1:length(cseg)
        rdist(kk) = mean(vecnorm(hoop - cseg(kk, :), 2, 2)) ;
    end
    [~, cid] = min(rdist) ;
    
    % Compute the radii
    radii(jj, :) = vecnorm(hoop - cseg(cid, :), 2, 2) ;
    cids(jj) = cid ;
    avgpts(jj, :) = mean(hoop) ; 
    % radii_from_mean(jj, :) = vecnorm(hoop - avgpts(jj, :), 2, 2) ;
end

% translate the centerline IDs into the pathlength values of the centerline
ssv = cseg_ss(cids) ;

end

