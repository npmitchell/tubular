function [tracks, trackM] = trackNearest(xyzs,maxdisp,param)
% NPMitchell 2022
% documentation adapted from Crocker-Grier codebase.
% 
% NAME:
% track
% PURPOSE:
% Constructs n-dimensional trajectories from a scrambled list of
% particle coordinates determined at discrete times (e.g. in
% consecutive video frames).
% CATEGORY:
% Image Processing
% CALLING SEQUENCE:
% result = track( positionlist, maxdisp, param )
%  set all keywords in the space below
% INPUTS:
% positionlist: an array listing the scrambled coordinates and data 
%     of the different particles at different times, such that:
%  positionlist(0:d-1,*): contains the d coordinates and
%     data for all the particles, at the different times. must be positve
%  positionlist(d,*): contains the time t that the position 
%     was determined, must be integers (e.g. frame number.  These values must 
%               be monotonically increasing and uniformly gridded in time.
% maxdisp: an estimate of the maximum distance that a particle 
%     would move in a single time interval.(see Restrictions)
%  OPTIONAL INPUT:
%   param:  a structure containing a few tracking parameters that are
%       needed for many applications.  If param is not included in the
%       function call, then default values are used. 
%         param.dim: if the user would like to unscramble non-coordinate data
%             for the particles (e.g. apparent radius of gyration for
%             the particle images), then positionlist should
%             contain the position data in positionlist(0:param.dim-1,*)
%             and the extra data in positionlist(param.dim:d-1,*). It is then
%             necessary to set dim equal to the dimensionality of the
%             coordinate data to so that the track knows to ignore the
%             non-coordinate data in the construction of the 
%             trajectories. The default value is two.
%
% OUTPUTS:
% result:  a list containing the original data rows sorted 
%     into a series of trajectories.  To the original input 
%     data structure there is appended an additional column 
%     containing a unique 'id number' for each identified 
%     particle trajectory.  The result array is sorted so 
%     rows with corresponding id numbers are in contiguous 
%     blocks, with the time variable a monotonically
%     increasing function inside each block.  For example:
%     
%     For the input data structure (positionlist):
%         (x)      (y)      (t)
%     pos = 3.60000      5.00000      0.00000
%           15.1000      22.6000      0.00000
%           4.10000      5.50000      1.00000 
%           15.9000      20.7000      2.00000
%           6.20000      4.30000      2.00000
% ;
%     >> res = track(pos,5,mem=2)
% ;
%     track will return the result 'res'
%         (x)      (y)      (t)          (id)
%     res = 3.60000      5.00000      0.00000      0.00000
%           4.10000      5.50000      1.00000      0.00000
%           6.20000      4.30000      2.00000      0.00000
%           15.1000      22.6000      0.00000      1.00000
%           15.9000      20.7000      2.00000      1.00000
% 
%     NB: for t=1 in the example above, one particle temporarily
%     vanished.  As a result, the trajectory id=1 has one time
%     missing, i.e. particle loss can cause time gaps to occur 
%     in the corresponding trajectory list. In contrast:
% 
%     >> res = track(pos,5)
% 
%     track will return the result 'res'
%         (x)      (y)      (t)          (id)
%     res = 15.1000      22.6000      0.00000      0.00000
%                   3.60000      5.00000      0.00000      1.00000
%               4.10000      5.50000      1.00000      1.00000
%               6.20000      4.30000      2.00000      1.00000
%               15.9000      20.7000      2.00000      2.00000
% 
%     where the reappeared 'particle' will be labelled as new
%     rather than as a continuation of an old particle since
%     mem=0.  It is up to the user to decide what setting of 
%     'mem' will yeild the highest fidelity .
% 
% SIDE EFFECTS:
% Produces informational messages.  
% RESTRICTIONS:
% maxdisp should be set to a value somewhat less than the mean 
% spacing between the particles. 
% PROCEDURE:
% Given the positions for n particles at time t(i), and m possible
% new positions at time t(i+1), this function considers each old position 
% independently and latches on to one new position. If multiple old
% positions latch on to the same new position, each track takes that new
% branch.
% 
% Example Usage
% -------------
%  xyzs = [0,0,0; 0,1,0; 1,1,0; 1,0,0; ...
%         0,0,1; 0.1,0,1; 1,1,1; 10,1,1] ;
% [tracks, trackM] = trackNearest(xyzs,1) ;
% >> tracks = [0 0 0 1;
%              0 1 0 2;
%              1 1 0 3;
%              1 0 0 4;
%              
%
% 

dim = 2 ;
if nargin < 3
    param = struct() ;
end
if isfield(param, 'dim')
    dim = param.dim ;
end

timePoints = unique(xyzs(:, end)) ;
  
% all possible timestamps
tids = xyzs(:, end) ;
% Initial positions of cells to track
xy0 = xyzs(find(tids == timePoints(1)), 1:dim) ;

% cell ids for tracks
cids = 1:size(xy0, 1) ;
ncells = length(cids) ;

% initialize the tracks to make
trackM = zeros(length(timePoints), size(xy0, 1), size(xyzs, 2)) ;
tracks = zeros(length(timePoints) * size(xy0, 1), size(xyzs, 2)) ;

% First track is starting segmentation
info1 =  xyzs( find(tids == timePoints(1)) , :) ;
trackM(1, :, :) = info1 ;  
tracks = [reshape(trackM(1, :, :), [size(trackM, 2), size(trackM, 3)]), cids' ] ;
    
for tidx = 2:length(timePoints)
    
    % positions for next timestamp
    info1 =  xyzs( find(tids == timePoints(tidx)) , :) ;
    xy1 = info1( :, 1:dim ) ;
    
    % next ids
    [nids, dist] = dsearchn(xy1, xy0)  ;
    % nids are indices of xy1 that are closest to xy0
    
    % keep only the closest to redundant points
    [uniqueNIds, posUnique] = unique( nids, 'stable' );
    duplicate_indices = setdiff( 1:numel(nids), posUnique ) ;
    if ~isempty(duplicate_indices)
        for qq = 1:length(duplicate_indices)
            % find the closest match of the duplicates
            thisDuplicate = duplicate_indices(qq) ;
            p2consider = find(nids == nids(thisDuplicate)) ;
            [~, minId] = min(dist(p2consider)) ;
            p2rm = setdiff(p2consider, minId) ;
            dist(p2rm) = NaN ;
        end
    end
    
    % plot it
    % clf
    % plot(xy0(:, 1), xy0(:, 2), '.'); hold on;
    % plot(xy1(:, 1), xy1(:, 2), 'o')
    % for qq = 1:ncells
    %     if dist(qq) < maxdisp
    %         plot( [xy0(qq, 1), xy1(nids(qq), 1)]', [xy0(qq, 2), xy1(nids(qq), 2)]', '-')
    %     end
    % end
    
    % Update current positions
    xy0 = xy1(nids, :) ;
    xy0(dist > maxdisp | isnan(dist), :) = NaN ;
    
    % add to trackM
    trackM(tidx, :, :) = info1(nids, :) ;   
    trackM(tidx, dist > maxdisp | isnan(dist), 1:end-1) = NaN ;
    
    % cat to tracks
    t2add = [reshape(trackM(tidx, :, :), [size(trackM, 2), size(trackM, 3)]), cids' ] ; 
    
    tracks = cat(1, tracks, t2add) ;
end


 
