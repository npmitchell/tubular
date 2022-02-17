function [cutMesh, adIDx, pdIDx, cutP, Tw] = generateCutMeshFixedTwist(mesh, ...
    adIDx, pdIDx, centerline, nsegs4path, prevTw, previousP, varargin)
% GENERATECUTMESH(mesh, adIDx, pdIDx, centerline, nsegs4path)
%   Generating cutMesh from a mesh and the start/endpoints a/pdIDx.
%   Here, ensure that the topology of the orbifold mapping of the
%   cylinder is consistent with the previous timepoint by using twist of
%   the cutpath around the centerline as a proxy for the coiling of the
%   cutpath. Only accept the cutpath that is within MaxTwChange of the
%   previous twist value prevTw (ie the twist of the previous cutpath 
%   around the previous centerline). If the twist does not match, then
%   perturb the previous cutpath by a random kick, find the closest
%   piecewise geodesic path on the current mesh and try that as the
%   cutpath. Iterate until a path is found or maximum number of iterations
%   is exhausted. Ramp up the randomization via jitter_growth_rate and ramp
%   up the number of previous cutpath's sampled points to point match to 
%   current mesh for piecewise geodesic definition via nsegs_growth_rate.
% 
% Parameters
% -----------
% mesh
% adIDx 
% pdIDx
% centerline : N x 3 float
%   approx centerline of the mesh in original pixel coordinates
% nsegs4path : int
%   how many segments to divide the path into initially to compute a nearby
%   cut
% prevTw : float
%   The twist value or winding number of previous mesh
% outcutfn : str
%   The output Base filename for the cutPath (evaluatable string with
%   sprintf)
% cylinderMeshCleanBase : str
%   The output Base filename for the cylinderMesh 
% t : int
%   time index of current mesh, whos path is to be compared to previous
%
% Optional Input Arguments
% ------------------------
% MaxTwChange : float
%   Maximum twist change allowed in the path
% MaxJitter : float 
%   Distance in mesh coordinate space units for how much to perturb the
%   guessed cut path if twist is too far from previous value
% prevCenterLine : Qx3 float array, optional
%   the previous timepoints centerline, the twist about which we will match
% centerlineIsErratic : bool
%   if true, compare the cutPath not only to the current centerline for
%   twist computation, but also to the previous centerline. If twist around
%   either is consistent with previous twist, then accept the cut path.
%
% Returns
% -------
% prevTw : 
% previousP
% compute_pullback = success
% prevTw = Tw : the current Twist value, to be stored for next time
%
% NPMitchell 2019 


% success = true ;  % mark as successful unless otherwise caught

max_Tw_change = 0.25 ;
max_jitter = 10 ;
max_nsegs4path = 10 ;
centerlineIsErratic = false ;
prevcntrline = [] ;
jitter_growth_rate = 10 ;  %linear growth rate
nsegs_growth_rate = 1 ; % linear growth rate
for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
    end
    if isa(varargin{i},'logical')
        continue;
    end
    
    if ~isempty(regexp(varargin{i},'^[Mm]ax[Tt]w[Cc]hange','match'))
        max_Tw_change = varargin{i+1} ;
    end
    if ~isempty(regexp(varargin{i},'^[Mm]ax[Jj]itter','match'))
        max_jitter = varargin{i+1} ;
    end
    if ~isempty(regexp(varargin{i},'^[Pp]rev[Cc]enter[Ll]ine','match'))
        prevcntrline = varargin{i+1} ;
    elseif ~isempty(regexp(varargin{i},'^[Pp]rev[Cc]ntr[Ll]ine','match'))
        prevcntrline = varargin{i+1} ;
    end
    if ~isempty(regexp(varargin{i}, '^[Cc]enterline[Ii]s[Ee]rratic', ...
            'match'))
        centerlineIsErratic = varargin{i+1} ;
    end
end

% First try geodesic approach
cutOptions.method = 'fastest' ;
disp(['Cutting mesh using method ' cutOptions.method])

% try fastest approach first
cutOptions.bbWeight = 100 ;
cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );
cutP = cutMesh.pathPairs(:, 1) ;
adIDx = cutP(1) ;
pdIDx = cutP(end) ;

% Assert that the cutMesh is a topological disk
eIDx = topologicalStructureTools(triangulation(cutMesh.f, cutMesh.v)) ;
eulerChar = size(cutMesh.v, 1) - size(eIDx, 1) + size(cutMesh.f, 1) ;
assert(eulerChar == 1)

% If twist about centerline has jumped, find nearest piecewise geodesic
Tw = twist(mesh.v(cutP, :), centerline) ;
disp(['Found Twist = ' num2str(Tw) '; previous twist = ' num2str(prevTw)])

% Initialize some book-keeping variables
trykk = 0 ;
jitter_amp = 0 ;
nsegs4path_kk = nsegs4path ;

% Determine if twist has changed too much, indicating a change in topology
% of the cutPath for the orbifold mapping
twist_changed = abs(Tw - prevTw) > max_Tw_change ;

% Allow leniency to above conditional 
if centerlineIsErratic
    if isempty(prevcntrline)
        error('Too allow erratic centerline, please pass previous centerline to generateCutMeshFixedTwist()')
    end
    Tw_wprevcline = twist(mesh.v(cutP, :), prevcntrline) ;
    twist_changed = twist_changed && ...
        abs(Tw_wprevcline - prevTw) > max_Tw_change ;
end

% Initialize the indices of the previous cutPath to use for nearest curve
inds2use = [1, length(previousP)] ;
    
while twist_changed
    disp(['Twist out of range. Forcing nearest curve: Tw=' num2str(Tw) ', Twprev=' num2str(prevTw)])
    % If this is the first timepoint we are considering
    % or if we have to iteratively refine the target 
    % curve, load previous cutP.

%     %%%%%%%%%%%%%%%%%%%
%     % OPTION 1
%     % subsample the previous path to get previousP_kk
%     % Note: avoid oversampling by taking min, with
%     % fudge factor to avoid perfect sampling (which
%     % is too dense also).
%     disp(['nsegs4path = ' num2str(min(max_nsegs4path, nsegs4path_kk))])
%     pstep = round(length(previousP) / min(max_nsegs4path, nsegs4path_kk)) ;
%     disp(['pstep = ' num2str(pstep)])
%     disp(['jitter_amp = ', num2str(min(jitter_amp, max_jitter))])
%     previousP_kk = previousP(1:pstep:end, :);
%     jitter = min(jitter_amp, max_jitter) * (rand(size(previousP_kk)) - 0.5);
%     previousP_kk = previousP_kk + jitter ;
%     %%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%
    % OPTION 2: Adjust points of maximum divergence
    % Use previous path to find point of maximum difference from current.
    % Use this to build previousP_kk
    % Note: avoid oversampling by taking min, with
    % fudge factor to avoid perfect sampling (which
    % is too dense also).
    disp(['nsegs4path = ' num2str(min(max_nsegs4path, nsegs4path_kk))])
    % Translate number of segments to cut previous path into number of
    % points along that path to identify. 
    % By default, inds2find will be computed to be one. 
    % This makes sure that each iteration, we add only one landmark each
    % iteration, so that we don't identify a bunch of points clustered
    % near one errant segment of the path.
    inds2find = max(1, min(max_nsegs4path, nsegs4path_kk-1)) - (length(inds2use)-2) ;
    % point-match previous curve with current
    dists = zeros(size(previousP,1), 1) ;
    thisP = mesh.v(cutP, :) ;
    for i = 1:size(previousP,1)
        dists(i) = min(sqrt((thisP(:,1) - previousP(i,1)).^2 +...
            (thisP(:,2) - previousP(i,2)).^2 + ...
            (thisP(:,3) - previousP(i,3)).^2));
    end
    % Find maximum deviation indices and use these to define the sampling
    % of the previous path
    [~, maxdevID] = maxk(dists, inds2find) ;    
    inds2use = sort([inds2use, maxdevID']) ;
    disp(['jitter_amp = ', num2str(min(jitter_amp, max_jitter))])
    previousP_kk = previousP(inds2use, :);
    jitter = min(jitter_amp, max_jitter) * (rand(size(previousP_kk)) - 0.5);
    previousP_kk = previousP_kk + jitter ;
    %%%%%%%%%%%%%%%%%%%%%

    cutOptions.method = 'nearest' ;
    cutOptions.path = previousP_kk;
    disp(['Cutting mesh using method ' cutOptions.method])
    try
        figure(1); clf
        set(gcf, 'visible', 'on')
        trisurf(cutMesh.f, cutMesh.v(:, 1), cutMesh.v(:, 2), ...
            cutMesh.v(:, 3), cutMesh.v(:, 3), ...
            'EdgeColor', 'none', 'FaceAlpha', 0.1)
        hold on;
        plot3(cutMesh.v(cutP, 1), cutMesh.v(cutP, 2), cutMesh.v(cutP, 3), 'k-')
        plot3(previousP_kk(:, 1), previousP_kk(:, 2), previousP_kk(:, 3), 'o')
        plot3(centerline(:, 1), centerline(:, 2), centerline(:, 3), '-')
        plot3(previousP(:, 1), previousP(:, 2), previousP(:, 3), 'b--')
        if ~isempty(prevcntrline)
            plot3(prevcntrline(:, 1), prevcntrline(:, 2), prevcntrline(:, 3), '--')
        end

        title(['update cutP sampling = ' num2str(length(previousP_kk))])
        hold off 
        axis equal
        pause(1e-4)
            
        cutMesh = cylinderCutMesh( mesh.f, mesh.v, mesh.vn, adIDx, pdIDx, cutOptions );  
        cutP = cutMesh.pathPairs(:,1) ;
        
        % Here we compute both twist of path around current centerline and
        % the previous centerline, to allow for the possibility that the
        % centerline jitters a bit too much
        Tw = twist(mesh.v(cutP, :), centerline) ;
        Tw_wprevcline = twist(mesh.v(cutP, :), prevcntrline) ;
        
        % If the twist is too far from the previous twist, resample the
        % curve that cuts the mesh
        disp(['Found new Twist = ', num2str(Tw), ...
            ', Tw with prev cline=', num2str(Tw_wprevcline)])
        twist_changed = abs(Tw - prevTw) > max_Tw_change ;
         
        % Note: could add leniency to above conditional via: 
        if centerlineIsErratic
            twist_changed = twist_changed && ...
                abs(Tw_wprevcline - prevTw) > max_Tw_change ;
        end
        
        % If still no good, increase number of segments for guided path and
        % increase jitter applied to path landmarks
        if twist_changed
            nsegs4path_kk = round(nsegs4path_kk + nsegs_growth_rate) ;
            jitter_amp = jitter_amp + jitter_growth_rate ;
            trykk = trykk + 1 ;

            % show progress
            figure(1); clf
            set(gcf, 'visible', 'on')
            trisurf(cutMesh.f, cutMesh.v(:, 1), cutMesh.v(:, 2), ...
                cutMesh.v(:, 3), cutMesh.v(:, 3), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.1)
            hold on;
            plot3(cutMesh.v(cutP, 1), cutMesh.v(cutP, 2), cutMesh.v(cutP, 3), 'k-')
            plot3(previousP_kk(:, 1), previousP_kk(:, 2), previousP_kk(:, 3), 'o')
            plot3(centerline(:, 1), centerline(:, 2), centerline(:, 3), '-')
            plot3(previousP(:, 1), previousP(:, 2), previousP(:, 3), 'b--')
            if ~isempty(prevcntrline)
                plot3(prevcntrline(:, 1), prevcntrline(:, 2), prevcntrline(:, 3), '--')
            end
            
            title(['sampling = ' num2str(length(previousP_kk))])
            hold off 
            axis equal
            pause(1e-4)
        end
    catch
        disp('Could not generate cutMesh, likely not a topological disk')
        
        eulerChar = eulerCharacteristic(cutMesh) ;
        disp(['Euler characteristic = ', num2str(eulerChar)])
        pause(0.0001)
        close all
        nsegs4path_kk = round(nsegs4path_kk + nsegs_growth_rate) ;
        jitter_amp = jitter_amp + jitter_growth_rate ;
        trykk = trykk + 1 ;
    end

    % Give up on perturbing the curve after 10 tries, instead construct
    % the solution by projecting onto an annulus
    % if trykk > 10
    %     % Glue the previous mesh back together again
    %     % prevcmesh = closeRectilinearCylMesh(prevmesh) ;
    %     % Glue the current mesh back together again
    %     % cmesh = closeRectilinearCylMesh(mesh) ;
    %     wN = annularPathWindingNumber(cmesh.f, cmesh.v, cutMesh.pathPairs) ;
    % end
end
disp('Twist within range. Continuing...')

% Redefine the endpoints
adIDx = cutP(1) ;
pdIDx = cutP(end) ;
