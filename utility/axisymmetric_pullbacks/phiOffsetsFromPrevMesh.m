function [phi0s, residuals] = phiOffsetsFromPrevMesh(TF, TV2D, TV3D, uspace, ...
    vspace, prev3d_sphi, lowerbound, upperbound, phiOptions)
%PHIOFFSETSFROMPREVMESH(TF, TV2D, TV3Drs, nU, vpsace, prev3d_sphi) 
%   Find the offset in phi (the y dimension of the 2d pullback) that
%   minimizes the difference in 3D of the positions of each DV hoop from
%   that of the previous timepoint.
%   NEW COORDINATES: phi = v - phi0
%
% Parameters
% ----------
% TF : nU*nV x 3 int array
%   The mesh connectivity list, indexing into the vertex arrays TV2D and
%   TV3Drs
% TV2D : nU*nV x 2 float array
%   The mesh vertex locations in 2d
% TV3Drs : nU*nV x 3 float array
%   The mesh vertex locations in 3d
% uspace : nU float array
%   The values of u for each line of constant v in pullback space
% vspace : nV float array OR nU x nV float array as grid
%   If nV x 1 float array, the values of v for each line of constant u in 
%   pullback space, otherwise the values for the whole grid
% prev3d_sphi : nU x nV x 3 float array
%   The 3D coordinates of the embedding for the reference timepoint
%   (previous timepoint, for ex) at the 2D locations given by uspace and 
%   vspace. Note that uspace is not used explicitly, only nU is used to 
%   extract the strips over which we iterate, minimizing for phi0 for each
%   strip.
% vargin : variable input arguments
%   preview : bool
%   optimization : struct, default=omtimset()
%       optimization options, as output of optimset()
%   resample_hoops : bool, default=true
% 
% Returns
% -------
% phi0s : nV x 1 float array
%   the offsets in V dimension that minimize distances of each DV hoop from
%   analogous hoop in previous xyz (from previous timepoint, for ex)
%
%
% NPMitchell 2019

%% Interpret Options
if nargin > 8
    %%%%%%%%%%%%%%%%%%%%
    % Unpack Options
    %%%%%%%%%%%%%%%%%%%%
    
    % boolean for visualization
    if isfield(phiOptions, 'preview')
        preview = phiOptions.preview ;
        if preview
            fig = figure('visible', 'on') ;
        end
    else
        preview = false ;
    end
    
    % options for optimization problem
    if isfield(phiOptions, 'optimization')
        options = phiOptions.optimization ;
    else
        % options = optimset('PlotFcns','optimplotfval','TolX',1e-7);
        options = optimset() ; 
    end
    
    % resample_hoops: boolean to resample each DV hoop uniformly for the 
    % geometric optimization (rotation)
    if isfield(phiOptions, 'resample_hoops')
        resample_hoops = phiOptions.resample_hoops ;
    else
        resample_hoops = true ; 
    end
    
    % 
    if isfield(phiOptions, 'avgpts') && isfield(phiOptions, 'prev_avgpts')
        avgpts = phiOptions.avgpts ;
        prev_avgpts = phiOptions.prev_avgpts ;
        optimize_hoop_distance = true ;
        % Get pathlengths along piecewise curve
        prev_ss = phiOptions.prev_avgpts_ss ;
        prev_ss = prev_ss / max(prev_ss) ;
    elseif isfield(phiOptions, 'avgpts') || isfield(phiOptions, 'prev_avgpts')
        error(['Hoop-averaged centerlines were passed to phiOptions, ', ...
            'but only for one mesh!'])
    else
        optimize_hoop_distance = false ;
        error('not optimizing hoop distance -- pass avgpts etc in phiOptions')
    end
else
    preview = true ;
    % options = optimset('PlotFcns','optimplotfval', 'TolX',1e-7); 
    options = optimset() ; % 'TolX',1e-7); 
    resample_hoops = true ;
    optimize_hoop_distance = false ;
end

%% Consider each value of u in turn
% Fit for phi0 such that v = phi - phi0
nU = length(uspace) ;
if any(size(vspace) == 1)
    % we find that vspace has been passed as a 1d array, as in linspace
    nV = length(vspace) ;
    input_v_is_function_of_u = false ;
else
    % we find that vspace has been passed as a grid, as in meshgrid
    nV = size(vspace, 1) ;
    input_v_is_function_of_u = true ;
end

%% Decide if we compute residuals based on output args
if nargout > 1
    compute_residual = true ;
    residuals = zeros(nU, nV - 1) ;
else
    compute_residual = false ;
end

%% Optimize phi(u) for each u value in 1:nU discretization 
phi0s = zeros(nU, 1) ;
prog = repmat('.', [1 floor(nU/10)]) ;
for qq = 1:nU
    tic 
    
    % Define vqq, the input v values for this u=const strip of the pullback
    if input_v_is_function_of_u
        vqq = vspace(qq, :)' ;
    else
        % The V values here are linspace as given, typically (0...1)
        vqq = vspace ;
    end
    
    % The previous 3d embedding values are stored in prev3d_sphi
    % Choose which previous mesh's hoops against which to compare the 
    % current hoop (qq).
    if optimize_hoop_distance 
        % two matching hoops are how far away?
        % NOTE: get the first hoop behind, and the first hoop ahead
        [~, ~, ss] = distance2curve(prev_avgpts, avgpts(qq, :), 'linear') ;
        ahead = find(ss < prev_ss, 1) ;
        behind = find(ss > prev_ss, 1, 'last') ; 
        inds = [ahead, behind] ;
        if length(inds) > 1
            % There are two nearby DV rings: one ahead, one behind. 
            % This should be the case every time except at endpts
            dists = vecnorm(avgpts(qq, :) - prev_avgpts(inds, :), 2, 2) ; 
            assert(any(dists)) 
            % Handle if one of the matches is exact (this should never happen)
            if dists(1) == 0
                prev3dvals = squeeze(prev3d_sphi(inds(1), :, :)) ;
            else
                % The usual case. There are two hoops nearby. Use both,
                % weighted by distance.
                weights = 1 - dists / (dists(1) + dists(2)) ;
                weights = weights / sum(weights) ;
                try
                    assert(all(weights < 1))
                    assert((sum(weights) == 1))
                catch
                    disp('weights not right!')
                end
                prev3dvals = weights(1) * squeeze(prev3d_sphi(inds(1), :, :)) + ...
                             weights(2) * squeeze(prev3d_sphi(inds(2), :, :)) ;
            end
        elseif ~isempty(inds)
            prev3dvals = squeeze(prev3d_sphi(inds, :, :)) ;
        else
            error('No matching hoops. Handle this case here')
        end
    else
        % Here we simply line up the hoop with the previous mesh's hoop 
        % matching the same index qq.
        prev3dvals = squeeze(prev3d_sphi(qq, :, :)) ;
    end
    
    % Check it
    % close all
    % disp('previous values are')
    % plot3(prev3dvals(:, 1), prev3dvals(:, 2), prev3dvals(:, 3), '.')
    % hold on;
    % tmpx = prev3d_sphi(:, :, 1) ;
    % tmpy = prev3d_sphi(:, :, 2) ;
    % tmpz = prev3d_sphi(:, :, 3) ;
    % scatter3(tmpx(:), tmpy(:), tmpz(:), 2, 'MarkerFaceAlpha', 0.4,...
    %     'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'c')
    % temp_uv = [uspace(qq) * ones(nV, 1), mod(vqq, 1)] ;
    % temp3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, temp_uv) ;
    % plot3(temp3d(:, 1), temp3d(:, 2), temp3d(:, 3), 'ro')
    % uspace(qq)
    % button_is_enter = false ;
    % while ~button_is_enter
    %     button = waitforbuttonpress() ;
    %     if button && ~strcmp(get(gcf, 'CurrentKey'), 'return')
    %         disp('stopping!')
    %         error('stop')
    %     else
    %         disp('continuing...')
    %         button_is_enter = true ;
    %     end
    % end
    % clf
        
    % Note: argument structure is
    %   --> interpolate2Dpts_3Dmesh(cutMeshrs.f, cutMeshrs.u, 
    %   -->                         cutMeshrs.v, uv) 
    
    % Used to do simple search fmin, now do constrained
    % phi0s(qq) = fminsearch(@(phi0)...
    %     sum(vecnorm(...
    %     interpolate2Dpts_3Dmesh(TF, TV2D, ...
    %         TV3Drs, [uspace(qq) * ones(nV, 1), mod(vqq + phi0(1), 1)]) ...
    %         - prev3dvals, 2, 2) .^ 2), [0.], options);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~resample_hoops
        %% Average over hoop to match position in hoop rotation
        % This works quite well but does have some drift
        phi0s(qq) = fminbnd(@(phi0)...
            sum(vecnorm(...
            interpolate2Dpts_3Dmesh(TF, TV2D, ...
                TV3D, [uspace(qq) * ones(nV, 1), mod(vqq + phi0(1), 1)] ) ...
                - prev3dvals, 2, 2) .^ 2), lowerbound, upperbound, options);
            
       if compute_residual
           residuals(qq, :) = vecnorm(...
            interpolate2Dpts_3Dmesh(TF, TV2D, ...
                TV3D, [uspace(qq) * ones(nV, 1), mod(vqq + phi0s(qq), 1)] ) ...
                - prev3dvals, 2, 2) ;
       end
    else        
        %% Try resampling each hoop every time
        % check that startpt is endpt for supplied prev3dvals 
        approx_dv = mean(vecnorm(diff(prev3dvals, 1), 2, 2)) ;
        ept_is_spt = all(abs(prev3dvals(1, :) - prev3dvals(end, :)) < 1e-6 * approx_dv) ;

        NN = size(prev3dvals, 1) ;
        closed_curv = true ;
        utmp = uspace(qq) * ones(nV, 1) ;
        if ept_is_spt
            phi0s(qq) = fminbnd(@(phi0)...
                sum(vecnorm(...
                openClosedCurve(resampleCurvReplaceNaNs( ...
                    interpolate2Dpts_3Dmesh(TF, TV2D, ...
                        TV3D, [utmp, mod(vqq + phi0(1), 1)] ), ...
                    NN, closed_curv)) ...
                - prev3dvals(1:end-1, :), 2, 2) .^ 2), lowerbound, upperbound, options);
        else
            error('Have not handled case when hoop is not closed curve.')
        end
        
        % if the residual is requested, compute it
        if compute_residual
            residuals(qq, :) = vecnorm(...
                openClosedCurve(resampleCurvReplaceNaNs( ...
                    interpolate2Dpts_3Dmesh(TF, TV2D, ...
                        TV3D, [utmp, mod(vqq + phi0s(qq), 1)] ), ...
                    NN, closed_curv)) ...
                - prev3dvals(1:end-1, :), 2, 2) ;
            
            % DEBUG
            % if qq == 50 && false
            %     uspace(qq)
            %     max(TV2D)
            %     disp('debugging here')
            %     new3d = openClosedCurve(resampleCurvReplaceNaNs( ...
            %         interpolate2Dpts_3Dmesh(TF, TV2D, ...
            %             TV3D, [utmp, mod(vqq + phi0s(qq), 1)] ), ...
            %         NN, closed_curv)) ;
            %     scatter3(prev3dvals(1:end-1, 1), ...
            %         prev3dvals(1:end-1, 2), prev3dvals(1:end-1, 3), 10, 1:99)
            %     hold on;
            %     scatter3(new3d(:, 1), new3d(:, 2), new3d(:, 3), 20, 1:99)
            %     hoop1 = squeeze(prev3d_sphi(inds(1), 1:end-1, :)) ;
            %     hoop2 = squeeze(prev3d_sphi(inds(2), 1:end-1, :)) ;
            %     hoop3 = squeeze(prev3d_sphi(max(inds)+1, 1:end-1, :)) ;
            %     scatter3(hoop1(:, 1), hoop1(:, 2), hoop1(:, 3), 5, 1:99)
            %     scatter3(hoop2(:, 1), hoop2(:, 2), hoop2(:, 3), 5, 1:99)
            %     scatter3(hoop3(:, 1), hoop3(:, 2), hoop3(:, 3), 5, 'r')
            %     pause(2)
            % 
            %     figure
            %     plot3(avgpts(:, 1), avgpts(:, 2), avgpts(:, 3), 'k.')
            %     hold on;
            %     plot3(prev_avgpts(:, 1), prev_avgpts(:, 2), prev_avgpts(:, 3), 'ro')
            % 
            % 
            %     figure
            %     plot(prev_ss, vecnorm(avgpts(qq, :) - prev_avgpts, 2, 2), '.-')
            %     hold on;
            %     plot(prev_ss(inds), vecnorm(avgpts(qq, :) - prev_avgpts(inds, :), 2, 2), 'o')
            %     plot([ss, ss], [0, 100], '--')
            %     pause(1)
            % end
            % % end debug
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Try to eliminate drift by using just branch cut? 
    % This is more erratic, since it is optimizing a single point instead
    % of all the points on a hoop
    % phi0s(qq) = fminbnd(@(phi0)...
    %     sum(vecnorm(...
    %     interpolate2Dpts_3Dmesh(TF, TV2D, ...
    %         TV3Drs, [uspace(qq), mod(phi0(1), 1)] ) ...
    %         - prev3dvals(1, :), 2, 2) .^ 2), lowerbound, upperbound, options);
        
    %% Visualize the minimization output values
    if preview
        % Plot the phi_0s computed thus far
        figure(1)
        plot(phi0s)
        title(['computed \phi_0: ' num2str(qq) ' / ' num2str(nU)])
        xlabel('index of uspace')
        ylabel('\phi_0')
        
        % Plot the points being adjusted
        % disp('phiOffsetsFromPrevMesh: casting into 3d')
        
        tmp = interpolate2Dpts_3Dmesh(TF, TV2D, ...
            TV3D, [uspace(qq) * ones(nV, 1), mod(vqq + phi0s(qq), 1)]) ;
        %% Try resampling each hoop every time
        tmp = resampleCurvReplaceNaNs(tmp, NN, closed_curv) ;
        
        figure(2)
        % colormap jet
        subs = 1:10:length(prev3dvals(:, 1)) ;
        scatter3(prev3dvals(subs, 1), prev3dvals(subs, 2), prev3dvals(subs, 3), ...
            6, linspace(0, 1, length(subs))) ;
        hold on;
        subs = 1:10:length(tmp(:, 1)) ;
        scatter3(tmp(subs, 1), tmp(subs, 2), tmp(subs, 3), ...
            4, vqq(subs), 'filled') 
        title('small/filled=new, big/open=old')
        axis equal
        pause(0.000000001)
    end 
    
    runtimeIter = toc ;
    % Display progress bar
    if mod(qq, round(nU * 0.2)) == 1
        prog(min(length(prog), max(1, floor(qq/10)))) = '*' ;
        fprintf([prog '(' num2str(runtimeIter) 's per u value)\n'])
    end
end
% if visualize
%     close all
% end
end

% returns phi0s
