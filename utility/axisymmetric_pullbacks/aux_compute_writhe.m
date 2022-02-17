function [Wr, Wr_density, dWr, Length_t, clines_resampled] = ...
    aux_compute_writhe(clineDVhoopBase, timepoints, filter_curve,...
    omit_endpts, flipy, preview)
%AUX_COMPUTE_WRITHE 
% auxiliary function for computing writhe in Axisymmetric Orbifold Pullback
% Pipeline.
% 
% Parameters
% ----------
% clineDVhoopBase : str
% timepoints : int array
% filter_curve : int
%   window size for polynomial filter (cubic) and fitting
% omit_endpts : int or two ints
%   how many points on either end of the centerline to omit in writhe
%   calculation. If length == 2, then interpret as [anterior endpoint 
%   omission, posterior endpoint omission].
% flipy : bool
%   invert the sign of the y coordinate of the curve
% preview : bool
%   
% Returns
% -------
% 
% 
% NPMitchell 2020

nTPs = length(timepoints) ;

% Iterate through each mesh
Wrp = zeros(nTPs, 1) ;
wrp_densities = cell(nTPs, 1);
Wrpseg1 = zeros(nTPs, 1) ;
WrL = zeros(nTPs, 1) ;
wrL_densities = cell(nTPs, 1);
WrG = zeros(nTPs, 1) ;
wrG_densities = cell(nTPs, 1);
lengths = zeros(nTPs, 1) ;
clines_resampled = cell(nTPs, 1) ;

if preview
    close all
    fig = figure('Visible', 'on') ;
    pausetime = 1;
end

for ii = 1:nTPs 
    t = timepoints(ii) ;
    % load the centerline for this timepoint
    fn = sprintf(clineDVhoopBase, t) ;
    disp(['t = ' num2str(t)])
    disp(['Loading DVhoop centerline from ' fn])
    load(fn, 'mss', 'mcline', 'avgpts')
    if flipy
        mcline(:, 2) = -mcline(:, 2) ;
    end

    % Load centerline 
    ds = diff(mss) ; 
    ds = [ds; ds(length(ds))] ;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute writhe of the centerline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if filter_curve < 1
        ssx = mss ;
        xpt = mcline(:, 1) ;
        ypt = mcline(:, 2) ;
        zpt = mcline(:, 3) ;
    else
        % Filter the result
        % Smoothing parameters
        framelen = filter_curve ;  % must be odd
        polyorder = 3 ;
        [ssx, xpt, ypt, zpt, coeffs] = ...
            smooth_curve_via_fit_3d(mss, mcline, polyorder, framelen, 21) ;

        % Check the smoothing
        if preview && mod(ii, 10) == 1
            % Save the difference
            % fig = figure;
            % set(fig, 'Visible', 'Off')
            clf
            xcoeffs = coeffs.xcoeffs ;
            ycoeffs = coeffs.ycoeffs ;
            zcoeffs = coeffs.zcoeffs ;
            Mux = coeffs.Mux ;
            Muy = coeffs.Muy ;
            Muz = coeffs.Muz ;
            ox = polyval(xcoeffs, (mss - Mux(1)) / Mux(2)) ;
            oy = polyval(ycoeffs, (mss - Muy(1)) / Muy(2)) ;
            oz = polyval(zcoeffs, (mss - Muz(1)) / Muz(2)) ;
            plot(mss, ox - mcline(:, 1))
            hold on
            plot(mss, oy - mcline(:, 2))
            plot(mss, oz - mcline(:, 3))
            title('\Delta = smoothing - original curve')
            ylabel(['\Delta [\mu' 'm]'])
            xlabel(['s [\mu' 'm]'])
            % set(fig, 'PaperUnits', 'centimeters');
            % set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
            % waitfor(gcf)
            pause(pausetime)
            
            % Save the polynomial fit
            % fig = figure;
            % set(fig, 'Visible', 'Off')
            clf
            plot3(mcline(:, 1), mcline(:, 2), mcline(:, 3), '-.')
            hold on 
            plot3(xpt, ypt, zpt, 's')
            title('Polynomial fit to the centerline')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            % % save the figure
            % xlim([xsmin xsmax])
            % ylim([ysmin ysmax])
            % zlim([zsmin zsmax])
            axis equal
            % set(fig, 'PaperUnits', 'centimeters');
            % set(fig, 'PaperPosition', [0 0 xwidth ywidth]);
            % waitfor(gcf)
            pause(pausetime)
        end
    end

    % %% Get tangent, normal, binormal
    % [tangent, normal, binormal] = frenetSeretFrame(ss, xpt, ypt, zpt) ;
    % 
    % % Now compute how the normal changes along the AP axis    
    % v1 = normal(1:end-1, :) ;
    % v2 = normal(2:end, :) ;
    % dphi = acos(sum(v1 .* v2, 2)) ;
    % signphi = sign(sum(tangent(1:end-1, :) .* cross(v1, v2), 2)) ;
    % dphi = dphi .* signphi ;
    % chirality = dphi ./ dsx(1:end-1)' ;
    % 
    % % Check the angle change along AP axis
    % if preview
    %     close all
    %     fig = figure ;
    %     set(fig, 'Visible', 'Off')
    %     plot(ssx(1:end-1), dphi)
    %     title('$d\phi/ds$', 'Interpreter', 'Latex')
    %     ylabel('$d\phi/ds$', 'Interpreter', 'Latex')
    %     xlabel('$s$ [$\mu$m]', 'Interpreter', 'Latex')
    %     saveas(fig, [figsmoutname '_dphids.png'])
    % end

    % Compute the writhe
    % Wr = 1/4pi \int \int T(s) xx T(s') \cdot [R(s) - R(s') / |R(s) - R(s')|^3]
    % Here use the polynomial fit to the curve to compute
    % xpt, ypt, zpt
    % first compute vec from each point to every other point
    if size(xpt, 2) == 1
        yzxp = [ypt, zpt, xpt] ;
        xyzp = [xpt, ypt, zpt] ;
    else
        yzxp = [ypt', zpt', xpt'] ;
        xyzp = [xpt', ypt', zpt'] ;
    end
    
    clines_resampled{ii}.ssx = ssx ;
    clines_resampled{ii}.xyzp = xyzp ;
    % clines_resampled{ii}.coeffs = coeffs ;
    
    % traditional writhe:
    % wr = zeros(length(ssx), 1) ;
    % for jj=1:length(ssx)
    %     oind = setdiff(1:length(ssx), jj) ;
    %     rmr = xyzp(jj, :) - xyzp(oind, :) ;
    %     rmrmag = vecnorm(rmr')' ;
    %     txt = cross(tangent(jj,:) .* ones(length(oind), 3), tangent(oind, :));
    %     % Take row-wise inner product
    %     integrand = sum(sum(txt .* rmr, 2) ./ (rmrmag.^3 .* ones(size(txt))), 2) ;
    %     % Writhe per unit length is wr
    %     wr(jj) = sum(integrand) ;
    % end
    
    if length(omit_endpts) == 1
        keep = omit_endpts:(length(yzxp)-omit_endpts) ;
    else
        keep = omit_endpts(1):(length(yzxp)-omit_endpts(2)) ;
    end
    ssxkeep = ssx(keep) ;
    [Wrp(ii), wrp_local, wrp_nl, turns, segments, seg_pairs] = polarWrithe(yzxp(keep, :), ssxkeep(:)) ;
    % Store the first segment writhe: sum of wrp_local of first segment
    if ~isempty(turns)
        Wrpseg1(ii) = nansum(wrp_local(1:turns(1))) ;
    else
        Wrpseg1(ii) = nansum(wrp_local) ;
    end
    wrp_densities{ii} = wrp_local ;

    % Segment cross product writhe
    [WrL(ii), wrL_densities{ii}] = writheLevitt(xyzp(keep, :), false) ;
    [WrG(ii), wrG_densities{ii}] = writheGaussIntegral(xyzp(keep, :), ssxkeep(:)) ;

    lengths(ii) = max(mss) ;
    if preview && mod(ii, 10) == 1
        % Plot the writhe 
        % fig = figure('Visible', 'Off');
        clf
        scatter3(yzxp(keep, 1), yzxp(keep, 2), yzxp(keep, 3), 55, wrp_local', 'filled')
        cb = colorbar ;
        caxis([-0.01, 0.01])
        colormap(blueblackred)
        title('Writhe density')
        xlabel('ap position [$\mu$m]', 'Interpreter', 'Latex') ;
        ylabel('lateral position [$\mu$m]', 'Interpreter', 'Latex') ;
        zlabel('dv position [$\mu$m]', 'Interpreter', 'Latex') ;
        ylabel(cb, 'writhe density, $wr$ [$\mu$m$^{-1}$]', 'Interpreter', 'Latex') ;
        axis equal
        %  waitfor(gcf)
        pause(pausetime)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Wr(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the writhe
windowSize = 7; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
WrPsm = smoothdata(Wrp, 'rlowess', 5) ;
wPsmooth = filter(b, a, WrPsm) ;
dwrP = gradient(wPsmooth) ;
% write_txt_with_header(wrpfn, ...
%     [timestamps, Wrp, dwrP], 'timestamp, polar Writhe, smoothed d(Wrp)/dt');
% save(wdpfn, 'wrp_densities') ;

% Save the Levitt writhe
windowSize = 7; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
WrLsm = smoothdata(WrL, 'rlowess', 5) ;
wsmooth = filter(b, a, WrLsm) ;
dwrL = gradient(wsmooth) ;
% write_txt_with_header(wrLfn, ...
%     [timestamps, WrL, dwrL], 'timestamp, Levitt Writhe, smoothed d(WrL)/dt');
% save(wdLfn, 'wrL_densities') ;

% Save the Gauss writhe
windowSize = 7; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
WrGsm = smoothdata(WrG, 'rlowess', 5) ;
wsmooth = filter(b, a, WrGsm) ;
dwrG = gradient(wsmooth) ;
% write_txt_with_header(wrGfn, ...
%     [timestamps, WrG, dwrG], 'timestamp, Gauss Writhe, smoothed d(WrG)/dt');
% save(wdGfn, 'wrG_densities') ;

% Save length vs time and dlength 
% Filter the data for derivative
windowSize = 7; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
lsm = smoothdata(lengths, 'rlowess', 5) ;
lsmooth = filter(b, a, lsm) ;
dl = gradient(lsmooth) ;
% save(fullfile(meshdir, 'lengths_over_time.mat'), 'lengths', 'dl') ;


% Dump output into structs
Wr.polar = Wrp ;
Wr.polar_seg1 = Wrpseg1 ;
Wr.Levitt = WrL ;
Wr.Gauss = WrG ;
Wr_density.polar = wrp_densities ;
Wr_density.Levitt = wrL_densities ;
Wr_density.Gauss = wrG_densities ;
dWr.polar = dwrP ;
dWr.Levitt = dwrL ;
dWr.Gauss = dwrG ;
Length_t.lengths = lengths ;
Length_t.dl = dl ;

disp('Done loading/computing chirality/writhe')

