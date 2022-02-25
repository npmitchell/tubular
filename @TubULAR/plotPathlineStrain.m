function plotPathlineStrain(QS, options)
% plotPathlineStrain(QS, options)
%   Plot the strain (from piv pathlines) along pathlines as kymographs and 
%   correlation plots. These are computed via finite differencing of 
%   pathline mesh metrics. That is, vertices are advected along piv flow 
%   and projected into 3D from pullback space to embedding, then the
%   metrics g_{ij} of those meshes are measured and differenced. 
% 
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   plot_kymographs         : bool
%   plot_kymographs_cumsum  : bool
%   plot_correlations       : bool 
%   plot_fold_strainRate    : bool
%   plot_lobe_strainRate    : bool
%   plot_fold_strain        : bool
%   plot_lobe_strain        : bool
% 
% Returns
% -------
% <none>
%
% NPMitchell 2020

%% Default options 
overwrite = false ;
plot_kymographs = true ;
plot_kymographs_strain = true ;
plot_fold_strainRate = true ;
plot_lobe_strainRate = true ;
plot_fold_strain = true ;
plot_lobe_strain = true ;
maxWFrac = 0.03 ;  % Maximum halfwidth of feature/fold regions as fraction of ap length
t0 = QS.t0set() ;
t0Pathline = t0 ;

strain_trace_label = '$\frac{1}{2}\mathrm{Tr} [\bf{g}^{-1}\varepsilon]$';
strain_trace_label_norm = '$||\mathrm{Tr} [\varepsilon]||$';
strain_deviator_label = ...
    '$||\varepsilon-\frac{1}{2}$Tr$\left[\mathbf{g}^{-1}\varepsilon\right]\bf{g}||$' ;
strain_deviator_label_short = ...
    '$||$Dev$\left[\varepsilon\right]||$' ;
% strain_theta_label = '$\theta_{\mathrm{Dev}[{\varepsilon}]}$'; 

%% Parameter options
lambda = QS.smoothing.lambda ;
lambda_mesh = QS.smoothing.lambda_mesh ;
nmodes = QS.smoothing.nmodes ;
zwidth = QS.smoothing.zwidth ;
climit = 0.2 ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
%% parameter options
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
end

% smoothing parameter options
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'nmodes')
    nmodes = options.nmodes ;
end
if isfield(options, 'zwidth')
    zwidth = options.zwidth ;
end

% plotting and sampling
if isfield(options, 'climit')
    climit = options.climit ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'climitWide')
    climitWide = options.climitWide ;
else
    climitWide = climit * 3; 
end
if isfield(options, 'maxWFrac')
    maxWFrac = options.maxWFrac ;
end

%% Operational options
if isfield(options, 'plot_kymographs')
    plot_kymographs = options.plot_kymographs ;
end
if isfield(options, 'plot_kymographs_strain')
    plot_kymographs_strain = options.plot_kymographs_strain ;
end
if isfield(options, 'plot_fold_strain')
    plot_fold_strain = options.plot_fold_strain ;
end
if isfield(options, 'plot_lobe_strain')
    plot_lobe_strain = options.plot_lobe_strain ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
else
    error("Could not parse samplingResolution: set to '1x' ")
end

%% Unpack QS
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
bbr256 = blueblackred(256) ;
% clf
% set(gcf, 'visible', 'off')
% imagesc([-1, 0, 1; -1, 0, 1])
% caxis([0, 1])
% pos256 = bluewhitered(256) ;
close all
pm256 = phasemap(256) ;
bluecolor = QS.plotting.colors(1, :) ;
orangecolor = QS.plotting.colors(2, :) ;
yellowcolor = QS.plotting.colors(3, :) ;
purplecolor = QS.plotting.colors(4, :) ;
greencolor = QS.plotting.colors(5, :) ;
graycolor = QS.plotting.colors(8, :) ;
browncolor = QS.plotting.colors(9, :) ;

%% Choose colors
trecolor = yellowcolor ;
Hposcolor = greencolor ;
Hnegcolor = purplecolor ;
Hsz = 3 ;  % size of scatter markers for mean curvature

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

%% Test incompressibility of the flow on the evolving surface
% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - t0 ;

% Output directory is inside metricKinematics dir
mKPDir = sprintf(QS.dir.strainRate.pathline.root, t0Pathline) ;
datdir = fullfile(mKPDir, 'measurements') ;
% Data for kinematics on meshes (defined on vertices) [not needed here]
% mdatdir = fullfile(mKDir, 'measurements') ;

% Unit definitions for axis labels
unitstr = [ '[unitless]' ];

%% Compute or load all timepoints
apSKymoFn = fullfile(datdir, 'apKymographPathlineStrain.mat') ;
lSKymoFn = fullfile(datdir, 'leftKymographPathlineStrain.mat') ;
rSKymoFn = fullfile(datdir, 'rightKymographPathlineStrain.mat') ;
dSKymoFn = fullfile(datdir, 'dorsalKymographPathlineStrain.mat') ;
vSKymoFn = fullfile(datdir, 'ventralKymographPathlineStrain.mat') ;
files_exist = exist(apSKymoFn, 'file') && ...
    exist(lSKymoFn, 'file') && exist(rSKymoFn, 'file') && ...
    exist(dSKymoFn, 'file') && exist(vSKymoFn, 'file') ;
if files_exist
    load(apSKymoFn, 'tr_apM', 'dv_apM', 'th_apM')
    load(lSKymoFn, 'tr_lM', 'dv_lM', 'th_lM')
    load(rSKymoFn, 'tr_rM', 'dv_rM', 'th_rM')
    load(dSKymoFn, 'tr_dM', 'dv_dM', 'th_dM')
    load(vSKymoFn, 'tr_vM', 'dv_vM', 'th_vM')
else
    % preallocate for kymos
    ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
    dv_apM = zeros(ntps, nU) ;   % dv averaged
    tr_apM = zeros(ntps, nU) ;   
    th_apM = zeros(ntps, nU) ;   
    dv_lM = zeros(ntps, nU) ;    % left averaged
    tr_lM = zeros(ntps, nU) ;
    th_lM = zeros(ntps, nU) ;
    dv_rM = zeros(ntps, nU) ;    % right averaged
    tr_rM = zeros(ntps, nU) ;
    th_rM = zeros(ntps, nU) ;
    dv_dM = zeros(ntps, nU) ;    % dorsal averaged
    tr_dM = zeros(ntps, nU) ;
    th_dM = zeros(ntps, nU) ;
    dv_vM = zeros(ntps, nU) ;    % ventral averaged
    tr_vM = zeros(ntps, nU) ;
    th_vM = zeros(ntps, nU) ;

    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        close all
        disp(['t = ' num2str(tp)])
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        srfn = fullfile(datdir, sprintf('strainRate_%06d.mat', tp))   ;

        % Load timeseries measurements
        load(srfn, 'tre_ap', 'tre_l', 'tre_r', 'tre_d', 'tre_v', ...
            'dev_ap', 'dev_l', 'dev_r', 'dev_d', 'dev_v', ...
            'theta_ap', 'theta_l', 'theta_r', 'theta_d', 'theta_v')
   
        %% Store in matrices
        % dv averaged
        tr_apM(tidx, :) = tre_ap ;
        dv_apM(tidx, :) = dev_ap ;
        th_apM(tidx, :) = theta_ap ;

        % left quarter
        tr_lM(tidx, :) = tre_l ;
        dv_lM(tidx, :) = dev_l ;
        th_lM(tidx, :) = theta_l ;

        % right quarter
        tr_rM(tidx, :) = tre_r ;
        dv_rM(tidx, :) = dev_r ;
        th_rM(tidx, :) = theta_r ;

        % dorsal quarter
        tr_dM(tidx, :) = tre_d ;
        dv_dM(tidx, :) = dev_d ;
        th_dM(tidx, :) = theta_d ;

        % ventral quarter
        tr_vM(tidx, :) = tre_v ;
        dv_vM(tidx, :) = dev_v ;
        th_vM(tidx, :) = theta_v ;
        
        %% Save kymograph data
        save(apSKymoFn, 'tr_apM', 'dv_apM', 'th_apM')
        save(lSKymoFn, 'tr_lM', 'dv_lM', 'th_lM')
        save(rSKymoFn, 'tr_rM', 'dv_rM', 'th_rM')
        save(dSKymoFn, 'tr_dM', 'dv_dM', 'th_dM')
        save(vSKymoFn, 'tr_vM', 'dv_vM', 'th_vM')
    end
end


%% Store kymograph data in cell arrays
trsK = {0.5*tr_apM, 0.5*tr_lM, 0.5*tr_rM, 0.5*tr_dM, 0.5*tr_vM} ;
dvsK = {dv_apM, dv_lM, dv_rM, dv_dM, dv_vM} ;
thsK = {th_apM, th_lM, th_rM, th_dM, th_vM} ;

%% Make kymographs averaged over dv, or left, right, dorsal, ventral 1/4
dvDir = fullfile(mKPDir, 'avgDV') ;
lDir = fullfile(mKPDir, 'avgLeft') ;
rDir = fullfile(mKPDir, 'avgRight') ;
dDir = fullfile(mKPDir, 'avgDorsal') ;
vDir = fullfile(mKPDir, 'avgVentral') ;
outdirs = {dvDir, lDir, rDir, dDir, vDir} ;

%% Now plot different measured quantities as kymographs
if plot_kymographs
    titleadd = {': circumferentially averaged', ...
        ': left side', ': right side', ': dorsal side', ': ventral side'} ;

    for qq = 1:length(outdirs)
        % Prep the output directory for this averaging
        odir = outdirs{qq} ;
        if ~exist(odir, 'dir')
            mkdir(odir)
        end

        % Unpack what to plot (averaged kymographs, vary averaging region)
        trK = trsK{qq} ;
        dvK = dvsK{qq} ;
        thK = thsK{qq} ;
        
        titles = {['dilation, ' strain_trace_label],...
            ['shear, ', strain_deviator_label]} ;
        labels = {[strain_trace_label ' ' unitstr], ...
            [strain_deviator_label ' ' unitstr]} ;
        names = {'strain_dilation', 'strain_deviator'} ;
        
        %% Plot strainRate DV-averaged Lagrangian pathline kymographs 
        % Check if images already exist on disk

        % Consider both wide color limit and narrow
        zoomstr = {'_wide', '', '_zoom'} ;
        climits = {climit, 0.5*climit, 0.25*climit} ;
        
        for pp = 1:length(climits)
            %% Plot STRAIN traceful DV-averaged pathline kymograph
            % Check if images already exist on disk
            fn = fullfile(odir, [ names{1} zoomstr{pp} '.png']) ;
            if ~exist(fn, 'file') || ~exist(fn_early, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                imagesc((1:nU)/nU, tps, trK)
                caxis([-climits{pp}, climits{pp}])
                colormap(bbr256)
                
                
                % Titles 
                title([titles{1}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                cb = colorbar() ;
                
                % title and save
                ylabel(cb, labels{1}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                
                clf
            end

            %% DEVIATOR -- strain
            fn = fullfile(odir, [ names{2} zoomstr{pp} '.png']) ;
            if ~exist(fn, 'file') || overwrite
                close all
                set(gcf, 'visible', 'off')
                % Map intensity from dev and color from the theta
                indx = max(1, round(mod(2*thK(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
                colors = pm256(indx, :) ;
                devKclipped = min(dvK / climits{pp}, 1) ;
                colorsM = devKclipped(:) .* colors ;
                colorsM = reshape(colorsM, [size(dvK, 1), size(dvK, 2), 3]) ;
                imagesc((1:nU)/nU, tps, colorsM)
                caxis([0, climits{pp}])

                % Add folds to plot
                hold on;

                % Colorbar and phasewheel
                colormap(gca, phasemap)
                phasebar('colormap', phasemap, ...
                    'location', [0.82, 0.7, 0.1, 0.135], 'style', 'nematic') ;
                ax = gca ;
                get(gca, 'position')
                cb = colorbar('location', 'eastOutside') ;
                drawnow
                axpos = get(ax, 'position') ;
                cbpos = get(cb, 'position') ;
                set(cb, 'position', [cbpos(1), cbpos(2), cbpos(3), cbpos(4)*0.6])
                set(ax, 'position', axpos) 
                hold on;
                caxis([0, climits{pp}])
                colormap(gca, gray)

                % title and save
                title([titles{2}, titleadd{qq}], 'Interpreter', 'Latex')
                ylabel(['time [' QS.timeUnits ']'], 'Interpreter', 'Latex')
                xlabel('ap position [$\zeta/L$]', 'Interpreter', 'Latex')
                ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
                    
                % title and save
                ylabel(cb, labels{2}, 'Interpreter', 'Latex')  
                disp(['saving ', fn])
                export_fig(fn, '-png', '-nocrop', '-r200')   
                clf
            end
        end
    end
end


disp('done')



