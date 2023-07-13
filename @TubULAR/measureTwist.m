function measureTwist(tubi, options)
%measureTwist(tubi, options)
% Measure twist of flow field from pullback: dv_phi / dzeta, where v_phi 
%   is the angular velocity (in radians / timeUnits, NOT in spaceUnits /
%   timeUnits). This is obtained by simply taking the derivative of the 
%   2d velocity field in the (s,phi) coordinate system wrt the longitudinal 
%   coordinate. 
%   
% We might consider projecting onto the instantaneous direction
%   perpendicular to the longitudinal axis to disambiguate orthogonal
%   contributions, but we have not done this as of 2023.
%   
% Parameters
% ----------
% tubi : TubULAR class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   writheStyle : str
%       style of writhe computation to compare to twist
%
% NPMitchell 2020-2023

%% Unpack options
% pivMethod = 'Lagrangian' ; % 'Lagranfian' or 'simpleAvg' ;
% This is a retired option, where  pivMethod: 'Lagrangian' or 'simpleAvg' 
%       whether to time-average along advected pathlines in some (small) 
%       time window at each point in spacetime ('Lagrangian') or average
%       the velocity at each (s,phi) coordinate along some (small) time
%       window with other timepoints' same (s, phi) coordinate
%       ('simpleAvg'). This is no longer used, as we always use Lagrangian
%       averaging. 

writheStyle = 'Levitt' ;
overwrite = false ;
foldW = 0.04 ;  % full width of fraction of zeta attributed to each fold

if nargin == 1
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite  ;
end
if isfield(options, 'writheStyle')
    writheStyle = options.writheStyle ;
end
if isfield(options, 'foldW')
    foldW = options.foldW ;
end

%% Unpack tubi
nU = tubi.nU ;
nV = tubi.nV ;
timePoints = tubi.xp.fileMeta.timePoints ;
outdir = fullfile(tubi.dir.piv.avg, 'twist') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
% piv3dfn = tubi.fullFileBase.piv3d ;
t0 = tubi.t0set() ;

%% Velocity -> twist
% (1) Load face-based velocity measurements
% (2) Resolve into tangential component
% (3) Compute deriv wrt proper length along longitudinal axis

nfaces = 2 * (nU - 1) * (nV - 1) ;
ntps = length(timePoints) - 1 ;

% Load velocity
vels = tubi.getVelocityAverage() ;
vfsmM = vels.vf ;

vv2dMfn = tubi.fileName.pivAvg.vv2d ;
% Used to be: fullfile(outdir, 'v0t2dvsmM_simpletimeavg.mat') ;
vf2dMfn = tubi.fileName.pivAvg.vf2d ;
% Used to be: fullfile(outdir, 'v0t2dfsmM_simpletimeavg.mat') ;
twistMfn = tubi.fileName.pivAvg.twist ;
% Used to be: fullfile(outdir, 'twistM_simpletimeavg.mat') ;
redo_calc = overwrite || ~exist(vv2dMfn, 'file') || ...
    ~exist(vf2dMfn, 'file') ;
redo_twist_only = ~exist(twistMfn, 'file') ;

if redo_calc 
    % preallocate v0t2dvM, v0t2dfM 
    v0t2dvM = zeros(ntps, nU * (nV-1), 2) ;
    v0t2dfM = zeros(ntps, nfaces, 2) ;
    twistM = zeros(ntps, nU) ;

    % Consider each timepoint
    for tidx = 1:ntps
        tp = timePoints(tidx) ;
        disp(['t = ', num2str(tp)])

        %% OLD METHOD
        % load(sprintf(piv3dfn, tp), 'piv3dstruct') ;
        % v0t2d = piv3dstruct.v0t2d ;
        % vvals = reshape(v0t2d, [nU, nV]) ;

        %% NEW METHOD
        % Load smoothed rotated scaled closed mesh
        load(sprintf(tubi.fullFileBase.spcutMeshSmRSC, tp), 'spcutMeshSmRSC') ;
        mesh = spcutMeshSmRSC ;
        fieldfaces = (1:nfaces)' ;          % one field for each face

        % Resolve into tangential component
        vf = squeeze(vfsmM(tidx, :, :)) ;
        [v0n, v0t, v0t2d_faces, jac] = ...
            resolveTangentNormalVelocities(mesh.f, mesh.v, vf, fieldfaces, mesh.u) ;

        %% Compute deriv wrt proper length along longitudinal axis (zeta)
        [V2F, F2V] = meshAveragingOperators(mesh.f, mesh.v) ;
        v0t2d_vertices = F2V * v0t2d_faces ;

        %% Save v0t2d_vertices and v0t2d_faces
        v0t2dvM(tidx, :, :) = v0t2d_vertices ;
        v0t2dfM(tidx, :, :) = v0t2d_faces ;

        %% Reshape and take gradients along zeta
        v0t2d_vertices = reshape(v0t2d_vertices, [nU, nV - 1, 2]) ;
        vphi = mean(v0t2d_vertices(:, 2), 2) ;
        twist = gradient(vphi) ;

        %% Save twist measurement
        twistM(tidx, :) = twist ;
    end

    %% Save all results
    save(vv2dMfn, 'v0t2dvM')
    save(vf2dMfn, 'v0t2dfM')
    save(twistMfn, 'twistM')
elseif redo_twist_only
    % The other calculations are already on disk, just recompute twist
    load(vv2dMfn, 'v0t2dvM')
    twistM = zeros(ntps, nU) ;
    
    % Consider each timepoint
    for tidx = 1:ntps
        v0t2d_vertices = v0t2dvM(tidx, :, :) ;
        
        % Reshape and take gradients along zeta
        v0t2d_vertices = reshape(v0t2d_vertices, [nU, nV - 1, 2]) ;
        vphi = mean(v0t2d_vertices(:, 2), 2) ;
        twist = gradient(vphi) ;
        
        % Save twist measurement
        twistM(tidx, :) = twist ;
    end
    save(twistMfn, 'twistM')
else
    load(twistMfn, 'twistM')
end

%% Prepare for plots
dtwLabel = '$\partial_{\zeta} \langle v_\phi \rangle_{\mathrm{dv}}$' ;
twLabel = '$\int_0^t \mathrm{d}t \, \partial_{\zeta} \langle v_\phi \rangle_{\mathrm{dv}}$';
dwrLabel = '$\partial_t \mathrm{Wr}$' ;
wrLabel = '$\mathrm{Wr}$' ;

% Load writhe
tmp = load(tubi.fileName.writhe) ;
if strcmp(writheStyle, 'Levitt')
    Wr = tmp.Wr.Levitt ;
    dWr = tmp.dWr.Levitt(1:ntps) ;
else
    error(['Handle this style here: ' writheStyle])
end


%% Save twist measurement figure
close all
xL = (1:nU) / nU ;
tps = timePoints(1:end-1) - t0 ; 
imagesc(xL, tps, twistM)
cb = colorbar() ;
title('Tissue twist rate')
ylabel(cb, dtwLabel, 'Interpreter', 'latex')
colormap bwr
caxis([-2, 2]) 
xlabel('ap position, $\zeta$', 'Interpreter', 'latex')
ylabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'latex')
saveas(gcf, fullfile(outdir, 'dTw_kymograph.png'))
close all

% %% Filter image
% selectIdx = xL > 0.1 & xL < 0.9 ;
% xLs = xL(selectIdx) ;
% hh = [1 2 3 4 5 6 7 8 7 6 5 4 3 2 1] ;
% hh = hh / sum(hh) ;
% im = imfilter(twistM(:, selectIdx), hh, 'circular') ;
% imagesc(xLs, tps, im) ;
% colormap bwr
% caxis([-0.5, 0.5]) 
% cb = colorbar() ;
% title('Tissue twist rate')
% ylabel(cb, dtwLabel, 'Interpreter', 'latex')
% xlabel('ap position, $\zeta$', 'Interpreter', 'latex')
% ylabel(['time [' tubi.timeUnits ']'], 'Interpreter', 'latex')
% saveas(gcf, fullfile(outdir, 'dTw_filt_kymograph.png'))
% close all


%% Plot averages in bins
clf

% Load folds
tubi.getFeatures()
folds = tubi.features.folds ;

% fold half width in integers sampling zeta axis
fhw = foldW * 0.5 * nU ;
if size(folds, 2) == 3
    bins = cat(2, 0*folds(:, 1), ...
                folds(:, 1) - fhw, folds(:, 1) + fhw, ...
                folds(:, 2) - fhw, folds(:, 2) + fhw, ...
                folds(:, 3) - fhw, folds(:, 3) + fhw, ...
                    nU * ones(size(folds(:, 1)))) ;
    labels = {'anterior lobe', 'anterior fold', ...
        'second lobe', 'middle fold', 'third lobe', ...
        'posterior fold', 'fourth lobe'} ;
else
    error('handle case of >3 or <3 folds here.')
end
bins = max(1, min(bins, nU)) ;
colors = tubi.plotting.colors ;
markers = tubi.plotting.markers ;
Twnets = cell(size(bins, 2) - 1, 1) ; 
alphaVal = 0.6 ;
dmyk = 0 ;
for qq = 1:(size(bins, 2) - 1)
    for tidx = 1:ntps
        % get first zeta index and last zeta index for this bin (a lobe or
        % a fold)
        start = bins(tidx, qq) ;
        finish = bins(tidx, qq + 1) ;
        dtwnet(tidx) = sum(twistM(tidx, start:finish), 2) / nU;
        Twnet(tidx) = sum(dtwnet(1:tidx)) ;
    end
    twfilt = savgol(dtwnet, 2, 21) ;
    if contains(labels{qq}, 'lobe')
        dmyk = dmyk + 1 ;
        ph(qq) = plot(tps, twfilt, 'linestyle', 'none', ...
            'marker', markers{dmyk}, 'color', colors(dmyk, :)) ;
    else
        ph(qq) = plot(tps, twfilt, '-', ...
            'color', colors(dmyk, :)) ;        
    end
    hold on;
    % Raw data
    % pskip(qq) = plot(tps, twnet, '--', 'color', colors(qq, :))  ;
    
    % Store sums 
    Twnets{qq} = Twnet ;    
end
legend(ph, labels, 'location', 'northeastoutside')
title(['Tissue twist rate, ' dtwLabel], 'Interpreter', 'Latex')
ylabel(dtwLabel, 'Interpreter', 'Latex')
xlabel(['time [' tubi.timeUnits ']'])
saveas(gcf, fullfile(outdir, 'dTw_vs_time.png'))

% Add dWrithe
ph(length(ph)+1) = plot(tps, dWr, 'k-', 'linewidth', 4) ;
labels{length(labels)+1} = dwrLabel ;
ylabel(['Tissue twist rate (', dtwLabel, '), or writhing (', dwrLabel, ')'], ...
    'Interpreter', 'Latex')
legend(ph, labels, 'Interpreter', 'Latex')
saveas(gcf, fullfile(outdir, 'dTw_dWr_vs_time_bins.png'))

%% Plot sums (accumulated Writhe)
close all
clearvars ph
hold on;
dmyk = 0 ;
for qq = 1:(size(bins, 2) - 1)
    Twnet = Twnets{qq} ;
    if contains(labels{qq}, 'lobe')
        dmyk = dmyk + 1 ;
        ph(qq) = plot(tps, Twnet, 'linestyle', 'none', ...
            'marker', markers{dmyk}, 'color', colors(dmyk, :)) ;
    else
        ph(qq) = plot(tps, Twnet, '-', ...
            'color', colors(dmyk, :)) ;   
    end
end
legend(ph, labels, 'location', 'northeastoutside')
title(['Tissue twist, ' twLabel], 'Interpreter', 'Latex')
ylabel(dtwLabel, 'Interpreter', 'Latex')
xlabel(['time [' tubi.timeUnits ']'])
saveas(gcf, fullfile(outdir, 'Tw_vs_time.png'))

%% Add Writhe
ph(length(ph)+1) = plot(tps, Wr(1:ntps), 'k-', 'linewidth', 4) ;
labels{length(labels)+1} = wrLabel ;
ylabel(['Tissue twist (', twLabel, '), and Writhe (', wrLabel, ')'], ...
    'Interpreter', 'Latex')
legend(ph, labels, 'Interpreter', 'Latex')
saveas(gcf, fullfile(outdir, 'Tw_Wr_vs_time_bins.png'))

%% Total twist rate
clf
clearvars ph labels
dtwnet = sum(im, 2) / size(im, 2) ;
twfilt = savgol(dtwnet, 2, 21) ;
ph(1) = plot(tps, dtwnet, '.', 'color', colors(1, :)) ;
hold on;
plot(tps, twfilt, '--', 'color', colors(1, :))  ;
title('Net tissue twist rate, $\langle\partial_{\zeta} \langle v_\phi \rangle_{\mathrm{dv}}\rangle_{\mathrm{ap}}$', ...
'Interpreter', 'Latex')
ylabel('$\langle\partial_{\zeta} \langle v_\phi \rangle_{\mathrm{dv}}\rangle_{\mathrm{ap}}$', ...
'Interpreter', 'Latex')
xlabel(['time [' tubi.timeUnits ']'])
saveas(gcf, fullfile(outdir, 'dTw_vs_time_total.png'))

%% Add Writhe
labels{1} = twLabel ;
ph(2) = plot(tps, Wr(1:ntps), 'k-', 'linewidth', 1) ;
labels{2} = wrLabel ;
ylabel(['Tissue twist (', twLabel, '), and Writhe (', wrLabel, ')'], ...
    'Interpreter', 'Latex')
legend(ph, labels, 'Interpreter', 'Latex', 'location', 'northeastoutside')
saveas(gcf, fullfile(outdir, 'Tw_Wr_vs_time_total.png'))

%% Take selective sums 
Twnet_lobes = Twnets{1} + Twnets{3} + Twnets{5} + Twnets{7} ;
Twnet_folds = Twnets{2} + Twnets{4} + Twnets{6} ;
close all
ph(1) = plot(tps, Twnet_lobes, '-', 'linewidth', 1, 'color', colors(1, :)) ;
hold on;
labels{1} = [twLabel ' in lobes'];
ph(2) = plot(tps, Twnet_folds, '-', 'linewidth', 1, 'color', colors(2, :)) ;
labels{2} = [twLabel ' in folds'];
ph(3) = plot(tps, Wr(1:ntps), 'k-', 'linewidth', 1) ;
labels{3} = wrLabel ;
ylabel(['Tissue twist (', twLabel, '), and Writhe (', wrLabel, ')'], ...
    'Interpreter', 'Latex')
legend(ph, labels, 'Interpreter', 'Latex', 'location', 'northeastoutside')
saveas(gcf, fullfile(outdir, 'Tw_Wr_vs_time_total_selective.png'))



