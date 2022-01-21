function computeRicciMeshes(QS, options)
% computeRicciMeshes(QS, options) 
%   Compute all Ricci meshes (if flow converges) and plot aspect ratio for isothermal PB over time
%
% Parameters
% ----------
% options : optional struct with fields
%   resample : bool (default=true) ;
%   
%
% NPMitchell 2021

% Default options
maxIter = 200 ;
resample = true ;

if nargin < 2
    options = struct();
end
if isfield(options, 'maxIter')
    maxIter = options.maxIter ;
end
if isfield(options, 'resample')
    resample = options.resample ;
end

% Compute each Ricci flow mesh
opts = struct() ;
opts.maxIter = maxIter ;
opts.resample = resample ;

t0 = QS.t0set() ;
timePoints = QS.xp.fileMeta.timePoints ;
startID = 87 ;
tidx2do0 = [QS.xp.tIdx(153), ...
    QS.xp.tIdx(183), QS.xp.tIdx(213),  ...
    QS.xp.tIdx(108), QS.xp.tIdx(168), QS.xp.tIdx(206)] ;
tidx2do1 = [tidx2do0, ...
    setdiff(startID:50:length(timePoints), tidx2do0)] ;
tidx2do1 = [tidx2do1, setdiff(startID:30:length(timePoints), tidx2do1)] ;
tidx2do2 = [tidx2do1, setdiff(startID:20:length(timePoints), tidx2do1)] ;
tidx2do3 = setdiff(setdiff(startID:10:length(timePoints), ...
    tidx2do1), tidx2do2) ;
tidx2do123 = [tidx2do1, tidx2do2, tidx2do3] ;
tidx2do1234 = [tidx2do123, ...
    setdiff(startID:2:length(timePoints), tidx2do123)] ;
tidx2do = [tidx2do1234, ...
    setdiff(1:length(timePoints), tidx2do1234)] ;
tidx0 = QS.xp.tIdx(t0) ;
tidx2do = [tidx0, setdiff(tidx2do, [tidx0], 'stable')] ;
% check that there are no repeats
assert(length(tidx2do) == length(unique(tidx2do)))
assert(all(tidx2do == unique(tidx2do, 'stable')))

for tiix = 1:length(tidx2do)
    tidx = tidx2do(tiix) ;
    tp = timePoints(tidx) ;
    disp(['t = ', num2str(tp)])
    QS.setTime(tp)
    if resample
        QS.generateRicciMeshTimePoint(tp, opts) 
    else
        try
            QS.generateRicciMeshTimePoint(tp, opts) 
        catch
            disp('could not generate Ricci mesh -- self intersections?')
        end
    end
end

% Plot aspect ratios compared to conformal factor denoting inside/outside
% radius of mapped annulus
aratio_r = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
aratio_a = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
aratio_i = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
rads = zeros(length(QS.xp.fileMeta.timePoints), 1) ;
for tidx = 1:length(QS.xp.fileMeta.timePoints)
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    disp(['t = ', num2str(tp)])
    QS.setTime(tp)
    try
        rmesh = QS.loadCurrentRicciMesh() ;
        aratio_r(tidx) = max(rmesh.rectangle.u(:,1)) / (2*pi) ;
        aratio_a(tidx) = max(vecnorm(rmesh.annulus.u, 2, 2)) / min(vecnorm(rmesh.annulus.u, 2, 2)) ;
    catch
        aratio_r(tidx) = NaN ;
        aratio_a(tidx) = NaN ;
    end
    [~, radius_cline] = QS.getRadii() ;
    rads(tidx) = mean(radius_cline) ;
    
    % Also load isoareal minimization aspect ratio
    thisSPCutMesh = QS.loadCurrentSPCutMesh ;
    aratio_i(tidx) = thisSPCutMesh.ar ;
end
Length_t = QS.measureLength() ;
lengs = Length_t.lengths ;

%% plot it
clf
t0 = QS.t0set() ;
tidx0 = QS.xp.tIdx(t0) ;
tps = (QS.xp.fileMeta.timePoints - t0) * QS.timeInterval ;
min2hr = strcmpi(QS.timeUnits, 'min') && (length(tps) > 90) ;
if min2hr
    tps = tps / 60 ;
end

plot(tps, aratio_r, '.-')
hold on
plot(tps, log10(aratio_a), '.-')
plot(tps, lengs / lengs(tidx0), '.-')
plot(tps, lengs ./ (2*pi*rads), '.-')
plot(tps, aratio_i, '.-')
legend({'Ricci $L_\zeta/L_\phi$', ...
    'Ricci $\log_{10}(\textrm{max}\rho/\textrm{min}\rho)$', ...
    '$L_{\Re^3}/L_{\Re^3}^{0}$', ...
    '$L_{\Re^3}/(2\pi\langle r \rangle))$', ...
    'aspect ratio of $uv$ coords'}, ...
    'interpreter', 'latex', 'location', 'best')
title('Aspect ratios over time', 'interpreter', 'latex')
if min2hr
    xlabel('time, $t$ [hr]', 'interpreter', 'latex')
else
    xlabel(['time, $t$ [' QS.timeUnits ']'], 'interpreter', 'latex')
end
saveas(gcf, fullfile(QS.dir.ricci.data, 'aspect_ratios.pdf'))

% PLot against each other
close all
fig = figure('position', [100, 100, 300, 300], 'units', 'centimeters') ;
xx = lengs ./ (2*pi*rads) ;
xx(isnan(aratio_r)) = NaN ;
scatter( xx, aratio_r) ;
hold on;
minval = min(min(aratio_r), min(xx)) ;
maxval = max(max(aratio_r), max(xx)) ;
plot([minval, maxval], [minval, maxval], 'k--')
xlabel('$L_{\Re^3}/(2\pi\langle r \rangle))$', 'interpreter', 'latex') ;
ylabel('Ricci $L_\zeta/L_\phi$', 'interpreter', 'latex') ;
saveas(gcf, fullfile(QS.dir.ricci.data, 'ricciAR_versus_geodesicAR.pdf'))


close all
fig = figure('position', [100, 100, 300, 300], 'units', 'centimeters') ;
xx = lengs / lengs(tidx0) ;
xx(isnan(aratio_r)) = NaN ;
scatter( xx, aratio_r) ;
hold on;
minval = min(min(aratio_r), min(xx)) ;
maxval = max(max(aratio_r), max(xx)) ;
plot([minval, maxval], [minval, maxval], 'k--')
xlabel('$L_{\Re^3}/L_{\Re^3}^{0}$', 'interpreter', 'latex') ;
ylabel('Ricci $L_\zeta/L_\phi$', 'interpreter', 'latex') ;
saveas(gcf, fullfile(QS.dir.ricci.data, 'ricciAR_versus_geodesicARDelta.pdf'))


close all
fig = figure('position', [100, 100, 300, 300], 'units', 'centimeters') ;
xx = aratio_i ;
xx(isnan(aratio_r)) = NaN ;
scatter( xx, aratio_r) ;
hold on;
minval = min(min(aratio_r), min(xx)) ;
maxval = max(max(aratio_r), max(xx)) ;
plot([minval, maxval], [minval, maxval], 'k--')
xlabel('aspect ratio of $uv$ coords', 'interpreter', 'latex') ;
ylabel('Ricci $L_\zeta/L_\phi$', 'interpreter', 'latex') ;

xkeep = xx(~isnan(xx) & tps' > -0.5 ) ;
ykeep = aratio_r(~isnan(aratio_r) & tps' > -0.5 ) ;
[pp, ss] = polyfit(xkeep, ykeep, 1) ;
uncs =  sqrt(diag(inv(ss.R)*inv(ss.R')).*ss.normr.^2./ss.df) ;
title(['$y=($', sprintf('%0.3f', pp(1)), '$\pm$'...
    sprintf('%0.3f', uncs(1)), '$)x+$', ...
    sprintf('%0.3f', pp(2)), '$\pm$', ...
    sprintf('%0.3f', uncs(2))], 'interpreter', 'latex')
saveas(gcf, fullfile(QS.dir.ricci.data, 'ricciAR_versus_uvAR.pdf'))


% Save aspect ratios
aratioFn = fullfile(QS.dir.ricci.data, 'aspect_ratios.mat') ;
save(aratioFn, 'aratio_r', 'aratio_a', 'lengs', 'rads', 'aratio_i')


