function measureRelativeMotion(QS, options)
% measureRelativeMotion(QS, options)
%
% Use velocity data of each layer to compare the relative motion
% 
% 
% NPMitchell 2021

%% Unpack options
channels = 1:length(QS.xp.expMeta.channelsUsed) ;
t0 = QS.t0set() ;
vtscale = 2 ;
overwrite = false ;

if nargin < 2
    options = struct() ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
else
    overwriteImages = overwrite ;
end
if isfield(options, 'channels')
    channels = options.channels ;
    if ~assert(length(channels) == 2)
        error(['Must compare only two channels at a time! '...
            'Supply channels as field in options'])
    end
end

doubleCovered = contains(QS.pivMultiChannel.imCoords, 'e') ;


% raw PIV
% piv = QS.getPIVMultiChannel() ;
ntps = length(QS.xp.fileMeta.timePoints)-1 ;
tidx2do = 1:20:ntps ;
tidx2do = [tidx2do, setdiff(1:ntps, tidx2do)] ;

% relative motion in tangent plane (to scale in um/min)
QS.setTime(QS.xp.fileMeta.timePoints(1)) 
mesh = QS.getCurrentSPCutMeshSmRS() ;
nFaces = length(mesh.f) ;
dvtT = zeros(length(tidx2do), nFaces, 2) ;
v1tT = zeros(length(tidx2do), nFaces, 3)  ;
v2tT = v1tT ;
v1trmsV = zeros(length(tidx2do), 1) ;
v2trmsV = zeros(length(tidx2do), 1) ;
dvtrmsV = zeros(length(tidx2do), 1) ;

% Compute/load relative motion for each timepoint
for tidx = tidx2do 
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    disp(['t = ' num2str(tp)])
    
    outfn = fullfile(QS.dir.pivMultiChannel.relativeMotion, ...
        'measurements', sprintf('relativeMotion_%06d.mat', tp)) ;
    if ~exist(outfn, 'file') || overwrite
        QS.setTime(tp) 

        % 3d PIV
        vel = QS.getCurrentVelocityMultiChannel() ;

        % Load mesh and image
        im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
        mesh = QS.getCurrentSPCutMeshSmRS() ;
        mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
        mesh.u = QS.uv2XY(im0, mesh.u, doubleCovered, 1, 1) ;

        v1 = vel.piv3d{channels(1)} ;
        v2 = vel.piv3d{channels(2)} ;

        % Get tangential component of face-based velocities
        v1f = v1.v3dfaces_rs ;
        v2f = v2.v3dfaces_rs ;

        % RESOLVE COMPONENTS
        [v1n, v1t, v12d] = ...
             resolveTangentNormalVelocities(mesh.f, mesh.v, v1f, ...
             1:size(mesh.f, 1), mesh.u) ;
        [v2n, v2t, v22d] = ......
             resolveTangentNormalVelocities(mesh.f, mesh.v, v2f, ...
             1:size(mesh.f, 1), mesh.u) ;
        [dvn, dvt, dv2d, ~, ~, ~, dilation] = ...
            resolveTangentNormalVelocities(mesh.f, mesh.v, v2f - v1f, ...
            1:size(mesh.f, 1), mesh.u) ;

        %% CHECK
        % Ensure that v1t is tangent to faces
        % assert that 
        % assert(all(sum(v1f .* facenormals, 2) == v1n))
        % for ff = 1:size(mesh.f, 1)
        %     assert(abs(sum(v1t(ff, :) .* facenormals(ff, :))) < 1e-7)
        % end
        % Note this is different than: sum(v1f .* facenormals, 2) 

        %% RMS and difference
        v1frms = sqrt(sum(vecnorm(v1f, 2, 2).^2) / size(v1t, 1)) ;
        v2frms = sqrt(sum(vecnorm(v2f, 2, 2).^2) / size(v2t, 1)) ;
        v1trms = sqrt(sum(vecnorm(v1t, 2, 2).^2) / size(v1t, 1)) ;
        v2trms = sqrt(sum(vecnorm(v2t, 2, 2).^2) / size(v2t, 1)) ;
        dvtrms = sqrt(sum(vecnorm(dvt, 2, 2).^2) / size(dvt, 1)) ;

        %% Save measures
        disp(['saving relativeMotion to ' outfn])
        save(outfn, 'v1t', 'v2t', 'v1n', 'v2n', ...
            'v12d', 'v22d', 'dvtrms', ...
            'v1trms', 'v2trms','v1frms', 'v2frms', ... 
            'dvt', 'dvn', 'dv2d', 'dilation')
    else
        load(outfn, 'dv2d', 'v1t', 'v2t', 'dilation', ...
            'v1trms', 'v2trms', 'dvtrms') ;
        %  'v12d', 'v22d', 
    end
    
    % Collate
    dvtT(tidx, :, :) = dv2d ./ dilation(:) ;
    v1tT(tidx, :, :) = v1t ;
    v2tT(tidx, :, :) = v2t ;
    v1trmsV(tidx) = v1trms ;
    v2trmsV(tidx) = v2trms ;
    dvtrmsV(tidx) = dvtrms ;
    
    %% Check visually
    % bc = barycenter(mesh.v, mesh.f) ;
    % trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none')
    % hold on;
    % quiver3(bc(:, 1), bc(:, 2), bc(:, 3), v2t(:, 1), v2t(:, 2), v2t(:, 3), 0)
    % axis equal
    
    %% Plot 2d velocity difference
    imDir = fullfile(QS.dir.pivMultiChannel.relativeMotion, ...
        'images_v2d') ;
    if ~exist(imDir, 'dir')
        mkdir(imDir)
    end
    
    dv2dfn = fullfile(QS.dir.pivMultiChannel.relativeMotion, ...
        'images_v2d', [sprintf(QS.fileBase.name, tp) '.png']) ;
    
    if ~exist(dv2dfn, 'file') || overwriteImages
        disp('Creating figure of difference in velocity')
        im = max(im0, [], 3) ;
        imRGB = cat(3, im, im, im) ;
        opts = struct() ;
        opts.outfn = dv2dfn ;
        washout2d = 0.5 ;
        % Load the image to put flow on top
        im = imRGB * washout2d + max(im0(:)) * (1-washout2d) ;
        opts.label = ['scaled relative tangential motion, ', ...
            '$|\Delta v_t|/||g^{-1}||$ [$\mu$m/', QS.timeUnits, ']'] ;
        opts.xlim = [0, size(im, 2)] ;
        if doubleCovered
            opts.ylim = size(im, 1) * [0.25, 0.75] ;
        else
            opts.ylim = size(im, 1) * [0, 1] ;
        end
        opts.title = ['relative motion, $t=$', ...
            sprintf('%03d', tp-t0), ' ', QS.timeUnits] ;

        % Quiver opts
        opts.qsubsample = 5 ;
        opts.qscale = 5 ;
        opts.overlay_quiver = false ;
        % mesh.nU = QS.nU ;
        % mesh2d = cutRectilinearCylMesh(mesh) ;
        % mesh2d.v = mesh2d.u ;
        mesh.v = mesh.u ;
        vectorFieldHeatPhaseOnImage(im, mesh, ...
                    dv2d(:,1) ./ dilation(:), ...
                    dv2d(:, 2) ./ dilation(:), vtscale, opts)
    end
end

%% Average velocity over time
dvtAvgTime = mean(dvtT, 1) ;
dvtAvgSpace = squeeze(mean(dvtT, 2)) ;


%% Autocorrelation of relative motion
plotResult = false ;
lags = 10 ;
pearsonFace = zeros(nFaces, lags) ;
pearsonFaceComplex = zeros(nFaces, lags) ;
pearsonFaceNoMean = zeros(nFaces, lags) ;
for face = 1:size(dvtT, 2)
    if mod(face, 250) == 0
        disp(['computing pearson for face ' num2str(face)])
    end
    % Could treat as one complex velocity
    dvtc = dvtT(:, face, 1) + 1j * dvtT(:, face, 2) ;
    pearsonFaceComplex(face, :) = acf(dvtc, lags, true, true, plotResult) ;
    
    dvtc = [dvtT(:, face, 1); dvtT(:, face, 2)] ;
    pearsonFace(face, :) = acf(dvtc, lags, true, true, plotResult) ;
    
    % No mean subtraction
    dvtc = [dvtT(:, face, 1); dvtT(:, face, 2)] ;
    pearsonFaceNoMean(face, :) = acf(dvtc, lags, true, false, plotResult) ;
end

%% Save pearson correlations
pearsonfn = fullfile(QS.dir.pivMultiChannel.relativeMotion, ...
    'pearsonCorrCoeffFaces') ;
disp(['Saving pearsonFace to ' pearsonfn])
save([pearsonfn '.mat'], 'pearsonFace', 'pearsonFaceComplex', ...
    'pearsonFaceNoMean', 'dvtAvgTime', 'dvtAvgSpace', 'dvtrmsV')

%% mean pearson coefficient across space (faces) -- NO MEAN
pearsonfn = fullfile(QS.dir.pivMultiChannel.relativeMotion, 'pearson_corrcoeff_autocov') ;
for noMean = [true, false]
    if noMean
        pearsonAvg = mean(pearsonFaceNoMean, 1) ;
        pearsonFileName = [ pearsonfn '_noMean' ];
    else
        pearsonAvg = mean(pearsonFace, 1) ;
        pearsonFileName = pearsonfn ;
    end

    clf
    if ~exist([pearsonFileName '.pdf'], 'file') || overwriteImages
        % Plot rejection region lines for test of individual autocorrelations
        % H_0: rho(tau) = 0 at alpha=.05
        bar(pearsonAvg)
        H0 = (1.96)*(1/sqrt(ntps)) ;
        line([0 lags+.5], H0*ones(1,2))
        line([0 lags+.5], -H0*ones(1,2))

        % Some figure properties
        line_hi = H0 + .05;
        line_lo = -H0 - .05;
        bar_abs = max(abs(pearsonAvg))+.05 ;

        % limits
        xlim([0 lags+.60])
        if any(abs([line_hi line_lo]) > abs(bar_abs)) 
            % if rejection lines might not appear on graph
            ylim([line_lo line_hi])
        else
            ylim([-bar_abs bar_abs])
        end
        if noMean
            ylabel(['Pearson coefficient, ', ...
                '$\rho=', ...
                '\Pi (\Delta v_t)(\Delta v_{t+\tau}) / \sigma^2$'], ...
                'interpreter', 'latex')
        else
            ylabel(['Pearson coefficient, ', ...
                '$\rho=', ...
                '\Pi (\Delta v_t-\langle \Delta v_t \rangle)', ...
                '(\Delta v_{t+\tau} - \langle \Delta v_t \rangle)', ...
                '/ \sigma^2$'], ...
                'interpreter', 'latex')
        end
        xlabel(['lag time $\tau$ [' QS.timeUnits ']'], 'interpreter', 'latex')
        legend({'$\rho$', '$H_0$: $\rho=0$ at $\alpha=0.05$'}, ...
            'interpreter', 'latex')
        saveas(gcf, [pearsonFileName '.pdf'])
        saveas(gcf, [pearsonFileName '.png'])
    end
end

%% mean Complex pearson coefficient
pearsonComplexMag = sqrt(pearsonFaceComplex .* conj(pearsonFaceComplex)) ;
pearsonComplexAvg = mean(pearsonComplexMag, 1) ;
pearsonfn = fullfile(QS.dir.pivMultiChannel.relativeMotion, 'pearsonComplex_corrcoeff_autocov') ;

clf
if ~exist([pearsonfn  '_Complex.pdf'], 'file') || overwriteImages
    % Plot rejection region lines for test of individual autocorrelations
    % H_0: rho(tau) = 0 at alpha=.05
    bar(pearsonComplexAvg)
    line([0 lags+.5], (1.96)*(1/sqrt(ntps))*ones(1,2))

    % Some figure properties
    line_hi = (1.96)*(1/sqrt(ntps))+.05;
    % line_lo = -(1.96)*(1/sqrt(ntps))-.05;
    bar_hi = max(pearsonComplexAvg)+.05 ;
    % bar_lo = -max(pearsonAvg)-.05 ;

    if (abs(line_hi) > abs(bar_hi)) % if rejection lines might not appear on graph
        axis([0 lags+.60 0 line_hi])
    else
        axis([0 lags+.60 0 bar_hi])
    end
    ylabel(['Pearson coefficient, ', ...
        '$\rho=', ...
        '\left|\Pi \overline{(v_t - \langle v_t \rangle)}(v_{t+\tau} - \langle v_t \rangle)', ...
        '/ \sigma^2\right|$'], ...
        'interpreter', 'latex')
    xlabel(['lag time $\tau$ [' QS.timeUnits ']'], 'interpreter', 'latex')
    legend({'$\rho$', '$H_0$: $\rho=0$ at $\alpha=0.05$'}, ...
        'interpreter', 'latex')
    saveas(gcf, [pearsonfn  '_Complex.pdf'])
    saveas(gcf, [pearsonfn  '_Complex.png'])
end

%% Compare to random walk in confined well
% res = zeros(20, lags) ;
% for pp = 1:20
%     nsteps = 1e6 ;
%     k = 1.001 ;
%     kicks = normrnd(0, 0.1, nsteps, 1) ;
%     xx = zeros(nsteps, 1) ;
%     xx(1) = normrnd(0, 1) ;
%     for q = 1:nsteps-1
%         xx(q+1) = xx(q) + kicks(q) - k*xx(q)  ;
%     end
%     res(pp, :) = acf(xx, lags, true, false, true) ;
% end
% 
% NN = 20 * nsteps ;
% bar(mean(res, 1))
% line([0 lags+.5], (1.96)*(1/sqrt(NN))*ones(1,2))
% line([0 lags+.5], (-1.96)*(1/sqrt(NN))*ones(1,2))



%% Spatial map of pearson's
% initial (reference) mesh
clf
if ~exist([pearsonfn '.png'], 'file') || overwriteImages
    pearsonMapDir = fullfile(QS.dir.pivMultiChannel.relativeMotion, ...
            'pearsonSpatialMaps') ;
    if ~exist(pearsonMapDir, 'dir')
        mkdir(pearsonMapDir)
    end
    alpha0p05 = (1.96)*(1/sqrt(ntps))+.05 ;
    clf
    QS.setTime(t0) ; 
    mesh = QS.getCurrentSPCutMeshSmRS() ;
    mesh.u(:,1) = mesh.u(:, 1) / max(mesh.u(:, 1));
    for lag = [1, 1:lags]
        % vangle = mod(angle(pearsonFace(:, lag)), 2*pi) ;

        pcoef = pearsonComplexMag(:, lag) ;
        % speed = vecnorm([real(pearsonFace(:, lag)), imag(pearsonFace(:, lag))], 2, 2) ;
        % all(abs(speed - speed2) < 1e-13)
        opts = struct() ;
        h2 = patch( 'Faces', mesh.f, 'Vertices', mesh.u, ...
            'FaceVertexCData', pcoef, 'FaceColor', 'flat', ...
            'EdgeColor', 'none') ;
            % , 'FaceVertexAlphaData', 1, 'faceAlpha', 'flat') ;
        xlabel('ap position, $u$', 'interpreter', 'latex')
        ylabel('circumferential position, $v$', 'interpreter', 'latex')
        axis equal
        maxp = 0.01 * round(100* max(pcoef)) ;
        caxis([0, max(alpha0p05, maxp)])
        cb = colorbar() ;
        if maxp < alpha0p05
            cb.YTick = [0, maxp, alpha0p05] ;
        else
            cb.YTick = [0, alpha0p05, maxp];
            cb.YTickLabel = {'0', '$H_0^{\alpha=0.05}$', num2str(max(pcoef))} ;
            cb.TickLabelInterpreter = 'latex' ;
        end
        colormap(cividis)
        ylabel(cb, '$\rho$', 'interpreter', 'latex')
        title({'Pearson coefficient for relative motion', ...
            ['$\rho(t, t+$' num2str(lag) ' ' QS.timeUnits ')']}, ...
            'interpreter', 'latex')
        axis off

        pearsonfn = fullfile(pearsonMapDir, ...
            sprintf('pearson_corrcoeff_autocov_map%03d', lag)) ;
        saveas(gcf, [pearsonfn '.pdf'])
        saveas(gcf, [pearsonfn '.png'])
    end
end

%% Correlation coefficient between the two channels
cc = corrcoef(v1tT(:), v2tT(:)) ;
cc = cc(1, 2) ;

clf
% Time-dependent cross-correlation
ct = zeros(ntps, 1) ;
for tidx = 1:ntps 
    % tp = QS.xp.fileMeta.timePoints(tidx) ;
    tmp = corrcoef(v1tT(tidx, :), v2tT(tidx, :)) ;
    ct(tidx) = tmp(1, 2) ;
end
plot(QS.xp.fileMeta.timePoints(1:end-1), ct, '.-')
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
ylabel('cross correlation', 'interpreter', 'latex')
yyaxis right

plot(QS.xp.fileMeta.timePoints(1:end-1) - t0, v1trmsV, '-')
hold on;
plot(QS.xp.fileMeta.timePoints(1:end-1) - t0, v2trmsV, '.')
ylabel('$\sqrt{(v_t)^2}$', 'interpreter', 'latex')
legend({'cross-correlation', '$v_{\textrm{rms}}^A$', ...
    '$v_{\textrm{rms}}^B$'}, 'interpreter', 'latex')

figfn = fullfile(QS.dir.pivMultiChannel.relativeMotion, ...
            'cross_correlation.pdf') ;
saveas(gcf, figfn)

save(fullfile(QS.dir.pivMultiChannel.relativeMotion, ...
            'cross_correlation.mat'), 'cc', 'ct', 'v1trmsV', 'v2trmsV') ;


% scatter(v1tT(:), v2tT(:), 3, 'filled', 'markerfacealpha', 0.3)

