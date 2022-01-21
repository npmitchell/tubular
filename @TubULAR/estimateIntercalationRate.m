function estimateIntercalationRate(QS, options)
% estimateIntercalationRate(QS, options)
% TODO: use lagrangian pathlines of cell semgentation to estimate
% contribution from flow only
%
% Parameters
% ----------
% options: struct with fields
%   
%
% NPMitchell 2021


% Unpack options
timePoints = QS.xp.fileMeta.timePoints ;
maxCellSize = 200 ;  % maximum size allowed to include a cell
overwrite = false ;
overwriteImages = false ;
debug = false ;
[~, ~, ~, xyzlims] = QS.getXYZLims() ;
t0 = QS.t0set() ;

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 'xyzlims')
    xyzlims = options.xyzlims ;
end

% Directory preparation
imDir = fullfile(QS.dir.segmentation, 'seg3d', 'intercalation_images') ;
if ~exist(imDir, 'dir')
    mkdir(imDir)
end

%% Plot each lobe
features = QS.getFeatures() ;
folds = features.folds ;
try
    tidxs = QS.xp.tIdx(timePoints) ;
catch
    for qq = 1:length(timePoints)
        tidxs(qq) = QS.xp.tIdx(timePoints(qq)) ;
    end
end
nLobes = size(folds, 2) + 1 ;

% Colors
colors = define_colors() ;

% preallocate
cos2thetaM = cell(length(folds), 1) ;
sin2thetaM = cos2thetaM ;
aspectM = cos2thetaM ;
nbins = 10 ;
for qq = 1:length(folds) + 1
    cos2thetaM{qq} = zeros(nbins, length(tidxs)) ;
    sin2thetaM{qq} = zeros(nbins, length(tidxs)) ;
    aspectM{qq} = zeros(nbins, length(tidxs)) ;
end

% preallocate
ecell = cell(length(tidxs), 1) ;
c2t_piv = ecell ;
s2t_piv = ecell ;
ar_cell = ecell ;
theta_cell = ecell ;
areas_cell = ecell ;
cntrd_cell = ecell ;

% Average cell shape measurements
earr = zeros(length(tidxs), 1) ;
ar_mcell = earr ;
theta_mcell = earr ;
c2t_mcell = earr ;
s2t_mcell = earr ;

% Ensure that pullback pathlines are loaded
QS.getPullbackPathlines(t0, 'vertexPathlines', 'vertexPathlines3d') ;

%% Create synthetic pathlines for cells from t=0
% synethetic cell vertex pathlines in pullback space
if ~exist(fullfile(QS.dir.segmentation, 'pathlines'), 'dir')
    mkdir(fullfile(QS.dir.segmentation, 'pathlines'))
end
cellVertexPathlineFn = fullfile(QS.dir.segmentation, 'pathlines', ...
    sprintf('cellVertexPathlines_%06t0.mat', t0)) ;
if ~exist(cellVertexPathlineFn, 'file') || overwrite 
    QS.setTime(t0) ;
    seg2d = getCurrentSegmentation2D(QS) ;
    cellV0 = seg2d.seg2d.vdat.v ;
    opts = struct('preview', true) ;
    [segVertexPathlines2D, segVertexPathlines3D] = ...
        QS.samplePullbackPathlines(cellV0, opts) ;
    save(cellVertexPathlineFn, 'segVertexPathlines2D', ...
        'segVertexPathlines3D')
else
    load(cellVertexPathlineFn, 'segVertexPathlines2D', ...
        'segVertexPathlines3D')
end

meanQLobeAspects_cell = zeros(nLobes, length(timePoints)) ;
meanQLobeAspectStds_cell = zeros(nLobes, length(timePoints)) ;
meanQLobeThetas_cell = zeros(nLobes, length(timePoints))  ;
meanQLobeAspects_piv = zeros(nLobes, length(timePoints)) ;
meanQLobeAspectStds_piv = zeros(nLobes, length(timePoints)) ;
meanQLobeThetas_piv = zeros(nLobes, length(timePoints))  ;
dmy = 1 ;
for tidx = tidxs 
    tp = timePoints(dmy) ;
    disp(['t= ', num2str(tp)])
    QS.setTime(tp) ;
    nU = QS.nU ;
    foldt = double(folds(tidx, :)) / double(nU) ;
    nLobes = length(foldt(:)) + 1 ;
    foldt = [foldt(:); 1] ;
    
    % Load segs
    seg2d = QS.getCurrentSegmentation2D() ;
    coordSys = seg2d.coordSys ;
    segIm = seg2d.segIm ;
    seg2d = seg2d.seg2d ;
    Lx = size(segIm, 2) ;
    seg3d = QS.getCurrentSegmentation3D() ;
    % ensure 3d and 2d segmentation are based on the same coordsys
    assert(strcmpi(seg3d.coordSys, coordSys))
    seg3d = seg3d.seg3d ;
    
    % Compare mean ellipse minus reference ellipse to strain 
    stat = seg3d.statistics ;
    keep = stat.keep ;
    ar_cell{dmy} = sqrt(seg3d.qualities.moment2(keep) ./ seg3d.qualities.moment1(keep)) ;
    theta_cell{dmy} = seg3d.qualities.ang1(keep) ;
    areas_cell{dmy} = seg3d.qualities.areas(keep) ;
    cntrd_cell{dmy} = seg3d.cdat.centroids_uv(keep, :) ;
    moi_cell{dmy} = seg3d.qualities.mInertia(keep, :, :) ;
    QQ_cell{dmy} = seg3d.qualities.nematicTensor(keep, :, :) ;
    meanQLobeAspects_cell(:, dmy) = seg3d.statistics.lobes.meanQLobeAspect ;
    meanQLobeThetas_cell(:, dmy) = seg3d.statistics.lobes.meanQLobeTheta ;
    meanQLobeAspectStds_cell(:, dmy) = seg3d.statistics.lobes.meanQLobeAspectStd ;
    
    % Average cell shape measurements
    % ar_mcell(dmy) = sqrt(stat.meanMoment2WeightBounded ./ stat.meanMoment1WeightBounded) ;
    % theta_mcell(dmy) = mod(stat.meanAng1WeightBounded, 2*pi) ;
    % c2t_mcell(dmy) = cos(2 * stat.meanAng1WeightBounded) ;
    % s2t_mcell(dmy) = sin(2 * stat.meanAng1WeightBounded) ;
    
    % Average nematic shape: 
    % strength is = norm(meanQ) * 2 + 1 =sqrt(abs(stat.meanQMoment2)) + 1
    ar_mcell(dmy) = (norm(stat.meanQWeighted) * 2) + 1 ;
    theta_mcell(dmy) = stat.meanQThetaWeighted ;
    c2t_mcell(dmy) = cos(2 * stat.meanQThetaWeighted) ;
    s2t_mcell(dmy) = sin(2 * stat.meanQThetaWeighted) ;
    
    % seg3d.apBins
    % seg3d.apCos2Theta
    % seg3d.apSin2Theta
    
    % Sample strain at faces of advected triangulation at cell barycenters    
    % Load strain
    strain = getCurrentPathlineStrain(QS, t0, 'beltrami') ;
    % tre = strain.strain.tre ;
    % dev = strain.strain.dev ;
    % thet = strain.strain.theta ;
    muv = strain.beltrami.mu_material_filtered ;
    
    % Evaluated at cell barycenters
    % Load advected lagrangian triangulation from pathlines
    vX = QS.pathlines.vertices.vX(tidx, :) ;
    vY = QS.pathlines.vertices.vY(tidx, :) ;
    vx = QS.pathlines.vertices3d.vXrs(tidx, :) ;
    vy = QS.pathlines.vertices3d.vYrs(tidx, :) ;
    vz = QS.pathlines.vertices3d.vZrs(tidx, :) ;
    v3d = [vx(:), vy(:), vz(:)] ;
    try
        assert(QS.pathlines.vertices.Lx(tidx) == Lx)
    catch
        error(['Pullback images for pathlines and segmentation are ' ...
            'different size. Could normalize each here, but are you ' ...
            'sure you are using the same images?'])
    end
    % Check that the mesh matches the QS mesh params
    assert(QS.pathlines.refMesh.nU == nU) ;
    nV = QS.nV ;
    assert(QS.pathlines.refMesh.nV == nV) ;
    pathlineMesh = struct('f', QS.pathlines.refMesh.f, ...
        'u', [vX(:), vY(:)], ...
        'v', v3d, ...
        'vn', [real(muv(:)), imag(muv(:)), 0*real(muv(:))], ...
        'pathPairs', [1:nU; nU*(nV-1)+1:nU*nV]') ;
    [TF, TU, TV3D, muV] = tileAnnularCutMesh(pathlineMesh, [-1, 1]) ;
    
    tri = triangulation(TF, TU(:, 1), TU(:, 2));
    fieldfaces = pointLocation(tri, seg2d.cdat.centroid);
    
    % Now sample mu values
    [V2F, ~] = meshAveragingOperators(TF, TV3D) ;
    muF = V2F * muV ;
    try
        muff = muF(fieldfaces(keep), 1) + 1j * muF(fieldfaces(keep), 2) ;
    catch
        muR = scatteredInterpolant(TU(:, 1), TU(:, 2), muV(:, 1), ...
            'linear', 'nearest') ;
        muI = scatteredInterpolant(TU(:, 1), TU(:, 2), muV(:, 2), ...
            'linear', 'nearest') ;
        muReval = muR(seg2d.cdat.centroid(keep, 1), ...
            seg2d.cdat.centroid(keep, 2)) ;
        muIeval = muI(seg2d.cdat.centroid(keep, 1), ...
            seg2d.cdat.centroid(keep, 2)) ;
        muff = muReval + 1j * muIeval ;
    end
    
    % check centroid sampling of the tiled mesh pullback
    if debug
        clf
        trisurf(TF, TU(:, 1), TU(:, 2), 0*TU(:, 1), muV(:, 1), ...
            'edgeColor', 'none', 'FaceColor', 'flat') ;
        hold on;
        scatter(seg2d.cdat.centroid(keep, 1), ...
            seg2d.cdat.centroid(keep, 2), 10, muff(:, 1))
        axis equal
        caxis([-0.5, 0.5])
        colormap(blueblackred)
        view(2)
        waitfor(gcf)
     end
    
    % mean of sampled strain faces
    if tp ~= t0
        mu_piv{dmy} = muff ;
        ar_piv{dmy} = (1 + abs(muff)) ./ (1 - abs(muff)) ;
        theta_piv{dmy} = 0.5 * angle(muff) ;
        c2t_piv{dmy} = cos(2 * theta_piv{dmy}) ;
        s2t_piv{dmy} = sin(2 * theta_piv{dmy}) ;
        
        % Weighted mean and unc
        weights = areas_cell{dmy} / sum(areas_cell{dmy}) ;
        
        % ap_x = dev .* cos(2*thet) ;
        % ap_y = dev .* sin(2*thet) ;
        % ap_xy = [mean(ap_x), mean(ap_y)] ;
        % theta_strain = 0.5 * mod(atan2(ap_xy(:, 2), ap_xy(:, 1)), 2*pi) ;
        % c2t_strain(dmy) = cos(2 * theta_strain) ;
        % s2t_strain(dmy) = cos(2 * theta_strain) ;
        % dev_strain(dmy) = sqrt(ap_xy(:, 1).^2 + ap_xy(:, 1).^2) ;
        % tr_strain(dmy) = mean(tre) ;

        %% LOBES: Obtain the strain rates in each chamber at centroid locations
        % from PIV
        nU = QS.nU ;
        fold0 = double(folds(tidx, :)) / double(nU) ;
        nLobes = length(fold0(:)) + 1 ;
        foldt = [0; fold0(:); 1] ;

        % Weights for means
        weights = areas_cell{dmy} / sum(areas_cell{dmy}) ;
        ap_pos = seg3d.cdat.centroids_uv(keep, 1) ;

        %% In terms of beltrami coefficients
        % [~, lobes_muff_re, lobes_std_muff_re] = binDataMeanStdWeighted(ap_pos, ...
        %     real(muff), foldt, weights) ;
        % [~, lobes_muff_im, lobes_std_muff_im] = binDataMeanStdWeighted(ap_pos, ...
        %     imag(muff), foldt, weights) ;
        % lobes_muff = lobes_muff_re + 1j * lobes_muff_im ;
        % lobes_std_muff = lobes_std_muff_re + 1j * lobes_std_muff_im ;

        %% In terms of nematic shape tensor
        strengths = ar_piv{dmy} - 1 ;
        nn = [cos(theta_piv{dmy}), sin(theta_piv{dmy})] ;
        [meanQ_lobes{dmy}, stdQ_lobes{dmy}, meanQLobeStrength, meanQLobeTheta, ...
            meanQLobeAspectStd] = ...
            binDataMeanStdWeightedNematic2D(nn, strengths, ap_pos, foldt, weights) ;
        meanQLobeAspect = meanQLobeStrength + 1 ;
    else
        mu_piv{dmy} = muff ;
        ar_piv{dmy} = (1 + abs(muff)) ./ (1 - abs(muff)) ;
        theta_piv{dmy} = NaN ;
        c2t_piv{dmy} = NaN ;
        s2t_piv{dmy} = NaN ;
        % dev_strain(dmy) = NaN ;
        % tr_strain(dmy) = NaN ;
        meanQLobeAspect = ones(nLobes, 1) ;
        meanQLobeAspectStd = zeros(nLobes, 1) ;
        meanQLobeTheta = nan(nLobes, 1) ;
    end
        
    % Store for later plotting
    meanQLobeAspects_piv(:, dmy) = meanQLobeAspect ;
    meanQLobeAspectStds_piv(:, dmy) = meanQLobeAspectStd ;
    meanQLobeThetas_piv(:, dmy) = meanQLobeTheta ;
    
    % Update dummy index for next timepoint
    dmy = dmy + 1 ;
end

%% Convert to relative shear compared to t0
dmy = 1 ;
dmy0 = find(timePoints == t0) ;

% Use average cell shape measurements to find simple shear that transforms
% average cell at t0 into unit circle
ar0 = ar_mcell(dmy0) ;
theta0 = mod(theta_mcell(dmy0), pi) ;
c2t0 = c2t_mcell(dmy0) ;
s2t0 = s2t_mcell(dmy0) ;
% matrix taking nematic to isotropic
% rot1 = [cos(theta0), -sin(theta0); sin(theta0), cos(theta0)] ;
% rot2 = [cos(theta0), sin(theta0); -sin(theta0), cos(theta0)] ;
% pure_shear = [1/sqrt(ar0), 0; 0, sqrt(ar0)] ; 
% pure_shear2 = [1/ar0, 0; 0, ar0] ; 
% RR = rot1 * pure_shear * rot2 ; % For application on tensor with dims of length
% RM = rot1 * pure_shear2 * rot2 ; % For vector ar*[cos(theta), sin(theta)]
if abs(theta0/pi - 0.5) < 1e-2
    pure_shear = [1/ar0, 0; 0, ar0] ;
else
    error('handle case here')
end

% Get pure shear for each lobe
for lobe = 1:nLobes
    ar0_lobe = meanQLobeAspects_cell(lobe, dmy0) ;
    pure_shear_lobe{lobe} = [1/ar0_lobe, 0; 0, ar0_lobe] ;
end

% We can convert to beltrami coefficient like this
dmy = 1 ;
for tidx = tidxs
    % Apply shear/rotation to unity
    nCells = length(ar_cell{dmy}) ;
    
    % Simple 1d picture here
    if sign(c2t_mcell(dmy)) == sign(c2t0)
        ar_mcell0_1d(dmy) = ar0 / ar_mcell(dmy) ;
    else
        ar_mcell0_1d(dmy) = ar_mcell(dmy) * ar0 ;
    end
    dmy = dmy + 1 ;
end

dmy = 1 ;
ar_cell0 = cell(length(tidxs), 1) ;
moment1_cell0 = cell(length(tidxs), 1) ;
moment2_cell0 = cell(length(tidxs), 1) ;
theta_cell0 = cell(length(tidxs), 1) ;
mu_cell0 = cell(length(tidxs), 1) ;
Q0 = cell(length(tidxs), 1) ;
Q0_lobe0 = cell(length(tidxs), 1) ;
meanQ0 = cell(length(tidxs), 1) ;
ar_mcell0 = zeros(length(tidxs), 1) ;
theta_mcell0 = zeros(length(tidxs), 1) ;
meanQLobeAspects_cell0 = zeros(nLobes, length(tidxs)) ;
meanQLobeAspectStds_cell0 = zeros(nLobes, length(tidxs)) ;
meanQLobeThetas_cell0 = zeros(nLobes, length(tidxs)) ;
for tidx = tidxs
    QS.setTime(timePoints(dmy)) ;
    seg3d = QS.getCurrentSegmentation3D() ;
    seg3d = seg3d.seg3d ;
    nCells = length(ar_cell{dmy}) ;
    keep = seg3d.statistics.keep ;
    ap_pos = seg3d.cdat.centroids_uv(keep, 1) ;
    
    % Apply shear to each cell here
    ar_cell0{dmy} = zeros(nCells, 1) ;
    moment1_cell0{dmy} = zeros(nCells, 1) ;
    moment2_cell0{dmy} = zeros(nCells, 1) ;
    theta_cell0{dmy} = zeros(nCells, 1) ;
    Q0{dmy} = zeros(nCells, 2, 2) ;
    for cc = 1:length(ar_cell{dmy})
        % ang = theta_cell{dmy}(cc) ;
        % vec2head(cc, :) = RM * ar_cell{dmy}(cc) * [cos(ang), sin(ang)]' ;
        
        moi = squeeze(moi_cell{dmy}(cc, :, :)) ;
        moi0 = pure_shear * moi ;
        [ eig_vec, eig_val ] = eig(moi0);
        bigMoment = max([eig_val(1,1), eig_val(2,2)]) ;
        smallMoment = min([eig_val(1,1), eig_val(2,2)]) ;
        moment1_cell0{dmy}(cc) = smallMoment ;
        moment2_cell0{dmy}(cc) = bigMoment ;
        theta_cell0{dmy}(cc) = atan2( eig_vec(2,1), eig_vec(1,1) );
        
        % nematic tensor 
        mratio = bigMoment / smallMoment ;
        tt = mod(theta_cell0{dmy}(cc), pi) ;
        nn = [cos(tt); sin(tt)] ;
        % Create traceless symmetric matrix using unit vec
        strength = abs(sqrt(mratio)) - 1 ;
        Q0{dmy}(cc, :, :) = strength * (nn' * nn - 0.5 * [1, 0; 0, 1]) ;
    end
    % aspect ratio of all cells
    ar_cell0{dmy} = sqrt(moment2_cell0{dmy} ./ moment1_cell0{dmy}) ;
    
    %% aspect of mean cell In terms of nematic shape tensor
    weights = seg3d.qualities.areas(seg3d.statistics.keep) ;
    weights = weights ./ sum(weights) ;
    meanQ0{dmy} = squeeze(sum(weights .* Q0{dmy}, 1) / sum(weights)) ;
    [ eig_vec, eig_val ] = eig(meanQ0{dmy}) ;
    meanQMoment1 = eig_val(1,1);
    meanQMoment2 = eig_val(2,2);
    % assert(norm(meanQ0) == abs(meanQMoment1))
    theta_mcell0(dmy) = atan2( eig_vec(2,2), eig_vec(1,2) );
    c2t_mcell0(dmy) = cos(2 * theta_mcell0(dmy)) ;
    s2t_mcell0(dmy) = sin(2 * theta_mcell0(dmy)) ;
    meanQ0{dmy} = squeeze(sum(Q0{dmy}, 1)) ;
    ar_mcell0(dmy) = norm(meanQ0{dmy})*2 + 1 ;
    
    %% Lobe classification
    [~, ~, loc] = histcounts(ap_pos, foldt);
    Q0_lobe0{dmy} = zeros(nCells, 2, 2) ;
    for lobe = 1:nLobes
        cellsInLobe = find(loc == lobe) ;
        for cidx = 1:length(cellsInLobe)
            cc = cellsInLobe(cidx) ;
            % ang = theta_cell{dmy}(cc) ;
            % vec2head(cc, :) = RM * ar_cell{dmy}(cc) * [cos(ang), sin(ang)]' ;

            moi = squeeze(moi_cell{dmy}(cc, :, :)) ;
            moi0 = pure_shear_lobe{lobe} * moi ;
            [ eig_vec, eig_val ] = eig(moi0);
            bigMoment = max([eig_val(1,1), eig_val(2,2)]) ;
            smallMoment = min([eig_val(1,1), eig_val(2,2)]) ;
            moment1_cell0_lobe0{dmy}(cc) = smallMoment ;
            moment2_cell0_lobe0{dmy}(cc) = bigMoment ;
            theta_cell0_lobe0{dmy}(cc) = atan2( eig_vec(2,1), eig_vec(1,1) );

            % nematic tensor 
            mratio = bigMoment / smallMoment ;
            tt = mod(theta_cell0_lobe0{dmy}(cc), pi) ;
            nn = [cos(tt), sin(tt)] ;
            % Create traceless symmetric matrix using unit vec
            strength = abs(sqrt(mratio)) - 1 ;
            Q0_lobe0{dmy}(cc, :, :) = strength * (nn' * nn - 0.5 * [1, 0; 0, 1]) ;
        end
    end
    
    % aspect ratio of all cells
    ar_cell0_lobe0{dmy} = sqrt(moment1_cell0_lobe0{dmy} ./ moment1_cell0_lobe0{dmy}) ;
    [meanQ0_lobe0{dmy}, stdQ0_lobe0{dmy}, meanQLobeStrength0, meanQLobeTheta0, ...
        meanQLobeAspectStd0] = ...
        binDataMeanStdWeightedNematic2D(Q0_lobe0{dmy}, ...
        ones(size(Q0_lobe0{dmy}, 1), 1), ap_pos, foldt, weights) ;
    meanQLobeAspects_cell0(:, dmy) = meanQLobeStrength0 + 1 ;
    meanQLobeAspectStds_cell0(:, dmy) = meanQLobeAspectStd0 ;
    meanQLobeThetas_cell0(:, dmy) = meanQLobeTheta0 ;
    
    %% Beltrami casting -- true cells
    % % kappa = aspect ratio = (1+|mu|) / (1-|mu|), 
    % % |mu| = (kappa - 1) / (1 + kappa)
    % abs_mu_cell = (ar_cell{dmy} - 1) ./ (ar_cell{dmy} + 1) ;
    % mu_cell{dmy} = exp(1j * 2 * theta_cell{dmy}) .* abs_mu_cell ;    
    
    %% Beltrami casting -- reference cells
    % % kappa = aspect ratio = (1+|mu|) / (1-|mu|), 
    % % |mu| = (kappa - 1) / (1 + kappa)
    abs_mu_cell = (ar_cell0{dmy} - 1) ./ (ar_cell0{dmy} + 1) ;
    mu_cell0{dmy} = exp(1j * 2 * theta_cell0{dmy}) .* abs_mu_cell ;    
        
    dmy = dmy + 1 ;
end

%% Compare the two shear deformations (cells versus tissue/PIV)
dmy = 1 ;
for tidx = tidxs 
    % Timing
    tp = timePoints(dmy) ;
    disp(['t= ', num2str(tp)])
    QS.setTime(tp) ;
    
    % cells
    ac2t_c = ar_cell0{dmy} .* cos(2*theta_cell0{dmy}) ;
    as2t_c = ar_cell0{dmy} .* sin(2*theta_cell0{dmy}) ;
    % areas_cell{dmy} ;
    % cntrd_cell{dmy}  ;
    
    % Plot cell against piv measurements --> ac2t_p = aspect*cos(2theta)_piv
    ac2t_p = ar_piv{dmy} .* c2t_piv{dmy} ;
    as2t_p = ar_piv{dmy} .* s2t_piv{dmy} ;
    
    % subplot(2, 1, 1)
    clf
    plot(ac2t_p, ac2t_c, '.')    
    % axis equal
    % ylim([-2, 2])
    % xlabel('piv $a\cos 2\theta$', 'interpreter', 'latex')
    % ylabel('segmentation $a\cos 2\theta$', 'interpreter', 'latex')
    % subplot(2, 1, 2)
    hold on;
    plot(as2t_p, as2t_c, '.') 
    axis equal
    xlabel('piv $a\cos 2\theta$ or $a\sin 2\theta$', 'interpreter', 'latex')
    ylabel('segmentation $a\cos 2\theta$ or $a\sin 2\theta$', 'interpreter', 'latex')
    titlestr = ['$t=$ ' num2str(tp - t0) ' ' QS.timeUnits] ;
    title(titlestr, 'interpreter', 'latex')
    saveas(gcf, fullfile(imDir, sprintf('piv_seg_%06d.png', tp)))
    
    % Beltrami scatter
    clf
    plot(real(mu_cell0{dmy}), real(mu_piv{dmy}), '.') ;
    hold on
    plot(imag(mu_cell0{dmy}), imag(mu_piv{dmy}), '.') ;
    legend({'$\Re\mu$', '$\Im\mu$'}, 'interpreter', 'latex')
    axis equal
    xlabel('piv $\mu$', 'interpreter', 'latex')
    ylabel('segmentation $\mu$', 'interpreter', 'latex')
    titlestr = ['$t=$ ' num2str(tp - t0) ' ' QS.timeUnits] ;
    title(titlestr, 'interpreter', 'latex')
    saveas(gcf, fullfile(imDir, sprintf('mus_%06d.png', tp)))
    
    % Average cell shape measurements
    seg3d = QS.getCurrentSegmentation3D() ;
    seg3d = seg3d.seg3d ;
    
    % Compare mean ellipse minus reference ellipse to strain 
    % error('use newly computed MOIs with subtracted ellipse')
    % stat = seg3d.statistics ;
    % ar_mcell(dmy) = sqrt(stat.meanMoment2 ./ stat.meanMoment1) ;
    % theta_mcell(dmy) = mod(stat.meanAng1, 2*pi) ;
    % c2t_mcell(dmy) = cos(2 * stat.meanAng1) ;
    % s2t_mcell(dmy) = sin(2 * stat.meanAng1) ;
    
    ar_mpiv(dmy) = (1+abs(mean(mu_piv{dmy}))) / (1-abs(mean(mu_piv{dmy}))) ;
    
    if tp == t0
        theta_mpiv(dmy) = NaN ;
        c2t_mpiv(dmy) = NaN ;
        s2t_mpiv(dmy) = NaN ;
    else
        theta_mpiv(dmy) = mod(0.5 * angle(mean(mu_piv{dmy})), 2*pi) ;
        c2t_mpiv(dmy) = cos(2 * theta_mpiv(dmy)) ;
        s2t_mpiv(dmy) = sin(2 * theta_mpiv(dmy)) ;
    end
    
    dmy = dmy + 1 ;
end


%% Plot nematic strength and direction for each lobe -- PIV
imfn = fullfile(QS.dir.segmentation, 'seg3d', 'piv_anisotropy_lobes.png') ;
clf
for lobe = 1:nLobes
    subplot(ceil(nLobes * 0.5), 2, lobe)
    midline = squeeze(meanQLobeAspects_piv(lobe, :)) ;
    uncs = squeeze(meanQLobeAspectStds_piv(lobe, :)) ;
    timestamps = timePoints - t0 ;
    if contains(QS.timeUnits, 'min')
        timestamps = timestamps / 60 ;
    end
    x2 = [timestamps, fliplr(timestamps)] ;
    fill(x2,[midline-uncs, fliplr(midline+uncs)], ...
        colors(1, :), 'facealpha', 0.3, 'edgecolor', 'none');
    hold on;
    plot(timestamps, midline, '.-', 'color', colors(1, :))
    ylim([1, Inf])
    
    yyaxis right
    % fill(x2, [c2t_low, fliplr(c2t_high)], redcol, 'facealpha', 0.3, 'edgecolor', 'none');
    hold on;
    plot(timestamps, mod(meanQLobeThetas_piv(lobe, :), pi)/pi, '.-')
    % plot(timestamps, sin(2*meanQLobeThetas_piv(lobe, :)), '.-')
    % 'color', redcol)
    ylim([0, 1])
    title(['chamber ' num2str(lobe)], 'interpreter', 'latex')
    
    if mod(lobe, 2) == 1 
        yyaxis left
        ylabel('aspect ratio, $a=2||Q|| + 1$',   'interpreter', 'latex')
    else
        yyaxis right
        ylabel('nematic orientation $\theta/\pi$',   'interpreter', 'latex')
    end
    
    % Time label
    if lobe > nLobes - 2
        if contains(QS.timeUnits, 'min')
            xlabel('time [hr]', 'interpreter', 'latex')  
        else
            xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
        end
    end
end
sgtitle('endoderm orientation over time', 'interpreter', 'latex')
saveas(gcf, imfn)


%% Plot nematic strength and direction for each lobe as ONE PLOT
imfn = fullfile(QS.dir.segmentation, 'seg3d', ...
    'piv_anisotropy_lobes_signed.png') ;
clf
for lobe = 1:nLobes
    % subplot(ceil(nLobes * 0.5), 2, lobe)
    c2t = cos(2*meanQLobeThetas_piv(lobe, :))  ;
    c2t = fillmissing(c2t, 'linear') ;
    midline = c2t .* (squeeze(meanQLobeAspects_piv(lobe, :)) - 1) ;
    uncs = c2t .* (squeeze(meanQLobeAspectStds_piv(lobe, :)) - 1) ;
    timestamps = timePoints - t0 ;
    if contains(QS.timeUnits, 'min')
        timestamps = timestamps / 60 ;
    end
    x2 = [timestamps, fliplr(timestamps)] ;
    fill(x2,[midline-abs(uncs), fliplr(midline+abs(uncs))], ...
        colors(lobe, :), 'facealpha', 0.3, 'edgecolor', 'none', ...
        'HandleVisibility', 'off');
    hold on;
    hs{lobe} = plot(timestamps, midline, '.-', 'color', colors(lobe, :)) ;
    
    legendentries{lobe} = ['chamber ' num2str(lobe)] ;
end
% ylims = ylim() ;
% ylim([-max(abs(ylims)), max(abs(ylims))])

% Mark zero line
plot(timestamps, 0*timestamps, 'k--', 'HandleVisibility','off')
% Labels
legend(legendentries, 'interpreter', 'latex', 'location', 'northwest')
ylabel('elongation, $2||Q|| \cos 2\theta$',   'interpreter', 'latex')
if contains(QS.timeUnits, 'min')
    xlabel('time [hr]', 'interpreter', 'latex')  
else
    xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')  
end
sgtitle('endoderm orientation over time', 'interpreter', 'latex')
saveas(gcf, imfn)
saveas(gcf, [imfn(1:end-3), 'pdf'])


%% Plot the shearing by PIV -- global average
clf 
scatter(timePoints - t0, ar_mpiv, 100, timePoints, 'filled')
ylabel('sheared aspect ratio from piv', 'interpreter', 'latex')
xlabel(['time [' QS.timeUnits ']'], 'interpreter', 'latex')
saveas(gcf, fullfile(QS.dir.segmentation, 'seg3d', ...
    'piv_aspect_ratio.png'))

%% COMPARE PIV TO CELL SHAPES
clf
% plot(ar_mpiv, ar_mcell0, '.-')
if contains(QS.timeUnits, 'min')
    timestamps = (timePoints - t0)/60 ;
end
scatter(ar_mpiv, ar_mcell0_1d, 100, timestamps, 'filled')
% TODO: ADD uncertainty here
xlabel('shear ratio from piv', 'interpreter', 'latex')
ylabel('cell shear ratio with respect to $t=0$', 'interpreter', 'latex')
hold on;
plot([1, 3], [1, 3], 'k--')
cb = colorbar() ;
colormap(magma)
if contains(QS.timeUnits, 'min')
    ylabel(cb, 'time [hrs]')
end
saveas(gcf, fullfile(QS.dir.segmentation, 'seg3d', ...
    'aspect_ratio_comparison.png'))
    

%% COMPARE LOBE BY LOBE
maxAR = -Inf ;
for raw = [true false]
    for shear_applied = [true false]
        if shear_applied
            shearStr = '_sheared' ;
        else
            shearStr = ''; 
        end
        for timeLimit = [76, Inf]
            clf
            tidx0 = find(timePoints == t0) ;
            if isempty(tidx0)
                error('please include t0 in the list of timePoints for now')
            end
            for lobe = 1:nLobes

                keepTime = find(timePoints - t0 < timeLimit) ;

                pivAR = fillmissing(squeeze(meanQLobeAspects_piv(lobe, keepTime)-1), 'linear') ;
                pivUnc = fillmissing(squeeze(meanQLobeAspectStds_piv(lobe, keepTime)), 'linear') ;
                pivTheta = fillmissing(squeeze(meanQLobeThetas_piv(lobe, keepTime)), 'linear') ;

                if shear_applied
                    cellAR = fillmissing(squeeze(meanQLobeAspects_cell0(lobe, keepTime)-1), 'linear') ;
                    cellTheta = fillmissing(squeeze(meanQLobeThetas_cell0(lobe, keepTime)), 'linear') ;
                    cellUnc = fillmissing(squeeze(meanQLobeAspectStds_cell0(lobe, keepTime)), 'linear');
                else
                    cellAR = fillmissing(squeeze(meanQLobeAspects_cell(lobe, keepTime)-1), 'linear') ;
                    cellTheta = fillmissing(squeeze(meanQLobeThetas_cell(lobe, keepTime)), 'linear') ;
                    cellUnc = fillmissing(squeeze(meanQLobeAspectStds_cell(lobe, keepTime)), 'linear');
                end

                % Cell elongation and piv elongation (subtracting off t0 for cells)
                cellE = cellAR .* cos(2*cellTheta) ;
                if ~raw
                    cellE = cellE - cellE(tidx0) ;
                end
                pivE = pivAR .* cos(2*pivTheta) ;

                % % ERRORBAR METHOD FOR VISUALIZATION
                % % Could include errors in both shear and cells
                h = errorbar(pivE, cellE, ...
                    cellUnc, cellUnc, 'color', colors(lobe, :), ...
                    'CapSize', 0) ; % , 'LineStyle', 'none') ;
                % h = errorbar(pivE, cellE, ...
                %     cellUnc, cellUnc, pivUnc, pivUnc, 'color', colors(lobe, :), ...
                %     'CapSize', 0, 'LineStyle', 'none') ;
                % Set transparency (undocumented)
                alpha = 0.3 ;
                set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', ...
                    'ColorData', [h.Line.ColorData(1:3); 255*alpha])
                hold on;
                scatter(pivE, cellE, 20, ...
                    'filled', 'markerfacecolor', colors(lobe, :), ...
                    'markeredgecolor', colors(lobe, :), 'handlevisibility', 'off')


                % % FILL METHOD
                % fill([pivE, fliplr(pivE)], ...
                %     [cellE-abs(cellUnc), fliplr(cellE+abs(cellUnc))], ...
                %      colors(lobe, :), 'facealpha', 0.1, 'edgecolor', 'none',...
                %     'handlevisibility', 'off') ;
                % hold on;
                % plot(pivE, cellE, '.-', 'color', colors(lobe, :), ...
                %     'markersize', 20)

                maxAR = max(maxAR, min(max(cellE), max(pivE))) ;
                legendEntries{lobe} = ['chamber ' num2str(lobe)] ;
            end
            axis equal
            axis square
            legend(legendEntries, 'interpreter', 'latex', 'location', 'southeast')
            plot([0, maxAR], [0, maxAR], 'k--', 'handleVisibility', 'off')
            xlabel('shear aspect ratio from piv', 'interpreter', 'latex')

            if raw
                rawStr = '_raw' ;
                ylabel('cell elongation, $E=2||Q|| \cos 2\theta$', 'interpreter', 'latex')
            else
                rawStr = '_modE0' ;
                ylabel('cell elongation, $\tilde{E}=2||Q|| \cos 2\theta - E_0$', 'interpreter', 'latex')
            end
            if timeLimit < Inf
                timeLimitStr = ['_timeLimit' num2str(timeLimit)] ;
            else
                timeLimitStr = '' ;
            end
            disp(['saving fig: ' shearStr rawStr timeLimitStr])
            saveas(gcf, fullfile(QS.dir.segmentation, 'seg3d', ...
                ['aspect_orientation_comparison_lobes' shearStr rawStr timeLimitStr '.png']))
            saveas(gcf, fullfile(QS.dir.segmentation, 'seg3d', ...
                ['aspect_orientation_comparison_lobes' shearStr rawStr timeLimitStr '.pdf']))
        end
    end
end

%% Extra direction info
clf
subplot(1, 2, 1)
scatter(c2t_mpiv, c2t_mcell0, 100, timePoints, 'filled')
axis equal
hold ;
xlim([-1.1, 1.1])
ylim([-1.1, 1.1])
ylabel('$\cos 2 \theta$ relative to initial cell shape', 'interpreter', 'latex')
xlabel('$\cos 2 \theta$ from piv', 'interpreter', 'latex')
subplot(1, 2, 2)
scatter(s2t_mpiv, s2t_mcell0, 100, timePoints, 'filled')
axis equal
xlim([-1.1, 1.1])
ylim([-1.1, 1.1])
ylabel('$\sin 2 \theta$ relative to initial cell shape', 'interpreter', 'latex')
xlabel('$\sin 2 \theta$ from piv', 'interpreter', 'latex')
saveas(gcf, fullfile(QS.dir.segmentation, 'seg3d', ...
    'aspect_orientation_comparison.png'))

