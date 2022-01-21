function [rawRicciMesh, ricciMu] = generateRawRicciMeshTimePoint(QS, tp, options)
% generateRawRicciMeshTimePoint(QS, tp, options)
% 
% Compute solution to Ricci flow for a given timepoint's cylinderMesh. 
% 
% See also
% --------
%  generateRicciMeshTimePoint(QS, tp, options) --> for Ricci mesh of
%   gridded coordinate meshes after re-parameterization of the mesh into 
%   spcutMesh or uv coords etc
% 
% Parameters
% ----------
% QS : quapSlap class instance
% tp : int (default=QS.t0)
%   timestamp of timepoint for which to generate Ricci flow mesh
% options : optional struct with optional fields
%   maxIter : int (default=100)
%   radiusTolerance : float (default=0.01)
%       maximum allowed fractional deviation of Ricci flow solution inner
%       and outer radius from a true circle with fixed radius (variable
%       inner radius, fixed outer radius of 1)
%   save_ims : bool (default = true)
%       save images of the ricci flow solution
%
% Returns
% -------
%
% Saves to disk
% -------------
% Ricci flow solution
%   ricciFn = sprintf(QS.fullFileBase.rawRicciSolution, maxIter, tp) 
% rawRicciMesh has fields annulus and rectangle
%   rawRicciMeshFn = sprintf(QS.fullFileBase.rawRicciMesh, maxIter, tp) ;
% Beltrami coeff mu for this #iterations
%   mufn = sprintf(QS.fullFileBase.rawRicciMu, maxIter, tp) ;
%
% Images saved
% ------------
% Plot Ricci result in annulus
%   fullfile(imDir, sprintf('%06d_RicciSolution.png', tp)) ;
% plot beltrami for Ricci flow
%   fullfile(imDir, sprintf('%06d_RicciFlowBeltrami.png', tp)) ;
% Plot just the bare triangulation
%   fullfile(imDir, sprintf('%06d_RicciFlowSolution.png', tp)) ;
% Histogram |mu| for each case
%   fullfile(imDir, sprintf('%06d_BeltramiCoefficients.png', tp)) ;
% Image of corrected vertices on inner and outer annulus
%   fullfile(imDir, sprintf('%06d_ricci_InnerCorrection.png', tp)) ;
%   fullfile(imDir, sprintf('%06d_ricci_OuterCorrection.png', tp)) ;
% Orientation and branch cut
%   fullfile(imDir, sprintf('%06d_DrhoDphi_flipped.png', tp)) ;
%   fullfile(imDir, sprintf('%06d_phiOrderInitial.png', tp)) ;
%   fullfile(imDir, sprintf('%06d_phiOrderFinal.png', tp)) ;
%
%
% NPMitchell 2020

%% Default parameters
overwrite = false ;
radiusTolerance = 0.01 ;
maxIter = 200 ;
save_ims = true ;
if isempty(QS.t0)
    t0 = QS.xp.fileMeta.timePoints(1) ;
else
    t0 = QS.t0 ;
end
%% Unpack parameters
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'maxIter')
    maxIter = options.maxIter ;
end
if isfield(options, 't0')
    t0 = options.t0 ;
end
if isfield(options, 'radiusTolerance')
    radiusTolerance = options.radiusTolerance ;
end
if isfield(options, 'save_ims')
    save_ims = options.save_ims  ;
elseif isfield(options, 'saveIms')
    save_ims = options.saveIms ;
end


%% Load raw mesh that we want to flow by Ricci flow
if isfield(options, 'mesh')
    mesh = options.mesh ;
    % check that supplied mesh is topological cylinder
    try
        assert(eulerCharacteristic(mesh) == 0)
    catch
        disp('input mesh is not topological cylinder. Attempting to glue seam...')
        mesh = glueCylinderCutMeshSeam(mesh) ;
    end
else
    disp('Loading 3d vertices in embedding space (original data frame)')
    mesh = read_ply_mod(sprintf(QS.fullFileBase.cylinderMeshClean, tp)) ;
    tidx = QS.xp.tIdx(tp) ;
end
imDir = fullfile(QS.dir.rawRicci.mesh, 'images', sprintf('%04diter', maxIter)) ;
if ~exist(imDir, 'dir') && save_ims
    mkdir(imDir)
end

%% Generate conformal parameterization in the unit disk

rawRicciMeshFn = sprintf(QS.fullFileBase.rawRicciMesh, maxIter, tp) ;

if ~exist(rawRicciMeshFn, 'file') || overwrite
    disp(['rawRicciMesh not on disk, computing: ' rawRicciMeshFn])
    
    % Load or compute ricci flow solution (initial mapping)
    
    ricciFn = sprintf(QS.fullFileBase.rawRicciSolution, maxIter, tp) ;
    try
        load(ricciFn, 'U')
    catch
        disp('Ricci solution not on disk, computing...')
        [~, U, ~] = DiscreteRicciFlow.EuclideanRicciFlow(mesh.f, mesh.v, ...
            'BoundaryType', 'Fixed', 'BoundaryShape', 'Circles', ...
            'MaxCircIter', maxIter);
        % [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(cutMesh) ;
        % Let maxIter >= 50
        F = mesh.f ;
        save(ricciFn, 'U', 'F') ;
    end

    % Plot Ricci result in annulus
    outfn = fullfile(imDir, sprintf('%06d_RicciSolution.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        triplot(triangulation(mesh.f, U), 'color', 'k')
        axis equal
        axis tight
        title('Ricci flow solution', 'interpreter', 'latex')
        xlabel('$\tilde{x}$', 'interpreter', 'latex')
        ylabel('$\tilde{y}$', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    %% Beltrami
    mu_annulus0 = bc_metric(mesh.f, U, mesh.v, 3) ;
    
    %% Aligned cylinder clean mesh for plotting etc
    alignCylMesh = mesh ;
    alignCylMesh.v = QS.xyz2APDV(mesh.v) ;
    

    %% plot beltrami for Ricci flow
    mesh2d = triangulation(mesh.f, [U(:, 1), U(:, 2)]) ; %  ; , 0*U(:, 1)]) ;
    options.view = [0, 90] ;
    options.axisOff = true ;
    options.visible = true ;
    options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
        '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;
    outfn = fullfile(imDir, sprintf('%06d_RicciFlowBeltrami.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        nFieldsOnSurface({mesh2d, mesh2d, mesh2d, ...
            alignCylMesh, alignCylMesh, alignCylMesh}, ...
            {real(mu_annulus0), imag(mu_annulus0), abs(mu_annulus0), ...
            real(mu_annulus0), imag(mu_annulus0), abs(mu_annulus0)}, ...
            options)
        sgtitle('Beltrami coefficients from Ricci flow', ...
            'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    % %% Compare to beltrami from UVprime cutmesh
    % uvMesh = load('./uvpcutMesh_000160.mat') ;
    % uvMesh = uvMesh.uvpcutMesh.resampled ;
    % m2duv = struct() ;
    % m2duv.f = uvMesh.f ;
    % m2duv.v = uvMesh.u ;
    % mu_uv = bc_metric(uvMesh.f, uvMesh.u, uvMesh.v, 3) ;
    % 
    % % plot beltrami
    % options.view = [0, 90] ;
    % options.axisOff = true ;
    % options.visible = true ;
    % options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
    %     '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;
    % outfn = fullfile(imDir, sprintf('%06d_DirichletBeltrami.png', tp)) ;
    % if ~exist(outfn, 'file') && save_ims
    %     clf
    %     nFieldsOnSurface({m2duv, m2duv, m2duv, ...
    %         mesh, mesh, mesh}, ...
    %         {real(mu_uv), imag(mu_uv), abs(mu_uv), ...
    %         real(mu_uv), imag(mu_uv), abs(mu_uv)}, ...
    %         options)
    %     sgtitle('Beltrami coefficients from Dirichlet energy minimization', ...
    %         'interpreter', 'latex')
    %     saveas(gcf, outfn)
    % end

    %% Plot just the bare triangulation
    outfn = fullfile(imDir, sprintf('%06d_RicciFlowSolution.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        triplot(triangulation(mesh.f, U), 'color', 'k')
        axis equal; axis tight
        hold on
        scatter(0,0,50, 'filled', 'b')
        title('off-center inner boundary', 'interpreter', 'latex')
        xlabel('$\tilde{x}$', 'interpreter', 'latex')
        ylabel('$\tilde{y}$', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end


    %% Push barycenter of inner hole to origin via Mobius transformation
    boundaries = DiscreteRicciFlow.compute_boundaries(mesh.f) ;
    % which boundary is the inner one? Find which one has shorter circumference
    LL = [0, 0] ;
    for qq = 1:length(boundaries)
        % make periodic 
        LL(qq) = sum(sqrt(sum((U(boundaries{qq},:)- ...
            circshift(U(boundaries{qq},:), 1, 1)).^2, 2))) ;
    end
    % which is smaller?
    [~, innerID] = min(LL) ;
    inner = boundaries{innerID} ;
    outerID = setdiff([1,2], innerID) ;
    outer = boundaries{outerID} ;

    % barycenter of inner -- THIS IS BIASED
    % baryc = mean(U(inner, :), 1) ;
    % baryc = complex(baryc(1), baryc(2)) ;

    % Take center (not barycenter) of circle
    % note: https://www.mathworks.com/matlabcentral/fileexchange/5557-circle-fit
    [xc,yc, innerR] = circfit(U(inner, 1), U(inner, 2)); 
    baryc = xc + 1j * yc ;

    % Covert U to complex
    zz = complex(U(:, 1), U(:, 2)) ;

    % Mobius transform to place center of annulus at origin
    zz = (zz - baryc) ./ (1 - conj(baryc) .* zz) ;
    UU = [real(zz), imag(zz) ] ;

    % inspect centered mesh
    clf
    triplot(triangulation(mesh.f, UU))
    hold on;
    plot(UU(boundaries{1}, 1), UU(boundaries{1}, 2), 'co-')
    plot(UU(boundaries{2}, 1), UU(boundaries{2}, 2), 'co-')
    scatter(0,0, 'r', 'filled')
    axis equal
    xlim([-2*innerR, 2*innerR])

    %% Enforce circularity in inner Boundary
    % push inner boundary onto circle with exact radius
    phaseInner = atan2(UU(inner, 2), UU(inner, 1)) ; 
    % check that this is a minor correction
    radii = vecnorm(UU(inner, :), 2, 2) ;
    disp(['Correcting radial coordinate by a maximum of ' ...
        num2str(max(abs(radii - innerR) / innerR)*100) '%'])
    try
        assert(max(abs(radii - innerR) / innerR) < radiusTolerance)
    catch
        msg = 'ERROR: adjustment change is larger than tolerance! ' ;
        msg = [msg '\n change=' num2str(max(abs(radii - innerR) / innerR))] ;
        msg = [msg '\n tol=' num2str(radiusTolerance) ] ;
        msg = [msg '\n Continue with large change? [N/y]:'] ;
        cont = input(msg, 's') ;
        if ~contains(lower(cont), 'y')
            error('ERROR: adjustment change is larger than tolerance!')
        end
    end
    UU(inner, 1) = innerR * cos(phaseInner) ;
    UU(inner, 2) = innerR * sin(phaseInner) ;

    % Enforce circularity in outer boundary ==> radius=1
    phaseOuter = atan2(UU(outer, 2), UU(outer, 1)) ; 
    % check that this is a minor correction
    radii = vecnorm(UU(outer, :), 2, 2) ;
    disp(['Correcting radial coordinate by a maximum of ' ...
        num2str(max(abs(radii - 1))*100) '%'])
    assert(max(abs(radii - 1)) < radiusTolerance)
    UU(outer, 1) = cos(phaseOuter) ;
    UU(outer, 2) = sin(phaseOuter) ;

    % inspect the adjusted mesh
    outfn1 = fullfile(imDir, sprintf('%06d_ricci_InnerCorrection.png', tp)) ;
    outfn2 = fullfile(imDir, sprintf('%06d_ricci_OuterCorrection.png', tp)) ;
    if ~exist(outfn1, 'file')
        triplot(triangulation(mesh.f, UU), 'color', 'k')
        hold on;
        plot(UU(boundaries{1}, 1), UU(boundaries{1}, 2), 'b.-')
        plot(UU(boundaries{2}, 1), UU(boundaries{2}, 2), 'b.-')
        scatter(0,0, 'r', 'filled')
        axis equal;
        xlim([-2*innerR, 2*innerR])
        xlabel('$\tilde{x}$', 'interpreter', 'latex')
        ylabel('$\tilde{y}$', 'interpreter', 'latex')
        sgtitle('Ricci flow result with circularity correction', ...
            'interpreter', 'latex')
        saveas(gcf, outfn1)
        ylim([-1,1])
        axis equal; axis tight
        saveas(gcf, outfn2)
    end
    clf

    %% Take log
    zz = complex(UU(:, 1), UU(:, 2)) ;
    rho = real(log(zz)) ;
    phi = imag(log(zz)) ;
    triplot(triangulation(mesh.f, [rho, phi])) ;
    axis equal

    %% Attempt to rectify orientation according to the following rules:
    % rho should increase with x of alignedMesh (ap axis), on average
    % phi should increase with angle of anterior boundary about its mean
    
    
    % Check direction
    if corr(alignCylMesh.v(:, 1), rho(:)) < 0
        rho = - rho ;
    end
    [sepB, hairpins] = separateFreeBoundaryCurves(freeBoundary(triangulation(mesh.f, mesh.v))) ;
    
    try
        assert(~any(isempty(hairpins)))
    catch
        error('Cleaned cylinder meshes ought not to have hairpin loops on boundaries')
    end
    
    % which endcap is which -- for now ask for x value in APDV coords
    x1 = mean(alignCylMesh.v(sepB{1}(:, 1), 1)) ;
    x2 = mean(alignCylMesh.v(sepB{2}(:, 1), 1)) ;
    if x2 < x1 
        sepBnew = {sepB{2}, sepB{1}} ;
        sepB = sepBnew ;
    end
    
    % Take phase around tangent plane of each boundary
    flipPhi = [false, false] ;
    for endcapId = 1:2
        endCap = alignCylMesh.v(sepB{endcapId}(:, 1), :)...
            - mean(alignCylMesh.v(sepB{endcapId}(:, 1), :)) ;

        dx = endCap(:, 1) - endCap(:, 1)' ;
        dy = endCap(:, 2) - endCap(:, 2)' ;
        dz = endCap(:, 3) - endCap(:, 3)' ;
        maxDist = max(max(sqrt(dx.^2+dy.^2+dz.^2))) ;
        ptCloud1 = pointCloud(endCap) ;
        Pmodel = pcfitplane(ptCloud1, maxDist) ;
        
        % Rotate so normal of plane of best fit is along x
        % rotM = rotate3dToAlignAxis(Pmodel.Normal) ;
        % rotEC = (rotM * endCap')' ;
        % if endcapId == 2
        %     % rotate by pi around z axis
        %     rotEC = ([-1 0 0; 0 -1 0; 0 0 1] * rotEC')' ;
        % end
        
        % Check 
        % plot3(rotEC(:, 1), rotEC(:, 2), rotEC(:, 3), '.')
        % axis equal; xlabel('x'); ylabel('y'); zlabel('z')
        
        % Phi estimate is the arctan in plane YZ.
        guessPhi = atan2(endCap(:, 3), endCap(:, 2)) ;
        
        % Check
        subplot(2, 1, 1)
        scatter3(endCap(:, 1), endCap(:, 2), endCap(:, 3), 5, guessPhi, 'filled')
        axis equal; xlabel('x'); ylabel('y'); zlabel('z')
        hold on;
        subplot(2, 1, 2)
        scatter3(endCap(:, 1), endCap(:, 2), endCap(:, 3), 5, phi(sepB{endcapId}(:, 1)), 'filled')
        axis equal; xlabel('x'); ylabel('y'); zlabel('z')
        hold on;
        % clf;
        % plot(mod(diff(phi(sepB{1}(:, 1))), 2*pi), mod(diff(guessPhi), 2*pi), '.')
        
        sign0 = mean(mod(diff(phi(sepB{endcapId}(:, 1))), 2*pi)) < pi ;
        signGuess = mean(mod(diff(guessPhi), 2*pi)) < pi ;
        
        if signGuess == sign0
            flipPhi(endcapId) = QS.flipy ;
        else
            flipPhi(endcapId) = ~QS.flipy ;
        end
        
    end
    
    % For now use only first boundary (anterior)
    if flipPhi(1)
        phi = -phi + pi ;
    % elseif any(flipPhi)
    %     error('Conflicting results for flipping phi from either endcap')
    end
    
    phi0 = phi ;
    
    minrho = min(rho) ;
    if abs(minrho) > 0 
        disp('Translating rho to origin')
        rho = rho - minrho ;
    end

    outfn = fullfile(imDir, sprintf('%06d_rhophi_rectified.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        subplot(2, 1, 1)
        trisurf(triangulation(alignCylMesh.f, ...
            alignCylMesh.v), rho, 'edgecolor', 'none')
        axis equal ; axis tight 
        colormap viridis
        cb = colorbar ;
        ylabel(cb, '$\rho$', 'interpreter', 'latex')
        subplot(2, 1, 2)
        trisurf(triangulation(alignCylMesh.f, ...
            alignCylMesh.v), phi, 'edgecolor', 'none')
        axis equal ; axis tight 
        colormap viridis
        cb = colorbar ;
        ylabel(cb, '$\phi$', 'interpreter', 'latex')
        sgtitle('reorienting pullback')
        saveas(gcf, outfn)
    end

    %% plot initial cutpath
    outfn = fullfile(imDir, sprintf('%06d_phiOrderInitial.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        trisurf(triangulation(mesh.f, alignCylMesh.v), phi, ...
            'edgecolor', 'none') ;
        axis equal; xlabel('x'); ylabel('y'); zlabel('z')
        title('initial branch cut', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    phi_recut = phi ;

    %% Minimize offset to all phi values based on previous mesh vertices 
    % in pullback space ---> ignore this for example script
    
    % Load apIDx for phi offset
    adIDx = h5read(QS.fileName.aBoundaryDorsalPtsClean,...
        ['/' sprintf('%06d', tp)]) ;
    phi_recut = mod(phi_recut - phi_recut(adIDx), 2*pi) ;


    %% plot final cutpath
    outfn = fullfile(imDir, sprintf('%06d_phiOrderFinal.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        trisurf(triangulation(alignCylMesh.f, alignCylMesh.v), phi_recut, ...
            'edgecolor', 'none') ;
        axis equal; xlabel('x'); ylabel('y'); zlabel('z')
        title('final branch cut', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    %% Unwrap from annulus into rectangle and save in struct
    rawRicciMesh = struct() ;
    rawRicciMesh.annulus = struct('f', mesh.f, 'u', UU, ...
        'v', mesh.v, 'vrs', alignCylMesh.v) ;
    riccicutMesh = rawRicciMesh.annulus ;
    riccicutMesh.u = [rho, phi_recut] ;
    rawRicciMesh.rectangle = struct('f', riccicutMesh.f, 'u', ...
        riccicutMesh.u, 'v', mesh.v, 'vrs', alignCylMesh.v) ;

    % Save rawRicciMesh
    % note: rawRicciMesh has fields annulus and rectangle
    disp(['Saving rawRicciMesh to ' rawRicciMeshFn])
    save(rawRicciMeshFn, 'rawRicciMesh')
    
    computed_mesh = true ;
        

    %% Preview final ordering
    outfn = fullfile(imDir, sprintf('%06d_ricciFinal.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        facecolors = (1:length(rawRicciMesh.rectangle.f)) ;
        trisurf2d(rawRicciMesh.rectangle.f, rawRicciMesh.rectangle.u,...
            'FaceVertexCData', facecolors(:), 'Facealpha',0.5)

        saveas(gcf, outfn)
    end
else
    disp('rawRicciMesh already on disk, loading...')
    load(rawRicciMeshFn, 'rawRicciMesh') 
    computed_mesh = false ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recompute mu in rectangular space 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mufn = sprintf(QS.fullFileBase.rawRicciMu, maxIter, tp) ;    
if ~exist(mufn, 'file') || computed_mesh || overwrite
    
    % Compute annular mu after all the alterations
    mu_annulus = bc_metric(mesh.f, rawRicciMesh.annulus.u, mesh.v, 3) ;
    
    % It should be conformal since this is a conformal map of a conformal map
    % into the disk
    rMesh = rawRicciMesh.rectangle ;
    mu_rectangle = bc_metric(rMesh.f, rMesh.u, rMesh.v, 3) ;

    %% plot beltrami for Ricci flow
    mesh2d = triangulation(rMesh.f, rMesh.u) ;
    options.view = [0, 90] ;
    options.axisOff = true ;
    options.visible = true ;
    options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
        '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;

    outfn = fullfile(imDir, sprintf('%06d_RicciFlowBeltrami_rectangle.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        nFieldsOnSurface({mesh2d, mesh2d, mesh2d, ...
            alignCylMesh, alignCylMesh, alignCylMesh}, ...
            {real(mu_rectangle), imag(mu_rectangle), abs(mu_rectangle), ...
            real(mu_rectangle), imag(mu_rectangle), abs(mu_rectangle)}, ...
            options)
        sgtitle('Beltrami coefficients from $\log$ of Ricci flow', ...
            'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    %% Histogram |mu| for each case
    outfn = fullfile(imDir, sprintf('%06d_BeltramiCoefficients.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        maxx = max([max(abs(mu_annulus)), max(abs(mu_rectangle))]) ;
        % subplot(3, 1, 1)
        % histogram(abs(mu_uv))
        % xlim([0, maxx])
        % xlabel('$|\mu|$ for Dirichlet minimization', 'interpreter', 'latex') 
        % ylabel('counts', 'interpreter', 'latex')
        subplot(2, 1, 1)
        histogram(abs(mu_annulus))
        xlim([0, maxx])
        xlabel('$|\mu|$ for Ricci flow to annulus', 'interpreter', 'latex')
        ylabel('counts', 'interpreter', 'latex')
        subplot(2, 1, 2)
        histogram(abs(mu_rectangle))
        xlim([0, maxx])
        xlabel('$|\mu|$ for rectilinear domain from Ricci flow', 'interpreter', 'latex')
        ylabel('counts', 'interpreter', 'latex')
        sgtitle('Conformality test', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end


    %% Save mu for this #iterations
    disp('saving mu for this computation')
    save(mufn, 'mu_annulus', 'mu_rectangle')
else
    load(mufn, 'mu_annulus', 'mu_rectangle')
end

% Output beltramis if requested
if nargout > 1
    ricciMu = struct();
    ricciMu.mu_annulus = mu_annulus ;
    ricciMu.mu_rectangle = mu_rectangle ;
end