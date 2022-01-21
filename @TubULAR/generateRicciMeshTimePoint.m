function [ricciMesh, ricciMu] = generateRicciMeshTimePoint(QS, tp, options)
% generateRicciMeshTimePoint(QS, tp, options)
% 
% Compute solution to Ricci flow for each timepoint. If there is only one
% timepoint, we by default flow the spcutMesh. Otherwise we by default flow
% the spcutMeshSm.
% 
% See also
% --------
%  generateRawRicciMeshTimePoint(QS, tp, options) --> for raw Ricci mesh
%  before re-parameterization of the mesh into spcutMesh or uv coords etc
% 
% Parameters
% ----------
% QS : quapSlap class instance
% tp : int (default=QS.t0)
%   timestamp of timepoint for which to generate Ricci flow mesh
% options : optional struct with optional fields
%   resample : bool (default = true)
%       perform isotropic resampling of the mesh before performing Ricci
%       flow to ensure vertex quality during operation
%   maxIter : int (default=100)
%   radiusTolerance : float (default=0.01)
%       maximum allowed fractional deviation of Ricci flow solution inner
%       and outer radius from a true circle with fixed radius (variable
%       inner radius, fixed outer radius of 1)
%   save_ims : bool (default = true)
%       save images of the ricci flow solution
%   pathline_computation : bool
%       compute the ricci flow for the pathline mesh of this timepoint
%   t0Pathlines
%       only used if pathline_computation, which pathline to look onto for
%       this timpoint (for example load the mesh for t=50 advected along
%       pathlines from t0=1).
%
% Returns
% -------
%
% Saves to disk
% -------------
% Ricci flow solution
%   ricciFn = sprintf(QS.fullFileBase.ricciSolution, maxIter, tp) 
% ricciMesh has fields annulus and rectangle
%   ricciMeshFn = sprintf(QS.fullFileBase.ricciMesh, maxIter, tp) ;
% Beltrami coeff mu for this #iterations
%   mufn = sprintf(QS.fullFileBase.ricciMu, maxIter, tp) ;
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
t0 = QS.t0set() ;
radiusTolerance = 0.01 ;
maxIter = 200 ;
save_ims = true ;
pathline_computation = false ;  % compute the ricci flow for the pathline mesh of this timepoint
t0Pathlines = t0 ;              % only used if pathline_computation 
if numel(QS.xp.fileMeta.timePoints) > 1
    coordSys = 'spsm' ;
else
    coordSys = 'sp' ;
end
resample = true ;

%% Unpack parameters
if nargin < 2
    tp = t0 ;
elseif isempty(tp)
    tp = t0 ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'maxIter')
    maxIter = options.maxIter ;
end
if isfield(options, 'radiusTolerance')
    radiusTolerance = options.radiusTolerance ;
end
if isfield(options, 'save_ims')
    save_ims = options.save_ims  ;
elseif isfield(options, 'saveIms')
    save_ims = options.saveIms ;
end
if isfield(options, 'resample')
    resample = options.resample ;
end

% In case we are computing for pathline, unpack relevant options
if isfield(options, 'pathline_computation')
    pathline_computation = options.pathline_computation ;
end
% Demand t0Pathlines is explicitly passed if we do pathlines so as to not
% confusedly write to a different set of pathlines
if pathline_computation
    if isfield(options, 't0Pathlines')
        t0Pathlines = options.t0Pathlines ;
    else
        error('Must pass t0Pathlines if pathline_computation')
    end
    
    % Ensure directories exist as well
    if resample
        dirs2make = {...
            sprintf(QS.dir.pathlines.ricci.dataWithResampling, t0Pathlines), ...
            sprintf(QS.dir.pathlines.ricci.meshWithResampling, t0Pathlines), ...
            sprintf(QS.dir.pathlines.ricci.quasiconformalWithResampling, t0Pathlines), ...
            sprintf(QS.dir.pathlines.ricci.solutionWithResampling, t0Pathlines), ...
            sprintf(QS.dir.pathlines.ricci.muWithResampling, t0Pathlines)} ;
    else
        dirs2make = {...
            sprintf(QS.dir.pathlines.ricci.data, t0Pathlines), ...
            sprintf(QS.dir.pathlines.ricci.mesh, t0Pathlines), ...
            sprintf(QS.dir.pathlines.ricci.quasiconformal, t0Pathlines), ...
            sprintf(QS.dir.pathlines.ricci.solution, t0Pathlines), ...
            sprintf(QS.dir.pathlines.ricci.mu, t0Pathlines)} ;
   
    end
    for qq = 1:length(dirs2make)
        if ~exist(dirs2make{qq}, 'dir')
            mkdir(dirs2make{qq})
        end
    end
end    

%% Load cutmesh to Ricci flow
if pathline_computation
    if isfield(options, 'cutMesh')
        cutMesh = options.cutMesh ;
        glueMesh = glueCylinderCutMeshSeam(cutMesh) ;
        try
            assert(isfield(cutMesh, 'pathPairs'))
        catch 
            refMesh = load(sprintf(QS.fileName.pathlines.refMesh. t0Pathlines)) ;
            cutMesh.pathPairs = refMesh.pathPairs ;
        end
    else
        disp('Loading 3d vertices in embedding RS space (rotated & scaled)')
        cutMesh = load(sprintf(QS.fileName.pathlines.refMesh. t0Pathlines)) ;
        tidx = QS.xp.tIdx(tp) ;
        cutMesh.v = zeros(size(vP3d.vX, 2) * size(vP3d.vX, 3), 3) ;
        cutMesh.v(:, 1) = reshape(squeeze(vP3d.vXrs(tidx, :, :)), [], 1) ;
        cutMesh.v(:, 2) = reshape(squeeze(vP3d.vYrs(tidx, :, :)), [], 1) ;
        cutMesh.v(:, 3) = reshape(squeeze(vP3d.vZrs(tidx, :, :)), [], 1) ;
    end
    nU = cutMesh.nU ;
    nV = cutMesh.nV ;
    if resample
        pathlineRicciDir = sprintf(QS.dir.pathlines.ricci.dataWithResampling, t0Pathlines) ;
    else
        pathlineRicciDir = sprintf(QS.dir.pathlines.ricci.data, t0Pathlines) ;
    end
    imDir = fullfile(pathlineRicciDir, 'images', sprintf('%04diter', maxIter)) ;
else
    if strcmpi(coordSys, strrep('spsm', '_', ''))
        cutMesh = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
        cutMesh = cutMesh.spcutMeshSmRS ;
    elseif strcmpi(coordSys, strrep('sp', '_', ''))
        cutMesh = load(sprintf(QS.fullFileBase.spcutMesh, tp)) ;
        tmp = cutMesh.spcutMesh ;
        cutMesh = struct() ;
        cutMesh.nU = tmp.nU ;
        cutMesh.nV = tmp.nV ;
        cutMesh.u = tmp.sphi ;
        cutMesh.pathPairs = tmp.pathPairs ;
        cutMesh.f = tmp.f ;
        cutMesh.v = QS.xyz2APDV( tmp.v ) ;
    else
        error('did not recognize coordSys upon which we perform Ricci flow')
    end
    glueMesh = glueCylinderCutMeshSeam(cutMesh) ;
    nU = cutMesh.nU ;
    nV = cutMesh.nV ;
    if resample
        imDir = fullfile(QS.dir.ricci.meshWithResampling, ...
            'images', sprintf('%04diter', maxIter)) ;
    else
        imDir = fullfile(QS.dir.ricci.mesh, 'images', sprintf('%04diter', maxIter)) ;
    end
end
if ~exist(imDir, 'dir') && save_ims
    mkdir(imDir)
end

%% Generate conformal parameterization in the unit disk
if pathline_computation
    if resample
        ricciMeshFn = sprintf(QS.fullFileBase.pathlines.ricciMeshWithResampling, t0Pathlines, maxIter, tp) ;
    else
        ricciMeshFn = sprintf(QS.fullFileBase.pathlines.ricciMesh, t0Pathlines, maxIter, tp) ;
    end
else
    if resample 
        ricciMeshFn = sprintf(QS.fullFileBase.ricciMeshWithResampling, maxIter, tp) ;
    else
        ricciMeshFn = sprintf(QS.fullFileBase.ricciMesh, maxIter, tp) ;
    end
end
if ~exist(ricciMeshFn, 'file') || overwrite
    disp(['ricciMesh not on disk, computing: ' ricciMeshFn])
    
    % Load or compute ricci flow solution (initial mapping)
    if pathline_computation
        if resample
            ricciFn = sprintf(...
                QS.fullFileBase.pathlines.ricciSolutionWithResampling, ...
                t0Pathlines, maxIter, tp) ;
        else
            ricciFn = sprintf(QS.fullFileBase.pathlines.ricciSolution, ...
                t0Pathlines, maxIter, tp) ;
        end
    else
        if resample
            ricciFn = sprintf(QS.fullFileBase.ricciSolutionWithResampling, maxIter, tp) ;
        else
            ricciFn = sprintf(QS.fullFileBase.ricciSolution, maxIter, tp) ;
        end
    end
    try
        load(ricciFn, 'U')
        
        % Check that the result obeys topological constraint
        if resample
            load(ricciFn, 'gMeshResample')
            load(ricciFn, 'glueMesh')
            % Note: eulerCharacteristic of an annulus is 0.
            % Check topological indices of resampled mesh
            assert(eulerCharacteristic(gMeshResample) == ...
                eulerCharacteristic(glueMesh))
        end
    catch
        disp('Ricci solution not on disk, computing...')
        
        if resample
            tarEdgePctile =  5 ;
            success_conv = false ;
            while ~success_conv
                try
                    numIterations = 10 ;
                    % by default choose edgelength target to be smaller than mean
                    eL = edge_lengths(glueMesh.v, glueMesh.f) ;
                    targetEdgeLength = prctile(eL(:), tarEdgePctile) ;
                    [ gMeshResample, bc0, qF0] = ...
                        isotropicRemeshAnnuluarCutMesh(...
                        glueMesh, targetEdgeLength, numIterations) ;
                    
                    % Note: eulerCharacteristic of an annulus is 0.
                    % Check topological indices of resampled mesh
                    assert(eulerCharacteristic(gMeshResample) == ...
                        eulerCharacteristic(glueMesh))
                
                    [~, U0, ~] = DiscreteRicciFlow.EuclideanRicciFlow(gMeshResample.f, gMeshResample.v, ...
                        'BoundaryType', 'Fixed', 'BoundaryShape', 'Circles', ...
                        'MaxCircIter', maxIter);
                    success_conv = true ;
                catch
                    fprintf(...
                        ['Could not make mesh: decreasing', ...
                        ' tarEdgePctile = %0.f by 20 percent...\n'],...
                        tarEdgePctile)
                    tarEdgePctile = tarEdgePctile * 0.8 ;
                    numIterations = numIterations + 1 ;
                end
            end
            % [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(cutMesh) ;
            % Let maxIter >= 50
            
            % -------------------------------------------------------------
            % This DOESN'T work since periodic -- must query faces and
            % barycenters in 3D, not 2D 
            % pointLocation(triangulation(gMeshResample.f, ...
            %    gMeshResample.u, glueMesh.u) ;
            % -------------------------------------------------------------
            
            samplingOptions = struct('epsVtx', 1) ;
            [U, qF, bc] = barycentricMap3DtoND(...
                gMeshResample.f, gMeshResample.v,...
                U0, glueMesh.v, samplingOptions) ;

            readme = ['resampled with numIterations, so now ', ...
                'resampled coords: ', ...
                'U0 = bc0(:, 1) .* U(glueMesh.f(qF0, 1), :) + ' , ...
                'bc0(:, 2) .* U(glueMesh.f(qF0, 2), :) + ', ...
                'bc0(:, 3) .* U(glueMesh.f(qF0, 3), :) ;', ...
                ' and regular coords:' , ...
                ' U = bc(:, 1) .* U0(gMeshResample.f(qF, 1), :) + ' , ...
                'bc(:, 2) .* U0(gMeshResample.f(qF, 2), :) + ', ...
                'bc(:, 3) .* U0(gMeshResample.f(qF, 3), :) ;'] ;
            save(ricciFn, 'U', 'U0', 'bc0', ...
                'qF0', 'readme', ...
                'bc', 'qF', ...
                'numIterations', 'targetEdgeLength', ...
                'gMeshResample', 'glueMesh') ;
            
            % -------------------------------------------------------------
            % Check visually
            % trisurf(triangulation(gMeshResample.f, gMeshResample.v), ...
            %      'edgecolor', 'none') ; 
            %  hold on ; 
            %  scatter3(glueMesh.v(:, 1), glueMesh.v(:, 2), ...
            %      glueMesh.v(:, 3), 10, 'filled')
            close all
             triplot(triangulation(gMeshResample.f, U0), 'color',...
                 [0    0.4470    0.7410]) ; 
             hold on ; 
             triplot(triangulation(glueMesh.f, U),'color',...
                 [0.8500    0.3250    0.0980]) ; 
             axis equal
             xlabel('u'); ylabel('v')
             legend('resampled', 'original cutMesh')
             saveas(gcf, [ricciFn(1:end-3) 'png'])
            % -------------------------------------------------------------
            
            
        else
            [~, U, ~] = DiscreteRicciFlow.EuclideanRicciFlow(glueMesh.f, glueMesh.v, ...
                'BoundaryType', 'Fixed', 'BoundaryShape', 'Circles', ...
                'MaxCircIter', maxIter);
            % [labels, dbonds, topStructTools] = labelRectilinearMeshBonds(cutMesh) ;
            % Let maxIter >= 50
            save(ricciFn, 'U') ;
        end
    end

    % Plot Ricci result in annulus
    outfn = fullfile(imDir, sprintf('%06d_RicciSolution.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        triplot(triangulation(glueMesh.f, U), 'color', 'k')
        axis equal
        axis tight
        title('Ricci flow solution', 'interpreter', 'latex')
        xlabel('$\tilde{x}$', 'interpreter', 'latex')
        ylabel('$\tilde{y}$', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    %% Beltrami
    mu_annulus0 = bc_metric(glueMesh.f, U, glueMesh.v, 3) ;

    %% plot beltrami for Ricci flow
    mesh2d = triangulation(glueMesh.f, [U(:, 1), U(:, 2)]) ; %  ; , 0*U(:, 1)]) ;
    options.view = [0, 90] ;
    options.axisOff = true ;
    options.visible = true ;
    options.labels = {'$\Re\mu$', '$\Im\mu$', '$|\mu|$', ...
        '$\Re\mu$', '$\Im\mu$', '$|\mu|$'} ;
    outfn = fullfile(imDir, sprintf('%06d_RicciFlowBeltrami.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        nFieldsOnSurface({mesh2d, mesh2d, mesh2d, ...
            glueMesh, glueMesh, glueMesh}, ...
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
    %         glueMesh, glueMesh, glueMesh}, ...
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
        triplot(triangulation(glueMesh.f, U), 'color', 'k')
        axis equal; axis tight
        hold on
        scatter(0,0,50, 'filled', 'b')
        title('off-center inner boundary', 'interpreter', 'latex')
        xlabel('$\tilde{x}$', 'interpreter', 'latex')
        ylabel('$\tilde{y}$', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end


    %% Push barycenter of inner hole to origin via Mobius transformation
    boundaries = DiscreteRicciFlow.compute_boundaries(glueMesh.f) ;
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
    triplot(triangulation(glueMesh.f, UU))
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
    try
        assert(max(abs(radii - 1)) < radiusTolerance)
    catch
        error('WARNING: RADIUS IS NOT SUFFICIENTLY CIRCULAR')
    end
    UU(outer, 1) = cos(phaseOuter) ;
    UU(outer, 2) = sin(phaseOuter) ;

    % inspect the adjusted mesh
    outfn1 = fullfile(imDir, sprintf('%06d_ricci_InnerCorrection.png', tp)) ;
    outfn2 = fullfile(imDir, sprintf('%06d_ricci_OuterCorrection.png', tp)) ;
    if ~exist(outfn1, 'file')
        triplot(triangulation(glueMesh.f, UU), 'color', 'k')
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
    cutMesh.pathPairs(:, 1)
    zz = complex(UU(:, 1), UU(:, 2)) ;
    rho = real(log(zz)) ;
    phi = imag(log(zz)) ;
    triplot(triangulation(glueMesh.f, [rho, phi])) ;
    axis equal
    rhoM = reshape(rho, [nU, nV-1]) ;
    phiM = reshape(phi, [nU, nV-1]) ;

    % Check direction
    grX = diff(rhoM, 1, 1) ;
    grY = diff(phiM, 1, 2) ;
    outfn = fullfile(imDir, sprintf('%06d_DrhoDphi_rectified.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        subplot(2, 1, 1)
        scatter(reshape(rhoM(1:nU-1, :), [], 1), ...
            reshape(phiM(1:nU-1,:), [], 1), 10, grX(:))
        title('$\Delta \tilde{u}$', 'interpreter', 'latex')
        xlabel('$\tilde{u}$', 'interpreter', 'latex')
        xlabel('$\tilde{v}$', 'interpreter', 'latex')
        caxis([-0.05, 0.05])
        axis equal ; axis tight 
        colormap blueblackred
        cb = colorbar ;
        ylabel(cb, '$\Delta \tilde{u}$', 'interpreter', 'latex')
        subplot(2, 1, 2)
        scatter(reshape(rhoM(:, 1:end-1), [], 1), ...
            reshape(phiM(:,1:end-1), [], 1), 10, grY(:))
        title('$\Delta \tilde{v}$', 'interpreter', 'latex')
        xlabel('$\tilde{u}$', 'interpreter', 'latex')
        xlabel('$\tilde{v}$', 'interpreter', 'latex')
        caxis([-0.05, 0.05])
        axis equal ; axis tight 
        colormap blueblackred
        cb = colorbar() ;
        sgtitle('initial pullback orientation')
        ylabel(cb, '$\Delta \tilde{v}$', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    phi0 = phi ;
    % Push y range to approx (0, 2pi), fliping if needed
    if all(mean(sign(grY), 2) < 0)
        disp('Flipping angles phi')
        % Push phi to approx (0, 2pi)
        phi = pi - phi ;
    elseif all(mean(sign(grY), 2) > 0)
        disp('phi angles properly ordered')
        % Push phi to approx (0, 2pi)
        phi = phi + pi ;
    else
        error('Could not determine direction of phi in meshgrid')
    end
    % Determine whether to flip in x direction
    if all(all(sign(grX) < 0))
        rho = -rho ;
    end
    minrho = min(rho) ;
    if abs(minrho) > 0 
        disp('Translating rho to origin')
        rho = rho - minrho ;
    end

    % Remake grids of rhoM and phiM
    rhoM = reshape(rho, [nU, nV-1]) ;
    phiM = reshape(phi, [nU, nV-1]) ;

    % Check direction
    grX = diff(rhoM, 1, 1) ;
    grY = diff(phiM, 1, 2) ;

    outfn = fullfile(imDir, sprintf('%06d_DrhoDphi_flipped.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        subplot(2, 1, 1)
        scatter(reshape(rhoM(1:nU-1, :), [], 1), ...
            reshape(phiM(1:nU-1,:), [], 1), 10, grX(:))
        title('$\Delta \tilde{u}$', 'interpreter', 'latex')
        xlabel('$\tilde{u}$', 'interpreter', 'latex')
        ylabel('$\tilde{v}$', 'interpreter', 'latex')
        caxis([-0.05, 0.05])
        axis equal ; axis tight 
        colormap blueblackred
        cb = colorbar ;
        ylabel(cb, '$\Delta \tilde{u}$', 'interpreter', 'latex')
        subplot(2, 1, 2)
        scatter(reshape(rhoM(:, 1:end-1), [], 1), ...
            reshape(phiM(:,1:end-1), [], 1), 10, grY(:))
        title('$\Delta \tilde{v}$', 'interpreter', 'latex')
        xlabel('$\tilde{u}$', 'interpreter', 'latex')
        ylabel('$\tilde{v}$', 'interpreter', 'latex')
        caxis([-0.05, 0.05])
        axis equal ; axis tight 
        colormap blueblackred
        cb = colorbar() ;
        ylabel(cb, '$\Delta \tilde{v}$', 'interpreter', 'latex')
        sgtitle('reorienting pullback')
        saveas(gcf, outfn)
    end

    %% plot initial cutpath
    outfn = fullfile(imDir, sprintf('%06d_phiOrderInitial.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        triplot(triangulation(glueMesh.f, [rho, phi]), 'color', 'k') ;
        hold on;
        scatter(rho(1:nU), phi(1:nU), 'filled', 'c')
        axis equal; axis tight
        title('initial cut path', 'interpreter', 'latex')
        xlabel('$\tilde{u}$', 'interpreter', 'latex')
        ylabel('$\tilde{v}$', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    %% Push all 1:nU vertices to near phi = 0 (cutPath leveling in PB space)
    % Make first phi in each row the lowest phi value
    clf 
    phi_recut = phi ;
    for qq = 1:nU
        % Consider this column of the rectilinear pullback structure
        phis = phi(qq:nU:nU*(nV-1)) ;
        rhos = rho(qq:nU:nU*(nV-1)) ;
        phis(phis < phis(1)) = phis(phis < phis(1)) + 2*pi ;

        % Correct all first phis to be in the same register 
        % Note: this means they might not be on the same branch cut!
        if qq == 1
            % For first column, make phi(1) zero
            phi0 = phis(1) ;
            phis = phis - phi0 ;
        else
            % For subsequent columns, find the right branch cut that connects
            % most closely to the adjacent column without allowing
            % a flip in the normal of the mesh. Note that this assumes nothing
            % totally insane is happening to the geometry of the mesh across
            % mesh sampling distance.
            phis = phis - phi0 ;

            % I think it is impossible to be off by more than one branch cut --
            % ie by more than 2pi, so lets find the best match of those
            % possibilities in adjacent branch cuts
            possible = [-2*pi, 0, 2*pi] ;
            [min_dphi, ind_min] = min(abs(phis(1) + possible - prevphi1)) ;
            phis = phis + possible(ind_min) ;
            disp(['Translating phi values by '...
                num2str(round(possible(ind_min)/pi)) 'pi'])

            % todo: perform check on the sign of the cross product of a face
            % including this vertex to check for no mesh intersections
            plot(rhos, phis, '.'); 
            hold on ;
        end

        phi_recut(qq:nU:nU*(nV-1)) = phis ;
        prevphi1 = phis(1) ;
    end

    %% Minimize offset to all phi values based on previous mesh vertices 
    % in pullback space ---> ignore this for example script
    tidx = QS.xp.tIdx(tp) ;
    tidx0 = QS.xp.tIdx(QS.t0set()) ;
    if tidx ~= tidx0 
        % load next timepoint phi values
        prevTP = QS.xp.fileMeta.timePoints(tidx0) ;
        try
            prevMesh = load(sprintf(QS.fullFileBase.ricciMesh, maxIter, prevTP)) ;
        catch
            msg = 'Could not load tidx0 mesh -- run t0 timepoint first!' ;
            msg = [msg '--> ' sprintf(QS.fullFileBase.ricciMesh, maxIter, prevTP)] ;
            error(msg)
        end
        phi2match = prevMesh.ricciMesh.rectangle.u(1:nU*(nV-1), 2) ;
        overall_offset = mean(phi_recut(:) - phi2match(:)) ; 
        phi_recut = phi_recut - overall_offset ;
    end

    %% plot final cutpath
    outfn = fullfile(imDir, sprintf('%06d_phiOrderFinal.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        clf
        triplot(triangulation(glueMesh.f, [rho, phi_recut]), 'color', 'k') ;
        hold on;
        scatter(rho(1:nU), phi_recut(1:nU), 'filled','c')
        axis equal; axis tight
        title('final cut path', 'interpreter', 'latex')
        xlabel('$\tilde{u}$', 'interpreter', 'latex')
        ylabel('$\tilde{v}$', 'interpreter', 'latex')
        saveas(gcf, outfn)
    end

    %% Preview final ordering
    outfn = fullfile(imDir, sprintf('%06d_phiOrderFinal.png', tp)) ;
    if ~exist(outfn, 'file') && save_ims
        facecolors = (1:length(ricciMesh.rectangle.f)) ;
        facecolors = mod(facecolors, 2*nV-2) ;
        facecolors(facecolors==0) = 2*nV - 2 ;
        trisurf2d(ricciMesh.rectangle.f, ricciMesh.rectangle.u, 'FaceVertexCData', facecolors(:), 'Facealpha',0.5)

        saveas(gcf, outfn)
    end

    %% Unwrap from annulus into rectangle and save in struct
    ricciMesh = struct() ;
    ricciMesh.annulus = struct('f', glueMesh.f, 'u', UU, 'v', glueMesh.v, ...
        'nU', nU, 'nV', nV) ;
    riccicutMesh = ricciMesh.annulus ;
    riccicutMesh.u = [rho, phi_recut] ;
    opts = struct('vmax', 2 * pi, 'ignoreRectangularConstraint', true) ;
    riccicutMesh = cutRectilinearCylMesh(riccicutMesh, opts) ;
    ricciMesh.rectangle = struct('f', riccicutMesh.f, 'u', ...
        riccicutMesh.u, 'v', riccicutMesh.v, ...
        'nU', nU, 'nV', nV) ;

    % Save ricciMesh
    % note: ricciMesh has fields annulus and rectangle
    disp(['Saving ricciMesh to ' ricciMeshFn])
    save(ricciMeshFn, 'ricciMesh')
    
    computed_mesh = true ;
else
    disp('ricciMesh already on disk, loading...')
    load(ricciMeshFn, 'ricciMesh') 
    computed_mesh = false ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recompute mu in rectangular space 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pathline_computation
    mufn = sprintf(QS.fullFileBase.pathlines.ricciMu, t0Pathlines, maxIter, tp) ;
else
    mufn = sprintf(QS.fullFileBase.ricciMu, maxIter, tp) ;    
end
if ~exist(mufn, 'file') || computed_mesh || overwrite
    
    % Compute annular mu after all the alterations
    mu_annulus = bc_metric(glueMesh.f, ricciMesh.annulus.u, glueMesh.v, 3) ;
    
    % It should be conformal since this is a conformal map of a conformal map
    % into the disk
    rMesh = ricciMesh.rectangle ;
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
            glueMesh, glueMesh, glueMesh}, ...
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