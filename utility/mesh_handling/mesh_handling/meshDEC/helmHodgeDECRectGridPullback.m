function [divs, rots, harms, lapvs, glueMesh] = ...
    helmHodgeDECRectGridPullback(cutM, facevf, Options, varargin)
% helmHodgeDECRectGridPullback(cutM, facevf, varargin)
%   Perform Hodge decomposition using discrete exterior calculus on
%   vector field defined on faces. Clip, denoise and/or smooth the results 
%   based on varargin.
%
% Parameters
% ----------
% cutM : struct with fields v,f
% facevf : #faces x 3 float array
%   vector field defined on faces of cutM?
% Options : struct with fields
%   lambda : smoothing diffusion constant
% varargin : keyword arguments (optional)
%   nSpectralFilterModes : int (default = 0)
%       number of circumferential modes to keep in spectral filter
%       (low-pass filtering)
%   spectralFilterWidth : int (default = 0)
%       longitudinal half-width along which to average modes in spectral
%       filtering as tripulse filter.
%   niterSmoothing : int or three ints 
%       how many smoothing steps to perform on U fields (div,rot,harm
%       scalar fields)
%       If two values given, applies to div, rot separately
%       If method is 'both', applies to 
%       div-denoise, div-smooth, rot-denoise, rot-smooth
%   filterMethod : 'smooth' 'denoise' 'both'
%       What smoothing method to apply to U fields (div,rot,harm scalar
%       fields)
%           denoise : average local extrema with neighboring face values
%           smooth : average all values with neighboring face values
%           both : first denoise, then smooth
%   epsilon : float
%       small value setting threshold for local extremum for method == denoise
%   clipDiv : list of two floats
%       values to clip the rotation field
%   clipRot : list of two floats
%       values to clip the curl field
%   preview : bool 
%       view intermediate results
%   outdir : str
%       output directory for preview results and smoothing analysis
%
% Returns
% -------
% divs : struct
% rots : struct 
% harms : struct
% glueMesh : struct
%   glued rectilinear mesh
%
% NPMitchell 2020


% Method options
max_niter_div = 1000 ;
max_niter_rot = 1000 ;
niterU2d_div = 0 ;
niterU2d_rot = 0 ;
niterU2d_harm = 0 ;
clipDiv = [-Inf, Inf] ;  
clipRot = [-Inf, Inf] ;  
lambda_smooth = 0 ;
lambda_mesh = 0 ;
nmodes = 0 ;
zwidth = 0 ;
method = 'smooth' ;     % options: smooth, denoise, both
eps = 1e-16 ;
preview = false ;
do_calibration = false ;
computeLaplacian = false ;

%% Unpack options
if isfield(Options, 'lambda')
    lambda_smooth = Options.lambda ;
end
if isfield(Options, 'lambda_mesh')
    lambda_mesh = Options.lambda_mesh ;
end
if isfield(Options, 'nSpectralFilterModes')
    nmodes = Options.nSpectralFilterModes ;
end
if isfield(Options, 'spectralFilterWidth')
    zwidth = Options.spectralFilterWidth;
end
if isfield(Options, 'outdir')
    outdir = Options.outdir;
else
    outdir = '' ;
end
if isfield(Options, 'do_calibration')
    do_calibration = Options.do_calibration;
end
if isfield(Options, 'computeLaplacian')
    computeLaplacian = Options.computeLaplacian;
end

%% varargin options
for i = 1:length(varargin)
    
    if isa(varargin{i}, 'double')
        continue;
    end
    if isa(varargin{i}, 'logical')
        continue;
    end
        
    if ~isempty(regexp(varargin{i}, '^[Nn]iter[Ss]moothing', 'match'))
        niter = varargin{i+1} ;
        % Allow for different niters for divergence and rotation
        if length(niter) > 1
            max_niter_mesh = niter(1) ;
            max_niter_div = niter(2) ;
            max_niter_rot = niter(3) ;
        else
            max_niter_mesh = niter ;
            max_niter_div = niter ;
            max_niter_rot = niter ;
        end
    end    
    
    if ~isempty(regexp(varargin{i}, '^[Nn]iter[Ss]moothing[Uu]2[Dd]', 'match'))
        niter = varargin{i+1} ;
        % Allow for different niters for divergence and rotation
        if length(niter) > 1
            niterU2d_div = niter(1) ;
            niterU2d_rot = niter(2) ;
            niterU2d_harm = niter(3) ;
        else
            niterU2d_div = niter ;
            niterU2d_rot = niter ;
            niterU2d_harm = niter ;
        end
    end    
    
    if ~isempty(regexp(varargin{i}, '^[Cc]lip[Dd]iv', 'match'))
        clipDiv = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Cc]lip[Rr]ot', 'match'))
        clipRot = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Mm]ethod', 'match'))
        method = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Ee]psilon', 'match'))
        eps = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Pp]review', 'match'))
        preview = varargin{i+1} ;
    end    
end

% Unpack the cutMesh
FF = cutM.f ;
V2D = cutM.u ;
vtx3drs = cutM.v ;
nU = cutM.nU ;
nV = cutM.nV ;

% % compute COM for each triangle in 2d and 3d --> note that we do this
% % before gluing so that the 2d case is not messed up, and the 3d case is
% % the same either way
% bc = cat( 3, v3drs(FF(:,1), :), v3drs(FF(:,2), :), v3drs(FF(:,3), :) );
% bc = mean( bc, 3 ) ;
% bc2d = cat( 3, V2D(FF(:,1), :), V2D(FF(:,2), :), V2D(FF(:,3), :) );
% bc2d = mean( bc2d, 3 ) ;

% num faces in each row, col is nfU, nfV
% faces are arranged as (nU-1)*(nV-1) * 2
% [~, faceIDgrid] = defineFacesRectilinearGrid(V2D, nU, nV) ;

% % Check that this makes sense
% trisurf(FF, v3drs(:, 1), v3drs(:, 2), v3drs(:, 3), divv, ...
%     'FaceAlpha', 0.2, 'EdgeColor', 'none')

% hold on;
% inds = 1:21:length(bc) ;
% quiver3(bc(inds, 1), bc(inds, 2), bc(inds, 3), ...
%     vfsm(inds, 1), vfsm(inds, 2), vfsm(inds, 3), 1) 
% % quiver3(bc(inds, 1), bc(inds, 2), bc(inds, 3), ...
% %     v0t(inds, 1), v0t(inds, 2), v0t(inds, 3), 1) 
% axis equal

% Take divergence and curl
        
%% Glue the mesh back together, FF will change
% cutMC is cutM that is Closed at the seam
[glueMesh, glue2cut] = glueCylinderCutMeshSeam(cutM) ;  

% Check smoothing
if preview
    clf
    subplot(1, 2, 1)
    trisurf(glueMesh.f, glueMesh.v(:, 1), glueMesh.v(:, 2), ...
        glueMesh.v(:, 3), 'edgecolor', 'none')
    axis equal
    view(0, 0)
    title('before smoothing')
end

% Laplacian smooth mesh vertices (lightly)
triglued = triangulation(glueMesh.f, glueMesh.v) ;
% Fix free boundaries of glued mesh, let others vary in space
fixed_verts = triglued.freeBoundary ; 
fixed_verts = fixed_verts(:, 1) ;
glueMesh.v = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
    lambda_mesh, 'explicit', glueMesh.v, max_niter_mesh) ;

if preview
    subplot(1, 2, 2)
    trisurf(glueMesh.f, glueMesh.v(:, 1), glueMesh.v(:, 2), ...
        glueMesh.v(:, 3), 'edgecolor', 'none')
    axis equal
    view(0, 0)
    title('after smoothing')
    pause(1)
    close all
end

%% Create DEC instance
DEC = DiscreteExteriorCalculus( glueMesh.f, glueMesh.v ) ;

% Now resolve the vector field for decomposition
[v0n, v0t, v0t2d, jac3d_to_2d, ~, ~, dilation] = ...
    resolveTangentNormalVelocities(FF, vtx3drs, facevf, 1:length(FF), V2D ) ;

divv = DEC.divergence(v0t) ;
rotv = DEC.curl(v0t) ;
[~, gF2V] = meshAveragingOperators(glueMesh.f, glueMesh.v) ;

% Optionally, compute laplacian too
if computeLaplacian
    lapv = DEC.laplacian(v0t) ;
else
    lapv = [] ;
end

% % Check that this gives the same thing as laplacian dim by dim
% checks = 0*lapv ;
% for dim = 1:3
%     checks(:, dim) = DEC.laplacian(gF2V * v0t(:, dim)) ;
% end

% Clip the fields using supplied bounds for valid values
divv(divv < clipDiv(1)) = clipDiv(1) ;
divv(divv > clipDiv(2)) = clipDiv(2) ;
rotv(rotv < clipRot(1)) = clipRot(1) ;
rotv(rotv > clipRot(2)) = clipRot(2) ;

% Perform Helmholtz-Hodge decomposition
[divU, rotU, harmU, scalarP, vectorP] = ...
    DEC.helmholtzHodgeDecomposition(facevf) ;

% Preview divergence
if preview 
    tmp = reshape(divv, [nU, (nV-1)]) ;
    imagesc(tmp)
    colormap bwr 
    colorbar
    caxis([-0.4, 0.4])
    title('Divergence, before denoising/smoothing')
    pause(1)
end

%% First pass: do calibration for good params
if do_calibration && false
    nmodes2explore = 1:10 ;
    lambdas2explore = linspace(0, 0.005, 10) ;
    errs = zeros(length(nmodes2explore), length(lambdas2explore)) ;
    errm = errs ;
    % cut off endcaps
    divvCut = reshape(divv, [nU, nV-1]) ;
    divvCut = divvCut(3:end-3, :) ;

    for ii = 1:length(nmodes2explore)
        for jj = 1:length(lambdas2explore)
            nmodesij = nmodes2explore(ii) ;
            lambda_smooth_ij = lambdas2explore(jj) ;
            if nmodesij > 0
                disp('Mode filtering divv on vertices')
                filterOpts.nmodesY = nmodesij ;
                filterOpts.widthX = 1 ;
                divvsm = modeFilterQuasi1D(reshape(divv, [nU, (nV-1)]), filterOpts) ;
            else
                divvsm = divv ;
            end

            if lambda_smooth_ij > 0 
                disp('Laplacian smoothing divv on vertices')
                fixed_verts = [] ;  % note: could use boundaries here, seems unnecessary
                divvsm = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
                    lambda_smooth_ij, 'explicit', divvsm(:), max_niter_div) ;
            end

            % cut off endcaps
            divvsm = reshape(divvsm, [nU, nV-1]) ;
            divvsm = divvsm(3:end-3, :) ;

            errs(ii, jj) = mean((divvsm(:) - divvCut(:)).^2) / var(divvCut(:)) ;
            errm(ii, jj) = max((divvsm(:) - divvCut(:)).^2) / max(abs(divvCut(:))) ;        
        end
    end
    subplot(1, 2, 1)
    imagesc(nmodes2explore, lambdas2explore, errs')
    xlabel('\#modes', 'interpreter', 'latex')
    ylabel('$\lambda$', 'interpreter', 'latex')
    set(gca, 'YDir', 'normal');
    colorbar('location', 'southoutside')
    title('$\frac{\langle\Delta\nabla \cdot v_t \rangle}{\sigma^2_{\nabla \cdot v_t}}$', 'interpreter', 'latex')
    subplot(1, 2, 2)
    imagesc(nmodes2explore, lambdas2explore, errm')
    xlabel('\#modes', 'interpreter', 'latex')
    ylabel('$\lambda$', 'interpreter', 'latex')
    set(gca, 'YDir', 'normal');
    colorbar('location', 'southoutside')
    title('$\frac{\mathrm{max}\Delta\nabla \cdot v_t}{\mathrm{max}\nabla \cdot v_t}$', 'interpreter', 'latex')
    sgtitle('Smoothing analysis for divergence', 'interpreter', 'latex')
    saveas(gcf, fullfile(outdir, 'smoothing_error_analysis.png'))

    %% Change width of filter in longitude
    nmodes2explore = 1:10 ;
    zw2explore = 0:5 ;
    errs = zeros(length(nmodes2explore), length(zw2explore)) ;
    errm = errs ;
    % cut off endcaps
    divvCut = reshape(divv, [nU, nV-1]) ;
    divvCut = divvCut(3:end-3, :) ;
    for ii = 1:length(nmodes2explore)
        for jj = 1:length(zw2explore)
            nmij = nmodes2explore(ii) ;
            zwij = zw2explore(jj) ; 
            if nmodes > 0
                disp('Mode filtering divv on vertices')
                filterOpts.nmodesY = nmij ;
                filterOpts.widthX = zwij ;
                divvsm = modeFilterQuasi1D(reshape(divv, [nU, (nV-1)]), filterOpts) ;
            else
                divvsm = divv ;
            end

            % cut off endcaps
            divvsm = reshape(divvsm, [nU, nV-1]) ;
            divvsm = divvsm(3:end-3, :) ;

            errs(ii, jj) = mean((divvsm(:) - divvCut(:)).^2) / var(divvCut(:)) ;
            errm(ii, jj) = max((divvsm(:) - divvCut(:)).^2) / max(abs(divvCut(:))) ;        
        end
    end
    subplot(1, 2, 1)
    imagesc(nmodes2explore, lambdas2explore, errs')
    xlabel('\#modes', 'interpreter', 'latex')
    ylabel('$w$', 'interpreter', 'latex')
    set(gca, 'YDir', 'normal');
    colorbar('location', 'southoutside')
    title('$\frac{\langle\Delta\nabla \cdot v_t \rangle}{\sigma^2_{\nabla \cdot v_t}}$', 'interpreter', 'latex')
    subplot(1, 2, 2)
    imagesc(nmodes2explore, lambdas2explore, errm')
    xlabel('\#modes', 'interpreter', 'latex')
    ylabel('$w$', 'interpreter', 'latex')
    set(gca, 'YDir', 'normal');
    colorbar('location', 'southoutside')
    title('$\frac{\mathrm{max}\Delta\nabla \cdot v_t}{\mathrm{max}\nabla \cdot v_t}$', 'interpreter', 'latex')
    sgtitle('Smoothing analysis for divergence', 'interpreter', 'latex')
    saveas(gcf, fullfile(outdir, 'smoothing_error_analysis2.png'))
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAPLACIAN SMOOTHING for divergence field (on vertices)
if nmodes > 0 || zwidth > 0
    disp('Mode filtering divv and lapv on vertices')
    filterOpts.nmodesY = nmodes ;
    filterOpts.widthX = zwidth ;
    divvm = modeFilterQuasi1D(reshape(divv, [nU, (nV-1)]), filterOpts) ;
    
    if computeLaplacian
        lapvm = zeros(nU*(nV-1), 3) ;
        tmpX = modeFilterQuasi1D(reshape(lapv(:, 1), [nU, (nV-1)]), filterOpts) ;
        tmpY = modeFilterQuasi1D(reshape(lapv(:, 2), [nU, (nV-1)]), filterOpts) ;
        tmpZ = modeFilterQuasi1D(reshape(lapv(:, 3), [nU, (nV-1)]), filterOpts) ;
        lapvm(:, 1) = tmpX(:) ;
        lapvm(:, 2) = tmpY(:) ;
        lapvm(:, 3) = tmpZ(:) ;
    else
        lapvm = lapv ;
    end
else
    divvm = divv ;
    lapvm = lapv ;
end

if lambda_smooth > 0 
    disp('Laplacian smoothing divv on vertices')
    fixed_verts = [] ;  % note: could use boundaries here, seems unnecessary
    divvsm = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
        lambda_smooth, 'explicit', divvm(:), max_niter_div) ;
    lapvsm = lapvm ;
    if computeLaplacian
        lapvsm(:, 1) = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
            lambda_smooth, 'explicit', lapvm(:, 1), max_niter_div) ;
        lapvsm(:, 2) = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
            lambda_smooth, 'explicit', lapvm(:, 2), max_niter_div) ;
        lapvsm(:, 3) = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
            lambda_smooth, 'explicit', lapvm(:, 3), max_niter_div) ;
    end
else
    divvsm = divvm(:) ;
    if computeLaplacian
        lapvsm = reshape(lapvm, [nU*(nV-1), 3]) ;
    else
        lapvsm = [] ;
    end
end

% View results on divergence
if preview 
    clf
    clim = max(abs(divv)) ;
    tmp = reshape(divv, [nU, (nV-1)]) ;
    subplot(2, 2, 1)
    imagesc(tmp')
    colormap bwr 
    colorbar
    caxis([-clim, clim]); axis off ;
    title('Divergence')
    tmp = reshape(divvm, [nU, (nV-1)]) ;
    subplot(2, 2, 2)
    imagesc(tmp')
    colormap bwr 
    colorbar
    caxis([-clim, clim]); axis off ;
    title('after filtering')
    subplot(2, 2, 3)
    tmp = reshape(divvsm, [nU, (nV-1)]) ;
    imagesc(tmp')
    colormap bwr 
    colorbar
    caxis([-clim, clim]); axis off ;
    title('after filtering + smoothing')
    % Compare to no spectral filter at all
    subplot(2, 2, 4)
    tmp = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...  
        lambda_smooth, 'explicit', divv(:), max_niter_div) ;
    tmp = reshape(tmp, [nU, (nV-1)]) ;
    imagesc(tmp')
    colormap bwr 
    colorbar
    caxis([-clim, clim]); axis off ;
    title('just smoothing')
    sgtitle(['\#modes=' num2str(nmodes) ', $w$=' num2str(zwidth) ', $\lambda$=' ...
        num2str(lambda_smooth) ' $\lambda_m$=' num2str(lambda_mesh)], ...
        'interpreter', 'latex') ;
    pause(1)
end

%% SMOOTHING for rotational field
% weight curl by face area
% Note: hodge dual (star) is applied at end of Curl(), so no area exists in
% the curl/rotational field, so just average with any weighting we like (we
% think).
% See retired code at end for old face smoothing. Now we put the field onto
% vertices.

% First inspect the field
if preview 
    tmp = reshape(rotv, [nU-1, (nV-1)*2]) ;
    imagesc(tmp)
    colormap bwr 
    colorbar
    caxis([-clim, clim])
    title('Curl, before denoising/smoothing')
    pause(1)
end

%% SMOOTHING on vertices (rotation field)
fixed_verts = [] ;  % note: could use boundaries here, seems unnecessary
rotvsm = gF2V * rotv ;
% Spectral filtering
if nmodes > 0 
    disp('Spectral filtering rotv on vertices')
    filterOpts.nmodesY = nmodes ;
    filterOpts.widthX = zwidth ;
    rotvsm = modeFilterQuasi1D(reshape(rotvsm, [nU, (nV-1)]), filterOpts) ;
end
% Average using laplacian smoothing
if lambda_smooth > 0 
    disp('Laplacian smoothing rotv on vertices')
    rotvsm = laplacian_smooth(glueMesh.v, glueMesh.f, 'uniform', fixed_verts, ...
        lambda_smooth, 'implicit', rotvsm(:), max_niter_rot) ;
end
rotvsm = rotvsm(:) ;

% add nU points back to divv from phi=0 by duplication --> convert back to
% cut mesh indices
divvCut = divvsm(glue2cut) ;
% note: rot was on faces, and faces are preserved, but now on vertices
rotvCut = rotvsm(glue2cut) ; 
if computeLaplacian
    lapvCut = lapvsm(glue2cut, :) ;
else
    lapvCut = lapvsm ;
end

% Note: rotU, divU, and harmU are on faces, and faces are preserved when
% converting from glueMesh to cutMesh

% Pullback Vector Fields to Domain of Parameterization ====================
divU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    divU2d(f,:) = jac3d_to_2d{FF(f)} * divU(f,:)' / dilation(f) ;
end
rotU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    rotU2d(f,:) = jac3d_to_2d{FF(f)} * rotU(f,:)'  / dilation(f) ;
end
harmU2d = zeros(size(FF, 1), 2);
for f = 1:size(FF,1)
    harmU2d(f,:) = jac3d_to_2d{FF(f)} * harmU(f,:)' / dilation(f) ;
end

% Smooth output 2d fields
for dim=1:2
    divU2d(:, dim) = laplacian_smooth_faces(divU2d(:, dim), glueMesh, ...
        'Weight', 'Area', 'Method', method, 'Epsilon', eps, ...
        'niter', niterU2d_div) ;
    rotU2d(:, dim) = laplacian_smooth_faces(rotU2d(:, dim), glueMesh, ...
        'Weight', 'Area', 'Method', method, 'Epsilon', eps, ...
        'niter', niterU2d_rot) ;
    harmU2d(:, dim) = laplacian_smooth_faces(harmU2d(:, dim), glueMesh, ...
        'Weight', 'Area', 'Method', method, 'Epsilon', eps, ...
        'niter', niterU2d_harm) ;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Store fields in structs
divs.raw = divv ;
divs.divv = divvCut ;
divs.divU = divU ;
divs.divU2d = divU2d ;
rots.raw = rotv ;
rots.rotv = rotvCut ;
rots.rotU = rotU ;
rots.rotU2d = rotU2d ;
harms.harmU = harmU ;
harms.harmU2d = harmU2d ;
lapvs.raw = lapv ;
lapvs.lapv = lapvCut ;

return


% Retired code
% 
% % SMOOTHING rot on faces -- DENOISE, SMOOTH, or ITERATE BOTH
% if strcmp(method, 'denoise')
%     % apply harmonic smoothing several times to extrema (DENOISE)
%     for qq = 1:niter_rot
%         rotvsm = laplacian_smooth_faces(rotvsm, glueMesh, ...
%             'Weight', 'Area', 'Method', 'denoise', 'Epsilon', eps) ;
% 
%         if preview 
%             tmp = reshape(rotvsm, [nU-1, (nV-1)*2]) ;
%             imagesc(tmp)
%             colormap bwr 
%             colorbar
%             caxis([-0.4, 0.4])
%             title(['Curl, denoising #' num2str(qq)])
%             pause(3)
%         end
%     end
% elseif strcmp(method, 'smooth')
%     % apply harmonic smoothing several times (SMOOTHING)
%     for qq = 1:niter_rot
%         rotvsm = laplacian_smooth_faces(rotvsm, glueMesh,...
%             'Weight', 'Area', 'Method', 'smooth', 'Epsilon', eps) ;
% 
%         % Check it
%         if preview
%             tmp = reshape(rotvsm, [nU-1, (nV-1)*2]) ;
%             imagesc(tmp)
%             colormap bwr 
%             colorbar
%             caxis([-0.4, 0.4])
%             title(['Curl, smoothing #' num2str(qq)])
%             pause(1)
%         end
%     end
% elseif strcmp(method, 'both')
%     % apply denoising and harmonic smoothing several times (BOTH)
%     for qq = 1:niter_rot
%         rotvsm = laplacian_smooth_faces(rotvsm, glueMesh,...
%             'Weight', 'Area', 'Method', 'denoise', 'Epsilon', eps) ;
% 
%         % Check it
%         if preview
%             tmp = reshape(rotvsm, [nU-1, (nV-1)*2]) ;
%             imagesc(tmp)
%             colormap bwr 
%             colorbar
%             caxis([-0.4, 0.4])
%             title(['Curl, denoising #' num2str(qq)])
%             pause(1)
%         end
%         
%         rotvsm = laplacian_smooth_faces(rotvsm, glueMesh,...
%             'Weight', 'Area', 'Method', 'smooth', 'Epsilon', eps) ;
%         
%         % Check it
%         if preview
%             tmp = reshape(rotvsm, [nU-1, (nV-1)*2]) ;
%             imagesc(tmp)
%             colormap bwr 
%             colorbar
%             caxis([-0.4, 0.4])
%             title(['Curl, smoothing #' num2str(qq)])
%             pause(1)
%         end
%         
%     end
% else
%     error('SmoothingMethod not recognized')
% end


