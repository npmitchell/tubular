function generateUVPrimeCutMeshes(QS, options)
% Generate pullbacks in the (u',v') coordinate system
%
% Parameters
% ----------
% QS : 
% options : struct with fields
%
%
% Returns
% -------
% none
%
% NPMitchell 2020


%% Unpack QS
timePoints = QS.xp.fileMeta.timePoints ;

%% Unpack options
overwrite = true ;
save_ims = false ;
preview = false ;
if isfield(options, 'overwrite')
    overwrite = options.overwrite ; 
end
if isfield(options, 'save_ims')
    save_ims = options.save_ims ; 
end
if isfield(options, 'preview')
    preview = options.preview ; 
end
if save_ims
    imdir = fullfile(QS.dir.uvpcutMesh, 'images') ; 
    mkdir(imdir) ;
end
t0 = QS.t0set() ;


% First load all phi offsets phi0(u=0,1) from spcutMesh
for tidx = 1:length(timePoints)
    tp = timePoints(tidx) ;
    fn = sprintf(QS.fullFileBase.uvpcutMesh, tp) ;
    
    if ~exist(fn, 'file') || overwrite
        QS.setTime(tp) ;    

        % % OPTION 1: Could load spcutMesh.phi0_fit, but this would be before 
        % % smoothing. 
        % % To do so, first load cylMesh for this timepoint to pull back via 
        % conformal map
        % QS.loadCurrentCylinderMeshClean() ;
        % cylMesh = QS.currentMesh.cylinderMeshClean ;
        % % Extract offsets phi0 and phi1 to use for boundaries
        % phi0 = QS.cutMesh.phi0_fit(1) ; phi1 = QS.cutMesh.phi0_fit(end) ;
        % % Could then average phi0s in time by hand.
        % % Prep filter
        % disp('Building tripulse filter equivalent to tripuls(-0.5:0.1:0.5)')
        % tripulse = 0:0.2:1 ;
        % tripulse = [tripulse, fliplr(tripulse(1:end-1))] ;
        % tripulse = tripulse ./ sum(tripulse(:)) ;
        % tripulse = reshape(tripulse, [length(tripulse), 1]) ;

        % % OPTION 2:  Compute phi0 anew
        % uvals_boundary = [1e-9, max(spcutMeshSm.u(:, 1))] ;
        % tileCount = [1,1] ;
        % [ TF, TV2D, TV3D ] = tileAnnularCutMesh( spcutMeshSm, tileCount );
        % phiv_kk = [0 0] ;
        % prev3d_sphi_dsdv = mesh0.v ;
        % save_ims = false ;
        % plotfn = 'none' ;
        % smoothingMethod = 'none' ;
        % phiOpts = struct() ;
        % patchOpts = struct() ;
        % [phi0_fit_kk, phi0s_kk] = fitPhiOffsetsFromPrevMesh(TF,...
        %     TV2D, TV3D, uvals_boundary, phiv_kk, ...
        %     prev3d_sphi_dsdv, -0.45, 0.45, ...
        %     save_ims, plotfn, smoothingMethod, ...
        %     phiOpts, patchOpts) ;

        %% OPTION 3: Load spcutmesh, which already has incorporated phi(u=0,1)
        % Then we only need to conformally map the triangles to a rectilinear
        % domain and voila!
        QS.loadCurrentSPCutMeshSm() ;
        spcutMeshSm = QS.currentMesh.spcutMeshSm ;
        
        % Conformally map to disk
        tileCount = [1 1] ;
        
        rawMesh = flattenAnnulus(spcutMeshSm) ;   
        % Relax Affine transformation along x
        rawMesh.ar = minimizeIsoarealAffineEnergy( rawMesh.f, rawMesh.v, rawMesh.u );
        
        if QS.flipy
            rawMesh.u = [ rawMesh.u(:,1), 1.0 - rawMesh.u(:,2) ];
        end
        
        % Compute Beltrami coefficient (quasiconformal)
        % affine_factor = rawMesh.ar ;
        rawMesh.mu = bc_metric(rawMesh.f, rawMesh.u, rawMesh.v, 3) ;
        rawMesh.ar_mu = (1 + mean(real(rawMesh.mu))) / (1 - mean(real(rawMesh.mu))) ;
        uvpcutMesh.raw = rawMesh ;
        
        % Check that aspect ratio ar relaxes mu
        vtxtmp = rawMesh.u ;
        vtxtmp(:, 1) = vtxtmp(:, 1) * rawMesh.ar_mu ;
        tmp = mean(real(bc_metric(rawMesh.f, vtxtmp, rawMesh.v, 3))) ;
        disp(['Residual quasiconformal component: ' num2str(tmp)])
        
        %% Check that ar minimizes mu -- IT DOES NOT!
        % dmyk = 1 ;
        % ars = linspace(rawMesh.ar-0.1, rawMesh.ar+0.5, 10) ;
        % meanmu = 0*ars ;
        % stdmu = 0*ars ;
        % for affine_factor = ars
        %     v2d_raw = rawMesh.u ;
        %     v2d_raw(:, 1) = affine_factor * v2d_raw(:, 1) ;
        %     mus = bc_metric(rawMesh.f, v2d_raw, rawMesh.v, 3) ;
        %     disp([num2str(mean(real(mus))) '+/-' num2str(std(real(mus)))])
        %     meanmu(dmyk) = mean(abs(mus)) ;
        %     stdmu(dmyk) = std(abs(mus)) ;
        %     dmyk = dmyk + 1 ;
        % end
        % clf
        % h1 = errorbar(ars, meanmu, stdmu); 
        % hold on;
        % h2 = errorbar(rawMesh.ar, mean(real(rawMesh.mu)), std(real(rawMesh.mu)), 'o')
        % xlabel('possible affine factors')
        % ylabel('$\langle \mu \rangle \pm \sigma_\mu$ across mesh', ...
        %     'interpreter', 'latex')
        % legend('test values', 'relaxed value')
        
        
        % Check against previously saved
        % tmp = load(sprintf(QS.fullFileBase.uvpcutMesh, tp)) ;
        % all(all(tmp.uvpcutMesh.raw.u == rawMesh.u))
        % all(all(tmp.uvpcutMesh.raw.v == rawMesh.v))

        %% Make avgpts in pixel space (not RS), from rawMesh
        fprintf('Resampling uvgrid3d curves in pix...\n')
        nU = rawMesh.nU ;
        nV = rawMesh.nV ;
        curves3d_pix = reshape(rawMesh.v, [nU, nV, 3]) ;
        c3d_dsv_pix = zeros(size(curves3d_pix)) ;  % in units of pix
        avgpts_pix = zeros(nU, 3) ;
        radius_pix = zeros(nU, nV) ;
        for i=1:nU
            % Note: no need to add the first point to the curve
            % since the endpoints already match exactly in 3d and
            % curvspace gives a curve with points on either
            % endpoint (corresponding to the same 3d location).
            c3d_dsv_pix(i, :, :) = resampleCurvReplaceNaNs(squeeze(curves3d_pix(i, :, :)), nV, true) ;
            if vecnorm(squeeze(c3d_dsv_pix(i, 1, :)) - squeeze(c3d_dsv_pix(i, end, :))) > 1e-7
                error('endpoints do not join! Exiting')
            end
            % Drop the final endpoint in the mean pt determination
            avgpts_pix(i, :) = mean(squeeze(c3d_dsv_pix(i, 1:end-1, :)), 1) ; 
            radius_pix(i, :) = vecnorm(squeeze(curves3d_pix(i, :, :)) - avgpts_pix(i, :), 2, 2) ;
        end
        
        % uvpcutMesh.raw.avgpts_pix = avgpts_pix ;
        % uvpcutMesh.raw.radius_pix = radius_pix ;
        uvpcutMesh.raw.avgpts_um = QS.xyz2APDV(avgpts_pix) ;
        uvpcutMesh.raw.radius_um = radius_pix * QS.APDV.resolution ;

        %% Compute rectified rectilinear grid in (u', v') pullback space   
        % First tile the raw mesh along v'
        tileCount = [1,1] ;
        [ TF, TV2D, TV3D ] = tileAnnularCutMesh( uvpcutMesh.raw, tileCount );
        tmp = uvpcutMesh.raw.radius_um(:) ;
        assert(length(tmp) == nU*nV)
        tiled_radius_um = [tmp(1:nU*(nV-1)); tmp; tmp(nU+1:end)] ;
        
        onesUV = ones(nU, nV) ;
        uspace = linspace(0, 1, nU) ;
        vspace = linspace(0, 1, nV) ;
        uu = (uspace .* onesUV)' ;
        vv = vspace .* onesUV' ;
        % uv is equally spaced in u' and v'
        uv = [uu(:), vv(:)] ;
        uvgrid3d = interpolate2Dpts_3Dmesh(TF, TV2D, TV3D, uv) ;
        resMesh.nU = nU ;
        resMesh.nV = nV ;
        resMesh.f = rawMesh.f ;
        resMesh.u = uv ;
        resMesh.v = uvgrid3d ;
        resMesh.ar = rawMesh.ar ;
        resMesh.pathPairs = rawMesh.pathPairs ;
        resMesh.mu = bc_metric(resMesh.f, resMesh.u, resMesh.v, 3) ;
        uvpcutMesh.resampled = resMesh ;
        assert(all(~any(isnan(resMesh.v))))
        assert(all(~any(isnan(resMesh.u))))

        % Interpolate radius onto resampled grid
        ugrid = reshape(uv(:, 1), [nU, nV]) ;
        radInterp = scatteredInterpolant(TV2D(:, 1), TV2D(:, 2),...
            tiled_radius_um(:), 'linear', 'nearest') ;
        resRadius_um = radInterp(resMesh.u(:, 1), resMesh.u(:, 2)) ;
        resRadius_um = reshape(resRadius_um, [nU, nV]) ;
        uvpcutMesh.resampled.radius_pix = resRadius_um ;
        
        if preview
            close all
            scatter(rawMesh.u(:, 1), rawMesh.u(:, 2), 10, radius_um(:)) 
            hold on;
            scatter(resMesh.u(:, 1), resMesh.u(:, 2), 10, resRadius_um(:))
            cb = colorbar() ;
            ylabel(cb, ['radius [' QS.spaceUnits ']'], 'interpreter', 'latex')
            pause(1)
        end
                
        %% Save the result
        disp(['Saving uvpcutMesh t=' num2str(tp) ': ' fn])
        save(fn, 'uvpcutMesh')
    else
        computed = false ;
    end

    %% Check orientation
    imfn = [sprintf(QS.fileBase.uvpcutMesh, tp) '.png'] ;
    imfn_full = fullfile(imdir, imfn) ;
    imfn_raw = ['mu_raw_' sprintf(QS.fileBase.uvpcutMesh, tp) '.png'] ;
    imfn_rawFull = fullfile(imdir, imfn_raw) ;
    imfn_raw3d = ['mu_raw3d_' sprintf(QS.fileBase.uvpcutMesh, tp) '.png'] ;
    imfn_raw3dFull = fullfile(imdir, imfn_raw3d) ;
    imfn_res = ['mu_res_' sprintf(QS.fileBase.uvpcutMesh, tp) '.png'] ;
    imfn_resFull = fullfile(imdir, imfn_res) ;
    not_on_disk = ~exist(imfn_full, 'file') || ...
        ~exist(imfn_rawFull, 'file') || ...
        ~exist(imfn_raw3dFull, 'file') || ...
        ~exist(imfn_resFull, 'file') ;
    
    if save_ims && not_on_disk
        if ~computed
            load(fn, 'uvpcutMesh')
            rawMesh = uvpcutMesh.raw ;
            resMesh = uvpcutMesh.resampled ;
            uspace = linspace(0, 1, rawMesh.nU) ;
            vspace = linspace(0, 1, rawMesh.nV) ;
            resRadius_um = uvpcutMesh.resampled.radius_pix ;
        end
        
        close all
        if ~exist(imfn_full, 'file')
            set(gcf, 'visible', 'off')
            imagesc(uspace, vspace, resRadius_um') ;
            title('Conformal $(u'',v'')$ cutMesh', 'interpreter', 'latex')
            axis equal
            axis tight
            xlabel('$u$', 'interpreter', 'latex')
            ylabel('$v$', 'interpreter', 'latex')
            h = gca() ;
            set(h,'Xcolor','none')
            h.XAxis.Label.Color=[0 0 0];
            h.XAxis.Label.Visible='on';
            set(h,'Ycolor','none')
            h.YAxis.Label.Color=[0 0 0];
            h.YAxis.Label.Visible='on';
            cb = colorbar() ;
            ylabel(cb, ['radius [' QS.spaceUnits ']'], 'interpreter', 'latex') 
            saveas(gcf, imfn_full)
            close all
        end

        % Save quasiconformal mesh -- raw2d
        if ~exist(imfn_rawFull, 'file')
            set(gcf, 'visible', 'off')
            subplot(1, 2, 1)
            options.labels = {'$\Re \mu$', '$\Im \mu$'} ;
            [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                twoScalarFieldsOnSurface({rawMesh.f, ...
                [rawMesh.u(:, 1), rawMesh.u(:, 2), 0*rawMesh.u(:, 1)]}, ...
                real(rawMesh.mu) - mean(real(rawMesh.mu)), ...
                imag(rawMesh.mu), options) ;
            sgtitle(['$\mu($embedding, pullback$)$, $t = $', ...
                sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                'interpreter', 'latex') ;
            set(gcf,'CurrentAxes', ax1)
            view(2)
            axis off
            set(gcf,'CurrentAxes', ax2)
            view(2)
            axis off
            disp(['Saving ' imfn_raw ': ' fullfile(imdir, imfn_raw)])
            saveas(gcf, imfn_rawFull)
            close all
        end

        % Save quasiconformal mesh -- raw3d
        if ~exist(imfn_raw3dFull, 'file')
            v3draw = QS.xyz2APDV(rawMesh.v) ;
            set(gcf, 'visible', 'off')
            subplot(1, 2, 1)
            options.labels = {'$\Re \mu$', '$\Im \mu$'} ;
            [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                twoScalarFieldsOnSurface({rawMesh.f, ...
                v3draw}, ...
                real(rawMesh.mu) - mean(real(rawMesh.mu)), ...
                imag(rawMesh.mu), options) ;
            sgtitle(['$\mu($embedding, pullback$)$, $t = $', ...
                sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                'interpreter', 'latex') ;
            set(gcf,'CurrentAxes', ax1)
            axis equal
            axis off
            set(gcf,'CurrentAxes', ax2)
            axis equal
            axis off
            disp(['Saving ' imfn_raw3d ': ' fullfile(imdir, imfn_raw3d)])
            saveas(gcf, imfn_raw3dFull)
            close all
        end

        % Save quasiconformal mesh -- resampled
        if ~exist(imfn_resFull, 'file')
            set(gcf, 'visible', 'off')
            subplot(1, 2, 1)
            options.labels = {'$\Re \mu$', '$\Im \mu$'} ;
            [ax1, ax2, cb1, cb2, mesh1, mesh2] = ...
                twoScalarFieldsOnSurface({resMesh.f, ...
                [resMesh.u(:, 1), resMesh.u(:, 2), 0*resMesh.u(:, 1)]}, ...
                real(resMesh.mu) - mean(real(resMesh.mu)), ...
                imag(resMesh.mu), options) ;
            sgtitle(['resampled mesh: $\mu($embedding, pullback$)$, $t = $', ...
                sprintf('%03d', tp-t0), ' ', QS.timeUnits], ...
                'interpreter', 'latex') ;
            set(gcf,'CurrentAxes', ax1)
            view(2)
            axis off
            set(gcf,'CurrentAxes', ax2)
            view(2)
            axis off
            saveas(gcf, imfn_resFull)
            close all
        end
    end
end
    
    
% Notes to self:
% uvpcutMesh = flattenAnnulusTiltedBoundaries(cutMesh, phi0, phi1, 'Dirichlet') ;
% ar = minimizeIsoarealAffineEnergy( cutMesh.f, cutMesh.v, cutMesh.u );
