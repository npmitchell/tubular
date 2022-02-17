function aux_smooth_avgptcline_radius_before_mesh_smoothing(...
    overwrite_spcutMesh_smoothradii, timePoints, ...
    spcutMeshBase, clineDVhoopBase, ...
    radiusImDir, rot, trans, resolution, xyzlim, nU, nV)
%AUX_SMOOTH_AVGPTCLINE_RADIUS_BEFORE_MESH_SMOOTHING
%
% Parameters
% ----------
% timePoints :
%   same as xp.fileMeta.timePoints
%
% Returns
% -------
% 
% NPMitchell 2020

redo_smoothed_centerline_and_radius = overwrite_spcutMesh_smoothradii ;
kk = 1;
while ~redo_smoothed_centerline_and_radius && kk < (length(timePoints) + 1)
    t = timePoints(kk) ;
    disp(['Checking if cline smoothing needs to be done for t=' num2str(t)])
    % Load mcline
    tmp = load(sprintf(clineDVhoopBase, t)) ;
    % Add radii from smoothed mcline to spcutMesh if not already present
    load(sprintf(spcutMeshBase, t), 'spcutMesh') ;

    radii_in_cutMeshsm = isfield(spcutMesh, 'radii_from_smoothed_mcline') ;
    avgpts_projected_in_cline = isfield(tmp, 'avgpts_projected') ;
    if ~radii_in_cutMeshsm  || ~avgpts_projected_in_cline
        redo_smoothed_centerline_and_radius = true;
    end
    kk = kk + 1 ;
end

% If any data is missing on disk, or if overwrite is true, redo
if redo_smoothed_centerline_and_radius
    mclineM = zeros(length(timePoints), 1000, 3) ;
    for kk=1:length(timePoints)
        t = timePoints(kk) ;
        % Load mcline
        load(sprintf(clineDVhoopBase, t), 'mcline') ;
        mclineM(kk, :, :) = mcline ;
    end

    % Time-average the centerline
    mcline_smM = movmean(mclineM,5,1) ;

    % Optional: polynomial fit to smooth

    % Re-compute radii from smoothed centerlines by pointmatching avgpts onto
    % smoothed centerline and plot them
    radImDir0 = fullfile(radiusImDir, 'view0') ;
    radImDirD = fullfile(radiusImDir, 'viewD') ;
    radImDirV = fullfile(radiusImDir, 'viewV') ;
    radImDirA = fullfile(radiusImDir, 'viewA') ;
    radImDirP = fullfile(radiusImDir, 'viewP') ;
    radImDirL = fullfile(radiusImDir, 'viewL') ;
    radImDirR = fullfile(radiusImDir, 'viewR') ;
    radImDirs = fullfile(radiusImDir, 'sphi') ;
    radImDiru = fullfile(radiusImDir, 'uphi') ;
    dirs2do = {radImDir0, radImDirD, radImDirV, ...
        radImDirA, radImDirP, radImDirL, radImDirR, ...
        radImDirs, radImDiru} ;
    for qq = 1:length(dirs2do)
        if ~exist(dirs2do{qq}, 'dir')
            mkdir(dirs2do{qq})
        end
    end
    % Now compute and plot
    for kk=1:length(timePoints)
        t = timePoints(kk) ;
        disp(['Computing radii from smoothed centerline for t=' num2str(t)])
        % Load mcline
        load(sprintf(clineDVhoopBase, t), 'avgpts', 'mcline', 'mss') ;
        mcline_sm = squeeze(mcline_smM(kk, :, :)) ;
        mss_sm = ss_from_xyz(mcline_sm) ;
        inds = pointMatch(avgpts, mcline_sm) ;
        avgpts_projected = mcline_sm(inds, :) ;

        % Save the new results with the old
        save(sprintf(clineDVhoopBase, t), 'avgpts', 'mcline', 'mss',...
            'avgpts_projected', 'mcline_sm', 'mss_sm') ;

        % Add radii from smoothed mcline to spcutMesh if not already present
        load(sprintf(spcutMeshBase, t), 'spcutMesh') ;
        vrs = ((rot * spcutMesh.v')' + trans) * resolution ;
        radii_in_cutMeshsm = isfield(spcutMesh, 'radii_from_smoothed_mcline') ;
        if ~radii_in_cutMeshsm || overwrite_spcutMesh_smoothradii
            curvesDV = reshape(vrs, [nU, nV, 3]) ;
            radii_from_avgpts_sm = zeros(size(curvesDV, 1), size(curvesDV, 2)) ;
            for jj = 1:size(curvesDV, 1)
                % Consider this hoop
                hoop = squeeze(curvesDV(jj, :, :)) ;
                radii_from_avgpts_sm(jj, :) = vecnorm(hoop - avgpts_projected(jj, :), 2, 2) ;
            end
            % Save it in spcutMesh
            spcutMesh.radii_from_smoothed_mcline = radii_from_avgpts_sm ;
            save(sprintf(spcutMeshBase, t), 'spcutMesh') ;
        else
            radii_from_avgpts_sm = spcutMesh.radii_from_smoothed_mcline ;
        end

        % Plot the radii as 3D image
        rad2dfigfn_s = fullfile(radImDirs, sprintf('radius_dvhoop_sphi_%06d.png', t)) ;
        % No longer use uphi as a map: commented out
        % rad2dfigfn_u = fullfile(radImDiru, sprintf('radius_dvhoop_uphi_%06d.png', t)) ;
        rad3dfigfn = fullfile(radImDir0, sprintf('radius_dvhoop_xyz_%06d.png', t)) ;
        if ~exist(rad3dfigfn, 'file') || ~exist(rad2dfigfn_s, 'file') || ...
                overwrite_spcutMesh_smoothradii
            if ~exist(rad3dfigfn, 'file') || ~exist(rad2dfigfn_s, 'file') 
                disp('Writing figures for radii from cline_sm')
            else
                disp('OVERWRITING figures for radii from cline_sm')
            end
            close all
            fig = figure('visible', 'off') ;
            trisurf(spcutMesh.f, vrs(:, 1), vrs(:, 2), ...
                vrs(:, 3), radii_from_avgpts_sm(:), ...
                'Edgecolor', 'none')
            c = colorbar ;
            c.Label.String = 'radius [\mum]' ;
            xlabel('x [\mum]')
            ylabel('y [\mum]')
            zlabel('z [\mum]')
            axis equal
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            title('Radius via DV curves and smoothed centerline')
            saveas(fig, rad3dfigfn)
            title(['Radius via DV curves, Dorsal view [' radImDirA ']'])
            view(0, 90)
            fn = sprintf('radius_dvhoop_xyD_%06d.png', t) ;
            disp(['Saving ' fn])
            saveas(gcf, fullfile(radImDirD, fn)) 
            title(['Radius via DV curves, Ventral view [' radImDirA ']'])
            view(0, -90)
            fn = sprintf('radius_dvhoop_xyV_%06d.png', t) ;
            disp(['Saving ' fn])
            saveas(gcf, fullfile(radImDirV, fn)) 
            title(['Radius via DV curves, Posterior view [' radImDirA ']'])
            view(90, 0)
            fn = sprintf('radius_dvhoop_yzP_%06d.png', t) ;
            disp(['Saving ' fn])
            saveas(gcf, fullfile(radImDirP, fn)) 
            title(['Radius via DV curves, Anterior view [' radImDirA ']'])
            view(-90, 0)
            fn = sprintf('radius_dvhoop_yzA_%06d.png', t) ;
            disp(['Saving ' fn])
            saveas(gcf, fullfile(radImDirA, fn)) 
            title(['Radius via DV curves, Lateral view [' radImDirA ']'])
            view(0, 0)
            fn = sprintf('radius_dvhoop_xzL_%06d.png', t) ;
            disp(['Saving ' fn])
            saveas(gcf, fullfile(radImDirL, fn)) 
            title(['Radius via DV curves, Lateral view [' radImDirA ']'])
            view(0, 180)
            fn = sprintf('radius_dvhoop_xzR_%06d.png', t) ;
            disp(['Saving ' fn])
            saveas(gcf, fullfile(radImDirR, fn)) 

            % Plot the radii as 2D image
            tmp = {spcutMesh.sphi} ;  %, spcutMesh.uphi} ;
            rad2dfigfns = {rad2dfigfn_s}; % , rad2dfigfn_u} ;
            xlabels = {'AP position, s/L', 'AP position, u/L'} ;
            for qq=1:length(tmp)
                uu = tmp{qq} ;
                if qq == 1
                    uu(:, 1) = uu(:, 1) / max(uu(:, 1)) ;
                else
                    uu(:, 1) = 0.5 * uu(:, 1) ;
                end
                close all
                fig = figure('visible', 'off') ;
                trisurf(spcutMesh.f, uu(:, 1), uu(:, 2), ...
                    radii_from_avgpts_sm(:), 'EdgeColor', 'none')
                caxis([0, 80]) ;
                % also plot tiled meshes above and below
                hold on;
                trisurf(spcutMesh.f, uu(:, 1), uu(:, 2) + 1, ...
                    radii_from_avgpts_sm(:), 'Edgecolor', 'none')
                trisurf(spcutMesh.f, uu(:, 1), uu(:, 2) - 1, ...
                    radii_from_avgpts_sm(:), 'Edgecolor', 'none')
                c = colorbar ;
                c.Label.String = 'radius [\mum]' ;
                xlabel(xlabels{qq})
                ylabel('\phi/2\pi')
                ylim([0, 1])
                view(2)
                disp(['Saving 3D image to ' rad2dfigfns{qq}])
                saveas(fig, rad2dfigfns{qq})
                close all
            end
        end
    end
end
disp('done')

% note: see extract_radius_from_DVhoops for old version of radius plotting
