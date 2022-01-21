function measurePIV3dMultiChannel(QS, options)
% measurePIV3dMultiChannel(QS, options)
%   Measure 3d flows from PIV results on smoothed (s,phi) MIP pullbacks
%   for each channel
%   
%   default options for created MIPS (done before this method) are:
%     pbOptions = struct() ;
%     pbOptions.overwrite = false ;
%     pbOptions.numLayers = [5, 5] ;
%     pbOptions.layerSpacing = 0.75 ;
%     pbOptions.generate_rsm = true ;
%     pbOptions.generate_spsm = true ;
%     pbOptions.generate_sphi = false ;
%     QS.data.adjustlow = 1.00 ;
%     QS.data.adjusthigh = 99.999 ;
%     QS.generateCurrentPullbacks([], [], [], pbOptions) ;
%
% 
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields 
%   overwrite : bool
%       overwrite previous results
%   preview : bool
%       view intermediate results
%   timePoints : numeric 1D array
%       the timepoints to consider for the measurement. For ex, could
%       choose subset of the QS experiment timePoints
%
%
% Saves to disk
% -------------
% m0XY : Px2 float: 2d mesh vertices in pullback image pixel space [pullback pix[" ;
% m0f : #facesx3 float: mesh connectivity list" ;
% m0v3d : Px3 float: 3d mesh coordinates in embedding pixel space [mesh pix / dt]" ;
% x0 : QxR float: x value of 2d velocity evaluation coordinates in pullback image pixel space [pullback pix]" ;
% y0 : QxR float: y value of 2d velocity evaluation coordinates in pullback image pixel space [pullback pix]" ;
% pt0 : Nx3 float: 3d location of evaluation points [mesh pix]" ;
% pt1 : Nx3 float: 3d location of advected point in next mesh [mesh pix]" ;
% v0 : Nx3 float: 3d velocity [mesh pix / dt]" ;
% v0n : Nx1 float: normal velocity [mesh pix / dt]" ;
% v0t : Nx3 float: tangential velocity in 3d [mesh pix / dt]" ;
% facenormals : #faces x 3 float: face normals for all faces [unitless, mesh pix / mesh pix]" ;
% v3dfaces : #faces x 3 float: face-centered velocities in embedding space [mesh pix / dt]" ;
% v3dvertices : #mesh vertices x 3 float: vertex-centered velocities in embedding space [mesh pix / dt]" ;
% v0_rs : Nx3 float: rotated/scaled 3d velocities [um/min]" ;
% v0n_rs : Nx3 float: scaled normal velocities [um/min]" ;
% v0t_rs : Nx3 float: rotated/scaled tangential 3d velocities [um/min]" ;
% normals_rs : Nx3 float: normal vector of field faces in 3d [unitless, um/um]" ;
% v0t2d : Nx2 float: in-plane velocity [pullback pix / dt]" ;
% g_ab : Nx2x2 float: pullback metric" ;
% dilation : Nx1 float: dilation of face from 3d to 2d (A_2d [pullback_pix^2] / A_3d [mesh_pix^2])" ;
% jac : #faces x 1 cell: jacobian of 3d->2d transformation" ;
% fieldfaces : Nx1 int: face indices into spcutMesh where v defined" ;
%
% NPMitchell 2020-2021

%% Default options
vnscale = 2 ;
vtscale = 5 ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;
end
if isfield(options, 'preview')
    preview = options.preview ;
else
    preview = false ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
else
    timePoints = QS.xp.fileMeta.timePoints ;
end
if isfield(options, 'doubleCovered')
    doubleCovered = options.doubleCovered ; 
else
    doubleCovered = true;
end
if isfield(options, 'pivimCoords')
    pivimCoords = options.pivimCoords ;
else
    pivimCoords = 'sp_sme' ;
end
if isfield(options, 'show_v3d_on_data')
    show_v3d_on_data = options.show_v3d_on_data ;
else
    show_v3d_on_data = true ;
end
if isfield(options, 'save_ims')
    save_ims = options.save_ims ;
else
    save_ims = true ;
end
if isfield(options, 'save_piv3d_im')
    save_piv3d_im = options.save_piv3d_im ;
else
    save_piv3d_im = false ;
end
if isfield(options, 'vtscale')
    vtscale = options.vtscale ;
end
if isfield(options, 'vnscale')
    vnscale = options.vnscale ;
end
if isfield(options, 'rgbChannels')
    % Indices of channels in pullback images to load
    rgbChannels = 1:length(QS.xp.expMeta.channelsUsed) ;
end

%% Unpack QS
piv3dfn = QS.fullFileBase.piv3dMultiChannel ;
ntps = length(timePoints) ;
[rot, ~] = QS.getRotTrans() ;
resolution = QS.APDV.resolution ; 
[~, ~, ~, xyzlim_APDV] = QS.getXYZLims() ;
axis_order = QS.data.axisOrder ;
pivOutDir = QS.dir.pivMultiChannel.v3d ;
blue = QS.plotting.colors(1, :) ;
red = QS.plotting.colors(2, :) ;
green = QS.plotting.colors(4, :) ;
t0 = QS.t0set() ;

%% Check if all timepoints' piv3d exist already
redo_piv3d = overwrite ; 
ii = 1 ;
while ~redo_piv3d && ii < length(timePoints)
    for channel = QS.xp.expMeta.channelsUsed
        redo_piv3d = redo_piv3d || ...
            ~exist(sprintf(piv3dfn, ch, timePoints(ii)), 'file') ;
        if redo_piv3d
            disp([sprintf(piv3dfn, ch, timePoints(ii))...
                ' has not been computed --> redo'])
        end
    end
    ii = ii + 1 ;
end

%% Compute or load results
if ~redo_piv3d
    disp(['Loading piv3d from ' piv3dfn])
    piv3d = cell(ntps - 1, 1) ;
    for ii=1:(ntps-1)
        if mod(timePoints(ii), 10) == 0
            disp(['loading piv3d for time t=' num2str(timePoints(ii))])
        end
        load(sprintf(piv3dfn, timePoints(ii)), 'piv3dstruct')
        piv3d{ii} = piv3dstruct ;
    end
else
    disp('Overwriting/computing piv3d data...')
    
    % preallocate piv3d
    piv3d = cell(length(QS.xp.expMeta.channelsUsed), 1) ;

    disp('Loading raw PIV results')    
    for ch = QS.xp.expMeta.channelsUsed
        piv{ch} = load(sprintf(QS.fileName.pivRawMultiChannel.raw, ch)) ;
    end
    
    % Iterate over all images with flow fields ----------------------------
    for ii=1:(ntps - 1)    
        for qq = 1:length(QS.xp.expMeta.channelsUsed)
            ch = QS.xp.expMeta.channelsUsed(qq) ;
            
            % Get timing
            timestr = sprintf('%03d', timePoints(ii) - t0) ;
            disp(['t = ' timestr])
            tp = timePoints(ii) ;
            dt = timePoints(ii + 1) - tp ;

            % Get scale of image
            if strcmp(pivimCoords, 'sp_sme') || strcmp(pivimCoords, 'spsme') 
                im0 = imread(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
                im1 = imread(sprintf(QS.fullFileBase.im_sp_sme, timePoints(ii+1))) ;
                doubleCovered = true ;
            else
                error(['Unrecognized pivimCoords: ' pivimCoords])
            end
            Xsz0 = size(im0, 2) ;
            Ysz0 = size(im0, 1) ;  
            Ysz1 = size(im1, 1) ;

            % Load spcutMesh 
            mesh0 = load(sprintf(QS.fullFileBase.spcutMeshSm, ...
                timePoints(ii)), 'spcutMeshSm') ;
            mesh0 = mesh0.spcutMeshSm ;
            umax0 = max(mesh0.u(:, 1)) ;
            vmax0 = max(mesh0.u(:, 2)) ;

            % Load next timepoint's spcutMesh
            mesh1 = load(sprintf(QS.fullFileBase.spcutMeshSm, ...
                timePoints(ii + 1)), 'spcutMeshSm') ;
            mesh1 = mesh1.spcutMeshSm ;
            umax1 = max(mesh1.u(:, 1)) ;
            vmax1 = max(mesh1.u(:, 2)) ;

            m0XY = QS.uv2XY(im0, mesh0.u, doubleCovered, umax0, vmax0) ;
            % Note: no need for m1XY since we only use tiled mesh (tm1XY)
            % m1XY = QS.uv2XY(im1, mesh1.u, doubleCovered, umax1, vmax1) ;
            m0x = m0XY(:, 1) ;
            m0y = m0XY(:, 2) ;

            % Load the positions of the velocity vectors in pixels
            x0 = piv{ch}.x{ii} ;
            y0 = piv{ch}.y{ii} ;
            uu = piv{ch}.u_filtered{ii} ;
            vv = piv{ch}.v_filtered{ii} ; 

            % Ensure no NaNs in uu and vv
            if any(isnan(uu(:))) || any(isnan(vv(:)))
               disp('inpainting NaNs in uu & vv')
               uu = inpaint_nans(uu) ;
               vv = inpaint_nans(vv) ;
            end

            % Get position in next timepoint in pixels (in XY pixel plane)
            % Clip the x position to the size of the image, and wrap the y position
            x1 = x0 + uu ;
            eps = 1e-8 ;
            x1 = max(x1, 1.0 + eps) ;  % the minimum value of tm1XY(:, 1) is 1.0
            x1 = min(x1, Xsz0 - eps) ;
            y1 = mod(y0 + vv, Ysz1) ;  

            % Two methods of obtaining the 3d vel evaluation pts. 
            % Barycenteric version:
            % [pt0, pt1, tr0, tr0_orig, tr1, tr1_orig, tria] = aux_barycentricInterpTiledMesh(mesh0, mesh1,...
            %     Ysc0, umax0, Ysc1, xsc1, use_shifty, shifty, i, x0, y0, x1, y1, x2Xpix, y2Ypix, dx2dX, dy2dY) ;
            % mf0c = tr0.ConnectivityList ;
            % m0xy = tr0.Points ;
            % m0x = tr0_orig.Points(:, 1) ;
            % m0y = tr0_orig.Points(:, 2) ;
            % mesh1x = tr1_orig.Points(:, 1) ;
            % mesh1y = tr1_orig.Points(:, 2) ;

            % % Instead, do interpolation:
            disp('Interpolating 3d vertices for tiled mesh 0')
            tileCount = [2, 2] ;
            if any(isnan(mesh0.v))
                % fill missing x,y,z values in mesh0.v
                disp('Filling missing x,y,z values in mesh0.v')
                for dim = 1:3
                    dgrid = reshape(mesh0.v(:, dim), [mesh0.nU, mesh0.nV]) ;
                    dgrid = fillmissing(dgrid, 'linear') ;
                    mesh0.v(:, dim) = dgrid(:) ;
                end
                clearvars dgrid
            end

            % Note: tm0f   will be [2*(nU-1)*(nV-1)*5] x 3
            % Note: tm0v2d will be [nU*(nV-1)*4 + nU*nV] x 2
            % Note: pt0    will be [#X0 * #Y0] x 3
            [tm0f, tm0v2d, tm0v3d, ~] = tileAnnularCutMesh(mesh0, tileCount);
            tm0XY = QS.uv2XY(im0, tm0v2d, doubleCovered, umax0, vmax0) ;
            [pt0, fieldfaces, tr0] = interpolate2Dpts_3Dmesh(tm0f, tm0XY,...
                                            tm0v3d, [x0(:), y0(:)]) ;

            % note: for interpolation of velocity onto face centers, use
            mesh_for_interp.pathPairs = mesh0.pathPairs ;
            mesh_for_interp.u = m0XY ;
            mesh_for_interp.f = mesh0.f ;

            % If any fieldfaces are NaN, jitter the points a bit
            if any(isnan(fieldfaces))
                disp('Field faces are nan. Jittering positions to fix (hack)...')
                % Fill in missing data. Will this work?
                facegrid = reshape(fieldfaces, size(x0)) ;
                facegrid = fillmissing(facegrid, 'linear') ;
                fieldfaces = uint8(facegrid(:)) ;

                if preview
                    % which are the bad ones?
                    badID = find(isnan(fieldfaces)) ;

                    % check it in a plot
                    triplot(tr0.ConnectivityList, tr0.Points(:, 1), tr0.Points(:, 2), 0*tr0.Points(:, 2), 'Edgecolor', 'none')
                    hold on;
                    bd = freeBoundary(tr0) ;
                    plot3(tr0.Points(bd(:, 1), 1), tr0.Points(bd(:, 1), 2), 0*tr0.Points(bd(:,1), 1), '.')
                    hold on;
                    plot3(x0(badID), y0(badID), 0*x0(badID), 'o')
                    xlim([min(x0(badID)) - 1, min(x0(badID)) + 1])
                    title('NaNs in fieldfaces')
                    waitfor(fig)
                end
            end

            % Interpolate for next timepoint
            disp('Interpolating 3d vertices for tiled mesh 1')
            tileCount = [2, 2] ;
            [tm1f, tm1v2d, tm1v3d, ~] = tileAnnularCutMesh(mesh1, tileCount);
            tm1XY = QS.uv2XY(im1, tm1v2d, doubleCovered, umax1, vmax1) ;
            pt1 = interpolate2Dpts_3Dmesh(tm1f, tm1XY, tm1v3d, [x1(:), y1(:)]) ;

            % Ensure no NaNs in pt0 and pt1
            if any(isnan(pt0(:))) || any(isnan(pt1(:)))
               % disp('inpainting NaNs in pt0 & pt1')
               error('why nans?')
               % pt0 = inpaint_nans(pt0) ;
               % pt1 = inpaint_nans(pt1) ;
               close all
               figure ;
               scatter(x1(:), y1(:), 10, pt1(:, 1))
               bad = find(isnan(pt1(:, 1))) ;
               hold on; 
               xx1 = x1(:) ;
               yy1 = y1(:) ;
               scatter(xx1(bad), yy1(bad), 10, 'k')
            end

            % Old version
            % Xbi = scatteredInterpolant(tm1X, tm1Y, tm1v3d(:, 1)) ;
            % Ybi = scatteredInterpolant(tm1X, tm1Y, tm1v3d(:, 2)) ;
            % Zbi = scatteredInterpolant(tm1X, tm1Y, tm1v3d(:, 3)) ;
            % pt1 = [Xbi(x1(:), y1(:)), Ybi(x1(:), y1(:)), Zbi(x1(:), y1(:))] ;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (preview || save_ims) && show_v3d_on_data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot the advected mesh
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Load 3D data for coloring mesh pullback
                disp(['loading timepoint data t=' num2str(tp)])
                QS.setTime(tp) ;
                QS.getCurrentData() ;
                IV0 = QS.currentData.IV ;

                % also plot the next timepoint
                t1 = timePoints(ii + 1) ;
                disp(['loading subsequent timepoint data t=' num2str(t1)])
                QS.setTime(t1);
                QS.getCurrentData() ;
                IV1 = QS.currentData.IV ;            

                disp('previewing')
                % evaluate flow field at (x, Ysz - y). Add (u,v) to (x, Ysz - y).
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot the advected mesh
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot the advected tissue on top of the next frame
                % interpolate velocities at (x0, y0) onto mesh (m0XY)
                uuinterp = griddedInterpolant(x0', y0', uu', 'linear', 'nearest') ; % in inverted Y
                vvinterp = griddedInterpolant(x0', y0', vv', 'linear', 'nearest') ; 
                addx = uuinterp(m0x, m0y) ; % interpolate in correct Y
                addy = vvinterp(m0x, m0y) ;
                mesh0adv_pix = [ m0x + addx, m0y + addy] ;

                % Texture image options
                Options.imSize = ceil( Xsz0 .* [(1+sum(tileCount))/QS.a_fixed, 1] );  % row, col pix
                Options.yLim = [0 Ysz0];
                Options.xLim = [0 Xsz0];
                % original image RED
                % im0 = texture_patch_to_image(tm0f, [tm0XY(:, 1), tm0XY(:, 2)], ...
                %     tm0f, tm0v3d(:, axis_order), IV0, Options );
                % advected image GREEN
                patchIma = texture_patch_to_image(mesh0.f, mesh0adv_pix, ...
                    mesh0.f, mesh0.v(:, axis_order), IV0, Options );
                % patchIma = adapthisteq(patchIma, 'NumTiles', [round(a_fixed * ntiles), round(2 * ntiles)]) ;
                % Next timepoint BLUE
                % im1 = texture_patch_to_image(tm1f, [tm1XY(:, 1), tm1XY(:, 2)], ...
                %     tm1f, tm1v3d(:, axis_order), IV1, Options );
                % Load next timepoint BLUE
                % im1 = imread(fullfile(fns(i+1).folder, fns(i+1).name)) ;
                % patchImRGB = cat(3, im0, uint8(patchIma * 255), im1) ;
                minrow = (tileCount(1) - 0.5)*Xsz0 / QS.a_fixed ;
                patchImaCrop = patchIma(minrow+1:minrow+Ysz0, :) ;
                patchImRGB = cat(3, squeeze(mean(im0, 3)), ...
                    patchImaCrop, im1) ;

                % Plot the image and advected image, in inverted Y space but
                % with increasing Y going down
                close all; 
                fig1 = figure(1) ; 
                h = imshow( patchImRGB );
                hold on;  
                quiver(m0x, m0y, addx, addy, 0, 'color', yellow, 'linewidth', 2)
                % plot(m0x, m0y, 'o')
                % plot(mesh0adv_pix(:, 1), mesh0adv_pix(:, 2), 's')
                triplot(mesh0.f, m0x, m0y, 'color', red, 'linewidth', 2)
                triplot(mesh0.f, m0x + addx, m0y + addy, 'color', green, 'linewidth', 2)
                % axis equal
                title('(r:t0) (g:advected t0) (b:t1) (y:PIV)')
                waitfor(gcf)

                % Check the displacement by toggling between frames
                % fig = figure(2) ;
                % pressed_enter = false ; 
                % kk = 0 ;
                % while ~pressed_enter
                %     if kk == 0
                %         imshow(flipud(im0))
                %         %set( gca, 'YDir', 'Normal' );
                %         titlestr = '0: <Enter> to exit, <-> to toggle' ;
                %     else
                %         imshow(flipud(im1))
                %         %set( gca, 'YDir', 'Normal' );
                %         titlestr = '1: <Enter> to exit, <-> to toggle' ;
                %     end
                %     title(titlestr)
                %     was_a_key = waitforbuttonpress;
                %     left = strcmp(get(fig, 'CurrentKey'), 'leftarrow');
                %     rght = strcmp(get(fig, 'CurrentKey'), 'rightarrow') ;
                %     if was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
                %         pressed_enter = true ;
                %     elseif was_a_key && left || rght 
                %         kk = mod(kk + 1, 2) ;
                %     end
                % end

                % Check differences
                d0 = mat2gray(im1, [0, double(max(im1(:)))]) - mat2gray(im0, [0, double(max(im0(:)))]);
                d0pos = 0*d0 ;
                d0pos(d0 > 0) = d0(d0 > 0) ;
                d0neg = 0*d0 ;
                d0neg(d0 < 0) = abs(d0(d0 < 0)) ;
                pd0 = cat(3, d0pos, 0*d0, d0neg) ; 

                % diff between next and advected 
                da = mat2gray(im1, [0, double(max(im1(:)))]) - patchIma ;
                dapos = 0*da ;
                dapos(da > 0) = da(da > 0) ;
                daneg = 0*da ;
                daneg(da < 0) = abs(da(da < 0)) ;
                pda = cat(3, dapos, 0*da, daneg) ;
                % plot both
                fig1 = figure(1) ;
                imshow(pd0) 
                title(['\langle|t1 - t0|\rangle = ' num2str(100 * mean(pd0(:))) '%'])
                waitfor(fig1)

                fig2 = figure(2) ;
                imshow(pda) 
                title(['\langle|t1 - advected t0|\rangle = ' num2str(100 * mean(pda(:))) '%'])
                waitfor(fig2)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Take the difference to get the velocity field ------------------
            v0 = (pt1 - pt0) / dt ;

            % Decompose/Resolve tangential and normal velocities ==============
            % Also obtain jacobian.
            [v0n, v0t, v0t2d, jac, facenormals, g_ab, dilation] = ...
                resolveTangentNormalVelocities(tm0f, tm0v3d, v0, fieldfaces, tm0XY) ;

            % Check for NaNs
            assert(~any(isnan(v0t2d(:)))) 

            % check this
            % imagesc(x0(:, 1), y0(1, :), reshape(v0n, size(x0)))
            % error('checking here -- erase this check if no debug')
            % trisurf(triangulation(tm0f, tm0v3d), 'edgecolor', 'none')
            % hold on;
            % plot3(pt1(:, 1), pt1(:, 2), pt1(:, 3)) ;

            % Ensure no NaNs in uu and vv
            if any(isnan(v0t(:))) || any(isnan(v0n(:)))
                error('why do we have NaNs in v0t?')
                disp('inpainting NaNs in v0t & v0n')
                v0t = inpaint_nans(reshape(v0t, [])) ;
                v0n = inpaint_nans(reshape(v0n, [])) ;
            end

            % I have checked that det(jac * jac') = det(jjac * jjac') for a
            % triangle. 

            % todo: check that the dilation = ratio of areas of triangles2d /
            % triangles3d

            if preview
                % Checking the dilation
                % plot(tm0XY(:, 1), tm0XY(:, 2), '.')
                clf
                hold on;
                scatter(x0(:), y0(:), 10, dilation, 'filled')
                set(gcf, 'visible', 'on')
                pause(2)
                clf

                % Now independently measure the areas of triangles in 3d
                % Evaluate for every face in tm0
                fa3d = doublearea(tm0v3d, tm0f) * 0.5 ;
                % also do 2d areas
                fa2d = doublearea(tm0XY, tm0f) * 0.5 ;
                arearatio = fa2d(fieldfaces) ./ fa3d(fieldfaces) ;
                figure;
                plot(dilation, arearatio, '.') ;
                hold on; 
                plot([0, 5], [0, 5], '--')
                xlabel('dilation from jacobian')
                ylabel('area ratio of triangles')
                title('Checking dilation from 3d->2d: should be y=x')
                set(gcf, 'visible', 'on')
                pause(2)
                clf
            end

            % If the velocity sampling is finer than 3 vectors on each face,
            % do one averaging pass before resolving onto faces
            [bincounts] = histc(fieldfaces, min(fieldfaces)-1:max(fieldfaces)+1) ;
            if mean(bincounts(bincounts > 0)) > 20
                msg = 'Velocity is much finer than mesh by factor of ';
                msg = [msg num2str(mean(bincounts(bincounts > 0))) ] ;
                msg = [msg '-- Do averaging step here. Have not needed to do this yet.'] ;
                error(msg) 
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Assign one velocity per face via interpolation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Obtain circumcenters for the current timepoint's 2D mesh
                % tr0 = triangulation(tm0f, tm0XY) ;  <-- already done in interpolate2dpts_3Dmesh
                % bc = circumcenter(tr0) ;
                % bc = barycenter(tm0XY, tm0f) ; % this is in gptoolbox

                % Calculate centroid of faces -- here evaluate only on untiled mesh
                bc = cat( 3, mesh_for_interp.u(mesh_for_interp.f(:,1), :), ...
                    mesh_for_interp.u(mesh_for_interp.f(:,2), :),...
                    mesh_for_interp.u(mesh_for_interp.f(:,3), :) );
                bc = mean( bc, 3 );

                % cast as NxMx3, then as N*M x 1 arrays for vx,vy,vz separately
                v3dgrid = reshape(v0, [size(piv{ch}.x{1}, 1), size(piv{ch}.x{1}, 2), 3]) ;
                xvel = squeeze(v3dgrid(:, :, 1)) ;
                yvel = squeeze(v3dgrid(:, :, 2)) ;
                zvel = squeeze(v3dgrid(:, :, 3)) ;
                % Find the velocity at the face barycenter using only 2d coords
                Fx = griddedInterpolant(x0', y0', xvel', 'linear', 'nearest') ;
                Fy = griddedInterpolant(x0', y0', yvel', 'linear', 'nearest') ;
                Fz = griddedInterpolant(x0', y0', zvel', 'linear', 'nearest') ;
                v3dfaces = [Fx(bc(:, 1), bc(:, 2)), ...
                    Fy(bc(:, 1), bc(:, 2)), Fz(bc(:, 1), bc(:, 2))] ;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Assign one tangential velocity per face via interpolation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                v3dvertices = [Fx(mesh_for_interp.u(:, 1), mesh_for_interp.u(:, 2)), ...
                    Fy(mesh_for_interp.u(:, 1), mesh_for_interp.u(:, 2)), ...
                    Fz(mesh_for_interp.u(:, 1), mesh_for_interp.u(:, 2))] ;
                % %% Convert pivSimAvg to evaluated at vertices of simple smoothed meshes
                % piv3dfn = fullfile(fullfile(pivDir, 'piv3d'), 'piv3d_%04d.mat') ;
                % % load piv3d as a cell array
                % disp(['Loading piv3d from ' piv3dfn])
                % piv3d = cell(ntps) ;
                % for i=1:(ntps-1)
                %     if mod(time(i), 10) == 0
                %         disp(['loading piv3d for time t=' num2str(time(i))])
                %     end
                %     load(sprintf(piv3dfn, time(i)), 'piv3dstruct')
                %     piv3d{i} = piv3dstruct ;
                % end
            end

            % Test validity of result
            v0_rs =  QS.dx2APDV(v0) / dt ;
            if any(isnan(v0_rs(:))) || any(isnan(v0_rs(:)))
               % disp('inpainting NaNs in pt0 & pt1')
               error('why nans?')
               % pt0 = inpaint_nans(pt0) ;
               % pt1 = inpaint_nans(pt1) ;
               close all
               figure ;
               scatter(x1(:), y1(:), 10, v0_rs(:, 1))
               bad = find(isnan(v0_rs(:, 1))) ;
               hold on; 
               xx1 = x1(:) ;
               yy1 = y1(:) ;
               scatter(xx1(bad), yy1(bad), 10, 'k')
            end

            %% Save the results in datstruct ----------------------------------
            % v0, v0n, v0t are in units of um/min,  
            % while v0t2d, g_ab, and jacobian are in pullback pixels
            datstruct.m0XY = m0XY ;
            datstruct.m0f = mesh0.f ;
            datstruct.m0v3d = mesh0.v ;
            datstruct.x0 = x0 ;
            datstruct.y0 = y0 ;
            datstruct.pt0 = pt0 ;
            datstruct.pt1 = pt1 ;
            datstruct.v0 = v0 / dt ;
            datstruct.v0n = v0n / dt ;
            datstruct.v0t = v0t / dt ;
            datstruct.facenormals = facenormals ;

            % Face-centered velocities
            datstruct.v3dvertices = v3dvertices ;
            datstruct.v3dfaces = v3dfaces ;
            datstruct.v3dfaces_rs = QS.dx2APDV(v3dfaces) / dt ; 
            assert(all(size(v3dfaces) == size(datstruct.m0f)))

            % rotated and scaled velocities
            datstruct.v0_rs = QS.dx2APDV(v0) / dt ;
            datstruct.v0t_rs = QS.dx2APDV(v0t) / dt ;
            if QS.flipy
                % Direct normals inward
                datstruct.v0n_rs = -v0n * resolution / dt ;
                normals_rs = (rot * facenormals')' ;            
                % Put the normal direction as inward -- revert X,Z
                datstruct.normals_rs = normals_rs ;
                datstruct.normals_rs(:, 1) = -normals_rs(:, 1) ;
                datstruct.normals_rs(:, 3) = -normals_rs(:, 3) ;            
            else
                datstruct.v0n_rs = v0n * resolution / dt ;
                datstruct.normals_rs = (rot * facenormals')' ;
            end

            inds = 1:10:length(x0(:)) ;
            xx = x0(:) ;
            yy = y0(:) ;
            h3 = quiver(xx(inds), yy(inds), ...
                v0t2d(inds, 1), v0t2d(inds,2), 1, 'k', 'LineWidth', 1.2) ;

            % 2d pullback velocities
            datstruct.v0t2d = v0t2d ;           % in pullback pix / dt
            datstruct.g_ab = g_ab ;             % pullback metric
            datstruct.dilation = dilation ;     % dilation of face from 3d to 2d (2d area/3d area)
            datstruct.dilation_rs = dilation / resolution ;     % dilation of face from 3d to 2d
            datstruct.jacobian = jac ;          % jacobian of 3d->2d transformation
            datstruct.fieldfaces = fieldfaces ; % faces where v defined from PIV

            % Check for nans        
            try
                assert(~any(isnan(v0(:)))) 
                assert(~any(isnan(v0n(:)))) 
                assert(~any(isnan(v0t(:)))) 
                assert(~any(isnan(datstruct.m0v3d(:)))) 
                assert(~any(isnan(datstruct.v0t_rs(:)))) 
                assert(~any(isnan(datstruct.v0t2d(:))))
            catch
                error('Velocities contain NaNs. Debug here.')
            end

            % package the struct with a readme
            readme = struct ;
            readme.m0XY = "Px2 float: 2d mesh vertices in pullback image pixel space [pullback pix[" ;
            readme.m0f = "#facesx3 float: mesh connectivity list" ;
            readme.m0v3d = "Px3 float: 3d mesh coordinates in embedding pixel space [mesh pix / dt]" ;
            readme.x0 = "QxR float: x value of 2d velocity evaluation coordinates in pullback image pixel space [pullback pix]" ;
            readme.y0 = "QxR float: y value of 2d velocity evaluation coordinates in pullback image pixel space [pullback pix]" ;
            readme.pt0 = "Nx3 float: 3d location of evaluation points [mesh pix]" ;
            readme.pt1 = "Nx3 float: 3d location of advected point in next mesh [mesh pix]" ;
            readme.v0 = "Nx3 float: 3d velocity [mesh pix / dt]" ;
            readme.v0n = "Nx1 float: normal velocity [mesh pix / dt]" ;
            readme.v0t = "Nx3 float: tangential velocity in 3d [mesh pix / dt]" ;
            readme.facenormals = "#faces x 3 float: face normals for all faces [unitless, mesh pix / mesh pix]" ;
            readme.v3dfaces = "#faces x 3 float: face-centered velocities in embedding space (eval at face barycenters) [mesh pix / dt]" ;
            readme.v3dvertices = "#mesh vertices x 3 float: vertex-centered velocities in embedding space [mesh pix / dt]" ;
            readme.v0_rs = "Nx3 float: rotated/scaled 3d velocities at PIV evaluation points [um/min]" ;
            readme.v0n_rs = "Nx3 float: scaled normal velocities at PIV evaluation points [um/min]" ;
            readme.v0t_rs = "Nx3 float: rotated/scaled tangential 3d velocities at PIV evaluation points [um/min]" ;
            readme.normals_rs = "Nx3 float: normal vector of field faces in 3d (at PIV evaluation points) [unitless, um/um]" ;
            readme.v0t2d = "Nx2 float: in-plane velocity (at PIV evaluation points) [pullback pix / dt]" ;
            readme.g_ab = "Nx2x2 float: pullback metric" ;
            readme.dilation = "Nx1 float: dilation of face from 3d to 2d (A_2d [pullback_pix^2] / A_3d [mesh_pix^2])" ;
            readme.jac = "#faces x 1 cell: jacobian of 3d->2d transformation" ;
            readme.fieldfaces = "Nx1 int: face indices into spcutMesh where v defined" ;
            datstruct.readme = readme ;
            piv3d{qq}{ii} = datstruct ;
        end

        %% Visualize the flow in 3d --------------------------------------- 
        piv3dimagefn = fullfile(pivOutDir, ['piv3d_' sprintf('%06d', tp) '.png']) ;
        if save_ims && save_piv3d_im && (overwrite || ~exist(piv3dimagefn, 'file'))
            disp('visualize flow in 3d')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            close all
            fig = figure('Visible', 'Off') ;
            hold on            
            % for checking purposes, grab the first few indices
            m0 = mesh0 ;
            % rmIDx = 600:length(m0.v) ;
            % [m0.f, m0.v, oldvIdx] = remove_vertex_from_mesh(m0.f, m0.v, rmIDx) ;
            % m0.vn = m0.vn(oldvIdx, :) ;
            m1 = mesh1 ;
            
            if show_v3d_on_data
                % Instead of using texture_patch_3d
                % create griddedInterpolant 
                IV0interp = griddedInterpolant(single(32768 + 0.5 * IV0), 'cubic') ;
                v0i = IV0interp(m0.v(:, 2), m0.v(:, 1), m0.v(:, 3)) ;
                patch('Faces', m0.f, 'Vertices', m0.v, ...
                    'FaceVertexCData', v0i, 'FaceColor', 'interp', ...
                    'EdgeColor', 'none') ; 
                hold on;
                IV1i = griddedInterpolant(single(32768 - 0.5 * IV1)) ;
                v1i = IV1i(m1.v(:, 2), m1.v(:, 1), m1.v(:, 3)) ;
                patch('Faces', m1.f, 'Vertices', m1.v, ...
                    'FaceVertexCData', v1i, 'FaceColor', 'interp', ...
                    'EdgeColor', 'none') ;
                axis equal
                colormap(bwr)
                % Option 2:
                % least squares vertex baking : intead of evaluating vertices
                % in image volume and interpolating values, baking will
                % minimize an energy functional that makes the intensities
                % closer to the proper values from texture mapping
                % 

                % Option 3 (expensive): texture mapping
                % clearvars Options
                % Options.PSize = 2;
                % Options.EdgeColor = 'none';
                % % Options.Rotation = rot ;
                % % Options.Translation = trans ;
                % texture_patch_3d( m0.f, m0.v, ...
                %     m0.f, m0.v(:, [2 1 3]), 32768 + 0.5 * IV0, Options );
                % hold on            
                % texture_patch_3d( m1.f, m1.v, ...
                %     m1.f, m1.v(:, [2 1 3]), 32768 - 0.5 * IV1, Options );
                % hold on            
                % axis equal
                % nearby = false(size(pt0(:, 1), 1), 1) ;
                % for qq = 1:length(pt0(:, 1))
                %     if min(vecnorm(pt0(qq, :) - m0.v, 2, 2)) < 10
                %         nearby(qq) = true ;
                %     end
                % end
                % nearby = find(nearby) ;
                % quiver3(pt0(nearby, 1), pt0(nearby, 2), pt0(nearby, 3), ...
                %     v0(nearby, 1), v0(nearby, 2), v0(nearby, 3), 0, 'color', green) 
            else
                m0vr = QS.xyz2APDV(m0.v) ;
                m1vr = QS.xyz2APDV(m1.v) ;
                trisurf(m0.f, m0vr(:, 1), m0vr(:, 2), m0vr(:, 3), ...
                    'FaceColor', blue, 'EdgeColor', 'none', 'FaceAlpha', 0.1)
                hold on;
                trisurf(m1.f, m1vr(:, 1), m1vr(:, 2), m1vr(:, 3), ...
                    'FaceColor', red, 'EdgeColor', 'none', 'FaceAlpha', 0.1)
            end
            hold on      
            pt0r = QS.xyz2APDV(pt0) ;
            for ch = QS.xp.expMeta.channelsUsed
                quiver3(pt0r(:, 1), pt0r(:, 2), pt0r(:, 3), ...
                    datstruct{ch}.v0_rs(:, 1), ...
                    datstruct{ch}.v0_rs(:, 2), ...
                    datstruct{ch}.v0_rs(:, 3), 0)
            end
            
            % scatter3(pt0(:, 1), pt0(:, 2), pt0(:, 3))
            % plot3(pt0(:, 1), pt0(:, 2), pt0(:, 3), 'o', 'color', yellow)
            % plot3(pt1(:, 1), pt1(:, 2), pt1(:, 3), 's', 'color', yellow)
            axis equal
            title(['t=' timestr]) 
            xlabel('x')
            ylabel('y')
            zlabel('z')
            xlim(xyzlim_APDV(1, :))
            ylim(xyzlim_APDV(2, :))
            zlim(xyzlim_APDV(3, :))
            % plot3(m1.v(:, 1), m1.v(:, 2), m1.v(:, 3), '.')
            
            saveas(gcf, piv3dimagefn)
            close all
        end
        
        %% Draw 2D flows
        vtdir2d = QS.dir.pivMultiChannel.vt2d ;
        vndir2d = QS.dir.pivMultiChannel.vn2d ;
        if ~exist(vtdir2d, 'dir')
            mkdir(vtdir2d)
        end
        if ~exist(vndir2d, 'dir')
            mkdir(vndir2d)
        end
        
        % check if figure output exists already
        outimfn_n = fullfile(vndir2d, sprintf([QS.fileBase.name '.png'], tp)) ;
        outimfn_t = fullfile(vtdir2d, sprintf([QS.fileBase.name '.png'], tp)) ;
        if save_ims && (~exist(outimfn_n, 'file') || overwrite)
            washout2d = 0.5 ;
            close all; clear alpha 
            alphaVal = 0.7 ;
            
            figure('visible', 'off')
            for qq = 1:length(QS.xp.expMeta.channelsUsed)
                % set up axis
                subplot(length(QS.xp.expMeta.channelsUsed), 1, qq, ...
                    'visible', 'off')
                
                ch = QS.xp.expMeta.channelsUsed(qq) ;
                disp(['plotting v2d for ch=', num2str(ch)])
                
                % Load the image to put flow on top
                assert(length(size(im0)) == 3)  % convert to rgb for no cmap change
                
                chIm = squeeze(im0(:, :, rgbChannels(qq))) ;
                im = cat(3, chIm, chIm, chIm) ;
                im = im * washout2d + max(im(:)) * (1-washout2d) ;
                options.label = ['scaled tangential velocity, ', ...
                    '$|v_t|/||g^{-1}||$ [$\mu$m/', QS.timeUnits, ']'] ;
                options.xlim = [0, size(im, 2)] ;
                options.axPosition = [] ;
                if qq == 1
                    options.pbPosition = [] ;
                    options.cbPosition = [] ;
                    options.title = ['tangential velocity, $t=$', ...
                        sprintf('%03d', tp-t0), ' ', QS.timeUnits] ;
                else
                    options.title = '' ;
                end
                
                if doubleCovered
                    options.ylim = size(im, 1) * [0.25, 0.75] ;
                else
                    options.ylim = size(im, 1) * [0, 1] ;
                end
                options.fig = gcf ;
                if length(x0(1, :)) > 200
                    options.qsubsample = 10 ;
                    options.qscale = 10 ;
                else
                    options.qsubsample = 5 ;
                    options.qscale = 5 ;
                end
                % plot the 2d tangential flow, adjusted for dilation
                % Control for dilation in quiver
                v0t2dsc = piv3d{qq}{ii}.v0t2d ./ dilation * resolution ;
                xyf.x = x0 ; % (1,:)  ;
                xyf.y = y0 ; %  (:,1)' ;
                vectorFieldHeatPhaseOnImage(im, xyf, ...
                    v0t2dsc(:,1), v0t2dsc(:, 2), vtscale, options) ;
            end
            disp(['Saving figure: ', outimfn_t])
            saveas(gcf, outimfn_t) ;   
            close all
                           
            % Create the figure
            figure('visible', 'off')
            for qq = 1:length(QS.xp.expMeta.channelsUsed)
                % set up axis
                subplot(length(QS.xp.expMeta.channelsUsed), 1, qq, ...
                    'visible', 'off')
                
                % xlims = xlim ;
                % ylims = ylim ;
                % hold on
                % quiver(x0(:), y0(:), v0t2dsc(:, 1), v0t2dsc(:, 2), 0) ;
                % axis equal
                % xlim(xlims) ;
                % ylim(ylims) ;
                % % Extract image from figure axes
                % patchIm = getframe(gca);
                % outimfn = fullfile(vtdir2d, [fileName '.png']) ;
                % imwrite( patchIm.cdata, outimfn );

                %% Now draw normal flow as heatmap
                % image = im * washout2d + max(im) * (1-washout2d) ;
                % xlims = xlim ;
                % ylims = ylim ;
                % hold on
                % % normalflow = pcolor(x0, y0, reshape(v0n, size(x0)));
                % normalflow = scatter(x0(:), y0(:), 10, v0n, 'filled') ;
                % axis equal
                % colormap( normalflow, bwr);
                % colorbar( normalflow);
                % set( normalflow, 'AlphaData', alpha );   
                options.caxis = [-10 10] ;
                options.cmap = bwr ;
                options.alpha = 0.5 ;
                v0ngrid = reshape(datstruct.v0n_rs, size(x0)) ;

                % Option 1: heatmap on alpha image
                % --------------------------------
                % xfield = x0(1, :)' ;
                % yfield = y0(:, 1) ;
                % image = mat2gray(im, [0, 256]) ;
                % c_handle = heatmap_on_alphaimage(image, yfield, xfield, v0grid, options) ;
                % % Extract image from figure axes
                % patchIm = getframe(gca);
                % % Write figure to file
                % outimfn = fullfile(vndir2d, [fileName '.png']) ;
                % % print('-dpng','-r300', outimfn)
                % imwrite( patchIm.cdata, outimfn );

                % Option 2: rgb + overlay
                % --------------------------------
                
                chIm = squeeze(im0(:, :, rgbChannels(qq))) ;
                im = cat(3, chIm, chIm, chIm) ; 
                imshow(im) ; hold on;
                pcolor(x0, y0, v0ngrid)
                shading interp
                colormap(bwr)
                alpha(alphaVal)
                set(gca, 'clim', [-vnscale vnscale])
                ylim([500 1500])
                c = colorbar();
                % Manually flush the event queue and force MATLAB to render the colorbar
                % necessary on some versions
                drawnow
                % Get the color data of the object that correponds to the colorbar
                cdata = c.Face.Texture.CData;
                % Change the 4th channel (alpha channel) to 10% of it's initial value (255)
                cdata(end,:) = uint8(alphaVal * cdata(end,:));
                % Ensure that the display respects the alpha channel
                c.Face.Texture.ColorType = 'truecoloralpha';
                % Update the color data with the new transparency information
                c.Face.Texture.CData = cdata;
                c.Label.String = 'normal velocity, $v_n$ [$\mu$m/min]' ;
                c.Label.Interpreter = 'Latex' ;
                if qq == 1
                title(['normal velocity, $t = $', sprintf('%03d', tp - t0),...
                    ' ', QS.timeUnits], 'Interpreter', 'Latex')
                end
            end
            % Write figure to file
            saveas(gcf, outimfn_n)
            
        end

        % Save dilation field as image
        dilfn = fullfile(QS.dir.pivMultiChannel.dilation, [sprintf('%04d', tp) '.png']) ;
        if save_ims && (~exist(dilfn, 'file') || overwrite)
            close all
            
            % fig = figure('units', 'normalized', ...
            %         'outerposition', [0 0 1 1], 'visible', 'off') ;
            fig = figure('visible', 'off') ;
            for qq = 1:length(QS.xp.expMeta.channelsUsed)
                subplot(length(QS.xp.expMeta.channelsUsed), 1, qq, ...
                    'visible', 'off')
                
                labelOptsDil.title = 'dilation, $\log_{10}||J||$' ;
                alphaVal = 0.5 ;
                scale = 0.5 ;
                
                chIm = im0(:, :, rgbChannels(qq)) ;
                im = cat(3, chIm, chIm, chIm) ;
                scalarFieldOnImage(im, [x0(1, :)', y0(:, 1)], ...
                    reshape(log10(piv3d{qq}{ii}.dilation), ...
                        [length(x0(1,:)), length(y0(:,1))]),...
                    alphaVal, scale,...
                    labelOptsDil, 'style', 'diverging') ;
                ylim([size(im0, 2) * 0.25, size(im0, 2) * 0.75])
            end
            saveas(gcf, dilfn)
            close all
        end
        
        clear x0 y0 uu vv x1 y1
        clear pt0 v0 m0x m0y 
        clear pt1 v1 mesh1x mesh1y
        clear meshxy meshabove meshbelow 
        clear m0xy mf0 mf0c tr0 fieldfaces baryc0
        clear m1xy mf1 mf1c tr1 t1_contain baryc1
        clear vxa vya vza x123a y123a z123a
        clear vxb vyb vzb x123b y123b z123b
        clear datstruct
        
        % Save this timepoint's piv3d
        disp(['Saving to file: ' sprintf(piv3dfn, timePoints(ii))])
        piv3dstruct = {} ;
        for qq = 1:length(QS.xp.expMeta.channelsUsed)
            piv3dstruct = cat(1, piv3dstruct, piv3d{qq}{ii}) ;
        end
        save(sprintf(piv3dfn, timePoints(ii)), 'piv3dstruct')
        
    end 
end
disp('done with piv3dstruct')
