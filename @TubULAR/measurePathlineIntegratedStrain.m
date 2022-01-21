function measurePathlineIntegratedStrain(QS, options)
% integratePathlineStrain(QS, options)
%   Integrate the metric strain rate along Lagrangian pathlines.
%   Allow for median filtering along Lagrangian pathlines to avoid 
%   spurious spikes in accumulated strain.
%   Plot results in 2d and/or 3d for each timepoint.
%   
%
% Parameters
% ----------
% QS : QuapSlap class instance
% options : struct with fields
%   overwrite : bool, overwrite data on disk
%   overwriteImages : bool, overwrite images of results on disk
%   plot_comparison : bool, plot a comparison with DEC traceful dilation
%   median_filter_strainRates : bool
% 
% NPMitchell 2020

%% Default options 
overwrite = false ;
overwriteImages = false ;
plot_dzdp = false ;
plot_comparison = false ;
median_filter_strainRates = false ;
climitInitial = 0.05 ;
climitRamp = 0.005 ;
climitRatio = 1 ;

%% Parameter options
lambda_mesh = 0.002 ;
lambda = 0.01 ; 
debug = false ;
% Sampling resolution: whether to use a double-density mesh
samplingResolution = '1x'; 
averagingStyle = "Lagrangian" ;
% Load time offset for first fold, t0 -- default pathline t0
QS.t0set() ;
t0 = QS.t0 ;
% By default, t0Pathline = t0 (see below)

%% Unpack options & assign defaults
if nargin < 2
    options = struct() ;
end
if isfield(options, 'median_filter_strainRates')
    median_filter_strainRates = options.median_filter_strainRates ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'overwriteImages')
    overwriteImages = options.overwriteImages ;
end
if isfield(options, 'plot_comparison')
    plot_comparison = options.plot_comparison ;
end
if isfield(options, 'plot_dzdp')
    plot_dzdp = options.plot_dzdp ;
end
%% parameter options
if isfield(options, 'climitInitial')
    climitInitial = options.climitInitial ;
end
if isfield(options, 'climitRamp')
    climitRamp = options.climitRamp ;
end
if isfield(options, 'climitRatio')
    climitRatio = options.climitRatio ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'samplingResolution')
    samplingResolution = options.samplingResolution ;
end
if isfield(options, 'averagingStyle')
    averagingStyle = options.averagingStyle ;
end
if isfield(options, 'debug')
    debug = options.debug ;
end
if isfield(options, 't0Pathline')
    t0Pathline = options.t0Pathline ;
else
    t0Pathline = t0 ;
end

%% Determine sampling Resolution from input -- either nUxnV or (2*nU-1)x(2*nV-1)
if strcmp(samplingResolution, '1x') || strcmp(samplingResolution, 'single')
    doubleResolution = false ;
    sresStr = '' ;
elseif strcmp(samplingResolution, '2x') || strcmp(samplingResolution, 'double')
    doubleResolution = true ;
    sresStr = 'doubleRes_' ;
else 
    error("Could not parse samplingResolution: set to '1x' or '2x'")
end

%% Unpack QS
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;
if strcmp(averagingStyle, 'Lagrangian')
    sKDir = fullfile(QS.dir.strainRate.root, ...
        strrep(sprintf([sresStr 'lambda%0.3f_lmesh%0.3f'], ...
        lambda, lambda_mesh), '.', 'p'));
else
    error('Have not implemented strain rate measurements based on simple averaging')
end
folds = load(QS.fileName.fold) ;
fons = folds.fold_onset - QS.xp.fileMeta.timePoints(1) ;

%% Colormap
close all
set(gcf, 'visible', 'off')
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;
bbr256 = blueblackred(256) ;
close all

%% load from QS
if doubleResolution
    nU = QS.nU * 2 - 1 ;
    nV = QS.nV * 2 - 1 ;
else
    nU = QS.nU ;
    nV = QS.nV ;    
end

% We relate the normal velocities to the divergence / 2 * H.
tps = QS.xp.fileMeta.timePoints(1:end-1) - t0;

% Unit definitions for axis labels
unitstr = [ '[1/' QS.timeUnits ']' ];
vunitstr = [ '[' QS.spaceUnits '/' QS.timeUnits ']' ];
    
% DONE WITH PREPARATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load pathlines to build Kymographs along pathlines
QS.loadPullbackPathlines(t0Pathline, 'vertexPathlines')
vP = QS.pathlines.vertices ;

% Output directory is inside StrainRate dir
sKPDir = fullfile(sKDir, sprintf('pathline_%04dt0', t0Pathline)) ;
outdir = fullfile(sKPDir, 'measurements') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% Load Lx, Ly by loadingPIV. 
QS.loadPIV()
Xpiv = QS.piv.raw.x ;
Ypiv = QS.piv.raw.y ;

% Also need velocities to advect mesh
% QS.loadVelocityAverage('vv')
% vPIV = QS.velocityAverage.vv ;

% Discern if piv measurements are done on a double covering or the meshes
if strcmp(QS.piv.imCoords(end), 'e')
    doubleCovered = true ;
end

%% INTEGRATE STRAINRATE INTO STRAIN ON PATHLINES
% Compute or load all timepoints
strain = zeros(size(vP.vX, 2) * size(vP.vX, 3), 4) ; 
strain(:, 2) = 1e4 ;
strain(:, 3) = 1e4 ;

ntps = length(QS.xp.fileMeta.timePoints(1:end-1)) ;
for tidx = 1:ntps
    % Identify current timepoint
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    disp(['t = ' num2str(tp)])
    QS.setTime(tp) ;
    
    % Check for timepoint measurement on disk, on mesh vertices 
    estrainFn = fullfile(outdir, sprintf('strain_%06d.mat', tp)) ;
    
    if ~exist(estrainFn, 'file') || overwrite
        % Load timeseries measurements defined on mesh vertices 
        msg = 'Run QS.measurePathlineStrainRate() ' ;
        msg = [msg 'with lambdas=(field, mesh)=('] ;
        msg = [msg num2str(lambda) ','] ;
        msg = [msg num2str(lambda_mesh)] ;
        msg = [msg ') before running ', ...
                'QS.measurePathlineStrain(). '] ;
        msg = [msg '--> no file ' estrainFn] ;

        % Load current timepoint strainrate as 'strainrate'
        srfnMesh = fullfile(outdir, sprintf('strainRate_%06d.mat', tp)) ;
        try
            load(srfnMesh, 'strainrate', 'gg', 'dz', 'dp') 
        catch
            error(msg)
        end

        %% Identify prev and next timepoints's strainRates
        if median_filter_strainRates
            if tidx > 1
                tp_prev = QS.xp.fileMeta.timePoints(tidx-1) ;
                srfnMesh0 = fullfile(outdir, sprintf('strainRate_%06d.mat', tp_prev)) ;
                try
                    tmp = load(srfnMesh0, 'strainrate') ;
                    strainRate0 = tmp.strainrate ;
                    clearvars tmp
                catch 
                    error(msg)
                end
            end
            if tidx < ntps
                tp_next = QS.xp.fileMeta.timePoints(tidx + 1) ;
                srfnMesh2 = fullfile(outdir, sprintf('strainRate_%06d.mat', tp_next)) ;
                try
                    tmp = load(srfnMesh2, 'strainrate') ;
                    strainRate2 = tmp.strainrate ;
                    clearvars tmp
                catch
                    error(msg)
                end
            end
        end
        
        %% Recall current pathlines (vP) in pullback pixel space (XY)
        xx = vP.vX(tidx, :, :) ;
        yy = vP.vY(tidx, :, :) ;
        XY = [xx(:), yy(:)] ;
        Lx = vP.Lx(tidx) ;
        Ly = vP.Ly(tidx) ;
        options.Lx = Lx ;
        options.Ly = Ly ;
        XY = QS.doubleToSingleCover(XY, Ly) ;

        %% Recall strain rate at gridded vertices
        % strainrate is on pathlines
        % Metric gg is on pathlines
        % dz and dp scalings are on pathlines
        
        %% Save image of dz and dp
        dzdpDir = fullfile(sprintf(QS.dir.strainRate.pathline.root, ...
                lambda, lambda_mesh), 'dzdp') ;
        if plot_dzdp 
            if ~exist(dzdpDir, 'dir')
                mkdir(dzdpDir)
            end
            fn = fullfile(dzdpDir, sprintf('dzdp_%06d.png', tidx)) ;
            if ~exist(fn, 'file')
                % Check the ratio of lengths compared to pullback lengths
                subplot(2, 2, 1)
                trisurf(triangulation(mesh1.f, mesh1.v), dz, 'edgecolor', 'none')
                title('d$\zeta^{3D}/$d$\zeta^{2D}$', 'interpreter', 'latex')
                axis equal
                subplot(2, 2, 2)
                trisurf(triangulation(mesh1.f, mesh1.v), dp, 'edgecolor', 'none')
                title('d$\phi^{3D}/$d$\phi^{2D}$', 'interpreter', 'latex')
                axis equal
                subplot(2, 2, 3)
                meandz = mean(reshape(dz, [nU,nV]), 2) ;
                plot((1:nU)/nU, meandz, '.')
                if all(meandz < 600)
                    ylim([0, 800])
                end
                ylabel(['$\langle$d$\zeta\rangle_{dv}$ [' QS.spaceUnits ']'], ...
                    'Interpreter', 'Latex')
                xlabel('ap position, $\zeta$', 'interpreter', 'latex')
                % save this test image
                saveas(gcf, fn) ;
            end
        end

        %% Accumulate strain rate into STRAIN        
        % Transform the previous strain rate into current basis
        if tidx > 1
            % Load previous mesh
            tp_prev = QS.xp.fileMeta.timePoints(tidx-1) ;

            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
            mesh1 = tmp.spcutMeshSmRS ;
            clearvars tmp

            % Time increment
            dt = QS.timeInterval * (tp - tp_prev) ;
            assert(dt == 1)

            % DEBUG
            % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
            % mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) * mesh1.ar ;
            mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) ;
            umax = max(mesh1.u(:, 1)) ;  % Note umax = mesh1.ar
            vmax = max(mesh1.u(:, 2)) ;  % vmax is typically 1.0
            
            %--------------------------------------------------------------------------
            % Find where mesh1.u would be at previous timepoint by tracing
            % back along the flow. 
            % If we want to define a transformation between the frame of time point 1
            % and into the frame of time point 1 that ends up recapitulating the mesh
            % at time point 1, we must ask "Where were the vertices of the mesh at time
            % point 1 back during time point 0?".  We COULD extract this information
            % using the motion of the material points, but it would be more accurate to
            % use the PIV fields from which the material points were built.
            % The PIV consists of velocities vPIV and evaluation positions (Xpiv,Ypiv).
            %--------------------------------------------------------------------------
            % qq is the previous timepoint
            qq = tidx-1 ;

            % The 'inverse map' 
            % 1. Interpolate velocities of time t=0 at their advected 
            %    locations in t=1. 
            %
            %   .---->*     .<----*
            %     qq=.       qq+1=*
            % 
            %    Since the advected x0+u_qq,y0+v_qq are spatially unstructured, we
            %    interpolate t=0 velocity at the advected positions,
            %    ie what velocities took XY0 to XY1, evaluate the velocities 
            %    that pushed to the next positions AT the next positions XY1,
            %    and pull back XY1 along those velocities to get displaced coordinates 
            %    DXY0 which will land on XY1 when moved along t0's flow.
            %
            uu = QS.piv.raw.u_filtered{qq} ;
            vv = QS.piv.raw.v_filtered{qq} ;
            % Force boundaries to have zero flux at AP cuts. This removes 
            % some spurious movement of the anterior and posterior cuts
            uu(:, 1) = 0 ;
            uu(:, end) = 0 ;
            
            % Gently smooth mesh vertex displacements -- pad with zeros
            uu = imgaussfilt(uu, 0.5, 'Padding', 0) ;
            vv = imgaussfilt(vv, 0.5, 'Padding', 0) ;
            
            % Load Lx, Ly for t=1, which are the extents of the domain (needed for
            % periodic boundary consitions and to keep all the points in the domain of
            % interpolation). 
            x01 = Xpiv{tidx-1}(:) + uu(:) ;
            y01 = Ypiv{tidx-1}(:) + vv(:) ;
            [xa, ya] = QS.clipXY(x01, y01, Lx, Ly) ;
            xya_uv = QS.XY2uv([Lx, Ly], [xa, ya], ...
                doubleCovered, umax, vmax) ;

            ui = scatteredInterpolant(xya_uv(:, 1), xya_uv(:, 2), ...
                uu(:)*(umax / Lx), 'natural', 'nearest') ;
            vi = scatteredInterpolant(xya_uv(:, 1), xya_uv(:, 2), ...
                vv(:)*(vmax / Ly), 'natural', 'nearest') ;

            % 2. Evaluate at mesh1 vertices to pull them back along velocity to t=0.
            dx = ui(mesh1.u(:, 1), mesh1.u(:, 2)) ;
            dy = vi(mesh1.u(:, 1), mesh1.u(:, 2)) ;

            % 3. Pull mesh vertices back to t=0
            Xqq = mesh1.u(:, 1) - dx(:) ;
            Yqq = mesh1.u(:, 2) - dy(:) ;

            % The locations of the rectilinear vertices at time 0
            % 'Deformed Vertices from t1 at t0'
            DXY10 = [Xqq, Yqq ] ;

            % % Check the deformed mesh that will advect into the gridded
            % % mesh at t+1
            % clf
            % % scatter(x01(:), y01(:), 10, uu(:), 'filled');
            % hold on;
            % quiver(Xqq, Yqq, dx(:), dy(:), 1, 'c')
            % triplot(triangulation(mesh1.f, DXY10));
            % axis equal            

            % Construct the Jacobian matrix on each mesh face
            J01 = jacobian2Dto2DMesh(mesh1.u, DXY10, mesh1.f);
            
            % % Check that jacobians are nearly unity
            % avgJ = [0 0 0 0]' ;
            % for tj = 1:length(J01)
            %     avgJ = avgJ + J01{tj}(:) ;
            % end
            % avgJ = avgJ ./ length(J01) ;
            
            % If any mesh faces are flipped, force their Jacobian 
            % transformation to be unity
            normals = faceNormal(triangulation(mesh1.f, ...
                [DXY10, zeros(size(DXY10,1), 1)]) );
            if (numel(unique(normals(:,3))) ~= 1)
                warning('Triangles are not consistently ordered: forcing bad Jacobians to be identity')
                badTris = find(normals(:,3) ~= mode(sign(normals(:,3)))) ;
                for badq = badTris
                    J01{badq} = [1, 0; 0, 1] ;
                end
            end
            
            %% ------------------------------------------------------------
            % Find which jacobians to use for each pathline point
            tri = triangulation(mesh1.f, mesh1.u) ;
            if strcmp(QS.piv.imCoords, 'sp_sme')
                im = imfinfo(sprintf(QS.fullFileBase.im_sp_sme, tp)) ;
                im = [im.Height, im.Width] ;
                if im(1) ~= im(2)
                    error('check that dimensions are in correct order here')
                end
            else
                error('Handle this case for imCoords here')
            end
            uv = QS.XY2uv(im, XY, doubleCovered, umax, vmax) ;
            uv(:, 2) = mod(uv(:, 2), vmax) ; 
            uv(:, 1) = max(uv(:, 1), 1e-14) ; 
            uv(:, 1) = min(uv(:, 1), umax-1e-14) ; 
            fieldfaces = pointLocation(tri, uv) ;

            % Check mapping to pixel space
            assert(all(min(fieldfaces) > 0))
            try
                assert(~any(isnan(fieldfaces)))
            catch
                error('Ensure that all uv lie in the mesh.u')
            end
            
            % % Check fieldfaces
            % trisurf(mesh1.f(fieldfaces, :), mesh1.u(:, 1), mesh1.u(:, 2), ...
            %     mesh1.u(:, 2) * 0)
            % view(2)
            
            %% Median filter the strainRates while accumulating them 
            if median_filter_strainRates
                % Transform the strainRate from t0 to t1 using 
                %   sR0 = (inv(J01f) * strainRate0 * inv(J01f).') 
                % Compare this quantity to
                %   sR1 = strainRate1
                % and also to 
                %   sR2 = (inv(J21f) * strainRate2 * inv(J21f).').
                % Take the median of each diagonal and the off-diagonal. 
                
                %% Build sR0 (prev timePoint)
                sR0 = zeros(size(fieldfaces, 1), 3) ;
                for qq = 1:size(fieldfaces,1)
                    sR0qf = [strainRate0(qq, 1), strainRate0(qq, 2); 
                        strainRate0(qq, 3), strainRate0(qq, 4)] ;  
                    J01f = J01{fieldfaces(qq)} ;
                    sR0q = inv(J01f) * sR0qf * (inv(J01f).') ;
                    % arrange as [zz, 0.5(zp+pz), pp]
                    sR0(qq, :) = [sR0q(1, 1), ...
                        0.5 * (sR0q(1, 2) + sR0q(2, 1)), ...
                        sR0q(2, 2) ] ;
                end
                
                %% Build sR1 (current timePoint)
                assert(all(strainrate(:, 2) == strainrate(:, 3)))
                sR1 = [strainrate(:, 1),  strainrate(:, 2), strainrate(:, 4) ];
                
                %% Build sR2 (next timePoint)
                % Create J21f, the Jacobian from the next timepoint's mesh
                % to the current mesh
                % 1. Recall piv "velocities" (optical flow in pullback
                % pixel space)
                u1 = QS.piv.raw.u_filtered{tidx} ;
                v1 = QS.piv.raw.v_filtered{tidx} ;
                % Gently smooth mesh vertex displacements -- pad with zeros
                u1 = imgaussfilt(u1, 0.5, 'Padding', 0) ;
                v1 = imgaussfilt(v1, 0.5, 'Padding', 0) ;

                % 2. Interpolate optical flow onto mesh vertices
                uvpiv = QS.XY2uv([Lx, Ly], [Xpiv{tidx}(:), Ypiv{tidx}(:)], ...
                    doubleCovered, umax, vmax) ;
                Upiv = reshape(uvpiv(:, 1), size(Xpiv{tidx})) ;
                Vpiv = reshape(uvpiv(:, 2), size(Ypiv{tidx})) ;
                ui = griddedInterpolant(Upiv', Vpiv', ...
                    u1' * (umax / Lx), 'linear', 'nearest') ;
                vi = griddedInterpolant(Upiv', Vpiv', ...
                    v1' * (vmax / Ly), 'linear', 'nearest') ;
                u1 = ui(mesh1.u(:, 1), mesh1.u(:, 2)) ;
                v1 = vi(mesh1.u(:, 1), mesh1.u(:, 2)) ;
                
                % 3. Push mesh vertices forward to next timepoint t=2
                Xqq = mesh1.u(:, 1) + u1(:) ;
                Yqq = mesh1.u(:, 2) + v1(:) ;

                % The locations of the rectilinear vertices at time 2
                % 'Deformed Vertices from t1 at t2'
                DXY12 = [Xqq, Yqq ] ;
                
                % 4. Construct the Jacobian matrix on each mesh face
                J21 = jacobian2Dto2DMesh(mesh1.u, DXY12, mesh1.f);
                
                % If any mesh faces are flipped, force their Jacobian 
                % transformation to be unity
                normals = faceNormal(triangulation(mesh1.f, ...
                    [DXY12, zeros(size(DXY12,1), 1)]) );
                if (numel(unique(normals(:,3))) ~= 1)
                    warning('Triangles are not consistently ordered: forcing bad Jacobians to be identity')
                    badTris = find(normals(:,3) ~= mode(sign(normals(:,3)))) ;
                    for badq = badTris
                        J21{badq} = [1, 0; 0, 1] ;
                    end
                end
                
                % 5. Transport the strainRate from the next timepoint to 
                % the current timepoint
                sR2 = zeros(size(fieldfaces, 1), 3) ;
                for qq = 1:size(fieldfaces,1)
                    sR2qf = [strainRate0(qq, 1), strainRate0(qq, 2); 
                        strainRate0(qq, 3), strainRate0(qq, 4)] ;  
                    
                    J21f = J21{fieldfaces(qq)} ;
                    sR2q = inv(J21f) * sR2qf * (inv(J21f).') ;
                    % arrange as [zz, 0.5(zp+pz), pp]
                    sR2(qq, :) = [sR2q(1, 1), ...
                        0.5 * (sR2q(1, 2) + sR2q(2, 1)), ...
                        sR2q(2, 2) ] ;
                end
                
                %% Take median of the three
                strainrate = nanmedian(cat(3, sR0, sR1, sR2), 3) ;
                
                % check it
                if debug
                    titlenames = {'$\varepsilon_{\zeta\zeta}$ filtering', ...
                        '$\varepsilon_{\zeta\phi}$ filtering', ...
                        '$\varepsilon_{\phi\phi}$ filtering'} ;
                    for dim = 1:3
                        subplot(2, 2, 1)
                        imagesc(reshape(strainrate(:, dim), [nU, nV]))
                        title('median filtered')
                        colorbar
                        subplot(2, 2, 2)
                        imagesc(reshape(sR1(:, dim), [nU, nV]))
                        title('original')
                        colorbar
                        subplot(2, 2, 3)
                        imagesc(reshape(abs(strainrate(:, dim)-sR1(:, dim)),...
                            [nU, nV]))
                        title('difference')
                        colorbar
                        subplot(2, 2, 4)
                        histogram(strainrate(:, dim)-sR1(:, dim))
                        title('difference')
                        sgtitle(titlenames{dim}, 'interpreter', 'latex')
                        waitfor(gcf)
                    end
                end
                % Store strainrate in #vP x 4 array
                ezz = strainrate(:, 1) ;
                ezp = strainrate(:, 2) ;
                epz = strainrate(:, 2) ;
                epp = strainrate(:, 3) ;
            else
                % Store strainrate in #vP x 4 array
                ezz = strainrate(:, 1) ;
                ezp = strainrate(:, 2) ;
                epz = strainrate(:, 3) ;
                epp = strainrate(:, 4) ;
            end
            %% ------------------------------------------------------------
            % Consider each pathline, add the strainRate * dt to the
            % accumulated strain from previous timepoint, which is called 
            % strain0. 
            % Accumulate the strain via implicit Euler (backward Euler
            % scheme)
            % Transform as a (0,2)-tensor (NOTICE THE MATRIX INVERSES)
            strainM = cell(size(fieldfaces)) ;
            for qq = 1:size(fieldfaces,1)
                strain0 = [strain(qq, 1), strain(qq, 2); ...
                           strain(qq, 3), strain(qq, 4)] ;
                strainrateP = [ezz(qq), ezp(qq); epz(qq), epp(qq)] ;
                try
                    J01f = J01{fieldfaces(qq)} ;
                    strainM{qq} = inv(J01f) * strain0 * (inv(J01f).') + ...
                        dt * strainrateP ;
                catch
                    error('Ensure that all uv lie in the mesh.u')
                end
            end

            % Convert strain from cell (needed for face Jacobians) to array
            % NOTE: This updates the variable strain from previous
            % accumulated strain to current (forward-looking) accumulated
            % strain.
            for qq = 1:size(ezz, 1)
                strain(qq, 1) = strainM{qq}(1, 1) ;
                strain(qq, 2) = strainM{qq}(1, 2) ;
                strain(qq, 3) = strainM{qq}(2, 1) ;
                strain(qq, 4) = strainM{qq}(2, 2) ; 
            end
        else
            % Load mesh for this timepoint
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
            mesh1 = tmp.spcutMeshSmRS ;
            clearvars tmp

            % DEBUG
            % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
            % mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) * mesh1.ar ;
            mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1))  ;

            dt = QS.timeInterval ;
            % Note: strainrate is already interpolated onto pathlines
            strain = dt * strainrate ;
        end

        %% Trace/Determinant of strain 
        strain_tr = zeros(size(strain, 1), 1) ;
        strain_dv = zeros(size(strain, 1), 1) ;
        strain_th = zeros(size(strain, 1), 1) ;
        for qq = 1:size(strain, 1)
            eq = [strain(qq, 1), strain(qq, 2); ...
                  strain(qq, 3), strain(qq, 4)] ;
            gq = [gg(qq, 1), gg(qq, 2); ...
                  gg(qq, 3), gg(qq, 4)] ;
            [strain_tr(qq), strain_dv(qq), strain_th(qq)] = ...
                traceDeviatorPullback(eq, gq, dz(qq), dp(qq));
        end

        %% CHECK integration 
        if debug
            % Map intensity from dev and color from the theta
            close all
            pm256 = phasemap(256) ;
            indx = max(1, round(mod(2*strain_th(:), 2*pi)*size(pm256, 1)/(2 * pi))) ;
            colors = pm256(indx, :) ;
            devKclipped = min(strain_dv / max(strain_dv), 1) ;
            colorsM = devKclipped(:) .* colors ;
            colorsM = reshape(colorsM, [nU, nV, 3]) ;
            imagesc(1:nU, 1:nV, permute(colorsM, [2, 1, 3]))
            caxis([0, 0.1])
            pause(1)
        end

        %% OPTION 1: simply reshape, tracing each XY pathline pt to its t0
        % % grid coordinate
        strain_tr = reshape(strain_tr, [nU, nV]) ;
        strain_dv = reshape(strain_dv, [nU, nV]) ;
        strain_th = reshape(strain_th, [nU, nV]) ;

        %% OPTION 2: the following regrids onto original XY coordinates,
        % rendering the process of following pathlines moot. 
        % Average into AP bins and take mean along 1/4 DV hoop arcs
        % if doubleCovered
        %     vminmax = [0.25 * Ly, 0.75 * Ly] ;
        % else
        %     vminmax = [1, Ly] ;
        % end
        %
        % Note the transposition: to plot as APDV, imshow(m')
        % HH = binData2dGrid([XY, HH], [1,Lx], vminmax, nU, nV) ;
        % gdot = binData2dGrid([XY, gdot], [1,Lx], vminmax, nU, nV) ;
        % divv = binData2dGrid([XY, divv], [1,Lx], vminmax, nU, nV) ;
        % veln = binData2dGrid([XY, veln], [1,Lx], vminmax, nU, nV) ;
        % H2vn = binData2dGrid([XY, H2vn], [1,Lx], vminmax, nU, nV) ;

        %% Average STRAIN (accumulated strain) along DV 
        % Average along DV -- ignore last redudant row at nV
        [strain_dv_ap, strain_th_ap] = ...
            QS.dvAverageNematic(strain_dv(:, 1:nV-1), strain_th(:, 1:nV-1)) ;
        strain_tr_ap = mean(strain_tr(:, 1:nV-1), 2) ;

        % quarter bounds
        q0 = round(nV * 0.125) ;
        q1 = round(nV * 0.375) ;
        q2 = round(nV * 0.625) ;
        q3 = round(nV * 0.875) ;
        left = q0:q1 ;
        ventral = q1:q2 ;
        right = q2:q3 ;
        dorsal = [q3:nV, 1:q1] ;

        % left quarter
        [strain_dv_l, strain_th_l] = ...
            QS.dvAverageNematic(strain_dv(:, left), strain_th(:, left)) ;
        strain_tr_l = mean(strain_tr(:, left), 2) ;
        % right quarter
        [strain_dv_r, strain_th_r] = ...
            QS.dvAverageNematic(strain_dv(:, right), strain_th(:, right)) ;
        strain_tr_r = mean(strain_tr(:, right), 2) ;
        % dorsal quarter
        [strain_dv_d, strain_th_d] = ...
            QS.dvAverageNematic(strain_dv(:, dorsal), strain_th(:, dorsal)) ;
        strain_tr_d = mean(strain_tr(:, dorsal), 2) ;
        % ventral quarter
        [strain_dv_v, strain_th_v] = ...
            QS.dvAverageNematic(strain_dv(:, ventral), strain_th(:, ventral)) ;
        strain_tr_v = mean(strain_tr(:, ventral), 2) ;

        % % Check the strain accumulation
        % set(gcf, 'visible', 'on')
        % subplot(2, 1, 1)
        % imagesc(strain_tr')
        % title(['t = ', num2str(tp)])
        % limit = std(strain_tr(:)) ;
        % caxis([-limit, limit])
        % colormap(bbr256)
        % colorbar()
        % subplot(2, 1, 2)
        % imagesc(strain_dv')
        % colormap parula
        % caxis([0, std(abs(strain_dv(:)))])
        % colorbar()
        % drawnow
        
        %% Save the strain 
        readme.strain = 'integrated strain tensor' ;
        readme.strain_tr = 'integrated strain trace Tr[g^{-1} epsilon]';        
        readme.strain_dv = 'integrated strain deviator magnitude sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) on mesh vertices';
        readme.strain_th = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector on mesh vertices';
        readme.strain_dv_ap = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged circumferentially';
        readme.strain_dv_l = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on left quarter, on vertices';
        readme.strain_dv_r = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on right quarter, on vertices';
        readme.strain_dv_d = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on dorsal quarter, on vertices';
        readme.strain_dv_v = 'integrated sqrt( Tr[g^{-1} epsilon g^{-1} epsilon] ) averaged on ventral quarter, on vertices';
        readme.strain_th_ap = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged circumferentially, on vertices';
        readme.strain_th_l = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on left quarter, on vertices';
        readme.strain_th_r = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on right quarter, on vertices';
        readme.strain_th_d = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on dorsal quarter, on vertices';
        readme.strain_th_v = 'integrated strain deviator angle -- arctan( e_phi / e_zeta ), where e is the positive-eigenvalue eigenvector, averaged on ventral quarter, on vertices';
        readme.strain_tr_ap = 'integrated strain trace Tr[g^{-1} epsilon], averaged circumferentially, on vertices';
        readme.strain_tr_l = 'integrated strain trace Tr[g^{-1} epsilon], averaged on left quarter, on vertices';
        readme.strain_tr_r = 'integrated strain trace Tr[g^{-1} epsilon], averaged on right quarter, on vertices';
        readme.strain_tr_d = 'integrated strain trace Tr[g^{-1} epsilon], averaged on dorsal quarter, on vertices';
        readme.strain_tr_v = 'integrated strain trace Tr[g^{-1} epsilon], averaged on ventral quarter, on vertices';

        readme.note = ['Integrated strain evaluated on Lagrangian paths. ', ...
            'The pullback space is taken to range from ', ...
            'zeta=[0, a_relaxed] and phi=[0, 1], where a_relaxed minimizes ', ...
            'Dirichlet energy of the map'] ; 
        disp(['saving ', estrainFn])
        save(estrainFn, 'gg', 'strain', 'dz', 'dp', 'readme', ...
            'strain_tr', 'strain_dv', 'strain_th', ...
            'strain_dv_ap', 'strain_dv_l', 'strain_dv_r', ...
            'strain_dv_d', 'strain_dv_v', ...
            'strain_tr_ap', 'strain_tr_l', 'strain_tr_r', ...
            'strain_tr_d', 'strain_tr_v', ...
            'strain_th_ap', 'strain_th_l', 'strain_th_r', ...
            'strain_th_d', 'strain_th_v')
    else
        load(estrainFn, 'strain', 'strain_tr_ap', 'strain_dv_ap')
        % Load mesh 
        tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp)) ;
        mesh1 = tmp.spcutMeshSmRS ;
        clearvars tmp
        % Normalize the zeta to fixed aspect ratio (ar=aspectratio relaxed)
        % DEBUG
        % mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) * mesh1.ar ;
        mesh1.u(:, 1) = mesh1.u(:, 1) / max(mesh1.u(:, 1)) ;
    end
    
    % Plot the result
    plotOpts.overwrite = overwriteImages ;
    plotOpts.cutMesh = mesh1 ;
    plotOpts.lambda = lambda ;
    plotOpts.lambda_mesh = lambda_mesh ;
    plotOpts.debug = debug ;
    plotOpts.t0Pathline = t0Pathline ;
    plotOpts.plot_comparison = plot_comparison ;
    bulk = round([0.2 * nU, 0.7 * nU]) ;
    plotOpts.clim_trace = 2 * max(abs(strain_tr_ap(bulk))) ; 
    % climitInitial + climitRamp * (tp - QS.xp.fileMeta.timePoints(1)) * QS.timeInterval ;
    plotOpts.clim_deviator = 5 * max(abs(strain_dv_ap(bulk))) ;
    QS.plotPathlineStrainTimePoint(tp, plotOpts)
end
disp('done with integrated pathline strain calculations')


%% Combine DV-averaged profiles into kymographs
apKymoFn = fullfile(outdir, 'apKymographPathlineStrain.mat') ;
lKymoFn = fullfile(outdir, 'leftKymographPathlineStrain.mat') ;
rKymoFn = fullfile(outdir, 'rightKymographPathlineStrain.mat') ;
dKymoFn = fullfile(outdir, 'dorsalKymographPathlineStrain.mat') ;
vKymoFn = fullfile(outdir, 'ventralKymographPathlineStrain.mat') ;
files_exist = exist(apKymoFn, 'file') && ...
    exist(lKymoFn, 'file') && exist(rKymoFn, 'file') && ...
    exist(dKymoFn, 'file') && exist(vKymoFn, 'file') ;
if ~files_exist || overwrite || true
    disp('Compiling kymograph data to save to disk...')
    for tp = QS.xp.fileMeta.timePoints(1:end-1)
        tidx = QS.xp.tIdx(tp) ;

        % Check for timepoint measurement on disk
        srfn = fullfile(outdir, sprintf('strain_%06d.mat', tp))   ;

        % Load timeseries measurements
        load(srfn, 'strain_tr_ap', 'strain_tr_l', 'strain_tr_r', ...
            'strain_tr_d', 'strain_tr_v', ...
            'strain_th_ap', 'strain_th_l', 'strain_th_r', ....
            'strain_th_d', 'strain_th_v', ...
            'strain_dv_ap', 'strain_dv_l', 'strain_dv_r', ...
            'strain_dv_d', 'strain_dv_v') ;

        %% Store accumulated strain in matrices
        % dv averaged
        str_apM(tidx, :) = strain_tr_ap ;
        sdv_apM(tidx, :) = strain_dv_ap ;
        sth_apM(tidx, :) = strain_th_ap ;

        % left quarter
        str_lM(tidx, :) = strain_tr_l ;
        sdv_lM(tidx, :) = strain_dv_l ;
        sth_lM(tidx, :) = strain_th_l ;

        % right quarter
        str_rM(tidx, :) = strain_tr_r ;
        sdv_rM(tidx, :) = strain_dv_r ;
        sth_rM(tidx, :) = strain_th_r ;

        % dorsal quarter
        str_dM(tidx, :) = strain_tr_d ;
        sdv_dM(tidx, :) = strain_dv_d ;
        sth_dM(tidx, :) = strain_th_d ;

        % ventral quarter
        str_vM(tidx, :) = strain_tr_v ;
        sdv_vM(tidx, :) = strain_dv_v ;
        sth_vM(tidx, :) = strain_th_v ;
    end
    
    %% Save kymographs
    disp('Saving kymograph data files for Lagrangian pathlines')
    save(apKymoFn, 'str_apM', 'sdv_apM', 'sth_apM')
    disp(['Saved kymograph data to: ' apKymoFn])
    save(lKymoFn, 'str_lM', 'sdv_lM', 'sth_lM')
    disp(['Saved kymograph data to: ' lKymoFn])
    save(rKymoFn, 'str_rM', 'sdv_rM', 'sth_rM')
    disp(['Saved kymograph data to: ' rKymoFn])
    save(dKymoFn, 'str_dM', 'sdv_dM', 'sth_dM')
    disp(['Saved kymograph data to: ' dKymoFn])
    save(vKymoFn, 'str_vM', 'sdv_vM', 'sth_vM')
    disp(['Saved kymograph data to: ' vKymoFn])
    disp('done with strain kymograph data saving')
else
    disp('strain kymograph data already on disk')    
end

