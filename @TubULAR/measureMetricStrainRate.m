function measureMetricStrainRate(QS, options)
%measureMetricStrainRate(QS, options)
%   
% Parameters
% ----------
%
% Returns 
% -------
%
% NPMitchell 2020

%% Default options
lambda = 0.01 ;
lambda_mesh = 0.002 ;
overwrite = false ;
preview = true ;
metric_style = 'strain' ;  % {'mesh', 'strain'} ;
clim_trgdot = 0.2 ;
clim_tg = 1 ;

%% Unpack options
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'preview')
    preview = options.preview ;
end
if isfield(options, 'lambda')
    lambda = options.lambda ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end

%% Unpack QS
if strcmp(metric_style, 'strain')
    egDir = QS.dir.gstrainRate ;
    egImDir = QS.dir.gstrainRateIm ;
elseif strcmp(metric_style, 'mesh')
    egDir = QS.dir.gstrainMesh ;
    egImDir = QS.dir.gstrainMeshIm ;
end
QS.getXYZLims ;
xyzlim = QS.plotting.xyzlim_um ;
buff = 10 ;
xyzlim = xyzlim + buff * [-1, 1; -1, 1; -1, 1] ;

%% Prepare for plots
colors = define_colors ;
blue = colors(1, :) ;
red = colors(2, :) ;
green = colors(3, :) ;

%% Colormap
close all
imagesc([-1, 0, 1; -1, 0, 1])
caxis([-1, 1])
bwr256 = bluewhitered(256) ;

%% Prepare both metric styles
dirs2make = { ...
    fullfile(egImDir, metric_style), ...
    fullfile(egImDir, ['trgdot_' metric_style]), ...
    fullfile(egImDir, ['gss_' metric_style]), ...
    fullfile(egImDir, ['gsphi_' metric_style]), ...
    fullfile(egImDir, ['gphiphi_' metric_style]), ...
    fullfile(egImDir, ['metricstrain_' metric_style '_2d'])} ;
for ii = 1:length(dirs2make)
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make)
    end
end

%% ------------------------------------------------------------------------
% Construct 2D mesh corresponding to the planar domain of
% parameterization
%--------------------------------------------------------------------------
% % Load from file
% path = fullfile(NESpath, 'NES_Examples') ;
% mesh = read_ply_mod(fullfile(path, 'tube_simple_h1p00_R0p00_w1p00.ply')) ;
% rmID = [length(mesh.v)-1, length(mesh.v)] ;
% [F, V] = remove_vertex_from_mesh(mesh.f, mesh.v, rmID) ;
% % Center the mesh around the x axis
% midx = 0.5 * (max(mesh.v(:, 1)) + min(mesh.v(:, 1))) ;
% V(:, 1) = V(:, 1) - midx ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Consider each kind of metric strain measurement
nU = QS.nU ;
nV = QS.nV ;
QS.t0set() ;
tfold = QS.t0 ;

% Load vertex-based velocity measurements
vvsmMfn = QS.fileName.pivAvg.vv  ;
tmp = load(vvsmMfn) ;
vertex_vels = tmp.vvsmM ;
vfsmMfn = QS.fileName.pivAvg.vf ;
tmp = load(vfsmMfn) ;
face_vels = tmp.vfsmM ;

% Pre-assign timepoints with velocities
tpts = QS.xp.fileMeta.timePoints(1:end-1) ;
tp2do = [60, tpts(1:10:end), setdiff(tpts, tpts(1:10:end))] ;

% Build metric from mesh
for tp = tp2do
    disp(['t = ' num2str(tp)])
    tidx = QS.xp.tIdx(tp) ;

    % Load current mesh
    tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp)) ;
    mesh = tmp.spcutMeshSmRSC ;
    clearvars tmp

    % Define metric strain filename        
    if strcmp(metric_style, 'mesh')
        gstrainFn = sprintf(QS.fullFileBase.gstrainMesh, tp) ;
    elseif strcmp(metric_style, 'strain')
        gstrainFn = sprintf(QS.fullFileBase.gstrainRate, tp) ;
    end

    % Compute the metric strain if not on disk
    if ~exist(gstrainFn, 'file') || overwrite
        if exist(gstrainFn, 'file')
            disp('Overwriting metric strain on disk')
        else
            disp('Computing metric strain anew')
        end
        
        % Advect the mesh along instantaneous velocity or load next tp mesh
        % (depending on which style to compute)
        if strcmp(metric_style, 'mesh')
            % Load next mesh also
            tmp = load(sprintf(QS.fullFileBase.spcutMeshSmRSC, tp + 1)) ;
            mesh2 = tmp.spcutMeshSmRSC ;
            clearvars tmp
        elseif strcmp(metric_style, 'strain')

            % Smooth the mesh vertices
            if lambda_mesh > 0 
                tri = triangulation(mesh.f, mesh.v) ;
                fbndy = tri.freeBoundary ;
                fbndy = fbndy(:, 1) ;
                mesh.v = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                    lambda_mesh, 'implicit', mesh.v) ;
                % check smoothed mesh
                % trisurf(triangulation(mesh.f, mesh.v), 'edgecolor', 'none')
            end

            % Advect the mesh to mesh2
            mesh2 = mesh ;

            % Smooth the velocities in space using gptoolbox
            vraw = squeeze(vertex_vels(tidx, 1:(nV-1)*nU, :)) ;
            vs = laplacian_smooth(mesh.v, mesh.f, 'cotan', fbndy, ...
                lambda, 'implicit', vraw) ;

            mesh2.v = mesh.v + vs ;

            % Check the vels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if preview
                meshes = [mesh, mesh2] ;
                for meshId = 1:2
                    mm = meshes(meshId) ;
                    clf
                    vlim = 0.8 ;
                    vxM = reshape(vraw(:, 1), [nU, nV-1]) ;
                    vyM = reshape(vraw(:, 2), [nU, nV-1]) ;
                    vzM = reshape(vraw(:, 3), [nU, nV-1]) ;
                    vxsM = reshape(vs(:, 1), [nU, nV-1]) ;
                    vysM = reshape(vs(:, 2), [nU, nV-1]) ;
                    vzsM = reshape(vs(:, 3), [nU, nV-1]) ;
                    subplot(2, 3, 1)
                    colormap bwr; 
                    trisurf(mm.f, mm.v(:, 1), mm.v(:, 2), mm.v(:, 3),...
                        vxM, 'edgecolor', 'none')
                    % colorbar(); 
                    caxis([-vlim, vlim]); 
                    axis equal; axis off ; 
                    % xlabel('x'); ylabel('y'); zlabel('z')
                    title('AP velocity snapshot -- pre-smoothing')
                    subplot(2, 3, 4)
                    trisurf(mm.f, mm.v(:, 1), mm.v(:, 2), mm.v(:, 3),...
                        vxsM, 'edgecolor', 'none')
                    % colorbar(); 
                    caxis([-vlim, vlim]); 
                    axis equal; axis off ; 
                    % xlabel('x'); ylabel('y'); zlabel('z')
                    title('AP velocity snapshot')
                    pause(1)

                    % Check other velocities
                    subplot(2, 3, 2)
                    colormap bwr; 
                    trisurf(mm.f, mm.v(:, 1), mm.v(:, 2), mm.v(:, 3),...
                        vyM, 'edgecolor', 'none')
                    % colorbar();
                    caxis([-vlim, vlim]); 
                    axis equal; axis off ; 
                    % xlabel('x'); ylabel('y'); zlabel('z')
                    title('lateral velocity snapshot -- pre-smoothing')
                    subplot(2, 3, 5)
                    trisurf(mm.f, mm.v(:, 1), mm.v(:, 2), mm.v(:, 3),...
                        vysM, 'edgecolor', 'none')
                    % colorbar(); 
                    caxis([-vlim, vlim]); 
                    axis equal; axis off ; 
                    % xlabel('x'); ylabel('y'); zlabel('z')
                    title('lateral velocity snapshot')
                    subplot(2, 3, 3)
                    colormap bwr; 
                    trisurf(mm.f, mm.v(:, 1), mm.v(:, 2), mm.v(:, 3),...
                        vzM, 'edgecolor', 'none')
                    % colorbar(); 
                    caxis([-vlim, vlim]); 
                    axis equal; axis off ; 
                    % xlabel('x'); ylabel('y'); zlabel('z')
                    title('DV velocity snapshot -- pre-smoothing')
                    subplot(2, 3, 6)
                    trisurf(mm.f, mm.v(:, 1), mm.v(:, 2), mm.v(:, 3),...
                        vzsM, 'edgecolor', 'none')
                    % colorbar(); 
                    caxis([-vlim, vlim]); 
                    axis equal; axis off ; 
                    % xlabel('x'); ylabel('y'); zlabel('z')
                    title('DV velocity snapshot')
                    set(gcf, 'visible', 'on')
                    pause(1)
                    clf
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            error('Bad strain rate style option')
        end
        
        %% Construct Topolgical Structure Tools ===============================
        % eg = metricStrainSPhiGridMesh(mesh, mesh2) ;

        % Cut the mesh into a cutMesh to open up the pullback mesh & also
        % grab du and dv for each face --  bond vecs along u and along v
        mesh.nU = nU ;
        tmp_options.preview = false ;
        [~, dbonds] = labelRectilinearMeshBonds(mesh, tmp_options) ;
        cutMesh = dbonds.cutMesh ; 
        
        % Cut the advected / next mesh
        mesh2.nU = nU ;
        cutMesh2 = cutRectilinearCylMesh(mesh2) ;
        
        % Metric for mesh
        g0cell = inducedMetric(cutMesh.f, cutMesh.v, cutMesh.u) ;
        % Convert metric tensor to #Faces x 2 x 2
        g0 = zeros(size(g0cell, 1), 2, 2) ;
        for qq = 1:size(g0cell, 1) 
            g0(qq, :, :) = g0cell{qq} ;
        end

        % Metric for mesh2
        g1cell = inducedMetric(cutMesh2.f, cutMesh2.v, cutMesh2.u) ;
        % Convert metric tensor to #Faces x 3
        g1 = zeros(size(g1cell, 1), 2, 2) ;
        for qq = 1:size(g1cell, 1)
            g1(qq, :, :) = g1cell{qq} ;
        end

        % Metric strain
        tre = zeros(size(g0, 1), 1) ;
        if strcmp(metric_style, 'strain')
            eg = zeros(size(g0)) ;
            for qq = 1:size(g0cell, 1)
                g0q = g0cell{qq} ;
                g1q = g1cell{qq} ;
                eg(qq, :, :) = 0.5 * (g1q - g0q) ;
                % tre(qq) = trace(0.5 * inv(g0q) * (g1q - g0q)) ;
                % tre(qq) = trace(0.5 * inv(g0q) * (g1q - g0q)) ;
            end
            a0 = 0.5 * doublearea(mesh.v, mesh.f) ;
            a1 = 0.5 * doublearea(mesh2.v, mesh2.f) ;
            dilation = (a1 - a0) ./ a0 ;
        elseif strcmp(metric_style, 'mesh')
            eg = g1 - g0 ;
        end

        % Build deformation tensor over mesh, whose elements are the change in
        % squared length in each direction (u.u,u.v,v.v) 
        tg = 0 * eg ;
        du = dbonds.baseSpace.u ;
        dv = dbonds.baseSpace.v ;
        for qq = 1:size(g1cell, 1)
            gq = squeeze(eg(qq, :, :)) ;

            % NORMALIZATION 
            % denom = [1 1 1 1] ;
            % % Frobenius norm of the metric
            % ggqq = squeeze(g0(qq, :, :)) ;
            % % denom = trace(gq .* ggqq') ;
            % 
            % denom(1) = du(qq, :) * ggqq * du(qq, :)' ;
            % denom(4) = dv(qq, :) * ggqq * dv(qq, :)' ;
            % % use fact that du and dv are orthogonal
            % denom(2) = sqrt(denom(1) * denom(4));
            % denom(3) = denom(2) ;
            % % Could consider norming by Euclidean norm of the three
            % % components so they are weighted equally

            tg(qq, 1, 1) = du(qq, :) * gq * du(qq, :)' ;
            tg(qq, 1, 2) = du(qq, :) * gq * dv(qq, :)' ;
            tg(qq, 2, 1) = dv(qq, :) * gq * du(qq, :)' ;
            tg(qq, 2, 2) = dv(qq, :) * gq * dv(qq, :)' ;
        end

        % save the metric strain
        readme.eg = ['strain of map from 2D to 3D measured ', ...
            'by ' metric_style ];
        readme.tg = ['strain of 3D bond lengths along s,phi ', ...
            'directions measured by ' metric_style];
        readme.lambda = 'Laplacian smoothing on velocities' ;
        readme.lambda_mesh = 'Laplacian smoothing on mesh vertices' ;
        disp(['saving ', gstrainFn])
        save(gstrainFn, 'eg', 'tg', 'tre', 'dilation', ...
            'lambda', 'lambda_mesh', 'readme')
        
        %% Check the result
        if preview 
            set(gcf, 'visible', 'on')
            bc = barycenter(dbonds.cutMesh.u, dbonds.cutMesh.f) ;
            scatter(bc(:, 1), bc(:, 2), 5, tg(:, 1, 1), 'filled')
            pause(1)
            
            clf
            set(gcf, 'visible', 'on')
            subplot(2, 2, 1); 
            trisurf(tri, eg(:, 1, 1), 'edgeAlpha', 0.01); 
            caxis(max(abs(eg(:, 1, 1))) * [-1, 1])
            zlim([0, 100]); view(0, -90); axis equal; colorbar
            title('$\dot{\epsilon}_{\zeta\zeta}$', 'Interpreter', 'Latex')
            subplot(2, 2, 2); 
            trisurf(tri, eg(:, 1, 2), 'edgeAlpha', 0.01); 
            caxis(max(abs(eg(:, 1, 2))) * [-1, 1])
            zlim([0, 100]); view(0, -90); axis equal; colorbar
            title('$\dot{\epsilon}_{\zeta\phi}$', 'Interpreter', 'Latex')
            subplot(2, 2, 3); 
            trisurf(tri, eg(:, 2, 1), 'edgeAlpha', 0.01); 
            caxis(max(abs(eg(:, 2, 1))) * [-1, 1])
            zlim([0, 100]); view(0, -90); axis equal; colorbar
            title('$\dot{\epsilon}_{\phi\zeta}$', 'Interpreter', 'Latex')
            subplot(2, 2, 4); 
            trisurf(tri, eg(:, 2, 2), 'edgeAlpha', 0.01); 
            caxis(max(abs(eg(:, 2, 2))) * [-1, 1])
            zlim([0, 100]); view(0, 090); axis equal; colorbar
            title('$\dot{\epsilon}_{\phi\phi}$', 'Interpreter', 'Latex')
            % saveas(gcf, fullfile(QS.dir.gstrainRate, 'debug2.png'))
            pause(1)
        end
    else
        % load the metric strain
        load(gstrainFn, 'eg', 'tg', 'tre', 'dilation')
    end        

    %% Plot the metric components on trisurf
    % denom = sqrt(tg(:, 1, 1) .* tg(:, 2, 2)) ;
    labels = {'$\partial_t g_{\zeta\zeta}$', ...
        '$\partial_t g_{\zeta\phi}$', ...
        '$\partial_t g_{\phi\phi}$', ...
        '$\frac{1}{2}\mathrm{Tr}[g^{-1} \dot{g}]$'} ;
    strainlabels = ...
        {'$\dot{\varepsilon}_{\zeta\zeta} \textrm{d}\zeta \textrm{d}\zeta$', ...
        '$\dot{\varepsilon}_{\zeta\phi}\textrm{d}\zeta \textrm{d}\phi$', ...
        '$\dot{\varepsilon}_{\phi\phi} \textrm{d}\phi \textrm{d}\phi $', ...
        '$\frac{1}{2}\mathrm{Tr}[g^{-1} \dot{g}]$'} ;
    glab = {'gss', 'gsphi', 'gphiphi', 'trgdot'} ;
    gelem = [ 1 2 4 ] ;
    time_in_units = (tp - tfold) * QS.timeInterval ;
    tstr = [': $t=$', sprintf('%03d', time_in_units), QS.timeUnits ];

    %% consider each metric element & plot in 3d
    fn = fullfile(egImDir, metric_style, ...
            sprintf([QS.fileBase.spcutMeshSmRSC '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        clf
        set(gcf, 'visible', 'off') ;
        for qq = 1:4
            % For each view (dorsal, ventral, left, right)
            % for pp = 1:4
            subplot(2, 2, qq)
            if qq < 4
                disp(['coloring by tg ' num2str(gelem(qq))])
                colors = tg(:, gelem(qq)) ;
            else
                disp('coloring by trace')
                colors = dilation ;
            end
            trisurf(mesh.f, mesh.v(:, 1), mesh.v(:, 2), mesh.v(:, 3), ...
                colors, 'edgecolor', 'none')
            axis equal
            cb = colorbar() ;

            if strcmp(metric_style, 'mesh')                 
                caxis([-0.5, 0.5])
                title(['surface deformation rate, ', labels{qq}, tstr], ...
                    'Interpreter', 'Latex')            
                ylabel(cb, labels{qq}, 'Interpreter', 'Latex')
            elseif strcmp(metric_style, 'strain')  
                if qq < 4
                    caxis([-clim_tg, clim_tg])
                else
                    caxis([-clim_trgdot, clim_trgdot])
                end
                title(strainlabels{qq}, 'Interpreter', 'Latex')            
                ylabel(cb, strainlabels{qq}, 'Interpreter', 'Latex')
            end
            % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
            % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
            % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
            colormap(bwr256)
            xlim(xyzlim(1, :))
            ylim(xyzlim(2, :))
            zlim(xyzlim(3, :))
            axis off

            % Save images;
            % dorsal
            % if pp == 1
            %     view(0, 90)
            % ventral
            % elseif pp == 2
            %     view(0, -90)
            % left
            % elseif pp == 3
            %     view(0, 0)
            % right
            % elseif pp == 4
            %     view(0, 180)
            % end
            % end
            
            % left view
            view(0, 0) 
        end
        sgtitle(['strain rate, ', tstr], 'Interpreter', 'latex') 

        % Save the image
        saveas(gcf, fn) ;
        clf

    end
    
    %% Now plot in 2d
    close all
    set(gcf, 'visible', 'off') ;
    fn = fullfile(egImDir, ['metricstrain_' metric_style '_2d'], ...
            sprintf([QS.fileBase.spcutMeshSm '.png'], tp));
    if ~exist(fn, 'file') || overwrite
        for qq = 1:4
            subplot(2, 2, qq) 
            if qq < 4
                disp(['coloring by tg ' num2str(gelem(qq))])
                colors = tg(:, gelem(qq)) ;
            else
                disp('coloring by trace')
                colors = dilation ;
            end
            trisurf(cutMesh.f, ...
                cutMesh.u(:, 1) / max(cutMesh.u(:, 1)), ...
                cutMesh.u(:, 2), 0 * cutMesh.u(:, 2), ...
                colors, 'edgecolor', 'none')
            daspect([1,1,1])
            cb = colorbar() ;

            if strcmp(metric_style, 'mesh')                 
                caxis([-0.5, 0.5])
                title(['surface deformation rate, ', labels{qq}, tstr], ...
                    'Interpreter', 'Latex')            
                ylabel(cb, labels{qq}, 'Interpreter', 'Latex')
            elseif strcmp(metric_style, 'strain')  
                if qq < 4
                    caxis([-clim_tg, clim_tg])
                else
                    caxis([-clim_trgdot, clim_trgdot])
                end
                title(['strain rate, ', strainlabels{qq}, tstr], ...
                    'Interpreter', 'Latex')            
                ylabel(cb, strainlabels{qq}, 'Interpreter', 'Latex')
            end
            % xlabel('AP position, [$\mu$m]', 'Interpreter', 'Latex')
            % ylabel('lateral position, [$\mu$m]', 'Interpreter', 'Latex')
            % zlabel('DV position, [$\mu$m]', 'Interpreter', 'Latex')
            colormap(bwr256)
            axis off
            view(2)
        end
        % Save the image
        saveas(gcf, fn) ;
        clf
    end    
    close all
end


%% Compare size of each contribution visually
% Consider each kind of metric strain measurement
for metric_style_kk = 3 
    metric_style = metric_styles{metric_style_kk} ;
    if strcmp(metric_style, 'mesh')
        gstrainBase = QS.fullFileBase.gstrainMesh ;
    elseif strcmp(metric_style, 'strain')
        gstrainBase = QS.fullFileBase.gstrainRate;
    end
    tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;
    
    % build arrays of average magnitude of g11, g12, g22
    g11avg = zeros(1, length(tps)) ;
    g12avg = zeros(1, length(tps)) ;
    g22avg = zeros(1, length(tps)) ;
    
    g11std = zeros(1, length(tps)) ;
    g12std = zeros(1, length(tps)) ;
    g22std = zeros(1, length(tps)) ;
    
    g11min = zeros(1, length(tps)) ;
    g12min = zeros(1, length(tps)) ;
    g22min = zeros(1, length(tps)) ;
    
    g11max = zeros(1, length(tps)) ;
    g12max = zeros(1, length(tps)) ;
    g22max = zeros(1, length(tps)) ;
    
    % Build metric from mesh
    for tp = xp.fileMeta.timePoints(1:end-1)
        disp(['t = ' num2str(tp)])
        tidx = xp.tIdx(tp) ;
        
        % define metric strain filename        
        gstrainFn = sprintf(gstrainBase, tp) ;
        
        % load it
        load(gstrainFn, 'tg')
        
        % Store in linear array for this timepoint
        g11avg(tidx) = mean(sqrt(tg(:, 1, 1).^2)) ;
        g12avg(tidx) = mean(sqrt(tg(:, 1, 2).^2)) ;
        g22avg(tidx) = mean(sqrt(tg(:, 2, 2).^2)) ;
        g11std(tidx) = std(sqrt(tg(:, 1, 1).^2)) ;
        g12std(tidx) = std(sqrt(tg(:, 1, 2).^2)) ;
        g22std(tidx) = std(sqrt(tg(:, 2, 2).^2)) ;
        g11min(tidx) = min(sqrt(tg(:, 1, 1).^2)) ;
        g12min(tidx) = min(sqrt(tg(:, 1, 2).^2)) ;
        g22min(tidx) = min(sqrt(tg(:, 2, 2).^2)) ;
        g11max(tidx) = max(sqrt(tg(:, 1, 1).^2)) ;
        g12max(tidx) = max(sqrt(tg(:, 1, 2).^2)) ;
        g22max(tidx) = max(sqrt(tg(:, 2, 2).^2)) ;
    end
    
    % Draw each magnitude with stdevs
    tps2 = [tps, fliplr(tps)] ; 
    alph = 0.1 ;
    close all 
    fill(tps2, [max(0, g11avg - g11std), fliplr(g11avg + g11std)], ...
        red, 'facealpha', alph, 'edgecolor', 'none') ;
    hold on;
    fill(tps2, [max(0, g12avg - g12std), fliplr(g12avg + g12std)], ...
        blue, 'facealpha', alph, 'edgecolor', 'none') ;
    fill(tps2, [max(0, g22avg - g22std), fliplr(g22avg + g22std)], ...
        green, 'facealpha', alph, 'edgecolor', 'none') ;
    plot(tps, g11avg, '-', 'color', red) ;
    plot(tps, g12avg, '-', 'color', blue) ;
    plot(tps, g22avg, '-', 'color', green) ;    
    legend({'$g_{\zeta\zeta}$', '$g_{\zeta\phi}$', '$g_{\phi\phi}$'},...
        'Interpreter', 'latex')
    xlabel('time [min]')
    ylabel('$|\partial_t g_{ij}|$', 'Interpreter', 'Latex')
    
    if strcmp(metric_style, 'mesh')
        title('Contribution of surface metric strain components')
        fn = fullfile(QS.dir.gstrain, ...
            ['gstrain_' metric_style '_magnitudes']) ;
    else
        title('Contribution of strain components')
        fn = fullfile(QS.dir.gstrain, ...
            ['gstrain_' metric_style '_magnitudes']) ;
    end
    saveas(gcf, [fn '.png']) ;
    ylim([0, 0.2])
    saveas(gcf, [fn '_zoom.png']) ;
    close all
    
end

%% Fit each hoop to Fourier series and add fit to saved .mat file
% Consider each kind of metric strain measurement
for metric_style_kk = 3 
    metric_style = metric_styles{metric_style_kk} ;
    if strcmp(metric_style, 'mesh')
        gstrainBase = QS.fullFileBase.gstrainMesh ;
    elseif strcmp(metric_style, 'strain')
        gstrainBase = QS.fullFileBase.gstrainRate;
    end
    tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;
    
    % Build metric from mesh
    for tp = xp.fileMeta.timePoints(1:end-1)
        disp(['t = ' num2str(tp)])
        tidx = xp.tIdx(tp) ;
        
        % define metric strain filename        
        gstrainFn = sprintf(gstrainBase, tp) ;
        
        % load it
        load(gstrainFn, 'tg')
        
        % Fit each hoop to series
        phi = 2*pi * (0:(nV-2)) / (nV - 1) ;
        
        fitres = struct() ;
        for tcomp = [1,2,4]
            % preallocate fit coefficients for saving
            fitcoeffs = zeros(nU, 3) ;
            for ihoop = 1:nU
                ids = ihoop:nV:nU*(nV-1) ;
                thoop = squeeze(tg(ids, tcomp))' ;

                % Option 1 is curvfitting toolbox
                % ft = fittype('a + b*sin((x - shift))', ...
                %     'coefficients', {'a', 'b', 'shift'}) ;
                % mdl = fit(X,Y,ft,'startpoint',[shiftguess,xscaleguess,yscaleguess]);

                % Instead use raw MATLAB 
                bguess = [0, 0.01, 0] ;
                % Function to fit
                fit = @(b,phi)  b(1) + b(2)*sin(phi + b(3)) ;  
                % Least-Squares cost function
                fcn = @(b) sum((fit(b,phi) - thoop).^2);
                % Minimise Least-Squares
                fitcoeffs(ihoop, :) = fminsearch(fcn, bguess) ;
            end
            if tcomp == 1
                fitres.ss = fitcoeffs ;
            elseif tcomp == 2
                fitres.sphi = fitcoeffs ;
            elseif tcomp == 4
                fitres.phiphi = fitcoeffs ;
            end
        end
        
        % Save fitres
        save(gstrainFn, 'tg', 'fitres') ;
                
    end
end

%% Plot the fits

tps = QS.xp.fileMeta.timePoints(1:end-1) - tfold;

% Build metric from mesh
fss = zeros(length(tps), 3) ;
fsv = zeros(length(tps), 3) ;
fvv = zeros(length(tps), 3) ;
for tp = xp.fileMeta.timePoints(1:end-1)
    disp(['t = ' num2str(tp)])
    tidx = xp.tIdx(tp) ;

    % define metric strain filename        
    gstrainFitFn = sprintf(QS.fullFileBase.gstrainRate, tp) ;
    
    % load it
    load(gstrainFitFn, 'fitres')
    fss = fitres.ss(:, 1:2) ;
    fsv = fitres.sphi(:, 1:2) ;
    fvv = fitres.phiphi(:, 1:2) ;
    
    % domain of parameterization 
    xx = (0:(nU-1))/(nU-1) ;
    
    fn = fullfile(QS.dir.gstrainRateIm, sprintf('ess_fit_%06d.png', tp)) ;
    plot(xx, fss(:, 1), '-'); hold on;
    plot(xx, fss(:, 2), '-')
    legend({'$\langle \varepsilon_{\zeta\zeta} \rangle$', ...
        '$M_1(\varepsilon_{\zeta\zeta})$'}, ...
        'Location', 'northeastoutside', 'Interpreter', 'Latex')
    title('$\varepsilon_{\zeta\zeta}$', 'Interpreter', 'latex')
    xlabel('AP position, $\zeta/L$', 'Interpreter', 'latex') 
    ylabel('moments', 'Interpreter', 'latex')
    ylim([-0.025, 0.025])
    saveas(gcf, fn)
    clf
    
    fn = fullfile(QS.dir.gstrainRateIm, sprintf('esphi_fit_%06d.png', tp)) ;
    plot(xx, fsv, '-'); hold on;
    plot(xx, fsv(:, 2), '-');
    legend({'$\langle \varepsilon_{\zeta\phi} \rangle$', ...
        '$M_1(\varepsilon_{\zeta\phi})$'}, ...
        'Location', 'northeastoutside', 'Interpreter', 'Latex')
    title('$\varepsilon_{\zeta\phi}$', 'Interpreter', 'latex')
    xlabel('AP position, $\zeta/L$', 'Interpreter', 'latex') 
    ylabel('moments', 'Interpreter', 'latex')
    ylim([-0.025, 0.025])
    saveas(gcf, fn)
    clf
        
    fn = fullfile(QS.dir.gstrainRateIm, sprintf('ephiphi_fit_%06d.png', tp)) ;
    plot(xx, fvv(:, 1), '-'); hold on;
    plot(xx, fvv(:, 2), '-');
    legend({'$\langle \varepsilon_{\phi\phi} \rangle$', ...
        '$M_1(\varepsilon_{\phi\phi})$'}, ...
        'Location', 'northeastoutside', 'Interpreter', 'Latex')
    title('$\varepsilon_{\phi\phi}$', 'Interpreter', 'latex')
    xlabel('AP position, $\zeta/L$', 'Interpreter', 'latex') 
    ylabel('moments', 'Interpreter', 'latex')
    ylim([-0.025, 0.025])
    saveas(gcf, fn)
    clf
    
end
