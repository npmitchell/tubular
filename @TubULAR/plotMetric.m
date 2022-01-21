function plotMetric(QS, options)
%plotMetric(QS, options)
%
% Parameters
% ----------
%
% Returns
% -------
%
% NPMitchell 2020

t0 = QS.t0set() ;
lambda_mesh = QS.smoothing.lambda_mesh ;
overwrite = false ;
makeRawMetricComponentFigures = true ;
coordSys = 'spsm_rs' ;
plot_director = false ;
ricciIter = 200 ;

% Unpack options
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
end
if isfield(options, 'ricciIter')
    ricciIter = options.ricciIter ;
end
if isfield(options, 'lambda_mesh')
    lambda_mesh = options.lambda_mesh ;
end
if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'plot_director')
    plot_director = options.plot_director ;
end
if isfield(options, 'makeRawMetricComponentFigures')
    makeRawMetricComponentFigures = options.makeRawMetricComponentFigures ;
end

ntp = length(QS.xp.fileMeta.timePoints) ;
tidx2do = 1:30:ntp ;
tidx2do = [tidx2do, setdiff(1:10:ntp, tidx2do)] ;
tidx2do = [tidx2do, setdiff(1:ntp, tidx2do)] ;
    
for tidx = tidx2do
    tp = QS.xp.fileMeta.timePoints(tidx) ;
    disp(['t = ', num2str(tp)])
    
    strClims = {'climVariable', 'climUniform'} ;
    outdirs = {sprintf(QS.dir.metric.data, coordSys, lambda_mesh), ...
        fullfile(sprintf(QS.dir.metric.g_images3d, coordSys, lambda_mesh), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.b_images3d, coordSys, lambda_mesh), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.g_images3d, coordSys, lambda_mesh), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.b_images3d, coordSys, lambda_mesh), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.g_images2d, coordSys, lambda_mesh), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), strClims{1}), ...
        fullfile(sprintf(QS.dir.metric.g_images2d, coordSys, lambda_mesh), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), strClims{2}), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), 'director'), ...
        fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), 'Hopf_differential')};
    for qq = 1:length(outdirs)
        if ~exist(outdirs{qq}, 'dir')
            mkdir(outdirs{qq}) ;
        end
    end
    
    % Load or compute the fundamental forms
    outfn = sprintf(QS.fullFileBase.metric, coordSys, lambda_mesh, tp) ;
    
    if ~exist(outfn, 'file') || overwrite || true
        % Load mesh
        if strcmpi(coordSys, 'spsm')
            mesh = load(sprintf(QS.fullFileBase.spcutMeshSm, tp), ...
                'spcutMeshSm') ;
            mesh = mesh.spcutMeshSm ;

            % Project to APDV coords if not already done
            try 
                mesh.v = mesh.vrs;
            catch
                mesh.v = QS.xyz2APDV(mesh.v) ;
            end
            mesh_exists = true ;
        elseif strcmpi(coordSys, 'spsm_rs')
            mesh = load(sprintf(QS.fullFileBase.spcutMeshSmRS, tp), ...
                'spcutMeshSmRS') ;
            mesh = mesh.spcutMeshSmRS ;
            mesh_exists = true ;
        elseif strcmpi(coordSys, 'ricci')
            try
                mesh = load(sprintf(QS.fullFileBase.ricciMesh, ricciIter, tp), ...
                    'ricciMesh') ;
                mesh = mesh.ricciMesh.rectangle ;
                mesh_exists = true ;
            catch
                mesh_exists = false ;
            end
        else
            error('handle this coordSys here')
        end
                
        if mesh_exists
            % Smooth the mesh with lambda_mesh
            if lambda_mesh > 0 
                disp('smoothing mesh vertices before computations')
                glueMesh = glueRectCylinderCutMeshSeam(mesh) ;
                tri = triangulation(glueMesh.f, glueMesh.v) ;
                fbndy = tri.freeBoundary ;
                fbndy = fbndy(:, 1) ;
                glueMesh.v = laplacian_smooth(glueMesh.v, glueMesh.f, 'cotan', fbndy, ...
                    lambda_mesh, 'implicit', glueMesh.v) ;
                opts = struct() ;
                opts.ignoreRectangularConstraint = true ;
                opts.vmax = 2 * pi ;
                mesh = cutRectilinearCylMesh(glueMesh, opts) ;
            else
                mesh.v = mesh.v ;
            end

            % Rescale u dimension for mu relaxation if not Ricci mesh
            if strcmpi(coordSys, 'ricci')
                mesh.mu = bc_metric(mesh.f, mesh.u, mesh.v, 3) ;
                aspectShear = 1.0 ;
            else
                uv = mesh.u ;
                uv(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) ;
                mesh.readme = 'mu is computed after rescaling refMesh.u(:, 1) to range from 0 to ar, as determined by mean pure shear in (s,phi)' ;
                mesh.mu = bc_metric(mesh.f, uv, mesh.v, 3) ;

                % Rescale and recompute mu
                aspectShear = (1 + mean(real(mesh.mu))) / (1 - mean(real(mesh.mu))) ;
                mesh.u(:, 1) = mesh.u(:, 1) / max(mesh.u(:, 1)) * aspectShear ; 
                mesh.mu = bc_metric(mesh.f, mesh.u, mesh.v, 3) ;
            end

            [gcell, bcell] = constructFundamentalForms(mesh.f, mesh.v, mesh.u) ;
            gg = zeros(length(gcell), 4) ;
            bb = zeros(length(gcell), 4) ;
            for qq = 1:length(gcell)
                gg(qq, 1) = gcell{qq}(1, 1) ;
                gg(qq, 2) = gcell{qq}(1, 2) ;
                gg(qq, 3) = gcell{qq}(2, 1) ;
                gg(qq, 4) = gcell{qq}(2, 2) ;
                bb(qq, 1) = bcell{qq}(1, 1) ;
                bb(qq, 2) = bcell{qq}(1, 2) ;
                bb(qq, 3) = bcell{qq}(2, 1) ;
                bb(qq, 4) = bcell{qq}(2, 2) ;
            end

            % glue mesh to topological cylinder
            glueMesh = glueCylinderCutMeshSeam(mesh) ;
            assert(eulerCharacteristic(glueMesh) == 0)
            [~, F2V] = meshAveragingOperators(glueMesh) ;

            % Get director as Hopf differential on faces
            LL = bb(:, 1) ;
            MM = bb(:, 2) * 0.5 ;
            NN = bb(:, 4) ;
            Hopf_faces = 0.25 * ((LL-NN) - 2*1j*MM ) ;

            % Check with bars
            % bc = barycenter(mesh.u, mesh.f) ;
            % bar1 = real(QQ) ;
            % bar1 = cat(2, bar1, imag(QQ)) ;
            % quiver(bc(:, 1), bc(:, 2), bar1(:, 1), bar1(:, 2), 1) ;

            % Check with scatter
            % subplot(2, 1, 1) ;
            % scatter(bc(:, 1), bc(:, 2), 10, real(QQ))
            % subplot(2, 1, 2) ;
            % scatter(bc(:, 1), bc(:, 2), 10, imag(QQ))

            % Get director on vertices
            nU = mesh.nU ;
            nV = mesh.nV ;
            QV = zeros(length(mesh.u), 2) ;
            evs = QV ;
            bbV = zeros(length(mesh.u), 4) ;
            bbV = F2V * bb ;
            bbV(nU*(nV-1)+1:nU*nV, :) = bbV(1:nU, :) ;

            LV = bbV(:, 1) ;
            MV = bbV(:, 2) * 0.5 ;
            NV = bbV(:, 4) ;
            QV = 0.25 * ((LV-NV) - 2*1j*MV ) ;

            % Check against moving QV directly -- checks out!
            % QV2 = F2V * QQ ;
            % QV2(nU*(nV-1)+1:nU*nV, :) = QV(1:nU, :) ;

            % Filter a bit
            QVs_r = imgaussfilt(reshape(real(QV), [nU, nV]), 1) ;
            QVs_i = imgaussfilt(reshape(imag(QV), [nU, nV]), 1) ;
            QVs = QVs_r + 1j * QVs_i ;

            mags = abs(QVs) ;
            thetas = atan2(imag(QVs), real(QVs)) ;
            uu = reshape(mesh.u(:, 1), [nU, nV]) ;
            vv = reshape(mesh.u(:, 2), [nU, nV]) ;
            Hopf_vertices = reshape(QVs, [nU, nV]) ;

            % Save fundamental forms with mesh
            save(outfn, 'gg', 'bb', 'mesh', 'aspectShear', ...
                'Hopf_faces', 'Hopf_vertices')
        end
    else
        load(outfn, 'gg', 'bb', 'mesh', 'aspectShear', ...
            'Hopf_faces', 'Hopf_vertices')
        mesh_exists = true ;
    end
    
    if mesh_exists
        %% Second fundamental form --> mean curvature and Hopf differential
        fnbb = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
            'Hopf_differential', sprintf(QS.fileBase.name, tp)) ;        
        fnbb = [fnbb '_b_fundForm_Hopf.png'] ;

        if ~exist(fnbb, 'file') || overwrite || true
            close all
            opts = struct() ;

            mags = abs(Hopf_vertices) ;
            thetas = atan2(imag(Hopf_vertices), real(Hopf_vertices)) ;
            nU = mesh.nU ;
            nV = mesh.nV ;
            uu = reshape(mesh.u(:, 1), [nU, nV]) ;
            vv = reshape(mesh.u(:, 2), [nU, nV]) ;

            % Check with polar field
            % opts = struct() ;
            % opts.clim_mag = rms1d(mags(:)) ;
            % plotPolarField(mags', thetas', opts)

            % Plot as polar in 2d AND 3d
            m2d = struct() ;
            m2d.f = mesh.f ;
            m2d.v = mesh.u ;
            [~,~,~,xyzlim] = QS.getXYZLims() ;
            m2d.v(:, 3) = 0 * mesh.u(:, 1) ;
            if strcmpi(coordSys, 'ricci')
                opts.xyzlims = {xyzlim, ...
                    [0,max(mesh.u(:, 1));0,max(mesh.u(:, 2));-0.1,0.1]} ;
            else
                opts.xyzlims = {xyzlim, [0,aspectShear;0,1;-0.1,0.1]} ;
            end
            opts.clim = rms1d(mags(:)) ;
            opts.axisOff = true ;
            opts.labels = {'', ''} ;
            opts.views = {[0, 0 ], [0, 90 ]};
            opts.cbarlabels = {['Hopf differential, $|Q|$ [' QS.spaceUnits '$^{-2}$]'], ...
                               ['Hopf differential, $|Q|$ [' QS.spaceUnits '$^{-2}$]']} ;
            opts.polarStyle = 'polar';
            close all
            nFieldsOnSurface({mesh, m2d}, ...
                {{mags, thetas}, {mags, thetas}}, opts) ;
            sgtitle({'Hopf differential', ...
                ['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                    ' ', QS.timeUnits]}, 'Interpreter', 'latex')
            disp(['Saving ' fnbb])
            saveas(gcf, fnbb)
            close all

            % SQUARE ROOT OF HOPF
            opts.clim = 1.2 * rms1d(sqrt(mags(:))) ;
            opts.labels = {'', ''} ;
            opts.cbarlabels = {['Hopf differential, $|\sqrt{Q}|$ [' QS.spaceUnits '$^{-1}$]'], ...
                               ['Hopf differential, $|\sqrt{Q}|$ [' QS.spaceUnits '$^{-1}$]']} ;
            fnbb2 = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
                'Hopf_differential', sprintf(QS.fileBase.name, tp)) ;        
            fnbb2 = [fnbb2 '_b_fundForm_Hopf_halfAngle.png'] ;
            opts.polarStyle = 'nematic';
            axs = nFieldsOnSurface({mesh, m2d}, ...
                {{sqrt(mags), thetas*0.5}, {sqrt(mags), thetas*0.5}}, opts) ;
            hold on; 

            % Plot orientation in 2d as well
            set(gcf, 'CurrentAxes', axs{2})
            cosD = imresize(reshape(cos((thetas) * 0.5), [nU, nV]), 0.25 ) ;
            sinD = imresize(reshape(sin((thetas) * 0.5), [nU, nV]), 0.25 ) ;
            magD = imresize(reshape(sqrt(mags), [nU, nV]), 0.25 ) ;
            xxD = imresize(reshape(m2d.v(:, 1), [nU, nV]), 0.25) ;
            yyD = imresize(reshape(m2d.v(:, 2), [nU, nV]), 0.25) ;
            quiver(xxD(:), yyD(:), magD(:) .* cosD(:), magD(:) .* sinD(:), 1, ...
                'w', 'ShowArrowHead', 'off')
            quiver(xxD(:), yyD(:), -magD(:) .* cosD(:), -magD(:) .* sinD(:), 1, ...
                'w', 'ShowArrowHead', 'off')
            if strcmpi(coordSys, 'ricci')
                xlim([0, max(mesh.u(:, 1))]) ;
                ylim([0, max(mesh.u(:, 2))]) ;
                zlim([-0.1, 0.1]) ;
            else
                xlim([0, aspectShear]) ;
                ylim([0, 1]) ;
                zlim([-0.1, 0.1]) ;
            end

            fieldfaces = pointLocation(triangulation(m2d.f, m2d.v(:, 1:2)), ...
                xxD(:), yyD(:)) ;
            n3d = pushVectorField2Dto3DMesh([cosD(:), sinD(:)], m2d.v, ...
                mesh.v, mesh.f, fieldfaces) ;
            bc3 = barycenter(mesh.v, mesh.f) ;

            % Plot orientation in 3d as well
            set(gcf, 'CurrentAxes', axs{1})
            quiver3(bc3(fieldfaces, 1), bc3(fieldfaces, 2), bc3(fieldfaces, 3), ...
                magD(:) .* n3d(:,1), magD(:) .* n3d(:,2), magD(:) .* n3d(:,3), ...
                1, 'w', 'ShowArrowHead', 'off')
            quiver3(bc3(fieldfaces, 1), bc3(fieldfaces, 2), bc3(fieldfaces, 3), ...
                -magD(:) .* n3d(:,1), -magD(:) .* n3d(:,2), -magD(:) .* n3d(:,3), ...
                1, 'w', 'ShowArrowHead', 'off')
            xlim(xyzlim(1, :)) ;
            ylim(xyzlim(2, :)) ;
            zlim(xyzlim(3, :)) ;


            sgtitle({'Square root of Hopf differential', ... % , rotated by $\pi/2$', ...
                ['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                    ' ', QS.timeUnits]}, 'Interpreter', 'latex')
            disp(['Saving ' fnbb2])
            set(gcf,'color','w');
            export_fig(fnbb2, '-nocrop', '-r150')
            close all

        end


        %% SAVE images of second fundamental form director if they dont exist
        fnbb = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
            'director', sprintf(QS.fileBase.name, tp)) ;        
        fnbb = [fnbb '_b_fundForm_director.png'] ;

        if (~exist(fnbb, 'file') || overwrite) && plot_director
            close all
            opts = struct() ;

            % glue mesh to topological cylinder
            glueMesh = glueCylinderCutMeshSeam(mesh) ;
            assert(eulerCharacteristic(glueMesh) == 0)
            [~, F2V] = meshAveragingOperators(glueMesh) ;

            % Get director
            QQ = zeros(length(bb), 2) ;
            evs = QQ ;
            for qq = 1:length(bb) 
                b2d = reshape(bb(qq, :), [2, 2]) ;
                [eigvect, eigval] = eig(b2d) ;
                % larger eigvect direction
                [~, id] = max([eigval(1, 1), eigval(2, 2)]) ;
                if sign(eigval(id, id)) > 0
                    QQ(qq, :) = eigvect(:, id) ;
                else
                    QQ(qq, :) = -eigvect(:, id) ;
                end
                % Note: big, then small
                otherId = setdiff([1,2], id) ;
                evs(qq, :) = [eigval(id, id) eigval(otherId, otherId)] ;
            end
            bc = barycenter(mesh.u, mesh.f) ;
            bar1 = evs(:, 1) .* QQ(:, 1) ;
            bar1 = cat(2, bar1, evs(:, 1) .* QQ(:, 2)) ;
            quiver(bc(:, 1), bc(:, 2), bar1(:, 1), bar1(:, 2), 1)

            % Get director on vertices
            nU = mesh.nU ;
            nV = mesh.nV ;
            QV = zeros(length(mesh.u), 2) ;
            evs = QV ;
            bbV = zeros(length(mesh.u), 4) ;
            bbV = F2V * bb ;
            bbV(nU*(nV-1)+1:nU*nV, :) = bbV(1:nU, :) ;
            for qq = 1:length(bbV) 
                b2d = reshape(bbV(qq, :), [2, 2]) ;
                [eigvect, eigval] = eig(b2d) ;

                % check that the tensor is still symmetric even though it is
                % averaged onto vertices
                assert(b2d(1, 2) == b2d(2, 1))

                % larger eigvect direction
                [~, id] = max([eigval(1, 1), eigval(2, 2)]) ;
                if sign(eigval(id, id)) > 0
                    QV(qq, :) = eigvect(:, id) ;
                else
                    QV(qq, :) = -eigvect(:, id) ;
                end
                % Note: big, then small
                otherId = setdiff([1,2], id) ;
                if sign(eigval(id, id)) > 0
                    evs(qq, :) = [eigval(id, id) eigval(otherId, otherId)] ;
                else
                    evs(qq, :) = -[eigval(id, id) eigval(otherId, otherId)] ;
                end
            end

            %% Bar c2t, s2t 
            thet = atan2(QV(:, 2), QV(:, 1)) ;
            c2t = evs(:, 1) .* cos(2 * thet) ;
            s2t = evs(:, 1) .* sin(2 * thet) ;
            c2t = imgaussfilt(c2t, 1) ;
            s2t = imgaussfilt(s2t, 1) ;
            uu = reshape(mesh.u(:, 1), [nU, nV]) ;
            vv = reshape(mesh.u(:, 2), [nU, nV]) ;

            assert(all(evs(:, 1) > 0))
            mags = vecnorm([c2t, s2t], 2, 2) ;
            thetas = mod(atan2(s2t, c2t), pi)  ;
            m2d = struct() ;
            m2d.f = mesh.f ;
            m2d.v = mesh.u ;
            m2d.v(:, 3) = 0 * mesh.u(:, 1) ;
            opts.clim = rms1d(evs(:, 1)) ;
            opts.axisOff = true ;
            opts.labels = {'', '2nd fundamental form as nematic'} ;
            opts.views = {[0, 0 ], [0, 90 ]};
            opts.cbarlabels = {'$b$ eigenvalue, $||b_1||$', ...
                               '$b$ eigenvalue, $||b_1||$'} ;
            nFieldsOnSurface({mesh, m2d}, ...
                {{mags, thetas}, {mags, thetas}}, opts) ;

            sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                    ' ', QS.timeUnits], 'Interpreter', 'latex')
            disp(['Saving ' fnbb])
            saveas(gcf, fnbb)
            close all

            %% Bar X,Y
            % bar1V = evs(:, 1) .* QV(:, 1) ;
            % bar1V = cat(2, bar1V, evs(:, 1) .* QV(:, 2)) ;
            % bQX = reshape(bar1V(:, 1), [nU, nV]) ;
            % bQY = reshape(bar1V(:, 2), [nU, nV]) ;
            % % bQXr = imresize(bQX, [0.5 * nU, NaN]) ;
            % % bQYr = imresize(bQY, [0.5 * nU, NaN]) ;
            % % uu = imresize(reshape(mesh.u(:, 1), [nU, nV]), [0.5 * nU, NaN]) ;
            % % vv = imresize(reshape(mesh.u(:, 2), [nU, nV]), [0.5 * nU, NaN]) ;
            % bQXr = imgaussfilt(bQX, 1) ;
            % bQYr = imgaussfilt(bQY, 1) ;
            % uu = reshape(mesh.u(:, 1), [nU, nV]) ;
            % vv = reshape(mesh.u(:, 2), [nU, nV]) ;
            % % quiver(mesh.u(:, 1), mesh.u(:, 2), bar1V(:, 1), bar1V(:, 2), 1)
            % quiver(uu, vv, bQXr, bQYr, 1) ;
            % hold off;

            % stream-like plot
            % flipYQ = bar1V(:, 2) < 0 ;
            % bar1V(flipYQ, :) = -bar1V(flipYQ, :) ;
            % streamline(uu', vv',  bQXr, bQYr, uu(1:10:end), vv(1:10:end)) ;
            % 
            % mags = vecnorm([bQXr(:), bQYr(:)], 2, 2) ;
            % thetas = mod(atan2(bQYr(:), bQXr(:)), pi) ;
            % m2d = struct() ;
            % m2d.f = mesh.f ;
            % m2d.v = mesh.u ;
            % m2d.v(:, 3) = 0 * mesh.u(:, 1) ;
            % opts.clim = 0.2 ;
            % opts.view = [0, 90 ];
            % nFieldsOnSurface({mesh, m2d}, ...
            %     {{mags, thetas}, {mags, thetas}}, opts) ;
        end

        % SAVE images of first and second fundamental forms if they dont exist
        files_exist = true ;
        for pp = 1:2
            fn1a = fullfile(sprintf(QS.dir.metric.g_images3d, coordSys, lambda_mesh), ...
                 strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
            fn1a = [fn1a '_g_fundForm.png'] ;

            fn1b = fullfile(sprintf(QS.dir.metric.g_images2d, coordSys, lambda_mesh), ...
                 strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
            fn1b = [fn1b '_g_fundForm.png'] ;


            fn2a = fullfile(sprintf(QS.dir.metric.b_images3d, coordSys, lambda_mesh), ...
                 strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
            fn2a = [fn2a '_b_fundForm.png'] ;
            fn2b = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
                 strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
            fn2b = [fn2b '_b_fundForm.png'] ;

            files_exist = files_exist && exist(fn1a, 'file') && ...
                exist(fn1b, 'file') && exist(fn2a, 'file') && ...
                exist(fn2b, 'file') ;
        end
        if ~files_exist && makeRawMetricComponentFigures
            for pp = 1:2
                opts = struct() ;
                if pp == 1
                    opts.clims = {max(abs(gg(:, 1))) * [-1, 1], ...
                        max(abs(gg(:, 2))) * [-1, 1], ...
                        max(abs(gg(:, 3))) * [-1, 1], ...
                        max(abs(gg(:, 4))) * [-1, 1]} ;
                else
                    opts.clim = max(abs(gg(:))) * [-1, 1] ;
                end

                % 3d case ---------------------------------------------------------
                labels = {'$\mathbf{g}_{\zeta\zeta}$', ...
                    '$\mathbf{g}_{\zeta\phi}$', ...
                    '$\mathbf{g}_{\phi\zeta}$', ...
                    '$\mathbf{g}_{\phi\phi}$'} ;
                m2view = mesh ;
                m2view.v = mesh.vrs ;
                opts.labels = labels ;
                [~, ~, ~, opts.xyzlims] = QS.getXYZLims() ;
                % gg 3d case ------------------------------------------------------
                close all
                [axs, cbs, meshHandles] = ...
                    nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
                    {gg(:, 1), gg(:, 2), gg(:, 3), gg(:, 4)}, opts) ;
                for qq = 1:4
                    set(gcf,'CurrentAxes', axs{qq})
                    view(0, 0)
                    axis off
                end
                sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                    ' ', QS.timeUnits], 'Interpreter', 'latex')
                fn = fullfile(sprintf(QS.dir.metric.g_images3d, coordSys, lambda_mesh), ...
                     strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
                fn = [fn '_g_fundForm.png'] ;
                saveas(gcf, fn)

                % gg 2d case ------------------------------------------------------
                close all
                m2d = mesh ;
                m2d.v = [mesh.u(:, 1) / aspectShear, mesh.u(:, 2), mesh.u(:, 1)*0] ;
                opts.labels = labels ;
                opts = rmfield(opts, 'xyzlims') ;
                % gg
                [axs, cbs, meshHandles] = ...
                    nFieldsOnSurface({m2d, m2d, m2d, m2d}, ...
                    {gg(:, 1), gg(:, 2), gg(:, 3), gg(:, 4)}, opts) ;
                for qq = 1:4
                    set(gcf,'CurrentAxes', axs{qq})
                    view(2)
                    axis off
                end
                sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                    ' ', QS.timeUnits], 'Interpreter', 'latex')
                fn = fullfile(sprintf(QS.dir.metric.g_images2d, coordSys, lambda_mesh), ...
                     strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
                fn = [fn '_g_fundForm.png'] ;
                saveas(gcf, fn)

                % bb 3d case ------------------------------------------------------
                close all
                opts = struct() ;
                if pp == 1
                    opts.clims = {max(abs(bb(:, 1))) * [-1, 1], ...
                        max(abs(bb(:, 2))) * [-1, 1], ...
                        max(abs(bb(:, 3))) * [-1, 1], ...
                        max(abs(bb(:, 4))) * [-1, 1]} ;
                else
                    opts.clim = max(abs(bb(:))) * [-1, 1] ;
                end
                labels = {'$\mathbf{b}_{\zeta\zeta}$', ...
                    '$\mathbf{b}_{\zeta\phi}$', ...
                    '$\mathbf{b}_{\phi\zeta}$', ...
                    '$\mathbf{b}_{\phi\phi}$'} ;
                opts.labels = labels ;
                [~, ~, ~, opts.xyzlims] = QS.getXYZLims() ;
                [axs, cbs, meshHandles] = ...
                    nFieldsOnSurface({m2view, m2view, m2view, m2view}, ...
                    {bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4)}, opts) ;
                % title
                sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                    ' ', QS.timeUnits], 'Interpreter', 'latex')
                fn = fullfile(sprintf(QS.dir.metric.b_images3d, coordSys, lambda_mesh), ...
                     strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
                fn = [fn '_b_fundForm.png'] ;
                for qq = 1:4
                    set(gcf,'CurrentAxes', axs{qq})
                    view(0, 0)
                    axis off
                end
                saveas(gcf, fn)

                % bb 2d case -------------------------------------------------
                close all
                opts = rmfield(opts, 'xyzlims') ;
                [axs, cbs, meshHandles] = ...
                    nFieldsOnSurface({m2d, m2d, m2d, m2d}, ...
                    {bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4)}, opts) ;
                for qq = 1:4
                    set(gcf,'CurrentAxes', axs{qq})
                    view(2)
                    axis off
                end
                % title
                sgtitle(['$t = $', sprintf('%03d', (tp - t0)*QS.timeInterval), ...
                    ' ', QS.timeUnits], 'Interpreter', 'latex')
                fn = fullfile(sprintf(QS.dir.metric.b_images2d, coordSys, lambda_mesh), ...
                     strClims{pp}, sprintf(QS.fileBase.name, tp)) ;        
                fn = [fn '_b_fundForm.png'] ;
                saveas(gcf, fn)
            end
        end
    end
    
end
