% Plot the cutMesh3D
            % Load the mesh
            disp(['NOW PROCESSING TIME POINT ', num2str(t)]);
            tidx = xp.tIdx(t);

            % Load the cylinder mesh
            cylmeshfn = sprintf( cylinderMeshCleanBase, t ) ;
            mesh = read_ply_mod( cylmeshfn );
            
            % Load the cut from h5 file
            cutP = dlmread(sprintf(outcutfn, t)) ;
            error('check that cutP was loaded correctly here')
            
            %% Plot the cutPath (cutP) in 3D
            disp('Plotting cut...')
            xyzrs = ((rot * mesh.v')' + trans) * resolution ;
            fig = figure('Visible', 'Off')  ;
            fig.PaperUnits = 'centimeters';
            fig.PaperPosition = [0 0 12 12];
            sh = trimesh(mesh.f, ...
                xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), xyzrs(:, 2), ...
                'FaceColor', 'interp', 'edgecolor', 'none', 'FaceAlpha', 0.3) ;
            hold on;
            ph = plot3(xyzrs(cutP, 1), xyzrs(cutP, 2), xyzrs(cutP, 3), 'k-', 'LineWidth', 3) ;
            axis equal
            xlim(xyzlim(1, :)); 
            ylim(xyzlim(2, :)); 
            zlim(xyzlim(3, :));
            title(['t=' sprintf('%04d', t)]) ;
            xlabel('x [$\mu$m]', 'Interpreter', 'Latex') ;
            ylabel('y [$\mu$m]', 'Interpreter', 'Latex') ;
            zlabel('z [$\mu$m]', 'Interpreter', 'Latex') ;
            cutfn = sprintf( fullfile(cutMeshImagesDir, [fileNameBase, '_cut.png']), t ) ;
            saveas(fig, cutfn)
            close all 