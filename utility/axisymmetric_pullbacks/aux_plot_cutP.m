disp('Plotting cut...')
xyzrs = ((rot * mesh.v')' + trans) * resolution ;
fig = figure('Visible', 'Off')  ;
fig.PaperUnits = 'centimeters';
% auxilliary function for Generate_Axisymmetric_Pullbacks_Orbifold.m

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
disp(['Saving cutpath figure to ' cutfn])
saveas(fig, cutfn)
close all 