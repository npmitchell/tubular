% Auxilliary function for plotting cleanCylinderMesh from
% Generate_Axisymmetric_Pullbacks_Orbifold.m script
% NOTE: THIS IS NOT USED IN THE MASTER PIPELINE, deprecated

close all
fig = figure('visible', 'off') ;
xyzrs = ((rot * mesh.v')' + trans) * resolution ;
trisurf( triangulation( mesh.f, xyzrs ), xyzrs(:, 2), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on
scatter3( xyzrs(adIDx,1), xyzrs(adIDx,2), xyzrs(adIDx,3), ...
    'filled', 'r' );
scatter3( xyzrs(pdIDx,1), xyzrs(pdIDx,2), xyzrs(pdIDx,3), ...
    'filled', 'c' );
hold off
axis equal
xlabel('x [\mum]')
ylabel('y [\mum]')
zlabel('z [\mum]')
xlim(xyzlim(1, :)); 
ylim(xyzlim(2, :)); 
zlim(xyzlim(3, :));
title(['cleaned cylinder mesh, t = ' num2str(t)])
saveas(fig, mesh3dfigfn)
close all