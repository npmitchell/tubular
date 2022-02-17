function aux_plot_clineDVhoop(QS, avgpts, avgpts_ss, cseg, cline, ...
    cseg_ss, curves3d, xyzlim, clineDVhoopFigBase, t)
% aux_plot_clineDVhoop: auxilliary function for
% Generate_Axisymmetric_Pullbacks_Orbifold script
%
% Plots the DVhoop as a scatter plot in 3d.
%
% cline has already been transformed to APDV
%
% NPMitchell 2019

% Should really transform avgpts back into XYZ coords, then into APDV.
% Otherwise we are off by a mirroring.
% if QS.flipy
%     avgpts(:, 2) = - avgpts(:, 2) ;
%     curves3d(:, 2) = - curves3d(:, 2) ;
% end

close all
fig = figure('visible', 'off') ;
scatter3(avgpts(:, 1), avgpts(:, 2), avgpts(:, 3), 10, avgpts_ss) ;
hold on;
scatter3(cseg(:, 1), cseg(:, 2), cseg(:, 3), 10, cseg_ss)
plot3(cline(:, 1), cline(:, 2), cline(:, 3), 'k');
% trisurf(triangulation(TF, TV3D), 'EdgeColor', 'none', 'FaceAlpha', 0.1); 
cmap = colormap ;
colorval = cmap(max(1, round(avgpts_ss/ max(avgpts_ss) * length(cmap))), :) ;
for jj = 1:size(curves3d, 1)
    s = scatter3(curves3d(jj, :, 1), curves3d(jj, :, 2), ...
        curves3d(jj, :, 3), 2, 'MarkerFaceColor', colorval(jj, :), 'MarkerEdgeColor', 'none') ;
    alpha(s,.2)
end
axis equal
xlim(xyzlim(1, :))
ylim(xyzlim(2, :))
zlim(xyzlim(3, :))
xlabel('x [\mum]')
ylabel('y [\mum]')
zlabel('z [\mum]')
fn = sprintf(clineDVhoopFigBase, t) ;
disp(['Saving image of new centerline to ' fn]) 
saveas(fig, fn);
close all