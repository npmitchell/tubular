function plotCutPath(QS, cutMesh, cutP)
%
% Plot the cut path of a cylinderCutMesh on the mesh in 3D and save image
%
% Parameters
% ----------
% QS: QuapSlap object
% cutMesh : cutMesh object with fields 
%   v : vertices
%   f : face connectivity list
% cutP : optional 3D path of the cut on the cylinderCutMesh
%
% NPMitchell 2019

% unpack options
if nargin < 2 || isempty(cutMesh)
    if isempty(QS.currentMesh.cutMesh)
        QS.loadCurrentCutMesh()
    end
    cutMesh = QS.currentMesh.cutMesh ;
end
if nargin < 3 || isempty(cutP)
    cutP = QS.currentMesh.cutPath ;
    if isempty(cutP)
        QS.loadCurrentCutMesh()
    end
    cutP = QS.currentMesh.cutPath ;
end

outdir = fullfile(QS.dir.cutMesh, 'images') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

tt = QS.currentTime;

% Now plot
disp('Plotting cut...')
xyzrs = QS.xyz2APDV(cutMesh.v) ;
fig = figure('Visible', 'Off')  ;
fig.PaperUnits = 'centimeters';
[~, ~, ~, xyzlim_um] = QS.getXYZLims() ;
xyzlim_um(:, 1) = xyzlim_um(:, 1) - QS.normalShift ;
xyzlim_um(:, 1) = xyzlim_um(:, 1) + QS.normalShift ;

% Figure generation
fig.PaperPosition = [0 0 12 12];
sh = trimesh(cutMesh.f, ...
    xyzrs(:, 1), xyzrs(:,2), xyzrs(:, 3), xyzrs(:, 2), ...
    'FaceColor', 'interp', 'edgecolor', 'none', 'FaceAlpha', 0.3) ;
hold on;
ph = plot3(xyzrs(cutP, 1), xyzrs(cutP, 2), xyzrs(cutP, 3), 'k-', 'LineWidth', 3) ;
axis equal
xlim(xyzlim_um(1, :)); 
ylim(xyzlim_um(2, :)); 
zlim(xyzlim_um(3, :));
title(['t=' sprintf('%04d', tt)]) ;
xlabel('x [$\mu$m]', 'Interpreter', 'Latex') ;
ylabel('y [$\mu$m]', 'Interpreter', 'Latex') ;
zlabel('z [$\mu$m]', 'Interpreter', 'Latex') ;
cutfn = sprintf( fullfile(outdir, [QS.fileBase.name, '_cut.png']), tt ) ;
disp(['Saving cutpath figure to ' cutfn])
saveas(fig, cutfn)
close all 
