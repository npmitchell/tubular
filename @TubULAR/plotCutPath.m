function plotCutPath(tubi, cutMesh, cutP)
%
% Plot the cut path of a cylinderCutMesh on the mesh in 3D and save image
%
% Parameters
% ----------
% tubi: TubULAR class instance object
% cutMesh : cutMesh object with fields 
%   v : vertices
%   f : face connectivity list
% cutP : optional 3D path of the cut on the cylinderCutMesh
%
% NPMitchell 2019

% unpack options
if nargin < 2 || isempty(cutMesh)
    if isempty(tubi.currentMesh.cutMesh)
        tubi.loadCurrentCutMesh()
    end
    cutMesh = tubi.currentMesh.cutMesh ;
end
if nargin < 3 || isempty(cutP)
    cutP = tubi.currentMesh.cutPath ;
    if isempty(cutP)
        tubi.loadCurrentCutMesh()
    end
    cutP = tubi.currentMesh.cutPath ;
end

outdir = fullfile(tubi.dir.cutMesh, 'images') ;
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

tt = tubi.currentTime;

% Now plot
disp('Plotting cut...')
xyzrs = tubi.xyz2APDV(cutMesh.v) ;
fig = figure('Visible', 'Off')  ;
fig.PaperUnits = 'centimeters';
[~, ~, ~, xyzlim_um] = tubi.getXYZLims() ;
xyzlim_um(:, 1) = xyzlim_um(:, 1) - tubi.normalShift ;
xyzlim_um(:, 1) = xyzlim_um(:, 1) + tubi.normalShift ;

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
cutfn = sprintfm( fullfile(outdir, [tubi.fileBase.name, '_cut.png']), tt ) ;
disp(['Saving cutpath figure to ' cutfn])
saveas(fig, cutfn)
close all 
