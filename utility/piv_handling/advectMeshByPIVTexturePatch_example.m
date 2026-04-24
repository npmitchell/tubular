% Tiny demo for advectMeshByPIVTexturePatch.m
% Run this script from anywhere after adding tubular to your MATLAB path.

addpath(genpath(fullfile('..', '..', 'TexturePatch')));
addpath(genpath(fullfile('..', '..', 'utility')));

clear; close all; clc;

% Ensure TexturePatch utilities are available on the MATLAB path.
if isempty(which('texture_patch_to_image'))
    thisDir = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(thisDir));
    tpDir = fullfile(repoRoot, 'TexturePatch');
    if exist(fullfile(tpDir, 'texture_patch_to_image.m'), 'file')
        addpath(tpDir);
    end
end

if isempty(which('texture_patch_to_image'))
    error(['Could not find texture_patch_to_image.m. Add ', ...
        '/TexturePatch to your MATLAB path.']);
end

% Use a MATLAB sample image and convert to grayscale.
Irgb = imread('peppers.png');
I = im2double(rgb2gray(Irgb));
[H, W] = size(I);

% Build a simple pullback mesh in pixel space.
nX = 40;
nY = 30;
[xv, yv] = meshgrid(linspace(1, W, nX), linspace(1, H, nY));
V2Dpix = [xv(:), yv(:)];
FF = delaunay(V2Dpix(:, 1), V2Dpix(:, 2));

% Build texture coordinates in [row, col, page] format.
% We use a 3D grayscale stack so TV3D can include a page index.
nPages = 5;
IV0 = repmat(I, 1, 1, nPages);
TV3D = [V2Dpix(:, 2), V2Dpix(:, 1), 3 * ones(size(V2Dpix, 1), 1)];

% Synthetic PIV field on a coarse grid (pixels/frame).
[xPiv, yPiv] = meshgrid(linspace(1, W, 50), linspace(1, H, 40));
uPiv = 8 * sin(2 * pi * yPiv / H);
vPiv = 5 * cos(2 * pi * xPiv / W);

% Optional comparison images for RGB QC overlay.
imCurrent = I;
imNext = circshift(I, [0, 8]);

% Texture patch options.
patchOpts = struct();
patchOpts.imSize = [H, W];
patchOpts.xLim = [1, W];
patchOpts.yLim = [1, H];

opts = struct();
opts.patchOpts = patchOpts;
opts.imCurrent = imCurrent;
opts.imNext = imNext;
opts.makeOverlay = true;

% Run advection + texture patch rendering.
[patchImAdv, meshAdvPix, overlayRGB, out] = advectMeshByPIVTexturePatch(...
    FF, V2Dpix, TV3D, xPiv, yPiv, uPiv, vPiv, IV0, opts);

% Visualize results.
figure('Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imshow(out.imCurrent, []);
title('Current image');

nexttile;
imshow(patchImAdv, []);
hold on;
triplot(FF, meshAdvPix(:, 1), meshAdvPix(:, 2), 'y');
title('Advected texture patch');

nexttile;
imshow(overlayRGB, []);
title('RGB overlay (R=current, G=advected, B=next)');
