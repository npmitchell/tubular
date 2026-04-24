function [patchImAdv, meshAdvPix, overlayRGB, out] = ...
    advectMeshByPIVTexturePatch(FF, V2Dpix, TV3D, xPiv, yPiv, uPiv, vPiv, IV0, opts)
%ADVECTMESHBYPIVTEXTUREPATCH Advect pullback mesh by PIV and render texture patch.
%   [patchImAdv, meshAdvPix, overlayRGB, out] = ...
%     advectMeshByPIVTexturePatch(FF, V2Dpix, TV3D, xPiv, yPiv, uPiv, vPiv, IV0, opts)
%
% Parameters
% ----------
% FF : #F x 3 int
%   Mesh face connectivity.
% V2Dpix : #V x 2 float
%   Mesh pullback vertex positions in pixel coordinates [x, y].
% TV3D : #V x 3 float
%   Texture-space coordinates in image indexing format [row, col, page].
% xPiv, yPiv : nRows x nCols float
%   Evaluation grid for PIV displacements. Must match uPiv, vPiv size.
% uPiv, vPiv : nRows x nCols float
%   PIV displacement field in pullback pixels.
% IV0 : image/volume
%   Source intensity volume used for texture mapping.
% opts : struct (optional)
%   patchOpts : struct
%       Passed directly to texture_patch_to_image().
%   interpMethod : char (default='linear')
%       Interpolation method for griddedInterpolant.
%   extrapMethod : char (default='nearest')
%       Extrapolation method for griddedInterpolant.
%   imCurrent : 2D numeric array (optional)
%       Current pullback image for RGB QC overlay (red channel).
%   imNext : 2D numeric array (optional)
%       Next pullback image for RGB QC overlay (blue channel).
%   makeOverlay : bool (default=true)
%       If true, return overlayRGB when a reference image is available.
%   clipToBounds : bool (default=false)
%       If true and patchOpts has xLim/yLim, advected vertices are clipped.
%
% Returns
% -------
% patchImAdv : 2D numeric array
%   Texture patch image generated from the advected mesh.
% meshAdvPix : #V x 2 float
%   Advected pullback mesh vertices [x + u, y + v].
% overlayRGB : H x W x 3 numeric array
%   QC overlay (R=current, G=advected, B=next). Empty if unavailable.
% out : struct
%   addx, addy, patchIm0, imCurrent, imNext.
%
% Notes
% -----
% - This function mirrors the workflow used in measurePIV3d and
%   fitPhiOffsetsFromPrevPullback for visual checks of alignment.
% - TV3D is expected in [row, col, page] indexing convention.

if nargin < 9
    opts = struct() ;
end

% Defaults
patchOpts = struct() ;
interpMethod = 'linear' ;
extrapMethod = 'nearest' ;
imCurrent = [] ;
imNext = [] ;
makeOverlay = true ;
clipToBounds = false ;

if isfield(opts, 'patchOpts')
    patchOpts = opts.patchOpts ;
end
if isfield(opts, 'interpMethod')
    interpMethod = opts.interpMethod ;
end
if isfield(opts, 'extrapMethod')
    extrapMethod = opts.extrapMethod ;
end
if isfield(opts, 'imCurrent')
    imCurrent = opts.imCurrent ;
end
if isfield(opts, 'imNext')
    imNext = opts.imNext ;
end
if isfield(opts, 'makeOverlay')
    makeOverlay = opts.makeOverlay ;
end
if isfield(opts, 'clipToBounds')
    clipToBounds = opts.clipToBounds ;
end

% Validate core array sizes
assert(size(V2Dpix, 2) == 2, 'V2Dpix must be #V x 2')
assert(size(TV3D, 2) == 3, 'TV3D must be #V x 3')
assert(all(size(xPiv) == size(yPiv)), 'xPiv and yPiv must be same size')
assert(all(size(xPiv) == size(uPiv)), 'xPiv and uPiv must be same size')
assert(all(size(xPiv) == size(vPiv)), 'xPiv and vPiv must be same size')

% Interpolate PIV vectors onto mesh pullback vertices
uInterp = griddedInterpolant(xPiv', yPiv', uPiv', interpMethod, extrapMethod) ;
vInterp = griddedInterpolant(xPiv', yPiv', vPiv', interpMethod, extrapMethod) ;

addx = uInterp(V2Dpix(:, 1), V2Dpix(:, 2)) ;
addy = vInterp(V2Dpix(:, 1), V2Dpix(:, 2)) ;
meshAdvPix = [V2Dpix(:, 1) + addx, V2Dpix(:, 2) + addy] ;

% Optional clipping into requested image bounds
if clipToBounds && isfield(patchOpts, 'xLim') && isfield(patchOpts, 'yLim')
    meshAdvPix(:, 1) = min(max(meshAdvPix(:, 1), patchOpts.xLim(1)), patchOpts.xLim(2)) ;
    meshAdvPix(:, 2) = min(max(meshAdvPix(:, 2), patchOpts.yLim(1)), patchOpts.yLim(2)) ;
end

% Render advected pullback image
patchImAdv = texture_patch_to_image(FF, meshAdvPix, FF, TV3D, IV0, patchOpts) ;

% Build optional RGB overlay
overlayRGB = [] ;
patchIm0 = [] ;
if makeOverlay
    if isempty(imCurrent)
        patchIm0 = texture_patch_to_image(FF, V2Dpix, FF, TV3D, IV0, patchOpts) ;
        imCurrent = patchIm0 ;
    end

    if ~isempty(imCurrent)
        % Use current image size as overlay target
        targetSize = size(imCurrent) ;
        if ndims(imCurrent) == 3
            targetSize = targetSize(1:2) ;
        end

        patchImAdv2 = ensureSize2D(patchImAdv, targetSize) ;
        imCurrent2 = ensureSize2D(imCurrent, targetSize) ;

        if isempty(imNext)
            imNext2 = zeros(targetSize, 'like', imCurrent2) ;
        else
            imNext2 = ensureSize2D(imNext, targetSize) ;
        end

        advCast = castLikeScaled(patchImAdv2, imCurrent2) ;
        nextCast = castLikeScaled(imNext2, imCurrent2) ;

        overlayRGB = cat(3, imCurrent2, advCast, nextCast) ;
    end
end

out = struct() ;
out.addx = addx ;
out.addy = addy ;
out.patchIm0 = patchIm0 ;
out.imCurrent = imCurrent ;
out.imNext = imNext ;
end


function im2 = ensureSize2D(im, targetSize)
% Ensure image is 2D and matches desired size.
if ndims(im) == 3
    im = squeeze(mean(im, 3)) ;
end

if ~isequal(size(im), targetSize)
    im2 = imresize(im, targetSize) ;
else
    im2 = im ;
end
end


function outim = castLikeScaled(inim, refim)
% Cast inim to class(refim), preserving dynamic range.
if isfloat(refim)
    outim = cast(inim, 'like', refim) ;
    return
end

imax = max(inim(:)) ;
if imax <= 0
    scaled = zeros(size(inim)) ;
else
    scaled = inim ./ imax ;
end

if isa(refim, 'uint8')
    outim = uint8(255 * scaled) ;
elseif isa(refim, 'uint16')
    outim = uint16((2^16 - 1) * scaled) ;
else
    outim = cast(inim, 'like', refim) ;
end
end