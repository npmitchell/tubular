function outim = mergeChannelsRGB(im1, im2, im1color, im2color, bitdepth)
% Combine multiple images as RGB colored layers and add them together.
%
% 
%
% Parameters
% ----------
% im1fn : str or 2d grayscale/BW image
%   full path to the first image or the image itself
% im2fn : str or 2d grayscale/BW image
%   full path to the second image or the image itself
%   todo: allow multiple additional image
% im1color : length(3) float (values between 0-1)
%   RGB values of the first image in merge
% im2color : length(3) float (values between 0-1) 
%   RGB values of the second image in merge
%
%
% NPMitchell 2021

if ischar(im1)
    im1 = imread(im1) ;
end

if ischar(im2)
    im2 = imread(im2) ;
end

outim = zeros(size(cat(3, im1, im1, im1)), class(im1)) ;
outim(:, :, 1) = im1 * im1color(1) ;
outim(:, :, 2) = im1 * im1color(2) ;
outim(:, :, 3) = im1 * im1color(3) ;

outim(:, :, 1) = outim(:, :, 1) + im2 * im2color(1) ;
outim(:, :, 2) = outim(:, :, 2) + im2 * im2color(2) ;
outim(:, :, 3) = outim(:, :, 3) + im2 * im2color(3) ;

% Clip at max intensity for this bitdepth
maxI = 2^bitdepth -1 ;
outim(outim > maxI) = maxI ;

