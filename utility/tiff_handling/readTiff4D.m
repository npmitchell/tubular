function im = readTiff4D(fn, ncolors, dimOrder)
% wrapper for fast loadtiff, with handling of multicolors. 
% This is a good one -NPM.
%
% NPMitchell 2020

if nargin < 3
    dimOrder = 'XYCZ' ;
end

if nargin < 2 
    ncolors = 1;
end

im = loadtiff(fn) ;

if ncolors > 1
    % Get the color axis to split that dimension
    cdim = find(lower(dimOrder) == 'c') ;
    splitID = size(im, cdim) / ncolors ;
    if mod(splitID, 1) ~= 0
        error('Color dimension does not split evenly between channels')
    end
    if strcmpi(dimOrder, 'xycz')
        for qq = 1:ncolors
            im2{qq} = im(:, :, qq:ncolors:end) ;
        end
        im = im2 ;
    elseif strcmpi(dimOrder, 'xyzc')
        for qq = 1:ncolors
            im2{qq} = im(:, :, max((qq-1)*splitID, 1):(qq*splitID-1)) ;
        end
        im = im2 ;
    else
        error('handle this dimension order here')
    end
    
    % Check that all channels are the same size
end