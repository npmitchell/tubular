function extendImages(directory, direc_e, fnbase, timePoints, options)
% EXTENDIMAGES(directory, direc_e, fileNameBase) Repeat an image above and
% below, and equilize histograms in ntiles in each dimension
%
% directory : str
%   path to the existing images
% direc_e : str
%   path to the place where extended images are to be saved
% fnbase : str
%   The file name of the images to load and save
% options : struct with fields (default is no histeq, overwrite==true)
%   outFnBase : str, default=fnbase
%       the output filename base for each extended image
%   histeq : bool
%       whether to equilize the LUT in tiles across the image
%   a_fixed : float
%       The aspect ratio of the pullback image: Lx / Ly
%   ntiles : int 
%       The number of bins in each dimension for histogram equilization for a
%       square original image. That is, the extended image will have (a_fixed *
%       ntiles, 2 * ntiles) bins in (x,y).
%   overwrite : bool
%       overwrite the existing extended image on disk
%
% NPMitchell 2019 

if nargin > 3
    if isfield(options, 'histeq')
        histeq = options.histeq ;
    else
        histeq = false ;
    end
    if isfield(options, 'overwrite')
        overwrite = options.overwrite ;
    else
        overwrite = false ;
    end
    if isfield(options, 'outFnBase')
        outFnBase = options.outFnBase  ;
    else
        outFnBase = fnbase ;
    end
else
    histeq = false ;
    outFnBase = fnbase ;
end

% disp(['Searching for files: ' fnbase])
% fns = dir(fullfile(directory, fnbase)) ;
% % Get original image size
% im = imread(fullfile(fns(1).folder, fns(1).name)) ;
% halfsize = round(0.5 * size(im, 1)) ;
% osize = size(im) ;

% Extend each timepoint's image
for tp=timePoints
    outfn = fullfile(direc_e, sprintf(outFnBase, tp)) ;
    if ~exist(outfn, 'file') || overwrite
        infn = fullfile(directory, sprintf(fnbase, tp)) ;
        
        % Declare if we are overwriting the file
        if exist(outfn, 'file') 
            disp(['Overwriting ' sprintf(outFnBase, tp) ': ' outfn])
        else
            disp(['Generating ' sprintf(outFnBase, tp) ': ' outfn])
        end
        disp(['Reading ' sprintf(fnbase, tp) ': ' infn])
        % fileName = split(fns(i).name, '.tif') ;
        % fileName = fileName{1} ;
        im = imread(infn) ;
        halfsize = round(0.5 * size(im, 1)) ;

        % im2 is as follows:
        % [ im(end-halfsize) ]
        % [     ...          ]
        % [    im(end)       ]
        % [     im(1)        ]
        % [     ...          ]
        % [    im(end)       ]
        % [     im(1)        ]
        % [     ...          ]
        % [  im(halfsize)    ]
        
        is2d = length(size(im)) == 2 || ...
            (length(size(im))==3 && size(im, 3)==1) ;
        
        if is2d
            im2 = uint8(zeros(size(im, 1) + 2 * halfsize, size(im, 2))) ;
            im2(1:halfsize, :) = im(end-halfsize + 1:end, :);
            im2(halfsize + 1:halfsize + size(im, 1), :) = im ;
            im2(halfsize + size(im, 1) + 1:end, :) = im(1:halfsize, :);
       
            % Histogram equilize (LUT tiled into bits)
            if histeq
                im2 = adapthisteq(im2, 'NumTiles', ...
                    [round(options.a_fixed * options.ntiles), round(2 * options.ntiles)]) ;
            end
        elseif length(size(im))==3
            im2 = uint8(zeros(size(im, 1) + 2 * halfsize, size(im, 2), 3)) ;
            im2(1:halfsize, :, :) = im(end-halfsize + 1:end, :, :);
            im2(halfsize + 1:halfsize + size(im, 1), :, :) = im ;
            im2(halfsize + size(im, 1) + 1:end, :, :) = im(1:halfsize, :, :);
       
            % Histogram equilize (LUT tiled into bits)
            if histeq
                error('check that adaptive histogram works on RGB images')
                im2 = adapthisteq(im2, 'NumTiles', ...
                    [round(options.a_fixed * options.ntiles), round(2 * options.ntiles)]) ;
            end
        else
            error('Could not identify image as RGB or grayscale. Is it 3d?')
        end
        
        imwrite( im2, outfn, 'TIFF' );
    else
        disp(['already exists: ' outfn])
    end

end