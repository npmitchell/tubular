function mip = tiffMIP(fileName, mipfn, slices, scale, dim, fileType)
% mip = tiffMIP(fileName, mipfn, slices, scale, dim)
%
% Parameters
% ----------
% fileName : str
%   tiff file to load
% mipfn : string or cell of strings  if dim is list
%   filename(s) to save the mip(s)
% slices : list of ints
%   slices of the tiff to take MIP of
% scale : optional numeric value
%   maximum intensity to cap mip at
% dim : int or list of ints
%   dimensions along which to take mip
%
% Example usage
% -------------
% 
% NPMitchell 2021


if nargin < 4
    scale = [] ;
    dim = 3 ;
elseif nargin < 5
    dim = 3 ;
end
if nargin < 6
    if any(contains(mipfn, 'tif'))
        fileType = 'tiff' ;
    elseif any(contains(mipfn, 'png'))
        fileType = 'png' ; 
    else
        fileType = 'tiff'; 
    end
end

if ~isempty(scale)
    use_scale = scale > 0 ;
else
    use_scale = false ;
end

% Load in a tiff, take mip of a chunk, save image 
disp([ 'reading ' fileName]) 
% data = readSingleTiff(fileName);
data = bfopen(fileName) ;
tmp = data{1} ;
dstack = zeros([size(tmp{1}), length(data{1})]) ;
for i = 1:length(data{1})
    dstack(:, :, i) = tmp{i};  
end

% Convert to grayscale at 16 bit depth
if ~use_scale
    scale = max(dstack(:)) ;
end
im2 = mat2gray(dstack, [0 scale]);
im2 = uint16(2^16 * im2);
imSize = size(im2);

% Creat MIPs (maximum intensity projections)
disp(['creating mips'])
% take max intensity projections of half-space volumes
if length(dim) == 1
    
    if isempty(slices)
        mip = max(im2(:,:,:),[],dim);
    else
        if dim == 3
            mip = max(im2(:,:,slices), [], dim);
        elseif dim == 2
            mip = max(im2(:,slices, :), [], dim);
        elseif dim == 1
            mip = max(im2(slices, :, :), [], dim);
        end
    end
    if ~isempty(mipfn)
        disp(['Saving MIP to: ' mipfn])
        imwrite(mip, mipfn, 'tiff','Compression','none');
    else
        disp('WARNING: No mipfn, returning mip')
    end
else
    disp('Saving one mip per dim in dims')
    mip = cell(length(dim), 1) ;
    for qq = 1:length(dim)
        dd = dim(qq) ;
        if isempty(slices)
            mip{qq} = squeeze(max(im2(:,:,:), [], dd)) ; 
        else
            if dd == 3
                if ischar(slices)
                    if contains(slices, 'mid')
                        mid = round(0.5*size(im2, 3)) + [-1,0,1] ;
                        mip{qq} = squeeze(max(im2(:,:,mid), [], dd)) ;
                    end
                else
                    mip{qq} = squeeze(max(im2(:,:,slices), [], dd)) ;
                end
            elseif dd == 2
                
                if ischar(slices)
                    if contains(slices, 'mid')
                        mid = round(0.5*size(im2, 2)) + [-1,0,1] ;
                        mip{qq} = squeeze(max(im2(:,mid,:), [], dd)) ;
                    end
                else
                    mip{qq} = squeeze(max(im2(:,slices, :), [], dd)) ;
                end
            elseif dd == 1
                
                if ischar(slices)
                    if contains(slices, 'mid')
                        mid = round(0.5*size(im2, 1)) + [-1,0,1] ;
                        mip{qq} = squeeze(max(im2(mid,:,:), [], dd)) ;
                    end
                else
                    mip{qq} = squeeze(max(im2(slices, :, :), [], dd)) ;
                end
            end
        end
    end
    
    if iscell(mipfn)
        for qq = 1:length(mipfn)
            disp(['Saving MIP to: ' mipfn{qq}])
            imwrite(mip{qq}, mipfn{qq}, fileType, 'Compression', 'none');
        end
    elseif ischar(mipfn) && length(dim) == 1
        disp('Writing single mip to disk')
        imwrite(mip{1}, mipfn, fileType, 'Compression', 'none')
    else
        disp('WARNING: No mipfns, returning mips')
    end
end

% declare done
clearvars im2
disp(['finished mip for ' fileName])