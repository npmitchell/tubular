function writeTiff5D(im, name_out, bitDepth)
% Write a TIFF file to disk with order XYCZT
% 
% Parameters
% ----------
% im : (nX x nY x nC x nZ x nT) 16-bit array 
% name_out : str
%   full path of output TIFF filename
% bitDepth : int (default=16)
%   bit depth of TIFF imageStack to write
%
% Returns
% -------
% none
%
% Saves to disk
% -------------
% name_out as image with bitDepth of 16 or supplied value
%
%
% NPMitchell 2020, adapted from https://www.mathworks.com/matlabcentral/answers/389765-how-can-i-save-an-image-with-four-channels-or-more-into-an-imagej-compatible-tiff-format#answer_438003

if nargin < 3
    bitDepth = 16 ;
    maxStr = '65535.0' ;
elseif bitDepth == 16 
    maxStr = '65535.0' ;
elseif bitDepth == 8 
    maxStr = '255.0' ;
else
    error(['Code for this bitDepth = ' num2str(bitDepth)])
end

fiji_descr = ['ImageJ=1.52p' newline ...
        'images=' num2str(size(im,3)*...
                          size(im,4)*...
                          size(im,5)) newline... 
        'channels=' num2str(size(im,3)) newline...
        'slices=' num2str(size(im,4)) newline...
        'frames=' num2str(size(im,5)) newline... 
        'hyperstack=true' newline...
        'mode=grayscale' newline...  
        'loop=false' newline...  
        'min=0.0' newline...      
        ['max=' maxStr]];  
    
        t = Tiff(name_out,'w') ;
        tagstruct.ImageLength = size(im,1);
        tagstruct.ImageWidth = size(im,2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = bitDepth;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.Compression = Tiff.Compression.LZW;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        tagstruct.ImageDescription = fiji_descr;
        for frame = 1:size(im,5)
            for slice = 1:size(im,4)
                for channel = 1:size(im,3)
                    % disp(['FSC = ', num2str(frame), '/', num2str(slice), '/', num2str(channel)])
                    t.setTag(tagstruct)
                    t.write(im(:,:,channel,slice,frame));
                    t.writeDirectory(); % saves a new page in the tiff file
                end
            end
        end
        t.close() 