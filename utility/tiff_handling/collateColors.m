function collateColors(fileNameIn, fileNameOut, timePoints, channels)
% collateColors(fileNameIn, fileNameOut, timePoints, channels)
%
% Parameters
% ----------
% fileNameIn : str
% fileNameOut : str
% timePoints : Nx1 int array
% channels #channels x 1 int array
%
% Saves to disk
% -------------
% single TIFF for each timepoint as XYZC, saved to sprintf(fileNameOut, tp)
%
% NPMitchell 2020

dimensionOrder = 'XYCZT' ;   % must be in {'XYZCT', 'XYZTC', 'XYCTZ', ...
                            %             'XYCZT', 'XYTCZ', 'XYTZC'}

%% Define stackSize, which is the number of frames in the z dimension
disp('defining stackSize (which is number of valid matching files)...')
done = false ;
stackSize = 0 ;
if isempty(timePoints)
    name_ref = sprintf(fileNameIn, channels(1)) ;
else
    name_ref = sprintf(fileNameIn, timePoints(1), channels(1)) ;
end
while ~done
    stackSize = stackSize + 1 ;
    try 
        tmp = imread(name_ref, stackSize);
        typename = class(tmp) ;
    catch
        done = true ;
    end
end    
stackSize = stackSize - 1 ;

%% Build image for each timepoint =========================================
disp('Running through timepoints to build ims...')
if isempty(timePoints)
    tidx_todo = 1 ;
else
    tidx_todoA = 1:10:length(timePoints) ;
    tidx_todoB = setdiff(1:length(timePoints), tidx_todoA) ;
    tidx_todo = [tidx_todoA, tidx_todoB] ;
end
for tid = tidx_todo    
    disp(['considering tidx = ' num2str(tid)])
    if ~isempty(timePoints)
        time = timePoints(tid);
        name_out = sprintf(fileNameOut, time) ;
    else
        time = 1 ;
        name_out = fileNameOut ;
    end
    
    tiff_exists = exist(name_out, 'file') ;
    if ~tiff_exists
        for cid = 1:length(channels)
            channel = channels(cid) ;
            if isempty(timePoints)
                imfn = sprintf(fileNameIn, channel);
            else
                imfn = sprintf(fileNameIn, time, channel);
            end
            if channel == channels(1)
                tmp = imread(imfn, 1) ;
                if strcmpi(dimensionOrder, 'xyzc')
                    im = zeros([size(tmp) stackSize length(channels)], typename) ;
                elseif strcmpi(dimensionOrder, 'xycz')
                    im = zeros([size(tmp) length(channels) stackSize], typename) ;
                end
                clear tmp
            end

            % Stuff the image for output TIFF
            disp(['Stuffing t=', num2str(time), ...
                ': channel=', num2str(channel), ...
                ' in ', dimensionOrder])
            for z = 1:stackSize
                if mod(z, 200) == 0
                    disp(['z = ', num2str(z)])
                end
                if strcmpi(dimensionOrder(1:4), 'xyzc')
                    im(:,:,z,cid)= imread(imfn,z);
                elseif strcmpi(dimensionOrder(1:4), 'xycz')
                    im(:,:,cid,z)= imread(imfn,z);
                end
            end
            disp('done stuffing into image')
        end

        % Write to disk. Note that if uint16 specified upon instantiation, then
        % imwrite will respect the class type.
        if ~tiff_exists || overwrite_tiffs
            disp(['Saving TIFF: t=' num2str(time) ' to ' name_out])
            bfsave(im, name_out, 'dimensionOrder', dimensionOrder)
            
            % Alternative to bfsave
            % fiji_descr = ['ImageJ=1.52p' newline ...
            %             'images=' num2str(size(MultiDimImg,3)*...
            %                               size(MultiDimImg,4)*...
            %                               size(MultiDimImg,5)) newline... 
            %             'channels=' num2str(size(MultiDimImg,3)) newline...
            %  ---> ignore? % 'slices=' num2str(size(MultiDimImg,4)) newline...
            %  ---> ignore? %  'frames=' num2str(size(MultiDimImg,5)) newline... 
            %             'hyperstack=true' newline...
            %             'mode=grayscale' newline...  
            %             'loop=false' newline...  
            %             'min=0.0' newline...      
            %             'max=65535.0'];  % change this to 256 if you use an 8bit image
            % 
            % t = Tiff(name_out, 'w')
            % tagstruct.ImageLength = size(MultiDimImg,1);
            % tagstruct.ImageWidth = size(MultiDimImg,2);
            % tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
            % tagstruct.BitsPerSample = 16;
            % tagstruct.SamplesPerPixel = 1;
            % tagstruct.Compression = Tiff.Compression.LZW;
            % tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            % tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
            % tagstruct.ImageDescription = fiji_descr;
            % % for frame = 1:size(MultiDimImg,5)
            % for channel = 1:size(im, 3)
            %     for slice = 1:size(im,4)
            %         t.setTag(tagstruct)
            %         t.write(im2uint16(im(:,:,slice,channel)));
            %         t.writeDirectory(); % saves a new page in the tiff file
            %     end
            % end
            % t.close() 
            
            disp(['saved image ' name_out ])
        end
    end
end

            