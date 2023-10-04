function convert32to16bit(timePoints, scales, dir32bit, dir16bit, ...
    file32Name, file16Name, channels)
% CONVERT32TO16BIT(timePoints, scale, dataDir, file32Name, file16Name)
% Read 32 bit file from deconvolution, save as 16 bit and save maximum
% intensity projections as we go. 
%
% Parameters
% ----------
% timePoints : 
% scales : #channels x 1 float array
%   intensity value of the maximum brightness of the TIFFs
% dir32bit : str
%   path to data where 32-bit TIFF files live
% dir16bit : path to data to store metadata for 16bit data
% file32Name : char
% file16Name : char
% channels : #channels x 1 int array (default=[1])
%   the channels for which to make mips, ex [1,2] or [0, 1] if 2color
%
% Outputs
% -------
% saves 16 bit files to fullfile(dataDir, sprintf(file16Name, t)) for each
% time t in timePoints.
%
% NPMitchell 2019-2020, based on a script by SJS 

% Input handling
if nargin < 7
    channels = [1] ;
end

%% Ensure output mips directory
mipsDir = fullfile(dir32bit, 'mips32bit') ;
dirs2make = {dir16bit, dir32bit, mipsDir, ...
    fullfile(mipsDir, 'view1'), ...
    fullfile(mipsDir, 'view2'), ...
    fullfile(mipsDir, 'view11'), ...
    fullfile(mipsDir, 'view21'), ...
    fullfile(mipsDir, 'view12'), ...
    fullfile(mipsDir, 'view22')} ;
for ii = 1:length(dirs2make) 
    dir2make = dirs2make{ii} ;
    if ~exist(dir2make, 'dir')
        mkdir(dir2make);
    end
end
logfn = fullfile(dir16bit, 'convert32to16_log.mat') ;

msgLevel = 1;
setpref('ImSAnE', 'msgLevel', msgLevel);

%% Cycle over timePoints for conversion
for time = timePoints
    for cqq = 1:length(channels)
        channel = channels(cqq) ;
        scale = scales(cqq) ;
        if sum(file32Name == '%') == 1
            fileName = fullfile(dir32bit, sprintf(file32Name, time)) ;
        else
            fileName = fullfile(dir32bit, sprintf(file32Name, time, channel)) ;
        end
        disp(['32 bit file ' fileName])

        % 16bit name
        if sum(file16Name == '%') == 1
            out16name = fullfile(dir16bit, sprintf(file16Name, time)) ;
        else
            out16name = fullfile(dir16bit, sprintf(file16Name, time, channel)) ;
        end
        disp(['Seeking ' out16name])
        if exist(out16name, 'file')
            disp(' --> 16bit file exists, skipping')
        else
            try
                % Read the 32 bit data and convert
                data = readSingleTiff(fileName);
                im2 = mat2gray(data, [0 scale]);
                im2 = uint16(2^16*im2);
                imSize = size(im2);

                % Could stabilize based on near half of the image 
                mip_1 = max(im2(:,:,1:round(imSize(end)/2)),[],3);
                mip_2 = max(im2(:,:,round(imSize(end)/2):end),[],3);
                mip_11 = squeeze(max(im2(1:round(imSize(1)/2),:,:),[],1));
                mip_21 = squeeze(max(im2(round(imSize(1)/2):end,:,:),[],1));
                mip_12 = squeeze(max(im2(:,1:round(imSize(2)/2),:),[],2));
                mip_22 = squeeze(max(im2(:,round(imSize(2)/2):end,:),[],2));

                m1fn = fullfile(mipsDir, 'view1', sprintf('mip_1_%03d_c%d.tif',  time, channel)) ;
                m2fn = fullfile(mipsDir, 'view2', sprintf('mip_2_%03d_c%d.tif',  time, channel)) ;
                m11fn = fullfile(mipsDir, 'view11', sprintf('mip_11_%03d_c%d.tif',time, channel)) ;
                m21fn = fullfile(mipsDir, 'view21', sprintf('mip_21_%03d_c%d.tif',time, channel)) ;
                m12fn = fullfile(mipsDir, 'view12', sprintf('mip_12_%03d_c%d.tif',time, channel)) ;
                m22fn = fullfile(mipsDir, 'view22', sprintf('mip_22_%03d_c%d.tif',time, channel)) ;

                imwrite(mip_1, m1fn, 'tiff', 'Compression', 'none');
                imwrite(mip_2, m2fn, 'tiff', 'Compression', 'none');
                imwrite(mip_11, m11fn, 'tiff', 'Compression', 'none');
                imwrite(mip_21, m21fn, 'tiff', 'Compression', 'none');
                imwrite(mip_12, m12fn, 'tiff', 'Compression', 'none');
                imwrite(mip_22, m22fn, 'tiff', 'Compression', 'none');

                for z = 1:imSize(3)
                    % write the first page as overwrite mode, then append
                    if z == 1
                        imwrite(im2(:,:,z), out16name, 'tiff', 'Compression','none');
                    else
                        imwrite(im2(:,:,z), out16name, 'tiff', 'Compression','none','WriteMode','append');    
                    end
                end

                % append the metadata used for this timepoint to scale.mat
                tconv = datestr(datetime, 'yyyy-mm-dd(HH:MM:SS)') ;
                tstamp = [sprintf('%06d: ', time) tconv ] ;
                fullstamp = [tstamp ' scale=' num2str(scale) '\n'] ;
                if exist(logfn, 'file')
                    load(logfn, 'metadata')
                    metadata = [metadata fullstamp] ;
                    save(logfn, 'metadata')
                else
                    metadata = tstamp ;
                    save(logfn, 'metadata')
                end
                sprintf(fullstamp)
            catch
                disp('WARNING! Could not convert timepoint!')
            end
        end
    end
end
disp('done')