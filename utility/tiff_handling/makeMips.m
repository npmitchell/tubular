function makeMips(timePoints, dir16bit, fileName, mipDir, Options)
% MAKEMIPS(timePoints, dir16bit, fileName, mipDir, Options)
% read 16 bit data images and output the mips 
%
% Parameters
% ----------
% Options : struct with fields
%   scale : float or int (default=timepoint-by-timepoint maximum)
%       intensity value of the maximum brightness of the TIFFs. If no scale
%       is supplied, code uses the max for each timepoint as the saturated
%       intensity value
%   overwrite_mips : bool (default=false)
%       overwrite the MIPS on file with newly computed
%   channels : #channels x 1 int array (default=[1])
%       the channels for which to make mips, ex [1,2] or [0, 1] if 2color
%   allowMissingFiles : bool
%       skip MIP file if 32 bit version is missing
%
% Outputs
% -------
%
% Example usage
% -------------
% makeMips(0:150, './', 'TP%d_Ch0_Ill0_Ang0,60,120,180,240,300.tif', './mips/')
% 
% NPMitchell 2019, based on original script by SJS
%
% NOTE:
% view1 = along third dimension, near half
% view2 = along third dimension, far half
% view11 = along first dimension, near half
% view12 = along first dimension, far half
% view21 = along second dimension, near half
% view22 = along second dimension, far half
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

overwrite_mips = true ;
use_scale = false ; 
channels = [1] ;
allowMissingFiles = true ;
if nargin < 5
    Options = struct() ;
end
if isfield(Options, 'overwrite_mips')
    overwrite_mips = Options.overwrite_mips ;
end
if isfield(Options, 'scale')
    if Options.scale > 0 
        use_scale = true ;
        scale = Options.scale ;
    end
end
if isfield(Options, 'channels')
    channels = Options.channels ;
end
if isfield(Options, 'allowMissingFiles')
    allowMissingFiles = Options.allowMissingFiles ;
end

msgLevel = 1;
setpref('ImSAnE', 'msgLevel', msgLevel);

%% Make the subdirectories for the mips if not already existing
for channel = channels 
    mipdirs = {mipDir, fullfile(mipDir, 'view1/'), ...
        fullfile(mipDir, 'view2/'), ...
        fullfile(mipDir, 'view11/'), ...
        fullfile(mipDir, 'view12/'),...
        fullfile(mipDir, 'view21/'),...
        fullfile(mipDir, 'view22/')} ;
    for i = 1:length(mipdirs)
        if ~exist(mipdirs{i},'dir')
            mkdir(mipdirs{i})
        end
    end
    % Define naming scheme for each half-volume MIP
    % along dim 3
    name1  = fullfile('view1', 'mip_1_%03d_c%d.tif');
    name2  = fullfile('view2', 'mip_2_%03d_c%d.tif');
    % along dim 1
    name11 = fullfile('view11', 'mip_11_%03d_c%d.tif');
    name21 = fullfile('view21', 'mip_21_%03d_c%d.tif');
    % along dim 2 
    name12 = fullfile('view12', 'mip_12_%03d_c%d.tif');
    name22 = fullfile('view22', 'mip_22_%03d_c%d.tif');

    %% Cycle through timepoints to make mips
    for time = timePoints
        disp(['Considering time=' num2str(time)])
        
        if count(fileName, '%') > 1
            fullFileName = fullfile(dir16bit, sprintf(fileName, time, channel)) ;
        else
            if length(channels) > 1
                error('Multiple channels specified but filename has only one variable (time)')
            end
            fullFileName = fullfile(dir16bit, sprintf(fileName, time)) ;
        end
        m1fn = fullfile(mipDir, sprintf(name1,  time, channel)) ;
        m2fn = fullfile(mipDir, sprintf(name2,  time, channel)) ;
        m11fn = fullfile(mipDir, sprintf(name11, time, channel)) ;
        m21fn = fullfile(mipDir, sprintf(name21, time, channel)) ;
        m12fn = fullfile(mipDir, sprintf(name12, time, channel)) ;
        m22fn = fullfile(mipDir, sprintf(name22, time, channel)) ;
        mexist = exist(m1fn, 'file') && exist(m2fn, 'file') && ...
            exist(m11fn, 'file') && exist(m21fn, 'file') && ...
            exist(m12fn, 'file') && exist(m22fn, 'file') ;

        if mexist
            disp(['MIPs already exist on disk for t=' num2str(time)])
            if exist(fullFileName, 'file') && overwrite_mips
                disp('Overwriting...')
            elseif ~exist(fullFileName, 'file') && overwrite_mips
                disp('cannot overwrite since tif does not exist')
            end
        end

        % Only consider this timepoint if mips don't exist, or overwrite
        if exist(fullFileName, 'file') && (~mexist || overwrite_mips)
            disp([ 'reading ' fullFileName]) 
            % data = readSingleTiff(fullFileName);
            data = bfopen(fullFileName) ;
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
            disp(['creating mips for timepoint=' num2str(time)])
            % take max intensity projections of half-space volumes
            mip_1 = max(im2(:,:,1:round(imSize(end)/2)),[],3);
            mip_2 = max(im2(:,:,round(imSize(end)/2):end),[],3);
            mip_11 = squeeze(max(im2(1:round(imSize(1)/2),:,:),[],1));
            mip_21 = squeeze(max(im2(round(imSize(1)/2):end,:,:),[],1));
            mip_12 = squeeze(max(im2(:,1:round(imSize(2)/2),:),[],2));
            mip_22 = squeeze(max(im2(:,round(imSize(2)/2):end,:),[],2));
            disp(['Saving MIP to: ' m1fn])
            imwrite(mip_1, m1fn,'tiff','Compression','none');
            imwrite(mip_2, m2fn,'tiff','Compression','none');
            imwrite(mip_11,m11fn,'tiff','Compression','none');
            imwrite(mip_21,m21fn,'tiff','Compression','none');
            imwrite(mip_12,m12fn,'tiff','Compression','none');
            imwrite(mip_22,m22fn,'tiff','Compression','none');

            % declare done
            clearvars im2
            disp(['finished mips for ' fullFileName])
        else
            if ~exist(fullFileName, 'file') && allowMissingFiles
                disp(['WARNING: file does not exist, skipping: ', fullFileName])
            elseif mexist
                disp(['MIPs skipped for ' fullFileName, ' -- skipping'])
            elseif ~exist(fullFileName, 'file')
                error(['File does not exist: ', fullFileName])
            else 
                error('Somehow the file exists and mips do not, but failed to execute. Investigate here.')
            end
        end
    end
    disp(['Done with channel ' num2str(channel)])
end
disp('done')
