function makeSubStackMips(timePoints, dir16bit, fileName, mipDir, Options)
% MAKEMIPS(timePoints, dir16bit, fileName, mipDir, Options)
% read 16 bit data images and output the mips 
%
% Parameters
% ----------
% Options : struct with fields
%   pages : Nx1 int or empty for default=(midplane+/-width pages)
%       pages of the TIFF over  which to take a MIP
%   width : int, ignored if Options.pages is supplied (default=2)
%       halfwidth of #pages of the TIFF over which to take a MIP, if pages
%       are not directly specified
%   scale : float or int (default=timepoint-by-timepoint maximum)
%       intensity value of the maximum brightness of the TIFFs. If no scale
%       is supplied, code uses the max for each timepoint as the saturated
%       intensity value
%   overwrite_mips : bool (default=false)
%       overwrite the MIPS on file with newly computed
%   channels : #channels x 1 int array (default=[1])
%       the channels for which to make mips, ex [1,2] or [0, 1] if 2color
%   dim : int
%       dimension along which to take the substack and project
%   colors : #channels x 3 RGB float or int array (optional)
%       each row is a color to use in overlay. 
%   adjust_high : float between 0 and 1
%       for multiColor overlay of all channels, what cdf value to use to 
%       clip the intensity 
%
% Outputs
% -------
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

% Default options
dim = 3 ;
overwrite_mips = true ;
use_scale = false ; 
channels = [1] ;
overlayChannels = false ;
overwrite_overlays = false ;
pages = [] ;
width = 2 ;
adjust_high = 0.95 ;

% Unpack options
if isfield(Options, 'overwrite_mips')
    overwrite_mips = Options.overwrite_mips ;
end
if isfield(Options, 'overwrite_overlays')
    overwrite_overlays = Options.overwrite_overlays ;
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
if length(channels) > 1
    if isfield(Options, 'overlayChannels')
        overlayChannels = Options.overlayChannels ;
    else
        overlayChannels = true ;
    end
end
if isfield(Options, 'dim')
    dim = Options.dim ;
elseif isfield(Options, 'dimension')
    dim = Options.dimension ;
end
if isfield(Options, 'pages')
    pages = Options.pages ;
else
    if isfield(Options, 'width')
        width = Options.width ;
    end
end
if isfield(Options, 'adjust_high')
    adjust_high = Options.adjust_high ;
end
if isfield(Options, 'colors')
    colors = Options.colors ;
else
    if length(channels) == 2
        colors = [1, 0.5, 0; 0, 0.5, 1] ;
    elseif length(channels) == 3
        colors = [1, 0, 0; 0, 1, 0; 0, 0, 1] ;
    else
        % pick an array of colors
        colors = [0,    0.4470,    0.7410;
            0.8500,    0.3250,    0.0980;
            0.9290,    0.6940,    0.1250;
            0.4940,    0.1840,    0.5560;
            0.4660,    0.6740,    0.1880;
            0.3010,    0.7450,    0.9330;
            0.6350,    0.0780,    0.1840;
            0.2000,    0.2000,    0.2000;
            0.5400,    0.2500,    0.0900;
                 0,    0.6500,    0.5200 ] ;
        colors = colors ./ vecnorm(colors, 2, 2) ;
    end
end

msgLevel = 1;
setpref('ImSAnE', 'msgLevel', msgLevel);

%% Make the subdirectories for the mips if not already existing
for channel = channels 
    
    mipcDir = fullfile(mipDir, 'mips') ;
    if ~exist(mipcDir,'dir')
        mkdir(mipcDir)
    end
    
    % Define naming scheme for each substack MIP
    % along dim sprintf(%d', dim)
    name  = fullfile('mip_dim%d_%03d_c%d.tif');

    %% Cycle through timepoints to make mips
    for time = timePoints
        disp(['Considering time=' num2str(time)])
        
        if count(fileName, '%') > 1
            fullFileName = fullfile(dir16bit, sprintf(fileName, time, channel)) ;
        else
            fullFileName = fullfile(dir16bit, sprintf(fileName, time)) ;
        end
        mfn = fullfile(mipcDir, sprintf(name, dim, time, channel)) ;
        mexist = exist(mfn, 'file') ;

        if mexist  
            disp(['MIP already exists on disk for t=' num2str(time)])
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
            % imSize = size(im2);

            % Creat MIPs (maximum intensity projections)
            disp(['creating mips for timepoint=' num2str(time)])
            
            % Identify substack
            if isempty(pages)
                npages = size(im2, dim) ;
                midplane = round(0.5 * npages) ;
                if isempty(pages)
                    pages = max(midplane-width, 1):min(midplane+width, npages) ;
                end
            end
            
            % take max intensity projections of substack volumes
            if dim == 1
                mip_data = max(im2(pages,:,:),[],dim) ;
            elseif dim == 2
                mip_data = max(im2(:,pages,:),[],dim) ;
            elseif dim == 3
                mip_data = max(im2(:,:,pages),[],dim) ;
            end
            disp(['Saving MIP to: ' mfn])
            imwrite(mip_data, mfn, 'tiff','Compression','none');

            % declare done
            clearvars im2
            disp(['finished mips for ' fullFileName])
        else
            if ~exist(fullFileName, 'file') 
                disp(['WARNING: file does not exist, skipping: ', fullFileName])
            elseif mexist
                disp(['MIPs skipped for ' fullFileName, ' -- skipping'])
            else 
                error('Somehow the file exists and mips do not, but failed to execute. Investigate here.')
            end
        end
    end
    disp(['Done with channel ' num2str(channel)])
end


%% Overaly all channels into one colored image
if overlayChannels
    cname  = fullfile('mip_overlay_dim%d_%03d.png');

    olDir = fullfile(mipDir, 'overlay') ;
    if ~exist(olDir, 'dir')
        mkdir(olDir)
    end
    
    %% Cycle through timepoints to make mips
    for time = timePoints
        disp(['overlay: considering time=' num2str(time)])
        
        cfn = fullfile(olDir, sprintf(cname, dim, time)) ;
        mexist = exist(cfn, 'file') ;
        if ~mexist || overwrite_mips || overwrite_overlays
            for qq = 1:length(channels)
                channel = channels(qq) ;
                mipcDir = fullfile(mipDir, 'mips') ;
                mfn = fullfile(mipcDir, sprintf(name, dim, time, channel)) ;
                imadd = double(imread(mfn)) ;
                imadd = mat2gray(imadd, [0, max(imadd(:))]) ;
                imadd = imadjustn(imadd, [0; adjust_high]) ;
                if qq == 1
                    % Initialize RGB image with first channel
                    im = zeros(size(imadd, 1), size(imadd, 2), 3) ;
                    for rgbi = 1:3
                        im(:,:,rgbi) = imadd*colors(qq,rgbi) ;
                    end
                else
                    % Add this channel to the RGB image
                    imadd = double(imread(mfn)) ;
                    imadd = mat2gray(imadd, [0, max(imadd(:))]) ;
                    imadd = imadjustn(imadd, [0; adjust_high]) ;
                    for rgbi = 1:3
                        im(:,:,rgbi) = im(:,:,rgbi)+imadd*colors(qq,rgbi) ;
                    end
                end
            end   
            
            % Save the image
            im = uint16(2^16 * im);
            im(im(:) > 65535) = 65535 ;
            imwrite(im, cfn) ; % 'png', 'Compression', 'none');
            % imwrite(im, cfn, 'tiff', 'Compression', 'none');
        end
    end
end


disp('done')
