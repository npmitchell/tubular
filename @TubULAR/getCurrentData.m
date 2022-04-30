function IV = getCurrentData(tubi, adjustIV, varargin)
    % IV = getCurrentData(tubi, adjustIV) 
    % Load/return volumetric intensity data for current timepoint
    % Note: axis order permutation is applied here upon loading and
    % assignment to self (current tubi instance). 
    % Note that the steps executed by this function depends on whether
    % xp is an ImSAnE experiment class or a simple struct belonging to the
    % TubULAR class instance.
    % 
    % Parameters
    % ----------
    % tubi : tubi class instance (self)
    % adjustIV : bool (default=true)
    %   apply the intensity adjustment stored in
    %   tubi.data.adjustlow/high
    %
    % Returns
    % -------
    % IV : #channels x 1 cell array of X*Y*Z intensities arrays
    %   volumetric intensity data of current timepoint
    %   Note this is also assigned to self as tubi.currentData.IV, so
    %   avoiding a duplicate copy in RAM is useful for large data
    %   by not returning an argument and instead calling
    %   tubi.currentData.IV when needed.
    
    
    if nargin < 2
        adjustIV = true ;
    end
    if isempty(tubi.currentTime)
        error('No currentTime set. Use tubi.setTime(timePoint)')
    end
    
    
    if isempty(tubi.currentData.IV)
        % Load 3D data for coloring mesh pullback
        if isa(tubi.xp, 'struct')
            disp('Treating xp as a simple metadata struct')

            % use bioformats if posssible
            hasbf = exist('bioformats_package.jar','file') ||...
                    exist('loci_tools.jar','file');
            if hasbf
                disp('loadStackBioformats() executing...')
                data = loadStackBioformats(tubi, varargin{:});
            else
                data = loadStackNoBioformats(tubi, varargin{:}) ;
            end
            IV = rescaleToUnitAspect(data, tubi.xp.fileMeta.stackResolution);
        else
            disp('Treating xp as an Experiment class instance')
            tubi.xp.loadTime(tubi.currentTime);
            tubi.xp.rescaleStackToUnitAspect();
            IV = tubi.xp.stack.image.apply() ;
        end

        if adjustIV
            adjustlow = tubi.data.adjustlow ;
            adjusthigh = tubi.data.adjusthigh ;
            if (any(adjustlow) ~= 0 || any(adjusthigh) ~=100) && ~any(isnan(adjustlow)) ...
                    && ~any(isnan(adjusthigh))
                IVtmp = tubi.adjustIV(IV, adjustlow, adjusthigh) ;
            else
                IVtmp = IV ;
            end
            if ~all(tubi.data.axisOrder == [1 2 3])
                for ch = 1:length(IVtmp)
                    tubi.currentData.IV{ch} = permute(IVtmp{ch}, tubi.data.axisOrder) ;
                end
            else
                tubi.currentData.IV = IVtmp ;
            end
        else
            if ~all(tubi.data.axisOrder == [1 2 3])
                for ch = 1:length(IV)
                    tubi.currentData.IV{ch} = permute(IV{ch}, tubi.data.axisOrder) ;
                end
            else                   
                tubi.currentData.IV = IV  ;
            end
        end
    end

    % create copy for output to function
    if nargout > 0
        IV = tubi.currentData.IV ;
    end
end

function stacks = rescaleToUnitAspect(data, resolution)
    % RESCALETOUNITASPECT Rescale the stack to unit aspect ratio
    %
    % newStack = rescaleToUnitAspect() 

    if all(resolution) == resolution(1)
        stacks = data;
        return;
    end

    % Right now we just handle the case of equal resolution in 2 dimensions
    % and unequal in the third dimension
    if length(unique(resolution)) > 2
        error('Handle the case of different resolution in all three dimensions here')
    else
        aspect = max(resolution)/min(resolution);
    end
    
    % prepare the perumutation of stacks, such that the low
    % resolution axis becomes the z axis;
    for channel = 1:size(data, 4)
        stacks{channel} = squeeze(data(:,:,:,channel)) ;
    end
    resolution = resolution([2,1,3]);
    [~,ii]     = sort(resolution);

    % Resample the 
    for i = 1:numel(stacks)

        fprintf(['Channel ' num2str(i) '\n']);

        curr = stacks{i};
        % permute;
        curr = permute(curr,ii);
        newnslices = round(size(curr,3)*aspect);
        scaled = zeros([size(curr,1) size(curr,2) newnslices], class(curr));
        for j=1:size(curr,1)
            % debugMsg(1, '.');
            % if rem(j,80) == 0
            %     debugMsg(1, '\n');
            % end
            scaled(j,:,:) = imresize(squeeze(curr(j,:,:)),...
                                        [size(curr,2) newnslices]);
        end
        % debugMsg(1,'\n');
        
        % permute back to original axis order
        scaled = ipermute(scaled,ii);
        stacks{i} = scaled;
    end
    
end

function data = loadStackBioformats(tubi, varargin)
    % Only used if tubi.xp is not ImSAnE Experiment class instance
    % Load stack from disc into project.
    %
    % loadStack()
    % loadStack(justMeta)
    %
    % Based on the metadata, load the stack at the current time
    % point using bioformats library. 
    % 
    % see also getCurrentData

    if nargin==2
        justMeta = varargin{1};
    else
        justMeta = 0;
    end

    fileName = sprintf(tubi.xp.fileMeta.filenameFormat, tubi.currentTime);
    fullFileName = fullfile(tubi.xp.fileMeta.dataDir, fileName);

    % load the Bio-Formats library into the MATLAB environment
    autoloadBioFormats = 1;
    status = bfCheckJavaPath(autoloadBioFormats);
    assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
        'to the static Java path or add it to the Matlab path.']);

    fprintf(['Using bioformats version ' char(loci.formats.FormatTools.VERSION) '\n']);

    % Get the channel filler
    fprintf(['loading file : ' fullFileName '\n']);
    r = bfGetReader(fullFileName);
    r.setSeries(0);

    % corrupted metadata can give the wrong stack resolution and
    % cause problems, we should set resolution by hand

    if tubi.xp.fileMeta.swapZT == 0
        stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeZ(), r.getSizeT()];
    else
        stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeT(), r.getSizeZ()];
    end
    fprintf(['stack size (xyzt) ' num2str(stackSize) '\n']);

    xSize = stackSize(1);
    ySize = stackSize(2);
    zSize = stackSize(3);

    % number of channels
    nChannels = r.getSizeC();

    nTimePts = stackSize(4);

    % update stack size in metadata 
    tubi.xp.fileMeta.stackSize = [xSize, ySize, zSize];
    assert(tubi.xp.fileMeta.nChannels == nChannels) ;

    % if only reading stack size, stop here-----------------------
    if justMeta
        return; 
    end

    nChannelsUsed = numel(tubi.xp.expMeta.channelsUsed);

    if isempty(nChannelsUsed)
        error('nChannelsUsed is empty');
    else
        disp(['Reading ' num2str(nChannelsUsed) ' channels \n'])
    end
    if any(tubi.xp.expMeta.channelsUsed > nChannels)
        error('channelsUsed specifies channel larger than number of channels');
    end

    fprintf(['Loading stack for time: ' num2str(tubi.currentTime) '\n']);

    % this is a stupid workaround
    % basically, even though you save each time point in a separate
    % file with the bioformats exporter, the metadata tells the
    % reader that it is part of a video and it will try to read all
    % the files as a single series
    % this is just a warning to the user, in the loop below is a
    % conditional to prevent it from reading all times
    if nTimePts > 1 
        fprintf('More than one time point in metadata!\n');
    end

    % read the data
    ticID = tic;
    data = zeros([ySize xSize zSize nChannelsUsed], 'uint16');

    for i = 1:r.getImageCount()

        ZCTidx = r.getZCTCoords(i-1) + 1;

        % in the fused embryo data coming out of the python script,
        % Z and T are swaped. In general this isn't the case, thus
        % introduce a file metaField swapZT
        if tubi.xp.fileMeta.swapZT == 0
            zidx = ZCTidx(1);
            tidx = ZCTidx(3);
        else 
            zidx = ZCTidx(3);
            tidx = ZCTidx(1);
        end
        cidx = ZCTidx(2);

        % see above: if there is only one timepoint all the planes
        % should be read, if there are multiple timepoints, only
        % the correct time should be read
        if nTimePts == 1 || (nTimePts > 1 && this.currentTime == tidx-1)
            if rem(i,80) == 0
                disp('...');
            end

            dataCidx = find(tubi.xp.expMeta.channelsUsed == cidx);
            if ~isempty(dataCidx)
                % bfGetPlane gives ALL Channels, and we stuff all
                % cidx into data at once.
                data(:,:, zidx, dataCidx) = bfGetPlane(r, i);
            else
                fprintf('skipping channel and z plane')
            end
        end
    end
    fprintf('\n');
    dt = toc(ticID);
    fprintf(['dt = ' num2str(dt) '\n']);

    % announce min and max
    fprintf('tubular.getCurrentData: loadStack() finished loading data volume')
end

function data = loadStackNoBioformats(tubi, varargin)
    % Only used if tubi.xp is not ImSAnE Experiment class instance
    %
    %
    fileName = sprintf(tubi.xp.fileMeta.filenameFormat, tubi.currentTime);
    fullFileName = fullfile(tubi.xp.fileMeta.dataDir, fileName);
    tmp = imfinfo(fullFileName);
    nImages = numel(tmp);
    tmp = tmp(1);
    isRGB = strcmp(tmp.ColorType,'truecolor');

    if nargin==2

        xSize = tmp.Width;
        ySize = tmp.Height;
        if isRGB
            zSize = nImages;
            assert(this.fileMeta.nChannels==3,...
                'your data is RGB, fileMeta.nChannels should equal 3');
        else
            zSize = nImages / this.fileMeta.nChannels;
        end
        this.fileMeta.stackSize = [xSize ySize zSize];

        return
    end

    nChannels = tubi.xp.fileMeta.nChannels;

    if strcmp(tmp.ColorType,'grayscale')
        assert(nImages == zSize*nChannels,...
             'fileMeta.nChannels is not consistent with data');
    elseif isRGB
        assert(nImages == zSize,...
             'fileMeta.nChannels is not consistent with data');
    end

    % number of channels used
    nChannelsUsed = numel(tubi.xp.expMeta.channelsUsed);

    % read the data
    ticID = tic;
    for i = 1:nImages
        % load image plane
        im = imread(fullFileName, i);

        % if first image, preallocate data array
        if i == 1
            ySize = size(im, 1) ;
            xSize = size(im, 2) ;
            zSize = nImages ;
            data = zeros([ySize xSize zSize nChannelsUsed], 'uint16');
        end

        if isRGB
            data(:,:, i, :) = im(:,:,tubi.xp.expMeta.channelsUsed);	
        else
            zidx = 1 + floor((i-1) / nChannels);
            cidx = 1 + mod(i-1, nChannels);
            if sum(tubi.xp.expMeta.channelsUsed == cidx)
                data(:,:, zidx, tubi.xp.expMeta.channelsUsed == cidx) = im;
            end
        end

        % progress indicator
        % debugMsg(1,'.');
        % if rem(i,80) == 0
        %     debugMsg(1,'\n');
        % end
    end
end