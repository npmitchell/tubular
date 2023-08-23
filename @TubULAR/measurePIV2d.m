function measurePIV2d(tubi, options)
% measurePIV2d(tubi, options)
%   Measure 2d PIV on pullbacks in some coordinate system (coordSys),
%   either using PIVLab (default is use_PIVLab=true) or using a simple
%   phasecorrelation with a single pass at a single fixed window size.
% 
%   Note that this code does NOT use the time interval tubi.timeInterval
%   between frames, since the 2D PIV connects the optical flow of adjacent 
%   frames and does not measure velocity in any physical sense.
%
% Parameters
% ----------
% tubi : TubULAR class instance
% options : struct with optional fields
%   use_PIVLab : bool (default=True)
%   edgeLength : size of interrogation window in pixels
%   histequilize : bool (default=True)
%       perform histogram equilization before piv calculation
%   timePoints : Nx1 numeric (default is tubi.xp.fileMeta.timePoints)
%       the timepoints for which to compute PIV
%   coordSys : string specifier (default is tubi.piv.imCoords)
%       the image coordinates in which PIV is to be performed
%   
% Returns
% -------
% none
% 
% Saves to disk
% -------------
% tubi.fileName.pivRaw.raw: 'x', 'y', 'u', 'v', 'u_filtered', 'v_filtered'
%
% NPMitchell 2022

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
else
    overwrite = false ;    
end
if isfield(options, 'usePIVLab')
    use_PIVLab = options.usePIVLab ;
else
    use_PIVLab = true ;    
end
if isfield(options, 'edgeLength')
    edgeLength = options.edgeLength ;
else
    edgeLength = 32 ;
end
if isfield(options, 'histequilize')
    histequilize = options.histequilize ;
else
    histequilize = true ;
end
if isfield(options, 'timePoints')
    timePoints = options.timePoints ;
else
    timePoints = tubi.xp.fileMeta.timePoints ;
end
if isfield(options, 'coordSys')
    coordSys = options.coordSys ;
else
    coordSys = tubi.piv.imCoords ;
end

if ~exist(tubi.fileName.pivRaw.raw, 'file') || overwrite 
    % Compute PIV for each timepoint requested
    for tidx = 1:(length(timePoints)-1)

        % Figure out which images to compare (First ask: which coordinate 
        % system are the images painted into? Perhaps spsme?)
        if strcmpi(strrep(coordSys, '_', ''), 'spsme')
            im1 = imread(sprintfm(tubi.fullFileBase.im_sp_sme, timePoints(tidx))) ;
            im2 = imread(sprintfm(tubi.fullFileBase.im_sp_sme, timePoints(tidx+1))) ;
        elseif strcmpi(strrep(coordSys, '_', ''), 'spsm')
            im1 = imread(sprintfm(tubi.fullFileBase.im_sp_sm, timePoints(tidx))) ;
            im2 = imread(sprintfm(tubi.fullFileBase.im_sp_sm, timePoints(tidx+1))) ;
        elseif strcmpi(strrep(coordSys, '_', ''), 'spe')
            im1 = imread(sprintfm(tubi.fullFileBase.im_spe, timePoints(tidx))) ;
            im2 = imread(sprintfm(tubi.fullFileBase.im_spe, timePoints(tidx+1))) ;
        elseif strcmpi(strrep(coordSys, '_', ''), 'uv')
            im1 = imread(sprintfm(tubi.fullFileBase.im_uv, timePoints(tidx))) ;
            im2 = imread(sprintfm(tubi.fullFileBase.im_uv, timePoints(tidx+1))) ;
        else
            error('did not recognize coordSys (image coordinate specifier for piv)')
        end

        % Now compare adjacent images either with a simple phase
        % correlation or using PIVLab 
        if ~use_PIVLab 
            [vv,uu,xx,yy] = GetPIV(im1, im2, [], [], edgeLength, histequilize) ;
    
            %% Filtering
            do_stdev_check  = 1 ;
            stdthresh       = 7  ;
            do_local_median = 1  ;
            neigh_thresh    = 3  ;
            
            if isfield(options, 'do_stdev_check') 
                do_stdev_check = options.do_stdev_check ;
            end
            if isfield(options, 'stdthres')
                stdthresh   = options.stdthres  ;
            end
            if isfield(options, 'do_local_median')
                do_local_median = options.do_local_median  ;
            end
            if isfield(options, 'neigh_thres')
                neigh_thresh = options.neigh_thres  ;
            end
            
            % local median check
            if do_local_median
                neigh_filt=medfilt2(u,[3,3],'symmetric');
                neigh_filt=inpaint_nans(neigh_filt);
                neigh_filt=abs(neigh_filt-u);
                u(neigh_filt>neigh_thresh)=nan;

                neigh_filt=medfilt2(v,[3,3],'symmetric');
                neigh_filt=inpaint_nans(neigh_filt);
                neigh_filt=abs(neigh_filt-v);
                v(neigh_filt>neigh_thresh)=nan;
            end
            
            % stddev check
            if do_stdev_check==1
                meanu=nanmean(u(:));
                meanv=nanmean(v(:));
                std2u=nanstd(reshape(u,size(u,1)*size(u,2),1));
                std2v=nanstd(reshape(v,size(v,1)*size(v,2),1));
                minvalu=meanu-stdthresh*std2u;
                maxvalu=meanu+stdthresh*std2u;
                minvalv=meanv-stdthresh*std2v;
                maxvalv=meanv+stdthresh*std2v;
                u(u<minvalu)=NaN;
                u(u>maxvalu)=NaN;
                v(v<minvalv)=NaN;
                v(v>maxvalv)=NaN;
            end
            
            u(isnan(v))=NaN;
            v(isnan(u))=NaN;
            u_filt=inpaint_nans(u,4);
            v_filt=inpaint_nans(v,4);
            
        else
            % Analyze the image with piv_FFTmulti
            disp('Performing PIV analysis with deforming windows and 4 passes...')

            % TYPICAL GUI PIPELINE
            % % Open PIVLab
            % % Select all frames in meshDir/PullbackImages_XXXstep_sphi/smoothed_extended/
            % % Select Sequencing style 1-2, 2-3, ... 
            % %
            % % Below are settings in PIVlab that work well for this dataset:
            % % Image Preprocessing (used to select all, but now:)
            % %  --> Enable CLAHE with 20 pix
            % %  --> DO NOT Enable highpass with 15 pix
            % %  --> DO NOT Enable Intensity capping
            % %  --> Wiener2 denoise filter with 3 pix
            % %  --> DO NOT Auto constrast stretch
            % % PIV settings: 
            % %  --> 128 (32 step), 64 (32 step), 32 (16 step), 16 (8 step)
            % %  --> Linear window deformation interpolator
            % %  --> 5x repeated correlation 
            % %  --> Disable auto-correlation
            % % Post-processing
            % %  --> Standard deviation filter: 7 stdev
            % %  --> Local median filter: thres=5, eps=0.1
            % %  --> Interpolate missing data
            % % Export 
            % %  --> File > Save > MAT file

            % Standard PIV Settings
            intArea1        = 128 ;
            step            = 32 ;
            subpixFindr     = 1  ;
            mask            = [] ;
            roi             = [] ;
            numPasses       = 4  ;
            intArea2        = 64 ;
            intArea3        = 32 ;
            intArea4        = 16 ;
            repeat          = 1  ;
            disAuto         = 0  ;
            % Image proc
            clahe           = 1  ;
            claheW          = 40 ;
            highPass        = 0  ;
            highPassSz      = 15 ;
            clipping        = 0  ;
            wiener          = 1  ;
            wienerW         = 3  ;
            % Post-processing        
            calu            = 1. ;
            calv            = 1  ;
            valid_vel       = [] ;
            do_stdev_check  = 1 ;
            stdthresh       = 7  ;
            do_local_median = 1  ;
            neigh_thresh    = 3  ;
            if isfield(options, 'intArea1')
                step = options.intArea1 ;
            end
            if isfield(options, 'subpixFindr')
                subpixFindr = options.subpixFindr ;
            end
            if isfield(options, 'mask')
                mask = options.mask ;
            end
            if isfield(options, 'roi')
                roi         = options.roi ;
            end
            if isfield(options, 'numPasses')
                numPasses   = options.numPasses  ;
            end
            if isfield(options, 'intArea2')
                intArea2    = options.intArea2 ;
            end
            if isfield(options, 'intArea3')
                intArea3    = options.intArea4 ;
            end
            if isfield(options, 'intArea4')
                intArea4    = options.intArea4 ;
            end
            if isfield(options, 'repeat')
                repeat      = options.repeat  ;
            end
            if isfield(options, 'disAuto')
                disAuto     = options.disAuto  ;
            end
            if isfield(options, 'clahe')
                clahe       = options.clahe  ;
            end
            if isfield(options, 'claheW')
                claheW      = options.claheW ;
            end
            if isfield(options, 'highPass')
                highPass    = options.highPass ;
            end
            if isfield(options, 'highPassSz')
                highPassSz  = options.highPassSz ;
            end
            if isfield(options, 'clipping')
                clipping    = options.clipping  ;
            end
            if isfield(options, 'wiener')
                wiener      = options.wiener  ;
            end
            if isfield(options, 'wienerW')
                wienerW     = options.wienerW  ;
            end
            if isfield(options, 'valid_vel')
                valid_vel   = options.valid_vel ;
            end
            if isfield(options, 'do_stdev_check') 
                do_stdev_check = options.do_stdev_check ;
            end
            if isfield(options, 'stdthres')
                stdthresh   = options.stdthres  ;
            end
            if isfield(options, 'do_local_median')
                do_local_median = options.do_local_median  ;
            end
            if isfield(options, 'neigh_thres')
                neigh_thresh = options.neigh_thres  ;
            end

            % To make it more readable, let's create a "settings table"
            s = cell(15,2); 
            %Parameter                       %Setting           %Options
            s{1,1}= 'Int. area 1';           s{1,2}=intArea1;   % window size of first pass
            s{2,1}= 'Step size 1';           s{2,2}=step;       % step of first pass
            s{3,1}= 'Subpix. finder';        s{3,2}=subpixFindr;% 1 = 3point Gauss, 2 = 2D Gauss
            s{4,1}= 'Mask';                  s{4,2}=mask;       % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
            s{5,1}= 'ROI';                   s{5,2}=roi;        % Region of interest: [x,y,width,height] in pixels, may be left empty
            s{6,1}= 'Nr. of passes';         s{6,2}=numPasses;  % 1-4 nr. of passes
            s{7,1}= 'Int. area 2';           s{7,2}=intArea2 ;  % second pass window size
            s{8,1}= 'Int. area 3';           s{8,2}=intArea3 ;  % third pass window size
            s{9,1}= 'Int. area 4';           s{9,2}=intArea4 ;  % fourth pass window size
            s{10,1}='Window deformation';    s{10,2}='*spline'; % '*spline' is more accurate, but slower
            s{11,1}='Repeated Correlation';  s{11,2}=repeat;    % 0 or 1 : Repeat the correlation four times and multiply the correlation matrices.
            s{12,1}='Disable Autocorrelation';  s{12,2}=disAuto;% 0 or 1 : Disable Autocorrelation in the first pass. 
            s{13,1}='Correlation style';     s{13,2}=0;         % 0 or 1 : Use circular correlation (0) or linear correlation (1).
            s{14,1}='Repeat last pass';      s{14,2}=repeat ;   % 0 or 1 : Repeat the last pass of a multipass analyis
            s{15,1}='Last pass quality slope';  s{15,2}=0.1;    % Repetitions of last pass will stop when the average difference to the previous pass is less than this number.

            % Standard image preprocessing settings
            p = cell(8,1);
            %Parameter                       %Setting           %Options
            p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
            p{2,1}= 'CLAHE';                 p{2,2}=clahe;      % 1 = enable CLAHE (contrast enhancement), 0 = disable
            p{3,1}= 'CLAHE size';            p{3,2}=claheW;     % CLAHE window size
            p{4,1}= 'Highpass';              p{4,2}=highPass;   % 1 = enable highpass, 0 = disable
            p{5,1}= 'Highpass size';         p{5,2}=highPassSz; % highpass size
            p{6,1}= 'Clipping';              p{6,2}=clipping;          % 1 = enable clipping, 0 = disable
            p{7,1}= 'Wiener';                p{7,2}=wiener;          % 1 = enable Wiener2 adaptive denaoise filter, 0 = disable
            p{8,1}= 'Wiener size';           p{8,2}=wienerW;          % Wiener2 window size
            p{9,1}= 'Minimum intensity';     p{9,2}=0.0;          % Minimum intensity of input image (0 = no change) 
            p{10,1}='Maximum intensity';     p{10,2}=1.0;         % Maximum intensity on input image (1 = no change)

            % Standard image postprocessing settings
            r = cell(6,1);
            %Parameter     %Setting                                     %Options
            r{1,1}= 'Calibration factor, 1 for uncalibrated data';      r{1,2}=calu;                % Calibration factor for u
            r{2,1}= 'Calibration factor, 1 for uncalibrated data';      r{2,2}=calv;                % Calibration factor for v
            r{3,1}= 'Valid velocities [u_min; u_max; v_min; v_max]';    r{3,2}=valid_vel ;          % Maximum allowed velocities, for uncalibrated data: maximum displacement in pixels
            r{4,1}= 'Stdev check?';                                     r{4,2}=do_stdev_check;      % 1 = enable global standard deviation test
            r{5,1}= 'Stdev threshold';                                  r{5,2}= stdthresh;          % Threshold for the stdev test
            r{6,1}= 'Local median check?';                              r{6,2}=do_local_median;     % 1 = enable local median test
            r{7,1}= 'Local median threshold';                           r{7,2}=neigh_thresh;        % Threshold for the local median test

            % PIV analysis:
            % for syntax see https://github.com/Shrediquette/PIVlab/blob/main/Accuracy.m
            try
                image1 = PIVlab_preproc (im1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2},p{9,2},p{10,2}); %preprocess images
            catch
                error('PIVlab image preprocessing failed -- is PIVLab installed?')
            end
            image2 = PIVlab_preproc (im2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2},p{9,2},p{10,2});
            tic % start timer for PIV analysis only
            cell2table(s)
            try
                disp(['computing PIV for tidx = ' num2str(tidx) '/' num2str(length(timePoints))])
                [xx, yy, uu, vv, ~] = piv_FFTmulti (image1, image2,...
                    s{1,2},s{2,2}, s{3,2}, s{4,2}, s{5,2}, s{6,2}, s{7,2} ,s{8,2},...
                    s{9,2},s{10,2},s{11,2},s{12,2},s{13,2},0,s{14,2},s{15,2});
            catch
                error('pivlab FFTmulti failed -- make sure PIVLab is up to date.')
            end

            %% Postprocessing
            [u_filt,v_filt] = PIVlab_postproc(uu,vv, r{1,2}, r{2,2}, r{3,2}, r{4,2},...
                r{5,2},	r{6,2},	r{7,2}) ;

            % typevector_filt = typevector; % initiate
            % typevector_filt(isnan(u_filt))=2;
            % typevector_filt(isnan(v_filt))=2;
            % typevector_filt(typevector==0)=0; %restores typevector for mask
            u_filt=inpaint_nans(u_filt,4);
            v_filt=inpaint_nans(v_filt,4);

        end
        u{tidx} = uu ;
        v{tidx} = vv ;
        x{tidx} = xx ;
        y{tidx} = yy ;
        u_filtered{tidx} = u_filt ;
        v_filtered{tidx} = v_filt ;
    end
    
    % Save results
    disp(['Saving 2d PIV to ' tubi.fileName.pivRaw.raw])
    save(tubi.fileName.pivRaw.raw, 'x', 'y', ...
        'u', 'v', 'u_filtered', 'v_filtered')
else
    disp('PIV already exists on file and overwrite = false, so skipping PIV.')
end

