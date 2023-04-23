function prepareIlastik(tubi)
    % Prepaere stack for Ilastik segmentation. This outputs an h5
    % file of subsampled intensity data at isotropic resolution
    % on which to train.
    %
    % prepareIlastik(stack)

    % Accoring to the specified options sub-sample the stack and
    % save it for analysis with ilastik.

    % Handle the two cases differently: if xp is an ImSAnE object or simply
    % a struct with metadata and options.
    if isa(tubi.xp, 'struct')

        xp = tubi.xp ;
        fn = xp.detectOptions.fileName ;

        for tt = tubi.xp.fileMeta.timePoints
            if ~exist(fullfile(tubi.dir.data, [sprintf(fn, tt) '.h5']), 'file')
                disp(['Did not find file: ', fullfile(tubi.dir.data, [sprintf(fn, tt) '.h5'])])

                tubi.setTime(tt) ;

                % make a copy of the detectOptions and change the fileName
                opts = xp.detectOptions ;
                opts.fileName = sprintf( fn, tubi.currentTime ) ;

                fileName = [opts.fileName,'.h5'];

                dsetName = '/inputData';

                % Determine the axis order for ilastik training data
                if strcmp(opts.preilastikaxisorder, 'xyzc') 
                    axperm = [1 2 3 4] ;
                elseif strcmp(opts.preilastikaxisorder, 'yxzc')
                    axperm = [2 1 3 4] ;
                elseif strcmp(opts.preilastikaxisorder, 'zxyc')
                    axperm = [3 1 2 4] ;
                elseif strcmp(opts.preilastikaxisorder, 'czxy')
                    axperm = [4 3 1 2] ;
                elseif strcmp(opts.preilastikaxisorder, 'czyx')
                    axperm = [4 3 2 1] ;
                elseif strcmp(opts.preilastikaxisorder, 'cxyz')
                    axperm = [4 1 2 3] ;
                else
                    error(['Have not coded for this axis permutation yet: ', ...
                        opts.preilastikaxisorder])
                end

                % Subsample the image to save for ilastik training
                adjustIV = false ;
                im = tubi.getCurrentData(adjustIV) ;

                for c = 1:length(im)
                    if c==1
                        tmp = im{c}(1:opts.ssfactor:end,1:opts.ssfactor:end,1:opts.ssfactor:end) ;
                        image = zeros(size(tmp, 1), size(tmp, 2), size(tmp, 3), length(im), class(im{c})) ;
                    end
                    image(:,:,:,c) = im{c}(1:opts.ssfactor:end,1:opts.ssfactor:end,1:opts.ssfactor:end);
                end

                % Now save the subsampled images to h5 using the axis order
                % specified by axperm
                if ndims(image)==4
                    disp(['Writing file: ' fileName])
                    if isa(image, 'uint8')
                        h5create(fileName,dsetName,[size(image,axperm(1)) size(image,axperm(2)), ...
                            size(image,axperm(3)) size(image,axperm(4))], ...
                            'datatype', 'uint8',...
                            'Chunksize', [size(image,axperm(1)) size(image,axperm(2)),...
                            size(image,axperm(3)) size(image,axperm(4))]);
                    elseif isa(image, 'uint16')
                        h5create(fileName,dsetName,[size(image,axperm(1)) size(image,axperm(2)),...
                            size(image,axperm(3)) size(image,axperm(4))], ...
                            'Datatype', 'uint16');
                    else
                        error('Did not recognize bitdepth of image. Add capability here')
                    end
                    h5write(fileName,dsetName,permute(image, axperm));
                else
                    % truncate the axis permutation to include just 3 dims
                    axperm = axperm(1:3) ;
                    disp(['Writing file: ' fileName])

                    if isa(image, 'uint8')
                        h5create(fileName,dsetName,[size(image,axperm(1)) size(image,axperm(2)),...
                            size(image,axperm(3))], ...
                            'datatype', 'uint8',...
                            'Chunksize', [size(image,axperm(1)) size(image,axperm(2)), size(image,axperm(3)) ]) ; 
                    elseif isa(image, 'uint16')
                        h5create(fileName,dsetName,...
                            [size(image,axperm(1)) size(image,axperm(2)), size(image,axperm(3))], ...
                            'Datatype', 'uint16')
                    else
                        error('Did not recognize bitdepth of image. Add capability here')
                    end
                    h5write(fileName,dsetName,permute(image, axperm));
                end

                % Tell user we are finished
                disp(['done outputting downsampled data h5: tp=' num2str(tt) ' for surface detection'])

            else
                disp(['h5 for timepoint ' num2str(tt) ' was already output, skipping...'])
            end
        end    
    else
        % xp is an Experiment class with its own methods
        disp('Detected that tubi.xp is an ImSAnE Experiment class instance')
        fn = tubi.xp.detect.options.fileName ;
        xp = tubi.xp ;
        % is this right??
        projectDir = xp.projectDir ;
        
        for tt = xp.fileMeta.timePoints
            if ~exist(fullfile(projectDir, [sprintf(fn, tt) '.h5']), 'file')
                disp(['Did not find file: ', fullfile(projectDir, [sprintf(fn, tt) '.h5'])])
                xp.loadTime(tt);
                xp.rescaleStackToUnitAspect();
                % make a copy of the detectOptions and change the fileName
                detectOpts2 = detectOptions ;
                detectOpts2.fileName = sprintf( fn, xp.currentTime ) ;
                xp.setDetectOptions( detectOpts2 );
                xp.detector.prepareIlastik(xp.stack);
                disp(['done outputting downsampled data h5: tp=' num2str(tt) ' for surface detection'])
            else
                disp(['h5 ' num2str(tt) ' was already output, skipping...'])
            end
        end    
        disp('Open with ilastik if not already done')
    end
end