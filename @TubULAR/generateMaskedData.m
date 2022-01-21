function generateMaskedData(QS, options)
% Given training on midgut cells and spurious (amnioserosa) cells as 
% h5 files in dir16bit/stabilized_h5s/, mask the original data and save 3d
% volumes.
% Save output h5s trained on stabilized h5s from iLastik as 
%   -->   <QS.fileBase.name>_Probabilities_mask3d.h5
% and 
%   -->   <QS.fileBase.name>_Probabilities_maskDorsal.h5
% Note: if QS.plotting.preview == True, displays intermediate results
%
% NPMitchell 2020

% Unpack QS/options
maskMask = true ;
overwrite = false ;
ssfactorMask = QS.ssfactor ;
maskedDir = QS.dir.maskedData ;
if nargin > 1
    if isfield(options, 'ssfactorMask')
        ssfactorMask = options.ssfactorMask ;
    end
    if isfield(options, 'maskMask')
        maskMask = options.maskMask ;
    end
    if isfield(options, 'overwrite')
        overwrite = options.overwrite ;
    end
end

% Parameters
cliplowDorsal = 0.3 ;
step = 10 ;

% Naming
maskBase = [QS.fileBase.name '_Probabilities_mask3d.h5'] ;
dorsalBase = [QS.fileBase.name '_Probabilities_maskDorsal.h5'] ;
maskDir = fullfile(QS.dir.data, 'stabilized_h5s') ;

% Change adjustment to allow for dim pixels (needed for MIPS)
old_adjusthigh = QS.data.adjusthigh ;
QS.data.adjusthigh = 99.999 ;


% Check how many channels
nChannels = QS.xp.expMeta.channelsUsed ;

% First create an average stab image to use for training
preview = QS.plotting.preview ;
for tt = QS.xp.fileMeta.timePoints
    disp(['t = ', num2str(tt)])
    
   % Check if already done
   redo_calc = false ;
    for ii = 1:nChannels
        if nChannels > 1
            out16name_ii = fullfile(maskedDir, ...
                sprintf([QS.fileBase.name '_masked.tif'], tt, ii)) ;
        else
            out16name_ii = fullfile(maskedDir, ...
                sprintf([QS.fileBase.name '_masked.tif'], tt)) ;
        end
        if ~exist(out16name_ii, 'file') || overwrite
            redo_calc = true ;
        end
    end
    
    if redo_calc
        maskFn = fullfile(maskDir, sprintf(maskBase, tt)) ;
        if maskMask
            dorsalFn = fullfile(maskDir, sprintf(dorsalBase, tt)) ;
        end

        % Mask data, rotate to lateral and ventral views
        try
            maskprob = h5read(maskFn, '/exported_data') ;
        catch
            msg = ['Create masked iLastik h5 output as ' maskFn] ;
            error(msg)
        end
        mask = squeeze(maskprob(1, :, :, :)) > 0.5 ;
        % Take largest connected component, erode and dilate.
        se2 = strel('sphere', 2);
        se3 = strel('sphere', 3);
        se4 = strel('sphere', 4);
        mask0 = imerode(mask, se3) ;
        mask0 = imdilate(mask0, se4) ;
        if preview
            for qq = 1:step:size(mask0, 3)
                title('eroded & dilated volume')
                imshow(squeeze(mask0(:, :, qq)))
                title('eroded & dilated volume')
            end
        end
        disp('erosion3 / dilation4')
        mask0 = imerode(mask0, se3) ;
        mask0 = imdilate(mask0, se4) ;
        if preview
            for qq = 1:step:size(mask0, 3)
                title('eroded & dilated volume -- pass 2')
                imshow(squeeze(mask0(:, :, qq)))
                title('eroded & dilated volume -- pass 2')
            end
        end
        disp('Isolating single component')
        CC = bwconncomp(mask0) ;
        numPixels = cellfun(@numel,CC.PixelIdxList) ;
        [~, bigID] = max(numPixels) ;
        biggest = zeros(size(mask0), 'logical') ;
        biggest(CC.PixelIdxList{bigID}) = 1 ;
        if preview
            for qq = 1:step:size(biggest, 3)
                title('biggest volume')
                imshow(squeeze(biggest(:, :, qq)))
                title('biggest volume')
            end
        end

        % % Now fill holes by keeping only the biggest outside region as zero
        % mask1 = (1-biggest) ;
        % CC = bwconncomp(mask1) ;
        % numPixels = cellfun(@numel,CC.PixelIdxList) ;
        % [vol, bigID] = max(numPixels) ;
        % bwsolid = ones(size(mask1), 'logical') ;
        % bwsolid(CC.PixelIdxList{bigID}) = 0 ;
        % if preview
        %     for qq = 1:step:size(bwsolid, 3)
        %         title('filled solid')
        %         imshow(squeeze(bwsolid(:, :, qq)))
        %     end
        % end

        % blur it & threshold to smoothen
        disp('blurring mask')
        % One final dilation before blurring and subtracting dorsal
        biggest = imdilate(biggest, se2) ;
        blurred = imgaussfilt3(single(biggest), 3) ;
        % bwsolid = blurred > 0.5 ;


        if preview
            for qq = 1:step:size(blurred, 3)
                imshow(squeeze(blurred(:, :, qq)))
            end
        end

        % active contour to smoothen more
        % bwsolid = activecontour(squeeze(maskprob(1, :, :, :)), bwsolid, ...
        %     10, 'Chan-Vese',...
        %     'ContractionBias', -0.1, 'SmoothFactor', 5);
        % disp('rescaling to full resolution')
        % mask_final = superkron(bwsolid, ones(ssfactorMask, ssfactorMask, ssfactorMask)) ;
        % % One final dilation
        % mask_final = imdilate(mask_final, se3) ;

        % Subtract off dorsal
        if maskMask
            disp('subtracting dorsal amnioserosa cells')
            dorsalprob = h5read(dorsalFn, '/exported_data') ;
            dorsalprob = squeeze(dorsalprob(1, :, :, :)) ;

            % Optional: smoothen the dorsal prob
            maskd = dorsalprob > cliplowDorsal ;
            % erode and dilate.
            maskd = imerode(maskd, se2) ;
            maskd = imdilate(maskd, se3) ;
            % convert back to uint16
            dorsalprob = imgaussfilt3(single(maskd), 5) ;
            % 
            % % clip the dorsal at low values and rescale LUT to (0,1)
            % disp('clipping dorsal')
            % dorsalprob = dorsalprob - cliplowDorsal ;
            % dorsalprob(dorsalprob < 0) = 0 ;
            % dorsalprob = single(double(dorsalprob) / double(max(dorsalprob(:)))) ;

            if preview
                for qq = 1:step:size(blurred, 3)
                    imshow(squeeze(dorsalprob(:, :, qq)))
                end
            end

            % multiply by probabilities
            prob_nondorsal = single(blurred) - dorsalprob ;
            prob_nondorsal(prob_nondorsal < 0) = 0 ;
            prob_nondorsal(prob_nondorsal > 1) = 1 ;

            if preview 
                for qq = 1:step:size(blurred, 3)
                    title('final mask')
                    imshow(squeeze(prob_nondorsal(:, :, qq)))
                    title('final mask')
                end
            end

            disp('rescaling to full resolution')
            mask_final = superkron(prob_nondorsal, ones(ssfactorMask, ssfactorMask, ssfactorMask)) ;
        else
            mask_final = superkron(blurred, ...
                ones(ssfactorMask, ssfactorMask, ssfactorMask)) ;
        end

        %% Apply the mask to the data
        disp(['Loading t = ' num2str(tt)])
        QS.setTime(tt)
        QS.getCurrentData() ;
        IV = QS.currentData.IV ;
        IVmasked = IV ;

        % clip the mask if larger than data
        if any(size(mask_final) > size(IV{1}))
            mask_final = mask_final(1:size(IV{1}, 1), ...
                1:size(IV{1}, 2), 1:size(IV{1}, 3)) ; 
        elseif any(size(mask_final) < size(IV{1}))
            error('handle this case here')
        end

        % inspect it
        nChannels = length(IV) ;
        for ii = 1:nChannels
            if nChannels > 1
                out16name = fullfile(maskedDir, ...
                    sprintf([QS.fileBase.name '_masked.tif'], tt, ii)) ;
            else
                out16name = fullfile(maskedDir, ...
                    sprintf([QS.fileBase.name '_masked.tif'], tt)) ;
            end
            IVii = IV{ii} ;
            IVmasked{ii} = uint16(single(IVii) .* single(mask_final)) ;

            % check it
            if preview
                clf
                for qq = 1:step:size(IVmasked{ii}, 3)
                    page = 0.5*squeeze(IVmasked{ii}(:, :, qq)) ;
                    bwpage = 0.5*squeeze(uint16(mask_final(:, :, qq)*2^16)) ;
                    imshow(cat(3, bwpage, page, bwpage + page))
                    % imshow(bwpage)
                end
            end

            % Save each image
            if ~exist(out16name, 'file') || overwrite
                disp(['writing data to ' out16name])
                out = IVmasked{ii} ;
                for z = 1:size(out, 3)
                    if mod(z, 100) == 0
                        disp(['page ' num2str(z) ' / ' num2str(size(out,3))])
                    end
                    % write the first page as overwrite mode, then append
                    if z == 1
                        imwrite(out(:,:,z), out16name, 'tiff',...
                            'Compression','none');
                    else
                        imwrite(out(:,:,z), out16name, 'tiff',...
                            'Compression','none', ...
                            'WriteMode','append');    
                    end
                end
                disp('done now')
            else
                disp(['file already exists: ' out16name])
            end
        end
    end
end

% Convert back to old adjust limit
QS.data.adjusthigh = old_adjusthigh ;
