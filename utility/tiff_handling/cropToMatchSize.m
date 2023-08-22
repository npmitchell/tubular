function dat2crop = cropToMatchSize(dat2crop, dataRef)
% dat2crop = cropToMatchSize(dat2crop, dataRef)
%   
% Parameters
% ----------
% dat2crop : data volume to crop/pad with zeros to match size of dataRef
% dataRef : data volume whose size we match in the crop/padded output
% 
% Returns 
% -------
% dat2crop : output data whose size matches dataRef
%

    if ~all(size(dat2crop) == size(data))

        % crop or pad dim 1
        aa = size(dat2crop, 1) ;
        bb = size(dataRef, 1) ;
        if aa > bb
            dat2crop = dat2crop(1:bb, :, :) ;
        elseif aa < bb
            dat2crop(aa+1:bb, :, :) = 0 ;
        end

        % crop or pad dim 2
        aa = size(dat2crop, 2) ;
        bb = size(data, 2) ;
        if aa > bb
            dat2crop = dat2crop(:, 1:bb, :) ;
        elseif aa < bb
            dat2crop(:, aa+1:bb, :) = 0 ;
        end

        % crop or pad dim 3
        aa = size(dat2crop, 3) ;
        bb = size(data, 3) ;
        if aa > bb
            dat2crop = dat2crop(:, :, 1:bb) ;
        elseif aa < bb
            dat2crop(:, :, aa+1:bb) = 0 ;
        end
    end

end