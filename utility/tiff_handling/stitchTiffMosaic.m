function [data, frame1] = stitchTiffMosaic(fns, nrow, ncol, datfn, oper,...
    dtype, maxI)
% stitchTiffMosaic(fns, nrow, ncol, datfn, oper)
%
%
% Parameters
% ----------
% fns : struct with fields
%   name
%   folder
% nrow : int
%   number of rows in the mosaic
% ncol : int 
%   number of columns in the mosaic
% datfn : str
%   path to h5 file to save whole stitched data
% oper : float (between 0-1)
%   overlap percent, as a float between 0-1
% dtype : string specifier ('uint8'=default)
% maxI : optional float (0-1)
%   clip intensity at this percentile
%
% Returns
% -------
% data
%
% NPMitchell 2020

if nargin < 6
    dtype = 'uint8' ;
end

first = true ;
for tile = 1:length(fns)
    disp(num2str(tile))
    fn = fullfile(fns(tile).folder, fns(tile).name) ;
    % Load the tiff stack
    dat = loadtiff(fn) ;
    if first
        WW = size(dat, 2) ;
        HH = size(dat, 1) ;
        Wtot = ceil(WW * ((ncol-1)*(1-oper) + 1)) ;
        Htot = ceil(HH * ((nrow-1)*(1-oper) + 1)) ;
        data = uint8(zeros(Htot, Wtot, size(dat, 3))) ;
        frame1 = uint8(zeros(Htot, Wtot)) ;
        first = false ;
    end

    % coordinates for placement
    row = floor((tile - 1) / ncol) ;
    col = mod(tile - 1, ncol);
    rpix = round((Htot - (row * HH * (1-oper))) - HH)  ;
    cpix = round(col * WW * (1 - oper));

    % check first frame
    rxx = 1 + rpix:(rpix + WW) ;
    cyy = 1 + cpix:(cpix + HH) ;

    % Allow for clipping
    % if any(rxx < 0)
    %     nclip = length(find(rxx < 0)) ;
    %     rxx = rxx(rxx > 0) ;
    % else
    %     nclip = 1 ;
    % end
    % frame1(rxx, cyy) = dat(nclip:end, :, 1) ;

    frame1(rxx, cyy) = dat(:, :, 1) ;
    data(rxx, cyy, :) = dat ;

    set(gcf, 'visible', 'on')
    imagesc(frame1)
    axis equal
    pause(0.000001)
end

% Normalize
for z = 1:size(data, 3)
    page = data(:, :, z) ;
    pagen = double(page - min(page(:))) / double(max(page(:)) - min(page(:))) ;
    if maxI > 0 
        pagen(pagen > maxI) = maxI ;
        pagen = pagen / maxI ;
    end
    if strcmp(dtype, 'uint8')
        data(:,:,z) = uint8(pagen * 255) ;
    else
        error('handle this dtype here')
    end
end

% Save the data to datfn
h5create(datfn, '/tileScan', size(data), 'datatype', dtype) 
h5write(datfn, '/tileScan', data)
