function expandSecondAxesRow(axs, ymove)
% expandSecondAxesRow(axs, ymove)
% For a 2-row set of axes, expand the axis positions of the second row to
% match those of the first row.
%
% Parameters
% ----------
% axs
% ymove
%
% Returns
% -------
%
% Example Usage
% -------------
% [axs, cbs] = ...
%     nFieldsOnSurface({m3, m3, m2, m2}, ...
%     {real(dmu), imag(dmu), real(dmu), imag(dmu)}, options) ;
% expandSecondAxesRow(axs, ymove)
%
% NPMitchell 2021

% Unpack options
if nargin < 2
    ymove = 0 ;
end

% Determine the widths and positions of axes in top row
pos = cell(ceil(length(axs)*0.5), 1) ;
for qq = 1:ceil(length(axs)*0.5)
    set(gcf, 'currentAxes', axs{qq})
    pos{qq} = get(gca, 'position') ;
end

% Apply the widths and x positions of axes to lower row
dmyk = 1 ;
for qq = (ceil(length(axs)*0.5) + 1):length(axs)
    set(gcf, 'currentAxes', axs{qq})
    pos0{qq} = get(gca, 'position') ;
    set(gca, 'position', [pos{dmyk}(1), pos0{qq}(2)+ymove, pos{dmyk}(3), pos{dmyk}(4)]) ;
    dmyk = dmyk + 1 ;
end

