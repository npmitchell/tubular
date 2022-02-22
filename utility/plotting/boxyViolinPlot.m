function boxyViolinPlot(xvals, yvals, boxWidths, ...
    boxHeights, style, faceColors, edgeColors, lineWidths)
%
% Parameters
% ----------
% xvals : length(N) numeric
%   x positions of each violin box
% yvals : length(N) numeric (default = [1])
%   y positions of each violin box
% boxWidths : N x 1 numeric or 1x1 numeric (default = diff(xvals) or 1)
%   width of each box in the plot
% boxHeights : N x 1 numeric or 1x1 numeric (default = diff(yvals) or 1)
%   height of each box in the plot
% style : string specifier ('left', 'center', default ='center')
%   if left, uses bottom left for (xvals, yvals)
%   if center, interprets (xvals, yvals) as center of box
% FaceColors
%   
% EdgeColors
%   
% LineWidths
%
% NPMitchell 2021


% Argument parsing
if nargin < 8 || isempty(lineWidths)
    lineWidths = 1 ;
end

if nargin < 7 || isempty(edgeColors)
    edgeColors = 'k' ;
end

if nargin < 6 || isempty(faceColors)
    faceColors = 'none' ;
end

if nargin < 5 || isempty(style)
    style = 'center' ;
end

if nargin < 4  || isempty(boxHeights)
    % try something reasonable for heights
    if numel(yvals) > 1
        boxHeights = diff(yvals) ;
        boxHeights = [boxHeights, boxHeights(end)] ;
    else
        boxHeights = 1 ;
    end
end

if nargin < 3  || isempty(boxWidths)
    % try something reasonable for widths
    if numel(xvals) > 1
        boxWidths = diff(xvals) ;
        boxWidths = [boxWidths, boxWidths(end)] ;
    else
        boxWidths = 1 ;
    end
end

if nargin < 2 
    yvals = 1 ;
end

% Input handling -- allow for single values in any of the supplied inputs
% if there are multiple xvals or multiple yvals
if numel(xvals) == 1 && numel(yvals) == 1
    assert(numel(boxWidths) == 1)
    assert(numel(boxHeights) == 1)
elseif numel(xvals) == 1
    xvals = xvals * ones(size(yvals)) ;
elseif numel(yvals) == 1
    yvals = yvals * ones(size(xvals)) ;
end

% If boxWidths/boxHeights is single value, match size of xvals
if numel(xvals) > 1
    % Now, if boxWidths is single value, match size of xvals
    if numel(boxWidths) == 1
        boxWidths = boxWidths * ones(size(xvals)) ;
    else
        try
            assert(numel(boxWidths) == numel(xvals)) ;
        catch
            error('Numel(boxWidths) must be 1 or same as # elements in xvals/yvals') 
        end
    end
    
    % Similarly, if boxHeightss is single value, match size of xvals
    if numel(boxHeights) == 1
        boxHeights = boxHeights * ones(size(xvals)) ;
    else
        try
            assert(numel(boxHeights) == numel(xvals)) ;
        catch
            error('Numel(boxWidths) must be 1 or same as # elements in xvals/yvals') 
        end
    end

    % Match size of xvals for faceColors, edgeColors, lineWidths
    if ischar(faceColors) && ~iscell(faceColors) 
        faceColors = repmat({faceColors}, numel(xvals),1) ;
    elseif numel(faceColors) == 3
        faceColors = faceColors .* ones(numel(xvals), 3) ;
    end
    if ischar(edgeColors) && ~iscell(edgeColors) 
        edgeColors = repmat({edgeColors}, numel(xvals),1) ;
    elseif numel(edgeColors) == 3
        edgeColors = edgeColors .* ones(numel(xvals), 3) ;
    end
    if numel(lineWidths) == 1
        lineWidths = lineWidths .* ones(size(xvals)) ;
    end
else
    % There should only be one boxWidth supplied
    try
        assert(numel(boxWidths) == numel(xvals)) ;
    catch
        error('Numel(boxWidths) must be 1 or same as # elements in xvals/yvals') 
    end
    % There should only be one boxHeight supplied
    try
        assert(numel(boxHeights) == numel(xvals)) ;
    catch
        error('Numel(boxWidths) must be 1 or same as # elements in xvals/yvals') 
    end
end
% done input handling/parsing


    
for kk = 1:length(xvals)
    % interpret current face color
    if iscell(faceColors)
        fColor = faceColors{kk} ;
    else
        fColor = faceColors(kk, :) ;
    end
    if iscell(edgeColors)
        eColor = edgeColors{kk} ;
    else
        eColor = edgeColors(kk, :) ;
    end
    
    % Plot the box
    switch style
        case 'left' 
            rectangle('Position',[xvals(kk), yvals(kk), ...
                boxWidths(kk), boxHeights(kk)],...
                'FaceColor', fColor,...
                'EdgeColor', eColor,...
                'LineWidth',lineWidths(kk))
        case 'center'
            xc = xvals(kk) -0.5*boxWidths(kk) ;
            yc =  yvals(kk) - 0.5*boxHeights(kk) ;
            rectangle('Position',...
                [xc, yc, ...
                boxWidths(kk), boxHeights(kk)],...
                'FaceColor', fColor,...
                'EdgeColor', eColor,...
                'LineWidth',lineWidths(kk))
    end
end

