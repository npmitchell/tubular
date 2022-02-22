function Fim = getFigureFrameCData(imW, imH, axisOrFigure)
% print the current fig or axis to an RGB image of determined width or size
%
% Parameters
% ----------
% imW : int
%   desired output RGB image width, in pixels
% imH : optional int
%   if not specified, figure image height will depend on current figure
%   appearance and imW
% axisOrFigure
%
% NPMitchell 2021


% Default options are:
% --> imH determined by figure appearance 
% --> whole figure printed
if nargin < 3
    axisOrFigure = 'fig' ;
end

% Capture the figure or axis
if nargin < 1 || (isempty(imW) && isempty(imH))
    % No size specified, capture the frame (fig or axis)
    if contains(lower(axisOrFigure), 'fig')
        FF = getframe(gcf) ;
    elseif contains(lower(axisOrFigure), 'ax')
        FF = getframe(gca) ;
    else
        error('Could not interpret passed axisOrFigure argument')
    end
    Fim = FF.cdata ;
else
    % First determine the resolution with which to print    
    if contains(lower(axisOrFigure), 'fig')
        set(gcf, 'units', 'inches')
        pos = get(gcf, 'position') ;
        dpi = ceil(imW / pos(3)) ;
        resolutionStr = ['-r' num2str(dpi)] ;
        Fim = print('-RGBImage',resolutionStr);
    elseif contains(lower(axisOrFigure), 'ax')
        set(gca, 'units', 'inches')
        pos = get(gca, 'position') ;
        dpi = ceil(imW / pos(3)) ;
        resolutionStr = ['-r' num2str(dpi)] ;
        Fim = print(gca, '-RGBImage',resolutionStr);
    else
        error('Could not interpret passed axisOrFigure argument')
    end
    
    % Fix overall size, to the pixel
    if nargin < 2
        % let the height be determined by the figure size and imW
        imH = round(imW * size(Fim, 2) / size(Fim, 1)) ;
    else
        if abs((size(Fim,1)/imH * imW - size(Fim, 2)) / size(Fim, 2)) > 0.05
            msg = 'Bad aspect ratio selection in imH/imW' ;
            msg = [msg '--> severe warping of image: fig=[' ] ;
            msg = [msg num2str(size(Fim,1)) ', ' num2str(size(Fim,2)) ']'] ;
            error(msg)
        end
    end
    
    % resize output image
    if ~all(size(Fim) == [imH,   imW, 3])
        Fim = imresize(Fim, [imH ,  imW]) ;
    end
    
end