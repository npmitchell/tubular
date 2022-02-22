function ax = phasebar(varargin) 
% phasebar places a circular or half-circular) donunt-shaped colorbar for 
% phase from -pi to pi, -180 degrees to 180 degrees, OR
%       from 0 to pi, 0 degrees to 180 degrees.
% If vargarin has 'style', styleStr,... and if lower(styleStr) contains 
%   'fill', then the inner radius of the donut is zero, so it is a half 
%   circle or a full circle.
%% Example Usage
% 
% cbs{1} = phasebar('colormap', phasemap, ...
%     'location', [0., 0.0, 1, 1], ...
%     'style', 'nematicFill') ;
% 
%% Syntax
% 
%  phasebar
%  phasebar(...,'location',Location) 
%  phasebar(...,'size',Size) 
%  phasebar('deg') 
%  phasebar('rad') 
%  ax = phasebar(...) 
% 
%% Description 
% 
% phasebar places a donut-shaped colorbar on the current axes. 
%
% phasebar(...,'location',Location) specifies the corner (e.g., 'northeast' or 'ne') 
% of the current axes in which to place the phasebar. Default location is the upper-right or 'ne' 
% corner. 
%
% phasebar(...,'size',Size) specifies a size fraction of the current axes.  Default is 0.3. 
%
% phasebar('deg') plots labels at every 90 degrees. 
%
% phasebar('rad') plots labels at every pi/2 radians. 
%
% phasebar('tickpos', 'inside') puts tick marks inside the phasebar annulus 
%
% phasebar('style', 'nematic') displays nematic order [0, pi]
%
% phasebar('colormap', my_cmap) uses custom colormap
%
% ax = phasebar(...) returns a handle ax of the axes in which the new axes 
%    are plotted. 
% 
%% Example
% 
% Z = 200*peaks(900); 
% Zw = phasewrap(Z,'degrees'); 
% imagesc(Zw) 
% phasemap(12)
% phasebar('location','se')
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas 
% at Austin's Institute for Geophysics (UTIG), May 2016. 
% This function includes Kelly Kearney's plotboxpos function as a subfunction. 
% NPMitchell added nematic and filled (half)circle functionality.
% 
% If the phasemap function is useful for you, please consider citing our 
% paper about it: 
% 
% Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 
% 2016. True colors of oceanography: Guidelines for effective and accurate 
% colormap selection. Oceanography 29(3):9?13. 
% http://dx.doi.org/10.5670/oceanog.2016.66
% 
% Edits and improvements by NPMitchell 2020 including
%   --> font size adjustment, placement of ticklabels
%   --> colormap options for non-default colormaps
%   --> TickPos = ('outside', 'inside')
%   --> style = ('polar', 'nematic')
%
% See also colorbar and phasemap. 


%% Set Defaults: 
usedegrees = false; 
axsize = 0.3; 
location = 'northeast'; 
cm = colormap; 
TickPos = 'outside'; 
style = 'polar' ;

% Try to automatically determine if current displayed data exist and are
% in radians or degrees: 
if max(abs(caxis))>(2 * pi + 1e-14)
   usedegrees = true;
else
   usedegrees = false; 
end

% If no data are already displayed use radians: 
if isequal(caxis,[0 1])
   usedegrees = false; 
end

%% Parse inputs: 

tmp = strncmpi(varargin,'location',3); 
if any(tmp) 
   location = varargin{find(tmp)+1}; 
end

tmp = strncmpi(varargin,'size',3); 
if any(tmp) 
   axsize = varargin{find(tmp)+1}; 
   assert(isscalar(axsize)==1,'Input error: axis size must be a scalar greater than zero and less than one.') 
   assert(axsize>0,'Input error: axis size must be a scalar greater than zero and less than one.') 
   assert(axsize<1,'Input error: axis size must be a scalar greater than zero and less than one.') 
end

if any(strncmpi(varargin,'radians',3))
   usedegrees = false; 
end

if any(strncmpi(varargin,'degrees',3)) 
   usedegrees = true; 
end

tmp = strncmpi(varargin,'colormap',3); 
if any(tmp)
    cm = varargin{find(tmp) + 1} ;
end

tmp = strncmpi(varargin, 'style',3); 
if any(tmp)
    style = varargin{find(tmp) + 1} ;
end

%% style: filled or not
if contains(lower(style), 'fill')
    innerRadius = 0 ;
else
    innerRadius = 10 ;
end

%% Starting settings: 

currentAx = gca; 
pos = plotboxpos(currentAx); 
xcol = get(currentAx,'XColor'); 

% Delete old phasebar if it exists: 
try
   oldphasebar = findobj(gcf,'tag','phasebar'); 
   delete(oldphasebar); 
catch
end

%% Created gridded surface: 
outerRadius = 10*1.618; 

[x,y] = meshgrid(linspace(-outerRadius,outerRadius,300));
[theta,rho] = cart2pol(x,y); 

% theta = rot90(-theta,3); 
theta(rho>outerRadius) = nan; 
theta(rho<innerRadius) = nan; 
if contains(lower(style), 'nematic')
    theta(y < 0) = nan ;
end

if usedegrees
   theta = theta*180/pi; 
end

%% Plot surface: 

ax = axes; 
if contains(lower(style), 'grad')
    alphaGrid = rho./outerRadius ;
    alphaGrid(alphaGrid > 1) = 1 ;
    ph = surf(x,y,theta,'FaceAlpha','flat',...
        'AlphaDataMapping','scaled',...
        'AlphaData',alphaGrid) ;
    view(2)
    axis equal
else
    ph = pcolor(x,y,theta) ;    
end
shading interp
hold on

if contains(lower(style), 'polar')
    % Plot a ring: 
    [xc1,yc1] = pol2cart(linspace(-pi,pi,360),innerRadius); 
    [xc2,yc2] = pol2cart(linspace(-pi,pi,360),outerRadius); 
    plot(xc1,yc1,'-','color',xcol,'linewidth',.2); 
    plot(xc2,yc2,'-','color',xcol,'linewidth',.2); 
elseif contains(lower(style), 'nematic')
    % Plot a ring: 
    [xc1,yc1] = pol2cart(linspace(0,pi,180),innerRadius); 
    [xc2,yc2] = pol2cart(linspace(0,pi,180),outerRadius); 
    plot(xc1,yc1,'-','color',xcol,'linewidth',.2); 
    plot(xc2,yc2,'-','color',xcol,'linewidth',.2);     
end

axis image off
colormap(gca, cm) 

if contains(lower(style), 'polar')
    if usedegrees
       caxis([-180 180]) 
    else
       caxis([-pi pi]) 
    end
elseif contains(lower(style), 'nematic')
    if usedegrees
       caxis([0 180]) 
    else
       caxis([0 pi]) 
    end
else
    error('Could not parse style definition')
end

%% Label: 
if contains(lower(style), 'polar')
    % Option 1: use inner radius 
    if strcmpi(TickPos, 'inside')
        [xt,yt] = pol2cart((-1:2)*pi/2+pi/2,innerRadius); 
        if usedegrees
           text(xt(1),yt(1),'0\circ','horiz','right','vert','middle'); 
           text(xt(2),yt(2),'90\circ','horiz','center','vert','top'); 
           text(xt(3),yt(3),'180\circ','horiz','left','vert','middle'); 
           text(xt(4),yt(4),'-90\circ','horiz','center','vert','bottom'); 
        else
           text(xt(1),yt(1),'0','horiz','right','vert','middle'); 
           text(xt(2),yt(2),'\pi/2','horiz','center','vert','top'); 
           text(xt(3),yt(3),'\pi','horiz','left','vert','middle'); 
           text(xt(4),yt(4),'-\pi/2','horiz','center','vert','bottom'); 
        end
    elseif strcmpi(TickPos, 'outside')
        % Option 2: use outer radius 
        [xt,yt] = pol2cart((-1:2)*pi/2+pi/2,outerRadius*1.05); 
        if usedegrees
           text(xt(1),yt(1),'0\circ','horiz','left','vert','middle'); 
           text(xt(2),yt(2),'90\circ','horiz','center','vert','bottom'); 
           text(xt(3),yt(3),'180\circ','horiz','right','vert','middle'); 
           text(xt(4),yt(4),'-90\circ','horiz','center','vert','top'); 
        else
           text(xt(1),yt(1),'0','horiz','left','vert','middle'); 
           text(xt(2),yt(2),'\pi/2','horiz','center','vert','bottom'); 
           text(xt(3),yt(3),'\pi','horiz','right','vert','middle'); 
           text(xt(4),yt(4),'-\pi/2','horiz','center','vert','top'); 
        end
    else
        error('Bad Tick Position definition')
    end
    
elseif contains(lower(style), 'nematic')
    % Option 1: use inner radius 
    if strcmpi(TickPos, 'inside')
        [xt,yt] = pol2cart((0:2)*pi/2,innerRadius); 
        if usedegrees
           text(xt(1),yt(1),'0\circ','horiz','right','vert','middle'); 
           text(xt(2),yt(2),'90\circ','horiz','center','vert','top'); 
           text(xt(3),yt(3),'180\circ','horiz','left','vert','middle'); 
        else
           text(xt(1),yt(1),'0','horiz','right','vert','middle'); 
           text(xt(2),yt(2),'\pi/2','horiz','center','vert','top'); 
           text(xt(3),yt(3),'\pi','horiz','left','vert','middle'); 
        end
    elseif strcmpi(TickPos, 'outside')
        % Option 2: use outer radius 
        [xt,yt] = pol2cart((0:2)*pi/2,outerRadius*1.05); 
        if usedegrees
           text(xt(1),yt(1),'0\circ','horiz','left','vert','middle'); 
           text(xt(2),yt(2),'90\circ','horiz','center','vert','bottom'); 
           text(xt(3),yt(3),'180\circ','horiz','right','vert','middle'); 
        else
           text(xt(1),yt(1),'0','horiz','left','vert','middle'); 
           text(xt(2),yt(2),'\pi/2','horiz','center','vert','bottom'); 
           text(xt(3),yt(3),'\pi','horiz','right','vert','middle'); 
        end
    else
        error('Bad Tick Position definition')
    end
else
    error('Style not recognized: must be polar or nematic')
end


%% Set position of colorwheel: 

try
    switch lower(location)
       case {'ne','northeast'} 
          set(ax,'position',...
              [pos(1)+(1-axsize)*pos(3), ...
              pos(2)+(1-axsize)*pos(4), ...
              axsize*pos(3), axsize*pos(4)]); 

       case {'se','southeast'} 
          set(ax,'position', ...
              [pos(1)+(1-axsize)*pos(3), pos(2),...
              axsize*pos(3), axsize*pos(4)]); 

       case {'nw','northwest'} 
          set(ax,'position',[pos(1), pos(2)+(1-axsize)*pos(4), ...
              axsize*pos(3), axsize*pos(4)]); 

       case {'sw','southwest'} 
          set(ax,'position',[pos(1), pos(2), axsize*pos(3), axsize*pos(4)]); 

       case {'neo', 'northeastoutside'} 
          set(ax,'position',[pos(1)+pos(3), pos(2)+(1-axsize)*pos(4),...
              axsize*pos(3), axsize*pos(4)]); 

       otherwise
          error('Unrecognized axis location.') 
    end
catch
    % location is given as 1 x 4 array
    set(ax,'position',location); 
end
      
%% Clean up 

set(ax,'tag','phasebar')

% Make starting axes current again: 
% axes(currentAx); 
set(gcf,'CurrentAxes',currentAx)

uistack(ax,'top'); 

% if nargout==0 
%    clear ax
% end

end




%% Kelly Kearney's plotboxpos: 

function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Copyright 2010 Kelly Kearney

% Check input

if nargin < 1
    h = gca;
end

if ~ishandle(h)  ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels

currunit = get(h, 'units');
set(h, 'units', 'pixels');
axisPos = get(h, 'Position');
set(h, 'Units', currunit);

% Calculate box position based axis limits and aspect ratios

darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    
    pos = axisPos;
    
else

    dx = diff(get(h, 'XLim'));
    dy = diff(get(h, 'YLim'));
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis
temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', get(h, 'parent'));
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);

end

