function [colors, names] = define_colors(varargin)
%DEFINE_COLORS Define a set of pleasant colors
%   Define a set of colors to use

if length(varargin) > 0
    ncolors = varargin{1} ;
    disp(['Defining ' num2str(ncolors) ' colors...'])
else
    ncolors = 10 ;
end

blue    = [0.0000, 0.4470, 0.7410] ; % 1
red     = [0.8500, 0.3250, 0.0980] ; % 2
yellow  = [0.9290, 0.6940, 0.1250] ; % 3
purple  = [0.4940, 0.1840, 0.5560] ; % 4 
green   = [0.4660, 0.6740, 0.1880] ; % 5 
sky     = [0.3010, 0.7450, 0.9330] ; % 6 
maroon  = [0.6350, 0.0780, 0.1840] ; % 7
gray    = [0.2000, 0.2000, 0.2000] ; % 8
brown   = [0.5400, 0.2500, 0.0900] ; % 9
teal    = [0.0000, 0.6500, 0.5200] ; % 10
pink    = [1.0000, 0.5137, 0.5137] ; % 11
dark_green =  [0.0392, 0.5059, 0.2745] ;   % 12

colors = [blue; red; yellow; purple; green; ...
          sky; maroon; gray; brown; teal; pink; dark_green] ;

if ncolors < 11
    colors = colors(1:ncolors, :) ;
elseif colors < 14
    next = [0.5 * (blue + red); 0.5 * (red+yellow); 0.5*(yellow+purple)] ;
    colors = [colors; next];
    colors = colors(1:ncolors, :) ;
else
    error('Need to define more colors for this')
end

if nargout > 1
    names = {'blue', 'red', 'yellow', 'purple', 'green', ...
        'sky', 'maroon', 'gray', 'brown', 'teal', 'pink', 'dark_green'};
    names = names(1:ncolors) ;
end

end

