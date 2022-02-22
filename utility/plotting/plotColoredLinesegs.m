function plotColoredLinesegs(lsegs, ecolors, varargin)
% PLOTCOLOREDLINESEGS(lsegs, ecolors, varargin)
%
% Parameters
% ----------
% lsegs : N x 4 or Nx6 float array
%   the linesegments, given as (x0,y0,x1,y1) or (x0,y0,z0,x1,y1,z1)
% ecolors : N x 3 int array
%   color specifications for each lineseg
% varargin : variable input arguments
%   passed to plot() or plot3()
%
% Returns
% -------
% 
%
% NPMitchell 2020

lw = 1 ;
verbose = false ;
for ii = 1:length(varargin)
    if isa(varargin{ii}, 'double')
        continue;
    end
    if isa(varargin{ii}, 'logical')
        continue;
    end
    
    % Target geometry parameters ------------------------------------------
    if ~isempty(regexp(varargin{ii},'^[Ll]ine[Ww]idth', 'match'))
        lw = varargin{ii+1};
    end    
   
    if ~isempty(regexp(varargin{ii},'^[Vv]erbose', 'match'))
        verbose = varargin{ii+1};
    end    
end

orighold = ishold ;
hold on; 

% Mark which linesegs/rows have already been plotted
not_done = ones(length(lsegs(:, 1)), 1) ;

dmyk = 1 ;
while any(not_done)
    % Bin each set of matching colors
    newcolor = ecolors(find(not_done, 1, 'first'), :) ;
    batch = find(all(ecolors == newcolor, 2));
    if verbose
        disp(['Plotting batch ' num2str(dmyk) ' with newcolor ' ])
        disp(ecolors(batch(1), :))
    end
    
    if dmyk > length(ecolors(:, 1))
        dmyk
        size(ecolors)
        length(ecolors(:, 1))
        error('Number of color batches has exceeded number of colors')
    end

    switch size(lsegs, 2)
        case 6
            plot3([lsegs(batch, 1)'; lsegs(batch, 4)'], ...
                [lsegs(batch, 2)'; lsegs(batch, 5)'], ...
                [lsegs(batch, 3)'; lsegs(batch, 6)'], ...
                'color', ecolors(batch(1), :), 'linewidth', lw)
        case 4
            plot([lsegs(batch, 1)'; lsegs(batch, 3)'], ...
                [lsegs(batch, 2)'; lsegs(batch, 4)'], ...
                'color', ecolors(batch(1), :), 'linewidth', lw)
    end
    
    % Mark the ones we've finished
    not_done(batch) = 0 ;
    dmyk = dmyk + 1 ;
end

if ~orighold
    hold off;
end