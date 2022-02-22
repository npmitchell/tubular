function markers = define_markers(nmarkers)
%DEFINE_MARKERS(nmarkers)
%   Return cell array of markerstyles (cyclic if nmarkers is large)
%   
% Parameters
% ----------
% nmarkers : int
%   The number of markers to hold in the array
%
% Returns
% -------
% markers : 1 x nmarkers cell array 
%   The symbols for linespec for each kind of marker
%
% NPMitchell 2019

markers_avail = {'o', 's', '^', '+', 'd', '*',  'v', 'x', '>', 'p', '<', 'h'} ;

% if nmarkers > number of markers available, set up a cycle here
if nmarkers < length(markers_avail)
    markers = markers_avail(1:nmarkers) ;
else
    disp('Warning: number of markers requested exceeds distinct built-in markers available, cycling...')
    % preallocate
    markers = cell(1, nmarkers) ;
    navail = length(markers_avail) ;
    for qq = 1:nmarkers
        idx = mod(qq, navail) ;
        if idx == 0
            idx = navail ;
        end
        markers{qq} = markers_avail{idx} ;
    end
end

end

