function [pieces, hairpins] = separateClosedCurveComponents(lsegs)
%[fb_pieces, hairpins] = separateFreeBoundaryCurves(fb)
%   Separate each topologically distinct closed curve of a free boundary 
%   or of vertex indices of curves in space into a separate list of free 
%   boundary vertex ids. 
%   If a closed curve contains a closed curve (ie a "hairpin" component)
%   that is a subset of the full path, then we return this information if 
%   desired in the second output.
%
%   For example, in the image below, there are two closed curves, one of
%   which contains two pieces, so the output is:
%      fb_pieces = {[ijkljmi], [npqn]}
%      hairpins = {{[ijmi],[jklj]}, {}}
%
%            k             p
%           o             o 
%   i    j / \           / \
%   o-----o---o l     n o---o q
%    \   /
%     \ /
%      o m
%
% Parameters
% ----------
% lsegs : #path points x 2 int
%   linesegment indices, so for the image above the input could be 
%   [i,j; j,k; k,l; l,j; j,m; n,p; p,q; pn] or equivalent. 
%
% NPMitchell 2020

% Input handling
if size(lsegs, 1) == 2
    lsegs = lsegs' ;
else
    try
        assert(size(lsegs, 2) ==2)
    catch
        error('Input linesegments must be Nx2 or 2xN int array of indices')
    end
end

% Now snip the linesegments into pieces
qq = 1;
continue_separating = true ;
while continue_separating
    [row, col] = find(lsegs == lsegs(1)) ;
    try
        assert(col(end) == 2)
    catch
        error('First startpoint in linesegment list is not repeated as endpoint')
    end
    pieces{qq} = lsegs(1:max(row), :) ;
    if max(row) == size(lsegs, 1)
        continue_separating = false ;
    else
       % Repeat the above for the remaining curve(s)
       lsegs = lsegs((max(row) + 1):end, :);
       qq = qq + 1 ;
    end
end

% We now have separate curves in pieces (a cell array). If desired, look
% for loops/hairpins
if nargout > 1
    disp('Determining hairpin loops')
    hairpins = cell(length(pieces), 1) ;
    for qq = 1:length(pieces)
        pq = pieces{qq} ;
        pt = pq' ;
        % Check that all linesegments connect in sequence
        dpt = diff(pt(:)) ;
        try
            assert(all(dpt(2:end:2) == 0))
        catch
            error(['Linesegments in piece ' num2str(qq) ' do not connect in sequence'])
        end
        
        % Find any repeats in column 1
        [uni, uID] = unique(pq(:, 1)) ;
        if length(uni) < size(pq, 1) 
            % Get the hairpin "pinching" node(s) (the repeated node(s))
            repeatID = setdiff(1:size(pq,2), uID) ;
            
            % For each hairpin, find first and second instance
            for pp = 1:length(repeatID)
                pp_loc = find(pq == repeatID(pp)) ;
                if length(pp_loc) > 2
                    error('handle multiple leaves of hairpin curve here')
                else
                    loop1 = [1:(pp_loc(1)-1), pp_loc(2):size(pq, 1)] ;
                    hairpins{qq}{1} = pq(loop1, :) ;
                    loop2 = pp_loc(1):(pp_loc(2)-1) ;
                    hairpins{qq}{2} = pq(loop2, :) ;
                end
            end
        else
            hairpins{qq} = {} ;
        end
    end
end

