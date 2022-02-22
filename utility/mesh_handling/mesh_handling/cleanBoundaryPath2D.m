function [cleantri] = cleanBoundaryPath2D(TRI, xy, cutPath, keep_ears, tri_is_ordered, allow_system_spanning_tris)
%CLEANBOUNDARYPATH(TRI, xy, cutPath) Prune boundary of triangulation 
%   Given triangulation, ensure that each boundary index specified in 
%   cutPath(ii) connects to no other points in cutPath other than 
%   cutPath(ii-1) and cutPath(ii+1) unless encloses a right handed "ear"
%   triangle that dangles at the boundary along an edge (if keep_ears ==
%   true).
%   TODO: Modify so that a face could span the whole system, from one
%   boundary edge to the other side. This will be easy: Count how many
%   edges are on the current (pre-pruned) boundary. If only one, keep. 
%   If two, remove. If three, then the triangle is isolated & disconnected.
%   
% Parameters
% ----------
% TRI : M x 3 int array
%   The indices into vertices of each triangular face
% vertices : Q x D float array
%   The spatial positions of the mesh vertices
% cutPath : N x 1 int array
%   Indices into vertices for the path to prune. Must be oriented
%   counterclockwise if closed or if ears (ii-1,ii,ii+1) are present unless
%   the ears are to be pruned.
% keep_ears : bool
%   Allow ears that are right-handed, like [cutPath(ii-1), cutPath(ii),
%   cutPath(ii+1)], when the right hand rule from r12->r23 is satisfied
% tri_is_ordered : bool
%   For speedup, can tell this function that the triangle faces have
%   already been ordered in a right handed fashion so (r12xr23)_3 > 0
% allow_system_spanning_tris : bool
%   Do not automatically kill 'nonlocal' tris with indices that differ by
%   more than 1 index into cutPath.
%
% Returns
% -------
% cleantri : (M - few) x 3 int array
%   The indices of vertices of each triangular face after pruning along
%   cutPath.
%
% NPMitchell 2019

if nargin < 4
    keep_ears = true;
    tri_is_ordered = false ;
    allow_system_spanning_tris = false ;
elseif nargin < 5
    tri_is_ordered = false ;
    allow_system_spanning_tris = false ;
elseif nargin < 6
    allow_system_spanning_tris = false ;
end

keep = true(size(TRI, 1), 1) ;
for ii = 1:size(TRI, 1)
    if length(find(ismember(TRI(ii, :), cutPath))) > 1
        % check for connections with other boundary pts that 
        % are not cutPath(ii-1) or cutPath(ii+1)
        tri_id = find(ismember(TRI(ii, :), cutPath)) ;
        if length(tri_id) == 3
            % Find if this is a right handed ear for counterclockwise path
            % cutPath
            if keep_ears
                % Determine if right or left handed not by ordering of TRI,
                % but by ordering of cutPath
                idt1 = find(cutPath == TRI(ii, 1)) ;
                idt2 = find(cutPath == TRI(ii, 2)) ;
                idt3 = find(cutPath == TRI(ii, 3)) ;
                idts = [idt1 idt2 idt3] ;
                % First consider nonlocal tris, such as...
                % [ii-1 ii ii+N] or [ii-N ii ii+1]
                if any(abs(diff(idts)) > 1)
                    if allow_system_spanning_tris
                        error('Have not written this case: count freeBoundary edges of this tri')
                        % Count how many edges are on the current (pre-pruned) boundary. If only one, keep. 
                        % If two, remove. If three, then the triangle is isolated & disconnected.
                    else
                        % Kill ears like [ii-1 ii ii+N]
                        disp(['pruning nonlocal tri [' num2str(TRI(ii, :)) ']'])
                        keep(ii) = false ;
                    end
                else
                    % Consider ears [ii-1 ii+1 ii] ... 
                    % or [ii-1 ii+N ii] for N > 1
                    if tri_is_ordered
                        % If all diffs > 0, then must be ii-1 ii ii+1
                        if any(diff(idts) < 0)
                            % Kill ears like [ii-1 ii+1 ii] or
                            % [ii-1 ii+N ii] for N > 1
                            disp(['pruning left-handed ear [' num2str(TRI(ii, :)) ']'])
                            keep(ii) = false;
                        end
                    else
                        % order according to cutPath
                        r12 = [xy(cutPath(idt2), :) - xy(cutPath(idt1), :), 0] ;
                        r23 = [xy(cutPath(idt3), :) - xy(cutPath(idt2), :), 0] ;
                        % r12 = [xy(TRI(ii, 2), :) - xy(TRI(ii, 1), :), 0] ;
                        % r23 = [xy(TRI(ii, 3), :) - xy(TRI(ii, 2), :), 0] ;
                        xprod = cross(r12, r23) ;
                        if xprod(3) < 0
                            % The ear is left handed. Remove it.
                            disp(['pruning left-handed ear [' num2str(TRI(ii, :)) ']'])
                            keep(ii) = false ;
                        end
                    end
                end
            else
                disp(['pruning ear [' num2str(TRI(ii, :)) ']'])
                keep(ii) = false ;
            end
         elseif length(tri_id) == 2
             % Kill the triangle if the two members are not adjacent in
             % cutPath list.
             pathid = find(ismember(cutPath, TRI(ii, :))) ;
             if abs(diff(pathid)) > 1
                 disp(['pruning triangle [' num2str(TRI(ii, :)) ']'])
                 keep(ii) = false;
             end
        else
            error('Only one match in TRI(ii, :). Should not end up here')
        end
    end
end

cleantri = TRI(keep, :) ;

% Check the results
triplot(cleantri, xy(:, 1), xy(:, 2))



% Old version uses NL
% NL = TRI2NL(TRI, vertices) ;
% for ii = 1:length(cutPath)
%     % check for connections with other boundary pts that 
%     % are not cutPath(ii-1) or cutPath(ii+1)
%     row = NL(cutPath(ii), :) ;
%     % Are elements of row in cutPath up to ii-1 or after ii+1?
%     totrim = [ cutPath(1:(ii-2)); cutPath((ii+2):end) ] ;
%     totrim = find(ismember(row, totrim)) ;
%     if length(totrim) > 1
%         disp(['trimming bonds of cutPath point ' num2str(cutPath(ii))])
%         NL(cutPath(ii), totrim) = 0 ;
%     end
% end
% 
% % Converting 
% BL = NL2BL(NL) ;
% cleantri = NL2TRI(NL, BL) ;

end

