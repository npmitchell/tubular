function doubleCoverToSingleCoverMesh(xy, umax, vmax)
% doubleCoverToSingleCoverMesh(xy, umax, vmax)
% This code has not been vetted/testing thoroughly but I think it
% should work.
%
%
% Parameters
% ----------
% xy : #vertices x 2 float array
%   vertex positions of a mesh in the pullback space 
% umax : float
%   extent of the X (u) direction in pullback space
% vmax : float 
%   extent of the Y (v) direction in pullback space
% 
% Returns
% -------
% 
% NPMitchell 2020

error('This code has not been vetted/testing thoroughly but I think it should work.')

% Check that mask limiting to singleCover does not introduce
% duplicates. First map coordinates to (0, 1), (-0.5, 1.5)
modxy = xy ;
modxy(:, 1) = modxy(:, 1) / umax ;  % normalizing u coord
modxy(:, 2) = modxy(:, 2) / vmax ;  % normalizing v coord
% Now map into the singleCover space
modxy(:, 2) = mod(modxy(:, 2), 1) ;
[~, mD] = knnsearch(modxy(singleCover, :),...
                       modxy(singleCover, :), 'K', 2) ;
if any(mD(:, 2) < 1e-7)
    error('duplicate points included in singleCover. Handle here')
end

% make singleCover version of connectivity (faces_singleCover),
% wrapping into a topological cylinder (perhaps with some holes, 
% but wraps around)
[mIDx, ~] = knnsearch(modxy, modxy, 'K', 2) ;
% Associate each non-singleCover index in faces 
% with a singleCover redundant point if one exists 
% (one may not exist near boundary of image).
% finc --> "face in cover"
finc = ismember(faces, find(singleCover)) ;

% remove all faces without any singleCover vertices
faces1 = faces(any(finc, 2), :) ;
finc = finc(any(finc, 2), :) ;

% Now replace indices in faces for each zero in new finc.
% redundant maps each particle to its nearest periodic image
% This should map to periodic images to the right particle in
% the singleCover space. 
% Grab the index of mIDx that is not 1:length(mIDx) -- that is, 
% find the redundant particle, not the singleCover particle itself.
% This is almost certainly the second column of mIDx.
% However, a more careful way is to set the self to zero and sum
% the row, adding back the offset which set the self to zero.
redundant = sum(mIDx - (1:length(mIDx))', 2) + (1:length(mIDx))' ;

% Check: this should equal mIDx(:, 2) if periodicity is not
% absolutely perfect (ie if mD(:, 2) > 0)
% assert(all(redundant == mIDx(:, 2)))

plot(xy(redundant(~singleCover), 1), xy(redundant(~singleCover), 2), '.')

% Check: redundant particles outside singleCover should map 
% into singleCover 
map_inside = ismember(redundant(find(~singleCover)), find(singleCover)) ;

% Tidy up faces_singleCover by replacing redundant points with
% their singleCover values
faces1(~finc) = redundant(faces1(~finc)) ;
[ faces_singleCover, xy_singleCover, ~] = ...
    remove_vertex_from_mesh( faces1, xy, find(~singleCover) ) ;
assert(all(all(xy_singleCover == xy(singleCover, :))))
