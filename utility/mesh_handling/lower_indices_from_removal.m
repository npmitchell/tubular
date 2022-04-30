function newList = lower_indices_from_removal(indices, oldVertices, newVertices, rmVIDx )
% Adjust indices in list to reflect removal of a subset of points from a
% collection.
%
% Parameters
% ----------
% newVertices (N-Q)xD float array
%   the set which has some subset removed
% rmVIDx : Qx1 int array
%   the indices removes in going from oldVertices to newVertices
% 
% Returns
% -------
% newList :  
%
% NPMitchell 2020

% Find the new index for each of the new vertices
[~, newVertexIDx] = ismember( oldVertices, newVertices, 'rows');

% Find any faces that contained the indicated vertices and remove them
newList = indices;

loc_v = false(size(face));
for i = 1:numel(rmVIDx)
    loc_v = loc_v | ( newList == rmVIDx(i) );
end
loc_v = find( loc_v );

[row_v, ~] = ind2sub(size(newList), loc_v);
row_v = unique(row_v);
newList(row_v, :) = [];

% Now update the vertex indices in the index list
newList = newVertexIDx(newList);

if length(newList) < length(indices)
    numremoved = length(indices) - length(newList) ;
    disp(['Removed ' num2str(numremoved) 'elements from input indices'])
end

