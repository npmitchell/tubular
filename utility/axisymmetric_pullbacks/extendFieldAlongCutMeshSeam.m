function field = extendFieldAlongCutMeshSeam(field_no_seam, cutMesh)
% For field defined on vertices, but not yet defined on the seam of
% cutMesh, extend the field from the lower pathPairs to the higher pathPair
% indices
%
%

field = field_no_seam ;
field(cutMesh.pathPairs(:, 2), :) = field(cutMesh.pathPairs(:, 1), :) ;
