function [vIdx, fIdx] = tiledAnnularCutMesh2SingleCover(TF, TV, tileCount, pathPairs) 
% [vIdx, fIdx] = tiledAnnularCutMesh2SingleCover(TF, TV, tileCount, pathPairs) 
% Note that this may seem a silly function, since our convention is that
%     vidx = 1:size(cutMesh.v, 1) ;
%     fidx = 1:size(cutMesh.f) ;
% would be reasonable, but if we are trying to extract something after
% smoothing of a periodic mesh, this is not the case.
%
% Parameters
% ----------
% TF : tiled face list
% TV : tiled vertex list
% tileCount : 2 x 1 int array
%   # tiles above and below the original single cover
% pathPairs : #vertices in seam x 2 int array
%   indices of single cover used to glue seam
%
% Returns
% -------
% vIdx : #vertices in single cover x 1 int array 
% fIdx : #faces in single cover x 1 int array
%
% See also
% --------
% tileAnnularCutMesh()
%
% NPMitchell 2021

if all(tileCount == [1 1])
    % vidx = 1:size(cutMesh.v, 1) ;
    % fidx = 1:size(cutMesh.f) ;
    
    
else
    error('handle this case here')
end
