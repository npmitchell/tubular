function sampleIDx = ...
    farthestPointVertexSampling(nPts, VV, FF, seedIDx, options)
%farthestPointVertexSampling()
%   Create an approximately uniform sampling of vertices by iterative
%   farthest point search. Requires gptoolbox for fast marching.
%
% Parameters
% ----------
% nPts : int
%   the number of vertices to identify in the mesh
% VV : #vertices x D float array
%   vertices of the mesh to subsample
% FF : #faces x 3 int array
%   faces of the mesh to subsample
% seedIDx : N x 1 int array
%   the vertices to include in initial sampling
% options : optional struct with fields
%   preview : bool
%       display progress and command line output
%
% Returns
% -------
%
% NPMitchell 2020

preview = false ;
if nargin > 4
    if isfield(options, 'preview')
        preview = options.preview ;
    end
end

% The vertex ID of each sample point
sampleIDx = zeros(nPts, 1);
sampleIDx(1:numel(seedIDx)) = seedIDx;

% Will hold the distances of all mesh vertices to the final sample set
% D(i,j) is the distance of the ith mesh vertex to the jth seed point
D = inf(size(VV, 1), nPts);

for i = 1:numel(seedIDx)
    if size(VV, 2) == 2
        VV(:, 3) = 0 * VV(:, 1) ;
    end
    if preview
        disp(['fast marching #' num2str(i) '/' num2str(numel(seedIDx))])
    end
    D(:,i) = perform_fast_marching_mesh(VV', FF, seedIDx(i));
end

% The minimum distance of each vertex in the mesh to any of the extant
% sample points
Dmin = min(D, [], 2);

if preview, progressbar(i, nPts); end
while (i < nPts)
    
    i = i+1;
    
    if preview, progressbar(i, nPts); end
    
    [~, newID] = max(Dmin);
    assert(~ismember(newID, sampleIDx), ...
        'Update vertex is already a sample point!');
    sampleIDx(i) = newID;
    
    % Update the minimum distances
    D(:,i) = perform_fast_marching_mesh(VV, FF, newID);
    Dmin = min(Dmin, D(:,i));
    
end
