function sfsm = laplacian_smooth_faces(sf, mesh, varargin) 
%LAPLACIAN_SMOOTH_FACES(sf, mesh, Options)
% Smooth a scalar field defined on faces by replacing each face value with 
% the weighted average of neighboring face values. Optional: select
% 'Method' == 'Denoise' to smooth only local extrema, up to a value epsilon
%  
%
% Parameters
% ----------
% sf : Nfaces x 1 float array
%   the scalar field to smooth, defined on faces of mesh (mesh.f)
% mesh : struct with fields f, v
%   the mesh on which to smooth the scalar field
% varargin : keyword arguments
%   options for how to perform smoothing
%       'Weight': 'areas'
%       'Method': 'smooth' or 'denoise'
%       'Epsilon': threshold for finding peaks (optional, default = 1e-16).
%		 A vertex must be within eps of a local maximum to be denoised 
% 
% Returns  
% -------
% sfsm : Nfaces x 1 float array
%   the smoothed scalar field defined on faces of mesh (mesh.f)
% 
% NPMitchell 2020


% Method options
weight = 'area' ;
method = 'smooth' ;  % options: smooth, denoise
eps = 1e-16 ;
for i = 1:length(varargin)
    
    if isa(varargin{i}, 'double')
        continue;
    end
    if isa(varargin{i}, 'logical')
        continue;
    end
    
    if ~isempty(regexp(varargin{i}, '^[Ww]eight', 'match'))
        weight = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Mm]ethod', 'match'))
        method = varargin{i+1} ;
    end    
    if ~isempty(regexp(varargin{i}, '^[Ee]psilon', 'match'))
        eps = varargin{i+1} ;
    end    
end

% Apply the weighting for neighboring face contributions
if  ~isempty(regexp(weight, '^[Aa]rea', 'match'))
    % Weighting is by area of neighboring faces
    weights = 0.5 * doublearea(mesh.v, mesh.f) ;
end

tri = triangulation(mesh.f, mesh.v) ;
% compute face neighbors 
faceNeighbors = [ (1:size(mesh.f,1)).', tri.neighbors() ];
% Handle the NaNs, which do not have four adjacent faces 
fNfix = faceNeighbors ;
[badr, badc] = find(isnan(faceNeighbors)) ;
 
% Add a new value to the weights (zero) and to the scalar field (large 
% index) to do sum efficiently
assert(size(mesh.f, 1) == length(sf)) 
newID = size(mesh.f, 1) + 1 ;
fNfix(badr, badc) = newID ;
sf(newID) = 0 ;
weights(newID) = 0 ;
sfsm = sum(sf(fNfix) .* weights(fNfix) ./ sum(weights(fNfix), 2), 2);

% If we denoise, then we only remove speckle (local extrema, up to eps)
if ~isempty(regexp(method, '^[Dd]enoise', 'match'))
    % Check to find which values are local extrema of their neighbors
    % Find local minima
    sf(newID) = max(sf(:)) + 1 ;
    isminimum = all(sf(1:end-1) < sf(fNfix) + eps, 2) ;
    % Find local maxima
    sf(newID) = min(sf(:)) - 1 ;
    ismaximum = all(sf(1:end-1) > sf(fNfix) - eps, 2) ;
    
    % Make sure we don't add a new face value to sfout
    isminimum(newID) = false ;
    ismaximum(newID) = false ;
    
    sfout = sf(1:end-1) ;
    sfout(ismaximum | isminimum) = sfsm(ismaximum | isminimum) ;
    sfsm = sfout ;
end    

return
