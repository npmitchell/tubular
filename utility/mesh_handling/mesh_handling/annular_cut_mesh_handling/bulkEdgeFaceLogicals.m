function [edgeIns, faceIns] = bulkEdgeFaceLogicals(faces, eIDx, uvLists, maxV)
% Keep only the edges/faces from middle swath for each rotation of annular
%   pullback in log space. For example we could have two uv maps, one with
%   branch cut at 0 and one rotated by pi, and then the logspace maps of
%   those (nonrotated and rotated) annuli would be uvLists{1} and
%   uvLists{2}.
%
% Parameters
% ----------
% faces : # faces x 3 int
%   face indices into uv
% eIDx : # edges x 2 int
%   edge indices into uv
% uvLists : #uMaps x 1 cell of #vertices x 1 complex arrays
%   complex 1d logspace coordinates for annular meshes
% maxV : float
%   cutoff for max value of Im(u) in log space pullback map of annular
%   pullback (now rectilinear)
% 
% Returns
% -------
% edgeIns : #uMaps x 1 cell of #edges x 1 booleans
%   for each map, is each edge inside the bulk range or not
% faceIns : #uMaps x 1 cell of #faces x 1 booleans
%   for each map, is each face inside the bulk range or not

if nargin < 4
    maxV = 2 ;
end
nVtx = size(uvLists{1}, 1) ;
if nVtx == 1
    nVtx = size(uvLists{1}', 1) ;
end

if length(uvLists) == 2
    v1 = imag(uvLists{1}) ;
    v2 = imag(uvLists{2}) ;
    
    all_accounted = false ;
    while ~all_accounted
        if maxV > 3.1
            msg = ['Could not find complete set of bonds in two ', ...
                'rotated frames. Decrease edge lengths or use three rotations'] ;
            error(msg)
        end
        keep1e = find(all(abs(v1(eIDx)) < maxV, 2)) ;
        keep2e = find(all(abs(v2(eIDx)) < maxV, 2)) ;
        keep1f = find(all(abs(v1(faces)) < maxV, 2)) ;
        keep2f = find(all(abs(v2(faces)) < maxV, 2)) ;
        all_accounted = isempty(setdiff(1:nVtx, unique([keep1e; keep2e]))) & ...
             isempty(setdiff(1:nVtx, unique([keep1f; keep2f]))) ;
        maxV = maxV + 1e-3 ;
    end
    disp('Accounted for all bonds and faces with two rotations')
    edgeIns = {keep1e, keep2e} ;
    faceIns = {keep1f, keep2f} ;
    
else
    error('Handle more than two rotations here')
end