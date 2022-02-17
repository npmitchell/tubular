function [FF, VV] = closeMeshEndcaps(mesh, nU, nV)
% Close the endcaps of a topological cylinder rectilinear mesh
%
% Parameters
% ----------
% mesh : struct with fields 
%   f : 
%   v : 
%   nU :
%   nV :
%
%

if nargin < 2
    nU = mesh.nU ;
end
if nargin < 3
    nV = mesh.nV ;
    try
        assert(size(mesh.v, 1) == (nV-1)*nU)
    catch
        error('Mesh is not nU x nV sized. This function assumes cylindrical mesh of nU x (nV-1) vertices')
    end
end

vtx = reshape(mesh.v, [mesh.nU, mesh.nV-1, 3]) ;
nvtx = size(mesh.v, 1) ;
endpt1 = mean(squeeze(vtx(1, :, :)), 1) ;
endpt2 = mean(squeeze(vtx(end, :, :)), 1) ;
end1 = 1:nU:nU*(nV-1) ;
endf1 = [];
for qq = 1:length(end1)
    if qq + 1 <= numel(end1) 
        endf1 = [endf1; [end1(qq), end1(qq+1), nvtx + 1]] ;
    else
        endf1 = [endf1; [end1(qq), end1(mod(qq+1, numel(end1))), nvtx + 1]] ;
    end
end
end2 = nU:nU:nU*(nV-1) ;
endf2 = [];
for qq = 1:length(end2)
    if qq + 1 <= numel(end2) 
        endf2 = [endf2; [end2(qq), nvtx+2, end2(qq+1)]] ;
    else
        endf2 = [endf2; [end2(qq), nvtx+2, end2(mod(qq+1, numel(end2)))]] ;
    end
end

% add endcap points to vtx and triangulation
VV = [mesh.v; endpt1; endpt2] ;
FF = [mesh.f; endf1; endf2] ;