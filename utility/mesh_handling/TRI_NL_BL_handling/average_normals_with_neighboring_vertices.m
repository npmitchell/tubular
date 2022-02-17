function vnx = average_normals_with_neighboring_vertices(mesh, onsite_weight)
%AVERAGE_NORMALS_WITH_NEIGHBORING_VERTICES(MESH, ONSITE_WEIGHT)
%   Average the normals of neighboring vertices with each vertex's normal.
%   Gove the N neighbors a weight of (1-onsite_weight) / N.
%   Give the vertex's original normal a weight of onsite_weight.
%
% Parameters
% ----------
% mesh : struct with fields v, f, vn (vn is optional)
%   The mesh on which to average normals
% onsite_weight : float (0-1)
%   The weight to give the original normal vs mean of neighbors' normals
%
% Returns
% -------
% vnx : N x 3 float array
%   The normal vectors at each vertex
%
% NPMitchell KITP 2019

% Get neighbors for each vertex
NLcell = TRI2NLcell(mesh.f, mesh.v) ;

if isfield(mesh, 'vn')
    vnx = mesh.vn ;
else
    % Note: must use angle weighting for mathematical rigor
    vnx = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
end

% Average each vertex normal with those of its neighbors
% Slow way
vnq = 0 * vnx ;
for vqq = 1:length(NLcell)
    vnq(vqq, :) = mean(vnx(NLcell{vqq}, :), 1) ;
end
vnx = onsite_weight * vnx + (1 - onsite_weight) * vnq ;
% Re-normalize the vertex array
vnx = vnx ./ vecnorm(vnx, 2, 2) ;

end
