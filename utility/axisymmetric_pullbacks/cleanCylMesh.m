function [mesh] = cleanCylMesh(mesh, options)
%CLEANCYLMESH(mesh) Orient faces, smooth normals, remove ears & hairpins
%   Given a cylindrical mesh, adjust the boundaries to not have vertices
%   with only two neighboring vertices. That is, all boundary vertices
%   have both bulk and boundary vertices as neighbors. Also, if there are
%   hairpin turns in the boundary, so that the indices vary as 
%   [i j k l m j n], then the sequence [klmj] is removed -> [ijn].
%   Additionally, faces are ensured to be properly orientated and normals
%   are smoothed.
%
% Parameters
% ----------
% mesh : cylinderMesh with fields v, f
% options : optional struct with fields
%   niter : int
%       max number of iterations to seek hairpin-free boundaries. If zero,
%       then we do not fix hairpin turns in the boundary           
% 
% Returns
% -------
% mesh : cleaned cylinderMesh with fields v, f, vn
% 
% NPMitchell 2019

%% Default options & unpacking
niter = 1000 ;  % max iterations for eradicating hairpin boundaries
if nargin > 1
    if isfield(options, 'niter')
        niter = options.niter;
    end
end

%% Consistently orient mesh faces
disp('orienting faces')
mesh.f = bfs_orient( mesh.f );

%% Compute vertex normals by weighting the incident triangles by
% their incident angles to the vertex
% This uses gptoolbox.
mesh.vn = per_vertex_normals(mesh.v, mesh.f, 'Weighting', 'angle') ;
% Average normals with neighboring normals
disp('Averaging normals with neighboring normals')
mesh.vn = average_normals_with_neighboring_vertices(mesh, 0.5) ;

%% Clip the ears of the triangulation and update the AD/PD points ---------
disp('Clipping ears...') ;
[ ff, ~, vv, newind] = clip_mesh_mod( mesh.f, mesh.v );
% Remove zeros from newind in order to assign the new vtx normals
newind_keep = newind(newind > 0) ;
if all(diff(newind_keep) > 0)
    mesh.vn = mesh.vn(newind > 0, :) ;
else
    error('Index array is non-monotonic, handle slightly fancier case here')
end

%% Clip hairpin turns in the boundary to remove dangling components -------
if niter > 0
    fb = freeBoundary(triangulation(ff, vv)) ;
    [pieces, hairpins] = separateFreeBoundaryCurves(fb) ;
    
    try
        assert(length(hairpins)==2)
    catch
        error(['There are ' num2str(length(hairpins)) ' boundary loops'])
    end
    hairpin_exists = false ;
    for qq = 1:2
        if ~isempty(hairpins{qq}) && ~hairpin_exists
            hairpin_exists = true ;
            % Keep only the largest loop on the boundary
            nloops = length(hairpins{qq}) ;
            loopL = zeros(nloops, 1) ;
            for pp = 1:nloops
                loopL(pp) = length(hairpins{qq}{pp}) ;
            end
            % which is the biggest loop?
            [~, bigID] = max(loopL) ;
            loops2rm = setdiff(1:nloops, bigID) ;
            
            pts2remove = [] ;
            for l2rm = loops2rm
                pts2remove = [pts2remove, hairpins{qq}{l2rm}] ; 
            end
        end
    end
    while hairpin_exists && kk < niter
        disp(['Removing hairpins: ' num2str(pts2remove) ' vertices'])
        [ff, vv, keepID] = remove_vertex_from_mesh(ff, vv, pts2remove) ;
        mesh.vn = mesh.vn(keepID, :) ;
    end
end

mesh.f = ff ;
mesh.v = vv ;

%% % Check it and the normals
% tmp = mesh ;
% gone = find(newind == 0) ;
% inds = 1:min(gone)
% ind = min(gone) ;
% tmp.v = tmp.v + tmp.vn * (130) ;
% trisurf(tmp.f, tmp.v(:, 1), tmp.v(:, 2), tmp.v(:, 3), 'EdgeColor', 'none')
% hold on
% axis equal
%
% Show which points were removed
% plot3(tmp0.v(gone, 1), tmp0.v(gone, 2), tmp0.v(gone, 3), 'ro')
% axis equal
% 
% tmp2 = mesh ;
% % plot3(tmp2.v(:, 1), tmp2.v(:, 2), tmp2.v(:, 3), '.')
% trisurf(tmp2.f, tmp2.v(:, 1), tmp2.v(:, 2), tmp2.v(:, 3), 'EdgeColor', 'none')
% hold on;
% quiver3(tmp2.v(:, 1), tmp2.v(:, 2), tmp2.v(:, 3), tmp2.vn(:, 1), tmp2.vn(:, 2), tmp2.vn(:, 3), 10)
% %quiver3(tmp2.v(inds, 1), tmp2.v(inds, 2), tmp2.v(inds, 3), tmp2.vn(inds, 1), tmp2.vn(inds, 2), tmp2.vn(inds, 3), 10)
% quiver3(tmp2.v(ind, 1), tmp2.v(ind, 2), tmp2.v(ind, 3), tmp2.vn(ind, 1), tmp2.vn(ind, 2), tmp2.vn(ind, 3), 10, 'color', 'r')
% axis equal
% mags = vecnorm(mesh.vn - mesh2.vn(1:length(mesh.vn), :), 2, 1)


end

