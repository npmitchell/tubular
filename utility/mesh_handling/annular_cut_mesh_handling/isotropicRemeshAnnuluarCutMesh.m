function [ glueMeshnew, bc, qF] = ...
    isotropicRemeshAnnuluarCutMesh(...
    glueMesh, targetEdgeLength, numIterations, options) 
%[Fnew, Vnew, Unew, Fnnew, Vnnew] = isotropicRemeshWithPeriodicSeam(...
%   FF, VV, UU, targetEdgeLength, numIterations, options) 
%
% Assumes CGAL_Code/isotropic_remeshing.cpp is in path and compiled.
% Assumes Y is periodic with translation of 2*pi as in a Ricci mesh /
% angle.
%
% Parameters
% ----------
% glueMesh : struct with fields
%   f : #faces x 3 int array
%       face connectivity list before remeshing
%   v : #vertices x 3 float array
%       vertices of cutMesh before remeshing. Cut seam has two copies of each
%       vertex. 
%   u : #vertices x 2 float array
%       pullback 2D coordinates of glueMesh (topological annulus)
%       before remeshing. Note that input state is topological annulus,
%       unlike in isotropicRemeshAnnuluarCutMeshWithLog()
% targetEdgeLength : float
%   target edge length to acheive uniformly in remeshing
% numIterations : int
%   number of iterations to acheive isotropic remeshing
% options : struct with fields
%   preview : bool 
%       preview results as we go
%   epsVtx : float
%       very small value, for assertion that projected resample points are
%       essentially already on the original mesh faces / surface. 
%
% Returns
% -------
% glueMeshnew : struct with fields
%   
% bc : barycentric coordinates such that, for ex: 
%     Unew = bc(:, 1) .* UU(FF(qF, 1), :) + ...
%            bc(:, 2) .* UU(FF(qF, 2), :) + ...
%            bc(:, 3) .* UU(FF(qF, 3), :) ;
% qF : faces on which query points reside, useful for ex for:
%     Unew = bc(:, 1) .* UU(FF(qF, 1), :) + ...
%            bc(:, 2) .* UU(FF(qF, 2), :) + ...
%            bc(:, 3) .* UU(FF(qF, 3), :) ;
%

epsVtx = 1e-5 ;
preview = false ;
if nargin > 3
    if isfield(options, 'epsVtx')
        epsVtx = options.epsVtx;
    end
    if isfield(options, 'preview')
        preview = options.preview;
    end
end

[Fnew, Vnew, Fnnew, Vnnew] = ...
       isotropic_remeshing(glueMesh.f, glueMesh.v, targetEdgeLength, numIterations) ;
[sqrD, qF, qV] = point_mesh_squared_distance(Vnew, glueMesh.v, glueMesh.f ) ;
% check that after remeshing, pushing the remeshed vertices onto the
% original mesh requires almost no motion in 3d
try
    assert(all(vecnorm(qV - Vnew, 2, 2) < epsVtx))
catch
    error(['Isotropic remeshing created vertices which left the original surface by more than ' num2str(epsVtx)])
end
Vnew = qV ;
glueMeshnew = struct('f', Fnew, 'v', Vnew, 'vn', Vnnew, 'fn', Fnnew) ;

% initial remeshing -- cut the seam
% normalsTemp = cat(2, 0*Vnew(:, 1), 0*Vnew(:, 1), ones(size(Vnew(:, 1)))) ;
normsTemp = per_vertex_normals(Vnew, Fnew) ;

% Get barycentric coordinates via glueMesh, but find Unew via cutMesh,
% which has the same number of faces, but the seam vertex ids are
% different. That doesn't affect the computation, however! 
Vg = glueMesh.v ;
Fg = glueMesh.f ;
Ug = glueMesh.u ;
bc = barycentric_coordinates(qV, Vg(Fg(qF, 1), :), ...
    Vg(Fg(qF, 2), :), Vg(Fg(qF, 3), :)) ;

if any(any(bc < 0))
    error('Some negative barycentric coordinates! Why?')
end
Unew = bc(:, 1) .* Ug(Fg(qF, 1), :) + ...
     bc(:, 2) .* Ug(Fg(qF, 2), :) + ...
     bc(:, 3) .* Ug(Fg(qF, 3), :) ;
 
% Find Unew for cutMesh version in pb space so we can find gg,bb
glueMeshnew.u = Unew ;

% check remeshing -- glueMesh
if preview
    plot(Unew(:, 1), Unew(:, 2), '.'); hold on
    scatter(Unew(:, 1), Unew(:, 2), 5, 1:size(Unew,1), 'filled')
    trisurf(triangulation(Fnew, cat(2, Unew, 0*Unew(:, 1))))
    pause(0.2)
end
