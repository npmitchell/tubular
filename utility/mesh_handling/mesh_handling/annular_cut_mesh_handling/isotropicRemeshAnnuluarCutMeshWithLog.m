function [ cutMeshnew, glueMeshnew, bc, qF, badVertices] = ...
    isotropicRemeshAnnuluarCutMeshWithLog(...
    cutMesh, targetEdgeLength, numIterations, options) 
%[Fnew, Vnew, Unew, Fnnew, Vnnew] = isotropicRemeshAnnuluarCutMeshWithLog(...
%   FF, VV, UU, targetEdgeLength, numIterations, options) 
%
% Assumes CGAL_Code/isotropic_remeshing.cpp is in path and compiled.
% Assumes Y is periodic with translation of 2*pi as in a Ricci mesh /
% angle.
% This does supply a cutMesh pullback, but not all the bonds will be
% correct. Some faces will connect across the periodic bond. To avoid
% headaches, use exponentiated map using glueMesh.
%
% Parameters
% ----------
% cutMesh : struct with fields
%   f : #faces x 3 int array
%       face connectivity list before remeshing
%   v : #vertices x 3 float array
%       vertices of cutMesh before remeshing. Cut seam has two copies of each
%       vertex. 
%   u : #vertices x 2 float array
%       pullback 2D coordinates of cutMesh (topological cylinder with cut seam)
%       before remeshing. Note that input state is topological disk (rectangle)
%       with implied identification along the cutPath.
%   pathPairs : #seamVertices x 2 int array
%       indices into cutMesh.v identifying cutMesh.v(pathPairs(:, 1), :)
%       with cutMesh.v(pathPairs(:, 2), :) along cut path
% targetEdgeLength : float
%   target edge length to acheive uniformly in remeshing
% numIterations : int
%   number of iterations to acheive isotropic remeshing
% options : struct with fields
%   preview : bool 
%       preview results as we go
%   epsVtx : 
%       threshold distance between vertices for considering them as
%       separate points between vertices before and after remeshing.
%       Vertices on free boundary can (and do!) end up on the wrong side
%       of the branch cut, so this is used to identify which ones are
%       candidates to be moved to the other side.
%
% Returns
% -------
% Fnew, Vnew, 
% Unew,
% Fnnew : face normals in remeshed mesh
% Vnnew : vertex normals in remeshed mesh
% pathPairs : 
% bc : barycentric coordinates such that, for ex: 
%     Unew = bc(:, 1) .* UU(FF(qF, 1), :) + ...
%            bc(:, 2) .* UU(FF(qF, 2), :) + ...
%            bc(:, 3) .* UU(FF(qF, 3), :) ;
% qF : faces on which query points reside, useful for ex for:
%     Unew = bc(:, 1) .* UU(FF(qF, 1), :) + ...
%            bc(:, 2) .* UU(FF(qF, 2), :) + ...
%            bc(:, 3) .* UU(FF(qF, 3), :) ;
%
% See also
% --------
% isotropicRemeshAnnuluarCutMesh
% isotropicRemeshWithPeriodicSeam
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

VV = cutMesh.v ;
FF = cutMesh.f ;
UU = cutMesh.u ;
glueMesh = glueCylinderCutMeshSeam(cutMesh) ;
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

% % If we need to find bc coords of a later re-meshing in terms of this
% % remeshing, the cut the mesh along pathPairs for handling seam faces
% % Find pathPairs for unravelling
% bnew = trinew.freeBoundary() ;
% cp1 = pointMatch(VV(1, :), qV(bnew(:, 1), :)) ;
% cp2 = pointMatch(VV(nU, :), qV(bnew(:, 1), :)) ;
% cp1 = bnew(cp1) ;
% cp2 = bnew(cp2) ;
% % check that these are good points
% clf;
% trisurf(trinew, 'edgecolor', 'none', 'faceVertexAlpha', 0.1); hold on;
% scatter3(qV(bnew, 1), qV(bnew, 2), qV(bnew, 3), 10) ;
% scatter3(qV([cp1 cp2], 1), qV([cp1 cp2], 2), qV([cp1 cp2], 3), 500, 'filled')
% scatter3(VV([1 nU], 1), VV([1 nU], 2), VV([1 nU], 3), 20)

% initial remeshing -- cut the seam
% normalsTemp = cat(2, 0*Vnew(:, 1), 0*Vnew(:, 1), ones(size(Vnew(:, 1)))) ;
normsTemp = per_vertex_normals(Vnew, Fnew) ;

% Find cutpath start/endpts
[~, cp1] = min(vecnorm(Vnew - VV(cutMesh.pathPairs(1, 1), :), 2, 2)) ;
[~, cp2] = min(vecnorm(Vnew - VV(cutMesh.pathPairs(end, 1), :), 2, 2)) ;
opts = struct() ;
opts.method = 'nearest' ;
opts.path = cutMesh.v(cutMesh.pathPairs(:, 1), :) ;
cutMeshnew = cylinderCutMesh(Fnew, Vnew, normsTemp, cp1, cp2, opts) ;

% Check it
cutTrinew = triangulation(cutMeshnew.f, cutMeshnew.v) ;
cutBndnew = cutTrinew.freeBoundary() ;
if preview
    trisurf(cutTrinew) ; hold on;
    plot3(cutMeshnew.v(cutBndnew, 1), cutMeshnew.v(cutBndnew, 2), ...
        cutMeshnew.v(cutBndnew, 3), 'r.-')
    axis equal
end

% Get barycentric coordinates via glueMesh, but find Unew via cutMesh,
% which has the same number of faces, but the seam vertex ids are
% different. That doesn't affect the computation, however! 
Vg = glueMesh.v ;
Fg = glueMesh.f ;
bc = barycentric_coordinates(qV, Vg(Fg(qF, 1), :), ...
    Vg(Fg(qF, 2), :), Vg(Fg(qF, 3), :)) ;
Unew = bc(:, 1) .* UU(FF(qF, 1), :) + ...
     bc(:, 2) .* UU(FF(qF, 2), :) + ...
     bc(:, 3) .* UU(FF(qF, 3), :) ;
 

%% By convention put cp1 and cp2 on lower cut -- this is a bit of a hack
if Unew(cp1, 2) > pi
    Unew(cp1, 2) = Unew(cp1, 2) - 2*pi ;
end
if Unew(cp2, 2) > pi
    Unew(cp2, 2) = Unew(cp2, 2) - 2*pi ;
end
 
%% Find Unew for cutMesh version in pb space so we can find gg,bb
glueMeshnew.u = Unew ;
nVtx_cut = size(cutMeshnew.v, 1) ;
nVtx_glue = size(Unew, 1) ;
try
    assert(nVtx_cut - nVtx_glue == size(cutMeshnew.pathPairs, 1))
catch
    error('The #new vertices should be the same as #vertices in pathPairs')
end
Unew_cut = Unew ;
Unew_cut((nVtx_glue + 1):nVtx_cut, :) = Unew(cutMeshnew.pathPairs(:, 1), :) ;

cutMeshnew.u = Unew_cut ;

% check remeshing -- glueMesh
if preview
    plot(Unew(:, 1), Unew(:, 2), '.'); hold on
    scatter(Unew(:, 1), Unew(:, 2), 5, 1:size(Unew,1), 'filled')
    trisurf(triangulation(Fnew, cat(2, Unew, 0*Unew(:, 1))))
    pause(0.2)
    
    clf
    plot(Unew_cut(:, 1), Unew_cut(:, 2), '.'); hold on
    scatter(Unew(:, 1), Unew(:, 2), 5, 1:size(Unew,1), 'filled')
    trisurf(triangulation(Fnew, cat(2, Unew, 0*Unew(:, 1))))
    plot(Unew_cut(cutBndnew, 1), Unew_cut(cutBndnew, 2), '.-')
    plot(Unew([cp1 cp2], 1), Unew([cp1, cp2], 2), 'o')
    pause(0.2)
end

% Check for indices that are not inside or on boundary 
% These will need to be pushed to their periodic image into the original
% bounding pullback polygon. 
tri = triangulation(FF, VV) ;
bnd0 = tri.freeBoundary() ;
bndPolygon = UU(bnd0, :) ;
inside = inpolygon(Unew(:, 1), Unew(:, 2), bndPolygon(:, 1), bndPolygon(:, 2)) ;
if any(~inside)
    error('handle here')
end
if preview 
    plot(bndPolygon(:, 1), bndPolygon(:, 2), '.-') ; hold on
    plot(Unew(:, 1), Unew(:, 2), '.')
end

%% Rename for ease
Fgnew = glueMeshnew.f ;
Vgnew = glueMeshnew.v ;
Ugnew = glueMeshnew.u ;
Fcnew = cutMeshnew.f ;
Vcnew = cutMeshnew.v ;
Ucnew = cutMeshnew.u ;

%% Start fixing triangles at seam which share periodic vertices
TR = triangulation(Fcnew, Vcnew);
edgeList = TR.edges;
[eL, ~] = calculateEdgeLengthsAndAngles(Fcnew, [0*Ucnew(:, 1), Ucnew(:, 2)]) ;
%%  Check 1: candidates for cutting have long bonds in theta direction but 
% also are very close to previous boundary points
cutCand = edgeList(eL > pi, :) ;
cutCand = unique(cutCand(:)) ;
cutCand = intersect(cutCand, cutBndnew) ;

% plot(Unew(cutCand(:), 1), Unew(cutCand(:), 2), '.')
% idx = pointMatch(UU(bnd0, :), Unew(cutCand(:), :)) ;
% XX = Unew_cut(cutBndnew, :) ;
% YY = Unew_cut(cutCand(:), :) ;
% dy = zeros(length(cutCand(:)), 1) ;
% yidx = zeros(length(cutCand(:)), 1) ;
% for i = 1:size(XX,1)
%     [dy(i), yidx(i)] = min(sqrt((YY(:,1) - XX(i,1)).^2 + (YY(:,2) - XX(i,2)).^2));
% end
% yidx = yidx(dy < 1e-4) ;
% cutCand = cutCand(yidx) ;
% plot(Unew(cutCand(:), 1), Unew(cutCand(:), 2), 'o')

%% CHECK 2: Are the reflected positions actually out of bounds?
cands = unique(cutCand(:)) ;
Ucnew2 = Ucnew ;
cUpper = Ucnew(cands, 2) > pi ;
cLower = Ucnew(cands, 2) < pi ;
cHi = cands(cUpper) ;
cLo = cands(cLower) ;

% Ensure that the reflected position isnt out of bounds
% inDn = inpolygon(Unew(cHi,1), Unew(cHi,2) - 2*pi, ...
%     bndPolygon(:, 1), bndPolygon(:, 2)) ;
% inUp = inpolygon(Unew(cLo,1), Unew(cLo,2) + 2*pi, ...
%     bndPolygon(:, 1), bndPolygon(:, 2)) ;
% cHi = cHi(inDn) ;
% cLo = cLo(inUp) ;
% if preview
%     plot(Unew(cHi, 1), Unew(cHi, 2), 's'); hold on;
%     plot(Unew(cLo, 1), Unew(cLo, 2), '^')
%     pause(1)
% end

%% CHECK 3: are the bulk neighbors of each candidate on the other
% side of the branch cut?
NLcell = TRI2NLcell(Fcnew, Ucnew) ;
bndVtx = TR.freeBoundary() ;

keep = false(size(cHi)) ;
for ii = 1:length(cHi)
    cand = cHi(ii) ;
    % obtain bulk neighbors for this vertex
    neighs = NLcell{cand} ;
    bulkN = setdiff(neighs, unique(bndVtx(:))) ;
    % are bulk neighbors far or near?
    uvTrans = Ucnew(cand, :) ;
    uvTrans(1, 2) = uvTrans(1, 2) - 2*pi ;
    dists = vecnorm(Ucnew(bulkN, :) - Ucnew(cand, :), 2, 2) ;
    dTrans = vecnorm(Ucnew(bulkN, :) - uvTrans, 2, 2) ;
    if all(dists > dTrans)
        keep(ii) = true ;
    elseif any(dists > dTrans)
        error('some bulk neighbors are near and some far! What to do?')
    end
end
if any(keep)
    disp('Some vertices on top branch cut should be on top. Moving...')
end
cHi2 = cHi(keep) ;
% Push particles to the correct side of the branch cut
Ucnew2(cHi2, 2) = Ucnew(cHi2, 2) - 2*pi ;

if preview
    clf
    plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-'); hold on
    trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1))))
    plot(Ucnew2(bndVtx, 1), Ucnew2(bndVtx, 2), '.-')
end

%% repeat for lo
keep = false(size(cLo)) ;
for ii = 1:length(cLo)
    cand = cLo(ii) ;
    % obtain bulk neighbors for this vertex
    neighs = NLcell{cand} ;
    bulkN = setdiff(neighs, unique(bndVtx(:))) ;
    % are bulk neighbors far or near?
    uvTrans = Ucnew2(cand, :) ;
    uvTrans(1, 2) = uvTrans(1, 2) + 2*pi ;
    dists = vecnorm(Ucnew2(bulkN, :) - Ucnew2(cand, :), 2, 2) ;
    dTrans = vecnorm(Ucnew2(bulkN, :) - uvTrans, 2, 2) ;
    if all(dists > dTrans)
        keep(ii) = true ;
        disp(['Flipping vertex: ' num2str(cand)])
    elseif any(dists > dTrans)
        error('some bulk neighbors are near and some far! What to do?')
    end

    % Check the polygon
    if preview
        clf; hold on
        plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-')
        trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1)))) 
        plot(Ucnew2(cand, 1), Ucnew2(cand, 2), 's')
        plot(Ucnew2(bulkN, 1), Ucnew2(bulkN, 2), '.')
        pause(0.1)
    end
end
if any(keep)
    disp('Some vertices on bottom branch cut should be on top. Moving...')
end
cLo2 = cLo(keep) ;
% Push particles to the correct side of the branch cut
Ucnew2(cLo2, 2) = Ucnew2(cLo2, 2) + 2 * pi ;

if preview
    clf
    hold on;
    trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1))))
    plot(Ucnew2(cHi2, 1), Ucnew2(cHi2, 2), 's')
    plot(Ucnew(cLo2, 1), Ucnew(cLo2, 2), '^')
    pause(1)
end

%% CHECK 4: Remaining edge vertices with long bonds may have only boundary 
% neighbors. Flip those if needed. 

[eL, ~] = calculateEdgeLengthsAndAngles(Fcnew, [0*Ucnew2(:, 1), Ucnew2(:, 2)]) ;
%%  Check 4a: candidates for cutting have long bonds in theta direction but 
% also are very close to previous boundary points
cutCand = edgeList(eL > pi, :) ;
cutCand = intersect(unique(cutCand(:)), unique(bndVtx(:))) ;

%% CHECK 4b: Are the reflected positions actually out of bounds? Which are upper and lower?
cands = unique(cutCand(:)) ;
cUpper = Ucnew2(cands, 2) > pi ;
cLower = Ucnew2(cands, 2) < pi ;
cHi = cands(cUpper) ;
cLo = cands(cLower) ;
% Ensure that the reflected position isnt out of bounds
inDn = inpolygon(Ucnew2(cHi,1), Ucnew2(cHi,2) - 2*pi, ...
    bndPolygon(:, 1), bndPolygon(:, 2)) ;
inUp = inpolygon(Ucnew2(cLo,1), Ucnew2(cLo,2) + 2*pi, ...
    bndPolygon(:, 1), bndPolygon(:, 2)) ;
cHi = cHi(inDn) ;
cLo = cLo(inUp) ;
if preview
    clf; hold on;
    trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1))))
    plot(Ucnew2(cHi, 1), Ucnew2(cHi, 2), 's')
    plot(Ucnew2(cLo, 1), Ucnew2(cLo, 2), '^')
    pause(1)
end

%% CHECK 4c: are the bulk neighbors of each candidate on the other
% side of the branch cut?
NLcell = TRI2NLcell(Fcnew, Ucnew) ;
bndVtx = TR.freeBoundary() ;

keep = false(size(cHi)) ;
for ii = 1:length(cHi)
    cand = cHi(ii) ;
    % obtain bulk neighbors for this vertex
    neighs = NLcell{cand} ;
    % bulkN = setdiff(neighs, unique(bndVtx(:))) ;

    % are edge neighbors far or near?
    uvTrans = Ucnew(cand, :) ;
    uvTrans(1, 2) = uvTrans(1, 2) - 2*pi ;
    dists = vecnorm(Ucnew(neighs, :) - Ucnew(cand, :), 2, 2) ;
    dTrans = vecnorm(Ucnew(neighs, :) - uvTrans, 2, 2) ;
    if all(dists > dTrans)
        keep(ii) = true ;
    end

    % Check the polygon
    if preview
        clf; hold on
        plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-')
        trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1)))) 
        plot(Ucnew2(cand, 1), Ucnew2(cand, 2), 's')
        plot(Ucnew2(neighs, 1), Ucnew2(neighs, 2), '.')
        pause(0.1)
    end
end
if any(keep)
    disp('Some vertices on top branch cut should be on bottom. Moving...')
end
cHi2 = cHi(keep) ;
% Push particles to the correct side of the branch cut
Ucnew2(cHi2, 2) = Ucnew(cHi2, 2) - 2*pi ;

if preview
    clf
    plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-'); hold on
    trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1))))
end

%% repeat for lo
keep = false(size(cLo)) ;
badVertices = [] ;
for ii = 1:length(cLo)
    cand = cLo(ii) ;
    % obtain bulk neighbors for this vertex
    neighs = NLcell{cand} ;
    % bulkN = setdiff(neighs, unique(bndVtx(:))) ;

    % are neighbors far or near?
    uvTrans = Ucnew2(cand, :) ;
    uvTrans(1, 2) = uvTrans(1, 2) + 2*pi ;
    dists = vecnorm(Ucnew2(neighs, :) - Ucnew2(cand, :), 2, 2) ;
    dTrans = vecnorm(Ucnew2(neighs, :) - uvTrans, 2, 2) ;
    if all(dists > dTrans)
        keep(ii) = true ;
        disp(['Flipping vertex: ' num2str(cand)])
    elseif any(dists > dTrans)
        badVertices = [badVertices, cand ] ;
        disp('WARNING: some bulk/edge neighbors are near and some far! Ignoring this tri.')
    end
    
    % Check the polygon
    if preview
        clf; hold on
        plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-')
        trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1)))) 
        plot(Ucnew2(cand, 1), Ucnew2(cand, 2), 's')
        plot(Ucnew2(neighs, 1), Ucnew2(neighs, 2), '.')
        pause(0.1)
    end
end
if any(keep)
    disp('Some vertices on bottom branch cut should be on top. Moving...')
end
cLo2 = cLo(keep) ;
% Push particles to the correct side of the branch cut
Ucnew2(cLo2, 2) = Ucnew2(cLo2, 2) + 2 * pi ;

if preview
    clf
    hold on;
    trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1))))
    plot(Ucnew2(cHi2, 1), Ucnew2(cHi2, 2), 's')
    pause(1)
end

%% DONE -- outputs
Ucnew = Ucnew2 ;
cutMeshnew.u = Ucnew ;

%% check remeshing
if preview
    clf
    plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-'); hold on
    trisurf(triangulation(Fcnew, cat(2, Ucnew2, 0*Ucnew2(:, 1))))
end

disp('done remeshing with seam')