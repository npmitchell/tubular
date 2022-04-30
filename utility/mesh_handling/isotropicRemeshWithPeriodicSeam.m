function [Fnew, Vnew, Unew, Fnnew, Vnnew, pathPairs, bc, qF] = ...
    isotropicRemeshWithPeriodicSeam(...
    FF, VV, UU, targetEdgeLength, numIterations, options) 
%[Fnew, Vnew, Unew, Fnnew, Vnnew] = isotropicRemeshWithPeriodicSeam(...
%   FF, VV, UU, targetEdgeLength, numIterations, options) 
%
% Assumes CGAL_Code/isotropic_remeshing.cpp is in path and compiled.
% Assumes Y is periodic with translation of 2*pi as in a Ricci mesh / angle
%
% Parameters
% ----------
% FF : #faces x 3 int array
%   face connectivity list before remeshing
% VV : #vertices x 3 float array
%   vertices of cutMesh before remeshing. Cut seam has two copies of each
%   vertex. 
% UU : #vertices x 2 float array
%   pullback 2D coordinates of cutMesh (topological cylinder with cut seam)
%   before remeshing. Note that input state is topological disk (rectangle)
%   with implied identification along the cutPath.
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
% [Fnew, Vnew, Unew, Fnnew, Vnnew, 
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

epsVtx = 1e-5 ;
preview = false ;
if nargin > 5
    if isfield(options, 'epsVtx')
        epsVtx = options.epsVtx;
    end
    if isfield(options, 'preview')
        preview = options.preview;
    end
end


[Fnew, Vnew, Fnnew, Vnnew] = ...
       isotropic_remeshing(FF, VV, targetEdgeLength, numIterations) ;
[sqrD, qF, qV] = point_mesh_squared_distance(Vnew, VV, FF ) ;
assert(all(vecnorm(qV - Vnew, 2, 2) < epsVtx))
Vnew = qV ;

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

% initial remeshing -- still has periodic boundary
bc = barycentric_coordinates(qV, VV(FF(qF, 1), :), ...
    VV(FF(qF, 2), :), VV(FF(qF, 3), :)) ;
Unew = bc(:, 1) .* UU(FF(qF, 1), :) + ...
     bc(:, 2) .* UU(FF(qF, 2), :) + ...
     bc(:, 3) .* UU(FF(qF, 3), :) ;
% check remeshing
if preview
    plot(Unew(:, 1), Unew(:, 2), '.'); hold on
    scatter(Unew(:, 1), Unew(:, 2), 5, 1:size(Unew,1), 'filled')
    trisurf(triangulation(Fnew, cat(2, Unew, 0*Unew(:, 1))))
    pause(0.2)
end

% Check for indices that are not inside or on boundary -- should be none
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

%% Start fixing triangles at seam which share periodic vertices
TR = triangulation(Fnew, qV);
edgeList = TR.edges;
[eL, ~] = calculateEdgeLengthsAndAngles(Fnew, [0*Unew(:, 1), Unew(:, 2)]) ;
%%  Check 1: candidates for cutting have long bonds in theta direction but 
% also are very close to previous boundary points
cutCand = edgeList(eL > pi, :) ;
% plot(Unew(cutCand(:), 1), Unew(cutCand(:), 2), '.')
% idx = pointMatch(UU(bnd0, :), Unew(cutCand(:), :)) ;
XX = UU(bnd0, :) ;
YY = Unew(cutCand(:), :) ;
dy = zeros(length(cutCand(:)), 1) ;
yidx = zeros(length(cutCand(:)), 1) ;
for i = 1:size(XX,1)
    [dy(i), yidx(i)] = min(sqrt((YY(:,1) - XX(i,1)).^2 + (YY(:,2) - XX(i,2)).^2));
end
yidx = yidx(dy < 1e-4) ;
cutCand = cutCand(yidx) ;
% plot(Unew(cutCand(:), 1), Unew(cutCand(:), 2), 'o')

%% CHECK 2: Are the reflected positions actually out of bounds?
cands = unique(cutCand(:)) ;
Unew2 = Unew ;
cUpper = Unew(cands, 2) > pi ;
cLower = Unew(cands, 2) < pi ;
cHi = cands(cUpper) ;
cLo = cands(cLower) ;
% Ensure that the reflected position isnt out of bounds
inDn = inpolygon(Unew(cHi,1), Unew(cHi,2) - 2*pi, ...
    bndPolygon(:, 1), bndPolygon(:, 2)) ;
inUp = inpolygon(Unew(cLo,1), Unew(cLo,2) + 2*pi, ...
    bndPolygon(:, 1), bndPolygon(:, 2)) ;
cHi = cHi(inDn) ;
cLo = cLo(inUp) ;
if preview
    plot(Unew(cHi, 1), Unew(cHi, 2), 's')
    plot(Unew(cLo, 1), Unew(cLo, 2), '^')
    pause(1)
end

%% CHECK 3: are the bulk neighbors of each candidate on the other
% side of the branch cut?
NLcell = TRI2NLcell(Fnew, Unew) ;
bndVtx = TR.freeBoundary() ;

keep = false(size(cHi)) ;
for ii = 1:length(cHi)
    cand = cHi(ii) ;
    % obtain bulk neighbors for this vertex
    neighs = NLcell{cand} ;
    bulkN = setdiff(neighs, unique(bndVtx(:))) ;
    % are bulk neighbors far or near?
    uvTrans = Unew(cand, :) ;
    uvTrans(1, 2) = uvTrans(1, 2) - 2*pi ;
    dists = vecnorm(Unew(bulkN, :) - Unew(cand, :), 2, 2) ;
    dTrans = vecnorm(Unew(bulkN, :) - uvTrans, 2, 2) ;
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
Unew2(cHi2, 2) = Unew(cHi2, 2) - 2*pi ;

if preview
    clf
    plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-'); hold on
    trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1))))
end

%% repeat for lo
keep = false(size(cLo)) ;
for ii = 1:length(cLo)
    cand = cLo(ii) ;
    % obtain bulk neighbors for this vertex
    neighs = NLcell{cand} ;
    bulkN = setdiff(neighs, unique(bndVtx(:))) ;
    % are bulk neighbors far or near?
    uvTrans = Unew2(cand, :) ;
    uvTrans(1, 2) = uvTrans(1, 2) + 2*pi ;
    dists = vecnorm(Unew2(bulkN, :) - Unew2(cand, :), 2, 2) ;
    dTrans = vecnorm(Unew2(bulkN, :) - uvTrans, 2, 2) ;
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
        trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1)))) 
        plot(Unew2(cand, 1), Unew2(cand, 2), 's')
        plot(Unew2(bulkN, 1), Unew2(bulkN, 2), '.')
        pause(0.1)
    end
end
if any(keep)
    disp('Some vertices on bottom branch cut should be on top. Moving...')
end
cLo2 = cLo(keep) ;
% Push particles to the correct side of the branch cut
Unew2(cLo2, 2) = Unew2(cLo2, 2) + 2 * pi ;

if preview
    clf
    hold on;
    trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1))))
    plot(Unew2(cHi2, 1), Unew2(cHi2, 2), 's')
    plot(Unew(cLo2, 1), Unew(cLo2, 2), '^')
    pause(1)
end

%% CHECK 4: Remaining edge vertices with long bonds may have only boundary 
% neighbors. Flip those if needed. 

[eL, ~] = calculateEdgeLengthsAndAngles(Fnew, [0*Unew2(:, 1), Unew2(:, 2)]) ;
%%  Check 4a: candidates for cutting have long bonds in theta direction but 
% also are very close to previous boundary points
cutCand = edgeList(eL > pi, :) ;
% plot(Unew(cutCand(:), 1), Unew(cutCand(:), 2), '.')
% idx = pointMatch(UU(bnd0, :), Unew(cutCand(:), :)) ;
XX = UU(bnd0, :) ;
YY = Unew2(cutCand(:), :) ;
dy = zeros(length(cutCand(:)), 1) ;
yidx = zeros(length(cutCand(:)), 1) ;
for i = 1:size(XX,1)
    [dy(i), yidx(i)] = min(sqrt((YY(:,1) - XX(i,1)).^2 + (YY(:,2) - XX(i,2)).^2));
end
yidx = yidx(dy < epsVtx) ;
cutCand = cutCand(yidx) ;
% plot(Unew(cutCand(:), 1), Unew(cutCand(:), 2), 'o')

%% CHECK 4b: Are the reflected positions actually out of bounds? Which are upper and lower?
cands = unique(cutCand(:)) ;
cUpper = Unew2(cands, 2) > pi ;
cLower = Unew2(cands, 2) < pi ;
cHi = cands(cUpper) ;
cLo = cands(cLower) ;
% Ensure that the reflected position isnt out of bounds
inDn = inpolygon(Unew2(cHi,1), Unew2(cHi,2) - 2*pi, ...
    bndPolygon(:, 1), bndPolygon(:, 2)) ;
inUp = inpolygon(Unew2(cLo,1), Unew2(cLo,2) + 2*pi, ...
    bndPolygon(:, 1), bndPolygon(:, 2)) ;
cHi = cHi(inDn) ;
cLo = cLo(inUp) ;
if preview
    clf; hold on;
    trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1))))
    plot(Unew2(cHi, 1), Unew2(cHi, 2), 's')
    plot(Unew2(cLo, 1), Unew2(cLo, 2), '^')
    pause(1)
end

%% CHECK 4c: are the bulk neighbors of each candidate on the other
% side of the branch cut?
NLcell = TRI2NLcell(Fnew, Unew) ;
bndVtx = TR.freeBoundary() ;

keep = false(size(cHi)) ;
for ii = 1:length(cHi)
    cand = cHi(ii) ;
    % obtain bulk neighbors for this vertex
    neighs = NLcell{cand} ;
    % bulkN = setdiff(neighs, unique(bndVtx(:))) ;

    % are edge neighbors far or near?
    uvTrans = Unew(cand, :) ;
    uvTrans(1, 2) = uvTrans(1, 2) - 2*pi ;
    dists = vecnorm(Unew(neighs, :) - Unew(cand, :), 2, 2) ;
    dTrans = vecnorm(Unew(neighs, :) - uvTrans, 2, 2) ;
    if all(dists > dTrans)
        keep(ii) = true ;
    end

    % Check the polygon
    if preview
        clf; hold on
        plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-')
        trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1)))) 
        plot(Unew2(cand, 1), Unew2(cand, 2), 's')
        plot(Unew2(neighs, 1), Unew2(neighs, 2), '.')
        pause(0.1)
    end
end
if any(keep)
    disp('Some vertices on top branch cut should be on bottom. Moving...')
end
cHi2 = cHi(keep) ;
% Push particles to the correct side of the branch cut
Unew2(cHi2, 2) = Unew(cHi2, 2) - 2*pi ;

if preview
    clf
    plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-'); hold on
    trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1))))
end

%% repeat for lo
keep = false(size(cLo)) ;
for ii = 1:length(cLo)
    cand = cLo(ii) ;
    % obtain bulk neighbors for this vertex
    neighs = NLcell{cand} ;
    % bulkN = setdiff(neighs, unique(bndVtx(:))) ;

    % are neighbors far or near?
    uvTrans = Unew2(cand, :) ;
    uvTrans(1, 2) = uvTrans(1, 2) + 2*pi ;
    dists = vecnorm(Unew2(neighs, :) - Unew2(cand, :), 2, 2) ;
    dTrans = vecnorm(Unew2(neighs, :) - uvTrans, 2, 2) ;
    if all(dists > dTrans)
        keep(ii) = true ;
        disp(['Flipping vertex: ' num2str(cand)])
    elseif any(dists > dTrans)
        error('some bulk/edge neighbors are near and some far! What to do?')
    end
    
    % Check the polygon
    if preview
        clf; hold on
        plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-')
        trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1)))) 
        plot(Unew2(cand, 1), Unew2(cand, 2), 's')
        plot(Unew2(neighs, 1), Unew2(neighs, 2), '.')
        pause(0.1)
    end
end
if any(keep)
    disp('Some vertices on bottom branch cut should be on top. Moving...')
end
cLo2 = cLo(keep) ;
% Push particles to the correct side of the branch cut
Unew2(cLo2, 2) = Unew2(cLo2, 2) + 2 * pi ;

if preview
    clf
    hold on;
    trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1))))
    plot(Unew2(cHi2, 1), Unew2(cHi2, 2), 's')
    pause(1)
end

%% DONE
Unew = Unew2 ;


%% check remeshing
if preview
    clf
    plot(bndPolygon(:, 1), bndPolygon(:, 2), 'o-'); hold on
    trisurf(triangulation(Fnew, cat(2, Unew2, 0*Unew2(:, 1))))
end

%% If asked for, supply cutMesh
if nargout > 5
    % For now assume that freeBoundary traverses counterclockwise and is
    % singly connected
    
    candidates = bndVtx ;
    left = find(Unew(candidates, 1) < (min(Unew(candidates, 1)) + epsVtx)) ;
    right = find(Unew(candidates, 1) > (max(Unew(candidates, 1)) - epsVtx)) ;
    [~, leftBottom ] = min(Unew(candidates(left), 2)) ;
    start1 = candidates(left(leftBottom)) ;
    [~, rightBottom ] = min(Unew(candidates(right), 2)) ;
    end1 = candidates(right(rightBottom)) ;
    [~, leftTop ] = max(Unew(candidates(left), 2)) ;
    start2 = candidates(left(leftTop)) ;
    [~, rightTop ] = max(Unew(candidates(right), 2)) ;
    end2 = candidates(right(rightTop)) ;
    
    % Start with s1, wind to e1 -- index into double cover of boundary
    fullPath = [bndVtx(:, 1); bndVtx(:, 1)] ; 
    s1ID = find(bndVtx(:, 1) == start1) ;
    e1ID = find(fullPath == end1) ;
    % end with first e1 that is larger than s1
    e1ID = min(e1ID(e1ID > s1ID)) ;
    
    % Repeat for the other side
    fullPT = flipud(fullPath) ;
    s2ID = find(fullPT == start2, 1) ;
    e2ID = find(fullPT == end2) ;
    % end with first e2 that is larger than s2
    e2ID = min(e2ID(e2ID > s2ID)) ;
    
    % Cat the two path pairs
    pathPairs = [fullPath(s1ID:e1ID), fullPT(s2ID:e2ID) ] ;
    
    if preview
        for ii = 1:length(bndVtx)
            plot(Unew(bndVtx(1:ii), 1), Unew(bndVtx(1:ii), 2), '.-')
            pause(0.00001)
        end
        hold on;
        plot(Unew(pathPairs(:, 1), 1), Unew(pathPairs(:, 1), 2), '.')
        plot(Unew(pathPairs(:, 2), 1), Unew(pathPairs(:, 2), 2), '.')
        plot(Unew(fullPath(s1ID), 1), Unew(fullPath(s1ID), 2), 'v')
        plot(Unew(fullPath(e1ID), 1), Unew(fullPath(e1ID), 2), 'o')
        plot(Unew(fullPT(s2ID), 1), Unew(fullPT(s2ID), 2), '^')
        plot(Unew(fullPT(e2ID), 1), Unew(fullPT(e2ID), 2), 's')
        pause(1e-5)
    end
end

disp('done remeshing with seam')