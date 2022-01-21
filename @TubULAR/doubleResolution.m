function [cutMesh2x, cutMesh2xC] = doubleResolution(cutMesh, preview)
% doubleResolution(cutMesh, preview)
%   Double the resolution of a cutMesh
% 
% Parameters
% ----------
% cutMesh : struct with fields 
%   nU :
%   nV :
%   f :
%   v :
%   u : 
%
% NPMitchell 2020

if nargin < 2
    preview = false;
end

% Double resolution in uv
uv0 = cutMesh.u ;
v3 = cutMesh.v ;
vn0 = cutMesh.vn ;
nU = cutMesh.nU ;
nV = cutMesh.nV ;

nU2 = nU * 2 - 1;
nV2 = nV * 2 - 1;

% Double resolution in V
unew = zeros(nU, nV2, 2) ;
v3new = zeros(nU, nV2, 3) ;
ugrid = reshape(uv0, [nU, nV, 2]) ;
v3grid = reshape(v3, [nU, nV, 3]) ;
for qq = 1:nV - 1
    % Duplicate column of u=const in 2d
    unew(:, 2*qq-1, :) = ugrid(:, qq, :) ;
    unew(:, 2 * qq, :) = 0.5 * (ugrid(:, qq, :) ...
                              + ugrid(:, qq + 1, :)) ;
    % Duplicate column of u=const in 3d
    v3new(:, 2*qq-1, :) = v3grid(:, qq, :) ;
    v3new(:, 2 * qq, :) = 0.5 * (v3grid(:, qq, :) ...
                              + v3grid(:, qq + 1, :)) ;
end
unew(:, nV2, :) = ugrid(:, nV, :) ;
v3new(:, nV2, :) = v3grid(:, nV, :) ;

% Double resolution in U
uv = zeros(nU2, nV2, 2) ;
v3d = zeros(nU2, nV2, 3) ;
for qq = 1:nU - 1
    uv(2*qq-1, :, :) = unew(qq, :, :) ;
    uv(2 * qq, :, :) = 0.5 * (unew(qq, :, :) ...
                               + unew(qq + 1, :, :)) ;
    % Duplicate column of v=const in 3d
    v3d(2*qq-1, :, :) = v3new(qq, :, :) ;
    v3d(2 * qq, :, :) = 0.5 * (v3new(qq, :, :) ...
                              + v3new(qq + 1, :, :)) ;
end
uv(nU2, :, :) = unew(nU, :, :) ;
v3d(nU2, :, :) = v3new(nU, :, :) ;

% Check it
if preview
    % check it
    imagesc(squeeze(v3new(:, :, 2)))
    waitfor(gcf)
    clf
    hold on;
    for qq = 1:length(v3d)
        plot3(v3d(qq, :, 1), v3d(qq, :, 2), v3d(qq, :, 3), '.')
        pause(0.0001)
    end
end

% Reshape uv and v3d
uv = reshape(uv, [nU2 * nV2, 2]) ;
v3d = reshape(v3d, [nU2 * nV2, 3]) ;

% Check it
if preview
    clf
    hold on;
    nn = 50 ;
    for qq = 1:nn:length(v3d)
        plot3(v3d(qq:qq+nn, 1), v3d(qq:qq+nn, 2), v3d(qq:qq+nn, 3), '.')
        pause(0.00001)
    end
end

% Output variables
faces = defineFacesRectilinearGrid(uv, nU2, nV2) ;
% Store in output cutMesh struct
cutMesh2x = struct() ;
cutMesh2x.v = v3d ;
cutMesh2x.u = uv ;
cutMesh2x.f = faces ;
cutMesh2x.nU = nU2 ;
cutMesh2x.nV = nV2 ;
cutMesh2x.pathPairs = [1:nU2; (nV2-1)*nU2 + (1:nU2)]' ;

cutMesh2xC = glueCylinderCutMeshSeam(cutMesh2x) ;
vn = zeros(size(cutMesh2x.v)) ;
vn(1:nU2*(nV2-1), :) = cutMesh2xC.vn ;
vn(nU2*(nV2-1)+1:nU2*nV2, :) = cutMesh2xC.vn(1:nU2, :) ;
cutMesh2x.vn = vn ;

%%% ALTERNATE APPROACH
% % Create vertex normals via interpolation
% % push boundaries inward slightly by epsilon
% leftId = uv(:, 1) < min(uv0(:, 1)) + eps ;
% rightId = uv(:, 1) > max(uv0(:, 1)) - eps ;
% uv(leftId, 1) = uv(leftId, 1) + eps ; 
% uv(rightId, 1) = uv(rightId, 1) + eps ;
% [ TF, TV2D, ~, TVN3D ] = tileAnnularCutMesh(cutMesh, [1, 1]) ;
% vn = interpolate2Dpts_3Dmesh(TF, TV2D, TVN3D, uv) ;
% vn = vn ./ vecnorm(vn, 2, 2) ;
% 
% dmyk = 0;
% while any(isnan(vn(:))) && dmyk < 1000 
%     % Fix bad normals
%     bad = find(isnan(vn(:, 1))) ;
%     disp(['found ', num2str(length(bad)), ...
%         ' bad vertex normals, attempting to fix...'])
%     jitter = 1e-14 * rand([size(bad, 1), 2]) ;
%     size(vn(bad, :))
%     vn(bad, :) = interpolate2Dpts_3Dmesh(TF, TV2D, TVN3D, uv(bad, :)+jitter) ;
%     dmyk = dmyk + 1 ;
% end
% if any(isnan(vn(:)))
%     error(['bad vertex normals. ', ...
%         'Cannot use cutMesh normals to define ', ...
%         'double Resolution cutMesh.'])
% end
% 
% % Check normals
% vn2 = per_vertex_normals(v3d, faces, 'Weighting', 'angle') ;
% distant = find(abs(vn(:) - vn2(:)) > 0.25) ;
% [distant, colind] = ind2sub(size(vn), distant) ;
% plot(vn(:), vn2(:), '.')
% waitfor(gcf)
% 
% % Check it
% % % bad indices have NaNs
% umax = max(uv(:, 1)) ;
% bad = find(isnan(vn(:, 1))) ;
if preview
    clf
    trisurf(cutMesh.f, cutMesh.v(:, 1), cutMesh.v(:, 2),...
        cutMesh.v(:, 3), 'edgecolor', 'none') ; 
    hold on; 
    % plot3(v3d(distant, 1), v3d(distant, 2), v3d(distant, 3), 'o')
    quiver3(v3d(:, 1), v3d(:, 2), v3d(:, 3), vn(:, 1), vn(:, 2), ...
        vn(:, 3), 1, 'r')
    % quiver3(v3d(:, 1), v3d(:, 2), v3d(:, 3), vn2(:, 1), vn2(:, 2), ...
    %     vn2(:, 3), 1, 'y')
    axis equal
    waitfor(gcf)
end
% % % plot(uv(:, 1)/umax, uv(:, 2), '.')
% % trisurf(TF, TV2D(:, 1)/umax, TV2D(:, 2), 0*TV2D(:, 1)) ; 
% % plot(uv(bad, 1)/umax, uv(bad, 2), 'o')

if preview
    trimesh(faces, v3d(:, 1), v3d(:, 2), v3d(:, 3), ...
        v3d(:, 1), 'Edgecolor', 'k', 'FaceColor', 'Interp')
    axis equal
    title('Preview double-resolution cutMesh')
end