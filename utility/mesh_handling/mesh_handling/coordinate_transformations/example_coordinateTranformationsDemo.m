%% Jacobian testing 
% to figure out if covariant or contravariant vectors are assumed as inputs
% This shows that jacobian functions like pushVectorField2Dto3DMesh act on
% contravariant vectors, NOT on covariant vectors

nU = 100 ;
nV = nU ;
ss = linspace(-2, 2, nV) ;
[yy, xx] = meshgrid(ss, ss) ;
sigma = 1/sqrt(2) ;
bump = exp(-(xx.^2 + yy.^2) / (2 * sigma.^2)) ;
x = xx(:) ;
y = yy(:) ;
vtx = [x, y, bump(:)] ;
faces = defineFacesRectilinearGrid(vtx, nU, nV) ;

% Check the mesh
trisurf(triangulation(faces, vtx), 'edgecolor', 'none')

%% Vector field
% tangential contravariant vector field
v2d = [x,y,0*x] ;
bc = barycenter(v2d, faces) ;
xf = bc(:, 1) ;
yf = bc(:, 2) ;
bc3d = barycenter(vtx, faces) ;
zf = bc3d(:, 3) ;

% Face-based vector field
vv = [xf.^2, yf.^2] ;

% The 3d vector field is then
vv_3d = pushVectorField2Dto3DMesh(vv, v2d, vtx, faces, 1:size(faces, 1)) ;

% Note that ex, ey are 
exz = - xf .* exp(-xf.^2 - yf.^2 / (2*sigma.^2)) / sigma.^2 ;
eyz = - yf .* exp(-xf.^2 - yf.^2 / (2*sigma.^2)) / sigma.^2 ;
ex = [ones(size(exz)), zeros(size(exz)), exz] ;
ey = [zeros(size(exz)), ones(size(exz)), eyz] ;

% Analytic push forward vector field
vv_3da = vv(:, 1) .* ex + vv(:, 2) .* ey ;

% Analytic metric
expfac = exp((-xf.^2 - yf.^2) / sigma.^2) / sigma.^4 ;
gg = [ 1 + expfac .* xf.^2, ...
    expfac .* xf .* yf, ...
    expfac .* xf .* yf, ...
    1 + expfac .* yf.^2 ] ;

% Compare the pushed covariant vector
vv_down = [gg(:, 1) .* vv(:, 1) + gg(:, 2) .* vv(:, 2),  ...
    gg(:, 3) .* vv(:, 1) + gg(:, 4) .* vv(:, 2)] ;

vv_down_3d = pushVectorField2Dto3DMesh(vv_down, v2d, vtx, faces, ...
    1:size(faces, 1)) ;


%% CONTRAVARIANT
%% Plot results -- raw difference -- CONTRAVARIANT
clf
h1 = plot(vv_3d(:, 1), vv_3d(:, 1) - vv_3da(:, 1), '.') ;
hold on;
h2 = plot(vv_3d(:, 2), vv_3d(:, 2) - vv_3da(:, 2), 'o') ;
h3 = plot(vv_3d(:, 3), vv_3d(:, 3) - vv_3da(:, 3), 's') ;
legend('\Delta v_x', '\Delta v_y', '\Delta v_z')
title('Raw difference in vectors in 3d - contravariant')

%% Plot results -- relative difference -- CONTRAVARIANT
clf
mags = vecnorm(vv_3d, 2, 2) ;
diffs = vv_3d - vv_3da ;
h1 = plot(vv_3d(:, 1), diffs(:, 1) ./ mags, '.') ;
hold on;
h2 = plot(vv_3d(:, 2), diffs(:, 2)./ mags, 'o') ;
h3 = plot(vv_3d(:, 3), diffs(:, 3)./ mags, 's') ;
legend('\Delta v_x/ |v|', '\Delta v_y/ |v|', '\Delta v_z / |v|')
title('Relative difference in vectors in 3d - contravariant')

%% Quiver in 3d -- CONTRAVARIANT
clf
quiver3(xf, yf, zf, vv_3d(:, 1), vv_3d(:, 2), vv_3d(:, 3), 1)
hold on;
quiver3(xf, yf, zf, vv_3da(:, 1), vv_3da(:, 2), vv_3da(:, 3), 1)
title('Relative difference in vectors in 3d - contravariant')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COVARIANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot results -- raw difference -- COVARIANT
clf
h1 = plot(vv_down_3d(:, 1), vv_down_3d(:, 1) - vv_3da(:, 1), '.') ;
hold on;
h2 = plot(vv_down_3d(:, 2), vv_down_3d(:, 2) - vv_3da(:, 2), 'o') ;
h3 = plot(vv_down_3d(:, 3), vv_down_3d(:, 3) - vv_3da(:, 3), 's') ;
legend('\Delta v_x', '\Delta v_y', '\Delta v_z')
title('Raw difference in vectors in 3d')

%% Plot results -- relative difference -- COVARIANT
clf
mags = vecnorm(vv_down_3d, 2, 2) ;
diffs = vv_down_3d - vv_3da ;
h1 = plot(vv_down_3d(:, 1), diffs(:, 1) ./ mags, '.') ;
hold on;
h2 = plot(vv_down_3d(:, 2), diffs(:, 2)./ mags, 'o') ;
h3 = plot(vv_down_3d(:, 3), diffs(:, 3)./ mags, 's') ;
legend('\Delta v_x/ |v|', '\Delta v_y/ |v|', '\Delta v_z / |v|')
title('Relative difference in vectors in 3d')

%% Quiver in 3d -- COVARIANT
clf
quiver3(xf, yf, zf, vv_down_3d(:, 1), vv_down_3d(:, 2), vv_down_3d(:, 3), 1)
hold on;
quiver3(xf, yf, zf, vv_3da(:, 1), vv_3da(:, 2), vv_3da(:, 3), 1)

