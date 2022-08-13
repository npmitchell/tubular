function [U, cornerIDx] = flattenSquare( F, V, cornerIDx, LType )
%FLATTENSQUARE This function constructs a conformal parameterization of a
%topological disk in the unit square [0 1] X [0 1]
%
%   INPUT PARAMETERS:
%
%       - F:            #Fx3 face connectivity list. NOTE: Assumes input
%                       mesh faces are CCW oriented
%
%       - V:            #Vx3 3D vertex coordinate list
%
%       - cornerIDx:    List of corner vertex IDs in CCW order from (0,1),
%                       i.e. [ (c1) (c2) (c3) (c4) ] (see below)
%
%       - LType:        The type of mesh Laplacian operator used to
%                       construct the linear system solved to find the
%                       embedding coordinates.  This method currently
%                       supports the 'Dirichlet' Laplacian and the 'Mean
%                       Value Coordinates' operator
%   
%   OUTPUT PARAMETERS:
%
%       - U:            #Vx2 2D vertex coordinate list
%
% by Dillon Cislo 2022/06/30

%==========================================================================
% THE GEOMETRY OF THE OUTPUT:
%
%           (4)
%    (c1)--------(c4)
%     |            |
%     |            |
% (1) |            | (3)
%     |            |
%     |            |
%    (c2)--------(c3)
%           (2)
%
% The output region in the (u,v) plane is the domain: [0 1] X [0 1]
%
% FINDING THE EMBEDDING:
%
% The embedding coordinates are found as the solution to the linear system
%
%   L u = 0 subject to A u = b
%
% Where L is some mesh Laplacian operator and u is a vector corresponding
% to the (u,v) coordinates of each vertex.  (A,b) is a linear set of
% equality constraints that specify the boundary conditions of the problem.
%
%==========================================================================

if (nargin < 1), error('Please supply face connectivity list'); end
if (nargin < 2), error('Please supply 3D vertex coordinates'); end
if (nargin < 3), cornerIDx = []; end
if (nargin < 4), LType = 'Dirichlet'; end

validateattributes(V, {'numeric'}, {'2d', 'ncols', 3, 'finite', 'real'});
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'positive', ...
    'finite', 'real', '<=', size(V,1)});
validateattributes(LType, {'char'}, {'vector'});

TR3D = triangulation(F, V);
E = TR3D.edges; % #Ex2 edge connectivity list

eulerChi = size(F,1) + size(V,1) - size(E,1);
assert(eulerChi == 1, 'Input mesh is NOT a topological disk');

bdyIDx = TR3D.freeBoundary; % #BEx2 sorted list of boundary edges

if ~isempty(cornerIDx)
    
    validateattributes(cornerIDx, {'numeric'}, ...
        {'vector', 'numel', 4, 'positive', 'integer', 'finite', 'real'});
    assert(all(ismember(cornerIDx, bdyIDx(:))), ...
        'Corner vertices must lie on the mesh boundary');
    
else
    
    % Boundary edge lengths
    bdyL = V(bdyIDx(:,2), :) - V(bdyIDx(:,1), :);
    bdyL = sqrt(sum(bdyL.^2, 2));
    
    cumBdyL = cumsum(bdyL); % Cumulative boundary length
    totalBdyL = sum(bdyL); % Total boundary length
    
    % First vertex is arbitrarily picked to be the first vertex in the
    % boundary vertex list. Subsquent vertices try to divide the boundary
    % into approximately equal length segments
    cornerIDx = [bdyIDx(1); ...
        bdyIDx(knnsearch(cumBdyL./totalBdyL, [0.25; 0.5; 0.75]), 2)];
    
end

%==========================================================================
% Partition into segments corresponding to the sides of the square domain
%==========================================================================

% Just retain the first vertex in each boundary edge
bdyIDx = bdyIDx(:,1);

% Shift the list so that the start point is in the first position
ind = find( bdyIDx == cornerIDx(1) );
bdyIDx = bdyIDx( [ ind:end, 1:ind-1 ] );

% Ensure that the list is consistently ordered (see diagram above)
if find( bdyIDx == cornerIDx(3) ) < find( bdyIDx == cornerIDx(2) )
    bdyIDx = [ bdyIDx(1); flipud( bdyIDx(2:end) ) ];
end

% Find the location of the corner vertices in the updated list
ind = [ find( bdyIDx == cornerIDx(1) ); find( bdyIDx == cornerIDx(2) );
    find( bdyIDx == cornerIDx(3) ); find( bdyIDx == cornerIDx(4) ) ];

goodCorners = (ind(1) == 1) && (ind(1) < ind(2)) && ...
    (ind(2) < ind(3)) && (ind(3) < ind(4));
assert(goodCorners, 'Invalid boundary segment division');

% Partition the boundary into segments
bdyPaths = cell(4,1);
for i = 1:length(ind)-1
    bdyPaths{i} = bdyIDx( ind(i):ind(i+1) );
end
bdyPaths{end} = [ bdyIDx(ind(end):end); bdyIDx(1) ];

%==========================================================================
% Construct the linear constraint sub-system
%==========================================================================

A = sparse(0, 2 * size(V,1));
b = [];

%--------------------------------------------------------------------------
% NOTE: I chose here to set hard positional constraints to force the true
% boundary vertices to lie on the lines supporting the vertical edges of
% the square at fixed positions.  Originally, I tried to allow the
% vertices to slide along these lines (see the 'freesquare' functionality
% of Noam Aigerman's 'euclidean_orbifold' package).  For particular meshes
% this led to issues with self-intersections in the embedding.  Perhaps a
% future version can figure out the pathology of these problems and improve
% this method
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Constrain boundary segment (1) to lie on the line (0,1)->(0,0)
%--------------------------------------------------------------------------

% The vertex IDs along segment (1)
segIDx = bdyPaths{1};

v1 = [0,1]; % The start coordinates of the target line
v2 = [0,0]; % The final coordinates of the target line

% Calculate distances between consecutive vertices
d = sqrt( sum( (V(segIDx(1:end-1),:) - V(segIDx(2:end),:)).^2, 2 ) );

% Convert these to a fractional distance along the length of the line
d = [ 0; cumsum(d)/sum(d) ];

% Interpolate to find the (u,v)-coordinates of each vertex
interpXY = [ v1(1) .* (1-d) + v2(1) .* d,  v1(2) .* (1-d) + v2(2) .* d ];

% Set constraints
for i = 1:(length(segIDx)-1)
    A(end+1, 2*segIDx(i)-1) = 1; b = [ b; interpXY(i,1) ];
    A(end+1, 2*segIDx(i)) = 1;   b = [ b; interpXY(i,2) ];
end

%--------------------------------------------------------------------------
% Constrain boundary segment (2) to lie on the line (0,0)->(1,0)
%--------------------------------------------------------------------------

% The vertex IDs along segment (2)
segIDx = bdyPaths{2};

v1 = [0,0]; % The start coordinates of the target line
v2 = [1,0]; % The final coordinates of the target line

% Calculate distances between consecutive vertices
d = sqrt( sum( (V(segIDx(1:end-1),:) - V(segIDx(2:end),:)).^2, 2 ) );

% Convert these to a fractional distance along the length of the line
d = [ 0; cumsum(d)/sum(d) ];

% Interpolate to find the (u,v)-coordinates of each vertex
interpXY = [ v1(1) .* (1-d) + v2(1) .* d,  v1(2) .* (1-d) + v2(2) .* d ];

% Set constraints
for i = 1:(length(segIDx)-1)
    A(end+1, 2*segIDx(i)-1) = 1; b = [ b; interpXY(i,1) ];
    A(end+1, 2*segIDx(i)) = 1;   b = [ b; interpXY(i,2) ];
end

%--------------------------------------------------------------------------
% Constrain boundary segment (3) to lie on the line (1,0)->(1,1)
%--------------------------------------------------------------------------

% The vertex IDs along segment (3)
segIDx = bdyPaths{3};

v1 = [1,0]; % The start coordinates of the target line
v2 = [1,1]; % The final coordinates of the target line

% Calculate distances between consecutive vertices
d = sqrt( sum( (V(segIDx(1:end-1),:) - V(segIDx(2:end),:)).^2, 2 ) );

% Convert these to a fractional distance along the length of the line
d = [ 0; cumsum(d)/sum(d) ];

% Interpolate to find the (u,v)-coordinates of each vertex
interpXY = [ v1(1) .* (1-d) + v2(1) .* d,  v1(2) .* (1-d) + v2(2) .* d ];

% Set constraints
for i = 1:(length(segIDx)-1)
    A(end+1, 2*segIDx(i)-1) = 1; b = [ b; interpXY(i,1) ];
    A(end+1, 2*segIDx(i)) = 1;   b = [ b; interpXY(i,2) ];
end

%--------------------------------------------------------------------------
% Constrain boundary segment (4) to lie on the line (1,1)->(0,1)
%--------------------------------------------------------------------------

% The vertex IDs along segment (4)
segIDx = bdyPaths{4};

v1 = [1,1]; % The start coordinates of the target line
v2 = [0,1]; % The final coordinates of the target line

% Calculate distances between consecutive vertices
d = sqrt( sum( (V(segIDx(1:end-1),:) - V(segIDx(2:end),:)).^2, 2 ) );

% Convert these to a fractional distance along the length of the line
d = [ 0; cumsum(d)/sum(d) ];

% Interpolate to find the (u,v)-coordinates of each vertex
interpXY = [ v1(1) .* (1-d) + v2(1) .* d,  v1(2) .* (1-d) + v2(2) .* d ];

% Set constraints
for i = 1:(length(segIDx)-1)
    A(end+1, 2*segIDx(i)-1) = 1; b = [ b; interpXY(i,1) ];
    A(end+1, 2*segIDx(i)) = 1;   b = [ b; interpXY(i,2) ];
end

%==========================================================================
% Construct the mesh Laplacian operator
%==========================================================================

if strcmpi( LType, 'Dirichlet' ) == 1
    
    L = cotmatrix(V, F);
    
    m = min( min( tril( L, -1 ) ) );
    if m < 0
        
        warning('Mesh is NOT Delaunay!');
        clamp = 1e-2;
        L(L<0) = clamp;
        
        inds = sub2ind( size(L), 1:length(L), 1:length(L) );
        L(inds) = 0;
        L(inds) = -sum(L);
        
    end
    
elseif strcmpi( LType, 'MVC' )
    
    
    L = mean_value_laplacian(V, F);
    
else
    
    error('Invalid Laplacian Type Supplied!');
    
end

%--------------------------------------------------------------------------
% Duplicate the Laplacian to work on both planar coordinates
%--------------------------------------------------------------------------

RealL = sparse( size(L,1)*2, size(L,2)*2 );
RealL(1:2:end, 1:2:end) = L;
RealL(2:2:end, 2:2:end) = L;

L = RealL;

%==========================================================================
% Solve the Full Constrained System
%==========================================================================

n_vars = size(A,2);
n_eq = size(A,1);

M = [ L A'; A sparse( n_eq, n_eq ) ];
rhs = [ zeros( n_vars, 1 ); b ];

u_lambda = M \ rhs;

e = max( abs( M*u_lambda -rhs ) );
if e > 1e-6
    error('Linear system solution failed!');
end

U = u_lambda(1:n_vars);
U = [ U(1:2:end), U(2:2:end) ];

end

