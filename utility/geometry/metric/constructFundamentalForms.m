function [g, b] = constructFundamentalForms(F, V, x)
%Constructs the second fundamental form of a surface
%represented by a mesh triangulation.
%
%   INPUT PARAMETERS
%
%       - F:      #Fx3 face connectivity list
%       - V:      #Vx3 3D vertex coordinate list
%       - x:      #Vx2 2D vertex coordinate list
%
%   OUTPUT PARAMETERS
%
%       - g:    #Fx1 cell array. Entries are a 2x2 matrix representing the
%               first fundamental form on each face
%
%       - b:    #Fx1 cell array. Entries are a 2x2 matrix representing the
%               second fundamental form on each face
%
% See also
% metricSPhiGridMesh.m
%
%   by Dillon Cislo 07/09/2020

%--------------------------------------------------------------------------
% Input Processing
%--------------------------------------------------------------------------

% Validate input triangulation
validateattributes( V, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real'} );
validateattributes( x, {'numeric'}, ...
    {'2d', 'ncols', 2, 'nrows', size(V,1), 'finite', 'real'} );
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'finite', ...
    'real', 'positive', '<=', size(V,1)} );

% A MATLAB-style representation of the input triangulation
TR = triangulation(F,x);

% The edge connectivity list
E = sort(edges(TR), 2);

%--------------------------------------------------------------------------
% Construct Mesh Topological Structure Tools
%--------------------------------------------------------------------------

% Construct face-edge correspondence tool ---------------------------------
efIDx = edgeAttachments(TR, E);

% Construct face-edge correspondence tool ---------------------------------
% Given a list of scalar edge quantities, 'EQ', the output of
% 'EQ(feIDx(f,i))' is that quantity corresponding to the edge opposite the
% ith vertex in face f

e1IDx = sort( [ F(:,3), F(:,2) ], 2 );
e2IDx = sort( [ F(:,1), F(:,3) ], 2 );
e3IDx = sort( [ F(:,2), F(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, E, 'rows' );
[~, e2IDx] = ismember( e2IDx, E, 'rows' );
[~, e3IDx] = ismember( e3IDx, E, 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

% Sparse matrix operator that averages vertex-based quantities onto faces -
V2F = sparse(repmat(1:size(F,1), 3, 1), F.', 1/3, size(F,1), size(V,1));

% Sparse matrix operator that averages face baced quantities onto vertices
% using angle weighting in 3D ---------------------------------------------

% Calculate edge lengths
L3D = V(E(:,2), :) - V(E(:,1), :);
L3D = sqrt(sum(L3D.^2, 2));

% Calculate internal angles
Gi = L3D(feIDx);
Gj = circshift(Gi, [0 -1]);
Gk = circshift(Gi, [0 -2]);
intAngles = ( Gj.^2 + Gk.^2 - Gi.^2 ) ./ ( 2 .* Gj .* Gk );
intAngles = acos(intAngles);
intAngles = intAngles(:);

% The sum of the internal around each vertex
angSumV = full(sparse( F(:), 1, intAngles(:), size(V,1), 1 ));

% The averaging weights
W = intAngles ./ angSumV(F(:));

F2V3D = sparse(F(:), repmat(1:size(F,1),1,3), W, size(V,1), size(F,1));

%--------------------------------------------------------------------------
% Construct Mesh Differential Operators
%--------------------------------------------------------------------------

% Row indices of sparse matrix entries
MI = reshape( repmat(1:size(F,1), 3, 1), [3 * size(F,1), 1] );

% Column indices of sparse matrix entries
MJ = reshape( F.', [3 * size(F,1), 1] );

% Extract edge vectors from faces
ei = x(F(:,3), :) - x(F(:,2), :);
ej = x(F(:,1), :) - x(F(:,3), :);
ek = x(F(:,2), :) - x(F(:,1), :);

% Extract edge vector components
eX = [ ei(:,1), ej(:,1), ek(:,1) ].';
eX = eX(:);

eY = [ ei(:,2), ej(:,2), ek(:,2) ].';
eY = eY(:);

% Extract signed double face areas
sdA = ( ei(:,1) .* ej(:,2) - ej(:,1) .* ei(:,2) );
sdA = repmat( sdA, 1, 3 ).';
sdA = sdA(:);

% Construct sparse operator for Dx
MX = -eY ./ sdA;
Dx = sparse( MI, MJ, MX, size(F,1), size(x,1) );

% Construct sparse operator for Dy
MY = eX ./ sdA;
Dy = sparse( MI, MJ, MY, size(F,1), size(x,1) );

%--------------------------------------------------------------------------
% 3D Geometry Processing
%--------------------------------------------------------------------------

% Calculate face unit normals
nF = faceNormal(triangulation(F,V));

% Calculate vertex unit normals
nV = F2V3D * nF;

%==========================================================================
% CALCULATE THE FUNDAMENTAL FORMS
%==========================================================================

% Calculate surface gradients
DFx = full(Dx * V);
DFy = full(Dy * V);

% Calculate the gradients of the normal vector
DNx = full(Dx * nV);
DNy = full(Dy * nV);

g = cell(size(F,1), 1);
b = cell(size(F,1), 1);

for f = 1:size(F,1)
    
    % The first fundamental form
    g{f} = [ dot( DFx(f,:), DFx(f,:), 2 ), ...
        dot( DFx(f,:), DFy(f,:), 2 ); ...
        dot( DFy(f,:), DFx(f,:), 2 ), ...
        dot( DFy(f,:), DFy(f,:), 2 ) ];
    
    % The second fundamental form
    % b{f} = -[ dot( DFx(f,:), DNx(f,:), 2), ...
    %     dot( DFx(f,:), DNy(f,:), 2); ...
    %     dot( DFy(f,:), DNx(f,:), 2), ...
    %     dot( DFy(f,:), DNy(f,:), 2) ];
    
    % The second fundamental form
    bxy = mean( [ dot( DFx(f,:), DNy(f,:), 2), ...
        dot( DFy(f,:), DNx(f,:), 2) ] );
    b{f} = -[ dot( DFx(f,:), DNx(f,:), 2), bxy; ...
        bxy, dot( DFy(f,:), DNy(f,:), 2) ];
    
end


end

