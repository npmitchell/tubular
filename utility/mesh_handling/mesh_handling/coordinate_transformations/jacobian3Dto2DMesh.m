function J_3D_To_2D = jacobian3Dto2DMesh(v2d, v3d, ff)
%JACOBIAN3DTO2DMESH(v2d, v3d, ff) Construct jacobian for mapping 2d->3d
%
% Parameters
% ----------
% v2d : #vertices x 2 float array
%   The vertex locations in 2d (domain space)
% v3d : #vertices x 3 float array
%   The vertex locations in 3d (image space)
% ff : #faces x 3 int array
%   connectivity list for triangulated meshes
% 
% Returns
% -------
% J_2D_To_3D : length(ff) x 1 cell array
%   A cell array containing the jacobian for each face as each element
% 
% Dillon Cislo, NPMitchell 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Construct Topological Structure Tools
%--------------------------------------------------------------------------
% Construct edge ID list --------------------------------------------------
tri3d = triangulation(ff, v3d) ;
EE = tri3d.edges;

% Construct face-edge correspondence tool ---------------------------------
% Given a list of scalar edge quantities, 'EQ', the output of
% 'EQ(feIDx(f,i))' is that quantity corresponding to the edge opposite the
% ith vertex in face f

% Look for the edge (32, 13, or 21) and be agnostic to the ordering (21=12)
e1IDx = sort( [ ff(:,3), ff(:,2) ], 2 );
e2IDx = sort( [ ff(:,1), ff(:,3) ], 2 );
e3IDx = sort( [ ff(:,2), ff(:,1) ], 2 );

[~, e1IDx] = ismember( e1IDx, EE, 'rows' );
[~, e2IDx] = ismember( e2IDx, EE, 'rows' );
[~, e3IDx] = ismember( e3IDx, EE, 'rows' );

feIDx = [ e1IDx e2IDx e3IDx ];

clear e1IDx e2IDx e3IDx

%% Calculate edge lengths in 3D -------------------------------------------
L_E = v3d(EE(:,2),:) - v3d(EE(:,1),:);
L_E = sqrt(sum(L_E.^2, 2));

L_F = L_E(feIDx);

%% Calculate internal angles in 3D ----------------------------------------

% Some convenience variables to vectorize the cosine law calculation
Gi = L_F; Gj = circshift(L_F, [0 -1]); Gk = circshift(L_F, [0 -2]);

% The cosine of the internal angles
cosAng = ( Gj.^2 + Gk.^2 - Gi.^2 ) ./ ( 2 .* Gj .* Gk );

% The sine of the internal angles
sinAng = sin(acos(cosAng));

%% Calculate the (e1Hat, t1Hat)-basis in 3D ===============================
% Calculate face unit normals  --------------------------------------------
normals = faceNormal( triangulation(ff, v3d) );

% Calculate the (e1Hat, t1Hat)-basis in 3D --------------------------------
e1Hat_3D = (v3d(ff(:,3),:) - v3d(ff(:,2),:)) ./ L_F(:,1);
t1Hat_3D = cross(normals, e1Hat_3D, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Jacobian Operator on Faces ===================================
% Calculate the components of the 2D Jacobian matrix ----------------------
% X-coordinates of triangles in 2D
x1 = v2d(ff(:,1), 1); x2 = v2d(ff(:,2), 1); x3 = v2d(ff(:,3), 1);
x23 = x3-x2; x21 = x1-x2;

% Y-coordinates of triangles in 2D
y1 = v2d(ff(:,1), 2); y2 = v2d(ff(:,2), 2); y3 = v2d(ff(:,3), 2);
y23 = y3-y2; y21 = y1-y2;

J11 = x23 ./ L_F(:,1);
J21 = y23 ./ L_F(:,1);

J12 = x21 ./ L_F(:,3) - x23 .* cosAng(:,2) ./ L_F(:,1);
J12 = J12 ./ sinAng(:,2);

J22 = y21 ./ L_F(:,3) - y23 .* cosAng(:,2) ./ L_F(:,1);
J22 = J22 ./ sinAng(:,2);

% clear x1 x2 x3 x23 x21 y1 y2 y3 y23 y21 cosTheta2 sinTheta2

% Construct full Jacobian from 3D to 2D on each face ----------------------
J_3D_To_2D = cell(size(ff,1),1);

for f = 1:size(ff,1)
    
    % The coordinate change operator into the (e1Hat, t1Hat)-basis
    CC = [ e1Hat_3D(f, :); t1Hat_3D(f, :) ];
    
    % The 2D Jacobian -matrix
    JF = [ J11(f) J12(f); J21(f) J22(f) ];
    
    % Combine to construct complete Jacobian operator
    J_3D_To_2D{f} = JF * CC;
    
end

% clear J11 J21 J12 J22 
% clear fN PP CC JF f