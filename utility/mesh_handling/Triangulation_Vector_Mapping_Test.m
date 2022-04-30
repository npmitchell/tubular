%% Triangulation Vector Mapping Test ======================================
% This script shows how to construct 'Jacobian' matrix operators that
% calculate the pushforward/pullback of tangent vector fields defined on
% mesh triangulation faces along discrete surface diffeomorphisms
%
% by Dillon Cislo 01/09/2020
%==========================================================================

clear; close all; clc;

%--------------------------------------------------------------------------
% Construct 3D/2D Mesh Triangulations
%--------------------------------------------------------------------------

% Gaussian surface --------------------------------------------------------
diskTri = diskTriangulation(15);

F = diskTri.ConnectivityList;
V2D = diskTri.Points;

AA = 1; % Maximum amplitude of the Gaussian
SS = 0.1; % Controls the width of the Gaussian
RR = 1; % Maximum radial distance in domain of parameterization

gauss3D = @(x,y) AA .* ( exp( -SS .* (x.^2 + y.^2) ) - ...
    exp( -SS .* RR.^2 ) ) ./ ( 1 - exp( -SS .* RR.^2 ) );

V3D = [V2D gauss3D(V2D(:,1), V2D(:,2))];

% Face centroids in 2D
COM2D = cat(3, V2D(F(:,1),:), V2D(F(:,2),:), V2D(F(:,3),:));
COM2D = mean(COM2D, 3);

% Face centroids in 3D
COM3D = cat(3, V3D(F(:,1),:), V3D(F(:,2),:), V3D(F(:,3),:));
COM3D = mean(COM3D, 3);

%--------------------------------------------------------------------------
% Construct Topological Structure Tools
%--------------------------------------------------------------------------

% Construct edge ID list --------------------------------------------------
E = diskTri.edges;

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

clear e1IDx e2IDx e3IDx

%--------------------------------------------------------------------------
% Geometry Processing
%--------------------------------------------------------------------------

% Calculate edge lengths in 2D --------------------------------------------
l_E = V2D(E(:,2),:) - V2D(E(:,1),:);
l_E = sqrt(sum(l_E.^2, 2));

l_F = l_E(feIDx);

clear l_E


% Calculate edge lengths in 3D --------------------------------------------
L_E = V3D(E(:,2),:) - V3D(E(:,1),:);
L_E = sqrt(sum(L_E.^2, 2));

L_F = L_E(feIDx);

error('break')
clear L_E

% Calculate internal angles in 3D -----------------------------------------

% Some convenience variables to vectorize the cosine law calculation
Gi = L_F; Gj = circshift(L_F, [0 -1]); Gk = circshift(L_F, [0 -2]);

% The cosine of the internal angles
cosAng = ( Gj.^2 + Gk.^2 - Gi.^2 ) ./ ( 2 .* Gj .* Gk );

% The sine of the internal angles
sinAng = sin(acos(cosAng));

clear Gi Gj Gk

% Calculate face unit normals  --------------------------------------------
n = faceNormal( triangulation(F, V3D) );

% Calculate the (e1Hat, t1Hat)-basis in 3D --------------------------------
e1Hat_3D = (V3D(F(:,3),:) - V3D(F(:,2),:)) ./ L_F(:,1);
t1Hat_3D = cross(n, e1Hat_3D, 2);

%--------------------------------------------------------------------------
% Construct Tangent Vector Fields in 2D and 3D
%--------------------------------------------------------------------------

% Gaussian surface --------------------------------------------------------

% The derivative of the Z-coordinate with respect to x
gaussDx = @(x,y) 2 .* AA .* SS .* x .* exp(-SS .* (x.^2 + y.^2)) ./ ...
    (exp( -SS .* RR.^2) - 1);

% The derivative of the Z-coordinate with respect to y
gaussDy = @(x,y) 2 .* AA .* SS .* y .* exp(-SS .* (x.^2 + y.^2)) ./ ...
    (exp(-SS .* RR.^2) - 1);

tvCoord = 'y';

switch tvCoord
    
    case 'x'
        
        TV2D = [ ones(size(F,1),1) zeros(size(F,1),1) ];
        TV3D = [ TV2D gaussDx(COM2D(:,1), COM2D(:,2)) ];
        
    case 'y'
        
        TV2D = [ zeros(size(F,1),1) ones(size(F,1),1) ];
        TV3D = [ TV2D gaussDy(COM2D(:,1), COM2D(:,2)) ];
        
    otherwise
        
        error('Invalid coordinate supplied');
        
end

clear tvCoord

%--------------------------------------------------------------------------
% View Results
%--------------------------------------------------------------------------
figure

subplot(1,2,1)
patch('Faces', F, 'Vertices', V3D, ...
    'FaceVertexCData', 0.8 .* ones(size(F)), 'FaceColor', 'flat');
hold on
quiver3( COM3D(:,1), COM3D(:,2), COM3D(:,3), ...
    TV3D(:,1), TV3D(:,2), TV3D(:,3), 1, ...
    'Color', 'r' );
hold off
axis equal
title('3D Surface and Tangent Vector Field');

subplot(1,2,2)
triplot(triangulation(F, V2D));
hold on
quiver(COM2D(:,1), COM2D(:,2), TV2D(:,1), TV2D(:,2), 1, ...
    'Color', 'r' );
hold off
axis equal
title('2D Parameterization and Tangent Vector Field');


%% ************************************************************************
% *************************************************************************
%  PULLBACK 3D VECTOR FIELDS TO 2D VECTOR FIELDS
% *************************************************************************
% *************************************************************************

%% Construct Jacobian Operator on Faces ===================================

% Calculate the components of the 2D Jacobian matrix ----------------------

% X-coordinates of triangles in 2D
x1 = V2D(F(:,1), 1); x2 = V2D(F(:,2), 1); x3 = V2D(F(:,3), 1);
x23 = x3-x2; x21 = x1-x2;

% Y-coordinates of triangles in 2D
y1 = V2D(F(:,1), 2); y2 = V2D(F(:,2), 2); y3 = V2D(F(:,3), 2);
y23 = y3-y2; y21 = y1-y2;

J11 = x23 ./ L_F(:,1);
J21 = y23 ./ L_F(:,1);

J12 = x21 ./ L_F(:,3) - x23 .* cosAng(:,2) ./ L_F(:,1);
J12 = J12 ./ sinAng(:,2);

J22 = y21 ./ L_F(:,3) - y23 .* cosAng(:,2) ./ L_F(:,1);
J22 = J22 ./ sinAng(:,2);

clear x1 x2 x3 x23 x21 y1 y2 y3 y23 y21 cosTheta2 sinTheta2

% Construct full Jacobian from 3D to 2D on each face ----------------------
J_3D_To_2D = cell(size(F,1),1);

for f = 1:size(F,1)
    
    % The coordinate change operator into the (e1Hat, t1Hat)-basis
    CC = [ e1Hat_3D(f, :); t1Hat_3D(f, :) ];
    
    % The 2D Jacobian -matrix
    JF = [ J11(f) J12(f); J21(f) J22(f) ];
    
    % Combine to construct complete Jacobian operator
    J_3D_To_2D{f} = JF * CC;
    
end

clear J11 J21 J12 J22 
clear fN PP CC JF f

%% TEST 1: Pullback Vector Fields to Domain of Parameterization ===========

TV2D_Test = zeros(size(TV2D));
for f = 1:size(F,1)
    TV2D_Test(f,:) = J_3D_To_2D{f} * TV3D(f,:)';
end

%--------------------------------------------------------------------------
% View Results
%--------------------------------------------------------------------------

% Determine the error in the calculation
vecErr = sqrt(sum((TV2D-TV2D_Test).^2,2));
rmsErr = sqrt(mean(vecErr.^2));
fprintf('RMS Error = %d\n', rmsErr);

figure

subplot(1,2,1)
triplot(triangulation(F,V2D), 'k');
hold on
quiver(COM2D(:,1), COM2D(:,2), TV2D(:,1), TV2D(:,2), ...
    1, 'Color', 'r', 'LineWidth', 2);
quiver(COM2D(:,1), COM2D(:,2), TV2D_Test(:,1), TV2D_Test(:,2), ...
    1, 'Color', 'g', 'LineWidth', 2);
hold off
axis equal
title('True Field (Red) and Calculated Field (Green)');

subplot(1,2,2)
patch('Faces', F, 'Vertices', V2D, 'FaceVertexCData', vecErr, ...
    'FaceColor', 'flat');
axis equal
colorbar
title('Calculation (Relative) Error');

clear vecErr rmsErr

%% TEST 2: Calculate 2D Edges From 3D Edges ===============================

vecErr = zeros(size(F,1),1);

for f = 1:size(F,1)
   
    % The edges of the current face in 3D
    e1_3D = V3D(F(f,3), :) - V3D(F(f,2), :);
    e2_3D = V3D(F(f,1), :) - V3D(F(f,3), :);
    e3_3D = V3D(F(f,2), :) - V3D(F(f,1), :);
    
    % The edges of the current face in 2D
    e1_2D = V2D(F(f,3), :) - V2D(F(f,2), :);
    e2_2D = V2D(F(f,1), :) - V2D(F(f,3), :);
    e3_2D = V2D(F(f,2), :) - V2D(F(f,1), :);
    
    % The Jacobian operator for the current face
    JF = J_3D_To_2D{f};
    
    % Pull back the edge vectors to the domain of parameterization
    u1_2D = (JF * e1_3D')';
    u2_2D = (JF * e2_3D')';
    u3_2D = (JF * e3_3D')';
    
    % Raw length of error vector for each edge
    fErr = [ sqrt(sum((e1_2D-u1_2D).^2, 2)), ...
        sqrt(sum((e2_2D-u2_2D).^2, 2)), ...
        sqrt(sum((e3_2D-u3_2D).^2, 2)) ];
    
    % Normalize by true edge length
    fErr = fErr ./ l_F(f,:);
    
    % Update relative error vector
    vecErr(f) = mean(fErr);
    
end

clear fErr JF f
clear e1_3D e2_3D e3_3D
clear e1_2D e2_2D e3_2D
clear u1_2D u2_2D u3_2D

rmsErr = sqrt(mean(vecErr.^2));
fprintf('RMS Error = %d\n', rmsErr);

%--------------------------------------------------------------------------
% View Results
%--------------------------------------------------------------------------
patch('Faces', F, 'Vertices', V2D, 'FaceVertexCData', vecErr, ...
    'FaceColor', 'flat');
axis equal
colorbar
title('Calculation (Relative) Error');

clear vecErr rmsErr


%% ************************************************************************
% *************************************************************************
%  PUSHFORWARD 2D VECTOR FIELDS TO 3D VECTOR FIELDS
% *************************************************************************
% *************************************************************************

%% Construct Jacobian Operator on Faces ===================================

% Some convenience variables
uv12 = V2D(F(:,2), :) - V2D(F(:,1), :);
uv13 = V2D(F(:,3), :) - V2D(F(:,1), :);

xyz12 = V3D(F(:,2), :) - V3D(F(:,1), :);
xyz13 = V3D(F(:,3), :) - V3D(F(:,1), :);

% This just returns the second element of an array
get2 = @(A) A(2);

% I don't really have a physical explanation for this operator, but I'm
% sure there is one. It returns the second element of the wedge product of
% two column vectors. Here {u,v} are each (2x1) column vectors
JOp = @(u,v) get2(u * v' - v * u');

J_2D_To_3D = cell(size(F,1), 1);

for f = 1:size(F,1)
    
    % Calculate the each component of the operator
    JF = zeros([3 2]);
    
    JF(1,1) = -JOp([uv12(f,2); xyz12(f,1)], [uv13(f,2); xyz13(f,1)]);
    JF(2,1) = -JOp([uv12(f,2); xyz12(f,2)], [uv13(f,2); xyz13(f,2)]);
    JF(3,1) = -JOp([uv12(f,2); xyz12(f,3)], [uv13(f,2); xyz13(f,3)]);
    JF(1,2) = JOp([uv12(f,1); xyz12(f,1)], [uv13(f,1); xyz13(f,1)]);
    JF(2,2) = JOp([uv12(f,1); xyz12(f,2)], [uv13(f,1); xyz13(f,2)]);
    JF(3,2) = JOp([uv12(f,1); xyz12(f,3)], [uv13(f,1); xyz13(f,3)]);
    
    JD = JOp(uv12(f,:)', uv13(f,:)');
    
    JF = JF ./ JD;
    
    J_2D_To_3D{f} = JF;
     
end

clear uv12 uv13 xyz12 xyz13 JOp get2

%% TEST 1: Push Vector Fields Forward to 3D ===============================

TV3D_Test = zeros(size(TV3D));
for f = 1:size(F,1)
    TV3D_Test(f,:) = J_2D_To_3D{f} * TV2D(f,:)';
end

%--------------------------------------------------------------------------
% View Results
%--------------------------------------------------------------------------

% Calculate relative error
vecErr = sqrt(sum((TV3D-TV3D_Test).^2,2)) ./ sqrt(sum(TV3D.^2, 2));
rmsErr = sqrt(mean(vecErr.^2));
fprintf('RMS Error = %d\n', rmsErr);

figure

subplot(1,2,1)
patch('Faces', F, 'Vertices', V3D, ...
    'FaceVertexCData', 0.8 .* ones(size(F)), 'FaceColor', 'flat');
hold on
quiver3( COM3D(:,1), COM3D(:,2), COM3D(:,3), ...
    TV3D(:,1), TV3D(:,2), TV3D(:,3), ...
    1, 'Color', 'r' );
quiver3( COM3D(:,1), COM3D(:,2), COM3D(:,3), ...
    TV3D_Test(:,1), TV3D_Test(:,2), TV3D_Test(:,3), ...
    1, 'Color', 'r' );
hold off
axis equal
title('True Field (Red) and Calculated Field (Green)');

subplot(1,2,2)
patch('Faces', F, 'Vertices', V2D, 'FaceVertexCData', vecErr, ...
    'FaceColor', 'flat');
axis equal
colorbar
title('Calculation (Relative) Error');

clear vecErr rmsErr


%% TEST 2: Calculate 3D Edges From 2D Edges ===============================

vecErr = zeros(size(F,1),1);

for f = 1:size(F,1)
   
    % The edges of the current face in 3D
    e1_3D = V3D(F(f,3), :) - V3D(F(f,2), :);
    e2_3D = V3D(F(f,1), :) - V3D(F(f,3), :);
    e3_3D = V3D(F(f,2), :) - V3D(F(f,1), :);
    
    % The edges of the current face in 2D
    e1_2D = V2D(F(f,3), :) - V2D(F(f,2), :);
    e2_2D = V2D(F(f,1), :) - V2D(F(f,3), :);
    e3_2D = V2D(F(f,2), :) - V2D(F(f,1), :);
    
    % The Jacobian operator for the current face
    JF = J_2D_To_3D{f};
    
    % Push fedge vectors forward to 3D
    u1_3D = (JF * e1_2D')';
    u2_3D = (JF * e2_2D')';
    u3_3D = (JF * e3_2D')';
    
    % Raw length of error vector for each edge
    fErr = [ sqrt(sum((e1_3D-u1_3D).^2, 2)), ...
        sqrt(sum((e2_3D-u2_3D).^2, 2)), ...
        sqrt(sum((e3_3D-u3_3D).^2, 2)) ];
    
    % Normalize by true edge length
    fErr = fErr ./ L_F(f,:);
    
    % Update relative error vector
    vecErr(f) = mean(fErr);
    
end

clear fErr JF f
clear e1_3D e2_3D e3_3D
clear e1_2D e2_2D e3_2D
clear u1_3D u2_3D u3_3D

rmsErr = sqrt(mean(vecErr.^2));
fprintf('RMS Error = %d\n', rmsErr);

%--------------------------------------------------------------------------
% View Results
%--------------------------------------------------------------------------
patch('Faces', F, 'Vertices', V2D, 'FaceVertexCData', vecErr, ...
    'FaceColor', 'flat');
axis equal
colorbar
title('Calculation (Relative) Error');

clear vecErr rmsErr





    








