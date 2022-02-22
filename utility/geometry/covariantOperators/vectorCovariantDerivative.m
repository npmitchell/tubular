function [DU3D, DU2D] = vectorCovariantDerivative(U, F, V3D, V2D)
%VECTORCOVARIANTDERIVATIVE Calculate the covariant derivative of a vector
%field on a mesh triangulation.
%
%   INPUT PARAMETERS:
%
%       - U:        A #Vx3 or #Fx3 vector field defined on the input mesh
%       - F:        #Fx3 face connectivity list
%       - V3D:      #Vx3 3D vertex coordinate list
%       - V2D:      #Vx2 2D vertex coordinate list
%
%   OUTPUT PARAMETERS:
%
%       - DU3D:     An #Fx1 cell array. Each entry is a 3x3 matrix. This
%                   matrix is a representation of the covariant derivative
%                   of U on each face.  For a (3x1) tangent vector W
%                   defined on a face f the action of DU3D{f} * W is to
%                   produce another tangent vector that is the covariant
%                   derivative of U in the direction of W
%
%       - DU2D:     An #Fx1 cell array. Each entry is a 2x2 matrix.  This
%                   this matrix is a representation of the covariant
%                   derivative in the (ex,ey)-basis induced by the
%                   Cartesian coordinates of the 2D parameterization of the
%                   3D surface mesh.
%
%   by Dillon Cislo 03/13/2020

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Check for proper number of inputs
if (nargin < 1), error('Please supply 3D vector field'); end
if (nargin < 2), error('Please supply face connectivity list'); end
if (nargin < 3), error('Please supply vertex coordinate list'); end
if (nargin < 4), V2D = []; end


% Validate input attributes
validateattributes( U, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan', 'real'} );
validateattributes( V3D, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan', 'real'} );
validateattributes( F, {'numeric'}, ...
    {'2d', 'ncols', 3', 'positive', ...
    'integer', 'real', '<=', size(V3D,1)} );

if ~isempty(V2D)
    validateattributes( V2D, {'numeric'}, ...
        {'2d', 'ncols', 2, 'nrows', size(V3D,1), ...
        'finite', 'nonnan', 'real'} );
end

% Re-define vector field on vertices if necessary
if ~isequal(size(U), size(V3D))
    
    % This is a failed attempt at creating a meshAveragingOperator
    % F2V = internalangles(V3D, F); % The internal angles of the mesh
    % F2V = sparse( F(:), repmat(1:size(F,1),1,3), F2V, ...
    %     size(V3D,1), size(F,1) );
    
    [~, F2V] = meshAveragingOperators(F, V3D) ;
    
    U = F2V * U;
    
end

%--------------------------------------------------------------------------
% CALCULATE 3D COVARIANT DERIVATIVES
%--------------------------------------------------------------------------

% Calculate face unit normal vectors
FN = faceNormal(triangulation(F,V3D));

% FEM gradient operator
G = grad(V3D, F);

% The directional derivatives
DUx = reshape( G * U(:,1), size(F) );
DUy = reshape( G * U(:,2), size(F) );
DUz = reshape( G * U(:,3), size(F) );

DU3D = cell(size(F,1), 1);
for i = 1:size(F,1)
    
    % Construct face projection operator
    P = eye(3) - FN(i,:).' * FN(i,:);
    
    % Calculate covariant derivative
    % extract the part tangent to the face, then ensure that result is
    % tangent to the face
    DU3D{i} = P * [ DUx(i,:); DUy(i,:); DUz(i,:) ] * P;
    
end

%--------------------------------------------------------------------------
% CALCULATE 2D COVARIANT DERIVATIVES
%--------------------------------------------------------------------------
% 2D results are returned as a type-(0,2) tensor

DU2D = [];
if ~isempty(V2D)
    
    % The Jacobian matrix dX/du
    J23 = jacobian2Dto3DMesh( V2D, V3D, F );
    
    % The Jacobian matrix du/dX
    J32 = jacobian3Dto2DMesh( V2D, V3D, F );
    
    % Transform the covariant derivative
    DU2D = cell(size(F,1), 1);
    for i = 1:size(F,1)
        
        % These actually produce the same results
        % The second line includes an explicit index lowering using the
        % discrete metric
        % DU2D{i} = J23{i}.' * DU3D{i}.' * J23{i};
        DU2D{i} = J23{i}.' * DU3D{i}.' * J32{i}.' * J23{i}.' * J23{i};
        
    end
    
end


end

