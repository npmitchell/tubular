function [phi, u, v] = rigidParameterizationAlignment( ...
    tarF, tarV3D, tarV2D, movF, movV3D, movV2D )
%RIGIDPARAMETERIZATIONALIGNMENT Finds the 2D rigid motion (rotation +
%translation) that most closely aligns the 2D parameterizations of a pair
%of 3D mesh triangulations
%
% Outputs are interpreted as follows:
% rotm = [cos(phi), sin(phi); -sin(phi), cos(phi)] ;
% newV2D = (movV2D + [u, v]) * rotm  ;
%
% Or equivalently up to machine precision:
% transV2D = [ movV2D(:,1) + u, movV2D(:,2) + v ];
% newV2D = [ cos(phi).*transV2D(:,1)-sin(phi).*transV2D(:,2), ...
%            sin(phi).*transV2D(:,1)+cos(phi).*transV2D(:,2) ];
% 
%
%   INPUT PARAMETERS:
%
%       - tarF:     #F1x3 connectivity list of the target triangulation
%       - tarV3D:   #V1x3 3D target vertex coordinate list
%       - tarV2D:   #V1x2 2D target vertex coordinate list
%       - movF:     #F2x3 connectivity list of the moving triangulation
%       - movV3D:   #V2x3 3D moving vertex coordinate list
%       - movV2D:   #V2x2 2D moving vertex coordinate list
%
%   OUTPUT PARAMETERS:
%
%       - phi:      The optimal rotation angle
%       - u:        The x-component of the optimal translation
%       - v:        The y-component of the optimal translation
%
%   by Dillon Cislo 12/27/2020

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Validate the target triangulation
validateattributes( tarF, {'numeric'}, ...
    {'2d', 'ncols', 3, 'positive', 'integer', 'finite', 'real'} );
validateattributes( tarV3D, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real'} );
validateattributes( tarV2D, {'numeric'}, ...
    {'2d', 'ncols', 2, 'finite', 'real'} );

assert( size(tarV3D,1) == size(tarV2D,1), ...
    'Target vertex lists are improperly sized');

% Validate the moving triangulation
validateattributes( movF, {'numeric'}, ...
    {'2d', 'ncols', 3, 'positive', 'integer', 'finite', 'real'} );
validateattributes( movV3D, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real'} );
validateattributes( movV2D, {'numeric'}, ...
    {'2d', 'ncols', 2, 'finite', 'real'} );

assert( size(tarV3D,1) == size(tarV2D,1), ...
    'Target vertex lists are improperly sized');

% Point match the 3D triangulations to establish a vertex correspondence
matchIDx = zeros(size(movV2D,1), 1);

for i = 1:size(movV2D,1)
    
    curID = tarV3D - repmat(movV3D(i,:), size(tarV3D,1), 1);
    curID = sqrt(sum(curID.^2, 2));
    
    [~, curID] = min(curID);
    
    matchIDx(i) = curID;
    
end

matchV2D = tarV2D(matchIDx, :);

%--------------------------------------------------------------------------
% MINIMIZATION PROCESSING
%-------------------------------------------------------------------------

fun = @rigidMotionEnergy;
x0 = [0 0 0];

options = optimoptions( 'fminunc', ...
    'Algorithm', 'trust-region', ...
    'SpecifyObjectiveGradient', true, ...
    'CheckGradients', false, ...
    'Display', 'none');

optimX = fminunc( fun, x0, options );
phi = optimX(1);
u = optimX(2);
v = optimX(3);

%--------------------------------------------------------------------------
% ENERGY FUNCTION
%--------------------------------------------------------------------------
    function [E, EG] = rigidMotionEnergy( x )
        % x = [phi u v]
        
        %------------------------------------------------------------------
        % Calculate the energy
        %------------------------------------------------------------------
        
        % Calculate the new vertex positions
        transV2D = [ movV2D(:,1) + x(2), movV2D(:,2) + x(3) ];
        
        newV2D = [ cos(x(1)).*transV2D(:,1)-sin(x(1)).*transV2D(:,2), ...
            sin(x(1)).*transV2D(:,1)+cos(x(1)).*transV2D(:,2) ];

        % The separation vector between the new vertex positions and their
        % targets
        diff2D = newV2D - matchV2D;
        
        E = sum( dot(diff2D, diff2D, 2) );
        
        %------------------------------------------------------------------
        % Calculate the energy gradient
        %------------------------------------------------------------------
        if ( nargout > 1 )
            
            % Calculate the gradient wrt phi
            gradPhiVec = ...
                [ -sin(x(1)).*transV2D(:,1)-cos(x(1)).*transV2D(:,2), ...
                cos(x(1)).*transV2D(:,1)-sin(x(1)).*transV2D(:,2) ];
            
            EG(1) = 2 .* sum( dot(diff2D, gradPhiVec, 2) );
            
            % Calculate the gradient wrt to u
            gradUVec = ...
                repmat([cos(x(1)), sin(x(1))], size(diff2D,1), 1);
            EG(2) = 2 .* sum( dot(diff2D, gradUVec, 2) );
            
            % Calculate the gradient wrt to v
            gradUVec = ...
                repmat([-sin(x(1)), cos(x(1))], size(diff2D,1), 1);
            EG(3) = 2 .* sum( dot(diff2D, gradUVec, 2) );
            
        end
        
    end

end

