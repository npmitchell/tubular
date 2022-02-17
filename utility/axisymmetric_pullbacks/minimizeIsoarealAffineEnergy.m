function a = minimizeIsoarealAffineEnergy( face, R3D, R2D )
%MINIMIZEISOAREALAFFINEENERGY Determines the affine transformation of the
%unit square { S = [0,1] X [0,1] }  of the form: f(x,y) = ( a x, y )',
%that most closely matches the (re-scaled) areas of the mesh triangles in
%the plane to their counterparts in 3D
%   Input Parameters:
%       - face               #Fx3 connectivity list
%       - R3D                #Vx3 3D vertex coordinate list
%       - R2D                #Vx2 2D vertex coordinate list.  Assumed to be
%                            a triangulation of the unit square.
%   Output Parameters:
%       - a                 The affine parameter of the best fit map


%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Calculate the area of the 3D mesh ---------------------------------------
e1_3D = R3D( face(:,3), : ) - R3D( face(:,2), : );
e2_3D = R3D( face(:,1), : ) - R3D( face(:,3), : );

A3D = cross( e1_3D, e2_3D, 2 );
A3D = sum( sqrt( sum( A3D.^2, 2 ) ) ./ 2 );

% Calculate the target edge lengths ---------------------------------------
tr3D = triangulation( face, R3D );
eIDx = tr3D.edges;

L3D = R3D( eIDx(:,2), : ) - R3D( eIDx(:,1), : );
L3D = sqrt( sum( L3D.^2, 2 ) );

%--------------------------------------------------------------------------
% MINIMIZATION PROCESSING
%--------------------------------------------------------------------------

fun = @isoarealAffineEnergy;
x0 = 1;
A = [];
b = [];
Aeq = [];
beq = [];
lb = 0 + eps;
ub = Inf;
nonlcon = [];

options = optimoptions( 'fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'none', ...
    'SpecifyObjectiveGradient', true, ...
    'CheckGradients', false );

a = fmincon( fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options );

%--------------------------------------------------------------------------
% ENERGY FUNCTION
%--------------------------------------------------------------------------

    function [ E, EG, EH ] = isoarealAffineEnergy( x )
        
        %------------------------------------------------------------------
        % Calculate the energy
        %------------------------------------------------------------------
        
        % Calculate the new vertex positions
        RR = sqrt(A3D) .* [ sqrt(x) .* R2D(:,1), R2D(:,2) ./ sqrt(x) ];
        
        % Calculate the new edge vectors
        e2D = RR(eIDx(:,2),:) - RR(eIDx(:,1),:);
        
        % Calculate the new edge lengths
        l2D =  sqrt( sum( e2D.^2, 2 ) );
        
        E = sum( ( l2D - L3D ).^2 );
        
        %------------------------------------------------------------------
        % Calculate the energy gradient
        %------------------------------------------------------------------
        if ( nargout > 1 ) % Gradient required
            
            % Calculate the edge unit vectors
            e2DHat = e2D ./ l2D;
            
            % Calculate the derivative of vertex positions
            dRR = sqrt(A3D) .* [ R2D(:,1), -R2D(:,2) ./ x ] ./ ...
                ( 2 .* sqrt(x) );
            
            % Calculate the derivative of edge vectors
            dEIJ = dRR(eIDx(:,2), :) - dRR(eIDx(:,1), :);
            
            % Calculate the derivative of edge lengths
            dl = dot( e2DHat, dEIJ, 2 );
            
            EG = 2 .* sum( ( l2D - L3D ) .* dl, 1 );
            
            %--------------------------------------------------------------
            % Calculate the energy Hessian
            %--------------------------------------------------------------
            if ( nargout > 2 ) % Hessian required
                
                % Calculate the second derivative of vertex positions
                ddRR = sqrt(A3D) .* [-R2D(:,1), 3 .* R2D(:,2) ./ x] ./ ...
                    ( 4 .* x.^(3/2) );
                
                % Calculate the second derivative of edge vectors
                ddEIJ = ddRR(eIDx(:,2), :) - ddRR(eIDx(:,1), :);
                
                % Calculate the second derivative of edge lengths
                ddl = ( dot( e2D, ddEIJ, 2 ) + ...
                    dot( dEIJ, dEIJ, 2 ) - dl.^2 ) ./ l2D;
                
                EH = 2 .* sum( dl.^2 + ( l2D - L3D ) .* ddl, 1 );
                
            end
            
        end
        
    end



end

