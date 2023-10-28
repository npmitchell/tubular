function [R, T, K, aa, tt, pp] = rigidAlignment3D(P, Q, varargin)
%RIGIDALIGNMENT3D Finds the rigid 3D transformation (rotation +
%translation) that best aligns two sets of 3D points. It is assumed that
%there is a direct correspondence between the points in each set. 
%
%   INPUT PARAMETERS:
%
%       - P:    #Px3 fixed point 3D coordinate list
%       - Q:    #Px3 moving point 3D coordinate list
%
%   OUTPUT PARAMETERS:
%
%       - R:    3x3 rotation matrix
%       - T:    3x1 translation vector
%       - K:    3x1 rotation unit vector
%       - aa:   Rotation angle
%       - tt:   Polar angle defining the rotation vector
%       - pp:   Azimuthal angle defining the rotation vector
%
%   by Dillon Cislo 02/26/2022

%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

validateattributes(P, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan', 'real'});
validateattributes(Q, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'nonnan', 'real'});

assert(size(P,1) == size(Q,1), ...
    'Point sets must contain the same number of points');

% The number of points in each point cloud
N = size(Q,1);

% Optional Input Processing -----------------------------------------------
% TODO: INITIAL GUESS PROCESSING

dispStyle = 'iter';
maxIter = 1000;
noTrans = false;
noRot = false;

for i = 1:length(varargin)
   if isa(varargin{i}, 'double')
       continue;
   end
   if isa(varargin{i}, 'logical')
       continue;
   end
   if strcmpi(varargin{i}, 'DisplayStyle')
       dispStyle = varargin{i+1};
       validateattributes( dispStyle, {'char'}, ...
           {'vector'} );
   end
   if strcmpi(varargin{i}, 'MaxIterations')
       maxIter = varargin{i+1};
       validateattributes( maxIter, {'numeric'}, ...
           {'scalar', 'positive', 'finite', 'real'} );
   end
   if strcmpi(varargin{i}, 'NoTranslation')
       noTrans = varargin{i+1};
       validateattributes( noTrans, {'logical'}, {'scalar'});
   end
   if strcmpi(varargin{i}, 'NoRotation')
       noRot = varargin{i+1};
       validateattributes( noRot, {'logical'}, {'scalar'});
   end
end

% Return a trivial identity map
if (noTrans && noRot)
    R = eye(3); T = zeros(3,1); K = [1 0 0];
    aa = 0; tt = pi/2; pp = 0;
    return
end

%--------------------------------------------------------------------------
% MINIMIZATION PROCESSING
%--------------------------------------------------------------------------

if noTrans
    
    % Energy function
    fun = @rigidMotionEnergyNoTrans;
    
    % Initial guess (default is identity)
    x0 = [0, pi/2, 0];
    
    % Upper and lower bounds
    lb = [ -pi, 0, -pi ];
    ub = [ pi, pi, pi ];
    
elseif noRot
    
    % Energy function
    fun = @rigidMotionEnergyNoRot;
    
    % Initial guess (default is identity)
    x0 = [0, 0, 0];
    
    % Upper and lower bounds
    lb = [ -Inf, -Inf, -Inf ];
    ub = [ Inf, Inf, Inf];
    
else
    
    % Energy function
    fun = @rigidMotionEnergy;
    
    % Initial guess (default is identity)
    x0 = [0, pi/2, 0, 0, 0, 0];
    
    % Upper and lower bounds
    lb = [ -pi, 0, -pi, -Inf, -Inf, -Inf ];
    ub = [ pi, pi, pi, Inf, Inf, Inf];
    
end

% Inequality constraints
A = [];
b = [];

% Equality constraints
Aeq = [];
beq = [];

% Nonlinear constraints
nonlcon = [];

options = optimoptions( 'fmincon', ...
    'Algorithm', 'interior-point', ...
    'ConstraintTolerance', 1e-6, ...
    'SpecifyObjectiveGradient', true, ...
    'CheckGradients', false, ...
    'MaxIterations', maxIter, ...
    'HessianApproximation', 'bfgs', ...
    'Display', dispStyle );

x = fmincon( fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options );

% Format transformation parameters for output -----------------------------

if noTrans
    
    aa = x(1); % Rotation angle
    tt = x(2); % Polar angle defining rotation vector
    pp = x(3); % Azimuthal angle defining rotation vector
    
    % The translation vector (stored as a column vector)
    T = zeros(3,1);
    
elseif noRot
    
    aa = 0; % Rotation angle
    tt = pi/2; % Polar angle defining rotation vector
    pp = 0; % Azimuthal angle defining rotation vector
    
    % The translation vector (stored as a column vector)
    T = x(1:3);
    if (size(T,2) ~= 1), T = T.'; end
    
else
    
    aa = x(1); % Rotation angle
    tt = x(2); % Polar angle defining rotation vector
    pp = x(3); % Azimuthal angle defining rotation vector
    
    % The translation vector (stored as a column vector)
    T = x(4:6);
    if (size(T,2) ~= 1), T = T.'; end
    
end

% Construct the rotation vector
K = [ sin(tt) * cos(pp); sin(tt) * sin(pp); cos(tt) ];

% The cross product matrix of the rotation vector
CKQ = [ 0, -K(3), K(2); K(3), 0, -K(1); -K(2), K(1), 0 ];

% Construct the rotation matrix from the axis-angle parameters
R = cos(aa) * eye(3) + sin(aa) * CKQ + (1-cos(aa)) * (K * K.');

%--------------------------------------------------------------------------
% ENERGY FUNCTIONS
%--------------------------------------------------------------------------

    function [E, EG] = rigidMotionEnergy( x )
        
        %------------------------------------------------------------------
        % Calculate the energy
        %------------------------------------------------------------------
        
        alpha = x(1); % Rotation angle
        theta = x(2); % Polar angle defining rotation vector
        phi = x(3); % Azimuthal angle defining rotation vector
        
        % The translation vector (stored as a row vector)
        t = x(4:6);
        if (size(t,1) ~= 1), t = t.'; end
        t = repmat(t, N, 1);
        
        % Construct the rotation vector
        k = [ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) ];
        k = repmat(k, N, 1);
        
        % Construct the rotated mpving point positions
        crossKQ = cross( k, Q, 2 );
        dotKQ = repmat(dot( k, Q, 2 ), 1, 3); 
        RQ = Q * cos(alpha) + crossKQ * sin(alpha) + ...
            k .* dotKQ * (1 - cos(alpha));
        
        % Construct aligned separation vectors
        D = P - (RQ + t);
        
        % Calculate energy
        E = sum( dot(D, D, 2) );
        
        % -----------------------------------------------------------------
        % Calculate the energy gradient
        %------------------------------------------------------------------
        
        if ( nargout > 1 )
            
            % Construct rotation vector derivatives
            k_T = [ cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) ];
            k_P = [ -sin(theta) * sin(phi), sin(theta) * cos(phi), 0 ];
            
            k_T = repmat(k_T, N, 1);
            k_P = repmat(k_P, N, 1);
            
            % Construct derivatives of rotated moving point positions
            RQ_A = -(Q - k .* dotKQ) * sin(alpha) + crossKQ * cos(alpha);
            
            crossKTQ = cross(k_T, Q, 2);
            dotKTQ = repmat(dot(k_T, Q, 2), 1, 3);
            RQ_T = crossKTQ * sin(alpha) + ...
                (k_T .* dotKQ + k .* dotKTQ) * (1-cos(alpha));
            
            crossKPQ = cross(k_P, Q, 2);
            dotKPQ = repmat(dot(k_P, Q, 2), 1, 3);
            RQ_P = crossKPQ * sin(alpha) + ...
                (k_P .* dotKQ + k .* dotKPQ) * (1-cos(alpha));
            
            % Calculate translation vector gradient
            EG_Trans = -2 * sum(D, 1);
            
            % Calculate rotation parameter gradients
            EG_Alpha = -2 * sum( dot(D, RQ_A, 2) );
            EG_Theta = -2 * sum( dot(D, RQ_T, 2) );
            EG_Phi = -2 * sum( dot(D, RQ_P, 2) );
            
            EG = [ EG_Alpha, EG_Theta, EG_Phi, EG_Trans ];
            
        end
        
    end

% NO TRANSLATION ----------------------------------------------------------

    function [E, EG] = rigidMotionEnergyNoTrans( x )
        
        %------------------------------------------------------------------
        % Calculate the energy
        %------------------------------------------------------------------
        
        alpha = x(1); % Rotation angle
        theta = x(2); % Polar angle defining rotation vector
        phi = x(3); % Azimuthal angle defining rotation vector
        
        % Construct the rotation vector
        k = [ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) ];
        k = repmat(k, N, 1);
        
        % Construct the rotated mpving point positions
        crossKQ = cross( k, Q, 2 );
        dotKQ = repmat(dot( k, Q, 2 ), 1, 3); 
        RQ = Q * cos(alpha) + crossKQ * sin(alpha) + ...
            k .* dotKQ * (1 - cos(alpha));
        
        % Construct aligned separation vectors
        D = P - RQ;
        
        % Calculate energy
        E = sum( dot(D, D, 2) );
        
        % -----------------------------------------------------------------
        % Calculate the energy gradient
        %------------------------------------------------------------------
        
        if ( nargout > 1 )
            
            % Construct rotation vector derivatives
            k_T = [ cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) ];
            k_P = [ -sin(theta) * sin(phi), sin(theta) * cos(phi), 0 ];
            
            k_T = repmat(k_T, N, 1);
            k_P = repmat(k_P, N, 1);
            
            % Construct derivatives of rotated moving point positions
            RQ_A = -(Q - k .* dotKQ) * sin(alpha) + crossKQ * cos(alpha);
            
            crossKTQ = cross(k_T, Q, 2);
            dotKTQ = repmat(dot(k_T, Q, 2), 1, 3);
            RQ_T = crossKTQ * sin(alpha) + ...
                (k_T .* dotKQ + k .* dotKTQ) * (1-cos(alpha));
            
            crossKPQ = cross(k_P, Q, 2);
            dotKPQ = repmat(dot(k_P, Q, 2), 1, 3);
            RQ_P = crossKPQ * sin(alpha) + ...
                (k_P .* dotKQ + k .* dotKPQ) * (1-cos(alpha));

            % Calculate rotation parameter gradients
            EG_Alpha = -2 * sum( dot(D, RQ_A, 2) );
            EG_Theta = -2 * sum( dot(D, RQ_T, 2) );
            EG_Phi = -2 * sum( dot(D, RQ_P, 2) );
            
            EG = [ EG_Alpha, EG_Theta, EG_Phi ];
            
        end
        
    end

% NO ROTATION -------------------------------------------------------------

function [E, EG] = rigidMotionEnergyNoRot( x )
        
        %------------------------------------------------------------------
        % Calculate the energy
        %------------------------------------------------------------------
        
        % The translation vector (stored as a row vector)
        t = x;
        if (size(t,1) ~= 1), t = t.'; end
        t = repmat(t, N, 1);
        
        % Construct aligned separation vectors
        D = P - (Q + t);
        
        % Calculate energy
        E = sum( dot(D, D, 2) );
        
        % -----------------------------------------------------------------
        % Calculate the energy gradient
        %------------------------------------------------------------------
        
        if ( nargout > 1 )

            % Calculate translation vector gradient
            EG = -2 * sum(D, 1);
            
        end
        
    end

end
