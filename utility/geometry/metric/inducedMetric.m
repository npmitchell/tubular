function g = inducedMetric(F, V3D, V2D)
%INDUCEDMETRIC 
%   Calculates the discrete induced metric of a surface
%   represented by a mesh triangulation. 
%   Unlike constructFundamentalForms.m, this one is not gradient-based, so 
%   it is cleaner. However, constructFundamentalForms.m also computes the 
%   second fundamental form.
%
%   INPUT PARAMETERS
%
%       - F:        #Fx3 face connectivity list
%       - V3D:      #Vx3 3D vertex coordinate list
%       - V3D:      #Vx2 2D vertex coordinate list
%
%   OUTPUT PARAMETERS
%
%       - g:    #Fx1 cell array. Each entry is a 2x2 matrix representing
%               the induced metric on a face
%
%   SEE ALSO
%   constructFundamentalForms.m : gives first and second fundamental form 
%
%   by Dillon Cislo 05/12/2020

% Validate input triangulation
validateattributes( V3D, {'numeric'}, ...
    {'2d', 'ncols', 3, 'finite', 'real'} );
validateattributes( V2D, {'numeric'}, ...
    {'2d', 'ncols', 2, 'nrows', size(V3D,1), 'finite', 'real'} );
validateattributes(F, {'numeric'}, ...
    {'2d', 'ncols', 3, 'integer', 'finite', ...
    'real', 'positive', '<=', size(V3D,1)} );

% The 3D directed edge vectors
e13D = V3D(F(:,3), :) - V3D(F(:,2), :);
e23D = V3D(F(:,1), :) - V3D(F(:,3), :);
e33D = V3D(F(:,2), :) - V3D(F(:,1), :);

% The 3D edge lengths
L1 = sqrt( sum( e13D.^2, 2 ) );
L2 = sqrt( sum( e23D.^2, 2 ) );
L3 = sqrt( sum( e33D.^2, 2 ) );

% The 2D directed edge vectors
e12D = V2D(F(:,3), :) - V2D(F(:,2), :);
e22D = V2D(F(:,1), :) - V2D(F(:,3), :);
e32D = V2D(F(:,2), :) - V2D(F(:,1), :);

% Augmented edge unit vectors (used for cross products)
e12D = [ e12D, zeros(size(e12D,1), 1) ];
e22D = [ e22D, zeros(size(e22D,1), 1) ];
e32D = [ e32D, zeros(size(e32D,1), 1) ];

% The face areas/face unit normals of the pullback mesh
% ( +Z if the face is CCW oriented, -Z if the face is CW oriented )
fN = cross( e12D, e22D, 2 );
fA = sqrt( sum( fN.^2, 2 ) ) ./ 2;

fN = fN ./ (2 .* fA);

% The outward pointing in-plane edge normals
t1 = cross(e12D, fN, 2 ); t1 = t1(:, [1 2]);
t2 = cross(e22D, fN, 2 ); t2 = t2(:, [1 2]);
t3 = cross(e32D, fN, 2 ); t3 = t3(:, [1 2]);

% Construct the induced metric
g = cell( size(F,1), 1 );
for f = 1:size(F,1)
    
    g{f} = -( ...
        ( L1(f).^2-L2(f).^2-L3(f).^2 ) .* kron( t1(f,:), t1(f,:)' ) + ...
        ( L2(f).^2-L3(f).^2-L1(f).^2 ) .* kron( t2(f,:), t2(f,:)' ) + ...
        ( L3(f).^2-L1(f).^2-L2(f).^2 ) .* kron( t3(f,:), t3(f,:)' ) ...
        ) ./ ( 8 .* fA(f).^2 );
    
end

end

