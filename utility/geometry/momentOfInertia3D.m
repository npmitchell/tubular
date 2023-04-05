function ig = momentOfInertia3D(pts) 
% MOMENTOFINERTA3D(pts)
%
% Parameters
% ----------
% pts : N x 3 float or int array
%   point cloud of which to take moment of inertia 
%
% Returns
% -------
% ig : 3 x 3 float array 
%   moment of inertia tensor in 3d
%
%
% NPMitchell 2020
% DJC 2023

validateattributes(pts, {'numeric'}, {'2d', 'ncols', 3, 'finite', 'real'});

% Translate point cloud so that the center of mass lies at the origin
com = mean(pts, 1);
pts = pts - repmat(com, size(pts,1), 1);

X = pts(:,1); Y = pts(:,2); Z = pts(:,3);
X2 = X.^2; Y2 = Y.^2; Z2 = Z.^2;

MOI11 = sum(Y2 + Z2); MOI12 = -sum(X .* Y); MOI13 = -sum(X .* Z);
MOI22 = sum(X2 + Z2); MOI23 = -sum(Y .* Z);
MOI33 = sum(Y2 + X2);

% Assemble moment of inertia
ig = [ MOI11, MOI12, MOI13;
       MOI12, MOI22, MOI23;
       MOI13, MOI23, MOI33 ];

end