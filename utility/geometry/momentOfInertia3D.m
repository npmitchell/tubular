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

com = mean(pts) ;
% Compute moment of inertia
ig = [0.0,0.0,0.0;
      0.0,0.0,0.0;
      0.0,0.0,0.0];
for ii = 1:size(pts, 1)
    x = pts(ii, 1) - com(1) ;
    y = pts(ii, 2) - com(2) ;
    z = pts(ii, 3) - com(3) ;
    ig = ig + [ (y^2 + z^2),    -x*y,         -x*z; ...
             -y*x,             (x^2 + z^2), -y*z; ...
             -x*z,            -y*z,          (y^2 + x^2)];
end