function tangent = tangent_from_curve(ss, xp, yp, zp)
%frenetSeretFrame Compute Frenet-Seret frame for a 1d curve in 3d
%   Detailed explanation goes here

dsx = gradient(ss) ;

% First calc rate of change of curve along curve
gradc_raw = [gradient(xp(:)), gradient(yp(:)), gradient(zp(:))] ; 
gradc = bsxfun(@rdivide, gradc_raw, dsx(:)) ;
gradc_ds = vecnorm(gradc, 2, 2) ;
% Compute the tangent to the curve
tangent = bsxfun(@rdivide, gradc, gradc_ds(:)) ;

end
