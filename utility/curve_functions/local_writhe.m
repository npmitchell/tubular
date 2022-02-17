function wr = local_writhe(ss, xyz)
%LOCAL_WRITHE Compute Wr contribution of each point in a curve 
%   This quantity (wr) is the amount of writhe stored in the segment 
%   (r, r + dr). Method based on Berger & Prior, The writhe of open and 
%   closed curves, 2006.
% 
% Note that in the paper, prime denotes d/dz.
dz = gradient(xyz(:, 3));
ds = gradient(ss) ;
lamb = dz ./ ds ;

% [tangent, ~, ~] = frenetSeretFrame(ss, xyz(:, 1), xyz(:, 2), xyz(:, 3)) ;

% Just compute the tangent
gradc_raw = [gradient(xyz(:, 1)), gradient(xyz(:, 2)), gradient(xyz(:, 3))] ; 
gradc = bsxfun(@rdivide, gradc_raw, ds(:)) ;
gradc_ds = vecnorm(gradc, 2, 2) ;
% Compute the tangent to the curve
tangent = bsxfun(@rdivide, gradc, gradc_ds(:)) ;

tangent_gradient = [gradient(tangent(:, 1)), ...
    gradient(tangent(:, 2)), ...
    gradient(tangent(:, 3))] ; 

tangentprime = tangent_gradient ./ dz ;

% take cross product of tangent with change of tangent
txtprime = cross(tangent, tangentprime) ;
wr = (1 / (2 * pi)) * (ones(length(lamb), 1) ./ (1 + abs(lamb))) .* txtprime(:, 3) ;
wr = wr .* abs(dz) ;

end

