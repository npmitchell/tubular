function [ringpath_ss, hoop_ss ] = ringpathsGridSampling(uspace, vspace, TF, TV2D, TV3Drs)
%RINGPATHSGRIDSAMPLING(uspace, vspace, nU, nV, TF, TV2D, TV3Drs)
% Compute the proper length from the left side of the 2d mesh image TV2D to
% each vertical (u=const) line/ring.
%
% Parameters
% ----------
% uspace : nU x 1 float array
%   sampling grid vector along u (first dim in 'mapping' space)
% vspace : nV x 1 float array
%   sampling grid vector along v (second dimension in 'mapping' space)
% TF : #triangles x 3 int array
%   connectivity list of the mesh
% TV2D : #vertices x 2 float array
%   image of the vertices in 'mapping' space coordinates (2d)
% TV3Drs : #vertices x 3 float array
%   embedding of the vertices in 'realspace' coordinates (3d, true length)
%
% Returns
% -------
% ringpath_ss : nU x 1 float array
%   the average proper pathlength from the left of the rectilinear 2d mesh
%   to the position at uspace(i) is given by ringpath_ss(i)
% hoop_ss : nV x 1 float array (optional)
%   the proper pathlength from the bottom of the rectilinear 2d mesh to the
%   position at vspace(j) and u=uspace(i) is hoop_ss(i, j)
%
% NPMitchell 2020

nU = length(uspace) ;
nV = length(vspace) ;

% NOTE: first dimension indexes u, second indexes v
curves3d = zeros(nU, nV, 3) ;
for kk = 1:nU
    if mod(kk, 20) == 0
        disp(['u = ' num2str(kk / nU)])
    end
    uv = [uspace(kk) * ones(size(vspace)), vspace] ;
    curves3d(kk, :, :) = interpolate2Dpts_3Dmesh(TF, TV2D, TV3Drs, uv) ;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute ds along the surface from each hoop to the next
% The distance from one hoop to another is the
% difference in position from (u_i, v_i) to (u_{i+1}, v_i).
dsuphi = reshape(vecnorm(diff(curves3d), 2, 3), [nU-1, nV]) ;
ringpath_ds = nanmean(dsuphi, 2) ;
ringpath_ss = cumsum([0; ringpath_ds]) ;

if nargout > 1 
    % also return the distances from the 'bottom' of the image (v=0) along
    % v (second dim in mapping space)
    hoop_ss = zeros(nU, nV) ;
    for kk = 1:nU
        error('have not debugged dimension of vecnorm and diff here')
        hoop_ds = reshape(vecnorm(diff(curves3d, 2), 2, 2), [nU-1, nV]) ;
        hoop_ss(kk, :) = cumsum([0; hoop_ds]) ;
    end
end
