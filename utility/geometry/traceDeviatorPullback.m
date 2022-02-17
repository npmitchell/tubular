function [tr, dev, theta, theta_pb] = traceDeviatorPullback(eq, gq, dx, dy)
%[tr, dev, theta] = trace_deviator(eq, gq)
% Compute traceful component (dilation), deviatoric magnitude and angle
% (deviator) from 2x2 tensor eq on metric gq.
% Eigenvector is first in units of (a,b) = a vec_zeta + b vec_phi, 
% where vec_zeta/phi are vectors of unit length in the pullback space, 
% but we want to have (a,b) = a hat_zeta + b hat_phi, 
% where hat_zeta/phi are vectors of unit length in the embedding space 
% along the directions of the transformed zetahat and phihat from the 
% pullback space. 
% To rescale here, scale a and b by the length in embedding space per unit
% length in the pullback space. a --> a * dl_embedding / dl_pullback for
% this face under consideration.
% 
% Parameters
% ----------
% eq : 2x2 numeric array
%   tensor to decompose
% gq : 2x2 numeric array
%   metric tensor, first fundamental form of surface
% dx : float
%   length in embedding space per unit length in the pullback space of zeta
%   This is used to obtain the proper deviatoric angle.
%   a --> a * dl_embedding / dl_pullback for this face under consideration.
% dy : float
%   length in embedding space per unit length in the pullback space of phi
%   This is used to obtain the proper deviatoric angle.
%   b --> b * dl_embedding / dl_pullback for this face under consideration.
%
%
% Returns
% -------
% tr : float
%   magnitude of strain dilation (traceful)
% dev: float
%   magnitude of strain deviator 
% theta: float 
%   angle of elongational component of strain deviator in embedding space 
%   relative to zetahat projected into embedding space. 
% theta_pb : float
%   angle of elongational component of strain deviator in pullback space 
%
%
% See also
% --------
% trace_deviator : does not scale deviator elongation axis by embedding
% 
% NPMitchell 2020

% traceful component -- 1/2 Tr[g^{-1} gdot] = Tr[g^{-1} eps] 
tr = trace(inv(gq) * (eq)) ;
% deviatoric component -- 
% || epsilon - 1/2 Tr[g^{-1} epsilon] g|| = sqrt(Tr[A A^T]),
% where A = epsilon - 1/2 Tr[g^{-1} epsilon] g.
AA = eq - 0.5 * tr * gq ;
dev = sqrt(trace(inv(gq) * (AA * (inv(gq) * AA)))) ;

%% angle of elongation -- first take eigvectors
[evec_dev, evals_dev] = eig( inv(gq) * AA ) ;
% [evec_dev, evals_dev] = eig( AA * inv(gq) ) ;
[~, idx] = sort(diag(evals_dev)) ;
evec_dev = evec_dev(:, idx) ;
pevec = evec_dev(:, end) ;
% pevec is in units of (a,b) = a vec_zeta + b vec_phi, 
% where vec_zeta/phi are vectors of unit length in the pullback space, 
% but we want to have (a,b) = a hat_zeta + b hat_phi, 
% where hat_zeta/phi are vectors of unit length in the embedding space 
% along the directions of the transformed zetahat and phihat from the 
% pullback space. 
% To rescale here, scale a and b by the length in embedding space per unit
% length in the pullback space. a --> a * dl_embedding / dl_pullback for
% this face under consideration.
try
    theta = mod(atan2(pevec(2) * dy, pevec(1) * dx), pi) ;
catch
    disp('WARNING: imaginary result for theta. Taking real component')
    theta = mod(atan2(real(pevec(2)) * dy, real(pevec(1)) * dx), pi) ;
end

if nargout > 3
    theta_pb = mod(atan2(pevec(2), pevec(1)), pi) ;
end

%% Test
% eq = [0., 0.1; 0.1, 0.] ;
% gq = [1, 0; 0, 2] ;
% tr = trace(inv(gq) * (eq)) ;
% % deviatoric component -- 
% % || epsilon - 1/2 Tr[g^{-1} epsilon] g|| = sqrt(Tr[A A^T]),
% % where A = epsilon - 1/2 Tr[g^{-1} epsilon] g.
% AA = eq - 0.5 * tr * gq ;
% dev = sqrt(trace(inv(gq) * (AA * (inv(gq) * AA)))) ;
% % angle of elongation -- first take eigvectors
% [evec_dev, evals_dev] = eig( AA * inv(gq) ) ;
% [~, idx] = sort(diag(evals_dev)) ;
% evec_dev = evec_dev(:, idx) ;
% pevec = evec_dev(:, end) ;
% theta = atan2(pevec(2), pevec(1)) ;


