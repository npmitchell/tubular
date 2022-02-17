function [tr, dev, theta] = trace_deviator(eq, gq)
%[tr, dev, theta] = trace_deviator(eq, gq)
% Compute traceful component (dilation), deviatoric magnitude and angle
% (deviator) from 2x2 tensor eq on metric gq
% 
% Parameters
% ----------
% eq : 2x2 numeric array
%   tensor to decompose
% gq : 2x2 numeric array
%   metric tensor, first fundamental form of surface
% 
% NPMitchell 2020

%% traceful component -- 1/2 Tr[g^{-1} gdot] = Tr[g^{-1} eps] 
tr = trace(inv(gq) * (eq)) ;
%% deviatoric component -- 
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
theta = atan2(pevec(2), pevec(1)) ;

%% Test to debug
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


