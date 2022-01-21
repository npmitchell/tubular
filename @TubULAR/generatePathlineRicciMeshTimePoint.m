function [ricciMesh, ricciMu] = generatePathlineRicciMeshTimePoint(QS, tp, options)
% generatePathlineRicciRicciMeshTimePoint(QS, tp, options)
%
%   Compute Ricci flow pullback for pathline-advected vertices in
%   QS.dir.piv/pathlines/t0_####/quasiconformal/pb2pb/
%   which is the same as
%   QS.dir.pathlines.quasiconformal/pb2pb
% 
% Parameters
% ----------
% QS : quapSlap class instance
% tp : int (default=QS.t0)
%   timestamp of timepoint for which to generate Ricci flow mesh
% options : optional struct with optional fields
%   maxIter : int (default=100)
%   radiusTolerance : float (default=0.01)
%       maximum allowed fractional deviation of Ricci flow solution inner
%       and outer radius from a true circle with fixed radius (variable
%       inner radius, fixed outer radius of 1)
%   save_ims : bool (default = true)
%       save images of the ricci flow solution
%
% Returns
% -------
%
% Saves to disk
% -------------
% Ricci flow solution
%   ricciFn = sprintf(QS.fullFileBase.ricciSolution, maxIter, tp) 
% ricciMesh has fields annulus and rectangle
%   ricciMeshFn = sprintf(QS.fullFileBase.ricciMesh, maxIter, tp) ;
% Beltrami coeff mu for this #iterations
%   mufn = sprintf(QS.fullFileBase.ricciMu, maxIter, tp) ;
%
% Images saved
% ------------
% Plot Ricci result in annulus
%   fullfile(imDir, sprintf('%06d_RicciSolution.png', tp)) ;
% plot beltrami for Ricci flow
%   fullfile(imDir, sprintf('%06d_RicciFlowBeltrami.png', tp)) ;
% Plot just the bare triangulation
%   fullfile(imDir, sprintf('%06d_RicciFlowSolution.png', tp)) ;
% Histogram |mu| for each case
%   fullfile(imDir, sprintf('%06d_BeltramiCoefficients.png', tp)) ;
% Image of corrected vertices on inner and outer annulus
%   fullfile(imDir, sprintf('%06d_ricci_InnerCorrection.png', tp)) ;
%   fullfile(imDir, sprintf('%06d_ricci_OuterCorrection.png', tp)) ;
% Orientation and branch cut
%   fullfile(imDir, sprintf('%06d_DrhoDphi_flipped.png', tp)) ;
%   fullfile(imDir, sprintf('%06d_phiOrderInitial.png', tp)) ;
%   fullfile(imDir, sprintf('%06d_phiOrderFinal.png', tp)) ;
%
%
% NPMitchell 2020

options.pathline_computation = true ;
[ricciMesh, ricciMu] = QS.generateRicciMeshTimePoint(tp, options) ;