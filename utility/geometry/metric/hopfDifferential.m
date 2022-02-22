function QQ = hopfDifferential(bb)
% Nematic representation of the traceless symmetric part of the second 
% fundamental form
%
% Parameters
% ----------
% bb : #forms x 4 array or #forms x 1 cell array of 2x2 matrices
%   second fundamental forms 
%
% Returns
% -------
% QQ : #forms x 1 complex array
%
% NPMitchell 2021

if iscell(bb)
    bbv = zeros(length(bb), 4) ;
    for qq = 1:length(bb)
        bbv(qq, :) = reshape(bb{qq}, [1, 4]) ;
    end
    bb = bbv ;
end

LV = bb(:, 1) ;
MV = bb(:, 2) * 0.5 ;
NV = bb(:, 4) ;
QQ = 0.25 * ((LV-NV) - 2*1j*MV ) ;
