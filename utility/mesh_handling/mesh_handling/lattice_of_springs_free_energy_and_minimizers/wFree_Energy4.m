% free energy of a lattice of springs
%
%
function [F] = wFree_Energy4(xyz,TRISm,BLm,KLm)


temp = (xyz(TRISm(:,1),:)-xyz(TRISm(:,2),:));
F = 1/2*sum(KLm.*(sqrt(dot(temp,temp,2))-BLm).^2);
