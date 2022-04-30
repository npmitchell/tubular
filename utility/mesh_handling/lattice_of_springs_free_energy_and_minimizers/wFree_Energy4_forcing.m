% free energy of a lattice of springs under forcing boundary conditions
% (von Neumann bcs --> not simply constant displacement)
%
%
function [F] = wFree_Energy4_forcing(xyz,TRISm,BLm,KLm,force,Fconst)


temp = (xyz(TRISm(:,1),:)-xyz(TRISm(:,2),:));
F = 1/2*sum(KLm.*(sqrt(dot(temp,temp,2))-BLm).^2);

%For each node where there is a boundary force, subtracting the work done
%by the boundary forces from the energy.
matrix = force .* xyz ;
F = F - sum(matrix(:)) + Fconst; 

return
