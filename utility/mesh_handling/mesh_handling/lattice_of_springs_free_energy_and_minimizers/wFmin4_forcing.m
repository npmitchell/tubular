% projectfcn should be a handle to a function that projects on the surface
% format should be projectfcn(x,y,z,varargin{:})

function [F] = wFmin4_forcing(eps,xyz,dF,TRISm,BLm,KLm,force,Fconst,projectfcn,varargin)

% Calculate discrete free energy and 
%
% %Debugging:
% fprintf(['size(xyz)=',num2str(size(xyz)),'\n'])
% fprintf(['size(dF)=' ,num2str(size(dF )),'\n'])

x = xyz(:,1)+eps*dF(:,1);
y = xyz(:,2)+eps*dF(:,2);
z = xyz(:,3)+eps*dF(:,3);

[x,y,z] = projectfcn(x,y,z,varargin{:});

% disp('wFmin4_forcing --> x(1:3) = ')
% z(1:3)
F = wFree_Energy4_forcing([x,y,z],TRISm,BLm,KLm,force,Fconst);

return;