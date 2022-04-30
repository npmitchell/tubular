% projectfcn should be a handle to a function that projects on the surface
% format should be projectfcn(x,y,z,varargin{:})

function [F] = wFmin4(eps,xyz,dF,TRISm,BLm,KLm,projectfcn,varargin)

% Calculate discrete free energy and 
% force on each particles
%
x=xyz(:,1)+eps*dF(:,1);
y=xyz(:,2)+eps*dF(:,2);
z=xyz(:,3)+eps*dF(:,3);

[x,y,z] = projectfcn(x,y,z,varargin{:});

F = wFree_Energy4([x,y,z],TRISm,BLm,KLm);

return;