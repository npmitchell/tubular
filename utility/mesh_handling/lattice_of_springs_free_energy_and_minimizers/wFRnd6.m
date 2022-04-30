% wFRnd5(tol,xyz,NL,KL,BL,TRISm,KLm,BLm,nd,projectfcn,varargin)
%
% wtmi 2011 - Chicago
% Fletcher-Reeves minimization of lattice energy
% Last changes 15 mar 2011 - Early code was based on Lucks' fr code, rewritten for 
% passing projection functions, protecting bonds, speedup and exit
% commented npm 2014-2015

function [xyz,val,iter,serie] = wFRnd6(tol,xyz,NL,BL,KL,TRISm,BLm,KLm,nd,projectfcn,varargin)

serie = [];

%NL(find(NL==0))=1;

[F,dF] = wFree_Energy4_dF(xyz,NL,BL,KL,TRISm,BLm,KLm);
% H is the negative gradient of the free energy
H=-dF; 
H(isnan(H))=0;
H(nd,:) = 0;

for i=1:3*size(xyz,1)
    FOld = F;
    HOld = H;
    %temp = varargin{:};
    
    %%Debugging:
    %fprintf(['size(xyz)=' ,num2str(size(xyz)),'\n'])
    %fprintf(['size(-dF)=' ,num2str(size(dF )),'\n'])
    
    %optimize amount of displacement in dir along dF to min FreeEngy
    eps = fminbnd(@wFmin4,0,1,[],xyz,H,TRISm,BLm,KLm,projectfcn,varargin{:});
    %eps;
    
    xyz = xyz + eps*H;
    [x,y,z] = projectfcn(xyz(:,1),xyz(:,2),xyz(:,3),varargin{:});
    xyz = [x,y,z];
    
    iter=i;
        
    [F,dF] = wFree_Energy4_dF(xyz,NL,BL,KL,TRISm,BLm,KLm);
    dF(isnan(dF))=0;
    dF(nd,:) = 0;
    H=-dF;
    
    serie = cat(1,serie,F);
    val = F;
    %F-FOld;
 
    if (2*abs(F-FOld) <= tol*(abs(F)+abs(FOld)))
        iter=i;
        val=F;
        return;
    else
        h_new = [H(:,1); H(:,2); H(:,3)];
        h_old = [HOld(:,1); HOld(:,2); HOld(:,3)];
        gamma = dot(h_new,h_new)/dot(h_old,h_old);
        
        if(dot(h_new,h_new) <= 0.)
          %'Zero Gradient Exit'
          iter=i;
          val=F;
          return;
        else
          H = H+gamma*HOld; 
          H(nd,:) = 0;

        end
    end

end
