function [pNL2] = bondLongerthan(xyz,pNL,pSL,pAL,len)
% WTMI & NPMitchell

pdr = double(0.0*pNL);

for ii=1:size(pNL,1)
    for jj=1:size(pNL,2)
        if (pNL(ii,jj) ~= 0)
            vec=xyz(ii,:)-xyz(pNL(ii,jj),:);
            dis=sqrt(vec*vec');
            pdr(ii,jj) = dis; %-pAL(ii,jj);
            %U(ii) = U(ii)+pSL(ii,jj)*(dis-pAL(ii,jj))^2;
            %%dF(i,:) = dF(i,:)+pSL(i,j)*vec*(dis-pAL(i,j))/dis;
            %zFac = 0;
            %dFac = pSL(ii,jj)*(dis-pAL(ii,jj))/dis;
            %dF(ii,1) = dF(ii,1) + dFac*(vec(1) - 0*vec(3)*xyz(ii,1)*zFac);
            %dF(ii,2) = dF(ii,2) + dFac*(vec(2) - 0*vec(3)*xyz(ii,2)*zFac);
            %dF(ii,3) = dF(ii,3) + dFac*(vec(3) - 0*vec(3)*xyz(ii,2)*zFac);
        end
    end
    %U(ii) = 0.5*U(ii);
end

%[tidx1,tidx2]
tidx = find((pdr<len));
pNL2 = pNL;
pNL2(tidx) = 0;

%pNL2 = pdr./pAL;

%F = 0.5*sum(U);
return;
    

