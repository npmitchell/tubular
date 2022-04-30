function [pdr] = bondLengths(xyz,NL)
% WTMI & NPMitchell

pdr = double(0.0*NL);

for ii=1:size(NL,1)
    for jj=1:size(NL,2)
        if (NL(ii,jj) ~= 0)
            vec=xyz(ii,:)-xyz(NL(ii,jj),:);
            dis=sqrt(vec*vec');
            pdr(ii,jj) = dis;
        end
    end
end

return;
    

