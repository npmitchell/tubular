% wtmi boston une 2011

function [TRI] = TRISm2TRI(TRISm)


nenemat=TRISm2NL(TRISm,9);
TRI = [];

for kk=1:size(nenemat,1)
    
    idxttemp = and(ismember(TRISm(:,1),nenemat(kk,:)),ismember(TRISm(:,2),nenemat(kk,:)));
    TRIStemp = TRISm(idxttemp,:);
    
    
    TRI = cat(1,TRI,cat(2,TRIStemp,kk*ones(size(TRIStemp,1),1)));
    
    
end

TRI = sort(TRI,2);
TRI = unique(TRI,'rows');