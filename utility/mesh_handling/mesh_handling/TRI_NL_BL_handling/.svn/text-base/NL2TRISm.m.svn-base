% wtmi 2010 modified:
% - aug 2010 (wtmi)
% 
%


function [TRISm] = NL2TRISm(NL)
TRISm = [];
%x = max(unique(TRIs));
for kk=1:length(NL)
    TRISm = cat(1,TRISm,[kk*ones(size(NL,2),1),NL(kk,:)']);
end

tempidx = find(TRISm(:,1).*TRISm(:,2));
tempidx = setdiff(1:size(TRISm,1),tempidx);
TRISm(tempidx,:) = [];
TRISm = sort(TRISm,2);
TRISm = unique(TRISm,'rows');


