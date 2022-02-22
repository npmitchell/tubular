function [BL] = NL2BL(NL)
%NL2BL(NL) Convert neighbor list to bond list BL
%
% Parameters
% ----------
% NL : array of dimension #pts x max(#neighbors)
%   The NL(i,:) contains indices for the neighbors for the ith point, 
%   buffered by zeros if a particle does not have the maximum # nearest 
%   neighbors.
% 
% Returns
% -------
% BL : array of dimension #bonds x 2
%     Each row is a bond and contains indices of connected points.  
%
%
% wtmi 2010, modified:
% - aug 2010 (wtmi)
% - nov 2019 (npmitchell)
%

BL = [];
%x = max(unique(TRIs));
for kk=1:length(NL)
    BL = cat(1,BL,[kk*ones(size(NL,2),1),NL(kk,:)']);
end

tempidx = find(BL(:,1).*BL(:,2));
tempidx = setdiff(1:size(BL,1),tempidx);
BL(tempidx,:) = [];
BL = sort(BL,2);
BL = unique(BL,'rows');
end

