function [TRI] = NL2TRI(NL, BL)
%NL2TRI(NL, BL) Convert neighbor list to connectivity list
%   
% Parameters
% ----------
% NL : array of dimension #pts x max(#neighbors)
%   The NL(i,:) contains indices for the neighbors for the ith point, 
%   buffered by zeros if a particle does not have the maximum # nearest 
%   neighbors.
% BL : array of dimension #bonds x 2
%     Each row is a bond and contains indices of connected points.
%
% Returns
% -------
% TRI : #faces x 3 int array
%   triangle array, with TRI(i, :) being indices into pts of ith face
%
% NPMitchell 2019

TRI = [];
for kk=1:size(NL,1)
    idxttemp = and(ismember(BL(:,1),NL(kk,:)),ismember(BL(:,2),NL(kk,:)));
    TRIStemp = BL(idxttemp,:);
    TRI = cat(1,TRI,cat(2,TRIStemp,kk*ones(size(TRIStemp,1),1)));
end

TRI = sort(TRI,2);
TRI = unique(TRI,'rows');
end

