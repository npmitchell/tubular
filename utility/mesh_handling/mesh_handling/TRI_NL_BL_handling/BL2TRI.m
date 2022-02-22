function [TRI] = BL2TRI(BL, varargin)
% BL2TRI(BL) Convert bond list to connectivity list TRI
%
% Parameters
% ----------
% BL : array of dimension #bonds x 2
%     Each row is a bond and contains indices of connected points.
%
% Returns
% ------- 
% TRI : #faces x 3 int array
%   triangle array, with TRI(i, :) being indices into pts of ith face
%
%
% wtmi boston une 2011
% NPMitchell 2019 

% Decide how many columns to have in NL
if size(varargin,2)==1
    n=varargin{1};
else
    n=15 ;
end

nenemat=BL2NL(BL,n);
TRI = [];
for kk=1:size(nenemat,1)
    idxttemp = and(ismember(BL(:,1),nenemat(kk,:)),ismember(BL(:,2),nenemat(kk,:)));
    TRIStemp = BL(idxttemp,:);
    TRI = cat(1,TRI,cat(2,TRIStemp,kk*ones(size(TRIStemp,1),1)));
end

TRI = sort(TRI,2);
TRI = unique(TRI,'rows');