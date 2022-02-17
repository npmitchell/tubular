function [BL] = TRI2BL(TRI)
% TRI2BL(TRI) Convert connectivity list to bond list
%
% wtmi boston june 2011
% npmitchell 2019

BL = TRI(:,[1,2]);
BL = cat(1,BL,TRI(:,[1,3]));
BL = cat(1,BL,TRI(:,[2,3]));

BL = sort(BL,2);
BL = unique(BL,'rows');