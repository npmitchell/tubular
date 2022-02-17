% wtmi boston une 2011

function [TRISm] = TRI2TRISm(TRI)

TRISm = TRI(:,[1,2]);
TRISm = cat(1,TRISm,TRI(:,[1,3]));
TRISm = cat(1,TRISm,TRI(:,[2,3]));


TRISm = sort(TRISm,2);
TRISm = unique(TRISm,'rows');