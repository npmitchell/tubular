% wtmi boston june 2011
% fixed tricomplete bug Oct 2011 chicago

function [TRI,TRISm,NL] = tricomplete(varargin)
    
    % go thro varargin and depending on number of columns call it 
    % TRI,TRISm or NL
    for kk=1:size(varargin,2)
        argin = varargin{kk};
        switch size(argin,2)
            case 2
                TRISm = argin;
            case 3
                TRI = argin;
            otherwise
                NL = argin;
        end
    end    
   
    fatti = [exist('TRI'),exist('TRISm'),exist('NL')];
    
    % make TRI
    if and(fatti(1)==0,fatti(2)==1)
        TRI = TRISm2TRI(TRISm);
    elseif and(fatti(1)==0,fatti(2)==0)
        TRISm = NL2TRISm(NL);
        TRI = TRISm2TRI(TRISm);
    end
    
    % make TRISm
    %TRISm = TRI2TRISm(TRI);
    TRISm = [TRI(:,[1,2]);TRI(:,[1,3]);TRI(:,[2,3])];
    TRISm = sort(TRISm,2);
    TRISm = unique(TRISm,'rows');
    
    
    fatti = [exist('TRI'),exist('TRISm'),exist('NL')];

    % make NL
    if ~exist('NL')
        NL = TRISm2NL(TRISm,60);
    end
        
    % Check for consistency ?
    
    
        
