function [NL] = BL2NL(BL,varargin)
% wtmi 2010 npmitchell 2013, 2019:
% - aug 2010 (wtmi)
% - june 2011 (wtmi) - varargin is size of NL
% - 2013, 2019 (npmitchell) automatically size NL based on neede rows
% - 2019 (npmitchell) renamed to BL2NL from TRISm2NL

if size(varargin,2)==1
    n=varargin{1};
else
    n=15 ;
end

NL = [];
x = max(unique(BL));
for kk=1:x
    temp = BL(or(BL(:,1)==kk,BL(:,2)==kk),:);
    temp = unique(temp(:));
    temp = setdiff(temp,kk)';
    if size(temp,1)>size(temp,2)
        temp = temp';
    end
    if ~isempty(temp)    
        temp =  cat(2,temp,zeros(1,(n-size(temp,2))));
        NL = cat(1,NL,temp);
    else
        temp = zeros(1,n);
        NL = cat(1,NL,temp);
    end
    %temp = temp(1:n);    
end

% Now truncate the unused columns is varargin is not supplied
if size(varargin,2)~=1
    NL = NL(:, 1:find(nonzeros(sum(NL, 1)), 1, 'last' )) ;
end

