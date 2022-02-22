function BL = Vdat2BL(Vdat, options)
%BL = Vdat2BL(Vdat, options)
% 
% Note: no functionality for periodic bonds
% 
% Parameters
% ----------
% Vdat : 1 x #vertices struct with fields
%     vertxcoord
%     vertycoord
%     ncells
%     nverts
%     fourfold
% options : optional struct with fields
%   verbose : bool (default = false)
%       verbose output to Command Window
%   
% 
% Returns
% -------
% BL : #bonds x 2 int array
%   bond list indexing into 
%
% NPMitchell

% unpack options
verbose = true ;
preallocationFactor = 1.6 ;   % 3/2=1.5=cell-like preallocation factor, 
                              % 6/2=3=mesh-like (ie triangulation) prellocation
                              % Set to average connectivity or slightly higher
if nargin < 2
    options = struct() ;
end
if isfield(options, 'verbose')
    verbose = options.verbose ;
end
if isfield(options, 'preallocationFactor')
    preallocationFactor = options.preallocationFactor ;
end
    

% Consider each vertex
% Preallocate a reasonably long array of bonds
BL = zeros(round(preallocationFactor * length(Vdat)), 2) ;
dmyk = 1 ;
for qq = 1:length(Vdat)
    if verbose && (mod(qq, 50) == 1) 
        disp(['configuring bonds for vtx ' num2str(qq) '/' num2str(length(Vdat))])
    end
    
    for nn = Vdat(qq).nverts
        if nn < qq
            new = [nn, qq] ;
        elseif nn > qq
            new = [qq, nn] ;
        else
            error('how can a vertex be its own neighbor?')
        end

        % add the bond to BL
        if ~ismember(BL, new, 'rows')
            BL(dmyk, :) = new ;
            dmyk = dmyk + 1 ;
            % BL(1:dmyk, :) 
        end
        % else
        %     if verbose
        %         disp(['bond already in BL: [' num2str(new) ']'])
        %     end
        % end
    end
end

% Truncate
BL = BL(1:dmyk-1, :) ;

end