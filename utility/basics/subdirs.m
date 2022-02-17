function [fullPaths, subNames] = subdirs(directory)

items = dir(directory) ;
dirFlags = [items.isdir] ;
% Extract only those that are directories.
out = items(dirFlags) ;

% Prepare for output as cell array(s)
dmyk = 1 ;
fullPaths = {} ;
if nargout > 1
    subNames = {} ;
end
% Populate cell arrays
for ii = 1:length(out)
    % Ingore hidden paths
    if ~strcmp(out(ii).name, '.') && ~strcmp(out(ii).name, '..')
        if nargout > 1
            subNames{dmyk} = out(ii).name ;
        end    
        fullPaths{dmyk} = fullfile(out(ii).folder, out(ii).name) ;
        dmyk = dmyk + 1 ;
    end
end