function save_all_pyplot_cmaps(outdir)
% 
% Example usage
% -------------
% save_all_pyplot_cmaps('/mnt/data/code/gut_matlab/plotting/PyColormap4Matlab/pyplot_colormaps/')

if nargin < 1
    outdir = './pyplot_colormaps/' ;
end

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

cmNames = getPyPlot_cMap('!GetNames');
for qq = 1:length(cmNames)
    cmap = getPyPlot_cMap(cmNames{qq}, 256);

    % Save to function
    outfn = fullfile(outdir, [lower(cmNames{qq}) '.m']) ;
    if exist(outfn, 'file')
        yn = input('file already exists! overwrite?', 's') ;
        if isempty(yn) || contains(yn, 'n') || contains(yn, 'N')
            error('exiting to avoid overwrite')
        end
    end
    
    fileID = fopen(outfn,'w');
    fprintf(fileID, ['function ' lower(cmNames{qq}) ' = ' lower(cmNames{qq}) '()\n']);
    fprintf(fileID, [ lower(cmNames{qq}) ' = [ ... \n']) ;
    for row = 1:length(cmap)-1
        fprintf(fileID,'       %d %d %d \n', cmap(row,1), cmap(row,2), cmap(row,3)) ;
    end
    fprintf(fileID, '       %d %d %d ] ; \n', cmap(end,1), cmap(end,2), cmap(end,3)) ;
    fclose(fileID);
    
    % Save to mat?
    % save(fullfile(outdir, cmNames{qq}), 'cmap') ;
end












