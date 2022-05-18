function save_all_scientific_cmaps(outdir)
%  FIRST DOWNLOAD ScientificColourMaps FROM INTERNET THEN RUN THIS
%
% Example usage
% -------------
% save_all_pyplot_cmaps('/mnt/data/code/gut_matlab/plotting/PyColormap4Matlab/pyplot_colormaps/')

if nargin < 1
    outdir = './' ;
end

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

cmNames ={'davos', 'devon', 'fes', 'grayC', ...
    'hawaii', 'imola', 'acton', 'lajolla', ...
    'bam', 'lapaz', ...
    'bamO', 'lisbon', ...
    'bamako', 'nuuk', ...
    'batlow', 'oleron', ...
    'batlowK', 'oslo', ...
    'batlowW', 'roma', ...
    'berlin', 'romaO', ...
    'bilbao', ...
    'broc', 'tofino', ...
    'brocO', 'tokyo', ...
    'buda', 'turku', ...
    'bukavu', 'vanimo', ...
    'cork', 'vik', ...
    'corkO', 'vikO' } ;
for qq = 1:length(cmNames)
    cmap = load(fullfile(cmNames{qq}, [cmNames{qq} '.mat'])) ;

    eval(['cmap = cmap.' cmNames{qq}])
    
    % Save to function
    outfn = fullfile(outdir, [lower(cmNames{qq}) '.m']) ;
    if exist(outfn, 'file')
        yn = input('file already exists! overwrite?; ...s') ;
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












