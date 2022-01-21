% function alignMaskedDataAPDV(QS)
% Given training on midgut cells and spurious (amnioserosa) cells as 
% h5 files in dir16bit/stabilized_h5s/, mask the original data and save 3d
% volumes.
% Note: if QS.plotting.preview == True, displays intermediate results
%
% THIS IS AT PRESENT UNFINISHED 2020
%
% NPMitchell 2020

% Parameters
cliplowDorsal = 0.3 ;
step = 10 ;

% Naming
maskBase = [QS.fileBase.name '_Probabilities_mask3d.h5'] ;
dorsalBase = [QS.fileBase.name '_Probabilities_maskDorsal.h5'] ;
maskDataDir = fullfile(QS.dir.dataDir, 'masked_data') ;


% First create an average stab image to use for training
preview = QS.plotting.preview ;
for tt = QS.xp.fileMeta.timePoints
    fullFileName = fullfile(maskDataDir, ...
        sprintf([QS.fileBase.name '_masked.tif'], tt)) ;
    dat = QS.loadBioFormats(fullFileName) ;
    ivii = griddedInterpolant(dat) ;
    
end