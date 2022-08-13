% make mips of data (for example, to see saggital plane of data over time)
% using tiffMIP
%
% See also makeSubstackMIPs()

addpath('/mnt/data/code/gut_matlab/tiff_handling')
slices = [] ;
outDir = './mipPreviews/' ;
for ch = [1,2]
    for tt = 0:10:110
        fileName = sprintf('./TP%d_Ch%d_Ill0_Ang0,60,120,180,240,300.tif', ...
            tt, ch) ;
        
        disp([ 'reading ' fileName]) 
        
        mipfns = {fullfile(outDir, sprintf('mip_t%04d_c%d_slice1.png', tt, ch)),...
            fullfile(outDir, sprintf('mip_t%04d_c%d_slice2.png', tt, ch)),...
            fullfile(outDir, sprintf('mip_t%04d_c%d_slice3.png', tt, ch))};
        
        if ~exist(mipfns{1}, 'file') && ~exist(mipfns{2}, 'file') ...
                && ~exist(mipfns{3}, 'file')
            mip = tiffMIP(fileName, mipfns, 'middle', 0.5, [1,2,3]) ;
        end
    end
end

%% 2 color overlay
addpath('/mnt/data/code/gut_matlab/tiff_handling')
% fileBase = './TP%d_Ch%d_Ill0_Ang0,60,120,180,240,300.tif';
fileBase = './Time_%06d_Angle_0_c%d_ls_1.ome.tif';
slices = [] ;
outDir = './mipPreviews/' ;
scales = [0.2, 0.005] ;
chs = [1, 2] ;
dims = [2,3]; % [1,2,3]

if ~exist(outDir, 'dir')
    mkdir(outDir)
end

for tt = 0:110
    mips = {} ;
    % Read each channel and store the mips along each dim
    for ch = chs
        fileName = sprintf(fileBase, ...
            tt, ch) ;
        mipfns = {} ;               
        mips{ch} = tiffMIP(fileName, mipfns, 'middle', scales(ch), dims) ;
        
    end

    % Overlay the two colors
    for qq = 1:length(mips{ch})
        mipfn = fullfile(outDir, sprintf('mip_t%04d_overlay_slice%d.png', tt, qq)) ;
        imRGB = cat(3, mips{chs(2)}{qq}, mips{chs(1)}{qq}, mips{chs(1)}{qq}) ;
        imwrite(imRGB, mipfn) 
    end
end


%% 2 color overlay -- RAW
addpath('/mnt/data/code/gut_matlab/tiff_handling')
% fileBase = './TP%d_Ch%d_Ill0_Ang0,60,120,180,240,300.tif';
fileBase = './Time_%06d_Angle_0_c%d_ls_1.ome.tif';
slices = [] ;
outDir = './mipPreviews/' ;
scales = [1000, 200] ;
chs = [1, 2] ;
dims = [3]; % [1,2,3]

if ~exist(outDir, 'dir')
    mkdir(outDir)
end

for tt = 0:99
    mips = {} ;
    % Read each channel and store the mips along each dim
    for ch = chs
        fileName = sprintf(fileBase, ...
            tt, ch) ;
        mipfns = {} ;               
        mips{ch} = tiffMIP(fileName, mipfns, 'middle', scales(ch), dims) ;
        
    end

    % Overlay the two colors
    % for qq = 1:length(mips{ch})
        mipfn = fullfile(outDir, sprintf('mip_t%04d_overlay_slice%d.png', tt, qq)) ;
        imRGB = cat(3, mips{chs(2)}, mips{chs(1)}, mips{chs(1)}) ;
        imwrite(imRGB, mipfn) 
    % end
end
