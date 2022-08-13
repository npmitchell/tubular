function subimages = readMultipageTiff(filename)
% NPM: I recommend using readTiff4D instead of this.
%
% Read a multipage tiff, assuming each file is the same size
    t = Tiff(filename, 'r');
    subimages(:,:,1) = t.read(); % Read the first image to get the array dimensions correct.
    if t.lastDirectory()
         return; % If the file only contains one page, we do not need to continue.
    end
    % Read all remaining pages (directories) in the file
    t.nextDirectory();
    while true
        subimages(:,:,end+1) = t.read();
        if t.lastDirectory()
            break;
        else
            t.nextDirectory();
        end
    end
end