function write_txt_with_header(fn, dat, header)
% WRITE_TXT_WITH_HEADER(FN, DAT, HEADER) Write data to a txt file with header
%
% Parameters
% ----------
% fn : char
%   the filename (full path) of the txt file
% dat : array
%   the data to write to disk
% header : char
%   the header string for the txt file
%
% NPMitchell 2019

fid = fopen(fn, 'wt');
% make header
fprintf(fid, [header '\n']);  
fclose(fid);
dlmwrite(fn, dat, '-append')

return