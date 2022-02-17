function dirpath = ensureDir(dirpath)
%ENSUREDIR(dirpath) create a directory if it does not exist already
%
% Parameters
% ----------
% dirpath : str
%   The path to create if it does not exist
%
% Returns
% -------
% dirpath : str
%   Returns the input, for convenience
% 
% NPMitchell 2019

if ~exist(dirpath, 'dir')
    mkdir(dirpath)
end