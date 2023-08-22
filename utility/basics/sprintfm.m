function out = sprintfm(inputString, varargin)
% Wrapper for sprintf()

if ispc
    inputString = regexprep(inputString, '(?<!\\)\\(?!\\)', '\\\\');
end

out = sprintf(inputString, varargin{:}) ;

end

