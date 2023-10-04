function out = sprintfm(inputString, varargin)
% Wrapper for sprintf()

if ispc
    % change \ to \\ but keep \\
    % inputString = regexprep(inputString, '(?<!\\)\\(?!\\)', '\\\\');

    % change \ to \\ and \\ to \\\\, etc
    inputString = replace(inputString, '\', '\\') ;
end

out = sprintf(inputString, varargin{:}) ;

end

