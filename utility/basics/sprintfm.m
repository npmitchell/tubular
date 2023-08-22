function out = sprintfm(inputString, number)

out = regexprep(inputString, '(?<!\\)\\(?!\\)', '\\\\');
out = sprintf(out, number) ;

end

