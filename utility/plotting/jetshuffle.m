function jetshuffle = jetshuffle(ncolors)
% Colormap of jet but with colors randomly shuffled
%
% NPMitchell 2021

if nargin < 1
    ncolors = 256 ;
end
jets = jet(ncolors) ;
jetshuffle = jets(randperm(ncolors), :) ;
    