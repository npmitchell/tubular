function ww = tukeyFunction(xx, wRise, wFall)
% For some scalar signal between 0 and 1, define a tukey window response to
% that signal, which rises up as a cosine, plateaus, and falls as cosine
% 
% Parameters
% ----------
% xx : signal to bound with window
%   data to bound, can be any shape 
% wRise : numeric / float 
%   width of rise
% wFall : optional numeric, default = wRise
%   width of fall
% 
%
% NPMitchell 2021

ww = ones(size(xx)) ;

% rise / increase
x1 = (xx < wRise*0.5) ;
ww(x1) = 0.5 * ( 1 + cos(2*pi / (wRise*2) * (xx(x1) + wRise)) ) ;
ww(xx < 0) = 0 ;

% plateau -- already satisfied
% ww(in_plateau) = 1 ;

% decrease
x2 = (xx > 1 - wFall*0.5) ;
ww(x2) = 0.5 * ( 1 + cos(2*pi / (wFall*2) * (xx(x2)-1 + wFall)) ) ;
ww(xx > 1) = 0 ;

