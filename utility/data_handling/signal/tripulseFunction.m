function tripKernel = tripulseFunction(halfwidth)
% tripulse(halfwidth)
% Build linear tripulse filter kernel with specified halfwidth.
%
%                .  1/integral(kernel)
%               / \
%              /   \
%             /     \
% 0 _________.       .___________ 0
%
% Parameters
% ----------
% halfwidth : int
%   half-width of the tripulse kernel. For ex, 2 gives a pulse of [0.33,
%   0.66, 1, 0.66, 0.33].
% 
% Returns
% -------
% tripKernel
%   tripulse filter kernal, rising from 0 to max to 0 with specified
%   halfwidth, and linear ramp. The max value is determined by the integral
%   so that the kernel is normalized.

if halfwidth == 0 
    tripKernel = 1 ;
elseif halfwidth > 0
    halfpulse = linspace(0, 1, halfwidth + 2) ;
    tripKernel = [ halfpulse(2:end), fliplr(halfpulse(2:end-1))] ;
    tripKernel = tripKernel ./ sum(tripKernel(:)) ;
else
    error('Specified halfwidth of tripulse kernel is negative. This is nonsensical.')
end
end