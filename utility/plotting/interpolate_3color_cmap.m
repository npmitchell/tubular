function[map] = interpolate_3color_cmap(s,rgb1, rgb2, rgb3, normalize)
%This function is based on Kenneth Moreland's code for greating Diverging
% Colormaps, modified by Andy Stein. Created by Noah P Mitchell 2019
% template accessed from: https://www.kennethmoreland.com/color-maps/
%s is a vector that goes between zero and one 
% Common rgb1, rgb2 pairs are (rgb1, rgb2), (rgb3, rgb4), etc with these
% each defined below
% 
% Parameters
% ----------
% s : array with values spanning [0, 1]
%   The values to be mapped to colors, spanning [0, 1]
% rgb1 : length(3) array or int
%   An RGB triplet for value 0 or an integer indexing common rgb triplets
% rgb2 : length(3) array or int
%   An RGB triplet for value 1 or an integer indexing common rgb triplets
%
% Output
% ------
% map : colormap
%   A mapping from [0, 1] to colors spanning rgb1 and rgb2, with length(s)
%   possible values.
% 
% Example Usage
% -------------
% % Blue to red:
% bwr = diverging_cmap([0:0.01:1], 1, 2)
% % Black to white:
% cmap2 = diverging_cmap([0:0.01:1], [0,0,0], [1,1,1])

% If the arguments rgb1 and rgb2 are ints, then use these arguments 
% to index from a list of commonly used rgb values
rgb_1 = [0.230, 0.299, 0.754] ;  % blue
rgb_2 = [0.706, 0.016, 0.150] ;  % red
rgb_3 = [0.436, 0.308, 0.631] ;  % dark purple
rgb_4 = [0.759, 0.334, 0.046] ;  % off red
rgb_5 = [0.217, 0.525, 0.910] ;  % light blue
rgbs = [rgb_1; rgb_2; rgb_3; rgb_4; rgb_5] ;
if all(size(rgb1) == 1)
    rgb1 = rgbs(rgb1, :) ;
end
if all(size(rgb2) == 1)
    rgb2 = rgbs(rgb2, :) ;
end
if all(size(rgb3) == 1)
    rgb3 = rgbs(rgb3, :) ;
end

map = zeros(length(s),3);
for i=1:length(s)
    map(i,:) = threecolor_map_1val(s(i),rgb1,rgb2, rgb3, normalize);
end
end

% Interpolate a diverging color map.
    function[result] = threecolor_map_1val(s, rgb1, rgb2, rgb3, normalize)
        % s1 is a number between 0 and 1 tuning through the colormap
        % rgb1,2,3 are RGB triplets
        
        % First color goes from 1 to zero in (0, 0.5)
        s1 = max(0, 1 - 2 * s) ;
        % First color goes from 1 to zero in (0, 0.5)
        s2 = 1 - abs(2 * (s - 0.5)) ;
        % Third color goes from zero to one in (0.5, 1)
        s3 = max(2 * s - 1, 0) ;
        
        result(1) = s1*rgb1(1) + s2*rgb2(1) + s3*rgb3(1);
        result(2) = s1*rgb1(2) + s2*rgb2(2) + s3*rgb3(2);
        result(3) = s1*rgb1(3) + s2*rgb2(3) + s3*rgb3(3);
        
        % Renormalize the result
        if normalize
            result = result / vecnorm(result) ;
        end
        xyz = RGBToXYZ(result) ;    
        result = XYZToRGB(xyz);
    end
