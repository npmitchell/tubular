function bnds = polygonsFromSegmentation(labelMat) 
% Convert a label matrix into a set of non-overlapping polygons.
% Tested on a labelMatrix made from a skeleton with connectivity=4.
%
% Parameters
% ----------
% labelMat : 2D int array
%   label matrix like that made from labelmatrix(skel)
% 
% Returns
% -------
% bnds : #polygons x 1 cell array of #vertices x 2 ints
%   the boundary locations of the polygons in the labelMatrix
% 
% NPMitchell 2021


se = strel('disk', 1) ;
tmp = se.Neighborhood ;
tmp(1, :) = 0 ;
tmp(2, 3) = 0 ;
tmp(3, 2) = 0 ;
se2 = strel('arbitrary', tmp)  ;
bnds = cell(max(labelMat(:)), 1) ;
for obj = 1:max(labelMat(:))
    if mod(obj, 10) == 0 
        disp(['boundaries for object ' num2str(obj)])
    end
    bw = labelMat == obj ;
    % bnd = bwboundaries(imdilate(imdilate(bw, se), se2), 4);
    bnd = bwboundaries(imdilate(bw, se), 8);
    if length(bnd) == 1
        bnds{obj} = bnd{1} ;
            
        % check it
        % plot(bnd{1}(:, 1), bnd{1}(:, 2), '-')
        % hold on;
        % pause(0.000001)
    end
    
end

