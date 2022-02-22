function [allVtx, allBonds] = nearestNeighborsAndBondsWithinKSteps(k, seeds, eIDx)
% Find all vertices and bonds that are within k steps of seed vertices.
%
% Parameters
% ----------
% k : int >= 0
%   How many steps to search
% seeds : ints
%   seed vertex indices for which to find Neighbors and Bonds within K
%   steps
% eIDx : #edges x 2 int array
%   edge pair indices, a 'bond list'
% 
% I'm sure this could be sped up with some matlab graph built-ins
% 
% NPMitchell 2021

allVtx = seeds(:) ;
allBonds = [] ;
if k == 0
    for qq = 1:length(seeds) 
        fixIdx = seeds(qq) ;
        [rows, ~ ] = find(eIDx == fixIdx) ;
        toadd = eIDx(rows, :) ;
        allBonds = unique([allBonds; rows]) ;
        allVtx = unique([allVtx; toadd(:) ]); 
    end
else
    for pp = 1:k
        for qq = 1:length(seeds) 
            fixIdx = seeds(qq) ;
            [rows, ~ ] = find(eIDx == fixIdx) ;
            toadd = eIDx(rows, :) ;
            allBonds = unique([allBonds; rows]) ;
            allVtx = unique([allVtx; toadd(:) ]); 
        end
        fixedIDx = unique(allVtx) ;

        % Prepare for next step by updating which seeds to cycle over
        seeds = fixedIDx ;
    end
end