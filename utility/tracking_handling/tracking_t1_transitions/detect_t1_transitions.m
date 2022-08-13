subdir = '/home/yuzhenglin/membraneproject';
labelDir = fullfile(subdir, 'labeled_groundTruth') ;

timePoints = 96:2:206;


if ~exist("neighbors.mat", 'File')
pairs_cell = cell(numel(timePoints),1);
for tidx = 1: numel(timePoints)
    tp = timePoints(tidx);
    load_label = load(fullfile(labelDir, sprintf("tracks_label_%06d.mat", tp)));
    label = load_label.imlabel;
    membrane = find(label == -1);
    pairs = zeros (numel(membrane),2);
    jj = 0;
    for i = 1:numel(membrane)
        [c1, c2] = ind2sub(size(label),membrane(i));
        empty = zeros(size(label));
        empty(c1,c2) = 1;
        not_empty = bwmorph(empty, 'dilate');
        mem = label(not_empty);
        [id1, id2] = cellPair(mem);
        if ~isnan(id1)
            jj = jj+1; % counting
            pairs(jj, 1) = id1;
            pairs(jj, 2) = id2;
        end
        pairs = pairs(1:jj, :); % get rid of zeros
    end
    pairs = o_g_t(sortrows(pairs, 1), 3);
    pairs_cell(tidx) = {pairs};
end
save("neighbors.mat", 'pairs_cell');
else
%% 
load("neighbors.mat");
union_pair = double.empty(0,2);
    for tidx = 1: numel(timePoints)
        tp = timePoints(tidx);
        pairs = cell2mat(pairs_cell(tidx));
        union_pair = union(pairs, union_pair, 'row'); 
        
    end
    
    roll = ones(size(union_pair,1),numel(timePoints));
    for tidx = 1: numel(timePoints)
        tp = timePoints(tidx);
        pairs = cell2mat(pairs_cell(tidx));
        [cnp, ind_cnp] = setdiff(union_pair, pairs, "row"); % cell not present. cells not in the current frame 
        roll(ind_cnp ,tidx) = 0;
    end
    [not_always_neighbors, not_a_n_ind] = setdiff(roll, ones(1, numel(timePoints)), 'row');
    not_always_neighbors_label = union_pair(not_a_n_ind, :);
    save('not_always_neighbors.mat', 'not_always_neighbors', 'not_always_neighbors_label');
end
function unique_list = o_g_t(sorted_list, n) % occurrence greater than
    i = 0;
    unique_list = zeros(size(sorted_list));
    for jj = 1: size(sorted_list, 1)-n+1
        if sorted_list(jj, :) == sorted_list(jj+ n-1, :)
           i = i + 1;
           unique_list(i, :)= sorted_list(jj, :);
        end
    end
    unique_list = unique_list(1:i,:);
    unique_list = unique(unique_list, "row");
end

function [id1, id2] = cellPair(mem) % returns a ordered pair of cells, mem is a 3x3 double matrix (or any list of numbers) where -1 is at the center.
    mem_mem = unique(mem);
    mem_mem = mem_mem(mem_mem ~= -1); %membrane
    mem_mem = mem_mem(mem_mem ~= 0); % empty space
    
    if numel(mem_mem) == 2
        id1 = min(mem_mem);
        id2 = max(mem_mem);
    else
        id1 = nan;
        id2 = nan;
    end
end