load("confirmed_t1.mat")
subdir = '/home/yuzhenglin/membraneproject';
labelDir = fullfile(subdir, 'labeled_groundTruth');
imOutDir = fullfile(subdir, 'potential_t1_transitions') ;
imgdir = '/caaxdataset/deconvolved_16bit/msls_output/gridCoords_nU0100_nV0100/PullbackImages_010step_sphi/smoothed_extended/';
rawImFileBase =  fullfile(subdir, imgdir, 'Time_%06d_c1_stab_pbspsme.tif');
timePoints = 96:2:206;
load('not_always_neighbors.mat');
load('neighbors.mat')

suspects = [];
for i = 1: size(split_events, 1)
    sus_merge = find(merge_events(:,2) < split_events(i, 2)+4 & merge_events(:,2) > split_events(i, 2)-4);
    for j = 1:length(sus_merge)
        if split_events(i,1) ~=  merge_events(sus_merge(j),1)
    suspects = [suspects; [split_events(i,1), merge_events(sus_merge(j),1), split_events(i,2)]];
        end
    %load_label = load(fullfile(labelDir, sprintf("tracks_label_%06d.mat", tp)));
    %label = load_label.imlabel;
    end
end

%%
for ii = 1: size(suspects, 1)
    pairid1 = suspects(ii, 1);
    pairid2 = suspects(ii, 2);
    cellpair1 = not_always_neighbors_label(pairid1,:);
    cellpair2 = not_always_neighbors_label(pairid2,:);
    tidx = suspects(ii, 3);
    neighbors = cell2mat(pairs_cell(tidx)); 
    if ismember([cellpair1(1), cellpair2(1)],neighbors, 'rows') ||...
            ismember([cellpair1(1), cellpair2(2)],neighbors, 'rows') ||...
            ismember([cellpair1(2), cellpair2(1)],neighbors, 'rows') ||...
            ismember([cellpair1(2), cellpair2(2)],neighbors, 'rows')
        
    else
        suspects(ii, :) = nan(1,3);
    end
    
end
%%
suspects(isnan(suspects(:,1)), :) = [];

%%
suspects(isnan(suspects(:,1)), :) = [];

for i = 1: size(suspects,1)
    pair_id1 = suspects(i, 1);
    pair_id2 = suspects(i,2);
    time_id = suspects(i,3)+2;
    flag = true;
    savefile = false;
    while flag
        display_potential_t1(pair_id1, pair_id2, time_id, savefile);
        k = waitforbuttonpress;
        %n: 110, y: 121, -> 29, <-28 ^ 30
        value = double(get(gcf,'CurrentCharacter'));
        switch value
            case 110 % no
                flag = false;
                close(gcf);
                suspects(i, :) = nan(1,3);
                suspects
            case 121 % yes
                flag = false;
                close(gcf);
            otherwise
        end
    end
end
save("t1_overlap.mat", "suspects")

%%

display_potential_t1(pair_id1, pair_id2, time_id, savefile);
%%
for i = 6%1: size(suspects,1)
    pair_id1 = suspects(i, 1);
    pair_id2 = suspects(i,2);
    time_id = suspects(i,3)-1:suspects(i,3)+4;
    savefile = false;
    display_potential_t1(pair_id1, pair_id2, time_id, savefile);
end