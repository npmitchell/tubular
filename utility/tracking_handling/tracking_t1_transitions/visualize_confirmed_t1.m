subdir = '/home/yuzhenglin/membraneproject';
labelDir = fullfile(subdir, 'labeled_groundTruth');
imOutDir = fullfile(subdir, 'potential_t1_transitions') ;
imgdir = '/caaxdataset/deconvolved_16bit/msls_output/gridCoords_nU0100_nV0100/PullbackImages_010step_sphi/smoothed_extended/';
rawImFileBase =  fullfile(subdir, imgdir, 'Time_%06d_c1_stab_pbspsme.tif');
timePoints = 96:2:206;
load('not_always_neighbors.mat');
load('confirmed_t1.mat');


mergers = unique(merge_events(:,1));
splitters = unique(split_events(:,1));

for i = 1: size(merge_events,1)
    display_potential_t1(merge_events(i, 1), merge_events(i,2)-4:merge_events(i, 2)+8);
end

