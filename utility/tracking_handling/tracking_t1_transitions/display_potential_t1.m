function display_potential_t1(pair_id1, pair_id2, time_id, savefile)

subdir = '/home/yuzhenglin/membraneproject';
labelDir = fullfile(subdir, 'labeled_groundTruth');
imOutDir = fullfile(subdir, 'confirmed_t1_transitions') ;
imgdir = '/caaxdataset/deconvolved_16bit/msls_output/gridCoords_nU0100_nV0100/PullbackImages_010step_sphi/smoothed_extended/';
rawImFileBase =  fullfile(subdir, imgdir, 'Time_%06d_c1_stab_pbspsme.tif');
timePoints = 96:2:206;
load('not_always_neighbors.mat');
zoom = true;


ii = pair_id1;
iii = pair_id2;
outName = strcat("pairId_", sprintf('%06d', ii), "_timeID_",sprintf('%06d', round((time_id(1) + time_id(end))/2)), ".tif");
disp(['considering i = ' num2str(ii)]);        
flag = true;
for tidx = intersect(time_id, 1:56)
    tp = timePoints(tidx);
    disp(['considering t = ' num2str(tp)]);
    load_label = load(fullfile(labelDir, sprintf("tracks_label_%06d.mat", tp)));
    label = load_label.imlabel;
    rawImg = imread(sprintf(rawImFileBase, timePoints(tidx)));
    %if ismember(not_always_neighbors_label(ii,:), label)
        cellpair1 = not_always_neighbors_label(ii,:);
        cellpair2 = not_always_neighbors_label(iii,:);
        cells = zeros(size(label));
        
        cells(label == cellpair1(1)) = 4;
        cells(label == cellpair1(2)) = 2;
        cells(label == cellpair2(1)) = 5;
        cells(label == cellpair2(2)) = 8;
        colored_im = labeloverlay( rawImg,cells,'Transparency',0.8);
        regp = regionprops(cells,'Centroid');
        ctd = vertcat(regp.Centroid);

        if flag == true % one off assignment
        window_y = max(floor(min(ctd(:,1)))-100, 1):min(floor(max(ctd(:,1)))+100, 2000);
        window_x = max(floor(min(ctd(:,2)))-100, 1):min(floor(max(ctd(:,2)))+100, 2000);
        flag = false;
        end
        %text(ctd(:,1),ctd(:,2), num2str(ii));
    %end
    if zoom == true
         colored_im = imresize(colored_im(window_x, window_y, :),5);
    end
        if savefile
            imwrite(colored_im, fullfile(imOutDir,outName), 'WriteMode', 'append', 'Compression','none');
        else
            imshow(colored_im);
            pause(0.1);
        end
end



% sort nanl such that the nanl is partitioned into disjoint parts, by which
% I mean each part has only unique elements. e.g. [1,2;3,4;5,6;1,3] rather than
% [1,2;1,3;3,4;5,6]
% 
% function reordered_list = disjoint_labels(nanl)
%     counter_list = [];
%     flag = false; 
%     counter1 = 1; % slow, traverse the list once
%     counter2 = 1; % fast, traverse the list multiple times.
%     while flag == false
%         if numel(setdiff(nanl(counter2, :), nanl(1:counter1, :))) == 2
%             % swap
%             counter1 = counter1 + 1;
%             counter2 = counter2 + 1;
%         else
%             counter2 = counter2+1;
%         end
%         
%     end
%         
% 
% end

end