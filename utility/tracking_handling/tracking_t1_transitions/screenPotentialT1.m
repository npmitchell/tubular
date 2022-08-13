subdir = '/home/yuzhenglin/membraneproject';
labelDir = fullfile(subdir, 'labeled_groundTruth');
imOutDir = fullfile(subdir, 'potential_t1_transitions') ;
imgdir = '/caaxdataset/deconvolved_16bit/msls_output/gridCoords_nU0100_nV0100/PullbackImages_010step_sphi/smoothed_extended/';
rawImFileBase =  fullfile(subdir, imgdir, 'Time_%06d_c1_stab_pbspsme.tif');
timePoints = 96:2:206;
load('not_always_neighbors.mat');
try
    load('continue_from.mat')
    start_from_id = pairid;
    start_from_tidx = tid;
catch
    start_from_id = 1;
    
    start_from_tidx = 1;
end
split = [1,1,0,0];
merge = [0,0,1,1];
try
    load('confirmed_t1.mat')
catch
    angle_split = [];
    angle_merge = [];
    split_events = [];
    merge_events = [];
end
for ii = start_from_id: size(not_always_neighbors_label,1)
    pairid = ii;
    %outName = strcat("pairId_", sprintf('%06d', ii), ".tif");
    disp(['considering i = ' num2str(ii)]);
    cellpair = not_always_neighbors_label(ii,:);
    for tidx = start_from_tidx:numel(timePoints)-(length(split)-1)
        tid = tidx;
        tp = timePoints(tidx);
        isSplit =  isequal(not_always_neighbors(ii, tidx:tidx + length(split)-1), split);
        isMerge =  isequal(not_always_neighbors(ii, tidx:tidx + length(split)-1), merge);
        if isSplit || isMerge
            load_label = load(fullfile(labelDir, sprintf("tracks_label_%06d.mat", tp)));
            label = load_label.imlabel;
            cells = zeros(size(label));
            cells(label == cellpair(1)) = 1;
            cells(label == cellpair(2)) = 2;
            regp = regionprops(cells,'Centroid');
            ctd = vertcat(regp.Centroid);
            window_y = max(floor(min(ctd(:,1)))-100, 1):min(floor(max(ctd(:,1)))+100, 2000);
            window_x = max(floor(min(ctd(:,2)))-100, 1):min(floor(max(ctd(:,2)))+100, 2000);
            flag = true;
            tt = 0;
            while flag == true
                load_label = load(fullfile(labelDir, sprintf("tracks_label_%06d.mat", timePoints(tidx + tt))));
                    rawImg = imread(sprintf(rawImFileBase, timePoints(tidx +tt)));
                    
                    label = load_label.imlabel;
                    cells = zeros(size(label));
                    cells(label == cellpair(1)) = 4;
                    cells(label == cellpair(2)) = 2;
                    composeimg = labeloverlay( rawImg,cells,'Transparency',0.8);
                    imshow(imresize(composeimg(window_x, window_y, :),5));
                    %title(strcat("current frame: ",num2str(tp)));
                k = waitforbuttonpress;
                %n: 110, y: 121, -> 29, <-28 ^ 30
                value = double(get(gcf,'CurrentCharacter'));
                try
                switch value
                    case 28
                        if tidx + tt >1
                            tt = tt-1;
                        end
                    case 29
                        if tidx + tt <numel(timePoints)
                            tt = tt + 1;
                        end
                    case 30
                        tt = 0;
                    case 110
                        flag = false;
                        close(gcf);
                    case 121
                        flag = false;
                        close(gcf);
                        if isSplit
                            angle_split = [angle_split, getOrientation(ctd)];
                            split_events = [split_events; ii, tidx]; % pair id, timepoint.
                        elseif isMerge
                            angle_merge = [angle_merge, getOrientation(ctd)];
                            merge_events = [merge_events; ii, tidx];
                        end
                    otherwise
                end
                catch
                    
                end
            end
        end
    save('confirmed_t1.mat', 'angle_split', 'angle_merge', 'split_events', 'merge_events');
    save('continue_from.mat', 'pairid', 'tid');
    end
    start_from_tidx = 1;
end

subplot(1,2,1)
histogram(angle_split,10)
axis([0, pi, 0, 30])
title('split angle')
set(gca,'XTick',0:pi/4:pi)
set(gca,'XTickLabel',{'0','pi/4','pi/2','3pi/4','pi'})
subplot(1,2,2)
histogram(angle_merge,10)
axis([0, pi, 0, 30])
title('merge angle')
set(gca,'XTick',0:pi/4:pi)
set(gca,'XTickLabel',{'0','pi/4','pi/2','3pi/4','pi'})

function angle = getOrientation(ctd)
try
    delta_x = ctd(1,1) - ctd(2,1);
    delta_y = ctd(1,2) - ctd(2,2);
    angle = atan(delta_y/delta_x);
    angle = mod(angle, pi);
catch
    angle = nan();
end
end
