seg_dir1 = '/first/pass/segmentation';
image_dir  = '/membrane/image/data/';
window = nan(2000,2000,5);
i = 0;
for t = 93:97
    
segname = strcat(seg_dir1, 'Time_', sprintf('%06d', t), '_c1_stab_segmentation_binary_map.png'); 
segimg = imread(segname);
window(:,:,5) = segimg;
if i < 4
else
    bigseg = sum(window, 3);
    bigseg(bigseg ~= 0) = 1;
    smoothed_seg = imgaussfilt(bigseg,4);
    img = imread(strcat(image_dir, 'Time_', sprintf('%06d', t-2), '_c1_stab_pbkspsme.tif'));
    %adjusted = smoothed_seg-im2double(imgaussfilt(img,5));
    %adjusted(adjusted<0)=0;
    %enhanced = adjusted + im2double(img);
    enhanced = smoothed_seg+im2double(img);
    imwrite(enhanced, strcat('/saved/data/dir/Time_', sprintf('%06d', t-2), '_c1_stab_pbspsme_enhanced.tif'))
    %imshowpair(enhanced, img, 'montage')
    %imshowpair(field, img, 'montage')
    %imshow(imoverlay(img, bigseg))
%imshow(imresize(bp(520:1000,200:650),2));
%pause(1)
end
window(:,:,1) = window(:,:,2);
window(:,:,2) = window(:,:,3);
window(:,:,3) = window(:,:,4);
window(:,:,4) = window(:,:,5);
i = i+1;

end
