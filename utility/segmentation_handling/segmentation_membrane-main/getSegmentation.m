function binary_map = getSegmentation(img, stepsize, windowsize)
%  getSegmentation returns the segmented image of a membrane image with default size 2000x2000, stepsize take an integer
%  and windowsize a 2-d array. Default values are 20 and [100,100]

if nargin == 1
    stepsize = 20;
    windowsize = [100,100];
end

canvas = zeros(2000,2000);
for i = 1: stepsize:2000-windowsize(2)-1
    for j = 500:stepsize:1500-windowsize(1)
        
        rectangle_position = [i,j,windowsize(1),windowsize(2)];
        cropped1 = imcrop(img, rectangle_position);
        
        
        
        cropped = adapthisteq(cropped1);
        I = max(cropped(:))-cropped;
        
        se = strel('disk',3);
        Ie = imerode(I,se);
        Iobr = imreconstruct(Ie,I);
        Iobrd = imdilate(Iobr,se);
        Iobrcbr = imreconstruct(Iobrd, Iobr);
        
        bw = imbinarize(Iobrcbr);
        DL = watershed(bwdist(bw));
        bgm = DL == 0;
        canvas(rectangle_position(2):rectangle_position(2)+rectangle_position(4), rectangle_position(1):rectangle_position(1)+rectangle_position(3))=canvas(rectangle_position(2):rectangle_position(2)+rectangle_position(4), rectangle_position(1):rectangle_position(1)+rectangle_position(3))+bgm;
        
        
    end
end

canvas(canvas ~=0) =1;
segimg = logical(canvas);
segimginverse = imcomplement(segimg);      % invert image
segimginverse = bwareaopen(segimginverse, 100, 4);
segimg = imcomplement(segimginverse);
segimg = bwskel(segimg);
segimg = bwmorph(segimg, 'spur', 1000);
segimg = bwareaopen(segimg, 2);
binary_map = segimg;

end