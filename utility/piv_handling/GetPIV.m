function [VX,VY, x, y, shiftx, shifty] = GetPIV(im1,im2,X1,Y1,EdgeLength, histequilize)
    %   GetPIV computes PIV flow field estimate based on image1 (im1) and
    %   image2 (im2) using the Phase Correlation method implemented in xcorr2fft. 
    %   If im1 & im2 do not have the same dimensions. 
    %   If supplied, the grid X1,Y1 determines interpolated evaluation
    %   points. Note that these are not the true PIV evaluation points.
    %   is assmued to be contained in the image domain with finite
    %   EdgeLength defining size of PIV box. 
    %   Output are components of the flow field on the grid in interpolated
    %   and raw coords.
    %   
    %   Written by: Sebastian J Streichan, KITP, February 01, 2013
    %   NPM added histeq step (optional) and modified I/O and argparsing
    %
    % Parameters
    % ----------
    % im1 : NxM numeric
    %   image to compare to im2, to find PIV taking v: im1->im2
    % im2 : PxQ numeric
    %   image compared to im1, to find PIV taking v: im1-> im2. This image
    %   is resized to MxN if P~=M and Q~=N
    % X1 : RxS numeric (optional, if empty use evaluation coords)
    %   grid of x coordinates onto which we interpolate the PIV, forming a
    %   meshgrid with Y1
    % Y1 : RxS numeric (optional, if empty use evaluation coords)
    %   grid of y coordinates onto which we interpolate the PIV, forming a
    %   meshgrid with X1
    % Edgelength : numeric
    %   edgelength of PIV interrogation window (extent in both x and y)
    % histequilize : bool
    %   perform histogram equilization on the image before computing PIV.
    %   This equilizes brightness and contrast across the image in 20x20
    %   subimage bins.
    % 
    % Returns
    % -------
    % VX, VY : RxS double arrays
    %   piv velocities interpolated onto X1 and Y1 supplied vals
    % x, y : U*V x 1 double arrays
    %   linear arrays of evaluation coordinates
    % shiftx, shifty : RxS 
    %   piv velocities at evaluation coords x and y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Histogram equilization
    if nargin > 5 
       if histequilize
           im1 = histeq(im1) ;
           im2 = histeq(im2) ;
       end
    end
    
    % Resize im2 if not equal to im1 size
    if ~all(size(im1) == size(im2)) 
        im2 = imresize(im2, size(im1),'bicubic'); % rescale image if desired
    end
    
    si      = size(im1);
    nb      = ceil(si/EdgeLength);
    counter = 0;
    shiftx  = zeros(1,nb(1)*nb(2));
    shifty  = shiftx;
    x       = shiftx;
    y       = shifty;
    for i = 1 : nb(1)
        for j = 1 : nb(2)

            counter = counter+1;
            temp1 = im1(max(1,((i-2)*EdgeLength+1)):min((i+1)*EdgeLength,si(1)),(max(1,(j-2)*EdgeLength+1)):min((j+1)*EdgeLength,si(2)));
            temp2 = im2(max(1,((i-2)*EdgeLength+1)):min((i+1)*EdgeLength,si(1)),(max(1,(j-2)*EdgeLength+1)):min((j+1)*EdgeLength,si(2)));    
            
            [shiftx(counter),shifty(counter)] = xcorr2fft(temp1,temp2);
            x(counter) = (i-.5)*EdgeLength;
            y(counter) = (j-.5)*EdgeLength;
        end
    end
    
    % Put the flow field on a grid 
    if ~isempty(X1) 
        VX  = griddata(x,y,shiftx,X1,Y1);
        VY  = griddata(x,y,shifty,X1,Y1);
    else
        VX = reshape(shiftx, si) ;
        VY = reshape(shifty, si) ;
    end
end
