% main program to demonstrate how NEDI works
% Author: Xin Li
% Reference: New edge-directed interpolation 
% Li, X and Orchard, M T 
% IEEE Transactions on Image Processing. Vol. 10, no. 10, pp. 1521-1527. Oct. 2001 

x=double(imread('fl_orig.ppm'));
imshow(x/255,[]);
title('original LR image');
for i=1:3
    y(:,:,i)=sri(x(:,:,i),1);
end

imshow(y/255,[]);
title('interpolated HR image');