function LRImg = downsample_avergae_pgm(HRImg)
[height,width,channel] = size(HRImg);
figure (1); imshow(HRImg/255);
title('HR image');

PSF = fspecial('gaussian',5,0.8);
Blurred = imfilter(HRImg,PSF,'symmetric','conv');
% figure (2); imshow(Blurred/255);title('Blurred Image');

if channel == 0
   LRImg = Blurred(1:2:height,1:2:width);
else
   LRImg = Blurred(1:2:height,1:2:width,channel);  
end
%%%%%%%%%%%%%%%%%low-pass filter%%%%%%%%%%%%%%
%  I = imread('board.tif');
% figure;imshow(I);title('Original Image');



% V = .02;
% BlurredNoisy = imnoise(Blurred,'gaussian',0,V);
% figure;imshow(BlurredNoisy);title('Blurred and Noisy Image');

figure (2); imshow(LRImg/255);
[height,width,channel] = size(LRImg);
% imwrite(uint8(LRImg), 'E:\xmliu\picture_base\graypictures\LR_average_baboon256.bmp', 'bmp' );
title('LR img');