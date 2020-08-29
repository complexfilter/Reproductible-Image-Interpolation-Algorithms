close all;clear;clc;
f = imread('man.tiff');
f = double(f(1:1024, 1:1024));
imH = f;

% Downsample the high-resolution by a factor of 2x2. 
imL = imH(1:2:end, 1:2:end);

% Super-resolution image zooming with mixing sparse estimators. 
% By defaut, no extrapolation.
option.extrapol = 0;
% Defaut parameter of orthogonal block matching pursuit.  
option.lambda = 0.6;
option.T = 0;
[imH_SR, imH_linear] = SME_ImageZoom(imL, option);

% Display the results. 
figure
imagesc(imL)
truesize
colormap('gray')
title(' LR Inputimage')

figure
imagesc(imH)
truesize
colormap('gray')
title('HR ground truth')

figure
imagesc(imH_linear)
truesize
colormap('gray')
if ( size(imH, 1) ~= size(imH_linear, 1) ) 
    PSNR = get_PSNR(imH(1:end-1, 1:end-1), imH_linear);
else
    PSNR = get_PSNR(imH, imH_linear);
end
tmp = sprintf('Cubic spline: %.2f dB', PSNR);
title(tmp)
disp(tmp)

figure
imagesc(imH_SR)
colormap('gray')
truesize
if ( size(imH, 1) ~= size(imH_linear, 1) ) 
    PSNR = get_PSNR(imH(1:end-1, 1:end-1), imH_SR);
else
    PSNR = get_PSNR(imH, imH_SR);
end
tmp = sprintf('SME supre-resolution %.2f dB', PSNR);
title(tmp)
disp(tmp)
imH_H=zeros(1024);
imH_H(1:1023,1:1023)=imH_SR(1:1023,1:1023);
load region_ind_edge_20180113_1024x1024.mat
region_ind=region_ind_edge_20180113_1024x1024;
l=length(region_ind);
vv4=sqrt(1/l*sum(sum((double(imH(region_ind)-imH_H(region_ind)).^2))));
PSNR4=20*log10(255/vv4);
l=numel(imH(234:324,449:539));
vv5=sqrt(1/l*sum(sum((double(imH(234:324,449:539)-imH_H(234:324,449:539)).^2))));
PSNR5=20*log10(255/vv5);



