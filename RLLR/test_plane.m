clc;
%addpath('/Users/complexfilter/Documents/paper code');
options = [];
options.lamda = 0;
options.deta = 0;
options.gama = 0.1;
options.th = 8;
options.win = 5;
options.swin = 3;

basic_number = 1;
options.deta =0.01;
original_HRImg = double(imread('boat.bmp'));

LRImg= downsample_direct(original_HRImg);
options.original_HRImg = original_HRImg;

t0 = cputime; 
kernel = make_kernel(options.swin);
options.kernel = kernel / sum(sum(kernel));%normalization
options.weight_fg = 1;%%
options.our_method = double(bicubic(LRImg));
our_HRImg = RMLS_Interpolation_my_understanding(LRImg,options);

diff_img = our_HRImg-original_HRImg;
mse=sum(diff_img(:).*diff_img(:))/numel((original_HRImg));
our_PSNR=PSNR(our_HRImg,original_HRImg);
ssim = ssim_index(original_HRImg,our_HRImg);
FSIM = FeatureSIM(original_HRImg,our_HRImg);

time1=cputime-t0;

edge_PSNR = EdgePsnr(original_HRImg,our_HRImg);

RLLR_result = [our_PSNR,edge_PSNR,ssim,FSIM,time1,mse];

figure (1);imshow(our_HRImg/255);
PSNR(our_HRImg,original_HRImg)
PSNR(our_HRImg(21:end-20,21:end-20),original_HRImg(21:end-20,21:end-20))