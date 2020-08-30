clc;

options = [];
options.lamda = 0;
options.deta = 0;
options.gama = 0.1;
options.th = 8;
options.win = 5;
options.swin = 3;

basic_number = 1;

for im_num = 1:10
    
% if im_num == 1
%    options.deta =0.001;
%    original_HRImg = double(imread('HR_airplane.pgm'));
% end
% 
% if im_num == 2
%    options.deta =0.001;
%    original_HRImg = double(imread('HR_lena.pgm'));
% end
% 
% 
% if im_num == 3
%    options.deta =0.01;
%    original_HRImg = double(imread('HR_flowers.pgm'));
% end
% 
% if im_num == 4
%    options.deta =0.01;
%    original_HRImg = double(imread('HR_girl.pgm'));
% end
% 
% if im_num == 5
%    options.deta =0.001;
%    original_HRImg = double(imread('HR_door.pgm'));
% end
% 
% if im_num == 6
%    options.deta =0.001;
%    original_HRImg = double(imread('HR_peppers.pgm'));
% end

if im_num == 7
   options.deta =0.01;
   original_HRImg = double(imread('HR_splash.pgm'));
end

if im_num == 8
   options.deta =0.1;
   original_HRImg = double(imread('baboon256.pgm'));
end

if im_num == 9
   options.deta =0.01;
   original_HRImg = double(imread('HR_butterfly.pgm'));
end
% %%%%%%%%%%%%%%%%%%%%%%%%

% original_HRImg = original_HRImg(:,:,channel);
% LRImg= downsample_average(original_HRImg);
LRImg= downsample_direct(original_HRImg);


options.original_HRImg = original_HRImg;


% % %%%%%%%%%%%%%%%%RLLR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = cputime; 
kernel = make_kernel(options.swin);
options.kernel = kernel / sum(sum(kernel));%normalization

options.weight_fg = 1;%%
options.our_method = double(bicubic(LRImg));
our_HRImg = RMLS_Interpolation(LRImg,options);

diff_img = our_HRImg-original_HRImg;
mse=sum(diff_img(:).*diff_img(:))/prod(size(original_HRImg));
our_PSNR=10*log10(255^2/mse);
ssim = ssim_index(original_HRImg,our_HRImg);
FSIM = FeatureSIM(original_HRImg,our_HRImg);

time1=cputime-t0;

edge_PSNR = EdgePsnr(original_HRImg,our_HRImg);

RLLR_result = [our_PSNR,edge_PSNR,ssim,FSIM,time1,mse];

figure (1);imshow(our_HRImg/255);

% % % % 

end