function all_method(LRImg,original_HRImg,im_num)
%% LAZA
method_name = 'LAZA_';
pic_name = 'pic_name';
name = [LAZA pic_name]
LAZA_HRImg = double(LAZA(LRImg));
diff_img = LAZA_HRImg-original_HRImg;

mse=sum(diff_img(:).*diff_img(:))/prod(size(original_HRImg));
LAZA_PSNR=10*log10(255^2/mse);
ssim = ssim_index(original_HRImg,LAZA_HRImg);

figure (1+(im_num-1)*5);imshow(LAZA_HRImg/255);
figure (2+(im_num-1)*5);imshow(diff_img/255);
title(sprintf('LAZA HR image, PSNR = %2.2f dB, SSIM = %2.2f', LAZA_PSNR,ssim));
imwrite(uint8(LAZA_HRImg), 'E:\xmliu\result\LAZA_butterfly.bmp', 'bmp' );
% % 
% % %%%%%%%%%%%%%%%%RLLR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.weight_fg = 0;%%´ú±íÓÃbilateral weights
options.our_method = double(bicubic(LRImg));
our_HRImg = RMLS_interpolation(LRImg,options);

diff_img = our_HRImg-original_HRImg;
mse=sum(diff_img(:).*diff_img(:))/prod(size(original_HRImg));
our_PSNR=10*log10(255^2/mse);
ssim = ssim_index(original_HRImg,our_HRImg);

figure (3+(im_num-1)*5);imshow(our_HRImg/255);
figure (4+(im_num-1)*5);imshow(diff_img/255);
figure (1);imshow(our_HRImg/255);
title(sprintf('RLLR HR image, PSNR = %2.2f dB, SSIM =%2.2f  ', our_PSNR,ssim));
imwrite(uint8(our_HRImg), 'E:\xmliu\result\RLLR_butterfly1.bmp', 'bmp' );

%%%%%%%
options.weight_fg = 0;
options.our_method = double(bicubic(LRImg));
our_HRImg = RMLS_interpolation(LRImg,options);

diff_img = our_HRImg-original_HRImg;
mse=sum(diff_img(:).*diff_img(:))/prod(size(original_HRImg));
our_PSNR=10*log10(255^2/mse);

% figure (2+(im_num-1)*5);imshow(our_HRImg/255);
figure (2);imshow(our_HRImg/255);
title(sprintf('RLLR HR image, PSNR = %2.2f dB', our_PSNR));
imwrite(uint8(our_HRImg), 'E:\xmliu\result\RLLR_butterfly2.bmp', 'bmp' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% the method of lei zhang
esintp_HRImg = esintp(LRImg);
diff_img = double(esintp_HRImg)-original_HRImg;
mse=sum(diff_img(:).*diff_img(:))/prod(size(original_HRImg));
PSNR=10*log10(255^2/mse);

figure (3+(im_num-1)*5);imshow(esintp_HRImg/255);
title(sprintf('esintp HR image, PSNR = %2.2f dB', PSNR));
imwrite(uint8(esintp_HRImg), 'E:\xmliu\result\FDI_butterfly.bmp', 'bmp' );
%% bilinear
bilinear_HRImg = bilinearup2(LRImg);
diff_img = bilinear_HRImg-original_HRImg;

mse=sum(diff_img(:).*diff_img(:))/prod(size(original_HRImg));
bilinear_PSNR=10*log10(255^2/mse);

figure (4+(im_num-1)*5);imshow(bilinear_HRImg/255);
title(sprintf('bilinear HR image, PSNR = %2.2f dB', bilinear_PSNR));
imwrite(uint8(bilinear_HRImg), 'E:\xmliu\result\bilinear_butterfly.bmp', 'bmp' );




%% bicubic
bicubic_HRImg = double(bicubic(LRImg));
diff_img = bicubic_HRImg-original_HRImg;
mse=sum(diff_img(:).*diff_img(:))/prod(size(original_HRImg));
bicubic_PSNR=10*log10(255^2/mse);

figure (5+(im_num-1)*5);imshow(bicubic_HRImg/255);
title(sprintf('bicubic HR image, PSNR = %2.2f dB', bicubic_PSNR));
imwrite(uint8(bicubicr_HRImg), 'E:\xmliu\result\bicubic_butterfly.bmp', 'bmp' );
% 
% Xiaolin Wu

SAI_HRImg = double(imread('SAI_splash.pgm')); 
original_HRImg_PGM = double(imread('HR_splash.pgm')); 
% 
diff_img = SAI_HRImg-original_HRImg_PGM;
mse=sum(diff_img(:).*diff_img(:))/prod(size(original_HRImg_PGM));
SAI_PSNR=10*log10(255^2/mse);
figure (6);imshow(SAI_HRImg/255);
title(sprintf('SAI HR image, PSNR = %2.2f dB', SAI_PSNR));

end